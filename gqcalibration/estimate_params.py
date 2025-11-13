#!/usr/bin/env python3
"""
Estimate beta-binomial mixture model parameters from genotyping truth set.
Fits parameters for three genotype classes: 0/0, 0/1, 1/1
"""
import sys
import json
import argparse
import kanpig
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from scipy.special import betaln, logsumexp

def beta_binomial_logpmf(k, n, mu, nu):
    """
    Log probability mass function for beta-binomial distribution.
    
    Parameters:
        k: number of successes (alt allele count)
        n: number of trials (total depth)
        mu: mean (alt allele fraction)
        nu: precision parameter (higher = less overdispersion)
    
    Returns beta-binomial with alpha = mu*nu, beta = (1-mu)*nu
    """
    alpha = mu * nu
    beta = (1 - mu) * nu
    
    return (betaln(k + alpha, n - k + beta) - 
            betaln(alpha, beta))

def mixture_loglikelihood(params, data):
    """
    Negative log-likelihood for beta-binomial mixture model.
    
    params: [frac0, frac1, mu0, mu1, mu2, nu0, nu1, nu2]
    data: dict with 'gt' (genotype 0,1,2), 'ad_alt', 'dp'
    """
    frac0, frac1 = params[0:2]
    frac2 = 1 - frac0 - frac1
    
    mu0, mu1, mu2 = params[2:5]
    nu0, nu1, nu2 = params[5:8]
    
    # Assign mixture component based on true genotype
    gt = data['gt']
    ad_alt = data['ad_alt']
    dp = data['dp']
    
    log_likes = []
    for i in range(len(gt)):
        if gt[i] == 0:
            ll = beta_binomial_logpmf(ad_alt[i], dp[i], mu0, nu0)
        elif gt[i] == 1:
            ll = beta_binomial_logpmf(ad_alt[i], dp[i], mu1, nu1)
        else:  # gt[i] == 2
            ll = beta_binomial_logpmf(ad_alt[i], dp[i], mu2, nu2)
        log_likes.append(ll)
    
    return -np.sum(log_likes)

def estimate_initial_params(df):
    """Estimate initial parameters from observed allele fractions per genotype."""
    initial = {}
    
    for gt in [0, 1, 2]:
        subset = df[df['Ogt'] == gt]
        if len(subset) == 0:
            continue
            
        # Calculate allele fraction
        af = subset['AD_alt'] / subset['DP']
        
        # Estimate mean
        mu = af.mean()
        
        # Estimate nu from variance (method of moments)
        # For beta-binomial: Var(p) = mu(1-mu)/(1+nu)
        var_p = af.var()
        if var_p > 0 and mu > 0 and mu < 1:
            nu = mu * (1 - mu) / var_p - 1
            nu = max(1.0, nu)  # ensure nu >= 1
        else:
            nu = 10.0
        
        initial[gt] = {'mu': mu, 'nu': nu, 'n': len(subset)}
    
    return initial

def fit_parameters(df, all_gts=False, min_dp=5, max_dp=60, all_hets=False):
    """
    Fit beta-binomial mixture parameters from truth set.
    """
    df = df[(df['DP'] >= min_dp) & (df['DP'] <= max_dp)].copy()
   
    df['Ogt'] = df['Ogt'].map({'REF':0,'HET':1, 'HOM':2})

    if all_hets:
        het_data = df[df['Ogt'] == 1]  # All hets
        ref_data = df[(df['Ogt'] == 0) & df['state']]  # Only correct refs
        hom_data = df[(df['Ogt'] == 2) & df['state']]  # Only correct homs
        df = pd.concat([het_data, ref_data, hom_data])
    elif all_gts:
        df = df[df['state']].copy()
 

    print(f"Using {len(df)} sites for parameter estimation")
    print(f"Genotype distribution: {df['Ogt'].value_counts().to_dict()}")
    
    # Get initial parameter estimates
    initial = estimate_initial_params(df)
    
    print("\nInitial parameter estimates:")
    for gt, params in initial.items():
        print(f"  GT {gt}: mu={params['mu']:.3f}, nu={params['nu']:.1f}, n={params['n']}")

    # Calculate empirical genotype fractions
    gt_counts = df['Ogt'].value_counts()
    total = len(df)
    frac = [gt_counts.get(i, 0) / total for i in [0, 1, 2]]
    
    # Prepare data for optimization
    data = {
        'gt': df['Ogt'].values,
        'ad_alt': df['AD_alt'].values,
        'dp': df['DP'].values
    }

    # Initial parameter vector
    x0 = [
        frac[0], frac[1],  # mixture fractions (frac2 = 1 - frac0 - frac1)
        initial.get(0, {'mu': 0.001})['mu'],
        initial.get(1, {'mu': 0.50})['mu'],
        initial.get(2, {'mu': 0.97})['mu'],
        initial.get(0, {'nu': 25.0})['nu'],
        initial.get(1, {'nu': 15.0})['nu'],
        initial.get(2, {'nu': 10.0})['nu']
    ]
    
    # Bounds: fractions in [0,1], mu in [0,1], nu > 0
    bounds = [
        (0.01, 0.98), (0.01, 0.98),  # frac0, frac1
        (0.0001, 0.05), (0.3, 0.7), (0.95, 0.9999),  # mu0, mu1, mu2
        (1.0, 100.0), (1.0, 100.0), (1.0, 100.0)  # nu0, nu1, nu2
    ]
    
    # Constraint: frac0 + frac1 <= 0.99
    constraints = {'type': 'ineq', 'fun': lambda x: 0.99 - x[0] - x[1]}
    
    print("\nOptimizing parameters...")
    
    result = minimize(
        mixture_loglikelihood,
        x0,
        args=(data,),
        method='SLSQP',
        bounds=bounds,
        constraints=constraints,
        options={'maxiter': 500}
    )
    
    if not result.success:
        print(f"Warning: optimization did not fully converge: {result.message}")
    
    # Extract fitted parameters
    frac0, frac1 = result.x[0:2]
    frac2 = 1 - frac0 - frac1
    mu0, mu1, mu2 = result.x[2:5]
    nu0, nu1, nu2 = result.x[5:8]
    
    fitted = {
        'frac': [frac0, frac1, frac2],
        'mu': [mu0, mu1, mu2],
        'nu': [nu0, nu1, nu2],
        'log_likelihood': -result.fun,
        'n_sites': len(df)
    }
    
    print("\nFitted parameters:")
    print(f"  Mixture fractions: [{frac0:.3f}, {frac1:.3f}, {frac2:.3f}]")
    print(f"  Means (mu):        [{mu0:.4f}, {mu1:.4f}, {mu2:.4f}]")
    print(f"  Precisions (nu):   [{nu0:.1f}, {nu1:.1f}, {nu2:.1f}]")
    print(f"  Log-likelihood:    {fitted['log_likelihood']:.1f}")
    
    return fitted


def save_config(fitted, filename, calibration=None, flat_priors=False):
    """Save parameters to JSON config file."""

    config = {
        'mixture_fractions': [0.33, 0.34, 0.33] if flat_priors else fitted['frac'],
        'means': fitted['mu'],
        'precisions': fitted['nu'],
        "calibration_table": [] if calibration is None else calibration,
        'metadata': {
            'n_sites': fitted['n_sites'],
            'log_likelihood': fitted['log_likelihood']
        }
    }
    
    with open(filename, 'w') as f:
        json.dump(config, f, indent=2)

    return config

def build_calibration(df):
    gqs = df['nGQ']
    correct = df['state']
    # Bin by GQ and calculate observed accuracy
    gq_bins = np.arange(0, 101, 5)
    observed_accuracies = []
    bin_midpoints = []

    for i in range(len(gq_bins) - 1):
        gq_min, gq_max = gq_bins[i], gq_bins[i+1]
        mask = (gqs >= gq_min) & (gqs < gq_max)

        if mask.sum() > 10:  # Need sufficient data
            accuracy = correct[mask].mean()
            observed_accuracies.append(accuracy)
            bin_midpoints.append((gq_min + gq_max) / 2)

    # Convert observed accuracy to calibrated GQ
    # P(correct) = accuracy, so P(error) = 1 - accuracy
    # GQ = -10 * log10(P(error))
    calibrated_gqs = [-10 * np.log10(max(1 - acc, 1e-10)) for acc in observed_accuracies]

    # Create interpolation function
    calibration_func = interp1d(bin_midpoints, calibrated_gqs,
                               kind='linear',
                               bounds_error=False,
                               fill_value=(calibrated_gqs[0], calibrated_gqs[-1]))
    return [[float(x), float(y)] for x, y in zip(bin_midpoints, calibrated_gqs)]

def parse_args(args):
    parser = argparse.ArgumentParser(prog="bench", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("IN", type=str,
                        help="Input CSV from `make_params_df.py`")
    parser.add_argument("OUT", type=str,
                        help="Output prefix file to write")
    parser.add_argument("--all", action="store_true",
                        help="Fit on all genotypes, not just correct")
    parser.add_argument("--min_dp", type=int, default=5,
                        help="Minimum GT depth to fit on (%(default)s)")
    parser.add_argument("--max_dp", type=int, default=60,
                        help="MaximumGT depth to fit on (%(default)s)")
    parser.add_argument("--all-hets", action="store_true",
                        help="Fit on all hets, but only correct REF/HOM")
    parser.add_argument("--flat-priors", action="store_true",
                        help="Use flat priors (0.33) instead of observed GT states")
    parser.add_argument("--no-calibrate", action="store_true",
                        help="Don't perform GQ calibration")
    parser.add_argument("--write-calib", action="store_true",
                        help="Write a csv of the calibrated GQs")
    parser.add_argument("--no-plots", action="store_true",
                        help="Skip plotting")
    args = parser.parse_args(args)
    if args.all and args.all-hets:
        print("Error! Can only fit either --all-hets XOR --all")
    return args

def make_plots(data, out_prefix):
    # Optional requirements
    import seaborn as sb
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score

    
    print("Making GT<->GQ Plot")
    # Sort data by GQ
    df_sorted = data.sort_values(by='GQ').reset_index(drop=True)
    roll = 1000
    # Calculate rolling averages
    state_rolling = df_sorted['state'].rolling(roll).mean()
    gq_rolling = (df_sorted['GQ']).rolling(roll).mean()
    gq_rolling = 1 - 10**(-gq_rolling / 10)
    # Create figure with twin y-axes
    fig, ax1 = plt.subplots(figsize=(8, 6), dpi=180)

    # First y-axis: State (accuracy)
    color1 = 'tab:blue'
    ax1.set_xlabel(f'Variants (sorted by GQ)', fontsize=12)
    ax1.set_ylabel('Accuracy (rolling avg)', color=color1, fontsize=12)
    ax1.plot(state_rolling, color=color1, linewidth=2, label='Observed')
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.set_ylim(0, 1)

    # Second y-axis: GQ
    ax2 = ax1#ax1.twinx()
    color2 = 'tab:orange'
    ax2.set_ylabel('Accuracy (rolling avg)', fontsize=12)
    ax2.plot(gq_rolling, color=color2, linewidth=2, label='GQ')
    ax2.tick_params(axis='y', labelcolor=color2)

    # Title and grid
    plt.title(f'Genotype Accuracy and Quality (200-sample rolling average)',
              fontsize=14, pad=20)
    ax1.grid(True, alpha=0.3)

    # Add legends
    ax1.legend(loc='lower right')

    plt.tight_layout()
    plt.savefig(out_prefix + '.GTGQ.png')

    print("Making ROC curves")
    fig, ax1 = plt.subplots(1, 1, figsize=(8, 6), dpi=180)
    # Colors for different scores
    colors = plt.cm.Set1(np.linspace(0, 1, 1))

    y_true = data['state'].astype(int)
    y_score = data['GQ']
    
    # ROC curve
    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)
    
    # Precision-Recall curve
    precision, recall, _ = precision_recall_curve(y_true, y_score)
    avg_precision = average_precision_score(y_true, y_score)
    
    # Plot ROC
    ax1.plot(fpr, tpr, color=colors[0], lw=2, 
             label=f'GQ (AUC = {roc_auc:.3f})')

    # Format ROC plot
    ax1.plot([0, 1], [0, 1], 'k--', lw=1, label='Random (AUC = 0.5)')
    ax1.set_xlim([0.0, 1.0])
    ax1.set_ylim([0.0, 1.05])
    ax1.set_xlabel('False Positive Rate', fontsize=12)
    ax1.set_ylabel('True Positive Rate', fontsize=12)
    ax1.set_title('Kanpig ROC Curve', fontsize=14, fontweight='bold')
    ax1.legend(loc="lower right")
    ax1.grid(alpha=0.3)

    plt.savefig(out_prefix + "ROC.png")

    print("Making Genotypes Plot")
    fig, ax = plt.subplots(3, 4, figsize=(12, 6), dpi=180)
    xlim=(0, data['GQ'].max() + 1)
    for i, m_ax in zip(['REF', 'HET', 'HOM'], ax):
        p = sb.histplot(data=data[data['Ogt'] == i],
                        x='GQ', hue='state', multiple='stack',
                        binwidth=1, ax=m_ax[0])
        p.set(title=f"Baseline", ylabel=i + ' Count', xlim=xlim)

        p = sb.histplot(data=data[data['Mgt'] == i],
                        x='GQ', hue='state', multiple='stack',
                        binwidth=1, ax=m_ax[1])
        p.set(title=f"Kanpig", ylabel= i + ' Count', xlim=xlim)

        subset = data[data['Ogt'] == i]
        af = subset['AD_alt'] / subset['DP']
        p = sb.histplot(af[subset['state']], bins=50, ax=m_ax[2], binwidth=0.02)
        p.set(xlabel='Allele Fraction',
              yscale='log',
              ylabel=i + ' Count (log)',
              xlim=(0,1),
              title='True GT')

        p = sb.histplot(af[~subset['state']], bins=50, ax=m_ax[3], binwidth=0.02)
        p.set(xlabel='Allele Fraction',
              yscale='log',
              xlim=(0,1),
              ylabel=i + ' Count (log)',
              title='False GT')

    plt.tight_layout()
    plt.savefig(out_prefix + '.STATE.png')


# Example usage
if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    
    # Load truth set
    df = pd.read_csv(args.IN)
    
    cnt = df.groupby(['Ogt', 'state']).size().unstack()
    cnt.loc['All'] = cnt.sum(axis=0)
    cnt['Acc'] = cnt[True] / cnt.sum(axis=1)
    print("Genotype Accuracy")
    print(cnt)
    print()

    print("GT Confusion Matrix")
    print(df.groupby(["Ogt", "Mgt"]).size().unstack())
    print()

    # Fit parameters
    fitted = fit_parameters(df,
                            all_gts=args.all,
                            min_dp=args.min_dp,
                            max_dp=args.max_dp,
                            all_hets=args.all_hets,
                            )
    
    # Now you need to go re-genotype everything and grab those GQs
    # Then you make the calibration table
    out_cfg = args.OUT + '.json'
    config = save_config(fitted, out_cfg, flat_priors=args.flat_priors)

    if not args.no_calibrate:
        print("Calibrating GQs")
        gt = kanpig.Genotyper(out_cfg)
        df['nGQ'] = df.apply((lambda x: gt.genotype(x['AD_ref'], x['AD_alt']).gq), axis=1)
        calibration = build_calibration(df)
        config = save_config(fitted, out_cfg, calibration, flat_priors=args.flat_priors)
        # And then we have to run again to actually get the calibrated GQs
        if args.write_calib or not args.no_plots:
            print("Regenotyping with config")
            gt = kanpig.Genotyper(out_cfg)
            df['nGQ'] = df.apply((lambda x: gt.genotype(x['AD_ref'], x['AD_alt']).gq), axis=1)
            df.to_csv(args.OUT + '.calibrated.csv')
    
    if not args.no_plots:
        print("Making original plots")
        df = df[(df['DP'] >= args.min_dp) & (df['DP'] <= args.max_dp)].copy()
        make_plots(df, args.OUT + '.original')
        if not args.no_calibrate:
            print("Making calibrated plots")
            df['GQ'] = df['nGQ']
            make_plots(df, args.OUT + '.calibrated')
    # Optional plotting here before/after GQs' ROC curves, accuracy curve, distribution plot, calibration plot
    print("Finished")
