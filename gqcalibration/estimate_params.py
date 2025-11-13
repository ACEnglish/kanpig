#!/usr/bin/env python3
"""
Estimate beta-binomial mixture model parameters from genotyping truth set.
Fits parameters for three genotype classes: 0/0, 0/1, 1/1
"""

import pandas as pd
import numpy as np
from scipy.optimize import minimize
from scipy.special import betaln, logsumexp
import kanpig
import json

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

def fit_parameters(df, verbose=True):
    """
    Fit beta-binomial mixture parameters from truth set.
    
    Args:
        df: DataFrame with columns 'Ogt' (0,1,2), 'AD_alt', 'DP'
        verbose: print progress
    
    Returns:
        dict with fitted parameters
    """
    # Filter out low depth sites
    df = df[(df['DP'] >= 5) & df['state']].copy()
    
    df['Ogt'] = df['Ogt'].map({'REF':0,'HET':1, 'HOM':2})
    if verbose:
        print(f"Using {len(df)} sites for parameter estimation")
        print(f"Genotype distribution: {df['Ogt'].value_counts().to_dict()}")
    
    # Get initial parameter estimates
    initial = estimate_initial_params(df)
    
    if verbose:
        print("\nInitial parameter estimates:")
        for gt, params in initial.items():
            print(f"  GT {gt}: mu={params['mu']:.3f}, nu={params['nu']:.1f}, n={params['n']}")
    
    # Prepare data for optimization
    data = {
        'gt': df['Ogt'].values,
        'ad_alt': df['AD_alt'].values,
        'dp': df['DP'].values
    }
    
    # Calculate empirical genotype fractions
    gt_counts = df['Ogt'].value_counts()
    total = len(df)
    frac = [gt_counts.get(i, 0) / total for i in [0, 1, 2]]
    
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
    
    if verbose:
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
    
    if verbose:
        print("\nFitted parameters:")
        print(f"  Mixture fractions: [{frac0:.3f}, {frac1:.3f}, {frac2:.3f}]")
        print(f"  Means (mu):        [{mu0:.4f}, {mu1:.4f}, {mu2:.4f}]")
        print(f"  Precisions (nu):   [{nu0:.1f}, {nu1:.1f}, {nu2:.1f}]")
        print(f"  Log-likelihood:    {fitted['log_likelihood']:.1f}")
    
    return fitted

def format_rust_output(fitted):
    """Format parameters for Rust code."""
    frac = fitted['frac']
    mu = fitted['mu']
    nu = fitted['nu']
    
    rust_code = f"""
// Fitted beta-binomial mixture parameters
let frac: &[f64] = &[{frac[0]:.4f}, {frac[1]:.4f}, {frac[2]:.4f}];
let mu = &[{mu[0]:.6f}, {mu[1]:.6f}, {mu[2]:.6f}];
let nu = &[{nu[0]:.2f}, {nu[1]:.2f}, {nu[2]:.2f}];
"""
    return rust_code

def save_config(fitted, filename, calibration=None):
    """Save parameters to JSON config file."""

    config = {
        'mixture_fractions': fitted['frac'],
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
    
    print(f"\nSaved parameters to {filename}")

def newgt(gt, row):
    """
    Rerun the genotyper on a row
    """
    return gt.genotype(row['AD_ref'], row['AD_alt']).gq

from scipy.interpolate import interp1d
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

# Example usage
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python estimate_params.py <truth_set.csv> [output.json]")
        sys.exit(1)
    
    # Load truth set
    df = pd.read_csv(sys.argv[1])
    
    # Fit parameters
    fitted = fit_parameters(df)
    
    print("Bestfit parameters:")
    print("="*60)
    print(format_rust_output(fitted))
    # Now you need to go re-genotype everything and grab those GQs
    # Then you make the calibration table
    # 
    # Save config if requested
    if len(sys.argv) >= 3:
        save_config(fitted, sys.argv[2])

    genotyper = kanpig.Genotyper(sys.argv[2])

    df['nGQ'] = df.apply((lambda x: newgt(genotyper, x)), axis=1)
    # Make calibrate
    calibration = build_calibration(df)
    if len(sys.argv) >= 3:
        save_config(fitted, sys.argv[2], calibration)
