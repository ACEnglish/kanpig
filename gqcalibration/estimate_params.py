#!/usr/bin/env python3
"""
Estimate beta-binomial mixture model parameters from genotyping truth set.
Fits parameters for three genotype classes: 0/0, 0/1, 1/1
"""
import os
import sys
import json
import logging
import argparse
import textwrap
from functools import partial
 
# pylint: disable=redefined-outer-name,no-name-in-module,no-member

import kanpig
import truvari
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from scipy.special import betaln
from scipy.optimize import minimize
from matplotlib.patches import Rectangle
from sklearn.metrics import roc_curve, auc

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
    """
    mu0, mu1, mu2 = params[2:5]
    nu0, nu1, nu2 = params[5:8]

    # Assign mixture component based on true genotype
    gt = data['gt']
    ad_alt = data['ad_alt']
    dp = data['dp']

    log_likes = []
    for i , _ in enumerate(gt):
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
        if var_p > 0 and 0 < mu < 1:
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

    df['Ogt'] = df['Ogt'].map({'REF': 0, 'HET': 1, 'HOM': 2})

    if all_hets:
        het_data = df[df['Ogt'] == 1]  # All hets
        ref_data = df[(df['Ogt'] == 0) & df['state']]  # Only correct refs
        hom_data = df[(df['Ogt'] == 2) & df['state']]  # Only correct homs
        df = pd.concat([het_data, ref_data, hom_data])
    elif all_gts:
        df = df[df['state']].copy()

    logging.info(f"Using {len(df)} sites for parameter estimation")
    logging.info(f"Genotype distribution:")
    logging.info(f"  {df['Ogt'].value_counts().to_dict()}")

    # Get initial parameter estimates
    initial = estimate_initial_params(df)

    logging.info("Initial parameter estimates:")
    for gt, params in initial.items():
        logging.info(
            f"  GT {gt}: mu={params['mu']:.3f}, nu={params['nu']:.1f}, n={params['n']}")

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

    logging.info("Optimizing parameters...")

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
        logging.warning(f"Optimization did not fully converge: {result.message}")

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

    logging.info("Fitted parameters:")
    logging.info(f"  Mixture fractions: [{frac0:.3f}, {frac1:.3f}, {frac2:.3f}]")
    logging.info(f"  Means (mu):        [{mu0:.4f}, {mu1:.4f}, {mu2:.4f}]")
    logging.info(f"  Precisions (nu):   [{nu0:.1f}, {nu1:.1f}, {nu2:.1f}]")
    logging.info(f"  Log-likelihood:    {fitted['log_likelihood']:.1f}")

    return fitted


def save_config(fitted, filename, calibration=None, flat_priors=False):
    """Save parameters to JSON config file."""

    config = {
        'mode': 'Beta',
        'mixture_fractions': [0.33, 0.34, 0.33] if flat_priors else fitted['frac'],
        'means': fitted['mu'],
        'precisions': fitted['nu'],
        "calibration_table": [] if calibration is None else calibration,
        'metadata': {
            'n_sites': fitted['n_sites'],
            'log_likelihood': fitted['log_likelihood'],
            'estimation_params': " ".join(sys.argv[1:])
        }
    }

    with open(filename, 'w') as f:
        json.dump(config, f, indent=2)

    return config


def build_calibration(df):
    """
    GQ Calibration
    """
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
    calibrated_gqs = [-10 * np.log10(max(1 - acc, 1e-10))
                      for acc in observed_accuracies]

    # Create interpolation function
    # I could be using this to create nGQ instead of calling kanpig again
    # But maybe better to actually test the kanpig code
    #calibration_func = interp1d(bin_midpoints, calibrated_gqs,
    #                            kind='linear',
    #                            bounds_error=False,
    #                            fill_value=(calibrated_gqs[0], calibrated_gqs[-1]))
    return [[float(x), float(y)] for x, y in zip(bin_midpoints, calibrated_gqs)]


def parse_args(args):
    """
    Parse arguments
    """
    parser = argparse.ArgumentParser(prog="bench", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("IN", type=str,
                        help="Input VCF file or `.genotypes.csv` from previous run")
    parser.add_argument("OUT", type=str,
                        help="Output prefix file to write")
    parser.add_argument("--bed", default=None, type=str,
                        help="Bed file for subsetting VCF entries to parse")
    parser.add_argument("--summary", action='store_true',
                        help="Only make summary report of genotypes")
    parser.add_argument("--leaveout", type=float, default=0,
                        help=("Leave out [0.0-1.0) from training and "
                              "make separate plots (%(default)s)"))
    parser.add_argument("--mindp", type=int, default=5,
                        help="Minimum GT depth to fit on (%(default)s)")
    parser.add_argument("--maxdp", type=int, default=60,
                        help="MaximumGT depth to fit on (%(default)s)")
    parser.add_argument("--sizemin", type=int, default=50,
                        help="Minimum SV size to parse from VCF (%(default)s)")
    parser.add_argument("--sizemax", type=int, default=10000,
                        help="Maximum SV size to parse from VCF (%(default)s)")
    parser.add_argument("--all", action="store_true",
                        help="Fit on all genotypes, not just correct ones")
    parser.add_argument("--all-hets", action="store_true",
                        help="Fit on all hets, but only correct REF/HOM")
    parser.add_argument("--flat-priors", action="store_true",
                        help="Use flat priors (0.33) instead of observed GT states")
    parser.add_argument("--no-calibrate", action="store_true",
                        help="Don't perform GQ calibration")
    parser.add_argument("--write-calib", action="store_true",
                        help="Write a csv of the calibrated GQs")
    args = parser.parse_args(args)
    return args


class HTMLReportBuilder:
    """Builder for creating HTML reports with multiple plot sections."""
    
    def __init__(self, title="QC Report", output_dir="qc_plots"):
        self.title = title
        self.output_dir = output_dir
        self.sections = []
        self.plot_counter = 0
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
    
    def add_section(self, section_title, plots):
        """
        Add a section with plots to the report.
        
        Args:
            section_title: Title for this section
            plots: List of tuples (plot_title, figure_object)
        """
        plot_files = []
        
        for plot_title, fig in plots:
            # Generate filename
            self.plot_counter += 1
            filename = f"plot_{self.plot_counter:03d}.png"
            filepath = os.path.join(self.output_dir, filename)
            
            # Save the figure
            fig.savefig(filepath, bbox_inches='tight')
            
            # Store relative path for HTML
            plot_files.append((plot_title, filename))
            
            # Close figure to free memory
            plt.close(fig)
        
        self.sections.append({
            'title': section_title,
            'type': 'plots',
            'plots': plot_files
        })
    
    def add_table_section(self, section_title, table1_title, table1_df, table2_title, table2_df):
        """
        Add a section with two tables displayed side-by-side.
        
        Args:
            section_title: Title for this section
            table1_title: Title for the first table
            table1_df: First pandas DataFrame
            table2_title: Title for the second table
            table2_df: Second pandas DataFrame
        """
        self.sections.append({
            'title': section_title,
            'type': 'tables',
            'table1': {'title': table1_title, 'df': table1_df},
            'table2': {'title': table2_title, 'df': table2_df}
        })
    
    def save(self, output_path):
        """Generate and save the HTML report."""
        # Get relative path from HTML to plots directory
        html_dir = os.path.dirname(os.path.abspath(output_path))
        plots_dir = os.path.abspath(self.output_dir)
        rel_path = os.path.relpath(plots_dir, html_dir)
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{self.title}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
            background-color: #F0F0F0;
        }}
        h1 {{
            color: #333;
            border-bottom: 3px solid #E8A5B5;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #555;
            margin-top: 40px;
            border-bottom: 2px solid #ddd;
            padding-bottom: 8px;
        }}
        .section {{
            background: white;
            padding: 20px;
            margin: 20px 0;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .plot-container {{
            margin: 30px 0;
        }}
        .plot-title {{
            font-size: 18px;
            font-weight: bold;
            color: #444;
            margin-bottom: 15px;
        }}
        .plot-image {{
            width: 100%;
            height: auto;
            border: 1px solid #ddd;
            border-radius: 4px;
            cursor: pointer;
        }}
        .plot-image:hover {{
            opacity: 0.9;
            box-shadow: 0 4px 8px rgba(0,0,0,0.2);
        }}
        .plot-link {{
            display: inline-block;
            margin-top: 5px;
            color: #E8A5B5;
            text-decoration: none;
            font-size: 14px;
        }}
        .plot-link:hover {{
            text-decoration: underline;
        }}
        .timestamp {{
            color: #888;
            font-size: 14px;
            text-align: right;
            margin-top: 20px;
        }}
        .info-box {{
            background: #FFE5EC;
            border-left: 4px solid #E8A5B5;
            padding: 15px;
            margin: 20px 0;
            border-radius: 4px;
        }}
        .table-container {{
            display: flex;
            gap: 20px;
            margin: 20px 0;
            flex-wrap: wrap;
        }}
        .table-wrapper {{
            flex: 1;
            min-width: 300px;
        }}
        .table-title {{
            font-size: 16px;
            font-weight: bold;
            color: #444;
            margin-bottom: 10px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 14px;
            background: white;
        }}
        th, td {{
            padding: 8px 12px;
            text-align: left;
            border: 1px solid #ddd;
        }}
        th {{
            background-color: #E8A5B5;
            color: black;
            font-weight: bold;
        }}
        tr:nth-child(even) {{
            background-color: #f9f9f9;
        }}
        tr:hover {{
            background-color: #f0f0f0;
        }}
        td:first-child {{
            font-weight: 500;
        }}
    </style>
</head>
<body>
    <h1>{self.title}</h1>
    <div class="info-box">
        <strong>Note:</strong> All plots are available as PNG files in the <code>{self.output_dir}/</code> directory.
    </div>
"""
        
        for section in self.sections:
            html_content += f"""
    <div class="section">
        <h2>{section['title']}</h2>
"""
            if section['type'] == 'plots':
                for plot_title, filename in section['plots']:
                    plot_path = os.path.join(rel_path, filename).replace('\\', '/')
                    html_content += f"""
        <div class="plot-container">
            <div class="plot-title">{plot_title}</div>
            <a href="{plot_path}" target="_blank">
                <img class="plot-image" src="{plot_path}" alt="{plot_title}">
            </a>
            <br>
            <a class="plot-link" href="{plot_path}" download>📥 Download {filename}</a>
        </div>
"""
            elif section['type'] == 'tables':
                html_content += """
        <div class="table-container">
"""
                for table_info in [section['table1'], section['table2']]:
                    html_content += f"""
            <div class="table-wrapper">
                <div class="table-title">{table_info['title']}</div>
                {table_info['df'].to_html(classes='', border=0, escape=False)}
            </div>
"""
                html_content += """
        </div>
"""
            html_content += """
    </div>
"""
        
        html_content += """
    <div class="timestamp">
        Report generated: <span id="timestamp"></span>
    </div>
    <script>
        document.getElementById('timestamp').textContent = new Date().toLocaleString();
    </script>
</body>
</html>
"""
        
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)

def make_plots(data, section_title, html_builder):
    """
    Generate QC plots and add them to the HTML report.
    
    Args:
        data: DataFrame with genotype data
        section_title: Title for this data section in the report
        html_builder: HTMLReportBuilder instance to add plots to
    """
   
    plots = []
    
    # GTGQ Plot
    df_sorted = data.sort_values(by='GQ').reset_index(drop=True)
    roll = max(len(df_sorted) // 50, 2)
    state_rolling = df_sorted['state'].rolling(roll).mean()
    gq_rolling = (df_sorted['GQ']).rolling(roll).mean()
    gq_rolling = 1 - 10**(-gq_rolling / 10)
    
    fig1, ax1 = plt.subplots(figsize=(8, 6), dpi=180)
    color1 = 'tab:blue'
    ax1.set_xlabel('Variants (sorted by GQ)', fontsize=12)
    ax1.set_ylabel('Accuracy (rolling avg)', fontsize=12)
    ax1.plot(state_rolling, color=color1, linewidth=2, label='Observed')
    ax1.tick_params(axis='y')
    ax1.set_ylim(0, 1)
    
    ax2 = ax1
    color2 = 'tab:orange'
    ax2.set_ylabel('Accuracy (rolling avg)', fontsize=12)
    ax2.plot(gq_rolling, color=color2, linewidth=2, label='GQ')
    ax2.tick_params(axis='y')
    plt.title(f'Genotype Accuracy and Quality ({roll}-sample rolling average)',
              fontsize=14, pad=20)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='lower right')
    plt.tight_layout()
    
    plots.append(("Genotype Accuracy vs Quality", fig1))
    
    # ROC Plot
    fig2, ax1 = plt.subplots(1, 1, figsize=(8, 6), dpi=180)
    colors = plt.cm.Set1(np.linspace(0, 1, 1))
    y_true = data['state'].astype(int)
    y_score = data['GQ']
    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)
    ax1.plot(fpr, tpr, color=colors[0], lw=2,
             label=f'GQ (AUC = {roc_auc:.3f})')
    ax1.plot([0, 1], [0, 1], 'k--', lw=1, label='Random (AUC = 0.5)')
    ax1.set_xlim([0.0, 1.0])
    ax1.set_ylim([0.0, 1.05])
    ax1.set_xlabel('False Positive Rate', fontsize=12)
    ax1.set_ylabel('True Positive Rate', fontsize=12)
    ax1.set_title('Kanpig ROC Curve', fontsize=14, fontweight='bold')
    ax1.legend(loc="lower right")
    ax1.grid(alpha=0.3)
    
    plots.append(("ROC Curve", fig2))
    
    # STATE Plot
    fig3, ax = plt.subplots(3, 4, figsize=(12, 6), dpi=180)
    xlim = (0, data['GQ'].max() + 1)
    for i, m_ax in zip(['REF', 'HET', 'HOM'], ax):
        if (data['Ogt'] == i).sum() == 0:
            logging.warning(f"No Ogt == {i} sites found. Skipping")
        else:
            p = sb.histplot(data=data[data['Ogt'] == i],
                            x='GQ', hue='state', multiple='stack',
                            hue_order=[False, True],
                            palette=sb.color_palette()[:2],
                            binwidth=1, ax=m_ax[0])
            p.set(title="Baseline", ylabel=i + ' Count', xlim=xlim)

            subset = data[data['Ogt'] == i]
            af = subset['AD_alt'] / subset['DP']
            p = sb.histplot(af[subset['state']],
                            color=sb.color_palette()[1],
                            ax=m_ax[2], binwidth=0.02)
            p.set(xlabel='Allele Fraction',
                  yscale='log',
                  ylabel=i + ' Count (log)',
                  xlim=(0, 1),
                  title='True GT')
            p = sb.histplot(af[~subset['state']],
                            color=sb.color_palette()[0],
                            ax=m_ax[3], binwidth=0.02)
            p.set(xlabel='Allele Fraction',
                  yscale='log',
                  xlim=(0, 1),
                  ylabel=i + ' Count (log)',
                  title='False GT')
        if (data['Kgt'] == i).sum() == 0:
            logging.warning(f"No Kgt == {i} sites found. Skipping")
        else:
            p = sb.histplot(data=data[data['Kgt'] == i],
                            x='GQ', hue='state', multiple='stack',
                            binwidth=1, ax=m_ax[1])
            p.set(title="Kanpig", ylabel=i + ' Count', xlim=xlim)

    plt.tight_layout()

    # Customize legend for the stacked histograms
    for i, m_ax in enumerate(ax):
        for j in [0, 1]:  # First two columns have the 'state' hue
            legend = m_ax[j].get_legend()
            if legend is not None:
                colors = sb.color_palette()[:2]
                legend.remove()
                
                # Get current title and add legend inline
                title = m_ax[j].get_title()
                m_ax[j].text(0.5, 1.08, title, transform=m_ax[j].transAxes,
                            fontsize=10, ha='center', va='center', fontweight='normal')
                
                # Add colored squares to the right of title
                m_ax[j].text(0.82, 1.08, 'F', transform=m_ax[j].transAxes,
                            bbox=dict(boxstyle='square,pad=0.2', facecolor=colors[0], 
                                    edgecolor='black', linewidth=0.5),
                            fontsize=7, ha='center', va='center')
                m_ax[j].text(0.75, 1.08, 'T', transform=m_ax[j].transAxes,
                            bbox=dict(boxstyle='square,pad=0.2', facecolor=colors[1], 
                                    edgecolor='black', linewidth=0.5),
                            fontsize=7, ha='center', va='center')
                
                # Remove the default title
                m_ax[j].set_title('')

    
    plots.append(("State Distribution", fig3))
    
    # Add all plots for this section to the HTML builder
    html_builder.add_section(section_title, plots)
    

def make_df(in_vcf, bed, sizemin=50, sizemax=10000):
    """
    Turn a VCF into the dataframe for parameter estimation
    """
    vcf = truvari.VariantFile(in_vcf)
    m_iter = vcf.fetch_bed(bed) if bed else vcf
    rows = []
    for entry in m_iter:
        if entry.chrom in ['chrX', 'chrY'] \
                or not (sizemin <= entry.var_size() <= sizemax) \
                or entry.is_monrefstar() \
                or None in entry.samples[1]['GT']:
            continue

        b_gt = truvari.get_gt(entry.gt(0))
        o_gt = truvari.get_gt(entry.gt(1))
        rows.append([b_gt == o_gt,
                     b_gt.name,
                     o_gt.name,
                     entry.samples[1]['DP'],
                     *entry.samples[1]['AD'],
                     entry.samples[1]['GQ'],
                     entry.samples[1]['FT'],
                     min(entry.samples[1]['KS']),
                     ])
    out = pd.DataFrame(rows, columns=['state', 'Ogt', 'Kgt',
                                      'DP', 'AD_ref', 'AD_alt', 'GQ', 'FT', 'KS'])
    return out


def regt(row, gt):
    """
    Runs kanpig genotyping on a row
    """
    result = gt.genotype(row['AD_ref'], 0, row['AD_alt'])
    return [result.state, result.state == row['Ogt'], int(round(result.gq))]


def calc_accuracy(df, prefix, html_builder=None):
    """
    Calculate and log genotype accuracy statistics.
    Optionally add tables to HTML report.

    Args:
        df: DataFrame with genotype data
        prefix: Prefix for log messages
        html_builder: Optional HTMLReportBuilder instance to add tables to
    """
    # Calculate accuracy table
    cnt = df.groupby(['Ogt', 'state']).size().unstack()
    cnt.loc['All'] = cnt.sum(axis=0)
    cnt['Acc'] = cnt[True] / cnt.sum(axis=1)
    table = textwrap.indent(cnt.to_string(), "    ")
    logging.info(f"{prefix} Genotype Accuracy:\n{table}")

    # Calculate confusion matrix
    cnt2 = df.groupby(["Ogt", "Kgt"]).size().unstack()
    table2 = textwrap.indent(cnt2.to_string(), "    ")
    logging.info(f"{prefix} GT Confusion Matrix:\n{table2}")

    # Add to HTML report if builder is provided
    if html_builder is not None:
        html_builder.add_table_section(
            section_title=f"{prefix} Statistics",
            table1_title="Genotype Accuracy",
            table1_df=cnt,
            table2_title="GT Confusion Matrix",
            table2_df=cnt2
        )

# Example usage
if __name__ == "__main__":
    args = parse_args(sys.argv[1:])

    truvari.setup_logging(stream=truvari.LogFileStderr(args.OUT + '.log'), show_version=False)

    if args.all and args.all_hets:
        logging.error("Can only fit either --all-hets XOR --all")
        sys.exit(1)
    
    report = HTMLReportBuilder(title="Kanpig GQ Calibration Report", output_dir=args.OUT + '_data')
    if args.IN.endswith("genotypes.csv"):
        df = pd.read_csv(args.IN)
    else:
        logging.info("Parsing VCF")
        df = make_df(args.IN, args.bed, args.sizemin, args.sizemax)

    calc_accuracy(df, "Full", report)

    if args.leaveout:
        leaveout = df.groupby(['Ogt']).sample(
            frac=args.leaveout, random_state=232)
        logging.info(f"Leaving out {args.leaveout} (N={len(leaveout)}) genotypes")
        calc_accuracy(leaveout, "Leaveout", report)
        df = df.drop(leaveout.index)

    logging.info("Summarizing Original Genotypes")
    make_plots(df, "Original", report)
    if args.summary:
        report.save(args.OUT + "_report.html")
        logging.info("Summary only - Finished")
        sys.exit(0)

    # Fit parameters
    fitted = fit_parameters(df,
                            all_gts=args.all,
                            min_dp=args.mindp,
                            max_dp=args.maxdp,
                            all_hets=args.all_hets,
                            )

    # Now you need to go re-genotype everything and grab those GQs
    # Then you make the calibration table
    out_cfg = args.OUT + '_gqconfig.json'
    config = save_config(fitted, out_cfg, flat_priors=args.flat_priors)

    if not args.no_calibrate:
        logging.info("Calibrating GQs")
        gt = kanpig.Genotyper.from_config_path(out_cfg)
        m_gtfunction = partial(regt, gt=gt)
        df[['nKgt', 'nState', 'nGQ']] = df.apply(
            m_gtfunction, axis=1, result_type='expand')
        calibration = build_calibration(df)
        config = save_config(fitted, out_cfg, calibration,
                             flat_priors=args.flat_priors)
        # And then we have to run again to actually get the calibrated GQs
        logging.info("Regenotyping with config")
        gt = kanpig.Genotyper.from_config_path(out_cfg)
        m_gtfunction = partial(regt, gt=gt)
        df[['nKgt', 'nState', 'nGQ']] = df.apply(
            m_gtfunction, axis=1, result_type='expand')

        if args.leaveout:
            leaveout[['nKgt', 'nState', 'nGQ']] = leaveout.apply(
                m_gtfunction, axis=1, result_type='expand')

    logging.info("Saving gqconfig")
    df.to_csv(args.OUT + '_data/genotypes.csv.gz', index=False, compression='gzip')
    if args.leaveout:
        leaveout.to_csv(args.OUT + '_data/leaveout_genotypes.csv.gz', index=False, compression='gzip')

    
    if not args.no_calibrate:
        logging.info("Making calibrated plots")
        df['GQ'] = df['nGQ']
        make_plots(df, 'Calibrated', report)
    if args.leaveout:
        logging.info("Making leaveout original plots",)
        make_plots(leaveout, "Leftout Original", report)
        if not args.no_calibrate:
            logging.info("Making leaveout calibrated plots")
            leaveout['GQ'] = leaveout['nGQ']
            make_plots(leaveout, "Leftout Calibrated", report)
    report.save(args.OUT + "_report.html")

    logging.info("Finished")
