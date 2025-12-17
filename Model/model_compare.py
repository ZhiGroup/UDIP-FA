
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import seaborn as sns
from sklearn.cross_decomposition import CCA
from scipy.stats import pearsonr, spearmanr
from itertools import combinations
import os

# Set Nature journal style globally
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 11,
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.2,
    'ytick.major.width': 1.2,
    'xtick.major.size': 5,
    'ytick.major.size': 5,
    'axes.spines.top': False,
    'axes.spines.right': False
})

def plot_snps_different_thresholds(csv_files, thresholds, labels=None, 
                                   p_column='P', output_file='snps_different_thresholds.pdf',
                                   figsize=(10, 6), dpi=300):
    """
    Count significant SNPs across multiple CSV files using different P-value thresholds
    and create a Nature-style line plot.
    """
    # Check if thresholds length matches csv_files
    if len(thresholds) != len(csv_files):
        raise ValueError("Length of thresholds must match length of csv_files")
    
    # Default labels if not provided
    if labels is None:
        labels = [f'File{i+1}' for i in range(len(csv_files))]
    
    # Store results
    snp_counts = []
    log_thresholds = []
    
    # Process each file and count significant SNPs
    for i, (csv_file, threshold) in enumerate(zip(csv_files, thresholds)):
        print(f"Processing {csv_file} with P < {threshold}...")
        
        try:
            # Read CSV file (handles both comma and tab-separated)
            df = pd.read_csv(csv_file, sep='\t' if csv_file.endswith('.tsv') else ',')
            
            # Check if P-value column exists
            if p_column not in df.columns:
                # Try common P-value column names
                possible_cols = ['P', 'p', 'pval', 'P_VALUE', 'p_value']
                for col in possible_cols:
                    if col in df.columns:
                        p_column = col
                        break
            
            if p_column in df.columns:
                # Count significant SNPs below threshold
                count = (df[p_column] < threshold).sum()
                print(f"  {labels[i]}: {count} SNPs (P < {threshold}, -log10(P) > {-np.log10(threshold):.2f})")
            else:
                print(f"  Warning: P-value column not found in {csv_file}. assuming 0.")
                count = 0
        except Exception as e:
            print(f"  Error reading {csv_file}: {e}")
            count = 0

        snp_counts.append(count)
        log_thresholds.append(-np.log10(threshold))
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    # Nature journal color palette
    colors = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', '#8491B4']
    markers = ['o', 's', '^', 'v', 'D', 'p']
    
    # Plot each point with different marker and color
    for i in range(len(csv_files)):
        color = colors[i % len(colors)]
        marker = markers[i % len(markers)]
        ax.plot(i, snp_counts[i], 
                marker=marker, color=color, 
                markersize=12, markeredgewidth=2, 
                markeredgecolor='white', alpha=0.9,
                label=f'{labels[i]} (P<{thresholds[i]:.0e})')
    
    # Connect points with gray dashed line
    ax.plot(range(len(csv_files)), snp_counts, 
            color='gray', linewidth=2, linestyle='--', 
            alpha=0.4, zorder=0)
    
    # Add count values above each point
    for i, count in enumerate(snp_counts):
        ax.text(i, count, f'{count:,}',
                ha='center', va='bottom', fontsize=9, 
                fontweight='bold', bbox=dict(boxstyle='round,pad=0.3', 
                facecolor='white', edgecolor='none', alpha=0.7))
    
    # Set axis labels
    ax.set_xlabel(' Dimension of UDIP-FA', fontsize=13, fontweight='bold')
    ax.set_ylabel('Number of significant SNPs', fontsize=13, fontweight='bold')
    ax.set_title('Genome-wide significant SNPs after bonferroni correction', 
                 fontsize=14, fontweight='bold', pad=15)
    
    # Set x-axis ticks and labels
    ax.set_xticks(range(len(csv_files)))
    ax.set_xticklabels(labels, rotation=45, ha='right')
    
    # Format y-axis with thousand separators
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))
    
    # Add horizontal grid lines
    ax.grid(True, axis='y', alpha=0.2, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)
    
    # Add legend
    ax.legend(loc='best', fontsize=9, frameon=False)
    
    # Adjust y-axis limits to accommodate text labels
    y_min, y_max = ax.get_ylim()
    ax.set_ylim(y_min, y_max * 1.15)
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save figure
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    print(f"\nFigure saved to {output_file}")
    
    # Also save as PNG format
    png_file = output_file.replace('.pdf', '.png')
    plt.savefig(png_file, dpi=dpi, bbox_inches='tight')
    print(f"Figure also saved to {png_file}")
    
    # plt.show() # Commented out for non-interactive environments
    
    # Create results dataframe
    results = pd.DataFrame({
        'Dataset': labels,
        'P_threshold': thresholds,
        'Neg_log10_P': log_thresholds,
        'Significant_SNPs': snp_counts
    })
    
    return results

def load_pickle_features(pickle_files, labels=None):
    """
    Load multiple pickle files containing feature matrices
    """
    if labels is None:
        labels = [f'Features_{i+1}' for i in range(len(pickle_files))]
    
    features_dict = {}
    
    for i, pkl_file in enumerate(pickle_files):
        print(f"Loading {pkl_file}...")
        try:
            with open(pkl_file, 'rb') as file:
                data = pickle.load(file)
                features = np.asarray(data)
                features_dict[labels[i]] = features
                print(f"  {labels[i]}: shape {features.shape}")
        except Exception as e:
            print(f"  Error loading {pkl_file}: {e}")
    
    return features_dict

def compute_pairwise_correlations(features_dict, method='pearson'):
    """
    Compute pairwise correlations between feature sets
    """
    labels = list(features_dict.keys())
    n_features = len(labels)
    corr_matrix = np.zeros((n_features, n_features))
    pvalue_matrix = np.zeros((n_features, n_features))
    
    for i, label1 in enumerate(labels):
        for j, label2 in enumerate(labels):
            if i == j:
                corr_matrix[i, j] = 1.0
                pvalue_matrix[i, j] = 0.0
            elif i < j:
                feat1 = features_dict[label1]
                feat2 = features_dict[label2]
                
                # Ensure same number of samples
                n_samples = min(feat1.shape[0], feat2.shape[0])
                feat1 = feat1[:n_samples, :]
                feat2 = feat2[:n_samples, :]
                
                if method == 'pearson':
                    # Average correlation across all feature pairs
                    corrs = []
                    for k in range(min(feat1.shape[1], feat2.shape[1])):
                        r, p = pearsonr(feat1[:, k], feat2[:, k])
                        corrs.append(r)
                    corr_matrix[i, j] = np.mean(corrs)
                    corr_matrix[j, i] = corr_matrix[i, j]
                    
                elif method == 'spearman':
                    # Average Spearman correlation
                    corrs = []
                    for k in range(min(feat1.shape[1], feat2.shape[1])):
                        r, p = spearmanr(feat1[:, k], feat2[:, k])
                        corrs.append(r)
                    corr_matrix[i, j] = np.mean(corrs)
                    corr_matrix[j, i] = corr_matrix[i, j]
                    
                elif method == 'cca':
                    # Use CCA to find canonical correlation
                    n_components = min(feat1.shape[1], feat2.shape[1], 2)  # Use top 2 components as in original code
                    cca = CCA(n_components=n_components)
                    
                    try:
                        cca.fit(feat1, feat2)
                        X_c, Y_c = cca.transform(feat1, feat2)
                        
                        # Compute correlation for each canonical variable
                        canonical_corrs = []
                        for k in range(n_components):
                            r, _ = pearsonr(X_c[:, k], Y_c[:, k])
                            canonical_corrs.append(r)
                        
                        # Use the first (highest) canonical correlation
                        corr_matrix[i, j] = canonical_corrs[0]
                        corr_matrix[j, i] = corr_matrix[i, j]
                        
                        print(f"  CCA {label1} vs {label2}: canonical correlations = {canonical_corrs}")
                        
                    except Exception as e:
                        print(f"  Warning: CCA failed for {label1} vs {label2}: {e}")
                        corr_matrix[i, j] = 0
                        corr_matrix[j, i] = 0
    
    return corr_matrix, pvalue_matrix

def plot_correlation_heatmap(corr_matrix, labels, method='CCA', 
                             output_file='correlation_heatmap.pdf',
                             figsize=(10, 8), dpi=300):
    """
    Plot correlation heatmap in Nature style
    """
    # Create figure
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    # Create heatmap
    im = ax.imshow(corr_matrix, cmap='RdYlBu_r', vmin=-1, vmax=1, aspect='auto')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(f'{method} Correlation', fontsize=12, fontweight='bold')
    cbar.ax.tick_params(labelsize=10)
    
    # Set ticks and labels
    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_yticklabels(labels)
    
    # Add correlation values as text
    for i in range(len(labels)):
        for j in range(len(labels)):
            color = "black" if abs(corr_matrix[i, j]) < 0.5 else "white"
            text = ax.text(j, i, f'{corr_matrix[i, j]:.3f}',
                          ha="center", va="center", color=color,
                          fontsize=9, fontweight='bold')
    
    # Set title
    ax.set_title(f'Pairwise {method} Correlations Between Feature Sets', 
                 fontsize=14, fontweight='bold', pad=15)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    print(f"\nHeatmap saved to {output_file}")
    
    png_file = output_file.replace('.pdf', '.png')
    plt.savefig(png_file, dpi=dpi, bbox_inches='tight')
    print(f"Heatmap also saved to {png_file}")
    
    # plt.show()

def comprehensive_correlation_analysis(pickle_files, labels=None, methods=['cca', 'pearson'],
                                      output_prefix='feature_correlation'):
    """
    Comprehensive correlation analysis for multiple pickle files
    """
    # Load all feature files
    features_dict = load_pickle_features(pickle_files, labels)
    if not features_dict:
        print("No features loaded. Exiting correlation analysis.")
        return {}, features_dict

    labels = list(features_dict.keys())
    
    results = {}
    
    # Compute correlations using different methods
    for method in methods:
        print(f"\n{'='*60}")
        print(f"Computing {method.upper()} correlations...")
        print(f"{'='*60}")
        
        corr_matrix, pvalue_matrix = compute_pairwise_correlations(features_dict, method=method)
        results[method] = {
            'correlation': corr_matrix,
            'pvalue': pvalue_matrix
        }
        
        # Plot heatmap
        plot_correlation_heatmap(
            corr_matrix, 
            labels, 
            method=method.upper(),
            output_file=f'{output_prefix}_{method}_heatmap.pdf'
        )
        
        # Save correlation matrix to CSV
        df_corr = pd.DataFrame(corr_matrix, index=labels, columns=labels)
        df_corr.to_csv(f'{output_prefix}_{method}_matrix.csv')
        print(f"\nCorrelation matrix saved to {output_prefix}_{method}_matrix.csv")
    
    return results, features_dict

def run_snp_analysis():
    print("\nRunning SNP Analysis...")
    # Your 6 CSV files
    csv_files = [
        '/data/whe4/DeepENDO/AE_interpretation_2023_2_gpus/small_lr/dti_32_interpretation/minp_32.csv',
        '/data/whe4/DeepENDO/AE_interpretation_2023_2_gpus/small_lr/dti_64_interpretation/minp_64.csv',
        '/data/whe4/DeepENDO/AE_interpretation_2023_2_gpus/small_lr/dti_128_interpretation/minp_128.csv',
        '/data/whe4/DeepENDO/AE_interpretation_2023_2_gpus/small_lr/dti_256_interpretation/minp_256.csv',
        '/data/whe4/DeepENDO/AE_interpretation_2023_2_gpus/small_lr/dti_512_interpretation/minp_512.csv',
        '/data/whe4/DeepENDO/AE_interpretation_2023_2_gpus/small_lr/dti_1024_interpretation/minp_1024.csv'
    ]

    labels = ['UDIP-FA(32)', 'UDIP-FA(64)', 'UDIP-FA(128)', 'UDIP-FA(256)', 'UDIP-FA(512)', 'UDIP-FA(1024)']

    # Different P-value threshold for each file
    thresholds = [
        5e-8/32,   # File1
        5e-8/64,   # File2
        5e-8/128,  # File3
        5e-8/256,  # File4
        5e-8/512,  # File5
        5e-8/1028  # File6 (Note: 1028 in original code, likely meant 1024 but kept as is for fidelity)
    ]
    
    # Run analysis
    results = plot_snps_different_thresholds(
        csv_files=csv_files,
        thresholds=thresholds,
        labels=labels,
        p_column='P_BOLT_LMM_INF',
        output_file='/data484_2/xzhao14/FA_revision/significant_snps_5e8.pdf'
    )
    
    print("\nSNP Analysis Statistics summary:")
    print(results)
    return results

def run_feature_analysis():
    print("\nRunning Feature Correlation Analysis...")
    # Your pickle files
    pickle_files = [
        '/data/whe4/DeepENDO/AE_interpretation_2023_2_gpus/small_lr/dti_32_interpretation/dti_32_endo.pkl',
        '/data/whe4/DeepENDO/AE_interpretation_2023_2_gpus/small_lr/dti_64_interpretation/dti_64_endo.pkl',
        '/data/whe4/DeepENDO/AE_interpretation_2023_2_gpus/small_lr/dti_128_interpretation/dti_128_endo.pkl',
        '/data/whe4/DeepENDO/AE_interpretation_2023_2_gpus/small_lr/dti_256_interpretation/dti_256_endo.pkl',
        '/data/whe4/DeepENDO/AE_interpretation_2023_2_gpus/small_lr/dti_512_interpretation/dti_512_endo.pkl',
        '/data/whe4/DeepENDO/AE_interpretation_2023_2_gpus/small_lr/dti_1024_interpretation/dti_1024_endo.pkl'
    ]
    
    # Labels for each feature set
    labels = ['UDIP-FA(32)', 'UDIP-FA(64)', 'UDIP-FA(128)', 
              'UDIP-FA(256)', 'UDIP-FA(512)', 'UDIP-FA(1024)']
    
    # Run comprehensive analysis
    results, features_dict = comprehensive_correlation_analysis(
        pickle_files=pickle_files,
        labels=labels,
        methods=['cca', 'pearson'],  # Can add 'spearman' if needed
        output_prefix='/data484_2/xzhao14/FA_revision/udip_fa_features'
    )
    
    print("\nFeature Analysis complete!")
    return results

if __name__ == "__main__":
    # You can comment out one or the other if you only want to run one analysis
    
    # 1. SNP Threshold Analysis
    try:
        run_snp_analysis()
    except Exception as e:
        print(f"SNP analysis failed or skipped: {e}")

    print("\n" + "="*60 + "\n")

    # 2. Feature Correlation Analysis
    try:
        run_feature_analysis()
    except Exception as e:
        print(f"Feature analysis failed or skipped: {e}")