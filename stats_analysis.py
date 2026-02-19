"""
Comprehensive Statistical Analysis Pipeline for Proteomics Data.

Features:
- Auto-detection of Peptide vs Protein data
- Preprocessing: Filter Missing -> Impute -> Log Transform
- Statistical Analysis: Mann-Whitney U, Fold Change, FDR correction
- Visualization: PCA, Volcano Plot, Heatmap
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.metrics import accuracy_score, roc_auc_score
import warnings

warnings.filterwarnings('ignore')

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

class Config:
    # Default Paths
    BASE_PATH = r
    DATA_FILE = f"{BASE_PATH}/"
    LABEL_FILE = f"{BASE_PATH}/"
    
    # Groups
    CASE_NAME = "Case"
    CONTROL_NAME = "Control"
    GROUP_COL = "Group"
    SAMPLE_ID_COL = "Sample_ID"
    
    # Preprocessing
    MAX_MISSING = 0.30
    IMPUTATION = "minimum"  # 'minimum', 'median', 'mean'
    LOG_TRANSFORM = True
    
    # Output
    OUTPUT_DIR = "stats_results"

def load_data(data_file, label_file):
    print(f"Loading data from: {data_file}")
    df = pd.read_csv(data_file)
    
    # Metadata columns
    meta_cols = ['Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 
                 'First.Protein.Description', 'Proteotypic', 'Stripped.Sequence', 
                 'Modified.Sequence', 'Precursor.Charge', 'Precursor.Id']
    meta_cols = [c for c in meta_cols if c in df.columns]
    sample_cols = [c for c in df.columns if c not in meta_cols]
    
    # Detect data type
    if 'Stripped.Sequence' in df.columns:
        print("Detected peptide-level data.")
    else:
        print("Detected protein-level data.")
    
    # Extract metadata and data
    metadata = df[meta_cols].copy()
    X = df[sample_cols].T
    X.columns = df.index
    
    # Load labels
    print(f"Loading labels from: {label_file}")
    labels = pd.read_excel(label_file)
    
    # Map samples
    sample_map = dict(zip(labels[Config.SAMPLE_ID_COL], labels[Config.GROUP_COL]))
    
    # Filter samples
    valid_samples = [s for s in X.index if s in sample_map]
    X = X.loc[valid_samples]
    y = np.array([1 if sample_map[s] == Config.CASE_NAME else 0 for s in X.index])
    
    print(f"Analysis set: {len(X)} samples ({np.sum(y==1)} {Config.CASE_NAME}, {np.sum(y==0)} {Config.CONTROL_NAME})")
    
    return X, y, metadata

def preprocess_data(X, config):
    print("\nPreprocessing...")
    
    # 1. Filter Missing
    missing_frac = X.isnull().mean()
    keep_cols = missing_frac[missing_frac <= config.MAX_MISSING].index
    X_filt = X[keep_cols]
    print(f"  Filtered features > {config.MAX_MISSING*100}% missing: {X.shape[1]} -> {X_filt.shape[1]}")
    
    # 2. Impute
    print(f"  Imputing missing values ({config.IMPUTATION})...")
    X_num = X_filt.apply(pd.to_numeric, errors='coerce')
    
    if config.IMPUTATION == 'minimum':
        fill_val = X_num.min()
        X_imp = X_num.fillna(fill_val)
    elif config.IMPUTATION == 'median':
        X_imp = X_num.fillna(X_num.median())
    elif config.IMPUTATION == 'mean':
        X_imp = X_num.fillna(X_num.mean())
    else:
        X_imp = X_num.fillna(0)
        
    # 3. Log Transform
    if config.LOG_TRANSFORM:
        print("  Applying Log2(x+1) transformation...")
        X_final = np.log2(X_imp + 1)
    else:
        X_final = X_imp
        
    return X_final

def run_statistics(X, y, metadata, X_raw=None):
    print("\nRunning Statistical Tests (Mann-Whitney U)...")
    results = []

    case_mask = (y == 1)
    ctrl_mask = (y == 0)

    for feat_id in X.columns:
        case_vals = X.loc[case_mask, feat_id]
        ctrl_vals = X.loc[ctrl_mask, feat_id]

        try:
            stat, p = mannwhitneyu(case_vals, ctrl_vals, alternative='two-sided')
        except:
            stat, p = 0, 1.0

        # Store log-transformed means for reference
        mean_case_log = case_vals.mean()
        mean_ctrl_log = ctrl_vals.mean()

        # Calculate non-log transformed means and fold changes from raw data
        if X_raw is not None and feat_id in X_raw.columns:
            case_vals_raw = X_raw.loc[case_mask, feat_id]
            ctrl_vals_raw = X_raw.loc[ctrl_mask, feat_id]
            mean_case_raw = case_vals_raw.mean()
            mean_ctrl_raw = ctrl_vals_raw.mean()
            # Fold change (not log-transformed)
            fc = mean_case_raw / mean_ctrl_raw if mean_ctrl_raw != 0 else np.nan
            log2fc = np.log2(fc) if not np.isnan(fc) and fc > 0 else np.nan
        else:
            mean_case_raw = np.nan
            mean_ctrl_raw = np.nan
            fc = np.nan
            log2fc = np.nan

        # Info
        info = {'feature_id': feat_id}
        if feat_id in metadata.index:
            row = metadata.loc[feat_id]
            if 'Genes' in row: info['Gene'] = row['Genes']
            if 'Stripped.Sequence' in row: info['Peptide'] = row['Stripped.Sequence']
            if 'Protein.Names' in row: info['Protein'] = row['Protein.Names']

        results.append({
            **info,
            'p_value': p,
            'log2fc': log2fc,
            'mean_case_log2': mean_case_log,
            'mean_ctrl_log2': mean_ctrl_log,
            'mean_case_raw': mean_case_raw,
            'mean_ctrl_raw': mean_ctrl_raw,
            'fold_change': fc
        })

    res_df = pd.DataFrame(results)

    # FDR Correction
    _, adj_p, _, _ = multipletests(res_df['p_value'], method='fdr_bh')
    res_df['adj_p_value'] = adj_p

    res_df = res_df.sort_values('p_value')
    print(f"  Significant features (FDR < 0.05): {sum(res_df['adj_p_value'] < 0.05)}")

    # Reorder columns: metadata first, then raw means/FC, then log2 means/FC, then p-values
    metadata_cols = ['feature_id', 'Gene', 'Peptide', 'Protein']
    stats_cols = ['mean_case_raw', 'mean_ctrl_raw', 'fold_change',
                  'mean_case_log2', 'mean_ctrl_log2', 'log2fc',
                  'p_value', 'adj_p_value']

    # Only include columns that exist in the DataFrame
    existing_metadata = [c for c in metadata_cols if c in res_df.columns]
    existing_stats = [c for c in stats_cols if c in res_df.columns]

    res_df = res_df[existing_metadata + existing_stats]

    return res_df

def plot_pca(X, y, output_dir):
    print("Generating PCA...")
    pca = PCA(n_components=2)
    X_scaled = StandardScaler().fit_transform(X)
    coords = pca.fit_transform(X_scaled)

    plt.figure(figsize=(8, 6))
    groups = [Config.CASE_NAME if label==1 else Config.CONTROL_NAME for label in y]

    # Add confidence ellipses (95%) - drawn first so they appear behind points
    from matplotlib.patches import Ellipse
    from scipy import stats

    unique_groups = [Config.CASE_NAME, Config.CONTROL_NAME]
    colors = {Config.CASE_NAME: '#440154', Config.CONTROL_NAME: '#fde724'}  # viridis palette colors

    # Track ellipse bounds for axis limits
    ellipse_bounds = []

    for group in unique_groups:
        # Get points for this group
        mask = np.array([g == group for g in groups])
        group_coords = coords[mask]

        if len(group_coords) > 2:  # Need at least 3 points for ellipse
            # Calculate mean and covariance
            mean = group_coords.mean(axis=0)
            cov = np.cov(group_coords.T)

            # Calculate eigenvalues and eigenvectors
            eigenvalues, eigenvectors = np.linalg.eigh(cov)

            # Calculate angle of ellipse
            angle = np.degrees(np.arctan2(eigenvectors[1, 1], eigenvectors[0, 1]))

            # Width and height are 2*sqrt(eigenvalue) * chi-square value for 95% confidence
            chi2_val = stats.chi2.ppf(0.95, df=2)
            # eigenvalues are sorted ascending, so eigenvalues[1] is the major axis
            width = 2 * np.sqrt(eigenvalues[1] * chi2_val)
            height = 2 * np.sqrt(eigenvalues[0] * chi2_val)

            # Create semi-transparent filled ellipse without border
            ellipse = Ellipse(xy=mean, width=width, height=height, angle=angle,
                            facecolor=colors.get(group, 'gray'), edgecolor='none',
                            alpha=0.2, zorder=0)
            plt.gca().add_patch(ellipse)

            # Track bounds for this ellipse
            ellipse_bounds.append([mean[0] - width/2, mean[0] + width/2,
                                 mean[1] - height/2, mean[1] + height/2])

    # Plot scatter points on top
    sns.scatterplot(x=coords[:,0], y=coords[:,1], hue=groups, palette='viridis', s=100, alpha=0.8, zorder=1)

    # Adjust axis limits to accommodate ellipses
    if ellipse_bounds:
        all_bounds = np.array(ellipse_bounds)
        x_min, x_max = all_bounds[:, 0].min(), all_bounds[:, 1].max()
        y_min, y_max = all_bounds[:, 2].min(), all_bounds[:, 3].max()

        # Add 10% padding
        x_range = x_max - x_min
        y_range = y_max - y_min
        plt.xlim(x_min - 0.1 * x_range, x_max + 0.1 * x_range)
        plt.ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)

    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
    plt.title("PCA Analysis with 95% Confidence Ellipses")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/pca_plot.png")
    plt.close()

def plot_plsda(X, y, metadata, output_dir, n_components=2):
    """
    Perform PLS-DA with cross-validation and VIP score extraction.
    """
    print("Generating PLS-DA...")

    # Scale the data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Fit PLS-DA
    plsda = PLSRegression(n_components=n_components)
    plsda.fit(X_scaled, y)

    # Get scores (sample projections)
    scores = plsda.x_scores_

    # Cross-validation to assess model performance
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    y_pred_cv = cross_val_predict(PLSRegression(n_components=n_components),
                                   X_scaled, y, cv=cv)
    y_pred_class = (y_pred_cv > 0.5).astype(int).ravel()

    accuracy = accuracy_score(y, y_pred_class)
    auc = roc_auc_score(y, y_pred_cv)

    # Calculate VIP scores
    # VIP = sqrt(p * sum(w^2 * SSY_explained) / SSY_total)
    t = plsda.x_scores_  # Scores
    w = plsda.x_weights_  # Weights
    q = plsda.y_loadings_  # Y loadings

    p = X_scaled.shape[1]  # Number of features
    h = n_components  # Number of components

    # Calculate VIP
    vip = np.zeros(p)
    s = np.diag(t.T @ t @ q.T @ q).reshape(h, -1)
    total_s = np.sum(s)

    for i in range(p):
        weight = np.array([(w[i, j] / np.linalg.norm(w[:, j]))**2 for j in range(h)])
        vip[i] = np.sqrt(p * (s.T @ weight) / total_s)

    vip_scores = vip.ravel()

    # Create VIP DataFrame
    vip_df = pd.DataFrame({
        'feature_id': X.columns,
        'VIP_score': vip_scores
    })

    # Add gene names if available
    gene_names = []
    for fid in X.columns:
        if fid in metadata.index and 'Genes' in metadata.columns:
            gene_names.append(metadata.loc[fid, 'Genes'])
        else:
            gene_names.append(fid)
    vip_df['Gene'] = gene_names

    vip_df = vip_df.sort_values('VIP_score', ascending=False)
    vip_df.to_csv(f"{output_dir}/plsda_vip_scores.csv", index=False)

    print(f"  PLS-DA 5-fold CV: Accuracy={accuracy:.2f}, AUC={auc:.2f}")
    print(f"  Features with VIP > 1: {sum(vip_scores > 1)}")

    # Plot 1: PLS-DA scores plot
    plt.figure(figsize=(8, 6))
    groups = [Config.CASE_NAME if label == 1 else Config.CONTROL_NAME for label in y]

    # Add confidence ellipses
    from matplotlib.patches import Ellipse
    from scipy import stats

    unique_groups = [Config.CASE_NAME, Config.CONTROL_NAME]
    colors = {Config.CASE_NAME: '#440154', Config.CONTROL_NAME: '#fde724'}

    for group in unique_groups:
        mask = np.array([g == group for g in groups])
        group_coords = scores[mask]

        if len(group_coords) > 2:
            mean = group_coords.mean(axis=0)
            cov = np.cov(group_coords.T)
            eigenvalues, eigenvectors = np.linalg.eigh(cov)
            angle = np.degrees(np.arctan2(eigenvectors[1, 1], eigenvectors[0, 1]))
            chi2_val = stats.chi2.ppf(0.95, df=2)
            # eigenvalues are sorted ascending, so eigenvalues[1] is the major axis
            width = 2 * np.sqrt(eigenvalues[1] * chi2_val)
            height = 2 * np.sqrt(eigenvalues[0] * chi2_val)

            ellipse = Ellipse(xy=mean, width=width, height=height, angle=angle,
                              facecolor=colors.get(group, 'gray'), edgecolor='none',
                              alpha=0.2, zorder=0)
            plt.gca().add_patch(ellipse)

    sns.scatterplot(x=scores[:, 0], y=scores[:, 1], hue=groups,
                    palette='viridis', s=100, alpha=0.8, zorder=1)

    plt.xlabel(f"PLS Component 1")
    plt.ylabel(f"PLS Component 2")
    plt.title(f"PLS-DA Analysis\n5-fold CV: Accuracy={accuracy:.2f}, AUC={auc:.2f}")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/plsda_plot.png")
    plt.close()

    # Plot 2: VIP scores bar plot (top 20)
    plt.figure(figsize=(10, 8))
    top_vip = vip_df.head(20)

    # Create bar colors based on VIP > 1 threshold
    bar_colors = ['#e74c3c' if v > 1 else '#3498db' for v in top_vip['VIP_score']]

    plt.barh(range(len(top_vip)), top_vip['VIP_score'].values, color=bar_colors)
    plt.yticks(range(len(top_vip)), top_vip['Gene'].values)
    plt.xlabel('VIP Score')
    plt.ylabel('Feature')
    plt.title('PLS-DA Variable Importance in Projection (VIP)\nTop 20 Features (Red = VIP > 1)')
    plt.axvline(x=1, color='k', linestyle='--', alpha=0.5, label='VIP = 1 threshold')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(f"{output_dir}/plsda_vip_plot.png")
    plt.close()

    return vip_df

def plot_volcano(res_df, output_dir):
    print("Generating Volcano Plot...")
    plt.figure(figsize=(10, 8))

    res_df['nlog10p'] = -np.log10(res_df['p_value'])

    # Define thresholds
    fdr_threshold = 0.05
    log2fc_threshold = 0.585

    # Color coding - only color features that pass BOTH thresholds
    colors = []
    for _, row in res_df.iterrows():
        if row['adj_p_value'] < fdr_threshold and abs(row['log2fc']) > log2fc_threshold:
            colors.append('red')
        else:
            colors.append('grey')

    plt.scatter(res_df['log2fc'], res_df['nlog10p'], c=colors, alpha=0.6, s=20)

    # Add threshold lines
    # Horizontal line for FDR threshold (convert to p-value line for visual reference)
    plt.axhline(-np.log10(0.05), color='k', linestyle='--', alpha=0.3, label='p=0.05')
    plt.axvline(log2fc_threshold, color='k', linestyle='--', alpha=0.3)
    plt.axvline(-log2fc_threshold, color='k', linestyle='--', alpha=0.3)

    # Label only features that pass both thresholds (with overlap detection)
    significant = res_df[(res_df['adj_p_value'] < fdr_threshold) & (abs(res_df['log2fc']) > log2fc_threshold)]

    # Limit number of labels to prevent overcrowding
    max_labels = 30

    # Sort by significance (lowest p-value first) and limit
    significant_sorted = significant.sort_values('p_value').head(max_labels)

    # Track label positions to avoid severe overlaps
    labeled_positions = []
    min_distance = 0.3  # Minimum distance between labels (in plot units)

    for _, row in significant_sorted.iterrows():
        x, y = row['log2fc'], row['nlog10p']

        # Check if position is too close to existing labels
        too_close = False
        for prev_x, prev_y in labeled_positions:
            distance = np.sqrt((x - prev_x)**2 + (y - prev_y)**2)
            if distance < min_distance:
                too_close = True
                break

        if not too_close:
            label = row.get('Gene', row.get('Peptide', str(row['feature_id'])))
            # Shorten long gene names with multiple entries
            if isinstance(label, str) and ';' in label:
                genes = label.split(';')
                if len(genes) > 2:
                    label = f"{genes[0]};{genes[1]}...+{len(genes)-2}"
            plt.text(x, y, str(label), fontsize=8, alpha=0.7, fontweight='bold')
            labeled_positions.append((x, y))

    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-Log10 P-value")
    plt.title("Volcano Plot")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/volcano_plot.png")
    plt.close()

def plot_heatmap(X, y, res_df, output_dir, top_n=50):
    print(f"Generating Heatmaps (Top {top_n})...")
    if len(res_df) < 2:
        print("Not enough features for heatmap.")
        return

    top_feats = res_df.head(top_n)['feature_id'].values
    X_sub = X[top_feats]

    # Scale for visualization (z-score)
    X_z = pd.DataFrame(StandardScaler().fit_transform(X_sub), index=X_sub.index, columns=X_sub.columns)

    # Map feature IDs to names
    labels = []
    for fid in X_sub.columns:
        row = res_df[res_df['feature_id'] == fid].iloc[0]
        labels.append(row.get('Gene', row.get('Peptide', str(fid))))
    X_z.columns = labels

    # Create group labels for coloring
    sample_labels = [Config.CASE_NAME if label==1 else Config.CONTROL_NAME for label in y]

    # Create color map for samples
    sample_colors = pd.Series(sample_labels, index=X_z.index)
    lut = {Config.CASE_NAME: "r", Config.CONTROL_NAME: "b"}
    row_colors = sample_colors.map(lut)

    try:
        # Heatmap 1: Group labels (original)
        X_z_grouped = X_z.copy()
        X_z_grouped.index = sample_labels

        g = sns.clustermap(X_z_grouped.T, col_colors=row_colors, cmap="vlag", center=0,
                           figsize=(12, 10), dendrogram_ratio=(0.1, 0.2),
                           cbar_pos=(0.02, 0.8, 0.03, 0.15))
        g.ax_heatmap.set_xlabel("Samples")
        g.ax_heatmap.set_ylabel("Features")
        plt.savefig(f"{output_dir}/heatmap_grouped.png")
        plt.close()

        # Heatmap 2: Individual sample IDs
        g2 = sns.clustermap(X_z.T, col_colors=row_colors, cmap="vlag", center=0,
                            figsize=(14, 10), dendrogram_ratio=(0.1, 0.2),
                            cbar_pos=(0.02, 0.8, 0.03, 0.15),
                            xticklabels=True, yticklabels=True)
        g2.ax_heatmap.set_xlabel("Samples")
        g2.ax_heatmap.set_ylabel("Features")

        # Force all sample labels to show with proper formatting
        g2.ax_heatmap.tick_params(axis='x', labelsize=6, rotation=90)
        g2.ax_heatmap.tick_params(axis='y', labelsize=8)

        plt.savefig(f"{output_dir}/heatmap_sample_ids.png", bbox_inches='tight')
        plt.close()

    except Exception as e:
        print(f"Heatmap generation failed: {e}")

def plot_boxplots(X, y, res_df, output_dir, top_n=10):
    print(f"Generating Boxplots (Top {top_n})...")
    if len(res_df) == 0:
        print("No features to plot.")
        return
        
    # Select top N features
    top_feats = res_df.head(top_n)
    
    # Setup grid
    n_cols = 5
    n_rows = int(np.ceil(len(top_feats) / n_cols))
    
    plt.figure(figsize=(15, 3 * n_rows))
    
    groups = [Config.CASE_NAME if label==1 else Config.CONTROL_NAME for label in y]
    
    for idx, (_, row) in enumerate(top_feats.iterrows()):
        feat_id = row['feature_id']
        label = row.get('Gene', row.get('Peptide', str(feat_id)))
        
        plt.subplot(n_rows, n_cols, idx + 1)
        
        # Data for plotting
        plot_data = pd.DataFrame({
            'Log2 Abundance': X[feat_id].values,
            'Group': groups
        })
        
        # Plot
        sns.boxplot(x='Group', y='Log2 Abundance', data=plot_data, 
                   palette={Config.CASE_NAME: "#e74c3c", Config.CONTROL_NAME: "#3498db"},
                   width=0.5, showfliers=False)
        sns.stripplot(x='Group', y='Log2 Abundance', data=plot_data, 
                     color='black', alpha=0.3, size=3, jitter=True)
                     
        plt.title(f"{label}\nFDR={row['adj_p_value']:.2e}", fontsize=9)
        plt.xlabel("")
        
    plt.tight_layout()
    plt.savefig(f"{output_dir}/top_features_boxplots.png")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Proteomics Statistical Analysis")
    parser.add_argument("--data", default=Config.DATA_FILE, help="Path to data matrix")
    parser.add_argument("--labels", default=Config.LABEL_FILE, help="Path to labels file")
    parser.add_argument("--output", default=Config.OUTPUT_DIR, help="Output directory")
    args = parser.parse_args()
    
    # Setup output
    out_path = Path(args.output)
    out_path.mkdir(exist_ok=True, parents=True)
    
    print("="*60)
    print("COMPREHENSIVE STATISTICAL ANALYSIS")
    print("="*60)
    
    # Pipeline
    try:
        X_raw, y, metadata = load_data(args.data, args.labels)
        X_proc = preprocess_data(X_raw, Config)

        # Statistics (pass both processed and raw data)
        res_df = run_statistics(X_proc, y, metadata, X_raw=X_raw)
        res_df.to_csv(out_path / "statistical_results.csv", index=False)
        
        # Plots
        plot_pca(X_proc, y, out_path)
        plot_plsda(X_proc, y, metadata, out_path)
        plot_volcano(res_df, out_path)
        plot_heatmap(X_proc, y, res_df, out_path)
        plot_boxplots(X_proc, y, res_df, out_path)

        print(f"\nAnalysis Complete. Results saved to {out_path.absolute()}")
        
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()