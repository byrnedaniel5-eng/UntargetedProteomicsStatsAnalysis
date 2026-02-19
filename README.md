# Untargeted Proteomics Stats Analysis

A one-stop statistical analysis pipeline for untargeted proteomics data. Takes DIA-NN output directly and runs a full analysis including preprocessing, differential expression testing, and visualisations.

---

## Features

- **Auto-detection** of peptide-level vs protein-level DIA-NN output
- **Preprocessing pipeline**: missingness filtering → imputation → log2 transformation
- **Statistical testing**: Mann-Whitney U test with Benjamini-Hochberg FDR correction
- **Fold change calculation** from raw (non-log-transformed) intensities
- **Visualisations**:
  - PCA scores plot with 95% confidence ellipses
  - PLS-DA scores plot with 95% confidence ellipses and cross-validated performance metrics
  - VIP (Variable Importance in Projection) scores bar chart
  - Volcano plot with gene/peptide labelling
  - Clustered heatmaps (group-labelled and sample ID-labelled)
  - Boxplots for top differentially expressed features

---

## Requirements

Install dependencies with pip:

```bash
pip install numpy pandas matplotlib seaborn scipy statsmodels scikit-learn openpyxl
```

---

## Input Files

### 1. Data matrix (CSV)
DIA-NN output in either format:
- **Protein-level**: `report.pg_matrix.csv`
- **Peptide-level**: `report.pr_matrix.csv`

The script auto-detects which type is provided based on the presence of a `Stripped.Sequence` column.

### 2. Group labels file (Excel `.xlsx`)
A spreadsheet with (at minimum) two columns:

| Sample_ID | Group   |
|-----------|---------|
| Sample1   | Case    |
| Sample2   | Case    |
| Sample3   | Control |
| Sample4   | Control |

Column names are configurable — see [Configuration](#configuration) below.

---

## Configuration

Edit the `Config` class at the top of `stats_analysis.py`:

```python
class Config:
    # Paths
    BASE_PATH = r"C:\path\to\your\data"
    DATA_FILE = f"{BASE_PATH}/report.pg_matrix.csv"   # or report.pr_matrix.csv
    LABEL_FILE = f"{BASE_PATH}/group_labels.xlsx"

    # Group names — must match values in your labels file exactly
    CASE_NAME = "Case"
    CONTROL_NAME = "Control"

    # Column names in the labels file
    GROUP_COL = "Group"
    SAMPLE_ID_COL = "Sample_ID"

    # Preprocessing
    MAX_MISSING = 0.30      # Exclude features with > 30% missing values
    IMPUTATION = "minimum"  # Options: 'minimum', 'median', 'mean'
    LOG_TRANSFORM = True    # Apply Log2(x+1) transformation

    # Output
    OUTPUT_DIR = "stats_results"
```

---

## Usage

Run with defaults configured in the script:

```bash
python stats_analysis.py
```

Or override paths from the command line:

```bash
python stats_analysis.py --data path/to/data.csv --labels path/to/labels.xlsx --output results_folder
```

### Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--data` | Path to DIA-NN data matrix (CSV) | `Config.DATA_FILE` |
| `--labels` | Path to group labels file (XLSX) | `Config.LABEL_FILE` |
| `--output` | Output directory | `Config.OUTPUT_DIR` |

---

## Output Files

All outputs are saved to the directory specified by `OUTPUT_DIR` (default: `stats_results/`).

| File | Description |
|------|-------------|
| `statistical_results.csv` | Full results table with p-values, FDR-adjusted p-values, fold changes, and group means |
| `pca_plot.png` | PCA scores plot (PC1 vs PC2) with 95% confidence ellipses |
| `plsda_plot.png` | PLS-DA scores plot with 5-fold CV accuracy and AUC |
| `plsda_vip_plot.png` | Horizontal bar chart of top 20 VIP scores |
| `plsda_vip_scores.csv` | Full VIP scores for all features |
| `volcano_plot.png` | Volcano plot (Log2FC vs -Log10 p-value); significant features highlighted in red |
| `heatmap_grouped.png` | Clustered heatmap with group labels on x-axis |
| `heatmap_sample_ids.png` | Clustered heatmap with individual sample IDs on x-axis |
| `top_features_boxplots.png` | Boxplots for the top 10 most significant features |

### Statistical Results Table Columns

| Column | Description |
|--------|-------------|
| `feature_id` | Row index from the data matrix |
| `Gene` | Gene name (from `Genes` column, if present) |
| `Peptide` | Peptide sequence (peptide-level data only) |
| `Protein` | Protein name (from `Protein.Names`, if present) |
| `mean_case_raw` | Mean intensity in cases (raw, non-transformed) |
| `mean_ctrl_raw` | Mean intensity in controls (raw, non-transformed) |
| `fold_change` | Raw fold change (case / control) |
| `mean_case_log2` | Mean Log2 intensity in cases |
| `mean_ctrl_log2` | Mean Log2 intensity in controls |
| `log2fc` | Log2 fold change |
| `p_value` | Mann-Whitney U p-value |
| `adj_p_value` | FDR-adjusted p-value (Benjamini-Hochberg) |

---

## Analysis Details

### Preprocessing
1. **Missingness filter**: Features with more than `MAX_MISSING` fraction of missing values are removed.
2. **Imputation**: Remaining missing values are replaced using the chosen strategy (`minimum` replaces with the per-feature minimum observed value).
3. **Log2 transformation**: `Log2(x + 1)` is applied to stabilise variance.

### Statistical Testing
- Non-parametric **Mann-Whitney U test** (two-sided) is used for each feature.
- **Fold change** is calculated from raw intensities (before log transformation) to give biologically interpretable values.
- Multiple testing correction uses the **Benjamini-Hochberg FDR** method.

### PLS-DA
- **5-fold stratified cross-validation** is used to evaluate model performance.
- Reports cross-validated **accuracy** and **AUC**.
- **VIP scores** (Variable Importance in Projection) identify the most discriminating features; VIP > 1 is the conventional significance threshold.

### Volcano Plot
- Points are coloured red if they pass **both** `adj_p_value < 0.05` **and** `|Log2FC| > 0.585` (≈1.5-fold change).
- Up to 30 significant features are labelled, with overlap avoidance applied.
