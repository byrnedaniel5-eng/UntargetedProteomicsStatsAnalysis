# Untargeted Proteomics Stats Analysis
One-stop data analysis pipeline for simple analysis of proteomics data - filtering of proteins for missingness, log transformation etc through to PCA / PLS-DA plots, boxplots for most statistically significant proteins, all in one place!

**REQUIRED INPUT FILES**

  DIA-NN output: report.pg_matrix.csv or report.pr_matrix.csv
  Group labels file (.xlsx): Two columns, one with group (case or control) and one with sample ID (e.g. Case1, Case2, Case3 etc)

**CONFIGURATION IN SCRIPT**

  BASE_PATH = r"path_to_data_folder"
  DATA_FILE = f{BASE_PATH}/"report.pg_or_pr_matrix.csv"
  LABEL_FILE = f{BASE_PATH}/"group_labels.xlsx"

  CASE_NAME = name matching case in group file
  CONTROL_NAME = name matching control in group file
  GROUP_COL = name of column in group file containing groups
  SAMPLE_ID_COL = name of column in group file containing sample IDs

  MAX_MISSING = maximum missingness (0.3 default, 30%)
  IMPUTATION = imputation method to be used for missingness below threshold (minimum, median, mean)
  LOG_TRANSFORM = whether to log transform data (True, False)

  OUTPUT_DIR = name of output directory ("stats_results" is default+)
