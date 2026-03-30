import os
import glob
import re
import pandas as pd
from scipy import stats
from collections import defaultdict
import numpy as np
import logging
import argparse

# ------------------------------------------------------------
# Hardcoded paths (możesz zmienić)
# ------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
GSP_AGGREGATED_FILE = os.path.join(BASE_DIR, "results", "plots", "gGSP_5.0.csv")
LDDT_GLOBAL_FILE = os.path.join(BASE_DIR, "results", "lddt", "finished", "global_report.csv")
OTHER_METRICS_DIR = os.path.join(BASE_DIR, "results", "other-metrics", "rnapuzzles.github.io")
OUTPUT_DIR = os.path.join(BASE_DIR, "results", "correlations")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ------------------------------------------------------------
# Logging
# ------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(os.path.join(BASE_DIR, 'gsp_correlation.log'), mode='w')
    ]
)
logger = logging.getLogger(__name__)

# ------------------------------------------------------------
# Metric orientation (True = higher is better)
# ------------------------------------------------------------
METRIC_ORIENTATION = {
    "RMSD": False,
    "DI_all": False,
    "INF_all": True,
    "INF_wc": True,
    "INF_nwc": True,
    "INF_stacking": True,
    "Clash_Score": False,
    "P-value": False,
    "mcq": False,
    "TM-score": True,
    "lDDT": True   # now included
}

# ------------------------------------------------------------
# Excluded targets
# ------------------------------------------------------------
EXCLUDED_TARGETS = ['PZ10tRNA', 'PZ26tRNA', 'PZ27tRNA', 'PZ28tRNA']

# ------------------------------------------------------------
# Load aggregated gGSP data
# ------------------------------------------------------------
def load_aggregated_gsp():
    """Loads gGSP data from the aggregated CSV file."""
    if not os.path.exists(GSP_AGGREGATED_FILE):
        raise FileNotFoundError(f"gGSP file not found: {GSP_AGGREGATED_FILE}")
    df = pd.read_csv(GSP_AGGREGATED_FILE)
    # Expected columns: Target, Type, Source, ModelGroup, gGSP
    if 'gGSP' not in df.columns:
        raise ValueError("gGSP column missing in aggregated file")
    df = df[~df['Target'].isin(EXCLUDED_TARGETS)]

    # Extract Lab and Num from ModelGroup (e.g., "TS029_1" -> lab="TS029", num="1")
    def split_modelgroup(mg):
        parts = mg.split('_')
        if len(parts) >= 2:
            lab = '_'.join(parts[:-1])
            num = parts[-1]
        else:
            lab = mg
            num = '0'
        return lab.lower(), num
    df[['Lab', 'Num']] = df['ModelGroup'].apply(lambda x: pd.Series(split_modelgroup(x)))
    logger.info(f"Loaded {len(df)} gGSP records from {GSP_AGGREGATED_FILE}")
    return df

# ------------------------------------------------------------
# Load aggregated lDDT data
# ------------------------------------------------------------
def load_aggregated_lddt():
    """Loads lDDT data from the global_report.csv file."""
    if not os.path.exists(LDDT_GLOBAL_FILE):
        logger.warning(f"lDDT file not found: {LDDT_GLOBAL_FILE}. lDDT will be omitted.")
        return None
    df = pd.read_csv(LDDT_GLOBAL_FILE)
    # Columns: Puzzle, ModelGroup, BestModel, BestSolution, lDDT
    if 'lDDT' not in df.columns:
        raise ValueError("lDDT column missing in global report")
    df.rename(columns={'Puzzle': 'Target'}, inplace=True)
    df = df[~df['Target'].isin(EXCLUDED_TARGETS)]

    # Extract Lab and Num from ModelGroup
    def split_modelgroup(mg):
        parts = mg.split('_')
        if len(parts) >= 2:
            lab = '_'.join(parts[:-1])
            num = parts[-1]
        else:
            lab = mg
            num = '0'
        return lab.lower(), num
    df[['Lab', 'Num']] = df['ModelGroup'].apply(lambda x: pd.Series(split_modelgroup(x)))
    logger.info(f"Loaded {len(df)} lDDT records from {LDDT_GLOBAL_FILE}")
    return df[['Target', 'Lab', 'Num', 'lDDT']]

# ------------------------------------------------------------
# Load other metrics from per‑target directories
# ------------------------------------------------------------
def load_other_metrics(target_dir):
    """
    Loads all metric CSV files (except gGSP/lDDT) from target_dir.
    Returns a dictionary {(lab, num): {metric: value}}.
    """
    metrics_files = glob.glob(os.path.join(target_dir, "*.csv"))
    all_metrics = defaultdict(dict)
    expected_metrics = set(METRIC_ORIENTATION.keys()) - {'gGSP', 'lDDT'}

    for file_path in metrics_files:
        metric_name = os.path.basename(file_path).replace(".csv", "")
        if metric_name not in expected_metrics:
            continue

        try:
            # Try reading with header
            df = pd.read_csv(file_path)
            # Normalize column names
            if 'Lab' not in df.columns or 'Num' not in df.columns:
                # Maybe first two columns are Lab and Num
                if len(df.columns) >= 3:
                    df.columns = ['Lab', 'Num', metric_name]
                else:
                    logger.warning(f"File {file_path} has insufficient columns, skipping")
                    continue
            df['Lab'] = df['Lab'].astype(str).str.lower().str.strip()
            df['Num'] = df['Num'].astype(str).str.strip()
            if metric_name not in df.columns:
                # Use third column as value
                df.rename(columns={df.columns[2]: metric_name}, inplace=True)

            for _, row in df.iterrows():
                try:
                    lab = row['Lab']
                    num = row['Num']
                    val = float(row[metric_name])
                    all_metrics[(lab, num)][metric_name] = val
                except Exception:
                    continue
        except Exception as e:
            logger.error(f"Error reading {file_path}: {e}")
    return all_metrics

# ------------------------------------------------------------
# Merge lDDT into metrics dictionary for a specific target
# ------------------------------------------------------------
def merge_lddt_into_metrics(lddt_df, target, metrics_dict):
    """Adds lDDT values from the aggregated lDDT DataFrame to metrics_dict."""
    if lddt_df is None:
        return metrics_dict
    target_lddt = lddt_df[lddt_df['Target'] == target]
    for _, row in target_lddt.iterrows():
        key = (row['Lab'], row['Num'])
        metrics_dict[key]['lDDT'] = row['lDDT']
    return metrics_dict

# ------------------------------------------------------------
# Display model data table
# ------------------------------------------------------------
def display_model_data(gsp_df, metrics_dict, target_name):
    """Prints a table of models with all available metrics."""
    rows = []
    for _, row in gsp_df.iterrows():
        key = (row['Lab'], row['Num'])
        model_id = f"{key[0]}_{key[1]}"
        base = {'Model': model_id, 'gGSP': row['gGSP']}
        if key in metrics_dict:
            base.update(metrics_dict[key])
        else:
            for metric in METRIC_ORIENTATION:
                if metric not in ('gGSP', 'lDDT'):
                    base[metric] = '---'
        rows.append(base)
    df_display = pd.DataFrame(rows)
    metric_cols = [m for m in METRIC_ORIENTATION if m not in ('gGSP', 'lDDT')]
    # Add lDDT if present
    if any('lDDT' in r for r in rows):
        metric_cols.insert(0, 'lDDT')
    cols = ['Model', 'gGSP'] + metric_cols
    df_display = df_display[cols]
    logger.info(f"\n--- Detailed data for target: {target_name} ---")
    logger.info("\n" + df_display.to_string(index=False))
    total = len(gsp_df)
    logger.info(f"Total models with gGSP: {total}")
    for metric in metric_cols:
        available = sum(1 for r in rows if r.get(metric, '---') != '---')
        logger.info(f"Metric {metric}: available for {available}/{total} models")

# ------------------------------------------------------------
# Calculate correlations
# ------------------------------------------------------------
def calculate_correlations(gsp_df, metrics_dict):
    """Computes Pearson, Spearman, Kendall correlations between gGSP and all other metrics."""
    results = []
    # Build combined data
    combined = []
    for _, row in gsp_df.iterrows():
        key = (row['Lab'], row['Num'])
        if key in metrics_dict:
            entry = {'gGSP': row['gGSP']}
            entry.update(metrics_dict[key])
            combined.append(entry)
    if not combined:
        logger.error("No overlapping models for this target.")
        return pd.DataFrame()
    combined_df = pd.DataFrame(combined)

    for metric in combined_df.columns:
        if metric == 'gGSP':
            continue
        clean = combined_df[['gGSP', metric]].dropna()
        if len(clean) < 3:
            logger.warning(f"Not enough data for {metric} (n={len(clean)}), skipping.")
            continue
        x = clean['gGSP'].values
        y = clean[metric].values
        # Reverse y if metric is "lower is better"
        if not METRIC_ORIENTATION.get(metric, True):
            y = -y
        pearson, p_pearson = stats.pearsonr(x, y)
        spearman, p_spearman = stats.spearmanr(x, y)
        kendall, p_kendall = stats.kendalltau(x, y)
        results.append({
            'metric': metric,
            'orientation': METRIC_ORIENTATION.get(metric, True),
            'pearson_coef': pearson,
            'pearson_pvalue': p_pearson,
            'spearman_coef': spearman,
            'spearman_pvalue': p_spearman,
            'kendall_coef': kendall,
            'kendall_pvalue': p_kendall,
            'n_models': len(clean)
        })
        logger.info(f"{metric}: Pearson={pearson:.3f}, Spearman={spearman:.3f}, Kendall={kendall:.3f} (n={len(clean)})")
    return pd.DataFrame(results)

# ------------------------------------------------------------
# Process a single target
# ------------------------------------------------------------
def process_target(target, gsp_df_target, lddt_df, output_dir):
    """Process one target: load metrics, merge lDDT, display, correlate, save."""
    logger.info(f"Processing target: {target}")

    target_metrics_dir = os.path.join(OTHER_METRICS_DIR, target)
    if not os.path.exists(target_metrics_dir):
        logger.error(f"Metrics directory missing for {target}")
        return None

    # Load per‑target metrics (excluding gGSP/lDDT)
    metrics_dict = load_other_metrics(target_metrics_dir)
    if not metrics_dict:
        logger.error(f"No metrics found for {target}")
        return None

    # Add lDDT from aggregated file
    metrics_dict = merge_lddt_into_metrics(lddt_df, target, metrics_dict)

    # Display data
    display_model_data(gsp_df_target, metrics_dict, target)

    # Calculate correlations
    corr_df = calculate_correlations(gsp_df_target, metrics_dict)
    if corr_df.empty:
        logger.error(f"No correlation results for {target}")
        return None

    # Save
    out_file = os.path.join(output_dir, f"{target}_gGSP_correlations.csv")
    corr_df.to_csv(out_file, index=False)
    logger.info(f"Saved correlation results to {out_file}")
    return corr_df

# ------------------------------------------------------------
# Interactive target selection
# ------------------------------------------------------------
def select_targets(gsp_df):
    """Let user choose which targets to process."""
    targets = sorted(gsp_df['Target'].unique())
    print("\n===== Target selection =====")
    for i, t in enumerate(targets, 1):
        print(f"{i}. {t}")
    print("\nEnter numbers (e.g., 1,3,5-7) or 'all':")
    selection = input("Selection: ").strip()
    if selection.lower() == 'all':
        return targets
    selected = []
    for part in selection.split(','):
        part = part.strip()
        if '-' in part:
            start, end = map(int, part.split('-'))
            for idx in range(start-1, end):
                if 0 <= idx < len(targets):
                    selected.append(targets[idx])
        else:
            try:
                idx = int(part) - 1
                if 0 <= idx < len(targets):
                    selected.append(targets[idx])
            except ValueError:
                continue
    return selected

# ------------------------------------------------------------
# Create summary report (with missing data analysis)
# ------------------------------------------------------------
def create_summary_report(correlation_dir, gsp_df, lddt_df):
    """
    Creates two summary reports:
    1. Simple summary (mean/std/count of coefficients) -> summary_report_simple.csv
    2. Detailed summary (orientation, n_targets, total_models, avg p-values, min/max) -> summary_report_detailed.csv
    Also produces detailed_missing_analysis.txt (same as before).
    """
    logger.info("Creating summary reports...")
    correlation_files = glob.glob(os.path.join(correlation_dir, "*_gGSP_correlations.csv"))
    if not correlation_files:
        logger.error("No correlation files found.")
        return

    all_dfs = []
    for f in correlation_files:
        df = pd.read_csv(f)
        target = os.path.basename(f).replace('_gGSP_correlations.csv', '')
        df['target'] = target
        all_dfs.append(df)
    combined_df = pd.concat(all_dfs, ignore_index=True)

    # ----- 1. Simple summary (mean, std, count) -----
    simple_summary = combined_df.groupby('metric').agg({
        'pearson_coef': ['mean', 'std', 'count'],
        'spearman_coef': ['mean', 'std'],
        'kendall_coef': ['mean', 'std'],
        'n_models': 'sum'
    }).round(4)
    simple_summary.columns = ['_'.join(col).strip() for col in simple_summary.columns.values]
    simple_summary = simple_summary.sort_values('spearman_coef_mean', ascending=False)
    simple_summary.to_csv(os.path.join(correlation_dir, "summary_report_simple.csv"))
    logger.info(f"Simple summary saved to {correlation_dir}/summary_report_simple.csv")

    # ----- 2. Detailed summary (like old version) -----
    detailed_rows = []
    for metric, group in combined_df.groupby('metric'):
        orientation = METRIC_ORIENTATION.get(metric, True)  # True = higher better
        n_targets = group['target'].nunique()
        total_models = group['n_models'].sum()

        # Pearson
        pearson_mean = group['pearson_coef'].mean()
        pearson_std = group['pearson_coef'].std()
        pearson_p_mean = group['pearson_pvalue'].mean()
        pearson_min = group['pearson_coef'].min()
        pearson_max = group['pearson_coef'].max()

        # Spearman
        spearman_mean = group['spearman_coef'].mean()
        spearman_std = group['spearman_coef'].std()
        spearman_p_mean = group['spearman_pvalue'].mean()
        spearman_min = group['spearman_coef'].min()
        spearman_max = group['spearman_coef'].max()

        # Kendall
        kendall_mean = group['kendall_coef'].mean()
        kendall_std = group['kendall_coef'].std()
        kendall_p_mean = group['kendall_pvalue'].mean()
        kendall_min = group['kendall_coef'].min()
        kendall_max = group['kendall_coef'].max()

        detailed_rows.append({
            'metric': metric,
            'orientation': orientation,
            'n_targets': n_targets,
            'total_models': total_models,
            'avg_pearson_coef': pearson_mean,
            'std_pearson_coef': pearson_std,
            'avg_pearson_pvalue': pearson_p_mean,
            'avg_spearman_coef': spearman_mean,
            'std_spearman_coef': spearman_std,
            'avg_spearman_pvalue': spearman_p_mean,
            'avg_kendall_coef': kendall_mean,
            'std_kendall_coef': kendall_std,
            'avg_kendall_pvalue': kendall_p_mean,
            'min_pearson': pearson_min,
            'max_pearson': pearson_max,
            'min_spearman': spearman_min,
            'max_spearman': spearman_max,
            'min_kendall': kendall_min,
            'max_kendall': kendall_max,
        })

    detailed_df = pd.DataFrame(detailed_rows)
    # Round float columns to 4 decimal places
    float_cols = detailed_df.select_dtypes(include=['float64']).columns
    detailed_df[float_cols] = detailed_df[float_cols].round(4)
    detailed_df = detailed_df.sort_values('avg_spearman_coef', ascending=False)
    detailed_df.to_csv(os.path.join(correlation_dir, "summary_report_detailed.csv"), index=False)
    logger.info(f"Detailed summary saved to {correlation_dir}/summary_report_detailed.csv")

    # ----- 3. Detailed missing data analysis (unchanged from original) -----
    logger.info("Creating detailed missing data analysis...")
    all_targets = set(gsp_df['Target'].unique())
    targets_with_data = {os.path.basename(f).replace('_gGSP_correlations.csv', '') for f in correlation_files}
    missing_targets = all_targets - targets_with_data

    metric_missing = {}

    # Process standard metrics (per-target files)
    for target in all_targets:
        # Get gGSP models for this target
        gsp_target = gsp_df[gsp_df['Target'] == target]
        gsp_models = set((row['Lab'], row['Num']) for _, row in gsp_target.iterrows())
        # Metrics directory
        metrics_dir = os.path.join(OTHER_METRICS_DIR, target)
        if not os.path.exists(metrics_dir):
            for metric in METRIC_ORIENTATION:
                if metric not in ('gGSP', 'lDDT'):
                    metric_missing.setdefault(metric, {})[target] = {
                        'missing': list(gsp_models),
                        'total_gsp': len(gsp_models),
                        'status': 'missing_dir'
                    }
            continue

        # Load each metric from file
        for metric in METRIC_ORIENTATION:
            if metric in ('gGSP', 'lDDT'):
                continue
            metric_file = os.path.join(metrics_dir, f"{metric}.csv")
            if not os.path.exists(metric_file):
                metric_missing.setdefault(metric, {})[target] = {
                    'missing': list(gsp_models),
                    'total_gsp': len(gsp_models),
                    'status': 'missing_file'
                }
                continue
            try:
                df_metric = pd.read_csv(metric_file)
                # Normalize
                if 'Lab' not in df_metric.columns or 'Num' not in df_metric.columns:
                    if len(df_metric.columns) >= 3:
                        df_metric.columns = ['Lab', 'Num', metric]
                df_metric['Lab'] = df_metric['Lab'].astype(str).str.lower()
                df_metric['Num'] = df_metric['Num'].astype(str)
                metric_models = set((lab, num) for lab, num in zip(df_metric['Lab'], df_metric['Num']))
                missing = gsp_models - metric_models
                if missing:
                    metric_missing.setdefault(metric, {})[target] = {
                        'missing': list(missing),
                        'total_gsp': len(gsp_models),
                        'total_metric': len(metric_models),
                        'status': 'partial'
                    }
                else:
                    metric_missing.setdefault(metric, {})[target] = {
                        'missing': [],
                        'total_gsp': len(gsp_models),
                        'total_metric': len(metric_models),
                        'status': 'complete'
                    }
            except Exception as e:
                metric_missing.setdefault(metric, {})[target] = {
                    'missing': list(gsp_models),
                    'total_gsp': len(gsp_models),
                    'status': f'error_reading: {e}'
                }

    # Handle lDDT separately (using aggregated lddt_df)
    if lddt_df is not None:
        for target in all_targets:
            gsp_target = gsp_df[gsp_df['Target'] == target]
            gsp_models = set((row['Lab'], row['Num']) for _, row in gsp_target.iterrows())
            lddt_target = lddt_df[lddt_df['Target'] == target]
            lddt_models = set((row['Lab'], row['Num']) for _, row in lddt_target.iterrows())
            missing = gsp_models - lddt_models
            if missing:
                metric_missing.setdefault('lDDT', {})[target] = {
                    'missing': list(missing),
                    'total_gsp': len(gsp_models),
                    'total_metric': len(lddt_models),
                    'status': 'partial'
                }
            else:
                metric_missing.setdefault('lDDT', {})[target] = {
                    'missing': [],
                    'total_gsp': len(gsp_models),
                    'total_metric': len(lddt_models),
                    'status': 'complete'
                }
    else:
        for target in all_targets:
            gsp_target = gsp_df[gsp_df['Target'] == target]
            gsp_models = set((row['Lab'], row['Num']) for _, row in gsp_target.iterrows())
            metric_missing.setdefault('lDDT', {})[target] = {
                'missing': list(gsp_models),
                'total_gsp': len(gsp_models),
                'status': 'no_lddt_file'
            }

    # Write missing data report
    missing_report_path = os.path.join(correlation_dir, "detailed_missing_analysis.txt")
    with open(missing_report_path, 'w') as f:
        f.write("="*120 + "\n")
        f.write("DETAILED MISSING DATA ANALYSIS\n")
        f.write("="*120 + "\n\n")
        f.write(f"Total targets: {len(all_targets)}\n")
        f.write(f"Targets with correlation results: {len(targets_with_data)}\n")
        f.write(f"Missing targets: {sorted(missing_targets)}\n\n")

        for metric in sorted(metric_missing.keys()):
            f.write(f"\n{'='*80}\n")
            f.write(f"METRIC: {metric}\n")
            f.write(f"{'='*80}\n")
            data = metric_missing[metric]
            if not data:
                f.write("No data for this metric.\n")
                continue
            for target, info in sorted(data.items()):
                f.write(f"\n{target}:\n")
                f.write(f"  gGSP models: {info.get('total_gsp', 'N/A')}\n")
                if info.get('status') == 'complete':
                    f.write("  Status: COMPLETE\n")
                elif info.get('status') == 'missing_file':
                    f.write("  Status: MISSING FILE\n")
                elif info.get('status') == 'missing_dir':
                    f.write("  Status: METRICS DIRECTORY MISSING\n")
                elif info.get('status') == 'no_lddt_file':
                    f.write("  Status: lDDT GLOBAL FILE NOT LOADED\n")
                elif info.get('status') == 'partial':
                    f.write(f"  Status: PARTIAL (models with metric: {info['total_metric']}/{info['total_gsp']})\n")
                    missing = info['missing']
                    if missing:
                        f.write("  Missing models:\n")
                        from collections import defaultdict
                        by_lab = defaultdict(list)
                        for lab, num in missing:
                            by_lab[lab].append(num)
                        for lab, nums in sorted(by_lab.items()):
                            f.write(f"    {lab}: {', '.join(sorted(nums, key=lambda x: int(x) if x.isdigit() else x))}\n")
                else:
                    f.write(f"  Status: ERROR ({info.get('status')})\n")
    logger.info(f"Missing data report saved to {missing_report_path}")

    return detailed_df

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
def main():
    # Load aggregated data
    gsp_df = load_aggregated_gsp()
    lddt_df = load_aggregated_lddt()

    # Select targets
    targets = select_targets(gsp_df)
    if not targets:
        logger.error("No targets selected.")
        return

    # Process each target
    for target in targets:
        gsp_target = gsp_df[gsp_df['Target'] == target]
        process_target(target, gsp_target, lddt_df, OUTPUT_DIR)

    # Ask for summary report
    print("\nCreate summary report? (y/n): ", end="")
    if input().strip().lower() == 'y':
        create_summary_report(OUTPUT_DIR, gsp_df, lddt_df)

    print("\nDone.")

if __name__ == "__main__":
    main()