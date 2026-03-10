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
# Path settings
# ------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_GSP_DIR = os.path.join(BASE_DIR, "results", "gsp", "finished")
BASE_METRICS_DIR = os.path.join(BASE_DIR, "results", "other-metrics", "rnapuzzles.github.io")
OUTPUT_DIR = os.path.join(BASE_DIR, "results", "correlations")

# ------------------------------------------------------------
# Logging configuration
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
# Metric orientation dictionary (TRUE = higher value is better)
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
    "gGSP": True
}

def parse_gsp_files(gsp_dir):
    """
    Robust reader for *C1'-gGSP.csv files.
    Handles CR targets (TS/TSR variants), PZ19 variants (with or without 'rpr'),
    and default pattern for other targets.
    Returns a DataFrame with columns: Lab (lowercase), Num (string), gGSP (float).
    """
    gsp_files = glob.glob(os.path.join(gsp_dir, "*_C1'-gGSP.csv"))
    gsp_data = []
    logger.info(f"Found {len(gsp_files)} gGSP files in {gsp_dir}")

    target_name = os.path.basename(gsp_dir)
    target_digits = ''.join(re.findall(r'\d+', target_name))  # e.g. "19" or "1107"

    for file_path in gsp_files:
        try:
            filename = os.path.basename(file_path)
            logger.debug(f"Processing gGSP file: {filename}")

            lab = None
            num = None

            # ---------- Special rules for CR targets ----------
            if target_name.startswith("CR"):
                # Try to extract core between _solution_ and _C1'-gGSP.csv
                m_core = re.search(r'_solution_(.*)_C1\'-gGSP\.csv$', filename)
                if m_core:
                    core = m_core.group(1)
                    logger.debug(f"CR core extracted: '{core}'")
                else:
                    # fallback: everything before _C1'-gGSP.csv minus leading target if present
                    alt = re.sub(r"_C1\'-gGSP\.csv$", "", filename)
                    if "solution_" in alt:
                        core = alt.split("solution_", 1)[-1]
                        logger.debug(f"CR alt core after 'solution_': '{core}'")
                    else:
                        core = alt
                        logger.debug(f"CR alt core used: '{core}'")

                # lab: prefer TSR/TS pattern
                lab_match = re.search(r'(TSR?\d+)', core, flags=re.IGNORECASE)
                if not lab_match:
                    lab_match = re.search(r'(TSR?\d+)', filename, flags=re.IGNORECASE)
                if lab_match:
                    lab = lab_match.group(1)
                    logger.debug(f"CR lab found: {lab}")
                else:
                    logger.warning(f"TS/TSR not found in CR filename: {filename}")

                # num: try specific patterns, then fallbacks
                num_match = re.search(r'_(\d+)(?:_refined|_rpr)?_C1', filename)
                if not num_match:
                    num_match = re.search(r'_(\d+)(?:_refined|_rpr)?$', core)
                if num_match:
                    num = num_match.group(1)
                    logger.debug(f"CR num found by pattern: {num}")
                else:
                    nums = re.findall(r'\d+', core)
                    nums_filtered = [n for n in nums if n != target_digits]
                    if nums_filtered:
                        num = nums_filtered[-1]
                        logger.debug(f"CR num fallback from core: {num}")
                    else:
                        nums_all = re.findall(r'\d+', filename)
                        nums_all_filtered = [n for n in nums_all if n != target_digits]
                        if nums_all_filtered:
                            num = nums_all_filtered[-1]
                            logger.debug(f"CR num fallback from filename: {num}")

                if not lab or not num:
                    logger.warning(f"Skipping CR file (couldn't extract lab/num): {filename}")
                    continue

            # ---------- Flexible rules for PZ19 ----------
            elif target_name.upper() in ("PZ19", "19"):
                # Look for '_solution_0_' or '_solution_'
                m = re.search(r'_solution_0?_(.*)_C1\'-gGSP\.csv$', filename)
                if m:
                    rest = m.group(1)
                    logger.debug(f"PZ19 rest extracted (via solution): '{rest}'")
                else:
                    # fallback: remove suffix and try to take part after first '_'
                    alt = re.sub(r"_C1\'-gGSP\.csv$", "", filename)
                    # Try to remove leading '19_' if present
                    if alt.startswith("19_solution_"):
                        rest = alt.split("19_solution_", 1)[-1]
                        logger.debug(f"PZ19 alt rest after '19_solution_': '{rest}'")
                    else:
                        # try take substring after first 'solution_' if any
                        if "solution_" in alt:
                            rest = alt.split("solution_", 1)[-1]
                            logger.debug(f"PZ19 alt rest after 'solution_': '{rest}'")
                        else:
                            # last resort: drop leading "19_" or "PZ19_"
                            rest = re.sub(r'^(?:19_|PZ19_)', '', alt)
                            logger.debug(f"PZ19 fallback rest: '{rest}'")

                # normalize rest: remove leading '19_' or 'PZ19_'
                rest = re.sub(r'^(?:19_|PZ19_)', '', rest, flags=re.IGNORECASE)

                # split by underscores to get tokens
                parts = [p for p in rest.split('_') if p != '']
                if not parts:
                    logger.warning(f"PZ19: couldn't parse tokens from rest '{rest}' in file {filename}")
                    continue

                # lab is first non-numeric token (or first token if all are non-numeric)
                lab_candidate = None
                for p in parts:
                    if not p.isdigit():
                        lab_candidate = p
                        break
                if lab_candidate is None:
                    lab_candidate = parts[0]
                lab = lab_candidate

                # num is the last numeric token in rest (or last token if numeric)
                nums = re.findall(r'\d+', rest)
                if nums:
                    num = nums[-1]
                else:
                    # sometimes last part is numeric-like but with letters — try last token
                    last_tok = parts[-1]
                    num_match = re.search(r'(\d+)', last_tok)
                    if num_match:
                        num = num_match.group(1)
                    else:
                        logger.warning(f"PZ19: no numeric model id found in '{rest}' for file {filename}")
                        continue

                logger.debug(f"PZ19 parsed lab='{lab}', num='{num}' from '{rest}'")

            # ---------- Default method for other targets ----------
            else:
                pattern = re.compile(r'_([^_]+)_(\d+)(?:_rpr)?(?:_refined)?_C1\'-gGSP\.csv$')
                match = pattern.search(filename)
                if not match:
                    logger.warning(f"Failed to match default pattern for: {filename}")
                    continue
                lab = match.group(1)
                num = match.group(2)

            # Final normalization
            lab = str(lab).lower()
            num = str(num)

            # Read gGSP numeric value from second line
            with open(file_path, 'r', encoding='utf-8') as f:
                _header = f.readline().strip()
                value_line = f.readline().strip()
                try:
                    value = float(value_line)
                except Exception:
                    logger.error(f"Cannot parse numeric gGSP in file {file_path}: '{value_line}'")
                    continue

            gsp_data.append({"Lab": lab, "Num": num, "gGSP": value})
            logger.debug(f"Appended gGSP: {lab}_{num} = {value}")

        except Exception as e:
            logger.error(f"Error processing {file_path}: {e}", exc_info=True)

    return pd.DataFrame(gsp_data)


def load_other_metrics(metrics_dir):
    """
    Loads all metric CSV files (excluding gGSP) from metrics_dir.
    Returns a dictionary {(lab, num): {metric: value}} with lab converted to lowercase.
    """
    metrics_files = glob.glob(os.path.join(metrics_dir, "*.csv"))
    all_metrics = defaultdict(dict)
    logger.info(f"Found {len(metrics_files)} metric files in {metrics_dir}")

    expected_metrics = list(METRIC_ORIENTATION.keys())

    for file_path in metrics_files:
        metric_name = os.path.basename(file_path).replace(".csv", "")

        if metric_name not in expected_metrics or metric_name == "gGSP":
            logger.debug(f"Skipping file: {os.path.basename(file_path)}")
            continue

        try:
            logger.debug(f"Processing metric file: {os.path.basename(file_path)}")

            # Attempt to read CSV
            try:
                df = pd.read_csv(file_path)
            except pd.errors.ParserError:
                logger.warning(f"Parser error for {file_path}, trying without header")
                df = pd.read_csv(file_path, header=None)
                if len(df.columns) >= 3:
                    df.columns = ["Lab", "Num", metric_name]

            if len(df.columns) < 3:
                logger.error(f"File {file_path} has only {len(df.columns)} columns, skipping")
                continue

            if metric_name not in df.columns:
                value_col = df.columns[2]
                df = df.rename(columns={value_col: metric_name})

            if "Lab" not in df.columns or "Num" not in df.columns:
                logger.error(f"Missing required columns in {file_path}, skipping")
                continue

            for _, row in df.iterrows():
                try:
                    # Convert lab to lowercase to match gGSP files
                    lab = str(row["Lab"]).lower()
                    num = str(row["Num"])
                    value = row[metric_name]

                    if pd.isna(value):
                        continue

                    try:
                        numeric_value = float(value)
                    except (ValueError, TypeError):
                        logger.warning(f"Invalid value in {metric_name} for {lab}_{num}: {value}")
                        continue

                    all_metrics[(lab, num)][metric_name] = numeric_value
                    logger.debug(f"Loaded {metric_name} for {lab}_{num} = {numeric_value}")

                except Exception as e:
                    logger.error(f"Error processing row in {file_path}: {e}")
                    continue

        except Exception as e:
            logger.error(f"Error processing file {file_path}: {str(e)}", exc_info=True)

    return all_metrics


def display_model_data(gsp_df, other_metrics, target_name):
    """
    Displays a table with models and available metrics.
    """
    rows = []
    for _, row in gsp_df.iterrows():
        key = (row["Lab"], str(row["Num"]))
        model_id = f"{key[0]}_{key[1]}"
        gsp_val = row["gGSP"]
        base = {"Model": model_id, "gGSP": gsp_val}

        if key in other_metrics:
            metrics_dict = other_metrics[key]
            for metric in METRIC_ORIENTATION:
                if metric == "gGSP":
                    continue
                base[metric] = metrics_dict.get(metric, "---")
        else:
            for metric in METRIC_ORIENTATION:
                if metric == "gGSP":
                    continue
                base[metric] = "---"
        rows.append(base)

    df_display = pd.DataFrame(rows)
    metric_cols = [m for m in METRIC_ORIENTATION if m != "gGSP"]
    cols = ["Model", "gGSP"] + metric_cols
    df_display = df_display[cols]

    logger.info(f"\n--- Detailed data for target: {target_name} ---")
    logger.info("\n" + df_display.to_string(index=False))

    total_models = len(gsp_df)
    logger.info(f"Total number of models with gGSP: {total_models}")
    for metric in metric_cols:
        available = sum(1 for r in rows if r[metric] != "---")
        logger.info(f"Metric {metric}: available for {available}/{total_models} models")


def calculate_correlations(gsp_df, other_metrics):
    """
    Calculates correlations between gGSP and other metrics.
    """
    results = []
    logger.info("Starting correlation calculations...")

    combined_data = []
    for _, row in gsp_df.iterrows():
        key = (row["Lab"], str(row["Num"]))
        if key in other_metrics:
            entry = {"Lab": key[0], "Num": key[1], "gGSP": row["gGSP"]}
            entry.update(other_metrics[key])
            combined_data.append(entry)
        else:
            logger.warning(f"No metrics for model: {key}")

    if not combined_data:
        logger.error("No data to calculate correlations!")
        return pd.DataFrame()

    combined_df = pd.DataFrame(combined_data)
    logger.info(f"Combined set contains {len(combined_df)} models")

    for metric in combined_df.columns:
        if metric in ["Lab", "Num"] or metric == "gGSP":
            continue

        clean_df = combined_df[["gGSP", metric]].copy()
        clean_df["gGSP"] = pd.to_numeric(clean_df["gGSP"], errors='coerce')
        clean_df[metric] = pd.to_numeric(clean_df[metric], errors='coerce')
        clean_df = clean_df.dropna()

        original_count = len(combined_df)
        clean_count = len(clean_df)
        if clean_count < original_count:
            logger.warning(f"Removed {original_count - clean_count} rows with invalid data for {metric}")

        if len(clean_df) < 3:
            logger.warning(f"Insufficient data for {metric} (only {len(clean_df)} pairs). Skipping...")
            continue

        try:
            x = clean_df["gGSP"].to_numpy()
            y = clean_df[metric].to_numpy()

            if not METRIC_ORIENTATION.get(metric, True):
                y = -y
                logger.debug(f"Reversed values for {metric} (lower is better)")

            logger.info(f"\n=== Data for metric: {metric} ===")
            for i, (xi, yi) in enumerate(zip(x, y), 1):
                logger.info(f"Model {i}: gGSP={xi:.4f}, {metric}={yi:.4f}")

            pearson_corr, pearson_p = stats.pearsonr(x, y)
            spearman_corr, spearman_p = stats.spearmanr(x, y)
            kendall_corr, kendall_p = stats.kendalltau(x, y)

            results.append({
                "metric": metric,
                "orientation": METRIC_ORIENTATION.get(metric, "unknown"),
                "pearson_coef": pearson_corr,
                "pearson_pvalue": pearson_p,
                "spearman_coef": spearman_corr,
                "spearman_pvalue": spearman_p,
                "kendall_coef": kendall_corr,
                "kendall_pvalue": kendall_p,
                "n_models": len(clean_df)
            })

            logger.info(f"[OK] Correlations for {metric}: Pearson={pearson_corr:.3f}, Spearman={spearman_corr:.3f}, Kendall={kendall_corr:.3f}")

        except Exception as e:
            logger.error(f"Error calculating correlations for {metric}: {str(e)}", exc_info=True)

    return pd.DataFrame(results)


def process_target(target_name, base_gsp_dir, base_metrics_dir, output_dir):
    """
    Processes a single target.
    """
    logger.info(f"Processing target: {target_name}")

    gsp_dir = os.path.join(base_gsp_dir, target_name)
    metrics_dir = os.path.join(base_metrics_dir, target_name)

    if not os.path.exists(gsp_dir):
        logger.error(f"GSP directory does not exist: {gsp_dir}")
        return None
    if not os.path.exists(metrics_dir):
        logger.error(f"Metrics directory does not exist: {metrics_dir}")
        return None

    gsp_df = parse_gsp_files(gsp_dir)
    if gsp_df.empty:
        logger.error(f"No gGSP data for {target_name}")
        return None

    other_metrics = load_other_metrics(metrics_dir)
    if not other_metrics:
        logger.error(f"No metrics for {target_name}")
        return None

    display_model_data(gsp_df, other_metrics, target_name)

    correlations_df = calculate_correlations(gsp_df, other_metrics)
    if correlations_df.empty:
        logger.error(f"No correlation results for {target_name}")
        return None

    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{target_name}_correlations.csv")
    correlations_df.to_csv(output_file, index=False)
    logger.info(f"Correlation results saved to: {output_file}")

    return correlations_df


def select_targets(base_dir):
    """
    Interactive target selection.
    """
    print("\n===== Target Selection =====")
    print("Please select datasets to process:")

    available_targets = [d for d in os.listdir(base_dir)
                         if os.path.isdir(os.path.join(base_dir, d))]

    if not available_targets:
        print("No targets found!")
        return []

    print("Available targets:")
    for i, target in enumerate(available_targets, 1):
        print(f"{i}. {target}")

    print("\nEnter target numbers (e.g., 1,3,5-7) or 'all'")
    selection = input("Selection: ").strip()

    if selection.lower() == 'all':
        return available_targets

    selected_targets = []
    for part in selection.split(','):
        part = part.strip()
        if '-' in part:
            start, end = map(int, part.split('-'))
            selected_targets.extend(available_targets[i-1] for i in range(start, end+1))
        else:
            try:
                index = int(part)
                if 1 <= index <= len(available_targets):
                    selected_targets.append(available_targets[index-1])
            except ValueError:
                continue

    return list(set(selected_targets))


def create_summary_report(correlation_dir, output_file=None):
    """
    Creates an overall summary from all CSV files in correlation_dir.
    """
    logger.info(f"Creating summary report from directory: {correlation_dir}")

    if not os.path.exists(correlation_dir):
        logger.error(f"Directory does not exist: {correlation_dir}")
        return None

    correlation_files = glob.glob(os.path.join(correlation_dir, "*_correlations.csv"))

    if not correlation_files:
        logger.error(f"No *_correlations.csv files found in {correlation_dir}")
        return None

    logger.info(f"Found {len(correlation_files)} correlation files")

    all_dataframes = []
    all_target_names = []

    for file_path in correlation_files:
        try:
            df = pd.read_csv(file_path)
            target_name = os.path.basename(file_path).replace('_correlations.csv', '')
            df['source_file'] = target_name
            all_dataframes.append(df)
            all_target_names.append(target_name)
            logger.debug(f"Loaded {os.path.basename(file_path)} with {len(df)} rows")
        except Exception as e:
            logger.error(f"Error loading {file_path}: {str(e)}")

    if not all_dataframes:
        logger.error("No data loaded from correlation files")
        return None

    combined_df = pd.concat(all_dataframes, ignore_index=True)
    logger.info(f"Total number of rows: {len(combined_df)}")

    required_columns = ['metric', 'orientation', 'pearson_coef', 'pearson_pvalue',
                        'spearman_coef', 'spearman_pvalue', 'kendall_coef',
                        'kendall_pvalue', 'n_models']
    missing_columns = [col for col in required_columns if col not in combined_df.columns]
    if missing_columns:
        logger.error(f"Missing columns: {missing_columns}")
        return None

    # Detailed missing data analysis (using global paths)
    detailed_missing_data = {}
    detailed_stats = {}

    for target in all_target_names:
        logger.info(f"Analyzing missing details for target: {target}")

        gsp_dir = os.path.join(BASE_GSP_DIR, target)
        metrics_dir = os.path.join(BASE_METRICS_DIR, target)

        if not os.path.exists(gsp_dir):
            logger.warning(f"Missing gGSP directory for target {target}")
            detailed_stats[target] = {"gsp_models": 0, "has_gsp": False}
            continue

        if not os.path.exists(metrics_dir):
            logger.warning(f"Missing metrics directory for target {target}")
            detailed_stats[target] = {"gsp_models": 0, "has_metrics": False}
            continue

        try:
            gsp_df = parse_gsp_files(gsp_dir)
            if gsp_df.empty:
                logger.warning(f"No gGSP data for target {target}")
                detailed_stats[target] = {"gsp_models": 0, "gsp_data": []}
                continue

            gsp_models = [(row["Lab"], str(row["Num"])) for _, row in gsp_df.iterrows()]
            detailed_stats[target] = {
                "gsp_models": len(gsp_models),
                "gsp_data": gsp_models,
                "has_gsp": True
            }
        except Exception as e:
            logger.error(f"Error loading gGSP for {target}: {str(e)}")
            detailed_stats[target] = {"gsp_models": 0, "has_gsp": False}
            continue

        for metric in METRIC_ORIENTATION.keys():
            if metric == "gGSP":
                continue

            metric_file = os.path.join(metrics_dir, f"{metric}.csv")
            if not os.path.exists(metric_file):
                if metric not in detailed_missing_data:
                    detailed_missing_data[metric] = {}
                detailed_missing_data[metric][target] = {
                    "missing": gsp_models,
                    "status": "missing_file"
                }
                continue

            try:
                metric_df = pd.read_csv(metric_file)
                metric_models = [(str(row["Lab"]).lower(), str(row["Num"])) for _, row in metric_df.iterrows()]
            except Exception as e:
                logger.error(f"Error loading {metric_file}: {str(e)}")
                if metric not in detailed_missing_data:
                    detailed_missing_data[metric] = {}
                detailed_missing_data[metric][target] = {
                    "missing": gsp_models,
                    "status": "error_reading"
                }
                continue

            missing_models = [m for m in gsp_models if m not in metric_models]
            if missing_models:
                if metric not in detailed_missing_data:
                    detailed_missing_data[metric] = {}
                detailed_missing_data[metric][target] = {
                    "missing": missing_models,
                    "total_gsp": len(gsp_models),
                    "total_metric": len(metric_models),
                    "status": "partial_data"
                }
            else:
                if metric not in detailed_missing_data:
                    detailed_missing_data[metric] = {}
                detailed_missing_data[metric][target] = {
                    "missing": [],
                    "total_gsp": len(gsp_models),
                    "total_metric": len(metric_models),
                    "status": "complete"
                }

    # Create statistical summary
    summary_data = []
    grouped = combined_df.groupby(['metric', 'orientation'])

    for (metric, orientation), group in grouped:
        n_targets = group['source_file'].nunique()
        total_models = group['n_models'].sum()

        summary_data.append({
            'metric': metric,
            'orientation': orientation,
            'n_targets': n_targets,
            'total_models': total_models,
            'avg_pearson_coef': group['pearson_coef'].mean(),
            'std_pearson_coef': group['pearson_coef'].std(),
            'avg_pearson_pvalue': group['pearson_pvalue'].mean(),
            'avg_spearman_coef': group['spearman_coef'].mean(),
            'std_spearman_coef': group['spearman_coef'].std(),
            'avg_spearman_pvalue': group['spearman_pvalue'].mean(),
            'avg_kendall_coef': group['kendall_coef'].mean(),
            'std_kendall_coef': group['kendall_coef'].std(),
            'avg_kendall_pvalue': group['kendall_pvalue'].mean(),
            'min_pearson': group['pearson_coef'].min(),
            'max_pearson': group['pearson_coef'].max(),
            'min_spearman': group['spearman_coef'].min(),
            'max_spearman': group['spearman_coef'].max(),
            'min_kendall': group['kendall_coef'].min(),
            'max_kendall': group['kendall_coef'].max()
        })

    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values('avg_spearman_coef', ascending=False)

    if output_file is None:
        output_file = os.path.join(correlation_dir, "summary_report.csv")
    summary_df.to_csv(output_file, index=False, float_format='%.6f')

    # Rounded version
    summary_rounded = summary_df.copy()
    numeric_cols = summary_rounded.select_dtypes(include=[np.number]).columns
    summary_rounded[numeric_cols] = summary_rounded[numeric_cols].round(4)
    summary_rounded_file = output_file.replace('.csv', '_rounded.csv')
    summary_rounded.to_csv(summary_rounded_file, index=False)

    # Generate detailed missing report
    report_lines = []
    report_lines.append("="*120)
    report_lines.append(f"OVERALL DATA ANALYSIS FOR {len(all_target_names)} TARGETS")
    report_lines.append("="*120)

    total_gsp_models = sum(stats.get("gsp_models", 0) for stats in detailed_stats.values())
    targets_with_gsp = sum(1 for stats in detailed_stats.values() if stats.get("has_gsp", False))

    report_lines.append(f"\n📊 OVERALL STATISTICS:")
    report_lines.append(f"   Targets with gGSP data: {targets_with_gsp}/{len(all_target_names)}")
    report_lines.append(f"   Total number of gGSP models: {total_gsp_models}")

    for metric in METRIC_ORIENTATION.keys():
        if metric == "gGSP":
            continue

        report_lines.append(f"\n{'='*80}")
        report_lines.append(f"ANALYSIS FOR METRIC: {metric}")
        report_lines.append(f"{'='*80}")

        if metric not in detailed_missing_data:
            report_lines.append(f"❌ NO DATA FOR THIS METRIC IN ANY TARGET")
            continue

        missing_info = detailed_missing_data[metric]

        targets_with_metric = 0
        targets_missing_file = 0
        targets_partial_data = 0
        total_missing_models = 0

        for target, info in missing_info.items():
            if info["status"] == "complete":
                targets_with_metric += 1
            elif info["status"] == "missing_file":
                targets_missing_file += 1
                total_missing_models += len(info["missing"])
            elif info["status"] == "partial_data":
                targets_partial_data += 1
                total_missing_models += len(info["missing"])

        total_targets_for_metric = targets_with_metric + targets_partial_data + targets_missing_file

        report_lines.append(f"\n📈 SUMMARY:")
        report_lines.append(f"   Targets with complete data: {targets_with_metric}/{total_targets_for_metric}")
        report_lines.append(f"   Targets with partial data: {targets_partial_data}/{total_targets_for_metric}")
        report_lines.append(f"   Targets missing file: {targets_missing_file}/{total_targets_for_metric}")
        report_lines.append(f"   Total number of missing models: {total_missing_models}")

        if targets_missing_file > 0:
            report_lines.append(f"\n📁 TARGETS MISSING FILE {metric}.csv:")
            missing_file_targets = [t for t, info in missing_info.items() if info["status"] == "missing_file"]
            for i, target in enumerate(missing_file_targets[:10], 1):
                gsp_count = detailed_stats.get(target, {}).get("gsp_models", 0)
                report_lines.append(f"   {i}. {target} (missing file, {gsp_count} gGSP models)")
            if len(missing_file_targets) > 10:
                report_lines.append(f"   ... and {len(missing_file_targets) - 10} more")

        if targets_partial_data > 0:
            report_lines.append(f"\n⚠️  TARGETS WITH PARTIAL DATA:")
            partial_targets = [(t, info) for t, info in missing_info.items() if info["status"] == "partial_data"]
            partial_targets.sort(key=lambda x: len(x[1]["missing"]), reverse=True)

            for i, (target, info) in enumerate(partial_targets[:5], 1):
                missing_count = len(info["missing"])
                total_gsp = info.get("total_gsp", 0)
                total_metric = info.get("total_metric", 0)

                report_lines.append(f"\n   {i}. {target}:")
                report_lines.append(f"       gGSP models: {total_gsp}")
                report_lines.append(f"       Models with {metric}: {total_metric}")
                report_lines.append(f"       Missing models: {missing_count} ({missing_count/total_gsp*100:.1f}%)")

                models_by_lab = {}
                for lab, num in info["missing"]:
                    if lab not in models_by_lab:
                        models_by_lab[lab] = []
                    models_by_lab[lab].append(num)

                for j, (lab, nums) in enumerate(sorted(models_by_lab.items())[:3]):
                    nums_str = ", ".join(sorted(nums, key=lambda x: int(x) if x.isdigit() else x))
                    report_lines.append(f"       {lab}: {nums_str}")
                if len(models_by_lab) > 3:
                    report_lines.append(f"       ... and {len(models_by_lab) - 3} other labs")

            if len(partial_targets) > 5:
                report_lines.append(f"\n   ... and {len(partial_targets) - 5} more targets with partial data")

    # Add summary from correlations.csv
    report_lines.append(f"\n{'='*120}")
    report_lines.append("SUMMARY FROM correlations.csv")
    report_lines.append(f"{'='*120}")

    for metric in METRIC_ORIENTATION.keys():
        if metric == "gGSP":
            continue
        metric_data = combined_df[combined_df['metric'] == metric]
        if not metric_data.empty:
            n_targets = metric_data['source_file'].nunique()
            avg_spearman = metric_data['spearman_coef'].mean()
            avg_pearson = metric_data['pearson_coef'].mean()
            status = "✅" if n_targets == len(all_target_names) else "⚠️"
            report_lines.append(f"{status} {metric}: {n_targets} targets, Spearman={avg_spearman:.3f}, Pearson={avg_pearson:.3f}")

    report_file_path = os.path.join(correlation_dir, "detailed_missing_analysis.txt")
    with open(report_file_path, 'w', encoding='utf-8') as f:
        f.write("\n".join(report_lines))

    print("="*120)
    print("SUMMARY REPORT - DATA MISSING ANALYSIS")
    print("="*120)
    print(f"Total number of targets: {len(all_target_names)}")
    print(f"Targets with gGSP data: {targets_with_gsp}")
    print(f"Total number of gGSP models: {total_gsp_models}")
    print("\n" + "="*120)
    print("METRICS STATUS (from correlations.csv):")
    for metric in METRIC_ORIENTATION.keys():
        if metric == "gGSP":
            continue
        metric_data = combined_df[combined_df['metric'] == metric]
        if not metric_data.empty:
            n_targets = metric_data['source_file'].nunique()
            avg_spearman = metric_data['spearman_coef'].mean()
            status = "COMPLETE" if n_targets == len(all_target_names) else f"MISSING {len(all_target_names) - n_targets}"
            print(f"  {metric}: {n_targets} targets, Spearman={avg_spearman:.3f} - {status}")

    print(f"\n📋 Full report saved to: {report_file_path}")
    logger.info(f"Summary report saved as: {output_file}")
    logger.info(f"Rounded report saved as: {summary_rounded_file}")
    logger.info(f"Detailed missing report saved as: {report_file_path}")

    return summary_df


def run_correlation_analysis():
    """
    Main function to run correlation analysis.
    """
    os.makedirs(BASE_GSP_DIR, exist_ok=True)
    os.makedirs(BASE_METRICS_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    logger.info(f"Main GSP directory: {BASE_GSP_DIR}")
    logger.info(f"Main metrics directory: {BASE_METRICS_DIR}")

    selected_targets = select_targets(BASE_GSP_DIR)

    if not selected_targets:
        logger.error("No targets selected!")
        return

    logger.info(f"Selected targets: {', '.join(selected_targets)}")

    all_results = {}
    for target in selected_targets:
        result = process_target(target, BASE_GSP_DIR, BASE_METRICS_DIR, OUTPUT_DIR)
        if result is not None:
            all_results[target] = result
            print(f"\nCorrelation summary for {target}:")
            print(result[["metric", "pearson_coef", "spearman_coef", "kendall_coef", "n_models"]])

    logger.info("Processing completed")
    print(f"\nAll results saved in: {OUTPUT_DIR}")
    print(f"Processed {len(all_results)} targets.")

    print("\n" + "="*80)
    create_summary = input("Do you want to create a summary report? (y/n): ").strip().lower()
    if create_summary in ('y', 'yes'):
        print(f"\nCreating summary report from: {OUTPUT_DIR}")
        summary = create_summary_report(OUTPUT_DIR)
        if summary is not None:
            print(f"\nSummary report created successfully!")

    print("\nDone!")


def main():
    parser = argparse.ArgumentParser(description='Correlation analysis between gGSP and other metrics')
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    analyze_parser = subparsers.add_parser('analyze', help='Run correlation analysis')
    summary_parser = subparsers.add_parser('summary', help='Create summary report from existing files')
    summary_parser.add_argument('--input-dir', required=True, help='Directory with correlation CSV files')
    summary_parser.add_argument('--output-file', help='Output file path (default: summary_report.csv in input directory)')

    args = parser.parse_args()

    if args.command == 'analyze':
        run_correlation_analysis()
    elif args.command == 'summary':
        create_summary_report(args.input_dir, args.output_file)
    else:
        print("No command provided. Running correlation analysis...")
        run_correlation_analysis()


if __name__ == "__main__":
    main()