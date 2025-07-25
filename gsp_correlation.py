import os
import glob
import re
import pandas as pd
from scipy import stats
from collections import defaultdict
import numpy as np
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('gsp_correlation.log', mode='w')
    ]
)
logger = logging.getLogger(__name__)

# Metric orientation dictionary
METRIC_ORIENTATION = {
    "RMSD": False,
    "DI_all": False,
    "INF_all": True,
    "INF_wc": True,
    "INF_stacking": True,
    "Clash_Score": False,
    "P-value": False,
    "mcq": False,
    "TM-score": True,
    "gGSP": True
}

def parse_gsp_files(gsp_dir):
    """Collect gGSP values from all files in the GSP directory"""
    gsp_files = glob.glob(os.path.join(gsp_dir, "*_C1'-gGSP.csv"))
    gsp_data = []
    logger.info(f"Found {len(gsp_files)} gGSP files in {gsp_dir}")
    
    # Pattern to extract model identifier
    pattern = re.compile(r'.*_(?P<lab>\w+)_(?P<num>\d+)_C1\'-gGSP\.csv$')
    
    for file_path in gsp_files:
        try:
            filename = os.path.basename(file_path)
            logger.debug(f"Processing gGSP file: {filename}")
            
            match = pattern.search(filename)
            if not match:
                logger.warning(f"Pattern not matched for: {filename}")
                continue
                
            lab = match.group('lab')
            num = match.group('num')
            
            # Read gGSP value
            with open(file_path, 'r') as f:
                # Skip header
                header = f.readline().strip()
                value = float(f.readline().strip())
            
            gsp_data.append({
                "Lab": lab,
                "Num": num,
                "gGSP": value
            })
            logger.debug(f"Extracted gGSP: {lab}_{num} = {value}")
        except Exception as e:
            logger.error(f"Error processing {file_path}: {str(e)}", exc_info=True)
    
    return pd.DataFrame(gsp_data)

def load_other_metrics(metrics_dir):
    """Load all metrics from the metrics directory"""
    metrics_files = glob.glob(os.path.join(metrics_dir, "*.csv"))
    all_metrics = defaultdict(dict)
    logger.info(f"Found {len(metrics_files)} metric files in {metrics_dir}")
    
    # List of expected metrics
    expected_metrics = list(METRIC_ORIENTATION.keys())
    
    for file_path in metrics_files:
        metric_name = os.path.basename(file_path).replace(".csv", "")
        
        # Skip unexpected metrics and gGSP
        if metric_name not in expected_metrics or metric_name == "gGSP":
            logger.debug(f"Skipping unexpected metric file: {os.path.basename(file_path)}")
            continue
            
        try:
            logger.debug(f"Processing metric file: {os.path.basename(file_path)}")
            df = pd.read_csv(file_path)
            
            for _, row in df.iterrows():
                lab = row["Lab"]
                num = str(row["Num"])
                value = row[df.columns[2]]  
                
                all_metrics[(lab, num)][metric_name] = value
                logger.debug(f"Loaded {metric_name} for {lab}_{num} = {value}")
        except Exception as e:
            logger.error(f"Error processing {file_path}: {str(e)}", exc_info=True)
    
    return all_metrics

def calculate_correlations(gsp_df, other_metrics):
    """Calculate correlations between gGSP and other metrics"""
    results = []
    logger.info("Starting correlation calculations...")
    
    # Prepare combined data
    combined_data = []
    for _, row in gsp_df.iterrows():
        key = (row["Lab"], str(row["Num"]))
        if key in other_metrics:
            entry = {"Lab": key[0], "Num": key[1], "gGSP": row["gGSP"]}
            entry.update(other_metrics[key])
            combined_data.append(entry)
            logger.debug(f"Combined data for {key}: {entry}")
        else:
            logger.warning(f"No metrics found for model: {key}")
    
    if not combined_data:
        logger.error("No data available for correlation calculations!")
        return pd.DataFrame()
    
    combined_df = pd.DataFrame(combined_data)
    logger.info(f"Combined dataset contains {len(combined_df)} models")
    logger.debug(f"Combined data columns: {combined_df.columns.tolist()}")
    
    # Calculate correlations for each metric
    for metric in combined_df.columns:
        if metric in ["Lab", "Num"] or metric == "gGSP":
            continue
            
        # Prepare data for analysis
        clean_df = combined_df[["gGSP", metric]].dropna()
        
        if len(clean_df) < 3:
            logger.warning(f"Insufficient data for {metric} (only {len(clean_df)} valid pairs). Skipping...")
            continue
        
        # Calculate correlations
        try:
            logger.debug(f"Calculating correlations for {metric} with {len(clean_df)} data points")
            
            # Pearson correlation
            pearson_corr, pearson_p = stats.pearsonr(clean_df["gGSP"], clean_df[metric])
            
            # Spearman correlation
            spearman_corr, spearman_p = stats.spearmanr(clean_df["gGSP"], clean_df[metric])
            
            # Kendall correlation
            kendall_corr, kendall_p = stats.kendalltau(clean_df["gGSP"], clean_df[metric])
            
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
            
            logger.info(f"Correlations for {metric}: "
                        f"Pearson={pearson_corr:.3f}, "
                        f"Spearman={spearman_corr:.3f}, "
                        f"Kendall={kendall_corr:.3f}")
        except Exception as e:
            logger.error(f"Error calculating correlations for {metric}: {str(e)}", exc_info=True)
    
    return pd.DataFrame(results)

def process_target(target_name, base_gsp_dir, base_metrics_dir, output_dir):
    """Process a single target and calculate correlations"""
    logger.info(f"Processing target: {target_name}")
    
    # Set up paths
    gsp_dir = os.path.join(base_gsp_dir, target_name)
    metrics_dir = os.path.join(base_metrics_dir, target_name)
    
    if not os.path.exists(gsp_dir):
        logger.error(f"GSP directory not found: {gsp_dir}")
        return
    if not os.path.exists(metrics_dir):
        logger.error(f"Metrics directory not found: {metrics_dir}")
        return
    
    # Load data
    gsp_df = parse_gsp_files(gsp_dir)
    if gsp_df.empty:
        logger.error(f"No gGSP data found for {target_name}")
        return
    
    other_metrics = load_other_metrics(metrics_dir)
    if not other_metrics:
        logger.error(f"No other metrics found for {target_name}")
        return
    
    # Calculate correlations
    correlations_df = calculate_correlations(gsp_df, other_metrics)
    if correlations_df.empty:
        logger.error(f"No correlation results for {target_name}")
        return
    
    # Save results
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{target_name}_correlations.csv")
    correlations_df.to_csv(output_file, index=False)
    logger.info(f"Saved correlation results to: {output_file}")
    
    # Print summary
    print(f"\nCorrelation summary for {target_name}:")
    print(correlations_df[["metric", "pearson_coef", "spearman_coef", "kendall_coef", "n_models"]])
    
    return correlations_df

def select_targets(base_dir):
    """Interactively select target directories"""
    print("\n===== Target Selection =====")
    print("Please select target datasets to process:")
    
    # Get available targets
    available_targets = [d for d in os.listdir(base_dir) 
                       if os.path.isdir(os.path.join(base_dir, d))]
    
    if not available_targets:
        print("No targets found!")
        return []
    
    print("Available targets:")
    for i, target in enumerate(available_targets, 1):
        print(f"{i}. {target}")
    
    # Allow multiple selection
    print("\nEnter target numbers (e.g. 1,3,5-7) or 'all'")
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
    
    return list(set(selected_targets))  # Remove duplicates

def main():
    # Configure base paths
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    BASE_GSP_DIR = os.path.join(BASE_DIR, "results", "gsp", "rnapuzzles.github.io")
    BASE_METRICS_DIR = os.path.join(BASE_DIR, "results", "other-metrics", "rnapuzzles.github.io")
    OUTPUT_DIR = os.path.join(BASE_DIR, "results", "correlations")
    
    # Create directories if needed
    os.makedirs(BASE_GSP_DIR, exist_ok=True)
    os.makedirs(BASE_METRICS_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    logger.info(f"Base GSP directory: {BASE_GSP_DIR}")
    logger.info(f"Base metrics directory: {BASE_METRICS_DIR}")
    
    # Select targets
    selected_targets = select_targets(BASE_GSP_DIR)
    
    if not selected_targets:
        logger.error("No targets selected!")
        return
    
    logger.info(f"Selected targets: {', '.join(selected_targets)}")
    
    # Process each selected target
    all_results = {}
    for target in selected_targets:
        result = process_target(
            target, 
            BASE_GSP_DIR, 
            BASE_METRICS_DIR, 
            OUTPUT_DIR
        )
        if result is not None:
            all_results[target] = result
    
    logger.info("Processing completed")
    print("\nAll processing completed. Results saved in:", OUTPUT_DIR)

if __name__ == "__main__":
    main()