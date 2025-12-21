import os
import glob
import re
import pandas as pd
from scipy import stats
from collections import defaultdict
import numpy as np
import logging
import argparse

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
    pattern = re.compile(r'.*_(?P<lab>\w+)_(?P<num>\d+)(?:_refined\d+)?_C1\'-gGSP\.csv$')
    
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
    """Load all metrics from the metrics directory with better error handling"""
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
            
            # Try to read CSV with different options
            try:
                df = pd.read_csv(file_path)
            except pd.errors.ParserError:
                logger.warning(f"Parser error for {file_path}, trying without header")
                df = pd.read_csv(file_path, header=None)
                if len(df.columns) >= 3:
                    df.columns = ["Lab", "Num", metric_name]
            
            # Validate columns
            if len(df.columns) < 3:
                logger.error(f"File {file_path} has only {len(df.columns)} columns, skipping")
                continue
            
            # Ensure we have the right column names
            if metric_name not in df.columns:
                # Assume the third column is the metric value
                value_col = df.columns[2]
                df = df.rename(columns={value_col: metric_name})
            
            # Check for required columns
            if "Lab" not in df.columns or "Num" not in df.columns:
                logger.error(f"Missing required columns in {file_path}, skipping")
                continue
            
            for _, row in df.iterrows():
                try:
                    lab = str(row["Lab"])
                    num = str(row["Num"])
                    value = row[metric_name]
                    
                    # Skip if value is None or NaN
                    if pd.isna(value):
                        continue
                    
                    # Try to convert to float, skip if not possible
                    try:
                        numeric_value = float(value)
                    except (ValueError, TypeError):
                        logger.warning(f"Non-numeric value in {metric_name} for {lab}_{num}: {value}")
                        continue
                    
                    all_metrics[(lab, num)][metric_name] = numeric_value
                    logger.debug(f"Loaded {metric_name} for {lab}_{num} = {numeric_value}")
                    
                except KeyError as e:
                    logger.error(f"Missing column in {file_path}: {e}")
                    break
                except Exception as e:
                    logger.error(f"Error processing row in {file_path}: {e}")
                    continue
                    
        except Exception as e:
            logger.error(f"Error processing {file_path}: {str(e)}", exc_info=True)
    
    return all_metrics

def calculate_correlations(gsp_df, other_metrics):
    """Calculate correlations between gGSP and other metrics with safe data handling"""
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
        
        # Prepare data for analysis with safe conversion
        clean_df = combined_df[["gGSP", metric]].copy()
        
        # Convert both columns to numeric, coercing errors to NaN
        clean_df["gGSP"] = pd.to_numeric(clean_df["gGSP"], errors='coerce')
        clean_df[metric] = pd.to_numeric(clean_df[metric], errors='coerce')
        
        # Remove rows with NaN values
        clean_df = clean_df.dropna()
        
        # Log data quality info
        original_count = len(combined_df)
        clean_count = len(clean_df)
        if clean_count < original_count:
            logger.warning(f"Removed {original_count - clean_count} rows with invalid data for {metric}")
        
        if len(clean_df) < 3:
            logger.warning(f"Insufficient data for {metric} (only {len(clean_df)} valid pairs). Skipping...")
            continue

        try:
            logger.debug(f"Calculating correlations for {metric} with {len(clean_df)} data points")

            # Extract x and y
            x = clean_df["gGSP"].to_numpy()
            y = clean_df[metric].to_numpy()

            # Adjust for orientation: invert if lower = better
            if not METRIC_ORIENTATION.get(metric, True):
                y = -y
                logger.debug(f"Inverted {metric} values for correlation (lower is better)")

            # Log xi, yi pairs with safe formatting
            logger.info(f"\n=== Data for metric: {metric} ===")
            for i, (xi, yi) in enumerate(zip(x, y), start=1):
                try:
                    # Try to format as floats
                    logger.info(f"Model {i}: gGSP={float(xi):.4f}, {metric}={float(yi):.4f}")
                except (ValueError, TypeError):
                    # Fallback to string representation
                    logger.info(f"Model {i}: gGSP={xi}, {metric}={yi}")

            # Log basic statistics
            logger.info(f"gGSP mean={np.mean(x):.4f}, std={np.std(x):.4f}")
            logger.info(f"{metric} mean={np.mean(y):.4f}, std={np.std(y):.4f}")

            # Calculate correlations
            try:
                # Pearson correlation
                pearson_corr, pearson_p = stats.pearsonr(x, y)
                logger.info(f"{metric} Pearson components: cov={np.cov(x, y)[0,1]:.4f}, "
                          f"var_x={np.var(x):.4f}, var_y={np.var(y):.4f}")
            except Exception as e:
                logger.error(f"Error calculating Pearson correlation for {metric}: {str(e)}")
                pearson_corr, pearson_p = np.nan, np.nan
            
            try:
                # Spearman correlation
                spearman_corr, spearman_p = stats.spearmanr(x, y)
                logger.info(f"{metric} Spearman ranks (x): {stats.rankdata(x)}")
                logger.info(f"{metric} Spearman ranks (y): {stats.rankdata(y)}")
            except Exception as e:
                logger.error(f"Error calculating Spearman correlation for {metric}: {str(e)}")
                spearman_corr, spearman_p = np.nan, np.nan
            
            try:
                # Kendall correlation
                kendall_corr, kendall_p = stats.kendalltau(x, y)
                logger.info(f"{metric} Kendall concordant pairs check: approx={kendall_corr:.3f}")
            except Exception as e:
                logger.error(f"Error calculating Kendall correlation for {metric}: {str(e)}")
                kendall_corr, kendall_p = np.nan, np.nan

            # Store results (even if some correlations failed)
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
    """Process a single target and calculate correlations"""
    logger.info(f"Processing target: {target_name}")
    
    # Set up paths
    gsp_dir = os.path.join(base_gsp_dir, target_name)
    metrics_dir = os.path.join(base_metrics_dir, target_name)
    
    if not os.path.exists(gsp_dir):
        logger.error(f"GSP directory not found: {gsp_dir}")
        return None
    if not os.path.exists(metrics_dir):
        logger.error(f"Metrics directory not found: {metrics_dir}")
        return None
    
    # Load data
    gsp_df = parse_gsp_files(gsp_dir)
    if gsp_df.empty:
        logger.error(f"No gGSP data found for {target_name}")
        return None
    
    other_metrics = load_other_metrics(metrics_dir)
    if not other_metrics:
        logger.error(f"No other metrics found for {target_name}")
        return None
    
    # Calculate correlations
    correlations_df = calculate_correlations(gsp_df, other_metrics)
    if correlations_df.empty:
        logger.error(f"No correlation results for {target_name}")
        return None
    
    # Save results
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{target_name}_correlations.csv")
    correlations_df.to_csv(output_file, index=False)
    logger.info(f"Saved correlation results to: {output_file}")
    
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

def create_summary_report(correlation_dir, output_file=None):
    """
    Tworzy ogólne podsumowanie z wszystkich plików CSV w katalogu.
    Sumuje wartości i oblicza średnie dla każdej metryki.
    """
    logger.info(f"Tworzenie raportu podsumowującego z katalogu: {correlation_dir}")
    
    if not os.path.exists(correlation_dir):
        logger.error(f"Katalog nie istnieje: {correlation_dir}")
        return None
    
    # Znajdź wszystkie pliki CSV z korelacjami
    correlation_files = glob.glob(os.path.join(correlation_dir, "*_correlations.csv"))
    
    if not correlation_files:
        logger.error(f"Nie znaleziono plików *_correlations.csv w {correlation_dir}")
        return None
    
    logger.info(f"Znaleziono {len(correlation_files)} plików z korelacjami")
    
    # Wczytaj wszystkie pliki
    all_dataframes = []
    all_target_names = []
    
    for file_path in correlation_files:
        try:
            df = pd.read_csv(file_path)
            target_name = os.path.basename(file_path).replace('_correlations.csv', '')
            df['source_file'] = target_name
            all_dataframes.append(df)
            all_target_names.append(target_name)
            logger.debug(f"Wczytano {os.path.basename(file_path)} z {len(df)} wierszami")
        except Exception as e:
            logger.error(f"Błąd wczytywania {file_path}: {str(e)}")
    
    if not all_dataframes:
        logger.error("Nie wczytano danych z plików korelacji")
        return None
    
    # Połącz wszystkie ramki danych
    combined_df = pd.concat(all_dataframes, ignore_index=True)
    logger.info(f"Łączna liczba wierszy w połączonych danych: {len(combined_df)}")
    
    # Sprawdź wymagane kolumny
    required_columns = ['metric', 'orientation', 'pearson_coef', 'pearson_pvalue',
                       'spearman_coef', 'spearman_pvalue', 'kendall_coef', 
                       'kendall_pvalue', 'n_models']
    
    missing_columns = [col for col in required_columns if col not in combined_df.columns]
    if missing_columns:
        logger.error(f"Missing required columns: {missing_columns}")
        logger.info(f"Available columns: {combined_df.columns.tolist()}")
        return None
    
    # PRZYGOTUJ SZCZEGÓŁOWE INFORMACJE O BRAKACH
    # Wczytaj oryginalne pliki z danymi dla każdego targetu
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    BASE_GSP_DIR = os.path.join(BASE_DIR, "results", "gsp", "rnapuzzles.github.io")
    BASE_METRICS_DIR = os.path.join(BASE_DIR, "results", "other-metrics", "rnapuzzles.github.io")
    
    detailed_missing_data = {}
    detailed_stats = {}
    
    # Najpierw zbiorczo przetwórz wszystkie targety
    for target in all_target_names:
        logger.info(f"Analizuję szczegóły braków dla targetu: {target}")
        
        # Ścieżki do danych dla tego targetu
        gsp_dir = os.path.join(BASE_GSP_DIR, target)
        metrics_dir = os.path.join(BASE_METRICS_DIR, target)
        
        if not os.path.exists(gsp_dir):
            logger.warning(f"Brak katalogu gGSP dla targetu {target}")
            detailed_stats[target] = {"gsp_models": 0, "has_gsp": False}
            continue
        
        if not os.path.exists(metrics_dir):
            logger.warning(f"Brak katalogu metryk dla targetu {target}")
            detailed_stats[target] = {"gsp_models": 0, "has_metrics": False}
            continue
        
        # Wczytaj modele gGSP
        try:
            gsp_df = parse_gsp_files(gsp_dir)
            if gsp_df.empty:
                logger.warning(f"Brak danych gGSP dla targetu {target}")
                detailed_stats[target] = {"gsp_models": 0, "gsp_data": []}
                continue
            
            # Zbuduj listę modeli z gGSP
            gsp_models = []
            for _, row in gsp_df.iterrows():
                lab = row["Lab"]
                num = str(row["Num"])
                gsp_models.append((lab, num))
            
            detailed_stats[target] = {
                "gsp_models": len(gsp_models),
                "gsp_data": gsp_models,
                "has_gsp": True
            }
            
        except Exception as e:
            logger.error(f"Błąd wczytywania gGSP dla {target}: {str(e)}")
            detailed_stats[target] = {"gsp_models": 0, "has_gsp": False}
            continue
        
        # Dla każdej metryki, sprawdź które modele z gGSP mają tę metrykę
        for metric in METRIC_ORIENTATION.keys():
            if metric == "gGSP":
                continue
            
            metric_file = os.path.join(metrics_dir, f"{metric}.csv")
            if not os.path.exists(metric_file):
                # Całkowity brak pliku metryki
                if metric not in detailed_missing_data:
                    detailed_missing_data[metric] = {}
                detailed_missing_data[metric][target] = {
                    "missing": gsp_models,  # Wszystkie modele brakują
                    "status": "missing_file"
                }
                continue
            
            # Wczytaj plik metryki
            try:
                metric_df = pd.read_csv(metric_file)
            except Exception as e:
                logger.error(f"Błąd wczytywania {metric_file}: {str(e)}")
                if metric not in detailed_missing_data:
                    detailed_missing_data[metric] = {}
                detailed_missing_data[metric][target] = {
                    "missing": gsp_models,
                    "status": "error_reading"
                }
                continue
            
            # Zbuduj listę modeli z metryki
            metric_models = []
            for _, row in metric_df.iterrows():
                lab = row["Lab"]
                num = str(row["Num"])
                metric_models.append((lab, num))
            
            # Znajdź modele, które są w gGSP ale nie w metryce
            missing_models = []
            for model in gsp_models:
                if model not in metric_models:
                    missing_models.append(model)
            
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
                # Wszystkie modele mają metrykę
                if metric not in detailed_missing_data:
                    detailed_missing_data[metric] = {}
                detailed_missing_data[metric][target] = {
                    "missing": [],
                    "total_gsp": len(gsp_models),
                    "total_metric": len(metric_models),
                    "status": "complete"
                }
    
    # PRZYGOTUJ DANE DO PODSUMOWANIA (jak wcześniej)
    summary_data = []
    grouped = combined_df.groupby(['metric', 'orientation'])
    
    for (metric, orientation), group in grouped:
        # Oblicz statystyki opisowe
        n_targets = group['source_file'].nunique()
        total_models = group['n_models'].sum()
        
        # Oblicz średnie współczynniki korelacji
        avg_pearson = group['pearson_coef'].mean()
        avg_spearman = group['spearman_coef'].mean()
        avg_kendall = group['kendall_coef'].mean()
        
        # Oblicz średnie p-wartości
        avg_pearson_p = group['pearson_pvalue'].mean()
        avg_spearman_p = group['spearman_pvalue'].mean()
        avg_kendall_p = group['kendall_pvalue'].mean()
        
        # Oblicz odchylenia standardowe
        std_pearson = group['pearson_coef'].std()
        std_spearman = group['spearman_coef'].std()
        std_kendall = group['kendall_coef'].std()
        
        summary_data.append({
            'metric': metric,
            'orientation': orientation,
            'n_targets': n_targets,
            'total_models': total_models,
            'avg_pearson_coef': avg_pearson,
            'std_pearson_coef': std_pearson,
            'avg_pearson_pvalue': avg_pearson_p,
            'avg_spearman_coef': avg_spearman,
            'std_spearman_coef': std_spearman,
            'avg_spearman_pvalue': avg_spearman_p,
            'avg_kendall_coef': avg_kendall,
            'std_kendall_coef': std_kendall,
            'avg_kendall_pvalue': avg_kendall_p,
            'min_pearson': group['pearson_coef'].min(),
            'max_pearson': group['pearson_coef'].max(),
            'min_spearman': group['spearman_coef'].min(),
            'max_spearman': group['spearman_coef'].max(),
            'min_kendall': group['kendall_coef'].min(),
            'max_kendall': group['kendall_coef'].max()
        })
    
    # Utwórz ramkę danych z podsumowaniem
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values('avg_spearman_coef', ascending=False)
    
    # Domyślna ścieżka wyjściowa
    if output_file is None:
        output_file = os.path.join(correlation_dir, "summary_report.csv")
    
    # Zapisz podsumowanie
    summary_df.to_csv(output_file, index=False, float_format='%.6f')
    
    # Utwórz również wersję z zaokrąglonymi wartościami
    summary_rounded = summary_df.copy()
    numeric_cols = summary_rounded.select_dtypes(include=[np.number]).columns
    summary_rounded[numeric_cols] = summary_rounded[numeric_cols].round(4)
    
    summary_rounded_file = output_file.replace('.csv', '_rounded.csv')
    summary_rounded.to_csv(summary_rounded_file, index=False)
    
    # GENERUJ SZCZEGÓŁOWY RAPORT O BRAKACH
    report_lines = []
    report_lines.append("="*120)
    report_lines.append(f"ŁĄCZNA ANALIZA DANYCH DLA {len(all_target_names)} TARGETÓW")
    report_lines.append("="*120)
    
    # Statystyki ogólne
    total_gsp_models = sum(stats.get("gsp_models", 0) for stats in detailed_stats.values())
    targets_with_gsp = sum(1 for stats in detailed_stats.values() if stats.get("has_gsp", False))
    
    report_lines.append(f"\n📊 STATYSTYKI OGÓLNE:")
    report_lines.append(f"   Targety z danymi gGSP: {targets_with_gsp}/{len(all_target_names)}")
    report_lines.append(f"   Łączna liczba modeli gGSP: {total_gsp_models}")
    
    # Analiza każdej metryki
    for metric in METRIC_ORIENTATION.keys():
        if metric == "gGSP":
            continue
        
        report_lines.append(f"\n{'='*80}")
        report_lines.append(f"ANALIZA METRYKI: {metric}")
        report_lines.append(f"{'='*80}")
        
        if metric not in detailed_missing_data:
            report_lines.append(f"❌ BRAK DANYCH DLA TEJ METRYKI WE WSZYSTKICH TARGETACH")
            continue
        
        missing_info = detailed_missing_data[metric]
        
        # Liczby
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
        
        report_lines.append(f"\n📈 PODSUMOWANIE:")
        report_lines.append(f"   Targety z kompletem danych: {targets_with_metric}/{total_targets_for_metric}")
        report_lines.append(f"   Targety z częściowymi danymi: {targets_partial_data}/{total_targets_for_metric}")
        report_lines.append(f"   Targety bez pliku: {targets_missing_file}/{total_targets_for_metric}")
        report_lines.append(f"   Łączna liczba brakujących modeli: {total_missing_models}")
        
        # Szczegółowe informacje o brakujących targetach
        if targets_missing_file > 0:
            report_lines.append(f"\n📁 TARGETY BEZ PLIKU {metric}.csv:")
            missing_file_targets = [t for t, info in missing_info.items() if info["status"] == "missing_file"]
            for i, target in enumerate(missing_file_targets[:10], 1):
                gsp_count = detailed_stats.get(target, {}).get("gsp_models", 0)
                report_lines.append(f"   {i}. {target} (brak pliku, {gsp_count} modeli gGSP)")
            if len(missing_file_targets) > 10:
                report_lines.append(f"   ... i {len(missing_file_targets) - 10} więcej")
        
        # Szczegółowe informacje o modelach z częściowymi danymi
        if targets_partial_data > 0:
            report_lines.append(f"\n⚠️  TARGETY Z CZĘŚCIOWYMI DANYMI:")
            partial_targets = [(t, info) for t, info in missing_info.items() if info["status"] == "partial_data"]
            partial_targets.sort(key=lambda x: len(x[1]["missing"]), reverse=True)
            
            for i, (target, info) in enumerate(partial_targets[:5], 1):
                missing_count = len(info["missing"])
                total_gsp = info.get("total_gsp", 0)
                total_metric = info.get("total_metric", 0)
                
                report_lines.append(f"\n   {i}. {target}:")
                report_lines.append(f"       Modele gGSP: {total_gsp}")
                report_lines.append(f"       Modele z {metric}: {total_metric}")
                report_lines.append(f"       Brakujące modele: {missing_count} ({missing_count/total_gsp*100:.1f}%)")
                
                # Grupuj brakujące modele po lab
                models_by_lab = {}
                for lab, num in info["missing"]:
                    if lab not in models_by_lab:
                        models_by_lab[lab] = []
                    models_by_lab[lab].append(num)
                
                # Pokaż 3 pierwsze laby
                for j, (lab, nums) in enumerate(sorted(models_by_lab.items())[:3]):
                    nums_str = ", ".join(sorted(nums, key=lambda x: int(x) if x.isdigit() else x))
                    report_lines.append(f"       {lab}: {nums_str}")
                
                if len(models_by_lab) > 3:
                    report_lines.append(f"       ... i {len(models_by_lab) - 3} innych labów")
            
            if len(partial_targets) > 5:
                report_lines.append(f"\n   ... i {len(partial_targets) - 5} więcej targetów z częściowymi danymi")
    
    # Dodaj podsumowanie z pliku correlations.csv
    report_lines.append(f"\n{'='*120}")
    report_lines.append("PODSUMOWANIE Z PLIKU correlations.csv")
    report_lines.append(f"{'='*120}")
    
    # Oblicz statystyki z correlations.csv
    for metric in METRIC_ORIENTATION.keys():
        if metric == "gGSP":
            continue
        
        metric_data = combined_df[combined_df['metric'] == metric]
        if not metric_data.empty:
            n_targets = metric_data['source_file'].nunique()
            avg_spearman = metric_data['spearman_coef'].mean()
            avg_pearson = metric_data['pearson_coef'].mean()
            
            status = "✅" if n_targets == len(all_target_names) else "⚠️"
            report_lines.append(f"{status} {metric}: {n_targets} targetów, Spearman={avg_spearman:.3f}, Pearson={avg_pearson:.3f}")
    
    # Zapisz raport do pliku
    full_report = "\n".join(report_lines)
    
    # Wyświetl krótkie podsumowanie w konsoli
    print("="*120)
    print("RAPORT PODSUMOWUJĄCY - ANALIZA BRAKÓW DANYCH")
    print("="*120)
    print(f"Łączna liczba targetów: {len(all_target_names)}")
    print(f"Targety z danymi gGSP: {targets_with_gsp}")
    print(f"Łączna liczba modeli gGSP: {total_gsp_models}")
    print("\n" + "="*120)
    print("STATUS METRYK (z pliku correlations.csv):")
    
    for metric in METRIC_ORIENTATION.keys():
        if metric == "gGSP":
            continue
        
        metric_data = combined_df[combined_df['metric'] == metric]
        if not metric_data.empty:
            n_targets = metric_data['source_file'].nunique()
            avg_spearman = metric_data['spearman_coef'].mean()
            status = "PEŁNE" if n_targets == len(all_target_names) else f"BRAKUJE {len(all_target_names) - n_targets}"
            print(f"  {metric}: {n_targets} targetów, Spearman={avg_spearman:.3f} - {status}")
    
    # Zapisz pełny raport do pliku
    report_file_path = os.path.join(correlation_dir, "detailed_missing_analysis.txt")
    with open(report_file_path, 'w', encoding='utf-8') as f:
        f.write(full_report)
    
    print(f"\n📋 Pełny raport zapisano do: {report_file_path}")
    
    logger.info(f"Raport podsumowujący zapisany jako: {output_file}")
    logger.info(f"Zaokrąglony raport zapisany jako: {summary_rounded_file}")
    logger.info(f"Szczegółowy raport o brakach zapisany jako: {report_file_path}")
    
    return summary_df

def run_correlation_analysis():
    """Główna funkcja do uruchomienia analizy korelacji"""
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
            # Print summary for this target
            print(f"\nCorrelation summary for {target}:")
            print(result[["metric", "pearson_coef", "spearman_coef", "kendall_coef", "n_models"]])
    
    logger.info("Processing completed")
    print(f"\nAll processing completed. Results saved in: {OUTPUT_DIR}")
    print(f"Processed {len(all_results)} targets.")
    
    # Ask if user wants to create summary report
    print("\n" + "="*80)
    create_summary = input("Do you want to create a summary report? (y/n): ").strip().lower()
    
    if create_summary == 'y':
        print(f"\nCreating summary report from: {OUTPUT_DIR}")
        summary = create_summary_report(OUTPUT_DIR)
        if summary is not None:
            print(f"\nSummary report created successfully!")
    
    print("\nDone!")

def main():
    """Główna funkcja z obsługą argumentów wiersza poleceń"""
    parser = argparse.ArgumentParser(description='Analyze correlations between gGSP and other metrics')
    
    # Subparsers dla różnych trybów działania
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Parser dla analizy korelacji
    analysis_parser = subparsers.add_parser('analyze', help='Run correlation analysis')
    
    # Parser dla tworzenia raportu podsumowującego
    summary_parser = subparsers.add_parser('summary', help='Create summary report from existing correlation files')
    summary_parser.add_argument('--input-dir', required=True, help='Directory with correlation CSV files')
    summary_parser.add_argument('--output-file', help='Output file path (default: summary_report.csv in input directory)')
    
    args = parser.parse_args()
    
    if args.command == 'analyze':
        # Uruchom analizę korelacji
        run_correlation_analysis()
    elif args.command == 'summary':
        # Utwórz raport podsumowujący
        create_summary_report(args.input_dir, args.output_file)
    else:
        # Domyślnie uruchom analizę korelacji
        print("No command specified. Running correlation analysis...")
        run_correlation_analysis()

if __name__ == "__main__":
    main()