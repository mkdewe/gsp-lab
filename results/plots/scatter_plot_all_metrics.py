import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from matplotlib.lines import Line2D

# --- KONFIGURACJA ŚCIEŻEK ---
GSP_AGGREGATED_CSV = r"D:\Studia\Projekty\gsp-lab\results\plots\gGSP_5.0.csv"
OTHER_METRICS_DIR = r"D:\Studia\Projekty\gsp-lab\results\other-metrics\rnapuzzles.github.io"
LDDT_GLOBAL_CSV = r"D:\Studia\Projekty\gsp-lab\results\lddt\finished\global_report.csv"
CSV_MAPPING_FILE = r"D:\Studia\Projekty\gsp-lab\results\plots\RNA-types.csv"
OUTPUT_DIR = r"D:\Studia\Projekty\gsp-lab\results\plots"

METRICS = ['RMSD', 'TM-score', 'INF_all', 'lDDT']
METRIC_COLORS = {
    'RMSD': '#66c2a5',
    'TM-score': '#fc8d62',
    'INF_all': '#8da0cb',
    'lDDT': '#ffd92f'
}
EXCLUDED_TARGETS = ['PZ10tRNA', 'PZ26tRNA', 'PZ27tRNA', 'PZ28tRNA']

# --- KONFIGURACJA ESTETYKI ---
sns.set_theme(style="whitegrid")
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'axes.labelweight': 'bold',
    'axes.titleweight': 'bold'
})

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ----------------------------------------------------------------------
# 1. FUNKCJE POMOCNICZE ŁADOWANIA DANYCH
# ----------------------------------------------------------------------
def split_modelgroup(mg):
    parts = mg.split('_')
    lab = '_'.join(parts[:-1]) if len(parts) >= 2 else mg
    num = parts[-1] if len(parts) >= 2 else '0'
    return lab.lower().strip(), str(num).strip()

def load_aggregated_gsp(csv_path):
    df = pd.read_csv(csv_path)
    df = df[~df['Target'].isin(EXCLUDED_TARGETS)]
    df[['Lab', 'Num']] = df['ModelGroup'].apply(lambda x: pd.Series(split_modelgroup(x)))
    return df

def load_aggregated_lddt(csv_path):
    if not os.path.exists(csv_path):
        return None
    df = pd.read_csv(csv_path)
    df.rename(columns={'Puzzle': 'Target'}, inplace=True)
    df = df[~df['Target'].isin(EXCLUDED_TARGETS)]
    df[['Lab', 'Num']] = df['ModelGroup'].apply(lambda x: pd.Series(split_modelgroup(x)))
    return df[['Target', 'Lab', 'Num', 'lDDT']]

def load_other_metrics(target_dir):
    result = None
    for metric in ['RMSD', 'TM-score', 'INF_all']:
        m_file = os.path.join(target_dir, f"{metric}.csv")
        if not os.path.exists(m_file):
            continue

        df_m = pd.read_csv(m_file)
        df_m.columns = [c.capitalize() if c.lower() in ['lab', 'num'] else c for c in df_m.columns]

        df_m['Lab'] = df_m['Lab'].astype(str).str.lower().str.strip()
        df_m['Num'] = df_m['Num'].astype(str).str.strip()

        val_col = [c for c in df_m.columns if c not in ['Lab', 'Num']][0]
        df_m = df_m[['Lab', 'Num', val_col]].rename(columns={val_col: metric})

        result = df_m if result is None else pd.merge(result, df_m, on=['Lab', 'Num'], how='outer')

    return result

def load_all_data():
    gsp_df = load_aggregated_gsp(GSP_AGGREGATED_CSV)
    lddt_df = load_aggregated_lddt(LDDT_GLOBAL_CSV)

    mapping_df = pd.read_csv(CSV_MAPPING_FILE)
    target_mapping = pd.Series(mapping_df.type.values, index=mapping_df.target_label).to_dict()

    merged_dfs = []

    for target in gsp_df['Target'].unique():
        gsp_target = gsp_df[gsp_df['Target'] == target].copy()

        other_path = os.path.join(OTHER_METRICS_DIR, target)
        other_metrics = load_other_metrics(other_path) if os.path.exists(other_path) else None

        merged = gsp_target

        if other_metrics is not None:
            merged = pd.merge(merged, other_metrics, on=['Lab', 'Num'], how='left')

        if lddt_df is not None:
            lddt_target = lddt_df[lddt_df['Target'] == target]
            merged = pd.merge(merged, lddt_target[['Lab', 'Num', 'lDDT']], on=['Lab', 'Num'], how='left')
        else:
            merged['lDDT'] = np.nan

        merged['Type'] = target_mapping.get(target, 'Other RNA Structures')
        merged_dfs.append(merged)

    return pd.concat(merged_dfs, ignore_index=True)

def safe_savefig(fig, path, **kwargs):
    try:
        fig.savefig(path, **kwargs)
        return path
    except PermissionError:
        base, ext = os.path.splitext(path)
        alt_path = f"{base}_new{ext}"
        fig.savefig(alt_path, **kwargs)
        return alt_path

# ----------------------------------------------------------------------
# 2. MAIN
# ----------------------------------------------------------------------
def main():
    full_df = load_all_data()

    # --- OBLICZANIE KORELACJI SPEARMANA ---
    spearman_data = []
    for (target, rna_type), group in full_df.groupby(['Target', 'Type']):
        for metric in METRICS:
            valid = group[['gGSP', metric]].dropna()
            if len(valid) > 3:
                rho, _ = spearmanr(valid['gGSP'], valid[metric])
                spearman_data.append({
                    'Target': target,
                    'Type': rna_type,
                    'Metric': metric,
                    'Rho': rho
                })

    df_rho = pd.DataFrame(spearman_data)

    if df_rho.empty:
        print("Brak danych do wykresu.")
        return

    # --- etykiety typów z n ---
    counts = df_rho.groupby('Type')['Target'].nunique()
    label_map = {t: f"{t}\n(n={counts[t]})" for t in counts.index}
    df_rho['Type_Label'] = df_rho['Type'].map(label_map)

    mapping_df = pd.read_csv(CSV_MAPPING_FILE)
    type_order = [label_map[t] for t in mapping_df['type'].unique() if t in label_map]
    type_to_x = {label: i for i, label in enumerate(type_order)}

    # --- pozycje dla metryk ---
    offsets = np.linspace(-0.30, 0.30, len(METRICS))

    # ------------------------------------------------------------------
    # 3. WYKRES 3B
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(14, 5))

    # Zmniejszenie grubości o ~35%
    violin_width = 0.16
    violin_edge_lw = 0.65
    stat_lw = 0.59
    stat_cap = 0.018
    mean_marker_size = 12
    single_marker_size = 40
    point_edge_lw = 0.45

    # najpierw rysujemy tylko grupy z n>1 jako violiny
    stats = (
        df_rho.groupby(['Type_Label', 'Metric'])['Rho']
        .agg(['mean', 'std', 'count'])
        .reset_index()
    )

    multi_df = stats[stats['count'] > 1].copy()

    for _, row in multi_df.iterrows():
        type_label = row['Type_Label']
        metric = row['Metric']
        x = type_to_x[type_label] + offsets[METRICS.index(metric)]

        subset = df_rho[
            (df_rho['Type_Label'] == type_label) &
            (df_rho['Metric'] == metric)
        ]['Rho'].dropna().values

        if len(subset) <= 1:
            continue

        parts = ax.violinplot(
            dataset=subset,
            positions=[x],
            widths=violin_width,
            showmeans=False,
            showmedians=False,
            showextrema=False,
            points=200
        )

        for body in parts['bodies']:
            body.set_facecolor(METRIC_COLORS[metric])
            body.set_edgecolor('black')
            body.set_linewidth(violin_edge_lw)
            body.set_alpha(0.85)

    # --- n=1: tylko kolorowa kropka z czarną obwódką ---
    single_stats = stats[stats['count'] == 1].copy()

    for _, row in single_stats.iterrows():
        type_label = row['Type_Label']
        metric = row['Metric']
        x = type_to_x[type_label] + offsets[METRICS.index(metric)]

        y = float(
            df_rho[
                (df_rho['Type_Label'] == type_label) &
                (df_rho['Metric'] == metric)
            ]['Rho'].iloc[0]
        )

        ax.scatter(
            x,
            y,
            s=single_marker_size,
            facecolors=METRIC_COLORS[metric],
            edgecolors='black',
            linewidths=point_edge_lw,
            marker='o',
            zorder=20
        )

    # --- n>1: mean + SD, węższe elementy ---
    for _, row in stats.iterrows():
        if int(row['count']) <= 1:
            continue

        type_label = row['Type_Label']
        metric = row['Metric']
        x = type_to_x[type_label] + offsets[METRICS.index(metric)]

        mean = float(row['mean'])
        sd = float(row['std']) if pd.notna(row['std']) else 0.0

        ax.plot(
            [x, x],
            [mean - sd, mean + sd],
            color='black',
            linewidth=stat_lw,
            zorder=15
        )

        ax.plot([x - stat_cap, x + stat_cap], [mean - sd, mean - sd], color='black', linewidth=stat_lw, zorder=15)
        ax.plot([x - stat_cap, x + stat_cap], [mean + sd, mean + sd], color='black', linewidth=stat_lw, zorder=15)

        ax.scatter(
            x,
            mean,
            marker='D',
            s=mean_marker_size,
            color='black',
            zorder=16
        )

    # --- estetyka osi ---
    ax.axhline(0, color="red", linestyle="--", alpha=0.4, linewidth=0.65)
    ax.set_title("Spearman Rank Correlation: gGSP vs Benchmarking Metrics", fontsize=13, pad=10)
    ax.set_ylabel("Spearman ρ", fontsize=11)
    ax.set_xlabel("")

    ax.set_xticks(range(len(type_order)))
    ax.set_xticklabels(type_order, rotation=0, ha="center")

    # --- legenda z kropkami i czarną obwódką ---
    legend_handles = [
        Line2D(
            [0], [0],
            marker='o',
            linestyle='None',
            markersize=7,
            markerfacecolor=METRIC_COLORS[m],
            markeredgecolor='black',
            markeredgewidth=0.7,
            label=m
        )
        for m in METRICS
    ]
    ax.legend(
        handles=legend_handles,
        title="Metric",
        loc="lower right",
        frameon=True
    )

    plt.tight_layout()

    pdf_path = os.path.join(OUTPUT_DIR, "Figure3b_FINAL.pdf")
    png_path = os.path.join(OUTPUT_DIR, "Figure3b_FINAL.png")

    saved_pdf = safe_savefig(fig, pdf_path, format='pdf', bbox_inches='tight')
    saved_png = safe_savefig(fig, png_path, dpi=300, bbox_inches='tight')

    plt.close(fig)

    print(f"Gotowe! Wygenerowano {len(df_rho)} punktów korelacji.")
    print(f"Zapisano PDF: {saved_pdf}")
    print(f"Zapisano PNG: {saved_png}")

if __name__ == "__main__":
    main()