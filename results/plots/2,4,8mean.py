import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

# ================= KONFIGURACJA ŚCIEŻEK =================
BASE_DIR = r"D:\Studia\Projekty\gsp-lab\results"
GSP_DIR = os.path.join(BASE_DIR, "plots")
OUTPUT_DIR = GSP_DIR
CSV_MAPPING_FILE = os.path.join(GSP_DIR, 'RNA-types.csv')
os.makedirs(OUTPUT_DIR, exist_ok=True)

EXCLUDED_TARGETS = ['PZ10tRNA', 'PZ26tRNA', 'PZ27tRNA', 'PZ28tRNA']

# ================= SPÓJNA KOLORYSTYKA (Project-wide) =================
# Kolory dla metryk porównawczych
METRIC_COLORS = {
    'GSP(5)': '#66c2a5',    # Ten sam kolor co Aptamer w drugim skrypcie (możesz zmienić jeśli wolisz unikalny)
    'GSP_avg': '#fc8d62',   # Ten sam kolor co Riboswitch
    'Difference': '#8da0cb' # Ten sam kolor co Ribozyme
}

# Predefiniowane kolory dla typów RNA (identyczne z Twoim drugim skryptem)
RNA_COLORS = {
    'Aptamer': '#66c2a5',
    'Riboswitch': '#fc8d62',
    'Ribozyme': '#8da0cb',
    'Viral element': '#e78ac3',
    'Regulatory element': '#ffd92f',
    'Synthetic construct': '#e5c494',
    'RNA-binding protein': '#b3b3b3',
    'Nano construct': '#8dd3c7',
    'Other RNA Structures': '#a6d854',
    'Unknown': '#d3d3d3'
}

# Estetyka (identyczna jak w Figure2/3)
sns.set_theme(style="whitegrid")
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'axes.labelweight': 'bold',
    'axes.titleweight': 'bold'
})

# ================= FUNKCJE POMOCNICZE =================
def split_modelgroup(mg):
    parts = mg.split('_')
    lab = '_'.join(parts[:-1]) if len(parts) >= 2 else mg
    num = parts[-1] if len(parts) >= 2 else '0'
    return lab.lower().strip(), str(num).strip()

def load_gsp(threshold):
    filename = f"gGSP_{threshold}.csv"
    filepath = os.path.join(GSP_DIR, filename)
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Brak pliku: {filepath}")
    df = pd.read_csv(filepath)
    df = df[~df['Target'].isin(EXCLUDED_TARGETS)]
    df[['Lab', 'Num']] = df['ModelGroup'].apply(lambda x: pd.Series(split_modelgroup(x)))
    return df[['Target', 'Lab', 'Num', 'gGSP']].rename(columns={'gGSP': f'gGSP_{threshold}'})

def safe_savefig(fig, filename):
    path = os.path.join(OUTPUT_DIR, filename)
    try:
        fig.savefig(path, bbox_inches='tight', dpi=300)
    except PermissionError:
        fig.savefig(path.replace(".pdf", "_new.pdf").replace(".png", "_new.png"), bbox_inches='tight', dpi=300)

# ================= GŁÓWNE PRZETWARZANIE =================
def main():
    thresholds = [2.0, 4.0, 5.0, 8.0]
    gsp_data = {}

    for th in thresholds:
        gsp_data[th] = load_gsp(th)

    # Połączenie danych
    merged = gsp_data[2.0]
    for th in [4.0, 5.0, 8.0]:
        merged = pd.merge(merged, gsp_data[th], on=['Target', 'Lab', 'Num'], how='inner')

    # Obliczenie średniej kroczącej (2,4,8)
    merged['gGSP_avg'] = (merged['gGSP_2.0'] + merged['gGSP_4.0'] + merged['gGSP_8.0']) / 3.0
    diff = merged['gGSP_5.0'] - merged['gGSP_avg']

    # 1. Scatter plot (GSP 5 vs Avg)
    fig1, ax1 = plt.subplots(figsize=(7, 7))
    ax1.scatter(merged['gGSP_5.0'], merged['gGSP_avg'],
                alpha=0.2, s=10, color=METRIC_COLORS['GSP(5)'], edgecolors='none')
    ax1.plot([0, 100], [0, 100], color='#e74c3c', linestyle='--', linewidth=1.5, label='y = x')
    ax1.set_xlabel('GSP Value (Threshold 5.0 Å)')
    ax1.set_ylabel('Mean GSP (Thresholds 2, 4, 8 Å)')
    ax1.set_title('Metric Robustness: GSP(5.0) vs. Multi-threshold Mean')
    ax1.legend()
    safe_savefig(fig1, 'gGSP5_vs_avg_scatter.png')
    plt.close()

    # 2. Histogramy
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    sns.kdeplot(merged['gGSP_5.0'], fill=True, color=METRIC_COLORS['GSP(5)'], label='GSP(5.0)', ax=ax2)
    sns.kdeplot(merged['gGSP_avg'], fill=True, color=METRIC_COLORS['GSP_avg'], label='Mean(2,4,8)', ax=ax2)
    ax2.set_title('Global Distribution Comparison')
    ax2.set_xlabel('gGSP Value')
    ax2.legend()
    safe_savefig(fig2, 'gGSP5_vs_avg_hist.png')
    plt.close()

    # 3. Wykres słupkowy per Target (Posortowany)
    target_means = merged.groupby('Target')[['gGSP_5.0', 'gGSP_avg']].mean().reset_index()
    target_means = target_means.sort_values('gGSP_5.0', ascending=False)

    fig4, ax4 = plt.subplots(figsize=(14, 6))
    x = np.arange(len(target_means))
    width = 0.4
    ax4.bar(x - width/2, target_means['gGSP_5.0'], width, label='GSP(5.0)', color=METRIC_COLORS['GSP(5)'])
    ax4.bar(x + width/2, target_means['gGSP_avg'], width, label='Mean(2,4,8)', color=METRIC_COLORS['GSP_avg'])
    ax4.set_xticks(x)
    ax4.set_xticklabels(target_means['Target'], rotation=90, fontsize=9)
    ax4.set_ylabel('Average gGSP')
    ax4.set_title('Target-wise Stability of GSP Metric')
    ax4.legend()
    safe_savefig(fig4, 'gGSP5_vs_avg_target_bars.png')
    plt.close()

    # 4. SCATTER PER TARGET (Kluczowy dla spójności kolorów RNA)
    mapping_df = pd.read_csv(CSV_MAPPING_FILE)
    target_type = dict(zip(mapping_df['target_label'], mapping_df['type']))
    target_means['Type'] = target_means['Target'].map(target_type).fillna('Unknown')

    fig5, ax5 = plt.subplots(figsize=(10, 7))
    
    # Sortujemy typy, aby legenda była w stałej kolejności
    unique_types_in_data = sorted(target_means['Type'].unique())
    
    for rna_type in unique_types_in_data:
        subset = target_means[target_means['Type'] == rna_type]
        # Pobieramy kolor ze słownika RNA_COLORS, jeśli brak - szary
        current_color = RNA_COLORS.get(rna_type, RNA_COLORS['Unknown'])
        
        ax5.scatter(subset['gGSP_5.0'], subset['gGSP_avg'],
                    s=100, alpha=0.8, label=rna_type, 
                    color=current_color, edgecolors='white', linewidth=1.5)

    ax5.plot([0, 100], [0, 100], 'k--', alpha=0.3)
    ax5.set_xlim(0, 105); ax5.set_ylim(0, 105)
    ax5.set_xlabel('GSP(5.0) Target Mean')
    ax5.set_ylabel('GSP(2,4,8) Multi-threshold Mean')
    ax5.set_title('Correlation of Metrics by RNA Functional Category')
    
    # Legenda poza wykresem
    ax5.legend(title="RNA Category", bbox_to_anchor=(1.05, 1), loc='upper left')
    
    safe_savefig(fig5, 'gGSP5_vs_avg_target_scatter_colored.pdf')
    safe_savefig(fig5, 'gGSP5_vs_avg_target_scatter_colored.png')
    plt.close()

    print(f"Analiza zakończona. Wykresy zapisano w: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
