import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations

# ================= KONFIGURACJA ŚCIEŻEK =================
BASE_DIR = r"D:\Studia\Projekty\gsp-lab\results"
GSP_FILE = os.path.join(BASE_DIR, "plots", "gGSP_5.0.csv")
RMSD_DIR = os.path.join(BASE_DIR, "other-metrics", "rnapuzzles.github.io")
LDDT_FILE = os.path.join(BASE_DIR, "lddt", "finished", "global_report.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "plots")
os.makedirs(OUTPUT_DIR, exist_ok=True)

EXCLUDED_TARGETS = ['PZ10tRNA', 'PZ26tRNA', 'PZ27tRNA', 'PZ28tRNA']

# ================= FUNKCJE POMOCNICZE =================
def split_modelgroup(mg):
    parts = str(mg).split('_')
    lab = '_'.join(parts[:-1]) if len(parts) >= 2 else str(mg)
    num = parts[-1] if len(parts) >= 2 else '0'
    return lab.lower().strip(), num.strip()

def normalize_series(s):
    if s.max() == s.min(): return s * 0.0
    return (s - s.min()) / (s.max() - s.min())

def load_all_data():
    gsp = pd.read_csv(GSP_FILE)
    gsp = gsp[~gsp['Target'].isin(EXCLUDED_TARGETS)]
    gsp[['Lab', 'Num']] = gsp['ModelGroup'].apply(lambda x: pd.Series(split_modelgroup(x)))
    
    all_data = []
    for target in gsp['Target'].unique():
        t_metrics = gsp[gsp['Target'] == target].copy()[['Lab', 'Num', 'gGSP', 'Target']]
        target_dir = os.path.join(RMSD_DIR, target)
        if not os.path.exists(target_dir): continue

        for m_name in ['RMSD', 'TM-score', 'INF_all']:
            fpath = os.path.join(target_dir, f"{m_name}.csv")
            if os.path.exists(fpath):
                df_m = pd.read_csv(fpath)
                df_m.columns = [c.capitalize() if c.lower() in ['lab', 'num'] else c for c in df_m.columns]
                df_m['Lab'] = df_m['Lab'].astype(str).str.lower().str.strip()
                df_m['Num'] = df_m['Num'].astype(str).str.strip()
                val_col = [c for c in df_m.columns if c not in ['Lab', 'Num']][0]
                t_metrics = pd.merge(t_metrics, df_m[['Lab', 'Num', val_col]].rename(columns={val_col: m_name}), on=['Lab', 'Num'], how='inner')

        if os.path.exists(LDDT_FILE):
            lddt = pd.read_csv(LDDT_FILE).rename(columns={'Puzzle': 'Target'})
            lddt = lddt[lddt['Target'] == target]
            lddt[['Lab', 'Num']] = lddt['ModelGroup'].apply(lambda x: pd.Series(split_modelgroup(x)))
            t_metrics = pd.merge(t_metrics, lddt[['Lab', 'Num', 'lDDT']], on=['Lab', 'Num'], how='inner')
        
        if len(t_metrics) > 1:
            all_data.append(t_metrics)
    return pd.concat(all_data, ignore_index=True)

# ================= GENEROWANIE WYKRESU I ANALIZA =================
def generate_anomaly_report(df):
    # Ustawienia stylu zgodne z Twoim kodem
    sns.set_theme(style="whitegrid")
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'pdf.fonttype': 42,
        'axes.labelweight': 'bold'
    })

    # Kolory (używam tych z Twojego słownika RNA_COLORS dla spójności)
    COLOR_M1 = '#66c2a5' # Aptamer green
    COLOR_M2 = '#fc8d62' # Riboswitch orange

    # Normalizacja dla wyboru najlepszych przykładów
    for col in ['gGSP', 'RMSD', 'TM-score', 'INF_all', 'lDDT']:
        df[f'{col}_norm'] = df.groupby('Target')[col].transform(normalize_series)

    all_anomalies = []
    for target, group in df.groupby('Target'):
        models = group.reset_index(drop=True)
        target_pairs = []
        for i, j in combinations(range(len(models)), 2):
            m1, m2 = models.iloc[i], models.iloc[j]
            d_gsp = abs(m1['gGSP_norm'] - m2['gGSP_norm'])
            # Średnia różnica metryk globalnych
            d_glob = np.mean([abs(m1[f'{c}_norm'] - m2[f'{c}_norm']) for c in ['RMSD', 'TM-score', 'INF_all', 'lDDT']])
            
            score = d_gsp / (d_glob + 0.001)
            target_pairs.append({
                'Target': target, 'M1': f"{m1['Lab']}_{m1['Num']}", 'M2': f"{m2['Lab']}_{m2['Num']}",
                'score': score, 'm1_data': m1, 'm2_data': m2
            })
        if target_pairs:
            all_anomalies.append(max(target_pairs, key=lambda x: x['score']))

    # Najlepsze 5 do konsoli
    top_5 = sorted(all_anomalies, key=lambda x: x['score'], reverse=True)[:5]
    print("\n" + "="*80)
    print("TOP 5 ANOMALIES (CONSISTENT FOLD VS LOCAL DIFFERENCE)")
    print("-" * 80)
    for i, ex in enumerate(top_5):
        print(f"{i+1}. {ex['Target']:<10} | {ex['M1']} vs {ex['M2']:<25} | Score: {ex['score']:.2f}")
    print("="*80)

    # --- GENEROWANIE WYKRESU DLA 1 NAJLEPSZEGO PRZYKŁADU ---
    best = top_5[0]
    metrics = ['gGSP', 'RMSD', 'TM-score', 'INF_all', 'lDDT']
    units = ['(0-100)', '[Å]', '(0-1)', '(0-1)', '(0-1)']
    
    fig, axes = plt.subplots(1, 5, figsize=(16, 5))
    fig.suptitle(f'Comparative Analysis: Target {best["Target"]} - {best["M1"]} vs {best["M2"]}', fontsize=16, fontweight='bold', y=1.08)

    for idx, m in enumerate(metrics):
        ax = axes[idx]
        val1 = best['m1_data'][m]
        val2 = best['m2_data'][m]
        
        bars = ax.bar(['Model 1', 'Model 2'], [val1, val2], color=[COLOR_M1, COLOR_M2], edgecolor='black', alpha=0.85, width=0.6)
        
        # OSOBNE OSI DLA KAŻDEGO SUBPLOTU
        ax.set_title(f"{m}\n{units[idx]}", fontsize=12, fontweight='bold', pad=10)
        
        # Dynamiczne granice osi Y dla czytelności
        if m == 'gGSP': 
            ax.set_ylim(0, 100)
        elif m == 'RMSD': 
            ax.set_ylim(0, max(val1, val2) * 1.3)
        else: 
            ax.set_ylim(0, 1.0)

        # Wartości nad słupkami
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., h + (ax.get_ylim()[1]*0.02),
                    f'{h:.2f}', ha='center', va='bottom', fontweight='bold', fontsize=11)
        
        # Etykiety osi X (skrócone dla estetyki)
        ax.set_xticklabels([best['M1'].split('_')[0], best['M2'].split('_')[0]], fontsize=9)

    plt.tight_layout()
    # Zapis zgodny z Twoim stylem
    plt.savefig(os.path.join(OUTPUT_DIR, f"Anomaly_Detail_{best['Target']}.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, f"Anomaly_Detail_{best['Target']}.pdf"), format='pdf', bbox_inches='tight')
    print(f"\nWykres szczegółowy zapisano w: {OUTPUT_DIR}")
    plt.show()

if __name__ == "__main__":
    try:
        data = load_all_data()
        generate_anomaly_report(data)
    except Exception as e:
        import traceback
        print(f"Błąd: {e}")
        traceback.print_exc()