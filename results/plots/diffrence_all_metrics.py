import os
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from matplotlib.patches import Patch

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
    if s.max() == s.min():
        return s * 0.0
    return (s - s.min()) / (s.max() - s.min())

def safe_savefig(fig, path, dpi=200):
    """
    Zapisuje figurę. Jeśli plik jest zablokowany albo zapis się nie uda,
    zapisuje do alternatywnej nazwy z timestampem.
    """
    try:
        fig.savefig(path, dpi=dpi, bbox_inches='tight')
        return path
    except PermissionError:
        base, ext = os.path.splitext(path)
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        alt_path = f"{base}_{ts}{ext}"
        fig.savefig(alt_path, dpi=dpi, bbox_inches='tight')
        return alt_path

# ================= WCZYTYWANIE DANYCH =================
def load_all_data():
    gsp = pd.read_csv(GSP_FILE)
    gsp = gsp[~gsp['Target'].isin(EXCLUDED_TARGETS)].copy()
    gsp[['Lab', 'Num']] = gsp['ModelGroup'].apply(lambda x: pd.Series(split_modelgroup(x)))

    all_data = []

    for target in gsp['Target'].unique():
        t_metrics = gsp[gsp['Target'] == target].copy()[['Lab', 'Num', 'gGSP', 'Target']]
        target_dir = os.path.join(RMSD_DIR, target)

        if not os.path.exists(target_dir):
            continue

        for m_name in ['RMSD', 'TM-score', 'INF_all']:
            fpath = os.path.join(target_dir, f"{m_name}.csv")
            if os.path.exists(fpath):
                df_m = pd.read_csv(fpath)
                df_m.columns = [c.capitalize() if c.lower() in ['lab', 'num'] else c for c in df_m.columns]

                if 'Lab' not in df_m.columns or 'Num' not in df_m.columns:
                    continue

                df_m['Lab'] = df_m['Lab'].astype(str).str.lower().str.strip()
                df_m['Num'] = df_m['Num'].astype(str).str.strip()

                val_candidates = [c for c in df_m.columns if c not in ['Lab', 'Num']]
                if not val_candidates:
                    continue
                val_col = val_candidates[0]

                t_metrics = pd.merge(
                    t_metrics,
                    df_m[['Lab', 'Num', val_col]].rename(columns={val_col: m_name}),
                    on=['Lab', 'Num'],
                    how='inner'
                )

        if os.path.exists(LDDT_FILE):
            lddt = pd.read_csv(LDDT_FILE).rename(columns={'Puzzle': 'Target'})
            lddt = lddt[lddt['Target'] == target].copy()

            if len(lddt) > 0 and 'ModelGroup' in lddt.columns and 'lDDT' in lddt.columns:
                lddt[['Lab', 'Num']] = lddt['ModelGroup'].apply(lambda x: pd.Series(split_modelgroup(x)))
                t_metrics = pd.merge(
                    t_metrics,
                    lddt[['Lab', 'Num', 'lDDT']],
                    on=['Lab', 'Num'],
                    how='inner'
                )

        if len(t_metrics) > 1:
            all_data.append(t_metrics)

    if not all_data:
        raise ValueError("Nie udało się wczytać i połączyć danych.")

    return pd.concat(all_data, ignore_index=True)

# ================= GŁÓWNA FUNKCJA =================
def plot_best_example(df):
    sns.set_theme(style="whitegrid")
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'pdf.fonttype': 42,
        'axes.labelweight': 'bold'
    })

    COLOR_M1 = '#66c2a5'
    COLOR_M2 = '#fc8d62'

    metrics = ['gGSP', 'RMSD', 'TM-score', 'INF_all', 'lDDT']

    # Normalizacja do wyboru najlepszego przykładu
    for col in metrics:
        if col in df.columns:
            df[f'{col}_norm'] = df.groupby('Target')[col].transform(normalize_series)

    best_examples = []

    for target, group in df.groupby('Target'):
        models = group.reset_index(drop=True)

        best_score = -1
        best_pair = None

        norm_cols_global = [c for c in ['RMSD', 'TM-score', 'INF_all', 'lDDT'] if f'{c}_norm' in models.columns]
        if not norm_cols_global:
            continue

        for i, j in combinations(range(len(models)), 2):
            m1, m2 = models.iloc[i], models.iloc[j]

            d_gsp = abs(m1['gGSP_norm'] - m2['gGSP_norm'])
            d_global = np.mean([abs(m1[f'{c}_norm'] - m2[f'{c}_norm']) for c in norm_cols_global])

            score = d_gsp / (d_global + 0.001)

            if score > best_score:
                best_score = score
                best_pair = (m1, m2, score, target)

        if best_pair:
            best_examples.append(best_pair)

    if not best_examples:
        raise ValueError("Nie znaleziono żadnego sensownego przykładu do narysowania.")

    m1, m2, score, target = sorted(best_examples, key=lambda x: x[2], reverse=True)[0]

    label1 = f"{m1['Lab']}_{m1['Num']}"
    label2 = f"{m2['Lab']}_{m2['Num']}"

    print("\nBEST EXAMPLE:")
    print(f"Target: {target}")
    print(f"{label1} vs {label2}")
    print(f"Score: {score:.2f}")

    units = ['(0–100)', '[Å]', '(0–1)', '(0–1)', '(0–1)']

    fig, axes = plt.subplots(1, 5, figsize=(17, 4.8))

    fig.suptitle(
        f'Comparative analysis for RNA-Puzzles target: {target}',
        fontsize=14,
        fontweight='bold',
        y=1.05
    )

    legend_handles = [
        Patch(facecolor=COLOR_M1, edgecolor='black', label=label1),
        Patch(facecolor=COLOR_M2, edgecolor='black', label=label2)
    ]

    x_labels = ['Model A', 'Model B']

    for i, m in enumerate(metrics):
        ax = axes[i]

        val1 = float(m1[m])
        val2 = float(m2[m])

        ax.bar(
            x_labels,
            [val1, val2],
            color=[COLOR_M1, COLOR_M2],
            edgecolor='black',
            alpha=0.9,
            width=0.6
        )

        ax.set_title(f"{m}\n{units[i]}", fontsize=11, fontweight='bold', pad=10)

        if m == 'gGSP':
            ax.set_ylim(0, 100)
        elif m == 'RMSD':
            ax.set_ylim(0, max(val1, val2) * 1.3 if max(val1, val2) > 0 else 1.0)
        else:
            ax.set_ylim(0, 1.0)

        ax.yaxis.set_ticks_position('left')
        ax.yaxis.set_label_position('left')
        ax.tick_params(axis='x', rotation=0, labelsize=9)

    fig.legend(
        handles=legend_handles,
        loc='upper center',
        ncol=2,
        frameon=True,
        bbox_to_anchor=(0.5, 0.98)
    )

    plt.tight_layout(rect=[0, 0, 1, 0.92])

    pdf_path = os.path.join(OUTPUT_DIR, f"Best_Example_{target}.pdf")
    png_path = os.path.join(OUTPUT_DIR, f"Best_Example_{target}.png")

    saved_pdf = safe_savefig(fig, pdf_path, dpi=200)
    saved_png = safe_savefig(fig, png_path, dpi=300)

    plt.close(fig)

    print("\nZapisano:")
    print(f"- {saved_pdf}")
    print(f"- {saved_png}")

# ================= MAIN =================
if __name__ == "__main__":
    try:
        data = load_all_data()
        plot_best_example(data)
    except Exception as e:
        import traceback
        print(f"Błąd: {e}")
        traceback.print_exc()