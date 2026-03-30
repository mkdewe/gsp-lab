import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.metrics import roc_curve, roc_auc_score, auc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from matplotlib.patches import Patch

# ================= KONFIGURACJA ŚCIEŻEK =================
BASE_DIR = r"D:\Studia\Projekty\gsp-lab\results"
GSP_FILE = os.path.join(BASE_DIR, "plots", "gGSP_5.0.csv")
RMSD_DIR = os.path.join(BASE_DIR, "other-metrics", "rnapuzzles.github.io")
OUTPUT_DIR = os.path.join(BASE_DIR, "plots")

# Wykluczone targety
EXCLUDED_TARGETS = ['PZ10tRNA', 'PZ26tRNA', 'PZ27tRNA', 'PZ28tRNA']

# Granica "dobrego" modelu
GROUND_TRUTH_THRESHOLD = 5.0

# Czy większy gGSP = lepszy model?
HIGHER_GGSP_IS_BETTER = True

# Debug
DEBUG = True
BOOTSTRAP_N = 1000
THRESHOLD_SWEEP = np.arange(3.0, 7.1, 0.5)

# Estetyka
sns.set_theme(style="whitegrid")
plt.rcParams.update({
    'font.family': 'sans-serif',
    'pdf.fonttype': 42
})


def split_modelgroup(mg):
    parts = str(mg).split('_')
    lab = '_'.join(parts[:-1]) if len(parts) >= 2 else str(mg)
    num = parts[-1] if len(parts) >= 2 else '0'
    return str(lab).lower().strip(), str(num).strip()


def load_and_merge():
    if not os.path.exists(GSP_FILE):
        raise FileNotFoundError(f"Brak pliku: {GSP_FILE}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df_gsp = pd.read_csv(GSP_FILE)

    required_cols = {'Target', 'ModelGroup', 'gGSP'}
    missing = required_cols - set(df_gsp.columns)
    if missing:
        raise ValueError(f"Brakuje kolumn: {missing}")

    df_gsp = df_gsp[~df_gsp['Target'].isin(EXCLUDED_TARGETS)].copy()

    df_gsp[['Lab', 'Num']] = df_gsp['ModelGroup'].apply(
        lambda x: pd.Series(split_modelgroup(x))
    )

    df_gsp['gGSP'] = pd.to_numeric(df_gsp['gGSP'], errors='coerce')
    df_gsp = df_gsp.dropna(subset=['gGSP']).copy()

    all_merged = []

    for target in df_gsp['Target'].dropna().unique():
        rmsd_path = os.path.join(RMSD_DIR, target, "RMSD.csv")
        if not os.path.exists(rmsd_path):
            if DEBUG:
                print(f"[DEBUG] Brak RMSD.csv dla targetu: {target}")
            continue

        df_rmsd = pd.read_csv(rmsd_path)

        df_rmsd.columns = [
            c.capitalize() if c.lower() in ['lab', 'num'] else c
            for c in df_rmsd.columns
        ]

        if 'Lab' not in df_rmsd.columns or 'Num' not in df_rmsd.columns:
            if DEBUG:
                print(f"[DEBUG] Pomijam {target} — brak kolumn Lab/Num")
            continue

        val_candidates = [c for c in df_rmsd.columns if c not in ['Lab', 'Num']]
        if not val_candidates:
            if DEBUG:
                print(f"[DEBUG] Pomijam {target} — brak kolumny z RMSD")
            continue

        val_col = val_candidates[0]

        df_rmsd['Lab'] = df_rmsd['Lab'].astype(str).str.lower().str.strip()
        df_rmsd['Num'] = df_rmsd['Num'].astype(str).str.strip()

        df_rmsd[val_col] = pd.to_numeric(df_rmsd[val_col], errors='coerce')
        df_rmsd = df_rmsd.dropna(subset=[val_col]).copy()

        target_rmsd = df_rmsd[['Lab', 'Num', val_col]].rename(columns={val_col: 'RMSD'})

        merged = pd.merge(
            df_gsp[df_gsp['Target'] == target],
            target_rmsd,
            on=['Lab', 'Num'],
            how='inner'
        )

        if DEBUG:
            print(f"[DEBUG] {target}: gGSP={len(df_gsp[df_gsp['Target'] == target])}, "
                  f"RMSD={len(target_rmsd)}, merged={len(merged)}")

        if not merged.empty:
            all_merged.append(merged)

    if not all_merged:
        raise ValueError("Brak połączonych danych.")

    df = pd.concat(all_merged, ignore_index=True)
    df = df.dropna(subset=['gGSP', 'RMSD']).copy()
    return df


def manual_auc_from_ranks(y_true, y_scores):
    y_true = np.asarray(y_true).astype(int)
    y_scores = np.asarray(y_scores, dtype=float)

    n_pos = int(y_true.sum())
    n_neg = int(len(y_true) - n_pos)

    if n_pos == 0 or n_neg == 0:
        return np.nan

    ranks = pd.Series(y_scores).rank(method='average').to_numpy()
    sum_pos_ranks = ranks[y_true == 1].sum()

    auc_manual = (sum_pos_ranks - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    return float(auc_manual)


def bootstrap_auc_ci(y_true, y_scores, n_boot=1000, seed=42):
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true)
    y_scores = np.asarray(y_scores)

    values = []
    n = len(y_true)

    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yt = y_true[idx]
        ys = y_scores[idx]

        if len(np.unique(yt)) < 2:
            continue

        values.append(roc_auc_score(yt, ys))

    if not values:
        return np.nan, np.nan

    lo = float(np.percentile(values, 2.5))
    hi = float(np.percentile(values, 97.5))
    return lo, hi


def debug_auc_report(df, y_true, y_scores, auc_value):
    fpr, tpr, thresholds = roc_curve(y_true, y_scores)
    auc_from_curve = auc(fpr, tpr)
    auc_manual = manual_auc_from_ranks(y_true, y_scores)
    ci_lo, ci_hi = bootstrap_auc_ci(y_true, y_scores, n_boot=BOOTSTRAP_N)

    pos_mask = y_true == 1
    neg_mask = y_true == 0

    print("\n================ AUC DEBUG ================")
    print(f"Próg RMSD: {GROUND_TRUTH_THRESHOLD}")
    print(f"Liczba próbek: {len(df)}")
    print(f"Pozytywne (RMSD < {GROUND_TRUTH_THRESHOLD}): {int(pos_mask.sum())}")
    print(f"Negatywne (RMSD >= {GROUND_TRUTH_THRESHOLD}): {int(neg_mask.sum())}")

    print("\n[gGSP statystyki]")
    print(f"Całość:   min={df['gGSP'].min():.6f}, max={df['gGSP'].max():.6f}, mean={df['gGSP'].mean():.6f}, std={df['gGSP'].std():.6f}")
    print(f"Pozytywne: min={df.loc[pos_mask, 'gGSP'].min():.6f}, max={df.loc[pos_mask, 'gGSP'].max():.6f}, mean={df.loc[pos_mask, 'gGSP'].mean():.6f}, std={df.loc[pos_mask, 'gGSP'].std():.6f}")
    print(f"Negatywne: min={df.loc[neg_mask, 'gGSP'].min():.6f}, max={df.loc[neg_mask, 'gGSP'].max():.6f}, mean={df.loc[neg_mask, 'gGSP'].mean():.6f}, std={df.loc[neg_mask, 'gGSP'].std():.6f}")

    print("\n[AUC porównanie metod]")
    print(f"roc_auc_score: {auc_value:.6f}")
    print(f"auc(fpr, tpr):  {auc_from_curve:.6f}")
    print(f"AUC ręczne:     {auc_manual:.6f}")
    print(f"Różnica 1-2:    {abs(auc_value - auc_from_curve):.12f}")
    print(f"Różnica 1-3:    {abs(auc_value - auc_manual):.12f}")

    print("\n[ROC próbki]")
    print("Pierwsze 10 punktów ROC:")
    for i in range(min(10, len(fpr))):
        thr = thresholds[i] if i < len(thresholds) else np.nan
        print(f"  i={i:02d}  FPR={fpr[i]:.6f}  TPR={tpr[i]:.6f}  threshold={thr:.6f}")

    if not np.isnan(ci_lo):
        print("\n[Bootstrap]")
        print(f"95% CI dla AUC: [{ci_lo:.6f}, {ci_hi:.6f}]")

    print("==========================================\n")


def threshold_sweep(df, thresholds):
    print("\n================ THRESHOLD SWEEP ================")
    print("Threshold | Positives | Negatives | AUC")
    print("-----------------------------------------------")

    y_scores = df['gGSP'].astype(float).to_numpy()
    if not HIGHER_GGSP_IS_BETTER:
        y_scores = -y_scores

    for thr in thresholds:
        y_true = (df['RMSD'] < thr).astype(int).to_numpy()
        n_pos = int(y_true.sum())
        n_neg = int(len(y_true) - n_pos)

        if n_pos == 0 or n_neg == 0:
            print(f"{thr:8.2f} | {n_pos:9d} | {n_neg:9d} |  nan")
            continue

        auc_val = roc_auc_score(y_true, y_scores)
        print(f"{thr:8.2f} | {n_pos:9d} | {n_neg:9d} | {auc_val:.6f}")

    print("=================================================\n")


def prepare_scores(df):
    df = df.copy()
    y_true = (df['RMSD'] < GROUND_TRUTH_THRESHOLD).astype(int).to_numpy()

    y_scores = df['gGSP'].astype(float).to_numpy()
    if not HIGHER_GGSP_IS_BETTER:
        y_scores = -y_scores

    auc_value = roc_auc_score(y_true, y_scores)

    df['Status'] = np.where(
        df['RMSD'] < GROUND_TRUTH_THRESHOLD,
        'RMSD < 5 Å',
        'RMSD ≥ 5 Å'
    )

    return df, y_true, y_scores, auc_value


def plot_boxplot(df, auc_value):
    plt.figure(figsize=(8, 5))

    palette = {
        'RMSD < 5 Å': '#66c2a5',
        'RMSD ≥ 5 Å': '#fc8d62'
    }

    ax = sns.boxplot(
        data=df,
        x='Status',
        y='gGSP',
        hue='Status',
        palette=palette,
        dodge=False,
        width=0.5,
        showfliers=False
    )

    sns.stripplot(
        data=df,
        x='Status',
        y='gGSP',
        color='black',
        size=2,
        alpha=0.2,
        jitter=True,
        ax=ax
    )

    if ax.legend_ is not None:
        ax.legend_.remove()

    legend_handles = [
        Patch(facecolor=palette['RMSD < 5 Å'], edgecolor='black', label='RMSD < 5 Å'),
        Patch(facecolor=palette['RMSD ≥ 5 Å'], edgecolor='black', label='RMSD ≥ 5 Å')
    ]
    ax.legend(
        handles=legend_handles,
        title='Model class',
        loc='upper left',
        frameon=True
    )

    plt.title(
        f'gGSP score distribution by model quality\nAUC = {auc_value:.3f}',
        fontsize=13,
        fontweight='bold'
    )
    plt.ylabel('gGSP score', fontsize=11)
    plt.xlabel('')

    save_path = os.path.join(OUTPUT_DIR, "gGSP_distribution_boxplot.pdf")
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()

    print(f"Boxplot zapisany: {save_path}")


def plot_roc_curve_zoom(df, y_true, y_scores, auc_value):
    fpr, tpr, _ = roc_curve(y_true, y_scores)

    n_pos = int(np.sum(y_true == 1))
    n_neg = int(np.sum(y_true == 0))

    fig, ax = plt.subplots(figsize=(6.5, 6.5))

    ax.plot(fpr, tpr, linewidth=2, label=f'AUC = {auc_value:.3f}')
    ax.plot([0, 1], [0, 1], linestyle='--', linewidth=1, color='gray')

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel('False Positive Rate', fontsize=11)
    ax.set_ylabel('True Positive Rate', fontsize=11)
    ax.set_title('ROC curve for gGSP score', fontsize=13, fontweight='bold')
    ax.legend(loc='lower right')

    # Liczby klas na wykresie
    stats_text = (
        f'Positive (RMSD < {GROUND_TRUTH_THRESHOLD:.1f} Å): {n_pos}\n'
        f'Negative (RMSD ≥ {GROUND_TRUTH_THRESHOLD:.1f} Å): {n_neg}'
    )
    ax.text(
        0.02, 0.18,
        stats_text,
        transform=ax.transAxes,
        fontsize=9,
        va='bottom',
        ha='left',
        bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.9)
    )

    # Wstawka z przybliżeniem w okolicy początku układu
    axins = inset_axes(ax, width="45%", height="45%", loc='lower right', borderpad=1.2)
    axins.plot(fpr, tpr, linewidth=2)
    axins.plot([0, 1], [0, 1], linestyle='--', linewidth=1, color='gray')

    # Przybliżenie tam, gdzie widać najwięcej różnicy
    axins.set_xlim(0.0, 0.05)
    axins.set_ylim(0.0, 0.60)
    axins.set_xticks([0.0, 0.02, 0.04])
    axins.set_yticks([0.0, 0.2, 0.4, 0.6])
    axins.tick_params(labelsize=8)

    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

    save_path = os.path.join(OUTPUT_DIR, "gGSP_ROC_curve_zoom.pdf")
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()

    print(f"ROC zapisany: {save_path}")


if __name__ == "__main__":
    try:
        data = load_and_merge()
        data, y_true, y_scores, auc_value = prepare_scores(data)

        print(f"MODELI: {len(data)}")
        print(f"AUC: {auc_value:.6f}")

        if DEBUG:
            debug_auc_report(data, y_true, y_scores, auc_value)
            threshold_sweep(data, THRESHOLD_SWEEP)

        plot_boxplot(data, auc_value)
        plot_roc_curve_zoom(data, y_true, y_scores, auc_value)

    except Exception as e:
        print(f"BŁĄD: {e}")