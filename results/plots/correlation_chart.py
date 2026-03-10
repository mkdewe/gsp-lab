import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, FixedFormatter

# === Load data ===
csv_path = r"D:\Studia\Projekty\gsp-lab\results\correlations\summary_report_rounded.csv"
df = pd.read_csv(csv_path)

metrics = df["metric"].tolist()
model_counts = df["total_models"].tolist()

# Prepare tick labels
tick_labels = [f"{m}\n({n})" for m, n in zip(metrics, model_counts)]

# Correlation stats
min_p = df["min_pearson"].values
max_p = df["max_pearson"].values
avg_p = df["avg_pearson_coef"].values

min_s = df["min_spearman"].values
max_s = df["max_spearman"].values
avg_s = df["avg_spearman_coef"].values

min_k = df["min_kendall"].values
max_k = df["max_kendall"].values
avg_k = df["avg_kendall_coef"].values

# x positions
spacing = 2.0
x_centers = np.arange(len(metrics)) * spacing
offset = 0.5

# colors
col_p = "#1f77b4"
col_s = "#ff7f0e"
col_k = "#2ca02c"

fig, ax = plt.subplots(figsize=(10, 4))

# plot ranges and means
for i, xc in enumerate(x_centers):
    # Pearson
    ax.vlines(xc - offset, min_p[i], max_p[i], color=col_p, linewidth=2)
    ax.scatter(xc - offset, avg_p[i], color=col_p, s=50)
    # Spearman
    ax.vlines(xc, min_s[i], max_s[i], color=col_s, linewidth=2)
    ax.scatter(xc, avg_s[i], color=col_s, s=50)
    # Kendall
    ax.vlines(xc + offset, min_k[i], max_k[i], color=col_k, linewidth=2)
    ax.scatter(xc + offset, avg_k[i], color=col_k, s=50)

# fixed x ticks
ax.xaxis.set_major_locator(FixedLocator(x_centers))
ax.xaxis.set_major_formatter(FixedFormatter(tick_labels))

# labels + grid
ax.set_ylabel("Correlation coefficient", fontsize=10)
ax.set_title("Range & Average of Correlation Types per Metric", fontsize=12)
ax.grid(axis='y', linestyle='--', alpha=0.3)

# rotate and center labels
for lbl in ax.get_xticklabels():
    lbl.set_rotation(45)
    lbl.set_ha('center')
    lbl.set_fontsize(8)

# legend in lower left
ax.legend(
    handles=[
        ax.scatter([], [], color=col_p, s=50),
        ax.scatter([], [], color=col_s, s=50),
        ax.scatter([], [], color=col_k, s=50),
    ],
    labels=["Pearson", "Spearman", "Kendall"],
    loc="lower left",
    fontsize=9,
    frameon=True
)

# adjust margins
fig.subplots_adjust(
    left=0.14,
    right=0.98,
    top=0.90,
    bottom=0.20
)

# save + show
plt.savefig("stable_correlation_ranges_leftlegend.png", dpi=200)
plt.show()

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# -----------------------------
# Ustawienia użytkownika
# -----------------------------
folder_path = r"D:\Studia\Projekty\gsp-lab\results\correlations"  # <- dostosuj jeśli trzeba
out_prefix = "results_pearson_spearman"  # prefix dla zapisanych plików
bootstrap_samples = 2000
random_seed = 0
# -----------------------------

# pomocniczna funkcja: znajdź pierwszą pasującą nazwę kolumny
def find_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

# kandydaci na nazwy kolumn (robi skrypt bardziej odpornym na różnice w CSV)
metric_cols = ["metric", "metric_name", "metric_id"]
pearson_cols = ["pearson", "pearson_coef", "pearson_coefficient", "pearson_correlation"]
spearman_cols = ["spearman", "spearman_coef", "spearman_coefficient", "spearman_correlation"]
kendall_cols = ["kendall", "kendall_coef", "kendall_coefficient", "kendall_correlation"]

# znajdź pliki
files = glob.glob(os.path.join(folder_path, "*_correlations.csv"))
if not files:
    raise SystemExit(f"Nie znaleziono plików *_correlations.csv w {folder_path}")

all_frames = []
for fpath in files:
    try:
        df = pd.read_csv(fpath)
    except Exception as e:
        print(f"Warning: nie można wczytać {fpath}: {e}")
        continue

    # wyciągnij nazwę targetu z pliku
    target_name = os.path.basename(fpath).replace("_correlations.csv", "")
    df["target"] = target_name

    # sprawdź podstawowe kolumny; jeśli brak, spróbuj zgadnąć (lub pominąć plik)
    col_metric = find_col(df, metric_cols)
    col_p = find_col(df, pearson_cols)
    col_s = find_col(df, spearman_cols)
    col_k = find_col(df, kendall_cols)

    if col_metric is None or col_p is None or col_s is None:
        print(f"Skipping {fpath} — brak wymaganych kolumn (metric/pearson/spearman).")
        continue

    # ujednolicenie nazw
    df = df.rename(columns={col_metric: "metric", col_p: "pearson", col_s: "spearman"})
    if col_k is not None:
        df = df.rename(columns={col_k: "kendall"})
    # keep only columns we need
    all_frames.append(df[["target", "metric", "pearson", "spearman"]].copy())

# concat
df_all = pd.concat(all_frames, ignore_index=True)
# konwersja na numeric i drop NaN
df_all["pearson"] = pd.to_numeric(df_all["pearson"], errors="coerce")
df_all["spearman"] = pd.to_numeric(df_all["spearman"], errors="coerce")
df_all = df_all.dropna(subset=["metric", "pearson", "spearman"]).reset_index(drop=True)

print(f"Wczytano {len(files)} plików, w DataFrame {df_all.shape[0]} wierszy (per-target).")
print("Unikalne metryki:", df_all["metric"].unique())

# -----------------------------
# Scatter: wszystkie punkty (pearson vs spearman), kolor = metric
# + globalna regresja, R^2, 95% CI (bootstrap)
# + centroidy metryk (średnie) z etykietami
# -----------------------------
x = df_all["pearson"].values
y = df_all["spearman"].values

# globalna regresja (scipy)
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
r2 = r_value**2
print(f"Global linear fit: slope={slope:.4f}, intercept={intercept:.4f}, R^2={r2:.3f}, p={p_value:.2e}")

# bootstrap predykcji by zbudować 95% CI dla linii regresji
rng = np.random.default_rng(random_seed)
x_line = np.linspace(x.min() - 0.01, x.max() + 0.01, 300)
pred_boot = np.empty((bootstrap_samples, x_line.size))

n = len(x)
for i in range(bootstrap_samples):
    idx = rng.integers(0, n, n)
    xb = x[idx]
    yb = y[idx]
    coef = np.polyfit(xb, yb, 1)
    pred_boot[i, :] = np.polyval(coef, x_line)

lower = np.percentile(pred_boot, 2.5, axis=0)
upper = np.percentile(pred_boot, 97.5, axis=0)
# punktowa przewidywana linia z oryginalnego dopasowania
y_line = intercept + slope * x_line

# przygotuj kolory dla metryk
metrics = sorted(df_all["metric"].unique())
cmap = plt.get_cmap("tab20")
colors = {m: cmap(i % 20) for i, m in enumerate(metrics)}

plt.figure(figsize=(8,8))
for m in metrics:
    sub = df_all[df_all["metric"] == m]
    plt.scatter(sub["pearson"], sub["spearman"],
                s=25, alpha=0.45, label=m, color=colors[m], edgecolor='none')

# centroidy metryk (średnie) — większe punkty i etykiety
centroids = df_all.groupby("metric")[["pearson", "spearman"]].mean().reset_index()
for _, row in centroids.iterrows():
    plt.scatter(row["pearson"], row["spearman"], s=120, facecolors='none', edgecolors='k', linewidth=0.9)
    plt.annotate(row["metric"], (row["pearson"], row["spearman"]),
                 xytext=(6, 4), textcoords='offset points', fontsize=8, weight='bold')

# linia regresji i CI
plt.plot(x_line, y_line, color='black', linewidth=1.8, label="global linear fit")
plt.fill_between(x_line, lower, upper, color='grey', alpha=0.25, label="95% CI (bootstrap)")

plt.xlabel("Pearson (per-target)")
plt.ylabel("Spearman (per-target)")
plt.title("Pearson vs Spearman — all targets")
plt.grid(alpha=0.3, linestyle='--')
plt.text(0.02, 0.98, f"R² = {r2:.3f}\nreg p = {p_value:.2e}\nN = {n}", transform=plt.gca().transAxes,
         va='top', fontsize=9, bbox=dict(boxstyle='round', facecolor='white', alpha=0.6))
plt.legend(loc='lower right', fontsize=7, framealpha=0.9, ncol=1)

plt.tight_layout()
plt.savefig(f"{out_prefix}_per_target_scatter.png", dpi=200)
plt.show()

# -----------------------------
# Boxplot: rozkłady Pearson per metric
# -----------------------------
plt.figure(figsize=(12,6))
data_for_box = [df_all[df_all["metric"] == m]["pearson"].values for m in metrics]
labels = [f"{m}\n(n={len(df_all[df_all['metric']==m])})" for m in metrics]

bp = plt.boxplot(data_for_box, labels=labels, patch_artist=True, showfliers=True)
# minimalistyczne style
for patch in bp['boxes']:
    patch.set_facecolor('white')
    patch.set_edgecolor('black')
for whisker in bp['whiskers']:
    whisker.set_color('black')
for cap in bp['caps']:
    cap.set_color('black')
for median in bp['medians']:
    median.set_color('red')

plt.xticks(rotation=45, ha='right', fontsize=8)
plt.ylabel("Pearson (per-target)")
plt.title("Rozkłady Pearson per metric (boxplots)")
plt.grid(axis='y', linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig(f"{out_prefix}_pearson_boxplots.png", dpi=200)
plt.show()

# -----------------------------
# Outliers: największe reszty względem globalnej linii
# -----------------------------
# oblicz residual = y - (intercept + slope*x)
df_all["residual_abs"] = np.abs(df_all["spearman"] - (intercept + slope * df_all["pearson"]))
top_outliers = df_all.sort_values("residual_abs", ascending=False).head(12)
print("\nTop 12 punktów odchylających się od globalnej regresji (largest absolute residuals):")
print(top_outliers[["target", "metric", "pearson", "spearman", "residual_abs"]].to_string(index=False))

# dodatkowo: statystyki per-metric (mediana, IQR, share above 0.5)
summary = df_all.groupby("metric")["pearson"].agg(['count','median','mean','std'])
summary["frac_above_0.5"] = df_all.groupby("metric").apply(lambda g: (g["pearson"]>0.5).mean())
print("\nSummary per metric (count, median, mean, std, frac_above_0.5):")
print(summary.sort_values("median", ascending=False).to_string())

# koniec skryptu
