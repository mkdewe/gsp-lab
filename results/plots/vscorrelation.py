import pandas as pd
import matplotlib.pyplot as plt

csv_path = r"D:\Studia\Projekty\gsp-lab\results\correlations\summary_report_rounded.csv"
df = pd.read_csv(csv_path)

# Rename columns to simple names
df = df.rename(columns={
    "avg_pearson_coef":"pearson",
    "avg_spearman_coef":"spearman",
    "avg_kendall_coef":"kendall"
})

fig, ax = plt.subplots(figsize=(8,6))

# Color by orientation
colors = df["orientation"].map({True:"red", False:"blue"})

# Plot
scatter = ax.scatter(df["pearson"], df["spearman"], c=colors, s=80)

for i, txt in enumerate(df["metric"]):
    ax.annotate(txt, (df["pearson"][i], df["spearman"][i]), fontsize=9, alpha=0.8)

ax.set_xlabel("Avg Pearson")
ax.set_ylabel("Avg Spearman")
ax.set_title("Pearson vs Spearman Correlation for Metrics")
ax.grid(True)

# custom legend
import matplotlib.patches as mpatches
red_patch = mpatches.Patch(color='red', label='orientation=True')
blue_patch = mpatches.Patch(color='blue', label='orientation=False')
ax.legend(handles=[red_patch, blue_patch])

plt.tight_layout()
plt.savefig("scatter_pearson_spearman.png", dpi=200)
plt.show()
