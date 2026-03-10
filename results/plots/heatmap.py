import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

# === Load data ===
csv_path = r"D:\Studia\Projekty\gsp-lab\results\zestawy_modele.csv"
df = pd.read_csv(csv_path)

# Rename if necessary
df = df.rename(columns={
    "nazwa_zestawu": "target",
    "liczba_modeli": "num_models"
})

# === Sort by number of models (malejąco) ===
df_sorted = df.sort_values("num_models", ascending=False).reset_index(drop=True)

# === Aesthetic settings ===
plt.style.use('seaborn-whitegrid')
fig, ax = plt.subplots(figsize=(16, 8))

# === Assign colors based on prefix ===
# PZ sets -> original #0099CC (RNA PUZZLE), others -> lighter shade #66CCFF (CASP-RNA)
colors = ['#0099CC' if name.startswith('PZ') else '#66CCFF' for name in df_sorted['target']]

# Create vertical bars with conditional colors
x_pos = np.arange(len(df_sorted))
bars = ax.bar(x_pos, df_sorted["num_models"], 
              color=colors,
              edgecolor='black', 
              linewidth=0.5)

# Set x-axis labels to target names, rotated
ax.set_xticks(x_pos)
ax.set_xticklabels(df_sorted["target"], rotation=90, fontsize=6)

# Titles and labels
ax.set_title("Number of Predicted Models per Target", fontsize=16, fontweight='bold')
ax.set_xlabel("Target", fontsize=14)
ax.set_ylabel("Number of Predicted Models", fontsize=14)

# === Mean line and label (integer) ===
mean_val = df_sorted["num_models"].mean()
ax.axhline(y=mean_val, color='red', linestyle='--', linewidth=2, label=f'Średnia: {int(mean_val)}')
ax.text(len(df_sorted)-0.5, mean_val + 0.5, f'{int(mean_val)}', color='red',
        ha='right', va='bottom', fontsize=10, fontweight='bold')

# === Legend for colors ===
# Create custom patches for legend
pz_patch = mpatches.Patch(color='#0099CC', label='RNA-Puzzles')
casp_patch = mpatches.Patch(color='#66CCFF', label='CASP-RNA')
# Add mean line to legend (optional)
# mean_patch = mpatches.Patch(color='red', linestyle='--', label=f'Średnia: {int(mean_val)}')
# Add legend in upper right corner
ax.legend(handles=[pz_patch, casp_patch], loc='upper right', fontsize=10)

# Grid formatting
ax.grid(axis='y', alpha=0.6)
ax.set_axisbelow(True)

# === Adjust layout to prevent cutting off labels and legend ===
plt.tight_layout()
plt.subplots_adjust(bottom=0.4)  # Większy dolny margines dla obróconych etykiet

# Save and show
plt.savefig("models_per_target_vertical.png", dpi=200, bbox_inches='tight')
plt.show()