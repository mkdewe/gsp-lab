# -*- coding: utf-8 -*-
"""
single_best_rank_final_palette_v2.py

Updated: loads lDDT values from aggregated global_report.csv (results/lddt/finished)
instead of gGSP data. Saves plots as Figure4a/4b/4c in PNG and PDF.

Modified: excluded targets (PZ10tRNA, PZ26tRNA, PZ27tRNA, PZ28tRNA) are skipped.
"""
import os
import re
import sys
import argparse
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import kendalltau

# === CONFIG (adjust paths if needed) ===
# Now wczytujemy globalny raport lDDT
LDDT_GLOBAL_CSV = r"D:\Studia\Projekty\gsp-lab\results\lddt\finished\global_report.csv"
OTHER_METRICS_DIR = r"D:\Studia\Projekty\gsp-lab\results\other-metrics\rnapuzzles.github.io"
CSV_MAPPING_FILE = r"D:\Studia\Projekty\gsp-lab\results\plots\RNA-types.csv"
OUTPUT_DIR = r"D:\Studia\Projekty\gsp-lab\results\plots"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# === Excluded targets (will be skipped during loading) ===
EXCLUDED_TARGETS = ['PZ10tRNA', 'PZ26tRNA', 'PZ27tRNA', 'PZ28tRNA']

# === Consistent palette for RNA categories ===
RNA_COLORS = {
    'Aptamer': '#66c2a5',               # green
    'Riboswitch': '#fc8d62',            # orange
    'Ribozyme': '#8da0cb',              # blue
    'Viral element': '#e78ac3',         # pink
    'Regulatory element': '#ffd92f',    # yellow
    'Synthetic construct': '#e5c494',   # beige
    'RNA-binding protein': '#b3b3b3',   # gray
    'Nano construct': '#8dd3c7',        # light turquoise
    'Other RNA Structures': '#a6d854',  # light green
    'Red category': '#e57373'           # soft red (fallback)
}

# === Defaults for selection / plots ===
MIN_MODELS_FOR_CANDIDATE = 6
TOP_K_BUMP = 15
TOP_N_CANDIDATES = 10

# === helper functions ===
def read_single_value_from_file(path):
    """Try to read a numeric value from a file: regex first, pandas fallback."""
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            txt = f.read()
    except Exception as e:
        warnings.warn(f"Cannot read file {path}: {e}")
        return None
    m = re.search(r'(-?\d+\.\d+|-?\d+)', txt)
    if m:
        try:
            return float(m.group(0))
        except:
            pass
    try:
        df = pd.read_csv(path, header=None)
        for col in df.columns:
            for v in df[col].dropna().astype(str).values:
                try:
                    return float(v)
                except:
                    continue
    except Exception:
        pass
    return None


def load_and_merge_all(lddt_csv_path, other_metrics_dir, mapping_csv):
    """
    Load lDDT values from global_report.csv, then merge with reference metrics.
    """
    mapping_df = pd.read_csv(mapping_csv)
    target_mapping = pd.Series(mapping_df.type.values, index=mapping_df.target_label).to_dict()

    # Wczytaj dane lDDT z pliku global_report.csv
    if not os.path.exists(lddt_csv_path):
        raise FileNotFoundError(f"lDDT global report not found: {lddt_csv_path}")
    lddt_df = pd.read_csv(lddt_csv_path)   # kolumny: Puzzle, ModelGroup, BestModel, BestSolution, lDDT

    # Usuń wykluczone targety
    lddt_df = lddt_df[~lddt_df['Puzzle'].isin(EXCLUDED_TARGETS)]

    # Przekształć na strukturę podobną do gGSP: Target, ModelGroup, lDDT (jako gGSP)
    lddt_df = lddt_df.rename(columns={'Puzzle': 'Target', 'lDDT': 'gGSP'})
    lddt_df['Source'] = lddt_df['Target'].apply(lambda x: 'CASP RNA' if x.startswith('CR') else 'RNA-Puzzles')
    lddt_df['Type'] = lddt_df['Target'].map(target_mapping).fillna('Other RNA Structures')

    # Dla każdego targetu przygotuj dane z podziałem na Lab i Num (z ModelGroup)
    target_lddt_dict = {}
    for target, group in lddt_df.groupby('Target'):
        # Wyciągnij Lab i Num z ModelGroup (zakładamy format "TS029_1")
        def split_modelgroup(mg):
            parts = mg.split('_')
            if len(parts) >= 2:
                # Pierwsza część to lab, druga to num
                return parts[0], parts[1]
            else:
                # fallback
                return mg, "0"
        group['Lab'] = group['ModelGroup'].apply(lambda x: split_modelgroup(x)[0])
        group['Num'] = group['ModelGroup'].apply(lambda x: split_modelgroup(x)[1])
        # Dodaj tymczasową kolumnę "file" – potrzebna później w rankingach
        group['file'] = group['ModelGroup']   # placeholder
        target_lddt_dict[target] = group[['Lab', 'Num', 'gGSP', 'file']].copy()

    merged = []
    print("Loading data...")
    for target, lddt_target_df in target_lddt_dict.items():
        other_dir = os.path.join(other_metrics_dir, target)
        if not os.path.exists(other_dir):
            warnings.warn(f"No reference-metrics folder for '{target}' -> skipping.")
            continue

        # find metric file: RMSD.csv, TM-score.csv, or any CSV with numeric column
        rmsd_fp = os.path.join(other_dir, "RMSD.csv")
        tm_fp = os.path.join(other_dir, "TM-score.csv")
        cand = glob.glob(os.path.join(other_dir, "*.csv"))
        df_m = None
        chosen_col = None
        if os.path.exists(rmsd_fp):
            df_m = pd.read_csv(rmsd_fp); chosen_col = "RMSD"
        elif os.path.exists(tm_fp):
            df_m = pd.read_csv(tm_fp); chosen_col = "TMscore"
        else:
            for c in sorted(cand):
                try:
                    tmp = pd.read_csv(c)
                except:
                    continue
                numeric_cols = tmp.select_dtypes(include=[np.number]).columns.tolist()
                numeric_cols = [cc for cc in numeric_cols if cc.lower() not in ['index','unnamed: 0']]
                if numeric_cols:
                    df_m = tmp; chosen_col = numeric_cols[0]; break
        if df_m is None:
            warnings.warn(f"No metric file found for '{target}' -> skipping.")
            continue

        # unify Lab/Num column names
        if 'Lab' not in df_m.columns or 'Num' not in df_m.columns:
            lab_col = num_col = None
            for c in df_m.columns:
                if c.lower() == 'lab': lab_col = c
                if c.lower() == 'num': num_col = c
            if lab_col and num_col:
                df_m.rename(columns={lab_col:'Lab', num_col:'Num'}, inplace=True)
            else:
                warnings.warn(f"Metric file for '{target}' has no 'Lab' and 'Num' -> skipping.")
                continue

        df_m['Lab'] = df_m['Lab'].astype(str).str.lower().str.strip()
        df_m['Num'] = df_m['Num'].astype(str).str.strip()

        # choose result column
        if chosen_col in df_m.columns:
            res_col = chosen_col
        else:
            numeric_cols = df_m.select_dtypes(include=[np.number]).columns.tolist()
            numeric_cols = [c for c in numeric_cols if c.lower() not in ['index','unnamed: 0']]
            if not numeric_cols:
                warnings.warn(f"No numeric columns in metric file for '{target}' -> skipping.")
                continue
            res_col = numeric_cols[0]

        df_m_sub = df_m[['Lab','Num',res_col]].rename(columns={res_col:'Metric'})
        # Wyrównaj format Lab i Num w danych lDDT
        lddt_target_df['Lab'] = lddt_target_df['Lab'].astype(str).str.lower().str.strip()
        lddt_target_df['Num'] = lddt_target_df['Num'].astype(str).str.strip()
        df_join = pd.merge(lddt_target_df, df_m_sub, on=['Lab','Num'], how='inner')
        if df_join.empty:
            warnings.warn(f"Merge lDDT<->Metric for '{target}' produced empty result -> skipping.")
            continue

        # Dodaj metadane
        df_join['Target'] = target
        df_join['Type'] = target_mapping.get(target, 'Other RNA Structures')
        df_join['MetricName'] = res_col
        df_join['MetricIsRMSD'] = ('rmsd' in res_col.lower())
        df_join['MetricIsTM'] = ('tm' in res_col.lower())
        merged.append(df_join)

    if not merged:
        raise RuntimeError("No merged data. Check paths and files.")
    full_df = pd.concat(merged, ignore_index=True)
    full_df['label'] = full_df['Lab'].astype(str) + "_" + full_df['Num'].astype(str)
    return full_df


def compute_ranks_and_tau(full_df):
    # rank lDDT (higher is better)
    def rank_lDDT_for_group(g):
        if g['gGSP'].nunique(dropna=True) <= 1:
            warnings.warn(f"Target '{g['Target'].iloc[0]}': identical lDDT values -> ranking by filename.")
            ordered = g.sort_values('file')
            ranks = pd.Series(np.arange(1, len(ordered)+1), index=ordered.index)
            return ranks.reindex(g.index)
        return g['gGSP'].rank(ascending=False, method='min')

    full_df['rank_gGSP'] = full_df.groupby('Target', group_keys=False).apply(lambda grp: rank_lDDT_for_group(grp))

    # rank Metric (RMSD smaller is better; TMscore bigger is better)
    def rank_metric_for_group(g):
        if g['MetricIsRMSD'].iloc[0]:
            return g['Metric'].rank(ascending=True, method='min')
        else:
            return g['Metric'].rank(ascending=False, method='min')

    full_df['rank_Metric'] = full_df.groupby('Target', group_keys=False).apply(lambda grp: rank_metric_for_group(grp))
    full_df['delta_rank_signed'] = full_df['rank_gGSP'] - full_df['rank_Metric']
    full_df['delta_rank_abs'] = full_df['delta_rank_signed'].abs()

    # compute Kendall tau per target and some summary stats
    rows = []
    for (t,tp), grp in full_df.groupby(['Target','Type']):
        valid = grp.dropna(subset=['rank_gGSP','rank_Metric'])
        n = len(valid)
        gstd = float(np.nanstd(valid['gGSP'])) if n>0 else np.nan
        mstd = float(np.nanstd(valid['Metric'])) if n>0 else np.nan
        gmean = float(np.nanmean(valid['gGSP'])) if n>0 else np.nan
        mmean = float(np.nanmean(valid['Metric'])) if n>0 else np.nan
        if n >= 4:
            tau, p = kendalltau(valid['rank_gGSP'], valid['rank_Metric'])
            tau = float(np.round(tau,4)) if not np.isnan(tau) else np.nan
        else:
            tau = np.nan
        rows.append({'Target':t, 'Type':tp, 'Tau':tau, 'Nmodels':n,
                     'gGSP_std':gstd, 'Metric_std':mstd, 'gGSP_mean':gmean, 'Metric_mean':mmean})
    df_tau = pd.DataFrame(rows)
    return full_df, df_tau


def choose_target_interactive(df_tau, full_df, args):
    # build candidates by score = Tau * std(gGSP) * std(Metric)
    candidates = []
    for _, r in df_tau.iterrows():
        n = int(r['Nmodels'])
        if n >= args.min_models and not np.isnan(r['Tau']):
            gstd = r['gGSP_std'] if not np.isnan(r['gGSP_std']) else 0.0
            mstd = r['Metric_std'] if not np.isnan(r['Metric_std']) else 0.0
            score = r['Tau'] * (gstd if gstd is not None else 0.0) * (mstd if mstd is not None else 0.0)
            candidates.append({'Target': r['Target'], 'Type': r['Type'], 'Nmodels': n, 'Tau': r['Tau'],
                               'gGSP_std': gstd, 'Metric_std': mstd, 'score': score})
    cand_df = pd.DataFrame(candidates).sort_values('score', ascending=False)

    if args.list_targets:
        print("\nAvailable targets (brief):")
        print(df_tau[['Target','Type','Nmodels','Tau']].sort_values(['Type','Target']).to_string(index=False))
        sys.exit(0)

    if args.target:
        if args.target not in full_df['Target'].unique():
            raise ValueError(f"Provided target {args.target} not found in data.")
        print(f"Selected target from argument: {args.target}")
        return args.target

    if not cand_df.empty:
        print("\nTop candidate targets (by score = Tau * std(lDDT) * std(Metric)):")
        print(cand_df[['Target','Type','Nmodels','Tau','score']].head(args.top_candidates).to_string(index=True))
    else:
        print("\nNo candidates meeting the minimum model count.")
        tby = df_tau.sort_values('Tau', ascending=False)
        print("Fallback - top targets by Tau:")
        print(tby[['Target','Type','Nmodels','Tau']].head(args.top_candidates).to_string(index=True))

    # interactive selection if possible
    try:
        if sys.stdin is None or not sys.stdin.isatty():
            # non-interactive: choose best candidate or fallback
            if not cand_df.empty:
                choice = cand_df.iloc[0]['Target']
                print(f"\nNo interactive terminal — automatically selecting: {choice}")
                return choice
            else:
                choice = df_tau.sort_values('Tau', ascending=False).iloc[0]['Target']
                print(f"\nNo interactive terminal — fallback selecting by Tau: {choice}")
                return choice

        prompt = ("\nType an index from the list above, or exact target name,\n"
                  "or press Enter to select the default (best): ")
        sel = input(prompt).strip()
        if sel == "":
            if not cand_df.empty:
                return cand_df.iloc[0]['Target']
            else:
                return df_tau.sort_values('Tau', ascending=False).iloc[0]['Target']
        if sel.isdigit():
            idx = int(sel)
            if not cand_df.empty and 0 <= idx < len(cand_df):
                return cand_df.reset_index(drop=True).iloc[idx]['Target']
            else:
                tby = df_tau.sort_values('Tau', ascending=False).reset_index(drop=True)
                if 0 <= idx < len(tby):
                    return tby.iloc[idx]['Target']
                raise ValueError("Invalid index.")
        if sel in full_df['Target'].unique():
            return sel
        else:
            raise ValueError("Typed target name not found.")
    except Exception as e:
        print(f"Interactive selection error: {e}. Selecting best candidate/fallback automatically.")
        if not cand_df.empty:
            return cand_df.iloc[0]['Target']
        return df_tau.sort_values('Tau', ascending=False).iloc[0]['Target']


def plot_and_save(target, full_df, df_tau, top_k_bump):
    # data for target
    target_df = full_df[full_df['Target'] == target].copy()
    if target_df.empty:
        raise RuntimeError(f"No data for selected target '{target}'")
    target_df = target_df.sort_values('rank_Metric').head(top_k_bump).reset_index(drop=True)

    print("\nDiagnostics (top K for selected target):")
    print(target_df[['label','file','gGSP','rank_gGSP','Metric','rank_Metric','delta_rank_signed']].to_string(index=False))

    sns.set_theme(style="whitegrid")

    # A) bump chart: left = model name + Metric, right = model name + lDDT
    fig_a, ax_a = plt.subplots(figsize=(8, max(6, 0.35*len(target_df))))
    x_left, x_right = 0, 1
    left_metric_offset = 0.22  # vertical offset for metric text under left label
    right_gsp_offset = 0.22    # vertical offset for lDDT text under right label
    right_label_x = x_right + 0.02

    for _, row in target_df.iterrows():
        y_metric = row['rank_Metric']
        y_gsp = row['rank_gGSP']
        diff = row['delta_rank_signed']
        # line color encodes whether lDDT ranks the model better (<0), worse (>0), or same
        line_color = '#a6d854' if diff < 0 else ('#e57373' if diff > 0 else '#b3b3b3')
        # marker face color according to RNA category (if available)
        marker_face = RNA_COLORS.get(row['Type'], '#bdbdbd')
        ax_a.plot([x_left, x_right], [y_metric, y_gsp],
                  marker='o', linewidth=2, color=line_color,
                  markerfacecolor=marker_face, markeredgecolor='k', markeredgewidth=0.5, alpha=0.95)

        # left: model name
        ax_a.text(x_left - 0.03, y_metric, row['label'], ha='right', va='center', fontsize=9)
        # left: metric value under the model name (gray)
        try:
            metric_val = float(row['Metric']); metric_text = f"{row['MetricName']} = {metric_val:.2f}"
        except:
            metric_text = f"{row['MetricName']} = {row['Metric']}"
        ax_a.text(x_left - 0.03, y_metric + left_metric_offset, metric_text, ha='right', va='top', fontsize=8, color='gray')

        # right: model name and lDDT value underneath (gray)
        ax_a.text(right_label_x, y_gsp, row['label'], ha='left', va='center', fontsize=9)
        try:
            lddt_val = float(row['gGSP']); gsp_text = f"lDDT = {lddt_val:.3f}"
        except:
            gsp_text = f"lDDT = {row['gGSP']}"
        ax_a.text(right_label_x, y_gsp + right_gsp_offset, gsp_text, ha='left', va='top', fontsize=8, color='gray')

    ax_a.set_xlim(-0.3, 1.35)
    ax_a.set_xticks([x_left, x_right])
    ax_a.set_xticklabels([f"{target_df['MetricName'].iloc[0]} ranking (reference)", "lDDT ranking"])
    ax_a.invert_yaxis()
    ax_a.set_ylabel("Rank position (1 = best model)", fontsize=10)
    ax_a.set_title(f"Ranking consistency — {target}", fontsize=14, weight='bold')
    out_a_png = os.path.join(OUTPUT_DIR, f"Figure4a_Bump_{target}.png")
    out_a_pdf = os.path.join(OUTPUT_DIR, f"Figure4a_Bump_{target}.pdf")
    fig_a.tight_layout()
    fig_a.savefig(out_a_png, dpi=300, bbox_inches='tight')
    fig_a.savefig(out_a_pdf, dpi=300, bbox_inches='tight')
    print(f"Saved bump chart: {out_a_png} and {out_a_pdf}")

    # B) histogram of absolute delta ranks
    fig_b, ax_b = plt.subplots(figsize=(8,5))
    max_bin = int(full_df['delta_rank_abs'].max()) if not full_df['delta_rank_abs'].isnull().all() else 1
    bins = np.arange(0, max_bin + 2) - 0.5
    hist_color = RNA_COLORS.get('Ribozyme', '#fc8d62')
    ax_b.hist(full_df['delta_rank_abs'].dropna(), bins=bins, color=hist_color, edgecolor='black', alpha=0.9)
    ax_b.set_xlabel("Absolute rank difference |rank_lDDT − rank_reference|", fontsize=10)
    ax_b.set_ylabel("Number of models", fontsize=10)
    ax_b.set_title("Distribution of absolute rank differences across all targets", fontsize=13)
    ax_b.set_xlim(left=-0.5)
    out_b_png = os.path.join(OUTPUT_DIR, "Figure4b_Hist_DeltaRank.png")
    out_b_pdf = os.path.join(OUTPUT_DIR, "Figure4b_Hist_DeltaRank.pdf")
    fig_b.tight_layout(); fig_b.savefig(out_b_png, dpi=300, bbox_inches='tight'); fig_b.savefig(out_b_pdf, dpi=300, bbox_inches='tight')
    print(f"Saved histogram of |Δrank|: {out_b_png} and {out_b_pdf}")

    # C) Kendall's tau per target (violin + mean±std). Single-target categories -> big colored dot
    fig_c, ax_c = plt.subplots(figsize=(14, 5))

    order = [
        'Aptamer', 'Riboswitch', 'Ribozyme', 'Viral element',
        'Regulatory element', 'Synthetic construct',
        'RNA-binding protein', 'Nano construct', 'Other RNA Structures'
    ]

    existing_order = [t for t in order if t in df_tau['Type'].unique()] + \
                    [t for t in sorted(df_tau['Type'].unique()) if t not in order]

    df_tau_plot = df_tau.copy()
    df_tau_nonan = df_tau_plot.dropna(subset=['Tau'])

    # violin only when >=2 targets
    violin_types = []
    for t in existing_order:
        n = len(df_tau_nonan[df_tau_nonan['Type'] == t])
        if n >= 2:
            violin_types.append(t)

    if violin_types:
        sns.violinplot(
            data=df_tau_nonan[df_tau_nonan['Type'].isin(violin_types)],
            x='Type',
            y='Tau',
            order=existing_order,
            palette={t: RNA_COLORS.get(t, '#bdbdbd') for t in existing_order},
            ax=ax_c,
            inner=None,
            cut=0
        )

    # mean ± std
    for i, t in enumerate(existing_order):
        vals = df_tau_nonan[df_tau_nonan['Type'] == t]['Tau'].values

        if len(vals) >= 2:
            mean = np.mean(vals)
            std = np.std(vals)
            ax_c.errorbar(
                i,
                mean,
                yerr=std,
                fmt='o',
                color='black',
                ecolor='black',
                elinewidth=2,
                capsize=5,
                markersize=7,
                zorder=5
            )
        elif len(vals) == 1:
            ax_c.scatter(
                i,
                vals[0],
                s=220,
                color=RNA_COLORS.get(t, '#777777'),
                edgecolors='black',
                linewidths=0.8,
                zorder=6
            )

    # Labels with target counts
    counts = df_tau_nonan['Type'].value_counts()
    labels = [
        f"{t}\n(n={counts[t]})" if t in counts else f"{t}\n(n=0)"
        for t in existing_order
    ]
    ax_c.set_xticks(range(len(existing_order)))
    ax_c.set_xticklabels(labels, rotation=0, ha='center', fontsize=10)
    ax_c.tick_params(axis='x', pad=10)

    ax_c.set_ylim(-0.2, 1.0)
    ax_c.axhline(0, color='black', linewidth=1, alpha=0.3)
    ax_c.axhline(0.5, color='red', linestyle='--', alpha=0.6, label="High agreement (τ = 0.5)")
    ax_c.set_title("Per-target ranking agreement (Kendall’s τ) by RNA type", fontsize=14, pad=20)
    ax_c.set_ylabel("Kendall’s τ (agreement between rankings)", fontsize=12)
    ax_c.set_xlabel("")
    ax_c.legend(loc='upper right')

    out_c_png = os.path.join(OUTPUT_DIR, "Figure4c_Kendall_Tau.png")
    out_c_pdf = os.path.join(OUTPUT_DIR, "Figure4c_Kendall_Tau.pdf")
    fig_c.tight_layout()
    fig_c.savefig(out_c_png, dpi=300, bbox_inches='tight')
    fig_c.savefig(out_c_pdf, dpi=300, bbox_inches='tight')
    print(f"Saved Kendall's Tau plot: {out_c_png} and {out_c_pdf}")

    # save diagnostics
    full_df.to_csv(os.path.join(OUTPUT_DIR, "full_joined_ranking_data.csv"), index=False)
    df_tau.to_csv(os.path.join(OUTPUT_DIR, "ranking_summary_per_target.csv"), index=False)
    print(f"Saved diagnostic CSVs to {OUTPUT_DIR}")


def main():
    parser = argparse.ArgumentParser(description="Compare lDDT ranking with reference metric (RMSD/TM).")
    parser.add_argument("--target", type=str, default=None, help="Target name (e.g. CR1138) to skip interactive selection.")
    parser.add_argument("--list-targets", action="store_true", help="List available targets and exit.")
    parser.add_argument("--min-models", type=int, default=MIN_MODELS_FOR_CANDIDATE, help="Minimum number of models to consider a target a candidate.")
    parser.add_argument("--top-k-bump", type=int, default=TOP_K_BUMP, help="How many models to show on bump chart.")
    parser.add_argument("--top-candidates", type=int, default=TOP_N_CANDIDATES, help="How many candidates to list.")
    args = parser.parse_args()

    # load and process – now from lDDT global report
    full_df = load_and_merge_all(LDDT_GLOBAL_CSV, OTHER_METRICS_DIR, CSV_MAPPING_FILE)
    full_df, df_tau = compute_ranks_and_tau(full_df)

    # choose
    chosen = choose_target_interactive(df_tau, full_df, args)
    plot_and_save(chosen, full_df, df_tau, args.top_k_bump)
    print("\nDone.")

if __name__ == "__main__":
    import glob   # brakowało wcześniej
    main()