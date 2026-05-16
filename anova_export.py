"""
ANOVA Analysis & CSV Export
============================
Reads two breast cancer datasets (分型.csv, 分期.csv), performs full
one-way ANOVA workflow, and exports every section of the results into
separate CSV files organised under output folders.

Output structure
----------------
  output_subtype/          (分型)
    1_descriptive.csv
    2_normality.csv
    3_levene.csv
    4_anova.csv
    5_kruskal_wallis.csv
    6a_pairwise_ttest.csv
    6b_pairwise_mannwhitney.csv
    7_summary.csv

  output_stage/            (分期)
    (same set of files)

Usage
-----
    python anova_export.py
"""

import os
import pandas as pd
import numpy as np
from scipy import stats
from itertools import combinations
import warnings

warnings.filterwarnings("ignore")

# ============================================================
# Paths  (put source CSVs next to this script or adjust paths)
# ============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_SUBTYPE = os.path.join(SCRIPT_DIR, "data", "subtype.csv")       # 分型
CSV_STAGE   = os.path.join(SCRIPT_DIR, "data", "stage.csv")         # 分期

# Fallback: also look in Downloads (for local convenience)
if not os.path.exists(CSV_SUBTYPE):
    CSV_SUBTYPE = r"C:\Users\jisizhe\Downloads\分型.csv"
if not os.path.exists(CSV_STAGE):
    CSV_STAGE = r"C:\Users\jisizhe\Downloads\分期.csv"


# ============================================================
# Helper: load & reshape
# ============================================================
def load_and_reshape(csv_path, groups):
    """Pool 3 columns per group into one flat array per group."""
    df = pd.read_csv(csv_path)
    group_data = {}
    for g_idx, group in enumerate(groups):
        cols = [g_idx * 3 + i for i in range(3)]
        vals = df.iloc[:, cols].values.flatten()
        vals = vals[~np.isnan(vals)]
        group_data[group] = vals
    return group_data


# ============================================================
# Core analysis → returns dict of DataFrames
# ============================================================
def run_anova_to_dfs(group_data, dataset_name):
    groups = list(group_data.keys())
    arrays = [group_data[g] for g in groups]
    dfs = {}

    # ---- 1. Descriptive ----
    rows = []
    for g in groups:
        d = group_data[g]
        rows.append({
            "Group": g,
            "N": len(d),
            "Mean": round(np.mean(d), 4),
            "SD": round(np.std(d, ddof=1), 4),
            "Median": round(np.median(d), 4),
            "Min": round(np.min(d), 4),
            "Max": round(np.max(d), 4),
            "Q1": round(np.percentile(d, 25), 4),
            "Q3": round(np.percentile(d, 75), 4),
        })
    dfs["1_descriptive"] = pd.DataFrame(rows)

    # ---- 2. Normality (Shapiro-Wilk) ----
    rows = []
    all_normal = True
    for g in groups:
        d = group_data[g]
        sample = d if len(d) <= 5000 else np.random.choice(d, 5000, replace=False)
        w, p = stats.shapiro(sample)
        normal = "Yes" if p > 0.05 else "No"
        if p <= 0.05:
            all_normal = False
        rows.append({
            "Group": g,
            "N": len(d),
            "W_statistic": round(w, 6),
            "p_value": f"{p:.6e}",
            "Normal_alpha_0.05": normal,
        })
    dfs["2_normality"] = pd.DataFrame(rows)

    # ---- 3. Levene's test ----
    lev_stat, lev_p = stats.levene(*arrays)
    homo = "Yes" if lev_p > 0.05 else "No"
    dfs["3_levene"] = pd.DataFrame([{
        "Test": "Levene",
        "Statistic": round(lev_stat, 4),
        "p_value": f"{lev_p:.6e}",
        "Equal_variance_alpha_0.05": homo,
    }])

    # ---- 4. One-Way ANOVA ----
    f_stat, anova_p = stats.f_oneway(*arrays)
    grand_mean = np.mean(np.concatenate(arrays))
    ss_between = sum(len(group_data[g]) * (np.mean(group_data[g]) - grand_mean) ** 2 for g in groups)
    ss_within  = sum(np.sum((group_data[g] - np.mean(group_data[g])) ** 2) for g in groups)
    ss_total   = ss_between + ss_within
    eta_sq = ss_between / ss_total if ss_total > 0 else 0
    df_between = len(groups) - 1
    df_within  = sum(len(group_data[g]) for g in groups) - len(groups)
    ms_between = ss_between / df_between
    ms_within  = ss_within / df_within
    sig = "***" if anova_p < 0.001 else ("**" if anova_p < 0.01 else ("*" if anova_p < 0.05 else "ns"))

    dfs["4_anova"] = pd.DataFrame([{
        "Source": "Between Groups",
        "SS": round(ss_between, 4),
        "df": df_between,
        "MS": round(ms_between, 4),
        "F": round(f_stat, 4),
        "p_value": f"{anova_p:.6e}",
        "Significance": sig,
        "Eta_squared": round(eta_sq, 4),
    }, {
        "Source": "Within Groups",
        "SS": round(ss_within, 4),
        "df": df_within,
        "MS": round(ms_within, 4),
        "F": "",
        "p_value": "",
        "Significance": "",
        "Eta_squared": "",
    }, {
        "Source": "Total",
        "SS": round(ss_total, 4),
        "df": df_between + df_within,
        "MS": "",
        "F": "",
        "p_value": "",
        "Significance": "",
        "Eta_squared": "",
    }])

    # ---- 5. Kruskal-Wallis ----
    kw_stat, kw_p = stats.kruskal(*arrays)
    kw_sig = "***" if kw_p < 0.001 else ("**" if kw_p < 0.01 else ("*" if kw_p < 0.05 else "ns"))
    dfs["5_kruskal_wallis"] = pd.DataFrame([{
        "Test": "Kruskal-Wallis",
        "H_statistic": round(kw_stat, 4),
        "p_value": f"{kw_p:.6e}",
        "Significance": kw_sig,
    }])

    # ---- 6a. Pairwise t-tests (Welch) + Bonferroni ----
    pairs = list(combinations(groups, 2))
    n_comp = len(pairs)
    rows = []
    for g1, g2 in pairs:
        t, p_raw = stats.ttest_ind(group_data[g1], group_data[g2], equal_var=False)
        p_adj = min(p_raw * n_comp, 1.0)
        s = "***" if p_adj < 0.001 else ("**" if p_adj < 0.01 else ("*" if p_adj < 0.05 else "ns"))
        rows.append({
            "Group_1": g1,
            "Group_2": g2,
            "t_statistic": round(t, 4),
            "p_raw": f"{p_raw:.6e}",
            "p_adjusted_Bonferroni": f"{p_adj:.6e}",
            "Significance": s,
        })
    dfs["6a_pairwise_ttest"] = pd.DataFrame(rows)

    # ---- 6b. Pairwise Mann-Whitney U + Bonferroni ----
    rows = []
    for g1, g2 in pairs:
        u, p_raw = stats.mannwhitneyu(group_data[g1], group_data[g2], alternative="two-sided")
        p_adj = min(p_raw * n_comp, 1.0)
        s = "***" if p_adj < 0.001 else ("**" if p_adj < 0.01 else ("*" if p_adj < 0.05 else "ns"))
        rows.append({
            "Group_1": g1,
            "Group_2": g2,
            "U_statistic": round(u, 1),
            "p_raw": f"{p_raw:.6e}",
            "p_adjusted_Bonferroni": f"{p_adj:.6e}",
            "Significance": s,
        })
    dfs["6b_pairwise_mannwhitney"] = pd.DataFrame(rows)

    # ---- 7. Summary ----
    if all_normal and homo == "Yes":
        rec = "One-Way ANOVA + Bonferroni post-hoc"
        key_result = f"F={f_stat:.4f}, p={anova_p:.6e} {sig}"
    elif all_normal:
        rec = "Welch ANOVA (unequal variance) + Games-Howell"
        key_result = f"F={f_stat:.4f}, p={anova_p:.6e} {sig}"
    else:
        rec = "Kruskal-Wallis + Mann-Whitney U (non-parametric)"
        key_result = f"H={kw_stat:.4f}, p={kw_p:.6e} {kw_sig}"

    dfs["7_summary"] = pd.DataFrame([
        {"Item": "Dataset", "Value": dataset_name},
        {"Item": "Number of groups", "Value": str(len(groups))},
        {"Item": "Total observations", "Value": str(sum(len(group_data[g]) for g in groups))},
        {"Item": "All groups normal (Shapiro-Wilk, a=0.05)?", "Value": "Yes" if all_normal else "No"},
        {"Item": "Equal variance (Levene, a=0.05)?", "Value": homo},
        {"Item": "Recommended method", "Value": rec},
        {"Item": "Key result", "Value": key_result},
        {"Item": "Effect size (Eta-squared)", "Value": str(round(eta_sq, 4))},
    ])

    return dfs


# ============================================================
# Export: one folder per dataset, one CSV per section
# ============================================================
def export_dfs(dfs, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for name, df in dfs.items():
        path = os.path.join(out_dir, f"{name}.csv")
        df.to_csv(path, index=False, encoding="utf-8-sig")
        print(f"  -> {path}")


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    # --- 1. Subtype (分型) ---
    groups_subtype = ["Normal", "Basal", "Her2", "LumA", "LumB"]
    data_subtype = load_and_reshape(CSV_SUBTYPE, groups_subtype)
    dfs_subtype = run_anova_to_dfs(data_subtype, "Breast Cancer Subtypes")

    out1 = os.path.join(SCRIPT_DIR, "output_subtype")
    print(f"\n[Subtype] Exporting to {out1}")
    export_dfs(dfs_subtype, out1)

    # --- 2. Stage (分期) ---
    groups_stage = ["Control", "Stage I", "Stage II", "Stage III", "Stage IV"]
    data_stage = load_and_reshape(CSV_STAGE, groups_stage)
    dfs_stage = run_anova_to_dfs(data_stage, "Cancer Staging")

    out2 = os.path.join(SCRIPT_DIR, "output_stage")
    print(f"\n[Stage] Exporting to {out2}")
    export_dfs(dfs_stage, out2)

    print("\nAll done.")
