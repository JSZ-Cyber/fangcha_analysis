"""
Merge CSV results into Excel workbooks (one per dataset, 7 sheets each).

Output:
  anova_results_subtype.xlsx   (分型)
  anova_results_stage.xlsx     (分期)
"""

import os
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Sheet name mapping: CSV filename -> friendly sheet name
SHEET_MAP = [
    ("1_descriptive.csv",            "1-Descriptive"),
    ("2_normality.csv",              "2-Normality"),
    ("3_levene.csv",                 "3-Levene"),
    ("4_anova.csv",                  "4-ANOVA"),
    ("5_kruskal_wallis.csv",         "5-Kruskal-Wallis"),
    ("6a_pairwise_ttest.csv",        "6a-Pairwise t-test"),
    ("6b_pairwise_mannwhitney.csv",  "6b-Pairwise MannWhitney"),
    ("7_summary.csv",                "7-Summary"),
]


def merge_csvs_to_excel(csv_dir, excel_path):
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        for csv_name, sheet_name in SHEET_MAP:
            csv_path = os.path.join(csv_dir, csv_name)
            if os.path.exists(csv_path):
                df = pd.read_csv(csv_path)
                df.to_excel(writer, sheet_name=sheet_name, index=False)
                print(f"  [{sheet_name}] <- {csv_name}  ({len(df)} rows)")
            else:
                print(f"  [SKIP] {csv_name} not found")
    print(f"  => {excel_path}\n")


if __name__ == "__main__":
    # Subtype
    print("[Subtype]")
    merge_csvs_to_excel(
        os.path.join(SCRIPT_DIR, "output_subtype"),
        os.path.join(SCRIPT_DIR, "anova_results_subtype.xlsx"),
    )

    # Stage
    print("[Stage]")
    merge_csvs_to_excel(
        os.path.join(SCRIPT_DIR, "output_stage"),
        os.path.join(SCRIPT_DIR, "anova_results_stage.xlsx"),
    )

    print("All done.")
