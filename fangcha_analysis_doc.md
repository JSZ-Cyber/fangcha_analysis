# ANOVA Analysis for Breast Cancer Data (方差分析)

## 1. Project Overview

This project performs **one-way ANOVA** (and non-parametric alternatives) on two breast cancer datasets to examine differences across groups:

| Dataset | File | Groups | Obs. per group |
|---------|------|--------|----------------|
| Subtype (分型) | `分型.csv` | Normal, Basal, Her2, LumA, LumB | 82–568 per column |
| Stage (分期)   | `分期.csv` | Control, Stage I, Stage II, Stage III, Stage IV | 19–619 per column |

Each dataset has 15 columns (3 repeated measures per group). All 3 columns are pooled for each group before analysis.

## 2. Analysis Workflow

The script `anova_export.py` executes the following steps for each dataset:

### Step 1 — Descriptive Statistics
Computes N, Mean, SD, Median, Min, Max, Q1, Q3 for each group.

### Step 2 — Normality Test (Shapiro-Wilk)
Tests whether each group's data follows a normal distribution (alpha = 0.05).

### Step 3 — Homogeneity of Variance (Levene's Test)
Tests whether all groups share equal variance (alpha = 0.05).

### Step 4 — One-Way ANOVA (F-test)
Computes between-group and within-group sum of squares, F-statistic, p-value, and effect size (Eta-squared).

### Step 5 — Kruskal-Wallis Test
Non-parametric alternative to ANOVA when normality / equal variance assumptions are violated.

### Step 6a — Pairwise t-tests (Welch's) with Bonferroni Correction
All 10 pairwise comparisons using Welch's t-test, p-values adjusted by Bonferroni method.

### Step 6b — Mann-Whitney U Tests with Bonferroni Correction
Non-parametric pairwise comparisons, also Bonferroni corrected.

### Step 7 — Summary
Consolidates key findings and recommends the appropriate statistical method.

## 3. Output Structure

```
fangcha_analysis/
├── anova_export.py                  # Main analysis script
├── requirements.txt                 # Python dependencies
├── fangcha_analysis_doc.md          # This documentation
│
├── output_subtype/                  # Results for 分型 (Subtype)
│   ├── 1_descriptive.csv
│   ├── 2_normality.csv
│   ├── 3_levene.csv
│   ├── 4_anova.csv
│   ├── 5_kruskal_wallis.csv
│   ├── 6a_pairwise_ttest.csv
│   ├── 6b_pairwise_mannwhitney.csv
│   └── 7_summary.csv
│
└── output_stage/                    # Results for 分期 (Stage)
    ├── 1_descriptive.csv
    ├── 2_normality.csv
    ├── 3_levene.csv
    ├── 4_anova.csv
    ├── 5_kruskal_wallis.csv
    ├── 6a_pairwise_ttest.csv
    ├── 6b_pairwise_mannwhitney.csv
    └── 7_summary.csv
```

## 4. How to Run

```bash
pip install -r requirements.txt
python anova_export.py
```

The script will read the source CSVs and write results into `output_subtype/` and `output_stage/`.

## 5. Key Results

### 5.1 Subtype (分型)

- **Normality**: All 5 groups fail the Shapiro-Wilk test (p < 0.05)
- **Variance homogeneity**: Levene p = 5.57e-21 (rejected)
- **Recommended method**: Kruskal-Wallis + Mann-Whitney U
- **Kruskal-Wallis**: H = 2920.66, **p < 0.001 \*\*\***, Eta-squared = 0.6336

Pairwise results (Mann-Whitney U, Bonferroni corrected):

| Comparison | p (adjusted) | Sig |
|------------|-------------:|:---:|
| Normal vs Basal | < 0.001 | \*\*\* |
| Normal vs Her2 | < 0.001 | \*\*\* |
| Normal vs LumA | < 0.001 | \*\*\* |
| Normal vs LumB | < 0.001 | \*\*\* |
| Basal vs Her2 | < 0.001 | \*\*\* |
| Basal vs LumA | < 0.001 | \*\*\* |
| Basal vs LumB | < 0.001 | \*\*\* |
| Her2 vs LumA | 1.000 | ns |
| Her2 vs LumB | 0.004 | \*\* |
| LumA vs LumB | < 0.001 | \*\*\* |

**Conclusion**: All groups differ significantly except **Her2 vs LumA** (ns).

### 5.2 Stage (分期)

- **Normality**: 4 of 5 groups fail Shapiro-Wilk (Control, Stage I-III)
- **Variance homogeneity**: Levene p = 0.017 (rejected)
- **Recommended method**: Kruskal-Wallis + Mann-Whitney U
- **Kruskal-Wallis**: H = 2544.71, **p < 0.001 \*\*\***, Eta-squared = 0.5750

Pairwise results (Mann-Whitney U, Bonferroni corrected):

| Comparison | p (adjusted) | Sig |
|------------|-------------:|:---:|
| Control vs Stage I | < 0.001 | \*\*\* |
| Control vs Stage II | < 0.001 | \*\*\* |
| Control vs Stage III | < 0.001 | \*\*\* |
| Control vs Stage IV | < 0.001 | \*\*\* |
| Stage I vs Stage II | 0.662 | ns |
| Stage I vs Stage III | 1.000 | ns |
| Stage I vs Stage IV | 1.000 | ns |
| Stage II vs Stage III | 0.062 | ns |
| Stage II vs Stage IV | 1.000 | ns |
| Stage III vs Stage IV | 1.000 | ns |

**Conclusion**: **Control differs significantly from all Stage groups**, but **Stage I through IV show no significant differences** among themselves.

## 6. Statistical Methods Reference

| Method | Purpose | Assumption |
|--------|---------|------------|
| Shapiro-Wilk | Normality test | — |
| Levene's test | Variance homogeneity | — |
| One-Way ANOVA | Compare means across groups | Normality + equal variance |
| Kruskal-Wallis | Non-parametric ANOVA alternative | Independent samples |
| Welch's t-test | Pairwise mean comparison | Does not assume equal variance |
| Mann-Whitney U | Non-parametric pairwise comparison | Independent samples |
| Bonferroni correction | Multiple comparison adjustment | Controls family-wise error rate |
