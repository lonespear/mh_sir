# NHANES Depression Severity: Feature Selection & Biological Correlates

_A JSM 2024 Poster Companion | Bayesian Modeling, SIR, & Interpretable Machine Learning in Population Mental Health_

---

## Overview

**Welcome!**  
This repository supports our Joint Statistical Meetings (JSM) 2025 poster:  
**“Feature Selection and Biological Correlates of Latent Depression Severity”**

We leverage modern statistical and machine learning approaches to uncover which **clinical, demographic, and behavioral variables** are most strongly associated with _depression severity_ in the U.S. population.

---

## Project Highlights

- **Data:** 2017–2020 National Health and Nutrition Examination Survey (NHANES)
- **Outcome:** _Latent depression severity_ scored via Bayesian Graded Response Model (BGRM) on PHQ-9 items
- **Methods:**  
    - Sliced Inverse Regression (SIR)
    - Principal Component Analysis (PCA)
    - LASSO Regression
    - XGBoost + SHAP values
- **Key Goals:**  
    - _Rank variables by predictive value_ for depression severity  
    - _Visualize and interpret_ the biological and social drivers of population-level depression

---

## Pipeline Diagram

![Analysis Pipeline Schematic](figures/pipeline_placeholder.png)

---

## Methods

### **1. Bayesian Graded Response Modeling**

- **Why:** Move beyond simple sum scores of PHQ-9 items; capture individual and item-level uncertainty in depression measurement.
- **How:** Fit a BGRM to estimate continuous latent depression (`\(\hat{\theta}_i\)`) for each participant, integrating item difficulties and response patterns.

### **2. Feature Selection & Dimension Reduction**

- **SIR (Sliced Inverse Regression):** Supervised technique to find linear combinations of predictors most strongly associated with depression severity.
- **PCA (Principal Component Analysis):** Unsupervised, finds principal axes of variation across all features.
- **LASSO:** Penalized regression for sparse, interpretable predictor selection.
- **XGBoost + SHAP:** Gradient-boosted trees, with Shapley-value-based importance to rank variable contributions.

### **3. Variable Grouping and Interpretation**

Variables were classified into **five clinically meaningful categories** for interpretation and visualization.

---

## Top Predictors Table

![Top Predictors Table](figures/top_predictors_table.png)

*Each cell is color-coded by variable group. See key below.*

---

## NHANES Variable Reference Table

| **Variable** | **Description**                        | **Group**              |
|--------------|----------------------------------------|------------------------|
| LBDSGBSI     | Globulin, serum                        | Blood Cell Indices     |
| LBXWBCSI     | White blood cell count                 | Blood Cell Indices     |
| LBXLYPCT     | Lymphocyte percent                     | Blood Cell Indices     |
| LBNEXPCT     | Neutrophil percent                     | Blood Cell Indices     |
| LBXHGB       | Hemoglobin                             | Blood Cell Indices     |
| LBDUIBSI     | CBC marker (unspecified)               | Blood Cell Indices     |
| LBXHCT       | Hematocrit                             | Blood Cell Indices     |
| LBDPCT       | Platelet count                         | Blood Cell Indices     |
| LBXMCHSI     | Mean corpuscular hemoglobin (MCH)      | Blood Cell Indices     |
| LBDIRNSI     | Iron saturation or CBC marker          | Blood Cell Indices     |
| LBXMCVSI     | Mean corpuscular volume (MCV)          | Blood Cell Indices     |
| LBXMC        | Monocyte percent/count                 | Blood Cell Indices     |
| LBXRDW       | Red cell distribution width            | Blood Cell Indices     |
| LBDTIB       | Total iron binding capacity            | Blood Cell Indices     |
| LBDBPBSI     | Bicarbonate (chemistry panel)          | Blood Chem/Metabolic   |
| LBDSALSI     | Albumin, serum                         | Blood Chem/Metabolic   |
| LBXHCOT      | Cotinine (smoking marker)              | Blood Chem/Metabolic   |
| LBXSCK       | Creatine kinase (muscle enzyme)        | Blood Chem/Metabolic   |
| BMXWAIST     | Waist circumference                    | Anthropometric         |
| BMXBMI       | Body mass index                        | Anthropometric         |
| BMXWT        | Weight                                 | Anthropometric         |
| BMXARMC      | Arm circumference                      | Anthropometric         |
| BMXHIP       | Hip circumference                      | Anthropometric         |
| BPAOCSZ      | Arm circumference (alternate code)     | Anthropometric         |
| RIDAGEYR     | Age                                    | Demographic/Socioecon  |
| INDFMPIR     | Income-to-poverty ratio                | Demographic/Socioecon  |
| DMDEDUC2     | Education level                        | Demographic/Socioecon  |
| MCQ160A      | Health insurance or chronic illness    | Demographic/Socioecon  |
| MCQ160L      | Physical function limitation           | Demographic/Socioecon  |
| PAQ620       | Physical activity                      | Demographic/Socioecon  |
| PAQ650       | Vigorous activity per week             | Demographic/Socioecon  |
| MCQ300A      | Hypertension diagnosis                 | Demographic/Socioecon  |
| MCQ160M      | Chronic illness (various)              | Demographic/Socioecon  |
| MCQ160P      | Chronic illness (e.g., depression)     | Survey/Medical History |
| DMDMARTZ     | Marital status                         | Survey/Medical History |
| BPQ020       | Told by doctor: hypertension           | Survey/Medical History |
| BPXOSY1/2/3  | Systolic blood pressure                | Survey/Medical History |
| BPXODI1/2    | Diastolic blood pressure               | Survey/Medical History |

---

## References

- Centers for Disease Control and Prevention (CDC), National Center for Health Statistics (NCHS). (2021). [National Health and Nutrition Examination Survey: 2017–2020 Continuous Data.](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?Cycle=2017-2020) Accessed May 30, 2025.
- Foley SF, Kirschbaum C, Grotegerd D, et al. *Peripheral blood cellular immunophenotype in depression: a systematic review and meta-analysis*. Molecular Psychiatry. 2023.
- Golub, G. H., & Van Loan, C. F. (2013). *Matrix Computations* (4th ed.). Johns Hopkins University Press.
- Hansen, P. C. (1987). The truncated SVD as a method for regularization. *BIT Numerical Mathematics*, 27(4), 534–553.
- Köhler-Forsberg O, Buttenschøn HN, Tansey KE, et al. *White blood cell count at first depression diagnosis as predictor for risk of subsequent hospitalization*. Neurology, Psychiatry and Brain Research. 2017.
- Kroenke, K., Spitzer, R. L., & Williams, J. B. W. (2001). [The PHQ-9: Validity of a Brief Depression Severity Measure. *Journal of General Internal Medicine*, 16(9), 606–613.](https://doi.org/10.1046/j.1525-1497.2001.016009606.x)
- Kuhn, M. (2023). [caret: Classification and Regression Training](https://cran.r-project.org/package=caret) (R package manual).
- Li, K.-C. (1991). [Sliced inverse regression for dimension reduction. *Journal of the American Statistical Association*, 86(414), 316–327.](https://www.tandfonline.com/doi/abs/10.1080/01621459.1991.10475035)
- Mazza MG, Lucchi S, Tringali A, et al. *Neutrophil/lymphocyte ratio and platelet/lymphocyte ratio in mood disorders: a meta-analysis*. Prog Neuropsychopharmacol Biol Psychiatry. 2018;84:229–236.
- McCullagh, P. (1980). [Regression models for ordinal data. *Journal of the Royal Statistical Society: Series B (Methodological)*, 42(2), 109–127.](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.2517-6161.1980.tb01109.x)
- Osimo EF, Baxter LJ, Lewis G, et al. *Prevalence of low-grade inflammation in depression: a systematic review and meta-analysis of CRP levels*. Psychological Medicine. 2019;49(11):1958–1970.
- Samejima, F. (1969). Estimation of latent ability using a response pattern of graded scores. *Psychometrika Monograph Supplement*, 34(4, Pt.2), 100.
- Sealock JM, Davis LK, Gymrek M, et al. *Depression polygenic scores and white blood cell count*. JAMA Psychiatry. 2021;78(12):1340–1350.
- Tikhonov, A. N., & Arsenin, V. Y. (1977). *Solutions of Ill-posed Problems*. Wiley.

*Additional references and code details available upon request.*

---

## About

**Author:** Jonathan Day 
**Contact:** jonathan.day@westpoint.edu
**Bio:** https://www.westpoint.edu/jonathan-l-day

For research/educational use only. Not medical advice.  
_Last updated: [2025-07-28]_

---

