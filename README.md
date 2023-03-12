# BC_subtypes_HSCBB22
Binary Classification of Breast Cancer subtype groups based on lncRNA expression.

# Dataset Preparation
The path: "./PrepareDatasets" contains four scripts. The first one combines patient and sample data from TCGA-BRCA cohort (cBioportal). The second collects information about hereditary mutations in the cohort and defines samples with sporadic mutations. The third one is responsible for the filtering of the clinical table by keeping only sporadic cases. The fourth and last script modifies the gene expression dataset and combines it with the filtered clinical table.

# Analysis
The analysis consists of multiple parts that involve Binary Classification, Feature Reduction, Enrichment Analysis, Survival Analysis. The path: "./Analysis/ML" contains the R script responsible for the definition of subtype groups, binary classification via training a LASSO model, hyperparameter tuning, feature selection and storage of informative features for downstream analysis.
