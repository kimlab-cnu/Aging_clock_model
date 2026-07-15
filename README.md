# Aging clock model
Aging research has largely focused on changes in protein expression levels, but aging cannot be explained by abundance changes alone. Here we focus on an under-explored dimension: structural changes in proteins and the dynamics of their interaction networks.
In this study, we integrated experimental proteomics with computational structural biology to comprehensively characterize the age-dependent structural dynamics of plasma proteins — that is, how plasma proteins change structurally over a healthy life course.

## 1. Data Preprocessing

Data_preprocessing_for_PAR.R
Computes the Peak Area Ratio (PAR) for each ion, using the 6 ions from the 1st, 2nd, and 3rd precursors of peptide-03 among the internal standards (IS), in order to select the IS used for normalization.
(Fig. 1B — code added 2025-10-19)

PAR_distribution_histogram.R
Applies a log10 transformation to the PAR1–6 values obtained above and plots their distribution as a histogram, together with the corresponding CV values.
(Fig. 1A, C — code added 2025-10-19)

PAR_input_preprocessing.R

Aging_log_norm_data_preprocessing.ipynb

PAR_data_preprocessing_Jhyh.R
Preprocesses the data required for the downstream analyses.

## 2. Statistical Analysis

Data_preprocessing_Jhyh+cosine similarity.R
Performs cosine similarity analysis on the preprocessed data.
(Fig. 2A, B)

limma_test_Jh_revised_added_fit3_final.R
Uses limma to extract only the ions satisfying the fold-change and p-value criteria for each pattern, and plots the mean PAR change per pattern.
(Fig. 3A, B)

overlap_pattern_fit3_final.R
Merges the cosine similarity and limma results and retains only the overlapping ions.

overlap_boxplot_fit3_final.R
Generates boxplots for each ion, grouped by pattern.
(Fig. 3C)

## 3. Structural Analysis

pymol_interaction_distance.txt
When inspecting AlphaFold3 multimer results in PyMOL, measures the minimum distance across the peptide– and protein–interface within each pattern's protein complex (1:1).

pymol_mean_plddt.txt
When inspecting AlphaFold3 multimer results in PyMOL, computes the mean pLDDT.
(Fig. 6B, C, D)


## Notes
Cytoscape: Full STRING network, confidence ≥ 0.7, Homo sapiens (see manual for details).
AF3 multimer web server: For PTMs, only N-linked glycosylation sites were added as NAG.
