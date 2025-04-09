# Code

## Preprocessing

Preprocessing, integration, celltype annotation, signature scoring and other data curation scripts.

1. `00_import_merge_seurat_v5.qmd`: Data preprocessing and quality control
2. `01_azimuth_liftover.qmd`: Cell type liftover using azimuth
3. `02_sketch_integration.qmd`: Seurat sketch integration using RPCA
4. `03_celltype_annotation.qmd`: Data clustering and cell type annotation
5. `04_specificity_analysis.qmd`: Cell type specificity analysis
6. `05_ucell_cohort_ranks.qmd`: UCell rank generation
7. `06_score_signatures_ucell.qmd`: UCell pathway scoring
8. `07_generate_ccle_sclc_signaure.qmd`: SCLC signature from CCLE data
9. `08_lung_lineage_genes.qmd`: Get lineage markers from the Human Lung Cell Atlas (HLCA)

## InferCNV Analysis

Code used to run inferCNV and assign malignant cells

9.  `10_run_infercnv.qmd`: Run InferCNV
10. `11_rebin_merge.qmd`: Rebin to 10Mb and merge across sampels
11. `12_infercnv_post.qmd`: InferCNV postprocessing and plotting
12. `13_assign_malignant.qmd`: Assign malignant cells per case

## Epithelial/Tumor Cell Analysis

Epithelial cell subset and histotime analysis on tumor cells

13.  `20_epithelial_object.qmd`: Subset, recluster, and annotate epithelial cell object
14.  `21_run_slingshot.qmd`: Run Slingshot
15.  `22_fit_gams.qmd`: Fit GAMs using TradeSeq
16.  `23_histotime.qmd`: Histotime analysis

## IMPACT Analysis

Code to load and plot IMPACT data

17. `30_load_impact.qmd`: Load IMPACT data
18. `31_impact_followup.qmd`: I think remove this...?

## External Datasets

Analysis code for [LuCA](https://pubmed.ncbi.nlm.nih.gov/36368318/) and TCGA datasets

20.  `40_luca_ranks.qmd`: UCell ranks for the LuCA Atlas data
21.  `41_luca_analysis.qmd`: LuCA Analysis
22.  `42_tcga_analysis.qmd`: TCGA Analysis

## Main Plotting

Core plotting script for figure production

23. `99_plots.qmd`: Main plotting script
