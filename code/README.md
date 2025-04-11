# Code

Many scripts here are meant to demonstrate how the data were processed, and are not intended to be rerunnable as is. Scripts that are possible to rerun are indicated in **`bold`**, and are meant to demonstrate figure reproduction.

## Download and convert objects

0. **`0_download_convert.qmd`**: Code to take H5AD objects downloaded from [CELLxGENE](https://cellxgene.cziscience.com) and convert into Seurat objects for [our study](https://cellxgene.cziscience.com/collections/df8ef04c-7a87-49d4-a127-f8d653d89d91), the [LuCA atlas](https://cellxgene.cziscience.com/collections/edb893ee-4066-4128-9aec-5eb2b03f8287). This script will also ensure other required data objects are available.

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

10.  `10_run_infercnv.qmd`: Run InferCNV
11. `11_rebin_merge.qmd`: Rebin to 10Mb and merge across sampels
12. `12_infercnv_post.qmd`: InferCNV postprocessing and plotting
13. `13_assign_malignant.qmd`: Assign malignant cells per case

## Epithelial/Tumor Cell Analysis

Epithelial cell subset and histotime analysis on tumor cells

14.  `20_epithelial_object.qmd`: Subset, recluster, and annotate epithelial cell object
15.  `21_run_slingshot.qmd`: Run Slingshot
16.  `22_fit_gams.qmd`: Fit GAMs using TradeSeq
17.  **`23_histotime.qmd`**: Histotime analysis

## IMPACT Analysis

Code to load and plot IMPACT data

18. **`30_load_impact.qmd`***: Load IMPACT data

## External Datasets

Analysis code for [LuCA](https://pubmed.ncbi.nlm.nih.gov/36368318/) and TCGA LUAD datasets.

1.   `40_luca_ranks.qmd`: UCell ranks for the LuCA Atlas data
2.   **`41_luca_analysis.qmd`**: LuCA Analysis
3.   **`42_tcga_cptac_analysis.qmd`**: TCGA and CPTAC bulk RNA dataset analysis

## Main Plotting

Core plotting script for figure production

22. **`99_plots.qmd`**: Main plotting script
