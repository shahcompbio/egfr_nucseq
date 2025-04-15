#' @export
load_project <- function(cohort_obj = here::here("paper_data/egfr_cohort.rds"),
                         epi_obj = here::here("paper_data/egfr_epi.rds"),
                         histo_obj = here::here("paper_data/egfr_histo.rds"),
                         infercnv_obj = here::here("paper_data/egfr_cohort_cnv.rds"),
                         # metapaper_data_sheet = "paper_data/supplemental_tables.xlsx",
                         # impact_obj = "impact_data.rds",
                         ct_spec = here::here("paper_data/cell_type_patient_specificity.txt.gz"),
                         tp_spec = here::here("paper_data/cell_type_timepoint_specificity.txt.gz")) {
  cli::cli_alert_info("Loading EGFR Project data")

  cli::cli_dl(c(
    "Cohort Object" = "{.path {cohort_obj}}",
    "Epithelial Object" = "{.path {epi_obj}}",
    "InferCNV Object" = "{.path {infercnv_obj}}"
    # "Sample Metadata" = "{.path {metadata_sheet}}",
    # "IMPACT Object" = "{.path {impact_obj}}"
  ))

  # cmap <- load_colors()

  srt <- readRDS(cohort_obj)

  # srt <- Seurat::UpdateSeuratObject(srt)

  # srt$time_point <- factor(srt$time_point, levels = names(cmap$time_point))

  # srt$cell_type <- factor(srt$cell_type, levels = names(cmap$cell_type))
  # srt$site_of_tissue_simple <- factor(srt$site_of_tissue_simple, levels = names(cmap$site_of_tissue_simple))

  # markers <- load_markers()

  # impact_data <- readRDS(impact_obj)

  ct_sp <- read.table(file = ct_spec, header = T, sep = "\t")
  tp_sp <- read.table(file = tp_spec, header = T, sep = "\t")
  # Batch corrected versions
  ct_b <- read.table(file.path(dirname(ct_spec), "cell_type_patient_specificity_batch.txt.gz"), header = T, sep = "\t")
  tp_b <- read.table(file.path(dirname(ct_spec), "cell_type_timepoint_specificity_batch.txt.gz"), header = T, sep = "\t")

  sp_tabs <- list(
    "cell_type" = ct_sp,
    "timepoint" = tp_sp,
    "cell_type_batch" = ct_b,
    "timepoint_batch" = tp_b
  )

  infercnv <- readRDS(infercnv_obj)

  if (!is.null(epi_obj)) {
    srt_epi <- readRDS(epi_obj) %>%
      Seurat::UpdateSeuratObject()
    # Temporary patch. The epi object needs to have this variable up front
    srt_epi$cell_type_epi <- srt$cell_type_epi[match(colnames(srt_epi), colnames(srt))]
  } else {
    srt_epi <- NULL
  }

  if (!is.null(epi_obj)) {
    srt_histo <- readRDS(histo_obj) %>%
      Seurat::UpdateSeuratObject()
  } else {
    srt_histo <- NULL
  }

  db <- list(
    srt = srt,
    srt_epi = srt_epi,
    srt_histo = srt_histo,
    # infercnv = infercnv,
    # metadata = metadata,
    # impact_data = impact_data,
    spec = sp_tabs
  )

  cli::cli_alert_success("EGFR Project Loaded Successfully!")
  cli::cli_text("DB objects: {.pkg {names(db)}}")

  return(db)
}
