#' @export
load_project <- function(cohort_obj = "/data1/shahs3/users/zatzmanm/work/repos/egfr/data/objects/egfr_cohort.rds",
                         epi_obj = "/data1/shahs3/users/zatzmanm/work/repos/egfr/data/objects/egfr_epi_v2.rds",
                         tumor_obj = "/data1/shahs3/users/zatzmanm/work/repos/egfr/data/objects/egfr_tumor_cells.rds",
                         ondisk_counts_dir = "/data1/shahs3/users/zatzmanm/work/repos/egfr/data/bp_cells/merged_counts/",
                         infercnv_calls = "/data1/shahs3/users/zatzmanm/work/repos/egfr/results/2024-12-15-egfr_infercnv_results.txt.gz",
                         infercnv_obj = "/data1/shahs3/users/zatzmanm/work/repos/egfr/cnv/infercnv/objects/egfr_cohort_cnv.rds",
                         metadata_sheet = "/data1/shahs3/users/zatzmanm/work/repos/egfr/metadata/metadata_update_mar_2025.txt",
                         impact_obj = "/data1/shahs3/users/zatzmanm/work/repos/egfr/results/impact_data.rds",
                         ct_spec = "/data1/shahs3/users/zatzmanm/work/repos/egfr/results/cell_type_patient_specificity.txt.gz",
                         tp_spec = "/data1/shahs3/users/zatzmanm/work/repos/egfr/results/cell_type_timepoint_specificity.txt.gz") {
  cli::cli_alert_info("Loading EGFR Project data")

  cli::cli_dl(c(
    "Cohort Object" = "{.path {cohort_obj}}",
    "Epithelial Object" = "{.path {epi_obj}}",
    "Tumor Object" = "{.path {tumor_obj}}",
    "On Disk Counts Location" = "{.path {ondisk_counts_dir}}",
    "InferCNV Calls" = "{.path {infercnv_calls}}",
    "InferCNV Object" = "{.path {infercnv_obj}}",
    "Sample Metadata" = "{.path {metadata_sheet}}",
    "IMPACT Object" = "{.path {impact_obj}}"
  ))

  cmap <- load_colors()

  srt <- readRDS(cohort_obj)

  # metadata <- load_metadata()
  metadata <- read.table(metadata_sheet, header = T, sep = "\t", quote = "\"")
  metadata$time_point <- factor(metadata$time_point, levels = c("TN", "MRD", "PD"))

  srt <- update_metadata(srt, metadata, match_by = "sample_id")

  # Adjust the directory
  srt[["RNA"]]$counts <- BPCells::open_matrix_dir(ondisk_counts_dir)

  srt <- Seurat::NormalizeData(srt)

  SeuratObject::Idents(srt) <- "time_point"

  srt <- Seurat::UpdateSeuratObject(srt)

  # Fix data structures
  # srt@graphs$RNA_nn <- srt@graphs$RNA_nn[colnames(srt), colnames(srt)]
  # srt@graphs$RNA_snn <- srt@graphs$RNA_snn[colnames(srt), colnames(srt)]
  # srt@reductions$RNA_rpca_integrated@cell.embeddings <- srt@reductions$RNA_rpca_integrated@cell.embeddings[colnames(srt), ]
  # srt@reductions$RNA_umap_rpca_integrated@cell.embeddings <- srt@reductions$RNA_umap_rpca_integrated@cell.embeddings[colnames(srt), ]

  tum_calls <- read.table(file = infercnv_calls, header = T, sep = "\t")
  rownames(tum_calls) <- tum_calls$cell_id

  srt <- Seurat::AddMetaData(srt, tum_calls)

  srt$cell_type_super <- srt@meta.data %>%
    dplyr::mutate(
      cell_type_super =
        case_when(RNA_cluster_labels_coarse == "Mast Cell" ~ "Myeloid", .default = RNA_cluster_labels_coarse)
    ) %>%
    dplyr::pull(cell_type_super) %>%
    factor()

  srt$time_point <- factor(srt$time_point, levels = names(cmap$time_point))

  srt$cell_type <- srt$RNA_cluster_labels_coarse

  srt$cell_type[srt$is_tumor_cell] <- "Tumor Cell"

  # Split out Plasma from B Cells
  srt$cell_type[srt$RNA_cluster_labels == "B-Cell (Plasma)"] <- "Plasma cell"

  srt$cell_type <- factor(srt$cell_type)

  new_meta <- srt@meta.data %>%
    dplyr::count(cell_type) %>%
    dplyr::mutate(cell_type_nlab = factor(glue::glue("{cell_type} (n={prettyNum(n, big.mark = ',')})"))) %>%
    dplyr::right_join(srt@meta.data)

  rownames(new_meta) <- new_meta$cell_id

  srt@meta.data <- new_meta

  SeuratObject::Idents(srt) <- "sample_id"

  srt <- SeuratObject::UpdateSeuratObject(srt)

  # srt[["umap"]] <- srt[["RNA_umap_rpca_integrated"]]
  # srt[["RNA_umap_rpca_integrated"]] <- NULL


  srt$cell_type <- factor(srt$cell_type, levels = names(cmap$cell_type))
  srt$site_of_tissue_simple <- factor(srt$site_of_tissue_simple, levels = names(cmap$site_of_tissue_simple))
  x <- table(distinct(srt@meta.data[, c("sample_id_new", "time_point")])$time_point)

  srt$sample_id_new <- factor(srt$sample_id_new,
    levels = c(
      paste0("TN", 1:x[["TN"]]),
      paste0("MRD", 1:x[["MRD"]]),
      paste0("PD", 1:x[["PD"]])
    )
  )

  # Update ordering and color_map with the labeled
  cmap$cell_type_nlab <- cmap$cell_type
  names(cmap$cell_type_nlab) <- unique(srt$cell_type_nlab)[pmatch(levels(srt$cell_type), unique(srt$cell_type_nlab))]

  srt$cell_type_nlab <- factor(srt$cell_type_nlab, levels = names(cmap$cell_type_nlab))

  markers <- load_markers()

  # Do this again...?
  # srt@graphs$RNA_nn <- srt@graphs$RNA_nn[colnames(srt), colnames(srt)]
  # srt@graphs$RNA_snn <- srt@graphs$RNA_snn[colnames(srt), colnames(srt)]

  impact_data <- readRDS(impact_obj)

  onco_genes <- impact_data$oncokb_dat$oncokb_summary_table %>%
    filter(gene %in% c("RB1", "TP53")) %>%
    select(dmp_sample, gene, variant) %>%
    pivot_wider(id_cols = dmp_sample, names_from = gene, values_from = variant, values_fn = function(x) paste(x, collapse = ", "), names_prefix = "variant_") %>%
    left_join(metadata[, c("sample_id", "closest_impact_to_sample")], by = c("dmp_sample" = "closest_impact_to_sample")) %>%
    dplyr::rename("closest_impact_to_sample" = "dmp_sample")

  metadata <- left_join(metadata, onco_genes) %>%
    mutate(
      TP53mut = !is.na(variant_TP53),
      RB1mut = !is.na(variant_RB1)
    )

  metadata <- dplyr::filter(metadata, sample_id %in% srt$sample_id)

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

    # histotime <- read.table(file = here("results", "histotime_cells.txt"), sep = "\t", header = T) %>%
    # column_to_rownames("cell_id")
    # srt_tum <- AddMetaData(srt_tum, histotime)

    # Update metadata on tumor cell object
    srt_epi <- update_metadata(srt_epi, metadata, match_by = "sample_id")

    srt_epi$sample_id_new <- factor(srt_epi$sample_id_new, levels = levels(srt$sample_id_new))
    srt_epi$site_of_tissue_simple <- factor(srt_epi$site_of_tissue_simple, levels = levels(srt$site_of_tissue_simple))
    srt_epi$cell_type_epi <- fct_recode(srt_epi$cell_type_epi, "Atypical" = "PD9/PD27 unk.")

    # Rotate UMAP
    # Since we resaved this is no longer needed
    # rot = -160
    # srt_epi[["umap"]]@cell.embeddings[,1:2] <- spdep::Rotation(srt_epi[["umap"]]@cell.embeddings[,1:2], angle = rot*(pi/180))
  } else {
    srt_epi <- NULL
  }

  if (!is.null(tumor_obj)) {
    srt_tum <- readRDS(tumor_obj) %>%
      Seurat::UpdateSeuratObject()

    histotime <- read.table(file = here("results", "histotime_cells.txt"), sep = "\t", header = T) %>%
      column_to_rownames("cell_id")
    srt_tum <- AddMetaData(srt_tum, histotime)

    # Update metadata on tumor cell object
    srt_tum <- update_metadata(srt_tum, metadata, match_by = "sample_id")
  } else {
    srt_tum <- NULL
  }

  db <- list(
    srt = srt,
    srt_tum = srt_tum,
    srt_epi = srt_epi,
    infercnv = infercnv,
    colors = cmap,
    markers = markers,
    metadata = metadata,
    impact_data = impact_data,
    spec = sp_tabs
  )


  cli::cli_alert_success("EGFR Project Loaded Successfully!")
  cli::cli_text("DB objects: {.pkg {names(db)}}")

  return(db)
}
