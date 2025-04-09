#' @export
load_metadata <- function(nucseq_excel = "/data1/shahs3/users/zatzmanm/work/repos/egfr/metadata/Nucseq_Shah_Lab_2023_5_24.xlsx", version = "V1", meta_dir = "/data1/shahs3/users/zatzmanm/work/repos/egfr/metadata/", save = FALSE) {
  # library(tidyverse)
  # library(readxl)
  # library(glue)
  # library(janitor)

  metadata <- suppressWarnings({
    readxl::read_excel(nucseq_excel,
      col_types = c(
        "text", "text", "date",
        "text", "skip", "skip", "text", "text",
        "date", "text", "date", "text", "text", "text", "text", "text", "text",
        "text", "text", "numeric", "numeric",
        "text", "text", "text", "text", "text",
        "text", "text", "text", "date", "text",
        "text", "text", "text", "text", "date",
        "text", "text", "text", "text", "text",
        "date", "date", "text", "text", "text",
        "text", "text", "text", "text", "text",
        "text", "text", "text", "text", "text",
        "text", "skip", "skip", "skip"
      ), na = c("", "#N/A", "N/A", "None", "na")
    ) %>%
      janitor::clean_names()
  })

  # Remove all date time variables
  dt_cols <- which(purrr::map_df(metadata, class)[1, ] == "POSIXct")
  metadata[, dt_cols] <- metadata[, dt_cols] %>%
    purrr::map_df(., \(x) as.character(as.Date(x)))

  metadata <- metadata %>%
    dplyr::mutate(
      sample_id = gsub("( |\\().*", "", stringr::str_to_title(project_id)),
      progressed = ifelse(grepl("Yes", patients_egfr_lc_progressed_ever_y_n), T, F),
      time_point = gsub(" .*", "", time_point),
      mapk_sustained_genomic = mapk_sustained_genomic == "y",
      histological_transformation_clinical_simple = histological_transformation_clinical != "",
      mech_resistance = ifelse(grepl(
        "Other*|Unknown",
        mech_of_resistance_category_by_alvaro
      ),
      "Other", mech_of_resistance_category_by_alvaro
      ),
      impact_pid = substr(metadata$impact_for_sample, 1, 9),
      impact_for_sample_clean = substr(metadata$impact_for_sample, 1, 17),
      .before = "project_id"
    )

  metadata$has_impact <- ifelse(!is.na(metadata$impact_for_sample_clean), T, F)
  metadata$has_nucseq <- T

  metadata$pass_scrna_qc <- T

  metadata$pass_scrna_qc[metadata$sample_id %in% c("Hpre12", "Hon7", "Hpost30", "Hpost39", "Ru1521b")] <- F

  metadata[metadata$patients_egfr_lc_progressed_ever_y_n == "Patient deceased from unknown cause", "progressed"] <- NA

  # If no impact and no annotation set to NA as this is potentially unevaluable
  metadata$mapk_sustained_genomic[is.na(metadata$mapk_sustained_genomic)] <- F

  metadata$mapk_sustained_genomic[metadata$has_impact == F & metadata$mapk_sustained_genomic == F] <- NA

  metadata <- metadata %>%
    dplyr::mutate(tissue_site = dplyr::case_when(
      grepl(x = site_of_tissue, "Primary", ignore.case = T) ~ "Primary",
      grepl(x = site_of_tissue, "bone|back|neck|sacrum", ignore.case = T) ~ "Bone/soft tissue",
      grepl(x = site_of_tissue, "LN", ignore.case = T) ~ "Lymph node",
      grepl(x = site_of_tissue, "Pleura", ignore.case = T) ~ "Pleura",
      grepl(x = site_of_tissue, "RLL|chest wall", ignore.case = T) ~ "Locoregional",
      grepl(x = site_of_tissue, "brain", ignore.case = T) ~ "Brain",
      grepl(x = site_of_tissue, "liver", ignore.case = T) ~ "Liver",
      .default = "other"
    ))

  # metadata$sample_id <- gsub("Av-1517_r", "R", metadata$sample_id)


  metadata <- metadata %>%
    dplyr::mutate_if(is.character, as.factor)

  metadata$time_point <- factor(metadata$time_point, levels = c("Pre", "On", "Progression"))

  metadata$sample_id <- factor(metadata$sample_id, levels = gtools::mixedsort(unique(as.character(metadata$sample_id))))

  if (save) {
    outfile <- file.path(meta_dir, glue::glue("egfr_metadata_{version}.txt"))

    cli::cli_alert_info("Saving to {.path {outfile}}")
    write.table(metadata, file = outfile, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }


  return(metadata)
}

#' @export
update_metadata <- function(obj, metadata, match_by) {
  new_meta <- obj@meta.data

  cols <- colnames(metadata)[!colnames(metadata) == match_by]
  new_meta[, cols] <- NULL

  new_meta <- dplyr::left_join(new_meta, metadata, by = match_by)
  rownames(new_meta) <- new_meta$cell_id

  obj@meta.data <- new_meta
  return(obj)
}
