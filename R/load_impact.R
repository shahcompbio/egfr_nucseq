#' Load IMPACT Data
#'
#' @param dmp_ids Vector of DMP IDs (e.g. 'P-0035733-T02-IM6')
#' @param impact_dir IMPACT directory location
#' @param save_to Optional. Path to save the impact object rds (e.g. 'data/impact_data.rds')
#'
#' @return A list containing: clin_dat_patient, clin_dat_sample, onco_meta, segments, data_cna, mutations, sv, onco_dat_full
#' @export
load_impact <- function(dmp_ids,
                        impact_dir = "/data1/shahs3/users/zatzmanm/work/repos/shared/msk-impact/",
                        onckb_dir = "/data1/shahs3/users/zatzmanm/work/repos/shared/oncokb-annotated-msk-impact",
                        facets_data = NULL,
                        save_to = NULL) {
  dmp_ids <- na.omit(as.character(dmp_ids))

  dmp_pids <- substr(dmp_ids, 1, 9) %>%
    unique()

  cli::cli_alert_info("Looking for IMPACT data for {length(dmp_ids)} sample{?s} in {length(dmp_pids)} patient{?s}. ")

  # Patient level clinical data
  clin_dat_patient <- read.table(file = file.path(impact_dir, "msk_solid_heme/data_clinical_patient.txt"), header = T, sep = "\t", quote = "\"", strip.white = T, blank.lines.skip = T) %>%
    dplyr::rename(dmp_pid = PATIENT_ID) %>%
    dplyr::filter(dmp_pid %in% dmp_pids)

  if (all(dmp_pids %in% clin_dat_patient$dmp_pid)) {
    cli::cli_alert_success("Found all patients")
  } else {
    missing_pts <- setdiff(dmp_pids, clin_dat_patient$dmp_pid)
    cli::cli_alert_warning("Missing {length(missing_pts)} patient{?s}: {.val {missing_pts}}")
  }

  # Sample level metadata
  clin_dat_sample <- data.table::fread(file.path(impact_dir, "msk_solid_heme/data_clinical_sample.txt"), skip = 4)
  clin_dat_sample <- clin_dat_sample[clin_dat_sample$PATIENT_ID %in% dmp_pids, ] %>%
    as.data.frame() %>%
    dplyr::rename(
      dmp_pid = PATIENT_ID,
      dmp_sample = SAMPLE_ID
    ) %>%
    dplyr::mutate(dmp_tumor = str_sub(dmp_sample, 1, 13)) %>%
    dplyr::relocate(dmp_pid, dmp_tumor, dmp_sample)

  if (all(dmp_ids %in% clin_dat_sample$dmp_sample)) {
    cli::cli_alert_success("Found all samples")
  } else {
    missing_samps <- setdiff(dmp_ids, clin_dat_sample$dmp_sample)
    cli::cli_alert_warning("Missing {length(missing_samps)} sample{?s}: {.val {missing_samps}}")
  }

  # Mutations
  mutations <- data.table::fread(file.path(impact_dir, "msk_solid_heme/data_mutations_extended.txt"), skip = 1, showProgress = F)
  mutations <- mutations[mutations$Tumor_Sample_Barcode %in% clin_dat_sample$dmp_sample, ] %>%
    as.data.frame() %>%
    dplyr::mutate(
      dmp_pid = str_sub(Tumor_Sample_Barcode, 1, 9),
      dmp_tumor = str_sub(Tumor_Sample_Barcode, 1, 13),
      dmp_sample = Tumor_Sample_Barcode,
      .before = 1
    ) %>%
    dplyr::mutate(source = ifelse(grepl("XS", dmp_sample), "access", "impact"))

  # Segment data
  segments <- readr::read_table(file.path(impact_dir, "msk_solid_heme/mskimpact_data_cna_hg19.seg"), col_types = readr::cols())
  segments$dmp_pid <- substr(segments$ID, 1, 9)
  segments <- filter(segments, dmp_pid %in% dmp_pids) %>%
    dplyr::rename(dmp_sample = ID) %>%
    dplyr::mutate(dmp_tumor = substr(dmp_sample, 1, 13)) %>%
    dplyr::relocate(dmp_pid, dmp_tumor, dmp_sample)
  segments$chrom <- paste0("chr", segments$chrom)

  # Gene level CNA data
  data_cna <- data.table::fread(file.path(impact_dir, "msk_solid_heme/data_CNA.txt"), showProgress = F)
  sample_cols <- purrr::map(dmp_pids, \(pid) grep(pid, colnames(data_cna))) %>%
    list_c()
  sample_cols <- sort(unique(sample_cols[!is.na(sample_cols)]))

  data_cna <- data_cna[, c(1, sample_cols), with = FALSE] %>%
    tidyr::pivot_longer(cols = !Hugo_Symbol, names_to = "dmp_sample") %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::rename(Gene = Hugo_Symbol) %>%
    dplyr::mutate(cna_type_detail = case_when(
      value >= 1.5 ~ "Amplification",
      value >= 1 ~ "Gain",
      value == 0 ~ "Neutral",
      value <= -1.5 ~ "Deep Deletion",
      value <= -1 ~ "Shallow Deletion"
    )) %>%
    dplyr::mutate(
      dmp_pid = substr(dmp_sample, 1, 9),
      dmp_tumor = substr(dmp_sample, 1, 13),
      source = "impact_cna",
    ) %>%
    dplyr::filter(value != 0) %>%
    dplyr::select(dmp_pid, dmp_tumor, dmp_sample, Gene, mutation_type = cna_type_detail, source, cna_value = value)


  # Facets CNA data
  if (!is.null(facets_data)) {
    cli::cli_alert_info("Loading facets data from {.path {facets_data}}")
    facets_df <- read.table(file = facets_data, header = T, sep = "\t")
  }


  # SV data
  sv <- data.table::fread(file.path(impact_dir, "msk_solid_heme/data_sv.txt"), showProgress = F)
  sv <- sv[sv$Sample_ID %in% clin_dat_sample$dmp_sample, ] %>%
    as.data.frame() %>%
    dplyr::mutate(
      dmp_pid = substr(Sample_ID, 1, 9),
      dmp_tumor = substr(Sample_ID, 1, 13),
      dmp_sample = Sample_ID,
      .before = 1
    ) %>%
    dplyr::filter(Site1_Hugo_Symbol != "" | Site2_Hugo_Symbol != "") %>%
    dplyr::mutate(
      source = ifelse(grepl("Archer", Event_Info), "archer",
        ifelse(grepl("XS", dmp_sample), "access", "impact")
      ),
      Class = ifelse(source == "archer", "RNA_FUSION", Class)
    ) %>%
    dplyr::mutate(
      Gene1 = HGNChelper::checkGeneSymbols(Site1_Hugo_Symbol)$Suggested.Symbol,
      Gene2 = HGNChelper::checkGeneSymbols(Site2_Hugo_Symbol)$Suggested.Symbol,
      .before = Site1_Hugo_Symbol
    )

  # Merge into a table
  mut_onco <- mutations %>%
    dplyr::select(dmp_pid, dmp_tumor, dmp_sample, source, Gene = Hugo_Symbol, mutation_type = Variant_Classification, t_ref_count, t_alt_count, n_ref_count, n_alt_count, HGVSp, HGVSp_Short) %>%
    dplyr::mutate(t_vaf = t_alt_count / (t_ref_count + t_alt_count)) %>%
    dplyr::relocate(t_vaf, .before = t_ref_count)

  sv_onco <- sv %>%
    dplyr::select(dmp_pid, dmp_tumor, dmp_sample, Gene1, Gene2, mutation_type = Class, source, sv_event_info = Event_Info) %>%
    dplyr::mutate(sv_single_gene_event = ifelse(is.na(Gene2), T, F)) %>%
    dplyr::mutate(
      sv_id = ifelse(sv_single_gene_event, Gene1, paste(Gene1, Gene2, sep = "_")),
    ) %>%
    tidyr::pivot_longer(cols = c("Gene1", "Gene2"), names_to = "sv_gene", values_to = "Gene") %>%
    dplyr::filter(Gene != "")

  # Bind all and add tumor sample identifier
  onco_dat_full <- purrr::list_rbind(list(mut_onco, data_cna, sv_onco))

  onco_dat_full$mutation_type <- janitor::make_clean_names(onco_dat_full$mutation_type, allow_dupes = TRUE, case = "title")


  # Map mutations
  onco_dat_full <- onco_dat_full %>%
    dplyr::mutate(mutation_type_simple = dplyr::case_when(
      grepl("In Frame", mutation_type) ~ "Inframe Indel",
      grepl("Frame Shift", mutation_type) ~ "Frameshift Indel",
      (!is.na(sv_id) & mutation_type != "Rna Fusion") ~ "Structural Variant",
      grepl("Translation Start|Splice Site|X5Flank", mutation_type) ~ "Other",
      grepl("In Frame|Nonframeshift Deletion", mutation_type) ~ "Inframe Indel",
      grepl("Nonsynonymous Snv", mutation_type) ~ "Missense Mutation",
      mutation_type == "Rna Fusion" ~ "RNA Fusion",
      .default = mutation_type
    ))

  # Load oncokb
  oncokb_dat <- load_oncokb(
    oncokb_dir = "/data1/shahs3/users/zatzmanm/work/repos/shared/oncokb-annotated-msk-impact",
    dmp_ids = dmp_ids
  )

  impact_obj <- list(onco_dat_full = onco_dat_full, clin_dat_patient = clin_dat_patient, clin_dat_sample = clin_dat_sample, segments = segments, data_cna = data_cna, mutations = mutations, sv = sv, facets_data = facets_df, oncokb_dat = oncokb_dat)

  if (!is.null(save_to)) {
    cli::cli_alert_info("Saving to {.path {save_to}}")
    saveRDS(impact_obj, file = save_to)
  }

  return(impact_obj)
}


load_facets <- function(facets_path = "/rtsess01/compute/juno/cmo/juno/work/ccs/shared/resources/impact/facets/",
                        dmp_ids,
                        save_to = NULL) {
  dmp_ids <- na.omit(as.character(dmp_ids))

  dmp_pids <- substr(dmp_ids, 1, 9) %>%
    unique()

  cli::cli_alert_info("Looking for FACETS data for {length(dmp_ids)} sample{?s} in {length(dmp_pids)} patient{?s}.")
  # Missing the shared facets data
  facets_results <- set_names(dmp_pids) %>%
    map(\(x) dir(glue("{facets_path}/all/{str_sub(x, 1, 7)}"), pattern = glue("{x}"), full.names = TRUE)) %>%
    list_c()

  res <- map(facets_results, \(f) {
    df <- read.table(file.path(f, "facets_qc.txt"), header = T, sep = "\t", stringsAsFactors = FALSE) %>%
      mutate(across(everything(), .fns = as.character))
  }) %>%
    list_rbind() %>%
    type_convert() %>%
    mutate(
      dmp_pid = str_sub(tumor_sample_id, 1, 9),
      dmp_tumor = str_sub(tumor_sample_id, 1, 13),
      dmp_sample = str_sub(tumor_sample_id, 1, 17),
      .before = 1
    ) %>%
    filter(is_best_fit) %>%
    filter(dmp_sample %in% dmp_ids) %>%
    distinct(dmp_sample, .keep_all = T)

  missing_samps <- setdiff(dmp_ids, res$dmp_sample)


  cli::cli_alert_warning("Missing {length(missing_samps)} sample{?s}: {.val {missing_samps}}")


  if (!is.null(save_to)) {
    write.table(res, file = save_to, quote = F, sep = "\t", col.names = T, row.names = F)
  }

  return(res)
}


load_oncokb <- function(oncokb_dir = "/data1/shahs3/users/zatzmanm/work/repos/shared/oncokb-annotated-msk-impact", dmp_ids) {
  # Get ONCOKB annotated mutations
  oncokb_clindat <- data.table::fread(file.path(oncokb_dir, "data_clinical_sample.oncokb.txt.gz")) %>%
    as.data.frame() %>%
    mutate(
      dmp_pid = str_sub(SAMPLE_ID, 1, 9),
      dmp_tumor = str_sub(SAMPLE_ID, 1, 13),
      dmp_sample = str_sub(SAMPLE_ID, 1, 17),
      .before = 1
    ) %>%
    filter(dmp_sample %in% dmp_ids) %>%
    janitor::clean_names()

  oncokb_cna <- data.table::fread(file.path(oncokb_dir, "data_CNA.oncokb.txt.gz")) %>%
    as.data.frame() %>%
    mutate(
      dmp_pid = str_sub(SAMPLE_ID, 1, 9),
      dmp_tumor = str_sub(SAMPLE_ID, 1, 13),
      dmp_sample = str_sub(SAMPLE_ID, 1, 17),
      .before = 1
    ) %>%
    filter(dmp_sample %in% dmp_ids) %>%
    janitor::clean_names()

  # Some are NAs but keeping as they indicate that archer was run
  oncokb_sv <- data.table::fread(file.path(oncokb_dir, "data_sv.oncokb.txt.gz")) %>%
    as.data.frame() %>%
    mutate(
      dmp_pid = str_sub(Sample_ID, 1, 9),
      dmp_tumor = str_sub(Sample_ID, 1, 13),
      dmp_sample = str_sub(Sample_ID, 1, 17),
      .before = 1
    ) %>%
    filter(dmp_sample %in% dmp_ids) %>%
    janitor::clean_names()

  oncokb_muts <- data.table::fread(file.path(oncokb_dir, "data_mutations_extended.oncokb.txt.gz")) %>%
    as.data.frame() %>%
    mutate(
      dmp_pid = str_sub(Tumor_Sample_Barcode, 1, 9),
      dmp_tumor = str_sub(Tumor_Sample_Barcode, 1, 13),
      dmp_sample = str_sub(Tumor_Sample_Barcode, 1, 17),
      .before = 1
    ) %>%
    filter(dmp_sample %in% dmp_ids) %>%
    janitor::clean_names()

  oncokb_summary_table <- oncokb_clindat %>%
    # resistance mutations are seperate
    mutate(mutations_full = ifelse(resistance_mutations != "", paste(oncogenic_mutations, resistance_mutations, sep = ";"), oncogenic_mutations)) %>%
    tidyr::separate_rows(mutations_full, sep = ";") %>%
    tidyr::separate(mutations_full, into = c("gene", "variant"), sep = " ", remove = F) %>%
    mutate(variant_type = case_when(
      variant %in% c("Amplification", "Deletion") ~ "CNA",
      variant == "fusion" ~ "Fusion",
      grepl("^p\\.", variant) ~ "SSM",
      .default = "Other"
    )) %>%
    dplyr::relocate(gene, variant, variant_type, .after = patient_id)

  oncokb_summary_table$sample_id <- NULL

  oncokb_obj <- list(oncokb_summary_table = oncokb_summary_table, oncokb_clindat = oncokb_clindat, oncokb_muts = oncokb_muts, oncokb_cna = oncokb_cna, oncokb_sv = oncokb_sv)

  return(oncokb_obj)
}
