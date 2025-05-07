#' @export
load_markers <- function() {
  markers <- list()
  # Based on HLCA (https://www.nature.com/articles/s41591-023-02327-2)
  markers$epi_markers_hlca <- list(
    "Epithelial" = c("FXYD3", "EPCAM", "ELF3"),
    "Alveolar Epithelium" = c("HOPX", "SFTA2", "SLC39A8"),
    "AT1" = c("CLIC5", "SPOCK2", "TIMP3"),
    "AT2" = c("MFSD2A", "TCIM", "C11orf96"),
    "AT2 proliferating" = c("CENPW", "CDKN3", "BIRC5"),
    "Airway Epithelium" = c("IGFBP2", "SERPINF1", "TSPAN1"),
    "Basal Resting" = c("KRT15", "KRT17", "TP63"),
    "Suprabasal" = c("KRT5", "SERPINB3"),
    "Hillock-like" = c("DSG3", "RAB38", "TNS4"), # KRT6A, KRT13, KRT14
    "Deuterosomal" = c("CDC20B", "KDELC2", "E2F7"),
    "Multiciliated (nasal)" = c("BEST4", "LYPD2"),
    "Multiciliated (non-nasal)" = c("RSPH1", "C20orf85", "C9orf24"),
    "Club (non-nasal)" = c("CYP2F1", "SCGB3A1", "BPIFB1"),
    "Club (nasal)" = c("SERPINB3", "TCN1", "ASRGL1"),
    "Goblet (nasal)" = c("CEACAM5", "LYPD2", "PSCA"),
    "Goblet (bronchial)" = c("GPX8", "ALDH1A3", "CEACAM5"),
    "Goblet (subsegmental)" = c("BPIFB1", "C16orf89", "NPDC1"),
    "AT0" = c("SFTPB", "SCGB3A2", "SFTA2"),
    "pre-TB secretory" = c("SFTPB", "RNASE1", "SFTA1P"),
    "Ionocyte" = c("BSND", "IGF1", "CLCNKB"),
    "Tuft" = c("DAB1"),
    "Neuroendocrine" = c("CELF3", "SLC6A17", "CDK5R2"),
    "Submucosal Gland" = c("TCN1", "FAM3D", "CCL28"),
    "SMG serous (nasal)" = c("STATH", "ODAM", "MUC7"),
    "SMG serous (bronchial)" = c("RP11-1143G9.4", "S100A1", "PRB3"),
    "SMG mucous" = c("BPIFB2", "MUC5B", "NKX3-1"),
    "SMG duct" = c("MGST1", "KRT19", "KLF5")
  )

  # Custom for main figure
  markers$epi_markers_final <- list(
    "Epithelial" = c("EPCAM", "ELF3", "FXYD3"),
    "Alveolar\nEpithelium" = c("NAPSA", "HOPX", "SLC39A8", "NKX2-1"),
    "AT1" = c("CLIC5", "TIMP3", "SPOCK2", "AGER"),
    "AT2" = c("SFTPB", "SFTPC", "SFTPD", "ETV5"), # MFSD2A also in hepatocytes
    # "Airway\nEpithelium" = c("IGFBP2", "SERPINF1", "TSPAN1"),
    "Basal" = c("KRT15", "KRT17", "TP63", "KRT5"),
    "Hillock" = c("KRT6A", "KRT13", "KRT14", "DSG3"), # Serpinb2 hillock luminal
    "Secretory" = c("SCGB1A1", "SCGB3A1", "SCGB3A2"),
    "AT0" = c("SFTPC", "SCGB3A2", "SFTA2"),
    "Multiciliated" = c("FOXJ1", "RSPH1"),
    "Neuroendocrine" = c("ASCL1", "NCAM1", "NEUROD1", "SYP"),
    "Cycling" = c("MKI67", "TOP2A", "POLA2", "CENPF"),
    "Hepatocyte" = c("ALB", "HP", "APOB") # ,
    # "PD9/PD27" = c("KDELC2", "ASRGL1", "CEACAM5", "DAB1", "MUC5B")
  )

  markers$cell_type_markers <- list(
    "Epithelial" =
      c(
        "EGFR",
        # "FXYD3",
        "EPCAM",
        "ELF3",
        "AGER",
        "NAPSA",
        "LAMP3",
        "SFTPA2"
      ),
    "T-Cell" =
      c(
        "CD3E",
        "CD3D",
        "CD247",
        # "FOXP3",
        # "CD8A",
        "IL7R"
      ),
    "B-Cell" =
      c(
        "MS4A1",
        "BANK1",
        "CD79A"
        # "CD19",
        # "CR2"
      ),
    "Plasma cell" =
      c(
        "JCHAIN",
        "IGKC",
        "PIM2"
      ),
    "Myeloid" =
      c(
        # "FCER1G",
        "CLEC7A",
        # "MARCO",
        # "CD68",
        "MRC1",
        "MSR1"
        # "FCGR3A"
      ),
    "Mast Cell" =
      c(
        "KIT",
        "CPA3",
        # "TPSAB1",
        "TPSB2"
      ),
    "Endothelial" =
      c(
        "VWF",
        "FLT1",
        "PECAM1",
        "PTPRB"
      ),
    "Fibroblast" =
      c(
        "COL1A2",
        "COL3A1",
        "DCN",
        "CACNA1C"
      )
  )

  markers$custom_markers <- list(
    "Epithelial" = c("KRT7", "PIGR", "ELF3", "CYB5A", "KRT8", "KRT19", "TACSTD2", "MUC1", "S100A14", "CXCL17"),
    "Basal" = c("KRT5", "TP63", "S100A2", "KRT6A", "TNS4"),
    "AT1" = c("AGER", "RTKN2", "CLIC5"),
    "AT2" = c("SFTPC", "SFTPA1", "SFTPA2", "WIF1", "HHIP", "CA2", "ETV5"),
    "Ciliated" = c("GSTA1", "DTHD1", "PIFO", "FOXJ1", "CCDC78"),
    "Neuroendrocrine" = c("CHGA", "CALCA", "ASCL1", "CHGB", "GRP", "BEX1"),
    "B-Cell" = c("MS4A1", "PAX5", "CD79A"),
    "Plasma" = c("JCHAIN", "IGKC", "IGHA1", "IGHG1", "MZB1", "ISG20"),
    "T-Cell" = c("IL7R", "CD3E", "CD3D"),
    "Endothelial" = c("VWF", "FLT1"),
    "Lymphatic EC" = c("MMRN1", "CCL21"),
    "Mast cell" = c("KIT", "CPA3"),
    "Fibroblast" = c("COL1A1", "DCN"),
    "Pericyte" = c("PDGFRB", "CALD1"),
    "Cycling cells" = c("TOP2A", "CENPF", "CENPP", "ASPM"),
    "PDC" = c("IL3RA", "TCF4", "LTB", "GZMB", "ITM2C", "IRF8", "PLD4", "IRF7")
  )
  return(markers)
}
