#' @export
split_integrate <- function(srt,
                            modality = c("RNA", "CITE"),
                            split_var,
                            sketch = FALSE,
                            ncells = 3000,
                            method = c("harmony", "rpca")) {
  # require("rlang")
  modality <- match.arg(modality)

  method <- match.arg(method)


  spl_var <- srt@meta.data[, split_var] # Weird workaround -- need named vector
  names(spl_var) <- colnames(srt)


  DefaultAssay(srt) <- modality
  plan("sequential")
  srt[[modality]] <- split(srt[[modality]], f = spl_var)

  if (method == "harmony") {
    int_method <- HarmonyIntegration
  } else if (method == "rpca") {
    int_method <- RPCAIntegration
  }

  if (modality == "RNA") {
    srt <- FindVariableFeatures(srt, verbose = TRUE)
  } else {
    VariableFeatures(srt) <- rownames(srt)
  }

  if (sketch == T) {
    sketch_assay_name <- glue::glue("sketch_{modality}")
    sketch_pca_name <- glue::glue("sketch_pca_{modality}")
    sketch_integrated_name <- glue::glue("sketch_integrated_{modality}")
    sketch_varname <- glue::glue("{modality}_leverage.score")
    sketch_umap_name <- glue::glue("{modality}_sketch_umap")

    srt <- SketchData(
      object = srt,
      ncells = ncells,
      method = "LeverageScore",
      sketched.assay = sketch_assay_name,
      var.name = sketch_varname
    )

    DefaultAssay(srt) <- sketch_assay_name

    if (modality == "RNA") {
      srt <- FindVariableFeatures(srt, verbose = TRUE)
    } else {
      VariableFeatures(srt) <- rownames(srt)
    }

    srt <- ScaleData(srt, verbose = T)
    srt <- RunPCA(srt, verbose = T, reduction.name = sketch_pca_name)

    srt <- IntegrateLayers(
      object = srt,
      method = int_method,
      assay = sketch_assay_name,
      orig.reduction = sketch_pca_name,
      new.reduction = sketch_integrated_name,
      verbose = TRUE
    )

    # Harmony doesn't correctly set new reduction name
    if (method == "harmony") {
      srt[[sketch_integrated_name]] <- srt[["harmony"]]
    }

    srt[[sketch_assay_name]] <- JoinLayers(srt[[sketch_assay_name]])

    # Error here with RNA data
    srt <- RunUMAP(srt,
      assay = sketch_assay_name,
      reduction = sketch_integrated_name,
      reduction.name = sketch_umap_name,
      dims = 1:30,
      return.model = T,
      verbose = T
    )

    # Reproject full dataset
    # resplit the sketched cell assay into layers this is required to project the integration onto
    # all cells
    spl_var <- srt@meta.data[, split_var] # Weird workaround -- need named vector
    names(spl_var) <- colnames(srt)
    srt[[sketch_assay_name]] <- split(srt[[sketch_assay_name]], f = spl_var)


    full_pca_name <- glue::glue("{modality}_{method}_integrated")
    full_umap_name <- glue::glue("{modality}_umap_{method}_integrated")

    srt <- ProjectIntegration(
      object = srt,
      sketched.assay = sketch_assay_name,
      assay = modality,
      reduction = sketch_integrated_name,
      reduction.name = full_pca_name
    )
    srt <- JoinLayers(srt)
  } else {
    pca_first_name <- glue::glue("{modality}_pca")
    full_pca_name <- glue::glue("{modality}_{method}_integrated")
    full_umap_name <- glue::glue("{modality}_umap_{method}_integrated")


    srt <- NormalizeData(srt)
    srt <- ScaleData(srt, verbose = T)
    srt <- RunPCA(srt, verbose = T, reduction.name = pca_first_name)

    srt <- IntegrateLayers(
      object = srt,
      method = int_method,
      assay = modality,
      orig.reduction = pca_first_name,
      new.reduction = full_pca_name,
      verbose = TRUE
    )

    # Harmony doesn't correctly set new reduction name
    if (method == "harmony") {
      srt[[sketch_integrated_name]] <- srt[["harmony"]]
    }
  }

  # Now we can remove the sketch assay
  # now that we have projected the full dataset, switch back to analyzing all cells
  DefaultAssay(srt) <- modality
  srt <- JoinLayers(srt)
  srt <- RunUMAP(srt, reduction = full_pca_name, dims = 1:30, reduction.name = full_umap_name)

  srt <- FindNeighbors(srt, reduction = full_pca_name, dims = 1:30)
  srt <- FindClusters(srt, resolution = 0.25, cluster.name = glue("{modality}_clusters"))

  return(srt)
}
