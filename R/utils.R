#' @export
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#' @export
standardize <- function(x, trim = c(0, 1)) {
  if (is.numeric(trim)) {
    x[x < quantile(x, prob = trim[1])] <- quantile(x, prob = trim[1])
    x[x > quantile(x, prob = trim[2])] <- quantile(x, prob = trim[2])
  }

  return((x - min(x)) / (max(x) - min(x)))
}

#' @export
memoSort <- function(M) {
  geneOrder <- sort(rowSums(M), decreasing = TRUE, index.return = TRUE)$ix
  scoreCol <- function(x) {
    score <- 0
    for (i in 1:length(x)) {
      if (x[i]) {
        score <- score + 2^(length(x) - i)
      }
    }
    return(score)
  }
  scores <- apply(M[geneOrder, ], 2, scoreCol)
  sampleOrder <- sort(scores, decreasing = TRUE, index.return = TRUE)$ix
  return(M[geneOrder, sampleOrder])
}


#' @export
DotPlotter <- function(obj,
                       features,
                       group.by,
                       size_limits = c(0, 100),
                       size_range = c(0.25, 5),
                       cluster_groups = F,
                       cluster_feats = F,
                       rotate_y_strip = F) {
  # Enable multiplotting of same feature
  # First uniqify features
  features_uniq <- unique(unlist(features))

  marker_dat <- Seurat::DotPlot(obj, features = features_uniq, group.by = group.by, cluster.idents = F)$data

  if (class(features) == "list") {
    x <- tibble::enframe(features, name = "feature.groups", value = "features.plot") %>%
      tidyr::unnest(cols = "features.plot") %>%
      dplyr::mutate(feature.groups.set = paste(feature.groups, features.plot, sep = "__"))
    marker_dat <- full_join(marker_dat, x, relationship = "many-to-many", by = "features.plot")
    marker_dat$feature.groups <- factor(marker_dat$feature.groups, levels = names(features))
    marker_dat$feature.groups.set <- factor(marker_dat$feature.groups.set, levels = unique(x$feature.groups.set))
    marker_dat <- marker_dat[!is.na(marker_dat$id), ]
    # Reorder?
    # marker_dat <- marker_dat %>%
    #   arrange(dplyr::desc(features.plot), dplyr::desc(feature.groups.set))
  }

  if (is.null(marker_dat[["feature.groups"]])) {
    marker_dat[["feature.groups"]] <- "Markers"
    marker_dat[["feature.groups.set"]] <- paste("Markers", marker_dat$features.plot, sep = "__")
  }



  marker_plot <- marker_dat %>%
    ggplot(aes(x = feature.groups.set, y = fct_rev(id))) +
    geom_point(pch = 21, aes(fill = avg.exp.scaled, size = pct.exp), color = "black", stroke = 0.5) +
    scale_size(limits = size_limits, range = size_range) +
    scale_x_discrete(breaks = marker_dat$feature.groups.set, labels = marker_dat$features.plot) +
    facet_grid(. ~ feature.groups, scales = "free", space = "free") +
    scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
    theme(
      axis.line = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      strip.background.x = element_blank(),
      strip.clip = "off"
    ) +
    guides(
      x = guide_axis(angle = 90),
      fill = guide_colorbar(keywidth = 0.4, keyheight = 3, order = 1),
      size = guide_legend(keyheight = 1, override.aes = list(pch = 19), order = 2)
    ) +
    labs(x = "Feature", y = "Cell Type", fill = "Scaled\nExpression", size = "Percent\nExpression")

  if (cluster_groups) {
    require(legendry)
    mat <- marker_dat %>%
      pivot_wider(id_cols = id, values_from = "avg.exp.scaled", names_from = "feature.groups.set") %>%
      column_to_rownames("id") %>%
      as.matrix()
    hc <- hclust(dist(mat), method = "ward.D2")
    marker_plot <- marker_plot +
      legendry::scale_y_dendro(clust = hc)
  }

  if (cluster_feats) {
    require(legendry)
    mat <- marker_dat %>%
      pivot_wider(id_cols = id, values_from = "avg.exp.scaled", names_from = "feature.groups.set") %>%
      column_to_rownames("id") %>%
      as.matrix() %>%
      t()
    hc2 <- hclust(dist(mat), method = "ward.D2")
    marker_plot <- marker_plot +
      legendry::scale_x_dendro(clust = hc2)
  }

  if (rotate_y_strip) {
    marker_plot <- marker_plot + theme(strip.text.x = element_text(angle = 90, hjust = 0))
  }
  return(marker_plot)
}


#' @export
bulk_exp <- function(obj, cell_type_col, sample_id_col, meta_cols = NULL, features, layer = "data", min_cells = 20) {
  df <- FetchData(
    object = obj,
    vars = c(sample_id_col, cell_type_col, features),
    layer = layer
  )

  df_test <- df %>%
    add_count(.data[[sample_id_col]], .data[[cell_type_col]], name = "n_type") %>%
    dplyr::filter(n_type >= min_cells) %>%
    # group_by(cell_type_l4, n_type_tp, across(all_of(meta_cols))) %>%
    group_by(.data[[sample_id_col]], .data[[cell_type_col]], n_type) %>%
    summarize(across(all_of(features), \(x) mean(x, na.rm = T))) %>%
    ungroup()

  if (!is.null(meta_cols)) {
    meta_df <- distinct(obj@meta.data[, c(sample_id_col, meta_cols)])
    df_test <- left_join(df_test, meta_df)
  }

  df_test$cell_type <- df_test[[cell_type_col]]

  return(df_test)
}


#' @export
umap_theme <- function(axis.title = element_text(hjust = 0), ...) {
  theme_add <- ggplot2::theme(
    axis.line = element_line(arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")), linewidth = 1),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = axis.title,
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_blank(),
    validate = T,
    ...
  )
  return(theme_add)
}

#' @export
label_num_generator <- function(x) {
  y <- table(x)
  y <- y[y > 0]
  paste(names(y), glue("({prettyNum(y, big.mark = ',')})"))
}


#' @export
umap_guides <- function(arrow_size = 3, units = "cm", ...) {
  umap_axis <- ggh4x::guide_axis_truncated(
    trunc_lower = unit(0, "npc"),
    trunc_upper = unit(arrow_size, units)
  )
  return(ggplot2::guides(
    x = umap_axis,
    y = umap_axis
  ))
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
