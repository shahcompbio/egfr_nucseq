#' @export
featureScore <- function(x) {
  s <- apply(x, 1, function(g) {
    g <- abs(g) + 2.2e-16
    p_i <- g / sum(g)
    crossprod(p_i, log2(p_i))
  })
  # scale, translate and return the result
  s <- 1 + s / log2(ncol(x))
  return(s)
}

#' @export
assignFeature <- function(x, n_mads = 3) {
  s <- featureScore(x)

  # filter for the genes whose score is greater than \mu + 3 \sigma
  th <- median(s) + n_mads * mad(s)
  sel <- s >= th

  # print( s[sel] )
  # print(sum(sel))

  # build a matrix with:
  #-> row#1=max column index, row#2=max value in row, row#3=row index
  temp <- 0
  g.mx <- apply(x, 1, base::which.max)
  max_val <- map(names(g.mx), .f = function(gene) {
    idx_max <- g.mx[gene]
    val_max <- x[gene, idx_max]
  }) %>%
    list_c() %>%
    setNames(names(g.mx))


  # test the second criteria
  med <- median(abs(x))
  sel2 <- max_val >= med

  res <- data.frame(Gene = names(g.mx), feature_score = s, factor = g.mx, max_value = max_val, mad_sig = sel, med_sig = sel2)

  cli::cli_alert_info("{sum(sel)} significant genes: ")
  cli::cli_alert_info("{paste(names(s[sel]), collapse = ', ')}")

  # return result
  return(res)
}
