#' Plot PCA across fitted models
#'
#' @param fits A list of fitted hypercubic inference models (output from hyperinf)
#' @param expt.names Optional vector of labels for each element of the fit list
#'
#' @return A ggplot object
#' @export
plot_hyperinf_pca <- function(fits,
                         expt.names = NULL) {
  df_list = list()
  i = 1
  for(this.fit in fits) {
    fit.type = hyperinf_gettype(this.fit)
    if(is.null(fit.type)) {
      return(NULL)
    }
    
    tmp = hyperinf_bubbles(this.fit)
    df_list[[i]] = tmp
    i = i+1
  }

    # ---- 1. collect all (feature, order) keys ----
    all_keys <- dplyr::distinct(
      dplyr::bind_rows(df_list),
      feature,
      order
    )
    
    all_keys <- all_keys[order(all_keys$feature, all_keys$order), ]
    
    # ---- 2. convert each df into aligned vector ----
    mat_list <- lapply(df_list, function(df) {
      
      df2 <- df[, c("feature", "order", "prob")]
      
      # join with full key grid
      merged <- dplyr::right_join(df2, all_keys, by = c("feature", "order"))
      
      # replace NA with 0
      merged$prob[is.na(merged$prob)] <- 0
      
      # order consistently
      merged <- merged[order(merged$feature, merged$order), ]
      
      merged$prob
    })
    
    # ---- 3. stack into matrix ----
    mat <- base::do.call(base::rbind, mat_list)
    
    base::colnames(mat) <- paste0(
      "f", all_keys$feature, "_o", all_keys$order
    )
    
    # ---- 4. PCA ----
    col_var <- apply(mat, 2, stats::var)
    keep <- is.finite(col_var) & col_var > 0
    
    mat <- mat[, keep, drop = FALSE]
    
    pca <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
    
    ### plot
    var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
    pc1_lab <- paste0("PC1 (", round(100 * var_expl[1], 1), "%)")
    pc2_lab <- paste0("PC2 (", round(100 * var_expl[2], 1), "%)")
    scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
    names(scores) <- c("PC1", "PC2")
    
    # add labels if provided
    if (!is.null(expt.names)) {
      scores$label <- expt.names
    } else {
      scores$label <- base::seq_len(nrow(scores))
    }
    
    # ---- plot ----
    p <- ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::labs(x = pc1_lab, y = pc2_lab) +
      ggplot2::theme_minimal()
    
    if (!is.null(expt.names)) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = label),
        vjust = -0.6,
        size = 3
      )
    }
    
   return(p)
}
