#' Propose a "phylogenetically"-linked dataset clustering observations based on similarity
#'
#' @param mat A required matrix or data.frame containing binary observations
#' @param method Character (default "clustering"): method to cluster observations. Either "clustering" (hierarchical clustering with Manhattan distance) or "nj" (neihgbour joining)
#'
#' @return A named list with a (possibly reduced) dataset and a putative tree
#' @examples
#' data = matrix(c(0,0,1, 1,0,0, 1,0,0, 1,0,0), ncol=3, nrow=4)
#' binary_phylogeny(data)
#' @export
binary_phylogeny <- function(data, method = "clustering") {
  
  root_tip = NULL
  mat = clean_data(data)
  mat <- mat[!duplicated(mat), ]
  rownames(mat) = paste("row-", 1:nrow(mat), sep="")
  d <- dist(mat, method = "manhattan")
  
  if(method == "clustering") {

  hc <- hclust(d, method = "complete")
  
  tree <- ape::as.phylo(hc)
  
  if(!is.null(root_tip))
    tree <- ape::root(tree, outgroup = root_tip, resolve.root = TRUE)
  } else {
    tree <- ape::nj(d)
    tree <- ape::root(tree, outgroup = tree$tip.label[1], resolve.root = TRUE)
  }
  return(list(data=mat,
              tree=tree))
}

#' Compare inference outputs assuming cross-sectional or minimally convergent dynamics
#'
#' @param mat A required matrix or data.frame containing binary observations
#' @param method Character (default "clustering"): method to cluster observations. Either "clustering" (hierarchical clustering with Manhattan distance) or "nj" (neihgbour joining)
#' @param ... Other arguments to pass to hyperinf
#'
#' @return A ggplot object
#' @examples
#' data = matrix(c(0,0,1, 1,0,0, 1,0,0, 1,0,0), ncol=3, nrow=4)
#' hyperinf_span_relatedness(data)
#' @export
hyperinf_span_relatedness = function(mat, method="clustering", ...) {
  fit.cs = hyperinf(mat, ...)
  clustered = binary_phylogeny(mat, method=method)
  fit.phy = hyperinf(clustered$data, clustered$tree, ...)
  plot_hyperinf_comparative(list(fit.cs, fit.phy), style="full")
}

