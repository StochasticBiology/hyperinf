#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # attach ggraph quietly so custom guides are registered
  if (!"ggraph" %in% loadedNamespaces()) {
    attachNamespace("ggraph")
  }
}

#' @import ggraph
NULL