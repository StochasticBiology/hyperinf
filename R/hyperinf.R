#' Curate data into a systematic format for inference
#'
#' Data can be flexibly supplied as a matrix or data frame, with or without a first labelling column, and with binary observations 0, 1 (with ? or 2 for uncertainty) as individual columns or concatenated strings
#'
#' @param data A required matrix or data.frame containing binary observations
#'
#' @return A matrix of observations in a standardised format
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' mat = clean_data(data)
#' @export
clean_data = function(data) {
  # Check input type
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("`data` must be a matrix or data.frame")
  }
  if(nrow(data) < 2 | ncol(data) < 2) {
    stop("There's not enough data here; nothing to do!")
  }
  
  x = data[1,ncol(data)]
  if(is.character(x)) {
    char.data = TRUE
  } else {
    char.data = FALSE
  }
  
  # If it's a data.frame, inspect the first column
  if (is.data.frame(data)) {
    first_col <- data[[1]]
    
    # Check if any entry is not 0, 1, or "?"
    # Coerce factors to character
    first_col_char <- as.character(unlist(strsplit(first_col, "")))
    invalid_entries <- !first_col_char %in% c("0", "1", "?", NA)
    
    if (any(invalid_entries)) {
      # Use columns 2 onward
      mat <- as.matrix(data[, -1, drop = FALSE])
    } else {
      # Use entire dataframe
      mat <- as.matrix(data)
    }
  } else {
    # Already a matrix
    mat <- data
  }
  
  if(char.data) {
    tmp = strsplit(mat, "")
    mat = matrix(unlist(tmp), byrow=TRUE, ncol=length(tmp[[1]]))
  }
  
  # Optionally coerce to numeric if needed
  # Convert "?" to 2
  mat[mat == "?"] <- 2
  mat[is.na(mat)] <- 2
  mat <- apply(mat, c(1,2), function(x) as.numeric(x))
  
  rownames(mat) = colnames(mat) = NULL
  
  return(mat)
}

#' Do hypercubic inference on some data
#'
#' Data can be flexibly supplied as a matrix or data frame, with or without an accompanying phylogeny
#' describing how observations are related
#'
#' Method can be specified: "hypermk", "hyperhmm", "hypertraps", "pli", or "hyperdags". If unspecified, the
#' most detailed approach compatible with the data will be chosen
#'
#' @param data A required matrix or data.frame containing binary observations
#' @param tree Optional tree object
#' @param losses Boolean (default FALSE) whether to consider losses rather than gains of features
#' @param method A character string, either empty (default) to allow automatic choice of algorithm, or one of the options in the text above
#' @param reversible Boolean (default FALSE) whether to allow reversible transitions
#' @param auto.cluster Boolean (default FALSE) whether to cluster the dataset by similarity to estimate relatedness
#' @param auto.cluster.method Character (default "clustering") method to use for clustering; see binary_phylogeny
#' @param boot.parallel Integer (default 0) number of bootstrap resamples to run in parallel (only meaningful for HyperHMM and HyperLAU)
#' @param ... other options to pass to the inference method used. For example, nboot for HyperHMM, length/kernel/walkers for HyperTraPS.
#'
#' @return A fitted hypercubic inference object
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit = hyperinf(data)
#' @export
hyperinf <- function(data,
                     tree = NULL,
                     losses = FALSE,
                     method = "",
                     reversible = FALSE,
                     auto.cluster = FALSE,
                     auto.cluster.method = "clustering",
                     boot.parallel = 0,
                     ...) {
  
  
  # TO DO 
  #       -- TIMINGS
  #       -- NOT SURE UNCERTAINTY'S HANDLED RIGHT
  
  mat = clean_data(data)
  
  if(auto.cluster == TRUE) {
    if(!is.null(tree)) {
      message("You provided a tree but also asked me to auto-cluster the data. I'm doing the latter and ignoring the tree.")
    }
    clustered = binary_phylogeny(mat, method=auto.cluster.method)
    data = clustered$data
    mat = clean_data(data)
    tree = clustered$tree
  }
  if(is.matrix(data) & !is.null(tree)) {
    if(is.null(rownames(data))) {
      message("You've given me a matrix without rownames and a tree, but to use tree-linked data I need a dataframe with IDs in the first column!")
      tree = NULL
      Sys.sleep(3)
    }
  }
  if(losses == TRUE) {
    mat = 1-mat
  }
  L = ncol(mat)
  feature.names = colnames(data[,(ncol(data)-L+1):ncol(data)])
  if(length(feature.names) != L) {
    feature.names = NULL
  }
  uncertainty = FALSE
  if(any(is.na(mat)) | any(mat == 2) | any(mat == -1)) {
    uncertainty = TRUE
  }
  
  if(!is.null(tree)) {
    if(is.matrix(data)) {
      df = cbind(data.frame(ID = rownames(data)), as.data.frame(mat))
    } else {
      df = cbind(data.frame(ID = data[,1]), as.data.frame(mat))
      colnames(df)[1] = "ID"
    }
    if(!setequal(df$ID, tree$tip.label)) {
      to.keep = which(tree$tip.label %in% df$ID)
      tree = ape::keep.tip(tree, to.keep)
      if(!setequal(df$ID, tree$tip.label)) {
        message("Warning: data IDs don't match tree tip labels. This won't work!")
        return("Error")
      }
      message("More tree tips than observations... dropping those without records.")
    }
    if(any(mat == 2) | any(mat == -1)) {
      it = TRUE
      dots <- list(...)
      if ("independent.transitions" %in% names(dots)) {
        it = dots$independent.transitions
      }
      df[df==2 | df == -1] = "?"
      c.tree = hyperlau::curate.uncertain.tree(tree, df, independent.transitions = it)
      if(method == "") {
        if(L < 10) {
          method = "hyperlau"
          message("Selecting HyperLAU...")
        } else {
          method = "pli"
          message("Selecting PLI...")
        }
      }
    } else {
      c.tree = hypertrapsct::curate.tree(tree, df)
    }
    if(nrow(c.tree$dests) < 2) {
      message("Warning: less than 2 independent transitions found! Results may not be what you expect.")
    }
  }
  
  if(method == "") {
    if(L <= 7 & reversible == TRUE) {
      method = "hypermk"
      message("Selecting HyperMk...")
    } else if(uncertainty == TRUE) {
      if(L <= 8) {
        method = "hyperlau"
        message("Selecting HyperLAU...")
      } else {
        method = "hypertraps"
        message("Selecting HyperTraPS...")
      }
    } else if(L <= 12) {
      method = "hyperhmm"
      message("Selecting HyperHMM...")
    } else if(L < 40) {
      method = "hypertraps"
      message("Selecting HyperTraPS...")
    } else {
      method = "hyperdags"
      message("Selecting HyperDAGs...")
    }
  } else {
    if(method == "hypermk" & L > 6) {
      message("L > 6 will be hard and unstable for HyperMk! Pausing in case you want to break...")
      Sys.sleep(3)
    }
    if(reversible == TRUE & method != "hypermk") {
      message("Only HyperMk can deal with reversibility. I'm turning off reversibility.")
      reversible = FALSE
    }
    if(method == "hyperhmm" & L > 18) {
      message("L > 18 for HyperHMM is untested and may cause memory errors. Consider HyperTraPS. Pausing in case you want to break...")
      Sys.sleep(3)
    }
    if(!(method %in% c("hypermk", "hyperhmm", "hyperlau", "hypertraps", "hyperdags", "pli"))) {
      message("I didn't recognise that method. Switching to HyperDAGs. Pausing in case you want to break...")
      Sys.sleep(3)
      method = "hyperdags"
    }
  }
  
  if(!is.null(tree)) {
    if(any(c.tree$srcs == 2) & !method %in% c("pli", "hyperlau")) {
      message("Only HyperLAU and phenotype landscape inference can deal with uncertain ancestors.")
      if(L < 8) {
        message("Switching to HyperLAU. Pausing in case you want to break...")
        method = "hyperlau"
      } else {
        message("Switching to PLI. Pausing in case you want to break...")
        method = "pli"
      }
      Sys.sleep(3)
    }
    if(!any(c.tree$srcs == 2) & any(c.tree$dests == 2) & !(method %in% c("pli", "hyperlau", "hypertraps"))) {
      message("Only HyperTraPS, HyperLAU, and PLI can deal with uncertain observations. Switching to HyperTraPS. Pausing in case you want to break...")
      Sys.sleep(3)
      method = "hypertraps"
    }
  } else {
    if(any(mat == 2) & !(method %in% c("pli", "hyperlau", "hypertraps"))) {
      message("Only HyperTraPS, HyperLAU, and PLI can deal with uncertain observations. Switching to HyperTraPS. Pausing in case you want to break...")
      Sys.sleep(3)
      method = "hypertraps"
    }
  }
  
  if(boot.parallel > 0 & !(method %in% c("hyperlau", "hyperhmm"))) {
    message("Warning: parallel bootstrapping only supported for HyperHMM and HyperLAU!")
  }
  
  if(method == "hypermk") {
    if (!is.null(tree)) {
      fit = hypermk::mk_infer_phylogenetic(mat, tree, reversible = reversible)
      this.data = list(obs=mat, tree=tree)
    } else {
      fit = hypermk::mk_infer_cross_sectional(mat, reversible = reversible)
      this.data = list(obs=mat)
    }
  } else if(method == "hyperhmm") {
    dots <- list(...)
    if (!"nboot" %in% names(dots)) {
      dots$nboot <- 0
    }
    if(!is.null(tree)) {
      identicals = which(apply(c.tree$dests == c.tree$srcs, 1, all))
      dests = c.tree$dests[-identicals,]
      srcs = c.tree$srcs[-identicals,]
      if(boot.parallel > 0) {
        dots$boot.parallel <- NULL
        b.dests = b.srcs = list()
        for(i in 1:(boot.parallel+1)) {
          if(i == 1) {
            refs = 1:nrow(dests)
          } else {
            refs = sample(1:nrow(dests), nrow(dests), replace=TRUE)
          }
          b.dests[[i]] = dests[refs,]
          b.srcs[[i]] = srcs[refs,]
        }
        fits = parallel::mclapply(1:length(b.dests),
                                  function(i) {
                                    do.call(hyperhmm::HyperHMM, c(list(obs = b.dests[[i]], 
                                                                       initialstates = b.srcs[[i]]),
                                                                  dots))
                                  }, mc.cores = parallel::detectCores())
        this.data = list(obs = b.dests[[1]], initialstates = b.srcs[[1]])
        fit = fits[[1]]
        fit$transitions$p.boot = 1
        for(i in 2:length(fits)) {
          tmp = fits[[i]]$transitions
          tmp$p.boot = i
          fit$transitions = rbind(fit$transitions, tmp)
        }
        fit$boots = fits
      } else {
        fit = do.call(hyperhmm::HyperHMM, c(list(obs = dests, initialstates = srcs), dots))
        this.data = list(obs = dests, initialstates = srcs)
      }
    } else {
      if(boot.parallel > 0) {
        dots$boot.parallel <- NULL
        b.mat = list()
        for(i in 1:(boot.parallel+1)) {
          if(i == 1) {
            refs = 1:nrow(mat)
          } else {
            refs = sample(1:nrow(mat), nrow(mat), replace=TRUE)
          }
          b.mat[[i]] = mat[refs,]
        }
        fits = parallel::mclapply(1:length(b.mat),
                                  function(i) {
                                    do.call(hyperhmm::HyperHMM, c(list(obs = b.mat[[i]]), 
                                                                  dots))
                                  }, mc.cores = parallel::detectCores())
        this.data = list(obs = b.mat[[1]])
        fit = fits[[1]]
        fit$transitions$p.boot = 1
        for(i in 2:length(fits)) {
          tmp = fits[[i]]$transitions
          tmp$p.boot = i
          fit$transitions = rbind(fit$transitions, tmp)
        }
        fit$boots = fits
      } else {
        fit = do.call(hyperhmm::HyperHMM, c(list(obs = mat), dots))
        this.data = list(obs = mat)
      }
    }
  } else if(method == "hypertraps" | method == "pli") {
    pli = 0
    dots <- list(...)
    if (!"length" %in% names(dots)) {
      dots$length <- 4
    }
    if(method == "pli") {
      pli = 1
    }
    if(!is.null(tree)) {
      fit = do.call(hypertrapsct::HyperTraPS, c(list(obs = c.tree$dests, initialstates = c.tree$srcs, pli=pli), dots))
      this.data = list(obs = c.tree$dests, initialstates = c.tree$srcs)
    } else {
      fit = do.call(hypertrapsct::HyperTraPS, c(list(obs = mat, pli=pli), dots))
      this.data = list(obs = mat)
    }
  } else if(method == "hyperlau") {
    dots <- list(...)
    if("model" %in% names(dots)) {
      this.model = dots[["model"]]
    } else {
      this.model = -1
    }
    dots = dots[names(dots) != "independent.transitions"]
    if(!is.null(tree)) {
      if(boot.parallel > 0) {
        dots$boot.parallel <- NULL
        b.dests = b.srcs = list()
        for(i in 1:(boot.parallel+1)) {
          if(i == 1) {
            refs = 1:nrow(c.tree$dests)
          } else {
            refs = sample(1:nrow(c.tree$dests), nrow(c.tree$dests), replace=TRUE)
          }
          b.dests[[i]] = c.tree$dests[refs,]
          b.srcs[[i]] = c.tree$srcs[refs,]
        }
        fits = parallel::mclapply(1:length(b.dests),
                                  function(i) {
                                    do.call(hyperlau::HyperLAU, c(list(obs = b.dests[[i]], 
                                                                       initialstates = b.srcs[[i]]),
                                                                  dots))
                                  }, mc.cores = parallel::detectCores())
        this.data = list(obs = b.dests[[1]], initialstates = b.srcs[[1]])
        fit = fits[[1]]
        fit$model = this.model
        fit$Dynamics$p.boot = 1
        for(i in 2:length(fits)) {
          tmp = fits[[i]]$Dynamics
          tmp$p.boot = i
          fit$Dynamics = rbind(fit$Dynamics, tmp)
        }
        fit$boots = fits
      } else {
        fit = do.call(hyperlau::HyperLAU, c(list(obs = c.tree$dests, initialstates = c.tree$srcs), dots))
        fit$model = this.model
        this.data = list(obs = c.tree$dests, initialstates = c.tree$srcs)
      }
    } else {
      if(boot.parallel > 0) {
        dots$boot.parallel <- NULL
        b.mat = list()
        for(i in 1:(boot.parallel+1)) {
          if(i == 1) {
            refs = 1:nrow(mat)
          } else {
            refs = sample(1:nrow(mat), nrow(mat), replace=TRUE)
          }
          b.mat[[i]] = mat[refs,]
        }
        fits = parallel::mclapply(1:length(b.mat),
                                  function(i) {
                                    do.call(hyperlau::HyperLAU, c(list(obs = b.mat[[i]]), 
                                                                  dots))
                                  }, mc.cores = parallel::detectCores())
        this.data = list(obs = b.mat[[1]])
        fit = fits[[1]]
        fit$Dynamics$p.boot = 1
        for(i in 2:length(fits)) {
          tmp = fits[[i]]$Dynamics
          tmp$p.boot = i
          fit$Dynamics = rbind(fit$Dynamics, tmp)
        }
        fit$boots = fits
      } else {
        fit = do.call(hyperlau::HyperLAU, c(list(obs = mat), dots))
        this.data = list(obs = mat)
      }
    }
  }
  else if(method == "hyperdags") {
    if(!is.null(tree)) {
      srcs = apply(c.tree$srcs, 1, paste0, collapse="")
      dests = apply(c.tree$dests, 1, paste0, collapse="")
      fit = hyperdags::simplest_DAG(srcs, dests)
      this.data = list(obs = c.tree$dests, initialstates = c.tree$srcs)
    } else {
      srcs = apply(matrix(0, nrow=nrow(mat), ncol=ncol(mat)), 1, paste0, collapse="")
      dests = apply(mat, 1, paste0, collapse="")
      fit = hyperdags::simplest_DAG(srcs, dests)
      this.data = list(obs = mat)
    }
    fit$L = ncol(mat)
  }
  fit$feature.names = feature.names
  fit$data = this.data
  
  return(fit)
}

# pulled from HyperTraPS-CT source code
# compute a flux dataframe given a set of transitions
# (PosteriorAnalysis provides the latter but not the former)
compute_hypertraps_fluxes = function(my.post,
                                     truncate = -1,
                                     max.samps = Inf,
                                     use.timediffs = FALSE,
                                     featurenames = FALSE,
                                     no.times = TRUE) {
  edge.from = edge.to = edge.time = edge.change = c()
  bigL = my.post$L
  if(truncate == -1 | truncate > bigL) { truncate = bigL }
  nsamps = min(max.samps, nrow(my.post$routes))
  for(i in 1:nsamps) {
    state = paste0(rep("0", bigL), collapse = "")
    for(j in 1:truncate) {
      edge.from = c(edge.from, state)
      locus = my.post$routes[i,j]+1
      substr(state, locus, locus) <- "1"
      edge.to = c(edge.to, state)
      edge.change = c(edge.change, my.post$routes[i,j])
      if(use.timediffs == TRUE) {
        edge.time = c(edge.time, my.post$timediffs[i,j])
      } else {
        edge.time = c(edge.time, my.post$times[i,j])
      }
    }
  }
  
  df = data.frame(From=sapply(edge.from, BinToDecS),
                  To=sapply(edge.to, BinToDecS),
                  Change=edge.change, 
                  Time=edge.time)
  dfu = unique(df[,1:3])
  if(length(featurenames) > 1) {
    dfu$Change = featurenames[dfu$Change+1]
  }
  dfu$Flux = dfu$MeanT = dfu$SDT = NA
  for(i in 1:nrow(dfu)) {
    this.set = which(df$From==dfu$From[i] & df$To==dfu$To[i])
    dfu$Flux[i] = length(this.set)
    dfu$MeanT[i] = mean(df$Time[this.set])
    if(length(this.set) > 1) {
      dfu$SDT[i] = sd(df$Time[this.set])
    }
    if(no.times == TRUE) {
      dfu$label[i] = paste(c("+", dfu$Change[i]), collapse="")
      dfu$tlabel[i] = paste(c(signif(dfu$MeanT[i], digits=2), " +- ", signif(dfu$SDT[i], digits=2)), collapse="") 
    } else {
      dfu$label[i] = paste(c("+", dfu$Change[i], ":\n", signif(dfu$MeanT[i], digits=2), " +- ", signif(dfu$SDT[i], digits=2)), collapse="") 
    }
    
  }
  dfu$Flux = dfu$Flux / nsamps
  return(dfu) 
}

#' Plot a fitted hypercubic inference model
#'
#' @param fit A fitted hypercubic inference model (output from hyperinf)
#' @param plot.type A character string, either empty (default) to allow standardised plot, or "native" to produce plots from the source algorithm
#' @param threshold Double (default 0.05), probability flux below which edges will not be plotted
#' @param uncertainty Boolean, whether to visualise uncertainty over bootstraps (only for HyperLAU and HyperHMM)
#' @param feature.names Boolean or character vector (default TRUE). If TRUE, use feature names from fit. If FALSE, use numerical labels. If vector of length L, use those labels.
#'
#' @return A ggplot object
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit = hyperinf(data)
#' plot_hyperinf(fit)
#' @export
plot_hyperinf = function(fit,
                         plot.type = "",
                         threshold = 0.05,
                         uncertainty = "",
                         feature.names = TRUE) {
  if("best.graph" %in% names(fit)) {
    fit.type = "DAG"
  } else if("raw.graph" %in% names(fit)) {
    fit.type = "arborescence"
  } else if("posterior.samples" %in% names(fit)) {
    fit.type = "hypertraps"
  } else if("Dynamics" %in% names(fit)) {
    fit.type = "hyperlau"
  } else if("viz" %in% names(fit)) {
    fit.type = "hyperhmm"
  } else if("fitted_mk" %in% names(fit)) {
    fit.type = "mk"
  } else {
    message("Didn't recognise this model")
    return(ggplot2::ggplot())
  }
  
  if(length(feature.names) == 1) {
    if(feature.names == TRUE) { 
      feature.names = fit$feature.names 
    } else {
      feature.names = NULL
    } 
  } else if(length(feature.names) != fit$L) {
      feature.names = NULL
    }
  
  if(uncertainty != FALSE) {
    if(fit.type == "hyperlau") {
      if(length(unique(fit$Dynamics$Bootstrap)) > 1 | 
         length(unique(fit$Dynamics$p.boot)) > 1) {
        uncertainty = TRUE
      }
    } else if(fit.type == "hyperhmm") {
      if(length(unique(fit$transitions$Bootstrap)) > 1 |
         length(unique(fit$transitions$p.boot)) > 1) {
        uncertainty = TRUE
      }
    } 
  }
  if(uncertainty != TRUE) {
    uncertainty = FALSE
  }
  
  reversible = FALSE
  if(plot.type == "native") {
    if(fit.type == "mk") {
      out.plot = hypermk::mk.inference.plot(fit)
    } else if(fit.type == "DAG") {
      out.plot = hyperdags::plot_stage_gen(fit$best.graph, label.size = 3)
    } else if(fit.type == "arborescence") {
      out.plot = hyperdags::plot_stage_gen(fit$rewired.graph, label.size = 3)
    } else if(fit.type == "hyperhmm") {
      out.plot = hyperhmm::plot_standard(fit)
    } else if(fit.type == "hypertraps") {
      out.plot = hypertrapsct::plotHypercube.sampledgraph2(fit, node.labels = FALSE,
                                                           no.times = TRUE, edge.label.size = 3)
    }
  } else {
    tmp1 = get_plot_graph(fit, fit.type, uncertainty, reversible, threshold, feature.names)
    plot.graph = tmp1[["plot.graph"]]
    layers = tmp1[["layers"]]
    fluxes = tmp1[["fluxes"]]
  }
  
  if(reversible) {
    out.plot=  ggraph::ggraph(plot.graph, layout="sugiyama", layers=layers) +
      ggraph::geom_edge_arc(ggplot2::aes(edge_width=Flux, edge_alpha=Flux, label=label, circular = FALSE),
                            strength = 0.05,
                            label_size = 3, label_colour="black", color="#AAAAFF",
                            label_parse = TRUE, angle_calc = "along", check_overlap = TRUE) +
      ggraph::scale_edge_width(limits=c(0,NA)) + ggraph::scale_edge_alpha(limits=c(0,NA)) +
      ggraph::theme_graph(base_family="sans")
  } else if(uncertainty) {
    library(ggraph)
    cvs = fluxes$FluxSD/fluxes$Flux
    maxcv = max( max(cvs), 0.5 )
    out.plot = ggraph::ggraph(plot.graph, layout="sugiyama", layers=layers) +
      ggraph::geom_edge_link(ggplot2::aes(edge_width=Flux, edge_color=FluxSD/Flux, label=label),
                             label_size = 3, label_colour="black",
                             label_parse = FALSE, angle_calc = "along", check_overlap = TRUE) +
      ggraph::scale_edge_width(limits=c(0,NA)) + 
      ggraph::scale_edge_color_gradient(name = "CV", low = "#AAAAFF", high = "#FFAAAA", na.value = "lightgrey", limits=c(0,maxcv)) +
      ggraph::theme_graph(base_family="sans")
  } else {
    out.plot = ggraph::ggraph(plot.graph, layout="sugiyama", layers=layers) +
      ggraph::geom_edge_link(ggplot2::aes(edge_width=Flux, edge_alpha=Flux, label=label),
                             label_size = 3, label_colour="black", color="#AAAAFF",
                             label_parse = FALSE, angle_calc = "along", check_overlap = TRUE) +
      ggraph::scale_edge_width(limits=c(0,NA)) + ggraph::scale_edge_alpha(limits=c(0,NA)) +
      ggraph::theme_graph(base_family="sans")
  }
  return(out.plot)
}

#' Plot some data for accumulation modelling
#'
#' @param data A required matrix or data.frame containing binary observations
#' @param tree Optional tree object
#' @param hjust Numeric (default 1), horizontal justification for rotated feature labels
#' @param bmargin Numeric (default 40), bottom margin size to prevent labels being truncated
#' @param auto.cluster Boolean (default FALSE) whether to cluster the dataset by similarity to estimate relatedness
#' @param auto.cluster.method Character (default "clustering") method to use for clustering; see binary_phylogeny
#' @param ... other options to pass to plotHypercube.curated.tree (if tree is provided)
#'
#' @return A ggplot object
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' plot_hyperinf_data(data)
#' @export
plot_hyperinf_data <- function(data,
                               tree = NULL,
                               feature.names = TRUE,
                               hjust = 1,
                               bmargin = 40,
                               auto.cluster = FALSE,
                               auto.cluster.method = "clustering",
                               ...) {
  mat = clean_data(data)
  
  L = ncol(mat)
  mat[mat == 2] = 0.5
  if(is.matrix(data) & !is.null(tree)) {
    if(is.null(rownames(data))) {
      message("You've given me a matrix without rownames and a tree, but to plot tree-linked data I need a dataframe with IDs in the first column!")
      tree = NULL
    }
  }
  if(length(feature.names) == 1) {
    if(feature.names == TRUE) {
      feature.names = colnames(data[,(ncol(data)-L+1):ncol(data)])
    } else {
      feature.names = NULL
    }
  } else if(length(feature.names) != L) {
    feature.names = NULL
  }
    
  if(auto.cluster == TRUE) {
    if(!is.null(tree)) {
      message("You provided a tree but also asked me to auto-cluster the data. I'm doing the latter and ignoring the tree.")
    }
    clustered = binary_phylogeny(mat, method=auto.cluster.method)
    data = clustered$data
    mat = clean_data(data)
    tree = clustered$tree
  }
  if(is.null(feature.names)) {
    feature.names = 1:L
  }
  if(!is.null(tree)) {
    if(is.matrix(data)) {
      df = cbind(data.frame(ID = rownames(data)), as.data.frame(mat))
    } else {
      df = cbind(data.frame(ID = data[,1]), as.data.frame(mat))
      colnames(df)[1] = "ID"
    }
    if(!setequal(df$ID, tree$tip.label)) {
      to.keep = which(tree$tip.label %in% df$ID)
      tree = ape::keep.tip(tree, to.keep)
      if(!setequal(df$ID, tree$tip.label)) {
        message("Warning: data IDs don't match tree tip labels. Unexpected behaviour will result.")
      } else {
        message("More tree tips than observations... dropping those without records.")
      }
    }
    df[df == 0.5] = "?"
    c.tree = hyperlau::curate.uncertain.tree(tree, df, independent.transitions = FALSE)
    colnames(c.tree$data)[1] = "label"
    colnames(c.tree$data)[2:ncol(c.tree$data)] = feature.names 
    out.plot = hypertrapsct::plotHypercube.curated.tree(c.tree, scale.fn = NULL, hjust=hjust, ...) +
      ggplot2::theme(plot.margin = ggplot2::margin(t = 5, r = 5, b = bmargin, l = 5)) +
      ggplot2::coord_cartesian(clip = "off")
  } else {
    df <- expand.grid(
      x = seq_len(ncol(mat)),
      y = seq_len(nrow(mat))
    )
    df$value <- as.vector(t(mat))
    
    df$color <- ifelse(
      df$value == 1, "one",
      ifelse(df$value == 0, "zero", "other")
    )
    
    out.plot = ggplot2::ggplot(df, ggplot2::aes(x, y, fill = color)) +
      ggplot2::geom_tile(width = 0.95, height = 0.95) +
      ggplot2::scale_fill_manual(
        values = c(
          one = "#888888",
          zero = "white",
          other = "red"
        )
      ) +
      ggplot2::coord_fixed() +
      ggplot2::scale_y_reverse() +
      ggplot2::theme_minimal() + 
      ggplot2::theme(legend.position="none",
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y  = ggplot2::element_blank()) +
      ggplot2::labs(x="Feature", y="Samples")
    
    if(length(feature.names) == L) {
      out.plot = out.plot + ggplot2::scale_x_continuous(breaks = 1:L,
                                               labels = feature.names) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust=hjust),
                       plot.margin = ggplot2::margin(t = 5, r = 5, b = bmargin, l = 5)) +
        ggplot2::coord_cartesian(clip = "off")
        
    }
  }
  
  return(out.plot)
}

#' Plot interactions from a fitted hypercubic inference model
#'
#' @param fit A fitted hypercubic inference model (output from hyperinf)
#' @param threshold Double (default 0.05), interaction strength below which edges will not be plotted
#' @param cv.threshold Double (default 1), coefficient of variation for interaction strength above which edges will not be plotted
#'
#' @return A ggplot object
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit = hyperinf(data)
#' plot_hyperinf_interactions(fit)
#' @export
plot_hyperinf_interactions = function(fit,
                                      threshold = 0.1,
                                      cv.threshold = 0.5) {
  if("best.graph" %in% names(fit)) {
    fit.type = "DAG"
    message("HyperDAGs not supported yet")
    return(ggplot2::ggplot())
  } else if("raw.graph" %in% names(fit)) {
    fit.type = "arborescence"
    message("HyperDAGs not supported yet")
    return(ggplot2::ggplot())
  } else if("posterior.samples" %in% names(fit)) {
    fit.type = "hypertraps"
    working.fit = fit
  } else if("Dynamics" %in% names(fit)) {
    fit.type = "hyperlau"
    working.fit = full_to_squared_fit(fit)
  } else if("viz" %in% names(fit)) {
    fit.type = "hyperhmm"
    working.fit = full_to_squared_fit(fit)
  } else if("fitted_mk" %in% names(fit)) {
    fit.type = "mk"
    working.fit = full_to_squared_fit(fit)
  } else {
    message("Didn't recognise this model")
    return(ggplot2::ggplot())
  }
  
  hypertrapsct::plotHypercube.influencegraph(working.fit,
                                             thresh = threshold,
                                             cv.thresh = cv.threshold)
}

