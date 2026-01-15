# simply returns a binary (character string) of length len from a decimal
DecToBin <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}

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
    invalid_entries <- !first_col_char %in% c("0", "1", "?")
    
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
  mat <- apply(mat, c(1,2), function(x) as.numeric(x))
  
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
                     ...) {
  
  
  # TO DO -- MAKE SURE DATA ALIGNS WITH TIPS
  #       -- TIMINGS
  #       -- NOT SURE UNCERTAINTY'S HANDLED RIGHT
  
  mat = clean_data(data)
  
  if(losses == TRUE) {
    mat = 1-mat
  }
  L = ncol(mat)
  
  uncertainty = FALSE
  if(any(mat == 2) | any(mat == -1)) {
    uncertainty = TRUE
  }
  
  #### XXX NEEDS FIXING
  
  if(!is.null(tree)) {
    df = cbind(data.frame(ID = data[,1]), as.data.frame(mat))
    if(any(mat == 2) | any(mat == -1)) {
      it = TRUE
      dots <- list(...)
      if ("independent.transitions" %in% names(dots)) {
        it = dots$independent.transitions
      }
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
    if(method == "hypermk" & L > 7) {
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
  
  if(method == "hypermk") {
    if (!is.null(tree)) {
      fit = hypermk::mk_infer_phylogenetic(mat, tree, reversible = reversible)
    } else {
      fit = hypermk::mk_infer_cross_sectional(mat, reversible = reversible)
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
      fit = do.call(hyperhmm::HyperHMM, c(list(obs = dests, initialstates = srcs), dots))
    } else {
      fit = do.call(hyperhmm::HyperHMM, c(list(obs = mat), dots))
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
    } else {
      fit = do.call(hypertrapsct::HyperTraPS, c(list(obs = mat, pli=pli), dots))
    }
  } else if(method == "hyperlau") {
    dots <- list(...)
    dots = dots[names(dots) != "independent.transitions"]
    if(!is.null(tree)) {
      fit = do.call(hyperlau::HyperLAU, c(list(obs = c.tree$dests, Xinitialstates = c.tree$srcs), dots))
    } else {
      fit = do.call(hyperlau::HyperLAU, c(list(obs = mat), dots))
    }
  }
  else if(method == "hyperdags") {
    if(!is.null(tree)) {
      srcs = apply(c.tree$srcs, 1, paste0, collapse="")
      dests = apply(c.tree$dests, 1, paste0, collapse="")
      fit = hyperdags::simplest_DAG(srcs, dests)
    } else {
      dests = apply(mat, 1, paste0, collapse="")
      fit = hyperdags::simplest_arborescence(dests)
    }
  }
  return(fit)
}

#' Plot a fitted hypercubic inference model
#'
#' @param fit A fitted hypercubic inference model (output from hyperinf)
#' @param plot.type A character string, either empty (default) to allow standardised plot, or "native" to produce plots from the source algorithm
#' @param threshold Double (default 0.05), probability flux below which edges will not be plotted
#' @param uncertainty Boolean, whether to visualise uncertainty over bootstraps (only for HyperLAU and HyperHMM)
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
                         uncertainty = "") {
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
    return(ggplot2::ggplot())
  }
  
  if(fit.type == "hyperlau") {
    if(length(unique(fit$Dynamics$Bootstrap)) > 1) {
      uncertainty = TRUE
    }
  } else if(fit.type == "hyperhmm") {
    if(length(unique(fit$transitions$Bootstrap)) > 1) {
      uncertainty = TRUE
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
    if(fit.type %in% c("mk", "hyperhmm", "hyperlau", "hypertraps")) {
      # our goal is now to get a From/To/Flux dataframe and eventually a graph to plot
      if(fit.type == "mk") {
        fluxes = fit$mk_fluxes
        fluxes$Flux = fluxes$Flux/sum(fluxes$Flux[fluxes$From == 0])
        if(any(fluxes$From > fluxes$To)) {
          reversible = TRUE
        }
        # decimal, 0-indexed labels
      } else if(fit.type == "hyperhmm") {
        if(uncertainty == FALSE) {
          fluxes = fit$transitions[fit$transitions$Bootstrap == 0, 2:ncol(fit$transitions)]
        } else {
          fluxes <- fit$transitions %>%
            # ensure all From–To pairs exist for every Bootstrap
            tidyr::complete(
              Bootstrap,
              From,
              To,
              fill = list(Flux = 0)
            ) %>%
            dplyr::group_by(From, To) %>%
            dplyr::summarise(
              mean_flux = mean(Flux),
              sd_flux = sd(Flux),
              .groups = "drop"
            )
          colnames(fluxes) = c("From", "To", "Flux", "FluxSD")
        }
        # decimal, 0-indexed labels
      } else if(fit.type == "hyperlau") {
        if(uncertainty == FALSE) {
          fluxes = fit$Dynamics[fit$Dynamics$Bootstrap == 0, 2:ncol(fit$Dynamics)]
        } else {
          fluxes = fit$Dynamics %>%
            # ensure all From–To pairs exist for every Bootstrap
            tidyr::complete(
              Bootstrap,
              From,
              To,
              fill = list(Flux = 0)
            ) %>%
            dplyr::group_by(From, To) %>%
            dplyr::summarise(
              mean_flux = mean(Flux),
              sd_flux = sd(Flux),
              .groups = "drop"
            )
          colnames(fluxes) = c("From", "To", "Flux", "FluxSD")
        }
        # decimal, 0-indexed labels
      } else if(fit.type == "hypertraps") {
        fluxes = fit$dynamics$trans
        # decimal, 0-indexed labels
      }
      
      L = fit$L
      fluxes$From.b = sapply(fluxes$From, DecToBin, len=L)
      fluxes$To.b = sapply(fluxes$To, DecToBin, len=L)
      fluxes$Change = L-log(abs(fluxes$From-fluxes$To), base=2)
      fluxes$label = paste0("+", fluxes$Change)
      if(reversible == TRUE) {
        fluxes$label[which(fluxes$From > fluxes$To)] = paste0("-", fluxes$Change[which(fluxes$From > fluxes$To)])
      }
      fluxes = fluxes[fluxes$Flux > threshold,]
      states = unique(c(fluxes$From, fluxes$To))
      states.b = unique(c(fluxes$From.b, fluxes$To.b))
      layers = sapply(states.b, stringr::str_count, pattern="1")
      names(layers) = states
      plot.graph = igraph::graph_from_data_frame(fluxes)
      
    } else if(fit.type == "DAG" | fit.type == "arborescence") {
      if(fit.type == "arborescence") {
        graphD = fit$rewired.graph
      } else {
        graphD = fit$best.graph
      }
      graphD.layers = sapply(igraph::V(graphD)$name, stringr::str_count, "1")
      L = stringr::str_length(igraph::V(graphD)$name[1])
      labels = as.character(1:L)
      graphD.size = igraph::neighborhood.size(graphD, L+1, mode="out")
      igraph::E(graphD)$Flux = as.numeric(graphD.size[igraph::ends(graphD, es = igraph::E(graphD), names = FALSE)[, 2]])
      igraph::E(graphD)$Flux = igraph::E(graphD)$Flux/max(igraph::E(graphD)$Flux)
      this.ends = igraph::ends(graphD, es=igraph::E(graphD))
      srcs = strsplit(this.ends[,1], split="")
      dests = strsplit(this.ends[,2], split="")
      for(i in 1:nrow(this.ends)) {
        igraph::E(graphD)$label[i] = paste0("+", paste0(labels[which(srcs[[i]]!=dests[[i]])], collapse="\n"), collapse="")
      }
      
      plot.graph = graphD
      layers = graphD.layers
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
      out.plot = ggraph::ggraph(plot.graph, layout="sugiyama", layers=layers) +
        ggraph::geom_edge_link(ggplot2::aes(edge_width=Flux, edge_color=FluxSD/Flux, label=label),
                               label_size = 3, label_colour="black",
                               label_parse = TRUE, angle_calc = "along", check_overlap = TRUE) +
        ggraph::scale_edge_width(limits=c(0,NA)) + 
        ggraph::scale_edge_color_gradient(name = "CV", low = "#AAAAFF", high = "#FFAAAA", na.value = "lightgrey", limits=c(0,NA)) +
        ggraph::theme_graph(base_family="sans")
    } else {
      out.plot = ggraph::ggraph(plot.graph, layout="sugiyama", layers=layers) +
        ggraph::geom_edge_link(ggplot2::aes(edge_width=Flux, edge_alpha=Flux, label=label),
                               label_size = 3, label_colour="black", color="#AAAAFF",
                               label_parse = TRUE, angle_calc = "along", check_overlap = TRUE) +
        ggraph::scale_edge_width(limits=c(0,NA)) + ggraph::scale_edge_alpha(limits=c(0,NA)) +
        ggraph::theme_graph(base_family="sans")
    }
  }
  return(out.plot)
}

#' Plot some data for accumulation modelling
#'
#' @param data A required matrix or data.frame containing binary observations
#' @param tree Optional tree object
#' @param ... other options to pass to plotHypercube.curated.tree (if tree is provided)
#'
#' @return A ggplot object
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' plot_hyperinf_data(data)
#' @export
plot_hyperinf_data <- function(data,
                               tree = NULL,
                               ...) {
  mat = clean_data(data)
  mat[mat == 2] = 0.5
  if(!is.null(tree)) {
    df = cbind(data.frame(ID = tree$tip.label), as.data.frame(mat))
    c.tree = hypertrapsct::curate.tree(tree, df)
    out.plot = hypertrapsct::plotHypercube.curated.tree(c.tree, scale.fn = NULL, ...)
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
  }
  
  return(out.plot)
}

