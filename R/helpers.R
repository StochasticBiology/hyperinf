#' @importFrom hyperdags fit_properties
#' @export
hyperdags::fit_properties

#' Parameteric bootstrap for a HyperMk model fit
#'
#' @param fit A fitted HyperMk model (output from hyperinf)
#' @param n_boot Numeric (default 10) number of bootstrap resamples to run
#' @param reversible Boolean (default TRUE) whether we're running a reversible HyperMk model
#' @param flux.samples Numeric (default 10000) number of samplers to run to characterise dynamics
#' 
#' @return The original fitted model with a new element boots, containing transition matrices and dynamic summaries across bootstraps
#' @export
bootstrap_mk_fit = function(fit, 
                            n_boot = 10,
                            reversible = TRUE,
                            flux.samples = 10000) {
  if(hyperinf_gettype(fit) != "hypermk") {
    message("This only works for HyperMk")
    return(NULL)
  }
  Q_hat <- fit$fitted_mk$transition_matrix
  
  fit$boots = list()
  reversible = TRUE
  L = fit$L
  
  for(i in 1:n_boot){
    
    # simulate new tip states under fitted model
    sim_states <- castor::simulate_mk_model(
      tree        = fit$data$tree,
      Q = Q_hat,
    )
    
    # refit
    if(FALSE) {
    fit_i <- castor::fit_mk(fit$data$tree, 2^L, tip_states = sim_states$tip_states,
                            rate_model = hypermk::mk_index_matrix(L, reversible=reversible))
    
    message("simulating fluxes")
    fit$boots[[i]] = list(Q_boot = fit_i$transition_matrix,
                          mk_df = hypermk::mk_pull_transitions(fit_i, reversible = reversible),
                          mk_fluxes = hypermk::mk_simulate_fluxes(fit_i, L, reversible = reversible,
                                                                  nwalker = flux.samples)
    )
    }
    fit$boots[[i]] = hypermk::mk.inference(fit$data$tree, L, use.priors=FALSE, sim_states$tip_states, reversible=reversible)
  }
  
  return(fit)
}

#' Get the fitting method used to produce a hyperinf model fit
#'
#' @param fit A fitted hypercubic inference model (output from hyperinf)
#' 
#' @return A character string describing the fit type
#' @export
hyperinf_gettype = function(fit) {
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
    if("mk2_fluxes" %in% names(fit)) {
      fit.type = "hypermk2"
    } else {
      fit.type = "hypermk"
    }
  } else {
    message("Didn't recognise this model")
    fit.type = NULL
  }
  return(fit.type)
}

#' Get a graph object reflecting a (plottable) transition network from a fitted model
#'
#' @param fit A fitted hypercubic inference model (output from hyperinf)
#' @param fit.type (default NULL) character vector giving type of fit (overridden by inferred type)
#' @param uncertainty Boolean, are we bootstrapping
#' @param reversible Boolean, are we consider reversible steps
#' @param threshold Numeric (default 0.05), threshold to place on fluxes
#' @param feature.names NULL or character vector (default NULL). If vector of length L, use those labels.
#' @param max.samps Numeric (default 1000) samples for HyperTraPS sampling
#' 
#' @return A graph object
#' @export
get_plot_graph = function(fit, fit.type = NULL, uncertainty = FALSE, 
                          reversible = FALSE, threshold = 0.05,
                          feature.names = NULL, max.samps = 1000) {
  fluxes = NULL
  if(is.null(fit.type)) {
    fit.type = hyperinf_gettype(fit)
  }
  if(length(feature.names) != fit$L) {
    if("feature.names" %in% names(fit)) {
      feature.names = fit$feature.names
    } else {
      feature.names = as.character(1:fit$L)
    }
  }
  if(fit.type %in% c("hypermk", "hyperhmm", "hyperlau", "hypertraps")) {
    # our goal is now to get a From/To/Flux dataframe and eventually a graph to plot
    if(fit.type == "hypermk") {
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
        fluxes <- fit$transitions |>
          # ensure all From–To pairs exist for every Bootstrap
          tidyr::complete(
            Bootstrap,
            From,
            To,
            fill = list(Flux = 0)
          ) |>
          dplyr::group_by(From, To) |>
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
        fluxes = fit$Dynamics |>
          # ensure all From–To pairs exist for every Bootstrap
          tidyr::complete(
            Bootstrap,
            From,
            To,
            fill = list(Flux = 0)
          ) |>
          dplyr::group_by(From, To) |>
          dplyr::summarise(
            mean_flux = mean(Flux),
            sd_flux = sd(Flux),
            .groups = "drop"
          )
        colnames(fluxes) = c("From", "To", "Flux", "FluxSD")
      }
      # decimal, 0-indexed labels
    } else if(fit.type == "hypertraps") {
      if("dynamics" %in% names(fit)) {
        fluxes = fit$dynamics$trans 
      } else {
        message("Computing fluxes...")
        fluxes = compute_hypertraps_fluxes(fit, max.samps = max.samps)
      }
      # decimal, 0-indexed labels
    }
    fluxes = fluxes[fluxes$Flux > 0,]
    L = fit$L
    fluxes$From.b = sapply(fluxes$From, DecToBinS, len=L)
    fluxes$To.b = sapply(fluxes$To, DecToBinS, len=L)
    fluxes$Change = L-log(abs(fluxes$From-fluxes$To), base=2)
    fluxes$Changelabel = feature.names[fluxes$Change]
    fluxes$label = paste0("+", fluxes$Changelabel)
    if(reversible == TRUE) {
      fluxes$label[which(fluxes$From > fluxes$To)] = paste0("-", fluxes$Changelabel[which(fluxes$From > fluxes$To)])
    }
    fluxes = fluxes[fluxes$Flux > threshold,]
    states = unique(c(fluxes$From, fluxes$To))
    states.b = unique(c(fluxes$From.b, fluxes$To.b))
    layers = sapply(states.b, stringr::str_count, pattern="1")
    names(layers) = states
    plot.graph = igraph::graph_from_data_frame(fluxes)
    
  } else if(fit.type == "hypermk2") {
      fluxes = fit$mk2_fluxes
      fluxes$Flux = fluxes$Flux/max(fluxes$Flux)
      if(any(fluxes$From > fluxes$To)) {
        reversible = TRUE
      }
      fluxes = fluxes[fluxes$Flux > 0,]
      fluxes$From.b = sapply(fluxes$From, DecToBinS, len=fit$L)
      fluxes$To.b = sapply(fluxes$To, DecToBinS, len=fit$L)
      fluxes$Changelabel = ""
      
      for(i in 1:nrow(fluxes)) {
        src = as.numeric(strsplit(fluxes$From.b[i], "")[[1]])
        dest = as.numeric(strsplit(fluxes$To.b[i], "")[[1]])
        adds = which(src<dest)
        losses = which(dest<src)
        for(add in adds) { fluxes$Changelabel[i] = paste0(fluxes$Changelabel[i], "+", feature.names[add], ",", collapse="")  }
        for(loss in losses) { fluxes$Changelabel[i] = paste0(fluxes$Changelabel[i], "-", feature.names[loss], ",", collapse="")  }
        
      }
      fluxes$label = fluxes$Changelabel

      fluxes = fluxes[fluxes$Flux > threshold,]
      states = unique(c(fluxes$From, fluxes$To))
      states.b = unique(c(fluxes$From.b, fluxes$To.b))
      layers = sapply(states.b, stringr::str_count, pattern="1")
      names(layers) = states
      plot.graph = igraph::graph_from_data_frame(fluxes)
      
    }
  else if(fit.type == "DAG" | fit.type == "arborescence") {
    if(fit.type == "arborescence") {
      graphD = fit$rewired.graph
    } else {
      graphD = fit$best.graph
    }
    graphD.layers = sapply(igraph::V(graphD)$name, stringr::str_count, "1")
    L = stringr::str_length(igraph::V(graphD)$name[1])
    this.ends = igraph::ends(graphD, es=igraph::E(graphD))
    srcs = strsplit(this.ends[,1], split="")
    dests = strsplit(this.ends[,2], split="")
    
    for(i in 1:length(igraph::V(graphD)$name)) {
      igraph::V(graphD)$name[i] = BinToDecS(igraph::V(graphD)$name[i])
    }
    names(graphD.layers) = igraph::V(graphD)$name
    
    labels = feature.names #as.character(1:L)
    graphD.size = igraph::neighborhood.size(graphD, L+1, mode="out")
    igraph::E(graphD)$Flux = as.numeric(graphD.size[igraph::ends(graphD, es = igraph::E(graphD), names = FALSE)[, 2]])
    igraph::E(graphD)$Flux = igraph::E(graphD)$Flux/max(igraph::E(graphD)$Flux)
    
    for(i in 1:nrow(this.ends)) {
      igraph::E(graphD)$label[i] = paste0(paste0("+", labels[which(srcs[[i]]!=dests[[i]])]), collapse="\n")
    }
    
    plot.graph = graphD
    layers = graphD.layers
  }
  return(list(plot.graph=plot.graph,
              layers=layers,
              fluxes=fluxes))
}
