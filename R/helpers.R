#' @importFrom hyperdags fit_properties
#' @export
hyperdags::fit_properties

# get a graph object reflecting a (plottable) transition network from a fitted model
get_plot_graph = function(fit, fit.type, uncertainty = FALSE, 
                          reversible = FALSE, threshold = 0.05) {
  fluxes = NULL
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
        fluxes = compute_hypertraps_fluxes(fit)
      }
      # decimal, 0-indexed labels
    }
    
    L = fit$L
    fluxes$From.b = sapply(fluxes$From, DecToBinS, len=L)
    fluxes$To.b = sapply(fluxes$To, DecToBinS, len=L)
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
      igraph::E(graphD)$label[i] = paste0(paste0("+", labels[which(srcs[[i]]!=dests[[i]])]), collapse="\n")
    }
    
    plot.graph = graphD
    layers = graphD.layers
  }
  return(list(plot.graph=plot.graph,
              layers=layers,
              fluxes=fluxes))
}
