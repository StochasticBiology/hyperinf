#' Produce dataframe of ordering probabilities
#'
#' @param fit A fitted hypercubic inference model (output from hyperinf)
#'
#' @return A dataframe with columns feature/order/prob
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit.1 = hyperinf(data)
#' hyperinf_bubbles(fit.1)
#' @export
hyperinf_bubbles = function(fit) {
  this.fit = fit
  if("best.graph" %in% names(this.fit)) {
    fit.type = "DAG"
  } else if("raw.graph" %in% names(this.fit)) {
    fit.type = "arborescence"
  } else if("posterior.samples" %in% names(this.fit)) {
    fit.type = "hypertraps"
  } else if("Dynamics" %in% names(this.fit)) {
    fit.type = "hyperlau"
  } else if("viz" %in% names(this.fit)) {
    fit.type = "hyperhmm"
  } else if("fitted_mk" %in% names(this.fit)) {
    fit.type = "mk"
  } else {
    message("Didn't recognise this model")
    return(ggplot2::ggplot())
  }
  if(!(fit.type %in% c("hyperhmm", "hypertraps", "mk", "hyperlau"))) {
    message("This fit type not yet supported!")
    stop()
  }
  if(fit.type == "mk") {
    c.df = old.df = this.fit$mk_fluxes
    c.df = c.df[c.df$To > c.df$From,]
    c.df$change = this.fit$L-log(c.df$To-c.df$From, base=2)
    bins = sapply(c.df$From, DecToBin, len=this.fit$L)
    c.df$level = colSums(bins)
    for(l in unique(c.df$level)) {
      l.sum = sum(c.df$Flux[c.df$level == l])
      c.df$Flux[c.df$level == l] = c.df$Flux[c.df$level == l]/l.sum
    } 
    bp = data.frame(feature=c.df$change,
                    order=c.df$level+1,
                    prob=c.df$Flux)
  }
  if(fit.type == "hyperlau") {
    c.df = this.fit$Dynamics
    c.df$change = this.fit$L-log(c.df$To-c.df$From, base=2)
    bins = sapply(c.df$From, DecToBin, len=this.fit$L)
    c.df$level = colSums(bins)
    bp = data.frame(feature=c.df$change,
                    order=c.df$level+1,
                    prob=c.df$Flux)
  }
  if(fit.type == "hyperhmm") {
    bp = this.fit$stats
    bp$prob = bp$mean
  }
  if(fit.type == "hypertraps") {
    bp = data.frame(feature=1+this.fit$bubbles$OriginalIndex,
                    order=1+this.fit$bubbles$Time,
                    prob=this.fit$bubbles$Probability)
  }
  return(bp[,c("feature","order","prob")])
}

#' Compare bubble plots for a set of fitted hypercubic inference models
#'
#' @param fits A list of fitted hypercubic inference models (output from hyperinf)
#' @param reorder Boolean (default FALSE), whether to reorder features by mean acquisition order
#' @param transpose Boolean (default FALSE), whether to transpose axes for the bubble plot
#' @param thetastep Integer (default 10), the number of angular steps to take when drawing polygons for each bubble segment. Lower this for many comparisons.
#' @param p.scale Numeric (default 1), the scaling factor to apply to the radius (probability) of bubbles
#' @param sqrt.trans Boolean (default FALSE), whether to sqrt-transform the probability to get the bubble radius
#' @param bins Integer (default 0, off), how many bins to divide the ordering axis into
#' @param expt.names Optional vector of labels for each element of the fit list
#' @param feature.names Optional vector of labels for the L features in each fit
#' @param fill.name Character (default "Experiment"), the label to use for the group type
#' 
#' @return A ggplot object
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit.1 = hyperinf(data)
#' fit.2 = hyperinf(1-data)
#' plot_hyperinf_bubbles(list(fit.1, fit.2))
#' @export
plot_hyperinf_bubbles = function(fits,
                                 reorder = FALSE, transpose = FALSE, 
                                 thetastep = 10, p.scale = 1, 
                                 sqrt.trans = FALSE, bins = 0, 
                                 expt.names = NULL, fill.name = "Experiment",
                                 feature.names = NULL) 
{
  if (bins == 0) {
    toplot = data.frame()
    for (i in 1:length(fits)) {
      tmp = hyperinf_bubbles(fits[[i]])
      tmp$expt = i
      if (length(expt.names) == 0) {
        tmp$exptname = i
      }
      else {
        tmp$exptname = expt.names[i]
      }
      toplot = rbind(toplot, tmp)
    }
  }
  else {
    toplot = data.frame()
    for (i in 1:length(fits)) {
      tmp = hyperinf_bubbles(fits[[i]])
      tmax = max(tmp$order)
      for (this.name in unique(tmp$feature)) {
        for (j in 0:(bins - 1)) {
          prob = sum(tmp$prob[tmp$feature == 
                                this.name & round((bins - 1) * tmp$order/tmax) == 
                                j])
          toplot = rbind(toplot, data.frame(feature = this.name, 
                                            order = j, prob = prob, expt = i))
        }
      }
    }
  }
  if (reorder == TRUE) {
    toplot$feature = factor(toplot$feature, levels = unique(toplot$feature))
  }
  if (transpose == TRUE) {
    toplot$x = toplot$feature
    toplot$y = toplot$order
  }
  else {
    toplot$x = toplot$order
    toplot$y = toplot$feature
  }
  if (sqrt.trans == TRUE) {
    toplot$prob = sqrt(toplot$prob)
  }
  polygons = data.frame()
  for (i in 1:nrow(toplot)) {
    thisx = as.numeric(toplot$x[i])
    thisy = as.numeric(toplot$feature[i])
    thisz = toplot$prob[i] * p.scale
    theta0 = 2 * pi * (toplot$expt[i] - 1)/max(toplot$expt)
    dtheta = (2 * pi/max(toplot$expt))/thetastep
    theta = theta0
    tmp = data.frame()
    tmp = rbind(tmp, data.frame(x1 = thisx, y1 = thisy, ref = i, 
                                expt = toplot$expt[i], exptname = toplot$exptname[i]))
    for (j in 0:(thetastep)) {
      theta = theta0 + j * dtheta
      tmp = rbind(tmp, data.frame(x1 = thisx + thisz * 
                                    cos(theta), y1 = thisy + thisz * sin(theta), 
                                  ref = i, expt = toplot$expt[i], exptname = toplot$exptname[i]))
    }
    polygons = rbind(polygons, tmp)
  }
  this.plot = ggplot2::ggplot(polygons, ggplot2::aes(x = x1, 
                                                     y = y1, group = ref, 
                                                     fill = factor(exptname))) + 
    ggplot2::geom_polygon() + ggplot2::theme_light() + ggplot2::labs(fill = fill.name, 
                                                                     x = "Ordinal Time", y = "")
  if(length(feature.names) != 0) {
    this.plot = this.plot + scale_y_continuous(breaks=1:fits[[1]]$L, labels = feature.names)
  }
  return(this.plot)
}

#' Compare transition graphs for a set of fitted hypercubic inference models
#'
#' @param fits A list of fitted hypercubic inference models (output from hyperinf)
#' @param threshold Double (default 0.05), probability flux below which edges will not be plotted
#' @param expt.names Optional vector of labels for each element of the fit list
#' @param feature.names Optional vector of labels for the L features in each fit
#' @param style Character (default "limited") giving plot style. "limited" condenses bootstrap resamples into single edges; others plot each as a different arc.
#' @param bend Numeric (default 0.5) the strength of the bend separating edges between the same pair of nodes
#' @param label_size Numeric (default 2) the size of edge labels
#'
#' @return A ggplot object
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit.1 = hyperinf(data)
#' fit.2 = hyperinf(1-data)
#' plot_hyperinf_comparative(list(fit.1, fit.2))
#' @export
plot_hyperinf_comparative = function(fits, threshold=0.05,
                                     expt.names=NULL,
                                     feature.names=NULL,
                                     style="limited",
                                     bend = 0.5,
                                     label_size = 2) {
  library(ggraph)
  plot.graphs = layers = es = list()
  i = 1
  for(this.fit in fits) {
    if("best.graph" %in% names(this.fit)) {
      fit.type = "DAG"
    } else if("raw.graph" %in% names(this.fit)) {
      fit.type = "arborescence"
    } else if("posterior.samples" %in% names(this.fit)) {
      fit.type = "hypertraps"
    } else if("Dynamics" %in% names(this.fit)) {
      fit.type = "hyperlau"
    } else if("viz" %in% names(this.fit)) {
      fit.type = "hyperhmm"
    } else if("fitted_mk" %in% names(this.fit)) {
      fit.type = "mk"
    } else {
      message("Didn't recognise this model")
      return(ggplot2::ggplot())
    }
    
    tmp = get_plot_graph(this.fit, fit.type)
    plot.graphs[[i]] = tmp[["plot.graph"]]
    layers[[i]] = tmp[["layers"]]
    if(length(expt.names) != length(fits)) { this.src = i }
    else { this.src = expt.names[i] }
    es[[i]] = igraph::as_data_frame(plot.graphs[[i]], "edges")
    es[[i]]$src = this.src
    i = i+1
  }
  edges = dplyr::bind_rows(es)
  edges = edges[edges$Flux > threshold,]
  if(length(feature.names) > 0) {
    edges$label = paste("+", feature.names[edges$Change], sep="")
  }
  
  if(style == "limited") {
  uniques = which(!duplicated(edges[,c("from", "to", "src")]))
  new.edges = data.frame()
  for(u in uniques) {
    tmp = edges[u,]
    tmp.set = which(edges$from==tmp$from &
                      edges$to==tmp$to &
                      edges$src==tmp$src)
    new.edges = rbind(new.edges, data.frame(from=tmp$from,
                                            to=tmp$to,
                                            src=tmp$src,
                                            label=tmp$label,
                                            nboot=length(tmp.set),
                                            mean=mean(edges$Flux[tmp.set])))
  }
  new.edges$label[new.edges$nboot < 2] = ""
  }
  
  u.layers = which(!duplicated(names(unlist(layers))))
  use.layers = unlist(layers)[u.layers]
  use.layers = use.layers[which(names(use.layers) %in% igraph::union(edges$from, edges$to))]
  if(style != "limited") {
    g_combined <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = NULL)
  } else {
    g_combined <- igraph::graph_from_data_frame(new.edges, directed = TRUE, vertices = NULL)
  }  
  
  # Ensure src is a factor
  igraph::E(g_combined)$Experiment <- factor(igraph::E(g_combined)$src)
  
  if(style == "limited") {
  g.plot = ggraph::ggraph(g_combined, layout = "sugiyama") +
    ggraph::geom_edge_fan(ggplot2::aes(
      edge_width = mean,
      edge_alpha = nboot,
      color = Experiment,
      label=label,
    ),
    strength = bend, check_overlap=TRUE, label_size=label_size
    ) +
    ggraph::scale_edge_width(range = c(0.1, 4)) +
    ggraph::scale_edge_alpha(range = c(0, 0.4)) +
    ggraph::theme_graph(base_family = "sans") 
  } else {
    g.plot = ggraph::ggraph(g_combined, layout = "sugiyama") +
      ggraph::geom_edge_fan(ggplot2::aes(
        edge_width = Flux,
        edge_alpha = Flux,
        color = Experiment,
        label=label,
      ),
      strength = bend, check_overlap=TRUE, label_size=label_size
      ) +
      ggraph::scale_edge_width(range = c(0.5, 2)) +
      ggraph::scale_edge_alpha(range = c(0.3, 1)) +
      ggraph::theme_graph(base_family = "sans") 
  }
  return(g.plot)
}

#' Comparative bubble plot for bootstrapped fitted models
#'
#' @param fit.1 A fitted hypercubic inference models with bootstrap info (output from hyperinf)
#' @param fit.2 A fitted hypercubic inference models with bootstrap info (output from hyperinf)
#'
#' @return A ggplot object
#' @export
plot_hyperinf_bootstrap = function(fit.1, fit.2) {
  if(!("boots" %in% names(fit.1) & "boots" %in% names(fit.2) )) {
    message("This doesn't look like a comparable pair of model fits")
    return(ggplot2::ggplot())
  }

  boots = list()
  boots = c(boots, fit.1$boots)
  boots = c(boots, fit.2$boots)

  bub.set = lapply(boots, hyperinf_bubbles)
  groups = c(rep(1, length(fit.1$boots)), rep(2, length(fit.2$boots)))
  sum.set = list()
  sum.df = data.frame()
  for(g in unique(groups)) {
    g.set = bub.set[which(groups == g)]
    result <- cbind(
      g.set[[1]][1:2],
      do.call(cbind, lapply(g.set, `[`, 3))
    )
    result$mean = apply(result[,3:ncol(result)], 1, mean)
    result$min = apply(result[,3:ncol(result)], 1, min)
    result$max = apply(result[,3:ncol(result)], 1, max)
    result$g = g
    sum.set[[g]] = result[,c("feature", "order", "mean", "min", "max")]
    sum.df = rbind(sum.df, result[,c("g", "feature", "order", "mean", "min", "max")])
  }
  m = unique(sum.df[,c("feature", "order")])
  m$sig = 0
  for(i in 1:nrow(m)) {
    g1 = sum.df[sum.df$g==1 & sum.df$feature==m$feature[i] & sum.df$order==m$order[i],]
    g2 = sum.df[sum.df$g==2 & sum.df$feature==m$feature[i] & sum.df$order==m$order[i],]
    if(g1$min > g2$max | g1$max < g2$min) { m$sig[i] = 1 }
  }

  plot_hyperinf_bubbles(boots, p.scale = 0.1) + 
    ggplot2::geom_point(data=m[m$sig==1,], ggplot2::aes(x=order, y=feature+0.15, group=1, fill="*"), shape="*", size=12) +
    ggplot2::theme(legend.position="none")
}

