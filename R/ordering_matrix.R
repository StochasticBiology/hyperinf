#' Compute a feature ordering matrix for a hypercubic inference model fit
#'
#' @param fit A fitted hypercubic inference model (output from hyperinf)
#' @param n.samples Integer (default 10k) number of samples to characterise dynamics in the reversible (HyperMk) case
#' @param type Character string (default "relative"). Either "relative", in which case P_ij gives the probability that feature i is acquired before feature j. Or "absolute", in which case P_ij is the probability that feature i is acquired after step j.
#' 
#' @return A matrix giving P_ij according to the type of comparison (see above)
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit.1 = hyperinf(data)
#' ordering_matrix(fit.1)
#' ordering_matrix(fit.1, type="absolute")
#' @export
ordering_matrix = function(fit, n.samples = 10000,
                           type = "relative") {
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
    fit.rev = this.fit
    ddf = fit.rev$mk_fluxes
    count = 0
    m = matrix(0, nrow=fit.rev$L, ncol=fit.rev$L)
    b4m = matrix(0, nrow=fit.rev$L, ncol=fit.rev$L)
    for(i in 1:n.samples) {
      state = 0
      for(j in 1:fit.rev$L) {
        outs = which(ddf$From == state)
        next.trans = sample(outs, 1, prob = ddf$Rate[outs])
        state.v = DecToBin(state, fit.rev$L)
        ones = which(state.v==1)
        level = sum(state.v)
        if(ddf$To[next.trans] > state) {
          change = fit.rev$L - log(ddf$To[next.trans]-state, base=2)
          if(level == 5 & change == 4) {
            cat(state.v, " ", next.trans, "\n")
          }
          m[level+1,change] = m[level+1,change]+1
          b4m[ones,change] = b4m[ones,change]+1
          count = count+1
        }
        state = ddf$To[next.trans]
      }
    }
  } else {
    fit.non.rev = this.fit
    #count
    #lvls = lvls/sum(lvls)
    
    # b4m[i,j] gets incremented if we have i when we get j
    
    if(fit.type == "hyperhmm") {
      ddf = fit.non.rev$transitions
    } else if(fit.type == "hypertraps") {
      if("routes" %in% names(fit.non.rev)) {
        routes_mat = fit.non.rev$routes
        m = matrix(0, nrow = ncol(routes_mat), ncol = ncol(routes_mat))
        for(i in 1:nrow(routes_mat)) {
          for(j in 1:ncol(routes_mat)) {
            if(type == "relative") {
              m[routes_mat[i,j]+1, 1+routes_mat[i,1:(j-1)]] = m[routes_mat[i,j]+1, 1+routes_mat[i,1:(j-1)]] + 1 
            } else {
              if(j > 1) {
                m[routes_mat[i,j]+1, 1:(j-1)] = m[routes_mat[i,j]+1, 1:(j-1)] + 1
              }
            }
          }
        }
        m = m/nrow(routes_mat)
        return(m)
      }
      ddf = fit.non.rev$dynamics$trans
    } else if(fit.type == "hyperlau") {
      ddf = fit.non.rev$Dynamics
    }
    if(is.null(ddf)) {
      message("Couldn't find transitions data frame in this model!")
      return(NULL)
    }
    if("p.boot" %in% colnames(ddf)) {
      message("Just looking at first resample")
      ddf = ddf[ddf$p.boot == 1,]
    }
    if("Bootstrap" %in% colnames(ddf)) {
      if(max(ddf$Bootstrap) > 0) {
        message("Just looking at first resample")
        ddf = ddf[ddf$Bootstrap == 0,]
      }
    }
    count = 0
    m = matrix(0, nrow=fit.non.rev$L, ncol=fit.non.rev$L)
    b4m = matrix(0, nrow=fit.non.rev$L, ncol=fit.non.rev$L)
    for(i in 1:nrow(ddf)) {
      state = ddf$From[i]
      state.v = DecToBin(state, fit.non.rev$L)
      next.state = ddf$To[i]
      change = fit.non.rev$L - log(next.state-state, base=2)
      level = sum(state.v)
      ones = which(state.v==1)
      m[level+1,change] = m[level+1,change]+1
      b4m[ones,change] =  b4m[ones,change] + ddf$Flux[i]
    }
  }
  
  for(i in 1:nrow(b4m)) {
    for(j in 1:ncol(b4m)) {
      if(i > j) {
        tot = b4m[i,j]+b4m[j,i]
        b4m[i,j] = b4m[i,j]/tot
        b4m[j,i] = b4m[j,i]/tot
        tot = m[i,j]+m[j,i]
        m[i,j] = m[i,j]/tot
        m[j,i] = m[j,i]/tot
      }
    }
  }
  if(type == "relative") {
    return(b4m)
  } else {
    return(m)
  }
}

#' Case multiple model fits into a form labelled as bootstrap samples (for compatibility) 
#'
#' @param fits A list of fitted hypercubic inference models (output from hyperinf)
#' 
#' @return A hypercubic inference model structure with the model fits stored as "bootstrap resamples"
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit.1 = hyperinf(data, method="hypertraps", seed=1)
#' fit.2 = hyperinf(data, method="hypertraps", seed=2)
#' fit = multiple_fits_to_booted_fit(list(fit.1, fit.2))
#' @export
multiple_fits_to_booted_fit = function(fits) {
  new.fit = list()
  new.fit$boots = list()
  for(i in 1:length(fits)) {
    new.fit$boots[[i]] = fits[[i]]
  }
  new.fit$L = fits[[i]]$L
  return(new.fit)
}


#' Compare ordering statistics from two fitted models
#'
#' @param fit.1 A fitted hypercubic inference model including bootstrap resamples (output from hyperinf)
#' @param fit.2 Another fitted hypercubic inference model including bootstrap resamples  (output from hyperinf)
#' @param n.samples Integer (default 10k) number of samples to characterise dynamics in the reversible (HyperMk) case
#' @param threshold Numeric (default 0). A threshold for differences. If zero, any separation of statistic (e.g. non-overlapping bootstrap distributions) is reported. If nonzero, we require one fit's statistic > 1-threshold and the other fit's statistic < threshold. 
#' @param type Character string (default "relative"). Either "relative", in which case P_ij gives the probability that feature i is acquired before feature j. Or "absolute", in which case P_ij is the probability that feature i is acquired after step j.
#' 
#' @return A hypercubic inference model structure with the model fits stored as "bootstrap resamples"
#' @export
compare_orderings = function(fit.1, fit.2, 
                             n.samples = 10000, 
                             threshold = 0,
                             type = "relative") {
  oms.1 = oms.2 = list()
  if(!("boots" %in% names(fit.1) & "boots" %in% names(fit.2))) {
    message("Didn't find bootstrap resamples in these fits")
    return(NULL)
  }
  message("Computing orderings 1")
  for(i in 1:length(fit.1$boots)) {
    oms.1[[i]] = ordering_matrix(fit.1$boots[[i]], n.samples = n.samples, type=type)
    if(is.null(oms.1[[i]])) {
      message("Missing ordering matrix!")
      return(NULL)
    }
  }
  message("Computing orderings 2")
  for(i in 1:length(fit.2$boots)) {
    oms.2[[i]] = ordering_matrix(fit.2$boots[[i]], n.samples = n.samples, type=type)
    if(is.null(oms.2[[i]])) {
      message("Missing ordering matrix!")
      return(NULL)
    }
  }
  # Stack into array
  arr.1 <- simplify2array(oms.1)
  arr.2 <- simplify2array(oms.2)
  
  # Min/max across 3rd dimension (list index)
  min_mat.1 <- apply(arr.1, c(1, 2), min)
  max_mat.1 <- apply(arr.1, c(1, 2), max)
  min_mat.2 <- apply(arr.2, c(1, 2), min)
  max_mat.2 <- apply(arr.2, c(1, 2), max)
  
  if(threshold == 0) {
    win.1 = which((min_mat.1 > max_mat.2) , arr.ind = TRUE)
    win.2 = which((min_mat.2 > max_mat.1) , arr.ind = TRUE)
    win.1 = data.frame(first=1, win.1)
    win.2 = data.frame(first=2, win.2)
    return(rbind(win.1,win.2))
  } else {
    return(which((min_mat.1 > 1-threshold & max_mat.2 < threshold) 
                 | (min_mat.2 > 1-threshold & max_mat.1 < threshold), arr.ind = TRUE))
  }
}

#' Plot comparative ordering statistics from two fitted models
#'
#' @param fit.1 A fitted hypercubic inference model including bootstrap resamples (output from hyperinf)
#' @param fit.2 Another fitted hypercubic inference model including bootstrap resamples  (output from hyperinf)
#' @param threshold Numeric (default 0). A threshold for differences. If zero, any separation of statistic (e.g. non-overlapping bootstrap distributions) is reported. If nonzero, we require one fit's statistic > 1-threshold and the other fit's statistic < threshold. 
#' @param ... additional parameters to pass to plot_hyperinf_bubbles
#' 
#' @return A ggplot object displaying pairs of features and timings for which separation is detected
#' @export
plot_hyperinf_compare_orderings = function(fit.1, fit.2, threshold=0.2, ...) {
  om.r = compare_orderings(fit.1, fit.2, threshold=threshold)
  om = compare_orderings(fit.1, fit.2, type="absolute", threshold=threshold)
  om.df = as.data.frame(om)
  om.r = om.r[om.r[,1] < om.r[,2],]
  if(nrow(om.r) >= 1) {
    om.r.df = data.frame(id=1:nrow(om.r), as.data.frame(om.r))
  } else {
    om.r.df = data.frame(id=0, row=0, col=0)
    om.r.df = om.r.df[NULL,]
  }
  b.plot = plot_hyperinf_bubbles(list(fit.1, fit.2), ...)
  b.plot + 
    ggplot2::geom_point(data = om.df, ggplot2::aes(x=col, y=row, fill=".", group=".")) +
    ggplot2::geom_segment(data = om.r.df, ggplot2::aes(x=-id/2, y=row, xend=-id/2, yend=col, fill=".", group="."))
}

