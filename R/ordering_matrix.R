#' Compute a feature ordering matrix for a hypercubic inference model fit
#'
#' There are several options here. 
#' "relative" returns a matrix where P_ij is the proportion of *precursor states* encountered with i and without j. 
#' "transitions" returns a matrix where P_ij gives the proportion of feature i acquisitions in which feature j is already present. 
#' "absolute" returns a matrix where P_ij is the probability, across *precursor states*, that feature i is acquired after at least j features are present.
#'
#' There are several ways that models can differ. With "absolute", if there is a >(1-q) probability 
#' that a feature is acquired before t in model 1, and a <q probability that a feature is acquired 
#' before t in model 2, the dynamics differ with some threshold q.
#' 
#' With "relative", if there is a >(1-q) probability that feature i is acquired before feature j in
#' model 1, and a <q probability that feature i is acquired before feature j in model 2, the dynamics
#' differ with some threshold q.
#' 
#' For bootstrapped model fits, if these conditions should hold across bootstrap resamples, the difference
#' can be described as statistically significant at a level that depends on the number of cases for which
#' the differences are observed.
#' 
#' To be interesting, q should be small (q << 0.5 means genuinely different dynamics are being observed),
#' and this should hold across all bootstraps.
#'
#' @param fit A fitted hypercubic inference model (output from hyperinf)
#' @param n.samples Integer (default 10k) number of samples to characterise dynamics in the reversible (HyperMk) case
#' @param type Character string (default "relative"). 
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
  fit.type = hyperinf_gettype(fit)
  if(is.null(fit.type)) {
    return(ggplot2::ggplot())
  } 
  if(!(fit.type %in% c("hyperhmm", "hypertraps", "hypermk", "hypermk2", "hyperlau"))) {
    message("This fit type not yet supported!")
    stop()
  }
  ddf = NULL
  if(fit.type == "hypermk") {
    fit.rev = fit
    ddf = fit.rev$mk_fluxes
    count = count2 = 0
    abs.mat = matrix(0, nrow=fit.rev$L, ncol=fit.rev$L)
    b4m = rel.mat = matrix(0, nrow=fit.rev$L, ncol=fit.rev$L)
    for(i in 1:n.samples) {
      state = 0
      for(j in 1:fit.rev$L) {
        outs = which(ddf$From == state)
        if(length(outs) == 0) {
          state = 0
          outs = which(ddf$From == state)
        } 
        if(length(outs) == 1) {
          next.trans = outs 
        } else {
          next.trans = sample(outs, 1, prob = ddf$Rate[outs])
        }
        nextstate = ddf$To[next.trans]
        state.v = DecToBin(state, fit.rev$L)
        nextstate.v = DecToBin(nextstate, fit.rev$L)
        ones = which(nextstate.v==1)
        zeroes = which(nextstate.v==0)
        if(length(ones) > 0 & length(zeroes) > 0) {
          rel.mat[ones, zeroes] = rel.mat[ones, zeroes] + 1
          count2 = count2 + 1
        }
        level = sum(state.v)
        if(ddf$To[next.trans] > state) {
          change = fit.rev$L - log(ddf$To[next.trans]-state, base=2)
          #if(level == 5 & change == 4) {
          #  cat(state.v, " ", next.trans, "\n")
          #}
            abs.mat[change,1:(level+1)] = abs.mat[change,1:(level+1)]+1
        
          b4m[ones,change] = b4m[ones,change]+1
          count = count+1
        }
        state = ddf$To[next.trans]
      }
    }
    b4m = b4m/count
    rel.mat = rel.mat/count2
    rs <- apply(abs.mat, 1, max)
    abs.mat <- abs.mat / ifelse(rs == 0, 1, rs)
  } else if(fit.type == "hypermk2") {
    fit.rev = fit
    ddf = fit.rev$mk2_fluxes
    count = count2 = 0
    abs.mat = matrix(0, nrow=fit.rev$L, ncol=fit.rev$L)
    b4m = matrix(0, nrow=fit.rev$L, ncol=fit.rev$L)
    rel.mat = matrix(0, nrow=fit.rev$L, ncol=fit.rev$L)
    for(i in 1:n.samples) {
      state = ddf$From[sample(1:nrow(ddf), 1)]
      for(j in 1:(2*fit.rev$L)) {
        outs = which(ddf$From == state)
        if(length(outs) == 0) {
          state = min(ddf$From)
          outs = which(ddf$From == state)
        }
        if(length(outs) == 1) {
          next.trans = outs
        } else {
          next.trans = sample(outs, 1, prob = ddf$Rate[outs])
        }
        nextstate = ddf$To[next.trans]
        state.v = DecToBin(state, fit.rev$L)
        nextstate.v = DecToBin(nextstate, fit.rev$L)
        ones = which(nextstate.v==1)
        zeroes = which(nextstate.v==0)
        if(length(ones) > 0 & length(zeroes) > 0) {
          rel.mat[ones, zeroes] = rel.mat[ones, zeroes] + 1
          count2 = count2 + 1
        }
        level = sum(state.v)

        adds = which(state.v == 0 & nextstate.v == 1)
        if(length(adds) != 0) {
          abs.mat[adds,1:(level+1)] = abs.mat[adds,1:(level+1)]+1
          if(FALSE & length(adds) > 1) {
            start.l = min(level+1+1, fit.rev$L)
            end.l = min(level+1+length(adds), fit.rev$L)
            abs.mat[adds,start.l:end.l] = abs.mat[adds,start.l:end.l]+1/length(adds)
            count = count + 1/length(adds)
          }
          b4m[ones,adds] = b4m[ones,adds]+1
          count = count+1
        }
        state = ddf$To[next.trans]
      }
    }
    b4m = b4m/count
    rel.mat = rel.mat/count2
    rs <- apply(abs.mat, 1, max)
    abs.mat <- abs.mat / ifelse(rs == 0, 1, rs)
  } else {
    fit.non.rev = fit
    process.ddf = TRUE
    #count
    #lvls = lvls/sum(lvls)
    
    # b4m[i,j] gets incremented if we have i when we get j
    
    if(fit.type == "hyperhmm") {
      ddf = fit.non.rev$transitions
    } else if(fit.type == "hypertraps") {
      if("dynamics" %in% names(fit.non.rev)) {
        ddf = fit.non.rev$dynamics$trans
      } else if("routes" %in% names(fit.non.rev)) {
        message("Processing via sampled routes")
        routes_mat = fit.non.rev$routes
        abs.mat = matrix(0, nrow = ncol(routes_mat), ncol = ncol(routes_mat))
        rel.mat = matrix(0, nrow = ncol(routes_mat), ncol = ncol(routes_mat))
        for(i in 1:nrow(routes_mat)) {
          for(j in 1:ncol(routes_mat)) {
            ones = routes_mat[i,1:j]+1
            zeroes = (1:fit.non.rev$L)[-ones]
            rel.mat[ones, zeroes] = rel.mat[ones, zeroes] + 1
            abs.mat[routes_mat[i,j]+1, 1:j] = abs.mat[routes_mat[i,j]+1,1:j] + 1
          }
        }
        rel.mat = rel.mat / (ncol(routes_mat)*nrow(routes_mat))
        #m = m/nrow(routes_mat)
        if(type != "transitions") {
          process.ddf = FALSE
        }
      }
    } else if(fit.type == "hyperlau") {
      ddf = fit.non.rev$Dynamics
    }
    if(process.ddf == TRUE & is.null(ddf)) {
      message("Couldn't find transitions data frame, or alternative, in this model!")
      return(NULL)
    }

    if(process.ddf == TRUE) {
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
      count = count2 = 0
      abs.mat = matrix(0, nrow=fit.non.rev$L, ncol=fit.non.rev$L)
      b4m = rel.mat = matrix(0, nrow=fit.non.rev$L, ncol=fit.non.rev$L)
      for(i in 1:nrow(ddf)) {
        state = ddf$From[i]
        nextstate = ddf$To[i]
        
        state.v = DecToBin(state, fit.non.rev$L)
        nextstate.v = DecToBin(nextstate, fit.non.rev$L)
        
        change = fit.non.rev$L - log(nextstate-state, base=2)
        level = sum(state.v)
        ones = which(nextstate.v==1)
        zeroes = which(nextstate.v==0)
        if(length(ones) > 0 & length(zeroes) > 0) {
          rel.mat[ones, zeroes] = rel.mat[ones, zeroes] + ddf$Flux[i]
          count2 = count2 + ddf$Flux[i]
        }
        abs.mat[change,1:(level+1)] = abs.mat[change,1:(level+1)]+ddf$Flux[i]
        b4m[ones,change] =  b4m[ones,change] + ddf$Flux[i]
      }
      rel.mat = rel.mat / count2
    }
    
    rs <- apply(abs.mat, 1, max)
    abs.mat <- abs.mat / ifelse(rs == 0, 1, rs)
    
  }
  if(type == "relative") {
    return(rel.mat)
  } else if(type == "transitions") {
    return(b4m)
  } else {
    return(abs.mat)
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
#' Comparisons of ordering statistics look at two components. First (a), is there a 
#' lack of overlap between bootstrap distributions for a given summary statistic?
#' Second (b), is the different in magnitudes in this summary statistic large? If (a)
#' but not (b), a precise statistic differs but the two models probably produce
#' similar dynamics. If (b) but not (a), the two models may produce different "mean"
#' dynamics, but there's enough noise that they overlap. If both (a) and (b), there
#' is a precisely specified, large difference in the dynamics.
#' 
#' The parameter "percentile" controls (a); "threshold" controls (b).
#'
#' @param fit.1 A fitted hypercubic inference model including bootstrap resamples (output from hyperinf)
#' @param fit.2 Another fitted hypercubic inference model including bootstrap resamples  (output from hyperinf)
#' @param n.samples Integer (default 10k) number of samples to characterise dynamics in the reversible (HyperMk) case
#' @param threshold Numeric or NULL (default 0.25). A threshold for differences. If NULL, any separation of statistic (e.g. non-overlapping bootstrap distributions) is reported. If nonzero, we require one fit's statistic > 1-threshold and the other fit's statistic < threshold. 
#' @param percentile Numeric (default 0.2). The percentile at which bootstrapped distributions should show no overlap.
#' @param type Character string (default "relative"). Either "relative", in which case P_ij gives the probability that feature i is acquired before feature j. Or "absolute", in which case P_ij is the probability that feature i is acquired after step j.
#' 
#' @return A hypercubic inference model structure with the model fits stored as "bootstrap resamples"
#' @export
compare_orderings = function(fit.1, fit.2, 
                             n.samples = 10000, 
                             threshold = 0.25,
                             percentile = 0.2,
                             type = "relative") {
  oms.1 = oms.2 = list()
  if(!("boots" %in% names(fit.1) & "boots" %in% names(fit.2))) {
    message("Didn't find bootstrap resamples in these fits")
    return(NULL)
  }
  if(is.null(threshold)) {
    message("Warning: comparing distributions without considering probabilities may not give a scientifically interesting result.")
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

  # pull percentiles across bootstraps (3rd dimension)
  min_mat.1 <- apply(arr.1, c(1, 2), quantile, probs=percentile) 
  max_mat.1 <- apply(arr.1, c(1, 2), quantile, probs=1-percentile) 
  min_mat.2 <- apply(arr.2, c(1, 2), quantile, probs=percentile) 
  max_mat.2 <- apply(arr.2, c(1, 2), quantile, probs=1-percentile) 
  
  if(is.null(threshold)) {
    win.1 = which((min_mat.1 > max_mat.2) , arr.ind = TRUE)
    win.2 = which((min_mat.2 > max_mat.1) , arr.ind = TRUE)
    wins = data.frame()
    if(length(win.1) > 0) {
      wins = rbind(wins, data.frame(first=1, win.1))
    } 
    if(length(win.2) > 0) {
      wins = rbind(wins, data.frame(first=2, win.2))
    }
    
    return(wins)
  } else {
    return(which((min_mat.1 > 1-threshold & max_mat.2 < threshold) 
                 | (min_mat.2 > 1-threshold & max_mat.1 < threshold), arr.ind = TRUE))
  }
}

#' Plot comparative ordering statistics from two fitted models
#'
#' See "compare_orderings" for a description of the "threshold" and "percentile"
#' parameters, which broadly control the size and the "significance" of reported
#' differences.
#'
#' @param fit.1 A fitted hypercubic inference model including bootstrap resamples (output from hyperinf)
#' @param fit.2 Another fitted hypercubic inference model including bootstrap resamples  (output from hyperinf)
#' @param threshold Numeric or NULL (default 0.25). A threshold for differences. If zero, any separation of statistic (e.g. non-overlapping bootstrap distributions) is reported. If nonzero, we require one fit's statistic > 1-threshold and the other fit's statistic < threshold. 
#' @param percentile Numeric (default 0.2). The percentile at which bootstrapped distributions should show no overlap.
#' @param ... additional parameters to pass to plot_hyperinf_bubbles
#' 
#' @return A ggplot object displaying pairs of features and timings for which separation is detected
#' @export
plot_hyperinf_compare_orderings = function(fit.1, fit.2, 
                                           threshold=0.25, 
                                           percentile=0.2,
                                           ...) {
  om.r = compare_orderings(fit.1, fit.2, threshold=threshold, percentile=percentile)
  om = compare_orderings(fit.1, fit.2, type="absolute", threshold=threshold, percentile=percentile)
  om.df = as.data.frame(om)
  om.r = om.r[om.r[,1] < om.r[,2], , drop=FALSE]
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


#' Plot ordering matrices from fitted models
#'
#' @param fits A list of fitted hypercubic inference models (or single model)
#' @param type Character string (default "absolute"): which type of ordering matrix to plot ("absolute", "relative"; see ordering_matrix)
#' @param expt.names Optional vector of labels for each element of the fit list
#' @param feature.names Boolean or character vector (default TRUE). If TRUE, use feature names from fit. If FALSE, use numerical labels. If vector of length L, use those labels.
#' @param thetastep Integer (default 5), the number of angular steps to take when drawing polygons for each bubble segment. Lower this for many comparisons.
#' @param p.scale Numeric (default 1), the scaling factor to apply to the radius (probability) of bubbles
#' @param sqrt.trans Boolean (default FALSE), whether to sqrt-transform the probability to get the bubble radius
#' 
#' @return A ggplot object comparing ordering matrices
#' @export
plot_hyperinf_ordering_matrices = function(fits,
                                           type = "absolute",
                                           expt.names = NULL,
                                           feature.names = TRUE,
                                           thetastep = 5,
                                           p.scale = 0.5,
                                           sqrt.trans = FALSE) {
  if("L" %in% names(fits)) {
    fits = list(fits)
  }
  if(!(type %in% c("absolute", "relative"))) {
    message("Only absolute and relative types supported")
    return(ggplot2::ggplot())
  }
  # convert matrix to dataframe
  mat_to_df <- function(mat, type) {
    expand.grid(
      y = seq_len(ncol(mat)),
      x = seq_len(nrow(mat))
    ) |>
      transform(
        value = as.vector(mat),
        type = type
      )
  }
  
  m = list()
  df = data.frame()
  for(i in 1:length(fits)) {
    m[[i]] = ordering_matrix(fits[[i]], type=type) 
   
    tmp = mat_to_df(m[[i]], i)
    df = rbind(df, tmp)
  }
  
  # semicircle generator
  semicircle <- function(x, y, r, index, nindex, n) {
    
    theta <- seq(2*pi*(index/nindex), 2*pi*((index+1)/nindex), length.out = thetastep+1)
    tmp = data.frame(x=x,y=y) 
    tmp = rbind(tmp, data.frame(
      x = x + r*cos(theta),
      y = y + r*sin(theta)
    ))
    tmp = rbind(tmp, data.frame(x=x,y=y))
  }
  
  # build polygon data
  poly_list <- lapply(seq_len(nrow(df)), function(i){
    
    row <- df[i,]
    if(sqrt.trans == TRUE) {
      val = sqrt(row$value)*p.scale
    } else {
      val = row$value*p.scale
    }
    sc <- semicircle(
      row$x,
      row$y,
      val,
      index = row$type,
      nindex = length(fits)
    )
    
    sc$value <- row$value
    sc$id <- i
    if(length(expt.names) == 0) {
      sc$type = row$type 
    } else {
      sc$type = expt.names[row$type]
    }
    sc
  })
  
  poly_df <- do.call(rbind, poly_list)
  if(type == "absolute") {
    poly_df$x = poly_df$x - 1
  }
  
  # plot
  this.plot = ggplot2::ggplot(poly_df, ggplot2::aes(x, y, group=id, fill=factor(type))) +
    ggplot2::geom_polygon() +
    #ggplot2::scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0.5) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() 
  
  if(type == "relative") {
    this.plot = this.plot + ggplot2::labs(fill = "Experiment", x = "Feature absent", y = "Feature present")
  } else {
    this.plot = this.plot + ggplot2::labs(fill = "Experiment", x = "At least n prior acquisitions?", y = "New acquisition")
  }
  
  yset = xset = FALSE
  if(length(feature.names) != 0) {
    if(length(feature.names) == 1) {
      if(feature.names == TRUE & !is.null(fits[[1]]$feature.names)) {
        this.plot = this.plot + ggplot2::scale_y_continuous(breaks=1:fits[[1]]$L, labels = fits[[1]]$feature.names)
        yset = TRUE
        if(type == "relative") {
          this.plot = this.plot + ggplot2::scale_x_continuous(breaks=1:fits[[1]]$L, labels = fits[[1]]$feature.names)
          xset = TRUE
        }
      }
    } else if(length(feature.names) == fits[[1]]$L) {
      this.plot = this.plot + ggplot2::scale_y_continuous(breaks=1:fits[[1]]$L, labels = feature.names)
      yset = TRUE
      if(type == "relative") {
        this.plot = this.plot + ggplot2::scale_x_continuous(breaks=1:fits[[1]]$L, labels = feature.names)
        xset = TRUE
      }
    }
  }
  if(xset == FALSE) {
    if(type == "relative") {
      this.plot = this.plot + ggplot2::scale_x_continuous(breaks=1:fits[[1]]$L)
    } else {
      this.plot = this.plot + ggplot2::scale_x_continuous(breaks=0:(fits[[1]]$L-1))
    }
  }
  if(yset == FALSE) {
    this.plot = this.plot + ggplot2::scale_y_continuous(breaks=1:fits[[1]]$L)
  }
  return(this.plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)))
}