#' Get the likelihood for a fitted hypercubic inference model
#'
#' @param fit A fitted hypercubic inference model (output from hyperinf)
#' @param emission Boolean (default FALSE), whether to consider emission probabilities P(s->t)Pemit(s, t) or just P(s->t), the probability that a random walker passes through s and t
#' @param expanded Boolean (default FALSE) for HyperTraPS, whether to report the range of sampled likelihoods or just the max of this sample
#'
#' @return A numerical value or (for expanded HyperTraPS) a vector of values
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit = hyperinf(data)
#' hyperinf_loglikelihood(fit)
#' @export
hyperinf_loglikelihood = function(fit, 
                                  emission = FALSE,
                                  expanded = FALSE) {
  L = fit$L
  fit.type = hyperinf_gettype(fit)
  if(fit.type == "DAG") {
    message("DAG models don't have associated likelihoods")
    return(NULL);
  } else if(fit.type == "arborescence") {
    message("DAG models don't have associated likelihoods")
    return(NULL);
  } else if(fit.type == "hypertraps") {
    likset = fit$lik.traces$CurrentLogLikelihood
    paramset = fit$lik.traces$nparam
    if(sd(paramset) > 0) {
      message("This inference was run with penalised likelihood; parameter counts change. Be aware that this likelihood will be a penalised one!")
    }
    loglik = fit$bestlik
    if(!("data" %in% names(fit))) {
      message("This model fit wasn't run with hyperinf, so I don't have the original data. All I can give you is the original, uncorrected HyperTraPS best-sampled loglik") 
      return(max(fit$lik.traces$CurrentLogLikelihood))
    }
    if("initialstates" %in% names(fit$data)) {
      diffs = rowSums(fit$data$obs) - rowSums(fit$data$initialstates)
    } else {
      diffs = rowSums(fit$data$obs)
    }
    emission.correction = sum(log(1/diffs))
    loglik = loglik - emission.correction 
    if(expanded == FALSE) {
      message("Reporting maximum sampled likelihood")
    } else {
      return(likset)
    }
  } else if(fit.type == "hyperlau") {
    loglik = max(fit$Likelihood$Likelihood[fit$Likelihood$Bootstrap == 0])
    if(!("data" %in% names(fit))) {
      message("This model fit wasn't run with hyperinf, so I don't have the original data. All I can give you is the original, uncorrected HyperLAU best-sampled loglik") 
      return(loglik)
    }
    logPemit = nrow(fit$data$obs) * log(2/((fit$L+1)*(fit$L+2)))
    loglik = loglik - logPemit
  } else if(fit.type == "hyperhmm") {
    loglik = fit$loglik
  } else if(fit.type == "hypermk" | fit.type == "hypermk2") {
    loglik = fit$fitted_mk$loglikelihood
  } else {
    return(NULL);
  }
  if(emission == TRUE) {
    if(!("data" %in% names(fit))) {
      message("This model fit wasn't run with hyperinf, so I don't have the original data. All I can give you is the original, uncorrected best likelihood") 
      return(loglik)
    }
    if("initialstates" %in% fit$data) {
      logPemit = nrow(fit$data$obs) * log((1/(fit$L+1))*(1/fit$L))
    } else {
      logPemit = nrow(fit$data$obs) * log(1/(fit$L+1))
    }
    if(fit.type %in% c("hyperhmm", "hypertraps", "hyperlau")) {
      loglik = loglik + logPemit
    }
  }
  return(loglik)
}

#' Get the AIC for a fitted hypercubic inference model
#'
#' @param fit A fitted hypercubic inference model (output from hyperinf)
#' @param ... other parameters to pass to hyperinf_loglikelihood
#'
#' @return A numerical value
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit = hyperinf(data)
#' hyperinf_AIC(fit)
#' @export
hyperinf_AIC = function(fit, ...) {
  L = fit$L
  fit.type = hyperinf_gettype(fit)
  if(fit.type == "DAG") {
    message("DAG models don't have associated likelihoods")
    return(NULL);
  } else if(fit.type == "arborescence") {
    message("DAG models don't have associated likelihoods")
    return(NULL);
  } else if(fit.type == "hyperlau") {
    if(!("model" %in% names(fit))) {
      message("This model fit wasn't run with hyperinf, so I can't be sure of the model structure. I'm assuming all-transitions-independent.")
      fit$model = -1
    }
  } else if(fit.type == "hyperhmm") {
    fit$model = -1
  } else if(fit.type == "hypermk" | fit.type == "hypermk2") {
    return(data.frame(loglik = fit$fitted_mk$loglikelihood,
                      nparam = (fit$fitted_mk$AIC + 2*fit$fitted_mk$loglikelihood)/2,
                      AIC = fit$fitted_mk$AIC))
  } else {
    return(NULL);
  }
  loglik = hyperinf_loglikelihood(fit, ...)
  if(fit.type %in% c("hyperhmm", "hypertraps", "hyperlau")) {
    if(fit$model == -1) {
      nparam = fit$L * (2**(fit$L-1))
    } else if(fit$model == 1) {
      nparam = fit$L
    } else if(fit$model == 2) {
      nparam = (fit$L**2)
    } else if(fit$model == 3) {
      nparam = (fit$L**2)*(fit$L+1)/2
    } else if(fit$model == 4) {
      nparam = (fit$L**2)*(fit$L+1)*(fit$L+2)/6
    }
    if(fit.type == "hypertraps" & "bestparams" %in% names(fit)) {
      if(fit$bestparams < nparam) {
        nparam = fit$bestparams
      }
    }
  }
  AIC = 2*nparam - 2*loglik
  return(data.frame(loglik = loglik, nparam = nparam, AIC = AIC))
}

#' Crude regularisation estimate of a hypercubic inference model
#'
#' @param fit A fitted hypercubic inference model (output from hyperinf)
#' @param threshold Numeric (default 1e-3), threshold of probability flux below which to remove a transition
#'
#' @return A list of regularised model and AIC statistics before and after "regularisation"
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit = hyperinf(data, method="hyperhmm")
#' hyperinf_estimate_regularised(fit)
#' @export
hyperinf_estimate_regularised = function(fit, threshold = 1e-3) {
  fit.type = hyperinf_gettype(fit)
  if(is.null(fit.type)) {
    return(NULL)
  } 

  pre.aic = hyperinf_AIC(fit)
  if(fit.type == "hyperhmm") {
    tdf = fit$transitions
    zeroes = which(tdf$Flux < threshold)
    fit$transitions$Flux[zeroes] = 0
  } 
  if(fit.type == "hyperlau") {
    if(fit$model == -1) {
      tdf = fit$Dynamics
      zeroes = which(tdf$Flux < threshold)
      fit$Dynamics$Flux[zeroes] = 0
    } else {
      message("Couldn't establish that this is an all-transitions-independent model!")
      return(NULL)
    }
  }
  if(fit.type == "hypermk") {
    tdf = fit$mk_fluxes
    zeroes = which(tdf$Flux < threshold*sum(tdf$Flux[tdf$From == 0]))
    fit$mk_fluxes$Flux[zeroes] = 0
  }
  if(fit.type == "hypertraps") {
    if(fit$model == -1) {
    tdf = fit$dynamics$trans
    zeroes = which(tdf$Flux < threshold)
    fit$dynamics$trans$Flux[zeroes] = 0
    } else {
      message("Couldn't establish that this is an all-transitions-independent model!")
      return(NULL)
    }
  }
  post.aic = pre.aic
  post.aic$nparam = post.aic$nparam - length(zeroes)
  post.aic$AIC = post.aic$AIC - 2*length(zeroes)
  return(list(regularised = fit,
              pre.aic = pre.aic,
              post.aic = post.aic))
}

#' Regularisation of a hypercubic inference model
#'
#' @param fit A fitted hypercubic inference model (output from hyperinf)
#' @param threshold Numeric (default 0), threshold of probability flux below which to remove a transition
#'
#' @return A list of regularised model and AIC statistics before and after "regularisation"
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit = hyperinf(data, method="hypermk")
#' hyperinf_regularise(fit)
#' @export
hyperinf_regularise = function(fit, threshold = 0) {
  fit.type = hyperinf_gettype(fit)
  if(is.null(fit.type)) {
    return(NULL)
  } 
  
  if(!(fit.type %in% c("hypermk"))) {
    if(threshold == 0) {
      message("Can't formally regularise this type yet. Running estimated regularisation.")
      return(hyperinf_estimate_regularised(fit))
    } else {
      message("Can't regularise this type yet!")
      return(NULL)
    }
  }
  pre.aic = hyperinf_AIC(fit)
  # AIC = 2k - 2 loglik
  nparams = (fit$fitted_mk$AIC + 2*fit$fitted_mk$loglikelihood)/2
  if(round(nparams) <= 2*fit$L) {
    message("Looks like this is already a first-order model!")
    post.aic = pre.aic
    new.fit = fit
    
    return(list(regularised = new.fit,
                pre.aic = pre.aic,
                post.aic = post.aic))
  } 

  new.fit = hypermk::mk_prune_model(fit, flux.threshold = threshold)
  new.fit$feature.names = fit$feature.names
  post.aic = hyperinf_AIC(new.fit)

  return(list(regularised = new.fit,
              pre.aic = pre.aic,
              post.aic = post.aic))
}

