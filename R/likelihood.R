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
  if("best.graph" %in% names(fit)) {
    fit.type = "DAG"
    message("DAG models don't have associated likelihoods")
    return(NULL);
  } else if("raw.graph" %in% names(fit)) {
    fit.type = "arborescence"
    message("DAG models don't have associated likelihoods")
    return(NULL);
  } else if("posterior.samples" %in% names(fit)) {
    fit.type = "hypertraps"
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
  } else if("Dynamics" %in% names(fit)) {
    fit.type = "hyperlau"
    loglik = max(fit$Likelihood$Likelihood[fit$Likelihood$Bootstrap == 0])
    if(!("data" %in% names(fit))) {
      message("This model fit wasn't run with hyperinf, so I don't have the original data. All I can give you is the original, uncorrected HyperLAU best-sampled loglik") 
      return(loglik)
    }
    logPemit = nrow(fit$data$obs) * log(2/((fit$L+1)*(fit$L+2)))
    loglik = loglik - logPemit
  } else if("viz" %in% names(fit)) {
    fit.type = "hyperhmm"
    loglik = fit$loglik
  } else if("fitted_mk" %in% names(fit)) {
    fit.type = "mk"
    loglik = fit$fitted_mk$loglikelihood
  } else {
    message("Didn't recognise this model")
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
  if("best.graph" %in% names(fit)) {
    fit.type = "DAG"
    message("DAG models don't have associated likelihoods")
    return(NULL);
  } else if("raw.graph" %in% names(fit)) {
    fit.type = "arborescence"
    message("DAG models don't have associated likelihoods")
    return(NULL);
  } else if("posterior.samples" %in% names(fit)) {
    fit.type = "hypertraps"
  } else if("Dynamics" %in% names(fit)) {
    fit.type = "hyperlau"
    if(!("model" %in% names(fit))) {
      message("This model fit wasn't run with hyperinf, so I can't be sure of the model structure. I'm assuming all-transitions-independent.")
      fit$model = -1
    }
  } else if("viz" %in% names(fit)) {
    fit.type = "hyperhmm"
    fit$model = -1
  } else if("fitted_mk" %in% names(fit)) {
    fit.type = "mk"
    return(fit.mk$fitted_mk$AIC)
  } else {
    message("Didn't recognise this model")
    return(NULL);
  }
  loglik = hyperinf_loglikelihood(fit, ...)
  if(fit.type == "hypertraps" & "bestparams" %in% names(fit)) {
     nparam = fit$bestparams
  } else if(fit.type %in% c("hyperhmm", "hypertraps", "hyperlau")) {
    if(fit.type == "hypertraps") {
      message("This fit used an old codebase, so I can't see the best parameterisation directly. I'll use the original model structure.")
    }
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
  }
  AIC = 2*nparam - 2*loglik
  return(AIC)
}