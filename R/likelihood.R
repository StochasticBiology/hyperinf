#' Get the likelihood for a fitted hypercubic inference model
#'
#' @param fit A fitted hypercubic inference model (output from hyperinf)
#' @param emission Boolean (default TRUE), whether to consider emission probabilities P(s->t)Pemit(s, t) or just P(s->t), the probability that a random walker passes through s and t
#' @param expanded Boolean (default FALSE) for HyperTraPS, whether to report the range of sampled likelihoods or just the max of this sample
#'
#' @return A numerical value or (for expanded HyperTraPS) a vector of values
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit = hyperinf(data)
#' hyperinf_loglikelihood(fit)
#' @export
hyperinf_loglikelihood = function(fit, 
                                  emission = TRUE,
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
  if("initialstates" %in% fit$data) {
    logPemit = nrow(fit$data$obs) * log((1/(fit$L+1))*(1/fit$L))
  } else {
    logPemit = nrow(fit$data$obs) * log(1/(fit$L+1))
  }
  if(emission == TRUE) {
    if(fit.type %in% c("hyperhmm", "hypertraps", "hyperlau")) {
      loglik = loglik + logPemit
    }
  }
  return(loglik)
}