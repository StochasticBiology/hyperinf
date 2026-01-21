#' Converts a fully parameterised model to a set of parameters for L^2 HyperTraPS
#'
#' @param fit A fitted model with a fully parameterised hypercube
#'
#' @return A vector of L^2 values corresponding to the parameterisation of a HyperTraPS model
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit = hyperinf(data)
#' l2params = full_to_squared(fit)
#' @export
full_to_squared = function(fit) {
  if("best.graph" %in% names(fit)) {
    fit.type = "DAG"
    message("HyperDAGs not yet supported")
    return(NA)
  } else if("raw.graph" %in% names(fit)) {
    fit.type = "arborescence"
    message("HyperDAGs not yet supported")
    return(NA)
  } else if("posterior.samples" %in% names(fit)) {
    fit.type = "hypertraps"
    df = fit$dynamics$trans
  } else if("Dynamics" %in% names(fit)) {
    fit.type = "hyperlau"
    df = fit$Dynamics
  } else if("viz" %in% names(fit)) {
    fit.type = "hyperhmm"
    df = fit$transitions
  } else if("fitted_mk" %in% names(fit)) {
    fit.type = "mk"
    df = fit$mk_fluxes
  } else {
    message("Didn't recognise this model")
    return(NA)
  }
  
  L = fit$L
  
  # covariates.mat stores presence of other features for each transition
  covariates.mat = replicate(
    L,
    matrix(numeric(0), nrow = 0, ncol = L),
    simplify = FALSE
  )
  # weights and response store those things for each transition
  weights = response = replicate(
    L,
    matrix(numeric(0), nrow = 0, ncol = 1),
    simplify = FALSE
  )
  
  # loop through rows in transition set
  for(i in 1:nrow(df)) {
    # characterise this change
    src = DecToBin(df$From[i], L)
    dest = DecToBin(df$To[i], L)
    change = which(dest != src)
    # which features are present?
    present = which(src != 0)
    coeffs = rep(0, L)
    coeffs[present] = 1
    # covariates.mat[changed feature] gets a new row containing these present features
    covariates.mat[[change]] = rbind(covariates.mat[[change]], coeffs)
    # response[change] gets the responseonse probability
    if(fit.type %in% c("hyperlau", "hyperhmm", "hypertraps")) {
      response[[change]] = rbind(response[[change]], df$Probability[i])
      weights[[change]] = rbind(weights[[change]], df$Flux[i])
    } else if(fit.type == "mk") {
      response[[change]] = rbind(response[[change]], df$Rate[i])
      weights[[change]] = rbind(weights[[change]], (1+df$Flux[i])/max(df$Flux))
    }
    
  }
  
  # l2 will store our coefficient estimates
  l2 = matrix(0, nrow=L, ncol=L)
  # for each feature
  for(i in 1:L) {
    # remove the feature itself from the covariate set
    this.covariates.mat = covariates.mat[[i]][,-i]
    # attempt to weighted-fit exp( beta_0 + sum beta_i covariate_i ) = response
    this.lm = lm(log(response[[i]]) ~ this.covariates.mat, weights=weights[[i]])
    # store modifiers and base rate
    l2[i,-i] = this.lm$coefficients[2:L]
    l2[i,i] = this.lm$coefficients[1]
  }
  # remove NAs and cast in form for HyperTraPS analysis
  l2[is.na(l2)] = 0
  l2v = as.vector(l2)
  l2tmp = as.matrix(l2v)
  return(l2tmp)
}

#' Produces a fitted HyperTraPS model from a parameterisation
#'
#' @param params A vector of L^2 values corresponding to a HyperTraPS parameterisation
#' @param reps An integer (default 100) of times to repeat this parameterisation to build up a model structure
#'
#' @return A HyperTraPS model structure
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit = hyperinf(data)
#' l2params = full_to_squared(fit)
#' fit.ht = hypertraps_from_params(l2params)
#' @export
hypertraps_from_params = function(params,
                                  reps = 100) {
  L = sqrt(length(params))
  l2m <- t(cbind(params, params[, 1, drop = FALSE][, rep(1, reps)]))
  est = list(posterior.samples=l2m,
             model=2,
             L=L)
  est.fit = hypertrapsct::PosteriorAnalysis(est)
  return(est.fit)
}

#' Estimates an L^2 HyperTraPS model from fully parameterised hypercubic model
#'
#' Combines full_to_squared and hypertraps_from_params
#' 
#' @param fit  A fitted model with a fully parameterised hypercube
#' @param reps An integer (default 100) of times to repeat this parameterisation to build up a model structure
#'
#' @return A HyperTraPS model structure
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit = hyperinf(data)
#' est.ht = full_to_squared_fit(fit)
#' @export
full_to_squared_fit = function(fit,
                               reps = 100) {
  params = full_to_squared(fit)
  est.fit = hypertraps_from_params(params, reps)
  return(est.fit)
}

