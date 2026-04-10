#' Run diagnostics on a model fit
#'
#' @param fit A fitted hypercubic inference model (output from hyperinf)
#' @param show.plot Boolean (default TRUE) whether to plot a diagnostic plot if appropriate
#'
#' @return A dataframe with descriptive statistics depending on fit type. Messages are printed to console if warning signs exist for model convergence.
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit.1 = hyperinf(data, method="hypertraps")
#' hyperinf_diagnostics(fit.1)
#' @export
hyperinf_diagnostics = function(fit, 
                                show.plot = TRUE) {
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
  }
  if(!(fit.type %in% c("hypertraps"))) {
    message("This fit type doesn't support diagnostics")
    return(NULL)
  }
  df = this.fit$lik.traces
  sds = apply(df[,5:7], 1, sd)
  cvs = abs(sds/rowMeans(df[,5:7]))
  kpss = tseries::kpss.test(df$CurrentLogLikelihood)
  static = length(which(duplicated(df$CurrentLogLikelihood)))
  if(length(static) == 0) { static = 0 }
  report = data.frame(max.cv = max(cvs),
                      mean.cv = mean(cvs),
                      kpss.p.val = kpss$p.value,
                      n.static = static,
                      n = length(cvs))
  problem = FALSE
  if(report$max.cv > 1e-2) {
    problem = TRUE
    message("WARNING: > ", round(report$max.cv*100), "% max discrepancy between likelihood estimates")
  } 
  if(report$mean.cv > 1e-3) {
    problem = TRUE
    message("WARNING: > ", round(report$max.cv*100, digits=1), "% mean discrepancy between likelihood estimates")
  }
  if(report$n.static > 1) {
    problem = TRUE
    message("WARNING: > ", round(100*report$n.static/report$n, digits=3), "% samples were stuck")
  }
  if(report$kpss.p.val < 0.05) {
    problem = TRUE
    message("WARNING: stationary p-value = ", report$kpss.p.val)
  }
  if(report$n < 10) {
    problem = TRUE
    message("WARNING: only ", report$n, " samples -- did the chain run long enough?")
  }
  if(problem == FALSE) {
    message("No immediate issues found.")
  }
  if(show.plot == TRUE) {
    print(hypertrapsct::plotHypercube.lik.trace(fit))
  }
  return(report)
}
