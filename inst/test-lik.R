library(ggpubr)
library(hyperinf)

try.mk = FALSE

# to do -- take a look at plugging HyperHMM fits into HyperTraPS -- are likelihoods preserved?
# can we just regularise HyperHMM internally?
# have hyperinf include native regularisation / pruning for HyperMk, HyperTraPS

obs = matrix(c(0,0,1, 1,0,0), nrow=2, ncol=3, byrow = TRUE)
obs = matrix(rbinom(60, 1, 0.5), nrow=10, ncol=6, byrow = TRUE)

fit.hm = hyperinf(obs)
fit.ht = hyperinf(obs, method="hypertraps", model=-1)
fit.ht2 = hyperinf(obs, method="hypertraps", model=2)
fit.htm = hyperinf(obs, method="hypertraps", model=-1, sa = 1)
fit.htp = hyperinf(obs, method="hypertraps", penalty = 1)
fit.hl = hyperinf(obs, method="hyperlau")
if(try.mk == TRUE) {
  fit.mk = hyperinf(obs, method="hypermk")
  fit.mki = hyperinf(obs, method="hypermk", reversible = FALSE)
}

hyperinf_AIC(fit.hm)
hyperinf_AIC(fit.hl)

hyperinf_AIC(fit.ht)
hyperinf_AIC(fit.ht2)
hyperinf_AIC(fit.htm)
hyperinf_AIC(fit.htp)

hyperinf_estimate_regularised(fit.hm)$post.aic
hyperinf_estimate_regularised(fit.hl)$post.aic

hyperinf_estimate_regularised(fit.ht)$post.aic
hyperinf_estimate_regularised(fit.ht2)$post.aic
hyperinf_estimate_regularised(fit.htm)$post.aic
hyperinf_estimate_regularised(fit.htp)$post.aic

hyperinf_estimate_regularised(fit.mk)$post.aic
hyperinf_estimate_regularised(fit.mki)$post.aic

if(try.mk == TRUE) {
  hyperinf_AIC(fit.mk)
  hyperinf_AIC(fit.mki)
  ggarrange(plot_hyperinf(fit.mk), plot_hyperinf(fit.mki),
            plot_hyperinf(fit.hm),
            plot_hyperinf(fit.ht), plot_hyperinf(fit.hl))
} else {
ggarrange(plot_hyperinf(fit.hm), plot_hyperinf(hyperinf_estimate_regularised(fit.hm)$regularised),
          plot_hyperinf(fit.ht), plot_hyperinf(fit.hl), plot_hyperinf(fit.htp))
}
