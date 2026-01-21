expt = 2
if(expt == 1) {
  m = matrix(rep(c(
    1,0,0,0,
    1,1,0,0,
    1,1,1,0), 5), byrow=TRUE, ncol=4)
} else if(expt == 2) {
  m = matrix(rep(c(0,0,0,1,
                   0,0,1,1,
                   0,1,1,1,
                   1,0,0,0,
                   1,1,0,0,
                   1,1,1,0), 5), byrow=TRUE, ncol=4)
} else {
  m = matrix(round(runif(80)), byrow=TRUE, ncol=4)
}
fitted.m = hyperinf::hyperinf(m)
fitted.ht = hyperinf::hyperinf(m, method="hypertraps")
est.fit = hyperinf::full_to_squared_fit(fitted.m)

# compare outputs
ggpubr::ggarrange(hyperinf::plot_hyperinf_data(m),
                  hyperinf::plot_hyperinf(fitted.m), 
                  hyperinf::plot_hyperinf(fitted.ht),
                  hyperinf::plot_hyperinf(est.fit),
                  nrow=2, ncol=2)
