set.seed(4)
my.tree = ape::rtree(4)
my.df = data.frame(label = my.tree$tip.label,
                   obs = c("001", "0??", "?10", "??1"))

c.utree = hyperlau::curate.uncertain.tree(my.tree, my.df)
fit = hyperlau::HyperLAU(c.utree$dests, c.utree$srcs)
fit$Dynamics

hyperinf::plot_hyperinf(fit)

# and we can hopefully run HyperLAU natively from hyperinf, including parallelised bootstrap uncertainty
fit.inf = hyperinf::hyperinf(my.df, my.tree, boot.parallel = 5)
hyperinf::plot_hyperinf(fit.inf)

# natively in hyperinf, with and without guaranteeing transition independence
fit.1 = hyperinf::hyperinf(my.df, my.tree)
fit.2 = hyperinf::hyperinf(my.df, my.tree, independent.transitions = FALSE)

ggpubr::ggarrange(hyperinf::plot_hyperinf_data(my.df, my.tree),
                  hyperinf::plot_hyperinf(fit.1),
                  hyperinf::plot_hyperinf(fit.2),
                  nrow=1)

# uncertainty and visualisation via the (serial) bootstrap
fit.1 = hyperinf::hyperinf(my.df, nboot = 5)
ggpubr::ggarrange(hyperinf::plot_hyperinf_data(my.df),
                  hyperinf::plot_hyperinf(fit.1, uncertainty = FALSE),
                  hyperinf::plot_hyperinf(fit.1, uncertainty = TRUE),
                  nrow=1)

# uncertainty and visualisation via the (serial) bootstrap with HyperHMM
m1 = matrix(rep(c(0,0,1,
                  0,1,1,
                  1,1,1), 1), byrow=TRUE, ncol=3)
m2 = matrix(rep(c(0,0,1,
                  0,1,1,
                  1,1,1), 5), byrow=TRUE, ncol=3)
fit.hmm1 = hyperinf::hyperinf(m1, nboot = 5)
fit.hmm2 = hyperinf::hyperinf(m2, nboot = 5)
ggpubr::ggarrange(hyperinf::plot_hyperinf_data(m2),
                  hyperinf::plot_hyperinf(fit.hmm1, uncertainty = TRUE),
                  hyperinf::plot_hyperinf(fit.hmm2, uncertainty = TRUE),
                  nrow=1)

# parallelised bootstrap in HyperHMM
fit.hmm.p = hyperinf::hyperinf(m1, boot.parallel = 50)
hyperinf::plot_hyperinf(fit.hmm.p)
