library(hyperinf)
library(phytools)

# tests in the below. (c)ross-sectional data; (p)hylogenetic data; (u)ncertain data; (b)ootstrapping
# HyperHMM: c,p,cb,pb. Not allowed: u
# HyperLAU: c,p,cb,pb,cu,pu,cub,pub.
# HyperTraPS: c,p,cu. Not allowed: b,pu
# PLI: c,p,cu,pu. Not allowed: b. 
# HyperDAGs: c,p. Not allowed: b,u

# Fail: including HyperDAGs output in comparative plots
# To do: automatically plot bootstrap outcomes in comparative bubble plots

# note: this is just testing for errors and plausible-looking outputs, not correct ones

#### synthetic data

# construct a simple dataset
set.seed(1)
n = 30
if(TRUE) {
  data = matrix(rep(c(0,0,1, 0,1,1, 1,1,1), n/3), byrow = TRUE, ncol=3, nrow=n)
} else {
  data = matrix(rep(c(0,0,1, 0,1,1, 1,1,1, 1,0,0, 1,1,0, 1,1,1), n/6), byrow = TRUE, ncol=3, nrow=n)
}
plot_hyperinf_data(data)
# and some nuances on top of it:
# ... a phylogeny connecting our observations
my.tree = rphylo(n, birth=1, death=1)
df = data.frame(ID = my.tree$tip.label, as.data.frame(data))
plot_hyperinf_data(df, my.tree)
# ... cross-sectional uncertain data
u.data = data
u.data[sample(1:(3*n), 20)] = NA
plot_hyperinf_data(u.data)
# ... and make this uncertain data phylogenetically linked
u.df = data.frame(ID = my.tree$tip.label, as.data.frame(u.data))
plot_hyperinf_data(u.df, my.tree) 

##### cross-sectional data

# do model fits with bootstrapping to this dataset and its "inverse"
fit.1 = hyperinf(data, boot.parallel = 50)
fit.2 = hyperinf(1-data, boot.parallel = 50)

# plot the data
plot_hyperinf_data(data)
# plot a single model fit
plot_hyperinf(fit.1)
# compare transition networks
plot_hyperinf_comparative(list(fit.1, fit.2))
# compare summary "bubble" plots
plot_hyperinf_bubbles(list(fit.1, fit.2), p.scale = 0.2)
# compare bootstrapped bubble plots
plot_hyperinf_bootstrap(fit.1, fit.2)

# try some different approaches
fit.1 = hyperinf(data)
fit.2 = hyperinf(data, reversible=TRUE)
fit.3 = hyperinf(data, method="hypertraps")
fit.4 = hyperinf(data, method="hyperlau")
fit.5 = hyperinf(data, method="pli")
fit.6 = hyperinf(data, method="hyperdags")

# BUG: HyperDAGs doesn't yet fit in comparative plots

plot_hyperinf_comparative(list(fit.1, fit.2, fit.3, fit.4, fit.5), 
                          style="full", 
                          expt.names = c("HyperHMM", "HyperMk", "HyperTraPS", "HyperLAU", "PLI"))

plot_hyperinf_bubbles(list(fit.1, fit.2, fit.3, fit.4, fit.5), 
                      p.scale = 0.2,
                      expt.names = c("HyperHMM", "HyperMk", "HyperTraPS", "HyperLAU", "PLI"))

##### phylogenetic data

# try some different approaches
fit.1 = hyperinf(df, my.tree)
fit.2 = hyperinf(df, my.tree, reversible=TRUE)
fit.3 = hyperinf(df, my.tree, method="hypertraps")
fit.4 = hyperinf(df, my.tree, method="hyperlau")
fit.5 = hyperinf(df, my.tree, method="pli")
fit.6 = hyperinf(df, my.tree, method="hyperdags")

plot_hyperinf_comparative(list(fit.1, fit.2, fit.3, fit.4, fit.5), 
                          style="full", 
                          expt.names = c("HyperHMM", "HyperMk", "HyperTraPS", "HyperLAU", "PLI"))

plot_hyperinf_bubbles(list(fit.1, fit.2, fit.3, fit.4, fit.5), 
                      p.scale = 0.2,
                      expt.names = c("HyperHMM", "HyperMk", "HyperTraPS", "HyperLAU", "PLI"))

##### uncertain data

# try some different methods
fit.u.1 = hyperinf(u.data)
fit.u.2 = hyperinf(u.data, method="hypertraps")
fit.u.3 = hyperinf(u.data, method="pli")
plot_hyperinf_comparative(list(fit.u.1, fit.u.2, fit.u.3), 
                          style="full",
                          expt.names = c("HyperLAU", "HyperTraPS", "PLI"))
plot_hyperinf_bubbles(list(fit.u.1, fit.u.2, fit.u.3), 
                      expt.names = c("HyperLAU", "HyperTraPS", "PLI"))

# phylogenetically-linked uncertain data
fit.ut.1 = hyperinf(u.df, my.tree)
fit.ut.2 = hyperinf(u.df, my.tree, method="pli")
plot_hyperinf_comparative(list(fit.ut.1, fit.ut.2), 
                          style="full",
                          expt.names = c("HyperLAU", "PLI"))
plot_hyperinf_bubbles(list(fit.u.1, fit.u.2), 
                      expt.names = c("HyperLAU", "PLI"))

##### bootstrap in these cases
# for phylogenetically-linked data, we bootstrap over transition set
fit.1.boot = hyperinf(df, my.tree, boot.parallel = 10)
fit.4.boot = hyperinf(df, my.tree, method="hyperlau", boot.parallel = 10)

# for uncertain data
fit.u.1.boot = hyperinf(u.data, boot.parallel = 10)
# need to enforce independent transitions here, otherwise we don't have enough data
fit.ut.1.boot = hyperinf(u.df, my.tree, 
                         independent.transitions=FALSE,  
                         boot.parallel = 10) 

plot_hyperinf_comparative(list(fit.1.boot, fit.4.boot, fit.u.1.boot, fit.ut.1.boot),
                          style="full",
                          expt.names = c("Phylo, HHMM", "Phylo, HLAU", "CS Unc, HHMM", "Phylo Unc, HLAU"))

plot_hyperinf_bubbles(list(fit.1.boot, fit.4.boot, fit.u.1.boot, fit.ut.1.boot),
                          expt.names = c("Phylo, HHMM", "Phylo, HLAU", "CS Unc, HHMM", "Phylo Unc, HLAU"))
