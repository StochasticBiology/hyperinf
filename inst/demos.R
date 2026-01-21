expt = "two.path"

set.seed(1)

if(expt == "one.path") {
  data = matrix(c(0,0,1,
                  0,1,1,
                  1,1,1), byrow = TRUE, ncol=3, nrow=3)
  # create random phylogeny
  tree = ape::rphylo(3, 1, 1)
} else if(expt == "two.path") {
  data = matrix(c(0,0,1,
                  0,1,1,
                  1,1,1,
                  1,0,0,
                  1,1,0,
                  1,1,1), byrow = TRUE, ncol=3, nrow=6)
  
  # create random phylogeny
  tree = ape::rphylo(6, 1, 1)
}

fit.mk = hyperinf(data, method="hypermk")
fit.hhmm = hyperinf(data, method="hyperhmm")
fit.dags = hyperinf(data, method="hyperdags")
fit.hct = hyperinf(data, method="hypertraps")
fit.pli = hyperinf(data, method="pli")

ggpubr::ggarrange(plot_hyperinf(fit.mk),plot_hyperinf(fit.hhmm),
                  plot_hyperinf(fit.dags), plot_hyperinf(fit.hct),
                  plot_hyperinf(fit.mk, plot.type="native"), plot_hyperinf(fit.hhmm, plot.type="native"),
                  plot_hyperinf(fit.dags, plot.type="native"), plot_hyperinf(fit.hct, plot.type="native"),
                  nrow=2, ncol=4)



fit.tree.mk = hyperinf(data, tree, method="hypermk")
fit.tree.hhmm = hyperinf(data, tree, method="hyperhmm")
fit.tree.dags = hyperinf(data, tree, method="hyperdags")
fit.tree.hct = hyperinf(data, tree, method="hypertraps")
ggpubr::ggarrange(plot_hyperinf(fit.tree.mk),plot_hyperinf(fit.tree.hhmm),
                  plot_hyperinf(fit.tree.dags), plot_hyperinf(fit.tree.hct),
                  plot_hyperinf(fit.tree.mk, plot.type="native"), plot_hyperinf(fit.tree.hhmm, plot.type="native"),
                  plot_hyperinf(fit.tree.dags, plot.type="native"), plot_hyperinf(fit.tree.hct, plot.type="native"),
                  nrow=2, ncol=4)

fit.tree.mk.rev = hyperinf(data, tree, method="hypermk", reversible = TRUE)
plot_hyperinf(fit.tree.mk)
plot_hyperinf(fit.tree.mk.rev)

data = matrix(rbinom(1000, 1, 0.5), ncol=50, nrow=20)
data = matrix(rbinom(9, 1, 0.5), ncol=3, nrow=3)
# create random phylogeny
tree = ape::rphylo(nrow(data), 1, 1)
fit.rnd = hyperinf(data, method="hyperdags")
plot_hyperinf(fit.rnd, plot.type="native")

data = matrix(c(0,0,1,
                0,2,1,
                1,1,1), ncol=3, nrow=3)

plot_hyperinf_data(data)
plot_hyperinf_data(data, tree)
test.missing = hyperinf(data)
test.missing.phylo = hyperinf(data, tree)
ggpubr::ggarrange(plot_hyperinf(test.missing), plot_hyperinf(test.missing.phylo))

