# hyperinf

An Ecosystem for Hypercubic Inference in Evolutionary Accumulation Modelling
---

This package brings together HyperTraPS(-CT), HyperHMM, HyperLAU, HyperDAGs and HyperMk(2), different algorithms for accumulation modelling. All are based around the idea of inferring hypercubic transition graphs that describe the accumulation of binary features.

Install with

`remotes::install_github("StochasticBiology/hyperinf")`

This includes content from:
* HyperTraPS [1] https://github.com/StochasticBiology/hypertrapsct
* HyperMk [2] https://github.com/StochasticBiology/hypermk
* HyperHMM [3] https://github.com/StochasticBiology/hyperhmm
* HyperDAGs [4] https://github.com/StochasticBiology/hyperdags
* HyperLAU [5] https://github.com/StochasticBiology/hyperlau
* ideas (but not code) from Cluster-HyperHMM [6]
* HyperMk2 https://github.com/StochasticBiology/hypermk2

Wrapper function `hyperinf` produces fitted hypercubic inference models from data; `plot_hyperinf` produces summary plots of the inferred transition graph. Arguments to `hyperinf` specify which fitting approach to use (default is to choose based on data structure). `full_to_squared_fit` converts a fully-parameterised hypercube (including those output from HyperHMM, HyperMK, HyperLAU) to a best estimate of an L<sup>2</sup>-parameterised HyperTraPS model (encoding pairwise interactions between features). `plot_hyperinf_comparative` and `plot_hyperinf_bubbles` produce comparative plots summarising dynamics across a list of model fits; `plot_hyperinf_bootstrap` compares bubble plots from two bootstrapped model fits. Functionality in `ordering_matrix`, `compare_orderings`, and `plot_hyperinf_compare_orderings` supports comparisons between matrices describing probabilities that feature i is acquired before j, where j is either another feature or an ordering. 

<p align="center">
<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/828cf552-875b-4219-9e3d-494b0c56d8c7" />
</p>

This will check and demo some of the functionality, based on HyperHMM. A more illustrative test bed is in `inst/test-bed.R`.

```
library(hyperinf)

# construct a simple dataset
data = matrix(rep(c(0,0,1, 0,1,1, 1,1,1), 10), byrow = TRUE, ncol=3, nrow=30)

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
```


References
---
[1] Aga, O.N., Brun, M., Dauda, K.A., Diaz-Uriarte, R., Giannakis, K. and Johnston, I.G., 2024. HyperTraPS-CT: Inference and prediction for accumulation pathways with flexible data and model structures. PLOS Computational Biology, 20(9), p.e1012393.

[2] Johnston, I.G. and Diaz-Uriarte, R., 2025. A hypercubic Mk model framework for capturing reversibility in disease, cancer, and evolutionary accumulation modelling. Bioinformatics, 41(1), p.btae737.

[3] Moen, M.T. and Johnston, I.G., 2023. HyperHMM: efficient inference of evolutionary and progressive dynamics on hypercubic transition graphs. Bioinformatics, 39(1), p.btac803.

[4] Giannakis, K., Aga, O.N., Moen, M.T., Drange, P.G. and Johnston, I.G., 2024. Identifying parsimonious pathways of accumulation and convergent evolution from binary data. bioRxiv, pp.2024-11.

[5] Renz, J., Brun, M. and Johnston, I.G., 2025. Flexible inference of evolutionary accumulation dynamics using uncertain observational data. arXiv preprint arXiv:2502.05872.

[6] Dauda, K.A., Aga, O.N. and Johnston, I.G., 2025. Clustering large-scale biomedical data to model dynamic accumulation processes in disease progression and anti-microbial resistance evolution. IEEE Access, 13, pp.13816-13831.
