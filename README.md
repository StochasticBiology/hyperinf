# hyperinf

An Ecosystem for Hypercubic Inference in Accumulation Modelling
---

This package brings together HyperTraPS(-CT), HyperHMM, HyperDAGs and HyperMk, different algorithms for accumulation modelling. All are based around the idea of inferring hypercubic transition graphs that describe the accumulation of binary features.

Install with

`remotes::install_github("StochasticBiology/hyperinf")`

This includes content from:
* HyperTraPS [1] https://github.com/StochasticBiology/hypertrapsct
* HyperMk [2] https://github.com/StochasticBiology/hypermk
* HyperHMM [3] https://github.com/StochasticBiology/hyperhmm
* HyperDAGs [4] https://github.com/StochasticBiology/hyperdags

Wrapper function `hyperinf` produces fitted hypercubic inference models from data; `plot_hyperinf` produces summary plots of the inferred transition graph.

<p align="center">
<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/828cf552-875b-4219-9e3d-494b0c56d8c7" />
</p>

References
---
[1] Aga, O.N., Brun, M., Dauda, K.A., Diaz-Uriarte, R., Giannakis, K. and Johnston, I.G., 2024. HyperTraPS-CT: Inference and prediction for accumulation pathways with flexible data and model structures. PLOS Computational Biology, 20(9), p.e1012393.
[2] Johnston, I.G. and Diaz-Uriarte, R., 2025. A hypercubic Mk model framework for capturing reversibility in disease, cancer, and evolutionary accumulation modelling. Bioinformatics, 41(1), p.btae737.
[3] Moen, M.T. and Johnston, I.G., 2023. HyperHMM: efficient inference of evolutionary and progressive dynamics on hypercubic transition graphs. Bioinformatics, 39(1), p.btac803.
[4] Giannakis, K., Aga, O.N., Moen, M.T., Drange, P.G. and Johnston, I.G., 2024. Identifying parsimonious pathways of accumulation and convergent evolution from binary data. bioRxiv, pp.2024-11.
