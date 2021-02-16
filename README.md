# An R package for estimating and forecasting financial data with MDSV model
The **MDSV** package implements the Multifractal Discrete Stochastic Volatility developed in Augustyniak, et al. (2021) for estimating and forecasting financial log-returns and realized variances (uniquely or jointly). It includes all functions for the replication of the results in Augustyniak, et al. (2021) and in my thesis. The exact replication codes are located in the folder [**test**](https://github.com/Abdoulhaki/MDSV/tree/master/test).

## Installation
### Requirements
- **MDSV** package needs [**R**](https://cran.r-project.org/) version 4.0.0 or later which can be installed on Windows. See [**CRAN**](https://cran.r-project.org/) for installation details.
- [**devtools**](https://cran.r-project.org/package=devtools) package should be installed on [**R**](https://cran.r-project.org/). If not already done, install [**devtools**](https://cran.r-project.org/package=devtools) using the code ` install.packages("devtools") `.
- [**Rtools**](https://cran.r-project.org/bin/windows/Rtools/) should be installed compatible with their [**R**](https://cran.r-project.org/) version.

### How to install
**MDSV** package can be installed from this GitHub repos using the `install_github` function of the [**devtools**](https://cran.r-project.org/package=devtools) package. All the dependencies will also be installed automatically.
```R
library(devtools)
install_github("Abdoulhaki/MDSV")
```
The option `build_vignettes = TRUE` can be added if one desires to install the vignettes (this will take several minutes).
### Load MDSV
Once the installation is done, **MDSV** can be loaded as a common package in [**R**](https://cran.r-project.org/).
```R
library(MDSV)
```
### How to use MDSV
I provide a [vignette](https://nbviewer.jupyter.org/github/Abdoulhaki/MDSV/blob/master/doc/MDSV_docx.pdf) that fully documents the package through simple and practical examples. Moreover, each function of the package has a help page accessible using `?name-of-the-function`.

