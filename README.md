Viral Evolution ReconStructiOn (VERSO)
================

| Branch | Status |
| --- | --- |
| master | [![R-CMD-check-bioc](https://github.com/BIMIB-DISCo/VERSO/actions/workflows/check-bioc.yml/badge.svg?branch=master)](https://github.com/BIMIB-DISCo/VERSO/actions/workflows/check-bioc.yml) |
| development | [![R-CMD-check-bioc](https://github.com/BIMIB-DISCo/VERSO/actions/workflows/check-bioc.yml/badge.svg?branch=development)](https://github.com/BIMIB-DISCo/VERSO/actions/workflows/check-bioc.yml) |

In this repository we provide an R implementation of *VERSO*. 

*VERSO* is an algorithmic framework that processes variants profiles from viral samples to produce phylogenetic 
models of viral evolution. The approach solves a Boolean Matrix Factorization problem with phylogenetic constraints, 
by maximizing a log-likelihood function. *VERSO* includes two separate and subsequent steps; in this repository we provide 
an R implementation of VERSO STEP #1. 

The R version of *VERSO* can be installed from Github. To do so, we need to install the R packages *VERSO* depends on and the devtools package. 

```r
# install VERSO dependencies
if (!require("ape")) install.packages("ape")
if (!require("Rfast")) install.packages("Rfast")

# install VERSO library
if (!require("devtools")) install.packages("devtools")
library("devtools")
install_github("BIMIB-DISCo/VERSO", ref = "master")

# load VERSO library
library("VERSO")
```

Please feel free to contact us if you have problems running our tool at daniele.ramazzotti1@gmail.com or d.maspero@campus.unimib.it. 
