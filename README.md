# FSPDM: Functional Supervised PCA by Symmetric Positive Definite Covariance Matrix Construction

The R package 'FSPDM'  is the implement of my current research about supervised functional principal component analysis(FPCA). The difference between classical FPCA and supervised FPCA is that the supervised FPCA can incorporate covariate information, often referred to as supervision information. Most supervised FPCA methods that do incorporate covariate information assume that only the scores are related to the covariates, which is often violated in practice. Our method proposes a framework to incorporate covariate information related to both the mean and covariance structure. To ensure the covariance matrix positive definite, we design a map from Euclidean space to the symmetric positive-definite matrix manifold. This package is developed to give an efficient estimation algorithm of our method and it includes the functions that fit our model on the training set and give the mean and eigenfunctions estimation. In addition, the function to give the prediction of each observation on the test set will be included as well.

## Installation

To install the latest version from Github, use
```s
library(devtools)
devtools::install_github("https://github.com/feiding333/FSPDM.git")
```
