# FSPDM: Functional Supervised PCA by Symmetric Positive Definite Covariance Matrix Construction

The R package 'FSPDM'  is the implement of my current research about supervised functional principal component analysis(FPCA). The difference between classical FPCA and supervised FPCA is that the supervised FPCA can incorporate covariate information, often referred to as supervision information. Most supervised FPCA methods that do incorporate covariate information assume that only the scores are related to the covariates, which is often violated in practice. Our method proposes a framework to incorporate covariate information related to both the mean and covariance structure. To ensure the covariance matrix positive definite, we design a map from Euclidean space to the symmetric positive-definite matrix manifold. This package is developed to give an efficient estimation algorithm of our method and it includes the functions that fit our model on the training set and give the mean and eigenfunctions estimation. In addition, the function to give the prediction of each observation on the test set will be included as well.

## Installation

To install the latest version of our package FSPDM from Github, use
```s
library(devtools)
devtools::install_github('https://github.com/feiding333/FSPDM')
library(FSPDM)
```

### Installation of another three packages
Our packages FSPDM depend on another three packages named FDABasics, mFPCA, and rOptManifold. These three packages provide basis construction functions and classical unsupervised PPC method which will be used in the initial steps of our own method. Please install the following three packages and this will take about 1-2 minutes.

```s
devtools::install_github('https://github.com/feiding333/mFPCA_a')
library(mFPCA)
devtools::install_github('https://github.com/feiding333/FDABasis')
library(FDABasics)
devtools::install_github('https://github.com/feiding333/rOptManifold_1.1')
library(rOptManifold)
```

## Usage
Library other packages, theses packages all can be installed from CRAN.
```s
## library require package
library(ggplot2)
library(reshape2)
library(rTensor)
library(MASS)
library(ggsci)
```

# Example
## Simulate Data
load simulated data, it is in the package data folder named data.
```s
load('data/data')
```
# construct the spline basis that used in our model
```s
tmin = 0 # the start point of the curve
tmax = 1 # the end point of the curve, i.e. the cuvre's regin is from tmin to tmax
ymin = 0 # the minimum of the covariates 
ymax = 1 # the maximum of the covariate 
num_bin = 20 # number of bin in initial steps of our algorithm
# spline basis for t
order = 4 # spline order
nknots =8 # number of knots
splineObj_t = new(orthoSpline,tmin,tmax,order,nknots)
# degree freedom of spline basis
M1 = splineObj_t$getDoF()
## basis with respect to y to get the tensor product basis
yknots = 3
splineObj_y = new(orthoSpline,ymin,ymax,order,(yknots))
# degree freedom of spline basis
M2 = splineObj_y$getDoF()
## basis with d that is in matrix C
# basis with respect y (d^T), in C matrix
dknots = 5
splineObj_d = new(orthoSpline,ymin,ymax,order,dknots)
# degree freedom of spline basis
M3 = splineObj_d$getDoF()
k = M3
Spline_func = c(splineObj_t,splineObj_d,splineObj_y)
```
## the the setting of algorithm
```s
## set the number of principle components
Eig_num = r =  3
## the the number of spline basis used
## spline basis for t 
M1 = 10
## spline basis for covariate in mean function
M2 = 5
## spline basis for covariate in  covariate structure
M3 = k = 7
## set the start point of parameters(option)
## the start point of theta
theta = rep(0,M1*M2)
## the start point of sigma2
sigma2 = 0.01
## the start point of beta
beta = rep(0,r*M1*M3)
```

## apply our FSPDM method to fit the model
```s
parameter_best = train_function (Data_generated = Data_generated,Eig_num = Eig_num,k = k, beta = beta,theta = theta,sigma2 = sigma2)
```
## use the parameters to get the estimation of mean function 
### construct true mean function and compare our method with true mean function
```s
mean_func = function(t,covariate){
  return(30*(t - covariate)^2)
}

get_mean_compare(mean_func = mean_func,Spline_func = Spline_func,theta_est = parameter_best$theta)
```

## use the parameters to get the estimation of eigenfunctions
```s
eigen_function_estimation = Estimate_Eigenfunction (est_beta = parameter_best$beta,r = r,M1 = M1,Spline_func = Spline_func)
```
## construct the true eigenfunction and compare our eigenfunctions estimation to true eigenfunction
```s
 phi_1 = function(t,covariate){
    return(cos(pi*(t+covariate))*sqrt(2))

  }
  phi_2 = function(t,covariate){
    return(sin(pi*(t+covariate))*sqrt(2))

  }
  phi_3 = function(t,covariate){
    return(cos(3*pi*(t-covariate))*sqrt(2))
  }
Eigen_func = c(phi_1,phi_2,phi_3)
```
# compare the eigenfunction estimation of our method with true eigenfunctions
```s
## comparation of first eigenfunctions
plotcompare(plot_eigenfunc = eigen_function_estimation, Eigen_func = Eigen_func,selK = 1)
## comparation of second eigenfunctions
plotcompare(plot_eigenfunc = eigen_function_estimation, Eigen_func = Eigen_func,selK = 2)
## comparation of third eigenfunctions
plotcompare(plot_eigenfunc = eigen_function_estimation, Eigen_func = Eigen_func,selK = 2)
```
