% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Fitting_algorithm.R
\name{Estimate_Eigenfunction}
\alias{Estimate_Eigenfunction}
\title{Get the estimation of Eigenfunctions using the estimators from the results of train_function}
\usage{
Estimate_Eigenfunction(est_beta, r, M1, Spline_func)
}
\arguments{
\item{est_beta}{the estimator of beta prameters, usually is the results from train_function}

\item{r}{the number of eigenfunctions}

\item{M1}{the dimension of basis with respect to time t in mean function.}

\item{Spline_func}{the spline basis function used approximate the mean and eigenfunction}
}
\value{
a list contains the estimation of eigenfunctions and the corresponding time sequence and covariate sequence
}
\description{
Get the estimation of Eigenfunctions using the estimators from the results of train_function
}
\examples{
tmin = 0 # the start point of the curve
tmax = 1 # the end point of the curve, i.e. the cuvre's regin is from tmin to tmax
ymin = 0 # the minimum of the covariates
ymax = 1 # the maximum of the covariate
num_bin = 20 # number of bin in initial steps of our algorithm
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
## set the number of principle components
Eig_num = r =  3
## the the number of spline basis used
## spline basis for t
M1 = 10
M2 = 5
M3 = k = 7
theta = rep(0,M1*M2)
sigma2 = 0.01
beta = rep(0,r*M1*M3)
parameter_best = train_function (Data_generated = Data_generated,Eig_num = Eig_num,k = k, beta = beta,theta = theta,sigma2 = sigma2)
mean_func = function(t,covariate){return(30*(t - covariate)^2)}
get_mean_compare(mean_func = mean_func,Spline_func = Spline_func,theta_est = parameter_best$theta)
eigen_function_estimation = Estimate_Eigenfunction (est_beta = parameter_best$beta,r = r,M1 = M1,Spline_func = Spline_func)
}
