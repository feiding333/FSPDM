% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Fitting_algorithm.R
\name{plotcompare}
\alias{plotcompare}
\title{plot function to compare the estimation of SFPDM eigen function and true eigen function}
\usage{
plotcompare(plot_eigenfunc, Eigen_func, Eigen_Gen = NULL, selK = NULL)
}
\arguments{
\item{plot_eigenfunc}{the estimation from the results of Estimate_Eigenfunction function}

\item{Eigen_func}{true eigenfunctions list}

\item{Eigen_Gen}{other model if users want to compare}

\item{selK}{the number of eigenfunction users want to compare}
}
\value{
the plot of comparison between eigenfunctions estimation and true eigenfunctions
}
\description{
plot function to compare the estimation of SFPDM eigen function and true eigen function
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
 phi_1 = function(t,covariate){return(cos(pi*(t+covariate))*sqrt(2))}
  phi_2 = function(t,covariate){return(sin(pi*(t+covariate))*sqrt(2))}
  phi_3 = function(t,covariate){return(cos(3*pi*(t-covariate))*sqrt(2))}
  Eigen_func = c(phi_1,phi_2,phi_3)
  plotcompare(plot_eigenfunc = eigen_function_estimation, Eigen_func = Eigen_func,selK = 1)
}
