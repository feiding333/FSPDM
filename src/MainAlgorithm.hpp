#ifndef F_ALGORITHM
#define F_ALGORITHM

#include "SPDMGradObj.hpp"
class MainAlgorithm: public ObjGradient{
public:

  MainAlgorithm(){
    // random initilization
    beta = beta0;
    theta = theta0;
    sigma2 = sigma20;
    upsigma2 = upsigma20;
  }

  // set data
  void set_data(vector<vec> X_,
                vector<vec> y_,
                vector<mat> B_,
                vector<mat> H_,
                vec ni_,
                mat Omega_eigen_t_,
                mat Omega_eigen_d_,
                mat Omega_mean_t_,
                mat Omega_mean_y_,
                double r_,
                double k_,double lambda_eigen_t_,double lambda_eigen_d_,double lambda_mean_t_,double lambda_mean_y_,vec beta0_,vec theta0_,double sigma20_,double upsigma20_,vector<vec> sigma2list_){
    set_data_raw (X_,y_,B_,H_,ni_,Omega_eigen_t_,Omega_eigen_d_,Omega_mean_t_,Omega_mean_y_,r_,k_,lambda_eigen_t_,lambda_eigen_d_,lambda_mean_t_,lambda_mean_y_,beta0_,theta0_,sigma20_,upsigma20_,sigma2list_);
  }

  // set parameters
  void set_parameters(double mainmaxinter_, double innermaxinter_, double tolerance_,double c1_,double c2_){
    mainmaxinter = mainmaxinter_;
    innermaxinter = innermaxinter_;
    tolerance = tolerance_;
    //lambda = lambda_;
    c1 = c1_;
    c2 = c2_;
  }

  // the objective function for all parameters
  double objfunc_with_all(vec cur_beta, vec cur_theta, double cur_sigma2);
  // gradient and objective function of beta
  vec grad_with_beta(vec beta_use);
  double objfunc_with_beta(vec mater_);
  // gradient and objective function of theta
  vec grad_with_theta(vec theta_use);
  double objfunc_with_theta(vec mater_);
  // gradient and objective function of sigma2
  double grad_with_sigma2(double sigma2_use);
  double objfunc_with_sigma2(double mater_);
  // gradient and objective function off exp(sigma2)
  double grad_with_upsigma2(double upsigma2_use);
  double objfunc_with_upsigma2(double upmater_);
  // compute function
  void compute();
  // initialize() function
  void initialize();
  // optimize() function
  arma::vec optimize();
  //protected:
protected:
  // set some parameters
  double mainmaxinter;
  double innermaxinter;
  double tolerance;
  //double lambda;
  double c1;
  double c2;


};


#endif

