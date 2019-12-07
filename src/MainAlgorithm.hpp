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

};


#endif

