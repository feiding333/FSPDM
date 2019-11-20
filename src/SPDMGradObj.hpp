#ifndef CLASS_F_Grad_Obj
#define CLASS_F_Grad_Obj


#include <RcppArmadillo.h>
#include <vector>
#define USE_RCPP_ARMADILLO
#include "optim.hpp"
using namespace std;
using namespace arma;


class ObjGradient{
public:

  // Set data
  // R user uses this function to specify data
  void set_data_raw(vector<vec> X_,
                    vector<vec> y_,
                    vector<mat> B_,
                    vector<mat> H_,
                    vec ni_,
                    mat Omega_eigen_t_,
                    mat Omega_eigen_d_,
                    mat Omega_mean_t_,
                    mat Omega_mean_y_,
                    double r_,
                    double  k_,double lambda_eigen_t_,double lambda_eigen_d_,double lambda_mean_t_,double lambda_mean_y_,vec beta0_, vec theta0_,double sigma20_,double upsigma20_,vector<vec> sigma2list_){
    X = X_;
    y = y_;
    B = B_;
    H = H_;
    ni = ni_;
    Omega_eigen_t = Omega_eigen_t_;
    Omega_eigen_d = Omega_eigen_d_;
    Omega_mean_t = Omega_mean_t_;
    Omega_mean_y = Omega_mean_y_;
    k = k_;
    r = r_;
    N = y.size();
    m = B.at(0).n_cols;
    u = H.at(0).n_cols;
    // intermediate variables
    A.resize(N);
    invA.resize(N);
    C.resize(N);
    partialC.resize(N);
    lambda_eigen_t = lambda_eigen_t_;
    lambda_eigen_d = lambda_eigen_d_;
    lambda_mean_t = lambda_mean_t_;
    lambda_mean_y = lambda_mean_y_;

    beta0 = beta0_;
    theta0 = theta0_;
    upsigma20 = upsigma20_;
    sigma20 = sigma20_;
    sigma2list = sigma2list_;

    BtB.resize(N);
    Btx.resize(N);
    HtH.resize(N);
    BtH.resize(N);
    Htx.resize(N);
    BtAinvB.resize(N);
    BtAinvx.resize(N);
    BtAinvH.resize(N);
    HtAinvH.resize(N);
    HtAinvx.resize(N);

  }


  void set_theta(vec theta_){
    theta = theta_;
  }

  vec get_theta(){
    return theta;
  }
  void set_beta(vec beta_){
    beta = beta_;
    update_C();
  }



  vec get_beta(){
    return beta;
  }
  // for sigma2
  void set_sigma2(double sigma2_){
    sigma2 = sigma2_;
    //update_C();
  }
  double get_sigma2(){
    return sigma2;
  }
  // upsigma2
  void set_upsigma2(double upsigma2_){
    upsigma2 = upsigma2_;
    sigma2 = exp(upsigma2);
    //update_C();
  }
  double get_upsigma2(){
    return upsigma2;
  }


  void update_C(){

    int startloc,endloc;
    vec tmp_beta;
    mat tmpC(m,r);

    for (int i = 0; i < N; i++) {
      for (int q = 0; q < r; q++) {
        for (int p = 0; p < m; p++) {
          startloc = q*m*k+ p*k ;
          endloc = startloc + k-1;
          tmp_beta = beta.subvec(startloc,endloc);

          tmpC(p,q) = as_scalar( y.at(i).t()*tmp_beta);

        }
      }
      C.at(i) = tmpC;
    }
  }

  double objfunc();
  vec gradient_theta();
  vec gradient_beta();
  double gradient_sigma2();
  double gradient_upsigma2();

  void store_intermediateData();


protected:
  // Raw Data
  vector<vec> X;
  vector<vec> y;
  vector<mat> B;
  vector<mat> H;
  vector<vec> sigma2list;
  vec ni;
  vec beta0;
  vec theta0;
  double sigma20;
  double upsigma20;
  mat Omega_eigen_t;
  mat Omega_eigen_d;
  mat Omega_mean_t;
  mat Omega_mean_y;



  double N, m,u,k,r,lambda_eigen_t,lambda_eigen_d,lambda_mean_t,lambda_mean_y;

  arma::vec theta;
  arma::vec beta;
  double sigma2;
  double upsigma2;

  // intermediate parameters
  vector<mat> A;
  vector<mat> C;
  vector<mat> invA;
  vector<mat> partialC;

  // invariant intermediate data
  vector<mat> BtB;
  vector<vec> Btx;
  vector<mat> HtH;
  vector<mat> BtH;
  vector<vec> Htx;
  vector<mat> BtAinvB;
  vector<mat> BtAinvH;
  vector<vec> BtAinvx;
  vector<vec> HtAinvx;
  vector<mat> HtAinvH;



};

#endif
