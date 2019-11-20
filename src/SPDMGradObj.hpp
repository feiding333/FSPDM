#ifndef CLASS_F_Grad_Obj
#define CLASS_F_Grad_Obj


#include <RcppArmadillo.h>
#include <vector>
#define USE_RCPP_ARMADILLO
#include "optim.hpp"
using namespace std;
using namespace arma;


class ObjGradient{



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
