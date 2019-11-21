#include "MainAlgorithm.hpp"

// write the store_store_intermediateData function to store the intermediate variables.
void ObjGradient::store_intermediateData(){
  mat SA_inv;
  if(sigma20!=100000){
    for(int i = 0;i < N; i++){
      BtB.at(i) = B.at(i).t()*B.at(i);
      HtH.at(i) = H.at(i).t()*H.at(i);
      BtH.at(i) = B.at(i).t()*H.at(i);
      Btx.at(i) = B.at(i).t()*X.at(i);
      Htx.at(i) = H.at(i).t()*X.at(i);
    }
  }else{
    for(int i = 0;i < N; i++){
      BtB.at(i) = B.at(i).t()*B.at(i);
      HtH.at(i) = H.at(i).t()*H.at(i);
      BtH.at(i) = B.at(i).t()*H.at(i);
      Btx.at(i) = B.at(i).t()*X.at(i);
      Htx.at(i) = H.at(i).t()*X.at(i);
      SA_inv = diagmat((1/sigma2list.at(i)));
      BtAinvB.at(i) = B.at(i).t() * SA_inv*B.at(i);
      BtAinvH.at(i) = B.at(i).t()*SA_inv*H.at(i);
      BtAinvx.at(i) = B.at(i).t()*SA_inv*X.at(i);
      HtAinvH.at(i) = H.at(i).t()*SA_inv*H.at(i);
      HtAinvx.at(i) = H.at(i).t()*SA_inv*X.at(i);
    }
  }
}

// Write the objective function.
double ObjGradient::objfunc(){
  double s =0;
  mat tmpFn;
  mat tmpxn;
  mat IdentityMr = eye<mat>(r,r);
  mat tmp_cholM;
  vec tmpgn;
  vec tmpdiagFn;
  vec SA;
  vec SA_sqrt;

  if(sigma20 != 100000){
    for (int i = 0; i < N ; i++) {
      tmp_cholM = IdentityMr + pow(sigma2,-1) * (C.at(i).t()*BtB.at(i)*C.at(i));
      //std::cout << "pow(sigma2,-1)="<< pow(sigma2,-1)  << std::endl;
      tmpgn = C.at(i).t()*(Btx.at(i) - BtH.at(i)*theta);

      // get F_n
      tmpFn = chol(tmp_cholM, "lower");

      // get x_n
      tmpxn = solve(tmpFn,tmpgn);
      tmpdiagFn = diagvec(tmpFn);
      s += 2*log(prod(tmpdiagFn)) - pow(sigma2,-2) * pow(norm(tmpxn,2),2)+ pow(sigma2,-1) * pow(norm((X.at(i)-H.at(i)*theta),2),2) + ni(i)*log(sigma2);


    }
  }else{
    for (int i = 0; i < N ; i++) {

      SA = sigma2list.at(i);
      SA_sqrt = sqrt( 1/sigma2list.at(i));

      tmp_cholM = IdentityMr + C.at(i).t()*BtAinvB.at(i)*C.at(i);
      tmpgn = (C.at(i).t())*(BtAinvx.at(i) - BtAinvH.at(i)*theta);

      // get F_n
      tmpFn = chol(tmp_cholM, "lower");
      // get x_n
      tmpxn = solve(tmpFn,tmpgn);
      tmpdiagFn = diagvec(tmpFn);
      s += 2*log(prod(tmpdiagFn)) -  pow(norm(tmpxn,2),2) + pow(norm(((X.at(i)-H.at(i)*theta)%SA_sqrt),2),2);

    }

  }

  // penalty
  s += lambda_eigen_t * as_scalar(beta.t()*Omega_eigen_t*beta)+lambda_eigen_d * as_scalar(beta.t()*Omega_eigen_d*beta)+
    lambda_mean_t * as_scalar(theta.t()*Omega_mean_t*theta)+lambda_mean_y * as_scalar(theta.t()*Omega_mean_y*theta);

  return s;
}
