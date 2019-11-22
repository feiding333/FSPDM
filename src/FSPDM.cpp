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

// gradient beta
vec ObjGradient::gradient_beta(){
  vec gbeta(m*r*k);
  vec vec_s;
  mat tmpLn;
  mat tmpDn;
  mat IdentityMm = eye<mat>(m,m);
  mat IdentityMr = eye<mat>(r,r);
  mat tmp_cholM_this;
  mat comM;
  mat tmpEn;

  if(sigma20 != 100000){

    for (int i=0; i < N; i++) {
      tmp_cholM_this = sigma2*IdentityMr + (C.at(i).t()*BtB.at(i)*C.at(i));
      tmpLn = chol(tmp_cholM_this, "lower");
      // get D_n
      tmpDn = solve(tmpLn,C.at(i).t());
      comM = IdentityMm - BtB.at(i)*tmpDn.t()*tmpDn;
      tmpEn = comM*(Btx.at(i) - BtH.at(i)*theta);
      partialC.at(i) = -2*pow(sigma2,-2)*tmpEn*tmpEn.t()*C.at(i)+ 2*pow(sigma2,-1)*comM*BtB.at(i)*C.at(i);
    }

  }else{
    // get the current \partial L \partial C_n
    for (int i=0; i < N; i++) {
      tmp_cholM_this = IdentityMr + (C.at(i).t()*BtAinvB.at(i)*C.at(i));
      tmpLn = chol(tmp_cholM_this, "lower");
      // get D_n
      tmpDn = solve(tmpLn,C.at(i).t());
      comM = IdentityMm - BtAinvB.at(i)*tmpDn.t()*tmpDn;
      tmpEn = comM*(BtAinvx.at(i) - BtAinvH.at(i)*theta);
      partialC.at(i) = -2*tmpEn*tmpEn.t()*C.at(i) + 2*comM*BtAinvB.at(i)*C.at(i);
    }
  }


  // get the gradient of whole beta
  uvec sub_loc = zeros<uvec>(k);

  for (int q = 0; q < r ; q++) {

    for (int p = 0; p < m; p++) {
      sub_loc =  regspace<uvec>(q*m*k+ p*k,q*m*k+ p*k +(k-1));
      vec_s = zeros<vec>(k);

      for (int i = 0; i < N; i++) {

        vec_s = vec_s + partialC.at(i)(p,q) * y.at(i);
      }
      gbeta.elem(sub_loc) = vec_s;
    }
  }
  gbeta += lambda_eigen_t*(Omega_eigen_t+Omega_eigen_t.t())*beta+lambda_eigen_d*(Omega_eigen_d+Omega_eigen_d.t())*beta;

  return gbeta;
}

// gtheta
// gradient theta
vec ObjGradient::gradient_theta(){
  vec gtheta(u);
  std::cout << "flag2"<< std::endl;
  vec vec_s = zeros<vec>(u);
  // *****
  mat tmpLn;
  mat tmpDn;
  mat IdentityMr = eye<mat>(r,r);
  mat tmp_cholM_this;

  if(sigma20 != 100000){
    // get the current \partial L \partial C_n

    for (int i=0; i < N; i++) {
      tmp_cholM_this = sigma2*IdentityMr + (C.at(i).t()*BtB.at(i)*C.at(i));
      tmpLn = chol(tmp_cholM_this, "lower");
      // get D_n
      tmpDn = solve(tmpLn,C.at(i).t());
      vec_s = vec_s+2*pow(sigma2,-1)*BtH.at(i).t()*tmpDn.t()*tmpDn*(Btx.at(i)-BtH.at(i)*theta);
      vec_s = vec_s-2*pow(sigma2,-1)*(Htx.at(i)-HtH.at(i)*theta);
    }

  }else{
    // get the current \partial L \partial C_n

    for (int i=0; i < N; i++) {
      tmp_cholM_this = IdentityMr + (C.at(i).t()*BtAinvB.at(i)*C.at(i));
      tmpLn = chol(tmp_cholM_this, "lower");
      // get D_n
      tmpDn = solve(tmpLn,C.at(i).t());
      vec_s = vec_s+2*BtAinvH.at(i).t()*tmpDn.t()*tmpDn*(BtAinvx.at(i)-BtAinvH.at(i)*theta);
      vec_s = vec_s-2*(HtAinvx.at(i)-HtAinvH.at(i)*theta);
    }

  }
  gtheta = vec_s + lambda_mean_t*(Omega_mean_t+Omega_mean_t.t())*theta+ lambda_mean_y*(Omega_mean_y+Omega_mean_y.t())*theta;
  // penalty
  //gtheta += lambda*(Omega+Omega.t())*beta;

  return gtheta;
}

// gsigma2
// gradient sigma2
double ObjGradient::gradient_sigma2(){
  double gsigma2;
  std::cout << "flag3"<< std::endl;
  double s = 0;
  mat tmpLn;
  mat tmpDn;
  vec Sigma_inv_y;
  mat IdentityMr = eye<mat>(r,r);
  mat tmp_cholM_this;
  // get the current \partial L \partial C_n

  for (int i=0; i < N; i++) {

    tmp_cholM_this = sigma2*IdentityMr + (C.at(i).t()*BtB.at(i)*C.at(i));

    tmpLn = chol(tmp_cholM_this, "lower");


    // get D_n
    tmpDn = solve(tmpLn,C.at(i).t());

    s = s+pow(sigma2,-1)*ni(i)-pow(sigma2,-1)*sum(diagvec(tmpDn.t()*tmpDn*BtB.at(i)));
    Sigma_inv_y = pow(sigma2,-1)*(X.at(i) - B.at(i)*tmpDn.t()*tmpDn*Btx.at(i));
    // caution
    s = s - as_scalar(Sigma_inv_y.t() * Sigma_inv_y);
  }

  gsigma2 = s;
  // penalty
  //gtheta += lambda*(Omega+Omega.t())*beta;


  return gsigma2;
}

// compute
void MainAlgorithm::compute(){
  initialize();
}

void MainAlgorithm::initialize(){
  theta = theta0;
  beta = beta0;
  upsigma2 = upsigma20;

  sigma2 = sigma20;

  store_intermediateData();

  update_C();


}
