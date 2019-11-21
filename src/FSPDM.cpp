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
