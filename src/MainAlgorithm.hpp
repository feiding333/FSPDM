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



};


#endif

