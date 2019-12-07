#include "MainAlgorithm.hpp"


RCPP_MODULE(outfun){
  Rcpp::class_<ObjGradient>("ObjGradient3")
  .constructor()
  .method("set_theta",&ObjGradient::set_theta)
  .method("set_beta",&ObjGradient::set_beta)
  .method("set_sigma2",&ObjGradient::set_sigma2)
  .method("set_upsigma2",&ObjGradient::set_upsigma2)
  ;
  Rcpp::class_<MainAlgorithm>("MainAlgorithm3")
    .derives<ObjGradient>("ObjGradient3")
    .constructor()
    .method("set_data",&MainAlgorithm::set_data)
    .method("set_parameters",&MainAlgorithm::set_parameters)
    .method("compute",&MainAlgorithm::compute)
    .method("initialize",&MainAlgorithm::initialize)
    .method("objfunc_with_all",&MainAlgorithm::objfunc_with_all)
    .method("objfunc_with_beta",&MainAlgorithm::objfunc_with_beta)
    .method("grad_with_beta",&MainAlgorithm::grad_with_beta)
    .method("objfunc_with_theta",&MainAlgorithm::objfunc_with_theta)
    .method("grad_with_theta",&MainAlgorithm::grad_with_theta)
    .method("objfunc_with_sigma2",&MainAlgorithm::objfunc_with_sigma2)
    .method("grad_with_sigma2",&MainAlgorithm::grad_with_sigma2)
    .method("objfunc_with_upsigma2",&MainAlgorithm::objfunc_with_upsigma2)
    .method("grad_with_upsigma2",&MainAlgorithm::grad_with_upsigma2)



  ;
}
