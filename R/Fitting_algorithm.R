## train function to fit FSPDM model
train_function = function(Data_generated,Eig_num,k, beta,theta,sigma2, lambda1 = 0, lambda2 = 0,lambda3 = 0, lambda4 = 0, testIndexes = NULL,sigma2_list = NULL,maxout = 1){
  list_beta = list()
  list_init_beta = list()
  list_theta = list()
  list_sigma2 = list()
  list_upsigma2 = list()
  list_responseAndpre_used = list()
  list_init_beta = list()
  set.seed(100)
  ## set the tolerance
  torrence_theta = 0.1
  torrence_beta = 0.1
  torrence_sigma2 = 0.1
  maxinner = 1
  maxout = 1
  # get the initial value of beta
  # use FSPDM to get the estimator of the beta
  # use the class function export from c code to get the estimation of our prameters
  SPDMEstimation = new(MainAlgorithm3)
  if(is.null(sigma2_list)){
    if(!is.null(testIndexes)){
      SPDMEstimation$set_data(Data_generated$X_in[-testIndexes],Data_generated$y_in[-testIndexes],Data_generated$B_in[-testIndexes],
                              Data_generated$H_in[-testIndexes],Data_generated$ni[-testIndexes],Data_generated$Omega_eig_t,Data_generated$Omega_eig_d,Data_generated$Omega_mean_t,Data_generated$Omega_mean_y,r=Eig_num,k,lambda1,lambda2,lambda3,lambda4,(beta),theta,sigma2,log(sigma2),list(1))

    }else{
      SPDMEstimation$set_data(Data_generated$X_in,Data_generated$y_in,Data_generated$B_in,
                              Data_generated$H_in,Data_generated$ni,Data_generated$Omega_eig_t,Data_generated$Omega_eig_d,Data_generated$Omega_mean_t,Data_generated$Omega_mean_y,r=Eig_num,k,lambda1,lambda2,lambda3,lambda4,(beta),theta,sigma2,log(sigma2),list(1))
    }
  }else{
    if(!is.null(testIndexes)){
      SPDMEstimation$set_data(Data_generated$X_in[-testIndexes],Data_generated$y_in[-testIndexes],Data_generated$B_in[-testIndexes],
                              Data_generated$H_in[-testIndexes],Data_generated$ni[-testIndexes],Data_generated$Omega_eig_t,Data_generated$Omega_eig_d,Data_generated$Omega_mean_t,Data_generated$Omega_mean_y,r=Eig_num,k,lambda1,lambda2,lambda3,lambda4,(beta),theta,100000,log(sigma2),sigma2_list)

    }else{
      SPDMEstimation$set_data(Data_generated$X_in,Data_generated$y_in,Data_generated$B_in,
                              Data_generated$H_in,Data_generated$ni,Data_generated$Omega_eig_t,Data_generated$Omega_eig_d,Data_generated$Omega_mean_t,Data_generated$Omega_mean_y,r=Eig_num,k,lambda1,lambda2,lambda3,lambda4,(beta),theta,100000,log(sigma2),sigma2_list)
    }

  }
  mainmaxinter = 1
  innermaxinter = 1
  tolerance = 1e-3
  c1 = 1e-4
  c2 = 1e-4
  SPDMEstimation$set_parameters(mainmaxinter,innermaxinter,tolerance,c1,c2)
  SPDMEstimation$compute()
  beta_cur = beta
  theta_cur = theta
  sigma2_cur = sigma2
  upsigma2_cur = log(sigma2_cur)

  for (outloop in 1:maxout) {
    flag_convergence = 0
    inner = 1
    theta_old = theta_cur
    theta_cur = update_theta(SPDMEstimation=SPDMEstimation,theta_cur)
    error_theta = max(abs(theta_cur - theta_old))
    beta_last = NULL
    while (flag_convergence == 0) {
      beta_old = beta_cur
      result_cur = update_beta(SPDMEstimation=SPDMEstimation,theta_cur = theta_cur,sigma2_cur = sigma2_cur,beta_last)
      beta_cur = result_cur$beta_cur
      init_beta=result_cur$init_beta
      responseAndpre_used = result_cur$responseAndpre_used
      error_beta = max(abs(beta_cur - beta_old))
      beta_last = beta_cur
      sigma2_old = sigma2_cur
      if(is.null(sigma2_list)){
        try({
          upsigma2_cur = update_upsigma2(SPDMEstimation=SPDMEstimation,upsigma2_cur)
          sigma2_cur = exp(upsigma2_cur)
        })
        if ('try-error' %in% class(upsigma2_cur)) {
          sigma2_cur = 0.01
        }
      }
      error_sigma2 = abs(sigma2_cur-sigma2_old)
      inner = inner + 1
      if((error_theta<torrence_theta)&&(error_beta<torrence_beta)&&(error_sigma2<torrence_sigma2)){
        flag_convergence = 1
      }
      if(inner > maxinner){
        flag_convergence = 1
      }
    }

    theta_old = theta_cur
    try({
      theta_cur = update_theta(SPDMEstimation=SPDMEstimation,theta_cur)
      error_theta = max(abs(theta_cur - theta_old))
    })
    if ('try-error' %in% class(theta_cur)) {
      print("theta updating meeting the problem of NA of INF")
    }

  }
  parameter_est = list()
  parameter_est$lambda_eigen_t = lambda1
  parameter_est$lambda_eigen_d = lambda2
  parameter_est$lambda_mean_t = lambda3
  parameter_est$lambda_mean_y = lambda4
  parameter_est$beta = beta_cur
  parameter_est$theta = theta_cur
  parameter_est$sigma2 = sigma2_cur
  parameter_est$init_beta = init_beta
  parameter_est$testIndexes = testIndexes
  parameter_est$responseAndpre_used = responseAndpre_used
  return(parameter_est)
}

### update theta using the function exports from C code.
update_theta = function(SPDMEstimation = SPDMEstimation,theta_cur){
  est_usepac = optim(theta_cur, SPDMEstimation$objfunc_with_theta,SPDMEstimation$grad_with_theta,
                     method = "BFGS",control = list(maxit = 1e6,reltol = 1e-35))
  theta_cur = est_usepac$par
  SPDMEstimation$set_theta(theta_cur)
  return(theta_cur)

}

# update functions
#### update beta parameters using the function exports from C code.
update_beta = function(SPDMEstimation = SPDMEstimation,theta_cur,sigma2_cur,beta_last = NULL){
  ## initial step when update beta
  if(!is.null(beta_last)){
    init_beta = beta_last
  }else{
    responseAndpre_used = get_resAndpre(Data_generated,splineObj_t,splineObj_d,num_bin,theta_cur,sigma2_cur)
    init_beta = Init_beta(responseAndpre_used)
    init_beta = init_beta$init_beta
  }
  ## using the function export from c code to update beta
  SPDMEstimation$set_beta(init_beta)
  est_usepac = optim(init_beta, SPDMEstimation$objfunc_with_beta,SPDMEstimation$grad_with_beta,
                     method = "BFGS",control = list(maxit = 1e6,reltol = 1e-35,abstol = 1e-30))
  beta_cur = est_usepac$par
  SPDMEstimation$set_beta(beta_cur)
  result_cur = list()
  result_cur$beta_cur = beta_cur
  result_cur$init_beta = init_beta
  result_cur$responseAndpre_used = responseAndpre_used
  return(result_cur)
}
### update sigma2 using the function exports from C code.
update_sigma2 = function(SPDMEstimation= SPDMEstimation,sigma2_cur){
  est_usepac = optim(sigma2_cur, SPDMEstimation$objfunc_with_sigma2,SPDMEstimation$grad_with_sigma2,
                     method = "BFGS",control = list(maxit = 1e6,reltol = 1e-35,abstol = 1e-30))
  sigma2_cur = est_usepac$par
  SPDMEstimation$set_sigma2(sigma2_cur)
  return(sigma2_cur)
}
### update upsigma2, i.e. exp(sigma2) ,using the function exports from C code.
update_upsigma2 = function(SPDMEstimation=SPDMEstimation,upsigma2_cur){
  est_usepac = optim(upsigma2_cur, SPDMEstimation$objfunc_with_upsigma2,SPDMEstimation$grad_with_upsigma2,
                     method = "BFGS",control = list(maxit = 1e6,reltol = 1e-35,abstol = 1e-30))
  upsigma2_cur = est_usepac$par
  SPDMEstimation$set_upsigma2(upsigma2_cur)
  return(upsigma2_cur)
}
