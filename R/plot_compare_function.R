
## test the initial estimation
Test_Initial = function(responseAndpre_used,splineObj_t,bin_num,ymin,ymax){
  seq_y = seq(ymin,ymax,length.out = (bin_num+1))
  tseq = seq(0,1,length.out = 100)
  basisMat = splineObj_t$evalSpline(tseq)
  num_eigen = dim(responseAndpre_used[[1]][[1]])[2]
  difference_matrix = matrix(0,bin_num,num_eigen)
  for(tmpj in 1:num_eigen){
    for( i in 1:bin_num){
      ploty = (seq_y[i]+seq_y[i+1])/2
      coef = responseAndpre_used[[1]][[i]][,tmpj]
      tmp_est_eigenvalue = sum(coef^2)
      tmp_true_eigenvalue = construc_value(ploty)
      tmp_true_eigenvalue = tmp_true_eigenvalue[tmpj]
      tmp_difference = tmp_est_eigenvalue - tmp_true_eigenvalue
      difference_matrix[i,tmpj] = tmp_difference
    }
  }
  return(difference_matrix)
}

# obj value
pred_function = function(parameter_est,testIndexes){
  SPDMEstimation = new(MainAlgorithm3)
  SPDMEstimation$set_data(Data_generated$X_in[testIndexes],Data_generated$y_in[testIndexes],Data_generated$B_in[testIndexes],
                          Data_generated$H_in[testIndexes],Data_generated$ni[testIndexes],Data_generated$Omega_eig_t,Data_generated$Omega_eig_d,
                          Data_generated$Omega_mean_t,Data_generated$Omega_mean_y,
                          r=Eig_num,k,parameter_est$lambda_eigen_t,parameter_est$lambda_eigen_d,parameter_est$lambda_mean_t,parameter_est$lambda_mean_y,
                          (parameter_est$beta),parameter_est$theta,parameter_est$sigma2,log(parameter_est$sigma2))

  mainmaxinter = 1
  innermaxinter = 1
  tolerance = 1e-3
  c1 = 1e-4
  c2 = 1e-4
  SPDMEstimation$set_parameters(mainmaxinter,innermaxinter,tolerance,c1,c2)
  SPDMEstimation$compute()
  obj_value = SPDMEstimation$objfunc_with_all(parameter_est$beta, parameter_est$theta,parameter_est$sigma2)
  return(obj_value)
}

# cv - tuning parameters
cv_tuning = function(Data_generated,Eig_num,k,beta,theta,sigma2,folds = 10,grid_lambda1,grid_lambda2,grid_lambda3,grid_lambda4,N,sigma_list_flag = NULL){
  min__nag_log = Inf
  min_store  = c()
  fold_split <- cut(seq(1,N),breaks=folds,labels=FALSE)
  tune_vec = c(2,2,2,2)
  for (tmpi in 1:length(grid_lambda1)) {
    lambda1 = grid_lambda1[tmpi]
    for (tmpj in 1:length(grid_lambda2)) {
      lambda2 = grid_lambda2[tmpj]
      for (tmpp in 1:length(grid_lambda3)) {
        lambda3 = grid_lambda3[tmpp]
        for (tmpq in 1:length(grid_lambda4)) {
          lambda4 = grid_lambda4[tmpq]
          obj_store = c()
          for (i in 1:folds) {
            testIndexes <- which(fold_split==i,arr.ind=TRUE)
            if(is.null(sigma_list_flag)){
              parameter_est = train_function(Data_generated= Data_generated,Eig_num = Eig_num,k = k, beta= beta,theta = theta, sigma2=sigma2,lambda1 = lambda1, lambda2 = lambda2,lambda3 = lambda3,lambda4=lambda4,testIndexes =  testIndexes)
            }else{
              parameter_est = train_function(Data_generated= Data_generated,Eig_num = Eig_num,k = k, beta= beta,theta = theta, sigma2=sigma2,lambda1 = lambda1, lambda2 = lambda2,lambda3 = lambda3,lambda4=lambda4,testIndexes =  testIndexes,sigma2_list =   Data_generated$sigma2list)

            }
            obj_store = c(obj_store,pred_function(parameter_est,testIndexes))
          }
          now_nag_log = mean(obj_store)
          min_store = c(min_store,now_nag_log)
          if(now_nag_log < min__nag_log){
            min__nag_log = now_nag_log
            tune_vec = c(lambda1,lambda2,lambda3,lambda4)
          }
        }
      }
    }
  }
  tune_result = c()
  tune_result$tune_vec = tune_vec
  tune_result$min_log = min__nag_log
  tune_result$min_store = min_store
  return(tune_result)
}

# get the estimation of score
get_est_score = function(theta_cur,beta_cur,sigma2_cur,r,M1,k,Data_test,obs_num = NULL,sigmaflag = NULL){
  N_in = length(Data_test$ni)
  if(is.null(obs_num)){
    seq_use = 1:N_in
  }else{
    seq_use = obs_num
  }

  score_est = c()
  for (i in seq_use) {
    tmpni = Data_test$ni[i]
    tmpBi = Data_test$B_in[[i]]
    tmp_zi = Data_test$y_in[[i]]
    tmpxi = Data_test$X_in[[i]]
    tmpHi = Data_test$H_in[[i]]
    if(!is.null(sigmaflag)){
      sigma2_cur = Data_test$sigma2list[[i]]
    }
    tmpBThetaC = get_CDTheta(covariate_spline =tmp_zi,beta_cur=beta_cur,r = r,M1=M1,k= k)
    C = tmpBThetaC$C
    if(Eig_num != 1){
      BThetaD = (tmpBi)%*%tmpBThetaC$Theta[,1:Eig_num]%*%diag(tmpBThetaC$value[1:Eig_num])
    }else{
      BThetaD = (tmpBi)%*%tmpBThetaC$Theta[,1:Eig_num]*(tmpBThetaC$value[1:Eig_num])
    }

    if(!is.null(sigmaflag)){
      cov_inv = solve(tmpBi%*%C%*%t(C)%*%t(tmpBi) + diag(sigma2_cur))
    }else{
      cov_inv = solve(tmpBi%*%C%*%t(C)%*%t(tmpBi) + sigma2_cur*diag(tmpni))
    }
    tmp_score = t(BThetaD)%*%cov_inv%*%(tmpxi - tmpHi%*%theta_cur)
    score_est = c(score_est,list(tmp_score))
  }
  return(score_est)
}

# get estimated covariance
get_est_cov = function(theta_cur,beta_cur,sigma2_cur,r,M1,k,Data_test,obs_num = NULL,sigmaflag = NULL){
  N_in = length(Data_test$ni)
  if(is.null(obs_num)){
    seq_use = 1:N_in
  }else{
    seq_use = obs_num
  }

  score_cov = c()
  for (i in seq_use) {
    tmpni = Data_test$ni[i]
    tmpBi = Data_test$B_in[[i]]
    tmp_zi = Data_test$y_in[[i]]
    tmpxi = Data_test$X_in[[i]]
    tmpHi = Data_test$H_in[[i]]
    if(!is.null(sigmaflag)){
      sigma2_cur = Data_test$sigma2list[[i]]
    }

    tmpBThetaC = get_CDTheta(covariate_spline =tmp_zi,beta_cur=beta_cur,r = r,M1=M1,k= k)
    C = tmpBThetaC$C
    tmpcol = tmpBThetaC$value
    if(Eig_num != 1){
      BThetaD = (tmpBi)%*%tmpBThetaC$Theta[,1:Eig_num]%*%diag(tmpBThetaC$value[1:Eig_num])
    }else{
      BThetaD = (tmpBi)%*%tmpBThetaC$Theta[,1:Eig_num]*(tmpBThetaC$value[1:Eig_num])
    }
    ###
    if(!is.null(sigmaflag)){
      cov_inv = solve(tmpBi%*%C%*%t(C)%*%t(tmpBi) + diag(sigma2_cur))
    }else{
      cov_inv = solve(tmpBi%*%C%*%t(C)%*%t(tmpBi) + sigma2_cur*diag(tmpni))
    }
    if(Eig_num != 1){
      tmp_covariance = diag((tmpcol)[1:Eig_num]) -  t(BThetaD)%*%cov_inv%*%BThetaD
    }else{
      tmp_covariance = diag((tmpcol)[1:Eig_num]) -  t(BThetaD)%*%cov_inv%*%BThetaD
    }

    score_cov = c(score_cov,list(tmp_covariance))
  }
  return(score_cov)
}

# get score estmation list version
get_est_score_list = function(parameter_best_list,r,M1,k,Data_test){
  N_in = length(Data_test$ni)
  score_est = c()
  for (i in 1:N_in) {
    tmpni = Data_test$ni[i]
    tmpBi = Data_test$B_in[[i]]
    tmp_zi = Data_test$y_in[[i]]
    tmpxi = Data_test$X_in[[i]]
    tmpHi = Data_test$H_in[[i]]
    tmpBThetaC = get_CDTheta(covariate_spline =tmp_zi,beta_cur=beta_cur,r = r,M1=M1,k= k)
    C = tmpBThetaC$C
    tmpcol = tmpBThetaC$value
    if(Eig_num != 1){
      BThetaD = (tmpBi)%*%tmpBThetaC$Theta[,1:Eig_num]%*%diag(tmpBThetaC$value[1:Eig_num])
    }else{
      BThetaD = (tmpBi)%*%tmpBThetaC$Theta[,1:Eig_num]*(tmpBThetaC$value[1:Eig_num])
    }
    ##
    cov_inv = solve(tmpBi%*%C%*%t(C)%*%t(tmpBi) + sigma2_cur*diag(tmpni))
    tmp_score = t(BThetaD)%*%cov_inv%*%(tmpxi - tmpHi%*%theta_cur)

    score_est = c(score_est,list(tmp_score))
  }
  return(score_est)
}

# divid the data into different sample for each curve
# test_percentage of each curve is testing data set and the other is training set
score_compare_divid = function(test_percentage){
  Data_generated1 = list()
  Data_generated2 = list()
  ni1 = rep(0,N)
  X_in1 = list()
  y_in1 = list()
  t_in1 = list()
  H_in1 =list()
  B_in1 = list()
  eigenValue_list1 = list()
  score_list1 = list()
  ni2 = rep(0,N)
  X_in2 = list()
  y_in2 = list()
  t_in2 = list()
  H_in2 =list()
  B_in2 = list()
  eigenvalue1 = list()
  eigenvalue2 = list()
  eigenValue_list2 = list()
  score_list2 = list()
  score1 = list()
  score2 = list()
  for (i in 1:N) {
    ni1[i] = floor(Data_generated$ni[i]*(1-test_percentage))
    ni2[i] = Data_generated$ni[i] - ni1[i]
    index_1 = sample(1:Data_generated$ni[i],size = ni1[i],replace = FALSE)
    index_2 = setdiff((1:Data_generated$ni[i]), index_1)
    X_in1 = c(X_in1,list(Data_generated$X_in[[i]][index_1]))
    X_in2 = c(X_in2,list(Data_generated$X_in[[i]][index_2]))
    B_in1 = c(B_in1,list(Data_generated$B_in[[i]][index_1,]))
    B_in2 = c(B_in2,list(Data_generated$B_in[[i]][index_2,]))
    H_in1 = c(H_in1,list(Data_generated$H_in[[i]][index_1,]))
    H_in2 = c(H_in2,list(Data_generated$H_in[[i]][index_2,]))
    t_in1 = c(t_in1,list(Data_generated$t_in[[i]][index_1]))
    t_in2 = c(t_in2,list(Data_generated$t_in[[i]][index_2]))
  }
  out_list1 = list()
  out_list1$X_in = X_in1
  out_list1$y_in = Data_generated$y_in
  out_list1$B_in = B_in1
  out_list1$H_in = H_in1
  out_list1$ni = ni1
  out_list1$t_in = t_in1
  out_list1$eigenvalue = Data_generated$eigenvalue
  out_list1$score = Data_generated$score
  out_list1$y = Data_generated$y
  out_list1$Omega_eig_t = Data_generated$Omega_eig_t
  out_list1$Omega_eig_d = Data_generated$Omega_eig_d
  out_list1$Omega_mean_t = Data_generated$Omega_mean_t
  out_list1$Omega_mean_y = Data_generated$Omega_mean_y

  out_list2 = list()
  out_list2$X_in = X_in2
  out_list2$y_in = Data_generated$y_in
  out_list2$B_in = B_in2
  out_list2$H_in = H_in2
  out_list2$ni = ni2
  out_list2$t_in = t_in2
  out_list2$eigenvalue = Data_generated$eigenvalue
  out_list2$score = Data_generated$score
  out_list2$y = Data_generated$y
  out_list2$Omega_eig_t = Data_generated$Omega_eig_t
  out_list2$Omega_eig_d = Data_generated$Omega_eig_d
  out_list2$Omega_mean_t = Data_generated$Omega_mean_t
  out_list2$Omega_mean_y = Data_generated$Omega_mean_y

  dividlist = list()
  dividlist$test_data = out_list2
  dividlist$train_data = out_list1
  return(dividlist)
}

# divid the data and the removed data is equally removed
# test_percentage of each curve is testing data set and the other is training set
Data_point_eually_removed = function(test_percentage,Data_generated,sigma2flag = NULL){
  N = length(Data_generated$X_in)
  Data_generated1 = list()
  Data_generated2 = list()
  ni1 = rep(0,N)
  X_in1 = list()
  y_in1 = list()
  t_in1 = list()
  H_in1 =list()
  B_in1 = list()
  eigenValue_list1 = list()
  score_list1 = list()
  ni2 = rep(0,N)
  X_in2 = list()
  y_in2 = list()
  t_in2 = list()
  H_in2 =list()
  B_in2 = list()
  sigma2list1 = list()
  sigma2list2 = list()
  eigenvalue1 = list()
  eigenvalue2 = list()
  eigenValue_list2 = list()
  score_list2 = list()
  score1 = list()
  score2 = list()
  for (i in 1:N) {
    num_by = round(1/(test_percentage))
    index_2 = seq(2,(Data_generated$ni[i] -1), by = num_by)
    index_1 = setdiff((1:Data_generated$ni[i]), index_2)
    ni1[i] = length(index_1)
    ni2[i] = length(index_2)
    X_in1 = c(X_in1,list(Data_generated$X_in[[i]][index_1]))
    X_in2 = c(X_in2,list(Data_generated$X_in[[i]][index_2]))
    if(!is.null(sigma2flag)){
      sigma2list1 = c(sigma2list1,list(Data_generated$sigma2list[[i]][index_1]))
      sigma2list2 = c(sigma2list2,list(Data_generated$sigma2list[[i]][index_2]))
    }
    B_in1 = c(B_in1,list(Data_generated$B_in[[i]][index_1,]))
    B_in2 = c(B_in2,list(Data_generated$B_in[[i]][index_2,]))
    H_in1 = c(H_in1,list(Data_generated$H_in[[i]][index_1,]))
    H_in2 = c(H_in2,list(Data_generated$H_in[[i]][index_2,]))
    t_in1 = c(t_in1,list(Data_generated$t_in[[i]][index_1]))
    t_in2 = c(t_in2,list(Data_generated$t_in[[i]][index_2]))
  }
  out_list1 = list()
  out_list1$X_in = X_in1
  out_list1$sigma2list = sigma2list1
  out_list1$y_in = Data_generated$y_in
  out_list1$B_in = B_in1
  out_list1$H_in = H_in1
  out_list1$ni = ni1
  out_list1$t_in = t_in1
  out_list1$eigenvalue = Data_generated$eigenvalue
  out_list1$score = Data_generated$score
  out_list1$y = Data_generated$y
  out_list1$Omega_eig_t = Data_generated$Omega_eig_t
  out_list1$Omega_eig_d = Data_generated$Omega_eig_d
  out_list1$Omega_mean_t = Data_generated$Omega_mean_t
  out_list1$Omega_mean_y = Data_generated$Omega_mean_y

  out_list2 = list()
  out_list2$X_in = X_in2
  out_list2$sigma2list = sigma2list2
  out_list2$y_in = Data_generated$y_in
  out_list2$B_in = B_in2
  out_list2$H_in = H_in2
  out_list2$ni = ni2
  out_list2$t_in = t_in2
  out_list2$eigenvalue = Data_generated$eigenvalue
  out_list2$score = Data_generated$score
  out_list2$y = Data_generated$y
  out_list2$Omega_eig_t = Data_generated$Omega_eig_t
  out_list2$Omega_eig_d = Data_generated$Omega_eig_d
  out_list2$Omega_mean_t = Data_generated$Omega_mean_t
  out_list2$Omega_mean_y = Data_generated$Omega_mean_y

  dividlist = list()
  dividlist$test_data = out_list2
  dividlist$train_data = out_list1
  return(dividlist)
}

## compare the estimate of the observation at the rmoved point, compare it to the true observation and give the square error
get_sq_error_obs = function(Data_generated2,parameter_best,score_est,obs_num = NULL,score_cov = NULL,sigmaflag = NULL){
  N_in = length(Data_generated2$ni)
  if(is.null(obs_num)){
    rep_seq_use = 1:N_in
  }else{
    rep_seq_use = obs_num
  }

  beta_cur = parameter_best$beta
  theta_cur = parameter_best$theta
  sigma2_cur = parameter_best$sigma2
  obs_true_list = c()
  obs_est_list = c()
  error_list = c()
  pre_var_list = c()
  for (i in rep_seq_use) {
    true_obs = Data_generated2$X_in[[i]]
    tmpni = Data_generated2$ni[i]
    tmpBi = Data_generated2$B_in[[i]]
    tmp_zi = Data_generated2$y_in[[i]]
    tmpHi = Data_generated2$H_in[[i]]
    if(!is.null(sigmaflag)){
      sigma2_cur = Data_generated2$sigma2list[[i]]
    }
    tmpBThetaC = get_CDTheta(covariate_spline =tmp_zi,beta_cur=beta_cur,r = r,M1=M1,k= k)
    C = tmpBThetaC$C
    tmpcol = tmpBThetaC$value
    if(Eig_num != 1){
      BThetaD = (tmpBi)%*%tmpBThetaC$Theta[,1:Eig_num]%*%diag(tmpBThetaC$value[1:Eig_num])
    }else{
      BThetaD = (tmpBi)%*%tmpBThetaC$Theta[,1:Eig_num]*(tmpBThetaC$value[1:Eig_num])
    }
    #
    if(is.null(obs_num)){
      score_in = score_est[[i]]
    }else{
      score_in = score_est[[1]]
    }

    if(Eig_num != 1){
      BThetascore = (tmpBi)%*%tmpBThetaC$Theta[,1:Eig_num]%*%score_in
    }else{
      BThetascore = (tmpBi)%*%tmpBThetaC$Theta[,1:Eig_num]*score_in
    }

    est_obs = tmpHi%*%theta_cur + BThetascore
    if(!is.null(score_cov)){
      if(is.null(obs_num)){
        score_cov_in = score_cov[[i]]
      }else{
        score_cov_in = score_cov[[1]]
      }
      pre_var = ((tmpBi)%*%(tmpBThetaC$Theta[,1:Eig_num]))%*%score_cov_in%*%(t(((tmpBi)%*%(tmpBThetaC$Theta[,1:Eig_num])))) + sigma2_cur*diag(length(true_obs))

    }
    obs_true_list = c(obs_true_list,list(true_obs))
    obs_est_list = c(obs_est_list,list(est_obs))
    if(!is.null(score_cov)){
      pre_var_list =c(pre_var_list,list(diag(pre_var)))
    }

    tmperror = sum((true_obs - est_obs)^2)
    error_list = c(error_list,tmperror)
  }
  list_return = list()
  list_return$obs_true = obs_true_list
  list_return$obs_est = obs_est_list
  list_return$error_square = error_list
  list_return$pre_var = pre_var_list
  return(list_return)
}

## get the prediction of removed points
## compare the estimate of the observation at the rmoved point
get_prediction_points = function(Data_generated2,parameter_best){
  N_in = length(Data_generated2$ni)
  beta_cur = parameter_best$beta
  theta_cur = parameter_best$theta
  sigma2_cur = parameter_best$sigma2
  score_est = ore(theta_cur,beta_cur,sigma2_cur,r,M1,k,Data_generated2)# here r, M1 and k are in the global spac
  obs_est_list = c()
  for (i in 1:N_in) {
    tmpni = Data_generated2$ni[i]
    tmpBi = Data_generated2$B_in[[i]]
    tmp_zi = Data_generated2$y_in[[i]]
    tmpHi = Data_generated2$H_in[[i]]
    tmpBThetaC = get_CDTheta(covariate_spline =tmp_zi,beta_cur=beta_cur,r = r,M1=M1,k= k)
    C = tmpBThetaC$C
    tmpcol = tmpBThetaC$value
    Theta = tmpBThetaC$Theta[,1:r]


    #
    BThetascore = ((tmpBi)%*%(Theta)%*%score_est[[i]])
    est_obs = tmpHi%*%theta_cur + BThetascore
    obs_est_list = c(obs_est_list,list(est_obs))
  }
  list_return = list()
  list_return$obs_est = obs_est_list
  return(list_return)
}

# function -- loss change with covariates
loss_change_with_covariates = function(Data_test,parameter_best_list,score_est,othermodel_rersult = NULL,mean_center = NULL){
  if(!is.null(mean_center)){
    mean_center = mean_center
  }else{
    mean_center = 0
  }
  repeat_num = length(parameter_best_list)
  if (repeat_num == 1) {
    rep_seq = 1
  }else{
    rep_seq = 1:repeat_num
  }
  loss_matrix = c()
  for (tmpi in rep_seq) {
    err_compare = get_sq_error_obs(Data_test,parameter_best_list[[tmpi]],score_est)
    n_test = length(Data_test$X_in)
    loss_use = c()
    loss_gen = c()
    for(i in 31:60){
      loss_use = c(loss_use,err_compare$error_square[i])
      if(!is.null(othermodel_rersult)){
        loss_gen = c(loss_gen, sum((othermodel_rersult[,i]+mean_center - Data_test$X_in[[i]])^2))
      }

    }
    loss_matrix = cbind(loss_matrix,loss_use)
  }
  mean_loss = apply(loss_matrix,1,mean)
  se_2 = ((apply(loss_matrix,1,sd))/sqrt(repeat_num))
  covariate_use = Data_test$y[1:30]

  plotData = data.frame(covariate = covariate_use, loss = mean_loss,se_2 = se_2,Curve = 'Our Model')
  if (!is.null(othermodel_rersult)) {
    tmp_gen = data.frame(covariate = covariate_use, loss = loss_gen,se_2 = 0,Curve = 'SupSFPC')
    plotData = rbind(plotData,tmp_gen)
  }
  colnames(plotData) = c("Covariate", "Loss",'se','Curve')

  p = ggplot(plotData, aes(Covariate, Loss,
                           group = Curve, color = Curve)) +
    geom_errorbar(aes(ymin=Loss-se, ymax=Loss+se),width=0.02)+
    geom_point()+geom_line()+theme_bw()

  p = p + scale_color_npg()

  return(p)
}

# loss bin form -- change with covariates
loss_change_with_covariates_bin = function(Data_test,parameter_best_list,num_bin,train_score_data,Index_score,Index_pred,othermodel_rersult = NULL,mean_center_list = NULL){
  num_pred_point_each = length(Data_test$X_in[[1]])
  repeat_num = length(parameter_best_list)
  if (repeat_num == 1) {
    rep_seq = 1
  }else{
    rep_seq = 1:repeat_num
  }
  loss_matrix = c()
  loss_matrix_gen = c()
  #
  y = Data_test$y
  # using get_index_interval to get the index in each bin
  group_interval = get_index_interval_loss(y = y,num_bin = (num_bin+5))

  #
  for (tmpi in rep_seq) {
    score_est = get_est_score(parameter_best_list[[tmpi]]$theta,parameter_best_list[[tmpi]]$beta,parameter_best_list[[tmpi]]$sigma2,r = Eig_num,M1 = M1,k = k,Data_test = train_score_data)
    err_compare = get_sq_error_obs(Data_test,parameter_best_list[[tmpi]],score_est)
    n_test = length(Data_test$X_in)
    loss_use = c()
    loss_gen = c()
    #
    covariate_bin = c()

    # get Gen's prediction
    if(!is.null(othermodel_rersult)){
      Gen_train_results = othermodel_rersult
      gen_score_est = get_gen_score(Gen_train_results = othermodel_rersult,Data_test_score = train_score_data,Index_score = Index_score,repeat_num = tmpi)
      Gen_pred = get_gen_pred(Gen_train_results,gen_score_est,Index_pred = Index_pred,repeat_num = tmpi)
    }

    for(i in 1:(num_bin-23)){
      index_used = group_interval[[i]]
      covariate_bin = c(covariate_bin,mean(y[index_used]))
      loss_use = c(loss_use,1/(num_pred_point_each)*mean(err_compare$error_square[index_used]))
      if(!is.null(othermodel_rersult)){
        bin_err = 0
        for (j in index_used) {
          tmp_bin_err = sum((Gen_pred[,j]+mean_center_list[[tmpi]][Index_pred] - Data_test$X_in[[j]])^2)
          bin_err = bin_err + tmp_bin_err
        }
        bin_err = bin_err/length(index_used)
        loss_gen = c(loss_gen, 1/(num_pred_point_each)*bin_err)
      }
    }
    loss_matrix = cbind(loss_matrix,loss_use)
    loss_matrix_gen = cbind(loss_matrix_gen,loss_gen)
  }

  mean_loss = apply(loss_matrix,1,mean)
  mean_loss_gen = apply(loss_matrix_gen,1,mean)
  se_2 = ((apply(loss_matrix,1,sd))/sqrt(repeat_num))
  se_2_gen = ((apply(loss_matrix_gen,1,sd))/sqrt(repeat_num))
  covariate_use = covariate_bin

  plotData = data.frame(covariate = covariate_use, loss = mean_loss,se_2 = se_2,Curve = 'SFPDM')
  if (!is.null(othermodel_rersult)) {
    tmp_gen = data.frame(covariate = covariate_use, loss = mean_loss_gen,se_2 = se_2_gen,Curve = 'SupSFPC')
    plotData = rbind(plotData,tmp_gen)
  }
  colnames(plotData) = c("Covariate", "Loss",'se','Curve')

  p = ggplot(plotData, aes(exp(1.8*4.483*(Covariate-8)), Loss,
                           group = Curve, color = Curve)) +
    geom_errorbar(aes(ymin=Loss-2*se, ymax=Loss+2*se),width=0.02)+
    geom_point()+geom_line(aes(linetype = Curve))+theme_bw() +scale_linetype_manual(values=c("twodash", "dotted")) +
    scale_colour_manual(values = c("blue","purple") )+
    theme(axis.title.x = element_text(size = 20),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14),legend.position = c(0.075,0.88), legend.title = element_blank())+xlab("covariate")+ylab("loss")

  #p = p + scale_color_npg()
  p = p + scale_color_d3()

  return(p)
}

loss_change_with_covariates_bin_raw = function(Data_test,parameter_best_list,num_bin,train_score_data){

  repeat_num = length(parameter_best_list)
  if (repeat_num == 1) {
    rep_seq = 1
  }else{
    rep_seq = 1:repeat_num
  }
  loss_matrix = c()
  #
  y = Data_test$y
  # using get_index_interval to get the index in each bin
  group_interval = get_index_interval_loss(y = y,num_bin = (num_bin+5))

  #
  for (tmpi in rep_seq) {
    score_est = get_est_score(parameter_best_list[[tmpi]]$theta,parameter_best_list[[tmpi]]$beta,parameter_best_list[[tmpi]]$sigma2,r = Eig_num,M1 = M1,k = k,Data_test = train_score_data)
    err_compare = get_sq_error_obs(Data_test,parameter_best_list[[tmpi]],score_est)
    n_test = length(Data_test$X_in)
    loss_use = c()
    #
    covariate_bin = c()


    for(i in 1:(num_bin-23)){
      index_used = group_interval[[i]]
      covariate_bin = c(covariate_bin,mean(y[index_used]))
      loss_use = c(loss_use,mean(err_compare$error_square[index_used]))
    }
    loss_matrix = cbind(loss_matrix,loss_use)
  }

  mean_loss = apply(loss_matrix,1,mean)
  se_2 = ((apply(loss_matrix,1,sd))/sqrt(repeat_num))
  covariate_use = covariate_bin

  plotData = data.frame(covariate = covariate_use, loss = mean_loss,se_2 = se_2,Curve = 'Our Model')
  colnames(plotData) = c("Covariate", "Loss",'se','Curve')

  p = ggplot(plotData, aes(Covariate, Loss,
                           group = Curve, color = Curve)) +
    geom_errorbar(aes(ymin=Loss-2*se, ymax=Loss+2*se),width=0.02)+
    geom_point()+geom_line()+theme_bw()
  return(p)
}

# loss change with each observation -- caculate the mean suqare of each observation in training data
# loss bin form -- change with covariates
loss_change_with_covariates_obs = function(Data_test,parameter_best_list,train_score_data,Index_score,Index_pred,othermodel_rersult = NULL,mean_center_list = NULL){

  repeat_num = length(parameter_best_list)
  N_use = length(Data_test$X_in)
  if (repeat_num == 1) {
    rep_seq = 1
  }else{
    rep_seq = 1:repeat_num
  }
  loss_matrix = c()
  loss_matrix_gen = c()
  #
  y = Data_test$y
  # using get_index_interval to get the index in each bin

  #
  for (tmpi in rep_seq) {
    score_est = get_est_score(parameter_best_list[[tmpi]]$theta,parameter_best_list[[tmpi]]$beta,parameter_best_list[[tmpi]]$sigma2,r = Eig_num,M1 = M1,k = k,Data_test = train_score_data)
    err_compare = get_sq_error_obs(Data_test,parameter_best_list[[tmpi]],score_est)
    n_test = length(Data_test$X_in)
    loss_use = c()
    loss_gen = c()
    #
    covariate_obs = c()

    # get Gen's prediction
    if(!is.null(othermodel_rersult)){
      gen_score_est = get_gen_score(Gen_train_results = othermodel_rersult,Data_test_score = train_score_data,Index_score = Index_score,repeat_num = tmpi)
      Gen_pred = get_gen_pred(Gen_train_results,gen_score_est,Index_pred = Index_pred,repeat_num = tmpi)
    }

    for(i in 1:N_use){
      loss_use = c(loss_use,(err_compare$error_square[i]))
      if(!is.null(othermodel_rersult)){
        bin_err = 0
        tmp_bin_err = sum((Gen_pred[,j]+mean_center_list[[tmpi]][Index_pred] - Data_test$X_in[[i]])^2)
        loss_gen = c(loss_gen, tmp_bin_err)
      }
    }
    loss_matrix = cbind(loss_matrix,loss_use)
    loss_matrix_gen = cbind(loss_matrix_gen,loss_gen)
  }

  mean_loss = apply(loss_matrix,1,mean)
  mean_loss = mean(mean_loss)
  mean_loss_gen = apply(loss_matrix_gen,1,mean)
  mean_loss_gen = mean(mean_loss_gen)

  se_2 = ((apply(loss_matrix,1,sd))/sqrt(repeat_num))
  se_2_gen = ((apply(loss_matrix_gen,1,sd))/sqrt(repeat_num))
  p = list()
  p$mean_loss_our_each_obs = mean_loss
  p$mean_loss_gen_each_obs = mean_loss_gen

  return(p)
}

#
variation_function = function(Data_test,parameter_best_list,num_com,train_score_data,Index_score,Index_pred,mean_center_list,Gen_train_results= NULL,index_gen_results = NULL){
  for_mean_matrix = matrix(0,nrow = length(Data_test$X_in[[1]]),ncol = length(Data_test$X_in))
  if(!is.null(Gen_train_results)){
    for(i in 1:length(Data_test$X_in)){
      for_mean_matrix[,i] = Data_test$X_in[[i]]
    }
    mean_use = rowMeans(for_mean_matrix)
  }

  n_use = length(Data_test$X_in)
  repeat_num = length(parameter_best_list)
  if (repeat_num == 1) {
    rep_seq = 1
  }else{
    rep_seq = 1:repeat_num
  }
  s_store = c()
  s_rate_store = c()
  s_rate_store1 = c()
  s_gen_store = c()
  s_gen_rate_store = c()
  s_gen_rate_store1 = c()
  for (tmpi in rep_seq) {
    if(is.null(index_gen_results)){
      index_gen_results = tmpi
    }
    s_total = 0
    s_pca = 0
    s_pca_other = 0
    beta_cur = parameter_best_list[[tmpi]]$beta
    theta_cur = parameter_best_list[[tmpi]]$theta
    score_est = get_est_score(parameter_best_list[[tmpi]]$theta,parameter_best_list[[tmpi]]$beta,parameter_best_list[[tmpi]]$sigma2,
                              r = num_com,M1 = M1,k = k,Data_test = train_score_data)

    if(!is.null(Gen_train_results)){
      gen_score_est = get_gen_score(Gen_train_results = Gen_train_results,Data_test_score = train_score_data,Index_score = Index_score,repeat_num =index_gen_results)

      if(num_com !=1){
        est_other = (Gen_train_results$V[[tmpi]][[1]][,1:num_com])%*%(t(gen_score_est)[1:num_com,])
      }else{
        tmpmid = (t(gen_score_est)[1,])
        tmpmid = matrix(tmpmid,nrow = 1)
        est_other =(Gen_train_results$V[[tmpi]][[1]][,1])%*%(tmpmid)
      }
    }
    for (i in 1:n_use) {

      s_total = s_total + sum((Data_test$X_in[[i]] - mean_use)^2)
      # for s_pca
      true_obs = Data_test$X_in[[i]]
      tmpni = Data_test$ni[i]
      tmpBi = Data_test$B_in[[i]]
      tmp_zi = Data_test$y_in[[i]]
      tmpHi = Data_test$H_in[[i]]
      tmpBThetaC = get_CDTheta(covariate_spline =tmp_zi,beta_cur=beta_cur,r = r,M1=M1,k= k)
      C = tmpBThetaC$C
      tmpcol = tmpBThetaC$value
      Theta = tmpBThetaC$Theta
      #
      if(num_com !=1){
        BTheta = ((tmpBi)%*%(Theta[,1:num_com]))
      }else{
        BTheta = (tmpBi)%*%(Theta[,1])
      }
      if(num_com != 1){
        BThetascore = BTheta%*%score_est[[i]][1:num_com]
      }else{
        BThetascore = BTheta*score_est[[i]][1]
      }


      est_obs = tmpHi%*%theta_cur + BThetascore
      tmperror = sum((true_obs - est_obs)^2)

      # gen results
      if(!is.null(Gen_train_results)){
        tmperror_other = sum((true_obs - mean_center_list[[tmpi]][Index_pred] -est_other[,i] )^2)
        s_pca_other = s_pca_other+tmperror_other
      }


      s_pca = s_pca + tmperror
    }
    s_store = c(s_store,s_pca)
    s_rate_store = c(s_rate_store,s_pca/s_total)
    s_rate_store1 = c(s_rate_store1,(1-s_pca/s_total))
    s_gen_store = c(s_gen_store,s_pca_other)
    s_gen_rate_store = c(s_gen_rate_store,s_pca_other/s_total)
    s_gen_rate_store1 = c(s_gen_rate_store1,(1- s_pca_other/s_total))

  }
  s_pca = mean(s_store)
  s_rate = mean(s_rate_store)
  s_rate1 = mean(s_rate_store1)
  se_pac = sd(s_store)/(sqrt(repeat_num))
  se_rate = sd(s_rate_store)/(sqrt(repeat_num))
  se_rate1 = sd(s_rate_store1)/(sqrt(repeat_num))

  total_points_num = sum(Data_test$ni)
  our_average_pca = s_pca/(sum(Data_test$ni))

  # gen
  if(!is.null(Gen_train_results)){
    s_pca_other = mean(s_gen_store)
    s_gen_rate = mean(s_gen_rate_store)
    s_gen_rate1 = mean(s_gen_rate_store1)
    se_gen_pac = sd(s_gen_store)/(sqrt(repeat_num))
    se_gen_rate = sd(s_gen_rate_store)/(sqrt(repeat_num))
    se_gen_rate1 = sd(s_gen_rate_store1)/(sqrt(repeat_num))
    gen_average_pca = s_pca_other/(sum(Data_test$ni))

  }
  out_list = list()
  out_list$total_variation = s_total
  out_list$pca_variation = s_pca
  out_list$rate_variation = s_rate
  out_list$rate_variation1 = s_rate1
  out_list$se_pac = se_pac
  out_list$se_rate = se_rate
  out_list$se_rate1 = se_rate1
  out_list$our_average_pca = our_average_pca
  out_list$total_points_num = total_points_num
  if(!is.null(Gen_train_results)){
    out_list$pca_variation_other = s_pca_other
    out_list$rate_variation_other = s_gen_rate
    out_list$rate_variation1_other = s_gen_rate1
    out_list$se_pac_other = se_gen_pac
    out_list$se_rate_other = se_gen_rate
    out_list$se_rate1_other = se_gen_rate1
    out_list$gen_average_pca = gen_average_pca
  }
  return(out_list)
}

# for non-grid data
variation_function_raw = function(Data_test,parameter_best_list,num_com,train_score_data){

  n_use = length(Data_test$X_in)
  #
  repeat_num = length(parameter_best_list)
  if (repeat_num == 1) {
    rep_seq = 1
  }else{
    rep_seq = 1:repeat_num
  }
  s_store = c()
  s_rate_store = c()
  s_rate_store1 = c()
  s_gen_store = c()
  s_gen_rate_store = c()
  s_gen_rate_store1 = c()
  for (tmpi in rep_seq) {
    s_total = 0
    s_pca = 0
    s_pca_other = 0
    beta_cur = parameter_best_list[[tmpi]]$beta
    theta_cur = parameter_best_list[[tmpi]]$theta
    score_est = get_est_score(parameter_best_list[[tmpi]]$theta,parameter_best_list[[tmpi]]$beta,parameter_best_list[[tmpi]]$sigma2,
                              r = num_com,M1 = M1,k = k,Data_test = train_score_data)


    for (i in 1:n_use) {

      # for s_pca
      true_obs = Data_test$X_in[[i]]
      tmpni = Data_test$ni[i]
      tmpBi = Data_test$B_in[[i]]
      tmp_zi = Data_test$y_in[[i]]
      tmpHi = Data_test$H_in[[i]]
      tmpBThetaC = get_CDTheta(covariate_spline =tmp_zi,beta_cur=beta_cur,r = r,M1=M1,k= k)
      C = tmpBThetaC$C
      tmpcol = tmpBThetaC$value
      Theta = tmpBThetaC$Theta
      #
      if(num_com !=1){
        BTheta = ((tmpBi)%*%(Theta[,1:num_com]))
      }else{
        BTheta = (tmpBi)%*%(Theta[,1])
      }
      if(num_com != 1){
        BThetascore = BTheta%*%score_est[[i]][1:num_com]
      }else{
        BThetascore = BTheta*score_est[[i]][1]
      }
      est_obs = tmpHi%*%theta_cur + BThetascore
      tmperror = sum((true_obs - est_obs)^2)

      s_pca = s_pca + tmperror
    }
    s_store = c(s_store,s_pca)

  }
  s_pca = mean(s_store)
  se_pac = sd(s_store)/(sqrt(repeat_num))

  total_points_num = sum(Data_test$ni)
  average_pca = s_pca/(sum(Data_test$ni))
  out_list = list()
  out_list$pca_variation = s_pca
  out_list$se_pac = se_pac
  out_list$total_points_num = total_points_num
  out_list$average_pca = average_pca

  return(out_list)
}

# get score small data
get_score_data = function(Data_test,Index_score,sigma2flag = NULL){
  N = length(Data_test$X_in)
  Data_generated1 = list()
  Data_generated2 = list()
  ni1 = rep(0,N)
  X_in1 = list()
  y_in1 = list()
  t_in1 = list()
  H_in1 =list()
  B_in1 = list()
  sigma2list1 = list()
  eigenValue_list1 = list()
  score_list1 = list()
  ni2 = rep(0,N)
  X_in2 = list()
  y_in2 = list()
  t_in2 = list()
  H_in2 =list()
  B_in2 = list()
  sigma2list2 = list()
  eigenvalue1 = list()
  eigenvalue2 = list()
  eigenValue_list2 = list()
  score_list2 = list()
  score1 = list()
  score2 = list()

  for (i in 1:N) {
    if(!is.list(Index_score)){
      index_1 = Index_score
    }else{
      index_1 = Index_score[[i]]
    }

    index_2 = setdiff((1:Data_test$ni[i]),index_1)
    ni1[i] = length(index_1)
    ni2[i] = length(index_2)
    X_in1 = c(X_in1,list(Data_test$X_in[[i]][index_1]))
    X_in2 = c(X_in2,list(Data_test$X_in[[i]][index_2]))
    B_in1 = c(B_in1,list(Data_test$B_in[[i]][index_1,]))
    B_in2 = c(B_in2,list(Data_test$B_in[[i]][index_2,]))
    H_in1 = c(H_in1,list(Data_test$H_in[[i]][index_1,]))
    H_in2 = c(H_in2,list(Data_test$H_in[[i]][index_2,]))
    t_in1 = c(t_in1,list(Data_test$t_in[[i]][index_1]))
    t_in2 = c(t_in2,list(Data_test$t_in[[i]][index_2]))
    if(!is.null(sigma2flag)){
      sigma2list1 = c(sigma2list1,list(Data_test$sigma2list[[i]][index_1]))
      sigma2list2 = c(sigma2list2,list(Data_test$sigma2list[[i]][index_2]))
    }


  }
  out_list1 = list()
  out_list1$X_in = X_in1
  out_list1$y_in = Data_test$y_in
  out_list1$B_in = B_in1
  out_list1$H_in = H_in1
  out_list1$ni = ni1
  out_list1$t_in = t_in1
  out_list1$eigenvalue = Data_test$eigenvalue
  out_list1$score = Data_test$score
  out_list1$y = Data_test$y
  out_list1$Omega_eig_t = Data_test$Omega_eig_t
  out_list1$Omega_eig_d = Data_test$Omega_eig_d
  out_list1$Omega_mean_t = Data_test$Omega_mean_t
  out_list1$Omega_mean_y = Data_test$Omega_mean_y
  if(!is.null(sigma2flag)){
    out_list1$sigma2list = sigma2list1
  }



  out_list2 = list()
  out_list2$X_in = X_in2
  out_list2$y_in = Data_test$y_in
  out_list2$B_in = B_in2
  out_list2$H_in = H_in2
  out_list2$ni = ni2
  out_list2$t_in = t_in2
  out_list2$eigenvalue = Data_test$eigenvalue
  out_list2$score = Data_test$score
  out_list2$y = Data_test$y
  out_list2$y_in = Data_test$y_in
  out_list2$Omega_eig_t = Data_test$Omega_eig_t
  out_list2$Omega_eig_d = Data_test$Omega_eig_d
  out_list2$Omega_mean_t = Data_test$Omega_mean_t
  out_list2$Omega_mean_y = Data_test$Omega_mean_y
  if(!is.null(sigma2flag)){
    out_list2$sigma2list = sigma2list2
  }
  dividlist = list()
  dividlist$train_score_data = out_list1
  dividlist$predict_data = out_list2
  return(dividlist)
}

# get Gen's model score using small data
get_gen_score = function(Gen_train_results,Data_test_score,Index_score,repeat_num = NULL){
  if(is.list(Gen_train_results$V)){
    r = dim(Gen_train_results$V[[repeat_num]][[1]])[2]
    se2 = Gen_train_results$se2[[repeat_num]][[1]]
    Sf =  Gen_train_results$Sf[[repeat_num]][[1]]
    B = Gen_train_results$B[[repeat_num]][[1]]
    V = Gen_train_results$V[[repeat_num]][[1]][Index_score,]
  }else{
    r = dim(Gen_train_results$V)[2]
    se2 = Gen_train_results$se2
    Sf =  Gen_train_results$Sf
    B = Gen_train_results$B
    V = Gen_train_results$V[Index_score,]
  }

  Sfinv = solve(Sf)
  weight = solve(diag(r) + se2[1]*Sfinv)
  N = length(Data_test_score$ni)
  U = matrix(0,nrow = N,ncol = r)
  for (i in 1:N) {
    y = Data_test_score$y[i]
    if(!is.list(Gen_train_results$V)){
      x = Gen_train_results$X.gens.trian.list[[repeat_num]][i,Index_score]
    }else{
      x = Data_test_score$X_in[[i]]
    }
    tmpscore = weight%*%((se2[1]*Sfinv)%*%(t(B))%*%y+t(V)%*%x)
    U[i,] = tmpscore
  }

  return(U)
}
# get gen's score cov
get_gen_score_cov = function(Gen_train_results,r,repeat_num = NULL){
  if(is.list(Gen_train_results$V)){
    se2 = Gen_train_results$se2[[repeat_num]][[1]]
    Sf =  Gen_train_results$Sf[[repeat_num]][[1]]
  }else{
    se2 = Gen_train_results$se2
    Sf =  Gen_train_results$Sf
  }

  Sfinv = solve(Sf)
  Sfinv = Sfinv[1:r,1:r]
  weight = solve(diag(r) + se2[1]*Sfinv)
  gen_score_cov = se2[1]*weight
  return(gen_score_cov)
}

# get gen's prediction cov by score cov
# get Gen's model score using small data
get_gen_pred_cov = function(Gen_train_results,gen_score_cov,Index_pred,repeat_num){
  Gen_U_cov = gen_score_cov
  if(is.list(Gen_train_results$V)){
    Gen_V = Gen_train_results$V[[repeat_num]][[1]][Index_pred,]
  }else{
    Gen_V = Gen_train_results$V[Index_pred,]
  }

  # Eigen_Gen = sqrt(dim(Gen_V)[1])*Gen_V
  Gen_pred_cov = (Gen_V[,1:Eig_num])%*%Gen_U_cov%*%t(Gen_V[,1:Eig_num])
  return(Gen_pred_cov)
}

# get gen model's predicition
get_gen_pred = function(Gen_train_results,gen_score_est,Index_pred,repeat_num){
  Gen_U = gen_score_est
  if(is.list(Gen_train_results$V)){
    Gen_V = Gen_train_results$V[[repeat_num]][[1]][Index_pred,]
  }else{
    Gen_V = Gen_train_results$V[Index_pred,]
  }

  # Eigen_Gen = sqrt(dim(Gen_V)[1])*Gen_V
  Gen_pred = Gen_V%*%t(Gen_U)
  return(Gen_pred)
}

# get C D Theta
get_CDTheta = function(covariate_spline,beta_cur,sigma2_cur,r,M1,k){
  C = matrix(0,M1,r)
  for (q in 1:r) {
    for (p in 1:M1) {
      num_subVec = (q-1)*M1*k+(p-1)*k + 1
      C[p,q] = sum(t(covariate_spline)*beta_cur[seq(num_subVec,(num_subVec+(k-1)))])
    }
  }

  Aplot = C%*%t(C)
  Plot_eig = eigen(Aplot)
  tmpeig_vec = Plot_eig$vectors
  tmpValue = Plot_eig$values
  outlist = list()
  outlist$Theta = tmpeig_vec
  outlist$C = C
  outlist$value = tmpValue
  return(outlist)
}

## compare with cafpca
## get the eigen function only
Estimate_Eigenfunction_only = function(est_beta,seq_t,seq_y,Num_Repr = 1,r,M1,Eigen_func,Spline_func,construc_value){
  Plot_y_num = length(seq_y)
  num_eigen = length(Eigen_func)
  B_plot = Spline_func[[1]]$evalSpline(seq_t)
  Plot_eigfunc_1 = list()
  Plot_eigfunc_2 = list()
  Plot_eigfunc_3 = list()
  est_eigenvalue_diff = c()
  est_eigen_mean_diff = c()
  cov_list = list()
  error_eigen_matrix = c()
  for (ploty in seq_y) {
    tmp_zi = Spline_func[[2]]$evalSpline(ploty)
    matrix_eigfunc_1 = c()
    matrix_eigfunc_2 = c()
    matrix_eigfunc_3 = c()
    C = matrix(0,M1,r)
    cov_dif = c()
    mean_diff = c()
    for (R in 1:Num_Repr) {
      for (q in 1:r) {
        for (p in 1:M1) {
          num_subVec = (q-1)*M1*k+(p-1)*k + 1
          if(Num_Repr == 1){
            C[p,q] =sum(t(tmp_zi)*est_beta[[R]][seq(num_subVec,(num_subVec+(k-1)))])
          }else{
            C[p,q] = sum(t(tmp_zi)*est_beta[[R]][seq(num_subVec,(num_subVec+(k-1)))])
          }
        }
      }
      tmpcol = colSums(C*C)
      tmpeig_vec = C%*%diag(1/sqrt(tmpcol))
      # get the true eigen value
      zzz = construc_value(ploty)
      # get the difference
      est_eigenvalue_diff = rbind(est_eigenvalue_diff,zzz - tmpcol)
      mean_diff = rbind(mean_diff,zzz - tmpcol)

      # get the difference of convariance matrix
      tmpPlot_eigfunc = t(B_plot)%*%tmpeig_vec
      eig_three_error = c()
      for (i in 1:num_eigen) {
        tmpPlot_ori_i = Eigen_func[[i]](seq_t,covariate = ploty)
        if(t(tmpPlot_ori_i)%*%tmpPlot_eigfunc[,i]<0){
          tmpPlot_eigfunc[,i] = -tmpPlot_eigfunc[,i]
        }
        eig_three_error = c(eig_three_error, mean((tmpPlot_eigfunc[,i] -tmpPlot_ori_i)^2) )

      }
      b = cbind(Eigen_func[[1]](seq_t,covariate = ploty),Eigen_func[[2]](seq_t,covariate = ploty),Eigen_func[[3]](seq_t,covariate = ploty))
      a = cbind(tmpPlot_eigfunc[,1],tmpPlot_eigfunc[,2],tmpPlot_eigfunc[,3])
      c = diag(tmpcol)
      d = diag(construc_value(ploty))
      dif_a = sum(abs(b%*%d%*%t(b) - a%*%c%*%t(a)))
      cov_dif = c(cov_dif,dif_a)

      # get the eigen function
      matrix_eigfunc_1 = cbind(matrix_eigfunc_1,tmpPlot_eigfunc[,1])
      matrix_eigfunc_2 = cbind(matrix_eigfunc_2,tmpPlot_eigfunc[,2])
      matrix_eigfunc_3 = cbind(matrix_eigfunc_3,tmpPlot_eigfunc[,3])
    }
    tmp_meandiff = colMeans(mean_diff)
    est_eigen_mean_diff = rbind(est_eigen_mean_diff,tmp_meandiff)
    Plot_eigfunc_1 = c(Plot_eigfunc_1,list(matrix_eigfunc_1))
    Plot_eigfunc_2 = c(Plot_eigfunc_2,list(matrix_eigfunc_2))
    Plot_eigfunc_3 = c(Plot_eigfunc_3,list(matrix_eigfunc_3))
    cov_list = c(cov_list,list(cov_dif))
    error_eigen_matrix = rbind(error_eigen_matrix,eig_three_error)
  }
  plot_eigenfunc = list()
  plot_eigenfunc$eigfun = c(list(Plot_eigfunc_1),list(Plot_eigfunc_2),list(Plot_eigfunc_3))
  plot_eigenfunc$num_eigen =num_eigen
  plot_eigenfunc$seq_t = seq_t
  plot_eigenfunc$seq_y = seq_y
  plot_eigenfunc$est_eigenvalue_diff = est_eigenvalue_diff
  plot_eigenfunc$cov_dif = cov_list
  plot_eigenfunc$est_eigen_mean_diff =est_eigen_mean_diff
  plot_eigenfunc$error_eigen_matrix = error_eigen_matrix
  return(plot_eigenfunc)
}
## get mean compare only
get_mean_out_only = function(mean_func,Spline_func,seq_t,seq_y,theta_est,othermodel = NULL){
  est_mean_list = list()
  mean_error = 0
  for (ploty in seq_y) {
    true_mean = mean_func(seq_t,ploty)
    Bi = Spline_func[[1]]$evalSpline(seq_t)
    Bi = t(Bi)
    tmpH_i = Spline_func[[3]]$evalSpline(ploty)
    tmp_tensor_list = list(Bi,t(tmpH_i));
    Hi = kronecker_list(tmp_tensor_list);
    rep_mean_num = 1
    if (rep_mean_num == 1) {
      rep_mean_num_seq = 1
    }else{
      rep_mean_num_seq = 1:rep_mean_num
    }
    est_mean = c()
    for (i in rep_mean_num_seq) {
      tmp_est_mean = Hi%*%theta_est
      est_mean = cbind(est_mean,tmp_est_mean)
    }
    ###
    if(t(est_mean)%*%true_mean<0){
      est_mean = -est_mean
    }
    est_mean_list = c(est_mean_list,list(est_mean))
    mean_error =  mean_error + mean((est_mean - true_mean)^2)
  }
  mean_error = mean_error/length(seq_y)
  mean_results_list = c()
  mean_results_list$est_mean_list = est_mean_list
  mean_results_list$mean_error = mean_error
  return(mean_results_list)
}
