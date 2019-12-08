# write the function get_index_interval to divid the covariate into small bins.
.get_index_interval = function(y,num_bin){
  cutinterval_inter = cut(y,breaks = num_bin,labels = FALSE)
  group_interval = list()
  # get the index of each bin
  for(i in 1:num_bin){
    index_interval = which(cutinterval_inter == i )
    group_interval = c(group_interval,list(index_interval))
  }
  return(group_interval)
}

# get the mediate variable C that is estimated by the mFPCA i.e. classical PCA method.
.get_resAndpre = function(Data_generated,splineObj_t,splineObj_d,num_bin,theta_cur,sigma2_cur){
  H_in = Data_generated$H_in
  X_in = Data_generated$X_in
  t_in = Data_generated$t_in
  ni=Data_generated$ni
  y = Data_generated$y
  k = splineObj_d$getDoF()
  # using get_index_interval to get the index in each bin
  group_interval = .get_index_interval(y = y,num_bin = num_bin)
  est_C = c()
  eigen_init = c()
  predictor_d = matrix(0,num_bin,k)
  eigenFunction_bin = c()
  covariate_bin = c()
  for(i in 1:num_bin){
    index_used = group_interval[[i]]
    X_used = X_in[index_used]
    H_used = H_in[index_used]
    covariate_bin = c(covariate_bin,mean(y[index_used]))
    tmpd = splineObj_d$evalSpline(mean(y[index_used]))
    predictor_d[i,] = tmpd
    n = length(X_used)
    tmpni = ni[index_used]
    # get the observation
    # get the time point
    t_used = t_in[index_used]
    obsID = c()
    obsT = c()
    obsY =c()
    for(tmpi in 1: n ){
      tmpobs = X_used[[tmpi]]- H_used[[tmpi]]%*%theta_cur
      obsY = c(obsY,tmpobs)
      tmpt = t_used[[tmpi]]
      obsT = c(obsT,tmpt)
      tmpID = rep(tmpi,length(tmpt))
      obsID = c(obsID,tmpID)
    }
    data_used = cbind(obsID,obsT,obsY)
    data_used = data.frame(data_used)

    # set the parameters in order to use the function: MFPCA_EstimateMLE
    controlList1 = list(alpha = 1, tol = 1e-8, iterMax = 100, sigma = 0.05)
    controlList2 = list(alpha = 1e-2, tol = 1e-4, sigma = 1e-3, beta = 0.618,
                        iterMax = 200, verbose = 1)
    optRank = Eig_num
    mfpcaModel = MFPCA_EstimateMLE(data_used, splineObj_t, optRank, mu2 = 0,
                                   controlList1, controlList2, sigmaSq = sigma2_cur)
    print('mfpcaModel')
    print(mfpcaModel)
    Sigma_est = mfpcaModel$SFinal[[1]]%*%mfpcaModel$SFinal[[2]]%*%t(mfpcaModel$SFinal[[1]])
    eig_Sigme = eigen(Sigma_est)
    tmp_est_C = eig_Sigme$vectors[,1:Eig_num] %*% diag(sqrt(eig_Sigme$values[1:Eig_num]))
    tmpeigenFunctions =  mfpcaModel$eigenFunctions
    eigenFunction_bin = c(eigenFunction_bin,list(tmpeigenFunctions))
    est_C = c(est_C,list(tmp_est_C))
    eigen_init = c(eigen_init,list(mfpcaModel$eigenValues))
  }
  resAndpre = list()
  resAndpre$est_C = est_C
  resAndpre$predictor_d = predictor_d
  resAndpre$eigenFunction_bin = eigenFunction_bin
  resAndpre$covariate_bin = covariate_bin
  resAndpre$eigen_init = eigen_init
  return(resAndpre)
}
# use the least square method to get the initial value of beta
# using least square to the the initial value of beta
.Init_beta = function(resAndpre_list,num_special){
  rownum = dim(resAndpre_list$est_C[[1]])[1]
  colnum = dim(resAndpre_list$est_C[[1]])[2]
  num_bin = length(resAndpre_list$est_C)
  for (i in 2:num_bin) {
    tmp_whole_C = resAndpre_list$est_C[[(i-1)]]
    C_bin_i = resAndpre_list$est_C[[i]]
    for (j in 1:colnum) {
      fix_com = tmp_whole_C[,j]
      com_eig = C_bin_i[,j]
      if(sum(fix_com*com_eig) < 0) com_eig = -com_eig
      resAndpre_list$est_C[[i]][,j] = com_eig
    }
  }
  init_beta = c()
  list_ly = list()
  list_beta_pq = list()
  max_maxtrix = matrix(0,rownum,colnum)
  for (q in 1:colnum) {
    for (p in 1:rownum) {
      lx = resAndpre_list$predictor_d
      ly = c()
      for (i in 1:num_bin) {
        tmpC = resAndpre_list$est_C[[i]][p,q]
        ly = c(ly,tmpC)
      }

      init_beta_pq = lm(formula =  ly ~ lx-1)


      init_beta_pq = as.vector( init_beta_pq$coefficients)
      init_beta_pq[is.na(init_beta_pq)] = 0
      init_beta = c(init_beta,init_beta_pq)
      list_ly = c(list_ly,list(ly))
      list_beta_pq = c(list_beta_pq,list(init_beta_pq))
    }
  }
  init_result = list()
  init_result$list_ly = list_ly
  init_result$list_beta_pq = list_beta_pq
  init_result$init_beta = init_beta
  return(init_result)
}
