# write the function get_index_interval to divid the covariate into small bins.
get_index_interval = function(y,num_bin){
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
get_resAndpre = function(Data_generated,splineObj_t,splineObj_d,num_bin){
  X_in = Data_generated$X_in
  t_in = Data_generated$t_in
  ni=Data_generated$ni
  y = Data_generated$y
  k = splineObj_d$getDoF()
  # using get_index_interval to get the index in each bin
  group_interval = get_index_interval(y = y,num_bin = num_bin)
  est_C = c()
  predictor_d = matrix(0,num_bin,k)
  for(i in 1:num_bin){
    index_used = group_interval[[i]]
    X_used = X_in[index_used]
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
      tmpobs = X_used[[tmpi]]
      obsY = c(obsY,tmpobs)
      tmpt = t_used[[tmpi]]
      obsT = c(obsT,tmpt)
      tmpID = rep(tmpi,length(tmpt))
      obsID = c(obsID,tmpID)
    }
    data_used = cbind(obsID,obsT,obsY)
    data_used = data.frame(data_used)

    # set the parameters inorder to use the function: MFPCA_EstimateMLE
    controlList1 = list(alpha = 2, tol = 1e-10, iterMax = 1000, sigma = 0.05)
    controlList2 = list(alpha = 2, tol = 1e-10, sigma = 1e-2, beta = 1e-5,
                        iterMax = 1000, verbose = 1)
    optRank = 3
    mfpcaModel = MFPCA_EstimateMLE(data_used, splineObj_t, optRank, mu2 = 0,
                                   controlList1, controlList2, sigmaSq = 0.01)

    Sigma_est = mfpcaModel$SFinal[[1]]%*%mfpcaModel$SFinal[[2]]%*%t(mfpcaModel$SFinal[[1]])
    eig_Sigme = eigen(Sigma_est)
    tmp_est_C = eig_Sigme$vectors[,1:3] %*% diag(sqrt(eig_Sigme$values[1:3]))
    est_C = c(est_C,list(tmp_est_C))
  }
  resAndpre = list()
  resAndpre$est_C = est_C
  resAndpre$predictor_d = predictor_d

  return(resAndpre)
}
# use the least square method to get the initial value of beta
Init_beta = function(resAndpre_list){
  # get the dimension of the C matrix.
  rownum = dim(resAndpre_list$est_C[[1]])[1]
  colnum = dim(resAndpre_list$est_C[[1]])[2]
  num_bin = length(resAndpre_list$est_C)
  init_beta = c()
  # construct the predictor and response of the least square method.
  for (q in 1:colnum) {
    for (p in 1:rownum) {
      ly = c()
      lx = resAndpre_list$predictor_d
      print(num_bin)
      for (i in 1:num_bin) {
        tmpC = resAndpre_list$est_C[[i]]
        ly = c(ly,tmpC[p,q])
      }
      print(dim(lx))
      init_beta_pq = lm(formula =  ly ~ lx-1)
      print(summary(init_beta_pq))
      init_beta_pq = as.vector( init_beta_pq$coefficients)
      init_beta = c(init_beta,init_beta_pq)
    }
  }

  return(init_beta)
}

