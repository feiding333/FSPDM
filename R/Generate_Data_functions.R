# generate the data
#Data_prameter = c(N,tmin,tmax,ymin,ymax)
# generate data function
generate_Data = function(N,tmin,tmax,ymin,ymax,sigma2,mean_func,Eigen_func,Spline_func,num_bin,nifix = NULL){
  # *******
  # penalty for t
  S_org_t = Spline_func[[1]]$get_Omega()
  seq_t_star = seq(tmin,tmax,length.out = Spline_func[[1]]$getDoF() )
  A_trans = t(Spline_func[[1]]$evalSpline(seq_t_star))
  A_inv = solve(A_trans)

  S_new_t = t(A_inv)%*%S_org_t%*%A_inv
  # penalty for d
  S_org_d = Spline_func[[2]]$get_Omega()
  seq_d_star = seq(ymin,ymax,length.out = Spline_func[[2]]$getDoF() )
  A_trans = t(Spline_func[[2]]$evalSpline(seq_d_star))
  # to be continute
  A_inv = solve(A_trans)
  S_new_d = t(A_inv)%*%S_org_d%*%A_inv
  # mean penalty
  S_org_y = Spline_func[[3]]$get_Omega()
  seq_y_star = seq(ymin,ymax,length.out = Spline_func[[3]]$getDoF() )
  A_trans = t(Spline_func[[3]]$evalSpline(seq_y_star))
  # to be continute
  A_inv = solve(A_trans)
  S_new_y = t(A_inv)%*%S_org_y%*%A_inv

  I_r = diag(1,num_eigen)
  I_t = diag(1,Spline_func[[1]]$getDoF())
  I_d = diag(1,Spline_func[[2]]$getDoF())
  I_y = diag(1,Spline_func[[3]]$getDoF())
  mid_kronecker_t = kronecker_list(list(S_new_t,I_d))
  mid_kronecker_d = kronecker_list(list(I_t,S_new_d))
  Omega_eig_t = kronecker_list(list(I_r,mid_kronecker_t))
  Omega_eig_d = kronecker_list(list(I_r,mid_kronecker_d))
  Omega_mean_t = kronecker_list(list(S_new_t,I_y))
  Omega_mean_y = kronecker_list(list(I_t,S_new_y))
  # ******
  num_eigen = length(Eigen_func)
  y = runif(N,ymin,ymax)
  ni = rep(0,N)
  X_in = list()
  y_in = list()
  t_in = list()
  H_in =list()
  B_in = list()
  eigenValue_list = list()
  score_list = list()
  for (i in 1:N) {
    if (is.null(nifix)) {
      ni[i] = sample(20:30, 1)
      tmpti =  runif(ni[i],min = tmin,max = tmax)
    }else{
      ni[i] = nifix
      tmpti = seq(tmin,tmax,length.out = nifix)
    }

    t_in = c(t_in, list(tmpti))

    # eigen value
    tmpeigen = c((y[[i]]+25)*2,(y[[i]]+20), y[[i]]+1)
    eigenValue_list = c(eigenValue_list,list(tmpeigen))

    # score
    tmpd = diag(tmpeigen)
    score = mvrnorm(n = 1,rep(0,num_eigen), tmpd)
    score_list = c(score_list,list(score))

    # mean function
    tmpmean = mean_func(tmpti,y[[i]])
    # combine the three eigenfunction
    tmpfirst = Eigen_func[[1]](tmpti,y[[i]])
    tmpsecond = Eigen_func[[2]](tmpti,y[[i]])
    tmpthird = Eigen_func[[3]](tmpti,y[[i]])
    tmpu = cbind(tmpfirst,tmpsecond,tmpthird)

    # error term
    tmpmesureerror = rnorm(0, ni[i], sigma2)

    # generate the observation
    # caution 0*mean
    sample_Sigma = tmpu%*%tmpd%*%t(tmpu)+sigma2*diag(ni[i])
    #tmpX = 0*tmpmean + tmpu%*%score + tmpmesureerror
    tmpX = mvrnorm(n = 1,tmpmean,sample_Sigma)
    X_in = c(X_in,list(tmpX))

    # spline B(t)
    Bi = Spline_func[[1]]$evalSpline(tmpti)
    Bi = t(Bi)
    B_in = c(B_in,list(Bi))

    # spline d(covariate) that in matrix C, writen as y_in
    tmp_d = t(Spline_func[[2]]$evalSpline(y[i]))
    y_in = c(y_in,list(tmp_d));

    # tensor product spline Hi
    tmpH_i = Spline_func[[3]]$evalSpline(y[i])
    tmp_tensor_list = list(Bi,t(tmpH_i));
    Hi = kronecker_list(tmp_tensor_list);
    H_in = c(H_in,list((Hi)));
  }



  #*******************##
  out_list = list()
  out_list$X_in = X_in
  out_list$y_in = y_in
  out_list$B_in = B_in
  out_list$H_in = H_in
  out_list$ni = ni
  out_list$t_in = t_in
  out_list$eigenvalue = eigenValue_list
  out_list$score = score_list
  out_list$y = y
  out_list$Omega_eig_t = Omega_eig_t
  out_list$Omega_eig_d = Omega_eig_d
  out_list$Omega_mean_t = Omega_mean_t
  out_list$Omega_mean_y = Omega_mean_y
  return(out_list)
  }
