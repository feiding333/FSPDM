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
    theta_cur = .update_theta(SPDMEstimation=SPDMEstimation,theta_cur)
    error_theta = max(abs(theta_cur - theta_old))
    beta_last = NULL
    while (flag_convergence == 0) {
      beta_old = beta_cur
      result_cur = .update_beta(SPDMEstimation=SPDMEstimation,theta_cur = theta_cur,sigma2_cur = sigma2_cur,beta_last)
      beta_cur = result_cur$beta_cur
      init_beta=result_cur$init_beta
      responseAndpre_used = result_cur$responseAndpre_used
      error_beta = max(abs(beta_cur - beta_old))
      beta_last = beta_cur
      sigma2_old = sigma2_cur
      if(is.null(sigma2_list)){
        try({
          upsigma2_cur = .update_upsigma2(SPDMEstimation=SPDMEstimation,upsigma2_cur)
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
      theta_cur = .update_theta(SPDMEstimation=SPDMEstimation,theta_cur)
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
.update_theta = function(SPDMEstimation = SPDMEstimation,theta_cur){
  est_usepac = optim(theta_cur, SPDMEstimation$objfunc_with_theta,SPDMEstimation$grad_with_theta,
                     method = "BFGS",control = list(maxit = 1e6,reltol = 1e-35))
  theta_cur = est_usepac$par
  SPDMEstimation$set_theta(theta_cur)
  return(theta_cur)

}

# update functions
#### update beta parameters using the function exports from C code.
.update_beta = function(SPDMEstimation = SPDMEstimation,theta_cur,sigma2_cur,beta_last = NULL){
  ## initial step when update beta
  if(!is.null(beta_last)){
    init_beta = beta_last
  }else{
    responseAndpre_used = .get_resAndpre(Data_generated,splineObj_t,splineObj_d,num_bin,theta_cur,sigma2_cur)
    init_beta = .Init_beta(responseAndpre_used)
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
.update_sigma2 = function(SPDMEstimation= SPDMEstimation,sigma2_cur){
  est_usepac = optim(sigma2_cur, SPDMEstimation$objfunc_with_sigma2,SPDMEstimation$grad_with_sigma2,
                     method = "BFGS",control = list(maxit = 1e6,reltol = 1e-35,abstol = 1e-30))
  sigma2_cur = est_usepac$par
  SPDMEstimation$set_sigma2(sigma2_cur)
  return(sigma2_cur)
}
### update upsigma2, i.e. exp(sigma2) ,using the function exports from C code.
.update_upsigma2 = function(SPDMEstimation=SPDMEstimation,upsigma2_cur){
  est_usepac = optim(upsigma2_cur, SPDMEstimation$objfunc_with_upsigma2,SPDMEstimation$grad_with_upsigma2,
                     method = "BFGS",control = list(maxit = 1e6,reltol = 1e-35,abstol = 1e-30))
  upsigma2_cur = est_usepac$par
  SPDMEstimation$set_upsigma2(upsigma2_cur)
  return(upsigma2_cur)
}

## get the estimation of eigenfunction.
# estimate eigenfunction
Estimate_Eigenfunction = function(est_beta,r,M1,Spline_func){
  est_beta = list(est_beta)
  phi_1 = function(t,covariate){
    return(cos(pi*(t+covariate))*sqrt(2))

  }
  phi_2 = function(t,covariate){
    return(sin(pi*(t+covariate))*sqrt(2))

  }
  phi_3 = function(t,covariate){
    return(cos(3*pi*(t-covariate))*sqrt(2))
  }
  Eigen_func = c(phi_1,phi_2,phi_3)
  construc_value = function(y){
    a = c((y+25)*2,(y+20), y+1)
    return(a)
  }
  seq_t = seq(tmin,tmax,length.out = 100)
  seq_y = seq(ymin,ymax,length.out = 12)
  Num_Repr = 1
  Plot_y_num = length(seq_y)
  num_eigen = length(Eigen_func)
  B_plot = Spline_func[[1]]$evalSpline(seq_t)
  Plot_eigfunc_1 = list()
  Plot_eigfunc_2 = list()
  Plot_eigfunc_3 = list()
  est_eigenvalue_diff = c()
  est_eigen_mean_diff = c()
  cov_list = list()
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
      # Aplot = C%*%t(C)
      # caculate the eigen value
      tmpcol = colSums(C*C)
      tmpeig_vec = C%*%diag(1/sqrt(tmpcol))
      #Plot_eig$values[1:num_eigen]
      # get the true eigen value
      zzz = construc_value(ploty)
      # get the difference
      est_eigenvalue_diff = rbind(est_eigenvalue_diff,zzz - tmpcol)
      mean_diff = rbind(mean_diff,zzz - tmpcol)

      # get the difference of convariance matrix
      tmpPlot_eigfunc = t(B_plot)%*%tmpeig_vec
      for (i in 1:num_eigen) {
        tmpPlot_ori_i = Eigen_func[[i]](seq_t,covariate = ploty)
        if(t(tmpPlot_ori_i)%*%tmpPlot_eigfunc[,i]<0){
          tmpPlot_eigfunc[,i] = -tmpPlot_eigfunc[,i]
        }
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
  }
  plot_eigenfunc = list()
  plot_eigenfunc$eigfun = c(list(Plot_eigfunc_1),list(Plot_eigfunc_2),list(Plot_eigfunc_3))
  plot_eigenfunc$num_eigen =num_eigen
  plot_eigenfunc$seq_t = seq_t
  plot_eigenfunc$seq_y = seq_y
  plot_eigenfunc$est_eigenvalue_diff = est_eigenvalue_diff
  plot_eigenfunc$cov_dif = cov_list
  plot_eigenfunc$est_eigen_mean_diff =est_eigen_mean_diff
  return(plot_eigenfunc)
}

## get the estimation of mean function
get_mean_compare = function(mean_func,Spline_func,theta_est,othermodel = NULL){
  theta_est = list(theta_est)
  plotData = data.frame();
  seq_t = seq(tmin,tmax,length.out = 100)
  seq_y = seq(ymin,ymax,length.out = 12)
  for (ploty in seq_y) {
    true_mean = mean_func(seq_t,ploty)
    Bi = Spline_func[[1]]$evalSpline(seq_t)
    Bi = t(Bi)
    tmpH_i = Spline_func[[3]]$evalSpline(ploty)
    tmp_tensor_list = list(Bi,t(tmpH_i));
    Hi = kronecker_list(tmp_tensor_list);
    rep_mean_num = length(theta_est)
    if (rep_mean_num == 1) {
      rep_mean_num_seq = 1
    }else{
      rep_mean_num_seq = 1:rep_mean_num
    }
    est_mean = c()
    for (i in rep_mean_num_seq) {
      tmp_est_mean = Hi%*%theta_est[[i]]
      est_mean = cbind(est_mean,tmp_est_mean)
    }
    ###
    est_mean = t(est_mean)
    est_mean = data.frame(est_mean)
    est_mean = melt(est_mean,measure.vars = colnames(est_mean))

    est_mean = .summarySE(est_mean, measurevar="value", groupvars="variable")
    fest_mean = est_mean$value
    ci_est = 2*est_mean$se
    ci_tru = 0
    tmp = data.frame(obsT = seq_t, obsY = true_mean, covariate_y = ploty,
                     curveID = "True Curve",ci = ci_tru, stringsAsFactors =  F)
    tmp2 = data.frame(obsT = seq_t, obsY = fest_mean, covariate_y = ploty,
                      curveID = "SFPDM",ci = ci_est,stringsAsFactors =  F)
    tmp = rbind(tmp, tmp2)
    plotData = rbind(plotData, tmp)
    ###


  }
  colnames(plotData) = c("Time", "X", "covariate","Curve","ci")
  plotData$covariate = paste('covariate:',round(plotData$covariate,2) )
  p = ggplot(plotData, aes(Time, X,
                           group = Curve, color = Curve)) +
    geom_line(aes(linetype = Curve))+geom_ribbon(aes(x = Time, y = X,ymin= X-ci,ymax=X+ci,linetype = Curve),alpha=I(1/7))
  p = p + facet_wrap(.~covariate)+
    theme_bw()+scale_linetype_manual(values=c("twodash","solid"))+
    scale_colour_manual(values = c("blue", "orange"))+
    theme(axis.title.x = element_text(size = 18),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14),legend.position = c(0.925,0.915), legend.title = element_blank())+xlab("t")+ labs(y=expression(mu~'(t,z)'))
  return(p)
}

# ***** compare the estimator of eigen functions and true eigen function*****#
plotcompare = function(plot_eigenfunc, Eigen_func,Eigen_Gen = NULL,selK = NULL){
  options(warn =-1)
  seq_t = plot_eigenfunc$seq_t
  seq_y = plot_eigenfunc$seq_y
  tmin = min(plot_eigenfunc$seq_t)
  tmax = max(plot_eigenfunc$seq_t)
  num_eigen = plot_eigenfunc$num_eigen
  nSeq = length(plot_eigenfunc$seq_t)
  plotData = data.frame();
  selR = 1:nSeq
  Plot_y_num = length(seq_y)
  # Plot_y_num is the number of covariates we want to plot
  for (tmpi in 1:Plot_y_num) {
    ploty = seq_y[tmpi]
    for(i in 1:num_eigen){
      options(warn =-1)
      fSeq0 = Eigen_func[[i]](seq_t,ploty)
      tmpfSeqHat = plot_eigenfunc$eigfun[[i]][[tmpi]]
      tmp_diff = (tmpfSeqHat - fSeq0)^2
      tmpfSeqHat = t(tmpfSeqHat)
      tmp_diff = t(tmp_diff)
      tmpfSeqHat = data.frame(tmpfSeqHat)
      tmp_diff = data.frame(tmp_diff)
      tmpfSeqHat = melt(tmpfSeqHat,measure.vars = colnames(tmpfSeqHat))
      tmp_diff = melt(tmp_diff,measure.vars = colnames(tmp_diff))

      tmpfSeqHat = .summarySE(tmpfSeqHat, measurevar="value", groupvars="variable")
      tmp_diff = .summarySE(tmp_diff, measurevar="value", groupvars="variable")
      fSeqHat = tmpfSeqHat$value
      if(!is.null(Eigen_Gen)){
        fSeqGen = Eigen_Gen[,i]

      }
      Diff = tmp_diff$value
      ci_est = 2*tmpfSeqHat$se
      Diff_ci = 2*tmp_diff$se
      ci_tru = 0
      if(sum(fSeq0*fSeqHat) < 0) fSeqHat = -fSeqHat
      if(!is.null(Eigen_Gen)){
        if(sum(fSeq0*fSeqGen) < 0) fSeqGen = -fSeqGen

      }
      tmpsquare_error = (fSeqHat - fSeq0)^2
      tmp = data.frame(obsT = seq_t, obsY = fSeq0, covariate_y = ploty,
                       pcaID =  i, curveID = "True curve",ci = ci_tru, square_error = 0,D_ci = 0, stringsAsFactors =  F)
      tmp2 = data.frame(obsT = seq_t, obsY = fSeqHat, covariate_y = ploty,
                        pcaID =  i,curveID = "SFPDM",ci = ci_est, square_error = Diff,D_ci = Diff_ci ,stringsAsFactors =  F)
      if(!is.null(Eigen_Gen)){
        tmp3 = data.frame(obsT = seq_t, obsY = fSeqGen, covariate_y = ploty,
                          pcaID =  i,curveID = "SupSFPC",ci = 0, square_error = 0,D_ci = 0 ,stringsAsFactors =  F)
        tmp = rbind(tmp, tmp2,tmp3)
      }else{
        tmp = rbind(tmp, tmp2)
      }

      plotData = rbind(plotData, tmp)

      selR = selR + nSeq
    }

  }

  colnames(plotData) = c("Time", "X", "covariate","pcaID", "Curve","ci","square_error","Dci")
  plotData$covariate = paste('covariate:',round(plotData$covariate,2) )
  if(!is.null(selK)){
    options(warn =-1)
    plotData = subset(plotData, plotData$pcaID == selK)
    p = ggplot(plotData, aes(Time, X,
                             group = Curve, color = Curve)) +
      geom_line(aes(linetype = Curve))+geom_ribbon(aes(x = Time, y = X,ymin= X-ci,ymax=X+ci,linetype = Curve ),alpha=I(1/7))
    p = p + facet_wrap(.~covariate) + theme_bw() + scale_linetype_manual(values=c("twodash", "dotted","solid"))

    p = p + scale_color_d3()+
      scale_colour_manual(values = c("blue","purple","orange") )+
      theme(axis.title.x = element_text(size = 18),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14),legend.position = c(0.925,0.88), legend.title = element_blank())+xlab("t")+ylab("f(t,z)")




  }else{
    options(warn =-1)
    p = ggplot(plotData, aes(Time, X,
                             group = Curve, color = Curve)) +
      geom_line()+geom_ribbon(aes(x = Time, y = X,ymin= X-ci,ymax=X+ci,group = "Curve"),alpha=I(1/7))
    p = p + facet_wrap(.~round(covariate,5))

  }
  return(p)
}

# # function to summary the data in order to get the error bar
.summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}
