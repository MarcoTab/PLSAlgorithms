rm(list=ls())
library(Matrix)
library(MASS)
library(Renvlp)
library(xtable)

source("our_functions_Feb2nd2022.R")
source("rrenv.R")
source("rrenvMU.R")
source("env.apweights.R")
source("rrenv.apweights.R")

# function to generate simulations at the ith repetition

sim_fun_a=function(i,p,r,dx,dy,N,M){
  set.seed(1)
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%   Envelopes model parameters          %
  #  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  trueG_temp=matrix(runif(dx*p,0,1),nrow=p)
  #print(trueG_temp)
  trueG = svd(trueG_temp)$u;   # X-envelope basis
  trueG0 =qr.Q(qr(trueG),complete=TRUE) [,-(1:dx)]
  
  trueH_temp=matrix(runif(dy*r,0,1),nrow=r)
  trueH = svd(trueH_temp)$u;   # Y-envelope basis
  trueH0 =qr.Q(qr(trueH),complete=TRUE) [,-(1:dy)]
  
  Eta = matrix(runif(dx*dy),nrow=dx)          # regression coefficient matrix of reduced Y on reduced X
  
 
  # case where simultaneous env is better than pls, and two block is not the best
  
  temp1=matrix(runif(dx*dx),nrow=dx)
  Omega =t(temp1)%*%temp1;    # material variation in X
  temp2=matrix(runif(dy*dy),nrow=dy)
  Phi = t(temp2)%*%temp2;      # material variation of Y|X
  temp3=matrix(runif((p-dx)^2),nrow=p-dx)
  Omega0 = t(temp3)%*%temp3;    # immaterial variation in X
  temp4=matrix(runif((r-dy)^2),nrow=r-dy)
  Phi0 = t(temp4)%*%temp4;       # immaterial variation of Y|X
  

  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #  %  Regression coefficient and covariance matrices      %
  #  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  beta = trueH%*%t(Eta)%*%t(trueG);
  sigmaXY = trueG%*%Omega%*%Eta%*%t(trueH);
  sigmaX = trueG%*%Omega%*%t(trueG) + trueG0%*%Omega0%*%t(trueG0);
  sigmaY = trueH%*%(Phi+t(Eta)%*%Omega%*%Eta)%*%t(trueH) + trueH0%*%Phi0%*%t(trueH0);
  sigmaC = rbind(cbind(sigmaX, sigmaXY),cbind(t(sigmaXY), sigmaY));
  sigmaD = bdiag(sigmaX,sigmaY);
  
  d_twoblock_pop=cal_d_twoblock_fun(SigmaX=sigmaX,SigmaY=sigmaY,SigmaXY=sigmaXY,mytol=1e-3)
  if(i==1){
    cat(sprintf("The dimension of population two block PLS is %d",d_twoblock_pop),"\n")
  }
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%   Data generation          %
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  set.seed(i)
  #print(i)
  datavec = mvrnorm(N+M, rep(0, p+r), sigmaC)
  X = datavec[,1:p];
  #print(X[1:3,])
  Y = datavec[,(p+1):(p+r)]
  Xtrain<-X[(1:N),] # extract training data
  Ytrain<-Y[(1:N),] # extract training data
  Xtest<-X[-(1:N),] # extract test data, the original paper only calculate the in sample MSE
  Ytest<-Y[-(1:N),]# extract test data
  if.Xscale=F
  if.Yscale=F
  
  d_pls=dx
  
  predMSE2_simul_pls_test=NA
  predMSE1_simul_pls_test=NA
  betaMSE_simul_pls=NA
  
  predMSE2_twoblock_test=NA
  predMSE1_twoblock_test=NA
  betaMSE_twoblock=NA
  
  
  predMSE2_pls_test=NA
  predMSE1_pls_test=NA
  betaMSE_pls=NA
  
  
  
  if (max(r,dx,dy,d_twoblock_pop)<=(N)){
    
    refit_result=refit_entire_training_and_get_testMSE_fun(Xtrain,Ytrain,Xtest,Ytest,
                                                           #d=min(dx,dy),
                                                           d=d_twoblock_pop,
                                                           d1=dx,
                                                           d2=dy,d_pls=c(dx,r),d_pls1=cbind(rep(dx,r),1),if.Xscale=if.Xscale,
                                                           if.Yscale=if.Yscale)
    
  }
  
  
  if (max(r,dx,dy,d_twoblock_pop)<=(N)){
    predMSE2_twoblock_test=mean(apply(refit_result$twoblock_Y_MSE.v, 2, mean))
    predMSE1_twoblock_test=mean(sqrt(apply(refit_result$twoblock_Y_MSE.v, 1, sum)))
    betaMSE_twoblock=mean(sqrt(apply((refit_result$twoblock_hat_beta-t(beta))^2,1,sum)))
    
    
    predMSE2_pls_test=mean(apply(refit_result$pls_Y_MSE.v, 2, mean))
    predMSE1_pls_test=mean(sqrt(apply(refit_result$pls_Y_MSE.v, 1, sum)))
    betaMSE_pls=mean(sqrt(apply((refit_result$pls_hat_beta-t(beta))^2,1,sum)))
    
    
    
    predMSE2_simul_pls_test=mean(apply(refit_result$simul_pls_Y_MSE.v, 2, mean))
    predMSE1_simul_pls_test=mean(sqrt(apply(refit_result$simul_pls_Y_MSE.v, 1, sum)))
    betaMSE_simul_pls=mean(sqrt(apply((refit_result$simul_pls_hat_beta-t(beta))^2,1,sum)))
    
  }
  
  
  ##only if.Xscale=FALSE and if.Yscale=FALSE
  
  if ((d_twoblock_pop<=(N)) & (max(r,dx,dy,d_twoblock_pop)>(N))){
    
    
    
    X=Xtrain
    Y=Ytrain
    S_XY=cov(X,Y)
    S_YX=cov(Y,X)
    S_X=var(X)
    S_Y=var(Y)
    mean_X.m=matrix(rep(colMeans(Xtrain),each=nrow(Xtest)),ncol=ncol(X))
    mean_Y.m=matrix(rep(colMeans(Ytrain),each=nrow(Ytest)),ncol=ncol(Y))
    
    twoblock_hat_beta=matrix(NA,nrow=ncol(X),ncol=ncol(Y))
    twoblock_hat_Y=Ytest*NA
    twoblock_Y_MSE.v=Ytest*NA
    tryCatch(
      {
        twoblock_result=twoblock_pls_fun(X,Y,d=d_twoblock_pop)
        twoblock_hat_beta=twoblock_result$u.m%*%solve(t(twoblock_result$u.m)%*%S_X%*%twoblock_result$u.m)%*%
          t(twoblock_result$u.m)%*%S_XY%*%twoblock_result$v.m%*%t(twoblock_result$v.m)
        twoblock_hat_Y=mean_Y.m+(Xtest-mean_X.m)%*%twoblock_hat_beta#### changed
        twoblock_Y_MSE.v=(Ytest-twoblock_hat_Y)^2  
      },
      error=function(cond) {
        # message("Here's the original error message:")
        message(cond)
      },
      warning=function(cond) {
        # message("Here's the original warning message:")
        message(cond)
      }
    )    
    predMSE2_twoblock_test=mean(apply(twoblock_Y_MSE.v, 2, mean))
    predMSE1_twoblock_test=mean(sqrt(apply(twoblock_Y_MSE.v, 1, sum)))
    betaMSE_twoblock=mean(sqrt(apply((twoblock_hat_beta-t(beta))^2,1,sum)))
    
    
  }
  
  ##simul pls when it is possible no scalling
  
  if ((dx<=(N) & dy<=(N)) & (max(r,dx,dy,d_twoblock_pop)>(N))){
    
    X=Xtrain
    Y=Ytrain
    S_XY=cov(X,Y)
    S_YX=cov(Y,X)
    S_X=var(X)
    S_Y=var(Y)
    mean_X.m=matrix(rep(colMeans(Xtrain),each=nrow(Xtest)),ncol=ncol(X))
    mean_Y.m=matrix(rep(colMeans(Ytrain),each=nrow(Ytest)),ncol=ncol(Y))
    
    d1=unlist(dx)
    d2=unlist(dy)
    
    simul_pls_hat_beta=matrix(NA,nrow=ncol(X),ncol=ncol(Y))
    simul_pls_hat_Y=Ytest*NA
    simul_pls_Y_MSE.v=Ytest*NA
    
    tryCatch(
      {
        simul_pls_result=simul_pls_fun(X,Y,d1=d1,d2=d2)
        simul_pls_hat_beta=simul_pls_result$u.m%*%solve(t(simul_pls_result$u.m)%*%S_X%*%simul_pls_result$u.m)%*%
          t(simul_pls_result$u.m)%*%S_XY%*%simul_pls_result$v.m%*%t(simul_pls_result$v.m)
        simul_pls_hat_Y=mean_Y.m+(Xtest-mean_X.m)%*%simul_pls_hat_beta
        simul_pls_Y_MSE.v=(Ytest-simul_pls_hat_Y)^2
      })
    
    
    predMSE2_simul_pls_test=mean(apply(simul_pls_Y_MSE.v, 2, mean))
    predMSE1_simul_pls_test=mean(sqrt(apply(simul_pls_Y_MSE.v, 1, sum)))
    betaMSE_simul_pls=mean(sqrt(apply((simul_pls_hat_beta-t(beta))^2,1,sum)))
  }
  
  
  
  
  ###X-pls when it is possible no scalling
  
  if ((dx<=(N)) & (max(r,dx,dy,d_twoblock_pop)>(N))){
    
    X=Xtrain
    Y=Ytrain
    S_XY=cov(X,Y)
    S_YX=cov(Y,X)
    S_X=var(X)
    S_Y=var(Y)
    mean_X.m=matrix(rep(colMeans(Xtrain),each=nrow(Xtest)),ncol=ncol(X))
    mean_Y.m=matrix(rep(colMeans(Ytrain),each=nrow(Ytest)),ncol=ncol(Y))
    
    d1=unlist(dx)
    d2=unlist(dy)
    
    simul_pls_hat_beta=matrix(NA,nrow=ncol(X),ncol=ncol(Y))
    simul_pls_hat_Y=Ytest*NA
    simul_pls_Y_MSE.v=Ytest*NA
    
    
    
    
    
    tryCatch(
      {
        pls_result=simul_pls_fun(X,Y,d1=dx,d2=r)
        pls_hat_beta=pls_result$u.m%*%solve(t(pls_result$u.m)%*%S_X%*%pls_result$u.m)%*%
          t(pls_result$u.m)%*%S_XY%*%pls_result$v.m%*%t(pls_result$v.m)
        pls_hat_Y=mean_Y.m+(Xtest-mean_X.m)%*%pls_hat_beta 
        pls_Y_MSE.v=(Ytest-pls_hat_Y)^2
      })
    
    
    predMSE2_pls_test=mean(apply(pls_Y_MSE.v, 2, mean))
    predMSE1_pls_test=mean(sqrt(apply(pls_Y_MSE.v, 1, sum)))
    betaMSE_pls=mean(sqrt(apply((pls_hat_beta-t(beta))^2,1,sum)))
  }
  
  
  
  
  
  
  # Just to use simultaneous PLS when no dimension reduction = OLS to get test OLS error
  
  predMSE2_ols_test=NA
  predMSE1_ols_test=NA
  betaMSE_ols=NA
  
  
  if (p<=(N) & r<=(N)){
    
    
    tryCatch(
      {        refit1_result=refit_entire_training_and_get_testMSE_fun(Xtrain,Ytrain,Xtest,Ytest,
                                                                       d=min(dx,dy),d1=p,
                                                                       d2=r,d_pls=c(dx,r),d_pls1=cbind(rep(dx,r),1),if.Xscale=if.Xscale,
                                                                       if.Yscale=if.Yscale)
      predMSE2_ols_test=mean(apply(refit1_result$simul_pls_Y_MSE.v, 2, mean))
      predMSE1_ols_test=mean(sqrt(apply(refit1_result$simul_pls_Y_MSE.v, 1, sum)))
      betaMSE_ols=mean(sqrt(apply((refit1_result$simul_pls_hat_beta-t(beta))^2,1,sum)))
      },
      error=function(cond) {
        # message("Here's the original error message:")
        message(cond)
      },
      warning=function(cond) {
        # message("Here's the original warning message:")
        message(cond)
      }
    )    
  }
  
  predMSE2_stenv_test=NA
  predMSE1_stenv_test=NA
  betaMSE_stenv=NA
  
  if (p<=(N) & r<=(N)){
    tryCatch(
      {
        fit <- stenv(Xtrain, Ytrain, dx,dy, asy = F)
        betahat_stenv <- fit$beta
        muhat_stenv <- fit$mu
        resi <- as.matrix(Ytest - matrix(1, M, 1) %*% 
                            t(muhat_stenv) - as.matrix(Xtest) %*% (betahat_stenv)) # stenv_Y_MSE
        
        
        predMSE2_stenv_test=mean(apply(resi^2, 2, mean))
        predMSE1_stenv_test=mean(sqrt(apply(resi^2,1, sum)))
        betaMSE_stenv=mean(sqrt(apply((betahat_stenv-t(beta))^2,1,sum)))
      },
      error=function(cond) {
        # message("Here's the original error message:")
        message(cond)
      },
      warning=function(cond) {
        # message("Here's the original warning message:")
        message(cond)
      }
    )    
  }
  
  
  
  
  
  
  
  
  
  result_to_save=c(predMSE2_ols_test = predMSE2_ols_test, predMSE1_ols_test = predMSE1_ols_test,
                   betaMSE_ols = betaMSE_ols, predMSE2_twoblock_test = predMSE2_twoblock_test, 
                   predMSE1_twoblock_test = predMSE1_twoblock_test,
                   betaMSE_twoblock = betaMSE_twoblock,predMSE2_stenv_test = predMSE2_stenv_test,
                   predMSE1_stenv_test = predMSE1_stenv_test, betaMSE_stenv = betaMSE_stenv,
                   predMSE2_pls_test = predMSE2_pls_test, predMSE1_pls_test = predMSE1_pls_test,
                   betaMSE_pls = betaMSE_pls, predMSE2_simul_pls_test = predMSE2_simul_pls_test,
                   predMSE1_simul_pls_test = predMSE1_simul_pls_test, betaMSE_simul_pls = betaMSE_simul_pls)
  return(result_to_save)
  
}
sim_num_a=1000





######################
######################
###TABLE 5.8
#####################
#####################

###FIRST PART Table 5.8


p=50;r=4;dx=10;dy=3;N=1000;M=1000



sim_result=apply(as.matrix(1:sim_num_a),1, sim_fun_a,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")
xtable(result.m)

###second PART Table 5.8

p=50;r=4;dx=10;dy=3;N=200;M=1000



sim_result=apply(as.matrix(1:sim_num_a),1, sim_fun_a,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")
xtable(result.m)


###third PART Table 5.8

p=50;r=4;dx=10;dy=3;N=100;M=1000


sim_result=apply(as.matrix(1:sim_num_a),1, sim_fun_a,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")
xtable(result.m)


###forth PART Table 5.8

p=50;r=4;dx=10;dy=3;N=70;M=1000



sim_result=apply(as.matrix(1:sim_num_a),1, sim_fun_a,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")
xtable(result.m)



###fifth PART Table 5.8

p=50;r=4;dx=10;dy=3;N=57;M=1000



sim_result=apply(as.matrix(1:sim_num_a),1, sim_fun_a,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")
xtable(result.m)


###sixth PART Table 5.8

p=50;r=4;dx=10;dy=3;N=50;M=1000

sim_result=apply(as.matrix(1:sim_num_a),1, sim_fun_a,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")
xtable(result.m)


###seven PART Table 5.8

p=50;r=4;dx=10;dy=3;N=25;M=1000


sim_result=apply(as.matrix(1:sim_num_a),1, sim_fun_a,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")
xtable(result.m)




