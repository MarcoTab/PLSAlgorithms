rm(list=ls())
library(Matrix)
library(MASS)
library(Renvlp)
library(xtable)
 
source("examples/chapter5/lib/ch5_misc_utils.R")
source("examples/chapter5/lib/rrenv.R")
source("examples/chapter5/lib/rrenvMU.R")
source("examples/chapter5/lib/env.apweights.R")
source("examples/chapter5/lib/rrenv.apweights.R")

# function to generate simulations at the ith repetition

sim_fun=function(i,p,r,dx,dy,N,M, simul_vers=1){
  if (simul_vers < 1 | simul_vers > 3) {
    return()
  }
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
  
  if (simul_vers == 1) {
    # case where two block go crazy and simul pls is better than simulenv at n=57 for case 1,2
    Omega = 500*diag(dx);    # material variation in X
    Phi = 0.05*diag(dy);      # material variation of Y|X
    Omega0 = diag(p-dx);    # immaterial variation in X
    Phi0 = 10*diag(r-dy);       # immaterial variation of Y|X
  } else if (simul_vers == 2) {
    # case where simultaneous env is better than pls, but two block is the best
    Omega = 1*diag(dx);    # material variation in X
    Phi = 0.1*diag(dy);      # material variation of Y|X
    Omega0 = diag(p-dx);    # immaterial variation in X
    Phi0 = 10*diag(r-dy);       # immaterial variation of Y|X
  } else {
    # case where simultaneous env is better than pls, and two block is not the best
    Omega = 50*diag(dx);    # material variation in X
    Phi = 0.01*diag(dy);      # material variation of Y|X
    Omega0 = diag(p-dx);    # immaterial variation in X
    Phi0 = 10*diag(r-dy);       # immaterial variation of Y|X
  }
  
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
sim_num=1000

 



######################
######################
###TABLE 5.3
#####################
#####################

###FIRST PART Table 5.3


p=50;r=4;dx=10;dy=3;N=1000;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)


###SECOND PART Table 5.3

p=50;r=4;dx=10;dy=3;N=200;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)




###THIRD PART Table 5.3

p=50;r=4;dx=10;dy=3;N=100;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)



###Fourth PART Table 5.3

p=50;r=4;dx=10;dy=3;N=70;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)



###FIFTH PART Table 5.3

p=50;r=4;dx=10;dy=3;N=57;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)


###sixt PART Table 5.3

p=50;r=4;dx=10;dy=3;N=50;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)





###seventh PART Table 5.3

p=50;r=4;dx=10;dy=3;N=25;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)







######################
######################
###TABLE 5.4
#####################
#####################

###FIRST PART Table 5.4


p=50;r=4;dx=40;dy=3;N=1000;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)


###SECOND PART Table 5.4

p=50;r=4;dx=40;dy=3;N=200;M=1000

 

sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)




###THIRD PART Table 5.4

p=50;r=4;dx=40;dy=3;N=100;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)



###Fourth PART Table 5.4

p=50;r=4;dx=40;dy=3;N=70;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)



###FIFTH PART Table 5.4

p=50;r=4;dx=40;dy=3;N=57;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)


###sixt PART Table 5.4

p=50;r=4;dx=40;dy=3;N=50;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)





###seventh PART Table 5.4

p=50;r=4;dx=40;dy=3;N=25;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)









######################
######################
###TABLE 5.5
#####################
#####################

###FIRST PART Table 5.5


p=30;r=50;dx=20;dy=2;N=1000;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)


###SECOND PART Table 5.5

p=30;r=50;dx=20;dy=2;N=200;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)




###THIRD PART Table 5.5

p=30;r=50;dx=20;dy=2;N=100;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)



###Fourth PART Table 5.5

p=30;r=50;dx=20;dy=2;N=85;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)



###FIFTH PART Table 5.5

p=30;r=50;dx=20;dy=2;N=70;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)


###sixt PART Table 5.5

p=30;r=50;dx=20;dy=2;N=50;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)





###seventh PART Table 5.5

p=30;r=50;dx=20;dy=2;N=25;M=1000



sim_result=apply(as.matrix(1:sim_num),1, sim_fun,p,r,dx,dy,N,M)
result.m=matrix(apply(sim_result,1,mean,na.rm=T),nrow=3)
row.names(result.m)=c("predMSE1","predMSE2","betaMSE")
colnames(result.m)=c("ols",'twoblock','stenv',"pls","simul_pls")

xtable(result.m)



