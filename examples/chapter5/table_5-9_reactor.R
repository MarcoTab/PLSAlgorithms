rm(list=ls())
library('pls')
 source("our_functions_Feb2nd2022.R")


###################
# Skagerberg data #
###################
source("our_functions_Feb2nd2022.R")
library(readr)
Skagerber_data <- read_csv("Process_data_Skagerberg1992 - Sheet1.csv")
 
 
X=as.matrix(Skagerber_data[,1:22])
Y=as.matrix(log(Skagerber_data[,23:28]))
Y=Y%*%diag(diag(cov(Y))^(-1/2))
a=sample(1:nrow(X),replace=F)
a=c(45 ,38 ,23 ,27, 35 ,52, 11 ,51 ,29 ,37 , 8 ,48,  7, 12, 32, 28 ,
    21 ,54 ,33 ,25, 34 ,22 ,50, 20 ,41 ,16, 43,
  1 , 3 ,44 ,24, 17 , 9 ,18 ,56 , 4 ,31, 10 ,26, 14, 53, 40, 15, 19, 42 ,39
  ,30 , 2 ,46, 47,  6 ,49,  5, 13,
    36 ,55)
X=X[a,]
Y=Y[a,]

Xtrain<-X # extract training data
Ytrain<-Y # extract training data
Xtest<-X # extract test data, the original paper only calculate the in sample MSE
Ytest<-Y # extract test data
CV_MSE_simul_pls_Y_test.v=NULL
CV_MSE_twoblock_Y_test.v=NULL
if.Xscale=F
if.Yscale=F
twoblock_CV_Y_MSE.m=NULL
twoblock_CV_Y_MSE.v=NULL
simul_pls_CV_Y_MSE.m=NULL
simul_pls_CV_Y_MSE.v=NULL
pls_CV_Y_MSE.m=NULL
pls_CV_Y_MSE.v=NULL
cor_twoblock_Y_Yhat=NULL
cor_simul_pls_Y_Yhat=NULL
cor_pls_Y_Yhat=NULL
 nfold=10
 
 
 #elige los d

cv_result=cv_fun(X=Xtrain,Y=Ytrain,nfold=nfold,if.Xscale=if.Xscale,if.Yscale=F)

result_LOO.list=apply(as.matrix(1:nrow(X)),1,cal_LOOMSE_onefold_fun,myX=X,myY=Y,d=cv_result$d,
                      d1=cv_result$d1,d2=cv_result$d2,
                      d_pls=cv_result$d_pls,d_pls1=cv_result$d_pls1, if.Xscale=if.Xscale,
                      if.Yscale=if.Yscale)
 
 
 
 



 

for(i in 1:nrow(X)){
  if(i==1){
    twoblock_MSE_sum=sqrt(sum(result_LOO.list[[i]]$twoblock_Y_MSE.v))
    simul_pls_MSE_sum=sqrt(sum(result_LOO.list[[i]]$simul_pls_Y_MSE.v))
    iter_simul_pls_MSE_sum=sqrt(sum(result_LOO.list[[i]]$iter_simul_pls_Y_MSE.v))
    pls_MSE_sum=sqrt(sum(result_LOO.list[[i]]$pls_Y_MSE.v))
    pls1_MSE_sum=sqrt(sum(result_LOO.list[[i]]$pls1_Y_MSE.v))
  }
  twoblock_MSE_sum = twoblock_MSE_sum + sqrt(sum(result_LOO.list[[i]]$twoblock_Y_MSE.v))
  simul_pls_MSE_sum = simul_pls_MSE_sum + sqrt(sum(result_LOO.list[[i]]$simul_pls_Y_MSE.v))
  pls_MSE_sum = pls_MSE_sum + sqrt(sum(result_LOO.list[[i]]$pls_Y_MSE.v))
  pls1_MSE_sum = pls_MSE_sum + sqrt(sum(result_LOO.list[[i]]$pls1_Y_MSE.v))
}

 

 
 
cal_ols_fun=function(i){
  ols_no1_Result=lm(Y[-i,]~X[-i,])
  mse_ols_1data=(Y[i,]-colMeans(Y[-i,])-(X[i,]-colMeans(X[-i,]))%*%ols_no1_Result$coefficients[-1,])^2
  return(mse_ols_1data)
}

mse_ols_1data=matrix(0,ncol=6,nrow=nrow(X))
for (i in 1:nrow(X)){
  ols_no1_Result=lm(Y[-i,]~X[-i,])
  mse_ols_1data[i,]=(Y[i,]-colMeans(Y[-i,])-(X[i,]-colMeans(X[-i,]))%*%ols_no1_Result$coefficients[-1,])^2

  }


 

 



  


#####simultaneous envelope
library(Renvlp)

a_s_env=u.stenv(X,Y)  
error_bic=matrix(0,ncol=ncol(Y),nrow=nrow(Y))


for (i in 1:nrow(X)){
  fit <- stenv(X[-i,],Y[-i,], as.numeric(a_s_env$u.bic[1]),as.numeric(a_s_env$u.bic[2]))
  
  betahat <- fit$beta
  
  
  error_bic[i,]=(Y[i,]-colMeans(Y[-i,])-(X[i,]-colMeans(X[-i,]))%*%fit$beta)^2
  
}


 
  ### todos lso errores juntos


A=c(mean(apply(sqrt(mse_ols_1data),1,sum)),twoblock_MSE_sum/nrow(X),mean(apply(sqrt(error_bic),1,sum)),pls_MSE_sum/nrow(X),simul_pls_MSE_sum/nrow(X),pls1_MSE_sum/nrow(X))
names(A)=c("OLS","Two blocks","Simultaneous ENVELOPE","PLS","Simultaneous PLS","PLS1")

 
print(round(A,2))
 

