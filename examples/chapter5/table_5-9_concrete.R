###################
# CONCRTE  #
###################
library(readr)
source("our_functions_Feb2nd2022.R")

 concrete <- read_table2("concrete.txt", col_names = FALSE)


Xtrain=as.matrix(concrete[1:78, 1:7])
Ytrain=as.matrix(concrete[1:78, 8:10])

#a=sample(1:nrow(Xtrain),replace=F)

#X=Xtrain[a,]
#Y=Ytrain[a,]

Xtrain<-X # extract training data
Ytrain<-Y # extract training data
Xtest<-as.matrix(concrete[79:103, 1:7]) # extract test data, the original paper only calculate the in sample MSE
Ytest<-as.matrix(concrete[79:103, 8:10]) # extract test data
CV_MSE_simul_pls_Y_test.v=NULL
CV_MSE_twoblock_Y_test.v=NULL
if.Xscale=F
if.Yscale=F
twoblock_CV_Y_MSE.m=NULL
twoblock_CV_Y_MSE.v=NULL
simul_pls_CV_Y_MSE.m=NULL
simul_pls_CV_Y_MSE.v=NULL
cor_twoblock_Y_Yhat=NULL
cor_simul_pls_Y_Yhat=NULL
nfold=10
set.seed(43)

### encuentro los d en el training con 10 cv
cv_result=cv_fun(X=Xtrain,Y=Ytrain,nfold=nfold,if.Xscale=if.Xscale,if.Yscale=if.Yscale)


 

##refiero con el trainnign con esos d y miro en el test
refit_result=refit_entire_training_and_get_testMSE_fun(Xtrain,Ytrain,Xtest,Ytest,
                                                       d=cv_result$d,d1=cv_result$d1,
                                                       d2=cv_result$d2,d_pls=cv_result$d_pls,d_pls1=cv_result$d_pls1,if.Xscale=if.Xscale,
                                                       if.Yscale=if.Yscale)



###error para TWO BLOCK
a1=mean(apply(refit_result$twoblock_Y_MSE.v, 2, mean))
 
a1=mean(sqrt(apply(refit_result$twoblock_Y_MSE.v, 1, sum)))
## error para X-pls

a2=mean(apply(refit_result$pls_Y_MSE.v, 2, mean))

a2=mean(sqrt(apply(refit_result$pls_Y_MSE.v, 1, sum)))


##error para XY-pls

a3=mean(apply(refit_result$simul_pls_Y_MSE.v, 2, mean))

a3=mean(sqrt(apply(refit_result$simul_pls_Y_MSE.v, 1, sum)))


##error para X-pls1

a4=mean(apply(refit_result$pls1_Y_MSE.v, 2, mean))

a4=mean(sqrt(apply(refit_result$pls1_Y_MSE.v, 1, sum)))






 
#envelope

library(Renvlp)
sim_env=u.stenv(Xtrain,Ytrain)
sim_env=stenv(Xtrain,Ytrain,as.numeric(sim_env$u.bic[1]),as.numeric(sim_env$u.bic[2]))




tmp0 <- (Ytest-(as.matrix(Xtest) %*% (sim_env$beta)+matrix(1, 25, 1) %*% t(sim_env$mu)))^2
 

a7=mean(apply(tmp0, 2, mean))

a7=mean(sqrt(apply(tmp0, 1, sum)))


###OLS

  
 



ols_result=refit_result=refit_entire_training_and_get_testMSE_fun(Xtrain,Ytrain,Xtest,Ytest,
                                                                  d=cv_result$d,d1=ncol(Xtrain),
                                                                  d2=ncol(Ytrain),d_pls=c(ncol(Xtrain),ncol(Ytrain)),d_pls1=cv_result$d_pls1,if.Xscale=if.Xscale,
                                                                  if.Yscale=if.Yscale)



#OLS

 
a6=mean(apply(ols_result$simul_pls_Y_MSE.v, 2, mean))

a6=mean(sqrt(apply(refit_result$pls_Y_MSE.v, 1, sum)))

M=matrix(c(a6,a1,a7,a2,a3,a4),nrow=1)
colnames(M)<-c("OLS", "two blocks","env","pls","sumult","pls1")

M

