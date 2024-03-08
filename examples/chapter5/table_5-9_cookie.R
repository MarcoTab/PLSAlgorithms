##########################################
##########################################
# cookie dough data third row of Table 5.9
##########################################
##########################################

library(groc)
library(readr)
library('pls')
#needed functions
source("examples/chapter5/lib/ch5_misc_utils.R")



##data
data(cookie)
library(Renvlp)

#read data
X<-(as.matrix(cookie[,1:700]) )# extract NIR spectra
Y<-as.matrix(cookie[,701:704]) # extract constituents

 
#training
Xtrain<-X[1:40,] # extract training data
Ytrain<-Y[1:40,] # extract training data

#test
Xtest<-X[41:72,] # extract test data
Ytest<-Y[41:72,] # extract test data

##take away outliers

Ytrain=as.matrix(Ytrain[-23,])

Ytest =as.matrix(Ytest[-21,])

Xtrain=Xtrain[-23,]


Xtest=Xtest[-21,]




##we use only 64 predictors

Xtrain1=Xtrain[,seq(141,651,8)]
Xtest1=Xtest[,seq(141,651,8)]



CV_MSE_simul_pls_Y_test.v=NULL
CV_MSE_twoblock_Y_test.v=NULL
if.Xscale=FALSE
if.Yscale=FALSE
twoblock_CV_Y_MSE.m=NULL
twoblock_CV_Y_MSE.v=NULL
simul_pls_CV_Y_MSE.m=NULL
simul_pls_CV_Y_MSE.v=NULL
cor_twoblock_Y_Yhat=NULL
cor_simul_pls_Y_Yhat=NULL
nfold=10
set.seed(10)
set.seed(8)

#choose d
cv_result=cv_fun(X=Xtrain1,Y=Ytrain,nfold=nfold,if.Xscale=if.Xscale,if.Yscale=if.Yscale)
 

 

probar=NULL

#fit the test data
refit_result=refit_entire_training_and_get_testMSE_fun(Xtrain1,Ytrain,Xtest1,Ytest,
                                                       d=cv_result$d,d1=cv_result$d1,
                                                       d2=cv_result$d2,d_pls=cv_result$d_pls,
                                                       d_pls1=cv_result$d_pls1,
                                                       if.Xscale=if.Xscale,
                                                       if.Yscale=if.Yscale)




 

#two blocks

a1=mean(sqrt(apply(refit_result$twoblock_Y_MSE.v, 1, sum)))

#pls

a2=mean(sqrt(apply(refit_result$pls_Y_MSE.v, 1, sum)))

#simul


a3= mean(sqrt(apply(refit_result$simul_pls_Y_MSE.v,1,sum)))



#pls1


a4=mean(sqrt(apply(refit_result$pls1_Y_MSE.v, 1, sum)))

#print of the third row of Table 5.9
M=matrix(c(a1,a2,a3,a4),nrow=1)
colnames(M)<-c("two blocks","pls","sumult","pls1")

M
