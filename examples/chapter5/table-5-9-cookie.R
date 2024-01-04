
#####################
# cookie dough data #
#####################

library(groc)
source("~/Dropbox/more-pls/pls-simultaneous/computing/our_functions_Feb2nd2022.R")
data(cookie)
library(Renvlp)

X<-(as.matrix(cookie[,1:700]) )# extract NIR spectra
Y<-as.matrix(cookie[,701:704]) # extract constituents

# X=X%*%eigen(cov(X))$vectors[,1:20]

#X=scale(X,center=FALSE)
#Y=scale(Y,center=FALSE)
#a=sample(1:72, replace=FALSE)
#X=X[a,]
#Y=Y[a,]

Xtrain<-X[1:40,] # extract training data

Ytrain<-Y[1:40,] # extract training data
Xtest<-X[41:72,] # extract test data
Ytest<-Y[41:72,] # extract test data

##saco unos outliers
 
Ytrain=as.matrix(Ytrain[-23,])

Ytest =as.matrix(Ytest[-21,])
 
Xtrain=Xtrain[-23,]

 
Xtest=Xtest[-21,]




##usamos solo 64 predictores

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
cv_result=cv_fun(X=Xtrain1,Y=Ytrain,nfold=nfold,if.Xscale=if.Xscale,if.Yscale=if.Yscale)
 #d=4
#d1=7
#d2=2

### el $d es el two block
### el $d_pls1 es pls 1 elige uno para cada Y
### el $d1 es el de simultaneo el primero
### el $d2 es el simultaneo segundo
### el $d_pls es el comun


#d_pls=7, 4

#cv_result$d_pls1
#dim_X dim_Y
#1    5     1
#2     12     1
#3     12     1
#4     5     1
#> 
#> 
#> 
#> 
#>
#>
#>

probar=NULL

refit_result=refit_entire_training_and_get_testMSE_fun(Xtrain1,Ytrain,Xtest1,Ytest,
                                                       d=cv_result$d,d1=cv_result$d1,
                                                       d2=cv_result$d2,d_pls=cv_result$d_pls,
                                                       d_pls1=cv_result$d_pls1,
                                                       if.Xscale=if.Xscale,
                                                       if.Yscale=if.Yscale)




#> 
#> 
#> 
#> 



#two blocks

a1=mean(sqrt(apply(refit_result$twoblock_Y_MSE.v, 1, sum)))

#pls

a2=mean(sqrt(apply(refit_result$pls_Y_MSE.v, 1, sum)))

#simul


a3= mean(sqrt(apply(refit_result$simul_pls_Y_MSE.v,1,sum)))



#pls1


a4=mean(sqrt(apply(refit_result$pls1_Y_MSE.v, 1, sum)))


M=matrix(c(a1,a2,a3,a4),nrow=1)
colnames(M)<-c("two blocks","pls","sumult","pls1")

M


