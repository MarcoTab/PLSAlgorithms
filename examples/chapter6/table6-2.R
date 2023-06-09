rm(list=ls())
library('pls')
library(EPPLS)
library(readr)
library(Renvlp) 

 
library(matrixcalc)

library(matrixcalc)


###################
# Skagerberg data #
###################

####PPLS PARTIAL OLS and ENV partial table 6.2 








####columns 2, 3, 4

aux <- read_csv("~/Dropbox/more-pls/pls-simultaneous/computing/datasets/Process_data_Skagerberg1992 - Sheet1.csv")



#Skagerber_data <- read_csv("Process_data_Skagerberg1992 - Sheet1.csv")

#para reducir
X1=as.matrix(aux[,1:20])

#sin reducir
X2=as.matrix(aux[,21:22])


Y=as.matrix(log(aux[,23:28]))
Y=Y%*%diag(diag(cov(Y))^(-1/2))

aux$Y=Y





 
Dmax=18
 
n=nrow(aux)
pls_predic=array(NaN,dim=c(n,Dmax,6) )
ols_predic=array(NaN,dim=c(n,Dmax,6) )
env_predic=array(NaN,dim=c(n,Dmax,6) )
 


####PPLS PARTIAL OLS and ENV partial
for (d in 1:Dmax){
  for (i in 1:n){
    
    #training
    aux1=aux[-i,]
    
    #testing
    aux2=aux[i,]
    
    ###matriz X que se reduce
    AA=as.matrix(aux1[,1:20],ncol=20)
    
     
    ###regresion de las continuas X en las que no reduzco
    A=lm(AA~Tw+S,data=aux1)
    
    #guardo los residuos
    aux1$SS=A$residuals 
    
    
    
    ###predigo en el nuevo dato
    A_1=predict(A,newdata=aux2)
    
    #cnetro la resputa en los de training
    aux1$ynew=scale(aux1$Y,scale=FALSE)
    
    ## hago pls en los residuos en el training
    m=plsr(ynew~SS,ncomp=d,data=aux1)
    
    #hago la regresion en las variables no reducidas
    j=lm(ynew~Tw+S,data=aux1)
    
    #predigo las no reducidas en el nuevo dato
    j_2=predict(j,newdata=aux2)
    
    #prediccion final en el nuevo datos
    yhat=apply(aux1$Y,2,mean)+j_2+as.numeric(as.matrix((aux2[,1:20]-A_1),nrow=1)%*%as.matrix(m$coefficients[,1:6,d],ncol=6))
    
    
    
    
    
    pls_predic[i,d,]=yhat
    
    
     
    
    
    
    a=lm(Y~T1+T2+T3+T4+T5+T6+T7+T8+T9+T10+T11+T12+T13+T14+T15+T16+T17+T18+T19+T20+Tw+S
         ,data=aux1)
    
    ols_predic[i,d,]=predict(a,newdata=aux2)
    
    
    
    
    
     M=eppls(X1[-i,],X2[-i,],Y[-i,],d)
    
   env_predic[i,d,]=t(M$beta1)%*%(X1[i,]-M$mu1)+t(M$beta2)%*%(X2[i,]-M$mu2)+M$muY
    
    
    
    
   
    
  }
}







 
 
 
 
 p.ols=NULL
 p.pls=NULL
 p.env=NULL
 
 for (d in 1:Dmax){
 
     
     p.ols[d]= mean(sqrt(apply((aux$Y-ols_predic[,d,])^2,1,mean)))
     p.pls[d]=mean(sqrt(apply((aux$Y-pls_predic[,d,])^2,1,mean)))
     p.env[d]=mean(sqrt(apply((aux$Y-env_predic[,d,])^2,1,mean)))
   
 }
 
 
 ##table 6.2 columns 2, 3,4
 
 p.ols
 d_pls=which.min(p.pls)
 d_env=which.min(p.env)
 
 d_pls
 d_env
 
 p.ols[1]
 p.pls[d_pls]
 p.env[d_env]
 
 
 
 #######
 
 #########COLUMNS 5,6 
 
 
 
 
 #para reducir
 X1=as.matrix(aux[,1:20])
 
 #sin reducir
 X2=as.matrix(aux[,21:22])
 
 
 Y=as.matrix(log(aux[,23:28]))
 Y=Y%*%diag(diag(cov(Y))^(-1/2))
 
 aux$Y=Y
 
 fit_penv=array(NaN, dim=c(nrow(X1),ncol(Y),ncol(X1)))
 error_penv=array(NaN, dim=c(nrow(X1),ncol(Y),ncol(X1)))
 
 fit_prpls=array(NaN, dim=c(nrow(X1),ncol(Y),ncol(X1)))
 error_prpls=array(NaN, dim=c(nrow(X1),ncol(Y),ncol(X1)))
 
 
 ##### voy a hacerlo para diferenes d el envelope
 for (d in 1:ncol(Y)){
   for (i in 1:nrow(X1)){
     fit <- penv(X1[-i,],X2[-i,],Y[-i,], d)
     
     fit_penv[i,,d]=colMeans(Y[-i,])+(X1[i,]-colMeans(X1[-i,]))%*%t(fit$beta1)+(X2[i,]-colMeans(X2[-i,]))%*%t(fit$beta2)
     error_penv[i,,d]=(Y[i,]-colMeans(Y[-i,])-(X1[i,]-colMeans(X1[-i,]))%*%t(fit$beta1)-(X2[i,]-colMeans(X2[-i,]))%*%t(fit$beta2))^2
     
   }
 }
 
 
 
 
 
 
 ####ppls lo hago para d =1 apartelo hago para d=1 aparte 
 
 for (i in 1:nrow(X1)){
   X1c=scale(X1[-i,],scale=FALSE)
   X2c=scale(X2[-i,],scale=FALSE)
   
   resY2=as.matrix(lm(Y[-i,]~X2[-i,])$residuals)
   res12=as.matrix(lm(X1[-i,]~X2[-i,])$residuals)
   
   
   
   m=plsr(res12~resY2,ncomp=1)
   
   
   A=Y[-i,]%*%t(t(m$projection))
   
   ss=lm(A~X1[-i,]+X2[-i,])
   
   
   auxx=Y[-i,]-X1[-i,]%*%ss$coefficients[2:21]%*%(t(m$projection))
   
   auxx1=lm(auxx~X2[-i,])$coefficients[2:3,]
   
   
   
   fit_prpls[i,,1]=colMeans(Y[-i,])+(X1[i,]-colMeans(X1[-i,]))%*%ss$coefficients[2:21]%*%(t(m$projection))+(X2[i,]-colMeans(X2[-i,]))%*%(auxx1)
   error_prpls[i,,1]=(Y[i,]-fit_prpls[i,,1])^2
   
 }
 
 
 ###ppls para el resto de los d
 
 
 
 
 
 for (d in 2:ncol(Y)){
   for (i in 1:nrow(X1)){
     
     resY2=as.matrix(lm(Y[-i,]~X2[-i,])$residuals)
     res12=as.matrix(lm(X1[-i,]~X2[-i,])$residuals)
     
     
     
     m=plsr(res12~resY2,ncomp=d)
     
     
     A=Y[-i,]%*%t(t(m$projection))
     
     ss=lm(A~X1[-i,]+X2[-i,])
     
     
     auxx=Y[-i,]-X1[-i,]%*%ss$coefficients[2:21,1:d]%*%(t(m$projection))
     
     auxx1=lm(auxx~X2[-i,])$coefficients[2:3,]
     
     
     
     fit_prpls[i,,d]=colMeans(Y[-i,])+(X1[i,]-colMeans(X1[-i,]))%*%ss$coefficients[2:21,1:d]%*%(t(m$projection))+(X2[i,]-colMeans(X2[-i,]))%*%(auxx1)
     error_prpls[i,,d]=(Y[i,]-fit_prpls[i,,d])^2
     
   }
 }
 
 
 p.ols=NULL
 p.pls=NULL
 p.env=NULL
 
 for (d in 1:ncol(Y)){
   
   
   p.env[d]=mean(sqrt(apply(error_penv[,,d],1,mean)))
   p.pls[d]=mean(sqrt(apply(error_prpls[,,d],1,mean)))
 }
 
 
 
 
 
 
 d_pls=which.min(p.pls)
 d_env=which.min(p.env)
 
 
 ##table 6.2  columns 5,6 
 
 d_pls
 d_env
 
 p.pls[d_pls]
 p.env[d_env]
 
 
 
 
  