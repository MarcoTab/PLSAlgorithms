rm(list=ls())
library('pls')
library(readr)
library(Renvlp) 
library(matrixcalc)

###################
# Skagerberg data #
###################

#### PPLS PARTIAL OLS and ENV partial table 6.2 

#### Columns 2, 3, 4

aux <- read_csv("data/chapter6/Process_data_Skagerberg1992 - Sheet1.csv")

# For reduction
X1=as.matrix(aux[,1:20])

# No reduction
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
    
    # training
    aux1=aux[-i,]
    
    # testing
    aux2=aux[i,]
    
    ### Matrix X to reduce
    AA=as.matrix(aux1[,1:20],ncol=20)
  
    ### Regression of the continuous X that we do not reduce
    A=lm(AA~Tw+S,data=aux1)
    
    # Save residuals
    aux1$SS=A$residuals 
    
    ### Predict on new data
    A_1=predict(A,newdata=aux2)
    
    # Center response in training 
    aux1$ynew=scale(aux1$Y,scale=FALSE)
    
    ## Do PLS on the training residuals
    m=plsr(ynew~SS,ncomp=d,data=aux1)
    
    # Perform regression on non-reduced variables
    j=lm(ynew~Tw+S,data=aux1)
    
    # Predict the new datapoint from non-reduced model
    j_2=predict(j,newdata=aux2)
    
    # Final prediction on new data
    yhat=apply(aux1$Y,2,mean)+j_2+as.numeric(as.matrix((aux2[,1:20]-A_1),nrow=1)%*%as.matrix(m$coefficients[,1:6,d],ncol=6))

    # Save PLS prediction
    pls_predic[i,d,]=yhat
    
    a=lm(Y~T1+T2+T3+T4+T5+T6+T7+T8+T9+T10+T11+T12+T13+T14+T15+T16+T17+T18+T19+T20+Tw+S
         ,data=aux1)
    
    # Save OLS prediction
    ols_predic[i,d,]=predict(a,newdata=aux2)
    
    M=eppls(X1[-i,],X2[-i,],Y[-i,],d)
    
    # Save Envelope prediction
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

table6_2 <- matrix(c(1:15), ncol=5, byrow=TRUE)

table6_2[1,] <- c('OLS', 'PLS (P.Predictor)', 'Envelope (P.Predictor)', 'PLS (P.Response)', 'Envelope (P.Response)')
rownames(table6_2) <- c('Method', 'No. components', 'predRMSE')

table6_2[3,1] = p.ols[1]
table6_2[2,1] = "p1 = 20"
table6_2[2,2] = which.min(p.pls)
table6_2[3,2] = p.pls[which.min(p.pls)]
table6_2[2,3] = which.min(p.env)
table6_2[3,3] = p.env[which.min(p.env)]

#######

#########COLUMNS 5,6 

# For reduction
X1=as.matrix(aux[,1:20])

# No reduction
X2=as.matrix(aux[,21:22])


Y=as.matrix(log(aux[,23:28]))
Y=Y%*%diag(diag(cov(Y))^(-1/2))

aux$Y=Y

fit_penv=array(NaN, dim=c(nrow(X1),ncol(Y),ncol(X1)))
error_penv=array(NaN, dim=c(nrow(X1),ncol(Y),ncol(X1)))

fit_prpls=array(NaN, dim=c(nrow(X1),ncol(Y),ncol(X1)))
error_prpls=array(NaN, dim=c(nrow(X1),ncol(Y),ncol(X1)))


##### Do envelope for different d
for (d in 1:ncol(Y)){
 for (i in 1:nrow(X1)){
   fit <- penv(X1[-i,],X2[-i,],Y[-i,], d)
   
   fit_penv[i,,d]=colMeans(Y[-i,])+(X1[i,]-colMeans(X1[-i,]))%*%t(fit$beta1)+(X2[i,]-colMeans(X2[-i,]))%*%t(fit$beta2)
   error_penv[i,,d]=(Y[i,]-colMeans(Y[-i,])-(X1[i,]-colMeans(X1[-i,]))%*%t(fit$beta1)-(X2[i,]-colMeans(X2[-i,]))%*%t(fit$beta2))^2
   
 }
}


#### PPLS, for d=1 we do it separately

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


### PPLS for the rest of the values for d

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




# Table 6.2 columns 5,6

table6_2[2,4] = which.min(p.pls)
table6_2[3,4] = p.pls[which.min(p.pls)]
table6_2[2,5] = which.min(p.env)
table6_2[3,5] = p.env[which.min(p.env)]

table6_2




