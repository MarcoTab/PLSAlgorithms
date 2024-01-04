#
#FIGURE 10.2 and 10.3
# This file produces all analyses from the article
#
#  
rm(list=ls())
library(pls)
library(Renvlp)
library(lavaan)
library(cSEM)
library(MASS)
library(chemometrics)
library(matrixpls)
library(CCA)
library(lavaan)
library(parallel)
library(pracma)
library(MASS)
library(psych)

source("examples/chapter10/rrenv.R")
source("examples/chapter10/rrr.R")
source("examples/chapter10/GE.R")
source("examples/chapter10/envMU.R")

 
set.seed(1)

#needs to be run with 100 and 1000
N <- 100

# Loadings are scaled so that errors have unit variances

Lrawx=c(.8,-.7,.6)
Lrawy=c(.8,-.7,.6)
 

Lambda <- matrix(c(Lrawx,rep(0,6),Lrawy), ncol = 2)

auxxx=Lrawx/as.numeric(sqrt(t(Lrawx)%*%(Lrawx)))
auxxy=Lrawy/as.numeric(sqrt(t(Lrawy)%*%(Lrawy)))
 


 
MODEL <- "
X =~ x1 + x2 + x3
Y =~ y1 + y2 + y3
Y ~ X
"


# Vary the population beta

BETAS <- seq(200/1500, 1000/1500,20/1500)

nrep=10
nose=array(NA,c(length(BETAS),7,nrep))


sigmacond=diag(rep(1,6))


for (h in 1:nrep){
  print("h")
  
  for (j in 1:length(BETAS)){
    Beta=BETAS[j]
    
    eta <- mvrnorm(N,c(0,0),matrix(c(1,Beta,Beta,1),2,2))
    
    data <- eta %*% t(Lambda) +  mvrnorm(N, rep(0,6),sigmacond )
     
    colnames(data) <- c("x1","x2","x3","y1","y2","y3")
    
    # CB | SEM estimates
     
    ml <- sem(MODEL,data)
    
     
    # Regression with PLS Mode A composites (what they called PLS)
    #MATRIXPLS
    
    pls_them1 <- matrixpls(cor(data), MODEL)
    
    # IDEAL TRUE Correlation between ideal composites (equivalent to standardized regression)
    ideal <- cor(lm(eta[,1] ~ data[,1:3])$fitted,
                 lm(eta[,2] ~ data[,4:6])$fitted) 
    
    
    
    ###TRUE PLS dimension 1
    X=data[,1:3]
    Y=data[,4:6]
    X1=X-(rep(1,dim(X)[1]))%*%t(apply(X,2,mean))
    Y1=X-(rep(1,dim(X)[1]))%*%t(apply(Y,2,mean))
    
    X= scale(X,center = TRUE, scale = FALSE)
    Y= scale(Y,center = TRUE, scale = FALSE)
    
    
    A=X%*%plsr(Y~X,ncomp=1)$projection[,1]
    B=Y%*%plsr(X~Y,ncomp=1)$projection[,1]
    
    
    
    A=A%*%sqrtm(cov(A))$Binv
    B=B%*%sqrtm(cov(B))$Binv
    
    au1=svd(cov(A,B))
    au2=svd(cov(B,A))
    
    true_pls_true=abs(cor(A%*%au1$u[,1],B%*%au2$u[,1]))
    
    
     
    
    
    
     
     
    ###PLS |SEM (pls|sem program by us...)
    X=data[,1:3]
    Y=data[,4:6]
    
    X= scale(X,center = TRUE, scale = FALSE)
    Y= scale(Y,center = TRUE, scale = FALSE)
    
    X= X%*%diag(diag(cov(X)*(N-1)/N)^(-1/2)) 
    Y= Y%*%diag(diag(cov(Y)*(N-1)/N)^(-1/2))
    
    au1=svd(cov(X,Y))
    au2=svd(cov(Y,X))
     
    m1=au1$u[,1]%*%cov(X)%*%au1$u[,1]
    m2=au2$u[,1]%*%cov(Y)%*%au2$u[,1]
    pls_SEM=abs(au1$d[1]/sqrt(m1*m2))
    
    
    
    
    
    ###MLE
    
    
    X=data[,1:3]
    Y=data[,4:6]
    
    
    X= scale(X,center = TRUE, scale = FALSE)
    Y= scale(Y,center = TRUE, scale = FALSE)
    
    
    r1<-rrr(X , Y, k=0,type = "cva",rrr.plot=FALSE)
    r2<-rrr(Y , X, k=0,type = "cva",rrr.plot=FALSE)
    
    abs(sum(diag(r1$C[[1]]%*%t(r2$C[[1]]))))^(1/2)
    
    CC=sqrtm(var(Y))$Binv%*%cov(Y,X)%*%sqrtm(var(X))$Binv
    auxiliar=svd(CC)
    betasiXY=auxiliar$d[1]*sqrtm(var(Y))$B%*%auxiliar$u[,1]%*%t(auxiliar$v[,1])%*%sqrtm(var(X))$Binv
    
    CC=sqrtm(var(X))$Binv%*%cov(X,Y)%*% sqrtm(var(Y))$Binv
    auxiliar=svd(CC)
    betasiYX=auxiliar$d[1]*sqrtm(var(X))$B%*%auxiliar$u[,1]%*%t(auxiliar$v[,1])%*%sqrtm(var(Y))$Binv
    
    ans=abs(sum(diag(betasiYX%*%(betasiXY))))^(1/2)
    
    
 
    
    
    
    
    #####ENV
    
    X=data[,1:3]
    Y=data[,4:6]
    
    
    X= scale(X,center = TRUE, scale = FALSE)
    Y= scale(Y,center = TRUE, scale = FALSE)
    
    GammaY=env(X, Y, 1, asy = TRUE)$Gamma
    
    
    GammaX=xenv(X, Y, 1, asy = TRUE)$Gamma
    
    GammaX=Lambda[1:3,1]
    GammaY=Lambda[1:3,1]
    nreps=10
    for (k in 1:nreps){
      Xnew=X%*%GammaX
      Ynew=Y%*%GammaY
      GammaX=xenv(X, Ynew, 1, asy = TRUE)$Gamma
      
      GammaY=env(Xnew, Y, 1, asy = TRUE)$Gamma
    }
    
    
    Ynew=Y%*%GammaY
    Xnew=X%*%GammaX 
    
    
    
    A=GammaY
    B=GammaX%*%cov(Xnew,Ynew)
    
    aux=xenv(X, Ynew, 1, asy = TRUE)
    SigmaX=xenv(X, Ynew, 1, asy = TRUE)$SigmaX
    
    SigmaXinv=GammaX%*%solve(aux$Omega)%*%t(GammaX)+aux$Gamma0%*%solve(aux$Omega0)%*%t(aux$Gamma0)
    aux=env(Xnew, Y, 1, asy = TRUE)
    
    SigmaY=env(Xnew, Y, 1, asy = TRUE)$Sigma+ GammaY%*%cov(Ynew,Xnew)%*%solve(xenv(X, Ynew, 1, asy = TRUE)$Omega)%*%cov(Xnew,Ynew)%*%t(GammaY)
    
    
    Meta=sqrt(t(A)%*%solve(SigmaY)%*%A)
    Mxi=sqrt(t(B)%*%(SigmaXinv)%*%B)
    
    res= Meta*Mxi
    
    
    
    
    
    
    
    nose[j,,h]=c(pls_them1[1], ideal,inspect(ml,"std")$beta[2],res,ans,
                pls_SEM,true_pls_true)
  }
}




estimates=apply(nose,1:2,function(x){mean(x,na.rm=TRUE)})


 


######FIGURE 10-3

 

 
 

par(mfrow=c(3,2))
library(latex2exp)

plot(estimates[,2],estimates[,5],xlab=TeX("$\\Psi $"),ylab="",cex.lab=2,col='black',cex.axis=1)
abline(0,1)
points(estimates[,2],estimates[,1],col="black")
title(ylab=TeX("$\\widehat{\\Psi }$"), line=2,cex.lab=2)
title("MLE")



plot(estimates[,2],estimates[,1],xlab=TeX("$\\Psi $"),ylab="",cex.lab=2,col='black',cex.axis=1)
abline(0,1)
title(ylab=TeX("$\\widehat{\\Psi }$"), line=2,cex.lab=2)
title("Matrixpls")


plot(estimates[,2],estimates[,4],xlab=TeX("$\\Psi$"),ylab="",cex.lab=2,col='black',cex.axis=1)
abline(0,1)
title(ylab=TeX("$\\widehat{\\Psi }$"), line=2,cex.lab=2)
title("ENV")


 

plot(BETAS,abs(estimates[,3]),xlab=TeX("$cor( xi, eta)$"),ylab="",cex.lab=2,col='black',cex.axis=1)
abline(0,1)
title(ylab=TeX("$\\widehat{cor}( xi, eta))$"), line=2,cex.lab=2)
title("CB | SEM")

 

 


plot(estimates[,2],estimates[,6],xlab=TeX("$\\Psi $"),ylab="",cex.lab=2,col='black',cex.axis=1)
abline(0,1)
title(ylab=TeX("$\\widehat{\\Psi }$"), line=2,cex.lab=2)
title("PLS | SEM")

plot(estimates[,2],estimates[,7],xlab=TeX("$\\Psi$"),ylab="",cex.lab=2,col='black',cex.axis=1)
abline(0,1)
title(ylab=TeX("$\\widehat{\\Psi}$"), line=2,cex.lab=2)
title("PLS")
 

