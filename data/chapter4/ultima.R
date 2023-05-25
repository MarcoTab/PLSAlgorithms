#we are doing n=p/2
library(pls)
library(spls)
#library(plsr)
library(MASS)
library(gglasso)

source('~/Dropbox/PLS-chemometrics/simus/solo_pls-simpls.R')
data(yeast)
ss=10

JJ=seq(50,1200,100)
nreps=40
PLSS=array(0,dim=c(nreps,length(JJ),6))
PLSS2=array(0,dim=c(nreps,length(JJ),6))
imp=matrix(0,ncol=length(JJ),nrow=nreps)

SPLSS=matrix(0,ncol=length(JJ),nrow =nreps)
SPLSS2=matrix(0,ncol=length(JJ),nrow=nreps)

OLS_PC=array(0,dim=c(nreps,length(JJ),6))
OLS_PC2=array(0,dim=c(nreps,length(JJ),6))


OLS_PLS=array(0,dim=c(nreps,length(JJ),6))
OLS_PLS2=array(0,dim=c(nreps,length(JJ),6))

PC=array(0,dim=c(nreps,length(JJ),6))
PC2=array(0,dim=c(nreps,length(JJ),6))

LASSO=matrix(0,ncol=length(JJ),nrow =nreps)
LASSO2=matrix(0,ncol=length(JJ),nrow=nreps)

ENV_PC=array(0,dim=c(nreps,length(JJ),6))
ENV_PC2=array(0,dim=c(nreps,length(JJ),6))


ENV_PLS=array(0,dim=c(nreps,length(JJ),6))
ENV_PLS2=array(0,dim=c(nreps,length(JJ),6))

 
error_env_pc_outsample=matrix(0,ncol=length(JJ),nrow=nreps)
error_env_pc_insample=matrix(0,ncol=length(JJ),nrow=nreps)

error_env_pls_outsample=matrix(0,ncol=length(JJ),nrow=nreps)
error_env_pls_insample=matrix(0,ncol=length(JJ),nrow=nreps)

error_ols_pc_outsample=matrix(0,ncol=length(JJ),nrow=nreps)
error_ols_pc_pc_insample=matrix(0,ncol=length(JJ),nrow=nreps)

error_pls_outsample=matrix(0,ncol=length(JJ),nrow=nreps)
error_pls_insample=matrix(0,ncol=length(JJ),nrow=nreps)

error_ols_pls_outsample=matrix(0,ncol=length(JJ),nrow=nreps)
error_ols_pls_insample=matrix(0,ncol=length(JJ),nrow=nreps)

 
for (i in 1:length(JJ)){
  
  p=JJ[i]
  qq=c(0,40,floor(p-p^(2/3)),floor(p-p^(1/2)),p-2,p)
  for (jj in 1:6){
    
    n=floor(p/3)
    q=qq[jj]
    
    for (hh in 1:nreps){
      
      cat("p = ", p, "\n")
      
      
      
      if (q!=p){
        H1=mvrnorm(1,rep(0,n),25*diag(rep(1,n)))
        H2=mvrnorm(1,rep(0,n),25*diag(rep(1,n)))
        if(q!=0){H3=mvrnorm(1,rep(0,n),25*diag(rep(1,n)))}
        
        ns=c(0,floor((p-q)/2),p-q,p)
        X=matrix(0,ncol=p,nrow=n)
        H1aux=matrix(rep(H1,ns[2]),ncol=ns[2])
        H2aux=matrix(rep(H2,ns[3]-ns[2]),ncol=ns[3]-ns[2])
        if(q!=0){
          H3aux=matrix(rep(H3,ns[4]-ns[3]),ncol=ns[4]-ns[3])}
        
        if(q!=0){X =cbind(H1aux+matrix(rnorm(ns[2]*n),ncol=ns[2]),H2aux+matrix(rnorm((ns[3]-ns[2])*n),ncol=ns[3]-ns[2]),H3aux+matrix(rnorm((ns[4]-ns[3])*n),ncol=ns[4]-ns[3]))}
        
        if(q==0){X =cbind(H1aux+matrix(rnorm(ns[2]*n),ncol=ns[2]),H2aux+matrix(rnorm((ns[3]-ns[2])*n),ncol=ns[3]-ns[2]))}
        
        
        y=rnorm(n,3*H1-4*H2,ss)
        
        
        if(q!=0){Xnew =cbind(H1aux+matrix(rnorm(ns[2]*n),ncol=ns[2]),H2aux+matrix(rnorm((ns[3]-ns[2])*n),ncol=ns[3]-ns[2]),H3aux+matrix(rnorm((ns[4]-ns[3])*n),ncol=ns[4]-ns[3]))}
        if(q==0){Xnew =cbind(H1aux+matrix(rnorm(ns[2]*n),ncol=ns[2]),H2aux+matrix(rnorm((ns[3]-ns[2])*n),ncol=ns[3]-ns[2]))}
        
        
        ynew=rnorm(n,3*H1-4*H2,ss)
        yeast$x=X
        yeast$y=y}
      
      
      
      if (q==p){
        
        
        X=matrix(0,ncol=p,nrow=n)
        H3aux=matrix(rep(H3,q),ncol=q)
        
        X =H3aux+matrix(rnorm((q)*n),ncol=q)
        
        
        
        y=rnorm(n,0,ss)
        
        
        Xnew =H3aux+matrix(rnorm((q)*n),ncol=q)
        
        
        ynew=rnorm(n,0,ss)
        yeast$x=X
        yeast$y=y}
      
      
      
      
      
      
      imp[hh,i]=sum(eigen(cov(yeast$x))$values[1:2])/sum(eigen(cov(yeast$x))$values)
     # yeast$y=as.matrix(yeast$y)
      j=1
      k=2
       
      
      a1=75/(1+25*(p-q)/2)
      a2=100/(1+25*(p-q)/2)
      beta.true=c(rep(a1,(p-q)/2),rep(-a2,(p-q)/2),rep(0,q))
      
      
      
      A1=estimate.beta.pls(yeast$x, yeast$y, 2,plsS=TRUE,plsN=TRUE, env_PLS=FALSE,quiet=TRUE,how.much=.975)
      # order(abs(t(A1[,3])%*%eigen(cov(yeast$x))$vectors),decreasing=TRUE)[1:10]
      
      
       
      ####pls
      
      
      xpred=(yeast$x%*%A1[,1])
      mirar=lm(yeast$y~xpred)
      error_pls_insample[hh,i]=mean(mirar$residuals^2)
      
      xnew=Xnew%*%A1[,1]
      error_pls_outsample[hh,i]=mean((ynew-mirar$coefficients[1]-mirar$coefficients[2]*xnew)^2)
      
      
      PLSS[hh,i,jj]=as.numeric(t(beta.true-A1[,1])%*%cov(yeast$x)%*%(beta.true-A1[,1]))
      
      PLSS2[hh,i,jj]=as.numeric(t(beta.true-A1[,1])%*%(beta.true-A1[,1]))
      
      
      #####env pls
      
      xpred=(yeast$x%*%A1[,1])
      mirar=lm(yeast$y~xpred)
      error_env_pls_insample[hh,i]=mean(mirar$residuals^2)
      
      xnew=Xnew%*%A1[,1]
      error_env_pls_outsample[hh,i]=mean((ynew-mirar$coefficients[1]-mirar$coefficients[2]*xnew)^2)
      
      
      ENV_PLS[hh,i,jj]=as.numeric(t(beta.true-A1[,3])%*%cov(yeast$x)%*%(beta.true-A1[,3]))
      
      ENV_PLS2[hh,i,jj]=as.numeric(t(beta.true-A1[,3])%*%(beta.true-A1[,3]))
      
      
      ###ols_pls
      
      xpred=(yeast$x%*%A1[,4])
      mirar=lm(yeast$y~xpred)
      error_ols_pls_insample[hh,i]=mean(mirar$residuals^2)
      
      xnew=Xnew%*%A1[,4]
      error_ols_pls_outsample[hh,i]=mean((ynew-mirar$coefficients[1]-mirar$coefficients[2]*xnew)^2)
      
      
      OLS_PLS[hh,i,jj]=t(beta.true-A1[,4])%*%cov(yeast$x)%*%(beta.true-A1[,4])
      
      OLS_PLS2[hh,i,jj]=t(beta.true-A1[,4])%*%(beta.true-A1[,4])
      
      
       #env pc 
       
      xpred=(yeast$x%*%A1[,5])
      mirar=lm(yeast$y~xpred)
      error_env_pc_insample[hh,i]=mean(mirar$residuals^2)
      
      xnew=Xnew%*%A1[,5]
      error_env_pc_outsample[hh,i]=mean((ynew-mirar$coefficients[1]-mirar$coefficients[2]*xnew)^2)
      
      
      ENV_PC[hh,i,jj]=t(beta.true-A1[,5])%*%cov(yeast$x)%*%(beta.true-A1[,5])
      
      ENV_PC2[hh,i,jj]=t(beta.true-A1[,5])%*%(beta.true-A1[,5])
      
      
      #ols_PC
      
      #env pc 
      
      xpred=(yeast$x%*%A1[,6])
      mirar=lm(yeast$y~xpred)
      error_env_pc_insample[hh,i]=mean(mirar$residuals^2)
      
      xnew=Xnew%*%A1[,6]
      error_env_pc_outsample[hh,i]=mean((ynew-mirar$coefficients[1]-mirar$coefficients[2]*xnew)^2)
      
      
      OLS_PC[hh,i,jj]=t(beta.true-A1[,6])%*%cov(yeast$x)%*%(beta.true-A1[,6])
      
      OLS_PC2[hh,i,jj]=t(beta.true-A1[,6])%*%(beta.true-A1[,6])
      
       
      
      
       
    }
    
  }
  
}


M=apply(PLSS,c(2,3),mean)

rownames(M)=JJ
colnames(M)=c(0,40,'floor(p-p^(2/3))','floor(p-p^(1/2))','p-2','p')

write.table(M,file='Sigmanormnnover2.txt')



N=apply(PLSS2,c(2,3),mean)
rownames(N)=JJ
colnames(N)=c(0,40,'floor(p-p^(2/3))','floor(p-p^(1/2))','p-2','p')

write.table(N,file='normnnover2.txt')

plot(Sigmanormnnover2[,1],Sigmanormnnover2[,2],type='l',ylim=c(min(Sigmanormnnover2[,2:7]),max(Sigmanormnnover2[,2:7])))
lines(Sigmanormnnover2[,1],Sigmanormnnover2[,3],type='l',ylim=c(min(Sigmanormnnover2[,2:7]),max(Sigmanormnnover2[,2:7])))
lines(Sigmanormnnover2[,1],Sigmanormnnover2[,4],type='l',ylim=c(min(Sigmanormnnover2[,2:7]),max(Sigmanormnnover2[,2:7])))
lines(Sigmanormnnover2[,1],Sigmanormnnover2[,5],type='l',ylim=c(min(Sigmanormnnover2[,2:7]),max(Sigmanormnnover2[,2:7])))
lines(Sigmanormnnover2[,1],Sigmanormnnover2[,5],type='l',ylim=c(min(Sigmanormnnover2[,2:7]),max(Sigmanormnnover2[,2:7])))
lines(Sigmanormnnover2[,1],Sigmanormnnover2[,6],type='l',ylim=c(min(Sigmanormnnover2[,2:7]),max(Sigmanormnnover2[,2:7])))
lines(Sigmanormnnover2[,1],Sigmanormnnover2[,7],type='l',ylim=c(min(Sigmanormnnover2[,2:7]),max(Sigmanormnnover2[,2:7])))

library(latex2exp)
plot(JJ,apply(error_spls_insample,2,mean),type='l',xlab='Number of predictors',ylab=TeX("$||Y - \\hat{Y}||$"))
lines(JJ, apply(error_pls_outsample,2,mean),type='l',col='red')
lines(JJ, apply(error_env_outsample,2,mean),type='l',col='blue')
title("n=p/2 and q=30")


plot(JJ ,apply(SPLSS ,2,mean),type='l',xlab='Number of predictors',ylab=TeX("$||\\hat{beta} - beta||_{Sigma}$"))
lines(JJ , apply(PLSS ,2,mean),type='l',col='red')
lines(JJ, apply(EPLSS ,2,mean),type='l',col='blue')
title("n=p/2 and q=30")
