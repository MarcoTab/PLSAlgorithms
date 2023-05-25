#FIGURE 4.7
library(pls)
library(spls)
#library(plsr)
library(MASS)
library(gglasso)

set.seed(100)
 
ss=10

JJ=seq(150,1150,100)
nreps=40
PLSS=array(0,dim=c(nreps,length(JJ),4))
 
yeast = NULL
 
for (i in 1:length(JJ)){
  
  p=JJ[i]
  
  ## qq is how many predictors do not make contribution
  ##qq=0 all make contribution means Case 1 or d=p
  ##qq=40 all except 40 make contribution (almost Case 1
  ##qq=floor(p-p^(2/3)) means p^(2/3) makes contribution Case 2
  ##qq=p-2 only two make contribution 

  qq=c(0,40,p-floor(p^(2/3)),p-2) 
  for (jj in 1:4){
    
    n=floor(p/3)
    q=qq[jj]
    
    for (hh in 1:nreps){
      
      cat("p = ", p, "q = ", q, "\n")
      
      
      
      if (q!=p) {
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
        
        # TODO Para que sirve esto??
        #if(q!=0){Xnew =cbind(H1aux+matrix(rnorm(ns[2]*n),ncol=ns[2]),H2aux+matrix(rnorm((ns[3]-ns[2])*n),ncol=ns[3]-ns[2]),H3aux+matrix(rnorm((ns[4]-ns[3])*n),ncol=ns[4]-ns[3]))}
        #if(q==0){Xnew =cbind(H1aux+matrix(rnorm(ns[2]*n),ncol=ns[2]),H2aux+matrix(rnorm((ns[3]-ns[2])*n),ncol=ns[3]-ns[2]))}
        
        
        #ynew=rnorm(n,3*H1-4*H2,ss)
        yeast$x=X
        yeast$y=y
      }
      
      
      
      if (q==p){
        
        
        X=matrix(0,ncol=p,nrow=n)
        H3=mvrnorm(1,rep(0,n),25*diag(rep(1,n)))
        H3aux=matrix(rep(H3,q),ncol=q)

        X =H3aux+matrix(rnorm((q)*n),ncol=q)
        
        
        
        y=rnorm(n,0,ss)
        
        # TODO Denuevo... pa'que?
        #Xnew =H3aux+matrix(rnorm((q)*n),ncol=q)
        
        #ynew=rnorm(n,0,ss)
        
        yeast$x=X
        yeast$y=y
      }

      j=1
      k=2
       
      a1=75/(1+25*floor((p-q)/2))
      a2=100/(1+25*ceiling((p-q)/2))
      beta.true=c(rep(a1,floor((p-q)/2)),rep(-a2,ceiling((p-q)/2)),rep(0,q))
        
      
      A1=plsr(yeast$y~yeast$x,ncomp=2,scale=FALSE)
      ####pls
      
      PLSS[hh,i,jj]=as.numeric(t(beta.true-A1$coefficients[,,2])%*%cov(yeast$x)%*%(beta.true-A1$coefficients[,,2]))
    }
    
  }
  
}


 

H=matrix(0,nrow=length(JJ),ncol=length(qq))
for (i in 1:nreps){
  H=H+PLSS[i,,]
}

H=H/nreps
## FIGURE 4.7
par(mar=c(5,6,4,1)+.1)
plot(JJ,H[,1],type='b' ,col='blue',ylim=c(0,25), lwd=2.5, xlim=c(100,1200), main="Figure 4.7", xlab=TeX("Number of predictors, $p$"), ylab=TeX("$||\\beta - \\hat{\\beta}||^2_{S_X}$"))
lines(JJ,H[,2],type='l' ,col='red', lwd=2.5)
 lines(JJ,H[,3],type='l' ,col='orange', lwd=2.5)
 lines(JJ,H[,4],type='l' ,col='dark green', lwd=2.5)
