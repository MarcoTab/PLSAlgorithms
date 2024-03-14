rm(list=ls())
library(fda)
library(fda.usc)
library(np)

 

X <- read.table("data/chapter9/X.txt", quote="\"", comment.char="")
Y <- read.table("data/chapter9/Y.txt", quote="\"", comment.char="")

j=1

 
library(chemometrics)





######choose d
##########
##########


#######
######choosing q, mutilivariate KOFFI AND DENNIS
#######
#######
set.seed(10)
dd=15



inv_g<-function(z,prex,yy,va){
  w=NULL
  Rf=NULL
  for (h in 1:dim(z)[1]){
    for (i in 1:dim(prex)[1]){
      w[i]=exp(-1/2*t(z[h,]-prex[i,])%*%solve(va)%*%(z[h,]-prex[i,]))
    }
    R=sum(w*yy)/sum(w)
    Rf[h]=R
  }
  Rf
}


inv_g_for_1<-function(z,prex,yy,va){
  w=NULL
  Rf=NULL
  for (i in 1:dim(prex)[1]){
    w[i]=exp(-1/2*(z-prex[i,])%*%solve(va)%*%t(z-prex[i,]))
  }
  R=sum(w*yy)/sum(w)
}



set.seed(10) 
 

K=nrow(X)
######### First part: chose how many proyections give the minimun error of prediction ########
in_sample=matrix(0,ncol=dd,nrow=K) 
out_sample=matrix(0,ncol=dd,nrow=K)
 
for (kk in 1:nrow(X)) {
  #cat(", i = ", i)
  #Randomly shuffle the data
  
  # Segement your data by fold using the which() function 
  testIndexes <- kk
  
  
  
  
  
  
  
  for (jjj in 2:dd){
    model=  pls2_nipals(X[-testIndexes,], Y[-testIndexes ,j],jjj,it = 2000,scale=FALSE)
    
    pred=Y[-testIndexes ,j]
    
    
    kk1_mult=t(t(X[-testIndexes,]))%*%model$P[,1:jjj]
    mod_inv=lm(kk1_mult~pred+I(pred^2)+I(pred^(1/2))+I(pred^3))
    
    pre=predict(mod_inv)
    
    
    
    
    va=var(mod_inv$residuals)
    
    res=inv_g(kk1_mult,mod_inv$fitted.values,Y[-testIndexes ,j],va)
    
    
    #in sample
    in_sample[kk,jjj]=sqrt(mean((res-Y[-testIndexes ,j])^2))
    
    
    
    kk1_mult_out=t(t(X[testIndexes,]))%*%model$P[,1:jjj]
    res=inv_g(kk1_mult_out,mod_inv$fitted.values,Y[-testIndexes ,j],va)
    
    
    #out of sample
    out_sample[kk,jjj]=sqrt(mean((res-Y[testIndexes ,j])^2))
    
    
    
     
     
    
    
    
    
  }
  
}



in_sample[,1]= in_sample[,2]
out_sample[,1]= out_sample[,2]

in_error=sqrt(apply(in_sample^2,2,mean))
out_error=sqrt(apply(out_sample^2,2,mean))

sss=which.min(out_error)


######## table 9.3 line (c)

c("dimnesion","CVRPE")
c(sss,mean(out_sample[,sss]))

  