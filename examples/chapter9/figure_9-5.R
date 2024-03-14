rm(list=ls())
library(readr)
solvents <- read_csv("data/chapter9/solvents.csv")




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





inv_g_1<-function(z,prex,yy,va){
  w=NULL
  Rf=NULL
  for (h in 1:length(z)){
    for (i in 1:length(prex)){
      w[i]=exp(-1/2*t(z[h]-prex[i])^2/va)
    }
    R=sum(w*yy)/sum(w)
    Rf[h]=R
  }
  Rf
}



dim(solvents)
X=solvents[,c(2:3,5:10)]
Y=solvents[,4]
YY=Y[Y<61]
XX=X[Y[,1]<61,]
err_n=NULL
error_p=NULL
XX=scale(XX,center=TRUE,scale=TRUE)




dd=8

K1=10

out_sample_p=matrix(0,ncol=dd,nrow=K1)

errro_pl=matrix(0,ncol=dd,nrow=K1)

 
set.seed(1)
folds1 <- cut(sample(seq(1, nrow(XX))), breaks=K1, labels=FALSE)


for (m in 1:8){
for (hh in 1:K1) {
  #cat(", i = ", i)
  #Randomly shuffle the data
  
  # Segement your data by fold using the which() function 
  testIndexes_a <- which(folds1 == hh, arr.ind=TRUE)
  
  x=XX[-testIndexes_a,]
  y=YY[-testIndexes_a]
  
  
  
  
  
   
  
  
  
   
    
    
    
    
   
  model=  pls2_nipals(x, y,m ,it = 2000,scale=FALSE)
  
  
  Xval=XX[testIndexes_a,]
  yval=YY[testIndexes_a]
  
  predict_2=scale(Xval,center=apply(x,2,mean),scale=FALSE)%*%model$B+mean(y)
  
  
  
  pred=y
  kk1_mult=x%*%model$P[,1:m]
  mod_inv=lm(kk1_mult~pred+I(pred^2)+log(pred)+I(pred^(1/3))+I(pred^(1/2))+I(pred^3))
  
  
  pre=predict(mod_inv)
  
  
  
  
  va=var(mod_inv$residuals)
  
  if(m==1){res=inv_g_1(kk1_mult,mod_inv$fitted.values,y,va)}else
  {res=inv_g(kk1_mult,mod_inv$fitted.values,y,va)}
  
  
  #in sample
  in_sample_p =sqrt(mean((res-y)^2))
  
  Xtest=Xval
  y1test=yval
  
  kk1_mult_out=Xtest%*%model$P[,1:m ]
  if(m==1){res=inv_g_1(kk1_mult_out,mod_inv$fitted.values,y,va)}else{
    res=inv_g(kk1_mult_out,mod_inv$fitted.values,y,va)}
  
  
  
  yhat =res
  
  #out of sample inverse regression 
  out_sample_p[hh,m] =sqrt(mean((res-y1test)^2))
  
  
  #out of sample linear pls
  errro_pl[hh,m]=sqrt(mean((predict_2-y1test)^2))
  
  
  
  
  
  
}

}

 


plot(1:8,apply(errro_pl,2,mean),ylim=c(4,7),type='l',ylab="PRMSE",xlab="Components")
lines(1:8,apply(out_sample_p,2,mean),col='red')
 
 
lines(1:8,c(4.55,4.51,4.6,4.28,4.45,4.32,4.38,4.26),col='blue')

title("Figure 9.5")














