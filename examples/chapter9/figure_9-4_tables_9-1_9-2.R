rm(list=ls())
library(fda)
library(fda.usc)
library(np)
data(tecator)
x=tecator$absorp.fdata$data
y=tecator$y
library(chemometrics)
names(y)

# PARA GENERAR TABLA
# table61 es un data.frame
print(knitr::kable(table61, "simple", caption="Table 6.1")) 

# Values of j correspond to
# j=1 FAT
# j=2 WATER
# j=3 PROTEIN
j=1


######choose d

#######
#######
#######
#######Table 9.1 line 5. and Table 9.2 line 5.
#######Inverse regression: NP-I-PLS with W TX
#######
#######
#######

set.seed(10)
#maximum d
dd=40


#####training
n=172
xx=x[1:n,]
yy=y[1:n,]


#definition of functions to be use (inverse regression)
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

##cross validation to choose d
set.seed(10)
K <- 5
folds <- cut(sample(seq(1, nrow(xx))), breaks=K, labels=FALSE)




######### First part: chose how many variables give the minimun error of prediction
######## For Line 5. from table 9.1

in_sample=matrix(0,ncol=dd,nrow=K) 
out_sample=matrix(0,ncol=dd,nrow=K)


for (kk in 1:K) {
  #cat(", i = ", i)
  #Randomly shuffle the data
  
  # Segement your data by fold using the which() function 
  testIndexes <- which(folds == kk, arr.ind=TRUE)
  
  
  
  
  
  
  
  for (jjj in 2:dd){
    model=  pls2_nipals(xx[-testIndexes,], yy[-testIndexes ,j],jjj,it = 2000,scale=FALSE)
    
    pred=yy[-testIndexes ,j]
    
   
    kk1_mult=xx[-testIndexes,]%*%model$P[,1:jjj]
    mod_inv=lm(kk1_mult~pred+I(pred^2)+I(pred^(1/2)))
    
    pre=predict(mod_inv)
    
    
    
    
    va=var(mod_inv$residuals)
    
    res=inv_g(kk1_mult,mod_inv$fitted.values,yy[-testIndexes ,j],va)
    
     
    #in sample
    in_sample[kk,jjj]=sqrt(mean((res-yy[-testIndexes ,j])^2))
    
    
    
    kk1_mult_out=xx[testIndexes,]%*%model$P[,1:jjj]
    res=inv_g(kk1_mult_out,mod_inv$fitted.values,yy[-testIndexes ,j],va)
    
    
    #out of sample
    out_sample[kk,jjj]=sqrt(mean((res-yy[testIndexes ,j])^2))
    
    
    
  }
 
}


 
in_sample[,1]= in_sample[,2]
out_sample[,1]= out_sample[,2]
in_error=sqrt(apply(in_sample^2,2,mean))
out_error=sqrt(apply(out_sample^2,2,mean))

##which one gives minimum error?
sss=which.min(out_error)





cat("Table 9.1 - Fat - Line 5")
cat("NP-I-PLS with W TX")
cat("dimension ", sss)  
 
 


#######Prediction with that error in testing with that dimension
testIndex=(n+1):dim(x)[1]
model=  pls2_nipals(x[-testIndexes,], y[-testIndexes ,j],sss,it = 2000,scale=FALSE)

pred=y[-testIndexes ,j]
kk1_mult=x[-testIndexes ,]%*%model$P[,1:sss]
mod_inv=lm(kk1_mult~pred+I(pred^2)+I(pred^(1/2)))

pre=predict(mod_inv)




va=var(mod_inv$residuals)
res=inv_g(kk1_mult,mod_inv$fitted.values,y[-testIndexes ,j],va)


#in sample
in_sample_p=sqrt(mean((res-y[-testIndexes ,j])^2))



kk1_mult_out=x[testIndexes,]%*%model$P[,1:sss]
res=inv_g(kk1_mult_out,mod_inv$fitted.values,y[-testIndexes ,j],va)


#out of sample
out_sample_p=sqrt(mean((res-y[testIndexes ,j])^2))


cat("Table 9.1 - Fat - Line 5","Table 9.2 - Fat - Line 5","Table 9.2 - Fat - Line 5")
cat("NP-I-PLS with W TX")
cat("dimennsion ", "In sample error", "Out of sample error")  
c(sss,in_sample_p,out_sample_p)

 







##################
##################
##################
##################
##################
##################
##################
##################
##################
##################names(y

###############


#######
######Table 9.1 and 9.2 line 1. to 4.
#######

inv_g_1<-function(z,prex,yy,va){
  w=NULL
 Rf=NULL
   for (h in 1:length(z)){
  for (i in 1:length(prex)){
    w[i]=exp(-1/2*(z[h]-prex[i])^2*va^(-1))
  }
  R=sum(w*yy)/sum(w)
  Rf[h]=R
  }
Rf
}




dd=40
######### First part: chose how many variables give the minimun error of prediction ########
in_sample_1=matrix(0,ncol=dd,nrow=K)
out_sample_1=matrix(0,ncol=dd,nrow=K)
pls_in=matrix(0,ncol=dd,nrow=K)
pls_out=matrix(0,ncol=dd,nrow=K)
non_param_out=matrix(0,ncol=dd,nrow=K)
non_param_in=matrix(0,ncol=dd,nrow=K)
cuad_in=matrix(0,ncol=dd,nrow=K)
cuad_out=matrix(0,ncol=dd,nrow=K)


K <- 5
set.seed(50)
folds <- cut(sample(seq(1, nrow(xx))), breaks=K, labels=FALSE)


for (kk in 1:K) {
  
  testIndexes <- which(folds == kk, arr.ind=TRUE)
  
  
  
  for (jjj in 1:dd){
    
    
    
    ###lineal pls
    model=  pls2_nipals(xx[-testIndexes,], yy[-testIndexes,j],jjj,it = 1000,scale=FALSE)
    
    
    predict_in=scale(xx[-testIndexes,],scale=FALSE)%*%model$B+mean(yy[-testIndexes,j])
    
    predict_out=scale(xx[testIndexes,],center=apply(xx[-testIndexes,],2,mean),scale=FALSE)%*%model$B+mean(yy[-testIndexes,j])
    
    
    ###mean square error for linear
    pls_in[kk,jjj]=sqrt(mean((yy[-testIndexes,j]-predict_in)^2))
    pls_out[kk,jjj]=sqrt(mean((yy[testIndexes,j]-predict_out)^2))
    
    
    
    
    
    
    kk1_uni=xx[-testIndexes,]%*%model$B
    
    
    kk1_uni_out=xx[testIndexes,]%*%model$B
    ########cuadratic
    
    
    modelo_2=lm(yy[-testIndexes,j]~kk1_uni+I(kk1_uni^2))
    
    
    
    
    
    betas=modelo_2$coefficients
    
    yfit=betas[1]+betas[2]*kk1_uni_out+betas[3]*kk1_uni_out^2
    
    
    h=cbind(yfit,yy[testIndexes,j],(yy[testIndexes,j]-yfit)^2)
    #out of sample error
    cuad_out[kk,jjj]=sqrt(mean(h[,3]))
    
    
    
    
    
    
     
    
    ######end cuadratci
    
    
    
    #######no parametric
    datos=data.frame(y=yy[-testIndexes,j], X=kk1_uni)
    
    
    model.np <- npreg(y ~ X,
                      regtype = "ll",
                      bwmethod = "cv.aic",
                      gradients = TRUE,
                      data = datos)
    
    non_param_in[kk,jjj]=sqrt(mean((fitted(model.np)-yy[-testIndexes,j])^2))
    
    
    datos.new=data.frame( X=kk1_uni_out)
    
    
    
    
    
    
    fit_new=predict(model.np,newdata=datos.new)
    
    #out non parametric
    non_param_out[kk,jjj]=sqrt(mean((fit_new-yy[testIndexes,j])^2))
    
    
    ####end non parametric
    
    ####3inverse regreession
    pred=yy[-testIndexes,j]
    mod_inv=lm(kk1_uni~pred+I(pred^2)+I(pred^(1/2)))
    
    pre=predict(mod_inv)
    
    
    
    
    va=var(mod_inv$residuals)
    
    
    res=inv_g_1(kk1_uni,mod_inv$fitted.values,yy[-testIndexes,j],va)
      
       
    
    #in sample
    in_sample_1[kk,jjj]=sqrt(mean((res-yy[-testIndexes,j])^2))
    
    res=inv_g_1(kk1_uni_out,mod_inv$fitted.values,yy[-testIndexes,j],va)
    
   
    
    #out of sample
    out_sample_1[kk,jjj]=sqrt(mean((res-yy[testIndexes,j])^2))
    
  }
}



inv_d=which.min(apply(out_sample_1^2,2,mean)[1:dd])
pls_d=which.min(apply(pls_out^2,2,mean)[2:dd])
cuad_d=which.min(apply(cuad_out^2,2,mean)[2:dd])
nonpar_d=which.min(apply(non_param_out^2,2,mean)[1:dd])




#######out of sample prediction
model=  pls2_nipals(x[1:n,], y[1:n,j],pls_d,it = 1000,scale=FALSE)

predict_in=scale(x[1:n,],scale=FALSE)%*%model$B+mean(y[1:n,j])

predict_out=scale(x[(n+1):dim(y)[1],],center=apply(x[1:n,],2,mean),scale=FALSE)%*%model$B+mean(y[1:n,j])


###mean square error for linear

in_sample_error=sqrt(mean((y[1:n,j]-predict_in)^2))
out_sample_error=sqrt(mean((y[(n+1):dim(y)[1],j]-predict_out)^2))

cat("Table 9.1 - Fat - Line 1","Table 9.2 - Fat - Line 1 - In sample","Table 9.2 - Fat - Line 1 - Out of sample")
cat("Pls lineal")  
c(pls_d,in_sample_error,out_sample_error)

 
 


#####quadratic PLS, get linear + quadratic fitting in the training
model=  pls2_nipals(x[1:n,], y[1:n,j],cuad_d,it = 1000,scale=FALSE)


kk1=scale(x[1:n,],center=rep(0,dim(x)[2]),scale=FALSE)%*%model$B 
kk1=x[1:n,]%*%model$B 


modelo=lm(y[1:n,j]~kk1+I(kk1^2))

 


kk1_test=x[(n+1):dim(x)[1],]%*%model$B
betas=modelo$coefficients
 
yfit=betas[1]+betas[2]*kk1_test+betas[3]*kk1_test^2

 
h=cbind(yfit,y[(n+1):dim(x)[1],j],(y[(n+1):dim(x)[1],j]-yfit)^2)
#out of sample error


cu_out=sqrt(mean(h[,3]))






yfit=betas[1]+betas[2]*kk1+betas[3]*kk1^2

h=cbind(yfit,y[1:n,j],(y[1:n,j]-yfit)^2)
#in sample error
cu_in=sqrt(mean(h[,3]))

cat("Table 9.1 - Fat - Line 2","Table 9.2 - Fat - Line 2 - In sample","Table 9.2 - Fat - Line 2 - Out of sample")
cat("Pls Quadratic")  
c(cuad_d,cu_in,cu_out)


 

######non parametric

model=  pls2_nipals(x[1:n,], y[1:n,j],nonpar_d,it = 1000,scale=FALSE)

kk1=x[1:n,]%*%model$B 


 



kk1_test=x[(n+1):dim(x)[1],]%*%model$B

datos=data.frame(y=y[1:n,j], X=kk1)

 
model.np <- npreg(y ~ X,
                   regtype = "ll",
                  bwmethod = "cv.aic",
                  gradients = TRUE,
                  data = datos)


 
no_in=sqrt(mean((fitted(model.np)-y[1:n,j])^2))


datos.new=data.frame( X=kk1_test)

 
 



fit_new=predict(model.np,newdata=datos.new)

#out non parametric
no_out= sqrt(mean((fit_new-y[(n+1):dim(x)[1],j])^2))



cat("Table 9.1 - Fat - Line 3","Table 9.2 - Fat - Line 3 - In sample","Table 9.2 - Fat - Line 3 - Out of sample")
cat("Non parametric PLS")  
c(nopar_d,no_in,no_out)
 


######
 
 
 # prediction with invere regression with beta
 
 model=  pls2_nipals(x[1:n,], y[1:n,j],inv_d,it = 1000,scale=FALSE)
 
 kk1=x[1:n,]%*%model$B 
 
 
 
 
 
 
 kk1_test=x[(n+1):dim(x)[1],]%*%model$B
 
 
 
 pred=y[1:n,j]
 #plot(kk1~pred)
 mod_inv=lm(kk1~pred+I(pred^2)+I(pred^(1/2)))
 inv_in=sqrt(mean(mod_inv$residuals^2))
 
 
 s=order(pred)
 pre=predict(mod_inv)
 
 

#lines(pred[s],pre[s])
va=var(mod_inv$residuals)

 res=inv_g_1(kk1,mod_inv$fitted.values,y[1:n,j],va)

 
#in sample
inv_in=sqrt(mean((res-y[1:n,j])^2))

res=inv_g_1(kk1_test,mod_inv$fitted.values,y[1:n,j],va)

 

#out of sample
inv_out=sqrt(mean((res-y[(n+1):dim(x)[1],j])^2))




cat("Table 9.1 - Fat - Line 4","Table 9.2 - Fat - Line 4 - In sample","Table 9.2 - Fat - Line 4 - Out of sample")
cat("Inverse PLS dimension 1")  
c(inv_d,inv_in, inv_out)

 




#####all results together

 c("d for linear pls","error in sample for linear pls","error out of sample for linear pls")
 c(pls_d,sqrt(mean((y[1:n,j]-predict_in)^2)),sqrt(mean((y[(n+1):dim(y)[1],j]-predict_out)^2)))
 
 
 c("d for quadratic pls","error in sample for quadratic pls","error out of sample for quadratic pls")
 c(cuad_d,cu_in,cu_out)

 c("d for NP pls","error in sample for NP pls","error out of sample for NP pls")
 c(nonpar_d,no_in, no_out)
 
 
 c("d for NP-I pls","error in sample for NP_I pls","error out of sample for NP_I pls")
 c(inv_d, inv_in, inv_out)
 
  
 ####Figure 9.4  (FAT)
 
 plot(4:25,sqrt(apply(pls_out^2,2,mean)[4:25]),ylim=c(1.5,4.5),type="l")
 lines(4:25,sqrt(apply(out_sample_1^2,2,mean)[4:25]))
 

 #############FINISH FAT
 
 