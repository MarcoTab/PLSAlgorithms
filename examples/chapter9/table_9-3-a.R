rm(list=ls())
library(fda)
library(fda.usc)
library(np)
library(chemometrics)
 
Xcal <- read.table("data/chapter9/Xcal.txt", quote="\"", comment.char="")
Xcal=t(Xcal)

Xval <- read.table("data/chapter9/Xval.txt", quote="\"", comment.char="")
Xval=t(Xval)

Ycal <- read.table("data/chapter9/Ycal.txt", quote="\"", comment.char="")

Yval <- read.table("data/chapter9/Yval.txt", quote="\"", comment.char="")
Xcal=Xcal[,3:9625]
Xval=Xval[,3:9625]

x=Xcal
y=Ycal


## first response
j=1
 

 
n=dim(Xcal)[1]
 


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
K <- 26
folds <- cut(sample(seq(1, nrow(x))), breaks=K, labels=FALSE)




######### how many dimension using W^T X for inverse and non parametri
######### non parametric is not here, for inversse is Table 9.3-5


in_sample=matrix(0,ncol=dd,nrow=K) 
out_sample=matrix(0,ncol=dd,nrow=K)
non_param_out_m=matrix(0,ncol=dd,nrow=K) 
non_param_in_m=matrix(0,ncol=dd,nrow=K) 

for (kk in 1:K) {
  #cat(", i = ", i)
  #Randomly shuffle the data
  
  # Segement your data by fold using the which() function 
  testIndexes <- which(folds == kk, arr.ind=TRUE)
  
  
  
  
  
  
  
  for (jjj in 2:dd){
    model=  pls2_nipals(x[-testIndexes,], y[-testIndexes ,j],jjj,it = 2000,scale=FALSE)
    
    pred=y[-testIndexes ,j]
    
    
    kk1_mult=x[-testIndexes,]%*%model$P[,1:jjj]
    mod_inv=lm(kk1_mult~pred+I(pred^2)+I(pred^(1/2))+I(pred^3))
    
    pre=predict(mod_inv)
    
    
    
    
    va=var(mod_inv$residuals)
    
    res=inv_g(kk1_mult,mod_inv$fitted.values,y[-testIndexes ,j],va)
    
    
    #in sample
    in_sample[kk,jjj]=sqrt(mean((res-y[-testIndexes ,j])^2))
    
    
    
    kk1_mult_out=x[testIndexes,]%*%model$P[,1:jjj]
    res=inv_g(kk1_mult_out,mod_inv$fitted.values,y[-testIndexes ,j],va)
    
    
    #out of sample
    out_sample[kk,jjj]=sqrt(mean((res-y[testIndexes ,j])^2))
    
    
    
    
    
    
    X<-data.frame(kk1_mult)
    y1<-c(y[-testIndexes,j])
    # npregbw(formula=y1~X, regtype="ll", bwmethod="cv.aic")
    
    
    # banda<-npregbw(xdat=X ,ydat=y1)
    
    bw <- npregbw(xdat=X, ydat=y1, regtype="ll", bwmethod="cv.aic")
    fit <- fitted(npreg(exdat=X , bws=bw))
    
    
    
    non_param_in_m[kk,jjj]=sqrt(mean((fit-y[-testIndexes,j])^2))
    
    
    X_new<-data.frame(kk1_mult_out)
    
    fit_new <- fitted(npreg(exdat=X_new , bws=bw))
    
    
    
    
    #fit_new=predict(model.np,newdata=X_new)
    
    #out non parametric
    non_param_out_m[kk,jjj]=sqrt(mean((fit_new-y[testIndexes,j])^2))
    
    
    
    
    
    
  }
  
}


non_param_out_m[,1]=non_param_out_m[,2]
in_sample[,1]= in_sample[,2]
out_sample[,1]= out_sample[,2]
in_error=sqrt(apply(in_sample^2,2,mean))
out_error=sqrt(apply(out_sample^2,2,mean))

n_p_out=sqrt(apply(non_param_out_m^2,2,mean))
ss=which.min(n_p_out)
sss=which.min(out_error)
sss
ss

 

 
###########
#######predecir con ese d (y luego vr en out of sample) Table 9.3. 5)
##########
 
model=  pls2_nipals(x, y[,1],sss,it = 2000,scale=FALSE)

pred=y[,1]
kk1_mult=x%*%model$P[,1:sss]

 
mod_inv=lm(kk1_mult~pred+I(pred^2)+I(abs(pred)^(1/2))+I(pred^3))

pre=predict(mod_inv)




va=var(mod_inv$residuals)
res=inv_g(kk1_mult,mod_inv$fitted.values,y[,1],va)


#in sample
in_sample_p=sqrt(mean((res-y[,1])^2))



kk1_mult_out=Xval%*%model$P[,1:sss]
res=inv_g(kk1_mult_out,mod_inv$fitted.values,y[,1],va)









#out of sample
out_sample_p=sqrt(mean((res-Yval[,1])^2))

c("Choosen d","Inverse mult out")
c(sss,out_sample_p)



########
##non-parametric
model=  pls2_nipals(x, y[,1],ss,it = 2000,scale=FALSE)

pred=y[,1]
kk1_mult=x%*%model$P[,1:ss]


X<-data.frame(kk1_mult)
y1<-c(y[,j])

bw <- npregbw(xdat=X, ydat=y1, regtype="ll", bwmethod="cv.aic")
fit <- fitted(npreg(exdat=X , bws=bw))

in_sample_nonp=sqrt(mean((fit-y[,1])^2))

kk1_mult_out=Xval%*%model$P[,1:ss]

X_new<-data.frame(kk1_mult_out)

fit_new <- fitted(npreg(exdat=X_new , bws=bw))

out_sample_nonp=sqrt(mean((fit_new-Yval[,1])^2))

c("choosen d","Inverse mult out")
c(ss,out_sample_nonp)





 
 









########


##################
##################Table 9.3 lines 1,2, 3,4

###############

 

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

set.seed(8)
K <- 26
folds <- cut(sample(seq(1, nrow(x))), breaks=K, labels=FALSE)



dd=15
######### First part: chose how many variables give the minimun error of prediction ########
in_sample_1=matrix(0,ncol=dd,nrow=K)
out_sample_1=matrix(0,ncol=dd,nrow=K)
pls_in=matrix(0,ncol=dd,nrow=K)
pls_out=matrix(0,ncol=dd,nrow=K)
non_param_out=matrix(0,ncol=dd,nrow=K)
non_param_in=matrix(0,ncol=dd,nrow=K)
cuad_in=matrix(0,ncol=dd,nrow=K)
cuad_out=matrix(0,ncol=dd,nrow=K)
for (kk in 1:K) {
  
  testIndexes <- which(folds == kk, arr.ind=TRUE)
  
  
  
  for (jjj in 1:dd){
    
    
    
    ###lineal pls
    model=  pls2_nipals(x[-testIndexes,], y[-testIndexes,j],jjj,it = 1000,scale=FALSE)
    
    
    predict_in=scale(x[-testIndexes,],scale=FALSE)%*%model$B+mean(y[-testIndexes,j])
    
    predict_out=(x[testIndexes,]-apply(x[-testIndexes,],2,mean))%*%model$B+mean(y[-testIndexes,j])
    
    
    ###mean square error for linear
    pls_in[kk,jjj]=sqrt(mean((y[-testIndexes,j]-predict_in)^2))
    pls_out[kk,jjj]=sqrt(mean((y[testIndexes,j]-predict_out)^2))
    
    
    
    
    
    
    
    
    kk1_uni=x[-testIndexes,]%*%model$B
    
    
    kk1_uni_out=x[testIndexes,]%*%model$B
    ########cuadratic
    
    
    modelo_2=lm(y[-testIndexes,j]~kk1_uni+I(kk1_uni^2))
    
    
    
    
    
    betas=modelo_2$coefficients
    
    yfit=betas[1]+betas[2]*kk1_uni_out+betas[3]*kk1_uni_out^2
    
    
   # h=cbind(yfit,y[testIndexes,j],(y[testIndexes,j]-yfit)^2)
    #out of sample error
    cuad_out[kk,jjj]=sqrt(mean((yfit-y[testIndexes,1])^2))
    
    
    
    
    
    
    yfit=betas[1]+betas[2]*kk1_uni+betas[3]*kk1_uni^2
    
    #h=cbind(yfit,y[1:n,j],(y[1:n,j]-yfit)^2)
    #in sample error
    cuad_in[kk,jjj]=sqrt(mean((y[-testIndexes,1]-yfit)^2))
    
    
    ######end cuadratci
    
    
    
    #######no parametric
    datos=data.frame(y=y[-testIndexes,j], X=kk1_uni)
    
    
    model.np <- npreg(y ~ X,
                      regtype = "ll",
                      bwmethod = "cv.aic",
                      gradients = TRUE,
                      data = datos)
    
    non_param_in[kk,jjj]=sqrt(mean((fitted(model.np)-y[-testIndexes,j])^2))
    
    
    datos.new=data.frame( X=kk1_uni_out)
    
    
    
    
    
    
    fit_new=predict(model.np,newdata=datos.new)
    
    #out non parametric
    non_param_out[kk,jjj]=sqrt(mean((fit_new-y[testIndexes,j])^2))
    
    
    ####end non parametric
    
    ####3inverse regreession
    pred=y[-testIndexes,j]
    mod_inv=lm(kk1_uni~pred+I(pred^2)+I(pred^(1/2)+I(pred^3)))
    
    pre=predict(mod_inv)
    
    
    
    
    va=var(mod_inv$residuals)
    
    
    res=inv_g_1(kk1_uni,mod_inv$fitted.values,y[-testIndexes,j],va)
    
    
    
    #in sample
    in_sample_1[kk,jjj]=sqrt(mean((res-y[-testIndexes,j])^2))
    
    res=inv_g_1(kk1_uni_out,mod_inv$fitted.values,y[-testIndexes,j],va)
    
    
    
    #out of sample
    out_sample_1[kk,jjj]=sqrt(mean((res-y[testIndexes,j])^2))
    
  }
}




pls_d=which.min(apply(pls_out^2,2,mean)[2:dd])
cuad_d=which.min(apply(cuad_out^2,2,mean)[2:dd])
nonpar_d=which.min(apply(non_param_out^2,2,mean)[1:dd])
inv_d=which.min(apply(out_sample_1^2,2,mean)[1:dd])


 
 



 
#######predecir con todos los metodos
model=  pls2_nipals(x , y[ ,1],pls_d,it = 1000,scale=FALSE)

predict_in=scale(x,scale=FALSE)%*%model$B+mean(y[,j])

predict_out=scale(Xval,center=apply(x,2,mean),scale=FALSE)%*%model$B+mean(y[,j])


###mean square error for linear
c("d for Pls", "Pls in", "Pls out")
c(pls_d,sqrt(mean((y[,j]-predict_in)^2)),sqrt(mean((Yval[,1]-predict_out)^2)))



#####quadratic PLS, get linear + quadratic fitting in the training
model=  pls2_nipals(x, y[,1],cuad_d,it = 1000,scale=FALSE)


kk1=scale(x,center=rep(0,dim(x)[2]),scale=FALSE)%*%model$B 
kk1=x%*%model$B 


modelo=lm(y[,1]~kk1+I(kk1^2))




kk1_test=Xval%*%model$B
betas=modelo$coefficients

yfit=betas[1]+betas[2]*kk1_test+betas[3]*kk1_test^2


h=cbind(yfit,Yval[,1],(Yval[,1]-yfit)^2)
#out of sample error


cu_out=sqrt(mean(h[,3]))






yfit=betas[1]+betas[2]*kk1+betas[3]*kk1^2

h=cbind(yfit,y[ ,j],(y[ ,j]-yfit)^2)
#in sample error
cu_in=sqrt(mean(h[,3]))
c("dimension quad pls","Quad in", "Quad out")
c(cuad_d,cu_in,cu_out)

######non parametric

model=  pls2_nipals(x[ ,], y[ ,j],nonpar_d,it = 1000,scale=FALSE)

kk1=x[ ,]%*%model$B 






kk1_test=Xval%*%model$B

datos=data.frame(y=y[ ,j], X=kk1)


model.np <- npreg(y ~ X,
                  regtype = "ll",
                  bwmethod = "cv.aic",
                  gradients = TRUE,
                  data = datos)


#model.np <- npreg(y ~ X,
#                  ckertype="epanechnikov", ckerorder=2,data= datos)
#in non parametric
no_in=sqrt(mean((fitted(model.np)-y[,j])^2))


datos.new=data.frame( X=kk1_test)






fit_new=predict(model.np,newdata=datos.new)

#out non parametric
no_out= sqrt(mean((fit_new-Yval[,1])^2))

c("dimension nonpar","No parametric in","No parametric out")
c(nonpar_d,no_in, no_out)



######


# prediction with invere regression with beta

model=  pls2_nipals(x[ ,], y[ ,j],inv_d,it = 1000,scale=FALSE)

kk1=x[,]%*%model$B 






kk1_test=Xval%*%model$B



pred=y[ ,j]
#plot(kk1~pred)
mod_inv=lm(kk1~pred+I(pred^2)+I(pred^(1/2)))
inv_in=sqrt(mean(mod_inv$residuals^2))


s=order(pred)
pre=predict(mod_inv)



#lines(pred[s],pre[s])
va=var(mod_inv$residuals)

res=inv_g_1(kk1,mod_inv$fitted.values,y[ ,j],va)


#in sample
inv_in=sqrt(mean((res-y[ ,j])^2))

res=inv_g_1(kk1_test,mod_inv$fitted.values,y[ ,j],va)



#out of sample
inv_out=sqrt(mean((res-Yval[,1])^2))


c("Dimension inv in beta","Inv 1 dimension in","Inv 1 dimension out")
c(inv_d,inv_in, inv_out)


 