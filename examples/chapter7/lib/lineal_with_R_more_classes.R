##quadratic and linear discimination with reduction

library(MASS)
library(nlme)
library(dr)
library(klaR)
library(chemometrics)
library(psych)
library(MASS)
library(mltools)
library(data.table)
library(devtools)

#  Xes is the data n times p, y is n by 1, K the number of folders
#  d de dimension and , y tiene que tener 0 y 1, solo permite eso
# estoy haciendo PLS 

lineal.junto.small.R.more.classes<-function(Xes,y,K,d){
  
  p.pls=y
  p.iso=y
  p.PFC=y
  #d should be small than the number of samples
  
  # library(torch)
  cv.prediction <- c()
  
  #making the folds
  
  folds <- cut(sample(seq(1, nrow(Xes))), breaks=K, labels=FALSE)
  
  

  
  for (kk in 1:K) {
    # Segement your data by fold using the which() function 
    testIndexes <- which(folds == kk, arr.ind=TRUE)
    
    auxxx=one_hot(as.data.table(as.factor(y[-testIndexes])))
    AUX=pls2_nipals(Xes[-testIndexes,],auxxx,d)$P[,1:d]
    
    datos_proy=Xes[-testIndexes,]%*%AUX
    datos_proy_test=Xes[testIndexes,]%*%AUX
    
    # PLS
    training.pls <- data.frame(X=datos_proy,Y=y[-testIndexes])
    testing.pls <- data.frame(X=datos_proy_test,Y=y[testIndexes])
    
    linear.pls <- lda(Y~., training.pls)
    
    p.pls[testIndexes] <- predict(linear.pls, testing.pls)$class
     
    
    ######iso
    reg=lm(Xes[-testIndexes,]~as.matrix(auxxx))
    vec=eigen(cov(reg$residuals))$vectors[,1:d]
    training.iso=data.frame(X=Xes[-testIndexes,]%*%vec,Y=y[-testIndexes])
    testing.iso=data.frame(X=Xes[testIndexes,]%*%vec,Y=y[testIndexes])
    
    
    linear.iso <- lda(Y~., training.iso)
    
    
    p.iso[testIndexes] <- predict(linear.iso, testing.iso)$class
    
    
    ######PFC
    
    
    if (d==1| d==2){
      training.PFC=data.frame(X=datos_proy,Y=y[-testIndexes])
      testing.PFC=data.frame(X=datos_proy_test,Y=y[testIndexes])}
    else{s0 <- dr(y[-testIndexes]~datos_proy)
    training.PFC=data.frame(X=datos_proy%*%s0$evectors[,1:2],Y=y[-testIndexes])
    testing.PFC=data.frame(X=datos_proy_test%*%s0$evectors[,1:2],Y=y[testIndexes])
    }
 
    
      linear.PFC <- lda(Y~., training.PFC)
      
      
      p.PFC[testIndexes] <- predict(linear.PFC, testing.PFC)$class
     
  
    }  
  
  list(predict_lineal_pls=p.pls,predict_lineal_PFC=p.PFC, predict_lineal_iso=p.iso)}
