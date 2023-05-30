##quadratic and linear discimination with reduction
 
library(MASS)
#library(nlme)
library(dr)
#library(klaR)
#library(psych)
#library(MASS)
 
library(devtools)
library(chemometrics)

#  Xes is the data n times p, y is n by 1, K the number of folds
lineal.junto.small.R <- function(Xes,y,K,d){
  
  p.pls=y
  p.iso=y
  p.PFC=y
  # d should be smaller than the number of samples
  
  # library(torch)
 
  cv.prediction <- c()
  
  #making the folds
  
  folds <- cut(sample(seq(1, nrow(Xes))), breaks=K, labels=FALSE)
  
  for (kk in 1:K) {
    #cat(", k = ", kk)
    #Randomly shuffle the data
    
    
    testIndexes <- which(folds == kk, arr.ind=TRUE)
    
    
    #divido los X en traing y test  
     
    XX=Xes[-testIndexes,]
    yy=y[-testIndexes]
    
    ## no se porque use esto
    AUX=chemometrics::pls2_nipals(XX,yy,d)$P[,1:d]
    
    datos_proy=Xes[-testIndexes,]%*%AUX
    datos_proy_test=Xes[testIndexes,]%*%AUX
    
    training.pls <- data.frame(X=datos_proy,Y=y[-testIndexes])
    testing.pls <- data.frame(X=datos_proy_test,Y=y[testIndexes])
    
    linear.pls <- lda(Y~., training.pls)
    
    
    p.pls[testIndexes] <- predict(linear.pls, testing.pls)$class
     
    
    ######iso
    reg=lm(Xes[-testIndexes,]~y[-testIndexes])
    
    vec=eigen(cov(reg$residuals))$vectors[,1:d]
    training.iso=data.frame(X=Xes[-testIndexes,]%*%vec,Y=y[-testIndexes])
    testing.iso=data.frame(X=Xes[testIndexes,]%*%vec,Y=y[testIndexes])
    
    
    linear.iso <- lda(Y~., training.iso)
    
    
    p.iso[testIndexes] <- predict(linear.iso, testing.iso)$class
    
    
    ######PFC
    if (d==1){
      training.PFC=data.frame(X=datos_proy,Y=y[-testIndexes])
      testing.PFC=data.frame(X=datos_proy_test,Y=y[testIndexes])}
    else{s0 <- dr(y[-testIndexes]~datos_proy)
    training.PFC=data.frame(X=datos_proy%*%s0$evectors[,1],Y=y[-testIndexes])
    testing.PFC=data.frame(X=datos_proy_test%*%s0$evectors[,1],Y=y[testIndexes])
    }
 
    
      linear.PFC <- lda(Y~., training.PFC)
      
      p.PFC[testIndexes] <- predict(linear.PFC, testing.PFC)$class
     
  
    }  
  
  list(predict_lineal_pls=p.pls,predict_lineal_PFC=p.PFC, predict_lineal_iso=p.iso)}
