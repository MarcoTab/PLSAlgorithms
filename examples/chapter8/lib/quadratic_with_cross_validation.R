##quadratic and linear discimination with reduction

library(MASS)
library(nlme)
library(tidyverse) 
source("examples/chapter8/pls_matrices.R")

#  Xes is the data n times p, y is n by 1, K the number of folders
#  d de dimension and , y tiene que tener 1,2,3, solo permite eso
# estoy haciendo PLS 

quadratic.discrimination<-function(Xes,y,K,d,para){

 # library(torch)
 # set.seed(10)
  cv.prediction <- c()
  
  #making the folds

folds <- cut(sample(seq(1, nrow(Xes))), breaks=K, labels=FALSE)


p.pls_1=NULL
p.pls_2=NULL

  
for (kk in 1:K) {
    #cat(", k = ", kk)
    #Randomly shuffle the data
    
    # Segement your data by fold using the which() function 
    testIndexes <- which(folds == kk, arr.ind=TRUE)
    
    
  #divido los X en traing y test  
  holdout <- as.matrix(Xes[-testIndexes,])
  X.test <- matrix(Xes[testIndexes,],ncol=ncol(Xes))
  
  
 
  
  d#training y
  y.holdhout<-y[-testIndexes]
  
  
  
  S_X <- cov(holdout)
  uY <- unique(y.holdhout) %>% set_names()
  S_k <- lapply(uY, function(k) cov(holdout[y.holdhout == k, , drop = FALSE])) 
  cov.mean <- lapply(uY, function(k) S_k[[k]] * mean(y.holdhout==k)) %>%
    Reduce('+', .)
  cov.ss1 <- lapply(uY, function(k) (S_X - S_k[[k]]) %*% t(S_X - S_k[[k]]) * mean(y.holdhout==k)) %>%
    Reduce('+', .)
  cov.ss2 <- lapply(uY, function(k) (cov.mean - S_k[[k]]) %*% t((cov.mean - S_k[[k]])) * mean(y.holdhout==k)) %>%
    Reduce('+', .)
  
  
  
  
   
  
  
  dif_d_delta2=(S_X-cov.mean)%*%t(S_X-cov.mean)
  
  
  
  A_al_para=para*dif_d_delta2+(1-para)*cov.ss2
  
  A_al1=cov.ss1
  
  AUX_para=pls_matrices(A_al_para, cov.mean, d)
  AUX=pls_matrices(A_al1, cov.mean, d)
   
   
    
  
  datos_proy_2=holdout%*%AUX_para$Gamma
  datos_proy_test_2=X.test%*%AUX_para$Gamma
  
  datos_proy_1=holdout%*%AUX$Gamma
  datos_proy_test_1=X.test%*%AUX$Gamma
  
  
  
  training.pls_1 <- data.frame(X=datos_proy_1,Y=y[-testIndexes])
  testing.pls_1 <- data.frame(X=datos_proy_test_1,Y=y[testIndexes])
  
  
  
  
  quad.pls_1 <- qda(Y~., training.pls_1)
  
  
  p.pls_1[testIndexes] <- predict(quad.pls_1, testing.pls_1)$class
  
  
  training.pls_2 <- data.frame(X=datos_proy_2,Y=y[-testIndexes])
  testing.pls_2 <- data.frame(X=datos_proy_test_2,Y=y[testIndexes])
  
  
  
  
  quad.pls_2 <- qda(Y~., training.pls_2)
  
  
  p.pls_2[testIndexes] <- predict(quad.pls_2, testing.pls_2)$class
  
  
}

list(predict_quad_1=p.pls_1,predict_quad_2=p.pls_2)}

