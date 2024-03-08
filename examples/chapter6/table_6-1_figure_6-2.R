rm(list=ls())
library(pls)
library(glmnet)
library(Rcpp)
library(Renvlp)
library(ggplot2)

library(readr)
library(foreign)

#Data Base without Exports, Imports and Government Consumption (because with these variables appear a lot of missing)
aux<-read.dta("data/chapter6/DepuredNewDataGovernance_withLagGDP_22.dta")


### Table 6.1 Line 1: n=161, 2003-2018
aux_2003<-aux[which(aux$year>=2003),]

### Table 6.1 Line 2: n=109, 2008-2018
aux_2008<-aux[which(aux$year>=2008),]

### Table 6.1 Line 3: n=86, 2010-2018
aux_2010<-aux[which(aux$year>=2010),]

### Table 6.1 Line 4: n=53, 2010-2014
aux_2010_2014<-aux[which(aux$year>=2010 & aux$year <= 2014),]

### Table 6.1 Line 5: n=51, 2013-2018
aux_2013<-aux[which(aux$year>=2013),]

### Table 6.1 Line 6: n=41, 2014-2018
aux_2014<-aux[which(aux$year>=2014),]

### Table 6.1 Line 7: n=32, 2015-2018
aux_2015<-aux[which(aux$year>=2015),]


table_lines = list(aux_2003, aux_2008, aux_2010, aux_2010_2014, aux_2013, aux_2014, aux_2015)

table61 <- data.frame(
  n = c(161, 109, 86, 53, 51, 41, 32),
  Period = c("2003-2018", "2008-2018", "2010-2018", "2010-2014", "2013-2018", "2014-2018", "2015-2018"),
  Countries = c(12, 12, 12, 12, 10, 9, 9),
  PPLS.q1 = c(rep("-", 7)),
  PPCA.q1 = c(rep("-", 7)),
  PPLS.1 = c(rep("-", 7)),
  PPCA.1 = c(rep("-", 7)),
  LASSO = c(rep("-", 7)),
  FULL = c(rep("-", 7)),
  PPENV.q1 = c(rep("-", 7)),
  PPENV.1 = c(rep("-", 7))
)


index = 1

cat("Generating Table 6.1\n")

for (this_aux in table_lines) {
  cat("Row", index, "\n")
  this_aux=this_aux[-2,]
  n=nrow(this_aux)
  
  # Using 24 indicators (to reduce)
  X1=this_aux[,c(4:9,11:13,15:17,19:21,23:25,27:29,31:33)]
  
  X2=this_aux[,c(34:39,41,43:53,55)] # with FDI, Domestic Credit, Age Dependency, Popul, Seconday School, Inflation and control countries 
  
  y=as.numeric(unlist(this_aux[,1]))
  
  mm=ncol(X1)-2
  
  #To Find optimal d for PLS
  
  pls_predic=matrix(0,ncol=mm,nrow=n)
  MSE_PLS=matrix(0,ncol=mm,nrow=n)
  
  ppls_predic=matrix(0,ncol=mm,nrow=n)
  MSE_PPLS=matrix(0,ncol=mm,nrow=n)
  
  pca_predic=matrix(0,ncol=mm,nrow=n)
  MSE_PCA=matrix(0,ncol=mm,nrow=n)
  
  ols_predic=matrix(0,ncol=mm,nrow=n)
  MSE_OLS=matrix(0,ncol=mm,nrow=n)
  
  lasso_predic=matrix(0,ncol=mm,nrow=n)
  MSE_LASSO=matrix(0,ncol=mm,nrow=n)
  
  env_predic=matrix(0,ncol=mm,nrow=n)
  MSE_EPLS=matrix(0,ncol=mm,nrow=n)
  
  for (d in 1:mm){
    for (i in 1:n){
      cat((d-1)*n+i, ", ", sep="")
      #training
      X11.train=X1[-i,]
      X22.train=X2[-i,]
      y.train=y[-i]
      
      #testing
      X11.test=X1[i,]
      X22.test=X2[i,]
      y.test=y[i]
      
      ### Matrix X to reduce
      X11.train=as.matrix(X11.train)
      
      M=data.frame(X22.train)
      
      ### Regression of the continuous X in which we do not reduce
      A=lm(X11.train~.,data=M)
      
      # Save residuals
      res.train=data.frame(RES=A$residuals)
      
      M1=data.frame(X22.test)
      
      ### Predict on new data
      A_1=predict(A,newdata=M1)
      
      # Center training response variable
      y.train.center=y.train-mean(y.train)
      
      ## Perform pls on residuals in training
      m=plsr(y.train.center~.,ncomp=d,data=res.train)
      
      #Perform regression on non-reduced variables
      j=lm(y.train.center~.,data=M)
      
      #Predict non-reduced on new data
      j_2=predict(j,newdata=M1)
      
      #Final prediction on new data
      yhat=mean(y.train)+j_2+as.numeric(as.numeric(X11.test-A_1,nrow=1)%*%(as.numeric(m$coefficients[,1,d],ncol=1)))
      
      pls_predic[i,d]=yhat
      MSE_PLS[i,d]=yhat-y.test

      ##########LASSO
      
      covariates=as.matrix(cbind(X11.train, X22.train))
      z=y.train
      cv_model <- cv.glmnet(covariates,z, alpha = 1)
      
      # Find optimal lambda value that minimizes test MSE
      best_lambda <- cv_model$lambda.min
      best_model <- glmnet(covariates, z, alpha = 1, lambda = best_lambda)
      # best_model$df
      covariates.test=as.matrix(cbind(X11.test,X22.test))
      # Use lasso regression model to predict response value
      lasso_predic[i]=predict(best_model, s = best_lambda, newx = covariates.test)
      
      yhat=lasso_predic[i]
      
      MSE_LASSO[i,d]=(yhat-y.test)
      
      
      ######## PARTIAL PCA 

      ## Perform PCR on training residuals
      m=pcr(y.train.center~.,ncomp=d,data=res.train)
      
      # Perform regression on non-reduced variables
      j=lm(y.train.center~.,data=M)
      
      #Predict non-reduced on new data
      j_2=predict(j,newdata=M1)
      
      #Final prediction on new data
      yhat=mean(y.train)+j_2+as.numeric(as.numeric(X11.test-A_1,nrow=1)%*%(as.numeric(m$coefficients[,1,d],ncol=1)))
      pca_predic[i,d]=yhat
      MSE_PCA[i,d]=yhat-y.test
      
      
      #####full
      
      if (nrow(X11.train)>ncol(X1)+ncol(X2)){
        train.predictors=data.frame(X11.train,X22.train)
        predictors=as.matrix(train.predictors)
        y.train=as.matrix(y.train)
        data.train= data.frame(xx=predictors,y=y.train)
     
        fit=lm(y~.,data=data.train)
        
        #Predicting
        scores = as.matrix(X11.test)
        test.predictors=data.frame(scores,X22.test)
        predictors=as.matrix(test.predictors)
        data.test= data.frame(xx=test.predictors)
        
        yhat= predict(fit, newdata = data.test)
  
        ols_predic[i,d]=yhat
        
        MSE_OLS[i,d]=yhat-y.test}
    
      if(nrow(X11.train)>ncol(aux)){
        
        M=eppls(X11.train,X22.train,y.train,d)
        
        env_predic[i,d]=t(M$beta1)%*%t(X11.test-M$mu1)+(t(M$beta2))%*%t(as.matrix(X22.test-as.matrix(M$mu2)))+M$muY
        
        yhat=env_predic[i,d]
        
        MSE_EPLS[i,d]=yhat-y.test}
  
    }
  }
  
  # Calculate best d for relevant models
  d_PCA=which.min(apply(MSE_PCA^2,2,mean))
  d_PLS=which.min(apply(MSE_PLS^2,2,mean))
  d_EPLS=which.min(apply(MSE_EPLS^2,2,mean))
  
  # Calculate values for the table
  
  table61$PPLS.q1[index] <- round(apply(MSE_PLS^2,2,mean)[d_PLS], 2)
  table61$PPCA.q1[index] <- round(apply(MSE_PCA^2,2,mean)[d_PCA], 2)
  table61$PPLS.1[index] <- round(apply(MSE_PLS^2,2,mean)[1], 2)
  table61$PPCA.1[index] <- round(apply(MSE_PCA^2,2,mean)[1], 2)
  table61$LASSO[index] <- round(mean(apply(MSE_LASSO^2,2,mean)), 2)
  table61$FULL[index] <- round(apply(MSE_OLS^2,2,mean)[1], 2)
  table61$PPENV.q1[index] <- round(apply(MSE_EPLS^2,2,mean)[d_EPLS], 2)
  table61$PPENV.1[index] <- round(apply(MSE_EPLS^2,2,mean)[1], 2)
  
  index = index + 1
  
}

# End table generation
print(knitr::kable(table61, "simple", caption="Table 6.1"))




####### Figure 6.2
aux_2003<-aux[which(aux$year>=2003),]
aux_2003=aux_2003[-2,]
n=nrow(aux_2003)

# Using 24 indicators (to reduce)
X1=aux_2003[,c(4:9,11:13,15:17,19:21,23:25,27:29,31:33)]

X2=aux_2003[,c(34:39,41,43:53,55)] # with FDI, Domestic Credit, Age Dependency, Popul, Seconday School, Inflation and control countries 

y=as.numeric(unlist(aux_2003[,1]))

lasso_prediction=rep(0,n)
pls_prediction=rep(0,n)

for (i in 1:n) {
  #training
  X11.train=X1[-i,]
  X22.train=X2[-i,]
  y.train=y[-i]
  
  #testing
  X11.test=X1[i,]
  X22.test=X2[i,]
  y.test=y[i]
  
  ### Matrix X to reduce
  X11.train=as.matrix(X11.train)
  
  M=data.frame(X22.train)
  
  ### Regression of the continuous X in which we do not reduce
  A=lm(X11.train~.,data=M)
  
  # Save residuals
  res.train=data.frame(RES=A$residuals)
  
  M1=data.frame(X22.test)
  
  ### Predict on new data
  A_1=predict(A,newdata=M1)
  
  # Center training response variable
  y.train.center=y.train-mean(y.train)
  
  ## Perform pls on residuals in training
  m=plsr(y.train.center~.,ncomp=8,data=res.train)
  
  #Perform regression on non-reduced variables
  j=lm(y.train.center~.,data=M)
  
  #Predict non-reduced on new data
  j_2=predict(j,newdata=M1)
  
  #Final prediction on new data
  yhat=mean(y.train)+j_2+as.numeric(as.numeric(X11.test-A_1,nrow=1)%*%(as.numeric(m$coefficients[,1,8],ncol=1)))
  
  pls_prediction[i]=yhat
  
  
  
  # LASSO
  
  covariates=as.matrix(cbind(X11.train, X22.train))
  z=y.train
  cv_model <- cv.glmnet(covariates,z, alpha = 1)
  
  # Find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  best_model <- glmnet(covariates, z, alpha = 1, lambda = best_lambda)
  # best_model$df
  covariates.test=as.matrix(cbind(X11.test,X22.test))
  # Use lasso regression model to predict response value
  lasso_prediction[i] =predict(best_model, s = best_lambda, newx = covariates.test)
  
  
}


grayPalette <- c("#AAAAAA",  "#000000FF")

dfp_proj <- as.data.frame(cbind(pls_prediction,y))
 
ggplot(dfp_proj, aes(x = pls_prediction, y = y) )+ 
  geom_point(alpha = 0.7, size = 2) +
  theme_bw() +
  labs(x = "Leave-one-out fitted values from PPLS-7",
       y = "Y") +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.9, 0.9), 
        legend.title = element_blank(),
        axis.title = element_text(size=16),
        legend.text = element_text(size = 13)
  ) 

dfp_proj <- as.data.frame(cbind(lasso_prediction,y))
 
ggplot(dfp_proj, aes(x = lasso_prediction, y = y) )+ 
  geom_point(alpha = 0.7, size = 2) +
  theme_bw() +
  labs(x = "Leave-one-out fitted values from LASSO",
       y = "Y") +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.9, 0.9), 
        legend.title = element_blank(),
        axis.title = element_text(size=16),
        legend.text = element_text(size = 13)
  ) 
