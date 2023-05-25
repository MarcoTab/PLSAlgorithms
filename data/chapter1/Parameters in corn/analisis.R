rm(list=ls())
library(pls)
setwd("~/Dropbox/pls-book/PLSbook/Data-analysis/mvc1_data_Frabricio/mvc1_data/Parameters in corn")
#parameters in corn
ycal <- read.table("y1cal.txt", quote="\"", comment.char="")

ytest <- read.table("y1test.txt", quote="\"", comment.char="")

Xcal <-  read.table("Xcal.txt", quote="\"", comment.char="")

Xtest <- read.table("Xtest.txt", quote="\"", comment.char="")




library(glmnet)



library(MASS)




x=t(Xcal)
xx=t(Xtest)
y=ycal$V1





xx_cent=scale(matrix(xx,ncol=ncol(xx)),center=apply(x,2,mean),scale=TRUE)
x_cent=scale(x,center=TRUE,scale=FALSE)
w = ginv(t(x_cent) %*% x_cent) %*% t(x_cent) %*% (y-mean(y))
pred_ols = xx_cent %*% w + mean(y)

set.seed(1)
pls_fit = plsr(y~x, scale = TRUE, validation = "CV")
summary(pls_fit)
validationplot(pls_fit, val.type = "MSEP")

#y1 9 compontnse
#y2 9 componetes
#y4 9 compontens
#y3 9 compontens
ncomp=9
gas <- plsr(y~x, ncomp = ncomp)

new=data.frame(x=xx)
new=t(data.frame(x=Xtest))
p=predict(gas,newdata=new)[,,ncomp]

pred_pls=p





#perform k-fold cross-validation to find optimal lambda value
cv_model <- cv.glmnet(x, y, alpha = 1)

#find optimal lambda value that minimizes test MSE
best_lambda <- cv_model$lambda.min



best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
best_model$df
new=t(Xtest)


#use lasso regression model to predict response value
pred_lasso=predict(best_model, s = best_lambda, newx = new)







pred_ols_in = x_cent %*% w + mean(y)
pred_pls_in=predict(gas,newdata=x)[,,ncomp]
pred_lasso_in=predict(best_model, s = best_lambda, newx = x)





ytest=ytest$V1




plot(ytest,pred_pls,xlab="Observed response, Y",ylab="Predicted response")

points(ytest,pred_lasso,col='red',pch='+')
#points(ytest,pred_ols,col='red')
title("Y vs fitted, green is lasso")
sqrt(mean((ytest-pred_pls)^2))
sqrt(mean((ytest-pred_lasso)^2))
sqrt(mean((ytest-pred_ols)^2))

mean((ycal$V1-pred_pls_in)^2)
mean((ycal$V1-pred_lasso_in)^2)
mean((ycal$V1-pred_ols_in)^2)

plot(pred_pls,pred_lasso)
title("Pred pls vs Pred lasso in test")
#####lasso is better
#####lasso is better
#####lasso is better


For y
> mean((ytest-pred_pls)^2)
[1] 0.0003359925
> mean((ytest-pred_lasso)^2)
[1] 0.01298456


fpr y1
[1] 0.0009054318
> mean((ytest-pred_lasso)^2)
[1] 0.01298456

for y2
mean((ytest-pred_pls)^2)
[1]  0.00192415
> mean((ytest-pred_lasso)^2)
[1] 0.007967773

for y3
> mean((ytest-pred_pls)^2)
[1] 0.009163168
> mean((ytest-pred_lasso)^2)
[1] 0.03831648

mean((ytest-pred_pls)^2)
[1]  0.03254397
> mean((ytest-pred_lasso)^2)
[1] 0.1649143