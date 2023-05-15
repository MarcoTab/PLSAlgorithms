rm(list=ls())
library(pls)
library(glmnet)
library(MASS)

# Load in moisture in corn data
ycal <- read.table("mvc1_data/Parameters in corn/y1cal.txt", quote="\"", comment.char="")
ytest <- read.table("mvc1_data/Parameters in corn/y1test.txt", quote="\"", comment.char="")
Xcal <-  read.table("mvc1_data/Parameters in corn/Xcal.txt", quote="\"", comment.char="")
Xtest <- read.table("mvc1_data/Parameters in corn/Xtest.txt", quote="\"", comment.char="")
x=t(Xcal)
xx=t(Xtest)
y=ycal$V1

# Cross validate for moisture in corn data
pls_fit = plsr(y~x, scale = FALSE, validation = "CV", segment.type="consecutive")
summary(pls_fit)

# Find best number of components to use acording to CV
best_ncomp <- which.min(pls_fit$validation$adj)

# Train model on best ncomps
gas <- plsr(y~x, ncomp=best_ncomp)

new=data.frame(x=xx)
new=t(data.frame(x=Xtest))
# Predict using best model
p=predict(gas,newdata=new)[,,best_ncomp]
pred_pls=p

ytest=ytest$V1

# Plot predicted vs observed responses regarding moisture in corn
plot(ytest,pred_pls,xlab="Observed response, Y",ylab="Predicted response")





# Load in protein in meat data
ycal <- read.table("mvc1_data/Parameters in meat/y3cal.txt", quote="\"", comment.char="")
ytest <- read.table("mvc1_data/Parameters in meat/y3test.txt", quote="\"", comment.char="")
Xcal <-  read.table("mvc1_data/Parameters in meat/Xcal.txt", quote="\"", comment.char="")
Xtest <- read.table("mvc1_data/Parameters in meat/Xtest.txt", quote="\"", comment.char="")
x=t(Xcal)
xx=t(Xtest)
y=ycal$V1

# Cross validate for protein in meat data
pls_fit = plsr(y~x, scale = FALSE, validation = "CV", segment.type="consecutive")
summary(pls_fit)

# Find best number of components to use acording to CV
best_ncomp <- which.min(pls_fit$validation$adj)

# Train model on best ncomps
gas <- plsr(y~x, ncomp=best_ncomp)

new=data.frame(x=xx)
new=t(data.frame(x=Xtest))

# Predict using best model
p=predict(gas,newdata=new)[,,best_ncomp]
pred_pls=p

ytest=ytest$V1

# Plot predicted vs observed responses regarding protein in meat
plot(ytest,pred_pls,xlab="Observed response, Y",ylab="Predicted response")





# Load in tetracycline data
ycal <- read.table("mvc1_data/Tetracycline in serum/ycal.txt", quote="\"", comment.char="")
ytest <- read.table("mvc1_data/Tetracycline in serum/ytest.txt", quote="\"", comment.char="")
Xcal <-  read.table("mvc1_data/Tetracycline in serum/Xcal.txt", quote="\"", comment.char="")
Xtest <- read.table("mvc1_data/Tetracycline in serum/Xtest.txt", quote="\"", comment.char="")
x=t(Xcal)
xx=t(Xtest)
y=ycal$V1

# Cross validate for tetracycline in serum data
pls_fit = plsr(y~x, scale = FALSE, validation = "CV", segment.type="consecutive")
summary(pls_fit)

# Find best number of components to use acording to CV
best_ncomp <- which.min(pls_fit$validation$adj)

# Train model on best ncomps
gas <- plsr(y~x, ncomp=best_ncomp)

new=data.frame(x=xx)
new=t(data.frame(x=Xtest))

# Predict using best model
p=predict(gas,newdata=new)[,,best_ncomp]
pred_pls=p

ytest=ytest$V1

# Plot predicted vs observed responses regarding tetracycline in serum
plot(ytest,pred_pls,xlab="Observed response, Y",ylab="Predicted response")


