rm(list=ls())
library(pls)
library(glmnet)
library(MASS)

# Load in moisture in corn data
ycal <- read.table("data/chapter1/Parameters in corn/y1cal.txt", quote="\"", comment.char="")
ytest <- read.table("data/chapter1/Parameters in corn/y1test.txt", quote="\"", comment.char="")
Xcal <-  read.table("data/chapter1/Parameters in corn/Xcal.txt", quote="\"", comment.char="")
Xtest <- read.table("data/chapter1/Parameters in corn/Xtest.txt", quote="\"", comment.char="")
x=t(Xcal)
xx=t(Xtest)
y=ycal$V1

# Cross validate for moisture in corn data in the training data
pls_fit = plsr(y~x, scale = FALSE, validation = "CV", segment.type="consecutive")
summary(pls_fit)

# Find best number of components to use acording to CV
best_ncomp <- which.min(pls_fit$validation$adj)

# Train model on best ncomps
gas <- plsr(y~x, ncomp=best_ncomp)


# Predicting in the test data
new=data.frame(x=xx)
new=t(data.frame(x=Xtest))
# Predict using best model
p=predict(gas,newdata=new)[,,best_ncomp]
pred_pls=p

ytest=ytest$V1

# Plot predicted vs observed responses regarding moisture in corn for PLS (part of Plot 1.1)
plot(ytest,pred_pls,xlab="Observed response, Y",ylab="Predicted response")



# Load in protein in meat data
ycal <- read.table("data/chapter1/Parameters in meat/y3cal.txt", quote="\"", comment.char="")
ytest <- read.table("data/chapter1/Parameters in meat/y3test.txt", quote="\"", comment.char="")
Xcal <-  read.table("data/chapter1/Parameters in meat/Xcal.txt", quote="\"", comment.char="")
Xtest <- read.table("data/chapter1/Parameters in meat/Xtest.txt", quote="\"", comment.char="")
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

# Plot predicted vs observed responses regarding protein in meat for PLS (part of Plot 1.1)
plot(ytest,pred_pls,xlab="Observed response, Y",ylab="Predicted response")





# Load in tetracycline data
ycal <- read.table("data/chapter1/Tetracycline in serum/ycal.txt", quote="\"", comment.char="")
ytest <- read.table("data/chapter1/Tetracycline in serum/ytest.txt", quote="\"", comment.char="")
Xcal <-  read.table("data/chapter1/Tetracycline in serum/Xcal.txt", quote="\"", comment.char="")
Xtest <- read.table("data/chapter1/Tetracycline in serum/Xtest.txt", quote="\"", comment.char="")
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

# Plot predicted vs observed responses regarding tetracycline in serum for PLS (part of Plot 1.1)
plot(ytest,pred_pls,xlab="Observed response, Y",ylab="Predicted response")


