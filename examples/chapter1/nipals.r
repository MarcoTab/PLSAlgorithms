rm(list=ls())
library(pls)
library(glmnet)
library(MASS)

best_comps = c(0,0,0)
rmses_11 = c(0,0,0)

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

# Find best number of components to use acording to CV
best_ncomp <- which.min(pls_fit$validation$adj)
best_comps[1] <- best_ncomp

cat("Found best number of components to be ", best_ncomp, "\n")
print(pls_fit$validation$adj)
cat("\n")

# Train model on best ncomps
gas <- plsr(y~x, ncomp=best_ncomp)


# Predicting in the test data
new=data.frame(x=xx)
new=t(data.frame(x=Xtest))
# Predict using best model
p=predict(gas,newdata=new)[,,best_ncomp]
pred_pls=p

ytest=ytest$V1

rmses_11[1] <- sqrt(mean((ytest-pred_pls)^2))

# Plot predicted vs observed responses regarding moisture in corn for PLS (part of Plot 1.1)
plot(ytest,pred_pls,xlab="Observed response, Y",ylab="Predicted response", main="Plot 1.1")



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

# Find best number of components to use according to CV
best_ncomp <- which.min(pls_fit$validation$adj)
best_comps[2] <- best_ncomp

cat("Found best number of components to be ", best_ncomp, "\n")
print(pls_fit$validation$adj)
cat("\n")

# Train model on best ncomps
gas <- plsr(y~x, ncomp=best_ncomp)

new=data.frame(x=xx)
new=t(data.frame(x=Xtest))

# Predict using best model
p=predict(gas,newdata=new)[,,best_ncomp]
pred_pls=p

ytest=ytest$V1

rmses_11[2] <- sqrt(mean((ytest-pred_pls)^2))

# Plot predicted vs observed responses regarding protein in meat for PLS (part of Plot 1.2)
plot(ytest,pred_pls,xlab="Observed response, Y",ylab="Predicted response", main="Plot 1.2")





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

# Find best number of components to use according to CV
best_ncomp <- which.min(pls_fit$validation$adj)
best_comps[3] <- best_ncomp

cat("Found best number of components to be ", best_ncomp, "\n")
print(pls_fit$validation$adj)
cat("\n")

# Train model on best ncomps
gas <- plsr(y~x, ncomp=best_ncomp)

new=data.frame(x=xx)
new=t(data.frame(x=Xtest))

# Predict using best model
p=predict(gas,newdata=new)[,,best_ncomp]
pred_pls=p

ytest=ytest$V1

rmses_11[3] <- sqrt(mean((ytest-pred_pls)^2))

# Plot predicted vs observed responses regarding tetracycline in serum for PLS (part of Plot 1.3)
plot(ytest,pred_pls,xlab="Observed response, Y",ylab="Predicted response", main="Plot 1.3")


titles = c("Moisture", "Protein", "Tetracycline")
cat("Table 1.1 (Only PLS)\n   Dataset | Components | RMSE\n")
for (x in 1:3) {
  cat(titles[x], " | ", best_comps[x], " | ", rmses_11[x], "\n")
}
