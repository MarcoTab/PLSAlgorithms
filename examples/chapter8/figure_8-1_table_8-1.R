rm(list=ls())
library(caret)

library(matrixcalc)

library(mosaic)
library(ggplot2)

marcewithout <- read.table("data/chapter8/marcewithout.txt", quote="\"", comment.char="")
source("examples/chapter8/lib/quadratic_with_cross_validation.R")
source("examples/chapter8/lib/pls_matrices.R")

Y=marcewithout[,1]+1

xx=marcewithout[,2:14]

X=matrix(0,ncol=ncol(xx),nrow=nrow(xx))
for (i in 1:ncol(xx)){
  for (j in 1:nrow(xx)){
    X[j,i]=as.numeric(xx[j,i])
  }
}

Dmax=11 

K=10
h=seq(0,1,.2) 

quad_1=matrix(0,ncol=length(h),nrow=Dmax)
quad_2=matrix(0,ncol=length(h),nrow=Dmax)
 
set.seed(2)

for (d in 1:Dmax){
  for (s in 1:length(h)){
    alp=h[s]

    aux=quadratic.discrimination(X,Y,K,d,alp)
  
    li_1=aux$predict_quad_1
    li_2=aux$predict_quad_2
   
    quad_1[d,s]=mean((li_1-Y)^2)
    quad_2[d,s]=mean((li_2-Y)^2)
  }
}

M1=which(quad_1 == min(quad_1), arr.ind = TRUE)

M2=which(quad_2 == min(quad_2), arr.ind = TRUE)

## Run data to make table 8.1

d=M1[1]
alp=h[M1[2]]
aux=quadratic.discrimination(X,Y,nrow(X),d,alp)
anq1 <- confusionMatrix(as.factor(aux$predict_quad_1),as.factor(Y))$overall[1]

d=M2[1]
alp=h[M2[2]]
aux=quadratic.discrimination(X,Y,nrow(X),d,alp)
anq2 <- confusionMatrix(as.factor(aux$predict_quad_2),as.factor(Y))$overall[1]

#without reduction
d=13
alp=1
aux=quadratic.discrimination(X,Y,nrow(X),d,alp)
qda_no_red <- confusionMatrix(as.factor(aux$predict_quad_1),as.factor(Y))$overall[1]


cat("Dim. Reduction Method |  u |   ϕqda | Accuracy (%)\n")
cat("----------------------|----|--------|-------------\n")
cat("AN-Q2 (λ=0.6)         | 6  | (8.15) |", 100*anq2, "\n")
cat("AN-Q1                 | 3  | (8.14) |", 100*anq1, "\n")
cat("QDA, no reduction     | 13 | (8.4)  |", 100*qda_no_red, "\n")

##### Figure 8.1

holdout=X
y.holdhout=Y

S_X <- cov(holdout)
uY <- unique(y.holdhout) %>% set_names()
S_k <- lapply(uY, function(k) cov(holdout[y.holdhout == k, , drop = FALSE])) 
cov.mean <- lapply(uY, function(k) S_k[[k]] * mean(y.holdhout==k)) %>%
  Reduce('+', .)
cov.ss1 <- lapply(uY, function(k) (S_X - S_k[[k]])^2 * mean(y.holdhout==k)) %>%
  Reduce('+', .)
cov.ss2 <- lapply(uY, function(k) (cov.mean - S_k[[k]])^2 * mean(y.holdhout==k)) %>%
  Reduce('+', .)


d=M1[1]
alp=h[M1[2]]

dif_d_delta2=(S_X-cov.mean)%*%t(S_X-cov.mean)

A_al1=cov.ss1

AUX=pls_matrices(A_al1, cov.mean, d)


d=M2[1]
alp=h[M2[2]]
A_al_para=alp*dif_d_delta2+(1-alp)*cov.ss2

AUX_para=pls_matrices(A_al_para, cov.mean, d)

datos_proy_2=X%*%AUX_para$Gamma
datos_proy_1=X%*%AUX$Gamma

Y=marcewithout$V1+1

dfp_proj <- as.data.frame(cbind(datos_proy_2))
origins <- c("Plane", "Car", "Bird")
dfp_proj$origins <- origins[Y] 

ggplot(dfp_proj, aes(x = V1, y = V2, shape = origins, color = origins)) + 
  geom_point(alpha = 0.7, size = 4) +
  theme_bw() +
  labs(x = "First projection using AN-Q2",
       y = "Second projection using AN-Q2")  +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.9, 0.9), 
        legend.title = element_blank(),
        axis.title = element_text(size=16),
        legend.text = element_text(size = 13)
  ) 

Y=marcewithout$V1+1

dfp_proj <- as.data.frame(cbind(datos_proy_1))
origins <- c("Plane", "Car", "Bird")
dfp_proj$origins <- origins[Y] 

ggplot(dfp_proj, aes(x = V1, y = V2, shape = origins, color = origins)) + 
  geom_point(alpha = 0.7, size = 4) +
  theme_bw() +
  labs(x = "First projection using AN-Q1",
       y = "Second projection using AN-Q1")  +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.9, 0.9), 
        legend.title = element_blank(),
        axis.title = element_text(size=16),
        legend.text = element_text(size = 13)
  ) 
