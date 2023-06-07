rm(list=ls())
library(caret)
library(matrixcalc)

source("examples/chapter8/pls_matrices.R")
source("examples/chapter8/quadratic_with_cross_validation.R")

 
MIR_Fruit_purees <- read.csv("data/chapter8/MIR_Fruit_purees.csv", header=FALSE)


Y=as.factor(t(MIR_Fruit_purees[1,2:984]))
Y=as.numeric(Y)
Y=matrix(Y,ncol=1)
xx=t(MIR_Fruit_purees[2:236,2:984])

X=matrix(0,ncol=ncol(xx),nrow=nrow(xx))
for (i in 1:ncol(xx)){
  for (j in 1:nrow(xx)){
    X[j,i]=as.numeric(xx[j,i])
  }
}

Dmax=20 

K=10

h=seq(0,1,.2) 

quad_1=matrix(0,ncol=length(h),nrow=Dmax)
quad_2=matrix(0,ncol=length(h),nrow=Dmax)

for (d in 1:Dmax){
  
  for (s in 1:length(h)){
    cat(d, " ")
    cat(s, "\n")
    alp=h[s]
    
    aux=quadratic.discrimination(X,Y,K,d,alp)
    
    li_1=aux$predict_quad_1
    li_2=aux$predict_quad_2
    
    
    quad_1[d,s]=mean((li_1-Y)^2)
    quad_2[d,s]=mean((li_2-Y)^2)
    
    
  }
  
}

###error leave one out

M1=which(quad_1 == min(quad_1), arr.ind = TRUE)
aux=quadratic.discrimination(X,Y,nrow(X),M1[1],h[M1[2]])
anq1 = confusionMatrix(as.factor(aux$predict_quad_1),as.factor(Y))$overall[1]

M2=which(quad_2 == min(quad_2), arr.ind = TRUE)
aux=quadratic.discrimination(X,Y,nrow(X),M2[1],h[M2[2]])
anq2 = confusionMatrix(as.factor(aux$predict_quad_2),as.factor(Y))$overall[1]

cat("Table 8.3")
cat("Dataset | AN-Q1 | AN-Q2\n")
cat("--------|-------|------\n")
cat("Fruit   |", 100*anq1,"|", 100*anq2,"\n")


## Figure 8.2
y.holdhout=Y
holdout=X

S_X <- cov(holdout)
uY <- unique(y.holdhout) %>% set_names()
S_k <- lapply(uY, function(k) cov(holdout[y.holdhout == k, , drop = FALSE])) 
cov.mean <- lapply(uY, function(k) S_k[[k]] * mean(y.holdhout==k)) %>%
  Reduce('+', .)
cov.ss1 <- lapply(uY, function(k) (S_X - S_k[[k]]) %*% t(S_X - S_k[[k]]) * mean(y.holdhout==k)) %>%
  Reduce('+', .)
cov.ss2 <- lapply(uY, function(k) (cov.mean - S_k[[k]]) %*% t(cov.mean - S_k[[k]]) * mean(y.holdhout==k)) %>%
  Reduce('+', .)

dif_d_delta2=(S_X-cov.mean)%*%t(S_X-cov.mean)

## for quad_1
A_al1=cov.ss1
alp=h[M1[2]]
d=M1[1]
AUX=pls_matrices(A_al1, cov.mean, d)

## for quad_2
alp=h[M2[2]]
d=M2[1]
A_al_para=alp*dif_d_delta2+(1-alp)*cov.ss2
AUX_para=pls_matrices(A_al_para, cov.mean, d)

datos_proy_2=X%*%AUX_para$Gamma
datos_proy_1=X%*%AUX$Gamma

library(mosaic)
library(ggplot2)

dfp_proj_2 <- as.data.frame(datos_proy_2)
origins <- c("NON-Strawberry", "Strawberry")
dfp_proj_2$origins <- origins[Y] 

ggplot(dfp_proj_2, aes(x = V1, y = V2, shape = origins, color = origins)) + 
  geom_point(alpha = 0.7, size = 2.5) +
  theme_bw() +
  labs(x = "First projection using AN-Q2",
       y = "Second projection using AN-Q2") +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.8, 0.9), 
        legend.title = element_blank(),
        axis.title = element_text(size=16),
        legend.text = element_text(size = 13)
  ) 


dfp_proj_1 <- as.data.frame(datos_proy_1)
origins <- c("NON-Strawberry", "Strawberry")
dfp_proj_1$origins <- origins[Y] 

ggplot(dfp_proj_1, aes(x = V1, y = V2, shape = origins, color = origins)) + 
  geom_point(alpha = 0.7, size = 2.5) +
  theme_bw() +
  labs(x = "First projection using AN-Q1",
       y = "Second projection using AN-Q1") +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.5, 0.9), 
        legend.title = element_blank(),
        axis.title = element_text(size=16),
        legend.text = element_text(size = 13)
  ) 

