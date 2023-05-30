
library(matrixcalc)

### Given a positive semi-definite matrix A and another one S and d, tries to find Gamma of dimensions p x d with p dimension of A and S,
### Atilde, Omega, Omega_0 such that A = Gamma Atilde Gamma^T and S = Gamma Omega Gamma^T + Gamma_0 Omega_0 Gamma_0^T with Gamma and Gamma_0 orthogonal.
### Returns all of the and also returns the estimated A, S and pls beta.

# Given a positive semi-definite matrix and a vector, finds the orthogonal projection wrt the dot product generated with that matrix
proy_ort<-function(M,v){
    m=dim(M)[1]
    #if (!is.positive.semi.definite(M)) stop("error M is not definite semi-positive ")
    #if (dim(v)[1]!=dim(M)[1]) stop("error M no es positiva definida")
    proy=diag(rep(1,m))-v%*%solve(t(v)%*%M%*%v)%*%t(v)%*%M
}


# Given positive semi-definite matrices A and S, calculates the d first PLS projections.
# A is Sigma_xy sigma_y^-1 Sigma_xy
# S is Sigma_x
# Returns Gamma, Omega, beta.
pls_matrices<-function(A, S, d){
   #if (!is.positive.semi.definite(A)) stop("error A is not definite semi-positive")
   #if (!is.positive.semi.definite(S)) stop("error S is not definite semi-positive")
   V=matrix(0,dim(A)[1],d)
   V[,1]=eigen(A)$vectors[,1]
   if (d>1){
     for (i in 2:d){
                  aux=proy_ort(S,V[,(1:(i-1))])
                   V[,i]=eigen(t(aux)%*%A%*%aux,symmetric=TRUE)$vectors[,1]
     }
   }
   
   invisible(list(Gamma=V, Omega=t(V)%*%A%*%V, Atilde=V%*%t(V)%*%A%*%V%*%t(V),beta =V%*%solve(t(V)%*%S%*%V)%*%t(V)%*%A))
}

  

          