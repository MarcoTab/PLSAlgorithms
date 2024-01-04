#######      Function rrr (Reduced-Rank Regression)


rrr <- function(x, y, k=0, type="Identity", rrr.plot=TRUE,
                missing="omit", title="") {
  
  # This file defines the function, rrr.  
  # rrr returns the the A, B, C, and mu matrices for the t values.
  #
  # Chuck Miller -- 6 Sept 2007
  #
  # NOTE: This Reduced-Rank Regression function is based on Alan Izenman's Text,
  # Modern Multivariate Statistical Techniques: Regression, Classification, 
  # and Manifold Learning, 
  # published by Springer (2008).
  #
  # This function removes missing values when missing is set to "omit" (default).
  #
  # This function can be used to perform Reduced-Rank Regression, PCA, CVA and LDA.
  # Set type="pca" for PCA, "cva" for CVA or "lda" for LDA
  #
  # k is an arbitrary value used to avoid singular matrices.  
  # It is used in a fashion similar to ridge regression.  
  # It is useful to present a rank trace plot that has a elbow shape.  
  
  #Put data into form where variables are in rows as in the Text
  x <- t(x)
  y <- t(y)
  
  #Remove missing values
  DF <- data.frame(t(y), t(x))
  if (missing=="omit"){
    y <- t(na.omit(DF)[1:dim(y)[1]])
    x <- t(na.omit(DF)[(dim(y)[1]+1):(dim(y)[1]+dim(x)[1])])
  }
  
  #Checking for missing values
  ok <- complete.cases(DF)
  if (!all(ok)) stop("missing values in data")
  
  #Check length of x and y
  n <- ncol(x)
  n2 <- min(ncol(y),length(y))
  if (n != n2)
    stop("lengths of x and y are different")
  
  #Linear Discriminant Values
  if (type=="LDA"||type=="lda") {
    yf <- as.factor(y)
    lev <- levels(yf)
    y <- matrix(0,n2, length(lev), dimnames=list(1:n2, lev))
    for(i in 1:(length(lev)-1)) {
      y[,i] <- yf==lev[i]
    }
    y <- t(y)
  }
  
  #Calculate mux and muy
  mux <- apply(x, 1, mean)
  muy <- apply(y, 1, mean)
  
  #Centering the x and y values
  xc <- apply(x, 1, function(x) (x-mean(x)))
  yc <- apply(y, 1, function(y) (y-mean(y)))
  
  #Place xc and yc back into form similar to book (6.82)
  xc <- t(xc); yc <- t(yc)
 
  #Compute sample sigma matrices
  Sxx <- (1/n)*tcrossprod(xc)+k*diag(nrow(xc))
  Syy <- (1/n)*tcrossprod(yc)+k*diag(nrow(yc))
  Syx <- (1/n)*tcrossprod(yc, xc)
  Sxy <- t(Syx)
 # print(Sxx)
 # print(Syx)
  #Calculate the max value of t
  tmax <- min(dim(x)[[1]], dim(y)[[1]])
  
  #Gamma Function (default value is identity matrix)
  
  if ((type=="identity")||(type=="Identity")) g <- diag(1, dim(y)[1]) else
    if ((type=="cva")||(type=="CVA")) g <- solve(Syy) else
      if ((type=="pca")||(type=="PCA")) g <- diag(1, dim(y)[1]) else
        if ((type=="lda")||(type=="LDA")) g <- solve(Syy) else
          stop("type is only valid for 'identity', 'cva', 'pca' or 'lda'")
  
  ev <- eigen(g)[[2]]; ed <- eigen(g)[[1]]
  g.5 <- ev%*%sqrt(diag(ed))%*%t(ev)
  
  #Compute eigenvalues and eigenvectors for (6.90)
  lambdaj <- eigen((g.5)%*%Syx%*%solve(Sxx)%*%Sxy%*%(g.5))[[1]]
  vj <- eigen((g.5)%*%Syx%*%solve(Sxx)%*%Sxy%*%(g.5))[[2]]
  
  #Calculate A (6.88), B (6.89), C (6.91), and mu (6.76)
  A <- list();B <- list();C <- list();mu <- list(); See <- list()
  for (j in 1:tmax) {
    A[[j]] <- solve(g.5)%*%vj[,1:j]
    B[[j]] <- t(vj[,1:j])%*%(g.5)%*%Syx%*%solve(Sxx)
    C[[j]] <- A[[j]] %*% B[[j]]
    mu[[j]] <-  muy - C[[j]] %*% mux
    See[[j]] <- (1/n)*(yc - (C[[j]]%*%xc))%*%t(yc - (C[[j]]%*%xc))
  }
  
  #Rank Trace -- used to determine value of t to use
  if (rrr.plot=="TRUE"){
    dC <- NULL; dE <- NULL
    for (j in 1:tmax){
      dC[j] <- (sqrt(sum(diag(tcrossprod(C[[tmax]]-C[[j]])))))/
        (sqrt(sum(diag(tcrossprod(C[[tmax]])))))
      
      dE[j] <- (sqrt(sum(diag(tcrossprod(See[[tmax]]-See[[j]])))))/
        (sqrt(sum(diag(tcrossprod(See[[tmax]]-Syy)))))
    }
    dC <- c(1, dC)
    dE <- c(1, dE)
    names(dC) <- c(0, 1:tmax)
    names(dE) <- c(0, 1:tmax)
    plot(dC, dE, col="red", type="b", main="Tobacco data, Gamma=Identity, k=0") 
    #	main=paste(title, " Rank Trace, Type=", type, ", k=", k))
    #  lines(dC, dE)
    text(dC, dE, names(dC), pos=3)
  }  #End of Rank Trace
  
  #Calculating principal component scores
  score <- B[[tmax]]%*%xc
  
  #Plot of principal component scores and scree plot for pca analysis
  if ((type=="PCA"||type=="pca")&&(rrr.plot=="TRUE")) {
    npcs <- min(10, tmax)
    pc.txt <- round(cumsum(lambdaj)/sum(lambdaj), digits=3)
    plot(1:npcs, lambdaj[1:npcs], type="b", bty="l", pch=16,
         main="Scree Plot", xlab="Components", ylab="Variance")
    lines(1:npcs, lambdaj[1:npcs])
    text(1:npcs, lambdaj[1:npcs], pc.txt[1:npcs], pos=3)
    
    plot(score[1,], score[2,], main="Principal Component Scores",
         xlab="Principal Component 1", ylab="Principal Component 2")
    
  }
  
  #Plots for Canonical Variate Analysis
  cva.yscore <- NULL
  cva.xscore <- NULL
  cva.corr <- NULL
  if (type=="CVA"||type=="cva"){
    cva.yscore <- solve(A[[tmax]])%*%y
    cva.xscore <- B[[tmax]]%*%x
    cva.corr <- lambdaj
    ##  if (rrr.plot=="TRUE"){
    ##  plot(cva.xscore[1,], cva.yscore[1,], main="First Canonical Variate Pair",
    ##       xlab="X Canonical Variate", ylab="Y Canonical Variate")
    ##  plot(cva.xscore[2,], cva.yscore[2,], main="Second Canonical Variate Pair",
    ##       xlab="X Canonical Variate", ylab="Y Canonical Variate")
    ##  }
  }
  
  #Plots for Linear Discriminant Analysis
 # if (type=="LDA"||type=="lda"){
 #   cva.yscore <- solve(A[[tmax]])%*%y
 #   cva.xscore <- B[[tmax]]%*%x
 #   cva.corr <- lambdaj
 #   if (rrr.plot=="TRUE"){
 #     plot(cva.xscore[1,], cva.xscore[2,], xlab="LD1", ylab="LD2", type="n",
 #          main="Linear Discriminant Analysis")
      #text(cva.xscore[1,], cva.xscore[2,], as.numeric(yf), col=as.numeric(yf))
 #     text(cva.xscore[1,], cva.xscore[2,], abbreviate(yf,1), col=as.numeric(yf))
      #text(cva.yscore[1,], cva.yscore[2,], "*", col=as.numeric(yf))
  #  }
  #}
  
  invisible(list (A=A, B=B, C=C, mu=mu, score=score,
                  cva.xscore=cva.xscore, cva.yscore=cva.yscore,
                  cva.corr=cva.corr))
} #End of rrr function

#  Note:  The values returned (eg. A, B, ...) are in text format
#         ie. variables are in the rows.
