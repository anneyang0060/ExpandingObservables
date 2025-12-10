truncation <- function(x, t){
  
  p <- length(x); y <- rep(0,p)
  
  for(j in 1:p){
    
    y[j] <- 0 + (x[j]-t[j])*(x[j]-t[j]>0) + (x[j]+t[j])*(x[j]+t[j]<0)
    
  }
  
  return(y)
}


proximal <- function(CXy, CXX, lambda, beta0, gamma = 0.001){
  
  beta <- truncation(beta0 + gamma * (CXy - as.vector(CXX%*%beta0)), gamma*lambda)
  flag <- 1
  
  while(mean((beta-beta0)^2)>1e-06 & flag<1000){
    
    beta0 <- beta
    beta <- truncation(beta0 + gamma * (CXy - as.vector(CXX%*%beta0)), gamma*lambda)
    flag <- flag + 1
    
  }
  
  return(beta)
}

scale <- function(Cy, CXy, CXX, sigma0, beta0, t){

  lambda0 <- rep(sqrt(2*log(length(CXy))/t),length(CXy))
  beta <- proximal(CXy, CXX, sigma0*lambda0, beta0)
  sigma <- as.vector((Cy-2*matrix(CXy,nrow=1)%*%matrix(beta,ncol=1)
            +matrix(beta,nrow=1)%*%CXX%*%matrix(beta,ncol=1))^(1/2))
  flag <- 1
  
  while(abs(sigma-sigma0)>1e-04 & flag<100){
    
    sigma0 <- sigma
    beta <- proximal(CXy, CXX, sigma0*lambda0, beta)
    sigma <- as.vector((Cy-2*matrix(CXy,nrow=1)%*%matrix(beta,ncol=1)
              +matrix(beta,nrow=1)%*%CXX%*%matrix(beta,ncol=1))^(1/2))
    flag <- flag + 1

  }
  
  return(sigma)
}

tunning <- function(sigma,CXX){
  
  # Xeps_max_approx <- c()
  p <- dim(CXX)[1]
  
  # for(nb in 1:Nb){
    
    zeta <- rmvn(1, rep(0,p), CXX+diag(1e-08,p))
    # Xeps_max_approx <- c(Xeps_max_approx, max(abs(zeta))*sigma)
    Xeps_max_approx <- max(abs(zeta))*sigma
    
  # }
  
  lambda <- Xeps_max_approx/sqrt(t*n) * 2
  
  return(lambda)
}


hardthresh <- function(beta, CXX, CX, sigma_hat, idx_thre, nb){
  
  beta1 <- beta[idx_thre]
  CXX1 <- CXX[idx_thre, idx_thre]
  CX1 <- CX[idx_thre]
  p <- length(beta1)
  thre_b <- c()
  for(i in 1:nb){
    e <- rnorm(p)
    thre_b <- cbind(thre_b, abs(ginv(CXX1)%*%(CX1*sigma_hat*e)))
  }
  thre <- c()
  for(i in 1:p){
    thre <- c(thre, sort(thre_b[i,])[nb*0.95])
  }
  
  for(i in 1:length(beta1)){
    # if(beta1[i] < thre[i]){beta1[i] <- 0}
    if(beta1[i] < thre[i]){beta1[i] <- 0}
  }
  beta[idx_thre] <- beta1
  return(beta)
}


proximal_OLL <- function(CXy, CXX, Obs, lambda, beta0, gamma = 0.001){
  
  beta_last <- beta0
  beta <- truncation(beta0 + gamma * (CXy - as.vector(CXX%*%beta_last)
                                      +1/n0*t(Obs[,-ncol(Obs)])%*%Obs[,-ncol(Obs)]%*%(beta_last-beta0)), gamma*lambda)
  flag <- 1
  
  while(mean((beta-beta0)^2)>1e-05 & flag<100){
    
    beta0 <- beta
    beta <- truncation(beta0 + gamma * (CXy - as.vector(CXX%*%beta_last)
                                        +1/n0*t(Obs[,-ncol(Obs)])%*%Obs[,-ncol(Obs)]%*%(beta_last-beta0)), gamma*lambda)
    flag <- flag + 1
    
  }
  
  return(beta)
}
