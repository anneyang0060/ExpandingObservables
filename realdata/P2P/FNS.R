ell_fun <- function(a) {
  log(log(exp(exp(1)) + a))
}

standardize_for_warmup <- function(X) {
  X <- as.matrix(X)
  cm <- colMeans(X)
  cs <- apply(X, 2, sd)
  cs[!is.finite(cs) | cs == 0] <- 1
  Xs <- sweep(X, 2, cm, "-")
  Xs <- sweep(Xs, 2, cs, "/")
  return(Xs)
}

cov_complexity_upper <- function(X_obs,
                                 grid_size,
                                 delta_sigma = 0.05,
                                 c_sigma = 2.5,
                                 q_cov = 0,
                                 R_sigma = NULL,
                                 radius_mult = 1,
                                 c_radius = NULL) {
  X_obs <- as.matrix(X_obs)
  N <- nrow(X_obs)
  a <- ncol(X_obs)
  
  S_hat <- crossprod(X_obs) / N
  
  X2 <- X_obs^2
  V_hat <- pmax(crossprod(X2) / N - S_hat^2, 0)
  
  log_fac <- log(2 * a^2 * grid_size / delta_sigma)
  Omega <- c_sigma * sqrt(V_hat * log_fac / N)
  
  S_tilde <- S_hat
  off_diag <- row(S_hat) != col(S_hat)
  S_tilde[off_diag & abs(S_hat) < Omega] <- 0
  diag(S_tilde) <- diag(S_hat)
  
  tau_max <- max(Omega[off_diag], na.rm = TRUE)
  if (!is.finite(tau_max)) tau_max <- max(Omega, na.rm = TRUE)
  
  if (is.null(R_sigma)) {
    R_sigma <- max(1, max(colSums(abs(S_tilde) > 0)))
  }
  if (is.null(c_radius)) {
    c_radius <- 2 + 3^(1 - q_cov)
  }
  
  r_sigma <- c_radius * R_sigma * tau_max^(1 - q_cov)
  M_hat <- as.numeric(norm(S_tilde, type = "1"))
  M_upper <- M_hat + radius_mult * r_sigma
  
  return(list(
    M_upper = M_upper,
    M_hat = M_hat,
    r_sigma = r_sigma,
    S_tilde = S_tilde,
    tau_max = tau_max,
    R_sigma = R_sigma
  ))
}

warmup_check <- function(X_obs,
                         grid_size,
                         delta_sigma = 0.05,
                         c_sigma = 2.5,
                         q_cov = 0,
                         R_sigma = NULL,
                         warm_const = 1,
                         radius_mult = 1,
                         standardize = TRUE) {
  X_obs <- as.matrix(X_obs)
  if (standardize) {
    X_obs <- standardize_for_warmup(X_obs)
  }
  
  N <- nrow(X_obs)
  a <- ncol(X_obs)
  
  cov_out <- cov_complexity_upper(
    X_obs = X_obs,
    grid_size = grid_size,
    delta_sigma = delta_sigma,
    c_sigma = c_sigma,
    q_cov = q_cov,
    R_sigma = R_sigma,
    radius_mult = radius_mult
  )
  
  rhs <- warm_const * ell_fun(a) * cov_out$M_upper^3 * (log(a))^2
  
  return(list(
    stop = (N >= rhs),
    N = N,
    a = a,
    rhs = rhs,
    M_upper = cov_out$M_upper,
    M_hat = cov_out$M_hat,
    r_sigma = cov_out$r_sigma,
    R_sigma = cov_out$R_sigma
  ))
}

truncation <- function(x, t){
  
  p <- length(x); y <- rep(0,p)
  
  for(j in 1:p){
    
    y[j] <- 0 + (x[j]-t[j])*(x[j]-t[j]>0) + (x[j]+t[j])*(x[j]+t[j]<0)
    
  }
  
  return(y)
}

# proximal(CXy_1, CXX_1, lam_vec, beta)
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
