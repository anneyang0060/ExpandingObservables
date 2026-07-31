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


ell_fun <- function(a) {
  log(log(exp(exp(1)) + a))
}

make_warm_grid <- function(min_blocks = 2, max_blocks, ratio = 1.25) {
  if (max_blocks <= min_blocks) return(max_blocks)
  grid <- unique(ceiling(min_blocks * ratio^(0:1000)))
  grid <- grid[grid <= max_blocks]
  grid <- unique(c(grid, max_blocks))
  return(grid)
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
                         radius_mult = 1) {
  X_obs <- as.matrix(X_obs)
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

make_chain_crossblock_cov <- function(idx_abs,
                                      rho,
                                      block_width = 500,
                                      m_dep = 10) {
  idx_abs <- as.integer(idx_abs)
  p_now <- length(idx_abs)
  
  if (p_now == 0) {
    return(matrix(0, 0, 0))
  }
  
  Sigma <- diag(1, p_now)
  
  # Feature-increment block id:
  # 1:500 -> 0, 501:1000 -> 1, etc.
  block_id <- floor((idx_abs - 1) / block_width)
  
  # Position inside each feature block:
  # 1, ..., block_width
  local_pos <- ((idx_abs - 1) %% block_width) + 1
  
  # Only the first m_dep positions in each block are allowed to have
  # cross-block dependence. This matches the construction of idx_dep.
  for (ell in seq_len(m_dep)) {
    ids <- which(local_pos == ell)
    
    if (length(ids) >= 2) {
      D <- abs(outer(block_id[ids], block_id[ids], "-"))
      Sigma[ids, ids] <- rho^D
    }
  }
  
  diag(Sigma) <- 1
  Sigma <- (Sigma + t(Sigma)) / 2
  
  return(Sigma)
}


cv_hard_const <- function(Cy, CXy, CXX, c_set, N_eff, eps = 1e-8) {
  
  K_cv <- length(Cy)
  p <- dim(CXX)[1]
  cv_err <- rep(0, length(c_set))
  
  for (k_cv in 1:K_cv) {
    
    ## training sufficient statistics
    Cy_train <- sum(Cy[-k_cv])
    CXy_train <- rowSums(CXy[, -k_cv, drop = FALSE])
    CXX_train <- apply(CXX[, , -k_cv, drop = FALSE], c(1, 2), sum)
    
    ## validation sufficient statistics
    Cy_val <- Cy[k_cv]
    CXy_val <- CXy[, k_cv]
    CXX_val <- CXX[, , k_cv]
    
    gram_train <- ginv(CXX_train)
    beta_train <- as.vector(gram_train %*% CXy_train)
    
    ## Since Cy_train, CXX_train, CXy_train are normalized by total N_eff,
    ## res_train is also normalized by total N_eff. Convert it to the
    ## training-fold residual variance by dividing by the training fraction.
    train_frac <- (K_cv - 1) / K_cv
    
    res_train <- as.numeric(
      Cy_train -
        2 * crossprod(CXy_train, beta_train) +
        crossprod(beta_train, CXX_train %*% beta_train)
    )
    
    sigma_train <- sqrt(max(res_train / train_frac, eps))
    
    se_train <- sqrt(pmax(diag(gram_train), 0)) *
      sigma_train / sqrt(N_eff)
    
    for (j in seq_along(c_set)) {
      
      beta_thr <- beta_train
      beta_thr[abs(beta_thr) <= c_set[j] * se_train] <- 0
      
      cv_err[j] <- cv_err[j] + as.numeric(
        Cy_val -
          2 * crossprod(CXy_val, beta_thr) +
          crossprod(beta_thr, CXX_val %*% beta_thr)
      )
    }
  }
  
  c_opt <- c_set[which.min(cv_err)]
  
  return(list(
    c_opt = c_opt,
    cv_err = cv_err,
    c_set = c_set
  ))
}
