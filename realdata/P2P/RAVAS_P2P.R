rm(list=ls())
gc()
setwd('C:/Users/annie/Downloads/accepted_2007_to_2018Q4.csv')


library(glmnet)
library(MASS)
library(fastDummies)
source('FNS.R')


K_train <- 2250 


####### training ############
{
  s <- 1; I_set <- list(); W_set <- c()
  Obs <- c(); lam_set <- c()
  sigma_o <- c(); p_old <- 0; K <- 0
  st <- 0; pt <- 0
  idx_imp <- c(); idx_abs <- c(); idx_del <- c()
  N_obs_old <- 0; p_obs_old <- 0
  err_test_o <- c(); k_o <- c()
  pt_o <- c(); err_y_o <- c(); time_o <- c(); resd <- c()
  beta_o <- list(); R2_o <- c(); y_hat <- c(); y_real <- c()
  sigma_hat <- 1; years <- c(); sigma_hat0 <- 0
  K_cv <- 2;  C_lam_opt0 <- 0; cv_num <- 0
}

## load data block
for(k in 1:K_train)
# k <- 1
{
  
  print(k)
  {
    # load
    file_cur <- paste0('P2Pdata_',k,'.csv')
    data_cur <- read.csv(file_cur)
    data_cur <- data_cur[,-c(1,2,3)]
    # response
    y <- data_cur$int_rate
    # covariate
    data_cur <- subset(data_cur, select = -c(int_rate))
    print(paste0('ncol data_cur =',ncol(data_cur)))
    X <- as.matrix(data_cur)
    p0 <- ncol(X)
    idx_del <- which(sapply(1:nrow(X), function(i){sum(is.na(X[i,]))})>0)
    if(length(idx_del)>0){
      X <- X[-idx_del, ]
      y <- y[-idx_del]
    }
    p_interact0 <- which(names(data_cur)=='home_ownership_MORTGAGE')
    j1_indices <- 1:(p_interact0-1)
    j2_indices <- p_interact0:p0
    combinations <- expand.grid(j1 = j1_indices, j2 = j2_indices)
    interaction_terms <- X[, combinations$j1] * X[, combinations$j2]
    X <- cbind(X, interaction_terms)    
    rm(data_cur)
    
    # data block features
    n_blk <- nrow(X)
    if(n_blk==0){
      print('n_blk=0!')
      next
    }
    p_obs <- ncol(X)
    print(paste0('p_obs=',p_obs))
    delta_p <- p_obs - p_obs_old
    N_obs <- N_obs_old + n_blk
    N_obs_old <- N_obs
  }
  print('data generation')
  
  t_start <- Sys.time()
  if(delta_p > 0){
    
    n0 <- n_blk
    print('delta_p>0')
    N_eff <- n_blk
    W <- n_blk/N_eff
    I <- (p_obs_old+1):p_obs
    idx_abs <- c(idx_abs, I)
    X <- X[,idx_abs]
    Cy <- rep(0, K_cv)
    CXX <- array(0,dim=c(length(idx_abs), length(idx_abs), K_cv))
    CXy <- matrix(0, length(idx_abs), K_cv)
    
    # update statistics
    {
      idx_imp <- 1:(length(idx_abs))
      Obs <- cbind(matrix(X,n_blk),y)
      for(k_cv in 1:K_cv){
        sub_obs <- (n_blk/K_cv*(k_cv-1)+1):(n_blk/K_cv*k_cv)
        Cy[k_cv] <- 1/N_eff*sum((y[sub_obs])^2)
        CXX[,,k_cv] <- 1/N_eff*t(X[sub_obs,])%*%X[sub_obs,]
        CXy[,k_cv] <- 1/N_eff*t(X[sub_obs,])%*%y[sub_obs]
      }
    }
    
    p_obs_old <- p_obs
    
  }else{
    
    N_eff <- N_eff + n_blk
    X <- X[,c(idx_abs)]
    
    
    if(N_eff<=n0){
      Obs <- rbind(Obs, cbind(matrix(X,n_blk),y))
    }
    
    X <- matrix(X, nrow = n_blk)
    for(k_cv in 1:K_cv){
      sub_obs <- (n_blk/K_cv*(k_cv-1)+1):(n_blk/K_cv*k_cv)
      Cy[k_cv] <- (N_eff-n_blk)/N_eff * Cy[k_cv] + n_blk/N_eff * mean((y[sub_obs])^2)
      CXX[,,k_cv] <- (N_eff-n_blk)/N_eff * CXX[,,k_cv] + 1/N_eff * t(X[sub_obs,])%*%X[sub_obs,]
      CXy[,k_cv] <- (N_eff-n_blk)/N_eff * CXy[,k_cv] + 1/N_eff * t(X[sub_obs,])%*%y[sub_obs]
    }
    
  }
  
  
  # initial estimate
  if(N_eff<=n0){
    if(N_eff==n0){
      set.seed(123) 
      cv_fit <- cv.glmnet(Obs[,-ncol(Obs)], Obs[,ncol(Obs)], nfolds = 5)  
      lam <- cv_fit$lambda.min
      modelfit <- glmnet(Obs[,-ncol(Obs)], Obs[,ncol(Obs)], lambda = lam)
      beta <- as.vector(modelfit$beta)
      sigma_hat <- min(sigma_hat,
                       sqrt(mean((Obs[,ncol(Obs)]-Obs[,-ncol(Obs)]%*%beta)^2)))
      tmpt_p <- length(intersect(I,idx_abs))
      C_lam <- lam/(sigma_hat * sqrt(log(tmpt_p) / (N_eff / 2)))
      C_lam_set <- seq(1, 2.8, 0.2)*C_lam
      print('initial estimation')
    }
    next
  }
  
  pt <- length(idx_imp)
  print(paste0('dim=',pt, ' sigma_hat=', round(sigma_hat,3)))
  
  if(N_eff <= length(I)){
    
    Cy_total <- sum(Cy); CXy_total <- rowSums(CXy); CXX_total <- apply(CXX, c(1,2), sum)
    tmpt_p <- pt # length(intersect(I,idx_abs))
    if(tmpt_p>1){
      delta <- min(1, log(N_eff)/log(tmpt_p))
      ## estimate noise level
      if(abs(sigma_hat-sigma_hat0)>0.01){
        sigma_hat0 <- sigma_hat
        sigma_hat <- scale(Cy_total, CXy_total, CXX_total, sigma_hat, beta, N_eff)
      }
      ## cross validation
      if(cv_num < 10){
        print('-------------- CV!!')
        err <- c()
        
        for(C_lam in C_lam_set){
          
          lam <- C_lam * sigma_hat * sqrt(log(tmpt_p)^(delta) / (N_eff / 2))            
          lam_vec <- rep(lam, tmpt_p)
          err1 <- numeric(K)
          
          for(k_cv in 1:K_cv){
            
            CXy_1 <- CXy_total - CXy[,k_cv]
            CXX_1 <- CXX_total - CXX[,,k_cv]
            beta1 <- proximal(CXy_1, CXX_1, lam_vec, beta)
            term1 <- Cy[k_cv]
            term2 <- sum(beta1 * (CXX[, , k_cv] %*% beta1))
            term3 <- 2 * crossprod(CXy[, k_cv], beta1)
            err1[k_cv] <- term1 + term2 - term3
            
          }
          
          err <- c(err, sum(err1))
          # print(C_lam)
        }
        C_lam_opt <- C_lam_set[which.min(err)]
        if(C_lam_opt==C_lam_opt0){cv_num <- cv_num+1}else{cv_num <- 0}
        C_lam_opt0 <- C_lam_opt
      }
    }
    if(tmpt_p==1){lam <-0}
    
    lam <- C_lam_opt*sigma_hat*sqrt(log(tmpt_p)^(delta)/(N_eff))
    lam_vec <- rep(lam, tmpt_p)
    beta <- proximal(CXy_total, CXX_total, lam_vec, beta)
    
  }else{
    print('-------------- hard threshold!!')
    Cy_total <- sum(Cy); CXy_total <- rowSums(CXy); CXX_total <- apply(CXX, c(1,2), sum)
    gram <- ginv(CXX_total)
    beta <- gram%*%matrix(CXy_total, pt, 1)
    sigma_hat <- sqrt(Cy_total - 2*matrix(beta,1)%*%CXy_total + matrix(beta,1)%*%CXX_total%*%beta)
    sigma_hat <- min(sigma_hat,10)
    beta1 <- rep(0.1/sqrt(n0)*as.vector(sigma_hat),pt)# C_lam_opt/sqrt(N_eff)*diag(gram)^(1/2)*as.vector(sigma_hat)
    beta <- beta * sapply(1:pt, function(j){abs(beta[j])>abs(beta1[j])})
  }
  
  t_end <- Sys.time()
  # compress statistics according to idx_imp
  {
    err_y <- Cy_total - 2*CXy_total%*%beta + matrix(beta,nrow=1)%*%CXX_total%*%beta
    err_y_o <- c(err_y_o, err_y)
    idx_imp <- unique(which(beta!=0))
    beta <- beta[idx_imp]
    idx_abs <- idx_abs[idx_imp]
    beta_expan <- rep(0, p_obs)
    beta_expan[idx_abs] <- beta
    CXX <- CXX[idx_imp, idx_imp,]
    CXy <- CXy[idx_imp,]
    pt_o <- c(pt_o, pt)
    time_o <- c(time_o, as.numeric(t_end-t_start))
    beta_o <-c(beta_o,list(beta_expan))
    R2 <- 1-err_y
    R2_o <- c(R2_o, R2)
    sigma_o <- c(sigma_o,as.vector(sigma_hat))
    k_o <- c(k_o, k)
  }
  
  if(k%%100==0){
    save.image(paste0('k=',k,'.Rdata'))
    file.remove(paste0('k=',k-100,'.Rdata'))
  }
  if(k==K_train){save.image(paste0('k=',k,'.Rdata'))}
}
save(pt_o,time_o,beta_o,R2_o,sigma_o,k_o,
     file = 'P2P_train_results.Rdata')



######### test ############
K <- 2261
y_test <- c(); X_test <- c()
for(k in (K_train+1):K)
{
  print(k)
  file_cur <- paste0('P2Pdata_',k,'.csv')
  data_cur <- read.csv(file_cur)
  data_cur <- data_cur[,-c(1,2,3)]
  y <- data_cur$int_rate
  # covariate
  data_cur <- subset(data_cur, select = -c(int_rate))
  X <- as.matrix(data_cur)
  p0 <- ncol(X)
  idx_del <- which(sapply(1:nrow(X), function(i){sum(is.na(X[i,]))})>0)
  if(length(idx_del)>0){
    X <- X[-idx_del, ]
    y <- y[-idx_del]
  }
  p_interact0 <- which(names(data_cur)=='home_ownership_MORTGAGE')
  j1_indices <- 1:(p_interact0-1)
  j2_indices <- p_interact0:p0
  combinations <- expand.grid(j1 = j1_indices, j2 = j2_indices)
  interaction_terms <- X[, combinations$j1] * X[, combinations$j2]
  X <- cbind(X, interaction_terms)    
  y_test <- c(y_test, y)
  X_test <- rbind(X_test, X)
  rm(data_cur)
}
# save(X_test, y_test, file = 'test_data.Rdata')
err_test_o <- c()
for(i in 1:length(beta_o)){
  if(i%%10==0){print(i)}
  beta_expan <- rep(0,ncol(X_test))
  beta_expan[1:length(beta_o[[i]])] <- beta_o[[i]]
  err_test_o <- c(err_test_o, mean((y_test-X_test%*%beta_expan)^2))
}

save(err_test_o, y_test, file = paste0('P2P_test_results.Rdata'))
