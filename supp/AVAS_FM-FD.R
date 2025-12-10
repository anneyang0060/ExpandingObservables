rm(list=ls())
# setwd('/home/yangy/OnlineHD')
setwd('C:/Users/annie/OneDrive/ImportantFiles/projects/1. InProgress/Online HD/fixed model/code')
library(mvnfast)
library(glmnet)
library(MASS)
source('FNS.R')
rho <- 0
p <- 1000
s <- 10
Tmax <- 10000
n_blk <- 10
sigma_eps<- 1
n0 <- 200
beta_true <- rep(0,p)
beta_true[1:10] <- rep(c(3,-3),s/2)
sigma_X <-matrix(0, p,p)
for(k in 1:p){
  sigma_X[k,] <- sapply(1:p, function(l){rho^(abs(k-l))})
}
K_cv <- 2
sd <- 1

############################# dynamic screening with lasso ##########
Mcl <- 50
cl <- makeCluster(Mcl)
registerDoParallel(cl)
res_c <- foreach(sd=sds, .packages = c('MASS' ,'mvnfast','glmnet')) %dopar%
{
  Cy <- rep(0, K_cv); CXX <- array(0,dim=c(p,p,K_cv)); CXy <- matrix(0,p,K_cv)
  Obs <- c(); C_lam_set <- 1:10
  sigma_hat <- 1; sigma_hat0 <- 0; p_old <- 0
  st <- 0; pt <- 0; sigma_o <- c()
  idx_imp <- c(); idx_abs <- 1:p
  err_beta_o <- c(); err_y_o <- c()
  time_o <- c(); pt_o <- c(); est_time <- c(); miss <- c()
  C_lam_opt0 <- 0; cv_num <- 0
  tmpt_p <- p
  
  set.seed(sd)
  
  for(t in seq(n_blk,Tmax,n_blk))
  {
    ### generate data
    {
      X <- rmvn(n_blk, rep(0,p), sigma_X)
      eps <- rnorm(n_blk,0,sigma_eps)
      y <- X %*% beta_true + eps
      t_start <- Sys.time()
    }
    
    # update statistics
    {
      if(t <= n0){
        Obs <- rbind(Obs, cbind(matrix(X,n_blk),y))
      }
      X <- X[,idx_abs]
      for(k in 1:K_cv){
        sub_obs <- (n_blk/K_cv*(k-1)+1):(n_blk/K_cv*k)
        Cy[k] <- ((t-n_blk)/t * Cy[k] + 1/t * sum(y[sub_obs]^2))
        CXX[,,k] <- (t-n_blk)/t * CXX[,,k] +1/t * t(X[sub_obs,])%*%X[sub_obs,]
        CXy[,k] <- (t-n_blk)/t * CXy[,k] + 1/t * t(X[sub_obs,])%*%y[sub_obs]  
      }
      
    }
    
    # initial estimate
    if(t <= n0){
      if(t ==n0){
        cvfit <- cv.glmnet(Obs[,-ncol(Obs)], Obs[,ncol(Obs)])
        modelfit <- glmnet(Obs[,-ncol(Obs)], Obs[,ncol(Obs)], lambda = cvfit$lambda.min)
        beta <- as.vector(modelfit$beta)
        beta_opt <- beta
        rm(cvfit, modelfit, Obs)
      }
      next
    }
    
    pt <- length(idx_abs)
    
    if(t <= 2*tmpt_p){
      
      Cy_total <- sum(Cy); CXy_total <- rowSums(CXy); CXX_total <- apply(CXX, c(1,2), sum)
      tmpt_p <- length(idx_abs)
      if(tmpt_p>1){
        delta <- min(1, log(t)/log(tmpt_p))
        ## estimate noise level
        if(abs(sigma_hat-sigma_hat0)>0.001){
          sigma_hat0 <- sigma_hat
          sigma_hat <- scale(Cy_total, CXy_total, CXX_total, sigma_hat, beta, t)
        }
        ## cross validation
        if(cv_num < 10){
          print('-------------- CV!!')
          err <- c()
          
          for(C_lam in C_lam_set){
            
            lam <- C_lam * sigma_hat * sqrt(log(tmpt_p)^(delta) / (t / 2))            
            lam_vec <- rep(lam, tmpt_p)
            err1 <- numeric(K_cv)
            
            for(k in 1:K_cv){
              
              CXy_1 <- CXy_total - CXy[,k]
              CXX_1 <- CXX_total - CXX[,,k]
              beta1 <- proximal(CXy_1, CXX_1, lam_vec, beta)
              term1 <- Cy[k]
              term2 <- sum(beta1 * (CXX[, , k] %*% beta1))
              term3 <- 2 * crossprod(CXy[, k], beta1)
              err1[k] <- term1 + term2 - term3
              
            }
            
            err <- c(err, sum(err1))
            print(C_lam)
          }
          C_lam_opt <- C_lam_set[which.min(err)]
          if(C_lam_opt==C_lam_opt0){cv_num <- cv_num+1}else{cv_num <- 0}
          C_lam_opt0 <- C_lam_opt
        }
      }
      lam <- C_lam_opt*sigma_hat*sqrt(log(tmpt_p)^(delta)/(t))
      lam_vec <- rep(lam, tmpt_p)
      beta <- proximal(CXy_total, CXX_total, lam_vec, beta)
      if(delta<1){
        lam <-  C_lam_opt*sigma_hat*sqrt(log(tmpt_p)/(t))
        lam_vec <- rep(lam, tmpt_p)
        beta_opt <- proximal(CXy_total, CXX_total, lam_vec, beta)
      }else{
        beta_opt <- beta
      }
      
    }else{
      print('-------------- hard threshold!!')
      Cy_total <- sum(Cy); CXy_total <- rowSums(CXy); CXX_total <- apply(CXX, c(1,2), sum)
      gram <- solve(CXX_total)
      beta <- gram%*%matrix(CXy_total, pt, 1)
      sigma_hat <- sqrt(Cy_total - 2*matrix(beta,1)%*%CXy_total + matrix(beta,1)%*%CXX_total%*%beta)
      beta1 <- 3/sqrt(t)*diag(gram)^(1/2)*as.vector(sigma_hat)
      beta_opt <- beta * sapply(1:pt, function(j){abs(beta[j])>abs(beta1[j])})
      beta <- beta_opt
    }
    
    # store expanded beta
    {
      err_beta <- sum((beta_true[idx_abs]-beta_opt)^2)
      idx_imp <- unique(which(beta!=0))
      beta <- beta[idx_imp]
      idx_abs <- idx_abs[idx_imp]
      idx_true <- 1:s
      miss1 <- 0
      for(i in idx_true){if(i%in%idx_abs==0){miss1<-miss1+1}}
      
      # compress statistics according to idx_imp
      CXX <- CXX[idx_imp, idx_imp,]
      CXy <- CXy[idx_imp,]
      t_end <- Sys.time()
      err_beta_o <- c(err_beta_o, err_beta)
      pt_o <- c(pt_o, pt)
      miss <- c(miss,miss1)
      sigma_o <- c(sigma_o,sigma_hat)
      time_o <- c(time_o, as.numeric(t_end-t_start))
      print(paste0('t=',t,' err_beta=',round(err_beta,3),
                   ' sigma_hat', round(sigma_hat,3),' pt=',pt, ' miss=',miss1))
    }
    # lam_o <- c(lam_o, list(lam_set))
  }
  if(!dir.exists(prefix)){dir.create(prefix)}
  save(err_beta_o, time_o, pt_o,  miss, sigma_o,
       file = paste0(prefix,'/','sd',sd,'.Rdata'))
  
}
stopCluster(cl)
