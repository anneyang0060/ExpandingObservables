# rm(list=ls())
# # setwd('/home/yangy/OnlineHD')
# setwd('C:/Users/annie/OneDrive/ImportantFiles/projects/1. InProgress/Online HD/code')
# library(mvnfast)
# library(glmnet)
# library(MASS)
# source('FNS.R')
# rho <- 0.8
# p <- 1000
# s <- 10
# Tmax <- 10000
# n_blk <- 10
# sigma_eps<- 1
# n0 <- 500
# sd <- 90

############################# dynamic screening with lasso ##########
Mcl <- 50
cl <- makeCluster(Mcl)
registerDoParallel(cl)
res_c <- foreach(sd=sds, .packages = c('MASS' ,'mvnfast','glmnet')) %dopar%
  {
    Cy <- 0; CX <- 0; CXX <- 0; CXy <- 0
    a1 <- 0; a2 <- 0; A1 <- 0
    Obs <- c(); lam_set <- c()
    sigma_hat <- 1; p_old <- 0; K <- 0
    err_beta_o <- c(); err_y_o <- c()
    time_o_bg <- c(); time_o_u <- c()
    err_beta_b_o <- c(); err_beta_g_o <- c()
    dt_o_b <- c(); miss_b <- c()
    dt_o_g <- c(); miss_g <- c()
    
    set.seed(sd)
    
    for(t in seq(n_blk,Tmax,n_blk))
    {
      gam <- 0.1
      ### generate data
      {
        beta_true <- rep(0,p) 
        beta_true[1:10] <- rep(c(3,-3),s/2)
        sigma_X <-matrix(0, p,p)
        for(k in 1:p){
          sigma_X[k,] <- sapply(1:p, function(l){rho^(abs(k-l))})
        }
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
        # select lambda
        # if(t>n0){
        #   beta0 <- beta_b
        #   if(rho<0.7){
        #     lam_set <-c(0.1, 0.15, 0.20, 0.25, 0.30)
        #   }else{
        #     lam_set <- c(0.05, 0.075, 0.1,0.125, 0.15,0.175, 0.20)
        #   }
        #   pre_err <- c()
        #   for(lam1 in lam_set){
        #     beta1 <- proximal(CXy, CXX, rep(lam1,p), beta0)
        #     pre_err <- c(pre_err,sum((y-X%*%beta1)^2))
        #   }
        #   lam <- lam_set[which.min(pre_err)]
        # }
        # lam <- sqrt(log(p)/t)
        CXmax0 <- max(abs(CXX))
        Cy <- (t-n_blk)/t * Cy + n_blk/t * mean(y^2)
        CX <- (t-n_blk)/t * CX + n_blk/t * colMeans(X)
        CXX <- (t-n_blk)/t * CXX + 1/t * t(X)%*%X
        CXy <- (t-n_blk)/t * CXy + 1/t * t(X)%*%y
        CXmax1 <- max(abs(CXX))
      }
      
      # initial estimate
      if(t <= n0){
        if(t ==n0){
          cvfit <- cv.glmnet(Obs[,-ncol(Obs)], Obs[,ncol(Obs)])
          lam <- cvfit$lambda.min
          lam_g <- sqrt(log(p)/n0)
          modelfit <- glmnet(Obs[,-ncol(Obs)], Obs[,ncol(Obs)], lambda = lam, gamma = gam)
          beta_b <- as.vector(modelfit$beta)
          beta_g <- beta_b
          gamma <- matrix(0,p-1,p)
          for(j in 1:p){
            modelfit <- glmnet(Obs[,-c(j,ncol(Obs))], Obs[,j], lambda = lam, gamma = gam)
            gamma[,j] <- as.vector(modelfit$beta)
          }
          a1 <- a1 + sapply(1:p, function(j){t(Obs[,j]-Obs[,-c(j,ncol(Obs))]%*%gamma[,j])%*%Obs[,j]})
          a2 <- a2 + sapply(1:p, function(j){t(Obs[,j]-Obs[,-c(j,ncol(Obs))]%*%gamma[,j])%*%Obs[,ncol(Obs)]})
          A1 <- A1 + sapply(1:p, function(j){t(Obs[,j]-Obs[,-c(j,ncol(Obs))]%*%gamma[,j])%*%(Obs[,1:p])%*%beta_b})
          rm(cvfit, modelfit)
        }
        next
      }
      
      # online estimate
      {
        # OR estimate
        lam_g <- sqrt((t-n_blk)/t*CXmax1/CXmax0)*lam_g
        beta_g1 <- proximal(CXy, CXX, rep(lam_g,p), beta_g, gamma = gam)
        if(sum((beta_g1-beta_g)^2)<=10){
          beta_g <- beta_g1
        }else{
          beta_g <- proximal(CXy, CXX, rep(lam_g,p), beta_g, gamma = 0.005)
        }
        t_end1 <- Sys.time()
        time_o_bg <- c(time_o_bg, as.numeric(t_end1-t_start))
        # ODL 
        lams <- rep(sqrt(log(p)/t), p)
        beta_b1 <- proximal(CXy, CXX, lams, beta_b, gamma = gam)
        if(sum((beta_b1-beta_b)^2)<=10){
          beta_b <- beta_b1
          for(j in 1:p){
            gamma[,j] <- proximal(CXX[-j,j], CXX[-j,-j], lams,gamma[,j], gamma = gam)
          }
          a1 <- a1 + sapply(1:p, function(j){t(X[,j]-X[,-j]%*%gamma[,j])%*%X[,j]})
          a2 <- a2 + sapply(1:p, function(j){t(X[,j]-X[,-j]%*%gamma[,j])%*%y})
          A1 <- A1 + sapply(1:p, function(j){t(X[,j]-X[,-j]%*%gamma[,j])%*%(X)%*%beta_b})
          beta <- beta_b+a1^(-1)*(a2-as.vector(A1))/n_blk
        }else{
          beta_b <- proximal(CXy, CXX, lams, beta_b, gamma = 0.005)
          for(j in 1:p){
            gamma[,j] <- proximal(CXX[-j,j], CXX[-j,-j], lams,gamma[,j])
          }
          a1 <- a1 + sapply(1:p, function(j){t(X[,j]-X[,-j]%*%gamma[,j])%*%X[,j]})
          a2 <- a2 + sapply(1:p, function(j){t(X[,j]-X[,-j]%*%gamma[,j])%*%y})
          A1 <- A1 + sapply(1:p, function(j){t(X[,j]-X[,-j]%*%gamma[,j])%*%(X)%*%beta_b})
          beta <- beta_b+a1^(-1)*(a2-as.vector(A1))/n_blk
        }
        t_end2 <- Sys.time()
        time_o_u <- c(time_o_u, as.numeric(t_end2-t_end1))
      }
      
      # store expanded beta
      {
        err_beta_b <- sum((beta_true-beta_b)^2)
        err_beta_g <- sum((beta_true-beta_g)^2)
        err_beta <- sum((beta_true-beta)^2)
        err_beta_o <- c(err_beta_o, err_beta)
        err_beta_b_o <- c(err_beta_b_o, err_beta_b)
        err_beta_g_o <- c(err_beta_g_o, err_beta_g)
        idx_true <- 1:s
        miss1 <- 0
        idx_abs <- which(beta_b!=0)
        for(i in idx_true){if(i%in%idx_abs==0){miss1<-miss1+1}}
        miss_b <- c(miss_b,miss1)
        dt_o_b <- c(dt_o_b, length(idx_abs))
        miss1 <- 0
        idx_abs <- which(beta_g!=0)
        for(i in idx_true){if(i%in%idx_abs==0){miss1<-miss1+1}}
        miss_g <- c(miss_g,miss1)
        dt_o_g <- c(dt_o_g, length(idx_abs))
        print(paste0('t=',t, ' err_g = ',round(err_beta_g,4), ' err_beta=',round(err_beta,4),
                     ' err_beta_b=',round(err_beta_b,4),' lam_g=', round(lam_g,3)))
      }
      # lam_o <- c(lam_o, list(lam_set))
    }
    if(!dir.exists(prefix)){dir.create(prefix)}
    save(err_beta_o,err_beta_b_o, time_o_bg, err_beta_g_o,time_o_u,
         dt_o_b,miss_b,dt_o_g,miss_g,
         file = paste0(prefix,'/','sd',sd,'.Rdata'))
    
    
  }
stopCluster(cl)
