# rm(list=ls())
# # setwd('/home/yangy/OnlineHD')
# setwd('C:/Users/annie/OneDrive/ImportantFiles/projects/1. InProgress/Online HD/code')
# library(mvnfast)
# library(glmnet)
# library(MASS)
# source('FNS.R')
# rho <- 0.5
# p <- 1000
# s <- 10
# Tmax <- 10000
# n_blk <- 100
# sigma_eps<- 1
# n0 <- 200

############################# dynamic screening with lasso ##########
Mcl <- 50
cl <- makeCluster(Mcl)
registerDoParallel(cl)
res_c <- foreach(sd=sds, .packages = c('MASS' ,'mvnfast','glmnet')) %dopar%
  {
    # sd <- 2
    Cy <- 0; CX <- 0; CXX <- 0; CXy <- 0
    Obs <- c(); sigma_hat <- 1; p_old <- 0
    err_beta_o <- c(); time_o <- c()
    dt_o <- c(); miss <- c()
    set.seed(sd)
    
    
    for(t in seq(n_blk,Tmax,n_blk))
    {
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
      
      # initial estimate
      {
        if(t <= n0){
          Obs <- rbind(Obs, cbind(matrix(X,n_blk),y))
          if(t ==n0){
            cvfit <- cv.glmnet(Obs[,-ncol(Obs)], Obs[,ncol(Obs)])
            lam <- cvfit$lambda.min
            modelfit <- glmnet(Obs[,-ncol(Obs)], Obs[,ncol(Obs)], lambda = lam, gamma = gam)
            beta0 <- as.vector(modelfit$beta)
            rm(cvfit, modelfit)
          }
          next
        }
      }
      
      # update statistics
      if(t>n0){
        gam <-  min((t-n0)/1000,1)*0.1
        Cy <- (t-n_blk-n0)/(t-n0) * Cy + n_blk/(t-n0) * mean(y^2)
        CX <- (t-n_blk-n0)/(t-n0) * CX + n_blk/(t-n0) * colMeans(X)
        CXX <- (t-n_blk-n0)/(t-n0) * CXX + 1/(t-n0) * t(X)%*%X
        CXy <- (t-n_blk-n0)/(t-n0) * CXy + 1/(t-n0) * t(X)%*%y
      }
      
      # estimate
      {
        lams <- rep(sqrt(log(p)/(t-n0)), p)
        beta1 <- proximal_OLL(CXy, CXX, Obs, lams, beta0, gamma = gam)
        if(sum((beta1-beta0)^2)<=10){
          beta <- beta1
        }else{ # id not converge, try a smaller step-size
          beta <- proximal_OLL(CXy, CXX, Obs, lams, beta0, gamma = 0.005)
        }
      }
      
      # store expanded beta
      {
        err_beta <- sum((beta_true-beta)^2)
        beta0 <- beta
        idx_true <- 1:s
        miss1 <- 0
        idx_abs <- which(beta!=0)
        for(i in idx_true){if(i%in%idx_abs==0){miss1<-miss1+1}}
        miss <- c(miss,miss1)
        dt_o <- c(dt_o, length(idx_abs))
        t_end <- Sys.time()
        err_beta_o <- c(err_beta_o, err_beta)
        time_o <- c(time_o, as.numeric(t_end-t_start))
        print(paste0('t=',t, ' err_beta=',round(err_beta,4)))
      }
      # lam_o <- c(lam_o, list(lam_set))
    }
    if(!dir.exists(prefix)){dir.create(prefix)}
    save(err_beta_o,time_o,dt_o,miss,
         file = paste0(prefix,'/','sd',sd,'_miss.Rdata'))
    
  }
stopCluster(cl)
