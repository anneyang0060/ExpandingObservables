# rm(list=ls())
# setwd('/home/yangy/OnlineHD')
# # setwd('C:/Users/annie/OneDrive/ImportantFiles/projects/1. InProgress/Online HD/code')
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
# beta_true <- rep(0,p)
# beta_true[1:10] <- rep(c(3,-3),s/2)
# sigma_X <-matrix(0, p,p)
# for(k in 1:p){
#   sigma_X[k,] <- sapply(1:p, function(l){rho^(abs(k-l))})
# }
# sd <- 90

############################# dynamic screening with lasso ##########
Mcl <- 50
cl <- makeCluster(Mcl)
registerDoParallel(cl)
res_c <- foreach(sd=sds, .packages = c('MASS' ,'mvnfast','glmnet')) %dopar%
  {
    Obs <- c(); idx_abs <- 1:p
    err_beta_o <- c(); time_o <- c(); dt_o <- c();miss <- c()
    
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
      
      # initial esimate
      {
        if(t < n0){
          Obs <- rbind(Obs, cbind(matrix(X,n_blk),y))
          next
        }
        if(t == n0){
          cvfit <- cv.glmnet(Obs[,-ncol(Obs)], Obs[,ncol(Obs)])
          modelfit <- glmnet(Obs[,-ncol(Obs)], Obs[,ncol(Obs)], lambda = cvfit$lambda.min)
          beta <- as.vector(modelfit$beta)
          rm(cvfit, modelfit)
        }
      }
      
      # TSGD
      if(t>n0){
        # etas <- seq(0.0005,0.01,0.0005)
        lam <- sqrt(log(p)/t)
        if(rho<=0.3){
          eta <- 0.1/t^(0.5) 
        }else{
          eta <- 0.1/t^0.25
        }
        # eta <- 0.1/t^0.5 when rho=0
        # eta <- 0.1/t^0.5 when rho = 0.3
        # eta <- 0.1/t^0.25 when rho=0.5
        # eta <- 0.1/t^0.25 when rho=0.7
        beta <- truncation(beta + as.vector(eta * t(y-X%*%beta)%*%X/n_blk),
                           rep(eta*lam,p))
        # beta <- truncation(beta, rep(eta*lam,p))
        # beta <- beta + 2*as.vector(eta * t(y-X%*%beta)%*%X)
      }
      
      # store expanded beta
      {
        err_beta <- sum((beta_true-beta)^2)
        t_end <- Sys.time()
        err_beta_o <- c(err_beta_o, err_beta)
        idx_true <- 1:s
        miss1 <- 0
        idx_abs <- which(beta!=0)
        for(i in idx_true){if(i%in%idx_abs==0){miss1<-miss1+1}}
        miss <- c(miss,miss1)
        dt_o <- c(dt_o, length(idx_abs))
        time_o <- c(time_o, as.numeric(t_end-t_start))
        print(paste0('t=',t,' err_beta=',round(err_beta,3)))
      }
      # lam_o <- c(lam_o, list(lam_set))
    }
    if(!dir.exists(prefix)){dir.create(prefix)}
    save(err_beta_o, time_o, miss, dt_o,
         file = paste0(prefix,'/','sd',sd,'_miss.Rdata'))
    
  }
stopCluster(cl)