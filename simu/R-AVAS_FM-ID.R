############################# dynamic screening with lasso ##########
Mcl <- 50
cl <- makeCluster(Mcl)
registerDoParallel(cl)
res_c <- foreach(sd=sds, .packages = c('MASS' ,'mvnfast','glmnet')) %dopar%
  {
    s <- 1; I_set <- list(); W_set <- c()
    Obs <- c(); lam_set <- c()
    sigma_hat <- 1; sigma_hat0<-0; p_old <- 0; K <- 0
    st <- 0; pt <- 0
    idx_imp <- c(); idx_abs <- c()
    err_beta_o <- c(); err_psue_o <- c()
    time_o <- c(); pt_o <- c(); est_time <- c(); miss <- c()
    C_lam_set <- 1:10; C_lam_opt0 <- 0; cv_num <- 0
    sigma_hat_o <- c()
    
    set.seed(sd)
    
    for(t in seq(n_blk,Tmax,n_blk))
    {
      ### generate data
      {
        k_tau <- max(which(tau<t))
        t_eff <- t-tau[k_tau]
        p_new <- p[k_tau]
        if(p_new>p_old){idx_abs <- c(idx_abs, (p_old+1):p_new)}
        if(t%in%(tau+n_blk)){
          p_obs <- length(idx_imp)+delta_p[k_tau]
        }else{
          p_obs <- length(idx_imp)
        }
        
        alpha <- beta_true[idx_abs]
        alpha_psue <- beta_psue[k_tau, idx_abs]
        gamma <- (beta_true[(p_new+1):length(beta_true)])[which(beta_true[(p_new+1):length(beta_true)]!=0)]
        Phi <- diag(0.4, p_obs)
        p_u <- sum(beta_true!=0)-sum(alpha!=0) 
        sigma_X <- matrix(0,p_obs,p_obs)
        for(k_dep in 0:(max(idx_abs)%/%m_dep)){
          idx_dep1 <- which(sapply(idx_abs, function(i){i%in%(idx_dep+k_dep*m_dep)}) == TRUE)
          sigma_X[idx_dep1, idx_dep1] <- rho
        }
        diag(sigma_X) <- 1
        if(p_u>1){
          sigma_U<-matrix(rho, p_u,p_u)
          diag(sigma_U) <- 1
          sigma_XU <- matrix(0, p_obs, p_u)
          idx_dep1 <- which(sapply(idx_abs, function(i){i%in%(idx_dep)}) == TRUE)
          sigma_XU[idx_dep1, ] <- rho
          Sigma <- rbind(cbind(sigma_X, sigma_XU), cbind(t(sigma_XU), sigma_U))
          XU <- rmvn(n_blk, rep(0,p_u+p_obs), Sigma)
          X <- XU[,1:p_obs]
          U <- XU[,(p_obs+1):(p_obs+p_u)]
          eps <- rnorm(n_blk,0,sigma_eps)
          y <- X %*% alpha + U %*% gamma + eps
        }else{
          X <- rmvn(n_blk, rep(0,p_obs), sigma_X)
          eps <- rnorm(n_blk,0,sigma_eps)
          y <- X %*% alpha + eps
          
        }
        
        
        t_start <- Sys.time()
      }
      
      if(p_new > p_old){
        
        Cy <- rep(0, K_cv); CXX <- array(0,dim=c(p_obs, p_obs, K_cv))
        CXy <- matrix(0, p_obs, K_cv)
        n_new <- n_blk
        W <- n_blk/t_eff
        I <- (p_old+1):p_new
        
        # update statistics
        {
          idx_imp <- 1:p_obs
          Obs <- cbind(matrix(X,n_blk),y)
          for(k in 1:K_cv){
            sub_obs <- (n_blk/K_cv*(k-1)+1):(n_blk/K_cv*k)
            Cy[k] <- 1/t_eff*sum((y[sub_obs])^2)
            CXX[,,k] <- 1/t_eff*t(X[sub_obs,])%*%X[sub_obs,]
            CXy[,k] <- 1/t_eff*t(X[sub_obs,])%*%y[sub_obs]
          }
        }
        
        p_old <- p_new
        cv_num <- 0
        
      }else{
        
        n_new <- n_new + n_blk
        
        if(n_new <= n0){
          Obs <- rbind(Obs, cbind(matrix(X,n_blk),y))
        }
        
        X <- matrix(X,n_blk,length(idx_abs))
        for(k in 1:K_cv){
          sub_obs <- (n_blk/K_cv*(k-1)+1):(n_blk/K_cv*k)
          Cy[k] <- (t_eff-n_blk)/t_eff * Cy[k] + n_blk/t_eff * mean((y[sub_obs])^2)
          CXX[,,k] <- (t_eff-n_blk)/t_eff * CXX[,,k] + 1/t_eff * t(X[sub_obs,])%*%X[sub_obs,]
          CXy[,k] <- (t_eff-n_blk)/t_eff * CXy[,k] + 1/t_eff * t(X[sub_obs,])%*%y[sub_obs]
        }
        
      }
      
      # initial estimate
      if(n_new <= n0){
        if(n_new==n0){
          cvfit <- cv.glmnet(Obs[,-ncol(Obs)], Obs[,ncol(Obs)])
          modelfit <- glmnet(Obs[,-ncol(Obs)], Obs[,ncol(Obs)], lambda = cvfit$lambda.min)
          beta <- as.vector(modelfit$beta)
          rm(cvfit, modelfit)
        }
        next
      }
      
      pt <- length(idx_imp)
      
      if(t_eff <= 2*length(I)){
        
        Cy_total <- sum(Cy); CXy_total <- rowSums(CXy); CXX_total <- apply(CXX, c(1,2), sum)
        tmpt_p <- pt; 
        if(tmpt_p>1){
          delta <- min(1, log(t_eff)/log(tmpt_p))
          ## estimate noise level
          if(abs(sigma_hat-sigma_hat0)>0.01){
            sigma_hat0 <- sigma_hat
            sigma_hat <- scale(Cy_total, CXy_total, CXX_total, sigma_hat, beta, t_eff)
            if(sigma_hat>10){sigma_hat <- 10; sigma_hat0 <- 0}
          }
          ## cross validation
          if(cv_num < 10){
            # print('-------------- CV!!')
            err <- c()
            
            for(C_lam in C_lam_set){
              
              lam <- C_lam * sigma_hat * sqrt(log(tmpt_p)^(delta) / (t_eff / 2))            
              lam_vec <- rep(lam, tmpt_p)
              err1 <- numeric(K)
              
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
              # print(C_lam)
            }
            C_lam_opt <- C_lam_set[which.min(err)]
            if(C_lam_opt==C_lam_opt0){cv_num <- cv_num+1}else{cv_num <- 0}
            C_lam_opt0 <- C_lam_opt
          }
        }
        if(tmpt_p==1){lam <-0}
        
        lam <- C_lam_opt*sigma_hat*sqrt(log(tmpt_p)^(delta)/(t_eff))
        lam_vec <- rep(lam, tmpt_p)
        beta <- proximal(CXy_total, CXX_total, lam_vec, beta)
        
      }else{
        # print('-------------- hard threshold!!')
        Cy_total <- sum(Cy); CXy_total <- rowSums(CXy); CXX_total <- apply(CXX, c(1,2), sum)
        gram <- ginv(CXX_total)
        beta <- gram%*%matrix(CXy_total, pt, 1)
        sigma_hat <- sqrt(Cy_total - 2*matrix(beta,1)%*%CXy_total + matrix(beta,1)%*%CXX_total%*%beta)
        sigma_hat <- min(sigma_hat,10)
        beta1 <- rep(0.5/sqrt(n0)*as.vector(sigma_hat), pt) # 1/sqrt(t_eff)*diag(gram)^(1/2)*as.vector(sigma_hat)
        beta <- beta * sapply(1:pt, function(j){abs(beta[j])>abs(beta1[j])})
      }
      
      # store expanded beta
      {
        err_beta <- sum((beta-alpha)^2)
        err_psue <- sum((beta-alpha_psue)^2)
        err_y <- Cy_total - 2*matrix(beta,1)%*%CXy_total + matrix(beta,1)%*%CXX_total%*%beta
        idx_imp <- unique(which(beta!=0))
        beta <- beta[idx_imp]
        idx_abs <- idx_abs[idx_imp]
        beta_expan <- rep(0, p_new)
        beta_expan[idx_abs] <- beta
        miss1 <- 0
        idx_true <- which(beta_true[1:p_new]!=0)
        for(i in idx_true){if(i%in%idx_abs==0){miss1<-miss1+1}}
        
        # compress statistics according to idx_imp
        CXX <- CXX[idx_imp, idx_imp,]
        CXy <- CXy[idx_imp,]
        t_end <- Sys.time()
        est_time <- c(est_time,t)
        err_beta_o <- c(err_beta_o, err_beta)
        err_psue_o <- c(err_psue_o, err_psue)
        pt_o <- c(pt_o, pt)
        miss <- c(miss,miss1)
        sigma_hat_o <- c(sigma_hat_o, sigma_hat)
        
        time_o <- c(time_o, as.numeric(t_end-t_start))
        # print(paste0('t=',t,' err_beta=',round(err_beta,3),
        #              ' err_y=',round(err_y,3),' pt=',pt, ' miss=',miss1))

      }
    }
    if(!dir.exists(prefix)){dir.create(prefix)}
    save(err_beta_o, err_psue_o, time_o, pt_o, est_time, miss, # idx_abs_ls,
         file = paste0(prefix,'/','sd',sd,'.Rdata'))
    
  }
stopCluster(cl)

