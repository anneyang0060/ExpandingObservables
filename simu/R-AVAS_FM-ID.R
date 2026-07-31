# rm(list=ls())
# setwd('')
# library(mvnfast)
# # library(scalreg)
# library(glmnet)
# library(MASS)
# source('FNS.R')
# rho <- 0.7
# source('simu/Setting_FM-ID.R')
# sd <- 1


############################# dynamic screening with lasso ##########
Mcl <- 50
cl <- makeCluster(Mcl)
registerDoParallel(cl)
res_c <- foreach(sd=sds, .packages = c('MASS' ,'mvnfast','glmnet')) %dopar%
  {
    s <- 1; I_set <- list(); W_set <- c()
    Obs <- c(); lam_set <- c()
    sigma_hat <- 1; sigma_hat0 <- 0; p_old <- 0; K <- 0
    st <- 0; pt <- 0
    idx_imp <- c(); idx_abs <- c()
    err_beta_o <- c(); err_psue_o <- c()
    time_o <- c(); pt_o <- c(); est_time <- c(); miss <- c()
    C_lam_set <- 1:10; 
    C_lam_opt0 <- 0; cv_num <- 0
    sigma_hat_o <- c()
    
    warm_done <- FALSE
    iota0_obs <- NA
    warm_grid <- NULL
    warm_grid_pos <- 1
    iota0_o <- c()
    M_upper_o <- c()
    warm_rhs_o <- c()
    
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
        
        # absolute indices of still-unobserved active variables
        idx_u_abs <- which(beta_true != 0 & seq_along(beta_true) > p_new)
        gamma <- beta_true[idx_u_abs]
        p_u <- length(idx_u_abs)
        
        # construct joint covariance for observed candidate variables and
        # unobserved active variables
        idx_joint <- c(idx_abs, idx_u_abs)
        
        Sigma_joint <- make_chain_crossblock_cov(
          idx_abs = idx_joint,
          rho = rho,
          block_width = block_width,
          m_dep = m_dep
        )
        
        p_obs <- length(idx_abs)
        
        sigma_X <- Sigma_joint[
          seq_len(p_obs),
          seq_len(p_obs),
          drop = FALSE
        ]
        
        eps <- rnorm(n_blk, 0, sigma_eps)
        
        if (p_u > 0) {
          
          XU <- rmvn(
            n_blk,
            rep(0, p_obs + p_u),
            Sigma_joint
          )
          
          X <- XU[, seq_len(p_obs), drop = FALSE]
          U <- XU[, p_obs + seq_len(p_u), drop = FALSE]
          
          y <- as.vector(X %*% alpha + U %*% gamma + eps)
          
        } else {
          
          X <- rmvn(
            n_blk,
            rep(0, p_obs),
            sigma_X
          )
          
          y <- as.vector(X %*% alpha + eps)
        }
        
        t_start <- Sys.time()
      }
      
      print(paste0('t=',t))
      
      if (p_new > p_old) {
        
        Cy <- rep(0, K_cv)
        CXX <- array(0, dim = c(p_obs, p_obs, K_cv))
        CXy <- matrix(0, p_obs, K_cv)
        
        n_new <- n_blk
        W <- n_blk / t_eff
        I <- (p_old + 1):p_new
        
        idx_imp <- 1:p_obs
        Obs <- cbind(matrix(X, n_blk), y)
        
        # initialize fold-specific sufficient statistics
        for (k in 1:K_cv) {
          sub_obs <- (n_blk / K_cv * (k - 1) + 1):(n_blk / K_cv * k)
          Cy[k] <- 1 / t_eff * sum((y[sub_obs])^2)
          CXX[, , k] <- 1 / t_eff * t(X[sub_obs, ]) %*% X[sub_obs, ]
          CXy[, k] <- 1 / t_eff * t(X[sub_obs, ]) %*% y[sub_obs]
        }
        
        # reset data-driven warm-up
        warm_done <- FALSE
        iota0_obs <- NA
        warm_grid_pos <- 1
        
        if (k_tau < length(tau)) {
          cycle_len_obs <- tau[k_tau + 1] - tau[k_tau]
        } else {
          cycle_len_obs <- Tmax - tau[k_tau]
        }
        max_warm_blocks <- max(1, floor(cycle_len_obs / n_blk))
        warm_grid <- make_warm_grid(
          min_blocks = warm_min_blocks,
          max_blocks = max_warm_blocks,
          ratio = warm_grid_ratio
        )
        
        p_old <- p_new
        cv_num <- 0
        
      }else{
        
        n_new <- n_new + n_blk
        
        if (!warm_done) {
          Obs <- rbind(Obs, cbind(matrix(X, n_blk), y))
        }
        
        X <- matrix(X, n_blk, length(idx_abs))
        
        for (k in 1:K_cv) {
          sub_obs <- (n_blk / K_cv * (k - 1) + 1):(n_blk / K_cv * k)
          Cy[k] <- (t_eff - n_blk) / t_eff * Cy[k] + n_blk / t_eff * mean((y[sub_obs])^2)
          CXX[, , k] <- (t_eff - n_blk) / t_eff * CXX[, , k] + 1 / t_eff * t(X[sub_obs, ]) %*% X[sub_obs, ]
          CXy[, k] <- (t_eff - n_blk) / t_eff * CXy[, k] + 1 / t_eff * t(X[sub_obs, ]) %*% y[sub_obs]
        }
        
      }
      
      # initial estimate
      # data-driven warm-up stage
      if (!warm_done) {
        
        m_now <- n_new / n_blk
        
        check_now <- FALSE
        if (warm_grid_pos <= length(warm_grid)) {
          check_now <- (m_now >= warm_grid[warm_grid_pos])
        }
        
        if (check_now) {
          
          X_warm <- Obs[, -ncol(Obs), drop = FALSE]
          
          warm_out <- warmup_check(
            X_obs = X_warm,
            grid_size = length(warm_grid),
            delta_sigma = delta_sigma,
            c_sigma = c_sigma,
            q_cov = q_cov,
            R_sigma = R_sigma,
            warm_const = warm_const,
            radius_mult = radius_mult
          )
          
          if (warm_out$stop) {
            
            warm_done <- TRUE
            iota0_obs <- n_new
            iota0_o <- c(iota0_o, iota0_obs)
            M_upper_o <- c(M_upper_o, warm_out$M_upper)
            warm_rhs_o <- c(warm_rhs_o, warm_out$rhs)
            
            cvfit <- cv.glmnet(
              Obs[, -ncol(Obs), drop = FALSE],
              Obs[, ncol(Obs)]
            )
            modelfit <- glmnet(
              Obs[, -ncol(Obs), drop = FALSE],
              Obs[, ncol(Obs)],
              lambda = cvfit$lambda.min
            )
            beta <- as.vector(modelfit$beta)
            rm(cvfit, modelfit)
            
            print(paste0(
              "-------------- warm-up done: N=", iota0_obs,
              " M_upper=", round(warm_out$M_upper, 3),
              " rhs=", round(warm_out$rhs, 1)
            ))
            
          } else {
            warm_grid_pos <- min(warm_grid_pos + 1, length(warm_grid))
          }
        }
        
        if (!warm_done) {
          next
        }
      }
      
      pt <- length(idx_imp)
      iota_star_obs <- ceiling(c_h * length(I))
      
      if (t_eff <= iota_star_obs || t_eff <= iota0_obs){
        
        Cy_total <- sum(Cy); CXy_total <- rowSums(CXy); CXX_total <- apply(CXX, c(1,2), sum)
        tmpt_p <- pt; #length(intersect(I,idx_abs))
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
            print('-------------- CV!!')
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
              print(C_lam)
            }
            C_lam_opt <- C_lam_set[which.min(err)]
            if(C_lam_opt==C_lam_opt0){cv_num <- cv_num+1}else{cv_num <- 0}
            C_lam_opt0 <- C_lam_opt
          }
        }
        if(tmpt_p==1){lam <-0}
        
        lam <- C_lam_opt*sigma_hat*sqrt(log(tmpt_p)^(delta)/(t_eff))
        # lam_vec <- rep(0,sum(idx_abs<min(I)))
        # if(tmpt_p>0){
        #   lam_vec <- c(lam_vec, rep(lam, tmpt_p))
        # }
        lam_vec <- rep(lam, tmpt_p)
        beta <- proximal(CXy_total, CXX_total, lam_vec, beta)
        
      }else{
        print('-------------- hard threshold!!')
        
        Cy_total <- sum(Cy)
        CXy_total <- rowSums(CXy)
        CXX_total <- apply(CXX, c(1, 2), sum)
        
        ## Candidate constants for hard thresholding
        C_hard_set <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3)
        
        cv_hard <- cv_hard_const(
          Cy = Cy,
          CXy = CXy,
          CXX = CXX,
          c_set = C_hard_set,
          N_eff = t_eff
        )
        
        C_hard_opt <- cv_hard$c_opt
        
        gram <- ginv(CXX_total)
        beta <- as.vector(gram %*% matrix(CXy_total, pt, 1))
        
        sigma_hat <- sqrt(
          Cy_total -
            2 * matrix(beta, 1) %*% CXy_total +
            matrix(beta, 1) %*% CXX_total %*% beta
        )
        sigma_hat <- min(sigma_hat, 10)
        
        beta1 <- C_hard_opt / sqrt(t_eff) *
          sqrt(pmax(diag(gram), 0)) *
          as.vector(sigma_hat)
        
        beta <- beta * sapply(
          1:pt,
          function(j) {
            abs(beta[j]) > abs(beta1[j])
          }
        )
        
      }
      
      # store expanded beta
      {
        err_beta <- sum((beta-alpha)^2)
        err_psue <- sum((beta-alpha_psue)^2)
        err_y <- Cy_total - 2*matrix(beta,1)%*%CXy_total + matrix(beta,1)%*%CXX_total%*%beta
        idx_imp <- unique(which(beta!=0))
        # if(t<=max(tau)){
        #   thre <- sort(abs(CXy)[-idx_imp], decreasing = TRUE)[min(K_cor,  length(idx_abs)-length(idx_imp))]
        #   idx_imp2 <- which(abs(CXy)>=thre)
        #   idx_imp <- sort(unique(c(idx_imp2, idx_imp)))
        # }
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
        print(paste0('t=',t,' err_beta=',round(err_beta,3),
                     ' err_y=',round(err_y,3),' pt=',pt, ' miss=',miss1))
        # if(t>eta[4] & t<eta[5]){print(alpha_psue)}
        
      }
    }
    
    if(!dir.exists(prefix)){dir.create(prefix)}
    save(err_beta_o, err_psue_o, time_o, pt_o, est_time, miss, # idx_abs_ls,
         file = paste0(prefix,'/','sd',sd,'.Rdata'))
    
  }
stopCluster(cl)

