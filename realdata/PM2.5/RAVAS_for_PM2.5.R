rm(list=ls())
setwd('C:/Users/annie/OneDrive/ImportantFiles/projects/1. InProgress/Online HD/fixed model/code')
# setwd('F:/datasets/PM2.5/PM2.5 regression')
# coordinate of PM2.5: 116.366,39.8673; bottom to top 8th; left to right 13th


###################### online train-test: data-driven warm-up #####################
library(glmnet)
library(MASS)

## Use the updated FNS.R corresponding to FNS(4).R.
source('FNS.R')

## Data-driven warm-up parameters.
## warm_target_obs is set to the old fixed warm-up length used in this script.
## Change it if you want a different finite-sample calibration.
warm_target_obs <- 700
warm_target_ratio <- 0.95
warm_const <- NA_real_

warm_grid_ratio <- 1.25
delta_sigma <- 0.10
c_sigma <- 1.5
q_cov <- 0
R_sigma <- NULL
radius_mult <- 1
warm_grid_size_default <- 2000

## If you want to force at least a few soft-selection updates after warm-up,
## set this to a positive integer. Setting it to 0 keeps the original transition
## rule as closely as possible.
soft_min_blocks <- 0

next_change_year_pm25 <- function(year) {
  if (year < 2000) {
    return(2000)
  } else if (year < 2008) {
    return(2008)
  } else if (year < 2013) {
    return(2013)
  } else {
    return(2018)
  }
}

{
  ####### training ############
  {
    s <- 1; I_set <- list(); W_set <- c()
    Obs <- c(); lam_set <- c()
    sigma_o <- c(); p_old <- 0; K <- 0
    st <- 0; pt <- 0
    idx_imp <- c(); idx_abs <- c(); idx_del <- c()
    N_obs_old <- 0; p_obs_old <- 0
    err_test_o <- c()
    pt_o <- c(); err_y_o <- c(); time_o <- c(); resd <- c()
    beta_o <- list(); R2_o <- c(); y_hat <- c(); y_real <- c()
    sigma_hat <- 1; years <- c(); sigma_hat0 <- 0
    K_cv <- 2; C_lam_set <- (1:10)/5; C_lam_opt0 <- 0; cv_num <- 0

    ## warm-up tracking
    warm_done <- FALSE
    iota0_obs <- NA_real_
    n0 <- NA_real_
    warm_grid <- NULL
    warm_grid_pos <- 1
    iota0_o <- c()
    warm_M_upper_o <- c()
    warm_M_hat_o <- c()
    warm_rhs_o <- c()
    warm_const_o <- c()
    warm_stop_year_o <- c()
    warm_stop_block_o <- c()
  }

  for (year in (1991:2017))
  {
    ## load yearly data
    {
      ## observed feature
      if(year < 2000){
        features <- c('lrad', 'prec', 'pres', 'shum', 'srad', 'temp', 'wind')
      }
      if(year >= 2000 & year < 2008){
        features <- c('lrad', 'prec', 'pres', 'shum', 'srad',
                      'temp', 'wind', 'O3', 'PM10')
      }
      if(year >= 2008 & year < 2013){
        features <- c('lrad', 'prec', 'pres', 'shum', 'srad',
                      'temp', 'wind', 'O3', 'PM10', 'NO2')
      }
      if(year >= 2013){
        features <- c('lrad', 'prec', 'pres', 'shum', 'srad', 'temp',
                      'wind', 'O3', 'PM10', 'NO2', 'CO', 'SO2')
      }

      ## feature
      X1 <- c()
      for(k_fea in 1:length(features))
      {
        feature <- features[k_fea]
        load(paste0('realdata/PM2.5/data/', feature,'_',year,'.Rdata'))
        if(k_fea<=7){X1 <- cbind(X1, X_year[1:364,])}else{X1 <- cbind(X1, X_year[2:365,])}
      }
      idx_del <- union(idx_del, which(is.na(X1[1,])))
      X1 <- X1[,-idx_del]
      X1 <- log(X1)
      X_year <- sapply(1:ncol(X1), function(i){(X1[,i]-mean(X1[,i]))/sd(X1[,i])})

      ## response
      load(paste0('realdata/PM2.5/data/PM2.5_',year,'_average.Rdata'))
      y_year <- y_year[2:365]
      y_year <- log(y_year)
      y_year <- (y_year - mean(y_year))/sd(y_year)
    }

    print(paste0('year=', year))

    for(k in 1:52){

      print(k)

      # read in data block
      {
        y <- y_year[((k-1)*7+1):(k*7)]
        X <- X_year[((k-1)*7+1):(k*7),]
        n_blk <- nrow(X)
        p_obs <- ncol(X)
        delta_p <- p_obs - p_obs_old
        N_obs <- N_obs_old + n_blk
        N_obs_old <- N_obs
      }

      t_start <- Sys.time()

      if(delta_p > 0){

        ## Start a new RAVAS cycle and reset the data-driven warm-up rule.
        warm_done <- FALSE
        iota0_obs <- NA_real_
        n0 <- NA_real_
        warm_grid_pos <- 1

        next_year <- next_change_year_pm25(year)
        max_warm_blocks <- max(1, (next_year - year) * 52)
        min_warm_blocks <- max(1, ceiling(warm_target_obs / n_blk))
        warm_grid <- make_warm_grid(
          min_blocks = min_warm_blocks,
          max_blocks = max_warm_blocks,
          ratio = warm_grid_ratio
        )

        N_eff <- n_blk
        W <- n_blk/N_eff
        I <- (p_obs_old+1):p_obs
        idx_abs <- c(idx_abs, I)
        X <- X[,idx_abs, drop = FALSE]
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
            CXX[,,k_cv] <- 1/N_eff*t(X[sub_obs,,drop=FALSE])%*%X[sub_obs,,drop=FALSE]
            CXy[,k_cv] <- 1/N_eff*t(X[sub_obs,,drop=FALSE])%*%y[sub_obs]
          }
        }

        p_obs_old <- p_obs

      }else{

        N_eff <- N_eff + n_blk
        X <- X[,c(idx_abs), drop = FALSE]

        ## Store raw data only during warm-up.
        if(!warm_done){
          Obs <- rbind(Obs, cbind(matrix(X,n_blk),y))
        }

        X <- matrix(X, nrow = n_blk)
        for(k_cv in 1:K_cv){
          sub_obs <- (n_blk/K_cv*(k_cv-1)+1):(n_blk/K_cv*k_cv)
          Cy[k_cv] <- (N_eff-n_blk)/N_eff * Cy[k_cv] + n_blk/N_eff * mean((y[sub_obs])^2)
          CXX[,,k_cv] <- (N_eff-n_blk)/N_eff * CXX[,,k_cv] + 1/N_eff * t(X[sub_obs,,drop=FALSE])%*%X[sub_obs,,drop=FALSE]
          CXy[,k_cv] <- (N_eff-n_blk)/N_eff * CXy[,k_cv] + 1/N_eff * t(X[sub_obs,,drop=FALSE])%*%y[sub_obs]
        }
      }

      ############################################################
      ## Data-driven warm-up and initial estimate
      ############################################################
      if(!warm_done){

        m_now <- ceiling(N_eff / n_blk)
        check_now <- FALSE
        if(!is.null(warm_grid) && warm_grid_pos <= length(warm_grid)){
          check_now <- (m_now >= warm_grid[warm_grid_pos])
        }

        if(check_now){

          X_warm <- Obs[, -ncol(Obs), drop = FALSE]
          grid_size_now <- ifelse(is.null(warm_grid), warm_grid_size_default, length(warm_grid))

          ## Calibrate warm_const once at the first eligible warm-up check so that
          ## the nominal stopping scale is comparable to the previous fixed n0.
          if(!is.finite(warm_const)){

            warm_tmp <- warmup_check(
              X_obs = X_warm,
              grid_size = grid_size_now,
              delta_sigma = delta_sigma,
              c_sigma = c_sigma,
              q_cov = q_cov,
              R_sigma = R_sigma,
              warm_const = 1,
              radius_mult = radius_mult
            )

            warm_const <- warm_target_ratio * warm_tmp$N / warm_tmp$rhs

            cat(
              "[warm-up calibration]",
              "year =", year,
              "block =", k,
              "N =", warm_tmp$N,
              "a =", warm_tmp$a,
              "M_upper =", round(warm_tmp$M_upper, 3),
              "rhs0 =", round(warm_tmp$rhs, 1),
              "warm_const =", signif(warm_const, 4),
              "\n"
            )
          }

          warm_out <- warmup_check(
            X_obs = X_warm,
            grid_size = grid_size_now,
            delta_sigma = delta_sigma,
            c_sigma = c_sigma,
            q_cov = q_cov,
            R_sigma = R_sigma,
            warm_const = warm_const,
            radius_mult = radius_mult
          )

          cat(
            "[warm-up check]",
            "year =", year,
            "block =", k,
            "N_eff =", N_eff,
            "a =", warm_out$a,
            "M_hat =", round(warm_out$M_hat, 3),
            "r_sigma =", round(warm_out$r_sigma, 3),
            "M_upper =", round(warm_out$M_upper, 3),
            "rhs =", round(warm_out$rhs, 1),
            "N/rhs =", round(warm_out$N / warm_out$rhs, 3),
            "\n"
          )

          if(warm_out$stop){

            warm_done <- TRUE
            iota0_obs <- N_eff
            n0 <- iota0_obs

            iota0_o <- c(iota0_o, iota0_obs)
            warm_M_upper_o <- c(warm_M_upper_o, warm_out$M_upper)
            warm_M_hat_o <- c(warm_M_hat_o, warm_out$M_hat)
            warm_rhs_o <- c(warm_rhs_o, warm_out$rhs)
            warm_const_o <- c(warm_const_o, warm_const)
            warm_stop_year_o <- c(warm_stop_year_o, year)
            warm_stop_block_o <- c(warm_stop_block_o, k)

            ## Keep the same initialization convention as the original PM2.5 code.
            lam <- (2e-3*(year<2000)+1e-3*(year >= 2000 & year < 2008)
                    +5e-4*(year >= 2008))
            modelfit <- glmnet(
              Obs[, -ncol(Obs), drop = FALSE],
              Obs[, ncol(Obs)],
              lambda = lam
            )
            beta <- as.vector(modelfit$beta)
            sigma_hat <- min(
              sigma_hat,
              sqrt(mean((Obs[, ncol(Obs)] -
                           Obs[, -ncol(Obs), drop = FALSE] %*% beta)^2))
            )
            rm(modelfit)

            print(paste0(
              "initial estimation by data-driven warm-up, iota0_obs=", iota0_obs
            ))

          }else{
            warm_grid_pos <- min(warm_grid_pos + 1, length(warm_grid))
          }
        }

        if(!warm_done){
          next
        }else{
          ## Same as the original code: after the initial estimate is computed,
          ## the algorithm starts selection from the next data block.
          next
        }
      }

      pt <- length(idx_imp)
      print(paste0('dim=',pt, ' sigma_hat=', round(sigma_hat,3)))

      iota_star_eff <- max(length(I), iota0_obs + soft_min_blocks * n_blk)

      if(N_eff <= iota_star_eff){

        Cy_total <- sum(Cy); CXy_total <- rowSums(CXy); CXX_total <- apply(CXX, c(1,2), sum)
        tmpt_p <- pt
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
              err1 <- numeric(K_cv)

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
        beta1 <- rep(0.5/sqrt(n0)*as.vector(sigma_hat),pt)
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
        CXX <- CXX[idx_imp, idx_imp, , drop = FALSE]
        CXy <- CXy[idx_imp, , drop = FALSE]
        pt_o <- c(pt_o, pt)
        time_o <- c(time_o, as.numeric(t_end-t_start))
        beta_o <- c(beta_o,list(beta_expan))
        R2 <- 1-err_y
        R2_o <- c(R2_o, R2)
        years <- c(years, year)
        sigma_o <- c(sigma_o,as.vector(sigma_hat))
      }
    }
  }

  ######### test ############
  y <- c(); X <- c()
  for(year in 2018)
  {
    features <- c('lrad', 'prec', 'pres', 'shum', 'srad', 'temp', 
                  'wind', 'O3', 'PM10', 'NO2', 'CO', 'SO2')

    ## feature
    X1 <- c()
    for(k_fea in 1:length(features))
    {
      feature <- features[k_fea]
      load(paste0('realdata/PM2.5/data/', feature,'_',year,'.Rdata'))
      if(k_fea<=7){X1 <- cbind(X1, X_year[1:364,])}else{X1 <- cbind(X1, X_year[2:365,])}
    }
    idx_del <- union(idx_del, which(is.na(X1[1,])))
    X1 <- X1[,-idx_del]
    X1 <- log(X1)
    X_year <- sapply(1:ncol(X1), function(i){(X1[,i]-mean(X1[,i]))/sd(X1[,i])})
    X <- rbind(X, X_year)

    ## response
    load(paste0('realdata/PM2.5/data/PM2.5_',year,'_average.Rdata'))
    y_year <- log(y_year[2:365])
    y_year <- (y_year-mean(y_year))/sd(y_year)
    y <- c(y,y_year)
  }

  for(i in 1:length(beta_o)){
    beta_expan <- rep(0,dim(X_year)[2])
    beta_expan[1:length(beta_o[[i]])] <- beta_o[[i]]
    err_test_o <- c(err_test_o, mean((y-X%*%beta_expan)^2))
  }

  save(idx_del, years, beta_o, err_test_o, err_y_o, time_o, pt_o, sigma_o,
       iota0_o, warm_M_upper_o, warm_M_hat_o, warm_rhs_o, warm_const_o,
       warm_stop_year_o, warm_stop_block_o,
       file = paste0('res/res_PM2.5.Rdata'))
}

####################### plot #####################

library(ggplot2)
library(ggpubr)
library(png)
library(pheatmap)

load('res/res_PM2.5.Rdata')

pdf('fig/test_err_PM2.5.pdf', width = 7, height=5)
plot(err_test_o[err_test_o<1.35], ylim = c(0.1, 1.35),xaxt='n',
     ylab='test error', xlab = 'data of year')
ts <-c(1)
for(year in c(1999,2007,2012)){ts <- c(ts, 1+max(which(years[err_test_o<1.35]==year)))}
axis(side=1, at=ts, labels = c(1991,2000,2008,2013))
for(t in ts){lines(rep(t,4), seq(0.1,1.2,length.out=4), lty=2)}
dev.off()

features <- c('lrad', 'prec', 'pres', 'shum', 'srad', 'temp', 
              'wind', 'O3', 'PM10', 'NO2', 'CO', 'SO2')
for(year in c(1999,2007,2012,2017)){
  k <- max(which(years==year))
  N_fea <- 7+2*(year>=2000)+1*(year>=2008)+2*(year>=2013)
  beta_expan <- rep(0,300*N_fea)
  idx_del_1 <- idx_del[idx_del<=length(beta_expan)]
  beta_expan[-idx_del_1] <- beta_o[[k]]
  beta_expan <- array(beta_expan, dim = c(15,20,N_fea))
  effect_var <- c()
  for(i in 1:N_fea){
    effect_var <- c(effect_var, sum(abs(beta_expan[,,i])))
  }
  print(year)
  print(sort(effect_var, decreasing = TRUE))
  print(features[order(effect_var, decreasing = TRUE)])
  print(order(effect_var, decreasing = TRUE))
}


## important grids
load(paste0('res/res_PM2.5.Rdata'))
tar_lat <- seq(37.2, 43, 0.4)
tar_lon <- seq(110.1, 120, 0.5)
coor <- cbind(rep(tar_lon,each = 15), rep(tar_lat, 20))
for(year in c(1999,2007,2012,2017))
{
  if(year==1999){thre <- 0.2}
  if(year==2007 | year==2012){thre <- 0.3}
  if(year==2017){thre <- 0.26}
  lats <- c(); lons <- c(); colors1 <- c(); sizes <- c()
  
  {
    if(year <= 2000){
      features <- c('lrad', 'prec', 'pres', 'shum', 'srad', 'temp', 'wind')
      ks <- c(1,3,7)
    }
    if(year > 2000 & year <= 2008){
      features <- c('lrad', 'prec', 'pres', 'shum', 'srad','temp', 'wind', 'O3', 'PM10')
      ks <- c(6,8,9)
    }
    if(year > 2008 & year <= 2013){
      features <- c('lrad', 'prec', 'pres', 'shum', 'srad', 'temp', 'wind', 'O3', 'PM10', 'NO2')
      ks <- c(8,9,10)
    }
    if(year > 2013){
      features <- c('lrad', 'prec', 'pres', 'shum', 'srad', 'temp', 
                    'wind', 'O3', 'PM10', 'NO2', 'CO', 'SO2')
      ks <- 10:12
    }
  }
  
  idx <- max(which(years==year))
  beta<-beta_o[[idx]]
  idx_del1 <- idx_del[idx_del<= 300*length(features)] 
  for(k in ks){
    beta1 <- rep(0, 300*length(features))
    beta1[-idx_del1] <- beta
    beta1 <- beta1[((k-1)*300+1):(k*300)]
    beta1 <- matrix(beta1,15,20)
    beta1[,20] <- 0
    if(k==1){color1 <- 'blue'}
    if(k==4){color1 <- 'lightblue'}
    if(k==6){color1 <- 'cornflowerblue'}
    if(k==7){color1 <- 'blue3'}
    if(k==8){color1 <- 'yellow'}
    if(k==9){color1 <- 'pink'}
    if(k==10){color1 <- 'darkorange'}
    if(k==11){color1 <- 'indianred1'}
    if(k==12){color1 <- 'red'}
    if(k<8){
      lats <- c(lats, coor[which(beta1< -thre),2])
      lons <- c(lons, coor[which(beta1< -thre),1])
      colors1 <- c(colors1, rep(color1, sum(beta1< -thre)))
      sizes <- c(sizes, abs(beta1[which(beta1< -thre)]))
    }else{
      lats <- c(lats, coor[which(beta1>thre),2])
      lons <- c(lons, coor[which(beta1>thre),1])
      colors1 <- c(colors1, rep(color1, sum(beta1>thre)))
      sizes <- c(sizes, abs(beta1[which(beta1>thre)]))
    }
  }
  sizes <- sizes*10
  if(year==1999){sizes <- sizes/2}
  image <- readPNG("realdata/PM2.5/map.PNG")
  data <- data.frame(x = lons, y = lats)
  if(year==1999){
    pp <- (ggplot(data,aes(x,y)) + ylim(37, 43.2) + xlim(110,120) 
           +  background_image(image) 
           + geom_point(color=colors1,size=(2.5*(sizes))))
  }
  if(year==2007 | year==2012){
    pp <- (ggplot(data,aes(x,y)) + ylim(37, 43.2) + xlim(110,120) 
           +  background_image(image) 
           + geom_point(color=colors1,size=((sizes)^1.2)))
  }
  if(year==2017){
    pp <- (ggplot(data,aes(x,y)) + ylim(37, 43.2) + xlim(110,120) 
           +  background_image(image) 
           + geom_point(color=colors1,size=(1.5*(sizes)^1.3)))
  }
  print(pp)
}

plot(c(),c(),ylim=c(0,1),xlim=c(0,1))
legend('topright',features[c(3,1,6:12)],
       col=c('lightblue','blue','cornflowerblue','blue3','yellow',
             'pink','darkorange','indianred1','red'), pch=16 )
