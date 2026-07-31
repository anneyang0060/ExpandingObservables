rm(list=ls())
setwd('/home/yangying/ExpandingObs')

library(foreach)
library(doParallel)
library(mvnfast)
library(glmnet)
library(MASS)

source('FNS.R')
sds <- 1:100

######## increasing dimension
radius_mult <- 1
for(rho in c(0.7, 0.5, 0.3, 0)){
  source('simu/Setting_FM-ID.R')
  source('simu/R-AVAS_FM-ID.R')
}

####### fixed dimension
p <- 1000 # model dimension
s <- 10 # sparsity
Tmax <- 10000 # maximal number of batches
sigma_eps<- 1 # noise level
n0 <- 200
beta_true <- rep(0,p)
beta_true[1:10] <- rep(c(3,-3),s/2)
radius_mult <- 0.1

for(rho in c(0,  0.3, 0.7)){
  
  n_blk <- 10
  sigma_X <-matrix(0, p,p)
  for(k in 1:p){
    sigma_X[k,] <- sapply(1:p, function(l){rho^(abs(k-l))})
  }
  
  n_blk <- 10
  K_cv <- 2
  prefix <- paste0('res/AVAS_rho',rho)
  source('simu/AVAS_FM-FD.R')
  
  prefix <- paste0('res/BR_fan_rho',rho)
  source('supp/BR_fan_FM-FD.R')


  prefix <- paste0('res/TSGD_rho',rho)
  source('supp/TSGD_FM-FD.R')

  n_blk <- 100

  prefix <- paste0('res/ODL_huang_rho',rho)
  source('supp/ODL_huang_FM-FD.R')

  prefix <- paste0('res/OLL_sun_rho',rho)
  source('supp/OLL_sun.R')
  
}

rm(p,s,Tmax,sigma_eps,n0)







