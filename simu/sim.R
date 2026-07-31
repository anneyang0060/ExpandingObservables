rm(list=ls())
setwd('.../ExpandingObservables')

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








