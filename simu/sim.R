rm(list=ls())
setwd('.../code')

library(foreach)
library(doParallel)
library(mvnfast)
library(glmnet)
library(MASS)

source('FNS.R')
sds <- 1:100

######## increasing dimension
for(rho in c(0.7, 0.5, 0.3, 0)){
  source('simu/Setting_FM-ID.R')
  source('simu/R-AVAS_FM-ID.R')
}



