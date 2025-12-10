rm(list=ls())
# setwd('.../code')
setwd('C:/Users/annie/OneDrive/ImportantFiles/projects/1. InProgress/Online HD/fixed model/code')
sds <- 1:100

############################### increasing dimension results #############################
tau <- c(0, seq(2500,50000,2500), seq(55000,70000,5000), 8e4, 1e5, 1.4e5, 1e7)
eta <- c(0, 2500,5000,12500,30000,80000, 1e7)
load(paste0('res/FM-ID_rho0/sd1.Rdata'))
idx <- c(seq(1,3801,20), which(pt_o>200),
         as.vector(sapply(2:28,function(i){which(est_time == tau[i])})), length(est_time))
idx <- sort(idx)
rhos <- c(0, 0.3, 0.5, 0.7)

#### iota0=6
err_full <- c(); time_full <- c(); pt_full <- c()
for(rho in rhos){
  err <- c(); err_psue <- c(); pt <- c(); runing_times <- c()
  for(sd in sds){
    load(paste0('res/FM-ID_rho',rho,'/sd',sd,'.Rdata'))
    if(max(miss)>0){print(paste0('rho',rho,'/sd',sd,'.Rdata'))}
    err_beta_o <- err_beta_o[idx] 
    err_psue_o <- err_psue_o[idx]
    time_o <- time_o[idx]
    pt_o <- pt_o[idx]
    est_time <- est_time[idx]
    err <- cbind(err, err_beta_o)
    err_psue <- cbind(err_psue, err_psue_o)
    pt <- cbind(pt, pt_o)
    runing_times <- cbind(runing_times, time_o)
}
  err_full <- rbind(err_full, rowMeans(err_psue))
  pt_full <- rbind(pt_full, rowMeans(pt))
  time_full <- rbind(time_full, rowMeans(runing_times))
}

cc <- 1.5
pdf('fig/err_id.pdf', height = 6, width = 8)
par(mfrow = c(4,1), mai = c(0.6,0.6,0.2,0.1))
for(i in 1:4){
  plot(est_time/50, err_full[i,], type='l', main = paste0('rho=', rhos[i]),
       xlab = 'data block', ylab = 'error',
       cex.axis = cc, cex.lab = cc, cex.main = cc)
}
dev.off()
pdf('fig/time_id.pdf', height = 6, width = 8)
par(mfrow = c(4,1), mai = c(0.6,0.6,0.2,0.1))
for(i in 1:4){
  plot(est_time[-c(1,2)]/50, time_full[i,-c(1,2)], type='l', main = paste0('rho=', rhos[i]),
       xlab = 'data block', ylab = 'time',
       cex.axis = cc, cex.lab = cc, cex.main = cc)
}
dev.off()
pdf('fig/dim_id.pdf', height = 6, width = 8)
par(mfrow = c(4,1), mai = c(0.6,0.6,0.2,0.1))
for(i in 1:4){
  plot(est_time/50, pt_full[i,], type='l', main = paste0('rho=', rhos[i]),
       xlab = 'data block', ylab = 'd_t',
       cex.axis = cc, cex.lab = cc, cex.main = cc)
}
dev.off()

err_full[,ncol(err_full)]
       
#### simu different iota0 (3, 6, 10)
rho <- 0.5; n0s <- c(150, 300, 500)
err_full <- c(); time_full <- c(); pt_full <- c(); est_time_full <- c()
for(n0 in n0s){
  err <- c(); err_psue <- c(); pt <- c(); runing_times <- c()
  for(sd in sds){
    if(n0==300){
      filename <-(paste0('res/FM-ID_rho',rho,'/sd',sd,'.Rdata'))
    }else{
      filename <-(paste0('res/FM-ID_rho',rho,'/n0_',n0,'_sd',sd,'.Rdata'))
    }
    load(filename)
    idx <- c(seq(1,3801,20), which(pt_o>200),
             as.vector(sapply(2:28,function(i){which(est_time == tau[i])})), length(est_time))
    idx <- sort(idx)
    if(max(miss)>0){print(filename)}
    err_beta_o <- err_beta_o[idx] 
    err_psue_o <- err_psue_o[idx]
    time_o <- time_o[idx]
    pt_o <- pt_o[idx]
    err <- cbind(err, err_beta_o)
    err_psue <- cbind(err_psue, err_psue_o)
    pt <- cbind(pt, pt_o)
    runing_times <- cbind(runing_times, time_o)
  }
  est_time <- est_time[idx]
  err_full <- rbind(err_full, rowMeans(err_psue))
  pt_full <- rbind(pt_full, rowMeans(pt))
  time_full <- rbind(time_full, rowMeans(runing_times))
  est_time_full <- rbind(est_time_full, est_time)
}

cc <- 1.5
pdf('fig/err_id_diff_iota0.pdf', height = 4.5, width = 8)
par(mfrow = c(3,1), mai = c(0.6,0.6,0.2,0.1))
for(i in 1:3){
  plot(est_time_full[i,]/50, err_full[i,], type='l', main = paste0('iota^0=', n0s[i]/50),
       xlab = 'data block', ylab = 'error',
       cex.axis = cc, cex.lab = cc, cex.main = cc)
}
dev.off()
pdf('fig/time_id_diff_iota0.pdf', height = 4.5, width = 8)
par(mfrow = c(3,1), mai = c(0.6,0.6,0.2,0.1))
for(i in 1:3){
  plot(est_time_full[i,][-c(1,2)]/50, time_full[i,-c(1,2)], type='l', main = paste0('iota^0=', n0s[i]/50),
       xlab = 'data block', ylab = 'time',
       cex.axis = cc, cex.lab = cc, cex.main = cc)
}
dev.off()
pdf('fig/dim_id_diff_iota0.pdf', height = 4.5, width = 8)
par(mfrow = c(3,1), mai = c(0.6,0.6,0.2,0.1))
for(i in 1:3){
  plot(est_time_full[i,]/50, pt_full[i,], type='l', main = paste0('iota^0=', n0s[i]/50),
       xlab = 'data block', ylab = 'd_t',
       cex.axis = cc, cex.lab = cc, cex.main = cc)
}
dev.off()

