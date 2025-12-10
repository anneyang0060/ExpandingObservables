rm(list=ls())
# setwd('.../code')
setwd('C:/Users/annie/OneDrive/ImportantFiles/projects/1. InProgress/Online HD/fixed model/code')
sds <- 1:100

############################### fixed dimension results #############################
for(rho in c(0, 0.3, 0.5, 0.7))
  # rho <- 0
{
  ## read in results of different methods
  {
    err <- c(); pt <- c(); runing_times <- c();sigmas <- c()
    for(sd in sds){
      filename <- paste0('res/AVAS_rho', rho,'/sd',sd,'.Rdata')
      load(filename)
      err <- cbind(err, err_beta_o)
      pt <- cbind(pt, pt_o)
      runing_times <- cbind(runing_times, time_o)
      sigmas <- cbind(sigmas, sigma_o)
    }
    err_AVAS <- rowMeans(err)
    time_AVAS <- rowMeans(runing_times)
    p_AVAS <- rowMeans(pt)
    
    err <- c(); runing_times <- c(); dt <- c()
    for(sd in sds){
      filename <- paste0('res/BR_fan_rho', rho,'/sd',sd,'.Rdata')
      load(filename)
      err <- cbind(err, err_beta_o)
      runing_times <- cbind(runing_times, time_o)
      dt <- cbind(dt, dt_o)
    }
    err_BR <- rowMeans(err)
    time_BR <- rowMeans(runing_times)
    p_BR <- rowMeans(dt)
    
    err <- c(); runing_times <- c(); dt <- c()
    for(sd in sds){
      filename <- paste0('res/TSGD_rho', rho,'/sd',sd,'.Rdata')
      load(filename)
      err <- cbind(err, err_beta_o)
      runing_times <- cbind(runing_times, time_o)
      dt <- cbind(dt, dt_o)
    }
    err_TSGD <- rowMeans(err)
    time_TSGD <- rowMeans(runing_times)
    p_TSGD <- rowMeans(dt)
    
    err <- c(); runing_times <- c(); dt <- c()
    for(sd in sds){
      filename <- paste0('res/OLL_sun_rho', rho,'/sd',sd,'.Rdata')
      if(file.exists(filename)){
        load(filename)
      }else{next}
      err <- cbind(err, err_beta_o)
      runing_times <- cbind(runing_times, time_o)
      dt <- cbind(dt, dt_o)
    }
    err_OLL <- rowMeans(err)
    time_OLL <- rowMeans(runing_times)
    p_OLL <- rowMeans(dt)
    
    err1 <- c(); err2<-c(); err3<-c()
    time1 <- c(); time2 <- c()
    dt1 <- c(); dt2 <- c()
    for(sd in sds){
      filename <- paste0('res/ODL_huang_rho', rho,'/sd',sd,'.Rdata')
      load(filename)
      err1 <- cbind(err1, err_beta_b_o)
      err2 <- cbind(err2, err_beta_o)
      err3 <- cbind(err3, err_beta_g_o)
      time1 <- cbind(time1, time_o_u)
      time2 <- cbind(time2, time_o_bg)
      dt1 <- cbind(dt1, dt_o_b)
      dt2 <- cbind(dt2, dt_o_g)
    }
    
    # err_b <- rowMeans(err1)
    err_ODL <- rowMeans(err1)
    err_OR <- rowMeans(err3)
    time_ODL <- rowMeans(time1)
    time_OR <- rowMeans(time2)
    p_ODL <- rowMeans(dt1)
    p_OR <- rowMeans(dt2)
  }  
  
  ## error fig
  {
    cc <- 1.6
    pdf(paste0('fig/err_fd_rho',rho,'.pdf'),height = 7, width = 5)
    ylim1 <- range(log(c(err_AVAS, err_BR[2:length(err_BR)],err_TSGD[2:length(err_TSGD)],
                         err_ODL, err_OLL, err_OR)))
    plot(seq(210,10000,20),log(err_AVAS[seq(1,length(err_AVAS),2)]),
         ylim =  ylim1, main=paste0('rho=',rho),
         type='l',col='red',xlab = 'data block',ylab = 'log(MSE)',
         cex.axis = cc, cex.lab=cc, cex.main=cc,
         mai = c(3,3,1,1), mar = c(4,4,1,1))
    lines(seq(210,10000,10),log(err_BR[2:length(err_BR)]),col='blue')
    lines(seq(210,10000,10),log(err_TSGD[2:length(err_TSGD)]),col='green')
    lines(seq(300,10000,100),log(err_ODL))
    lines(seq(300,10000,100),log(err_OR),col='orange')
    lines(seq(300,10000,100),log(err_OLL), col='purple')
    if(rho==0.7){legend(6000,0,c('AVAS','BR','TSGD','ODL','OR','OLL'),
                        col=c('red','blue','green','black','orange','purple'),lty=1,
                        cex=(cc-0.3))}
    dev.off()
  }
  
  ## time fig
  {
    cc <- 1.6
    png(paste0('fig/time_fd_rho',rho,'.png'), height = 700*0.6, width = 500*0.6, res = 80*0.8)
    plot(seq(220,10000,10),(time_AVAS[-c(1)]), yaxt='n',main=paste0('rho=',rho),
         xlim = c(50,10100), ylim = c(-0.01,0.15),xaxs='i', yaxs='i',
         type='l',col='red',xlab = 'data block',ylab = 'running time',
         cex.axis = cc, cex.lab=cc, cex.main=cc)
    axis(side=2,at=c(0, 0.025, 0.05, 0.075, 0.09, 0.11, 0.13, 0.15),
         labels = c(0,  0.025, 0.05, 0.075, 0.09,  9,     12,  15), cex.axis=cc)
    lines(seq(210,10000,10),(time_BR[-1]),col='blue')
    lines(seq(230,10000,10),(time_TSGD[-c(1,2,3)]),col='green')
    lines(seq(400,10000,100),(time_ODL)[-1]/150+0.05)
    lines(seq(300,10000,100),(time_OLL), col='purple')
    lines(seq(500,10000,100),(time_OR)[-c(1,2)],col='orange')
    segments(50,0.09,50,0.1,col="white",lwd=1)
    segments(50,0.09,300,0.095,col='black',lwd=1)
    segments(50,0.1,300,0.105,col='black',lwd=1)
    if(rho==0.7){legend('topright',c('AVAS','BR','TSGD','ODL','OLL','OR'),
                        col=c('red','blue','green','black','purple','orange'),lty=1,
                        cex = (cc-0.3))}
    dev.off()
    
  }
  
  ## dim fig
  {
    cc <- 1.6
    pdf(paste0('fig/dim_fd_rho',rho,'.pdf'),height = 7, width = 5)
    plot(seq(220,10000,10),(p_AVAS[-c(1)]), main=paste0('rho=',rho),
         xlim = c(50,10100), ylim = c(0,1000),
         type='l',col='red',xlab = 'data block',ylab = 'No. selected variables',
         cex.axis = cc, cex.lab=cc, cex.main=cc,
         mai = c(3,3,1,1), mar = c(4,4,1,1))
    lines(seq(200,10000,10),(p_BR),col='blue')
    lines(seq(200,10000,10),(p_TSGD),col='green')
    lines(seq(300,10000,100),(p_ODL))
    lines(seq(300,10000,100),(p_OLL), col='purple')
    lines(seq(300,10000,100),(p_OR),col='orange')
    if(rho==0.7){legend(6000,800,c('AVAS','BR','TSGD','ODL','OR','OLL'),
                        col=c('red','blue','green','black','orange','purple'),lty=1,
                        cex=(cc-0.3))}
    dev.off()
  }
  
}
