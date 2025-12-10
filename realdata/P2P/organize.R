rm(list=ls())
gc()
setwd('C:/Users/annie/Downloads/accepted_2007_to_2018Q4.csv')
load('pk.Rdata')
plot(pk)

pk_unique <- unique(pk)
new_batch <- c()
for(p in pk_unique){
  new_batch <- c(new_batch, max(which(pk==p)))
}

features <- list()
for(k in new_batch){
  file_cur <- paste0('P2Pdata_',k,'.csv')
  data_cur <- read.csv(file_cur)
  features <- c(features, list(names(data_cur)))
}


load('P2P_train_results_1.Rdata')
p_k <- c()
for(i in 1:length(beta_o)){
  p_k <- c(p_k, length(beta_o[[i]]))
}

final_b <- c()
for(i in new_batch[1:4]){
  final_b <- c(final_b, which(k_o==(i-1)))
}
final_b <- c(final_b, length(k_o))

for(b in final_b){
  print(paste0('################## features of beta_',b))
  file_cur <- paste0('P2Pdata_',b,'.csv')
  data_cur <- read.csv(file_cur)
  data_cur <- data_cur[,-c(1,2,3)]
  # covariate
  data_cur <- subset(data_cur, select = -c(int_rate))
  X <- as.matrix(data_cur)
  p0 <- ncol(X)
  idx_del <- which(sapply(1:nrow(X), function(i){sum(is.na(X[i,]))})>0)
  if(length(idx_del)>0){
    X <- X[-idx_del, ]
  }
  p_interact0 <- which(names(data_cur)=='home_ownership_MORTGAGE')
  j1_indices <- 1:(p_interact0-1)
  j2_indices <- p_interact0:p0
  combinations <- expand.grid(j1 = j1_indices, j2 = j2_indices)
  
  imp_idx <- which(beta_o[[b]]!=0)
  imp_idx1 <- imp_idx[which(imp_idx <= p0)]
  imp_idx2 <- imp_idx[which(imp_idx > p0)]
  for(k in imp_idx1){
    print(paste0(names(data_cur)[k]))
  }
  for(k in imp_idx2){
    i <- k - p0
    j1_value <- combinations[i, "j1"]
    j2_value <- combinations[i, "j2"]
    print(paste0(names(data_cur)[j1_value], '*', names(data_cur)[j2_value]))
  }
}


load('P2P_test_results_1.Rdata')
blk <- c(100,500,800,1300,2000,2250)
idx <- c()
for (b in blk) {
  idx <- c(idx, which(k_o == b))
}
pt_o[idx]
round(time_o[idx],4)
round(err_test_o[idx],4)



