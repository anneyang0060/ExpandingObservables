setwd('.../PM2.5')
lats <- c(); lons <- c()
for(i in 1001:2557){
  data <- read.table(paste0('China_',i,'.txt'),header = F,sep=",")
  lons <- c(lons, as.numeric(unique(data$V3)[2]))
  lats <- c(lats, as.numeric(unique(data$V4)[2]))
}
idx <- which(lats>=39.4 & lats <= 41.6 & lons>=115.7 & lons <= 117.4)

### average of multi-spots
y <- 0
for(i in c(1001:2557)[idx]){
  data <- read.table(paste0('China_',i,'.txt'),header = F,sep=",")
  y1 <- data$V2[which(data$V1=="1991-01-01"):which(data$V1=="2018-12-31")]
  y1 <- as.numeric(y1)
  y <- y + y1
}
dates <- data$V1[which(data$V1=="1991-01-01"):which(data$V1=="2018-12-31")]
y <- y/length(idx)

######## align dates #############
full_dates <- seq(as.Date('1991-01-01'), as.Date('2018-12-31'), by = 'day')
del_dates <- c(as.Date('1992-12-31'), as.Date('1996-12-31'),as.Date('2000-12-31'),
               as.Date('2004-12-31'),as.Date('2008-12-31'),as.Date('2012-12-31'),
               as.Date('2016-12-31'))
del_idx <- c()
for(d in del_dates){
  del_idx <- c(del_idx, which(full_dates==d))
}
full_dates <- full_dates[-del_idx]
y <- y[-del_idx]
miss_dates <- as.Date(setdiff(as.Date(full_dates), as.Date(dates)))
miss_idx <- c()
for(d in miss_dates){
  miss_idx <- c(miss_idx, which(full_dates==d))
}
y_full <- rep(0,length(full_dates))
y_full[-miss_idx] <- y
for(i in miss_idx){
  y_full[i] <- (y_full[i-1] + y_full[i+1])/2
}
setwd('D:\\datasets\\PM2.5\\PM2.5 regression\\data')
for(i in 1:28){
  year <- 1990 + i
  y_year <- y_full[((i-1)*365+1):(i*365)]
  save(y_year, file = paste0('PM2.5_', year, '_average.Rdata'))
}

######### scale #################
setwd()
features <- c('lrad', 'prec', 'pres', 'shum', 'srad', 'temp', 
              'wind', 'O3', 'PM10', 'NO2', 'CO', 'SO2')
M_X <- c()
SD_X <- c()
for(feature in features){
  if(feature %in% c('lrad', 'prec', 'pres', 'shum', 'srad', 'temp', 'wind')){year1 <- 1991}
  if(feature %in% c('O3', 'PM10')){year1 <- 2000}
  if(feature %in% c('NO2')){year1 <- 2008}
  if(feature %in% c('CO', 'SO2')){year1 <- 2013}
  X <- c()
  for(year in year1:2018){
    load(paste0('realdata/data/PM2.5/', feature,'_',year,'.Rdata'))
    X1 <- X_year[1:365,]
    X <- rbind(X,X1)
  }
  M_X <- c(M_X, colMeans(X))
  SD_X <- c(SD_X, sapply(1:ncol(X), function(i){sd(X[,i])}))
}
y <- c()
for(year in 1991:2018){
  load(paste0('realdata/data/PM2.5/PM2.5_',year,'_average.Rdata'))
  y <- c(y,y_year)
}
M_y <- mean(y)
SD_y <- sd(y)
save(M_X, SD_X, M_y, SD_y, file = 'realdata/data/PM2.5/statistics_log.Rdata')





