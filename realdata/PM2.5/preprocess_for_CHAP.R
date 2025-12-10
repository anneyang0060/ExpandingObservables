rm(list=ls())
setwd('.../CHAP')
# library(raster)
library(ncdf4)

pollutions <- c('O3','PM2.5','PM10','NO2','CO','SO2')

for(k in 1:6)
# k <- 1
{
  pollution <- pollutions[k]
  
  if(pollution %in% c('O3','PM2.5','PM10')){year1 <- 2000}
  if(pollution %in% c('NO2')){year1 <- 2008}
  if(pollution %in% c('CO','SO2')){year1 <- 2013}
  
  for(year in year1:2020)
  # year <- year1
  {
    save_file <- paste0('F:/datasets/air pollution/',pollution,'_',year,'.Rdata')
    if(!file.exists(save_file)){
      setwd(paste0('F:/datasets/air pollution/',pollution,'(',year1,'-2022)/',year))
      dates <- seq(as.Date(paste0(year,'-01-01')),as.Date(paste0(year,'-12-31')), by = 'day')
      X_year <- c()
      
      for(dd in dates)
        # dd <- dates[1]
      {
        dd1 <- as.Date(dd)
        dd1 <- gsub('-','',dd1)
        print(dd1)
        input_name <- paste0('CHAP_',pollution,'_D1K_',dd1,'.nc')
        if(!file.exists(input_name)){input_name <- paste0('CHAP_',pollution,'_D10K_',dd1,'.nc')      }
        input_file <- nc_open(input_name)#, readunlim = FALSE, suppress_dimvals = TRUE)
        obs_lon <- ncvar_get(input_file, 'lon')
        obs_lat <- ncvar_get(input_file, 'lat')
        obs_var <- ncvar_get(input_file, pollution)
        nc_close(input_file)
        rm(input_file)
        gc()
        # Sys.sleep(2)
        # tar_lat <- seq(36.5,43.7,0.4)
        # tar_lon <- seq(112.4,117.6,0.4)
        tar_lat <- seq(37, 43.2, 0.4)
        tar_lon <- seq(110, 120, 0.5)
        # BEIJING: tar_lon <- c(116, 116.4); tar_lat <- c(39.3, 39.7)
        idx_lon <- which(obs_lon>=min(tar_lon) & obs_lon<max(tar_lon))
        idx_lat <- which(obs_lat>=min(tar_lat) & obs_lat<max(tar_lat))
        # obs_var <- obs_var[idx_lon, idx_lat]
        x <- c()
        for(i in 2:length(tar_lon))
          for(j in 2:length(tar_lat)){
            idx_lon <- which(obs_lon>=tar_lon[i-1] & obs_lon<tar_lon[i])
            idx_lat <- which(obs_lat>=tar_lat[j-1] & obs_lat<tar_lat[j])
            x1 <- as.vector(obs_var[idx_lon, idx_lat])
            x1 <- x1[!is.na(x1)]
            x <- c(x, mean(x1))
          }
        X_year <- rbind(X_year, x)
      }
      
      save(X_year, file = save_file)
      print(paste0('-----------------',pollution,'_',year))
      rm(obs_lon,obs_lat,obs_var,X_year)
      gc()
    }
    # Sys.sleep(2)
  }
  
  gc()
  # Sys.sleep(10)
  
}





