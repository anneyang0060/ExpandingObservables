rm(list=ls())
setwd('.../CMFD')
# library(raster)
library(ncdf4)

# coordinate of PM2.5: 116.366,39.8673
vars <- c('wind','temp','srad','shum','pres','prec','lrad')

for(k in 1:7)
{
  var_name <- vars[k]

  for(year in 1991:2018)
  {
    input_name <- paste0(var_name,'_ITPCAS-CMFD_V0106_B-01_01dy_010deg_',year,'.nc')
    save_file <- paste0(var_name,'_',year,'.Rdata')

    X_year <- c()

    input_file <- nc_open(input_name)#, readunlim = FALSE, suppress_dimvals = TRUE)
    obs_lon <- ncvar_get(input_file, 'lon')
    obs_lat <- ncvar_get(input_file, 'lat')
    obs_var <- ncvar_get(input_file, var_name)
    nc_close(input_file)
    rm(input_file)
    gc()
    tar_lat <- seq(37, 43.2, 0.4)
    tar_lon <- seq(110, 120, 0.5)
    idx_lon <- which(obs_lon>=min(tar_lon) & obs_lon<max(tar_lon))
    idx_lat <- which(obs_lat>=min(tar_lat) & obs_lat<max(tar_lat))
    for(dd in 1:365){
      x <- c()
      for(i in 2:length(tar_lon))
        for(j in 2:length(tar_lat)){
          idx_lon <- which(obs_lon>=tar_lon[i-1] & obs_lon<tar_lon[i])
          idx_lat <- which(obs_lat>=tar_lat[j-1] & obs_lat<tar_lat[j])
          x1 <- as.vector(obs_var[idx_lon, idx_lat,dd])
          x1 <- x1[!is.na(x1)]
          x <- c(x, mean(x1))
        }
      X_year <- rbind(X_year, x)
    }
      

    save(X_year, file = save_file)
    print(paste0('-----------------',var_name,'_',year))
    rm(obs_lon,obs_lat,obs_var,X_year)
    gc()
    }

  gc()

}





