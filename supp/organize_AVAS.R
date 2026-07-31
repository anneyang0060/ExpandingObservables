rm(list=ls())
# setwd('.../code')
sds <- 1:100

if (!dir.exists("fig")) dir.create("fig")

############################################################
## Helper functions
############################################################

get_common_est_time_avas <- function(rho, sds, base_dir = "res") {
  
  time_list <- list()
  
  for (sd in sds) {
    filename <- paste0(base_dir, "/AVAS_rho", rho, "/sd", sd, ".Rdata")
    
    if (!file.exists(filename)) {
      stop("Missing file: ", filename)
    }
    
    e <- new.env()
    load(filename, envir = e)
    
    if (!exists("est_time", envir = e)) {
      stop("est_time is not saved in ", filename)
    }
    
    time_list[[length(time_list) + 1]] <- as.numeric(e$est_time)
  }
  
  common_time <- Reduce(intersect, time_list)
  common_time <- sort(unique(common_time))
  
  if (length(common_time) == 0) {
    stop("No common est_time found for rho = ", rho)
  }
  
  return(common_time)
}

read_avas_fd_intersection <- function(rho,
                                      sds,
                                      base_dir = "res",
                                      enforce_dim_monotone = TRUE) {
  
  common_time <- get_common_est_time_avas(
    rho = rho,
    sds = sds,
    base_dir = base_dir
  )
  
  err_mat <- matrix(NA_real_, nrow = length(common_time), ncol = length(sds))
  pt_mat <- matrix(NA_real_, nrow = length(common_time), ncol = length(sds))
  time_mat <- matrix(NA_real_, nrow = length(common_time), ncol = length(sds))
  sigma_mat <- matrix(NA_real_, nrow = length(common_time), ncol = length(sds))
  iota0_vec <- rep(NA_real_, length(sds))
  
  for (j in seq_along(sds)) {
    
    sd <- sds[j]
    filename <- paste0(base_dir, "/AVAS_rho", rho, "/sd", sd, ".Rdata")
    
    e <- new.env()
    load(filename, envir = e)
    
    id <- match(common_time, e$est_time)
    
    if (any(is.na(id))) {
      stop("Internal error: common_time not found in ", filename)
    }
    
    err_mat[, j] <- e$err_beta_o[id]
    
    pt_vec <- e$pt_o[id]
    if (enforce_dim_monotone) {
      pt_vec <- cummin(pt_vec)
    }
    pt_mat[, j] <- pt_vec
    
    time_mat[, j] <- e$time_o[id]
    
    if (exists("sigma_o", envir = e)) {
      sigma_mat[, j] <- e$sigma_o[id]
    }
    
    if (exists("iota0_o", envir = e) && length(e$iota0_o) > 0) {
      iota0_vec[j] <- e$iota0_o[1]
    }
  }
  
  return(list(
    est_time = common_time,
    err = rowMeans(err_mat, na.rm = FALSE),
    pt = rowMeans(pt_mat, na.rm = FALSE),
    time = rowMeans(time_mat, na.rm = FALSE),
    sigma = rowMeans(sigma_mat, na.rm = TRUE),
    iota0 = iota0_vec
  ))
}

safe_log_range <- function(...) {
  x <- unlist(list(...))
  x <- x[is.finite(x) & x > 0]
  range(log(x), na.rm = TRUE)
}

safe_range <- function(...) {
  x <- unlist(list(...))
  x <- x[is.finite(x)]
  rg <- range(x, na.rm = TRUE)
  pad <- 0.05 * diff(rg)
  if (!is.finite(pad) || pad == 0) pad <- 0.01
  c(rg[1] - pad, rg[2] + pad)
}

add_line_safe <- function(x, y, col, lty = 1, lwd = 1) {
  n <- min(length(x), length(y))
  x <- x[seq_len(n)]
  y <- y[seq_len(n)]
  ok <- is.finite(x) & is.finite(y)
  lines(x[ok], y[ok], col = col, lty = lty, lwd = lwd)
}

############################################################
## Fixed dimension results
############################################################

for (rho in c(0, 0.3, 0.7)) {
  
  ##########################################################
  ## Read AVAS/RAVAS using intersection of est_time
  ##########################################################
  
  avas_out <- read_avas_fd_intersection(
    rho = rho,
    sds = sds,
    base_dir = "res",
    enforce_dim_monotone = TRUE
  )
  
  x_AVAS <- avas_out$est_time
  err_AVAS <- avas_out$err
  p_AVAS <- avas_out$pt
  time_AVAS <- avas_out$time
  
  cat("rho =", rho,
      "number of common AVAS est_time =", length(x_AVAS),
      "first =", min(x_AVAS),
      "last =", max(x_AVAS),
      "\n")
  
  cat("rho =", rho,
      "mean iota0 =", round(mean(avas_out$iota0, na.rm = TRUE), 1),
      "median iota0 =", round(median(avas_out$iota0, na.rm = TRUE), 1),
      "\n")
  
  ##########################################################
  ## Read baseline methods
  ##########################################################
  
  ## BR
  err <- c(); runing_times <- c(); dt <- c()
  for (sd in sds) {
    filename <- paste0("res/BR_fan_rho", rho, "/sd", sd, ".Rdata")
    load(filename)
    err <- cbind(err, err_beta_o)
    runing_times <- cbind(runing_times, time_o)
    dt <- cbind(dt, dt_o)
  }
  err_BR <- rowMeans(err)
  time_BR <- rowMeans(runing_times)
  p_BR <- rowMeans(dt)
  
  ## TSGD
  err <- c(); runing_times <- c(); dt <- c()
  for (sd in sds) {
    filename <- paste0("res/TSGD_rho", rho, "/sd", sd, ".Rdata")
    load(filename)
    err <- cbind(err, err_beta_o)
    runing_times <- cbind(runing_times, time_o)
    dt <- cbind(dt, dt_o)
  }
  err_TSGD <- rowMeans(err)
  time_TSGD <- rowMeans(runing_times)
  p_TSGD <- rowMeans(dt)
  
  ## OLL
  err <- c(); runing_times <- c(); dt <- c()
  for (sd in sds) {
    filename <- paste0("res/OLL_sun_rho", rho, "/sd", sd, ".Rdata")
    if (file.exists(filename)) {
      load(filename)
    } else {
      next
    }
    err <- cbind(err, err_beta_o)
    runing_times <- cbind(runing_times, time_o)
    dt <- cbind(dt, dt_o)
  }
  err_OLL <- rowMeans(err)
  time_OLL <- rowMeans(runing_times)
  p_OLL <- rowMeans(dt)
  
  ## ODL and OR
  err1 <- c(); err2 <- c(); err3 <- c()
  time1 <- c(); time2 <- c()
  dt1 <- c(); dt2 <- c()
  
  for (sd in sds) {
    filename <- paste0("res/ODL_huang_rho", rho, "/sd", sd, ".Rdata")
    load(filename)
    err1 <- cbind(err1, err_beta_b_o)
    err2 <- cbind(err2, err_beta_o)
    err3 <- cbind(err3, err_beta_g_o)
    time1 <- cbind(time1, time_o_u)
    time2 <- cbind(time2, time_o_bg)
    dt1 <- cbind(dt1, dt_o_b)
    dt2 <- cbind(dt2, dt_o_g)
  }
  
  err_ODL <- rowMeans(err1)
  err_OR <- rowMeans(err3)
  time_ODL <- rowMeans(time1)
  time_OR <- rowMeans(time2)
  p_ODL <- rowMeans(dt1)
  p_OR <- rowMeans(dt2)
  
  ##########################################################
  ## Fixed x axes for baseline methods
  ##########################################################
  
  x_BR_err <- seq(210, 10000, 10)
  x_TSGD_err <- seq(210, 10000, 10)
  x_BR_dim <- seq(200, 10000, 10)
  x_TSGD_dim <- seq(200, 10000, 10)
  
  x_BR_time <- seq(210, 10000, 10)
  x_TSGD_time <- seq(230, 10000, 10)
  
  x_ODL <- seq(300, 10000, 100)
  x_OLL <- seq(300, 10000, 100)
  x_OR <- seq(300, 10000, 100)
  
  ##########################################################
  ## Error figure
  ##########################################################
  
  {
    cc <- 1.6
    
    pdf(paste0("fig/err_fd_rho", rho, ".pdf"), height = 6, width = 6)
    
    ylim1 <- safe_log_range(
      err_AVAS,
      err_BR[2:length(err_BR)],
      err_TSGD[2:length(err_TSGD)],
      err_ODL,
      err_OLL,
      err_OR
    )
    
    ok_AVAS_err <- is.finite(err_AVAS) & err_AVAS > 0
    
    plot(x_AVAS[ok_AVAS_err],
         log(err_AVAS[ok_AVAS_err]),
         ylim = ylim1,
         main = paste0("rho=", rho),
         type = "l",
         col = "red",
         xlab = "data block",
         ylab = "log(MSE)",
         cex.axis = cc,
         cex.lab = cc,
         cex.main = cc)
    
    add_line_safe(x_BR_err, log(err_BR[2:length(err_BR)]), col = "blue")
    add_line_safe(x_TSGD_err, log(err_TSGD[2:length(err_TSGD)]), col = "green")
    add_line_safe(x_ODL, log(err_ODL), col = "black")
    add_line_safe(x_OR, log(err_OR), col = "orange")
    add_line_safe(x_OLL, log(err_OLL), col = "purple")
    
    if (rho == 0.7) {
      legend("topright",
             c("AVAS", "BR", "TSGD", "ODL", "OR", "OLL"),
             col = c("red", "blue", "green", "black", "orange", "purple"),
             lty = 1,
             cex = cc - 0.3)
    }
    
    dev.off()
  }
  
  ##########################################################
  ## Running time figure
  ## ODL is excluded because it includes debiasing.
  ## No artificial y-axis truncation or rescaling.
  ##########################################################
  
  {
    cc <- 1.6
    
    png(paste0("fig/time_fd_rho", rho, ".png"),
        height = 600 * 0.6,
        width = 600 * 0.6,
        res = 80 * 0.8)

    
    ylim_time <- c(0,0.1)
    
    ok_AVAS_time <- which(time_AVAS<0.1)
    
    plot(x_AVAS[ok_AVAS_time],
         time_AVAS[ok_AVAS_time],
         main = paste0("rho=", rho),
         xlim = c(50, 10100),
         ylim = ylim_time,
         xaxs = "i",
         yaxs = "i",
         type = "l",
         col = "red",
         xlab = "data block",
         ylab = "running time",
         cex.axis = cc,
         cex.lab = cc,
         cex.main = cc)
    
    add_line_safe(x_BR_time, time_BR[-1], col = "blue")
    add_line_safe(x_TSGD_time, time_TSGD[-c(1, 2, 3)], col = "green")
    add_line_safe(x_OLL, time_OLL, col = "purple")
    add_line_safe(seq(500, 10000, 100), time_OR[-c(1, 2)], col = "orange")
    
    if (rho == 0.7) {
      legend("topright",
             c("AVAS", "BR", "TSGD", "OLL", "OR"),
             col = c("red", "blue", "green", "purple", "orange"),
             lty = 1,
             cex = cc - 0.3)
    }
    
    dev.off()
  }
  
  ##########################################################
  ## Dimension figure
  ##########################################################
  
  {
    cc <- 1.6
    
    pdf(paste0("fig/dim_fd_rho", rho, ".pdf"), height = 6, width = 6)
    
    ok_AVAS_dim <- is.finite(p_AVAS)
    
    plot(x_AVAS[ok_AVAS_dim],
         p_AVAS[ok_AVAS_dim],
         main = paste0("rho=", rho),
         xlim = c(50, 10100),
         ylim = c(0, 1000),
         type = "l",
         col = "red",
         xlab = "data block",
         ylab = "No. selected variables",
         cex.axis = cc,
         cex.lab = cc,
         cex.main = cc)
    
    add_line_safe(x_BR_dim, p_BR, col = "blue")
    add_line_safe(x_TSGD_dim, p_TSGD, col = "green")
    add_line_safe(x_ODL, p_ODL, col = "black")
    add_line_safe(x_OLL, p_OLL, col = "purple")
    add_line_safe(x_OR, p_OR, col = "orange")
    
    if (rho == 0.7) {
      legend("topright",
             c("AVAS", "BR", "TSGD", "ODL", "OR", "OLL"),
             col = c("red", "blue", "green", "black", "orange", "purple"),
             lty = 1,
             cex = cc - 0.3)
    }
    
    dev.off()
  }
}
