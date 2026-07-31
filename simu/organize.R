rm(list=ls())
# setwd('.../code')

############################### increasing dimension results #############################

tau <- c(0, seq(2500, 50000, 2500),
         seq(55000, 70000, 5000), 8e4, 1e5, 1.4e5, 1e7)

rhos <- c(0, 0.3, 0.5, 0.7)
sds <- 1:100

n_blk_id <- 50
Tmax_id <- 2e5

## In the simulation, each feature increment adds 500 variables.
delta_p_const <- 500
delta_p_by_cycle <- rep(delta_p_const, length(tau) - 1)

## Soft/hard transition used in the algorithm:
## old code used t_eff <= 2 * length(I)
c_h <- 2

if(!dir.exists("fig")) {
  dir.create("fig")
}

############################################################
## Helper functions
############################################################

cycle_id_alg <- function(time, tau) {
  ## Match the algorithm definition:
  ## k_tau = max(which(tau < t)).
  ## Therefore t = tau[k] still belongs to the previous cycle.
  id <- findInterval(time - 1e-8, tau)
  id <- pmax(id, 1)
  id <- pmin(id, length(tau) - 1)
  return(id)
}

align_locf <- function(time, value, grid_time) {
  ## Last observation carried forward.
  time <- as.numeric(time)
  value <- as.numeric(value)
  
  ord <- order(time)
  time <- time[ord]
  value <- value[ord]
  
  keep <- !duplicated(time, fromLast = TRUE)
  time <- time[keep]
  value <- value[keep]
  
  id <- findInterval(grid_time, time)
  out <- rep(NA_real_, length(grid_time))
  ok <- id > 0
  out[ok] <- value[id[ok]]
  return(out)
}

align_dim_cycle_fill_alg <- function(time, value, grid_time, tau) {
  ## Dimension is a maintained state. During warm-up, no permanent deletion occurs.
  ## So before the first recorded point within a cycle, fill by the first recorded
  ## dimension of that cycle.
  
  time <- as.numeric(time)
  value <- as.numeric(value)
  
  ord <- order(time)
  time <- time[ord]
  value <- value[ord]
  
  keep <- !duplicated(time, fromLast = TRUE)
  time <- time[keep]
  value <- value[keep]
  
  out <- rep(NA_real_, length(grid_time))
  
  obs_cycle <- cycle_id_alg(time, tau)
  grid_cycle <- cycle_id_alg(grid_time, tau)
  
  for(cc in intersect(unique(obs_cycle), unique(grid_cycle))) {
    
    id_obs <- which(obs_cycle == cc)
    id_grid <- which(grid_cycle == cc)
    
    if(length(id_obs) == 0 || length(id_grid) == 0) next
    
    tt <- time[id_obs]
    vv <- value[id_obs]
    
    oo <- order(tt)
    tt <- tt[oo]
    vv <- vv[oo]
    
    gg <- grid_time[id_grid]
    loc <- findInterval(gg, tt)
    
    tmp <- rep(NA_real_, length(gg))
    
    ## After first observed point: ordinary LOCF.
    ok_after <- loc > 0
    tmp[ok_after] <- vv[loc[ok_after]]
    
    ## Before first observed point in the same cycle:
    ## warm-up stage keeps candidate dimension unchanged.
    ok_before <- loc == 0
    tmp[ok_before] <- vv[1]
    
    out[id_grid] <- tmp
  }
  
  return(out)
}

align_bin_mean <- function(time, value, grid_time) {
  ## For running time: average all actual recorded times falling into
  ## the bin represented by each grid point.
  
  time <- as.numeric(time)
  value <- as.numeric(value)
  
  ok <- is.finite(time) & is.finite(value)
  time <- time[ok]
  value <- value[ok]
  
  if(length(time) == 0) {
    return(rep(NA_real_, length(grid_time)))
  }
  
  ord <- order(time)
  time <- time[ord]
  value <- value[ord]
  
  mid <- (head(grid_time, -1) + tail(grid_time, -1)) / 2
  breaks <- c(-Inf, mid, Inf)
  
  bin_id <- cut(
    time,
    breaks = breaks,
    include.lowest = TRUE,
    labels = FALSE
  )
  
  out <- rep(NA_real_, length(grid_time))
  
  for(g in seq_along(grid_time)) {
    vv <- value[bin_id == g]
    if(length(vv) > 0) {
      out[g] <- mean(vv, na.rm = TRUE)
    }
  }
  
  return(out)
}

row_mean_with_min_frac <- function(mat, min_frac = 0.8) {
  n_eff <- rowSums(!is.na(mat))
  out <- rowMeans(mat, na.rm = TRUE)
  out[is.nan(out)] <- NA_real_
  out[n_eff < ceiling(min_frac * ncol(mat))] <- NA_real_
  return(out)
}

enforce_cycle_monotone <- function(y, grid_time, tau) {
  ## This is only for displaying d_t.
  ## Within a cycle, candidate dimension should be non-increasing.
  
  out <- y
  grid_cycle <- cycle_id_alg(grid_time, tau)
  
  for(cc in unique(grid_cycle)) {
    id <- which(grid_cycle == cc & !is.na(out))
    if(length(id) > 1) {
      out[id] <- cummin(out[id])
    }
  }
  
  return(out)
}

############################################################
## Construct stage-aware grid_time
############################################################

make_stage_dense_grid_time <- function(rhos,
                                       sds,
                                       tau,
                                       delta_p_by_cycle,
                                       n_blk_id = 50,
                                       Tmax_id = 2e5,
                                       c_h = 2,
                                       base_step_blocks = 20,
                                       soft_extra_blocks = 5,
                                       pt_threshold = 200,
                                       high_buffer_blocks = 2,
                                       base_dir = "res") {
  
  ## 1. Sparse grid for global trend.
  base_grid <- seq(n_blk_id, Tmax_id, by = base_step_blocks * n_blk_id)
  
  ## 2. Feature-increment times.
  tau_grid <- tau[tau > 0 & tau <= Tmax_id]
  
  ## 3. Dense deterministic grid within theoretical soft-selection windows.
  soft_grid <- c()
  K_eff <- min(length(delta_p_by_cycle), length(tau) - 1)
  
  for(k in seq_len(K_eff)) {
    
    start_k <- tau[k] + n_blk_id
    
    soft_end_k <- tau[k] +
      c_h * delta_p_by_cycle[k] +
      soft_extra_blocks * n_blk_id
    
    end_k <- min(
      tau[k + 1] - n_blk_id,
      soft_end_k,
      Tmax_id
    )
    
    if(end_k >= start_k) {
      soft_grid <- c(
        soft_grid,
        seq(start_k, end_k, by = n_blk_id)
      )
    }
  }
  
  ## 4. Add actual high-dimensional records and peak-dimension records.
  high_times <- c()
  peak_times <- c()
  
  for(rho in rhos) {
    for(sd in sds) {
      
      filename <- paste0(base_dir, "/FM-ID_rho", rho, 'radius', 0.1, "/sd", sd, ".Rdata")
      if(!file.exists(filename)) next
      
      e <- new.env()
      load(filename, envir = e)
      
      if(!all(c("est_time", "pt_o") %in% ls(e))) next
      
      time_vec <- as.numeric(e$est_time)
      pt_vec <- as.numeric(e$pt_o)
      
      ok <- is.finite(time_vec) & is.finite(pt_vec)
      time_vec <- time_vec[ok]
      pt_vec <- pt_vec[ok]
      
      if(length(time_vec) == 0) next
      
      cyc_id <- cycle_id_alg(time_vec, tau)
      rel_time <- time_vec - tau[cyc_id]
      
      soft_limit <- c_h * delta_p_by_cycle[cyc_id] +
        soft_extra_blocks * n_blk_id
      
      ## Include all actual records in soft-selection windows.
      high_times <- c(high_times, time_vec[rel_time > 0 & rel_time <= soft_limit])
      
      ## Include all records with large selected dimension.
      high_times <- c(high_times, time_vec[pt_vec > pt_threshold])
      
      ## Include peak dimension time in each cycle.
      for(cyc in unique(cyc_id)) {
        idx_cyc <- which(cyc_id == cyc)
        if(length(idx_cyc) > 0) {
          idx_peak <- idx_cyc[which.max(pt_vec[idx_cyc])]
          peak_times <- c(peak_times, time_vec[idx_peak])
        }
      }
    }
  }
  
  key_times <- sort(unique(c(high_times, peak_times)))
  key_times <- key_times[is.finite(key_times)]
  
  if(length(key_times) > 0) {
    offsets <- seq(-high_buffer_blocks, high_buffer_blocks) * n_blk_id
    key_grid <- sort(unique(as.vector(
      sapply(key_times, function(tt) tt + offsets)
    )))
  } else {
    key_grid <- c()
  }
  
  key_grid <- key_grid[key_grid >= n_blk_id & key_grid <= Tmax_id]
  
  grid_time <- sort(unique(c(
    base_grid,
    tau_grid,
    soft_grid,
    key_grid,
    Tmax_id
  )))
  
  grid_time <- grid_time[grid_time >= n_blk_id & grid_time <= Tmax_id]
  
  info <- data.frame(
    n_base_grid = length(base_grid),
    n_tau_grid = length(tau_grid),
    n_soft_grid = length(unique(soft_grid)),
    n_key_grid = length(unique(key_grid)),
    n_final_grid = length(grid_time)
  )
  
  return(list(
    grid_time = grid_time,
    info = info,
    soft_grid = sort(unique(soft_grid)),
    key_grid = sort(unique(key_grid))
  ))
}

grid_obj <- make_stage_dense_grid_time(
  rhos = rhos,
  sds = sds,
  tau = tau,
  delta_p_by_cycle = delta_p_by_cycle,
  n_blk_id = n_blk_id,
  Tmax_id = Tmax_id,
  c_h = c_h,
  base_step_blocks = 20,
  soft_extra_blocks = 5,
  pt_threshold = 200,
  high_buffer_blocks = 2,
  base_dir = "res"
)

grid_time <- grid_obj$grid_time
print(grid_obj$info)

############################################################
## Organize results on the common grid_time
############################################################

err_full <- matrix(NA_real_, nrow = length(rhos), ncol = length(grid_time))
pt_full <- matrix(NA_real_, nrow = length(rhos), ncol = length(grid_time))
time_full <- matrix(NA_real_, nrow = length(rhos), ncol = length(grid_time))

n_eff_err <- matrix(0, nrow = length(rhos), ncol = length(grid_time))
n_eff_pt <- matrix(0, nrow = length(rhos), ncol = length(grid_time))
n_eff_time <- matrix(0, nrow = length(rhos), ncol = length(grid_time))

for(i_rho in seq_along(rhos)) {
  
  rho <- rhos[i_rho]
  
  err_mat <- matrix(NA_real_, nrow = length(grid_time), ncol = length(sds))
  pt_mat <- matrix(NA_real_, nrow = length(grid_time), ncol = length(sds))
  time_mat <- matrix(NA_real_, nrow = length(grid_time), ncol = length(sds))
  
  for(j_sd in seq_along(sds)) {
    
    sd <- sds[j_sd]
    filename <- paste0("res/FM-ID_rho", rho, "_radius1_hardcv/sd", sd, ".Rdata")
    
    if(!file.exists(filename)) {
      print(paste0("Missing file: ", filename))
      next
    }
    
    e <- new.env()
    load(filename, envir = e)
    
    if(exists("miss", envir = e)) {
      if(max(e$miss, na.rm = TRUE) > 0) {
        print(filename)
      }
    }
    
    ## Error: use current available estimate.
    err_mat[, j_sd] <- align_locf(
      time = e$est_time,
      value = e$err_psue_o,
      grid_time = grid_time
    )
    
    ## Dimension: use cycle-wise fill to avoid artificial spikes.
    ## If you later save pt_after_o, replace e$pt_o by e$pt_after_o here.
    pt_mat[, j_sd] <- align_dim_cycle_fill_alg(
      time = e$est_time,
      value = e$pt_o,
      grid_time = grid_time,
      tau = tau
    )
    
    ## Running time: use bin average.
    time_mat[, j_sd] <- align_bin_mean(
      time = e$est_time,
      value = e$time_o,
      grid_time = grid_time
    )
  }
  
  err_full[i_rho, ] <- row_mean_with_min_frac(err_mat, min_frac = 0.8)
  pt_full[i_rho, ] <- row_mean_with_min_frac(pt_mat, min_frac = 0.8)
  
  ## Running time is noisier and not a state variable, so use a lower threshold.
  time_full[i_rho, ] <- row_mean_with_min_frac(time_mat, min_frac = 0.3)
  
  n_eff_err[i_rho, ] <- rowSums(!is.na(err_mat))
  n_eff_pt[i_rho, ] <- rowSums(!is.na(pt_mat))
  n_eff_time[i_rho, ] <- rowSums(!is.na(time_mat))
}

## Optional monotone display correction for dimension curves.
pt_full_plot <- pt_full
for(i in seq_len(nrow(pt_full_plot))) {
  pt_full_plot[i, ] <- enforce_cycle_monotone(
    y = pt_full_plot[i, ],
    grid_time = grid_time,
    tau = tau
  )
}

############################################################
## Plot 1: estimation error
############################################################

cc <- 1.5

pdf("fig/err_id_radius1_hadrcv.pdf", height = 9, width = 13)
par(mfrow = c(4, 1), mai = c(0.6, 0.6, 0.2, 0.1))

for(i in seq_along(rhos)) {
  
  ok <- is.finite(err_full[i, ])
  
  plot(
    grid_time[ok] / n_blk_id,
    err_full[i, ok],
    type = "l",
    main = paste0("rho=", rhos[i]),
    xlab = "data block",
    ylab = "error",
    cex.axis = cc,
    cex.lab = cc,
    cex.main = cc
  )
  
  abline(
    v = tau[tau > 0 & tau <= Tmax_id] / n_blk_id,
    lty = 3,
    col = "gray"
  )
}

dev.off()

############################################################
## Plot 2: running time
############################################################

## Remove grid points with no running-time value in any panel.
keep_time <- colSums(is.na(time_full)) == 0

pdf("fig/time_id_radius1_hadrcv.pdf", height = 9, width = 13)
par(mfrow = c(4, 1), mai = c(0.6, 0.6, 0.2, 0.1))

for(i in seq_along(rhos)) {
  
  ok <- keep_time & is.finite(time_full[i, ])
  
  plot(
    grid_time[ok] / n_blk_id,
    time_full[i, ok],
    type = "l",
    main = paste0("rho=", rhos[i]),
    xlab = "data block",
    ylab = "time",
    cex.axis = cc,
    cex.lab = cc,
    cex.main = cc
  )
  
  abline(
    v = tau[tau > 0 & tau <= Tmax_id] / n_blk_id,
    lty = 3,
    col = "gray"
  )
}

dev.off()

############################################################
## Plot 3: selected dimension
############################################################

pdf("fig/dim_id_radius1_hadrcv.pdf", height = 9, width = 13)
par(mfrow = c(4, 1), mai = c(0.6, 0.6, 0.2, 0.1))

for(i in seq_along(rhos)) {
  
  ok <- is.finite(pt_full_plot[i, ])
  
  plot(
    grid_time[ok] / n_blk_id,
    pt_full_plot[i, ok],
    type = "l",
    main = paste0("rho=", rhos[i]),
    xlab = "data block",
    ylab = "d_t",
    cex.axis = cc,
    cex.lab = cc,
    cex.main = cc
  )
  
  abline(
    v = tau[tau > 0 & tau <= Tmax_id] / n_blk_id,
    lty = 3,
    col = "gray"
  )
}

dev.off()

