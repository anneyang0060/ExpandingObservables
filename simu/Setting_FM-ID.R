##### Setting : increasing dimension model
{
  L <- 5 # number of model changes
  Tmax <- 2e+05
  tau <- c(seq(2500,50000,2500), seq(55000,70000,5000), 8e4, 1e5, 1.4e5)
  K <- length(tau)
  eta <- c(2500,5000,12500,30000,80000)
  
  # dynamic dimension and sparsity
  M <- 5 # max No. penalties
  p0 <- 500
  delta_p <- rep(500, K) # increased dimensions
}

n_blk <- 50

# design and noise parameters
# rho <- 0.3
sigma_eps <- 1
tau <- c(0, tau)
eta <- c(0, eta)
p <- c(p0, p0 + cumsum(delta_p))
# p[K+1] <- 1e+04
s <- 5:10
s0 <- s[1]

# data-driven warm-up parameters
warm_grid_ratio <- 1.25
warm_min_blocks <- 2
delta_sigma <- 0.05
c_sigma <- 2.5
q_cov <- 0
R_sigma <- NULL
warm_const <- 0.01;5e-4

# hard-selection parameter
c_h <- 2

n_new <- 0 # number of obs of the new dimension
r <- 0.01

# true coefficient
idx <- 1:s0
for(i in 2:(K+1)){
  if(tau[i] %in% eta){
    idx <- c(idx, p[i-1]+1)
  }
}
beta_true <- rep(0, p[K+1])

set.seed(12)
beta_true[idx] <- sample(c((-10):(-1),1:10), max(s))
beta_true[idx] <- (beta_true[idx])[order(abs(beta_true[idx]),decreasing = TRUE)]

delta_p <- c(p0, delta_p)
m_dep <- 10
idx_dep <- as.vector(sapply(cumsum(delta_p)-500, function(i){(i+1):(i+m_dep)}))

block_width <- p0

beta_psue <- matrix(0, length(tau), max(p))

for (i in 1:(length(tau) - 1)) {
  
  p_cur <- p[i]
  
  idx_x <- idx_dep[idx_dep <= p_cur]
  idx_u <- idx_dep[idx_dep > p_cur]
  
  if (length(idx_x) > 0) {
    
    if (length(idx_u) > 0) {
      
      idx_joint <- c(idx_x, idx_u)
      
      Sigma_joint <- make_chain_crossblock_cov(
        idx_abs = idx_joint,
        rho = rho,
        block_width = block_width,
        m_dep = m_dep
      )
      
      p1 <- length(idx_x)
      q1 <- length(idx_u)
      
      SX <- Sigma_joint[seq_len(p1), seq_len(p1), drop = FALSE]
      SXU <- Sigma_joint[seq_len(p1), p1 + seq_len(q1), drop = FALSE]
      
      gamma <- beta_true[idx_u]
      
      alpha_psue <- beta_true[idx_x] +
        as.vector(ginv(SX) %*% SXU %*% gamma)
      
    } else {
      
      alpha_psue <- beta_true[idx_x]
      
    }
    
    beta_psue[i, idx_x] <- alpha_psue
  }
}

beta_psue[length(tau), ] <- beta_true

K_cv <- 2

prefix <- paste0('res/FM-ID_rho',rho,'_radius',radius_mult,'_hardcv')
