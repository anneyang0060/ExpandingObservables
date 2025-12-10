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
n0 <- 300
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

beta_psue <- matrix(0, length(tau), max(p))
eff_coeff <- beta_true[idx_dep]
for(i in 1:(length(tau)-1)){
  p_cur <- p[i]
  p1 <- sum(idx_dep <= p_cur)
  q1 <- sum(idx_dep > p_cur)
  SX <- matrix(rho,p1,p1)
  diag(SX) <- 1
  SXU <- matrix(rho, p1, q1)
  gamma <- eff_coeff[(p1+1):(p1+q1)]
  alpha_psue <- as.vector(eff_coeff[1:p1] + ginv(SX) %*% SXU %*% gamma)
  beta_psue[i, idx_dep[1:p1]] <- alpha_psue
  # print(sum((alpha_psue-eff_coeff[1:p1])^2))
}
beta_psue[length(tau), ] <- beta_true
K_cv <- 2

prefix <- paste0('res/FM-ID_rho',rho)
