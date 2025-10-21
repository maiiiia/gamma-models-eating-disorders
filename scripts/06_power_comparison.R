
set.seed(1)

lambda1 <- params_low["lambda"]     
beta1   <- params_low["beta"]

lambda2 <- params_high["lambda"]    
beta2   <- params_high["beta"]

sample_sizes <- seq(250, 4000, by = 250)

n_sim <- 100

power_results <- data.frame()

for (n in sample_sizes) {
  beta_t_sig <- 0
  beta_w_sig <- 0
  beta_ks_sig <- 0
  
  for (i in 1:n_sim) {
    x1 <- rgamma(n, shape = lambda1, scale = beta1)
    x2 <- rgamma(n, shape = lambda2, scale = beta2)
    
    zz_1 <- x1 * log(x1) - mean(x1) * mean(log(x1))
    zz_2 <- x2 * log(x2) - mean(x2) * mean(log(x2))
    
    beta_t <- t.test(zz_1, zz_2)$p.value
    beta_w <- wilcox.test(zz_1, zz_2)$p.value
    beta_ks <- ks.test(zz_1, zz_2)$p.value
    
    beta_t_sig <- beta_t_sig + as.integer(beta_t < 0.05)
    beta_w_sig <- beta_w_sig + as.integer(beta_w < 0.05)
    beta_ks_sig <- beta_ks_sig + as.integer(beta_ks < 0.05)
  }
  
  power_results <- rbind(power_results, data.frame(
    n = n,
    
    power_beta_t = beta_t_sig / n_sim,
    power_beta_w = beta_w_sig / n_sim,
    power_beta_ks = beta_ks_sig / n_sim
  ))
}




sample_sizes <- seq(20, 80, by = 5)

power_results_lambda <- data.frame()

for (n in sample_sizes) {
  lam_t_sig <- 0
  lam_w_sig <- 0
  lam_ks_sig <- 0
  
  
  for (i in 1:n_sim) {
    x1 <- rgamma(n, shape = lambda1, scale = beta1)
    x2 <- rgamma(n, shape = lambda2, scale = beta2)
    
    lam_t <- t.test(x1 / beta1, x2 / beta2)$p.value
    lam_w <- wilcox.test(x1 / beta1, x2 / beta2)$p.value
    lam_ks <- ks.test(x1 / beta1, x2 / beta2)$p.value
    
    lam_t_sig <- lam_t_sig + as.integer(lam_t < 0.05)
    lam_w_sig <- lam_w_sig + as.integer(lam_w < 0.05)
    lam_ks_sig <- lam_ks_sig + as.integer(lam_ks < 0.05)
    
  }
  
  power_results_lambda <- rbind(power_results_lambda, data.frame(
    n = n,
    power_lam_t = lam_t_sig / n_sim,
    power_lam_w = lam_w_sig / n_sim,
    power_lam_ks = lam_ks_sig / n_sim
    
  ))
}

