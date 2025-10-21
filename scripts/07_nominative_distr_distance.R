library(nleqslv)

set.seed(7)

sample_sizes <- seq(100, 50000, by = 300)
results_nom_gamma <- data.frame()

lambda1 <- params_low["lambda"]     
beta1   <- params_low["beta"]

kappa_test <- seq(0.75, 1, by = 0.01)

for (n in sample_sizes) {
  x <- rgamma(n, shape = lambda1, scale = beta1)
  
  est <- estimate_pars_optim(x)
  par_est <- c(kappa = 1, lambda = est[1], beta = est[2])
  
  syn_list_temp <- lapply(kappa_test, function(k2) {
    
    par2 <- find_synonimous(par_est, k2)
    H22 <- H_12(par2, par2)
    list(kappa2 = k2,  par2 = par2, H22 = H22)
    
  })
  
  syn_table_temp <- do.call(rbind, lapply(syn_list_temp, function(x) {
    data.frame(kappa2 = x$kappa2,  beta = x$par2[2], lambda = x$par2[3], H22 = x$H22)
  }))
  
  syn_table_temp <- syn_table_temp[order(syn_table_temp$H22), ]
  syn_table_temp$nominee <- syn_table_temp$H22 == min(syn_table_temp$H22)
  
  par_nom <- c(kappa = syn_table_temp[1,1], lambda = syn_table_temp[1,3], beta = syn_table_temp[1,2])
  
  J_val <- H_12(par_est, par_nom) - H_12(par_est, par_est) +
    H_12(par_nom, par_est) - syn_table_temp[1,4]
  
  results_nom_gamma <- rbind(results_nom_gamma, data.frame(n = n,
                                                           J = J_val,
                                                           lambda_hat = est[1],
                                                           beta_hat = est[2],
                                                           
                                                           kappa_nom = par_nom["kappa"],
                                                           lambda_nom = par_nom["lambda"],
                                                           beta_nom = par_nom["beta"]
  ))
  
}


