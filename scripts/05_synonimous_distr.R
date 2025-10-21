set.seed(7)
library(pracma)  

E_lnX <- function(kappa, lambda, alpha) {
  return((digamma(lambda) - log(alpha)) / kappa)
}

E_Xk <- function(k1, lambda1, alpha1, k2) {
  theta <- k2 / k1
  return(gamma(lambda1 + theta) / (gamma(lambda1) * alpha1^theta))
}

H_12 <- function(par1, par2) {
  k1 <- par1["kappa"]; beta1 <- par1["beta"]; lambda1 <- par1["lambda"]; 
  alpha1 <- 1 / beta1
  
  k2 <- par2["kappa"]; beta2 <- par2["beta"]; lambda2 <- par2["lambda"]; 
  alpha2 <- 1 / beta2
  
  elogx <- E_lnX(k1, lambda1, alpha1)
  exk <- E_Xk(k1, lambda1, alpha1, k2)
  
  h12 <- - (log(k2) + lambda2 * log(alpha2) - lgamma(lambda2)) -
    (k2 * lambda2 - 1) * elogx + alpha2 * exk
  return(h12)
}


H12 <- H_12(params_low, params_high)
H11 <- H_12(params_low, params_low)
delta1 <- H12 - H11

H21 <- H_12(params_high, params_low)
H22 <- H_12(params_high, params_high)
delta2 <- H21 - H22
J_low_high <- delta1+delta2

find_synonimous <- function(par1, kappa2) {
  kappa1 <- unname(par1["kappa"])
  beta1  <- unname(par1["beta"])
  lambda1 <- unname(par1["lambda"])
  
  alpha1 <- 1 / beta1
  theta <- kappa2 / kappa1
  
  psi_term <- digamma(lambda1 + theta) - digamma(lambda1)
  lambda2 <- 1 / (theta * psi_term)
  
  alpha2 <- lambda2 * alpha1^theta * gamma(lambda1) / gamma(lambda1 + theta)
  beta2 <- 1 / alpha2
  
  return(c(kappa = kappa2, beta = beta2, lambda = lambda2))
}

kappa_test <- seq(0.75, 1, by = 0.01)

syn_list_low <- lapply(kappa_test, function(k2) {
  par2 <- find_synonimous(params_low, k2)
  delta <- H_12(params_low, par2) - H_12(params_low, params_low)
  H22 <- H_12(par2, par2)
  list(kappa2 = k2, delta = unname(delta), par2 = par2, H22 = H22)
})

syn_table_low <- do.call(rbind, lapply(syn_list_low, function(x) {
  data.frame(kappa2 = x$kappa2, delta = x$delta, beta = x$par2[2], lambda = x$par2[3], H22 = x$H22)
}))

syn_table_low <- syn_table_low[order(syn_table_low$H22), ]
syn_table_low$nominee <- syn_table_low$H22 == min(syn_table_low$H22)

syn_list_high <- lapply(kappa_test, function(k2) {
  par2 <- find_synonimous(params_high, k2)
  delta <- H_12(params_high, par2) - H_12(params_high, params_high)
  H22 <- H_12(par2, par2)
  list(kappa2 = k2, delta = unname(delta), par2 = par2, H22 = H22)
})

syn_table_high <- do.call(rbind, lapply(syn_list_high, function(x) {
  data.frame(kappa2 = x$kappa2, delta = x$delta, beta = x$par2[2], lambda = x$par2[3], H22 = x$H22)
}))

syn_table_high <- syn_table_high[order(syn_table_high$H22), ]
syn_table_high$nominee <- syn_table_high$H22 == min(syn_table_high$H22)


lambda_low_nom <- syn_table_low[1, 4]
beta_low_nom   <- syn_table_low[1, 3]
kappa_low_nom <- syn_table_low[1, 1]
params_low_nom <- c(kappa = kappa_low_nom, lambda = lambda_low_nom, beta = beta_low_nom)

lambda_high_nom <- syn_table_high[1, 4]
beta_high_nom   <- syn_table_high[1, 3]
kappa_high_nom <- syn_table_high[1, 1]
params_high_nom <- c(kappa = kappa_high_nom, lambda = lambda_high_nom, beta = beta_high_nom)


J_low <- syn_table_low[1, 2] + H_12(params_low_nom, params_low) - H_12(params_low_nom, params_low_nom)
J_high <- syn_table_high[1, 2] + H_12(params_high_nom, params_high) - H_12(params_high_nom, params_high_nom)


# Проверка для параметра масштаба:
zz_low_nom  <- resid_low^kappa_low_nom * log(resid_low^kappa_low_nom) - 
  mean(resid_low^kappa_low_nom) * mean(log(resid_low^kappa_low_nom))

zz_high_nom <- resid_high^kappa_high_nom * log(resid_high^kappa_high_nom) - 
  mean(resid_high^kappa_high_nom) * mean(log(resid_high^kappa_high_nom))

t.test(zz_low_nom, zz_high_nom)
ks.test(zz_low_nom, zz_high_nom)
wilcox.test(zz_low_nom, zz_high_nom)


# Проверка для параметра формы:
x_low_norm  <- resid_low^kappa_low_nom / beta_low_nom
x_high_norm <- resid_high^kappa_high_nom / beta_high_nom

t.test(x_low_norm, x_high_norm)
ks.test(x_low_norm, x_high_norm)
wilcox.test(x_low_norm, x_high_norm)






