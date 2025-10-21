set.seed(7)

residuals_all <- unlist(residuals_list)
debq_all <- unlist(lapply(subsamples, function(df1) df1$DEBQ))
resid_low <- residuals_all[debq_all <= 7]
resid_high <- residuals_all[debq_all > 7]


library(nleqslv)

estimate_param_resid_high <- estimate_pars_optim(resid_high)
params_high <- c(kappa = 1, lambda = estimate_param_resid_high[1], beta = estimate_param_resid_high[2])

estimate_param_resid_low <- estimate_pars_optim(resid_low)
params_low <- c(kappa = 1, lambda = estimate_param_resid_low[1],  beta = estimate_param_resid_low[2])


library(stats)

ks.test(resid_low, resid_high)
wilcox.test(resid_high, resid_low)
t.test(resid_low, resid_high)

# Проверка для параметра формы:
resid_high_modif <- resid_high / params_high["beta"]
resid_low_modif <- resid_low / params_low["beta"]

t.test(resid_low_modif, resid_high_modif)
wilcox.test(resid_high_modif, resid_low_modif)
ks.test(resid_high_modif, resid_low_modif)

# Проверка для параметра масштаба:
zz_low  <- resid_low * log(resid_low)  - mean(resid_low) * mean(log(resid_low)) 
zz_high <- resid_high * log(resid_high)  - mean(resid_high) * mean(log(resid_high))

t.test(zz_low, zz_high)   
ks.test(zz_low, zz_high)  
wilcox.test(zz_low, zz_high)







