df2 <- data.frame(ИМТ = data_women$ИМТ, 
                  Возраст = data_women$Возраст, 
                  DEBQ = data_women$DEBQsumm, 
                  EDE = data_women$EDEРезультат, 
                  TAS = data_women$TASРезультат)

df2 <- drop_na(df2)

set.seed(7)
n_subsamples <- 20
subsample_size <- 200


subsamples <- lapply(1:n_subsamples, function(i) {
  df_i <- df2[sample(nrow(df2), subsample_size), ]
  
  preds <- names(df_i) != "ИМТ"
  df_i[ , preds] <- scale(df_i[ , preds], center = TRUE, scale = FALSE)
  
  df_i
})


fit_gamma_model_vars <- function(df, predictors, beta_scale = 0.6) {
  Y <- df$ИМТ
  X <- as.matrix(df[, predictors, drop = FALSE])
  
  M <- mean(Y)
  V <- var(Y)
  alpha0 <- M / V
  lambda0 <- M * alpha0
  
  beta0 <- coef(lm(as.formula(paste("ИМТ ~", paste(predictors, collapse = " + "), "-1")), data = df)) * beta_scale
  
  
  logLikGamma <- function(theta) {
    alpha <- theta[1]
    lambda <- theta[2]
    beta <- theta[3:length(theta)]
    mu <- X %*% beta
    residuals <- Y - mu
    if (any(residuals <= 0 | !is.finite(residuals))) return(1e6)
    
    ll <- -sum(dgamma(residuals, shape=lambda, rate=alpha, log=TRUE))
    
    if (!is.finite(ll)) return(1e6)
    return(ll)
  }
  
  theta_start <- c(alpha0, lambda0, beta0)
  
  fit <- tryCatch(
    optim(par = theta_start, fn = logLikGamma, method = "Nelder-Mead", control = list(maxit = 1000)),
    error = function(e) NULL
  )
  
  result <- c(alpha = fit$par[1], lambda = fit$par[2], fit$par[3:length(fit$par)], logLik = fit$value)
  names(result)[3:(length(result)-1)] <- predictors
  return(result)
  
}


get_gamma_residuals_generic <- function(df, predictors, beta_vec) {
  Y <- df$ИМТ
  X <- as.matrix(df[, predictors, drop = FALSE])
  mu <- X %*% beta_vec
  residuals <- Y - mu
  return(residuals)
}

extract_betas <- function(res, predictors) as.numeric(res[predictors])

predictors4 <- c("Возраст", "DEBQ", "EDE", "TAS")
results_4 <- lapply(subsamples, fit_gamma_model_vars, predictors = predictors4,  beta_scale = 0.8)

betas_list4 <- lapply(results_4, extract_betas, predictors = predictors4)

residuals_list_4 <- mapply(get_gamma_residuals_generic, subsamples, 
                           MoreArgs = list(predictors = predictors4), 
                           beta_vec = betas_list4, SIMPLIFY = FALSE)
lambda_vec4 <- sapply(results_4, function(res) res[["lambda"]])
alpha_vec4  <- sapply(results_4, function(res) res[["alpha"]])
beta_vec4 <- 1/alpha_vec4


compare_models_lambda_ttest <- function(resid_1, resid_2, beta_1, beta_2) {
  
  resid_1_m <- resid_1 / beta_1
  resid_2_m  <- resid_2 / beta_2
  
  t.test(resid_1_m, resid_2_m)$p.value
}

compare_beta_ks <- function(resid_1, resid_2) {
  zz1 <- resid_1 * log(resid_1) - mean(resid_1) * mean(log(resid_1))
  zz2 <- resid_2 * log(resid_2) - mean(resid_2) * mean(log(resid_2))
  ks.test(zz1, zz2)$p.value
}

predictors3 <- c("EDE", "TAS", "Возраст")
results_3 <- lapply(subsamples, fit_gamma_model_vars, predictors = predictors3,  beta_scale = 0.8)

betas_list3 <- lapply(results_3, extract_betas, predictors = predictors3)

residuals_list_3 <- mapply(get_gamma_residuals_generic, subsamples, 
                           MoreArgs = list(predictors = predictors3), 
                           beta_vec = betas_list3, SIMPLIFY = FALSE)


lambda_vec3 <- sapply(results_3, function(res) res[["lambda"]])
alpha_vec3  <- sapply(results_3, function(res) res[["alpha"]])
beta_vec3 <- 1/alpha_vec3

p_vals_lambda_comp <- mapply(
  compare_models_lambda_ttest,
  resid_1 = residuals_list_4,
  resid_2  = residuals_list_3,
  beta_1 = beta_vec3,
  beta_2  = beta_vec4
)

p_vals_beta_comp <- mapply(
  compare_beta_ks,
  resid_1 = residuals_list_4,
  resid_2 = residuals_list_3
)

sum(p_vals_lambda_comp <0.05) 
sum(p_vals_beta_comp <0.05)	


corr_age_resid <- mapply(
  function(resid, df) cor(resid, df$Возраст),
  resid   = residuals_list_4,
  df      = subsamples
)
 



