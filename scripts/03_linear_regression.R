
library(tidyr)
set.seed(7)

df1 <- data.frame(ИМТ = data_women$ИМТ, 
                  Возраст = data_women$Возраст, 
                  DEBQ = data_women$DEBQsumm)
df1 <- drop_na(df1)

n_samples <- 20
sample_size <- 200

subsamples <- lapply(1:n_samples, function(i) {
  df1[sample(1:nrow(df1), sample_size), ]
})


fit_glm <- function(data) {
  lm(ИМТ ~ Возраст + DEBQ, data = data)
}

glm_results <- lapply(subsamples, fit_glm)
glm_coefs <- lapply(glm_results, coef)
lapply(glm_results, summary)


set.seed(7)

library(nleqslv)
library(zoo)

fit_glm_age <- function(data) {
  lm(ИМТ ~ Возраст, data = data)
}

glm_results_age <- lapply(subsamples, fit_glm_age)
glm_coefs_age <- lapply(glm_results_age, coef)

shift_resid <- function(residuals)
{
  min_residual <- min(residuals)
  if (min_residual <= 0) {
    residuals_shifted <- residuals + abs(min_residual) + 0.0001
  } else {
    residuals_shifted <- residuals
  }
  return(residuals_shifted)
}

residuals_list <- lapply(glm_results_age, function(model) 
  residuals(model) |> as.vector() |> shift_resid())

estimate_param_res <- lapply(1:n_samples, function(i) 
  estimate_pars_optim(residuals_list[[i]]))

L <- sapply(1:n_samples, function(i) estimate_param_res[[i]][1])
B <- sapply(1:n_samples, function(i) estimate_param_res[[i]][2])

plot_list <- lapply(1:n_samples, function(i)
  hist(residuals_list[[i]], breaks=40))

breaks_distrF_list <- lapply(1:n_samples, function(i)
  pgamma(plot_list[[i]]$breaks, shape=L[i], scale=B[i]))

expected <- lapply(1:n_samples, function(i)
  rollapply(breaks_distrF_list[[i]], 2, function(x) x[2]-x[1]))

chisq_test_result <- lapply(1:n_samples, function(i)
  chisq.test(plot_list[[i]]$counts, p=expected[[i]],
             rescale.p=TRUE, simulate.p.value=FALSE))

m <-  2
n_bins <-  40
df_chi <- (n_bins - 1) - m

p_val_resid <- lapply(1:n_samples, function(i) 
  pchisq(chisq_test_result[[i]]$statistic, 
         df = df_chi, lower.tail = FALSE))

p_val_resid > 0.05









