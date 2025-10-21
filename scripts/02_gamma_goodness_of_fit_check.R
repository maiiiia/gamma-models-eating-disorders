# Оценка параметров гамма-распределения ИМТ, проверка согласия 
# по критерию хи-квадрат

gammaDistrDens <- function(x, s, r, coef)
{
  return (dgamma(x, shape = s, rate = 1/r)*coef)
}

count_lambda <- function(l)
{
  m <- mean(l)
  v <- var(l)
  b <- v/m
  lambda <- m/b
  return (lambda)
}

count_b <- function(l)
{
  m <- mean(l)
  v <- var(l)
  b <- v/m
  return (b)
}


estimate_parameters <- function(x_sample, lambda_approx, b_approx)
{
  eq_system <- function(params, x, n) {
    
    f1 <- n*log(1/params[2]) + sum(log(x)) - n * digamma(params[1])
    f2 <- n*params[1] * params[2] - sum(x)
    
    return(c(f1, f2))
    
  }
  
  fn <- function(params)
  {
    return(eq_system(params, x = x_sample, n=length(x_sample)))
  }
  
  result <- nleqslv(c(lambda_approx, b_approx), fn, method = "Newton")
  return(result$x)
}

estimate_pars_optim <- function(x_sample)
{
  lambda <- count_lambda(x_sample)
  beta <- count_b(x_sample)
  par <- estimate_parameters(x_sample, lambda, beta)
  return(par)
}

library(zoo)

run_chisq_test_batch <- function(df, n_samples = 20, 
                                 sample_size = 200, seed = 7,
                                 m = 2, n_bins = 50) 
{
  set.seed(seed)
  
  sample_list<- lapply(1:n_samples, function(i) 
    sample(df, size = sample_size, replace = F))
  
  estimate_param_res <- lapply(1:n_samples, function(i) estimate_pars_optim(sample_list[[i]]))
  
  L <- sapply(1:n_samples, function(i) estimate_param_res[[i]][1])
  B <- sapply(1:n_samples, function(i) estimate_param_res[[i]][2])
  
  plot_list <- lapply(1:n_samples, function(i)
    hist(sample_list[[i]], breaks=n_bins))
  
  breaks_distrF_list <- lapply(1:n_samples, function(i)
    pgamma(plot_list[[i]]$breaks, shape=L[i], scale=B[i]))
  
  expected <- lapply(1:n_samples, function(i)
    rollapply(breaks_distrF_list[[i]], 2, function(x) x[2]-x[1]))
  
  chisq_test_result <- lapply(1:n_samples, function(i)
    chisq.test(plot_list[[i]]$counts, p=expected[[i]],
               rescale.p=TRUE, simulate.p.value=FALSE))
  
  
  df_chi <- (n_bins - 1) - m
  
  p_val <- lapply(1:n_samples, function(i) 
    pchisq(chisq_test_result[[i]]$statistic, 
           df = df_chi, lower.tail = F))
  
  return(list(
    p_values = p_val,
    passed = which(p_val > 0.05)
  ))
}

df1 <- data.frame(ИМТ = data_women$ИМТ, 
                  Возраст = data_women$Возраст, 
                  DEBQ = data_women$DEBQsumm)
library(tidyr)
df1 <- drop_na(df1)

result_chisq_women <- run_chisq_test_batch(df1$ИМТ)

result_chisq_DEBQ1 <- run_chisq_test_batch(data_DEBQ1$ИМТ)
result_chisq_DEBQ2 <- run_chisq_test_batch(data_DEBQ2$ИМТ)
result_chisq_DEBQ3 <- run_chisq_test_batch(data_DEBQ3$ИМТ)
result_chisq_DEBQ4 <- run_chisq_test_batch(data_DEBQ4$ИМТ)
result_chisq_DEBQ5 <- run_chisq_test_batch(data_DEBQ5$ИМТ)

result_chisq_age0 <- run_chisq_test_batch(data_age0$ИМТ)
result_chisq_age1 <- run_chisq_test_batch(data_age1$ИМТ)
result_chisq_age2 <- run_chisq_test_batch(data_age2$ИМТ)
result_chisq_age3 <- run_chisq_test_batch(data_age3$ИМТ)


