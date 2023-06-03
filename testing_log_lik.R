
w <- data.frame(summary(fit_m, par = c("alpha_mu_gene", "alpha_sigma_gene", 
                                       "beta_mu_gene", "beta_sigma_gene"))$summary)

model_log <- rstan::stan_model(file = "tmp.stan")

e_a <- (w$mean[1:20]+w$mean[21])-(w$mean[22:41]+w$mean[42])
e_a <- exp(e_a)/sum(exp(e_a))

e_b <- (w$mean[1:20]+w$mean[21])+(w$mean[22:41]+w$mean[42])
e_b <- exp(e_b)/sum(exp(e_b))

o <- rstan::sampling(object = model_log,
                     data = list(N = 20,
                                 y_a = u$Y_1[,2],
                                 y_b = u$Y_2[,2],
                                 a = e_a,
                                 b = e_b),
                     chains = 1,
                     cores = 1,
                     iter = 3500,
                     warmup = 1500,
                     algorithm = "Fixed_param")
summary(o)$summary

summary(fit_m, par = "log_lik_test")$summary
