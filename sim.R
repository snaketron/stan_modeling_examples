require(rstan)
m <- rstan::stan_model(file = "dm_sim/dm_sim.stan")

set.seed(seed = 12345)
k <- 20
a <- rexp(n = k, rate = 0.1)
a <- rnorm(n = k, mean = 0, sd = 1)
b <- rnorm(n = k, mean = 0, sd = 0.5)
xi <- 100
n <- 10^4
plot(exp(a+b)/sum(exp(a+b)))

sim <- rstan::sampling(object = m,
                       data = list(k = k, a = a, b = b, xi = xi, n = n),
                       chain = 4,
                       cores = 4,
                       iter = 500,
                       warmup = 100,
                       algorithm = "Fixed_param",
                       seed = 12345)

plot(sim, par = "y_1")
plot(sim, par = "y_2")

y_1_summary <- summary(sim, par = "y_1")$summary
y_2_summary <- summary(sim, par = "y_2")$summary

e_y_1 <- rstan::extract(sim, par = "y_1")$y_1
e_y_2 <- rstan::extract(sim, par = "y_2")$y_2

y_1 <- data.frame(reshape::melt(e_y_1[1:6, ]))
colnames(y_1) <- c("sample_id", "gene_name", "gene_usage_count")
y_1$condition <- "A"

y_2 <- data.frame(reshape::melt(e_y_2[1:6, ]))
colnames(y_2) <- c("sample_id", "gene_name", "gene_usage_count")
y_2$condition <- "B"

gene_usage <- rbind(y_1, y_2)
rm(y_1, y_2)

source("../IgGeneUsage/R/Util.R")
gene_usage <- get_paired_usage(u = gene_usage)
save(gene_usage, file = "dm_sim/gene_usage.RData")



model_dm <- rstan::stan_model(file = "../IgGeneUsage/inst/extdata/pair_dm_analytical.stan")
model_mn <- rstan::stan_model(file = "../IgGeneUsage/inst/extdata/pair_mn.stan")


fit_dm <- rstan::sampling(object = model_dm,
                          data = gene_usage,
                          chains = 5,
                          cores = 5,
                          iter = 3500,
                          warmup = 1500,
                          control = list(adapt_delta = 0.99, 
                                         max_treedepth = 13),
                          algorithm = "NUTS")
save(fit_dm, file = "dm_sim/fit_dm.RData")


fit_mn <- rstan::sampling(object = model_mn,
                          data = gene_usage,
                          chains = 5,
                          cores = 5,
                          iter = 3500,
                          warmup = 1500,
                          control = list(adapt_delta = 0.99, 
                                         max_treedepth = 13),
                          algorithm = "NUTS")
save(fit_mn, file = "dm_sim/fit_mn.RData")



require(rstan)
plot(fit_mn, par = "Y_hat_1")

# DM
summary_y_hat_1 <- data.frame(summary(fit_dm, par = "Y_hat_1")$summary)
summary_y_hat_1$Y <-as.vector(t(gene_usage$Y_1))
x <- do.call(rbind, strsplit(x = gsub(pattern = "Y_hat_1|\\[|\\]", 
                                      replacement = '', 
                                      x = rownames(summary_y_hat_1)),
                             split = ','))
summary_y_hat_1$gene_name <- x[, 1]
summary_y_hat_1$sample_id <- x[, 2]
summary_y_hat_1$condition <- "A"

summary_y_hat_2 <- data.frame(summary(fit_dm, par = "Y_hat_2")$summary)
summary_y_hat_2$Y <-as.vector(t(gene_usage$Y_2))
x <- do.call(rbind, strsplit(x = gsub(pattern = "Y_hat_2|\\[|\\]", 
                                      replacement = '', 
                                      x = rownames(summary_y_hat_2)),
                             split = ','))
summary_y_hat_2$gene_name <- x[, 1]
summary_y_hat_2$sample_id <- x[, 2]
summary_y_hat_2$condition <- "B"

summary_y_dm <- rbind(summary_y_hat_1, summary_y_hat_2)
summary_y_dm$model <- "DM"


# MN
summary_y_hat_1 <- data.frame(summary(fit_mn, par = "Y_hat_1")$summary)
summary_y_hat_1$Y <-as.vector(t(gene_usage$Y_1))
x <- do.call(rbind, strsplit(x = gsub(pattern = "Y_hat_1|\\[|\\]", 
                                      replacement = '', 
                                      x = rownames(summary_y_hat_1)),
                             split = ','))
summary_y_hat_1$gene_name <- x[, 1]
summary_y_hat_1$sample_id <- x[, 2]
summary_y_hat_1$condition <- "A"

summary_y_hat_2 <- data.frame(summary(fit_mn, par = "Y_hat_2")$summary)
summary_y_hat_2$Y <-as.vector(t(gene_usage$Y_2))
x <- do.call(rbind, strsplit(x = gsub(pattern = "Y_hat_2|\\[|\\]", 
                                      replacement = '', 
                                      x = rownames(summary_y_hat_2)),
                             split = ','))
summary_y_hat_2$gene_name <- x[, 1]
summary_y_hat_2$sample_id <- x[, 2]
summary_y_hat_2$condition <- "B"

summary_y_mn <- rbind(summary_y_hat_1, summary_y_hat_2)
summary_y_mn$model <- "MN"

summary_y <- rbind(summary_y_dm, summary_y_mn)
rm(summary_y_hat_1, summary_y_hat_2, summary_y_dm, summary_y_mn, x)

ggplot(data = summary_y)+
  geom_abline(slope = 1, intercept = 0)+
  facet_grid(facets = sample_id~model)+
  geom_point(aes(y = mean, x = Y, col = condition))+
  geom_errorbar(aes(y = mean, x = Y, ymin = X2.5., ymax = X97.5., col = condition))+
  theme_bw()


