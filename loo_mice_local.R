require(rstan)


fs <- list.files(path = "dm_mice/", 
                 pattern = "_ds", 
                 full.names = TRUE)
fs_short <- list.files(path = "dm_mice/", 
                       pattern = "_ds", 
                       full.names = FALSE)



model_m <- rstan::stan_model(file = "dm_mice/m.stan")
model_dm <- rstan::stan_model(file = "dm_mice/dm.stan")


get_loo <- function(file_path, 
                    file_key,
                    out_path, 
                    model_m, 
                    model_dm) {
  
  ds <- get(load(file = file_path))
  for(u in 1:length(ds)) {
    d <- ds[[u]]
    
    ts <- vector(mode = "list", length = ncol(d$Y_1))
    fs_dm <- ts
    fs_m <- ts
    
    for(i in 1:ncol(d$Y_1)) {
      cat(i, "\n")
      ts[[i]] <- list(Y_1 = d$Y_1[,-i],
                      Y_2 = d$Y_2[,-i],
                      N = d$N[-i,],
                      N_sample = d$N_sample-1,
                      N_gene = d$N_gene,
                      Y_1_test = matrix(data = d$Y_1[,i], ncol = 1),
                      Y_2_test = matrix(data = d$Y_2[,i], ncol = 1),
                      N_test = matrix(data = d$N[i,], nrow = 1),
                      N_sample_test = 1)
      
      fs_dm[[i]] <- rstan::sampling(object = model_dm,
                                    data = ts[[i]],
                                    chains = 5,
                                    cores = 5,
                                    iter = 6500,
                                    warmup = 1500,
                                    control = list(adapt_delta = 0.99, 
                                                   max_treedepth = 14),
                                    algorithm = "NUTS")
      
      fs_m[[i]] <- rstan::sampling(object = model_m,
                                   data = ts[[i]],
                                   chains = 5,
                                   cores = 5,
                                   iter = 6500,
                                   warmup = 1500,
                                   control = list(adapt_delta = 0.99, 
                                                  max_treedepth = 14),
                                   algorithm = "NUTS")
      
    }
    loo <- list(fs_dm = fs_dm, fs_m = fs_m, ts = ts)
    rm(ts)
    out_name <- paste0(out_path, "/loo_", file_key, names(ds)[u])
    save(loo, out_name)
  }
  return("Done!")
}

get_loo(file_path = fs[1], 
        file_key = fs_short[1], 
        out_path = "dm_mice/", 
        model_m = model_m, 
        model_dm = model_dm)

get_loo(file_path = fs[1], 
        file_key = fs_short[1], 
        out_path = "dm_mice/", 
        model_m = model_m, 
        model_dm = model_dm)

get_loo(file_path = fs[1], 
        file_key = fs_short[1], 
        out_path = "dm_mice/", 
        model_m = model_m, 
        model_dm = model_dm)

get_loo(file_path = fs[1], 
        file_key = fs_short[1], 
        out_path = "dm_mice/", 
        model_m = model_m, 
        model_dm = model_dm)
