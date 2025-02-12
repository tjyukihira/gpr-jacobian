### Scripts for fitting S-map models to 5 species food web model data

# load required packages
library(macam)
library(magrittr)
library(foreach)

# Load defined functions
source("functions/optim_smap.R")

# time series length and the number of datasets
datasize <- 100
nmodel <- 50

# sets of the candidate alpha values for elastic-net models
vec_alpha <- seq(0.1, 0.9, 0.1)

# number of threads for parallel computing
nthreads <- parallelly::availableCores(logical = FALSE) - 1
# On Linux, you can run the following to detect physical cores:
#nthreads <- nrow(unique(matrix(grep("^core id|^physical id", readLines("/proc/cpuinfo"), value = TRUE), ncol=2, byrow=TRUE))) - 1

# --------------------------------------------------------------------
# Data with process noise
cat("Beginning fitting to data with process noise \n")

data_fw5_noise <- list()
data_fw5_noise_tr <- list()

# load each synthetic time series data
for (i in 1:nmodel) {
  data_fw5_noise[[i]] <- paste0("data/data_fw5_noise/data_fw5_noise_model", i, ".csv") %>% 
    read.delim(header = FALSE, sep = "\t") %>% 
    as.matrix()
  
  # standardise data
  data_fw5_noise_tr[[i]] <- data_fw5_noise[[i]][1:datasize,] %>% scale()
}

# fitting each S-map model to data
# standard S-map
res_smap_noise <- optim_smap(data_fw5_noise_tr, regularised = FALSE, nthreads = nthreads)
# LASSO
res_lasso_noise <- optim_smap(data_fw5_noise_tr, regularised = TRUE, alpha = 1, nthreads=nthreads)
# Ridge
res_ridge_noise <- optim_smap(data_fw5_noise_tr, regularised = TRUE, alpha = 0, nthreads=nthreads)

# elastic net models with different values of alpha
for(iter in 1:length(vec_alpha)){
  res_elastic <- optim_smap(data_fw5_noise_tr, regularised = TRUE, alpha = vec_alpha[iter], nthreads = nthreads, cl_type = "PSOCK")
  assign(paste0("res_elastic_alpha0", iter, "_noise"), res_elastic)
  cat("alpha =", vec_alpha[iter], " finished \n")
}

# save the results
write.csv(res_smap_noise$opt_param, "results/smap_fw5_noise/opt_param_smap.csv")
write.csv(res_ridge_noise$opt_param, "results/smap_fw5_noise/opt_param_ridge.csv")
write.csv(res_lasso_noise$opt_param, "results/smap_fw5_noise/opt_param_lasso.csv")
for (iter in 1:9){
  paste0("res_elastic_alpha0", iter, "_noise") %>% get() %>% getElement("opt_param") %>% write.csv(file = paste0("results/smap_fw5_noise/opt_param_elastic_alpha0", iter, ".csv"))
  paste0("res_elastic_alpha0", iter, "_noise") %>% get() %>% getElement("jmat_opt") %>% saveRDS(file=paste0("results/smap_fw5_noise/jmat_elastic_alpha0", iter, ".rds"))
}

saveRDS(res_smap_noise$jmat_opt, file="results/smap_fw5_noise/jmat_smap.rds")
saveRDS(res_ridge_noise$jmat_opt, file="results/smap_fw5_noise/jmat_ridge.rds")
saveRDS(res_lasso_noise$jmat_opt, file="results/smap_fw5_noise/jmat_lasso.rds")

# --------------------------------------------------------------------
# Data with high level of process noise
cat("Beginning fitting to data with high level of process noise \n")

data_fw5_highnoise <- list()
data_fw5_highnoise_tr <- list()

# load each synthetic time series data
for (i in 1:nmodel) {
  data_fw5_highnoise[[i]] <- paste0("data/data_fw5_highnoise/data_fw5_highnoise_model", i, ".csv") %>% 
    read.delim(header = FALSE, sep = "\t") %>% 
    as.matrix()
  
  # standardise data
  data_fw5_highnoise_tr[[i]] <- data_fw5_highnoise[[i]][1:datasize,] %>% scale()
}

# fitting each S-map model to data
# standard S-map
res_smap_highnoise <- optim_smap(data_fw5_highnoise_tr, regularised = FALSE, nthreads = nthreads)
# LASSO
res_lasso_highnoise <- optim_smap(data_fw5_highnoise_tr, regularised = TRUE, alpha = 1, nthreads=nthreads)
# Ridge
res_ridge_highnoise <- optim_smap(data_fw5_highnoise_tr, regularised = TRUE, alpha = 0, nthreads=nthreads)

# elastic net models with different values of alpha
for(iter in 1:length(vec_alpha)){
  res_elastic <- optim_smap(data_fw5_highnoise_tr, regularised = TRUE, alpha = vec_alpha[iter], nthreads = nthreads, cl_type = "PSOCK")
  assign(paste0("res_elastic_alpha0", iter, "_highnoise"), res_elastic)
  cat("alpha =", vec_alpha[iter], " finished \n")
}

# save the results
write.csv(res_smap_highnoise$opt_param, "results/smap_fw5_highnoise/opt_param_smap.csv")
write.csv(res_ridge_highnoise$opt_param, "results/smap_fw5_highnoise/opt_param_ridge.csv")
write.csv(res_lasso_highnoise$opt_param, "results/smap_fw5_highnoise/opt_param_lasso.csv")
for (iter in 1:9){
  paste0("res_elastic_alpha0", iter, "_highnoise") %>% get() %>% getElement("opt_param") %>% write.csv(file = paste0("results/smap_fw5_highnoise/opt_param_elastic_alpha0", iter, ".csv"))
  paste0("res_elastic_alpha0", iter, "_highnoise") %>% get() %>% getElement("jmat_opt") %>% saveRDS(file=paste0("results/smap_fw5_highnoise/jmat_elastic_alpha0", iter, ".rds"))
}

saveRDS(res_smap_highnoise$jmat_opt, file="results/smap_fw5_highnoise/jmat_smap.rds")
saveRDS(res_ridge_highnoise$jmat_opt, file="results/smap_fw5_highnoise/jmat_ridge.rds")
saveRDS(res_lasso_highnoise$jmat_opt, file="results/smap_fw5_highnoise/jmat_lasso.rds")

# --------------------------------------------------------------------
# Data with process and observational noise
cat("Beginning fitting to data with process and observational noise (high and modest level) \n")

data_fw5_highnoise_obs <- list()
data_fw5_highnoise_obs_tr <- list()

# load each synthetic time series data
for (i in 1:nmodel) {
  data_fw5_highnoise_obs[[i]] <- paste0("data/data_fw5_highnoise_obs/data_fw5_highnoise_obs_model", i, ".csv") %>% 
    read.delim(header = FALSE, sep = "\t") %>% 
    as.matrix()
  
  # standardise data
  data_fw5_highnoise_obs_tr[[i]] <- data_fw5_highnoise_obs[[i]][1:datasize,] %>% scale()
}

# fitting each S-map model to data
# standard S-map
res_smap_highnoise_obs <- optim_smap(data_fw5_highnoise_obs_tr, regularised = FALSE, nthreads = nthreads)
# LASSO
res_lasso_highnoise_obs <- optim_smap(data_fw5_highnoise_obs_tr, regularised = TRUE, alpha = 1, nthreads=nthreads)
# Ridge
res_ridge_highnoise_obs <- optim_smap(data_fw5_highnoise_obs_tr, regularised = TRUE, alpha = 0, nthreads=nthreads)

# elastic net models with different values of alpha
for(iter in 1:length(vec_alpha)){
  res_elastic <- optim_smap(data_fw5_highnoise_obs_tr, regularised = TRUE, alpha = vec_alpha[iter], nthreads = nthreads, cl_type = "PSOCK")
  assign(paste0("res_elastic_alpha0", iter, "_highnoise_obs"), res_elastic)
  cat("alpha =", vec_alpha[iter], " finished \n")
}

# save the results
write.csv(res_smap_highnoise_obs$opt_param, "results/smap_fw5_highnoise_obs/opt_param_smap.csv")
write.csv(res_ridge_highnoise_obs$opt_param, "results/smap_fw5_highnoise_obs/opt_param_ridge.csv")
write.csv(res_lasso_highnoise_obs$opt_param, "results/smap_fw5_highnoise_obs/opt_param_lasso.csv")
for (iter in 1:9){
  paste0("res_elastic_alpha0", iter, "_highnoise_obs") %>% get() %>% getElement("opt_param") %>% write.csv(file = paste0("results/smap_fw5_highnoise_obs/opt_param_elastic_alpha0", iter, ".csv"))
  paste0("res_elastic_alpha0", iter, "_highnoise_obs") %>% get() %>% getElement("jmat_opt") %>% saveRDS(file=paste0("results/smap_fw5_highnoise_obs/jmat_elastic_alpha0", iter, ".rds"))
}

saveRDS(res_smap_highnoise_obs$jmat_opt, file="results/smap_fw5_highnoise_obs/jmat_smap.rds")
saveRDS(res_ridge_highnoise_obs$jmat_opt, file="results/smap_fw5_highnoise_obs/jmat_ridge.rds")
saveRDS(res_lasso_highnoise_obs$jmat_opt, file="results/smap_fw5_highnoise_obs/jmat_lasso.rds")

# --------------------------------------------------------------------
# Data with high level of process and observational noise
cat("Beginning fitting to data with high level of process and observational noise \n")

data_fw5_highnoise_highobs <- list()
data_fw5_highnoise_highobs_tr <- list()

# load each synthetic time series data
for (i in 1:nmodel) {
  data_fw5_highnoise_highobs[[i]] <- paste0("data/data_fw5_highnoise_highobs/data_fw5_highnoise_highobs_model", i, ".csv") %>% 
    read.delim(header = FALSE, sep = "\t") %>% 
    as.matrix()
  
  # standardise data
  data_fw5_highnoise_highobs_tr[[i]] <- data_fw5_highnoise_highobs[[i]][1:datasize,] %>% scale()
}

# fitting each S-map model to data
# standard S-map
res_smap_highnoise_highobs <- optim_smap(data_fw5_highnoise_highobs_tr, regularised = FALSE, nthreads = nthreads)
# LASSO
res_lasso_highnoise_highobs <- optim_smap(data_fw5_highnoise_highobs_tr, regularised = TRUE, alpha = 1, nthreads=nthreads)
# Ridge
res_ridge_highnoise_highobs <- optim_smap(data_fw5_highnoise_highobs_tr, regularised = TRUE, alpha = 0, nthreads=nthreads)

# fitting elastic net models with different values of alpha
for(iter in 1:length(vec_alpha)){
  res_elastic <- optim_smap(data_fw5_highnoise_highobs_tr, regularised = TRUE, alpha = vec_alpha[iter], nthreads = nthreads, cl_type = "PSOCK")
  assign(paste0("res_elastic_alpha0", iter, "_highnoise_highobs"), res_elastic)
  cat("alpha =", vec_alpha[iter], " finished \n")
}

# save the results
write.csv(res_smap_highnoise_highobs$opt_param, "results/smap_fw5_highnoise_highobs/opt_param_smap.csv")
write.csv(res_ridge_highnoise_highobs$opt_param, "results/smap_fw5_highnoise_highobs/opt_param_ridge.csv")
write.csv(res_lasso_highnoise_highobs$opt_param, "results/smap_fw5_highnoise_highobs/opt_param_lasso.csv")
for (iter in 1:9){
  paste0("res_elastic_alpha0", iter, "_highnoise_highobs") %>% get() %>% getElement("opt_param") %>% write.csv(file = paste0("results/smap_fw5_highnoise_highobs/opt_param_elastic_alpha0", iter, ".csv"))
  paste0("res_elastic_alpha0", iter, "_highnoise_highobs") %>% get() %>% getElement("jmat_opt") %>% saveRDS(file=paste0("results/smap_fw5_highnoise_highobs/jmat_elastic_alpha0", iter, ".rds"))
}

saveRDS(res_smap_highnoise_highobs$jmat_opt, file="results/smap_fw5_highnoise_highobs/jmat_smap.rds")
saveRDS(res_ridge_highnoise_highobs$jmat_opt, file="results/smap_fw5_highnoise_highobs/jmat_ridge.rds")
saveRDS(res_lasso_highnoise_highobs$jmat_opt, file="results/smap_fw5_highnoise_highobs/jmat_lasso.rds")
