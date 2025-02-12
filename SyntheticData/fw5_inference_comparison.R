### scatter plot of jacobian coefficients
library(dplyr)
library(magrittr)
library(ggplot2)
library(patchwork)

# load functions to calculate summary statistics and plot the results
source("functions/stats_jmat.R")

npop <- 5
nmodel <- 50
datasize <- 100

## -----------------------------------------------------------------------------
## results with process noise
jmat_smap <- readRDS(file="./results/smap_fw5_noise/jmat_smap.rds")
jmat_gpr <- readRDS(file="./results/fw5_noise/jmat_gpr.rds")

# create Reg. S-map's Jacobian matrix from the models with the optimal alpha values which gives the minimum in-sample error for each data set and target species
res_jmat_rsmap <- get_jmat_rsmap(modeloption = "fw5", noisecond = "noise")
jmat_rsmap <- res_jmat_rsmap$jmat_rsmap

list_jmat_scaled_fw5_temp <- list()
list_jmat_scaled_fw5 <- list()

for (i in 1:nmodel) {
  list_jmat_scaled_fw5_temp[[i]] <- paste0("./data/data_fw5_noise/jmat_true_fw5_noise_model", i, ".csv") %>% read.delim(header = FALSE) %>% as.matrix()
  list_jmat_scaled_fw5[[i]] <- array(data=matrix(rep(0,5^2), nrow = 5,ncol = 5), dim=c(5,5,datasize))
  
  for (j in 1:datasize) {
    list_jmat_scaled_fw5[[i]][,,j] <- list_jmat_scaled_fw5_temp[[i]][(1+5*(j-1)):(5*j),]
  }
}

inference_jmat_smap <- get_stats_ver2(jmat_smap, list_jmat_scaled_fw5)
inference_jmat_rsmap <- get_stats_ver2(jmat_rsmap, list_jmat_scaled_fw5)
inference_jmat_gpr <- get_stats_ver2(jmat_gpr, list_jmat_scaled_fw5)

# save boxplot
plot_rmse <- stats_boxplot(inference_jmat_smap = inference_jmat_smap, inference_jmat_rsmap_best = inference_jmat_rsmap, inference_jmat_gpr = inference_jmat_gpr,
              crit = "rmse", ylim = c(0,2.4), fontsize = 22)
ggsave("./fig/fw5/fig2a_rmse_fw5_noise.png", plot_rmse, dpi = 300, width = 8, height = 4)

plot_rho <- stats_boxplot(inference_jmat_smap = inference_jmat_smap, inference_jmat_rsmap_best = inference_jmat_rsmap, inference_jmat_gpr = inference_jmat_gpr,
              crit = "rho", ylim = c(0,1), fontsize = 22)
ggsave("./fig/fw5/fig2g_rho_fw5_noise.png", plot_rho, dpi = 300, width = 8, height = 4)


## -----------------------------------------------------------------------------
## results with high process noise
jmat_smap <- readRDS(file="./results/smap_fw5_highnoise/jmat_smap.rds")
jmat_gpr <- readRDS(file="./results/fw5_highnoise/jmat_gpr.rds")

# create Reg. S-map's Jacobian matrix from the models with the optimal alpha values which gives the minimum in-sample error for each data set and target species
res_jmat_rsmap <- get_jmat_rsmap(modeloption = "fw5", noisecond = "highnoise")
jmat_rsmap <- res_jmat_rsmap$jmat_rsmap

list_jmat_scaled_fw5_temp <- list()
list_jmat_scaled_fw5 <- list()

for (i in 1:nmodel) {
  list_jmat_scaled_fw5_temp[[i]] <- paste0("./data/data_fw5_highnoise/jmat_true_fw5_highnoise_model", i, ".csv") %>% read.delim(header = FALSE) %>% as.matrix()
  list_jmat_scaled_fw5[[i]] <- array(data=matrix(rep(0,npop^2), nrow = npop,ncol = npop), dim=c(npop,npop,datasize))
  
  for (j in 1:datasize) {
    list_jmat_scaled_fw5[[i]][,,j] <- list_jmat_scaled_fw5_temp[[i]][(1+npop*(j-1)):(npop*j),]
  }
}

# calculate performance criteria
inference_jmat_smap <- get_stats_ver2(jmat_smap, list_jmat_scaled_fw5)
inference_jmat_rsmap <- get_stats_ver2(jmat_rsmap, list_jmat_scaled_fw5)
inference_jmat_gpr <- get_stats_ver2(jmat_gpr, list_jmat_scaled_fw5)

# plot and save violin plots for RMSE and rho
plot_rmse <- stats_boxplot(inference_jmat_smap = inference_jmat_smap, inference_jmat_rsmap_best = inference_jmat_rsmap, inference_jmat_gpr = inference_jmat_gpr,
                           crit = "rmse", ylim = c(0,2.4),fontsize = 22)
ggsave("./fig/fw5/fig2d_rmse_fw5_highnoise.png", plot_rmse, dpi = 300, width = 8, height = 4)

plot_rho <- stats_boxplot(inference_jmat_smap = inference_jmat_smap, inference_jmat_rsmap_best = inference_jmat_rsmap, inference_jmat_gpr = inference_jmat_gpr,
                          crit = "rho", ylim = c(0,1), fontsize = 22)
ggsave("./fig/fw5/fig2j_rho_fw5_highnoise.png", plot_rho, dpi = 300, width = 8, height = 4)


## -----------------------------------------------------------------------------
## results with process and observational noise (high and modest level)
jmat_smap <- readRDS(file="./results/smap_fw5_highnoise_obs/jmat_smap.rds")
jmat_gpr <- readRDS(file="./results/fw5_highnoise_obs/jmat_gpr.rds")

# create Reg. S-map's Jacobian matrix from the models with the optimal alpha values which gives the minimum in-sample error for each data set and target species
res_jmat_rsmap <- get_jmat_rsmap(modeloption = "fw5", noisecond = "highnoise_obs")
jmat_rsmap <- res_jmat_rsmap$jmat_rsmap

list_jmat_scaled_fw5_temp <- list()
list_jmat_scaled_fw5 <- list()

for (i in 1:nmodel) {
  list_jmat_scaled_fw5_temp[[i]] <- paste0("./data/data_fw5_highnoise_obs/jmat_true_fw5_highnoise_obs_model", i, ".csv") %>% read.delim(header = FALSE) %>% as.matrix()
  list_jmat_scaled_fw5[[i]] <- array(data=matrix(rep(0,5^2), nrow = 5,ncol = 5), dim=c(5,5,datasize))
  
  for (j in 1:datasize) {
    list_jmat_scaled_fw5[[i]][,,j] <- list_jmat_scaled_fw5_temp[[i]][(1+5*(j-1)):(5*j),]
  }
}

# calculate performance criteria
inference_jmat_smap <- get_stats_ver2(jmat_smap, list_jmat_scaled_fw5)
inference_jmat_rsmap <- get_stats_ver2(jmat_rsmap, list_jmat_scaled_fw5)
inference_jmat_gpr <- get_stats_ver2(jmat_gpr, list_jmat_scaled_fw5)

# plot and save violin plots for RMSE and rho
plot_rmse <- stats_boxplot(inference_jmat_smap = inference_jmat_smap, inference_jmat_rsmap_best = inference_jmat_rsmap, inference_jmat_gpr = inference_jmat_gpr,
                           crit = "rmse", ylim = c(0,2.4),fontsize = 22)
ggsave("./fig/fw5/S4a_rmse_fw5_highnoise_obs.png", plot_rmse, dpi = 300, width = 8, height = 4)

plot_rho <- stats_boxplot(inference_jmat_smap = inference_jmat_smap, inference_jmat_rsmap_best = inference_jmat_rsmap, inference_jmat_gpr = inference_jmat_gpr,
                          crit = "rho", ylim = c(0,1), fontsize = 22)
ggsave("./fig/fw5/S4g_rho_fw5_highnoise_obs.png", plot_rho, dpi = 300, width = 8, height = 4)


## -----------------------------------------------------------------------------
## results with high level of process and observational noise
jmat_smap <- readRDS(file="./results/smap_fw5_highnoise_highobs/jmat_smap.rds")
jmat_gpr <- readRDS(file="./results/fw5_highnoise_highobs/jmat_gpr.rds")

# create Reg. S-map's Jacobian matrix from the models with the optimal alpha values which gives the minimum in-sample error for each data set and target species
res_jmat_rsmap <- get_jmat_rsmap(modeloption = "fw5", noisecond = "highnoise_highobs")
jmat_rsmap <- res_jmat_rsmap$jmat_rsmap

list_jmat_scaled_fw5_temp <- list()
list_jmat_scaled_fw5 <- list()

for (i in 1:nmodel) {
  list_jmat_scaled_fw5_temp[[i]] <- paste0("./data/data_fw5_highnoise_highobs/jmat_true_fw5_highnoise_highobs_model", i, ".csv") %>% read.delim(header = FALSE) %>% as.matrix()
  list_jmat_scaled_fw5[[i]] <- array(data=matrix(rep(0,5^2), nrow = 5,ncol = 5), dim=c(5,5,datasize))
  
  for (j in 1:datasize) {
    list_jmat_scaled_fw5[[i]][,,j] <- list_jmat_scaled_fw5_temp[[i]][(1+5*(j-1)):(5*j),]
  }
}

# calculate performance criteria
inference_jmat_smap <- get_stats_ver2(jmat_smap, list_jmat_scaled_fw5)
inference_jmat_rsmap <- get_stats_ver2(jmat_rsmap, list_jmat_scaled_fw5)
inference_jmat_gpr <- get_stats_ver2(jmat_gpr, list_jmat_scaled_fw5)

# plot and save violin plots for RMSE and rho
plot_rmse <- stats_boxplot(inference_jmat_smap = inference_jmat_smap, inference_jmat_rsmap_best = inference_jmat_rsmap, inference_jmat_gpr = inference_jmat_gpr,
                           crit = "rmse", ylim = c(0,2.4), fontsize = 22)
ggsave("./fig/fw5/S4d_rmse_fw5_highnoise_highobs.png", plot_rmse, dpi = 300, width = 8, height = 4)

plot_rho <- stats_boxplot(inference_jmat_smap = inference_jmat_smap, inference_jmat_rsmap_best = inference_jmat_rsmap, inference_jmat_gpr = inference_jmat_gpr,
                          crit = "rho", ylim = c(0,1), fontsize = 22)
ggsave("./fig/fw5/S4j_rho_fw5_highnoise_highobs.png", plot_rho, dpi = 300, width = 8, height = 4)
