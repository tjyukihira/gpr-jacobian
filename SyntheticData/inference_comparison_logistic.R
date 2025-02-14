### scatter plot of jacobian coefficients
library(dplyr)
library(magrittr)
library(ggplot2)
library(patchwork)

# load functions to calculate summary statistics and plot the results
source("functions/stats_jmat.R")

npop <- 10
nmodel <- 50
datasize <- 100

## -----------------------------------------------------------------------------
## results with process noise
jmat_smap <- readRDS(file="./results/smap_logistic_10_noise/jmat_smap.rds")
jmat_gpr <- readRDS(file="./results/logistic_10_noise/jmat_gpr.rds")

# create Reg. S-map's Jacobian matrix from the models with the optimal alpha values which gives the minimum in-sample error for each data set and target species
res_jmat_rsmap <- get_jmat_rsmap(modeloption = "logistic_10", noisecond = "noise")
jmat_rsmap <- res_jmat_rsmap$jmat_rsmap

list_jmat_scaled_logistic_10_temp <- list()
list_jmat_scaled_logistic_10 <- list()

for (i in 1:nmodel) {
  list_jmat_scaled_logistic_10_temp[[i]] <- paste0("./data/data_logistic_10_noise/jmat_logistic_10_noise_model", i, ".csv") %>% read.delim(header = FALSE) %>% as.matrix()
  list_jmat_scaled_logistic_10[[i]] <- array(data=matrix(rep(0,npop^2), nrow = npop,ncol = npop), dim=c(npop,npop,datasize))
  
  for (j in 1:datasize) {
    list_jmat_scaled_logistic_10[[i]][,,j] <- list_jmat_scaled_logistic_10_temp[[i]][(1+npop*(j-1)):(npop*j),]
  }
}

inference_jmat_smap <- get_stats(jmat_smap, list_jmat_scaled_logistic_10)
inference_jmat_rsmap <- get_stats(jmat_rsmap, list_jmat_scaled_logistic_10)
inference_jmat_gpr <- get_stats(jmat_gpr, list_jmat_scaled_logistic_10)

# plot and save violin plots for RMSE and rho
plot_rmse <- stats_boxplot(inference_jmat_smap = inference_jmat_smap, inference_jmat_rsmap_best = inference_jmat_rsmap, inference_jmat_gpr = inference_jmat_gpr,
              crit = "rmse", fontsize = 22, ylim = c(0,1))
ggsave("./fig/logistic_10/fig2c_rmse_logistic_10_noise.png", plot_rmse, dpi = 300, width = 8, height = 4)

plot_rho <- stats_boxplot(inference_jmat_smap = inference_jmat_smap, inference_jmat_rsmap_best = inference_jmat_rsmap, inference_jmat_gpr = inference_jmat_gpr,
              crit = "rho", fontsize = 22, ylim=c(0, 1.0))
ggsave("./fig/logistic_10/fig2i_rho_logistic_10_noise.png", plot_rho, dpi = 300, width = 8, height = 4)


## -----------------------------------------------------------------------------
## results with high process noise
jmat_smap <- readRDS(file="./results/smap_logistic_10_highnoise/jmat_smap.rds")
jmat_gpr <- readRDS(file="./results/logistic_10_highnoise/jmat_gpr.rds")

# create Reg. S-map's Jacobian matrix from the models with the optimal alpha values which gives the minimum in-sample error for each data set and target species
res_jmat_rsmap <- get_jmat_rsmap(modeloption = "logistic_10", noisecond = "highnoise")
jmat_rsmap <- res_jmat_rsmap$jmat_rsmap

list_jmat_scaled_logistic_10_temp <- list()
list_jmat_scaled_logistic_10 <- list()

for (i in 1:nmodel) {
  list_jmat_scaled_logistic_10_temp[[i]] <- paste0("./data/data_logistic_10_highnoise/jmat_logistic_10_highnoise_model", i, ".csv") %>% read.delim(header = FALSE) %>% as.matrix()
  list_jmat_scaled_logistic_10[[i]] <- array(data=matrix(rep(0,npop^2), nrow = npop,ncol = npop), dim=c(npop,npop,datasize))
  
  for (j in 1:datasize) {
    list_jmat_scaled_logistic_10[[i]][,,j] <- list_jmat_scaled_logistic_10_temp[[i]][(1+npop*(j-1)):(npop*j),]
  }
}

# calculate performance criteria
inference_jmat_smap <- get_stats(jmat_smap, list_jmat_scaled_logistic_10)
inference_jmat_rsmap <- get_stats(jmat_rsmap, list_jmat_scaled_logistic_10)
inference_jmat_gpr <- get_stats(jmat_gpr, list_jmat_scaled_logistic_10)

# plot and save violin plots for RMSE and rho
plot_rmse <- stats_boxplot(inference_jmat_smap = inference_jmat_smap, inference_jmat_rsmap_best = inference_jmat_rsmap, inference_jmat_gpr = inference_jmat_gpr,
                           crit = "rmse", fontsize = 22, ylim=c(0,1))
ggsave("./fig/logistic_10/fig2f_rmse_logistic_10_highnoise.png", plot_rmse, dpi = 300, width = 8, height = 4)

plot_rho <- stats_boxplot(inference_jmat_smap = inference_jmat_smap, inference_jmat_rsmap_best = inference_jmat_rsmap, inference_jmat_gpr = inference_jmat_gpr,
                          crit = "rho", fontsize = 22, ylim=c(0, 1.0))
ggsave("./fig/logistic_10/fig2l_rho_logistic_10_highnoise.png", plot_rho, dpi = 300, width = 8, height = 4)


## -----------------------------------------------------------------------------
## results with process and observational noise (high and modest level)
jmat_smap <- readRDS(file="./results/smap_logistic_10_highnoise_obs/jmat_smap.rds")
jmat_gpr <- readRDS(file="./results/logistic_10_highnoise_obs/jmat_gpr.rds")

# create Reg. S-map's Jacobian matrix from the models with the optimal alpha values which gives the minimum in-sample error for each data set and target species
res_jmat_rsmap <- get_jmat_rsmap(modeloption = "logistic_10", noisecond = "highnoise_obs")
jmat_rsmap <- res_jmat_rsmap$jmat_rsmap

list_jmat_scaled_logistic_10_temp <- list()
list_jmat_scaled_logistic_10 <- list()

for (i in 1:nmodel) {
  list_jmat_scaled_logistic_10_temp[[i]] <- paste0("./data/data_logistic_10_highnoise_obs/jmat_logistic_10_highnoise_obs_scaled_model", i, ".csv") %>% read.delim(header = FALSE) %>% as.matrix()
  list_jmat_scaled_logistic_10[[i]] <- array(data=matrix(rep(0,npop^2), nrow = npop,ncol = npop), dim=c(npop,npop,datasize))
  
  for (j in 1:datasize) {
    list_jmat_scaled_logistic_10[[i]][,,j] <- list_jmat_scaled_logistic_10_temp[[i]][(1+npop*(j-1)):(npop*j),]
  }
}

# calculate performance criteria
inference_jmat_smap <- get_stats(jmat_smap, list_jmat_scaled_logistic_10)
inference_jmat_rsmap <- get_stats(jmat_rsmap, list_jmat_scaled_logistic_10)
inference_jmat_gpr <- get_stats(jmat_gpr, list_jmat_scaled_logistic_10)

# plot and save violin plots for RMSE and rho
plot_rmse <- stats_boxplot(inference_jmat_smap = inference_jmat_smap, inference_jmat_rsmap_best = inference_jmat_rsmap, inference_jmat_gpr = inference_jmat_gpr,
                           crit = "rmse", fontsize = 22, ylim=c(0,1))
ggsave("./fig/logistic_10/S4c_rmse_logistic_10_highnoise_obs.png", plot_rmse, dpi = 300, width = 8, height = 4)

plot_rho <- stats_boxplot(inference_jmat_smap = inference_jmat_smap, inference_jmat_rsmap_best = inference_jmat_rsmap, inference_jmat_gpr = inference_jmat_gpr,
                          crit = "rho", fontsize = 22, ylim=c(0,1))
ggsave("./fig/logistic_10/S4i_rho_logistic_10_highnoise_obs.png", plot_rho, dpi = 300, width = 8, height = 4)


## -----------------------------------------------------------------------------
## results with high level of process and observational noise
jmat_smap <- readRDS(file="./results/smap_logistic_10_highnoise_highobs/jmat_smap.rds")
jmat_gpr <- readRDS(file="./results/logistic_10_highnoise_highobs/jmat_gpr.rds")

# create Reg. S-map's Jacobian matrix from the models with the optimal alpha values which gives the minimum in-sample error for each data set and target species
res_jmat_rsmap <- get_jmat_rsmap(modeloption = "logistic_10", noisecond = "highnoise_highobs")
jmat_rsmap <- res_jmat_rsmap$jmat_rsmap

list_jmat_scaled_logistic_10_temp <- list()
list_jmat_scaled_logistic_10 <- list()

for (i in 1:nmodel) {
  list_jmat_scaled_logistic_10_temp[[i]] <- paste0("./data/data_logistic_10_highnoise_highobs/jmat_logistic_10_highnoise_highobs_scaled_model", i, ".csv") %>% read.delim(header = FALSE) %>% as.matrix()
  list_jmat_scaled_logistic_10[[i]] <- array(data=matrix(rep(0,npop^2), nrow = npop,ncol = npop), dim=c(npop,npop,datasize))
  
  for (j in 1:datasize) {
    list_jmat_scaled_logistic_10[[i]][,,j] <- list_jmat_scaled_logistic_10_temp[[i]][(1+npop*(j-1)):(npop*j),]
  }
}

# calculate performance criteria
inference_jmat_smap <- get_stats(jmat_smap, list_jmat_scaled_logistic_10)
inference_jmat_rsmap <- get_stats(jmat_rsmap, list_jmat_scaled_logistic_10)
inference_jmat_gpr <- get_stats(jmat_gpr, list_jmat_scaled_logistic_10)

# plot and save violin plots for RMSE and rho
plot_rmse <- stats_boxplot(inference_jmat_smap = inference_jmat_smap, inference_jmat_rsmap_best = inference_jmat_rsmap, inference_jmat_gpr = inference_jmat_gpr,
                           crit = "rmse", fontsize = 22, ylim=c(0,1))
ggsave("./fig/logistic_10/S4f_rmse_logistic_10_highnoise_highobs.png", plot_rmse, dpi = 300, width = 8, height = 4)

plot_rho <- stats_boxplot(inference_jmat_smap = inference_jmat_smap, inference_jmat_rsmap_best = inference_jmat_rsmap, inference_jmat_gpr = inference_jmat_gpr,
                          crit = "rho", fontsize = 22, ylim=c(0,1))
ggsave("./fig/logistic_10/S4l_rho_logistic_10_highnoise_highobs.png", plot_rho, dpi = 300, width = 8, height = 4)
