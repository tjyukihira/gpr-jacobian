### fitting GPR models to data under each noise condition

# load packages
## Use JuliaCall to call Julia functions for GPR
library(JuliaCall)
library(magrittr)

# load Julia functions for GPR inference
julia_source("functions/GPR.jl")
# load pre-defined functions 
source("functions/fit_gpr.R")

# -----------------------------------------------------------
# fitting GPR to data with process noise
res_gpr_noise <- fit_gprmodels(modeloption="logistic_10", Obs.noise=F, Processnoiselevel="low", save_juliadata=T)

# Jacobian
list_jmat_gpr_noise <- res_gpr_noise$list_jmat_gpr

# sd of Jacobian 
list_jmat_sd_noise <- res_gpr_noise$list_jmat_sd

saveRDS(list_jmat_gpr_noise, "results/logistic_10_noise/jmat_gpr.rds")
saveRDS(list_jmat_sd_noise, "results/logistic_10_noise/jmat_sd.rds")

# -----------------------------------------------------------
# fitting GPR to data with strong process noise
res_gpr_highnoise <- fit_gprmodels(modeloption="logistic_10", Obs.noise=F, Processnoiselevel="high", save_juliadata=T)

# Jacobian
list_jmat_gpr_highnoise <- res_gpr_highnoise$list_jmat_gpr

# sd of Jacobian 
list_jmat_sd_highnoise <- res_gpr_highnoise$list_jmat_sd

saveRDS(list_jmat_gpr_highnoise, "results/logistic_10_highnoise/jmat_gpr.rds")
saveRDS(list_jmat_sd_highnoise, "results/logistic_10_highnoise/jmat_sd.rds")

# -----------------------------------------------------------
# fitting GPR to data with process and observational noise
res_gpr_noise_obs <- fit_gprmodels(modeloption="logistic_10", Obs.noise=T, Processnoiselevel="high", Obs.noiselevel="low", save_juliadata=T)

# Jacobian
list_jmat_gpr_noise_obs <- res_gpr_noise_obs$list_jmat_gpr

# sd of Jacobian
list_jmat_sd_noise_obs <- res_gpr_noise_obs$list_jmat_sd

saveRDS(list_jmat_gpr_noise_obs, "results/logistic_10_highnoise_obs/jmat_gpr.rds")
saveRDS(list_jmat_sd_noise_obs, "results/logistic_10_highnoise_obs/jmat_sd.rds")

# -----------------------------------------------------------
# fitting GPR to data with process and observational noise
res_gpr_noise_highobs <- fit_gprmodels(modeloption="logistic_10", Obs.noise=T, Processnoiselevel="high", Obs.noiselevel="high", save_juliadata=T)

# Jacobian
list_jmat_gpr_noise_highobs <- res_gpr_noise_highobs$list_jmat_gpr

# sd of Jacobian
list_jmat_sd_noise_highobs <- res_gpr_noise_highobs$list_jmat_sd

saveRDS(list_jmat_gpr_noise_highobs, "results/logistic_10_highnoise_highobs/jmat_gpr.rds")
saveRDS(list_jmat_sd_noise_highobs, "results/logistic_10_highnoise_highobs/jmat_sd.rds")
