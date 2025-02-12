### functions for calculating  summary statistics of Jacobian matrix inference comparisons

# get Reg. S-map's Jacobian matrix according to the LOOCV results
get_jmat_rsmap <- function(modeloption, noisecond){
  jmat_ridge <- readRDS(file=paste0("results/smap_", modeloption, "_", noisecond, "/jmat_ridge.rds"))
  jmat_lasso <- readRDS(file=paste0("results/smap_", modeloption, "_", noisecond, "/jmat_lasso.rds"))
  
  for (iter in 1:9){
    jmat_elastic <- readRDS(file=paste0("results/smap_", modeloption, "_", noisecond, "/jmat_elastic_alpha0", iter, ".rds"))
    assign(paste0("jmat_elastic_alpha0", iter), jmat_elastic)
  }
  
  # loocv results for regularised S-map to decide the optimal alpha value which gives the minimum in-sample error (RMSE)
  loocv_ridge <- read.csv(file=paste0("results/smap_", modeloption, "_", noisecond, "/opt_param_ridge.csv"))
  loocv_lasso <- read.csv(file=paste0("results/smap_", modeloption, "_", noisecond, "/opt_param_lasso.csv"))
  for (iter in 1:9){
    loocv_elastic <- read.csv(file=paste0("results/smap_", modeloption, "_", noisecond, "/opt_param_elastic_alpha0", iter, ".csv"))
    assign(paste0("loocv_elastic_alpha0", iter), loocv_elastic)
  }
  
  # create Reg. S-map's Jacobian matrix from the models with the optimal alpha values for each data set and target
  df_alpha_rsmap <- data.frame(model=NA, target=NA, alpha=NA, rmse=NA)
  jmat_rsmap <- list()
  for (i in 1:nmodel) {
    jmat_rsmap[[i]] <- array(data=matrix(rep(0,npop^2), nrow = npop,ncol = npop), dim=c(npop,npop,datasize))
    for (j in 1:npop) {
      rmse_rsmap <- rep(0, 11)
      rmse_rsmap[1] <- loocv_ridge %>% filter(model == i & target == j) %>% select(rmse)
      rmse_rsmap[11] <- loocv_lasso %>% filter(model == i & target == j) %>% select(rmse)
      for (iter in 1:9){
        rmse_rsmap[iter+1] <- paste0("loocv_elastic_alpha0", iter) %>% get() %>% filter(model == i & target == j) %>% select(rmse)
      }
      alpha_opt <- (which.min(rmse_rsmap) - 1)/10
      rmse_min <- rmse_rsmap[which.min(rmse_rsmap)]
      df_alpha_rsmap[(i-1)*npop+j,] <- data.frame(model=i, target=j, alpha=alpha_opt, rmse=rmse_min)
      if(alpha_opt == 0.0){
        jmat_rsmap[[i]][j,,] <- jmat_ridge[[i]][j,,]
      }else if(alpha_opt == 1.0){
        jmat_rsmap[[i]][j,,] <- jmat_lasso[[i]][j,,]
      }else{
        jmat_rsmap[[i]][j,,] <- paste0("jmat_elastic_alpha0", which.min(rmse_rsmap) - 1) %>% get() %>% .[[i]] %>% .[j,,]
      }
    }
  }
  
  return(list(jmat_rsmap=jmat_rsmap, df_alpha_rsmap=df_alpha_rsmap))
}

# function to calculate summary statistics of the inference
get_stats <- function(list_jmat_smap, list_jmat_scaled, datasize=100){
  nmodel <- length(list_jmat_scaled)
  nsp <- ncol(list_jmat_scaled[[1]][,,1])
  
  vec_rmse <- vec_rho <- c()
  list_mat_rho <- list_mat_rmse <- list()
  
  for (model in 1:nmodel) {
    vec_rho_temp <- vec_rmse_temp <- c()
    list_mat_rho[[model]] <- list_mat_rmse[[model]] <- matrix(NA, nrow = nsp, ncol = nsp)
    for(i_row in 1:nsp){
      for (i_col in 1:nsp) {
        rho_temp <- suppressWarnings(cor(list_jmat_smap[[model]][i_row,i_col,1:(datasize-1)], list_jmat_scaled[[model]][i_row,i_col,1:(datasize-1)]))
        rmse_temp <- (list_jmat_smap[[model]][i_row,i_col,1:(datasize-1)] - list_jmat_scaled[[model]][i_row,i_col,1:(datasize-1)])^2 %>% mean %>% sqrt
        
        vec_rho_temp <- c(vec_rho_temp, rho_temp)
        vec_rmse_temp <- c(vec_rmse_temp, rmse_temp)
        list_mat_rho[[model]][i_row, i_col] <- rho_temp
        list_mat_rmse[[model]][i_row, i_col] <- rmse_temp
      }
    }
    
    vec_rmse[model] <- mean(vec_rmse_temp)
    vec_rho[model] <- mean(vec_rho_temp, na.rm = TRUE)
  }
  res_inference <- list(vec_rmse=vec_rmse, vec_rho=vec_rho, list_mat_rmse=list_mat_rmse, list_mat_rho=list_mat_rho)
  return(res_inference)
}

get_stats_ver2 <- function(list_jmat_est, list_jmat_true, datasize=100){
  nmodel <- length(list_jmat_true)
  nsp <- ncol(list_jmat_true[[1]][,,1])
  
  vec_rmse <- vec_rho <- c()
  
  for (model in 1:nmodel) {
    vec_rho[model] <- cor(c(list_jmat_est[[model]][,,1:(datasize-1)]), c(list_jmat_true[[model]][,,1:(datasize-1)]))
    vec_rmse[model] <- sqrt(mean((c(list_jmat_est[[model]][,,1:(datasize-1)]) - c(list_jmat_true[[model]][,,1:(datasize-1)]))^2))
  }
  res_inference <- list(vec_rmse=vec_rmse, vec_rho=vec_rho)
  return(res_inference)
}

# function for scatter plot of interaction strengths
scatter_jmat <- function(inference_jmat, list_jmat, model=1, ylim=c(-2,2), title=""){
  datasize <- dim(list_jmat[[1]])[3]
  npop <- ncol(list_jmat[[1]][,,1])
  vec_coeffs_true <- c()
  for(t in 1:(datasize-1)){
    vec_coeffs_true[(1+(t-1)*npop^2):(t*npop^2)] <- c(list_jmat[[model]][,,t])
  }
  df_sct <- data.frame(Estimated = inference_jmat$mat_coeffs[, model], Theoretical = vec_coeffs_true)
  ggplot(data = df_sct, aes(x = Theoretical, y = Estimated)) +
    geom_point(colour = "darkgreen", alpha = 0.3) +
    geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    geom_abline(slope = inference_jmat$vec_bhat[model], intercept = 0) +
    annotate("text", x=min(vec_coeffs_true)+0.28, y=ylim[2], label=paste("RMSE =", round(inference_jmat$vec_rmse[model], digits = 2))) +
    annotate("text", x=min(vec_coeffs_true)+0.2, y=ylim[2]-0.25, label=paste("beta ==", round(inference_jmat$vec_bhat[model], digits = 2)), parse = T) +
    annotate("text", x=min(vec_coeffs_true)+0.2, y=ylim[2]-0.5, label=paste("rho ==", round(inference_jmat$vec_rho[model], digits = 2)), parse = T) +
    ylim(ylim) +
    labs(title = title) +
    theme_light()
}

# function for boxplot of summary statistics of the inference
stats_boxplot <- function(inference_jmat_smap, inference_jmat_rsmap_best, inference_jmat_gpr, ylim=c(0,1), crit="rmse", fontsize=16){
  nmodel = length(inference_jmat_smap$vec_rho)
  df <- data.frame(rho = NA, rmse = NA, model = NA)
  #df[1:(nmodel*3), 1] <- c(inference_jmat_smap$vec_bhat, inference_jmat_rsmap_best$vec_bhat, inference_jmat_gpr$vec_bhat)
  df[1:(nmodel*3), 1] <- c(inference_jmat_smap$vec_rho, inference_jmat_rsmap_best$vec_rho, inference_jmat_gpr$vec_rho)
  df[1:(nmodel*3), 2] <- c(inference_jmat_smap$vec_rmse, inference_jmat_rsmap_best$vec_rmse, inference_jmat_gpr$vec_rmse)
  df[1:(nmodel*3), 3] <- c(rep("S-map", nmodel), rep("Regularised S-map" , nmodel), rep("GPR", nmodel))
  
  plt_rmse <- ggplot(df, aes(x = model, y = rmse, colour = model)) + 
    geom_violin() +
    geom_boxplot(width=.1) +
    geom_hline(yintercept = median(inference_jmat_gpr$vec_rmse) , lty = "dashed") +
    ylim(ylim) +
    theme_light() +
    theme(legend.position = "none", text = element_text(size = fontsize)) +
    ylab("RMSE") +
    xlab("")
  
  plt_rho <- ggplot(df, aes(x = model, y = rho, colour = model)) + 
    geom_violin() +
    geom_boxplot(width=.1) +
    geom_hline(yintercept = median(inference_jmat_gpr$vec_rho) , lty = "dashed") +
    theme_light() +
    ylim(ylim) +
    theme(legend.position = "none", text = element_text(size = fontsize)) +
    ylab(expression(rho)) +
    xlab("")
  
  if(crit=="rmse"){
    plt_rmse
  }else if(crit=="rho"){
    plt_rho
  }else if(crit=="all"){
    plt_rho - plt_rmse + plot_layout(nrow = 2)
  }
}