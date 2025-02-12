### function for GPR inference by calling julia functions

fit_gpr <- function(lib, save_juliadata=FALSE, file=NULL){
  if(save_juliadata){
    res_gpr <- julia_call("fit_gpr_synthetic", lib, filename=file, save_file=save_juliadata)
  }else{
    res_gpr <- julia_call("fit_gpr_synthetic", lib)
  }
  return(res_gpr)
}

fit_gprmodels <- function(modeloption="fw5", Obs.noise=FALSE, Processnoiselevel="low", Obs.noiselevel="low", save_juliadata=TRUE){
  
  data <- list()
  data_tr <- list()
  
  datasize <- 100
  nmodel <- 50
  
  if(modeloption=="fw5"){
    dataname <- "data_fw5_"
    
    npop <- 5
  }else if(modeloption=="switching"){
    dataname <- "data_switching_"
    
    npop <- 5
  }else if(modeloption=="logistic_10"){
    dataname <- "data_logistic_10_"
    
    npop <- 10
  }
  
  if(!Obs.noise){
    if(Processnoiselevel=="low"){
      noisecond <- "noise"
    }else if(Processnoiselevel=="high"){
      noisecond <- "highnoise" 
    }
  }else if(Obs.noise){
    if(Obs.noiselevel=="low"){
      noisecond <- "highnoise_obs"
    }else if(Obs.noiselevel=="high"){
      noisecond <- "highnoise_highobs"
    }
  }
  # data directory
  datadir <- paste0("data/", dataname, noisecond)
  
  # standardise data
  for (i in 1:nmodel) {
    data[[i]] <- paste0(datadir, "/" , dataname, noisecond, "_model", i, ".csv") %>% 
      read.delim(header = FALSE, sep = "\t") %>% 
      as.matrix()
    data_tr[[i]] <- data[[i]][1:datasize,] %>% scale()
  }
  
  res_gpr <- list()
  
  for(model in 1:nmodel) {
    filename <- paste0("results/", modeloption, "_", noisecond, "/model", model, ".jld2")
    res_gpr[[model]] <- fit_gpr(data_tr[[model]], save_juliadata = save_juliadata, file = filename)
    cat("model", model, "finished \n")
  }
  
  list_jmat_gpr <- list()
  list_jmat_sd <- list()
  
  for(i in 1:nmodel){
    list_jmat_gpr[[i]] <- res_gpr[[i]]$jmat[,,1:(datasize-1)]
    list_jmat_sd[[i]] <- res_gpr[[i]]$jmat_sd[,,1:(datasize-1)]
  }
  
  res_gprmodels <- list(list_jmat_gpr=list_jmat_gpr, list_jmat_sd=list_jmat_sd)
  return(res_gprmodels)
}