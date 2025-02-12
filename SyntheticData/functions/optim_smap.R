### Optimise Regularised S-map model

optim_smap <- function(data, regularised=TRUE, alpha=0, parallel=TRUE, nthreads=2, cl_type="PSOCK"){
  nmodel <- length(data)
  npop <- ncol(data[[1]])
  
  # candidate values of theta and lambda
  theta_test <- c(0, 0.1, 0.5, 1, 2, 3, 4, 6, 8)
  lambda_test <- c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 8)
  
  list_jmat_smap <- list()
  list_smap_stats <- list()
  
  len_grid = length(lambda_test)
  
  # function for grid-search of optimal theta and lambda
  fit_smap <- function(block, theta_test, lambda_test, alpha, regularised){
    npop <- ncol(block)
    list_jmat_smap_target <- list_smap_stats_target <- list()
    for (target in 1:npop) {#start of target
      y1_smap_stats <- data.frame(theta = NA, lambda = NA, N = NA, rho = NA, mae = NA, rmse = NA)
      list_smap_theta <- list()
      for(i in 1:length(theta_test))
      {#start of theta
        list_smap_lambda <- list()
        if(!regularised){
          # NOTE: glmnet_parallel = FALSE runs built-in parallel processing of extended_lnlp() from "macam", which may lead to slower fitting.
          y1_smap <- extended_lnlp(block, target_column = target, theta = theta_test[i], lambda = 0,
                                   regularized = FALSE, glmnet_parallel = TRUE, save_smap_coefficients = TRUE)
          for(j in 1:length(lambda_test))
          {#start of lambda
            list_smap_lambda[[j]] <- y1_smap$smap_coefficients[,-(npop+2)]
            
            # Summarise results
            y1_smap_stats[((i-1)*len_grid+j),] <- data.frame(theta = theta_test[i], lambda = 0, y1_smap$stats)
          }#end of lambda
        }else{
          for(j in 1:length(lambda_test))
          {#start of lambda
            y1_smap <- extended_lnlp(block, target_column = target, theta = theta_test[i], lambda = lambda_test[j],
                                     regularized = TRUE, alpha = alpha, glmnet_parallel = TRUE, save_smap_coefficients = TRUE)
            
            list_smap_lambda[[j]] <- y1_smap$smap_coefficients[,-(npop+2)]
            
            # Summarise results for each lambda and theta
            y1_smap_stats[((i-1)*len_grid+j),] <- data.frame(theta = theta_test[i], lambda = lambda_test[j], y1_smap$stats)
          }#end of lambda
        }#end of if(!regularised){}else{}
        
        list_smap_theta[[i]] <- list_smap_lambda
      }#end of theta
      list_smap_stats_target[[target]] <- y1_smap_stats
      
      list_jmat_smap_target[[target]] <- list_smap_theta
    }#end of target
    return(list(list_smap_stats_target=list_smap_stats_target, list_jmat_smap_target=list_jmat_smap_target))
  }
  
  nthreads <- ifelse(parallel, nthreads, 1)
  cl <- parallel::makeCluster(nthreads, type = cl_type, outfile="")
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  
  # Running parallel computing of fitting S-map models
  list_res_smap <- foreach(model=1:nmodel, .packages = "macam") %dopar% {
    cat("model", model, "\n")
    fit_smap(block = data[[model]], theta_test = theta_test, lambda_test = lambda_test, alpha = alpha, regularised = regularised)
  }
  
  opt_param_smap <- data.frame(model=NA, target=NA, theta=NA, lambda=NA, rmse=NA)
  
  jmat_smap <- list()
  
  # making estimated Jacobian matrices from S-map coefficients
  for (model in 1:nmodel) {
    list_jmat_smap <- c(list_jmat_smap, list(list_res_smap[[model]]$list_jmat_smap_target))
    list_smap_stats <- c(list_smap_stats, list(list_res_smap[[model]]$list_smap_stats_target))
    
    jmat_smap[[model]] <- array(data=matrix(nrow=npop,ncol=npop), dim = c(npop, npop, datasize))
  
    for (target in 1:npop) {
      list_smap_stats[[model]][[target]]$rmse
      opt_param_smap[((model-1)*npop+target),] <- data.frame(model=model, target=target, theta=list_smap_stats[[model]][[target]]$theta[which.min(list_smap_stats[[model]][[target]]$rmse)], lambda=list_smap_stats[[model]][[target]]$lambda[which.min(list_smap_stats[[model]][[target]]$rmse)], rmse=list_smap_stats[[model]][[target]]$rmse[which.min(list_smap_stats[[model]][[target]]$rmse)])
      theta_smap <- which(theta_test==opt_param_smap$theta[(opt_param_smap$model==model)&(opt_param_smap$target==target)], theta_test)
      lambda_smap <- which(lambda_test==opt_param_smap$lambda[(opt_param_smap$model==model)&(opt_param_smap$target==target)], lambda_test)
      
      jmat_smap[[model]][target,,1:datasize] <- list_jmat_smap[[model]][[target]][[theta_smap]][[lambda_smap]][,-1] %>% as.matrix() %>% t()
      }
  }
  
  return(list(jmat_opt=jmat_smap, opt_param=opt_param_smap, stats=list_smap_stats, jmat_all=list_jmat_smap))
}