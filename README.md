# Quantifying and exploring state-dependent ecological interactions from time series data using Gaussian process regression

This data set includes R and Julia scripts used in the study.

## Description of the scripts
### Synthetic data
In ./SyntheticData, data and scripts for the analyses of synsthetic data are located.

Running "sim_fw5.jl" and "sim_prey_switching.jl" creates synthetic time series for food web model 1 and 2 respectively.
Running "sim_logistic_10.jl" creates synthetic time series for logistic map.

To prepare theoretical Jacobian matrix and data with observational noise for each theoretical model, run "generate_data_*modelname*.jl" files  ("*modelname*"= "fw5", "switching", "logistic").
To prepare theoretical discretised Jacobian matrix for continuous time models, run "generate_jmat_discretised_*modelname*.jl"  ("*modelname*"= "fw5", "switching").

To fit GPR and S-map models, run "fit_GPR_*modelname*.R" and "fit_smap_*modelname*.R" files.

To create figures used to compare the performances in the article, run "*modelname*_inference_comparison.R" files.

To create figures used to test the performance of GPR for scenario exploration, run "sim_jmat_*modelname*.jl" files.

### Empirical data
In ./Beninca2009, data and scripts for the analyses of experimental mesocosm data are located.

To fit GPR models, run "Beninca2009.jl".
To create figure in the article, run "Beninca2009_figs.jl".

## Code/Software
Analyses for GPR were performed using Julia version 1.9.3 (Bezanson et. al. 2017) and DifferentialEquations.jl (version 7.10.0 for Solving Stochastic Differential Equations). In the model fitting of GPR, we used R package JuliaCall (version 0.17.5) to run Julia scripts on R version 4.1.3 (R Core Team 2022). Analyses for S-map and Regularised S-map were performed using R and the R packages macam (version 0.1.4, for the Regularised S-map and S-map), rEDM (version 1.14.3, for the S-map analysis), glmnet (version 4.1-7, for the regularised S-map analysis).
