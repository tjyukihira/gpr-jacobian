# Quantifying and exploring state-dependent ecological interactions from time series data using Gaussian process regression

This data set includes R and Julia scripts used in the study.

## Description of the scripts
### Synthetic data
In ./SyntheticData, data and scripts for the analyses of synsthetic data are located.

To create all synthetic time series for three artificial model systems, run "all_ts.jl" on Julia.
```julia:all_ts.jl
include("all_ts.jl")
```

To prepare all standaridised time series data for model fitting procedures and theoretical Jacobian matrices for performance validation, run "all_generate_data_jmat.jl".
```julia:all_ts.jl
include("all_generate_data_jmat.jl")
```

On R, running "all_fit.R " fits all models (GPR, S-map and regularised S-map) to time series data.
```R
source("all_fit.R")
```

Running "all_inference_comparison.R" creates figures used in the articles to compare inference performances.
```R
source("all_inference_comparison.R")
```

Running "all_sim_jmat.jl" on Julia creates figures for performance validation of scenario exploration for species interactions.
```julia
include("all_sim_jmat.jl")
```

#### All scripts used above
Running "ts_fw5.jl" and "ts_switching.jl" creates synthetic time series for food web model 1 and 2 respectively.
Running "ts_logistic.jl" creates synthetic time series for logistic map.

To prepare theoretical Jacobian matrix and data with observational noise for each theoretical model, run "generate_data_*modelname*.jl" files  ("*modelname*" = ["fw5", "switching", "logistic"]).
To prepare theoretical discretised Jacobian matrix for continuous time models, run "generate_jmat_discretised_*modelname*.jl"  ("*modelname*" = ["fw5", "switching"]).

To fit GPR and S-map models, run "fit_GPR_*modelname*.R" and "fit_smap_*modelname*.R" files.

To create figures used to compare the performances in the article (Figure 2 and S2), run "inference_comparison_*modelname*.R" files.

To create figures used to test the performance of GPR for scenario exploration (Figure 3 and S3), run "sim_jmat_*modelname*.jl" files.

"load_data_*modelname*.jl" are used for loading inferred Jacobian matrices for each model systems on Julia.

### Empirical data
In ./Beninca2009, data and scripts for the analyses of experimental mesocosm data are located.

To fit GPR models, run "Beninca2009.jl".
```julia
include("Beninca2009.jl")
```

To create figures in the article (Figure 4, 5 and S4), run "Beninca2009_figs.jl".
```julia
include("Beninca2009_figs.jl")
```

## Code/Software
Analyses for GPR were performed using Julia version 1.9.3 (Bezanson et. al. 2017) and DifferentialEquations.jl (version 7.10.0 for Solving Stochastic Differential Equations). In the model fitting of GPR, we used R package JuliaCall (version 0.17.5) to run Julia scripts on R version 4.1.3 (R Core Team 2022). Analyses for S-map and Regularised S-map were performed using R and the R packages macam (version 0.1.4, for the Regularised S-map and S-map), rEDM (version 1.14.3, for the S-map analysis), glmnet (version 4.1-7, for the regularised S-map analysis).

For installation of macam, please see its [https://github.com/ong8181/macam](GitHub) page. You may be required to install packages from Bioconductor. If so, install them  using BiocManager package. For example,
```R
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("phyloseq")
```
