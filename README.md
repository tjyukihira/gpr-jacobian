# Quantifying and exploring state-dependent ecological interactions from time series data using Gaussian process regression

This data set includes R and Julia scripts used in the study.

## Description of the scripts
### Synthetic data
In ./SyntheticData, data and scripts for the analyses of synsthetic data are located.

To create all synthetic time series for three artificial model systems, run "all_ts.jl" on Julia.
```julia
include("all_ts.jl")
```

To prepare all standardised time series data for model fitting procedures and theoretical Jacobian matrices for performance validation, run "all_generate_data_jmat.jl".
```julia
include("all_generate_data_jmat.jl")
```

On R, running "all_fit.R " fits all models (GPR, S-map and regularised S-map) to time series data.
```R
source("all_fit.R")
```

Running "fig2_figS2.R" creates figures used in the articles to compare inference performances (figure 2 and S2).
```R
source("fig2_figS2.R")
```

Running "fig3_figS3.jl" on Julia creates figures for performance validation of scenario exploration for species interactions (figure 3 and S3).
```julia
include("fig3_figS3.jl")
```

#### Each script used above
Running "ts_fw5.jl" and "ts_switching.jl" creates synthetic time series for food web model 1 and 2 respectively.
Running "ts_logistic.jl" creates synthetic time series for logistic map.

To prepare theoretical Jacobian matrix and data with observational noise for each theoretical model, run "generate_data_*modelname*.jl" files  ("*modelname*" = ["fw5", "switching", "logistic"]).
To prepare theoretical discretised Jacobian matrix for continuous time models, run "generate_jmat_discretised_*modelname*.jl"  ("*modelname*" = ["fw5", "switching"]).

To fit GPR and S-map models, run "fit_GPR_*modelname*.R" and "fit_smap_*modelname*.R" files.

To create figures used to compare the performances in the article (Figure 2 and S2), run "inference_comparison_*modelname*.R" files.

To create figures used to test the performance of GPR for scenario exploration (Figure 3 and S3), run "sim_jmat_*modelname*.jl" files.

"load_data_*modelname*.jl" are used for loading inferred Jacobian matrices for each model systems on Julia.

### Empirical data
In ./Beninca2009, data and scripts for the analyses of experimental mesocosm data (Beninca et. al. 2009, *Ecology Letters*) are located.

To fit GPR models, run "Beninca2009.jl".
```julia
include("Beninca2009.jl")
```

To create figure 4, 5 and S4, run "fig4_fig5_figS4.jl".
```julia
include("fig4_fig5_figS4.jl")
```

## Code/Software
Rprop algorithms for GPR were implmented using Julia version 1.9.3 (Bezanson et. al. 2017) and the packages StatsBase.jl (version 0.33.21), Distances.jl (0.10.8), DataFrames.jl (version 1.5.0) and JLD2.jl (version 0.4.31). We used DifferentialEquations.jl (version 7.10.0) and Distributions.jl (version 0.25.98) to create sysnthetic time series data. 

In the model fitting of GPR, we used R package JuliaCall (version 0.17.5) to run Julia scripts on R version 4.1.3 (R Core Team 2022). Analyses for S-map and Regularised S-map were performed using R and the R packages macam (version 0.1.4), rEDM (version 1.14.3), glmnet (version 4.1-7), foreach (version 1.5.2), doParallel (version 1.0.17), and parallelly (version 1.36.0). 

To create figures 3, 4, 5, S3, and S4, we used Julia packages Plots.jl (version 1.38.16), StatsPlots.jl (version 0.15.6) and LaTeXStrings.jl (version 1.3.0). In addition, R packages dplyr (version 1.1.2), ggplot2 (version 3.4.2), magrittr (version 2.0.3), and patchwork (version 1.1.3) were used for figures 2 and S2. 

For installation of macam, please see its [GitHub](https://github.com/ong8181/macam) page. You may be required to install packages from Bioconductor. If so, install them  using BiocManager package. For example,
```R
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("phyloseq")
```
