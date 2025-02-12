### Beninca et. al. (2009) Ecol. Lett.
using DelimitedFiles
using Plots
using StatsBase

include("functions/GPR.jl")

# load data
d_baltic = readdlm("data/Beninca_2009.csv", ',', header = true)

# use data after April, 1991 to ensure stationarity
ts_plankton = d_baltic[1][d_baltic[1][:, 1] .> 263, 2:end]

ts_lib = ts_plankton[1:end, :]

# the data is not scaled since the original data is already scaled
res_gpr = fit_gpr_synthetic(ts_lib; filename = "results/res_gpr.jld2", save_file = true)


