###Synthetic data of logistic map
###

using LinearAlgebra
using StatsBase
import Random.seed!
import Plots.plot
import DelimitedFiles.writedlm
import Dates.today
import Dates.format
import DataFrames.Not
include("functions/fun_logistic_map.jl")

# # of populations
npop = 10

# # of models
nmodel = 50

# intrinsic growth rates are sampled from Uniform(r_growth - 1.0, r_growth + 1.0)
r_growth = 3.2

# 0.05 <= abs(interaction strengths) <= 0.5
inter = 0.5

# probality with which each species interacts with other species
connectance = 1/3

# process noise intensity is 0.05 and 0.1.
sd_noise = 0.05

datasize = 100

# # of steps
iter = 1000

# set seed to replicate the results
seed = 3333

generate_logistic(nmodel = nmodel, r_growth = r_growth, inter = inter, connectance = connectance, sd_noise = sd_noise, npop= npop, datasize = datasize, iter = iter, seed = seed)
