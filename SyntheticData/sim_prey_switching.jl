### generate synthetic time series data from 5 species food web model with prey switching

using DifferentialEquations
using Plots
using DelimitedFiles
import StatsBase.mean
import Distributions.Uniform
import Random.seed!

include("functions/fun_prey_switching.jl")

seed!(8888)

# Parameter range is determined so that the resulting dynamics become nonstatic according to Post et. al. (2000).
pref = rand(Uniform(0.36, 0.64), 50)

tspan = (0.0, 599.999)

nmodel = 50

## data with process noise
## discard the results with low abundance
data_switching_noise = sim_prey_swtching_noise(nmodel = nmodel, pref = pref, tspan = tspan, intense = false)

plt_sde = plot(data_switching_noise[49], legend=false, title="Process noise", titlefontsize=9)


## data with stronger process noise
## discard the results with low abundance
data_switching_intensenoise = sim_prey_swtching_noise(nmodel = nmodel, pref = pref, tspan = tspan, intense = true)

plt_sde_intense = plot(data_switching_intensenoise[49], legend=false, title="Strong process noise", titlefontsize=9)

for i in 1:nmodel
    writedlm("data/data_switching_noise/data_switching_noise_model$(i).csv", data_switching_noise[i])
    writedlm("data/data_switching_highnoise/data_switching_highnoise_model$(i).csv", data_switching_intensenoise[i])
end

writedlm("data/data_switching/prey_preference.csv", pref)

plot(plt_sde, plt_sde_intense, layout=(1,:), title = "π = $(pref[49])")