using DifferentialEquations
using Plots
using DelimitedFiles
import Distributions.Uniform
import StatsBase.mean
import Random.seed!

include("functions/fun_fw5.jl")

seed!(4141)

nmodel = 50

#parameters
p = (0.1, 0.07, 3.2, 2.9, 0.5, 0.5, 0.15, 0.15, 2.5, 2.0, 0.3, 1.2)

tspan = (0.0, 599.999)

## data with process noise

data_fw5_noise = sim_fw5_noise(nmodel = nmodel, p = p, tspan = tspan, intense = false)

plt_sde = plot(data_fw5_noise[49], legend=false, title="Process noise", titlefontsize=9)


## data with higher level of process noise
## discard the results with low predator abundance (mean predator abundance < 0.5)

data_fw5_highnoise = sim_fw5_noise(nmodel = nmodel, p = p, tspan = tspan, intense = true)

plt_sde_high = plot(data_fw5_highnoise[49], legend=false, title="Intense process noise", titlefontsize=9)

for i in 1:nmodel
    writedlm("data/data_fw5_noise/data_fw5_noise_model$(i).csv", data_fw5_noise[i])
    writedlm("data/data_fw5_highnoise/data_fw5_highnoise_model$(i).csv", data_fw5_highnoise[i])
end

plot(plt_sde, plt_sde_high, layout=(1,:))
