#load synthetic time series data 
using Distributions
using StatsBase: transform
#load synthetic data of discrete LV model
using StatsBase
using Plots
using DelimitedFiles

include("functions/fun_fw5.jl")

p = (0.1, 0.07, 3.2, 2.9, 0.5, 0.5, 0.15, 0.15, 2.5, 2.0, 0.3, 1.2)

#library size
datasize = 100
nmodel = 50

data_fw5_noise = Vector{Matrix{Float64}}(undef, nmodel)
data_fw5_highnoise = copy(data_fw5_noise)
data_fw5_highnoise_obs = copy(data_fw5_noise)
data_fw5_highnoise_highobs = copy(data_fw5_noise)
data_fw5_noise_tr = copy(data_fw5_noise)
data_fw5_highnoise_tr = copy(data_fw5_noise)
data_fw5_highnoise_obs_tr = copy(data_fw5_noise)
data_fw5_highnoise_highobs_tr = copy(data_fw5_noise)

sd_fw5_noise = Vector{Vector{Float64}}(undef, nmodel)
sd_fw5_highnoise = copy(sd_fw5_noise)
sd_fw5_highnoise_obs = copy(sd_fw5_noise)
sd_fw5_highnoise_highobs = copy(sd_fw5_noise)

mu_fw5_noise = copy(sd_fw5_noise)
mu_fw5_highnoise = copy(sd_fw5_noise)
mu_fw5_highnoise_obs = copy(sd_fw5_noise)
mu_fw5_highnoise_highobs = copy(sd_fw5_noise)

vec_jmat_unscaled_fw5_noise = Vector{Array{Float64}}(undef, nmodel)
vec_jmat_scaled_fw5_noise = copy(vec_jmat_unscaled_fw5_noise)
vec_jmat_unscaled_fw5_highnoise = copy(vec_jmat_unscaled_fw5_noise)
vec_jmat_scaled_fw5_highnoise = copy(vec_jmat_unscaled_fw5_noise)
vec_jmat_unscaled_fw5_highnoise_obs = copy(vec_jmat_unscaled_fw5_noise)
vec_jmat_scaled_fw5_highnoise_obs = copy(vec_jmat_unscaled_fw5_noise)
vec_jmat_unscaled_fw5_highnoise_highobs = copy(vec_jmat_unscaled_fw5_noise)
vec_jmat_scaled_fw5_highnoise_highobs = copy(vec_jmat_unscaled_fw5_noise)

for i in 1:nmodel
    data_fw5_noise[i] = readdlm("./data/data_fw5_noise/data_fw5_noise_model$(i).csv")
    data_fw5_highnoise[i] = readdlm("./data/data_fw5_highnoise/data_fw5_highnoise_model$(i).csv")
    
    sd_fw5_noise[i] = std(data_fw5_noise[i], dims=1)[1,:]
    sd_fw5_highnoise[i] = std(data_fw5_highnoise[i], dims=1)[1,:]
    mu_fw5_noise[i] = mean(data_fw5_noise[i], dims=1)[1,:]
    mu_fw5_highnoise[i] = mean(data_fw5_highnoise[i], dims=1)[1,:]
    
    #scaling training data
    data_fw5_noise_tr[i] = standardize(ZScoreTransform, data_fw5_noise[i], dims=1)
    data_fw5_highnoise_tr[i] = standardize(ZScoreTransform, data_fw5_highnoise[i], dims=1)
    
    #adding observational noise
    data_fw5_highnoise_obs[i] = readdlm("./data/data_fw5_highnoise_obs/data_fw5_highnoise_obs_model$(i).csv")
    data_fw5_highnoise_highobs[i] = readdlm("./data/data_fw5_highnoise_highobs/data_fw5_highnoise_highobs_model$(i).csv")

    sd_fw5_highnoise_obs[i] = std(data_fw5_highnoise_obs[i], dims=1)[1,:]
    sd_fw5_highnoise_highobs[i] = std(data_fw5_highnoise_highobs[i], dims=1)[1,:]

    mu_fw5_highnoise_obs[i] = mean(data_fw5_highnoise_obs[i], dims=1)[1,:]
    mu_fw5_highnoise_highobs[i] = mean(data_fw5_highnoise_highobs[i], dims=1)[1,:]
    
    data_fw5_highnoise_obs_tr[i] = standardize(ZScoreTransform, data_fw5_highnoise_obs[i], dims=1)
    data_fw5_highnoise_highobs_tr[i] = standardize(ZScoreTransform, data_fw5_highnoise_highobs[i], dims=1)
    
    # theoretical jacobian matrix
    vec_jmat_unscaled_fw5_noise[i] = zeros(5,5,datasize)
    vec_jmat_scaled_fw5_noise[i] = zeros(5,5,datasize)
    vec_jmat_unscaled_fw5_highnoise[i] = zeros(5,5,datasize)
    vec_jmat_scaled_fw5_highnoise[i] = zeros(5,5,datasize)
    vec_jmat_unscaled_fw5_highnoise_obs[i] = zeros(5,5,datasize)
    vec_jmat_scaled_fw5_highnoise_obs[i] = zeros(5,5,datasize)
    vec_jmat_unscaled_fw5_highnoise_highobs[i] = zeros(5,5,datasize)
    vec_jmat_scaled_fw5_highnoise_highobs[i] = zeros(5,5,datasize)
    
    vec_jmat_unscaled_fw5_noise[i] = jmat_fw5(data_fw5_noise[i], p)
    vec_jmat_unscaled_fw5_highnoise[i] = jmat_fw5(data_fw5_highnoise[i], p)

    #scale jacobian matrix
    for j in 1:datasize
        vec_jmat_scaled_fw5_noise[i][:,:,j] = (vec_jmat_unscaled_fw5_noise[i][:,:,j].*(1 ./sd_fw5_noise[i])).*sd_fw5_noise[i]'
        vec_jmat_scaled_fw5_highnoise[i][:,:,j] = (vec_jmat_unscaled_fw5_highnoise[i][:,:,j].*(1 ./sd_fw5_highnoise[i])).*sd_fw5_highnoise[i]'
        vec_jmat_scaled_fw5_highnoise_obs[i][:,:,j] = (vec_jmat_unscaled_fw5_highnoise[i][:,:,j].*(1 ./sd_fw5_highnoise_obs[i])).*sd_fw5_highnoise_obs[i]'
        vec_jmat_scaled_fw5_highnoise_highobs[i][:,:,j] = (vec_jmat_unscaled_fw5_highnoise[i][:,:,j].*(1 ./sd_fw5_highnoise_highobs[i])).*sd_fw5_highnoise_highobs[i]'
    end
end

#plt=plot(data_fw5_tr[1])
pltnoise=plot(data_fw5_noise[1], label=permutedims(["P1", "P2", "C1", "C2", "R"]), legendfontsize=5, title="process noise", titlefontsize=7)
pltintnoise=plot(data_fw5_highnoise[1], label=permutedims(["P1", "P2", "C1", "C2", "R"]), legendfontsize=5, title="high process noise", titlefontsize=7)
pltintnoise_obs=plot(data_fw5_highnoise_obs[1], label=permutedims(["P1", "P2", "C1", "C2", "R"]), legendfontsize=5, title="high process and observational noise", titlefontsize=7)
pltintnoise_intobs=plot(data_fw5_highnoise_highobs[1], label=permutedims(["P1", "P2", "C1", "C2", "R"]), legendfontsize=5, title="high process and high observational noise", titlefontsize=7)

plot(pltnoise, pltintnoise, pltintnoise_obs, pltintnoise_intobs, layout=(2,:))
