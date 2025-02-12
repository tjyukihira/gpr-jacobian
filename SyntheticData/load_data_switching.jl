#load synthetic time series data 
using Distributions
using StatsBase: transform
#load synthetic data of discrete LV model
using StatsBase
using Plots
using DelimitedFiles

include("functions/fun_prey_switching.jl")

pref = vec(readdlm("./data/data_switching/prey_preference.csv"))

#library size
datasize = 100
nmodel = 50

data_switching_noise = Vector{Matrix{Float64}}(undef, nmodel)
data_switching_highnoise = copy(data_switching_noise)
data_switching_highnoise_obs = copy(data_switching_noise)
data_switching_highnoise_highobs = copy(data_switching_noise)
data_switching_noise_tr = copy(data_switching_noise)
data_switching_highnoise_tr = copy(data_switching_noise)
data_switching_highnoise_obs_tr = copy(data_switching_noise)
data_switching_highnoise_highobs_tr = copy(data_switching_noise)

sd_switching_noise = Vector{Vector{Float64}}(undef, nmodel)
sd_switching_highnoise = copy(sd_switching_noise)
sd_switching_highnoise_obs = copy(sd_switching_noise)
sd_switching_highnoise_highobs = copy(sd_switching_noise)

mu_switching_noise = copy(sd_switching_noise)
mu_switching_highnoise = copy(sd_switching_noise)
mu_switching_highnoise_obs = copy(sd_switching_noise)
mu_switching_highnoise_highobs = copy(sd_switching_noise)

vec_jmat_unscaled_switching_noise = Vector{Array{Float64}}(undef, nmodel)
vec_jmat_scaled_switching_noise = copy(vec_jmat_unscaled_switching_noise)
vec_jmat_unscaled_switching_highnoise = copy(vec_jmat_unscaled_switching_noise)
vec_jmat_scaled_switching_highnoise = copy(vec_jmat_unscaled_switching_noise)
vec_jmat_unscaled_switching_highnoise_obs = copy(vec_jmat_unscaled_switching_noise)
vec_jmat_scaled_switching_highnoise_obs = copy(vec_jmat_unscaled_switching_noise)
vec_jmat_unscaled_switching_highnoise_highobs = copy(vec_jmat_unscaled_switching_noise)
vec_jmat_scaled_switching_highnoise_highobs = copy(vec_jmat_unscaled_switching_noise)

for i in 1:nmodel
    data_switching_noise[i] = readdlm("./data/data_switching_noise/data_switching_noise_model$(i).csv")
    data_switching_highnoise[i] = readdlm("./data/data_switching_highnoise/data_switching_highnoise_model$(i).csv")
    
    sd_switching_noise[i] = std(data_switching_noise[i], dims=1)[1,:]
    sd_switching_highnoise[i] = std(data_switching_highnoise[i], dims=1)[1,:]
    mu_switching_noise[i] = mean(data_switching_noise[i], dims=1)[1,:]
    mu_switching_highnoise[i] = mean(data_switching_highnoise[i], dims=1)[1,:]
    
    #scaling training data
    data_switching_noise_tr[i] = standardize(ZScoreTransform, data_switching_noise[i], dims=1)
    data_switching_highnoise_tr[i] = standardize(ZScoreTransform, data_switching_highnoise[i], dims=1)
    
    #adding observational noise
    data_switching_highnoise_obs[i] = readdlm("./data/data_switching_highnoise_obs/data_switching_highnoise_obs_model$(i).csv")
    data_switching_highnoise_highobs[i] = readdlm("./data/data_switching_highnoise_highobs/data_switching_highnoise_highobs_model$(i).csv")

    sd_switching_highnoise_obs[i] = std(data_switching_highnoise_obs[i], dims=1)[1,:]
    sd_switching_highnoise_highobs[i] = std(data_switching_highnoise_highobs[i], dims=1)[1,:]

    mu_switching_highnoise_obs[i] = mean(data_switching_highnoise_obs[i], dims=1)[1,:]
    mu_switching_highnoise_highobs[i] = mean(data_switching_highnoise_highobs[i], dims=1)[1,:]
    
    data_switching_highnoise_obs_tr[i] = standardize(ZScoreTransform, data_switching_highnoise_obs[i], dims=1)
    data_switching_highnoise_highobs_tr[i] = standardize(ZScoreTransform, data_switching_highnoise_highobs[i], dims=1)
    
    # theoretical jacobian matrix
    vec_jmat_unscaled_switching_noise[i] = zeros(5,5,datasize)
    vec_jmat_scaled_switching_noise[i] = zeros(5,5,datasize)
    vec_jmat_unscaled_switching_highnoise[i] = zeros(5,5,datasize)
    vec_jmat_scaled_switching_highnoise[i] = zeros(5,5,datasize)
    vec_jmat_unscaled_switching_highnoise_obs[i] = zeros(5,5,datasize)
    vec_jmat_scaled_switching_highnoise_obs[i] = zeros(5,5,datasize)
    vec_jmat_unscaled_switching_highnoise_highobs[i] = zeros(5,5,datasize)
    vec_jmat_scaled_switching_highnoise_highobs[i] = zeros(5,5,datasize)

    #parameters (xp, yp, pref, xc, yc, C0, R0)
    p = (0.08, 1.7, pref[i], 0.15, 2.3, 0.5, 0.25)
    
    vec_jmat_unscaled_switching_noise[i] = jmat_switching(data_switching_noise[i], p)
    vec_jmat_unscaled_switching_highnoise[i] = jmat_switching(data_switching_highnoise[i], p)

    #scale jacobian matrix
    for j in 1:datasize
        vec_jmat_scaled_switching_noise[i][:,:,j] = (vec_jmat_unscaled_switching_noise[i][:,:,j].*(1 ./sd_switching_noise[i])).*sd_switching_noise[i]'
        vec_jmat_scaled_switching_highnoise[i][:,:,j] = (vec_jmat_unscaled_switching_highnoise[i][:,:,j].*(1 ./sd_switching_highnoise[i])).*sd_switching_highnoise[i]'
        vec_jmat_scaled_switching_highnoise_obs[i][:,:,j] = (vec_jmat_unscaled_switching_highnoise[i][:,:,j].*(1 ./sd_switching_highnoise_obs[i])).*sd_switching_highnoise_obs[i]'
        vec_jmat_scaled_switching_highnoise_highobs[i][:,:,j] = (vec_jmat_unscaled_switching_highnoise[i][:,:,j].*(1 ./sd_switching_highnoise_highobs[i])).*sd_switching_highnoise_highobs[i]'
    end
end

#plt=plot(data_switching_tr[1])
pltnoise=plot(data_switching_noise[1], label=permutedims(["P1", "P2", "C1", "C2", "R"]), legendfontsize=5, title="process noise", titlefontsize=7)
pltintnoise=plot(data_switching_highnoise[1], label=permutedims(["P1", "P2", "C1", "C2", "R"]), legendfontsize=5, title="high process noise", titlefontsize=7)
pltintnoise_obs=plot(data_switching_highnoise_obs[1], label=permutedims(["P1", "P2", "C1", "C2", "R"]), legendfontsize=5, title="high process and observational noise", titlefontsize=7)
pltintnoise_intobs=plot(data_switching_highnoise_highobs[1], label=permutedims(["P1", "P2", "C1", "C2", "R"]), legendfontsize=5, title="high process and high observational noise", titlefontsize=7)

plot(pltnoise, pltintnoise, pltintnoise_obs, pltintnoise_intobs, layout=(2,:))
