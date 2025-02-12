### Scripts for making theoretical Jacobian and data with observational noise
using Distributions
using StatsBase: transform
using StatsBase
using Plots
using DelimitedFiles
import Random.shuffle
import Random.seed!

include("functions/fun_prey_switching.jl")

seed!(8888)

pref = vec(readdlm("data/data_switching/prey_preference.csv"))

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
    data_switching_highnoise_obs[i] = data_switching_highnoise[i] + data_switching_highnoise[i].*rand(Normal(0, 0.1), datasize, 5)
    data_switching_highnoise_highobs[i] = data_switching_highnoise[i] + data_switching_highnoise[i].*rand(Normal(0, 0.2), datasize, 5)

    sd_switching_highnoise_obs[i] = std(data_switching_highnoise_obs[i], dims=1)[1,:]
    sd_switching_highnoise_highobs[i] = std(data_switching_highnoise_highobs[i], dims=1)[1,:]
    mu_switching_highnoise_obs[i] = mean(data_switching_highnoise_obs[i], dims=1)[1,:]
    mu_switching_highnoise_highobs[i] = mean(data_switching_highnoise_highobs[i], dims=1)[1,:]

    data_switching_highnoise_obs_tr[i] = standardize(ZScoreTransform, data_switching_highnoise_obs[i], dims=1)
    data_switching_highnoise_highobs_tr[i] = standardize(ZScoreTransform, data_switching_highnoise_highobs[i], dims=1)
    
    #theoretical jacobian matrix
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

    writedlm("./data/data_switching_highnoise_obs/data_switching_highnoise_obs_model$(i).csv", data_switching_highnoise_obs[i])
    writedlm("./data/data_switching_highnoise_highobs/data_switching_highnoise_highobs_model$(i).csv", data_switching_highnoise_highobs[i])
end

for i in 1:nmodel
    jmat_switching_noise_temp = zeros(5*datasize, 5)
    jmat_switching_highnoise_temp = zeros(5*datasize, 5)
    jmat_switching_highnoise_obs_temp = zeros(5*datasize, 5)
    jmat_switching_highnoise_highobs_temp = zeros(5*datasize, 5)
    for j in 1:datasize
        jmat_switching_noise_temp[(1+5*(j-1)):5*j, :] = vec_jmat_scaled_switching_noise[i][:,:,j]
        jmat_switching_highnoise_temp[(1+5*(j-1)):5*j, :] = vec_jmat_scaled_switching_highnoise[i][:,:,j]
        jmat_switching_highnoise_obs_temp[(1+5*(j-1)):5*j, :] = vec_jmat_scaled_switching_highnoise_obs[i][:,:,j]
        jmat_switching_highnoise_highobs_temp[(1+5*(j-1)):5*j, :] = vec_jmat_scaled_switching_highnoise_highobs[i][:,:,j]
    end
    writedlm("./data/data_switching_noise/jmat_switching_noise_model$i.csv", jmat_switching_noise_temp)
    writedlm("./data/data_switching_highnoise/jmat_switching_highnoise_model$i.csv", jmat_switching_highnoise_temp)
    writedlm("./data/data_switching_highnoise_obs/jmat_switching_highnoise_obs_model$i.csv", jmat_switching_highnoise_obs_temp)
    writedlm("./data/data_switching_highnoise_highobs/jmat_switching_highnoise_highobs_model$i.csv", jmat_switching_highnoise_highobs_temp)
end

# plot examplary time series for each condition
pltnoise=plot(data_switching_noise_tr[1], title="process noise", titlefontsize=7)
pltintnoise=plot(data_switching_highnoise_tr[1], title="high process noise", titlefontsize=7)
pltintnoise_obs=plot(data_switching_highnoise_obs_tr[1], title="high process and observational noise", titlefontsize=7)
pltintnoise_intobs=plot(data_switching_highnoise_highobs_tr[1], title="high process and high observational noise", titlefontsize=7)

plot(pltnoise, pltintnoise, pltintnoise_obs, pltintnoise_intobs, layout=(2,:))
