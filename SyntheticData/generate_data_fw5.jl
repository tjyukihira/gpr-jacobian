### Scripts for making theoretical Jacobian and data with observational data

# load packages
using Distributions
using StatsBase: transform
using StatsBase
using Plots
using DelimitedFiles
import Random.shuffle
import Random.seed!

include("functions/fun_fw5.jl")

seed!(1414)

p = (0.1, 0.07, 3.2, 2.9, 0.5, 0.5, 0.15, 0.15, 2.5, 2.0, 0.3, 1.2)

#library size
datasize = 100
nmodel = 50

data_fw5_noise = Vector{Matrix{Float64}}(undef, nmodel)
data_fw5_highnoise = copy(data_fw5_noise)
data_fw5_noise_obs = copy(data_fw5_noise)
data_fw5_noise_highobs = copy(data_fw5_noise)
data_fw5_highnoise_obs = copy(data_fw5_noise)
data_fw5_highnoise_highobs = copy(data_fw5_noise)
data_fw5_noise_tr = copy(data_fw5_noise)
data_fw5_highnoise_tr = copy(data_fw5_noise)
#data_fw5_noise_obs_tr = copy(data_fw5_noise)
#data_fw5_noise_highobs_tr = copy(data_fw5_noise)
data_fw5_highnoise_obs_tr = copy(data_fw5_noise)
data_fw5_highnoise_highobs_tr = copy(data_fw5_noise)

sd_fw5_noise = Vector{Vector{Float64}}(undef, nmodel)
sd_fw5_highnoise = copy(sd_fw5_noise)
#sd_fw5_noise_obs = copy(sd_fw5_noise)
#sd_fw5_noise_highobs = copy(sd_fw5_noise)
sd_fw5_highnoise_obs = copy(sd_fw5_noise)
sd_fw5_highnoise_highobs = copy(sd_fw5_noise)

mu_fw5_noise = copy(sd_fw5_noise)
mu_fw5_highnoise = copy(sd_fw5_noise)
#mu_fw5_noise_obs = copy(sd_fw5_noise)
#mu_fw5_noise_highobs = copy(sd_fw5_noise)
mu_fw5_highnoise_obs = copy(sd_fw5_noise)
mu_fw5_highnoise_highobs = copy(sd_fw5_noise)

vec_jmat_unscaled_fw5_noise = Vector{Array{Float64}}(undef, nmodel)
vec_jmat_scaled_fw5_noise = copy(vec_jmat_unscaled_fw5_noise)
vec_jmat_unscaled_fw5_highnoise = copy(vec_jmat_unscaled_fw5_noise)
vec_jmat_scaled_fw5_highnoise = copy(vec_jmat_unscaled_fw5_noise)

#vec_jmat_unscaled_fw5_noise_obs = copy(vec_jmat_unscaled_fw5_noise)
#vec_jmat_scaled_fw5_noise_obs = copy(vec_jmat_unscaled_fw5_noise)
#vec_jmat_unscaled_fw5_noise_highobs = copy(vec_jmat_unscaled_fw5_noise)
#vec_jmat_scaled_fw5_noise_highobs = copy(vec_jmat_unscaled_fw5_noise)
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
    #data_fw5_noise_obs[i] = data_fw5_noise[i] + data_fw5_noise[i].*rand(Normal(0, 0.1), datasize, 5)
    #data_fw5_noise_highobs[i] = data_fw5_noise[i] + data_fw5_noise[i].*rand(Normal(0, 0.2), datasize, 5)
    data_fw5_highnoise_obs[i] = data_fw5_highnoise[i] + data_fw5_highnoise[i].*rand(Normal(0, 0.1), datasize, 5)
    data_fw5_highnoise_highobs[i] = data_fw5_highnoise[i] + data_fw5_highnoise[i].*rand(Normal(0, 0.2), datasize, 5)

    #mean and sd for data with observational noise
    #sd_fw5_noise_obs[i] = std(data_fw5_noise_obs[i], dims=1)[1,:]
    #sd_fw5_noise_highobs[i] = std(data_fw5_noise_highobs[i], dims=1)[1,:]
    sd_fw5_highnoise_obs[i] = std(data_fw5_highnoise_obs[i], dims=1)[1,:]
    sd_fw5_highnoise_highobs[i] = std(data_fw5_highnoise_highobs[i], dims=1)[1,:]
    #mu_fw5_noise_obs[i] = mean(data_fw5_noise_obs[i], dims=1)[1,:]
    mu_fw5_highnoise_obs[i] = mean(data_fw5_highnoise_obs[i], dims=1)[1,:]
    mu_fw5_highnoise_highobs[i] = mean(data_fw5_highnoise_highobs[i], dims=1)[1,:]

    #data_fw5_noise_obs_tr[i] = standardize(ZScoreTransform, data_fw5_noise_obs[i], dims=1)
    #data_fw5_noise_highobs_tr[i] = standardize(ZScoreTransform, data_fw5_noise_highobs[i], dims=1)
    data_fw5_highnoise_obs_tr[i] = standardize(ZScoreTransform, data_fw5_highnoise_obs[i], dims=1)
    data_fw5_highnoise_highobs_tr[i] = standardize(ZScoreTransform, data_fw5_highnoise_highobs[i], dims=1)
    
    #theoretical jacobian matrix
    vec_jmat_unscaled_fw5_noise[i] = zeros(5,5,datasize)
    vec_jmat_scaled_fw5_noise[i] = zeros(5,5,datasize)
    vec_jmat_unscaled_fw5_highnoise[i] = zeros(5,5,datasize)
    vec_jmat_scaled_fw5_highnoise[i] = zeros(5,5,datasize)
    
    #vec_jmat_unscaled_fw5_noise_obs[i] = zeros(5,5,datasize)
    #vec_jmat_scaled_fw5_noise_obs[i] = zeros(5,5,datasize)
    #vec_jmat_unscaled_fw5_noise_highobs[i] = zeros(5,5,datasize)
    #vec_jmat_scaled_fw5_noise_highobs[i] = zeros(5,5,datasize)
    vec_jmat_unscaled_fw5_highnoise_obs[i] = zeros(5,5,datasize)
    vec_jmat_scaled_fw5_highnoise_obs[i] = zeros(5,5,datasize)
    vec_jmat_unscaled_fw5_highnoise_highobs[i] = zeros(5,5,datasize)
    vec_jmat_scaled_fw5_highnoise_highobs[i] = zeros(5,5,datasize)
    
    # compute theoretical Jacobian
    vec_jmat_unscaled_fw5_noise[i] = jmat_fw5(data_fw5_noise[i], p)
    vec_jmat_unscaled_fw5_highnoise[i] = jmat_fw5(data_fw5_highnoise[i], p)
    
    #scale jacobian matrix
    for j in 1:datasize
        vec_jmat_scaled_fw5_noise[i][:,:,j] = (vec_jmat_unscaled_fw5_noise[i][:,:,j].*(1 ./sd_fw5_noise[i])).*sd_fw5_noise[i]'
        vec_jmat_scaled_fw5_highnoise[i][:,:,j] = (vec_jmat_unscaled_fw5_highnoise[i][:,:,j].*(1 ./sd_fw5_highnoise[i])).*sd_fw5_highnoise[i]'
        
        #vec_jmat_scaled_fw5_noise_obs[i][:,:,j] = (vec_jmat_unscaled_fw5_noise[i][:,:,j].*(1 ./sd_fw5_noise_obs[i])).*sd_fw5_noise_obs[i]'
        #vec_jmat_scaled_fw5_noise_highobs[i][:,:,j] = (vec_jmat_unscaled_fw5_noise[i][:,:,j].*(1 ./sd_fw5_noise_highobs[i])).*sd_fw5_noise_highobs[i]'
        vec_jmat_scaled_fw5_highnoise_obs[i][:,:,j] = (vec_jmat_unscaled_fw5_highnoise[i][:,:,j].*(1 ./sd_fw5_highnoise_obs[i])).*sd_fw5_highnoise_obs[i]'
        vec_jmat_scaled_fw5_highnoise_highobs[i][:,:,j] = (vec_jmat_unscaled_fw5_highnoise[i][:,:,j].*(1 ./sd_fw5_highnoise_highobs[i])).*sd_fw5_highnoise_highobs[i]'
    end

    #writedlm("./data/data_fw5_noise/data_fw5_noise_model$(i).csv", data_fw5_noise[i])
    #writedlm("./data/data_fw5_highnoise/data_fw5_highnoise_model$(i).csv", data_fw5_highnoise[i])
    #writedlm("./data/data_fw5_noise_obs/data_fw5_noise_obs_model$(i).csv", data_fw5_noise_obs[i])
    #writedlm("./data/data_fw5_noise_highobs/data_fw5_noise_highobs_model$(i).csv", data_fw5_noise_highobs[i])
    writedlm("./data/data_fw5_highnoise_obs/data_fw5_highnoise_obs_model$(i).csv", data_fw5_highnoise_obs[i])
    writedlm("./data/data_fw5_highnoise_highobs/data_fw5_highnoise_highobs_model$(i).csv", data_fw5_highnoise_highobs[i])
end

# save theoretical Jacobian data
for i in 1:nmodel
    jmat_fw5_noise_temp = zeros(5*datasize, 5)
    jmat_fw5_highnoise_temp = zeros(5*datasize, 5)
    
    #jmat_fw5_noise_obs_temp = zeros(5*datasize, 5)
    #jmat_fw5_noise_highobs_temp = zeros(5*datasize, 5)
    jmat_fw5_highnoise_obs_temp = zeros(5*datasize, 5)
    jmat_fw5_highnoise_highobs_temp = zeros(5*datasize, 5)
    for j in 1:datasize
        jmat_fw5_noise_temp[(1+5*(j-1)):5*j, :] = vec_jmat_scaled_fw5_noise[i][:,:,j]
        jmat_fw5_highnoise_temp[(1+5*(j-1)):5*j, :] = vec_jmat_scaled_fw5_highnoise[i][:,:,j]
        
        #jmat_fw5_noise_obs_temp[(1+5*(j-1)):5*j, :] = vec_jmat_scaled_fw5_noise_obs[i][:,:,j]
        #jmat_fw5_noise_highobs_temp[(1+5*(j-1)):5*j, :] = vec_jmat_scaled_fw5_noise_highobs[i][:,:,j]
        jmat_fw5_highnoise_obs_temp[(1+5*(j-1)):5*j, :] = vec_jmat_scaled_fw5_highnoise_obs[i][:,:,j]
        jmat_fw5_highnoise_highobs_temp[(1+5*(j-1)):5*j, :] = vec_jmat_scaled_fw5_highnoise_highobs[i][:,:,j]
    end
    writedlm("./data/data_fw5_noise/jmat_fw5_noise_model$i.csv", jmat_fw5_noise_temp)
    writedlm("./data/data_fw5_highnoise/jmat_fw5_highnoise_model$i.csv", jmat_fw5_highnoise_temp)
    #writedlm("./data/data_fw5_noise_obs/jmat_fw5_noise_obs_model$i.csv", jmat_fw5_noise_obs_temp)
    #writedlm("./data/data_fw5_noise_highobs/jmat_fw5_noise_highobs_model$i.csv", jmat_fw5_noise_highobs_temp)
    writedlm("./data/data_fw5_highnoise_obs/jmat_fw5_highnoise_obs_model$i.csv", jmat_fw5_highnoise_obs_temp)
    writedlm("./data/data_fw5_highnoise_highobs/jmat_fw5_highnoise_highobs_model$i.csv", jmat_fw5_highnoise_highobs_temp)
end

#plot examplary time series for each condition
pltnoise=plot(data_fw5_noise[1], label=permutedims(["P1", "P2", "C1", "C2", "R"]), legendfontsize=5, title="process noise", titlefontsize=7)
pltintnoise=plot(data_fw5_highnoise[1], label=permutedims(["P1", "P2", "C1", "C2", "R"]), legendfontsize=5, title="high process noise", titlefontsize=7)
#pltnoise_obs=plot(data_fw5_noise_obs[1], label=permutedims(["P1", "P2", "C1", "C2", "R"]), legendfontsize=5, title="process and observational noise", titlefontsize=7)
#pltnoise_intobs=plot(data_fw5_noise_highobs_tr[1], title="process and high observational noise", titlefontsize=7)
pltintnoise_obs=plot(data_fw5_highnoise_obs[1], label=permutedims(["P1", "P2", "C1", "C2", "R"]), title="high process and observational noise", titlefontsize=7)
pltintnoise_intobs=plot(data_fw5_highnoise_highobs[1], label=permutedims(["P1", "P2", "C1", "C2", "R"]), legendfontsize=5, title="high process and high observational noise", titlefontsize=7)

plot(pltnoise, pltintnoise, pltintnoise_obs, pltintnoise_intobs, layout=(2,:))
