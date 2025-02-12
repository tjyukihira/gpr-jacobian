### Scripts for making theoretical Jacobian matrix and data with observational noise

# load packages
using Distributions
using StatsBase: transform
using StatsBase
using Plots
using DelimitedFiles
import Random.shuffle
import Random.seed!

# set seed for replication
seed!(2323)

#library size
npop = 10
datasize = 100
nmodel = 50

data_logistic_10_noise = Vector{Matrix{Float64}}(undef, nmodel)
data_logistic_10_highnoise = copy(data_logistic_10_noise)

data_logistic_10_highnoise_obs = copy(data_logistic_10_noise)
data_logistic_10_highnoise_highobs = copy(data_logistic_10_noise)

data_logistic_10_noise_tr = copy(data_logistic_10_noise)
data_logistic_10_highnoise_tr = copy(data_logistic_10_noise)

data_logistic_10_highnoise_obs_tr = copy(data_logistic_10_noise)
data_logistic_10_highnoise_highobs_tr = copy(data_logistic_10_noise)

sd_logistic_10_noise = Vector{Vector{Float64}}(undef, nmodel)
sd_logistic_10_highnoise = copy(sd_logistic_10_noise)
sd_logistic_10_highnoise_obs = copy(sd_logistic_10_noise)
sd_logistic_10_highnoise_highobs = copy(sd_logistic_10_noise)

mu_logistic_10_noise = copy(sd_logistic_10_noise)
mu_logistic_10_highnoise = copy(sd_logistic_10_noise)
mu_logistic_10_highnoise_obs = copy(sd_logistic_10_noise)
mu_logistic_10_highnoise_highobs = copy(sd_logistic_10_noise)

vec_jmat_unscaled_logistic_10_noise = Vector{Array{Float64}}(undef, nmodel)
vec_jmat_scaled_logistic_10_noise = copy(vec_jmat_unscaled_logistic_10_noise)
vec_jmat_unscaled_logistic_10_highnoise = copy(vec_jmat_unscaled_logistic_10_noise)
vec_jmat_scaled_logistic_10_highnoise = copy(vec_jmat_unscaled_logistic_10_noise)

vec_jmat_unscaled_logistic_10_highnoise_obs = copy(vec_jmat_unscaled_logistic_10_noise)
vec_jmat_scaled_logistic_10_highnoise_obs = copy(vec_jmat_unscaled_logistic_10_noise)
vec_jmat_unscaled_logistic_10_highnoise_highobs = copy(vec_jmat_unscaled_logistic_10_noise)
vec_jmat_scaled_logistic_10_highnoise_highobs = copy(vec_jmat_unscaled_logistic_10_noise)

for i in 1:nmodel
    data_logistic_10_noise[i] = readdlm("./data/data_logistic_10_noise/data_logistic_10_noise_model$(i).csv")
    data_logistic_10_highnoise[i] = readdlm("./data/data_logistic_10_highnoise/data_logistic_10_highnoise_model$(i).csv")
    
    sd_logistic_10_noise[i] = std(data_logistic_10_noise[i], dims=1)[1,:]
    sd_logistic_10_highnoise[i] = std(data_logistic_10_highnoise[i], dims=1)[1,:]
    
    mu_logistic_10_noise[i] = mean(data_logistic_10_noise[i], dims=1)[1,:]
    mu_logistic_10_highnoise[i] = mean(data_logistic_10_highnoise[i], dims=1)[1,:]
    
    #scaling training data
    data_logistic_10_noise_tr[i] = standardize(ZScoreTransform, data_logistic_10_noise[i], dims=1)
    data_logistic_10_highnoise_tr[i] = standardize(ZScoreTransform, data_logistic_10_highnoise[i], dims=1)
    
    #adding observational noise
    data_logistic_10_highnoise_obs[i] = data_logistic_10_highnoise[i] + data_logistic_10_highnoise[i].*rand(Normal(0, 0.1), datasize, npop)
    data_logistic_10_highnoise_highobs[i] = data_logistic_10_highnoise[i] + data_logistic_10_highnoise[i].*rand(Normal(0, 0.1), datasize, npop)

    #mean and sd for data with observational noise
    sd_logistic_10_highnoise_obs[i] = std(data_logistic_10_highnoise_obs[i], dims=1)[1,:]
    sd_logistic_10_highnoise_highobs[i] = std(data_logistic_10_highnoise_highobs[i], dims=1)[1,:]
    mu_logistic_10_highnoise_obs[i] = mean(data_logistic_10_highnoise_obs[i], dims=1)[1,:]
    mu_logistic_10_highnoise_highobs[i] = mean(data_logistic_10_highnoise_highobs[i], dims=1)[1,:]
    
    data_logistic_10_highnoise_obs_tr[i] = standardize(ZScoreTransform, data_logistic_10_highnoise_obs[i], dims=1)
    data_logistic_10_highnoise_highobs_tr[i] = standardize(ZScoreTransform, data_logistic_10_highnoise_highobs[i], dims=1)
    
    #theoretical jacobian matrix
    jmat_logistic_10_noise = readdlm("./data/data_logistic_10_noise/jmat_logistic_10_noise_model$(i).csv")
    jmat_logistic_10_highnoise = readdlm("./data/data_logistic_10_highnoise/jmat_logistic_10_highnoise_model$(i).csv")

    vec_jmat_unscaled_logistic_10_noise[i] = zeros(npop,npop,datasize)
    vec_jmat_scaled_logistic_10_noise[i] = zeros(npop,npop,datasize)
    vec_jmat_unscaled_logistic_10_highnoise[i] = zeros(npop,npop,datasize)
    vec_jmat_scaled_logistic_10_highnoise[i] = zeros(npop,npop,datasize)
    vec_jmat_unscaled_logistic_10_highnoise_obs[i] = zeros(npop,npop,datasize)
    vec_jmat_scaled_logistic_10_highnoise_obs[i] = zeros(npop,npop,datasize)
    vec_jmat_unscaled_logistic_10_highnoise_highobs[i] = zeros(npop,npop,datasize)
    vec_jmat_scaled_logistic_10_highnoise_highobs[i] = zeros(npop,npop,datasize)

    for j in 1:(datasize)
        vec_jmat_unscaled_logistic_10_noise[i][:,:,j] = jmat_logistic_10_noise[(1+(j-1)*npop):(j*npop), 1:npop]
        vec_jmat_unscaled_logistic_10_highnoise[i][:,:,j] = jmat_logistic_10_highnoise[(1+(j-1)*npop):(j*npop), 1:npop]
        
        vec_jmat_scaled_logistic_10_noise[i][:,:,j] = (vec_jmat_unscaled_logistic_10_noise[i][:,:,j].*(1 ./sd_logistic_10_noise[i])).*sd_logistic_10_noise[i]'
        vec_jmat_scaled_logistic_10_highnoise[i][:,:,j] = (vec_jmat_unscaled_logistic_10_highnoise[i][:,:,j].*(1 ./sd_logistic_10_highnoise[i])).*sd_logistic_10_highnoise[i]'

        vec_jmat_scaled_logistic_10_highnoise_obs[i][:,:,j] = (vec_jmat_unscaled_logistic_10_highnoise[i][:,:,j].*(1 ./sd_logistic_10_highnoise_obs[i])).*sd_logistic_10_highnoise_obs[i]'
        vec_jmat_scaled_logistic_10_highnoise_highobs[i][:,:,j] = (vec_jmat_unscaled_logistic_10_highnoise[i][:,:,j].*(1 ./sd_logistic_10_highnoise_highobs[i])).*sd_logistic_10_highnoise_highobs[i]'
    end

    writedlm("./data/data_logistic_10_highnoise_obs/data_logistic_10_highnoise_obs_model$(i).csv", data_logistic_10_highnoise_obs[i])
    writedlm("./data/data_logistic_10_highnoise_highobs/data_logistic_10_highnoise_highobs_model$(i).csv", data_logistic_10_highnoise_highobs[i])
end

# save theoretical jacobian matrix
for i in 1:nmodel
    jmat_logistic_10_noise_temp = zeros(npop*datasize, npop)
    jmat_logistic_10_highnoise_temp = zeros(npop*datasize, npop)
    jmat_logistic_10_highnoise_obs_temp = zeros(npop*datasize, npop)
    jmat_logistic_10_highnoise_highobs_temp = zeros(npop*datasize, npop)
    for j in 1:datasize
        jmat_logistic_10_noise_temp[(1+npop*(j-1)):npop*j, :] = vec_jmat_scaled_logistic_10_noise[i][:,:,j]
        jmat_logistic_10_highnoise_temp[(1+npop*(j-1)):npop*j, :] = vec_jmat_scaled_logistic_10_highnoise[i][:,:,j]
        jmat_logistic_10_highnoise_obs_temp[(1+npop*(j-1)):npop*j, :] = vec_jmat_scaled_logistic_10_highnoise_obs[i][:,:,j]
        jmat_logistic_10_highnoise_highobs_temp[(1+npop*(j-1)):npop*j, :] = vec_jmat_scaled_logistic_10_highnoise_highobs[i][:,:,j]
    end
    
    writedlm("./data/data_logistic_10_noise/jmat_logistic_10_noise_scaled_model$i.csv", jmat_logistic_10_noise_temp)
    writedlm("./data/data_logistic_10_highnoise/jmat_logistic_10_highnoise_scaled_model$i.csv", jmat_logistic_10_highnoise_temp)
    writedlm("./data/data_logistic_10_highnoise_obs/jmat_logistic_10_highnoise_obs_scaled_model$i.csv", jmat_logistic_10_highnoise_obs_temp)
    writedlm("./data/data_logistic_10_highnoise_highobs/jmat_logistic_10_highnoise_highobs_scaled_model$i.csv", jmat_logistic_10_highnoise_highobs_temp)
end

# plot an example time series for each noise condition
pltnoise=plot(data_logistic_10_noise[1], title="process noise", titlefontsize=7)
pltintnoise=plot(data_logistic_10_highnoise[1], title="high process noise", titlefontsize=7)
pltintnoise_obs=plot(data_logistic_10_highnoise_obs[1], title="high process and observational noise", titlefontsize=7)
pltintnoise_intobs=plot(data_logistic_10_highnoise_highobs[1], title="high process and high observational noise", titlefontsize=7)

plot(pltnoise, pltintnoise, pltintnoise_obs, pltintnoise_intobs, layout=(2,:))
