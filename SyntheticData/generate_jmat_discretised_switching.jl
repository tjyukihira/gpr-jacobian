### Scripts for making discretised theoretical Jacobian of continuous time models
using LinearAlgebra

include("functions/get_jmat_discrete.jl")
include("load_data_switching.jl")

vec_jmat_true_switching_noise = copy(vec_jmat_scaled_switching_noise)
vec_jmat_true_switching_highnoise = copy(vec_jmat_scaled_switching_highnoise)
vec_jmat_true_switching_highnoise_obs = copy(vec_jmat_scaled_switching_highnoise_obs)
vec_jmat_true_switching_highnoise_highobs = copy(vec_jmat_scaled_switching_highnoise_highobs)

for model in 1:nmodel
    jmat_switching_noise_temp = zeros(5*datasize, 5)
    jmat_switching_highnoise_temp = zeros(5*datasize, 5)
    jmat_switching_highnoise_obs_temp = zeros(5*datasize, 5)
    jmat_switching_highnoise_highobs_temp = zeros(5*datasize, 5)
    for i in 1:datasize
        vec_jmat_true_switching_noise[model][:,:,i] = get_jmat_discrete(vec_jmat_scaled_switching_noise[model][:,:,i])
        vec_jmat_true_switching_highnoise[model][:,:,i] = get_jmat_discrete(vec_jmat_scaled_switching_highnoise[model][:,:,i])
        vec_jmat_true_switching_highnoise_obs[model][:,:,i] = get_jmat_discrete(vec_jmat_scaled_switching_highnoise_obs[model][:,:,i])
        vec_jmat_true_switching_highnoise_highobs[model][:,:,i] = get_jmat_discrete(vec_jmat_scaled_switching_highnoise_highobs[model][:,:,i])

        jmat_switching_noise_temp[(1+5*(i-1)):5*i, :] = vec_jmat_true_switching_noise[model][:,:,i]
        jmat_switching_highnoise_temp[(1+5*(i-1)):5*i, :] = vec_jmat_true_switching_highnoise[model][:,:,i]
        jmat_switching_highnoise_obs_temp[(1+5*(i-1)):5*i, :] = vec_jmat_true_switching_highnoise_obs[model][:,:,i]
        jmat_switching_highnoise_highobs_temp[(1+5*(i-1)):5*i, :] = vec_jmat_true_switching_highnoise_highobs[model][:,:,i]
    end

    writedlm("data/data_switching_noise/jmat_true_switching_noise_model$model.csv", jmat_switching_noise_temp)
    writedlm("data/data_switching_highnoise/jmat_true_switching_highnoise_model$model.csv", jmat_switching_highnoise_temp)
    writedlm("data/data_switching_highnoise_obs/jmat_true_switching_highnoise_obs_model$model.csv", jmat_switching_highnoise_obs_temp)
    writedlm("data/data_switching_highnoise_highobs/jmat_true_switching_highnoise_highobs_model$model.csv", jmat_switching_highnoise_highobs_temp)
end

