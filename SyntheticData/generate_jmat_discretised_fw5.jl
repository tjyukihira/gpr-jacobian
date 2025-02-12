### Scripts for making discretised theoretical Jacobian of continuous time models
include("functions/get_jmat_discrete.jl")
include("load_data_fw5.jl")

vec_jmat_true_fw5_noise = copy(vec_jmat_scaled_fw5_noise)
vec_jmat_true_fw5_highnoise = copy(vec_jmat_scaled_fw5_highnoise)
vec_jmat_true_fw5_highnoise_obs = copy(vec_jmat_scaled_fw5_highnoise_obs)
vec_jmat_true_fw5_highnoise_highobs = copy(vec_jmat_scaled_fw5_highnoise_highobs)

for model in 1:nmodel
    jmat_fw5_noise_temp = zeros(5*datasize, 5)
    jmat_fw5_highnoise_temp = zeros(5*datasize, 5)
    jmat_fw5_highnoise_obs_temp = zeros(5*datasize, 5)
    jmat_fw5_highnoise_highobs_temp = zeros(5*datasize, 5)
    for i in 1:datasize
        vec_jmat_true_fw5_noise[model][:,:,i] = get_jmat_discrete(vec_jmat_scaled_fw5_noise[model][:,:,i])
        vec_jmat_true_fw5_highnoise[model][:,:,i] = get_jmat_discrete(vec_jmat_scaled_fw5_highnoise[model][:,:,i])
        vec_jmat_true_fw5_highnoise_obs[model][:,:,i] = get_jmat_discrete(vec_jmat_scaled_fw5_highnoise_obs[model][:,:,i])
        vec_jmat_true_fw5_highnoise_highobs[model][:,:,i] = get_jmat_discrete(vec_jmat_scaled_fw5_highnoise_highobs[model][:,:,i])

        jmat_fw5_noise_temp[(1+5*(i-1)):5*i, :] = vec_jmat_true_fw5_noise[model][:,:,i]
        jmat_fw5_highnoise_temp[(1+5*(i-1)):5*i, :] = vec_jmat_true_fw5_highnoise[model][:,:,i]
        jmat_fw5_highnoise_obs_temp[(1+5*(i-1)):5*i, :] = vec_jmat_true_fw5_highnoise_obs[model][:,:,i]
        jmat_fw5_highnoise_highobs_temp[(1+5*(i-1)):5*i, :] = vec_jmat_true_fw5_highnoise_highobs[model][:,:,i]
    end

    writedlm("data/data_fw5_noise/jmat_true_fw5_noise_model$model.csv", jmat_fw5_noise_temp)
    writedlm("data/data_fw5_highnoise/jmat_true_fw5_highnoise_model$model.csv", jmat_fw5_highnoise_temp)
    writedlm("data/data_fw5_highnoise_obs/jmat_true_fw5_highnoise_obs_model$model.csv", jmat_fw5_highnoise_obs_temp)
    writedlm("data/data_fw5_highnoise_highobs/jmat_true_fw5_highnoise_highobs_model$model.csv", jmat_fw5_highnoise_highobs_temp)
end

