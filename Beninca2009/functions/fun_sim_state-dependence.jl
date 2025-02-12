### functions for analyses of state-dependence
function sim_mean_strength_switching(vec_data, vec_res, pref; mu, sd)
    nmodel = length(vec_data)
    datasize = size(vec_data[1], 1)
    npop = size(vec_data[1], 2)

    vec_mean_strength = Vector{Vector{Float64}}(undef, nmodel)
    vec_mean_strength_plus2sd = copy(vec_mean_strength)
    vec_mean_strength_min2sd = copy(vec_mean_strength)
    vec_mean_strength_true = copy(vec_mean_strength)

    vec_jmat_gpr = Vector{Array{Float64}}(undef, nmodel)
    vec_jmat_true = copy(vec_jmat_gpr)

    for model in 1:nmodel

        mean_strength = zeros(datasize - 1)
        mean_strength_plus2sd = copy(mean_strength)
        mean_strength_min2sd = copy(mean_strength)
        mean_strength_true = copy(mean_strength)

        p = (0.08, 1.7, pref[model], 0.15, 2.3, 0.5, 0.25)

        Xlib = vec_data[model][1:(datasize-1), :]

        vec_jmat_gpr[model] = zeros(npop, npop, datasize-1)
        vec_jmat_true[model] = zeros(npop, npop, datasize-1)

        for i in 1:datasize-1
            jmat_mean = zeros(npop, npop)
            jmat_sd = copy(jmat_mean)
            jmat_true = copy(jmat_mean)

            xtarget = Xlib[i, :]
            # rescle data to the original scale
            xtarget_unscaled = xtarget.*sd[model] .+ mu[model]

            # perturb each state with Gaussian noise with mu = 0, sd = 0.1
            perturb = rand(Normal(0, 0.1), npop)
            xtarget_unscaled = xtarget_unscaled + xtarget_unscaled .* perturb

            #standardise pertrubed state
            xtarget = (xtarget_unscaled - mu[model]) ./ sd[model]

            for target in 1:npop

                Ylib = vec_data[model][2:end, target]
                vec_res_gpr_target = vec_res[model].vec_res_gpr[target]

                pars = vec_res_gpr_target[1]
                invL = vec_res_gpr_target[3]
                dmax = vec_res_gpr_target[4]


                jmat_mean[target, :] = deriv_mean(pars, Xlib, xtarget, Ylib, invL, dmax)
                jmat_sd[target, :] = deriv_sd(pars, Xlib, xtarget, invL, dmax)
            end

            jmat_mean = jmat_mean - diagm(ones(npop))

            # true jacobian matrix
            jmat_true = jmat_temp(xtarget_unscaled, p)
            # scale jacobian elements
            jmat_true_scaled = (jmat_true.*sd[model]')./sd[model]

            # mean value of mean interaction strength and +/- 2*sd
            mean_strength[i] = mean(jmat_mean)
            mean_strength_plus2sd[i] = mean(jmat_mean + 2*jmat_sd)
            mean_strength_min2sd[i] = mean(jmat_mean - 2*jmat_sd)

            # true mean strength
            mean_strength_true[i] = mean(jmat_true_scaled)

            vec_jmat_gpr[model][:,:,i] .= jmat_mean
            vec_jmat_true[model][:,:,i] .= jmat_true_scaled
        end


        vec_mean_strength[model] = mean_strength
        vec_mean_strength_plus2sd[model] = mean_strength_plus2sd
        vec_mean_strength_min2sd[model] = mean_strength_min2sd
        vec_mean_strength_true[model] = mean_strength_true
    end

    return (vec_mean_strength=vec_mean_strength,
    vec_mean_strength_plus2sd=vec_mean_strength_plus2sd,
    vec_mean_strength_min2sd=vec_mean_strength_min2sd,
    vec_mean_strength_true=vec_mean_strength_true,
    vec_jmat_gpr,
    vec_jmat_true)
end

function sim_mean_strength(vec_data, vec_res; mu, sd)
    nmodel = length(vec_data)
    datasize = size(vec_data[1], 1)
    npop = size(vec_data[1], 2)

    vec_mean_strength = Vector{Vector{Float64}}(undef, nmodel)
    vec_mean_strength_plus2sd = copy(vec_mean_strength)
    vec_mean_strength_min2sd = copy(vec_mean_strength)
    vec_mean_strength_true = copy(vec_mean_strength)

    vec_jmat_gpr = Vector{Array{Float64}}(undef, nmodel)
    vec_jmat_true = copy(vec_jmat_gpr)

    # food web model parameters
    p = (0.1, 0.07, 3.2, 2.9, 0.5, 0.5, 0.15, 0.15, 2.5, 2.0, 0.3, 1.2)

    for model in 1:nmodel

        mean_strength = zeros(datasize - 1)
        mean_strength_plus2sd = copy(mean_strength)
        mean_strength_min2sd = copy(mean_strength)
        mean_strength_true = copy(mean_strength)

        Xlib = vec_data[model][1:(datasize-1), :]

        vec_jmat_gpr[model] = zeros(npop, npop, datasize-1)
        vec_jmat_true[model] = zeros(npop, npop, datasize-1)

        for i in 1:datasize-1
            jmat_mean = zeros(npop, npop)
            jmat_sd = copy(jmat_mean)
            jmat_true = copy(jmat_mean)

            xtarget = Xlib[i, :]
            # rescle data to the original scale
            xtarget_unscaled = xtarget.*sd[model] .+ mu[model]

            # perturb each state x(t) with Gaussian noise with mu = 0, sd = 0.1*x(t)
            perturb = rand(Normal(0, 0.1), npop)
            xtarget_unscaled = xtarget_unscaled + xtarget_unscaled .* perturb

            #standardise pertrubed state
            xtarget = (xtarget_unscaled - mu[model]) ./ sd[model]

            for target in 1:npop

                Ylib = vec_data[model][2:end, target]
                vec_res_gpr_target = vec_res[model].vec_res_gpr[target]

                pars = vec_res_gpr_target[1]
                invL = vec_res_gpr_target[3]
                dmax = vec_res_gpr_target[4]


                jmat_mean[target, :] = deriv_mean(pars, Xlib, xtarget, Ylib, invL, dmax)
                jmat_sd[target, :] = deriv_sd(pars, Xlib, xtarget, invL, dmax)
            end

            jmat_mean = jmat_mean - diagm(ones(npop))

            # true jacobian matrix
            jmat_true = jmat_fw5_temp(xtarget_unscaled, p)
            # scale jacobian elements
            jmat_true_scaled = (jmat_true.*sd[model]')./sd[model]

            # mean value of mean interaction strength and +/- 2*sd
            mean_strength[i] = mean(jmat_mean)
            mean_strength_plus2sd[i] = mean(jmat_mean + 2*jmat_sd)
            mean_strength_min2sd[i] = mean(jmat_mean - 2*jmat_sd)

            # true mean strength
            mean_strength_true[i] = mean(jmat_true_scaled)

            vec_jmat_gpr[model][:,:,i] .= jmat_mean
            vec_jmat_true[model][:,:,i] .= jmat_true_scaled
        end


        vec_mean_strength[model] = mean_strength
        vec_mean_strength_plus2sd[model] = mean_strength_plus2sd
        vec_mean_strength_min2sd[model] = mean_strength_min2sd
        vec_mean_strength_true[model] = mean_strength_true
    end

    return (vec_mean_strength=vec_mean_strength,
    vec_mean_strength_plus2sd=vec_mean_strength_plus2sd,
    vec_mean_strength_min2sd=vec_mean_strength_min2sd,
    vec_mean_strength_true=vec_mean_strength_true,
    vec_jmat_gpr,
    vec_jmat_true)
end

function sim_mean_strength_logistic_10(vec_data, vec_res, vec_growth, intmat; mu, sd)
    nmodel = length(vec_data)
    datasize = size(vec_data[1], 1)
    npop = size(vec_data[1], 2)

    vec_mean_strength = Vector{Vector{Float64}}(undef, nmodel)
    vec_mean_strength_plus2sd = copy(vec_mean_strength)
    vec_mean_strength_min2sd = copy(vec_mean_strength)
    vec_mean_strength_true = copy(vec_mean_strength)

    vec_jmat_gpr = Vector{Array{Float64}}(undef, nmodel)
    vec_jmat_true = copy(vec_jmat_gpr)

    for model in 1:nmodel

        mean_strength = zeros(datasize - 1)
        mean_strength_plus2sd = copy(mean_strength)
        mean_strength_min2sd = copy(mean_strength)
        mean_strength_true = copy(mean_strength)

        Xlib = vec_data[model][1:(datasize-1), :]

        vec_jmat_gpr[model] = zeros(npop, npop, datasize-1)
        vec_jmat_true[model] = zeros(npop, npop, datasize-1)

        vec_growth_temp = vec_growth[model]
        intmat_temp = intmat[model]

        for i in 1:datasize-1
            jmat_mean = zeros(npop, npop)
            jmat_sd = copy(jmat_mean)
            jmat_true = copy(jmat_mean)

            xtarget = Xlib[i, :]
            # rescle data to the original scale
            xtarget_unscaled = xtarget.*sd[model] .+ mu[model]

            # perturb each state x(t) with Gaussian noise with mu = 0, sd = 0.1*x(t)
            perturb = rand(Normal(0, 0.1), npop)
            xtarget_unscaled = xtarget_unscaled + xtarget_unscaled .* perturb

            #standardise pertrubed state
            xtarget = (xtarget_unscaled - mu[model]) ./ sd[model]

            for target in 1:npop

                Ylib = vec_data[model][2:end, target]
                vec_res_gpr_target = vec_res[model].vec_res_gpr[target]

                pars = vec_res_gpr_target[1]
                invL = vec_res_gpr_target[3]
                dmax = vec_res_gpr_target[4]


                jmat_mean[target, :] = deriv_mean(pars, Xlib, xtarget, Ylib, invL, dmax)
                jmat_sd[target, :] = deriv_sd(pars, Xlib, xtarget, invL, dmax)
            end

            # true jacobian matrix
            jmat_true = jacobian_logistic(xtarget_unscaled, vec_growth_temp, intmat_temp)
            # scale jacobian elements
            jmat_true_scaled = (jmat_true.*sd[model]')./sd[model]

            # mean value of mean interaction strength and +/- 2*sd
            mean_strength[i] = mean(jmat_mean)
            mean_strength_plus2sd[i] = mean(jmat_mean + 2*jmat_sd)
            mean_strength_min2sd[i] = mean(jmat_mean - 2*jmat_sd)

            # true mean strength
            mean_strength_true[i] = mean(jmat_true_scaled)

            vec_jmat_gpr[model][:,:,i] .= jmat_mean
            vec_jmat_true[model][:,:,i] .= jmat_true_scaled
        end


        vec_mean_strength[model] = mean_strength
        vec_mean_strength_plus2sd[model] = mean_strength_plus2sd
        vec_mean_strength_min2sd[model] = mean_strength_min2sd
        vec_mean_strength_true[model] = mean_strength_true
    end

    return (vec_mean_strength=vec_mean_strength,
    vec_mean_strength_plus2sd=vec_mean_strength_plus2sd,
    vec_mean_strength_min2sd=vec_mean_strength_min2sd,
    vec_mean_strength_true=vec_mean_strength_true,
    vec_jmat_gpr,
    vec_jmat_true)
end

# function for simulate state-dependent interactions with data containing observational noise
function sim_mean_strength_obs(vec_data, # data without observational noise to calculate true Jacobian
     vec_data_obs, # training data
     vec_res; # gpr results
     mu_obs, # mean of training data
     sd_obs # sd of training data
     )
    nmodel = length(vec_data)
    datasize = size(vec_data[1], 1)
    npop = size(vec_data[1], 2)

    vec_mean_strength = Vector{Vector{Float64}}(undef, nmodel)
    vec_mean_strength_plus2sd = copy(vec_mean_strength)
    vec_mean_strength_min2sd = copy(vec_mean_strength)
    vec_mean_strength_true = copy(vec_mean_strength)

    vec_jmat_gpr = Vector{Array{Float64}}(undef, nmodel)
    vec_jmat_true = copy(vec_jmat_gpr)

    # food web model parameters
    p = (0.1, 0.07, 3.2, 2.9, 0.5, 0.5, 0.15, 0.15, 2.5, 2.0, 0.3, 1.2)

    for model in 1:nmodel

        mean_strength = zeros(datasize - 1)
        mean_strength_plus2sd = copy(mean_strength)
        mean_strength_min2sd = copy(mean_strength)
        mean_strength_true = copy(mean_strength)

        Xlib_obs = vec_data_obs[model][1:(datasize-1), :]
        Xlib = vec_data[model][1:(datasize-1), :]

        vec_jmat_gpr[model] = zeros(npop, npop, datasize-1)
        vec_jmat_true[model] = zeros(npop, npop, datasize-1)

        for i in 1:datasize-1
            jmat_mean = zeros(npop, npop)
            jmat_sd = copy(jmat_mean)
            jmat_true = copy(jmat_mean)

            xtarget_obs = Xlib_obs[i, :]
            xtarget = Xlib[i, :]
            
            # rescle data for GPR model to the original scale
            xtarget_obs_unscaled = xtarget_obs.*sd_obs[model] .+ mu_obs[model]

            # perturb each state x(t) with Gaussian noise with mu = 0, sd = 0.1*x(t)
            perturb = rand(Normal(0, 0.1), npop)
            xtarget_obs_unscaled = xtarget_obs_unscaled + xtarget_obs_unscaled .* perturb
            # true state with perturbations
            xtarget = xtarget + xtarget .* perturb

            #standardise pertrubed state
            xtarget_obs_scaled = (xtarget_obs_unscaled - mu_obs[model]) ./ sd_obs[model]

            for target in 1:npop

                Ylib = vec_data_obs[model][2:end, target]
                vec_res_gpr_target = vec_res[model].vec_res_gpr[target]

                pars = vec_res_gpr_target[1]
                invL = vec_res_gpr_target[3]
                dmax = vec_res_gpr_target[4]


                jmat_mean[target, :] = deriv_mean(pars, Xlib_obs, xtarget_obs_scaled, Ylib, invL, dmax)
                jmat_sd[target, :] = deriv_sd(pars, Xlib_obs, xtarget_obs_scaled, invL, dmax)
            end

            jmat_mean = jmat_mean - diagm(ones(npop))

            # true jacobian matrix
            jmat_true = jmat_fw5_temp(xtarget, p)
            # scale jacobian elements
            jmat_true_scaled = (jmat_true.*sd_obs[model]')./sd_obs[model]

            # mean value of mean interaction strength and +/- 2*sd
            mean_strength[i] = mean(jmat_mean)
            mean_strength_plus2sd[i] = mean(jmat_mean + 2*jmat_sd)
            mean_strength_min2sd[i] = mean(jmat_mean - 2*jmat_sd)

            # true mean strength
            mean_strength_true[i] = mean(jmat_true_scaled)

            vec_jmat_gpr[model][:,:,i] .= jmat_mean
            vec_jmat_true[model][:,:,i] .= jmat_true_scaled
        end


        vec_mean_strength[model] = mean_strength
        vec_mean_strength_plus2sd[model] = mean_strength_plus2sd
        vec_mean_strength_min2sd[model] = mean_strength_min2sd
        vec_mean_strength_true[model] = mean_strength_true
    end

    return (vec_mean_strength=vec_mean_strength,
    vec_mean_strength_plus2sd=vec_mean_strength_plus2sd,
    vec_mean_strength_min2sd=vec_mean_strength_min2sd,
    vec_mean_strength_true=vec_mean_strength_true,
    vec_jmat_gpr=vec_jmat_gpr,
    vec_jmat_true=vec_jmat_true)
end

function sim_mean_strength_switching_obs(vec_data, vec_data_obs, vec_res, pref; mu_obs, sd_obs)
    nmodel = length(vec_data)
    datasize = size(vec_data[1], 1)
    npop = size(vec_data[1], 2)

    vec_mean_strength = Vector{Vector{Float64}}(undef, nmodel)
    vec_mean_strength_plus2sd = copy(vec_mean_strength)
    vec_mean_strength_min2sd = copy(vec_mean_strength)
    vec_mean_strength_true = copy(vec_mean_strength)

    vec_jmat_gpr = Vector{Array{Float64}}(undef, nmodel)
    vec_jmat_true = copy(vec_jmat_gpr)

    for model in 1:nmodel

        mean_strength = zeros(datasize - 1)
        mean_strength_plus2sd = copy(mean_strength)
        mean_strength_min2sd = copy(mean_strength)
        mean_strength_true = copy(mean_strength)

        p = (0.08, 1.7, pref[model], 0.15, 2.3, 0.5, 0.25)

        Xlib_obs = vec_data_obs[model][1:(datasize-1), :]
        Xlib = vec_data[model][1:(datasize-1), :]

        vec_jmat_gpr[model] = zeros(npop, npop, datasize-1)
        vec_jmat_true[model] = zeros(npop, npop, datasize-1)

        for i in 1:datasize-1
            jmat_mean = zeros(npop, npop)
            jmat_sd = copy(jmat_mean)
            jmat_true = copy(jmat_mean)

            xtarget_obs = Xlib_obs[i, :]
            xtarget = Xlib[i, :]

            # rescle training data to the original scale
            xtarget_obs_unscaled = xtarget_obs.*sd_obs[model] .+ mu_obs[model]

            # perturb each state with Gaussian noise with mu = 0, sd = 0.1
            perturb = rand(Normal(0, 0.1), npop)
            xtarget_unscaled_perturbed = xtarget_obs_unscaled + xtarget_obs_unscaled .* perturb
            xtarget_perturbed = xtarget + xtarget .* perturb

            #standardise pertrubed state
            xtarget_scaled_perturbed = (xtarget_unscaled_perturbed - mu_obs[model]) ./ sd_obs[model]

            for target in 1:npop

                Ylib = vec_data_obs[model][2:end, target]
                vec_res_gpr_target = vec_res[model].vec_res_gpr[target]

                pars = vec_res_gpr_target[1]
                invL = vec_res_gpr_target[3]
                dmax = vec_res_gpr_target[4]


                jmat_mean[target, :] = deriv_mean(pars, Xlib_obs, xtarget_scaled_perturbed, Ylib, invL, dmax)
                jmat_sd[target, :] = deriv_sd(pars, Xlib_obs, xtarget_scaled_perturbed, invL, dmax)
            end

            jmat_mean = jmat_mean - diagm(ones(npop))

            # true jacobian matrix
            jmat_true = jmat_temp(xtarget_perturbed, p)
            # scale jacobian elements
            jmat_true_scaled = (jmat_true.*sd_obs[model]')./sd_obs[model]

            # mean value of mean interaction strength and +/- 2*sd
            mean_strength[i] = mean(jmat_mean)
            mean_strength_plus2sd[i] = mean(jmat_mean + 2*jmat_sd)
            mean_strength_min2sd[i] = mean(jmat_mean - 2*jmat_sd)

            # true mean strength
            mean_strength_true[i] = mean(jmat_true_scaled)

            vec_jmat_gpr[model][:,:,i] .= jmat_mean
            vec_jmat_true[model][:,:,i] .= jmat_true_scaled
        end


        vec_mean_strength[model] = mean_strength
        vec_mean_strength_plus2sd[model] = mean_strength_plus2sd
        vec_mean_strength_min2sd[model] = mean_strength_min2sd
        vec_mean_strength_true[model] = mean_strength_true
    end

    return (vec_mean_strength=vec_mean_strength,
    vec_mean_strength_plus2sd=vec_mean_strength_plus2sd,
    vec_mean_strength_min2sd=vec_mean_strength_min2sd,
    vec_mean_strength_true=vec_mean_strength_true,
    vec_jmat_gpr,
    vec_jmat_true)
end

function sim_mean_strength_logistic_10_obs(vec_data, vec_data_obs, vec_res, vec_growth, intmat; mu_obs, sd_obs)
    nmodel = length(vec_data)
    datasize = size(vec_data[1], 1)
    npop = size(vec_data[1], 2)

    vec_mean_strength = Vector{Vector{Float64}}(undef, nmodel)
    vec_mean_strength_plus2sd = copy(vec_mean_strength)
    vec_mean_strength_min2sd = copy(vec_mean_strength)
    vec_mean_strength_true = copy(vec_mean_strength)

    vec_jmat_gpr = Vector{Array{Float64}}(undef, nmodel)
    vec_jmat_true = copy(vec_jmat_gpr)

    for model in 1:nmodel

        mean_strength = zeros(datasize - 1)
        mean_strength_plus2sd = copy(mean_strength)
        mean_strength_min2sd = copy(mean_strength)
        mean_strength_true = copy(mean_strength)

        Xlib_obs = vec_data_obs[model][1:(datasize-1), :]
        Xlib = vec_data[model][1:(datasize-1), :]

        vec_jmat_gpr[model] = zeros(npop, npop, datasize-1)
        vec_jmat_true[model] = zeros(npop, npop, datasize-1)

        vec_growth_temp = vec_growth[model]
        intmat_temp = intmat[model]

        for i in 1:datasize-1
            jmat_mean = zeros(npop, npop)
            jmat_sd = copy(jmat_mean)
            jmat_true = copy(jmat_mean)

            xtarget_obs = Xlib_obs[i, :]
            xtarget = Xlib[i, :]

            # rescle training data to the original scale
            xtarget_obs_unscaled = xtarget_obs.*sd_obs[model] .+ mu_obs[model]

            # perturb each state with Gaussian noise with mu = 0, sd = 0.1
            perturb = rand(Normal(0, 0.1), npop)
            xtarget_unscaled_perturbed = xtarget_obs_unscaled + xtarget_obs_unscaled .* perturb
            xtarget_perturbed = xtarget + xtarget .* perturb

            #standardise pertrubed state
            xtarget_scaled_perturbed = (xtarget_unscaled_perturbed - mu_obs[model]) ./ sd_obs[model]

            for target in 1:npop

                Ylib = vec_data_obs[model][2:end, target]
                vec_res_gpr_target = vec_res[model].vec_res_gpr[target]

                pars = vec_res_gpr_target[1]
                invL = vec_res_gpr_target[3]
                dmax = vec_res_gpr_target[4]


                jmat_mean[target, :] = deriv_mean(pars, Xlib_obs, xtarget_scaled_perturbed, Ylib, invL, dmax)
                jmat_sd[target, :] = deriv_sd(pars, Xlib_obs, xtarget_scaled_perturbed, invL, dmax)
            end

            # true jacobian matrix
            jmat_true = jacobian_logistic(xtarget_perturbed, vec_growth_temp, intmat_temp)
            # scale jacobian elements
            jmat_true_scaled = (jmat_true.*sd_obs[model]')./sd_obs[model]

            # mean value of mean interaction strength and +/- 2*sd
            mean_strength[i] = mean(jmat_mean)
            mean_strength_plus2sd[i] = mean(jmat_mean + 2*jmat_sd)
            mean_strength_min2sd[i] = mean(jmat_mean - 2*jmat_sd)

            # true mean strength
            mean_strength_true[i] = mean(jmat_true_scaled)

            vec_jmat_gpr[model][:,:,i] .= jmat_mean
            vec_jmat_true[model][:,:,i] .= jmat_true_scaled
        end


        vec_mean_strength[model] = mean_strength
        vec_mean_strength_plus2sd[model] = mean_strength_plus2sd
        vec_mean_strength_min2sd[model] = mean_strength_min2sd
        vec_mean_strength_true[model] = mean_strength_true
    end

    return (vec_mean_strength=vec_mean_strength,
    vec_mean_strength_plus2sd=vec_mean_strength_plus2sd,
    vec_mean_strength_min2sd=vec_mean_strength_min2sd,
    vec_mean_strength_true=vec_mean_strength_true,
    vec_jmat_gpr,
    vec_jmat_true)
end

# calculating the simulated change of mean strength and stats of the inferences
function stats_sim_jmat(res_mean_strength, vec_res, vec_jmat_scaled; diag_correction=false)
    nmodel = length(vec_res)
    datasize = size(vec_jmat_scaled[1], 3)
    npop = length(vec_res[1].vec_res_gpr)

    vec_mean_strength_gpr = res_mean_strength.vec_mean_strength
    vec_mean_strength_true = res_mean_strength.vec_mean_strength_true
    vec_jmat_gpr = res_mean_strength.vec_jmat_gpr
    vec_jmat_true = res_mean_strength.vec_jmat_true

    vec_cor = zeros(nmodel)
    vec_cor_jmat = copy(vec_cor)
    vec_mean_strength_before_true = Vector{Vector{Float64}}(undef, nmodel)
    vec_mean_strength_before_gpr = copy(vec_mean_strength_before_true)

    vec_diff_mean_strength_gpr = copy(vec_mean_strength_before_gpr)
    vec_diff_mean_strength_true = copy(vec_mean_strength_before_gpr)

    vec_sign_diff = copy(vec_cor)

    for model in 1:nmodel
        vec_cor[model] = cor(vec_mean_strength_gpr[model], vec_mean_strength_true[model])
        vec_cor_jmat[model] = cor(vec(vec_jmat_gpr[model]), vec(vec_jmat_true[model]))

        mean_strength_before_gpr = zeros(datasize-1)
        mean_strength_before_true = copy(mean_strength_before_gpr)

        for iter in 1:datasize-1
            mean_strength_before_gpr[iter] = mean(vec_res[model].jmat[:,:,iter] .- ifelse(diag_correction, diagm(ones(npop)), 0))
            mean_strength_before_true[iter] = mean(vec_jmat_scaled[model][:,:,iter])
        end

        vec_mean_strength_before_gpr[model] = mean_strength_before_gpr
        vec_mean_strength_before_true[model] = mean_strength_before_true

        vec_diff_mean_strength_gpr[model] = vec_mean_strength_gpr[model] - mean_strength_before_gpr
        vec_diff_mean_strength_true[model] = vec_mean_strength_true[model] - mean_strength_before_true
        vec_sign_diff[model] = sum(sign.(vec_diff_mean_strength_gpr[model]) .== sign.(vec_diff_mean_strength_true[model]))/(datasize-1)
    end

    return (vec_diff_mean_strength_gpr=vec_diff_mean_strength_gpr, vec_diff_mean_strength_true=vec_diff_mean_strength_true, vec_cor_jmat=vec_cor_jmat, vec_sign_diff=vec_sign_diff)
end

function scatter_jmat_sim(vec_jmat_gpr, vec_jmat_scaled; model=1)
	datasize = size(vec_jmat_gpr[1], 3)
    nsp = size(vec_jmat_gpr[1], 1)
	
    res_jmat = copy(vec_jmat_gpr[model])
	
    vec_coeffs = vec(res_jmat[:,:,1])
    vec_coeffs_true = vec(vec_jmat_scaled[model][:,:,1])
    for t in 2:datasize
        vec_coeffs = vcat(vec_coeffs, vec(res_jmat[:,:,t]))
        vec_coeffs_true = vcat(vec_coeffs_true, vec(vec_jmat_scaled[model][:,:,t]))
	end

    coeff_min = minimum(vec_coeffs_true)
    coeff_max = maximum(vec_coeffs_true)
    rho = cor(vec_coeffs, vec_coeffs_true)

    bhat = [vec_coeffs_true ones((datasize)nsp^2)]\vec_coeffs
    scatter(vec_coeffs_true, vec_coeffs, markersize=2, msw=0, xlab="theoretical", ylab="estimated", label="ρ = $(round(rho, digits=2))", legendfontsize=19, margin=0.3Plots.cm)
    #rho = cor(vec(jmat_scaled_fw1080[nozero,nozero,timepoint]), vec(jmat_RDE_fw10[:,:,timepoint]))
    plot!(coeff_min:0.01:coeff_max, x -> x, c=:red, label=false)
    Plots.abline!(bhat..., label=false)
end
