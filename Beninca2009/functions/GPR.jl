### Functions for GPR inference

include("Rprop.jl")

using JLD2
import StatsBase.mean
import StatsBase.std
import DataFrames.Not

# load saved GPR results from .jld2 file
function load_gpresults_jld2(filepath; npop, jacobian=true)
    res_gpr = jldopen(filepath, "r") do file
        vec_res_gpr = Vector{NamedTuple}(undef, npop)
        loocv = Vector{NamedTuple}(undef, npop)
        for i in 1:npop
            vec_res_gpr[i] = file["gpr$i/res_gpr_target$i"]
            loocv[i] = file["loocv$i/loocv_target$i"]
        end

        if jacobian
            jmat = file["jmat/res_jmat"]
            jmat_sd = file["jmat/res_jmat_sd"]
        end
        
        (vec_res_gpr = vec_res_gpr, jmat = jmat, jmat_sd = jmat_sd, loocv = loocv)
    end
    return res_gpr
end

# LOOCV to check the fitting with in-sample errors
function error_loocv(pars, Ytr, Xtr, Ktr, dmax)
    datasize = length(Ytr)
    sqerror_loo = 0.0
    mu_loo = zeros(datasize)
    var_loo = similar(mu_loo)

    for i in 1:datasize
        noti = Not(i)

        cholK_woi = cholesky(Symmetric(Ktr[noti, noti]))
        invL_woi = inv(cholK_woi.L)
        invL_Ytr = invL_woi * Ytr[noti]

        ks = cov_pred_onestep(pars, Xtr[noti, :], Xtr[i, :], dmax)
        invL_ks = invL_woi * ks
        ktes = pars[end-1] + pars[end] #eta + sigma
        mu_loo[i] = invL_ks' * invL_Ytr
        var_loo[i] = ktes - invL_ks' * invL_ks
        sqerror_loo += (Ytr[i] - mu_loo[i])^2 + var_loo[i]
    end

    rmse_loo = sqrt(sqerror_loo / datasize)

    return (rmse_loo=rmse_loo, mu_loo=mu_loo, sd_loo=sqrt.(var_loo))
end


# Posterior mean of partial derivatives (Jacobian elements)
function deriv_mean(pars, Xtr, xtes, Ytr, invL, dmax, scale_dist=true)
    E = size(Xtr, 2) #embedding dimension
    theta = pars[1:E]
    eta = pars[end-1]

    deriv_mean = zeros(E)

    cholalpha = invL' * invL * Ytr

    ks = cov_pred_onestep(pars, Xtr, xtes, dmax)

    for i in 1:E
        dks = -2 * theta[i] * ((xtes[i] .- Xtr[:, i]) / dmax^2) .* ks
        deriv_mean[i] = dks' * cholalpha
    end

    return deriv_mean
end


#(2nd-order partial derivatives of se kernel)
function kernel_deriv2(pars, xtarget, x2, theta_i, dmax)
    E = length(xtarget)
    theta = pars[1:E]
    eta = pars[end-1]
    d2k = 2 * (theta[theta_i] / dmax^2) * (1 - 2 * (theta[theta_i] / dmax^2) * (xtarget[theta_i] - x2[theta_i])^2) * sekernel(eta, theta, xtarget, x2, dmax)
    return d2k
end

#(1st-order partial derivatives of se kernel)
function kernel_deriv(pars, xtarget, x2, theta_i, dmax)
    E = length(xtarget)
    theta = pars[1:E]
    eta = pars[end-1]
    dk = -2 * theta[theta_i] * ((xtarget[theta_i] - x2[theta_i]) / dmax^2) * sekernel(eta, theta, xtarget, x2, dmax)
    return dk
end


# Posterior sd of partial derivatives
function deriv_sd(pars, Xtr, xtarget, invL, dmax)
    E = length(xtarget) #embedding dimension

    deriv_sd = zeros(E)

    for i in 1:E
        ktarget = kernel_deriv2(pars, xtarget, xtarget, i, dmax)
        kcross(X2) = kernel_deriv(pars, xtarget, X2, i, dmax)
        ktarget_in = kcross.(eachrow(Xtr))

        deriv_sd[i] = sqrt(ktarget - (invL * ktarget_in)' * invL * ktarget_in)
    end

    return deriv_sd
end


# Main function for fitting GPR model to data
function fit_gpr_synthetic(data_tr;
    filename = "./result_gpr.jld2",
    save_file=false,
    jacobian=true,
    scale_dist=true
)

    npop = size(data_tr, 2)
    size_tr = size(data_tr, 1)

    vec_res_gpr = Vector{NamedTuple}(undef, npop)
    jmat_tr = Array{Float64}(undef, npop, npop, size_tr)
    jmat_sd = Array{Float64}(undef, npop, npop, size_tr)
    loocv = Vector{NamedTuple}(undef, npop)

    if jacobian
        for target in 1:npop
            Xtr = data_tr[1:(size_tr-1), :]
            Ytr = data_tr[2:end, target]
            vec_res_gpr[target] = Rprop(Xtr, Ytr; maxiter = 100, scale_dist = scale_dist)
            pars = vec_res_gpr[target][1]
            invL = vec_res_gpr[target][3]
            dmax = vec_res_gpr[target][4]
            Ktr = vec_res_gpr[target][6]
            loocv[target] = error_loocv(pars, Ytr, Xtr, Ktr, dmax)

            for i in 1:size_tr
                jmat_tr[target, :, i] = deriv_mean(pars, Xtr, data_tr[i , :], Ytr, invL, dmax)
                jmat_sd[target, :, i] = deriv_sd(pars, Xtr, data_tr[i , :], invL, dmax)
            end
            println("variable ", target, " finished")
        end
    else
        for target in 1:npop
            Xtr = data_tr[1:(size_tr-1), :]
            Ytr = data_tr[2:end, target]
            vec_res_gpr[target] = Rprop(Xtr, Ytr; maxiter = 100, scale_dist = scale_dist)
            pars = vec_res_gpr[target][1]
            invL = vec_res_gpr[target][3]
            dmax = vec_res_gpr[target][4]
            Ktr = vec_res_gpr[target][6]
            loocv[target] = error_loocv(pars, Ytr, Xtr, Ktr, dmax)
        end
    end

    if save_file
        jldopen("$filename", "w") do file
            for i in 1:npop
                file["gpr$i/res_gpr_target$i"] = vec_res_gpr[i]
                file["loocv$i/loocv_target$i"] = loocv[i]
            end
            if jacobian
                file["jmat/res_jmat"] = jmat_tr
                file["jmat/res_jmat_sd"] = jmat_sd
            end
        end
    end

    if jacobian
        return (vec_res_gpr = vec_res_gpr, jmat = jmat_tr, jmat_sd = jmat_sd, loocv = loocv)
    else
        return (vec_res_gpr = vec_res_gpr, loocv = loocv)
    end
end

