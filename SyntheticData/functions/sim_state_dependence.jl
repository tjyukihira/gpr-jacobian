### function to simulate state-dependence of interaction strengths
function sim_state_dependence(res_gpr, data, mu, sd; target_inc=1, target_dec=false, kernel="se")
    Xtr = data[1:end-1, :]

    nsp = size(Xtr, 2)
    datasize = size(Xtr, 1)

    jmat_sim = zeros(nsp, nsp, datasize)

    #data with 10% more abundance added to target species
    input_sim = copy(Xtr)
    if target_inc != false
        input_sim[:, target_inc] .= 1.2 * Xtr[:, target_inc] .+ 0.2 * mu[target_inc] / sd[target_inc]
    end
    if target_dec != false
        input_sim[:, target_dec] .= 0.8 * Xtr[:, target_dec] .- 0.2 * mu[target_dec] / sd[target_dec]
    end

    for i_sp in 1:nsp
        res_target = res_gpr[1][i_sp]

        pars = res_target[1]
        invL = res_target[3]
        dmax = res_target[4]

        Ytr = data[2:end, i_sp]

        for i in 1:datasize
            jmat_sim[i_sp, :, i] = deriv_mean(pars, Xtr, input_sim[i, :], Ytr, invL, dmax, kernel, true)
        end
    end

    return (jmat_sim=jmat_sim)
end