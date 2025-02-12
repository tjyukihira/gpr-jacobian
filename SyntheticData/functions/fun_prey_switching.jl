### simulation of 5-species food web model
## u[1]: P, u[2]: C1, u[3]: C2, u[4]: R1, u[5]: R2
function food_web5_switching!(du, u, p, t)
    xp, yp, pref, xc, yc, C0, R0 = p
    delta1 = pref * u[2] / (pref * u[2] + (1 - pref) * u[3])
    delta2 = (1 - pref) * u[3] / (pref * u[2] + (1 - pref) * u[3])
    du[1] = xp * u[1] * (delta1 * yp * u[2] / (u[2] + C0) + delta2 * yp * u[3] / (u[3] + C0) - 1)
    du[2] = u[2] * (xc * yc * u[4] / (u[4] + R0) - delta1 * xp * yp * u[1] / (u[2] + C0) - xc)
    du[3] = u[3] * (xc * yc * u[5] / (u[5] + R0) - delta2 * xp * yp * u[1] / (u[3] + C0) - xc)
    du[4] = u[4] * (1 - u[4] - xc * yc * u[2] / (u[4] + R0))
    du[5] = u[5] * (1 - u[5] - xc * yc * u[3] / (u[5] + R0))
end

### process noise
function g_food_web_switching(du, u, p, t)
    du[1] = u[1] * 0.05
    du[2] = u[2] * 0.05
    du[3] = u[3] * 0.05
    du[4] = u[4] * 0.05
    du[5] = u[5] * 0.05
end

### intense process noise
function g_food_web_switching_intense(du, u, p, t)
    du[1] = u[1] * 0.1
    du[2] = u[2] * 0.1
    du[3] = u[3] * 0.1
    du[4] = u[4] * 0.1
    du[5] = u[5] * 0.1
end

function sim_prey_swtching_noise(; nmodel=50, pref, tspan, intense=false)
    data_switching_noise = Vector{Array{Float64}}(undef, nmodel)
    iter = 1
    while iter <= nmodel
        p = (0.08, 1.7, pref[iter], 0.15, 2.3, 0.5, 0.25)
        #initial values
        u0 = vcat(rand(Uniform(0.8, 1.0)), rand(Uniform(0.2, 0.5), 4))

        g = ifelse(intense, g_food_web_switching_intense, g_food_web_switching)

        probsde = SDEProblem(food_web5_switching!, g, u0, tspan, p)

        solsde = solve(probsde, dt=0.05, adaptive=false)
        #plot(solode[10002:end], idxs=(1,2,4), xlab="P", ylab="C1", zlab="R1")
        #plot(solode[10002:end], idxs=(4, 3, 5))

        if mean(solsde[10002:20:end]'[1:end, 1]) > 0.5
            data_switching_noise[iter] = solsde[10002:20:end]'
            iter += 1
        end
    end

    return (data_switching_noise)
end

#function to calculate theoretical jacobian matrix
function jmat_switching(data, p)
    jmat = zeros(5, 5, size(data, 1))
    for t in 1:size(data, 1)
        u = data[t, :]
        jmat[:,:,t] .= jmat_temp(u, p)
    end

    return jmat
end

function jmat_temp(u, p)
    jmat = zeros(5, 5)

    xp, yp, pref, xc, yc, C0, R0 = p

    delta1 = pref * u[2] / (pref * u[2] + (1 - pref) * u[3])
    delta2 = (1 - pref) * u[3] / (pref * u[2] + (1 - pref) * u[3])
    ddelta1_c1 = (1 - pref) * pref * u[3] / (pref * u[2] + (1 - pref) * u[3])^2
    ddelta2_c2 = (1 - pref) * pref * u[2] / (pref * u[2] + (1 - pref) * u[3])^2
    ddelta1_c2 = -(1 - pref) * pref * u[2] / (pref * u[2] + (1 - pref) * u[3])^2
    ddelta2_c1 = -(1 - pref) * pref * u[3] / (pref * u[2] + (1 - pref) * u[3])^2

    jmat[1, 1] = delta1 * xp * yp * u[2] / (u[2] + C0) + delta2 * xp * yp * u[3] / (u[3] + C0) - xp
    jmat[1, 2] = xp * yp * u[1] * (delta1 * C0 / (u[2] + C0)^2 + ddelta1_c1 * u[2] / (u[2] + C0) + ddelta2_c1 * u[3] / (u[3] + C0))
    jmat[1, 3] = xp * yp * u[1] * (delta2 * C0 / (u[3] + C0)^2 + ddelta2_c2 * u[3] / (u[3] + C0) + ddelta1_c2 * u[2] / (u[2] + C0))

    jmat[2, 1] = -delta1 * xp * yp * u[2] / (u[2] + C0)
    jmat[2, 2] = xc * yc * u[4] / (u[4] + R0) - delta1 * xp * yp * u[1] * C0 / (u[2] + C0)^2 - ddelta1_c1 * xp * yp * u[2] * u[1] / (u[2] + C0) - xc
    jmat[2, 3] = -ddelta1_c2 * xp * yp * u[2] * u[1] / (u[2] + C0)
    jmat[2, 4] = xc * yc * u[2] * R0 / (u[4] + R0)^2

    jmat[3, 1] = -delta2 * xp * yp * u[3] / (u[3] + C0)
    jmat[3, 2] = -ddelta2_c1 * xp * yp * u[3] * u[1] / (u[3] + C0)
    jmat[3, 3] = xc * yc * u[5] / (u[5] + R0) - delta2 * xp * yp * u[1] * C0 / (u[3] + C0)^2 - ddelta2_c2 * xp * yp * u[3] * u[1] / (u[3] + C0) - xc
    jmat[3, 5] = xc * yc * u[3] * R0 / (u[5] + R0)^2

    jmat[4, 2] = -xc * yc * u[4] / (u[4] + R0)
    jmat[4, 4] = 1 - 2 * u[4] - xc * yc * u[2] * R0 / (u[4] + R0)^2

    jmat[5, 3] = -xc * yc * u[5] / (u[5] + R0)
    jmat[5, 5] = 1 - 2 * u[5] - xc * yc * u[3] * R0 / (u[5] + R0)^2

    return jmat
end