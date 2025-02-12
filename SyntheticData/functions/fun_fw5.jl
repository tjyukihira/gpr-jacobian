### simulation of 5-species food web model

function food_web5!(du, u, p, t)
    nu1, nu2, lambda1, lambda2, p3, p4, mu1, mu2, kappa1, kappa2, c5, K = p
    du[1] = nu1*u[1]*(lambda1*u[3]/(u[3] + p3) - 1)
    du[2] = nu2*u[2]*(lambda2*u[4]/(u[4] + p4) - 1)
    du[3] = u[3]*(mu1*kappa1*u[5]/(u[5] + c5) - nu1*lambda1*u[1]/(u[3] + p3) - mu1)
    du[4] = u[4]*(mu2*kappa2*u[5]/(u[5] + c5) - nu2*lambda2*u[2]/(u[4] + p4) - mu2)
    du[5] = u[5]*(1 - u[5]/K - mu1*kappa1*u[3]/(u[5] + c5) - mu2*kappa2*u[4]/(u[5] + c5))
end

### process noise
function g_food_web(du, u, p, t)
    du[1] = u[1]*0.05
    du[2] = u[2]*0.05
    du[3] = u[3]*0.05
    du[4] = u[4]*0.05
    du[5] = u[5]*0.05
end

### intense process noise
function g_food_web_intense(du, u, p, t)
    du[1] = u[1]*0.1
    du[2] = u[2]*0.1
    du[3] = u[3]*0.1
    du[4] = u[4]*0.1
    du[5] = u[5]*0.1
end

function sim_fw5_noise(;nmodel=50, p, tspan, intense=false)
    data_fw5_noise = Vector{Array{Float64}}(undef, nmodel)
    iter = 1
    while iter <= nmodel
        #initial values
        u0 = rand(Uniform(0.1,0.5),5)

        g = ifelse(intense, g_food_web_intense, g_food_web)
    
        probsde = SDEProblem(food_web5!, g, u0, tspan, p)
    
        solsde = solve(probsde, dt = 0.05, adaptive = false)
        #plot(solode[10002:end], idxs=(1,2,4), xlab="P", ylab="C1", zlab="R1")
        #plot(solode[10002:end], idxs=(4, 3, 5))
    
        if sum(mean(solsde[10002:20:end], dims = 2)[:, 1] .> 0.1) == 5
            data_fw5_noise[iter] = solsde[10002:20:end]'
            iter += 1
        end
    end

    return(data_fw5_noise)
end

#function to calculate theoretical jacobian matrix
function jmat_fw5(data, p)
    jmat = zeros(5, 5, size(data, 1))
    for t in 1:size(data, 1)
        u = data[t,:]
        jmat[:,:,t] = jmat_fw5_temp(u, p)
    end
    return jmat
end

function jmat_fw5_temp(u, p)
    nu1, nu2, lambda1, lambda2, p3, p4, mu1, mu2, kappa1, kappa2, c5, K = p
    jmat = [nu1*lambda1*u[3]/(u[3]+p3)-nu1 0.0 nu1*lambda1*u[1]*p3/(u[3]+p3)^2 0.0 0.0;
        0.0 nu2*lambda2*u[4]/(u[4]+p4)-nu2 0.0 nu2*lambda2*u[2]*p4/(u[4]+p4)^2 0.0;
        -nu1*lambda1*u[3]/(u[3]+p3) 0.0 mu1*kappa1*u[5]/(u[5]+c5)-nu1*lambda1*u[1]*p3/(u[3]+p3)^2-mu1 0.0 mu1*kappa1*u[3]*c5/(u[5]+c5)^2;
        0.0 -nu2*lambda2*u[4]/(u[4]+p4) 0.0 mu2*kappa2*u[5]/(u[5]+c5)-nu2*lambda2*u[2]*p4/(u[4]+p4)^2-mu2 mu2*kappa2*u[4]*c5/(u[5]+c5)^2;
        0.0 0.0 -mu1*kappa1*u[5]/(u[5]+c5) -mu2*kappa2*u[5]/(u[5]+c5) 1-2*u[5]/K-mu1*kappa1*u[3]*c5/(u[5]+c5)^2-mu2*kappa2*u[4]*c5/(u[5]+c5)^2]
    return jmat
end