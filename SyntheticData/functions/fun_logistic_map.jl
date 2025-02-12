# simulation of logistic map
using Distributions
import LinearAlgebra.diagm

function logisticmap(popvec, vec_growth, intmat; npop=npop, proc_noise=true, sd_noise=0.05)
    #popnew = vec_growth.*popvec.*(1.0 .+ intmat*popvec .+ ifelse(proc_noise, rand(Normal(0, 0.01), npop), zeros(npop)))
    popnew = vec_growth.*popvec.*(1.0 .+ intmat*popvec) .+ ifelse(proc_noise, popvec.*rand(Normal(0, sd_noise), npop), 0.0)
    popnew
end

function jacobian_logistic(popvec, vec_growth, intmat)
    npop = length(popvec)
    jacobian = vec_growth.*popvec.*intmat + vec_growth.*(1.0 .+ intmat*popvec).*diagm(ones(npop))
    jacobian
end

function sim_logistic(npop; initpop, vec_growth, intmat, datasize=100,  iter=2000, proc_noise=true, sd_noise = 0.05)
    #vec_growth = r_growth .+ rand(Uniform(-1.0, 0.5), npop)
    
    #intmat = rand(Normal(0., sd_int), npop, npop)
    #intmat = ((abs.(intmat) .< sd_int) .& (abs.(intmat) .> 0.01)).*intmat
    #intmat = (abs.(intmat) .< 0.5*sd_int).*intmat
    #intmat -= (1.0 .+ intmat).*diagm(ones(npop))
    
    popdata = zeros(iter, npop)
    jmat = zeros(npop, npop, datasize)
    #initpop = rand(Uniform(0.1, 0.2), npop)
    popvec = copy(initpop)

    for i in 1:iter
        popnew = logisticmap(popvec, vec_growth, intmat; npop = npop, proc_noise = proc_noise, sd_noise = sd_noise)
        popdata[i, :] = popnew
        if i > iter - datasize
            jmat[:,:,i-iter+datasize] = jacobian_logistic(popnew, vec_growth, intmat)
        end
        popvec = popnew
    end
    (popdata, intmat, jmat, vec_growth, initpop)
end

# function to generate time series data from logistic map
function generate_logistic(;nmodel, r_growth, inter = 0.5, connectance = 0.3, sd_noise = 0.05, npop, datasize, iter, seed = 1234)
    seed!(seed)
    model = 1
    while model <= nmodel
        initpop = rand(Uniform(0.1, 0.5), npop)

        connect = zeros(Int, npop, npop)
        pos_neg = copy(connect)

        for sp in 1:npop-1
            connect[sp+1:npop, sp] .= sample([0,1], Weights([1 - connectance, connectance]), npop - sp)
            pos_neg[:, sp] .= sample([-1,1], Weights([0.5, 0.5]), npop)
            #connect[sp, sp] = 1
        end

        pos_neg[:, npop] .= sample([-1,1], Weights([0.5, 0.5]), npop)

        # assuming thath the interactions are bidiretional
        connect += LowerTriangular(connect)'
        
        vec_growth = r_growth .+ rand(Uniform(-0.5, 0.5), npop)
        intmat = connect.*pos_neg.*rand(Uniform(0.05, inter), npop, npop)
        # The class of interaction can be conpetitive, predatory or mutualistic
        intmat -= (1.0 .+ intmat).*diagm(ones(npop))

        data_logistic_noise = sim_logistic(npop; initpop=initpop, vec_growth=vec_growth, intmat=intmat, datasize=datasize, iter=iter, proc_noise=true, sd_noise = sd_noise)
        data_logistic_highnoise = sim_logistic(npop; initpop=initpop, vec_growth=vec_growth, intmat=intmat, datasize=datasize, iter=iter, proc_noise=true, sd_noise = 2*sd_noise)
        ts_logistic_noise = data_logistic_noise[1]
        intmat_logistic_noise = data_logistic_noise[2]
        ts_logistic_highnoise = data_logistic_highnoise[1]
        intmat_logistic_highnoise = data_logistic_highnoise[2]

        jmat_logistic_noise = data_logistic_noise[3]
        jmat_logistic_highnoise = data_logistic_highnoise[3]

        if (sum(mean(ts_logistic_noise, dims=1) .>0.05) == sum(mean(ts_logistic_highnoise, dims=1) .>0.05) == npop)
            writedlm("data/data_logistic_$(npop)_noise/data_logistic_$(npop)_noise_model$(model).csv",ts_logistic_noise[iter-datasize+1:end,:])
            writedlm("data/data_logistic_$(npop)_highnoise/data_logistic_$(npop)_highnoise_model$(model).csv",ts_logistic_highnoise[iter-datasize+1:end,:])
            jmat_data_noise = zeros(npop*datasize, npop)
            jmat_data_highnoise = zeros(npop*datasize, npop)
            for i in 1:datasize
                jmat_data_noise[(npop*(i-1)+1):npop*i, :] = jmat_logistic_noise[:, :, i]
                jmat_data_highnoise[(npop*(i-1)+1):npop*i, :] = jmat_logistic_highnoise[:, :, i]
            end
            
            writedlm("data/data_logistic_$(npop)_noise/jmat_logistic_$(npop)_noise_model$(model).csv", jmat_data_noise)
            writedlm("data/data_logistic_$(npop)_noise/intmat_logistic_$(npop)_noise_model$(model).csv", intmat_logistic_noise)
            writedlm("data/data_logistic_$(npop)_noise/initpop_$(npop)_noise_model$(model).csv", initpop)
            writedlm("data/data_logistic_$(npop)_noise/vec_growth_$(npop)_noise_model$(model).csv", vec_growth)

            writedlm("data/data_logistic_$(npop)_highnoise/jmat_logistic_$(npop)_highnoise_model$(model).csv", jmat_data_highnoise)
            writedlm("data/data_logistic_$(npop)_highnoise/intmat_logistic_$(npop)_highnoise_model$(model).csv", intmat_logistic_highnoise)
            writedlm("data/data_logistic_$(npop)_highnoise/initpop_$(npop)_highnoise_model$(model).csv", initpop)
            writedlm("data/data_logistic_$(npop)_highnoise/vec_growth_$(npop)_highnoise_model$(model).csv", vec_growth)

            println("model $model finished")
            model = model + 1
        end
    end
end
