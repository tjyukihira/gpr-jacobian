### Simulation of state-dependent Jacobian matrix
using Plots;
import Random.seed!
import StatsPlots.boxplot
using LaTeXStrings
using StatsBase

include("functions/GPR.jl")
include("functions/fun_fw5.jl")
include("functions/fun_sim_state-dependence.jl")
include("load_data_fw5.jl")
include("functions/get_jmat_discrete.jl")

# set seed for replication
seed!(2345)

vec_jmat_true_fw5_noise = copy(vec_jmat_scaled_fw5_noise)
vec_jmat_true_fw5_highnoise = copy(vec_jmat_scaled_fw5_highnoise)
vec_jmat_true_fw5_highnoise_obs = copy(vec_jmat_scaled_fw5_highnoise_obs)
vec_jmat_true_fw5_highnoise_highobs = copy(vec_jmat_scaled_fw5_highnoise_highobs)

for model in 1:nmodel
    for i in 1:datasize
        vec_jmat_true_fw5_noise[model][:,:,i] = get_jmat_discrete(vec_jmat_scaled_fw5_noise[model][:,:,i])
        vec_jmat_true_fw5_highnoise[model][:,:,i] = get_jmat_discrete(vec_jmat_scaled_fw5_highnoise[model][:,:,i])
        vec_jmat_true_fw5_highnoise_obs[model][:,:,i] = get_jmat_discrete(vec_jmat_scaled_fw5_highnoise_obs[model][:,:,i])
        vec_jmat_true_fw5_highnoise_highobs[model][:,:,i] = get_jmat_discrete(vec_jmat_scaled_fw5_highnoise_highobs[model][:,:,i])
    end
end

# ------------------------------------------------------------------------------------------------
# Process noise
vec_res_fw5_noise = Vector{NamedTuple}(undef, 50)

for i in 1:50
    vec_res_fw5_noise[i] = load_gpresults_jld2("results/fw5_noise/model$(i).jld2";npop=5, jacobian=true)
end

vec_data = data_fw5_noise_tr
vec_jmat_scaled = vec_jmat_true_fw5_noise
mu = mu_fw5_noise
sd = sd_fw5_noise

npop = size(vec_data[1], 2)

res_mean_strength = sim_mean_strength(vec_data, vec_res_fw5_noise; mu=mu, sd=sd)

vec_mean_strength_gpr = res_mean_strength.vec_mean_strength
vec_mean_strength_true = res_mean_strength.vec_mean_strength_true

vec_jmat_gpr = res_mean_strength.vec_jmat_gpr
vec_jmat_true = res_mean_strength.vec_jmat_true

stats_sim_noise = stats_sim_jmat(res_mean_strength, vec_res_fw5_noise, vec_jmat_scaled)

vec_diff_mean_strength_gpr = stats_sim_noise.vec_diff_mean_strength_gpr
vec_diff_mean_strength_true = stats_sim_noise.vec_diff_mean_strength_true
vec_diff_gpr = stats_sim_noise.vec_diff_gpr
vec_diff_true = stats_sim_noise.vec_diff_true

vec_sign_diff_mean = stats_sim_noise.vec_sign_diff_mean
vec_cor_diff = stats_sim_noise.vec_cor_diff

println("High level of process noise")
println("Accuracy (mean strength): ", round(median(vec_sign_diff_mean)*100, digits=1), "\n")
println("Accuracy (strengths): ", round(median(vec_cor_diff), digits=2), "\n")

i_med = argmin(abs.(vec_sign_diff_mean .- median(vec_sign_diff_mean)))
i_med_diff = argmin(abs.(vec_cor_diff .- median(vec_cor_diff)))

plt_diff = plot(vec_diff_mean_strength_gpr[i_med], linewidth=4, xlabel="Time", ylabel="Change of mean strength", label="GPR", legendfontsize=18, legend=:bottomright)
plt_diff = plot!(plt_diff, vec_diff_mean_strength_true[i_med], linewidth=4, label="Ground truth")
plt_diff = annotate!(plt_diff, [75], [0.02], text("Accuracy: $(round(vec_sign_diff_mean[i_med]*100, digits=1))%", 20))
plt_diff = hline!(plt_diff, [0], linestyle=:dash, colour=:black, linewidth=3, label=false)

sct_diff = scatter_diff(vec_diff_gpr, vec_diff_true; model=i_med_diff, nsp=5, datasize=100)

box_sign_diff_mean = boxplot(vec_sign_diff_mean*100, xticks=false, label="Accuracy (median: $(round(median(vec_sign_diff_mean)*100, digits=1))%)", ylim=(minimum(vec_sign_diff_mean*100)-5,100), legend=:outertop, legendfontsize=19)
box_sign_diff_mean = annotate!(box_sign_diff_mean, [1], [93], text("$(sum(vec_sign_diff_mean.>0.5))/50", 20))

box_cor_diff = boxplot(vec_cor_diff, xticks=false, label="ρ (median: $(round(median(vec_cor_diff), digits=2)))", ylim=(0,1), legend=:outertop, legendfontsize=19)

plt1 = plot(plt_diff, box_sign_diff_mean, layout=grid(1,2, widths=(7/10,3/10)), size=(1600, 900), tickfontsize = 16, labelfontsize = 19, margin=0.8Plots.cm)

savefig("fig/figureS3/figS3a_sim_jmat_fw5_noise.png")

plt2 = plot(sct_diff, box_cor_diff, layout=grid(1,2, widths=(7/10,3/10)), size=(1600, 900), tickfontsize = 16, labelfontsize = 19, margin=0.8Plots.cm)
savefig("fig/figureS3/figS3b_diff_jmat_fw5_noise.png")

# ------------------------------------------------------------------------------------------------
# High level of process noise
vec_res_fw5_highnoise = Vector{NamedTuple}(undef, 50)

for i in 1:50
    vec_res_fw5_highnoise[i] = load_gpresults_jld2("results/fw5_highnoise/model$(i).jld2";npop=5, jacobian=true)
end

vec_data = data_fw5_highnoise_tr
vec_jmat_scaled = vec_jmat_true_fw5_highnoise
mu = mu_fw5_highnoise
sd = sd_fw5_highnoise

npop = size(vec_data[1], 2)

res_mean_strength = sim_mean_strength(vec_data, vec_res_fw5_highnoise; mu=mu, sd=sd)

vec_mean_strength_gpr = res_mean_strength.vec_mean_strength
vec_mean_strength_true = res_mean_strength.vec_mean_strength_true

vec_jmat_gpr = res_mean_strength.vec_jmat_gpr
vec_jmat_true = res_mean_strength.vec_jmat_true

stats_sim_highnoise = stats_sim_jmat(res_mean_strength, vec_res_fw5_highnoise, vec_jmat_scaled)

vec_diff_mean_strength_gpr = stats_sim_highnoise.vec_diff_mean_strength_gpr
vec_diff_mean_strength_true = stats_sim_highnoise.vec_diff_mean_strength_true
vec_diff_gpr = stats_sim_highnoise.vec_diff_gpr
vec_diff_true = stats_sim_highnoise.vec_diff_true

vec_sign_diff_mean = stats_sim_highnoise.vec_sign_diff_mean
vec_cor_diff = stats_sim_highnoise.vec_cor_diff

println("High level of process noise")
println("Accuracy (mean strength): ", round(median(vec_sign_diff_mean)*100, digits=1), "\n")
println("Accuracy (strengths): ", round(median(vec_cor_diff), digits=2), "\n")

i_med = argmin(abs.(vec_sign_diff_mean .- median(vec_sign_diff_mean)))
i_med_diff = argmin(abs.(vec_cor_diff .- median(vec_cor_diff)))

plt_diff = plot(vec_diff_mean_strength_gpr[i_med], linewidth=4, xlabel="Time", ylabel="Change of mean strength", label="GPR", legendfontsize=18, legend=:bottomright)
plt_diff = plot!(plt_diff, vec_diff_mean_strength_true[i_med], linewidth=4, label="Ground truth")
plt_diff = annotate!(plt_diff, [25], [0.02], text("Accuracy: $(round(vec_sign_diff_mean[i_med]*100, digits=1))%", 20))
plt_diff = hline!(plt_diff, [0], linestyle=:dash, colour=:black, linewidth=3, label=false)

sct_diff = scatter_diff(vec_diff_gpr, vec_diff_true; model=i_med_diff, nsp=5, datasize=100)

box_sign_diff_mean = boxplot(vec_sign_diff_mean*100, xticks=false, label="Accuracy (median: $(round(median(vec_sign_diff_mean)*100, digits=1))%)", ylim=(minimum(vec_sign_diff_mean*100)-5,100), legend=:outertop, legendfontsize=19)
box_sign_diff_mean = annotate!(box_sign_diff_mean, [1], [90], text("$(sum(vec_sign_diff_mean.>0.5))/50", 20))

box_cor_diff = boxplot(vec_cor_diff, xticks=false, label="ρ (median: $(round(median(vec_cor_diff), digits=2)))", ylim=(0,1), legend=:outertop, legendfontsize=19)

plt1 = plot(plt_diff, box_sign_diff_mean, layout=grid(1,2, widths=(7/10,3/10)), size=(1600, 900), tickfontsize = 16, labelfontsize = 19, margin=0.8Plots.cm)

savefig("fig/figureS3/figS3c_sim_jmat_fw5_highnoise.png")

plt2 = plot(sct_diff, box_cor_diff, layout=grid(1,2, widths=(7/10,3/10)), size=(1600, 900), tickfontsize = 16, labelfontsize = 19, margin=0.8Plots.cm)
savefig("fig/figureS3/figS3d_diff_jmat_fw5_highnoise.png")

# ------------------------------------------------------------------------------------------------
# Process and observational noise (high and modest)
vec_res_fw5_highnoise_obs = Vector{NamedTuple}(undef, 50)

for i in 1:50
    vec_res_fw5_highnoise_obs[i] = load_gpresults_jld2("results/fw5_highnoise_obs/model$(i).jld2";npop=5, jacobian=true)
end

vec_data_obs = data_fw5_highnoise_obs_tr
vec_data_unscaled = data_fw5_highnoise
vec_jmat_scaled = vec_jmat_true_fw5_highnoise_obs
mu = mu_fw5_highnoise_obs
sd = sd_fw5_highnoise_obs

npop = size(vec_data_obs[1], 2)

res_mean_strength = sim_mean_strength_obs(vec_data_unscaled, vec_data_obs, vec_res_fw5_highnoise_obs; mu_obs=mu, sd_obs=sd)

vec_mean_strength_gpr = res_mean_strength.vec_mean_strength
vec_mean_strength_true = res_mean_strength.vec_mean_strength_true

vec_jmat_gpr = res_mean_strength.vec_jmat_gpr
vec_jmat_true = res_mean_strength.vec_jmat_true

stats_sim_highnoise_obs = stats_sim_jmat(res_mean_strength, vec_res_fw5_highnoise_obs, vec_jmat_scaled)

vec_diff_mean_strength_gpr = stats_sim_highnoise_obs.vec_diff_mean_strength_gpr
vec_diff_mean_strength_true = stats_sim_highnoise_obs.vec_diff_mean_strength_true
vec_diff_gpr = stats_sim_highnoise_obs.vec_diff_gpr
vec_diff_true = stats_sim_highnoise_obs.vec_diff_true

vec_sign_diff_mean = stats_sim_highnoise_obs.vec_sign_diff_mean
vec_cor_diff = stats_sim_highnoise_obs.vec_cor_diff

println("Process and observational noise (high and modest)")
println("Accuracy (mean strength): ", round(median(vec_sign_diff_mean)*100, digits=1), "\n")
println("Accuracy (strengths): ", round(median(vec_cor_diff), digits=2), "\n")

i_med = argmin(abs.(vec_sign_diff_mean .- median(vec_sign_diff_mean)))
i_med_diff = argmin(abs.(vec_cor_diff .- median(vec_cor_diff)))

plt_diff = plot(vec_diff_mean_strength_gpr[i_med], linewidth=4, xlabel="Time", ylabel="Change of mean strength", label="GPR", legendfontsize=18, legend=:bottomright)
plt_diff = plot!(plt_diff, vec_diff_mean_strength_true[i_med], linewidth=4, label="Ground truth")
plt_diff = annotate!(plt_diff, [75], [0.022], text("Accuracy: $(round(vec_sign_diff_mean[i_med]*100, digits=1))%", 20))
plt_diff = hline!(plt_diff, [0], linestyle=:dash, colour=:black, linewidth=3, label=false)

sct_diff = scatter_diff(vec_diff_gpr, vec_diff_true; model=i_med_diff, nsp=5, datasize=100)

box_sign_diff_mean = boxplot(vec_sign_diff_mean*100, xticks=false, label="Accuracy (median: $(round(median(vec_sign_diff_mean)*100, digits=1))%)", ylim=(minimum(vec_sign_diff_mean*100)-5,100), legend=:outertop, legendfontsize=19)
box_sign_diff_mean = annotate!(box_sign_diff_mean, [1], [90], text("$(sum(vec_sign_diff_mean.>0.5))/50", 20))

box_cor_diff = boxplot(vec_cor_diff, xticks=false, label="ρ (median: $(round(median(vec_cor_diff), digits=2)))", ylim=(0,1), legend=:outertop, legendfontsize=19)

plt1 = plot(plt_diff, box_sign_diff_mean, layout=grid(1,2, widths=(7/10,3/10)), size=(1600, 900), tickfontsize = 16, labelfontsize = 19, margin=0.8Plots.cm)

savefig("fig/figureS3/figS3e_sim_jmat_fw5_highnoise_obs.png")

plt2 = plot(sct_diff, box_cor_diff, layout=grid(1,2, widths=(7/10,3/10)), size=(1600, 900), tickfontsize = 16, labelfontsize = 19, margin=0.8Plots.cm)
savefig("fig/figureS3/figS3f_diff_jmat_fw5_highnoise_obs.png")

# ------------------------------------------------------------------------------------------------
# High level of process and observational noise
vec_res_fw5_highnoise_highobs = Vector{NamedTuple}(undef, 50)

for i in 1:50
    vec_res_fw5_highnoise_highobs[i] = load_gpresults_jld2("results/fw5_highnoise_highobs/model$(i).jld2";npop=5, jacobian=true)
end

vec_data_highobs = data_fw5_highnoise_highobs_tr
vec_data_unscaled = data_fw5_highnoise
vec_jmat_scaled = vec_jmat_true_fw5_highnoise_highobs
mu = mu_fw5_highnoise_highobs
sd = sd_fw5_highnoise_highobs

npop = size(vec_data_obs[1], 2)

res_mean_strength = sim_mean_strength_obs(vec_data_unscaled, vec_data_highobs, vec_res_fw5_highnoise_highobs; mu_obs=mu, sd_obs=sd)

vec_mean_strength_gpr = res_mean_strength.vec_mean_strength
vec_mean_strength_true = res_mean_strength.vec_mean_strength_true

vec_jmat_gpr = res_mean_strength.vec_jmat_gpr
vec_jmat_true = res_mean_strength.vec_jmat_true

stats_sim_highnoise_highobs = stats_sim_jmat(res_mean_strength, vec_res_fw5_highnoise_highobs, vec_jmat_scaled)

vec_diff_mean_strength_gpr = stats_sim_highnoise_highobs.vec_diff_mean_strength_gpr
vec_diff_mean_strength_true = stats_sim_highnoise_highobs.vec_diff_mean_strength_true
vec_diff_gpr = stats_sim_highnoise_highobs.vec_diff_gpr
vec_diff_true = stats_sim_highnoise_highobs.vec_diff_true

vec_sign_diff_mean = stats_sim_highnoise_highobs.vec_sign_diff_mean
vec_cor_diff = stats_sim_highnoise_highobs.vec_cor_diff

println("high level of process and observational noise")
println("Accuracy (mean strength): ", round(median(vec_sign_diff_mean)*100, digits=1), "\n")
println("Accuracy (strengths): ", round(median(vec_cor_diff), digits=2), "\n")

i_med = argmin(abs.(vec_sign_diff_mean .- median(vec_sign_diff_mean)))
i_med_diff = argmin(abs.(vec_cor_diff .- median(vec_cor_diff)))

plt_diff = plot(vec_diff_mean_strength_gpr[i_med], linewidth=4, xlabel="Time", ylabel="Change of mean strength", label="GPR", legendfontsize=18, legend=:bottomright)
plt_diff = plot!(plt_diff, vec_diff_mean_strength_true[i_med], linewidth=4, label="Ground truth")
plt_diff = annotate!(plt_diff, [75], [0.02], text("Accuracy: $(round(vec_sign_diff_mean[i_med]*100, digits=1))%", 20))
plt_diff = hline!(plt_diff, [0], linestyle=:dash, colour=:black, linewidth=3, label=false)

sct_diff = scatter_diff(vec_diff_gpr, vec_diff_true; model=i_med_diff, nsp=5, datasize=100)

box_sign_diff_mean = boxplot(vec_sign_diff_mean*100, xticks=false, label="Accuracy (median: $(round(median(vec_sign_diff_mean)*100, digits=1))%)", ylim=(minimum(vec_sign_diff_mean*100)-5,100), legend=:outertop, legendfontsize=19)
box_sign_diff_mean = annotate!(box_sign_diff_mean, [1], [90], text("$(sum(vec_sign_diff_mean.>0.5))/50", 20))

box_cor_diff = boxplot(vec_cor_diff, xticks=false, label="ρ (median: $(round(median(vec_cor_diff), digits=2)))", ylim=(0,1), legend=:outertop, legendfontsize=19)

plt1 = plot(plt_diff, box_sign_diff_mean, layout=grid(1,2, widths=(7/10,3/10)), size=(1600, 900), tickfontsize = 16, labelfontsize = 19, margin=0.8Plots.cm)

savefig("fig/figureS3/figS3g_sim_jmat_fw5_highnoise_highobs.png")

plt2 = plot(sct_diff, box_cor_diff, layout=grid(1,2, widths=(7/10,3/10)), size=(1600, 900), tickfontsize = 16, labelfontsize = 19, margin=0.8Plots.cm)
savefig("fig/figureS3/figS3h_diff_jmat_fw5_highnoise_highobs.png")
