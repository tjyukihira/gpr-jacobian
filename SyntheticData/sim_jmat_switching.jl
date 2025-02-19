### Simulation of state-dependent Jacobian matrix
using Plots;
import Random.seed!
import StatsPlots.boxplot
using LaTeXStrings
using StatsBase

include("functions/GPR.jl")
include("functions/fun_prey_switching.jl")
include("functions/fun_sim_state-dependence.jl")
include("functions/get_jmat_discrete.jl")
include("load_data_switching.jl")


# set seed for replication
seed!(2345)

vec_jmat_true_switching_noise = copy(vec_jmat_scaled_switching_noise)
vec_jmat_true_switching_highnoise = copy(vec_jmat_scaled_switching_highnoise)
vec_jmat_true_switching_highnoise_obs = copy(vec_jmat_scaled_switching_highnoise_obs)
vec_jmat_true_switching_highnoise_highobs = copy(vec_jmat_scaled_switching_highnoise_highobs)

for model in 1:nmodel
    for i in 1:datasize
        vec_jmat_true_switching_noise[model][:,:,i] = get_jmat_discrete(vec_jmat_scaled_switching_noise[model][:,:,i])
        vec_jmat_true_switching_highnoise[model][:,:,i] = get_jmat_discrete(vec_jmat_scaled_switching_highnoise[model][:,:,i])
        vec_jmat_true_switching_highnoise_obs[model][:,:,i] = get_jmat_discrete(vec_jmat_scaled_switching_highnoise_obs[model][:,:,i])
        vec_jmat_true_switching_highnoise_highobs[model][:,:,i] = get_jmat_discrete(vec_jmat_scaled_switching_highnoise_highobs[model][:,:,i])
    end
end

# ------------------------------------------------------------------------------------------------
# Process noise
vec_res_switching_noise = Vector{NamedTuple}(undef, 50)

for i in 1:50
    vec_res_switching_noise[i] = load_gpresults_jld2("results/switching_noise/model$(i).jld2";npop=5, jacobian=true)
end

vec_data = data_switching_noise_tr
vec_jmat_scaled = vec_jmat_true_switching_noise
mu = mu_switching_noise
sd = sd_switching_noise

npop = size(vec_data[1], 2)

res_mean_strength = sim_mean_strength_switching(vec_data, vec_res_switching_noise, pref; mu=mu, sd=sd)

vec_mean_strength_gpr = res_mean_strength.vec_mean_strength
vec_mean_strength_true = res_mean_strength.vec_mean_strength_true

vec_jmat_gpr = res_mean_strength.vec_jmat_gpr
vec_jmat_true = res_mean_strength.vec_jmat_true

stats_sim_noise = stats_sim_jmat(res_mean_strength, vec_res_switching_noise, vec_jmat_scaled)

vec_diff_mean_strength_gpr = stats_sim_noise.vec_diff_mean_strength_gpr
vec_diff_mean_strength_true = stats_sim_noise.vec_diff_mean_strength_true
vec_diff_gpr = stats_sim_noise.vec_diff_gpr
vec_diff_true = stats_sim_noise.vec_diff_true

vec_sign_diff_mean = stats_sim_noise.vec_sign_diff_mean
vec_cor_diff = stats_sim_noise.vec_cor_diff


println("Process noise")
println("Accuracy (mean strength): ", round(median(vec_sign_diff_mean)*100, digits=1), "\n")
println("Accuracy (strengths): ", round(median(vec_cor_diff), digits=2), "\n")

i_med = argmin(abs.(vec_sign_diff_mean .- median(vec_sign_diff_mean)))
i_med_diff = argmin(abs.(vec_cor_diff .- median(vec_cor_diff)))

plt_diff = plot(vec_diff_mean_strength_gpr[i_med], linewidth=4, xlabel="Time", ylabel="Change of mean strength", label="GPR", legendfontsize=18, legend=:bottomright)
plt_diff = plot!(plt_diff, vec_diff_mean_strength_true[i_med], linewidth=4, label="Ground truth")
plt_diff = annotate!(plt_diff, [75], [0.02], text("Accuracy: $(round(vec_sign_diff_mean[i_med]*100, digits=1))%", 20))
plt_diff = hline!(plt_diff, [0], linestyle=:dash, colour=:black, linewidth=3, label=false)

sct_diff = scatter_diff(vec_diff_gpr, vec_diff_true; model=i_med_diff, nsp=5, datasize=100)

box_sign_diff_mean = boxplot(vec_sign_diff_mean*100, xticks=false, label="Accuracy (median: $(round(median(vec_sign_diff_mean)*100, digits=1))%)", ylim=(50,100), legend=:outertop, legendfontsize=19)
box_sign_diff_mean = annotate!(box_sign_diff_mean, [1], [95], text("$(sum(vec_sign_diff_mean.>0.5))/50", 20))

box_cor_diff = boxplot(vec_cor_diff, xticks=false, label="ρ (median: $(round(median(vec_cor_diff), digits=2)))", ylim=(0,1), legend=:outertop, legendfontsize=19)

plt1 = plot(plt_diff, box_sign_diff_mean, layout=grid(1,2, widths=(7/10,3/10)), size=(1600, 900), tickfontsize = 16, labelfontsize = 19, margin=0.8Plots.cm)

savefig("fig/figureS3/figS3i_sim_jmat_switching_noise.png")

plt2 = plot(sct_diff, box_cor_diff, layout=grid(1,2, widths=(7/10,3/10)), size=(1600, 900), tickfontsize = 16, labelfontsize = 19, margin=0.8Plots.cm)
savefig("fig/figureS3/figS3j_diff_jmat_switching_noise.png")

# ------------------------------------------------------------------------------------------------
# High level of process noise
vec_res_switching_highnoise = Vector{NamedTuple}(undef, 50)

for i in 1:50
    vec_res_switching_highnoise[i] = load_gpresults_jld2("results/switching_highnoise/model$(i).jld2";npop=5, jacobian=true)
end

vec_data = data_switching_highnoise_tr
vec_jmat_scaled = vec_jmat_true_switching_highnoise
mu = mu_switching_highnoise
sd = sd_switching_highnoise

npop = size(vec_data[1], 2)

res_mean_strength = sim_mean_strength_switching(vec_data, vec_res_switching_highnoise, pref; mu=mu, sd=sd)

vec_mean_strength_gpr = res_mean_strength.vec_mean_strength
vec_mean_strength_true = res_mean_strength.vec_mean_strength_true

vec_jmat_gpr = res_mean_strength.vec_jmat_gpr
vec_jmat_true = res_mean_strength.vec_jmat_true

stats_sim_highnoise = stats_sim_jmat(res_mean_strength, vec_res_switching_highnoise, vec_jmat_scaled)

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

plt_diff = plot(vec_diff_mean_strength_gpr[i_med], linewidth=4, xlabel="Time", ylabel="Change of mean strength", label="GPR", legendfontsize=18, legend=:topleft)
plt_diff = plot!(plt_diff, vec_diff_mean_strength_true[i_med], linewidth=4, label="Ground truth")
plt_diff = annotate!(plt_diff, [50], [0.02], text("Accuracy: $(round(vec_sign_diff_mean[i_med]*100, digits=1))%", 20))
plt_diff = hline!(plt_diff, [0], linestyle=:dash, colour=:black, linewidth=3, label=false)

sct_diff = scatter_diff(vec_diff_gpr, vec_diff_true; model=i_med_diff, nsp=5, datasize=100)

box_sign_diff_mean = boxplot(vec_sign_diff_mean*100, xticks=false, label="Accuracy (median: $(round(median(vec_sign_diff_mean)*100, digits=1))%)", ylim=(50,100), legend=:outertop, legendfontsize=19)
box_sign_diff_mean = annotate!(box_sign_diff_mean, [1], [95], text("$(sum(vec_sign_diff_mean.>0.5))/50", 20))

box_cor_diff = boxplot(vec_cor_diff, xticks=false, label="ρ (median: $(round(median(vec_cor_diff), digits=2)))", ylim=(0,1), legend=:outertop, legendfontsize=19)

plt1 = plot(plt_diff, box_sign_diff_mean, layout=grid(1,2, widths=(7/10,3/10)), size=(1600, 900), tickfontsize = 15, labelfontsize = 18, margin=0.8Plots.cm)
savefig("fig/figureS3/figS3k_sim_jmat_switching_highnoise.png")

plt_fig2cd = annotate!(plt1, -20, 0.027, text(L"(\textit{c})\hspace{10}"*"Food web 2 (κ = 0.1)", :left, 22), subplot=1)
plt_fig2cd = annotate!(plt1, 105, 0.027, text(L"(\textit{d})", :left, 22), subplot=1)

plt2 = plot(sct_diff, box_cor_diff, layout=grid(1,2, widths=(7/10,3/10)), size=(1600, 900), tickfontsize = 15, labelfontsize = 18, margin=0.8Plots.cm)
savefig("fig/figureS3/figS3l_diff_jmat_switching_highnoise.png")

plt_fig2ab = annotate!(plt2, -0.24, 0.46, text(L"(\textit{a})\hspace{10}"*"Food web 2 (κ = 0.1)", :left, 22), subplot=1)
plt_fig2ab = annotate!(plt2, 0.34, 0.46, text(L"(\textit{b})", :left, 22), subplot=1)

plt_fig2 = plot(plt_fig2ab, plt_fig2cd, layout=(2,1), size=(1500,1650), margin=0.8Plots.cm, top_margin=0.9Plots.cm, bottom_margin=1.0Plots.cm)
savefig("fig/figure3")

# ------------------------------------------------------------------------------------------------
# Process and observational noise (high and modest)
vec_res_switching_highnoise_obs = Vector{NamedTuple}(undef, 50)

for i in 1:50
    vec_res_switching_highnoise_obs[i] = load_gpresults_jld2("results/switching_highnoise_obs/model$(i).jld2";npop=5, jacobian=true)
end

vec_data_obs = data_switching_highnoise_obs_tr
vec_data_unscaled = data_switching_highnoise
vec_jmat_scaled = vec_jmat_true_switching_highnoise_obs
mu = mu_switching_highnoise_obs
sd = sd_switching_highnoise_obs

npop = size(vec_data_obs[1], 2)

res_mean_strength = sim_mean_strength_switching_obs(vec_data_unscaled, vec_data_obs, vec_res_switching_highnoise_obs, pref; mu_obs=mu, sd_obs=sd)

vec_mean_strength_gpr = res_mean_strength.vec_mean_strength
vec_mean_strength_true = res_mean_strength.vec_mean_strength_true

vec_jmat_gpr = res_mean_strength.vec_jmat_gpr
vec_jmat_true = res_mean_strength.vec_jmat_true

stats_sim_highnoise_obs = stats_sim_jmat(res_mean_strength, vec_res_switching_highnoise_obs, vec_jmat_scaled)

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

plt_diff = plot(vec_diff_mean_strength_gpr[i_med], linewidth=4, xlabel="Time", ylabel="Change of mean strength", label="GPR", legendfontsize=18, legend=:bottom)
plt_diff = plot!(plt_diff, vec_diff_mean_strength_true[i_med], linewidth=4, label="Ground truth")
plt_diff = annotate!(plt_diff, [75], [0.02], text("Accuracy: $(round(vec_sign_diff_mean[i_med]*100, digits=1))%", 20))
plt_diff = hline!(plt_diff, [0], linestyle=:dash, colour=:black, linewidth=3, label=false)

sct_diff = scatter_diff(vec_diff_gpr, vec_diff_true; model=i_med_diff, nsp=5, datasize=100)

box_sign_diff_mean = boxplot(vec_sign_diff_mean*100, xticks=false, label="Accuracy (median: $(round(median(vec_sign_diff_mean)*100, digits=1))%)", ylim=(30,100), legend=:outertop, legendfontsize=19)
box_sign_diff_mean = annotate!(box_sign_diff_mean, [1], [95], text("$(sum(vec_sign_diff_mean.>0.5))/50", 20))

box_cor_diff = boxplot(vec_cor_diff, xticks=false, label="ρ (median: $(round(median(vec_cor_diff), digits=2)))", ylim=(0,1), legend=:outertop, legendfontsize=19)

plt1 = plot(plt_diff, box_sign_diff_mean, layout=grid(1,2, widths=(7/10,3/10)), size=(1600, 900), tickfontsize = 16, labelfontsize = 19, margin=0.8Plots.cm)

savefig("fig/figureS3/figS3m_sim_jmat_switching_highnoise_obs.png")

plt2 = plot(sct_diff, box_cor_diff, layout=grid(1,2, widths=(7/10,3/10)), size=(1600, 900), tickfontsize = 16, labelfontsize = 19, margin=0.8Plots.cm)
savefig("fig/figureS3/figS3n_diff_jmat_switching_highnoise_obs.png")

# ------------------------------------------------------------------------------------------------
# Process and observational noise (Both high)
vec_res_switching_highnoise_highobs = Vector{NamedTuple}(undef, 50)

for i in 1:50
    vec_res_switching_highnoise_highobs[i] = load_gpresults_jld2("results/switching_highnoise_highobs/model$(i).jld2";npop=5, jacobian=true)
end

vec_data_highobs = data_switching_highnoise_highobs_tr
vec_data_unscaled = data_switching_highnoise
vec_jmat_scaled = vec_jmat_true_switching_highnoise_highobs
mu = mu_switching_highnoise_highobs
sd = sd_switching_highnoise_highobs

npop = size(vec_data_obs[1], 2)

res_mean_strength = sim_mean_strength_switching_obs(vec_data_unscaled, vec_data_highobs, vec_res_switching_highnoise_highobs, pref; mu_obs=mu, sd_obs=sd)

vec_mean_strength_gpr = res_mean_strength.vec_mean_strength
vec_mean_strength_true = res_mean_strength.vec_mean_strength_true

vec_jmat_gpr = res_mean_strength.vec_jmat_gpr
vec_jmat_true = res_mean_strength.vec_jmat_true

stats_sim_highnoise_highobs = stats_sim_jmat(res_mean_strength, vec_res_switching_highnoise_highobs, vec_jmat_scaled)

vec_diff_mean_strength_gpr = stats_sim_highnoise_highobs.vec_diff_mean_strength_gpr
vec_diff_mean_strength_true = stats_sim_highnoise_highobs.vec_diff_mean_strength_true
vec_diff_gpr = stats_sim_highnoise_highobs.vec_diff_gpr
vec_diff_true = stats_sim_highnoise_highobs.vec_diff_true

vec_sign_diff_mean = stats_sim_highnoise_highobs.vec_sign_diff_mean
vec_cor_diff = stats_sim_highnoise_highobs.vec_cor_diff

println("High level of process and observational noise")
println("Accuracy (mean strength): ", round(median(vec_sign_diff_mean)*100, digits=1), "\n")
println("Accuracy (strengths): ", round(median(vec_cor_diff), digits=2), "\n")

i_med = argmin(abs.(vec_sign_diff_mean .- median(vec_sign_diff_mean)))
i_med_diff = argmin(abs.(vec_cor_diff .- median(vec_cor_diff)))

plt_diff = plot(vec_diff_mean_strength_gpr[i_med], linewidth=4, xlabel="Time", ylabel="Change of mean strength", label="GPR", legendfontsize=18, legend=:bottom)
plt_diff = plot!(plt_diff, vec_diff_mean_strength_true[i_med], linewidth=4, label="Ground truth")
plt_diff = annotate!(plt_diff, [75], [0.02], text("Accuracy: $(round(vec_sign_diff_mean[i_med]*100, digits=1))%", 20))
plt_diff = hline!(plt_diff, [0], linestyle=:dash, colour=:black, linewidth=3, label=false)

sct_diff = scatter_diff(vec_diff_gpr, vec_diff_true; model=i_med_diff, nsp=5, datasize=100)

box_sign_diff_mean = boxplot(vec_sign_diff_mean*100, xticks=false, label="Accuracy (median: $(round(median(vec_sign_diff_mean)*100, digits=1))%)", ylim=(50,100), legend=:outertop, legendfontsize=19)
box_sign_diff_mean = annotate!(box_sign_diff_mean, [1], [95], text("$(sum(vec_sign_diff_mean.>0.5))/50", 20))

box_cor_diff = boxplot(vec_cor_diff, xticks=false, label="Accuracy (median: $(round(median(vec_cor_diff), digits=2)))", ylim=(0,1), legend=:outertop, legendfontsize=19)

plt1 = plot(plt_diff, box_sign_diff_mean, layout=grid(1,2, widths=(7/10,3/10)), size=(1600, 900), tickfontsize = 16, labelfontsize = 19, margin=0.8Plots.cm)

savefig("fig/figureS3/figS3o_sim_jmat_switching_highnoise_highobs.png")

plt2 = plot(sct_diff, box_cor_diff, layout=grid(1,2, widths=(7/10,3/10)), size=(1600, 900), tickfontsize = 16, labelfontsize = 19, margin=0.8Plots.cm)

savefig("fig/figureS3/figS3p_diff_jmat_switching_highnoise_highobs.png")
