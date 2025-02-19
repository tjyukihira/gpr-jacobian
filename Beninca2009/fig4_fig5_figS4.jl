### Beninca et. al. (2009) Ecol. Lett.
using DelimitedFiles
using Plots
using StatsBase
using Distributions
using LaTeXStrings

include("functions/GPR.jl")
include("functions/Auxiliar.jl")
include("functions/sim_state_dependence.jl")

d_baltic = readdlm("data/Beninca_2009.csv", ',', header = true)

# use data after April, 1991 to ensure stationarity
ts_plankton = d_baltic[1][d_baltic[1][:, 1] .> 263, 2:end]
v_days = d_baltic[1][d_baltic[1][:, 1] .> 263, 1]

ts_lib = ts_plankton[1:end, :] # the data is not scaled since the original data is already scaled

# load the GPR results
res_gpr = load_gpresults_jld2("results/res_gpr.jld2", npop = 4)

# in-sample correlation to check model fitting
cor1 = cor(ts_lib[2:end, 1], res_gpr.loocv[1].mu_loo)
cor2 = cor(ts_lib[2:end, 2], res_gpr.loocv[2].mu_loo)
cor3 = cor(ts_lib[2:end, 3], res_gpr.loocv[3].mu_loo)
cor4 = cor(ts_lib[2:end, 4], res_gpr.loocv[4].mu_loo)

# plot in-sample prediction with LOOCV
plt1 = plot(ts_lib[2:end, 1], title = "ρ = $cor1")
plt1 = plot!(plt1, res_gpr.loocv[1].mu_loo)

plt2 = plot(ts_lib[2:end, 2], title = "ρ = $cor2")
plt2 = plot!(plt2, res_gpr.loocv[2].mu_loo)

plt3 = plot(ts_lib[2:end, 3])
plt3 = plot!(plt3, res_gpr.loocv[3].mu_loo, title = "ρ = $cor3")

plt4 = plot(ts_lib[2:end, 4])
plt4 = plot!(plt4, res_gpr.loocv[4].mu_loo, title = "ρ = $cor4")

plot(plt1, plt2, plt3, plt4, layout=(2,2))

# plot the results on each interaction strengths
p12 = plot_strengths(res_gpr, v_days[1:end-1]; row = 1, column = 2, title = "Calanoids ⬅ Rotifers", digits = 3)

p13 = plot_strengths(res_gpr, v_days[1:end-1]; row = 1, column = 3, title = "Calanoids ⬅ Nano", digits = 3) 

p14 = plot_strengths(res_gpr, v_days[1:end-1]; row = 1, column = 4, title = "Calanoids ⬅ Pico", digits = 3) 

p21 = plot_strengths(res_gpr, v_days[1:end-1]; row = 2, column = 1, title = "Rotifers ⬅ Calanoids", digits = 3) 

p23 = plot_strengths(res_gpr, v_days[1:end-1]; row = 2, column = 3, title = "Rotifers ⬅ Nano", digits = 3) 

p24 = plot_strengths(res_gpr, v_days[1:end-1]; row = 2, column = 4, title = "Rotifers ⬅ Pico", digits = 3) 

p31 = plot_strengths(res_gpr, v_days[1:end-1]; row = 3, column = 1, title = "Nano ⬅ Calanoids", digits = 3) 

p32 = plot_strengths(res_gpr, v_days[1:end-1]; row = 3, column = 2, title = "Nano ⬅ Rotifers", digits = 3) 

p34 = plot_strengths(res_gpr, v_days[1:end-1]; row = 3, column = 4, title = "Nano ⬅ Pico", digits = 3) 

p41 = plot_strengths(res_gpr, v_days[1:end-1]; row = 4, column = 1, title = "Pico ⬅ Calanoids", digits = 3) 

p42 = plot_strengths(res_gpr, v_days[1:end-1]; row = 4, column = 2, title = "Pico ⬅ Rotifers", digits = 3) 

p43 = plot_strengths(res_gpr, v_days[1:end-1]; row = 4, column = 3, title = "Pico ⬅ Nano", digits = 3) 

# plot all of the results and save them
#plot(p34, p43, layout = (1,2), title = [L"(\textit{a})" L"(\textit{b})"], size = (1400, 500), titlelocation=:left, legendfontsize=9, titlefontsize=11, tickfontsize=8, labelfontsize=10, legend=:outertop, margin = 1.5Plots.cm, top_margin = 2Plots.px)
plot(p34, p43, layout = (1,2), size = (1500, 600), plot_title = L"(\textit{a}) \hspace{29} (\textit{b})", plot_titlelocation=(0.21,0.3), plot_titilevspan=0.22, plot_titlefontsize=21, legendfontsize=16, titlefontsize=21, tickfontsize=13, labelfontsize=15, legend=:outertop, margin = 1.8Plots.cm, top_margin = 25Plots.px)
savefig("fig/fig4_compet.png")

plot(p12, p13, p14, p21, p23, p24, p31, p32, p34, p41, p42, p43, layout = (4,3), size = (3000, 3000), legendfontsize=20, titlefontsize=27, tickfontsize=17, labelfontsize=19, legend=:outertop, margin = 2.5Plots.cm, top_margin = 8Plots.px)
savefig("fig/figS4_strengths.png")

# mean and sd of time series for scaling
mu = mean(ts_lib, dims=1)[1,:]
sd = std(ts_lib, dims=1)[1,:]

# simulate the depence of competition on Rotifers
res_sim_jmat_inc_roti = sim_state_dependence(res_gpr, ts_lib, mu, sd; target_inc = 2)

mean34_init = mean(res_gpr.jmat[3, 4, 1:end-1]); mean34_sim = mean(res_sim_jmat_inc_roti[3, 4, :])
plt_j34_inc_roti = plot(v_days[1:end-1], res_gpr.jmat[3, 4, 1:end-1], label = "In-sample "*replace("(mean $(round(mean34_init, digits=3)))", "-" => L"-"), ylim=(-0.24, 0), xlabel = "Time (days)", ylabel = "Nano ⬅ Pico", title = "Nano ⬅ Pico (+20% Rotifers)", titlefontsize = 14, legendfontsize=10)
plt_j34_inc_roti = plot!(plt_j34_inc_roti, v_days[1:end-1],  res_sim_jmat_inc_roti[3,4,:], label = replace("+20% Rotifers (mean $(round(mean34_sim, digits=3)))", "-" => L"-"))
plt_j34_inc_roti = annotate!(plt_j34_inc_roti, [500], [-0.2], text("ΔIS: $(round(-((mean34_sim-mean34_init)/mean34_init)*100, digits=1))%", 20))
plt_j34_inc_roti = annotate!(plt_j34_inc_roti, [0], [0.019], text(L"(\textit{a})", 20))

mean43_init = mean(res_gpr.jmat[4, 3, 1:end-1]); mean43_sim = mean(res_sim_jmat_inc_roti[4, 3, :])
plt_j43_inc_roti = plot(v_days[1:end-1], res_gpr.jmat[4, 3, 1:end-1], label = "In-sample "*replace("(mean $(round(mean43_init, digits=3)))", "-" => L"-"), xlabel = "Time (days)", ylabel = "Pico ⬅ Nano", title = "Pico ⬅ Nano (+20% Rotifers)", titlefontsize = 14, legendfontsize=10)
plt_j43_inc_roti = plot!(plt_j43_inc_roti, v_days[1:end-1], res_sim_jmat_inc_roti[4,3,:], label = replace("+20% Rotifers (mean $(round(mean43_sim, digits=3)))", "-" => L"-"))
plt_j43_inc_roti = annotate!(plt_j43_inc_roti, [500], [-0.2], text("ΔIS: $(round(-((mean43_sim-mean43_init)/mean43_init)*100, digits=1))%", 20))
plt_j43_inc_roti = annotate!(plt_j43_inc_roti, [0], [0.024], text(L"(\textit{b})", 20))

plot(plt_j34_inc_roti, plt_j43_inc_roti, layout=(2,1), size=(1500, 1500), legendfontsize=17, labelfontsize=13, tickfontsize=12, legend=:bottomright, linewidths=2, margin=1.5Plots.cm, top_margin = 0.65Plots.cm, bottom_margin = 1.3Plots.cm)
savefig("fig/fig5_sim_compet.png")



