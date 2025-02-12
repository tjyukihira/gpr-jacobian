import StatsPlots.boxplot
import StatsPlots.boxplot!

function plot_strengths(res_gpr, v_days; row = 1, column = 2, title = "R ⬅ G", CI = [0.05, 0.95], digits = 2)
    # showing 90% posterior CI on default
    percent = round(Int64, (CI[2]-CI[1])*100)
    int_low = quantile.(Normal.(res_gpr.jmat[row, column, 1:end-1], res_gpr.jmat_sd[row, column, 1:end-1]), CI[1]); int_up = quantile.(Normal.(res_gpr.jmat[row, column, 1:end-1], res_gpr.jmat_sd[row, column, 1:end-1]), CI[2])
    ylim = (minimum(int_low) - 0.05, maximum(int_up) + 0.05)
    mean_low = round(mean(int_low), digits = digits); mean_up = round(mean(int_up), digits = digits)
    p = plot(v_days, res_gpr.jmat[row, column, 1:end-1], ylim = ylim, ylabel = "IS", xlabel = "Time (days)", title = title, label = replace("Mean: $(round(mean(res_gpr.jmat[row, column, 1:end-1]), digits = digits)) $(percent)%CI[$mean_low, $mean_up]", "-" => L"-")) #?
    #int_low = res_gpr.jmat[row, column, 1:end-1] - 2*res_gpr.jmat_sd[row, column, 1:end-1]; int_up = res_gpr.jmat[row, column, 1:end-1] + 2*res_gpr.jmat_sd[row, column, 1:end-1]
    p = plot!(v_days, [res_gpr.jmat[row, column, 1:end-1] res_gpr.jmat[row, column, 1:end-1]], fillrange=[int_low int_up], fillalpha=0.2, c=:blue, label=false)

    return p
end
