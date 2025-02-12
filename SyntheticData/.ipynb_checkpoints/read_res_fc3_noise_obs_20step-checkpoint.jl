vec_res_fc3 = Vector{Tuple}(undef, 30)
for i in 1:30
    vec_res_fc3[i] = load_gpresults_jld2("./results/fc3_noise_obs_20step_dist/0809_fc3_noise_obs_20step_model$(i)_se.jld2";npop=3, jacobian=true)
end