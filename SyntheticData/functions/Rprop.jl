### Functions for GPR model fitting with Resillient Backpropagation

using StatsBase
using LinearAlgebra
using Plots
using Distances
import DataFrames.Not

# SE kernel function for a single data point
function sekernel(eta, theta, x1, x2, dmax)
    eta * exp(-sum(theta .* (((x1 - x2) ./ dmax) .^ 2)))
end

# Covariance matrix with SE kernel
function sekernel_dist(eta, dist)
    eta * exp(-dist^2)
end

# Distance matrix with ARD
function distm_ard(theta, X)
    theta_ard = sqrt.(theta)
    distm = pairwise(Euclidean(), X .* theta_ard', dims=1)
    distm
end

# covariance matrix with SE kernel
function covmat(eta, distm_scaled, sigma)
    n = size(distm_scaled, 1) # ncol
    K = zeros(n, n)

    @. K = sekernel_dist(eta, distm_scaled)
    for i in 1:n
        K[i, i] = eta + sigma
    end

    return K
end

# Covariance vector between target point and data point to calculate posterior predictive distribution
function cov_pred_onestep(pars, Xtr, xtes, dmax)
    E = size(Xtr, 2)
    theta = pars[1:E]
    eta = pars[E+1]
    theta_ard = sqrt.(theta)
    dist = sqrt.(sum(((Xtr .- xtes') .* theta_ard') .^ 2, dims=2)) / dmax

    # covariance vector bewtween target point and training data
    ks = sekernel_dist.(eta, dist)

    #need slicing since ks is a N x 1 matrix
    return ks[:, 1]
end

# Log priors in log-likelihood
function logpriors(theta, eta, sigma, maxvar, minvar=0.0001, par_beta=[2, 2])
    #half-normal(0,1^2) 
    lprior_theta = -0.5 * sum(theta .^ 2)
    dlprior_theta = -theta

    #default prior means of sigma and eta are 0.5*sigma_Y
    a_beta = par_beta[1]
    b_beta = par_beta[2]
    lprior_sigma = (a_beta - 1) * log((sigma - minvar) / (maxvar - minvar)) + (b_beta - 1) * log(1 - (sigma - minvar) / (maxvar - minvar))
    lprior_eta = (a_beta - 1) * log((eta - minvar) / (maxvar - minvar)) + (b_beta - 1) * log(1 - (eta - minvar) / (maxvar - minvar))

    #derivatives of priors
    dlprior_sigma = (a_beta - 1) / (sigma - minvar) - (b_beta - 1) / (maxvar - sigma)
    dlprior_eta = (a_beta - 1) / (eta - minvar) - (b_beta - 1) / (maxvar - eta)

    return ([lprior_theta, lprior_eta, lprior_sigma], vcat(dlprior_theta, dlprior_eta, dlprior_sigma))
end

# Logit transformation of the variance parameters (eta and sigma)
function logit(par_bound)
    log(par_bound / (1 - par_bound))
end

# setting boundaries for optimised parameters 
function setbound_par(par, parmax, parmin)
    (par - parmin) / (parmax - parmin)
end

# rescale parameters into the original scale
function rescale_etasig(par_scaled, parmax, parmin)
    parmin + (parmax - parmin) / (1 + exp(-par_scaled))
end

# negative log-likelihood function
function nllgp(par, X, Y, ldetK, invcholK_L, maxvar, minvar)
    n = size(X, 2)
    theta = par[1:n]
    eta = par[n+1]
    sigma = par[end]

    invL_Y = invcholK_L * Y
    lpriors = logpriors(theta, eta, sigma, maxvar, minvar, [2, 2])
    nlprior = -sum(lpriors[1])

    negllgp = 0.5 * invL_Y' * invL_Y + 0.5 * ldetK + nlprior

    return negllgp
end

# gradients of nll with SE kernel
function gradnllgp(par, distm, X, Y, invcholK_L, invcholK_U, etasigmax, etasigmin, dmax, scale=true)
    n = size(X, 2) #embedding dimension
    theta = par[1:n] # inverse length-scale parameters
    eta = par[n+1]
    sigma = par[end]

    invK_chol = invcholK_U * invcholK_L
    invK_Y = invcholK_U * invcholK_L * Y
    datasize = size(X, 1)
    lpriors = logpriors(theta, eta, sigma, etasigmax, etasigmin, [2, 2])
    dnlpriors = -lpriors[2]

    grad_nll = zeros(length(par))

    dK_theta = zeros(datasize, datasize)
    dK_eta = Matrix{Float64}(I, datasize, datasize)

    # gradients w.r.t each inverse length-scale parameters (theta)
    for i_theta in 1:n
        dXiXj = zeros(datasize, datasize)
        for i in 1:(datasize-1)
            dXiXj[i, (i+1):end] .= (X[i, i_theta] .- X[(i+1):datasize, i_theta]) / dmax
            dXiXj[(i+1):end, i] = view(dXiXj, i, (i+1):datasize)
        end

        @. dK_theta = -(dXiXj^2) * sekernel_dist(eta, distm)

        grad_nll[i_theta] = 0.5 * tr(invK_chol * dK_theta) - 0.5 * (invK_Y)' * dK_theta * invK_Y + dnlpriors[i_theta]
    end

    # derivatives of covariances w.r.t eta
    @. dK_eta = sekernel_dist(1.0, distm)

    # gradients w.r.t eta and sigma
    grad_nll[n+1] = 0.5 * tr(invK_chol * dK_eta) - 0.5 * (invK_Y)' * dK_eta * invK_Y + dnlpriors[n+1]
    grad_nll[end] = 0.5 * tr(invK_chol) - 0.5 * (invK_Y)' * invK_Y + dnlpriors[end]

    #scale gradients for Rprop optimization
    scale_grad = vcat(theta, (eta - etasigmin) * (1 - (eta - etasigmin) / (etasigmax - etasigmin)), (sigma - etasigmin) * (1 - (sigma - etasigmin) / (etasigmax - etasigmin)))

    ifelse(scale, grad_nll .* scale_grad, grad_nll)
end

# Main function for Rprop optimisation
function Rprop(X, Y; maxiter=100, scale_dist=true, print_res=false)
    npars = size(X, 2) + 2

    #optimisation parameters for Rprop
    delta0 = 0.1 * ones(npars)
    deltamin = repeat([0.000001], npars)
    deltamax = repeat([50.0], npars)
    eta_minus = -0.4
    eta_plus = 0.1

    #initial values
    pars = repeat([0.01], npars)
    iter = 0
    delta = delta0
    df = 1.0

    etasig_max = ifelse(var(Y) < 1.0, 1.0, var(Y))
    etasig_min = 0.0001

    theta = pars[1:(npars-2)]
    eta = pars[npars-1]
    sigma = pars[npars]

    distm_unscaled = pairwise(Euclidean(), X, dims=1)
    dmax = ifelse(scale_dist, maximum(vec(distm_unscaled)), 1.0)
    distm = distm_ard(theta, X)
    distm_scaled = distm / dmax

    K = covmat(eta, distm_scaled, sigma)
    cholK = cholesky(K)
    ldetK = 2.0 * sum(log.(diag(cholK.L)))
    invcholK_L = inv(cholK.L)
    invcholK_U = inv(cholK.L')

    nll = [nllgp(pars, X, Y, ldetK, invcholK_L, etasig_max, etasig_min)]
    grad = gradnllgp(pars, distm_scaled, X, Y, invcholK_L, invcholK_U, etasig_max, etasig_min, dmax)
    tau = sqrt(grad' * grad)

    #parameter transformation
    ##eta and sigma <= var(Y)
    pars[(npars-1):end] = setbound_par.(pars[(npars-1):end], etasig_max, etasig_min)
    pars_scaled = vcat(log.(pars[1:(npars-2)]), logit.(pars[(npars-1):end]))

    while iter < maxiter && df > 0.000001 && tau > 0.0001
        pars_scaled += -sign.(grad) .* delta

        #rescale parameters into the original scale for nll and grad calculation
        pars_resc = vcat(exp.(pars_scaled[1:(npars-2)]), rescale_etasig.(pars_scaled[(npars-1):end], etasig_max, etasig_min))
        theta, eta, sigma = pars_resc[1:(npars-2)], pars_resc[npars-1], pars_resc[npars]

        distm = distm_ard(theta, X)

        distm_scaled = distm / dmax

        K = covmat(eta, distm_scaled, sigma)
        cholK = cholesky(Symmetric(K))
        ldetK = 2.0 * sum(log.(diag(cholK.L)))
        invcholK_L = inv(cholK.L)
        invcholK_U = inv(cholK.L')
        nll_new = nllgp(pars_resc, X, Y, ldetK, invcholK_L, etasig_max, etasig_min)

        grad_new = gradnllgp(pars_resc, distm_scaled, X, Y, invcholK_L, invcholK_U, etasig_max, etasig_min, dmax)

        #convergence
        df = abs(nll_new / nll[iter+1] - 1.0)
        tau = sqrt(grad_new' * grad_new)

        gradsign = grad .* grad_new
        delta_new = delta .* (1.0 .+ eta_plus * (gradsign .> 0.0) .+ eta_minus * (gradsign .< 0.0))
        delta = ifelse.(delta_new .> deltamax, deltamax, ifelse.(delta_new .< deltamin, deltamin, delta_new))

        grad = grad_new
        nll = vcat(nll, nll_new)
        iter += 1
        pars = pars_resc
    end

    if print_res
        print("Iteration: ", iter, " ll_ratio: ", df, " |grad|: ", tau)
    end

    return (pars=pars, nll=nll, invL=invcholK_L, dmax=dmax, distm_scaled=distm_scaled, K=K)
end



