# Panel Econometrics: Problem Set 3
# Author: Joe Wilske
# Date: February 14, 2024

using Statistics, Random, Distributions, Optim, LinearAlgebra


##########################################################################################
### --- Problem 1 --- ####################################################################
##########################################################################################


# Define model
function model(α0, β0)
    x_i, ϵ_i = rand(Normal(0, 1), 2)

    y_i = α0 + x_i*β0 + ϵ_i

    (y_i > 0.0) ? (return [x_i, y_i]') : (return [x_i, 0.0]')
end

# Generate data
function generate_data(α0, β0, obs)
    α0 = fill(α0, obs)
    return vcat(broadcast(model, α0, β0)...)
end

# Generate data with α0 = 0.5, β0 = 1.0, and 100 0bservations
data = generate_data(0.5, 1.0, 100)


### --- Part (A) --- #######################################################################


# Define log-PDF
function f(α, β, x, y)
    if y > 0.0
        return log(pdf(Normal(0, 1), y - α - x*β))
    else
        return log(cdf(Normal(0, 1), - α - x*β))
    end
end

# Define negative log-likelihood function
function G(θ, data)
    log_likes = broadcast(f, θ[1], θ[2], data[:, 1], data[:, 2])
    return -sum(log_likes)
end

# Define minimizer function
function minimize(α_initial, β_initial, data)

    θ_initial = [α_initial, β_initial]

    # Minimize the negative log-likelihood function.
    res = optimize(θ -> G(θ, data), θ_initial)
    return Optim.minimizer(res)'
end

# Find minimum of negative log-likelihood
θ_hat = minimize(0.0, 0.0, data)

# Print
println("α estimate:  ", θ_hat[1], "\nβ estimate:  ", θ_hat[2])


### --- Part (B) --- #########################################################################


# Takes initial parameter guesses as inputs, as well as number of estimations to perform.
# Generates new data each time to avoid getting the same estimates every run.
function minimize_more(α_initial, β_initial, attempts, obs)

    α0 = fill(0.5, attempts)
    all_data = broadcast(generate_data, α0, 1.0, obs)

    return vcat(broadcast(minimize, α_initial, β_initial, all_data)...)
end

# Do MLE estimates 400 times, each with 100 observations.
estimates_100 = minimize_more(0.0, 0.0, 400, 100)

# Define mean bias function
function mean_bias(estimates)

    α_bias = mean(estimates[:, 1] .- 0.5)
    β_bias = mean(estimates[:, 2] .- 1.0)

    return α_bias, β_bias
end

# Define mean squared error function
function MSE(estimates)

    α_mse = mean((estimates[:, 1] .- 0.5) .^ 2)
    β_mse = mean((estimates[:, 2] .- 1.0) .^ 2)

    return α_mse, β_mse
end

# Find mean bias and MSE for α and β estimates for 100 observation data
bias_100 = mean_bias(estimates_100)
mse_100 = MSE(estimates_100)

# Print
println(
    "α mean-bias:  ", bias_100[1], "\nβ mean-bias:  ", bias_100[2],
    "\nα MSE:        ", mse_100[1], "\nβ MSE:        ", mse_100[2]
)


### --- Part (C) --- ##########################################################################


# Do MLE estimates 400 times, each with 200 observations, and then each with 400 obs.
estimates_200 = minimize_more(0.0, 0.0, 400, 200)
estimates_400 = minimize_more(0.0, 0.0, 400, 400)

# Find mean bias and MSE for α and β estimates for 200 observation data, and then 400 obs data.
bias_200 = mean_bias(estimates_200)
mse_200 = MSE(estimates_200)
bias_400 = mean_bias(estimates_400)
mse_400 = MSE(estimates_400)

# Print all
println(
    "\n100 OBSERVATIONS:", "\nα mean bias:  ", bias_100[1], "\nβ mean bias:  ", bias_100[2],
    "\nα MSE:        ", mse_100[1], "\nβ MSE:        ", mse_100[2], "\n",
    "\n200 OBSERVATIONS:", "\nα mean bias:  ", bias_200[1], "\nβ mean bias:  ", bias_200[2],
    "\nα MSE:        ", mse_200[1], "\nβ MSE:        ", mse_200[2], "\n", 
    "\n400 OBSERVATIONS:", "\nα mean bias:  ", bias_400[1], "\nβ mean bias:  ", bias_400[2],
    "\nα MSE:        ", mse_400[1], "\nβ MSE:        ", mse_400[2]
)

