# Panel Econometrics: Problem Set 3
# Author: Joe Wilske
# Date: February 14, 2024

using Statistics, Random, Distributions, Optim, LinearAlgebra


##########################################################################################
### --- Problem 1 --- ####################################################################
##########################################################################################


# Define model
function model(α0, β0, σ0)
    x_i = rand(Normal(0, 1))
    ϵ_i = rand(Normal(0, σ0))

    y_i = α0 + x_i*β0 + ϵ_i

    y_i > 0.0  ?  [x_i, y_i]'  :  [x_i, 0.0]'
end

# Define data generating function
generate_data(α0, β0, σ0, obs) = vcat(broadcast(model, fill(α0, obs), β0, σ0)...)

# Generate data
α0, β0, σ0, obs = (0.5, 1.0, 1.0, 100)
data = generate_data(α0, β0, σ0, obs)


### --- Part (A) --- #######################################################################


# Define log-PDF
function f(α, β, σ, x, y)
    σ = max(σ, nextfloat(0.0))

    if y > 0.0
        log(pdf(Normal(0, σ), y - α - x*β))
    else
        log(cdf(Normal(0, σ), - α - x*β))
    end

end

# Define negative log-likelihood function
log_likelihood(θ, data) = -sum(broadcast(f, θ[1], θ[2], θ[3], data[:, 1], data[:, 2]))

# Define minimizer function
function minimize(α_initial, β_initial, σ_initial, data)

    θ_initial = [α_initial, β_initial, σ_initial]

    # Minimize the negative log-likelihood function.
    Optim.minimizer(
        optimize(θ -> log_likelihood(θ, data), θ_initial)
    )'
end

# Find minimum of negative log-likelihood
α_initial, β_initial, σ_initial = (0.0, 0.0, 2.0)
θ_hat = minimize(α_initial, β_initial, σ_initial, data)

# Print
println("α estimate:  ", θ_hat[1], "\nβ estimate:  ", θ_hat[2])


### --- Part (B) --- #########################################################################


# Takes initial parameter guesses as inputs, as well as number of estimations to perform.
# Generates new data each time to avoid getting the same estimates every run.
function minimize_more(α_initial, β_initial, σ_initial, attempts, obs)

    vcat(
        broadcast(
            minimize, α_initial, β_initial, σ_initial, 
            broadcast(generate_data, fill(0.5, attempts), 1.0, 1.0, obs)
        )...
    )
end

# Do MLE estimates 400 times, each with 100 observations.
attempts, obs = (400, 100)
estimates_100 = minimize_more(α_initial, β_initial, σ_initial, attempts, obs)

# Define mean bias function. Returns in order of (α, β)
mean_bias(estimates) = ( mean(estimates[:, 1] .- 0.5), mean(estimates[:, 2] .- 1.0) )

# Define mean squared error function. Returns in order of (α, β)
MSE(estimates) = ( mean((estimates[:, 1] .- 0.5) .^ 2), mean((estimates[:, 2] .- 1.0) .^ 2) )

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
obs = 200
estimates_200 = minimize_more(α_initial, β_initial, σ_initial, attempts, obs)
obs = 400
estimates_400 = minimize_more(α_initial, β_initial, σ_initial, attempts, obs)

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

