# Panel Econometrics: Problem Set 5
# Author: Joe Wilske
# Date: April 2, 2024

using Statistics, Random, Distributions, Optim, LinearAlgebra, CSV, DataFrames, GLM


#################################################################################################
### --- Problem 1: Local Linear Estimator --- ###################################################
#################################################################################################


# Returns n × 2 matrix of x's and y's
function generate_data(n, γ0)

    x = rand(Uniform(-1, 1), n)
    ϵ = rand(Normal(0, 1), n)

    y = γ0 .+ x .+ ϵ

    return hcat(x, y)

end

# Takes (x, y) matrix as data and the central_x around which to estimate y.
# Uses Kernel bandwidth h = n^(-1/5)
# Returns value of the LLE objective function
function LLE_objective(θ, data, central_x)

    α, β = θ
    (x, y) = (data[:, 1], data[:, 2])
    h = length(data[:, 1])^(- 0.2)

    K = ifelse.(abs.(x .- central_x) .< h/2, 1, 0)

    return (
        sum( K.*(y .- α .- β.*(x .- central_x)).^2 ) / length(data[:, 1])
    )

end

# Given θ_initial = (α_initial, β_initial), returns the (α, β)
# that minimize the objective function as a vector
function local_linear_estimator(θ_initial, data, central_x)

    Optim.minimizer(
        optimize(θ -> LLE_objective(θ, data, central_x), θ_initial)
    )'

end

# Takes number of simulations to run, observations per simulation,
# and the central x around which to estimate. Returns a num_sims × 2 matrix
# of (α, β) estimates, where α = E[y | x = central_x]
function more_LLE(num_sims, n, central_x, γ0, θ_initial)

    simulated_data = generate_data.(fill(n, num_sims), γ0)

    estimates = vcat(
        local_linear_estimator.(Ref(θ_initial), simulated_data, central_x)...
    )

    return estimates

end

# Set true γ
γ0 = 0.5

# Set initial (α, β) for minimization
θ_initial = [1.0, 1.0]

# Estimate 401 simulations, with 400 observations, around x = 0
estimates_0 = more_LLE(401, 400, 0.0, γ0, θ_initial)
# Estimate 401 simulations, with 400 observations, around x = 1
estimates_1 = more_LLE(401, 400, 1.0, γ0, θ_initial)

# Define mean bias function for α estimates
mean_bias_LLE(estimates, central_x, γ0) = mean(estimates[:, 1] .- (central_x + γ0))
# Define root mean squared error (RMSE) function for α estimates
RMSE_LLE(estimates, central_x, γ0) = mean((estimates[:, 1] .- (central_x + γ0)) .^ 2) ^ 0.5

# Find mean-bias and RMSE of α estimates for central_x = 0 and central_x = 1
bias_0 = mean_bias_LLE(estimates_0, 0.0, γ0)
rmse_0 = RMSE_LLE(estimates_0, 0.0, γ0)
bias_1 = mean_bias_LLE(estimates_1, 1.0, γ0)
rmse_1 = RMSE_LLE(estimates_1, 1.0, γ0)

# Print all
println(
    "\n-----------------------",
    "\nLLE with x = 0",
    "\nmean bias:    ", round(bias_0, digits = 3),
    "\nRMSE:         ", round(rmse_0, digits = 3),
    "\n=======================",
    "\nLLE with x = 1",
    "\nmean bias:    ", round(bias_1, digits = 3),
    "\nRMSE:         ", round(rmse_1, digits = 3),
    "\n------------------------"
)


#################################################################################################
### --- Wooldridge, Problem 18.3 --- ############################################################
#################################################################################################


# Get path to directory
dir = @__DIR__

# Set path to CSV
csv_file = joinpath(dir, "jtrain2.csv")

# Read-in dataset
df = CSV.read(csv_file, DataFrame)

# Given data frame, returns a vector of individuals' calculated propensity score_estimator
# as estimated by a probit regression of train on some variables
function prop_score(df)

    data_mat = hcat(
        ones(length(df[:, :black])),
        df[:, :black],
        df[:, :hisp],
        df[:, :married],
        df[:, :nodegree],
        df[:, :re74].^2,
        df[:, :re75].^2
    )

    probit_model = glm(@formula(train ~ black + hisp + married + nodegree + re74^2 + re75^2), df, Binomial(), ProbitLink())

    probit_fit = coef(probit_model)
    
    return data_mat * probit_fit
end

# Execute above function. Propensity is a vector of individual propensity scores.
propensity = prop_score(df)

# Add propensity to data frame
df[!, :propensity] = propensity

# Run OLS regression 18.23 from Wooldridge
OLS_1823 = lm(@formula(unem78 ~ train + propensity), df)

# Add propensity deviation from mean to data frame
df[!, :propdeviation] = df[:, :train] .* (propensity .- mean(propensity))

# Run OLS regression 18.24 from Wooldridge
OLS_1824 = lm(@formula(unem78 ~ train + propensity + propdeviation), df)

# Run a third OLS regression, this time with no propensity score controls
OLS_3 = lm(@formula(unem78 ~ train + black + hisp + married + nodegree + re74^2 + re75^2), df)

