# Panel Econometrics: Problem Set 4, Problem 1
# Author: Joe Wilske
# Date: March 28, 2024

using Statistics, Random, Distributions, Optim, LinearAlgebra

#################################################################################################
### --- Problem 1 --- ###########################################################################
#################################################################################################


# Define model
function model(α0, β0, σ0)
    x_i = rand(Normal(0, 1))
    ϵ_i = rand(Normal(0, σ0))

    y_i = α0 + x_i*β0 + ϵ_i

    (y_i > 0.0) ? (d_i = 1) : (d_i = 0.0)

    (y_i > 0.0) ? [x_i, y_i, d_i]' : [x_i, 0.0, d_i]'
end

# Define data generating function
generate_data(α0, β0, σ0, obs) = vcat(model.(fill(α0, obs), β0, σ0)...)

# Set true values of variables
(α0, β0, σ0) = (0.5, 1.0, 1.0)
θ0 = [α0, β0, σ0]


### --- Part (A) --- ############################################################################

### --- MAXIMUM SCORE --- #######################################################################


# Define the score function
function s(β, data)

    α = sqrt( 1 - β^2 )

    (x, d) = (data[:,1], data[:,3])

    indicator = ifelse.( α .+ x.*β .> 0, 1, 0 )
    
    return (
        mean( (d .- indicator).^2 )
    )

end

# Define score estimator 
function score_estimator(data, precision)

    value = 1.0
    β_hat = 0.0

    for β in 0:(10.0^(-precision)):1
        
        min = s(β, data)

        if min < value

            value = min
            β_hat = β

        end

    end

    α_hat = round(sqrt( 1 - β_hat^2 ), digits = precision)

    return [α_hat, β_hat]'

end

# Takes number of simluation estimations to perform, observations per simulation,
# precision, and true parameter values in θ0.
function more_score_estimation(attempts, obs, precision, θ0)

    (α0, β0, σ0) = (θ0[1], θ0[2], θ0[3])

    vcat(
        score_estimator.( 
            generate_data.(α0, β0, σ0, fill(obs, attempts)),
            precision
        )...
    )

end

# Set how many decimal places over which to grid search.
precision = 3

# Do score estimation 400 times, first with 100 observations each, then with 200 obs, then with 400 obs.
s_estimates_100 = more_score_estimation(400, 100, precision, θ0)
s_estimates_200 = more_score_estimation(400, 200, precision, θ0)
s_estimates_400 = more_score_estimation(400, 400, precision, θ0)

# Score estimator can only give estimates to scale, so we'll compare the estimates to normalized α0 and β0.
(α_norm, β_norm) = ( normalize([α0, β0])[1], normalize([α0, β0])[2] )

# Define mean bias function. Returns in order of (α, β)
mean_bias(estimates) = ( mean(estimates[:, 1] .- α_norm), mean(estimates[:, 2] .- β_norm) )

# Define mean squared error function. Returns in order of (α, β)
MSE(estimates) = ( mean((estimates[:, 1] .- α_norm) .^ 2), mean((estimates[:, 2] .- β_norm) .^ 2) )

# Find mean bias and MSE for α and β estimates for 100, 200, and 400 observation data.
bias_100 = mean_bias(s_estimates_100)
mse_100 = MSE(s_estimates_100)
bias_200 = mean_bias(s_estimates_200)
mse_200 = MSE(s_estimates_200)
bias_400 = mean_bias(s_estimates_400)
mse_400 = MSE(s_estimates_400)

# Print all
println(
    "\nMAX SCORE, 100 OBSERVATIONS", "\nα mean bias:  ", round(bias_100[1], digits = precision),
    "\nβ mean bias:  ", round(bias_100[2], digits = precision),
    "\nα MSE:        ", round(mse_100[1], digits = precision), 
    "\nβ MSE:        ", round(mse_100[2], digits = precision), "\n",
    "\nMAX SCORE, 200 OBSERVATIONS", "\nα mean bias:  ", round(bias_200[1], digits = precision), 
    "\nβ mean bias:  ", round(bias_200[2], digits = precision),
    "\nα MSE:        ", round(mse_200[1], digits = precision), 
    "\nβ MSE:        ", round(mse_200[2], digits = precision), "\n", 
    "\nMAX SCORE, 400 OBSERVATIONS", "\nα mean bias:  ", round(bias_400[1], digits = precision), 
    "\nβ mean bias:  ", round(bias_400[2], digits = precision),
    "\nα MSE:        ", round(mse_400[1], digits = precision), 
    "\nβ MSE:        ", round(mse_400[2], digits = precision)
)


### --- MAXIMUM RANK --- ########################################################################


# Define the rank function
function r(β, data)

    (x, d) = (data[:,1], data[:,3])

    rank = 0.0

    for i in 1:length(d) 
        
        for j in 1:length(d)

            if i ≠ j

                (d[i] > d[j]) ? (y_ind = 1.0) : (y_ind = 0.0)

                (x[i]*β > x[j]*β) ? (x_ind = 1.0) : (x_ind = 0.0)

                rank = rank + (y_ind * x_ind)

            end
        end
    end
    return (rank)
end


# Define rank estimator 
function rank_estimator(data)

    value = 1.0
    β_hat = 0.0

    for β in -1:1
        
        max = r(β, data)

        if max > value

            value = max
            β_hat = β

        end
    end
    return β_hat
end

# Takes number of simluation estimations to perform, observations per simulation,
# precision, and true parameter values in θ0.
function more_rank_estimation(attempts, obs, θ0)

    (α0, β0, σ0) = (θ0[1], θ0[2], θ0[3])

    rank_estimator.(
        generate_data.(α0, β0, σ0, fill(obs, attempts))
    )

end

# Do rank estimation 400 times, first with 100 observations each, then with 200 obs, then with 400 obs.
r_estimates_100 = more_rank_estimation(400, 100, θ0)
r_estimates_200 = more_rank_estimation(400, 200, θ0)
r_estimates_400 = more_rank_estimation(400, 400, θ0)

# Rank estimation can only get you the sign of the variable β at best, so there's no point in comparing 
# bias and MSE. The estimator works very well at finding the sign, and we'll leave it at that.
println("\nMAX RANK, 100 OBSERVATIONS", "\nmean β estimate:  ", mean(r_estimates_100), "\n",
"\nMAX RANK, 200 OBSERVATIONS", "\nmean β estimate:  ", mean(r_estimates_200), "\n",
"\nMAX RANK, 400 OBSERVATIONS",  "\nmean β estimate:  ", mean(r_estimates_400))


### --- Part (B) --- ############################################################################

### --- CENSORED LEAST ABSOLUTE DEVIATIONS (CLAD) --- ###########################################


# Define the CLAD function
function clad(θ, data)

    mean(
        abs.(
            data[:,2] .- max.(θ[1] .+ data[:,1] .* θ[2], 0)
        )
    )

end

# Define CLAD minimizer
function clad_estimator(θ_initial, data)

    Optim.minimizer(
        optimize(θ -> clad(θ, data), θ_initial)
    )'

end

# Takes number of simluation estimations to perform, observations per simulation,
# θ_initial, and true parameter values in θ0.
function more_clad_estimation(attempts, obs, θ_initial, θ0)

    (α0, β0, σ0) = (θ0[1], θ0[2], θ0[3])

    vcat(
        clad_estimator.(
            Ref(θ_initial),
            generate_data.(α0, β0, σ0, fill(obs, attempts))
        )...
    )

end

# Set initial guess for α, β
θ_initial = [0.75, 0.75]

# Do CLAD estimation 400 times, first with 100 observations each, then with 200 obs, then with 400 obs.
c_estimates_100 = more_clad_estimation(400, 100, θ_initial, θ0)
c_estimates_200 = more_clad_estimation(400, 200, θ_initial, θ0)
c_estimates_400 = more_clad_estimation(400, 400, θ_initial, θ0)

# Find mean bias and MSE for α and β estimates for 100, 200, and 400 observation data.
bias_100 = mean_bias(c_estimates_100)
mse_100 = MSE(c_estimates_100)
bias_200 = mean_bias(c_estimates_200)
mse_200 = MSE(c_estimates_200)
bias_400 = mean_bias(c_estimates_400)
mse_400 = MSE(c_estimates_400)

# Print all
println(
    "\nCLAD, 100 OBSERVATIONS", "\nα mean bias:  ", round(bias_100[1], digits = 4),
    "\nβ mean bias:  ", round(bias_100[2], digits = 4),
    "\nα MSE:        ", round(mse_100[1], digits = 4), 
    "\nβ MSE:        ", round(mse_100[2], digits = 4), "\n",
    "\nCLAD, 200 OBSERVATIONS", "\nα mean bias:  ", round(bias_200[1], digits = 4), 
    "\nβ mean bias:  ", round(bias_200[2], digits = 4),
    "\nα MSE:        ", round(mse_200[1], digits = 4), 
    "\nβ MSE:        ", round(mse_200[2], digits = 4), "\n", 
    "\nCLAD, 400 OBSERVATIONS", "\nα mean bias:  ", round(bias_400[1], digits = 4), 
    "\nβ mean bias:  ", round(bias_400[2], digits = 4),
    "\nα MSE:        ", round(mse_400[1], digits = 4), 
    "\nβ MSE:        ", round(mse_400[2], digits = 4)
)


### --- PAIRWISE DIFFERENCE (PWD) --- ###########################################################


# Define transformation function
e(data, θ, i, j) = max( (data[i,2] - data[i,1]*θ[2] - θ[1]), (- data[j,1]*θ[2] - θ[1]) )

# Define PWD value function
function pwd(θ, data)

    value = 0.0

    for i in 1:length(data[:,1])

        for j in 1:length(data[:,1])

            if i ≠ j

                value = value + e(data, θ, i, j) - e(data, θ, j, i)

            end
        end
    end

    return (
        value / ( length(data[:,1]) * (1 - length(data[:,1])) )
    )
end

# Define PWD minimizer
function pwd_estimator(θ_initial, data)

    Optim.minimizer(
        optimize(θ -> pwd(θ, data), θ_initial)
    )'

end

# Takes number of simluation estimations to perform, observations per simulation,
# θ_initial, and true parameter values in θ0.
function more_pwd_estimation(attempts, obs, θ_initial, θ0)

    (α0, β0, σ0) = (θ0[1], θ0[2], θ0[3])

    vcat(
        pwd_estimator.(
            Ref(θ_initial),
            generate_data.(α0, β0, σ0, fill(obs, attempts))
        )...
    )

end

# Set initial guess for α, β
θ_initial = [0.75, 0.75]

# Do PWD estimation 400 times, first with 100 observations each, then with 200 obs, then with 400 obs.
p_estimates_100 = more_pwd_estimation(400, 100, θ_initial, θ0)
p_estimates_200 = more_pwd_estimation(400, 200, θ_initial, θ0)
p_estimates_400 = more_pwd_estimation(400, 400, θ_initial, θ0)

# Find mean bias and MSE for α and β estimates for 100, 200, and 400 observation data.
bias_100 = mean_bias(p_estimates_100)
mse_100 = MSE(p_estimates_100)
bias_200 = mean_bias(p_estimates_200)
mse_200 = MSE(p_estimates_200)
bias_400 = mean_bias(p_estimates_400)
mse_400 = MSE(p_estimates_400)

# Print all
println(
    "\nPWD, 100 OBSERVATIONS", "\nα mean bias:  ", round(bias_100[1], digits = 4),
    "\nβ mean bias:  ", round(bias_100[2], digits = 4),
    "\nα MSE:        ", round(mse_100[1], digits = 4), 
    "\nβ MSE:        ", round(mse_100[2], digits = 4), "\n",
    "\nPWD, 200 OBSERVATIONS", "\nα mean bias:  ", round(bias_200[1], digits = 4), 
    "\nβ mean bias:  ", round(bias_200[2], digits = 4),
    "\nα MSE:        ", round(mse_200[1], digits = 4), 
    "\nβ MSE:        ", round(mse_200[2], digits = 4), "\n", 
    "\nPWD, 400 OBSERVATIONS", "\nα mean bias:  ", round(bias_400[1], digits = 4), 
    "\nβ mean bias:  ", round(bias_400[2], digits = 4),
    "\nα MSE:        ", round(mse_400[1], digits = 4), 
    "\nβ MSE:        ", round(mse_400[2], digits = 4)
)

