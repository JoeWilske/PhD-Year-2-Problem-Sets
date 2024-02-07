### Industrial Organization 1: Problem Set 2: Problem 1
### Topic: Entry
### Author: Joe Wilske
### Date: October 24, 2023

# This is the final version of my response to PS1 Problem 1, using Optim.jl

using CSV, DataFrames, Statistics, Random, Distributions, Optim

# Load dataset
path = "C:\\Gitdir\\Projects\\PhD-Year-2-Problem-Sets\\IO1_PS2\\entryData.csv"
entry_ref = CSV.read(path, DataFrame)
entry = copy(entry_ref)

# Define model
ϕ(Z, α, u) = Z*α + u
π(X, β, δ, N, Z, α, u) = X*β .- δ*log(N) .- ϕ(Z, α, u)
F = 3
M = 100
(α, β, δ) = (1, 1, 1)

X = entry.MarketChar
Z = hcat(entry.CrossChar1, entry.CrossChar2, entry.CrossChar3)
Strats = hcat(entry.Strat1, entry.Strat2, entry.Strat3)

data = (X, β, δ, Z, α)

# Given data, finds the number of firms N that can enter a market m 
# for a particular u_m. Returns N.
function find_N(data, u_m, m)
    (X, β, δ, Z, α) = data

    result_1 = 0; result_2 = 0; result_3 = 0

    for f in 1:F
        if π(X, β, δ, 1, Z[m, f], α, u_m[f])[m] > 0
            result_1 += 1
        end
        if π(X, β, δ, 2, Z[m, f], α, u_m[f])[m] > 0
            result_2 += 1
        end
        if π(X, β, δ, 3, Z[m, f], α, u_m[f])[m] > 0
            result_3 += 1
        end
    end

    if result_3 == 3
        N = 3
    elseif result_2 >= 2
        N = 2
    elseif result_1 >= 1
        N = 1
    else
        N = 0
    end

    return N
end

# Given data and N, finds which firms enter a market m for a particular u, u_m.
# Returns the market strategy profile s_m.
function determine_market_strategies(data, u_m, m, N)
    (X, β, δ, Z, α) = data

    if N == 0
        s_m = [0,0,0]'
    elseif N == 1
        profitabilities = [π(X, β, δ, 1, Z[m, 1], α, u_m[1])[m],
                           π(X, β, δ, 1, Z[m, 2], α, u_m[2])[m],
                           π(X, β, δ, 1, Z[m, 3], α, u_m[3])[m]]
        first_firm = findmax(profitabilities)[2]
        if first_firm == 1
            s_m = [1, 0, 0]'
        elseif first_firm == 2
            s_m = [0, 1, 0]'
        else
            s_m = [0, 0, 1]'
        end
    elseif N == 2
        profitabilities = [π(X, β, δ, 2, Z[m, 1], α, u_m[1])[m],
                           π(X, β, δ, 2, Z[m, 2], α, u_m[2])[m],
                           π(X, β, δ, 2, Z[m, 3], α, u_m[3])[m]]
        last_firm = findmin(profitabilities)[2]
        if last_firm == 1
            s_m = [0, 1, 1]'
        elseif last_firm == 2
            s_m = [1, 0, 1]'
        else
            s_m = [1, 1, 0]'
        end
    else
        s_m = [1,1,1]'
    end
    return s_m
end

# Given data, finds the strategy profile of all firms in all markets for a partiuclar u.
# Returns the strategy profile s.
function determine_all_strategies(data, u)
    m = 1
    u_m = [u[m], u[m + M], u[m + 2*M]]
    N = find_N(data, u_m, m)
    s = determine_market_strategies(data, u_m, m, N)

    for m in 2:M
        u_m = [u[m], u[m + M], u[m + 2*M]]
        N = find_N(data, u_m, m)
        s_m = determine_market_strategies(data, u_m, m, N)
        s = vcat(s, s_m)
    end

    return s
end

# Given data, runs simulations_per_guess number of simulations for a 
# partiuclar μ, σ guess. Returns a vector of all strategy profiles generated
# in the simulations, s.
function run_simulation(μ, σ, data, simulations_per_guess)
    s = []
    for draw in 1:simulations_per_guess
        if σ < 0
            σ = 0
        end
        u = rand(Normal(μ, σ), F*M)
        strat_matrix = determine_all_strategies(data, u)
        s = vcat(s, [strat_matrix])
    end

    return s
end

# Compares simulated strategies to true strategies. Builds a vector that contains 
# information on how often each simulated market matched the true market, called frequencies.
# Returns frequencies.
function compare_simulation_to_data(s, Strats)
    frequencies = zeros(M)
    for simulation in s
        for m in 1:M
            if simulation[m, :] == Strats[m, :]
                frequencies[m] += 1
            end
        end
    end
    for i in eachindex(frequencies)
        if frequencies[i] == 0
            frequencies[i] = 1
        end
    end
    return frequencies
end

# Objective function. Runs simulations and generates frequency. 
# Returns negative product of frequency vector, which is proportional 
# to the probability of the simulated outcome matching the true outcome.
function G(θ)
    (μ, σ) = θ
    display("Attempt")
    s = run_simulation(μ, σ, data, simulations_per_guess)
    frequencies = compare_simulation_to_data(s, Strats)
    
    return -prod(frequencies ./ (simulations_per_guess / 10))
end

# Switches between objective functions G or H depending on use_G argument.
# Switches between optimization types depending on use_boundaries argument.
function minimization_type(use_boundaries, lower, upper, θ_initial)
    if use_boundaries == true
        res = optimize(G, lower, upper, θ_initial)
    else
        res = optimize(G, θ_initial)
    end
    return res
end


### --- Control Panel: -----------------------------------------------------------------

# Set use_boundaries to true if you want to use "box minimization" with boundaries. 
# Set use_boundaries to false if you want to use Nelder-Mead minimization without boundaries.
use_boundaries = false

# Lower and upper bounds for μ and σ.
lower = [-10, 0]
upper = [10, 10]

# Initial guesses for μ and σ.
(μ_hat, σ_hat) = (1.5, 1.5)

# How many simulations per (μ, σ) guess.
simulations_per_guess = 1000

# This is where the magic happens.
θ_initial = copy([μ_hat, σ_hat])
res = minimization_type(use_boundaries, lower, upper, θ_initial)
θ_hat = Optim.minimizer(res)

# Print.
println("μ_hat = ", θ_hat[1])
if θ_hat[2] < 0; θ_hat[2] = 0; end
println("σ_hat = ", θ_hat[2])
