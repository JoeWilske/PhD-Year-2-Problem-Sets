### Industrial Organization 1: Problem Set 2: Problem 1
### Topic: Entry
### Author: Joe Wilske
### Date: October 24, 2023

# This is the first version of my response to PS2 Problem 1, using a gridsearch instead of 
# an optimization algorithm.

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
        u = rand(Normal(μ, σ), F*M)
        strat_matrix = determine_all_strategies(data, u)
        s = vcat(s, [strat_matrix])
    end

    return s
end

# Compares simulated strategies to true strategies. Builds a vector that contains 
# information on how many markets in the true strategies matched the simulated 
# strategies. Returns this vector, frequenencies.
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
    return frequencies'
end

# Builds a matrix Θ, with height equal to the number of (μ, σ) pairs tested,
# and width equal to the number of markets. Each element of Θ is the number of
# times the (μ, σ) pair was able to match the corresponding number of markets' strategies.
# Displays number of test pairs remaining. Returns Θ.
function search_parameters(test_vector, data, Strats, simulations_per_guess)
    Θ = zeros(M)'
    for test_number in eachindex(test_vector[:,1])
        display(length(test_vector[:,1]) - test_number + 1)
        (μ, σ) = test_vector[test_number, :]
        s = run_simulation(μ, σ, data, simulations_per_guess)
        frequencies = compare_simulation_to_data(s, Strats)
        Θ = vcat(Θ, frequencies)
    end

    return Θ[2:end, :]
end

# Checks to see which row in Θ has the highest product, which is proportional to
# the probability of matching a simulated market to the data. Returns the number
# of the row with the highest probability.
function interpret_simulations(Θ, simulations_per_guess)
    freq_vec = []
    for rownum in eachindex(Θ[:, 1])
        freq_vec = vcat(freq_vec, prod(Θ[rownum, :] ./ (simulations_per_guess / 10)))
    end

    return findmax(freq_vec)[2]
end

# Given specifications on which (μ, σ) pairs to test, creates the vector (technically a matrix)
# of pairs to test. Returns test_vector.
function build_test_vector(decimal_granularity, μ_lower_bound, μ_upper_bound, σ_lower_bound, σ_upper_bound)
    test_vector = [0,0]'
    for μ in (μ_lower_bound * 10^(decimal_granularity) : μ_upper_bound * 10^(decimal_granularity)) / (10^decimal_granularity)
        for σ in (σ_lower_bound * 10^(decimal_granularity) : σ_upper_bound * 10^(decimal_granularity)) / (10^decimal_granularity)
            test_vector = vcat(test_vector, [μ, σ]')
        end
    end
    return test_vector[2:end, :]
end

# Reads the results of interpret_simulations() and pulls out the corresponding (μ, σ) from test_vector.
# Returns the estimates of μ and σ as a vector.
function spit_out_estimates(Θ, simulations_per_guess, test_vector)

    test_vector = build_test_vector(decimal_granularity,
                                    μ_lower_bound, μ_upper_bound, σ_lower_bound, σ_upper_bound)
    winners = interpret_simulations(Θ, simulations_per_guess)

    μ_hat = []
    σ_hat = []

    for index in winners
        μ_hat = vcat(μ_hat, test_vector[index, 1])
        σ_hat = vcat(σ_hat, test_vector[index, 2])
    end

    return (mean(μ_hat), mean(σ_hat))
end


### --- Control Panel: -------------------------------------------------------

# Controls accuracy.
# Multiplying simulations_per_guess by 10 will multiply time-to-compute by 10.
simulations_per_guess = 100

# Controls precision (number of decimal places).
# Adding 1 to decimal_granularity will multiply time-to-compute by 10.
decimal_granularity = 1

μ_lower_bound = 1
μ_upper_bound = 3

σ_lower_bound = 0
σ_upper_bound = 2

# Create the test vector.
test_vector = build_test_vector(decimal_granularity, μ_lower_bound, μ_upper_bound, σ_lower_bound, σ_upper_bound)

# Create Θ
Θ = search_parameters(test_vector, data, Strats, simulations_per_guess)

# Recover results from Θ.
(μ_hat, σ_hat) = spit_out_estimates(Θ, simulations_per_guess, test_vector)                                   

# Print.
println("μ_hat = ", μ_hat)
println("σ_hat = ", σ_hat)
