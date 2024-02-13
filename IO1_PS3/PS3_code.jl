### Industrial Organization 1: Problem Set 3
### Topic: Dynamics
### Author: Joe Wilske
### Date: December 21, 2023


using Statistics, Random, Distributions, Optim, LinearAlgebra


####################################################################################
### --- (3) Contraction Mapping --- ################################################
####################################################################################


# Profit function, not including ε shocks.
π(a, i, μ, R) = (i == 0) ? (a*μ) : (R)

# Log-sum formula of expected value. EV_0 is a 5-vector of expected values when i = 0,
# and EV_1 is a 5-vector of expected values when i = 1.
function log_sum(a, i, μ, R, EV_0, EV_1, β)

    (i == 0) ? (a = min(a + 1, 5)) : (a = 1)

    log( exp( π(a, 0, μ, R) + β*EV_0[a]) + exp( π(a, 1, μ, R) + β*EV_1[a]) )
end

# Present value function, taking all future expected values into account.
value_function(a, i, μ, R, EV_0, EV_1, β) = π(a, i, μ, R) + β*log_sum(a, i, μ, R, EV_0, EV_1, β)

# Performs contraction mapping. Returns a 5×2 matrix, where column 1 is
# the values for i = 0 and column 2 is the values for i = 1.
function contraction_map(μ, R, β)
    
    EV_0, EV_1, V_0, V_1 = [fill(1.0, 5) for _ in 1:4]
            
    α = 1:5

    while true

        EV_0 = log_sum.(α, 0, μ, R, Ref(EV_0), Ref(EV_1), β)
        EV_1 = log_sum.(α, 1, μ, R, Ref(EV_0), Ref(EV_1), β)

        prev_V_0 = deepcopy(V_0)
        prev_V_1 = deepcopy(V_1)

        V_0 = value_function.(α, 0, μ, R, Ref(EV_0), Ref(EV_1), β)
        V_1 = value_function.(α, 1, μ, R, Ref(EV_0), Ref(EV_1), β)
        
        cond_0 = abs.(V_0 - prev_V_0) .== 0.0
        cond_1 = abs.(V_1 - prev_V_1) .== 0.0
        
        if all(x -> x == 1, cond_0) && all(y -> y == 1, cond_1)
            break
        end

    end

    hcat(V_0, V_1)
end


### --- (3 a) --- ######################################################################


# Specify parameters
(μ, R, β) = (-1, -3, 0.9)

# Find values for all possible (a, i).
values = contraction_map(μ, R, β)

# Firm is indifferent between replacement or not if ε0 - ε1 = indif, below
indif = values[2, 2] - values[2, 1]

# Difference between two Gumbel(0,1) variables follows a
# Logistic(0,1) distribution, so P(firm replaces when a = 2) is
probability = cdf(Logistic(0, 1), indif)

# The value of a firm when a = 4, ε0 = 1, and ε1 = 1.5 is
val = max(values[4, 1] + 1, values[4, 2] + 1.5)


########################################################################################
### --- (4) Simulate Data --- ##########################################################
########################################################################################


# Given an a, takes ε draws and compares resultant firm values. 
# Returns the optimal i decision.
function i_from_a(a, μ, R, V, β)
    
    V_0 = V[a, 1] + rand(Gumbel(0, 1))
    V_1 = V[a, 2] + rand(Gumbel(0, 1))
    
    (V_0 > V_1)  ?  0  :  1
end

# Given μ, R, β, and number of observations T, generates a simulated dataset.
# Returns an T×2 matrix of a's and i's.
function generate_data(μ, R, β, T)
    a_vec = zeros(Int, T); a_vec[1] = 1
    i_vec = zeros(Int, T)

    V = contraction_map(μ, R, β)
    
    for index in 2:T

        i_vec[index - 1] = i_from_a(a_vec[index - 1], μ, R, V, β)

        if i_vec[index - 1] == 0

            if a_vec[index - 1] < 5
                a_vec[index] = a_vec[index - 1] + 1
            else 
                a_vec[index] = 5
            end

        else 
            a_vec[index] = 1
        end
    end

    i_vec[T] = i_from_a(a_vec[T], μ, R, V, β)

    hcat(a_vec, i_vec)
end

# Generate simulated data
data = generate_data(-1, -3, 0.9, 20000)


#######################################################################################
### --- (5) Estimate θ Using Rust's NFPA --- #############################################
######################################################################################


# Given (a, i) and firm values, returns the evaluated logit formula.
P(a, i, values) = exp(values[a, i + 1]) / (exp(values[a, 1]) + exp(values[a, 2]))

# Given target (a, i) and data (a_point, i_point),
# returns 1 if (a, i) match (a_point, i_point). Returns 0 otherwise.
count_combos(a, i, a_point, i_point) = (a_point == a) && (i_point == i) ? 1 : 0

# Given θ and data, returns the value of the joint negative log-likelihood function.
function loglikelihood_1(θ, data)
    (μ, R) = θ

    a_vec = data[:, 1]
    i_vec = data[:, 2]

    values = contraction_map(μ, R, 0.9)

    P_vec = zeros(10)
    Q_vec = zeros(10)

    for index in 1:5
        P_vec[index] = P(index, 0, values)
        Q_vec[index] = sum(count_combos.(index, 0, a_vec, i_vec))
    end

    for index in 6:10
        P_vec[index] = P(index - 5, 1, values)
        Q_vec[index] = sum(count_combos.(index - 5, 1, a_vec, i_vec))
    end

    -sum( log.(P_vec) .* Q_vec )
end

# Initial guesses for μ and R.
μ_hat, R_hat = (0.0, 0.0)
θ_initial = copy([μ_hat, R_hat])

# Minimize the negative log-likelihood function.
res = @time optimize(θ -> loglikelihood_1(θ, data), θ_initial)
θ_hat = Optim.minimizer(res)

# Print.
println("μ estimate: ", θ_hat[1], "\nR estimate: ", θ_hat[2])


###################################################################################
### --- (6) Estimate θ Using Hotz & Miller's CCP Approach --- #####################
###################################################################################


# Given data, returns 5×2 matrix of probalities:
# prob[a, i] = probability of choice i at age a
function get_replacement_probs(data)

    prob = zeros(5, 2)
    a_vec = data[:, 1]
    i_vec = data[:, 2]

    for a in 1:5

        a_occurrences = count(x -> x == a, a_vec)
        
        prob[a, 1] = sum(count_combos.(a, 0, a_vec, i_vec)) / a_occurrences
        prob[a, 2] = 1 - prob[a, 1]

    end

    prob
end

# Given θ, replacement probability matrix, number of periods per simulation, 
# number of simulations, and β, returns 5×2 matrix of simulated firm values:
# values[a, i] = value of firm who begins with choice i at age a.
function simulate_forward_values(θ, replace_probs, periods, sims, β)
    (μ, R) = θ
    values = zeros(5, 2, sims)
    γ = MathConstants.eulergamma

    for sim in 1:sims
        for a in 1:5
            for i in 0:1

                values[a, i + 1, sim] += π(a, i, μ, R)

                α = a
                ι = i

                for period in 1:periods
    
                    if ι == 1
                        α = 1
                    elseif ι == 0 && α < 5
                        α += 1
                    else
                        α = 5
                    end
                        
                    ι = rand(Binomial( 1, replace_probs[α, 2] ))

                    values[a, i + 1, sim] += (β^period)*(π(α, ι, μ, R) + γ - log(replace_probs[α, ι + 1]))
                end
            end
        end
    end
    mean(values, dims = 3)[:, :, 1]
end

# Given θ, data, and probabilities of replacement,
# returns the value of the negative joint log-likelihood function.
function loglikelihood_2(θ, data, replace_probs)
    (μ, R) = θ

    a_vec = data[:, 1]
    i_vec = data[:, 2]

    values = simulate_forward_values(θ, replace_probs, 500, 500, 0.9)

    P_vec = zeros(10)
    Q_vec = zeros(10)

    for index in 1:5
        P_vec[index] = P(index, 0, values)
        Q_vec[index] = sum(count_combos.(index, 0, a_vec, i_vec))
    end

    for index in 6:10
        P_vec[index] = P(index - 5, 1, values)
        Q_vec[index] = sum(count_combos.(index - 5, 1, a_vec, i_vec))
    end

    -sum( log.(P_vec) .* Q_vec )
end

# Get replacement probabilities.
replace_probs = get_replacement_probs(data)

# Initial guesses for μ and R.
μ_hat, R_hat = (-5.0, -5.0)
θ_initial = copy([μ_hat, R_hat])

# Minimize the negative log-likelihood function.
res = @time optimize(
    θ -> loglikelihood_2(θ, data, replace_probs), θ_initial,
    Optim.Options(show_trace = true, iterations = 100)
)

θ_hat = Optim.minimizer(res)

# Print.
println("μ estimate: ", θ_hat[1], "\nR estimate: ", θ_hat[2])

