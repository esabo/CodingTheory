# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
    _minimum_distance_GA(C::AbstractLinearCode; max_gens::Int=1000, pop_size::Int=200, verbose::Bool=false)

Standalone Genetic Algorithm solver for the minimum distance problem based on Askali et al. (2013).
Uses Tournament Selection, 2-Point Crossover, and Elitism.
"""
function _minimum_distance_GA(C::AbstractLinearCode; max_gens::Int=1000, pop_size::Int=200, verbose::Bool=false)
    G = _convert_binary_to_int_matrix(generator_matrix(C, true))
    k, n = size(G)
    
    # Initialize Population
    population = Vector{Vector{Int}}(undef, pop_size)
    for i in 1:pop_size
        msg = zeros(Int, k)
        target_wt = rand(1:k)
        indices = randperm(k)[1:target_wt]
        msg[indices] .= 1
        population[i] = msg
    end
    
    best_overall_wt = n + 1
    best_overall_cw = zeros(Int, n)
    
    function evaluate_fitness(m::Vector{Int})
        cw = vec(Array((m' * G) .% 2))
        wt = count(!iszero, cw)
        return wt == 0 ? (n + 1) : wt, cw
    end

    for gen in 1:max_gens
        fitnesses = [evaluate_fitness(m) for m in population]
        
        # Track Global Best
        for (i, (wt, cw)) in enumerate(fitnesses)
            if wt < best_overall_wt
                best_overall_wt = wt
                copyto!(best_overall_cw, cw)
                verbose && println("Generation $gen: New minimum distance found -> $best_overall_wt")
            end
        end
        
        # Sort for Elitism (Top 10%)
        sort_idx = sortperm([f[1] for f in fitnesses])
        next_gen = Vector{Vector{Int}}(undef, pop_size)
        elite_count = div(pop_size, 10)
        
        for i in 1:elite_count
            next_gen[i] = copy(population[sort_idx[i]])
        end
        
        # Tournament Selection & 2-Point Crossover
        for i in (elite_count + 1):pop_size
            p1_idx = sort_idx[rand(1:3)]
            p2_idx = sort_idx[rand(1:3)]
            parent1 = population[p1_idx]
            parent2 = population[p2_idx]
            
            pt1, pt2 = sort(rand(1:k, 2))
            child = copy(parent1)
            child[pt1:pt2] .= parent2[pt1:pt2]
            
            # Classic Mutation
            for b in 1:k
                if rand() < 0.02 # 2% mutation rate as per Askali
                    child[b] ⊻= 1
                end
            end
            
            next_gen[i] = child
        end
        
        population = next_gen
    end
    
    witness = matrix(C.F, 1, n, [c == 1 ? one(C.F) : zero(C.F) for c in best_overall_cw])
    return best_overall_wt, witness * C.P_stand
end

"""
    _minimum_distance_ACO(C::AbstractLinearCode; y_max::Int=500, m_ants::Int=100, verbose::Bool=false)

Standalone Ant Colony Optimization solver based on Bouzkraoui et al. (2018).
Ants iteratively construct message vectors guided by pheromone trails and visibility.
"""
function _minimum_distance_ACO(C::AbstractLinearCode; y_max::Int=500, m_ants::Int=0, verbose::Bool=false)
    G = _convert_binary_to_int_matrix(generator_matrix(C, true))
    k, n = size(G)
    
    # Optimal ACO Parameters (Bouzkraoui et al.)
    α = 1.0; β = 2.0; ρ = 0.5; Q = 10.0
    m_ants = m_ants == 0 ? k : m_ants # Recommendation: m ≈ k
    
    # Pheromone matrix τ and Visibility matrix η 
    # Dimensions: k bits x 2 choices (0 or 1)
    τ = fill(1.0, k, 2)
    η = fill(0.5, k, 2) 
    
    best_overall_wt = n + 1
    best_overall_cw = zeros(Int, n)
    
    for y in 1:y_max
        # Pheromone tracking for this iteration
        Δτ = zeros(Float64, k, 2)
        
        for ant in 1:m_ants
            msg = zeros(Int, k)
            
            # Ant constructs the message bit by bit
            for i in 1:k
                # Bouzkraoui Weight Culling Rule: w(c) >= w(m)
                if count(!iszero, msg) >= best_overall_wt
                    msg[i] = 0
                    continue
                end
                
                # Probabilistic Choice
                prob_0 = (τ[i, 1]^α) * (η[i, 1]^β)
                prob_1 = (τ[i, 2]^α) * (η[i, 2]^β)
                total_prob = prob_0 + prob_1
                
                if rand() < (prob_1 / total_prob)
                    msg[i] = 1
                end
            end
            
            # Evaluate Code Weight
            cw = vec(Array((msg' * G) .% 2))
            wt = count(!iszero, cw)
            
            if wt > 0
                if wt < best_overall_wt
                    best_overall_wt = wt
                    copyto!(best_overall_cw, cw)
                    verbose && println("ACO Iteration $y: New minimum distance found -> $best_overall_wt")
                end
                
                # Ant deposits pheromone based on solution quality
                deposit = Q / wt
                for i in 1:k
                    choice_idx = msg[i] + 1
                    Δτ[i, choice_idx] += deposit
                end
            end
        end
        
        # Evaporation and Global Pheromone Update
        for i in 1:k
            for j in 1:2
                τ[i, j] = (1 - ρ) * τ[i, j] + Δτ[i, j]
            end
        end
    end
    
    witness = matrix(C.F, 1, n, [c == 1 ? one(C.F) : zero(C.F) for c in best_overall_cw])
    return best_overall_wt, witness * C.P_stand
end

"""
$(TYPEDSIGNATURES)

Executes a standalone Genetic Algorithm (Askali et al., 2013) to estimate the minimum distance of a linear code.
Return `(d_upper_bound, witness)`. Note that this is a probabilistic metaheuristic and does NOT 
mathematically guarantee the true minimum distance.
"""
function heuristic_minimum_distance_ga(C::AbstractLinearCode; max_gens::Int=1000, pop_size::Int=200, verbose::Bool=false)
    # Early exit if the exact distance is already known
    if !ismissing(C.d)
        verbose && println("Minimum distance already known: $(C.d). Bypassing GA heuristic.")
        return C.d, (isdefined(C, :witness) ? C.witness : zero_matrix(C.F, 1, C.n))
    end
    
    # Execute the internal GA algorithm
    d_est, witness = _minimum_distance_GA(C; max_gens=max_gens, pop_size=pop_size, verbose=verbose)
    
    # We strictly update the upper bound, NEVER the exact distance C.d
    if ismissing(C.u_bound) || d_est < C.u_bound
        C.u_bound = d_est
    end
    
    return d_est, witness
end

"""
$(TYPEDSIGNATURES)

Executes a standalone Ant Colony Optimization algorithm (Bouzkraoui et al., 2018) to estimate 
the minimum distance of a linear code.
Return `(d_upper_bound, witness)`. Note that this is a probabilistic metaheuristic and does NOT 
mathematically guarantee the true minimum distance.
"""
function heuristic_minimum_distance_aco(C::AbstractLinearCode; y_max::Int=500, m_ants::Int=0, verbose::Bool=false)
    # Early exit if the exact distance is already known
    if !ismissing(C.d)
        verbose && println("Minimum distance already known: $(C.d). Bypassing ACO heuristic.")
        return C.d, (isdefined(C, :witness) ? C.witness : zero_matrix(C.F, 1, C.n))
    end

    # Execute the internal ACO algorithm
    d_est, witness = _minimum_distance_ACO(C; y_max=y_max, m_ants=m_ants, verbose=verbose)
    
    # We strictly update the upper bound, NEVER the exact distance C.d
    if ismissing(C.u_bound) || d_est < C.u_bound
        C.u_bound = d_est
    end
    
    return d_est, witness
end

"""
    _minimum_distance_GGA_order(C::AbstractLinearCode; max_gens::Int=500, pop_size::Int=50, verbose::Bool=false)

Implements the GGA-Order metaheuristic from Cuéllar et al. (2020).
Searches the permutation space `S_n` to find a column permutation whose RREF 
contains a minimum-weight row. Native to any finite field F_q.
"""
function _minimum_distance_GGA_order(C::AbstractLinearCode; max_gens::Int=500, pop_size::Int=50, verbose::Bool=false)
    G = generator_matrix(C, true)
    k, n = size(G)
    F = base_ring(G)
    
    # 1. Initialize Population (Permutations of 1:n)
    population = [randperm(n) for _ in 1:pop_size]
    
    best_overall_wt = n + 1
    best_overall_cw = zeros(F, n)
    
    # Fitness Evaluator: Apply permutation, RREF, and find the lightest row
    function evaluate_fitness(perm::Vector{Int})
        Gp = G[:, perm]
        rnk, R = rref(Gp)
        
        min_wt = n + 1
        best_row_idx = -1
        
        for i in 1:rnk
            row_wt = count(!iszero, view(R, i, :))
            if row_wt > 0 && row_wt < min_wt
                min_wt = row_wt
                best_row_idx = i
            end
        end
        
        # Reconstruct the codeword in the ORIGINAL coordinate space
        cw_original = zeros(F, n)
        if best_row_idx != -1
            cw_permuted = vec(Array(R[best_row_idx, :]))
            # If `perm` maps original col `i` to permuted col `perm[i]`, 
            # we place the permuted value back into its original slot.
            for i in 1:n
                cw_original[perm[i]] = cw_permuted[i]
            end
        end
        
        return min_wt, cw_original
    end

    for gen in 1:max_gens
        fitnesses = [evaluate_fitness(p) for p in population]
        
        # Track Global Best
        for (i, (wt, cw)) in enumerate(fitnesses)
            if wt < best_overall_wt
                best_overall_wt = wt
                copyto!(best_overall_cw, cw)
                verbose && println("GGA-Order Gen $gen: New upper bound found -> $best_overall_wt")
            end
        end
        
        # Sort for Elitism
        sort_idx = sortperm([f[1] for f in fitnesses])
        next_gen = Vector{Vector{Int}}(undef, pop_size)
        
        # Keep the top 10% (Elitism)
        elite_count = max(2, div(pop_size, 10))
        for i in 1:elite_count
            next_gen[i] = copy(population[sort_idx[i]])
        end
        
        # Algebraic Crossover (AX_2) & 2-Swap Mutation
        for i in (elite_count + 1):2:pop_size
            p1 = population[sort_idx[rand(1:pop_size)]]
            p2 = population[sort_idx[rand(1:pop_size)]]
            
            # Crossover: Composition of permutations (p1 ∘ p2 and p2 ∘ p1)
            child1 = [p1[p2[j]] for j in 1:n]
            child2 = [p2[p1[j]] for j in 1:n]
            
            # 2-Swap Mutation: Swap a pivot (1:k) with a non-pivot (k+1:n)
            for child in (child1, child2)
                if rand() < 0.10 # 10% chance to mutate
                    idx1 = rand(1:k)
                    idx2 = rand((k+1):n)
                    child[idx1], child[idx2] = child[idx2], child[idx1]
                end
            end
            
            next_gen[i] = child1
            if i + 1 <= pop_size
                next_gen[i+1] = child2
            end
        end
        
        population = next_gen
    end
    
    witness = matrix(F, 1, n, best_overall_cw)
    return best_overall_wt, witness
end

"""
$(TYPEDSIGNATURES)

Executes the GGA-Order metaheuristic (Cuéllar et al., 2020) to estimate the minimum distance.
Explores the permutation space `S_n` rather than the discrete message space `F_q^k`.
Return `(d_upper_bound, witness)`. Does NOT mathematically guarantee the true minimum distance.
"""
function heuristic_minimum_distance_gga_order(C::AbstractLinearCode; max_gens::Int=500, pop_size::Int=50, verbose::Bool=false)
    if !ismissing(C.d)
        verbose && println("Minimum distance already known: $(C.d). Bypassing GGA-Order heuristic.")
        return C.d, (isdefined(C, :witness) ? C.witness : zero_matrix(C.F, 1, C.n))
    end
    
    d_est, witness = _minimum_distance_GGA_order(C; max_gens=max_gens, pop_size=pop_size, verbose=verbose)
    
    if ismissing(C.u_bound) || d_est < C.u_bound
        C.u_bound = d_est
    end
    
    return d_est, witness
end

"""
$(TYPEDSIGNATURES)

Implements the non-binary probabilistic algorithm (Algorithm 5) from Irons (2005).
Randomly selects combinations of up to `p_max` rows of the generator matrix and evaluates their weights.
Extremely fast Monte Carlo estimation for generic non-binary codes.
"""
function heuristic_minimum_distance_irons(C::AbstractLinearCode; num_iters::Int=1000, p_max::Int=3, verbose::Bool=false)
    if !ismissing(C.d) return C.d, (isdefined(C, :witness) ? C.witness : zero_matrix(C.F, 1, C.n)) end
    
    G = generator_matrix(C, true)
    k, n = size(G)
    F = base_ring(G)
    nonzero_scalars = filter(!iszero, collect(F))
    
    best_wt = n + 1
    best_cw = zeros(F, n)
    
    for _ in 1:num_iters
        # Pick a random linear combination size up to p_max
        p = rand(1:min(p_max, k))
        row_indices = randperm(k)[1:p]
        
        # Assign random non-zero scalars to those rows
        scalars = [rand(nonzero_scalars) for _ in 1:p]
        
        cw = zeros(F, n)
        for (idx, r) in enumerate(row_indices)
            sc = scalars[idx]
            for c in 1:n
                cw[c] += sc * G[r, c]
            end
        end
        
        wt = count(!iszero, cw)
        if wt > 0 && wt < best_wt
            best_wt = wt
            copyto!(best_cw, cw)
            verbose && println("Irons Randomized Search: New minimum found -> $best_wt")
        end
    end
    
    if ismissing(C.u_bound) || best_wt < C.u_bound C.u_bound = best_wt end
    return best_wt, matrix(F, 1, n, best_cw)
end

"""
$(TYPEDSIGNATURES)

Implements the Nearest Nonzero Codeword Search (NNCS) using the Bit-Reversing noise 
pattern from Hu et al. (2004). Uses an OSD-0 (Order Statistic Decoding) hard-decision 
erasure strategy tailored for LDPC codes.
"""
function heuristic_minimum_distance_nncs(C::AbstractLinearCode; verbose::Bool=false)
    if !ismissing(C.d) return C.d, (isdefined(C, :witness) ? C.witness : zero_matrix(C.F, 1, C.n)) end
    
    G = generator_matrix(C, true)
    k, n = size(G)
    F = base_ring(G)
    
    best_wt = n + 1
    best_cw = zeros(F, n)
    
    # NNCS Bit-Reversing: We systematically force each bit i to be part of the basis.
    for i in 1:n
        perm = collect(1:n)
        perm[1], perm[i] = perm[i], perm[1]
        
        Gp = G[:, perm]
        rnk, G_rref = rref(Gp)
        
        if rnk < k continue end
        
        # The first row of the RREF corresponds to a message m = (1, 0, ..., 0).
        # This perfectly simulates Bit-Reversing (bit i = 1) with an OSD-0 assumption 
        # that the remaining basis bits are uncorrupted (0).
        cw_permuted = vec(Array(G_rref[1, :]))
        wt = count(!iszero, cw_permuted)
        
        if wt > 0 && wt < best_wt
            best_wt = wt
            cw_original = zeros(F, n)
            for j in 1:n cw_original[perm[j]] = cw_permuted[j] end
            copyto!(best_cw, cw_original)
            verbose && println("NNCS Bit-Reversing: New minimum found -> $best_wt")
        end
    end
    
    if ismissing(C.u_bound) || best_wt < C.u_bound C.u_bound = best_wt end
    return best_wt, matrix(F, 1, n, best_cw)
end

"""
$(TYPEDSIGNATURES)

Return a `Dict` of weight to estimated count for the low-weight part of the weight
distribution ``A_w`` of an LDPC code, using the method of Hirotomo et al. (2005).
Stern's algorithm is run `t_iters` times to populate the frequency spectrum ``B_w``,
and the multiplicity is estimated by dividing by Stern's success probability
``\\pi_w``.
"""
function heuristic_weight_distribution_hirotomo(C::AbstractLinearCode, target_w::Int; t_iters::Int=1000)
    k, n = C.k, C.n
    p = 2 # Standard Stern parameter for LDPC
    l = min(16, n - k)
    
    B_w = Dict{Int, Int}()
    
    # Run standard Stern attack `t_iters` times to populate B_w
    for _ in 1:t_iters
        found_vecs = Stern_attack(C, target_w; p=p, l=l, num_find=0, max_iters=1, unroll=true)
        for v in found_vecs
            wt = count(!iszero, v)
            B_w[wt] = get(B_w, wt, 0) + 1
        end
    end
    
    A_w_estimate = Dict{Int, Float64}()
    
    # Calculate pi_w and the Hirotomo approximation A_w = B_w / (pi_w * t)
    for (w, count) in B_w
        # Hirotomo Equation (3) calculation
        num = binomial(big(w), big(p)) * binomial(big(n - w), big(div(k, 2) - p)) * binomial(big(w - p), big(p)) * binomial(big(n - w - div(k, 2) + p), big(div(k, 2) - p))
        den = binomial(big(n), big(div(k, 2))) * binomial(big(n - div(k, 2)), big(div(k, 2)))
        
        prob_base = Float64(num) / Float64(den)
        
        prob_tail = Float64(binomial(big(n - k - w + 2*p), big(l))) / Float64(binomial(big(n - k), big(l)))
        
        pi_w = prob_base * prob_tail
        
        if pi_w > 0
            A_w_estimate[w] = count / (pi_w * t_iters)
        end
    end
    
    return A_w_estimate
end