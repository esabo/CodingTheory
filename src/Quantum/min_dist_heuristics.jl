# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

function _css_heuristic_generator(H::_CSSBinaryMatrix, logical_checks::_CSSBinaryMatrix)
    F2 = Oscar.Nemo.Native.GF(2)
    H_int = _convert_binary_to_int_matrix(H)
    L_int = _convert_binary_to_int_matrix(logical_checks)
    C = LinearCode(matrix(F2, H_int), true)
    return _convert_binary_to_int_matrix(generator_matrix(C)), L_int
end

function _css_evaluate_message(
    message::AbstractVector{Bool}, packed_rows::Vector{Vector{UInt64}},
    packed_labels::Vector{Vector{UInt64}}, n::Int
)
    word = zeros(UInt64, cld(n, 64))
    label = isempty(packed_labels) ? UInt64[] :
        zeros(UInt64, length(first(packed_labels)))
    @inbounds for row in eachindex(message)
        message[row] || continue
        _xor_packed!(word, packed_rows[row])
        _xor_packed!(label, packed_labels[row])
    end
    any(!iszero, label) || return n + 1, word
    return _packed_weight(word), word
end

function _css_evaluate_information_set(
    G::Matrix{Int}, logical_checks::Matrix{Int}, permutation::Vector{Int}
)
    k, n = size(G)
    σ = copy(permutation)
    G_local = G[:, σ]
    try
        _make_systematic_gf!(G_local, σ, k)
    catch
        return n + 1, zeros(Int, n)
    end
    labels = (G_local * transpose(logical_checks[:, σ])) .% 2
    packed_rows = _pack_binary_rows(G_local)
    best_weight = n + 1
    best_local = zeros(Int, n)
    for row in 1:k
        any(!iszero, @view labels[row, :]) || continue
        weight = _packed_weight(packed_rows[row])
        if weight < best_weight
            best_weight = weight
            best_local .= @view G_local[row, :]
        end
    end
    best_weight == n + 1 && return best_weight, best_local

    witness = zeros(Int, n)
    @inbounds for local_index in 1:n
        witness[σ[local_index]] = best_local[local_index]
    end
    return best_weight, witness
end

function _css_heuristic_ga(
    G::Matrix{Int}, logical_checks::Matrix{Int}, lower_bound::Int,
    upper_bound::Int; max_gens::Int=1_000, pop_size::Int=200,
    mutation_rate::Float64=0.02, seed::Integer=0, verbose::Bool=false
)
    k, n = size(G)
    packed_rows = _pack_binary_rows(G)
    packed_labels = _pack_binary_rows((G * transpose(logical_checks)) .% 2)
    population = Vector{BitVector}(undef, pop_size)
    for individual in 1:pop_size
        rng = Random.Xoshiro(UInt64(seed) + UInt64(individual))
        message = falses(k)
        support_size = rand(rng, 1:min(k, max(1, upper_bound)))
        message[randperm(rng, k)[1:support_size]] .= true
        population[individual] = message
    end

    best_weight = upper_bound + 1
    best_word = zeros(UInt64, cld(n, 64))
    for generation in 1:max_gens
        fitness = fill(n + 1, pop_size)
        words = [zeros(UInt64, cld(n, 64)) for _ in 1:pop_size]
        Threads.@threads for individual in 1:pop_size
            fitness[individual], words[individual] = _css_evaluate_message(
                population[individual], packed_rows, packed_labels, n)
        end
        generation_best = argmin(fitness)
        if fitness[generation_best] < best_weight
            best_weight = fitness[generation_best]
            best_word .= words[generation_best]
            verbose && println("CSS GA generation $generation: upper bound $best_weight.")
            best_weight == lower_bound && break
        end

        order = sortperm(fitness)
        elite_count = max(1, pop_size ÷ 10)
        parent_pool = max(2, pop_size ÷ 3)
        next_population = Vector{BitVector}(undef, pop_size)
        for individual in 1:elite_count
            next_population[individual] = copy(population[order[individual]])
        end
        for individual in elite_count + 1:pop_size
            rng = Random.Xoshiro(
                UInt64(seed) + UInt64(generation) * UInt64(pop_size) +
                UInt64(individual))
            parent1 = population[order[rand(rng, 1:parent_pool)]]
            parent2 = population[order[rand(rng, 1:parent_pool)]]
            point1, point2 = sort(rand(rng, 1:k, 2))
            child = copy(parent1)
            child[point1:point2] .= parent2[point1:point2]
            @inbounds for bit in 1:k
                rand(rng) < mutation_rate && (child[bit] = !child[bit])
            end
            any(child) || (child[rand(rng, 1:k)] = true)
            next_population[individual] = child
        end
        population = next_population
    end
    return best_weight <= upper_bound ?
        (best_weight, _unpack_binary_vector(best_word, n)) :
        (-1, zeros(Int, n))
end

function _css_heuristic_aco(
    G::Matrix{Int}, logical_checks::Matrix{Int}, lower_bound::Int,
    upper_bound::Int; max_iters::Int=500, num_ants::Int=0,
    seed::Integer=0, verbose::Bool=false
)
    k, n = size(G)
    num_ants = iszero(num_ants) ? k : num_ants
    packed_rows = _pack_binary_rows(G)
    packed_labels = _pack_binary_rows((G * transpose(logical_checks)) .% 2)
    pheromone = ones(Float64, k, 2)
    visibility = fill(0.5, k, 2)
    best_weight = upper_bound + 1
    best_word = zeros(UInt64, cld(n, 64))

    for iteration in 1:max_iters
        messages = Vector{BitVector}(undef, num_ants)
        fitness = fill(n + 1, num_ants)
        words = [zeros(UInt64, cld(n, 64)) for _ in 1:num_ants]
        Threads.@threads for ant in 1:num_ants
            rng = Random.Xoshiro(
                UInt64(seed) + UInt64(iteration) * UInt64(num_ants) +
                UInt64(ant))
            message = falses(k)
            for bit in 1:k
                count(message) >= min(best_weight, upper_bound) && break
                p0 = pheromone[bit, 1] * visibility[bit, 1]^2
                p1 = pheromone[bit, 2] * visibility[bit, 2]^2
                message[bit] = rand(rng) < p1 / (p0 + p1)
            end
            any(message) || (message[rand(rng, 1:k)] = true)
            messages[ant] = message
            fitness[ant], words[ant] = _css_evaluate_message(
                message, packed_rows, packed_labels, n)
        end

        deposits = zeros(Float64, k, 2)
        for ant in 1:num_ants
            weight = fitness[ant]
            weight <= upper_bound || continue
            if weight < best_weight
                best_weight = weight
                best_word .= words[ant]
                verbose && println("CSS ACO iteration $iteration: upper bound $best_weight.")
            end
            deposit = 10.0 / weight
            @inbounds for bit in 1:k
                deposits[bit, messages[ant][bit] ? 2 : 1] += deposit
            end
        end
        pheromone .= 0.5 .* pheromone .+ deposits
        best_weight == lower_bound && break
    end
    return best_weight <= upper_bound ?
        (best_weight, _unpack_binary_vector(best_word, n)) :
        (-1, zeros(Int, n))
end

function _css_heuristic_gga_order(
    G::Matrix{Int}, logical_checks::Matrix{Int}, lower_bound::Int,
    upper_bound::Int; max_gens::Int=500, pop_size::Int=50,
    automorphisms::Vector{<:AbstractVector{<:Integer}}=Vector{Vector{Int}}(),
    seed::Integer=0, verbose::Bool=false
)
    n = size(G, 2)
    population = Vector{Vector{Int}}()
    append!(population, [Int.(σ) for σ in automorphisms])
    while length(population) < pop_size
        rng = Random.Xoshiro(UInt64(seed) + UInt64(length(population) + 1))
        push!(population, randperm(rng, n))
    end
    resize!(population, pop_size)

    best_weight = upper_bound + 1
    best_witness = zeros(Int, n)
    for generation in 1:max_gens
        fitness = fill(n + 1, pop_size)
        witnesses = [zeros(Int, n) for _ in 1:pop_size]
        Threads.@threads for individual in 1:pop_size
            fitness[individual], witnesses[individual] =
                _css_evaluate_information_set(
                    G, logical_checks, population[individual])
        end
        generation_best = argmin(fitness)
        if fitness[generation_best] < best_weight
            best_weight = fitness[generation_best]
            best_witness .= witnesses[generation_best]
            verbose && println(
                "CSS GGA-Order generation $generation: upper bound $best_weight.")
            best_weight == lower_bound && break
        end

        order = sortperm(fitness)
        elite_count = max(2, pop_size ÷ 10)
        next_population = Vector{Vector{Int}}(undef, pop_size)
        for individual in 1:elite_count
            next_population[individual] = copy(population[order[individual]])
        end
        for individual in elite_count + 1:pop_size
            rng = Random.Xoshiro(
                UInt64(seed) + UInt64(generation) * UInt64(pop_size) +
                UInt64(individual))
            parent1 = population[order[rand(rng, 1:pop_size)]]
            parent2 = population[order[rand(rng, 1:pop_size)]]
            child = parent1[parent2]
            if rand(rng) < 0.1
                i, j = rand(rng, 1:n, 2)
                child[i], child[j] = child[j], child[i]
            end
            next_population[individual] = child
        end
        population = next_population
    end
    return best_weight <= upper_bound ?
        (best_weight, best_witness) : (-1, zeros(Int, n))
end

function _css_heuristic_nncs(
    G::Matrix{Int}, logical_checks::Matrix{Int}, lower_bound::Int,
    upper_bound::Int; automorphisms::Vector{<:AbstractVector{<:Integer}}=
        Vector{Vector{Int}}(), verbose::Bool=false
)
    n = size(G, 2)
    permutations = [collect(1:n)]
    append!(permutations, [Int.(σ) for σ in automorphisms])
    for qubit in 2:n
        permutation = collect(1:n)
        permutation[1], permutation[qubit] =
            permutation[qubit], permutation[1]
        push!(permutations, permutation)
    end

    fitness = fill(n + 1, length(permutations))
    witnesses = [zeros(Int, n) for _ in eachindex(permutations)]
    Threads.@threads for index in eachindex(permutations)
        fitness[index], witnesses[index] = _css_evaluate_information_set(
            G, logical_checks, permutations[index])
    end
    best = argmin(fitness)
    verbose && fitness[best] <= upper_bound &&
        println("CSS NNCS upper bound: $(fitness[best]).")
    return fitness[best] <= upper_bound ?
        (fitness[best], witnesses[best]) : (-1, zeros(Int, n))
end

function _minimum_distance_css_heuristic_binary(
    H::_CSSBinaryMatrix, logical_checks::_CSSBinaryMatrix;
    alg::Symbol=:GGAOrder, lower_bound::Int=1,
    upper_bound::Int=ncols(H), max_iters::Int=500,
    pop_size::Int=50, mutation_rate::Float64=0.02,
    automorphisms::Vector{<:AbstractVector{<:Integer}}=Vector{Vector{Int}}(),
    seed::Integer=0, verbose::Bool=false
)
    alg ∈ (:GA, :ACO, :GGAOrder, :NNCS) ||
        throw(ArgumentError("CSS heuristics support `:GA`, `:ACO`, `:GGAOrder`, and `:NNCS`."))
    1 <= lower_bound <= upper_bound <= ncols(H) ||
        throw(DomainError((lower_bound, upper_bound),
            "CSS heuristic bounds must satisfy 1 ≤ lower ≤ upper ≤ n."))
    max_iters > 0 ||
        throw(DomainError(max_iters, "`max_iters` must be positive."))
    pop_size >= 2 ||
        throw(DomainError(pop_size, "`pop_size` must be at least two."))
    0.0 <= mutation_rate <= 1.0 ||
        throw(DomainError(mutation_rate,
            "`mutation_rate` must lie between zero and one."))
    G, L = _css_heuristic_generator(H, logical_checks)
    if alg == :GA
        return _css_heuristic_ga(
            G, L, lower_bound, upper_bound; max_gens=max_iters,
            pop_size=pop_size, mutation_rate=mutation_rate,
            seed=seed, verbose=verbose)
    elseif alg == :ACO
        return _css_heuristic_aco(
            G, L, lower_bound, upper_bound; max_iters=max_iters,
            num_ants=pop_size, seed=seed, verbose=verbose)
    elseif alg == :GGAOrder
        return _css_heuristic_gga_order(
            G, L, lower_bound, upper_bound; max_gens=max_iters,
            pop_size=pop_size, automorphisms=automorphisms,
            seed=seed, verbose=verbose)
    end
    return _css_heuristic_nncs(
        G, L, lower_bound, upper_bound;
        automorphisms=automorphisms, verbose=verbose)
end

"""
$(TYPEDSIGNATURES)

Return a CSS distance upper bound found using a quotient-aware metaheuristic.
These algorithms never certify a lower bound.
"""
function heuristic_minimum_distance(
    S::AbstractStabilizerCodeCSS; which::Symbol=:full,
    alg::Symbol=:GGAOrder, max_iters::Int=500, pop_size::Int=50,
    mutation_rate::Float64=0.02, automorphisms=nothing,
    seed::Integer=0, verbose::Bool=false
)
    which ∈ (:full, :X, :Z) ||
        throw(ArgumentError("Expected `which` to be `:full`, `:X`, or `:Z`."))
    Int(order(field(S))) == 2 ||
        throw(ArgumentError("CSS heuristics are currently implemented only over GF(2)."))
    S.k > 0 || throw(ArgumentError(
        "Logical minimum distance is undefined for a CSS stabilizer state with k = 0."))
    resolved_automorphisms = isnothing(automorphisms) ?
        distance_automorphisms(S) :
        [_css_isd_permutation_vector(σ, S.n) for σ in automorphisms]

    if which == :full
        dx, wx = heuristic_minimum_distance(
            S; which=:X, alg=alg, max_iters=max_iters,
            pop_size=pop_size, mutation_rate=mutation_rate,
            automorphisms=resolved_automorphisms, seed=seed,
            verbose=verbose)
        dz, wz = heuristic_minimum_distance(
            S; which=:Z, alg=alg, max_iters=max_iters,
            pop_size=pop_size, mutation_rate=mutation_rate,
            automorphisms=resolved_automorphisms, seed=seed + 1,
            verbose=verbose)
        return dz == -1 || (dx != -1 && dx <= dz) ? (dx, wx) : (dz, wz)
    end

    keys = _css_distance_cache_keys(which)
    haskey(S.cache, keys.exact) &&
        return _css_cached_distance_result(S, which)
    lower = get(S.cache, keys.lower, 1)
    upper = get(S.cache, keys.upper, S.n)
    cached_witness = get(S.cache, keys.upper_witness, nothing)
    search_upper = isnothing(cached_witness) ? upper : upper - 1
    search_upper < lower &&
        return isnothing(cached_witness) ?
            (-1, zero_matrix(S.F, 1, 2 * S.n)) :
            (upper, cached_witness)

    H, logical_checks = _css_distance_problem(S, which)
    d, vector_witness = _minimum_distance_css_heuristic_binary(
        H, logical_checks; alg=alg, lower_bound=lower,
        upper_bound=search_upper, max_iters=max_iters,
        pop_size=pop_size, mutation_rate=mutation_rate,
        automorphisms=resolved_automorphisms, seed=seed,
        verbose=verbose)
    if d == -1
        return isnothing(cached_witness) ?
            (-1, zero_matrix(S.F, 1, 2 * S.n)) :
            (upper, cached_witness)
    end
    witness = _css_quantum_witness(S, which, vector_witness)
    set_minimum_distance_upper_bound!(S, d, witness; which=which)
    return d, witness
end

heuristic_minimum_distance_ga(
    S::AbstractStabilizerCodeCSS; kwargs...
) = heuristic_minimum_distance(S; alg=:GA, kwargs...)
heuristic_minimum_distance_aco(
    S::AbstractStabilizerCodeCSS; kwargs...
) = heuristic_minimum_distance(S; alg=:ACO, kwargs...)
heuristic_minimum_distance_gga_order(
    S::AbstractStabilizerCodeCSS; kwargs...
) = heuristic_minimum_distance(S; alg=:GGAOrder, kwargs...)
heuristic_minimum_distance_nncs(
    S::AbstractStabilizerCodeCSS; kwargs...
) = heuristic_minimum_distance(S; alg=:NNCS, kwargs...)
