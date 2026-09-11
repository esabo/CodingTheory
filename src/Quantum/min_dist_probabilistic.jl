# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

function _css_isd_permutation_vector(perm, n::Int)
    result = perm isa AbstractVector ? Int.(perm) : _matrix_to_perm_vector(perm)
    length(result) == n && sort(result) == collect(1:n) ||
        throw(ArgumentError("Every ISD preprocessor permutation must permute 1:$n."))
    return result
end

function _css_isd_preprocessor_permutations(
    G::Matrix{Int}, info_set_alg::Symbol,
    automorphisms::Vector{<:AbstractVector{<:Integer}}
)
    n = size(G, 2)
    automorphism_permutations =
        [_css_isd_permutation_vector(σ, n) for σ in automorphisms]
    permutations = copy(automorphism_permutations)
    info_set_alg == :random && return permutations

    info_set_alg ∈ (:Brouwer, :Zimmermann, :White, :Chen, :Bouyuklieva, :Edmonds) ||
        throw(ArgumentError("Unknown CSS ISD information-set preprocessor `$info_set_alg`."))
    F = Oscar.Nemo.Native.GF(2)
    _, raw_permutations, ranks = information_sets(
        matrix(F, G), info_set_alg; permute=false)
    k = size(G, 1)
    for (perm, rnk) in zip(raw_permutations, ranks)
        rnk == k || continue
        base_permutation = _css_isd_permutation_vector(perm, n)
        push!(permutations, base_permutation)
        for automorphism in automorphism_permutations
            push!(permutations, automorphism[base_permutation])
        end
    end
    return unique(permutations)
end

function _css_isd_prange_candidates!(
    G::Matrix{Int}, row_labels::Matrix{Int}, max_weight::Threads.Atomic{Int},
    report!::F
) where F
    packed_rows = _pack_binary_rows(G)
    @inbounds for row in axes(G, 1)
        any(!iszero, @view row_labels[row, :]) || continue
        weight = _packed_weight(packed_rows[row])
        weight < max_weight[] || continue
        report!(vec(copy(@view G[row, :])), weight)
    end
    return nothing
end

function _css_isd_lee_brickell_candidates!(
    G::Matrix{Int}, row_labels::Matrix{Int}, p::Int,
    max_weight::Threads.Atomic{Int}, report!::F
) where F
    k, n = size(G)
    1 <= p <= k || return nothing
    tail_length = n - k
    tail_rows = _pack_binary_rows(@view G[:, k + 1:n])
    packed_labels = _pack_binary_rows(row_labels)
    current_tail = zeros(UInt64, cld(tail_length, 64))
    current_label = zeros(UInt64, cld(size(row_labels, 2), 64))
    support = Int[]

    function visit!(next_row::Int, remaining::Int)
        if iszero(remaining)
            any(!iszero, current_label) || return
            weight = p + _packed_weight(current_tail)
            weight < max_weight[] || return
            candidate = zeros(Int, n)
            candidate[support] .= 1
            candidate[k + 1:n] .= _unpack_binary_vector(current_tail, tail_length)
            report!(candidate, weight)
            return
        end

        final_row = k - remaining + 1
        for row in next_row:final_row
            push!(support, row)
            _xor_packed!(current_tail, tail_rows[row])
            _xor_packed!(current_label, packed_labels[row])
            visit!(row + 1, remaining - 1)
            _xor_packed!(current_tail, tail_rows[row])
            _xor_packed!(current_label, packed_labels[row])
            pop!(support)
        end
    end
    visit!(1, p)
    return nothing
end

function _css_isd_stern_candidates!(
    G::Matrix{Int}, row_labels::Matrix{Int}, p::Int, l::Int,
    max_weight::Threads.Atomic{Int}, report!::F
) where F
    k, n = size(G)
    half_k = k ÷ 2
    half_k >= 1 || return _css_isd_prange_candidates!(
        G, row_labels, max_weight, report!)
    p = min(max(1, p), half_k, k - half_k)
    0 <= l <= n - k || throw(DomainError(
        l, "The Stern collision window must not exceed the parity tail."))
    l <= 64 || throw(DomainError(
        l, "The Stern collision window must fit in one UInt64."))

    X_rows = 1:half_k
    Y_rows = half_k + 1:k
    tail_length = n - k - l
    tail_rows = _pack_binary_rows(@view G[:, k + l + 1:n])
    packed_labels = _pack_binary_rows(row_labels)
    window_rows = zeros(UInt64, k)
    @inbounds for row in 1:k, col in 1:l
        iszero(G[row, k + col]) ||
            (window_rows[row] |= UInt64(1) << (col - 1))
    end

    Label = typeof(Tuple(zeros(UInt64, cld(size(row_labels, 2), 64))))
    Tail = typeof(Tuple(zeros(UInt64, cld(tail_length, 64))))
    Entry = Tuple{Label, Tail, Vector{Int}}
    table = Dict{UInt64, Vector{Entry}}()
    x_tail = zeros(UInt64, cld(tail_length, 64))
    x_label = zeros(UInt64, cld(size(row_labels, 2), 64))
    x_support = Int[]

    function build_x!(offset::Int, remaining::Int, window::UInt64)
        if iszero(remaining)
            entry = (Tuple(x_label), Tuple(x_tail), copy(x_support))
            push!(get!(table, window, Entry[]), entry)
            return
        end
        final_offset = length(X_rows) - remaining + 1
        for offset_i in offset:final_offset
            row = X_rows[offset_i]
            push!(x_support, row)
            _xor_packed!(x_tail, tail_rows[row])
            _xor_packed!(x_label, packed_labels[row])
            build_x!(offset_i + 1, remaining - 1, window ⊻ window_rows[row])
            _xor_packed!(x_tail, tail_rows[row])
            _xor_packed!(x_label, packed_labels[row])
            pop!(x_support)
        end
    end
    build_x!(1, p, UInt64(0))

    y_tail = zeros(UInt64, cld(tail_length, 64))
    y_label = zeros(UInt64, cld(size(row_labels, 2), 64))
    y_support = Int[]

    function probe_y!(offset::Int, remaining::Int, window::UInt64)
        if iszero(remaining)
            entries = get(table, window, nothing)
            entries === nothing && return
            for (stored_label, stored_tail, stored_support) in entries
                logical = false
                @inbounds for chunk in eachindex(y_label)
                    if y_label[chunk] ⊻ stored_label[chunk] != UInt64(0)
                        logical = true
                        break
                    end
                end
                logical || continue

                candidate_tail = UInt64[stored_tail...]
                _xor_packed!(candidate_tail, y_tail)
                weight = 2p + _packed_weight(candidate_tail)
                weight < max_weight[] || continue
                candidate = zeros(Int, n)
                candidate[stored_support] .= 1
                candidate[y_support] .= 1
                candidate[k + l + 1:n] .=
                    _unpack_binary_vector(candidate_tail, tail_length)
                report!(candidate, weight)
            end
            return
        end

        final_offset = length(Y_rows) - remaining + 1
        for offset_i in offset:final_offset
            row = Y_rows[offset_i]
            push!(y_support, row)
            _xor_packed!(y_tail, tail_rows[row])
            _xor_packed!(y_label, packed_labels[row])
            probe_y!(offset_i + 1, remaining - 1, window ⊻ window_rows[row])
            _xor_packed!(y_tail, tail_rows[row])
            _xor_packed!(y_label, packed_labels[row])
            pop!(y_support)
        end
    end
    probe_y!(1, p, UInt64(0))
    return nothing
end

"""
$(TYPEDSIGNATURES)

Run a threaded, quotient-aware ISD search for a low-weight vector in `ker(H)`
with a nonzero logical label. This is a probabilistic upper-bound algorithm;
failure to find a vector does not certify a lower bound.
"""
function _minimum_distance_css_isd_binary(
    H::_CSSBinaryMatrix, logical_checks::_CSSBinaryMatrix;
    alg::Symbol=:Prange, lower_bound::Int=1, max_weight::Int=ncols(H),
    max_iters::Int=10_000, p::Int=2, l::Int=12,
    info_set_alg::Symbol=:random,
    automorphisms::Vector{<:AbstractVector{<:Integer}}=Vector{Vector{Int}}(),
    seed::Union{Nothing, Integer}=nothing, verbose::Bool=false
)
    alg ∈ (:Prange, :LeeBrickell, :Stern) ||
        throw(ArgumentError("CSS ISD supports `:Prange`, `:LeeBrickell`, and `:Stern`."))
    ncols(H) == ncols(logical_checks) ||
        throw(ArgumentError("The parity-check and logical-check matrices must have the same number of columns."))
    1 <= max_weight <= ncols(H) ||
        throw(DomainError(max_weight, "`max_weight` must lie between one and the code length."))
    1 <= lower_bound <= max_weight ||
        throw(DomainError(lower_bound,
            "`lower_bound` must lie between one and `max_weight`."))
    max_iters > 0 ||
        throw(DomainError(max_iters, "`max_iters` must be positive."))

    F2 = Oscar.Nemo.Native.GF(2)
    H_int = _convert_binary_to_int_matrix(H)
    L_int = _convert_binary_to_int_matrix(logical_checks)
    normalizer = LinearCode(matrix(F2, H_int), true)
    G = _convert_binary_to_int_matrix(generator_matrix(normalizer))
    k, n = size(G)
    permutations = _css_isd_preprocessor_permutations(
        G, info_set_alg, automorphisms)

    best_weight = Threads.Atomic{Int}(max_weight + 1)
    best_witness = zeros(Int, n)
    result_lock = ReentrantLock()
    completed = Threads.Atomic{Int}(0)
    report_interval = max(1, max_iters ÷ 20)
    master_seed = isnothing(seed) ? rand(UInt64) : UInt64(seed)

    Threads.@threads for iteration in 1:max_iters
        best_weight[] == lower_bound && continue
        rng = Random.Xoshiro(master_seed + UInt64(iteration))
        σ = if iteration <= length(permutations)
            copy(permutations[iteration])
        else
            randperm(rng, n)
        end
        G_local = G[:, σ]
        try
            _make_systematic_gf!(G_local, σ, k)
        catch
            continue
        end
        L_local = L_int[:, σ]
        row_labels = (G_local * transpose(L_local)) .% 2

        function report!(candidate::Vector{Int}, weight::Int)
            lock(result_lock) do
                if weight < best_weight[]
                    witness = zeros(Int, n)
                    @inbounds for local_index in 1:n
                        witness[σ[local_index]] = candidate[local_index]
                    end
                    best_witness .= witness
                    best_weight[] = weight
                    verbose && println(
                        "$alg found a logical of weight $weight at iteration $iteration.")
                end
            end
            return nothing
        end

        if alg == :Prange
            _css_isd_prange_candidates!(
                G_local, row_labels, best_weight, report!)
        elseif alg == :LeeBrickell
            _css_isd_lee_brickell_candidates!(
                G_local, row_labels, p, best_weight, report!)
        else
            _css_isd_stern_candidates!(
                G_local, row_labels, p, l, best_weight, report!)
        end

        done = Threads.atomic_add!(completed, 1) + 1
        if verbose && (done % report_interval == 0 || done == max_iters)
            current = best_weight[] <= max_weight ? string(best_weight[]) : "none"
            println("$alg progress: $done/$max_iters information sets; best = $current.")
        end
    end

    return best_weight[] <= max_weight ?
        (best_weight[], best_witness) : (-1, zeros(Int, n))
end

function _record_css_isd_upper_bound!(
    S::AbstractStabilizerCodeCSS, which::Symbol, d::Int, witness
)
    d == -1 && return d, witness
    set_minimum_distance_upper_bound!(S, d, witness; which=which)
    return d, witness
end

"""
$(TYPEDSIGNATURES)

Return an `X`, `Z`, or full-distance upper bound found using quotient-aware
binary ISD. Returned witnesses are symplectic `[X | Z]` rows.
"""
function probabilistic_minimum_distance(
    S::AbstractStabilizerCodeCSS; which::Symbol=:full, alg::Symbol=:Prange,
    max_weight::Int=S.n, max_iters::Int=10_000, p::Int=2, l::Int=12,
    info_set_alg::Symbol=:random,
    automorphisms=nothing,
    seed::Union{Nothing, Integer}=nothing, verbose::Bool=false
)
    which ∈ (:full, :X, :Z) ||
        throw(ArgumentError("Expected `which` to be `:full`, `:X`, or `:Z`."))
    Int(order(field(S))) == 2 ||
        throw(ArgumentError("CSS ISD is currently implemented only over GF(2)."))
    S.k > 0 || throw(ArgumentError(
        "Logical minimum distance is undefined for a CSS stabilizer state with k = 0."))
    resolved_automorphisms = isnothing(automorphisms) ?
        distance_automorphisms(S) : automorphisms

    if which == :full
        dx, wx = probabilistic_minimum_distance(
            S; which=:X, alg=alg, max_weight=max_weight,
            max_iters=max_iters, p=p, l=l, info_set_alg=info_set_alg,
            automorphisms=resolved_automorphisms, seed=seed, verbose=verbose)
        z_max = dx == -1 ? max_weight : min(max_weight, dx)
        dz, wz = probabilistic_minimum_distance(
            S; which=:Z, alg=alg, max_weight=z_max,
            max_iters=max_iters, p=p, l=l, info_set_alg=info_set_alg,
            automorphisms=resolved_automorphisms,
            seed=isnothing(seed) ? nothing : seed + 1, verbose=verbose)
        return dz == -1 || (dx != -1 && dx <= dz) ? (dx, wx) : (dz, wz)
    end

    keys = _css_distance_cache_keys(which)
    if haskey(S.cache, keys.exact)
        return _css_cached_distance_result(S, which)
    end
    lower = get(S.cache, keys.lower, 1)
    upper = get(S.cache, keys.upper, S.n)
    cached_witness = get(S.cache, keys.upper_witness, nothing)
    search_max = min(max_weight, isnothing(cached_witness) ? upper : upper - 1)
    if search_max < lower
        return isnothing(cached_witness) ?
            (-1, zero_matrix(S.F, 1, 2 * S.n)) :
            (upper, cached_witness)
    end

    H, logical_checks = _css_distance_problem(S, which)
    d, vector_witness = _minimum_distance_css_isd_binary(
        H, logical_checks; alg=alg, lower_bound=lower,
        max_weight=search_max,
        max_iters=max_iters, p=p, l=l, info_set_alg=info_set_alg,
        automorphisms=resolved_automorphisms, seed=seed, verbose=verbose)
    if d == -1 && !isnothing(cached_witness)
        return upper, cached_witness
    end
    witness = d == -1 ?
        zero_matrix(S.F, 1, 2 * S.n) :
        _css_quantum_witness(S, which, vector_witness)
    return _record_css_isd_upper_bound!(S, which, d, witness)
end

probabilistic_minimum_distance_prange(
    S::AbstractStabilizerCodeCSS; kwargs...
) = probabilistic_minimum_distance(S; alg=:Prange, kwargs...)

probabilistic_minimum_distance_lee_brickell(
    S::AbstractStabilizerCodeCSS; kwargs...
) = probabilistic_minimum_distance(S; alg=:LeeBrickell, kwargs...)

probabilistic_minimum_distance_stern(
    S::AbstractStabilizerCodeCSS; kwargs...
) = probabilistic_minimum_distance(S; alg=:Stern, kwargs...)
