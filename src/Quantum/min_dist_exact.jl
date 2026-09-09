# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

const _CSSBinaryMatrix = Union{CTMatrixTypes, AbstractMatrix}

"""
    _minimum_distance_css_ILP(H, logical_checks; kwargs...)

Extension point for the quotient-aware CSS integer-programming solver. Load
JuMP and GLPK to activate it.
"""
function _minimum_distance_css_ILP end

function _pack_binary_columns(A::_CSSBinaryMatrix)
    nr, nc = size(A)
    chunks = cld(nr, 64)
    packed = [zeros(UInt64, chunks) for _ in 1:nc]

    if A isa SparseMatrixCSC
        vals = SparseArrays.nonzeros(A)
        rows = SparseArrays.rowvals(A)
        @inbounds for c in 1:nc, ptr in SparseArrays.nzrange(A, c)
            if !iszero(vals[ptr])
                r = rows[ptr]
                packed[c][(r - 1) ÷ 64 + 1] |= UInt64(1) << ((r - 1) % 64)
            end
        end
    else
        @inbounds for c in 1:nc, r in 1:nr
            if !iszero(A[r, c])
                packed[c][(r - 1) ÷ 64 + 1] |= UInt64(1) << ((r - 1) % 64)
            end
        end
    end
    return packed
end

function _pack_binary_rows(A::_CSSBinaryMatrix)
    nr, nc = size(A)
    chunks = cld(nc, 64)
    packed = [zeros(UInt64, chunks) for _ in 1:nr]

    if A isa SparseMatrixCSC
        vals = SparseArrays.nonzeros(A)
        rows = SparseArrays.rowvals(A)
        @inbounds for c in 1:nc, ptr in SparseArrays.nzrange(A, c)
            iszero(vals[ptr]) && continue
            packed[rows[ptr]][(c - 1) ÷ 64 + 1] |=
                UInt64(1) << ((c - 1) % 64)
        end
    else
        @inbounds for r in 1:nr, c in 1:nc
            iszero(A[r, c]) && continue
            packed[r][(c - 1) ÷ 64 + 1] |=
                UInt64(1) << ((c - 1) % 64)
        end
    end
    return packed
end

@inline function _xor_packed!(target::Vector{UInt64}, source::Vector{UInt64})
    @inbounds @simd for i in eachindex(target, source)
        target[i] ⊻= source[i]
    end
    return target
end

@inline _packed_weight(v::Vector{UInt64}) =
    sum(count_ones, v; init=0)

function _unpack_binary_vector(v::Vector{UInt64}, n::Int)
    result = zeros(Int, n)
    @inbounds for c in 1:n
        result[c] = Int((v[(c - 1) ÷ 64 + 1] >> ((c - 1) % 64)) & UInt64(1))
    end
    return result
end

"""
    _minimum_distance_css_gray_binary(H, logical_checks)

Exactly compute the minimum weight of a binary vector `v` satisfying
`H * v' == 0` and `logical_checks * v' != 0` by Gray-code enumeration of
`ker(H)`. Returns `(distance, witness::Vector{Int})`.
"""
function _minimum_distance_css_gray_binary(
    H::_CSSBinaryMatrix, logical_checks::_CSSBinaryMatrix
)
    ncols(H) == ncols(logical_checks) ||
        throw(ArgumentError("The parity-check and logical-check matrices must have the same number of columns."))

    F = Oscar.Nemo.Native.GF(2)
    H_int = _convert_binary_to_int_matrix(H)
    L_int = _convert_binary_to_int_matrix(logical_checks)
    C = LinearCode(matrix(F, H_int), true)
    G = _convert_binary_to_int_matrix(generator_matrix(C))
    k_normalizer, n = size(G)

    iszero(k_normalizer) && return -1, zeros(Int, n)
    k_normalizer <= 62 || throw(ArgumentError(
        "Gray enumeration supports normalizer dimension at most 62; use `alg=:Wagner` for this code."))

    row_labels = (G * transpose(L_int)) .% 2
    packed_rows = _pack_binary_rows(G)
    packed_labels = _pack_binary_rows(row_labels)
    stop = UInt64(1) << k_normalizer
    candidate_count = stop - 1
    task_count = min(Threads.nthreads(), Int(candidate_count))
    local_best = fill((n + 1, typemax(UInt64)), task_count)

    Threads.@threads for task in 1:task_count
        first_i = UInt64(1) +
            (UInt64(task - 1) * candidate_count) ÷ UInt64(task_count)
        last_i =
            (UInt64(task) * candidate_count) ÷ UInt64(task_count)
        current = zeros(UInt64, cld(n, 64))
        current_label = zeros(UInt64, cld(size(L_int, 1), 64))
        gray = first_i ⊻ (first_i >> 1)

        bits = gray
        while !iszero(bits)
            row = trailing_zeros(bits) + 1
            _xor_packed!(current, packed_rows[row])
            _xor_packed!(current_label, packed_labels[row])
            bits &= bits - UInt64(1)
        end

        best_weight = n + 1
        best_index = typemax(UInt64)
        for i in first_i:last_i
            if any(!iszero, current_label)
                weight = _packed_weight(current)
                if weight < best_weight
                    best_weight = weight
                    best_index = i
                end
            end
            i == last_i && break
            next_gray = (i + UInt64(1)) ⊻ ((i + UInt64(1)) >> 1)
            changed = trailing_zeros(gray ⊻ next_gray) + 1
            _xor_packed!(current, packed_rows[changed])
            _xor_packed!(current_label, packed_labels[changed])
            gray = next_gray
        end
        local_best[task] = (best_weight, best_index)
    end

    best_weight, best_index = minimum(local_best)
    best_weight == n + 1 && return -1, zeros(Int, n)

    best = zeros(UInt64, cld(n, 64))
    bits = best_index ⊻ (best_index >> 1)
    while !iszero(bits)
        row = trailing_zeros(bits) + 1
        _xor_packed!(best, packed_rows[row])
        bits &= bits - UInt64(1)
    end
    return best_weight, _unpack_binary_vector(best, n)
end

"""
    _minimum_distance_css_wagner_binary(H, logical_checks; max_d=ncols(H), verbose=false)

Search in increasing physical weight for a vector with zero syndrome under `H`
and nonzero logical label under `logical_checks`. The meet-in-the-middle table
retains up to two distinct logical labels for each syndrome; this is sufficient
because a right label can forbid only one left label.
"""
function _minimum_distance_css_wagner_binary(
    H::_CSSBinaryMatrix, logical_checks::_CSSBinaryMatrix;
    max_d::Int=ncols(H), verbose::Bool=false
)
    ncols(H) == ncols(logical_checks) ||
        throw(ArgumentError("The parity-check and logical-check matrices must have the same number of columns."))
    0 <= max_d <= ncols(H) ||
        throw(DomainError(max_d, "`max_d` must lie between zero and the code length."))

    n = ncols(H)
    mid = n ÷ 2
    H_columns = _pack_binary_columns(H)
    L_columns = _pack_binary_columns(logical_checks)
    H_chunks = cld(nrows(H), 64)
    L_chunks = cld(nrows(logical_checks), 64)
    Syndrome = typeof(Tuple(zeros(UInt64, H_chunks)))
    Label = typeof(Tuple(zeros(UInt64, L_chunks)))
    Entry = Tuple{Label, Vector{Int}}

    for w in 1:max_d
        verbose && println("Checking CSS logical operators of weight $w...")
        splits = collect(max(0, w - (n - mid)):min(w, mid))
        witnesses = Vector{Union{Nothing, Vector{Int}}}(nothing, length(splits))
        found = Threads.Atomic{Bool}(false)

        Threads.@threads for split_index in eachindex(splits)
            found[] && continue
            w_left = splits[split_index]
            w_right = w - w_left
            table = Dict{Syndrome, Vector{Entry}}()
            left_syndrome = zeros(UInt64, H_chunks)
            left_label = zeros(UInt64, L_chunks)
            left_support = Int[]

            function build_left!(next_column::Int, remaining::Int)
                found[] && return
                if iszero(remaining)
                    syndrome = Tuple(left_syndrome)
                    label = Tuple(left_label)
                    bucket = get!(table, syndrome, Entry[])
                    if all(entry -> entry[1] != label, bucket) && length(bucket) < 2
                        push!(bucket, (label, copy(left_support)))
                    end
                    return
                end

                final_column = mid - remaining + 1
                for c in next_column:final_column
                    push!(left_support, c)
                    @inbounds @simd for i in 1:H_chunks
                        left_syndrome[i] ⊻= H_columns[c][i]
                    end
                    @inbounds @simd for i in 1:L_chunks
                        left_label[i] ⊻= L_columns[c][i]
                    end
                    build_left!(c + 1, remaining - 1)
                    @inbounds @simd for i in 1:H_chunks
                        left_syndrome[i] ⊻= H_columns[c][i]
                    end
                    @inbounds @simd for i in 1:L_chunks
                        left_label[i] ⊻= L_columns[c][i]
                    end
                    pop!(left_support)
                end
            end
            build_left!(1, w_left)

            right_syndrome = zeros(UInt64, H_chunks)
            right_label = zeros(UInt64, L_chunks)
            right_support = Int[]

            function probe_right!(next_column::Int, remaining::Int)
                found[] && return nothing
                if iszero(remaining)
                    syndrome = Tuple(right_syndrome)
                    haskey(table, syndrome) || return nothing
                    label = Tuple(right_label)
                    for (stored_label, stored_support) in table[syndrome]
                        if stored_label != label
                            witness = zeros(Int, n)
                            witness[stored_support] .= 1
                            witness[right_support] .= 1
                            return witness
                        end
                    end
                    return nothing
                end

                final_column = n - remaining + 1
                for c in next_column:final_column
                    push!(right_support, c)
                    @inbounds @simd for i in 1:H_chunks
                        right_syndrome[i] ⊻= H_columns[c][i]
                    end
                    @inbounds @simd for i in 1:L_chunks
                        right_label[i] ⊻= L_columns[c][i]
                    end
                    witness = probe_right!(c + 1, remaining - 1)
                    witness === nothing || return witness
                    @inbounds @simd for i in 1:H_chunks
                        right_syndrome[i] ⊻= H_columns[c][i]
                    end
                    @inbounds @simd for i in 1:L_chunks
                        right_label[i] ⊻= L_columns[c][i]
                    end
                    pop!(right_support)
                end
                return nothing
            end

            witness = probe_right!(mid + 1, w_right)
            if witness !== nothing
                witnesses[split_index] = witness
                found[] = true
            end
        end
        for witness in witnesses
            witness === nothing || return w, witness
        end
    end
    return -1, zeros(Int, n)
end

"""
    _minimum_distance(H, logical_checks; alg=:auto, max_d=ncols(H), verbose=false)

Compute the minimum weight in `ker(H)` having a nonzero logical label. Both
matrices are binary; `H` may be a sparse parity-check matrix.
"""
function _minimum_distance(
    H::_CSSBinaryMatrix, logical_checks::_CSSBinaryMatrix;
    alg::Symbol=:auto, max_d::Int=ncols(H), verbose::Bool=false,
    time_limit_sec::Union{Nothing, Float64}=nothing,
    ilp_parity_cut_max_degree::Int=10
)
    alg ∈ (:auto, :Gray, :Wagner, :ILP) ||
        throw(ArgumentError("CSS matrix minimum distance supports `:auto`, `:Gray`, `:Wagner`, and `:ILP`."))
    0 <= max_d <= ncols(H) ||
        throw(DomainError(max_d, "`max_d` must lie between zero and the code length."))

    ilp_available = !isempty(methods(_minimum_distance_css_ILP))
    chosen = alg
    if chosen == :auto
        normalizer_dimension = if H isa SparseMatrixCSC
            ncols(H) - _binary_sparse_rank(H)
        else
            H_dense = matrix(Oscar.Nemo.Native.GF(2), _convert_binary_to_int_matrix(H))
            ncols(H) - rank(H_dense)
        end
        if normalizer_dimension <= 24 && ncols(H) <= 256
            chosen = :Gray
        elseif ilp_available && (H isa SparseMatrixCSC || ncols(H) >= 96)
            chosen = :ILP
        else
            chosen = :Wagner
        end
    end

    if chosen == :Gray
        d, witness = _minimum_distance_css_gray_binary(H, logical_checks)
        return d > max_d ? (-1, zeros(Int, ncols(H))) : (d, witness)
    elseif chosen == :ILP
        ilp_available || throw(ArgumentError(
            "Load JuMP and GLPK before selecting the CSS `:ILP` solver."))
        return _minimum_distance_css_ILP(H, logical_checks;
            max_d=max_d, verbose=verbose, time_limit_sec=time_limit_sec,
            parity_cut_max_degree=ilp_parity_cut_max_degree)
    end
    return _minimum_distance_css_wagner_binary(
        H, logical_checks; max_d=max_d, verbose=verbose)
end

function _css_distance_problem(S::AbstractStabilizerCodeCSS, which::Symbol)
    logs = logicals_matrix(S)
    n = S.n
    if which == :X
        return Z_stabilizers(S), _remove_empty(logs[:, n + 1:2n], :rows)
    elseif which == :Z
        return X_stabilizers(S), _remove_empty(logs[:, 1:n], :rows)
    end
    throw(ArgumentError("Expected `which` to be `:X` or `:Z`."))
end

function _css_quantum_witness(S::AbstractStabilizerCodeCSS, which::Symbol, v::Vector{Int})
    z = zeros(Int, S.n)
    data = which == :X ? vcat(v, z) : vcat(z, v)
    return matrix(S.F, 1, 2 * S.n, data)
end

function _css_cached_distance_result(S::AbstractStabilizerCodeCSS, which::Symbol)
    distance_key = which == :X ? :dx : which == :Z ? :dz : :d
    witness_key = which == :X ? :X_minimum_distance_witness :
        which == :Z ? :Z_minimum_distance_witness : :minimum_distance_witness
    return S.cache[distance_key], get(
        S.cache, witness_key, zero_matrix(S.F, 1, 2 * S.n))
end

function _compute_css_distance(
    S::AbstractStabilizerCodeCSS, which::Symbol;
    alg::Symbol=:auto, max_d::Int=S.n, verbose::Bool=false,
    time_limit_sec::Union{Nothing, Float64}=nothing,
    ilp_parity_cut_max_degree::Int=10
)
    distance_key = which == :X ? :dx : :dz
    haskey(S.cache, distance_key) && return _css_cached_distance_result(S, which)

    H, logical_checks = _css_distance_problem(S, which)
    d, v = _minimum_distance(H, logical_checks;
        alg=alg, max_d=max_d, verbose=verbose,
        time_limit_sec=time_limit_sec,
        ilp_parity_cut_max_degree=ilp_parity_cut_max_degree)
    d == -1 && return d, zero_matrix(S.F, 1, 2 * S.n)

    witness = _css_quantum_witness(S, which, v)
    if which == :X
        set_X_minimum_distance!(S, d)
        S.cache[:X_minimum_distance_witness] = witness
    else
        set_Z_minimum_distance!(S, d)
        S.cache[:Z_minimum_distance_witness] = witness
    end
    return d, witness
end

"""
    minimum_distance(S::AbstractStabilizerCodeCSS; which=:full, alg=:auto,
                     max_d=S.n, verbose=false)

Compute the `X`, `Z`, or full minimum distance of a binary CSS stabilizer code.
The witness is returned in symplectic `[X | Z]` form.
"""
function minimum_distance(
    S::AbstractStabilizerCodeCSS; which::Symbol=:full, alg::Symbol=:auto,
    max_d::Int=S.n, verbose::Bool=false,
    time_limit_sec::Union{Nothing, Float64}=nothing,
    ilp_parity_cut_max_degree::Int=10
)
    which ∈ (:full, :X, :Z) ||
        throw(ArgumentError("Expected `which` to be `:full`, `:X`, or `:Z`."))
    Int(order(field(S))) == 2 ||
        throw(ArgumentError("CSS minimum-distance solvers are currently implemented only over GF(2)."))
    S.k > 0 || throw(ArgumentError(
        "Logical minimum distance is undefined for a CSS stabilizer state with k = 0."))

    if which != :full
        return _compute_css_distance(
            S, which; alg=alg, max_d=max_d, verbose=verbose,
            time_limit_sec=time_limit_sec,
            ilp_parity_cut_max_degree=ilp_parity_cut_max_degree)
    end
    haskey(S.cache, :d) && return _css_cached_distance_result(S, :full)

    dx, wx = _compute_css_distance(
        S, :X; alg=alg, max_d=max_d, verbose=verbose,
        time_limit_sec=time_limit_sec,
        ilp_parity_cut_max_degree=ilp_parity_cut_max_degree)
    z_limit = dx == -1 ? max_d : min(max_d, dx)
    dz, wz = _compute_css_distance(
        S, :Z; alg=alg, max_d=z_limit, verbose=verbose,
        time_limit_sec=time_limit_sec,
        ilp_parity_cut_max_degree=ilp_parity_cut_max_degree)
    (dx == -1 && dz == -1) &&
        return -1, zero_matrix(S.F, 1, 2 * S.n)

    if dz == -1 || (dx != -1 && dx <= dz)
        set_minimum_distance!(S, dx)
        S.cache[:minimum_distance_witness] = wx
        return dx, wx
    end
    set_minimum_distance!(S, dz)
    S.cache[:minimum_distance_witness] = wz
    return dz, wz
end

X_minimum_distance(S::AbstractStabilizerCodeCSS; kwargs...) =
    minimum_distance(S; which=:X, kwargs...)
Z_minimum_distance(S::AbstractStabilizerCodeCSS; kwargs...) =
    minimum_distance(S; which=:Z, kwargs...)
XZ_minimum_distance(S::AbstractStabilizerCodeCSS; kwargs...) =
    (X_minimum_distance(S; kwargs...)[1], Z_minimum_distance(S; kwargs...)[1])
