# Copyright (c) 2024 - 2025 Benjamin Ide, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

function _thickened_cone(S::AbstractStabilizerCodeCSS, HX_A::CTMatrixTypes, HZ_A::CTMatrixTypes,
        f1::CTMatrixTypes, f0::CTMatrixTypes, type::Symbol, r::Int = 1)

    r > 0 || (r == 0 && (return S;)) || throw(DomainError(r, "Must be a positive integer."))
    type in (:X, :Z) || throw(DomainError(type, "Must choose `type` to be `:X` or `:Z`."))
    
    HX_C = type == :X ? X_stabilizers(S) : Z_stabilizers(S)
    HZ_C = type == :X ? Z_stabilizers(S) : X_stabilizers(S)
    F = field(S)
    HX_C = _dense_code_matrix(HX_C, F)
    HZ_C = _dense_code_matrix(HZ_C, F)
    HX_A = _dense_code_matrix(HX_A, F)
    HZ_A = _dense_code_matrix(HZ_A, F)
    f1 = _dense_code_matrix(f1, F)
    f0 = _dense_code_matrix(f0, F)
    all(==(F), (base_ring(HX_A), base_ring(HZ_A), base_ring(f1), base_ring(f0))) ||
        throw(ArgumentError("All cone data must use the code field."))
    
    nC = length(S)
    nA = size(HX_A, 2)
    size(HZ_A, 2) == nA ||
        throw(ArgumentError("The auxiliary X and Z checks must have the same width."))
    size(f1) == (nC, size(HX_A, 1)) ||
        throw(ArgumentError("f1 has incompatible dimensions."))
    size(f0) == (size(HZ_C, 1), nA) ||
        throw(ArgumentError("f0 has incompatible dimensions."))
    HZ_C * f1 == f0 * transpose(HX_A) ||
        throw(ArgumentError("The cone maps do not satisfy the chain-map equation."))
    n = nC + r * nA + (r - 1) * size(HX_A, 1)
    
    I_r = identity_matrix(F, r)
    I_rm1 = identity_matrix(F, r - 1)
    I_nA = identity_matrix(F, nA)
    I_szHXA = identity_matrix(F, size(HX_A, 1))
    rep_tr = _rep_pcm_tr(F, r)
    rep_pcm = _rep_pcm(F, r)
    Z_blocks1 = zero_matrix(F, size(HX_C, 1), n - nC)
    Z_blocks2 = zero_matrix(F, (r - 1) * size(HX_A, 1), nC)
    
    HX = vcat(hcat(HX_C, Z_blocks1),
              hcat(vcat(transpose(f1), Z_blocks2), kronecker_product(I_r, HX_A), kronecker_product(rep_tr, I_szHXA)))
              
    # ... Similar Sparse-Aware logic for HZ ...
    HZ = vcat(hcat(HZ_C, f0, zero_matrix(F, size(HZ_C, 1), (r - 1) * (nA + size(HX_A, 1)))),
              hcat(zero_matrix(F, (r - 1) * nA, nC), kronecker_product(rep_pcm, I_nA), kronecker_product(I_rm1, transpose(HX_A))), 
              hcat(zero_matrix(F, size(HZ_A, 1), nC + (r - 1) * nA), HZ_A, zero_matrix(F, size(HZ_A, 1), (r - 1) * size(HX_A, 1))))
              
    stabs = type == :X ? direct_sum(HX, HZ) : direct_sum(HZ, HX)

    implied_stabs = zero_matrix(F, 0, nC)
    N = kernel(transpose(HX_A), side = :right)
    for i in axes(N, 2)
        implied_stabs = vcat(implied_stabs, transpose(f1 * N[:, i:i]))
    end
    temp = CSSCode(vcat(HX_C, implied_stabs), HZ_C)
    X_logs = dimension(temp) > 0 ?
        _remove_empty(logicals_matrix(temp)[:, 1:nC], :rows) :
        zero_matrix(F, 0, nC)
    X_logs = hcat(X_logs, zero_matrix(F, size(X_logs, 1), n - nC))

    implied_stabs_Z = zero_matrix(F, 0, nA)
    N_Z = transpose(kernel(transpose(HZ_C), side = :right))
    for i in axes(N_Z, 1)
        implied_stabs_Z = vcat(implied_stabs_Z, N_Z[i:i, :] * f0)
    end
    Z_gauges = if size(HZ_A, 1) + size(implied_stabs_Z, 1) == 0
        zero_matrix(F, 0, nA)
    else
        temp = CSSCode(HX_A, vcat(HZ_A, implied_stabs_Z))
        dimension(temp) > 0 ?
            _remove_empty(logicals_matrix(temp)[:, nA + 1:2nA], :rows) :
            zero_matrix(F, 0, nA)
    end
    Z_gauges = hcat(zero_matrix(F, size(Z_gauges, 1), nC + (r - 1) * nA), Z_gauges, zero_matrix(F, size(Z_gauges, 1), (r - 1) * size(HX_A, 1)))

    new_X, new_Z, new_mixed, _ = _complete_pairs(stabs, type == :X ? direct_sum(X_logs, Z_gauges) : direct_sum(Z_gauges, X_logs))
    isempty(new_mixed) ||
        error("The cone construction produced mixed logical partners.")
    logs = type == :X ?
        direct_sum(X_logs, new_Z[:, n + 1:2n]) :
        direct_sum(new_X[:, 1:n], X_logs)
    gauges = type == :X ? direct_sum(new_X[:, 1:n], Z_gauges) : direct_sum(Z_gauges, new_Z[:, n + 1:2n])

    return isempty(gauges) ? StabilizerCode(stabs) :
        SubsystemCode(stabs, logs, gauges)
end
_thickened_cone(S::AbstractStabilizerCodeCSS, A::AbstractStabilizerCodeCSS, f1::CTMatrixTypes,
    f0::CTMatrixTypes, type::Symbol, r::Int = 1) = _thickened_cone(
    S, X_stabilizers(A), Z_stabilizers(A), f1, f0, type, r)

"""
$(TYPEDSIGNATURES)

Return the (mapping) cone code associated with measuring the logical(s) `L` of the CSS stabilizer
code `S`.

# Optional Arguments
All paramaters are aligned with their respective papers.
- `style` - `:Xanadu`, `:IBM`, or `:Cohen`
- `max_iters` - used for `:Xanadu` and `:IBM`
- `r` - used for `:Cohen`
- `improve_cycles` - used for `IBM`
- `remove_and_improve_cycles` - used for `IBM`, supersedes previous parameter
- `log_checking` - set to false to intentionally use an `L` that isn't a logical
"""
function homological_measurement(S::AbstractStabilizerCodeCSS, L::CTMatrixTypes; style::Symbol =
    :Xanadu, r::Int = 1, max_iters::Int = 50000, cellulate::Bool = false, improve_cycles::Bool =
    true, remove_and_improve_cycles::Bool = false, log_checking::Bool = true,
    rng::AbstractRNG = Random.default_rng())

    is_positive(r) || throw(DomainError(r, "Must be a positive integer."))
    is_positive(max_iters) || throw(DomainError(max_iters, "Must be a positive integer."))
    cellulate &&
        throw(ArgumentError("The `cellulate` option is not implemented."))
    L_red = _remove_empty(L, :rows)
    nrows(L_red) == 1 || throw(ArgumentError("Requires a single logical of the code."))
    is_logical(S, L_red) || !log_checking || throw(ArgumentError("The input matrix is not a logical of the code."))
    F = field(S)
    Fone = F(1)
    Int(order(F)) == 2 || throw(ArgumentError("Only defined for binary codes."))
    n = length(S)
    # k = dimension(S)

    type, stabs, log = if iszero(L_red[1:1, 1 + n:2n])
        :X, Z_stabilizers(S), L_red[1:1, 1:n]
    elseif iszero(L_red[1:1, 1:n])
        :Z, X_stabilizers(S), L_red[1:1, 1 + n:2n]
    else
        throw(DomainError(L, "Only defined for pure X or Z logicals `L`."))
    end

    Q = getindex.(findall(!iszero, log), 2)
    f1 = matrix(F, Int[Q[j] == i for i in 1:n, j in 1:length(Q)])
    nonzero = findall(!iszero(stabs[i, Q]) for i in 1:size(stabs, 1))
    f0 = zero_matrix(F, size(stabs, 1), length(nonzero))
    for (i, idx) in enumerate(nonzero)
        f0[idx, i] = Fone
    end
    HX = transpose(stabs[nonzero, Q])

    if style == :Xanadu
        temp = size(HX, 2)
        HX = matrix(F, _add_edges(_Flint_matrix_to_Julia_int_matrix(HX); rng=rng))
        f0 = hcat(f0, zero_matrix(F, size(f0, 1), size(HX, 2) - temp))
        HZ = _remove_empty(rref(transpose(kernel(HX, side = :right)))[2], :rows)
        a = transpose(kernel(transpose(stabs), side = :right))
        b = _remove_empty(rref(a * f0)[2], :rows)
        if isempty(b)
            HZ = _find_low_weights_rand(HZ, max_iters; rng=rng)
        else
            HZ = _find_low_weight_cycle_subspace(HZ, b, max_iters; rng=rng)
        end
        return _thickened_cone(S, HX, HZ, f1, f0, type)
    elseif style == :IBM
        HZ = _remove_empty(rref(transpose(kernel(HX, side = :right)))[2], :rows)

        if remove_and_improve_cycles
            a = transpose(kernel(transpose(stabs), side = :right))
            b = _remove_empty(rref(a * f0)[2], :rows)
            if isempty(b)
                HZ = _find_low_weights_rand(HZ, max_iters; rng=rng)
            else
                HZ = _find_low_weight_cycle_subspace(HZ, b, max_iters; rng=rng)
            end
        elseif improve_cycles
            HZ = _find_low_weights_rand(HZ, max_iters; rng=rng)
        end

        r = ceil(Int, 1 / Cheeger_constant(_Flint_matrix_to_Julia_int_matrix(HX)))
        return _thickened_cone(S, HX, HZ, f1, f0, type, r)
    elseif style == :Cohen
        HZ = zero_matrix(F, 0, size(HX, 2))
        return _thickened_cone(S, HX, HZ, f1, f0, type, r)
    else
        throw(ArgumentError("Unknown `style` parameter $style"))
    end
end

#############################
     # general functions
#############################

_rep_pcm_tr(F::CTFieldTypes, d::Int) = matrix(F, diagm(d, d - 1, 0 => ones(Int, d - 1), -1 =>
    ones(Int, d - 1)))
_rep_pcm(F::CTFieldTypes, d::Int) = matrix(F, diagm(d - 1, d, 0 => ones(Int, d - 1), 1 => ones(Int,
    d - 1)))

function _complete_pairs(stabs::CTMatrixTypes, logs::CTMatrixTypes)
    F = base_ring(stabs)
    @assert iseven(size(stabs, 2))
    @assert size(stabs, 2) == size(logs, 2)
    @assert rank(logs) == size(logs, 1)
    n = div(size(stabs, 2), 2)
    
    # 1. Fast symplectic dot product without allocating Ω
    X_logs = logs[:, 1:n]
    Z_logs = logs[:, n+1:end]
    sign_factor = Int(characteristic(F)) == 2 ? 1 : -1
    find_pairs = X_logs * transpose(Z_logs) + sign_factor * Z_logs * transpose(X_logs)
    
    @assert all(count(!iszero, find_pairs[i:i, :]) in (0, 1) for i in 1:size(find_pairs, 1))
    needs_pair = findall(iszero(find_pairs[i:i, :]) for i in 1:size(find_pairs, 1))

    new_X = zero_matrix(F, 0, 2n)
    new_Z = zero_matrix(F, 0, 2n)
    new_mixed = zero_matrix(F, 0, 2n)

    # 2. Block-wise LHS construction to avoid Ω multiplication
    X_stabs = stabs[:, 1:n]
    Z_stabs = stabs[:, n+1:end]
    
    for i in needs_pair
        RHS = zero_matrix(F, size(logs, 1) + size(stabs, 1) + n, 1)
        RHS[i, 1] = 1

        # Force the X block to zero, producing a pure Z partner.
        LHS_X_sym = vcat(Z_logs, Z_stabs, identity_matrix(F, n))
        LHS_Z_sym = vcat(sign_factor * X_logs, sign_factor * X_stabs, zero_matrix(F, n, n))
        LHS = hcat(LHS_X_sym, LHS_Z_sym)
        
        flag, sol = can_solve_with_solution(LHS, RHS, side = :right)
        if flag
            logs = vcat(logs, transpose(sol))
            new_Z = vcat(new_Z, transpose(sol))
            continue
        end

        # Force the Z block to zero, producing a pure X partner.
        LHS_X_sym[size(logs, 1) + size(stabs, 1) + 1:end, :] = zero_matrix(F, n, n)
        LHS_Z_sym[size(logs, 1) + size(stabs, 1) + 1:end, :] = identity_matrix(F, n)
        LHS = hcat(LHS_X_sym, LHS_Z_sym)
        
        flag, sol = can_solve_with_solution(LHS, RHS, side = :right)
        if flag
            logs = vcat(logs, transpose(sol))
            new_X = vcat(new_X, transpose(sol))
            continue
        end

        # try mixed
        LHS_Z_sym[size(logs, 1) + size(stabs, 1) + 1:end, :] = zero_matrix(F, n, n)
        LHS = hcat(LHS_X_sym, LHS_Z_sym)
        
        flag, sol = can_solve_with_solution(LHS, RHS, side = :right)
        if flag
            logs = vcat(logs, transpose(sol))
            new_mixed = vcat(new_mixed, transpose(sol))
            continue
        end
    end

    return new_X, new_Z, new_mixed, logs
end

"""
$(TYPEDSIGNATURES)

Return the Cheeger constant of the matrix `M` assuming `M` is a vertex-edge incidence matrix.
"""
function Cheeger_constant(M::Matrix{T}) where T <: Integer
    m, n = size(M)
    2 <= m < 64 ||
        throw(ArgumentError("The exact Cheeger constant requires 2 to 63 vertices."))
    edge_masks = UInt64[
        sum(UInt64(1) << (i - 1) for i in 1:m if isodd(M[i, j]))
        for j in 1:n
    ]
    h = Inf
    for mask in UInt64(1):((UInt64(1) << (m - 1)) - 1)
        subset_size = count_ones(mask)
        denominator = min(subset_size, m - subset_size)
        boundary = count(isodd(count_ones(mask & edge)) for edge in edge_masks)
        h = min(h, boundary / denominator)
    end
    return h
end

"""
    Cheeger_constant(S::AbstractSubsystemCodeCSS, L::CTMatrixTypes)

Return the Cheeger constant of the incidence graph induced by the pure `X` or
`Z` logical operator `L` and the opposite-type stabilizers of `S`.
"""
function Cheeger_constant(S::AbstractSubsystemCodeCSS, L::CTMatrixTypes)
    L_red = _remove_empty(L, :rows)
    nrows(L_red) == 1 || throw(ArgumentError("Requires a single logical of the code."))
    is_logical(S, L_red) || throw(ArgumentError("The input matrix is not a logical of the code."))

    n = length(S)
    stabs, log = if iszero(L_red[1:1, 1 + n:2n])
        Z_stabilizers(S), L_red[1:1, 1:n]
    elseif iszero(L_red[1:1, 1:n])
        X_stabilizers(S), L_red[1:1, 1 + n:2n]
    else
        throw(DomainError(L_red, "Only defined for pure `X` or `Z` logicals `L`."))
    end
    Q = getindex.(findall(!iszero, log), 2)
    nonzero = findall(!iszero(stabs[i:i, Q]) for i in 1:size(stabs, 1))
    graph = transpose(stabs[nonzero, Q])
    return Cheeger_constant(_Flint_matrix_to_Julia_int_matrix(graph))
end

function _sparsest_cut(M::AbstractMatrix{T}; rng::AbstractRNG = Random.default_rng()) where T <: Integer
    m, n = size(M)
    2 <= m < 64 ||
        throw(ArgumentError("The exact sparsest cut requires 2 to 63 vertices."))
    edge_masks = UInt64[
        sum(UInt64(1) << (i - 1) for i in 1:m if isodd(M[i, j]))
        for j in 1:n
    ]
    h = Inf
    sparse_cut = zeros(T, m)
    ties = 0
    for mask in UInt64(1):((UInt64(1) << (m - 1)) - 1)
        subset_size = count_ones(mask)
        denominator = min(subset_size, m - subset_size)
        boundary = count(isodd(count_ones(mask & edge)) for edge in edge_masks)
        temp = boundary / denominator
        if temp < h
            h = temp
            ties = 1
            for i in 1:m
                sparse_cut[i] = T((mask >> (i - 1)) & 1)
            end
        elseif temp == h
            ties += 1
            if rand(rng, 1:ties) == 1
                for i in 1:m
                    sparse_cut[i] = T((mask >> (i - 1)) & 1)
                end
            end
        end
    end

    return h, sparse_cut
end

function _add_edges(M::Matrix{T}; rng = Xoshiro()) where T <: Integer
    m, n_orig = size(M)
    
    # Pre-allocate a buffer to avoid hcat inside the loop
    capacity = n_orig + 50
    M_new = zeros(T, m, capacity)
    M_new[:, 1:n_orig] .= M
    
    M_temp = zeros(T, m, capacity)
    curr_cols = n_orig
    
    Tzero = T(0)
    Tone = T(1)
    
    while true
        # Operate only on the active view
        active_view = view(M_new, :, 1:curr_cols)
        h, v = _sparsest_cut(active_view, rng = rng)
        h < 1 || break
        
        # Expand capacity if we hit the limit
        if curr_cols == capacity
            capacity += 50
            M_new_expand = zeros(T, m, capacity)
            M_new_expand[:, 1:curr_cols] .= M_new[:, 1:curr_cols]
            M_new = M_new_expand
            M_temp = zeros(T, m, capacity)
        end
        
        curr_cols += 1
        M_sum = vec(sum(view(M_new, :, 1:curr_cols-1), dims = 2))
        
        S_idx = findall(!iszero, v)
        Sc_idx = findall(iszero, v)
        S_min = minimum(M_sum[S_idx])
        Sc_min = minimum(M_sum[Sc_idx])
        
        is = S_idx[findall(x -> M_sum[x] == S_min, S_idx)]
        js = Sc_idx[findall(x -> M_sum[x] == Sc_min, Sc_idx)]
        
        h_best = -Inf
        best_i, best_j = -1, -1
        
        # Copy current state to temp buffer exactly once per expansion
        M_temp[:, 1:curr_cols] .= M_new[:, 1:curr_cols]
        
        for (i, j) in Iterators.product(shuffle(rng, is), shuffle(rng, js))
            M_temp[i, curr_cols] = Tone
            M_temp[j, curr_cols] = Tone
            
            temp_h = Cheeger_constant(M_temp[:, 1:curr_cols])
            if temp_h > h_best
                h_best = temp_h
                best_i, best_j = i, j
            end
            
            M_temp[i, curr_cols] = Tzero
            M_temp[j, curr_cols] = Tzero
        end
        
        M_new[best_i, curr_cols] = Tone
        M_new[best_j, curr_cols] = Tone
    end
    
    return M_new[:, 1:curr_cols]
end

function _random_matrix!(A::CTMatrixTypes, F::CTFieldTypes;
    rng::AbstractRNG = Random.default_rng())
    # Operates in-place on pre-allocated A to avoid GC overhead
    nr, nc = size(A)
    for i in 1:nr
        for j in 1:nc
            A[i, j] = rand(rng, F)
        end
    end
    return A
end

function _find_low_weight_cycle_subspace(all_cycles::CTMatrixTypes,
    already_covered_cycles::CTMatrixTypes, max_iters::Int, f::T = maximum;
    rng::AbstractRNG = Random.default_rng()) where T <: Function

    @assert size(all_cycles, 2) == size(already_covered_cycles, 2)
    F = base_ring(all_cycles)
    
    A = _remove_empty(rref(already_covered_cycles)[2], :rows)
    pivots = findall(count(!iszero(A[i, j]) for i in axes(A, 1)) == 1 for j in axes(A, 2))
    B = _remove_empty(rref(all_cycles)[2], :rows)
    
    for i in axes(B, 1)
        for j in axes(A, 1)
            if !iszero(B[i, pivots[j]])
                x_inv = inv(B[i, pivots[j]])
                for k in axes(B, 2)
                    B[i, k] = x_inv * A[j, k] + B[i, k]
                end
            end
        end
    end
    B = _remove_empty(rref(B)[2], :rows)
    isempty(B) && return B

    w = size(A, 2) + 1
    C = zero_matrix(F, size(B, 1), size(B, 2))
    
    # Pre-allocations for the hot loop
    n_B = size(B, 1)
    n_A = size(A, 1)
    x = zero_matrix(F, n_B, n_B)
    y = zero_matrix(F, n_B, n_A)
    
    for i in 1:max_iters
        _random_invertible_matrix!(x, F, n_B; rng=rng)
        _random_matrix!(y, F; rng=rng)
        
        # Branch 1: temp = x * B + y * A
        temp = x * B + y * A
        temp_w = f(count(!iszero, temp[row, :]) for row in 1:size(temp, 1))
        if temp_w < w
            w = temp_w
            C = temp
        end
        
        # Branch 2: temp = x * B
        temp = x * B
        temp_w = f(count(!iszero, temp[row, :]) for row in 1:size(temp, 1))
        if temp_w < w
            w = temp_w
            C = temp
        end
    end
    
    return C
end

function _find_low_weights_rand(M::CTMatrixTypes, max_iters::Int, f::T = maximum;
    rng::AbstractRNG = Random.default_rng()) where T <: Function
    isempty(M) && return M
    initial_w = f(count(!iszero, M[i, :]) for i in 1:size(M, 1))
    A = _remove_empty(rref(M)[2], :rows)
    isempty(A) && return A
    
    F = base_ring(A)
    w = size(A, 2) + 1
    n_A = size(A, 1)
    
    # Preallocate once
    x = zero_matrix(F, n_A, n_A)
    
    for i in 1:max_iters
        _random_invertible_matrix!(x, F, n_A; rng=rng)
        temp = x * A
        
        # Fast weight check
        temp_w = f(count(!iszero, temp[j, :]) for j in 1:n_A)
        if temp_w < w
            w = temp_w
            A = temp
        end
    end
    return initial_w < w ? M : A
end

function _random_invertible_matrix!(A::CTMatrixTypes, F::CTFieldTypes, n::Int;
    rng::AbstractRNG = Random.default_rng())
    # Operates in-place on pre-allocated A to avoid GC overhead
    for i in 1:n
        for j in 1:n
            A[i, j] = F(0)
        end
    end
    
    # Generate a random lower triangular matrix with 1s on diagonal (always invertible)
    for i in 1:n
        A[i, i] = F(1)
        for j in 1:i-1
            A[i, j] = rand(rng, F)
        end
    end
    
    # Randomly permute rows to spread the entropy
    for i in n:-1:2
        swap_idx = rand(rng, 1:i)
        if swap_idx != i
            # Swap rows i and swap_idx
            for j in 1:n
                temp = A[i, j]
                A[i, j] = A[swap_idx, j]
                A[swap_idx, j] = temp
            end
        end
    end
    return A
end
