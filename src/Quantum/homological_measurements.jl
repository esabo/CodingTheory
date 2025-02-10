# Copyright (c) 2024 - 2025 Benjamin Ide, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

# TODO Check logicals are being used correctly throughout

function _thickened_cone(S::AbstractStabilizerCodeCSS, HX_A::CTMatrixTypes, HZ_A::CTMatrixTypes,
        f1::CTMatrixTypes, f0::CTMatrixTypes, type::Symbol, r::Int = 1)

    r > 0 || (r == 0 && (return S;)) || throw(DomainError(r, "Must be a positive integer."))
    type in (:X, :Z) || throw(DomainError(type, "Must choose `type` to be `:X` or `:Z`."))
    HX_C = type == :X ? X_stabilizers(S) : Z_stabilizers(S)
    HZ_C = type == :X ? Z_stabilizers(S) : X_stabilizers(S)
    # HX_A = type == :X ? X_stabilizers(A) : Z_stabilizers(A)
    # HZ_A = type == :X ? Z_stabilizers(A) : X_stabilizers(A)
    F = field(S)
    @assert F == base_ring(HX_A) == base_ring(HZ_A) == base_ring(f1) == base_ring(f0)
    @assert size(HX_A, 2) == size(HZ_A, 2)
    nC = length(S)
    nA = size(HX_A, 2)
    n = nC + r * nA + (r - 1) * size(HX_A, 1)
    @assert size(f1) == (nC, size(HX_A, 1))
    @assert size(f0) == (size(HZ_C, 1), nA)
    @assert HZ_C * f1 == f0 * transpose(HX_A)

    # build stabilizers
    HX = vcat(hcat(HX_C, zero_matrix(F, size(HX_C, 1), n - nC)),
        hcat(vcat(transpose(f1), zero_matrix(F, (r - 1) * size(HX_A, 1), nC)),
        identity_matrix(F, r) ⊗ HX_A, _rep_pcm_tr(F, r) ⊗ identity_matrix(F, size(HX_A, 1))))
    HZ = vcat(hcat(HZ_C, f0, zero_matrix(F, size(HZ_C, 1), (r - 1) * (nA + size(HX_A, 1)))),
        hcat(zero_matrix(F, (r - 1) * nA, nC), _rep_pcm(F, r) ⊗ identity_matrix(F, nA),
        identity_matrix(F, r - 1) ⊗ transpose(HX_A)), hcat(zero_matrix(F, size(HZ_A, 1), nC +
        (r - 1) * nA), HZ_A, zero_matrix(F, size(HZ_A, 1), (r - 1) * size(HX_A, 1))))
    stabs = type == :X ? direct_sum(HX, HZ) : direct_sum(HZ, HX)

    # build X logicals
    implied_stabs = zero_matrix(F, 0, nC)
    N = kernel(transpose(HX_A), side = :right)
    for i in axes(N, 2)
        implied_stabs = vcat(implied_stabs, transpose(f1 * N[:, i:i]))
    end
    temp = CSSCode(vcat(HX_C, implied_stabs), HZ_C)
    X_logs = if dimension(temp) > 0
        # TODO could reduce allocations by preallocating and then doing row = logicals(temp)[i][1][1:nc]
        _remove_empty(logicals_matrix(temp)[:, 1:nC], :rows)
    else
        zero_matrix(F, 0, nC)
    end
    # this can have 0 rows
    X_logs = hcat(X_logs, zero_matrix(F, size(X_logs, 1), n - nC))

    # build Z gauges
    implied_stabs = zero_matrix(F, 0, nA)
    N = transpose(kernel(transpose(HZ_C), side = :right))
    for i in axes(N, 1)
        implied_stabs = vcat(implied_stabs, N[i:i, :] * f0)
    end
    Z_gauges = if size(HZ_A, 1) + size(implied_stabs, 1) == 0
        zero_matrix(F, 0, nA)
    else
        temp = CSSCode(HX_A, vcat(HZ_A, implied_stabs))
        if dimension(temp) > 0
            # TODO can preallocate as above
            _remove_empty(logicals_matrix(temp)[:, 1 + nA:2nA], :rows)
        else
            zero_matrix(F, 0, nA)
        end
    end
    # Z_gauges = hcat(zero_matrix(F, size(Z_gauges, 1), nC), Z_gauges)
    Z_gauges = hcat(zero_matrix(F, size(Z_gauges, 1), nC + (r - 1) * nA), Z_gauges,
        zero_matrix(F, size(Z_gauges, 1), (r - 1) * size(HX_A, 1)))

    # build full logs and gauges
    new_X, new_Z, new_mixed, _ = _complete_pairs(stabs, type == :X ? direct_sum(X_logs, Z_gauges) :
        direct_sum(Z_gauges, X_logs))
    isempty(new_mixed) || error() # TODO
    logs = type == :X ? direct_sum(X_logs, new_Z[:, n + 1:2n]) : direct_sum(new_X[:, 1:n], X_logs)
    gauges = type == :X ? direct_sum(new_X[:, 1:n], Z_gauges) : direct_sum(Z_gauges, new_Z[:, n +
        1:2n])

    # TODO why are logs not also passed in here?
    # TODO can we remove the type instability here?
    return isempty(gauges) ? StabilizerCode(stabs) : SubsystemCode(stabs, logs, gauges)
end
_thickened_cone(S::AbstractStabilizerCodeCSS, A::AbstractStabilizerCodeCSS, f1::CTMatrixTypes,
    f0::CTMatrixTypes, type::Symbol, r::Int = 1) = _thickened_cone(S, A.X_stabs, A.Z_stabs, f1, f0,
    type, r)

"""
$TYPEDSIGNATURES

Return the (mapping) cone code associated with measuring the logical(s) `L` of the CSS stabilizer
code `S`.

# Optional Arguments
All paramaters are aligned with their respective papers.
- `style` - `:Xanadu`, `:IBM`, or `:Cohen`
- `max_iters` - used for `:Xanadu` and `:IBM`
- `r` - used for `:Cohen`
- `improve_cycles` - used for `IBM`
- `remove_and_improve_cycles` - used for `IBM`, supersedes previous parameter
"""
function homological_measurement(S::AbstractStabilizerCodeCSS, L::CTMatrixTypes; style::Symbol =
    :Xanadu, r::Int = 1, max_iters::Int = 50000, cellulate::Bool = false, improve_cycles::Bool =
    true, remove_and_improve_cycles::Bool = false)

    is_positive(r) || throw(DomainError(r, "Must be a positive integer."))
    is_positive(max_iters) || throw(DomainError(max_iters, "Must be a positive integer."))
    L_red = _remove_empty(L, :rows)
    nrows(L_red) == 1 || throw(ArgumentError("Requires a single logical of the code."))
    is_logical(S, L_red) || throw(ArgumentError("The input matrix is not a logical of the code."))
    F = field(S)
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
    f0 = identity_matrix(F, size(stabs, 1))[:, nonzero]
    HX = transpose(stabs[nonzero, Q])

    if style == :Xanadu
        temp = size(HX, 2)
        HX = matrix(F, _add_edges(_Flint_matrix_to_Julia_int_matrix(HX)))
        # TODO: cellulate
        f0 = hcat(f0, zero_matrix(F, size(f0, 1), size(HX, 2) - temp))
        HZ = _remove_empty(rref(transpose(kernel(HX, side = :right)))[2], :rows)
        a = transpose(kernel(transpose(stabs), side = :right))
        b = _remove_empty(rref(a * f0)[2], :rows)
        if isempty(b)
            HZ = _find_low_weights_rand(HZ, max_iters)
        else
            HZ = _find_low_weight_cycle_subspace(HZ, b, max_iters)
        end
        return _thickened_cone(S, HX, HZ, f1, f0, type)
    elseif style == :IBM
        HZ = _remove_empty(rref(transpose(kernel(HX, side = :right)))[2], :rows)

        if remove_and_improve_cycles
            a = transpose(kernel(transpose(stabs), side = :right))
            b = _remove_empty(rref(a * f0)[2], :rows)
            if isempty(b)
                HZ = _find_low_weights_rand(HZ, max_iters)
            else
                HZ = _find_low_weight_cycle_subspace(HZ, b, max_iters)
            end
        elseif improve_cycles
            HZ = _find_low_weights_rand(HZ, max_iters)
        end

        r = ceil(Int, 1 / Cheeger_constant(_Flint_matrix_to_Julia_int_matrix(HX)))

        if type == :X
            return _thickened_cone(S, HX, HZ, f1, f0, type, r)
        else
            return _thickened_cone(S, HZ, HX, f1, f0, type, r)
        end
    elseif style == :Cohen
        HZ = zero_matrix(F, 0, size(HX, 2))
        if type == :X
            return _thickened_cone(S, HX, HZ, f1, f0, type, r)
        else
            return _thickened_cone(S, HZ, HX, f1, f0, type, r)
        end
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
    # we can remove some of these assertions since it's a private function and we control the input
    # @assert base_ring(stabs) == base_ring(logs)
    F = base_ring(stabs)
    @assert iseven(size(stabs, 2))
    @assert size(stabs, 2) == size(logs, 2)
    @assert rank(logs) == size(logs, 1)
    n = div(size(stabs, 2), 2)
    Ω = vcat(hcat(zero_matrix(F, n, n), identity_matrix(F, n)), hcat(identity_matrix(F, n),
        zero_matrix(F, n, n)))
    # @assert iszero(vcat(stabs, logs) * Ω * transpose(stabs))
    find_pairs = logs * Ω * transpose(logs)
    @assert all(count(!iszero, find_pairs[i:i, :]) in (0, 1) for i in 1:size(find_pairs, 1))
    needs_pair = findall(iszero(find_pairs[i:i, :]) for i in 1:size(find_pairs, 1))

    new_X = zero_matrix(F, 0, 2n)
    new_Z = zero_matrix(F, 0, 2n)
    new_mixed = zero_matrix(F, 0, 2n)

    for i in needs_pair
        RHS = zero_matrix(F, size(logs, 1) + size(stabs, 1) + n, 1)
        RHS[i, 1] = 1

        # try pure X
        # TODO don't understand this line
        LHS = vcat(logs, stabs, Ω[n + 1:2n, :])
        flag, sol = can_solve_with_solution(LHS * Ω, RHS, side = :right)
        if flag
            logs = vcat(logs, transpose(sol))
            new_X = vcat(new_X, transpose(sol))
            continue
        end

        # try pure Z
        LHS[size(logs, 1) + size(stabs, 1) + 1:end, :] = Ω[1:n, :]
        flag, sol = can_solve_with_solution(LHS * Ω, RHS, side = :right)
        if flag
            logs = vcat(logs, transpose(sol))
            new_Z = vcat(new_Z, transpose(sol))
            continue
        end

        # try mixed
        LHS[size(logs, 1) + size(stabs, 1) + 1:end, :] = zero_matrix(F, n, 2n)
        flag, sol = can_solve_with_solution(LHS * Ω, RHS)
        if flag
            logs = vcat(logs, transpose(sol))
            new_mixed = vcat(new_mixed, transpose(sol))
            continue
        end
    end

    return new_X, new_Z, new_mixed, logs
end

"""
$TYPEDSIGNATURES

Return the Cheeger constant of the matrix `M` assuming `M` is a vertex-edge incidence matrix.
"""
function Cheeger_constant(M::Matrix{T}) where T <: Integer
    m, n = size(M)
    # get one more bit in using an unsigned integer...
    U = UInt64
    if m >= 64
        error("Not implemented for more than 64 vertices.")
    end
    r = div(m, 2)
    h = Inf
    for x in U(1):U(2)^U(m - 1)
        v = digits(T, x, base = 2, pad = m)
        s = sum(v)
        s > r && continue
        h = min(h, count(isodd, dot(v, M[:, c]) for c in 1:n) / s)
    end
    return h
end

# TODO finish
"""
$TYPEDSIGNATURES

Return 
"""
function Cheeger_constant(S::AbstractSubsystemCode, L::CTMatrixTypes)
    # TODO why subsystem if never using gauges?
    # TODO do we require it not be a graph state here?

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

function _sparsest_cut(M::Matrix{T}; rng = Xoshiro()) where T <: Integer
    m, n = size(M)
    r = div(m, 2)
    h = Inf
    sparse_cut = zeros(T, m)
    # get one more bit in using an unsigned integer...
    U = UInt64
    if m >= 64
        error("Not implemented for more than 64 vertices.")
    end

    # the shuffle allows different choices of the sparsest cuts to be chosen
    for x in shuffle(rng, U(1):U(2)^U(m - 1))
        # v corresponds to a subset of the vertices
        v = digits(T, x, base = 2, pad = m)
        s = sum(v)
        s > r && continue
        temp = count(isodd, dot(v, M[:, c]) for c in 1:n) / s
        if temp < h
            h = temp
            sparse_cut .= v
        end
    end

    return h, sparse_cut
end

function _add_edges(M::Matrix{T}; rng = Xoshiro()) where T <: Integer
    M_new = copy(M)
    while true

        # check the Cheeger constant, quit if it's already good
        h, v = _sparsest_cut(M_new, rng = rng)
        h < 1 || break

        # M_sum is a vector of the vertex weights
        M_sum = vec(sum(M_new, dims = 2))

        # # i is a vertex in supp(v) with minimal degree, j is a vertex in V \ supp(v) with minimal degree
        # i = argmin(i -> iszero(v[i]) ? typemax(T) : M_sum[i], axes(M_new, 1))
        # j = argmin(j -> !iszero(v[j]) ? typemax(T) : M_sum[j], axes(M_new, 1))

        # # add an edge that connects these two vertices
        # M_new = [M_new zeros(T, size(M_new, 1), 1)]
        # M_new[i, end] = M_new[j, end] = T(1)

        S = findall(!iszero, v)
        Sc = findall(iszero, v)
        S_min = minimum(M_sum[S])
        Sc_min = minimum(M_sum[Sc])
        is = S[findall(M_sum[S[i]] == S_min for i in eachindex(S))]
        js = Sc[findall(M_sum[Sc[i]] == Sc_min for i in eachindex(Sc))]
        h = -Inf
        M_new = [M_new zeros(T, size(M_new, 1), 1)]
        M_temp = copy(M_new)
        Tzero = T(0)
        Tone = T(1)
        for (i, j) in Iterators.product(shuffle(rng, is), shuffle(rng, js))
            M_temp[i, end] = M_temp[j, end] = Tone
            temp = Cheeger_constant(M_temp)
            if temp > h
                h = temp
                M_new[:, end] .= Tzero
                M_new[i, end] = M_new[j, end] = Tone
            end
            M_temp[i, end] = M_temp[j, end] = Tzero
        end
    end

    return M_new
end

function _find_low_weight_cycle_subspace(all_cycles::CTMatrixTypes,
    already_covered_cycles::CTMatrixTypes, max_iters::Int, f::T = maximum) where T <: Function

    @assert size(all_cycles, 2) == size(already_covered_cycles, 2)
    F = base_ring(all_cycles)
    A = _remove_empty(rref(already_covered_cycles)[2], :rows)
    pivots = findall(count(!iszero(A[i, j]) for i in axes(A, 1)) == 1 for j in axes(A, 2))
    B = _remove_empty(rref(all_cycles)[2], :rows)
    for i in axes(B, 1)
        for j in axes(A, 1)
            if !iszero(B[i, pivots[j]])
                x = inv(B[i, pivots[j]])
                for k in axes(B, 2)
                    B[i, k] = x * A[j, k] + B[i, k]
                end
            end
        end
    end
    B = _remove_empty(rref(B)[2], :rows)
    isempty(B) && (return B;)

    w = size(A, 2) + 1
    C = zero_matrix(F, size(B, 1), size(B, 2))
    for i in 1:max_iters
        x = _random_invertible_matrix(F, size(B, 1))
        y = matrix(F, rand(F, size(B, 1), size(A, 1)))
        temp = x * B + y * A
        temp_w = f(count(!iszero(temp[i, j]) for j in size(temp, 2)) for i in size(temp, 1))
        if temp_w < w
            w = temp_w
            C = temp
        end
        temp = x * B
        temp_w = f(count(!iszero(temp[i, j]) for j in size(temp, 2)) for i in size(temp, 1))
        if temp_w < w
            w = temp_w
            C = temp
        end
    end
    return C
end

function _find_low_weights_rand(M::CTMatrixTypes, max_iters::Int, f::T = maximum) where
    T <: Function

    initial_w = f(count(!iszero(M[i, j]) for j in size(M, 2)) for i in size(M, 1))
    A = _remove_empty(rref(M)[2], :rows)
    isempty(A) && (return A;)
    F = base_ring(A)
    w = size(A, 2) + 1
    for i in 1:max_iters
        x = _random_invertible_matrix(F, size(A, 1))
        temp = x * A
        temp_w = f(count(!iszero(temp[i, j]) for j in size(temp, 2)) for i in size(temp, 1))
        if temp_w < w
            w = temp_w
            A = temp
        end
    end
    return initial_w < w ? M : A
end

function _random_invertible_matrix(n::Int)
    # @assert n > 0
    inds = collect(1:n)
    A = zeros(UInt8, n, n)
    T = zeros(UInt8, n, n)
    for k in 1:n
        v = rand(0x00:0x01, n - k + 1)
        while iszero(v)
            v .= rand(0x00:0x01, n - k + 1)
        end
        r = findfirst(!iszero, v)
        A[k, inds[r]] = 0x01
        A[k + 1:end, inds[r]] .= rand(0x00:0x01, n - k)
        T[inds[r], inds] .= v
        deleteat!(inds, r)
    end
    return A * T .% 0x02
end

function _random_invertible_matrix(F::CTFieldTypes, n::Int)
    # @assert n > 0
    inds = collect(1:n)
    A = zero_matrix(F, n, n)
    T = zero_matrix(F, n, n)
    Fone = F(1)
    for k in 1:n
        v = matrix(F, rand(F, 1, n - k + 1))
        while iszero(v)
            v = matrix(F, rand(F, 1, n - k + 1))
        end
        r = findfirst(!iszero, v)[2]
        A[k, inds[r]] = Fone
        for i in k + 1:n
            A[i, inds[r]] = rand(F)
        end
        for i in eachindex(inds)
            T[inds[r], inds[i]] = v[1, i]
        end
        deleteat!(inds, r)
    end
    return A * T
end
