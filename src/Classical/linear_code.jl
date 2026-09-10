# Copyright (c) 2021 - 2026 Eric Sabo, Benjamin Ide, David Marquis
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

# CodingTheory also exports `nonzeros` for cyclic/Goppa codes, so SparseArrays
# accessors must be fully qualified on CSC matrices.
Oscar.nrows(A::SparseMatrixCSC) = size(A, 1)
Oscar.ncols(A::SparseMatrixCSC) = size(A, 2)

function _code_matrix_base_ring(A::SparseMatrixCSC)
    vals = SparseArrays.nonzeros(A)
    if eltype(A) <: Integer
        all(x -> iszero(x) || isone(x), vals) ||
            throw(ArgumentError("Integer sparse code matrices are interpreted over GF(2) and must contain only zeros and ones."))
        return Oscar.Nemo.Native.GF(2)
    end

    isempty(vals) && throw(ArgumentError(
        "Cannot infer the finite field of an all-zero SparseMatrixCSC. Use an Oscar matrix with an explicit base ring."))
    F = parent(first(vals))
    all(x -> parent(x) == F, vals) ||
        throw(ArgumentError("All sparse matrix entries must belong to the same finite field."))
    return F
end
_code_matrix_base_ring(A::CTMatrixTypes) = base_ring(A)

function _binary_sparse_rank(A::SparseMatrixCSC)
    nr, nc = size(A)
    chunks = cld(nc, 64)
    rows = [zeros(UInt64, chunks) for _ in 1:nr]
    vals = SparseArrays.nonzeros(A)
    rows_idx = SparseArrays.rowvals(A)

    @inbounds for c in 1:nc
        chunk = (c - 1) ÷ 64 + 1
        bit = UInt64(1) << ((c - 1) % 64)
        for ptr in SparseArrays.nzrange(A, c)
            iszero(vals[ptr]) || (rows[rows_idx[ptr]][chunk] |= bit)
        end
    end

    pivots = Dict{Int, Vector{UInt64}}()
    rnk = 0
    for row in rows
        while true
            pivot = 0
            for chunk in eachindex(row)
                if !iszero(row[chunk])
                    pivot = 64 * (chunk - 1) + trailing_zeros(row[chunk]) + 1
                    break
                end
            end
            iszero(pivot) && break

            if haskey(pivots, pivot)
                @inbounds @simd for chunk in eachindex(row)
                    row[chunk] ⊻= pivots[pivot][chunk]
                end
            else
                pivots[pivot] = row
                rnk += 1
                break
            end
        end
    end
    return rnk
end

function _code_matrix_rank(A::SparseMatrixCSC, F::CTFieldTypes)
    Int(order(F)) == 2 && return _binary_sparse_rank(A)
    return rank(_dense_code_matrix(A, F))
end
_code_matrix_rank(A::SMat, F::CTFieldTypes) = rank(matrix(A))
_code_matrix_rank(A::CTMatrixTypes, F::CTFieldTypes) = rank(A)

function _dense_code_matrix(A::SparseMatrixCSC, F::CTFieldTypes)
    B = zero_matrix(F, size(A)...)
    vals = SparseArrays.nonzeros(A)
    rows = SparseArrays.rowvals(A)
    @inbounds for c in axes(A, 2), ptr in SparseArrays.nzrange(A, c)
        B[rows[ptr], c] = F(vals[ptr])
    end
    return B
end
_dense_code_matrix(A::SMat, F::CTFieldTypes) = matrix(A)
_dense_code_matrix(A::CTMatrixTypes, F::CTFieldTypes) = A

_is_sparse_code_matrix(A) = A isa Union{SparseMatrixCSC, SMat}
_sparse_code_matrix(A::Union{SparseMatrixCSC, SMat}) = A
function _sparse_code_matrix(A::CTMatrixTypes)
    F = base_ring(A)
    if Int(order(F)) == 2
        values = [iszero(A[r, c]) ? 0 : 1
                  for r in 1:nrows(A), c in 1:ncols(A)]
        return sparse(values)
    end
    iszero(A) && return sparse_matrix(F, nrows(A), ncols(A))
    return sparse_matrix(F, Matrix(A))
end
_normalize_quantum_matrix(A::SMat) =
    matrix(A)
_normalize_quantum_matrix(A::CTMatrixTypes) = A

function _prime_subfield(F::CTFieldTypes)
    degree(F) == 1 && return F
    return GF(Int(characteristic(F)))
end

function _additive_expansion(A::CTMatrixTypes, F::CTFieldTypes=base_ring(A))
    dense = _dense_code_matrix(A, F)
    degree(F) == 1 && return dense
    prime_field = _prime_subfield(F)
    basis, _ = primitive_basis(F, prime_field)
    return expand_matrix(dense, prime_field, basis)
end

function _additive_rank(A::CTMatrixTypes, F::CTFieldTypes=base_ring(A))
    degree(F) == 1 && return _code_matrix_rank(A, F)
    return rank(_additive_expansion(A, F))
end

"""
$(TYPEDSIGNATURES)

Return the linear code constructed with matrix `G`. If `parity` is `true`, `G` is 
treated as the parity-check matrix. Safely handles sparse matrices and uses lazy 
evaluation to avoid eager rank and dual computations.
"""
function LinearCode(G::CTMatrixTypes, parity::Bool = false)
    is_sparse = G isa SparseMatrixCSC || G isa SMat
    F = _code_matrix_base_ring(G)
    iszero(G) && return parity ? IdentityCode(F, ncols(G)) : ZeroCode(F, ncols(G))

    G_clean = G isa SMat ? copy(G) : _remove_empty(deepcopy(G), :rows)
    n = ncols(G_clean)
    
    # Sparse parity checks are often overcomplete, so nrows(G) is not a valid
    # substitute for rank. Binary CSC matrices use packed GF(2) elimination.
    actual_rank = _code_matrix_rank(G_clean, F)
    cache = Dict{Symbol, Any}()

    if parity
        k = n - actual_rank
        k == n && return IdentityCode(F, n)
        
        ub = is_sparse ? n : _min_wt_row(G_clean)[1]
        cache[:H] = G_clean
        return LinearCode(F, n, k, missing, 1, ub, cache)
    else
        k = actual_rank
        k == n && return IdentityCode(F, n)
        
        ub = is_sparse ? n : _min_wt_row(G_clean)[1]
        cache[:G] = G_clean
        return LinearCode(F, n, k, missing, 1, ub, cache)
    end
end

"""
$(TYPEDSIGNATURES)

Construct a linear code explicitly providing both `G` and `H`. 
Safely handles mixed sparse and dense matrix types.
"""
function LinearCode(G::CTMatrixTypes, H::CTMatrixTypes; check_orthogonality::Bool = false)
    ncols(G) == ncols(H) || throw(ArgumentError("G and H must have the same number of columns (received ncols(G) = $(ncols(G)), ncols(H) = $(ncols(H)))."))
    
    is_sparse = (G isa SparseMatrixCSC || G isa SMat) || (H isa SparseMatrixCSC || H isa SMat)
    
    F = _code_matrix_base_ring(G)
    F == _code_matrix_base_ring(H) || throw(ArgumentError("G and H must be over the same field."))
    G_clean = G isa SMat ? copy(G) : _remove_empty(deepcopy(G), :rows)
    H_clean = H isa SMat ? copy(H) : _remove_empty(deepcopy(H), :rows)
    
    if check_orthogonality
        G_temp = _dense_code_matrix(G_clean, F)
        H_temp = _dense_code_matrix(H_clean, F)
        iszero(G_temp * transpose(H_temp)) || throw(ArgumentError("H is not orthogonal to G."))
    end
    
    k = _code_matrix_rank(G_clean, F)
    _code_matrix_rank(H_clean, F) == ncols(G_clean) - k ||
        throw(ArgumentError("H does not have the correct rank to be a parity-check matrix for G."))

    ub1, _ = is_sparse ? (ncols(G_clean), 0) : _min_wt_row(G_clean)
    ub2, _ = is_sparse ? (ncols(H_clean), 0) : _min_wt_row(H_clean)
    ub = min(ub1, ub2)
    
    cache = Dict{Symbol, Any}(:G => G_clean, :H => H_clean)
    return LinearCode(F, ncols(G_clean), k, missing, 1, ub, cache)
end

"""
$(TYPEDSIGNATURES)

Construct a linear code from a Julia `Matrix{Int}` over the finite field of order `q`.
"""
function LinearCode(G::Matrix{Int}, q::Int, parity::Bool = false)
    factors = Nemo.factor(q)
    (length(factors) == 1 && q > 1) || throw(ArgumentError("There is no finite field of order $q."))

    p, m = first(factors)
    F = m == 1 ? GF(p) : GF(p, m, :ω)
    G1 = matrix(F, G)
    
    # REMOVED: rref!(G1) to preserve the user's input topology!
    return LinearCode(G1, parity)
end

"""
$(TYPEDSIGNATURES)

Construct a linear code by vertically concatenating a vector of matrices.
"""
function LinearCode(Gs::Vector{<:CTMatrixTypes})
    s = size(Gs[1])
    all(s == size(Gs[i]) for i in 2:length(Gs)) || throw(ArgumentError("Not all matrices in `Gs` have the same dimensions."))

    G = reduce(vcat, Gs)
    # REMOVED: rref!(G) to preserve the user's input topology!
    return LinearCode(G)
end

"""
$(TYPEDSIGNATURES)

Construct a linear code from a vector of Julia `Vector{Int}` over the finite field of order `q`.
"""
function LinearCode(Gs::Vector{Vector{Int}}, q::Int, parity::Bool = false)
    s = length(Gs[1])
    all(s == length(Gs[i]) for i in 2:length(Gs)) || throw(ArgumentError("Not all vectors in `Gs` are the same length."))

    # Safely convert Vector{Vector{Int}} to Matrix{Int} where each vector is a row
    G_mat = permutedims(reduce(hcat, Gs))
    return LinearCode(G_mat, q, parity)
end

"""
$(TYPEDSIGNATURES)

Return a random `[n, k]` linear code over `F`. 
Bypasses standard constructor checks for massive speedups on large codes.
"""
function random_linear_code(F::CTFieldTypes, n::Int, k::Int; rng::AbstractRNG = Random.seed!())
    rand_mat = zero_matrix(F, k, n - k)
    for r in 1:k, c in 1:(n - k) 
        rand_mat[r, c] = rand(rng, F) 
    end
    
    # The matrix is inherently in standard form [I_k | A]
    G_stand = hcat(identity_matrix(F, k), rand_mat)
    
    # Bypass the O(n^3) rank check since we know rank == k
    ub, _ = _min_wt_row(G_stand)
    
    # Inject directly into the cache
    cache = Dict{Symbol, Any}(:G => G_stand, :G_stand => G_stand)
    return LinearCode(F, n, k, missing, 1, ub, cache)
end 

"""
$(TYPEDSIGNATURES)

Return a random `[n, k]` linear code over `GF(q)`.
"""
function random_linear_code(q::Int, n::Int, k::Int; rng::AbstractRNG = Random.seed!())
    _, e, p = is_prime_power_with_data(q)
    if e == 1 
        field = Oscar.Nemo.Native.GF(p)
    elseif e > 1
        field = GF(p, e, :x)
    else
        throw(ArgumentError("q must be a prime power"))
    end
    return random_linear_code(field, n, k, rng = rng)
end

#############################
      # getter functions
#############################

"""
$(TYPEDSIGNATURES)

Return the base ring of the generator matrix of `C`.
"""
field(C::AbstractLinearCode) = C.F

"""
$(TYPEDSIGNATURES)

Return the length of `C`.
"""
length(C::AbstractLinearCode) = C.n

"""
$(TYPEDSIGNATURES)

Return the dimension of `C`.
"""
dimension(C::AbstractLinearCode) = C.k

"""
$(TYPEDSIGNATURES)

Return the cardinality of `C`.
"""
cardinality(C::AbstractLinearCode) = BigInt(order(C.F))^C.k

"""
$(TYPEDSIGNATURES)

Return the rate of `C`.
"""
rate(C::AbstractLinearCode) = C.k / C.n

"""
$(TYPEDSIGNATURES)

Return the generator matrix of `C`. Evaluates lazily and caches the result.
If the optional parameter `stand_form` is set to `true`, the standard form is returned.
"""
function generator_matrix(C::AbstractLinearCode, stand_form::Bool = false)
    cache = getfield(C, :cache)
    
    if !haskey(cache, :G)
        if isa(C, AbstractQuasiCyclicCode)
            if C.A_type == :G
                cache[:G] = lift(C.A)
            else
                G_raw = kernel(lift(C.A), side = :right)
                rnk_G = rank(G_raw)
                if ncols(G_raw) == rnk_G
                    cache[:G] = transpose(G_raw)
                else
                    nr = nrows(G_raw)
                    G_tr = zero_matrix(base_ring(G_raw), rnk_G, nr)
                    for r in 1:nr, c in 1:rnk_G
                        !iszero(G_raw[r, c]) && (G_tr[c, r] = G_raw[r, c])
                    end
                    cache[:G] = G_tr
                end
            end
        # TODO: Add `elseif isa(C, AbstractCyclicCode)` here to build from C.g
        else
            !haskey(cache, :H) && error("Fatal: Neither G nor H found in cache for a standard LinearCode.")
            H_mat = cache[:H]
            
            if H_mat isa SparseMatrixCSC || H_mat isa SMat
                @warn "Computing the generator matrix from a sparse parity-check matrix requires computing the kernel. This will likely result in a dense matrix and may cause severe memory explosion for large codes (n > 10,000)."
                H_mat = _dense_code_matrix(H_mat, C.F)
            end
            
            G_raw = kernel(H_mat, side = :right)
            rnk_G = rank(G_raw)
            
            if ncols(G_raw) == rnk_G
                cache[:G] = transpose(G_raw)
            else
                nr = nrows(G_raw)
                G_tr = zero_matrix(base_ring(G_raw), rnk_G, nr)
                for r in 1:nr, c in 1:rnk_G
                    !iszero(G_raw[r, c]) && (G_tr[c, r] = G_raw[r, c])
                end
                cache[:G] = G_tr
            end
        end
    end
    
    if stand_form
        if !haskey(cache, :G_stand)
            G_for_standard_form = _dense_code_matrix(cache[:G], C.F)
            G_stand, H_stand, P_stand, _ = _standard_form(G_for_standard_form)
            cache[:G_stand] = G_stand
            cache[:H_stand] = H_stand
            cache[:P_stand] = P_stand
        end
        return cache[:G_stand]
    end
    
    return cache[:G]
end

"""
$(TYPEDSIGNATURES)

Return the parity-check matrix of `C`. Evaluates lazily and caches the result.
If the optional parameter `stand_form` is set to `true`, the standard form is returned.
"""
function parity_check_matrix(C::AbstractLinearCode, stand_form::Bool = false)
    cache = getfield(C, :cache)
    
    if !haskey(cache, :H)
        if isa(C, AbstractQuasiCyclicCode)
            if C.A_type == :H
                cache[:H] = lift(C.A)
            else
                H_raw = kernel(lift(C.A), side = :right)
                cache[:H] = transpose(H_raw)
            end
        # TODO: Add `elseif isa(C, AbstractCyclicCode)` here to build from C.h
        else
            !haskey(cache, :G) && error("Fatal: Neither G nor H found in cache for a standard LinearCode.")
            G_mat = cache[:G]
            
            if G_mat isa SparseMatrixCSC || G_mat isa SMat
                @warn "Computing the parity-check matrix from a sparse generator matrix requires computing the kernel. This will likely result in a dense matrix and may cause severe memory explosion for large codes (n > 10,000)."
                G_mat = _dense_code_matrix(G_mat, C.F)
            end
            
            H_raw = kernel(G_mat, side = :right)
            rnk_H = rank(H_raw)
            
            if ncols(H_raw) == rnk_H
                cache[:H] = transpose(H_raw)
            else
                nr = nrows(H_raw)
                H_tr = zero_matrix(base_ring(H_raw), rnk_H, nr)
                for r in 1:nr, c in 1:rnk_H
                    !iszero(H_raw[r, c]) && (H_tr[c, r] = H_raw[r, c])
                end
                cache[:H] = H_tr
            end
        end
    end
    
    if stand_form
        if !haskey(cache, :H_stand)
            G_stand, H_stand, P_stand, _ = _standard_form(generator_matrix(C)) 
            cache[:G_stand] = G_stand
            cache[:H_stand] = H_stand
            cache[:P_stand] = P_stand
        end
        return cache[:H_stand]
    end
    
    return cache[:H]
end

"""
$(TYPEDSIGNATURES)

Return the permutation matrix required to permute the columns of the code matrices 
to have the same row space as the matrices in standard form. 
Evaluates lazily if not already computed.
"""
function standard_form_permutation(C::AbstractLinearCode)
    cache = getfield(C, :cache)
    if !haskey(cache, :P_stand)
        # Forcing the standard form getter naturally populates P_stand in the cache
        generator_matrix(C, true)
    end
    return cache[:P_stand]
end

"""
$(TYPEDSIGNATURES)

Return the relative minimum distance of `C` if ``d`` is known;
otherwise return `missing`.
"""
relative_distance(C::AbstractLinearCode) = ismissing(C.d) ? missing : C.d / C.n

"""
$(TYPEDSIGNATURES)

Return the genus, ``n + 1 - k - d``, of the code.
"""
genus(C::AbstractLinearCode) = ismissing(C.d) ? missing : C.n + 1 - C.k - C.d

"""
$(TYPEDSIGNATURES)

Return the current lower bound on the minimum distance of `C`.
"""
minimum_distance_lower_bound(C::AbstractLinearCode) = C.l_bound

"""
$(TYPEDSIGNATURES)

Return the current upper bound on the minimum distance of `C`.
"""
minimum_distance_upper_bound(C::AbstractLinearCode) = C.u_bound

"""
$(TYPEDSIGNATURES)

Return `true` if code is maximum distance separable (MDS).
Calculates the exact minimum distance if currently missing.
"""
function is_MDS(C::AbstractLinearCode)
    # BUG FIX: Was previously calling minimum_distance_Gray(C)
    if ismissing(C.d)
        minimum_distance(C) # Forces calculation and updates C.d
    end
    return C.d == Singleton_bound(C.n, C.k)
end

"""
$(TYPEDSIGNATURES)

Return the number of correctable errors for the code.

# Notes
* The number of correctable errors is ``t = \\floor{(d - 1) / 2}``.
"""
number_correctable_errors(C::AbstractLinearCode) = ismissing(C.d) ? missing : Int(fld(C.d - 1, 2))

"""
$(TYPEDSIGNATURES)

Return `true` if the generator matrix is overcomplete. If the optional parameter is
set to `:H`, returns `true` if the parity-check matrix is overcomplete.
"""
function is_overcomplete(C::AbstractLinearCode, which::Symbol=:G)
    if which == :G
        # Must use getter instead of C.G to ensure it evaluates if missing
        return nrows(generator_matrix(C)) > C.k
    elseif which == :H
        # Must use getter instead of C.H to ensure it evaluates if missing
        return nrows(parity_check_matrix(C)) > C.n - C.k
    end
    throw(ArgumentError("Received symbol $which, expected :G or :H"))
end

#############################
      # setter functions
#############################

"""
$(TYPEDSIGNATURES)

Set the lower bound on the minimum distance of `C`, if `l` is better than the current bound.
"""
function set_distance_lower_bound!(C::AbstractLinearCode, l::Int)
    1 <= l <= C.u_bound || throw(DomainError(l, "The lower bound must be between 1 and the upper bound (currently $(C.u_bound))."))
    C.l_bound < l && (C.l_bound = l)
    if C.l_bound == C.u_bound
        @info "The new lower bound is equal to the upper bound; setting the minimum distance."
        C.d = C.l_bound
    end
end

"""
$(TYPEDSIGNATURES)

Set the upper bound on the minimum distance of `C`, if `u` is better than the current bound.
"""
function set_distance_upper_bound!(C::AbstractLinearCode, u::Int)
    C.l_bound <= u <= C.n || throw(DomainError(u, "The upper bound must be between the lower bound (currently $(C.l_bound)) and the code length."))
    u < C.u_bound && (C.u_bound = u)
    if C.l_bound == C.u_bound
        @info "The new upper bound is equal to the lower bound; setting the minimum distance."
        C.d = C.l_bound
    end
end

"""
$(TYPEDSIGNATURES)

Set the exact minimum distance of the code to `d` and update the bounds.
"""
function set_minimum_distance!(C::AbstractLinearCode, d::Int)
    0 < d <= C.n || throw(DomainError(d, "The minimum distance of a code must be between 1 and n; received: d = $d."))
    C.d = d
    C.l_bound = d
    C.u_bound = d
end

"""
$(TYPEDSIGNATURES)

In-place change of the base ring of `C` to `F`. Safely updates all cached matrices.
"""
function change_field!(C::T, F::CTFieldTypes) where T <: AbstractLinearCode
    # TODO: implementation for cyclic codes
    T <: AbstractCyclicCode && @error "Not implemented for cyclic codes yet."

    cache = getfield(C, :cache)
    
    # Safely convert all cached matrices to the new field if they exist
    haskey(cache, :G) && (cache[:G] = change_base_ring(F, cache[:G]))
    haskey(cache, :H) && (cache[:H] = change_base_ring(F, cache[:H]))
    haskey(cache, :G_stand) && (cache[:G_stand] = change_base_ring(F, cache[:G_stand]))
    haskey(cache, :H_stand) && (cache[:H_stand] = change_base_ring(F, cache[:H_stand]))
    haskey(cache, :P_stand) && (cache[:P_stand] = change_base_ring(F, cache[:P_stand]))

    # Wipe cached enumerators/distributions if the field order changes
    if order(F) != order(C.F)
        delete!(cache, :weight_enum)
        delete!(cache, :weight_dist)
    end

    C.F = F
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return a new linear code which is `C` changed to the base ring `F`.
"""
function change_field(C::AbstractLinearCode, F::CTFieldTypes)
    C2 = deepcopy(C)
    change_field!(C2, F)
    return C2
end

#############################
     # general functions
#############################

# # Arguments
# - `C` - a linear code 
# ==============================================================================
# 8. INFORMATION SETS
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return a set of column indices corresponding to an information set of `C`.
"""
function information_set(C::AbstractLinearCode)
    # Safely use the getter to ensure G is populated
    G = generator_matrix(C)
    _, _, perm = _rref_col_swap_perm(G)
    return on_sets(collect(1:C.k), inv(perm)) 
end

"""
$(TYPEDSIGNATURES)

Return a set of column indices corresponding to a random information set of `C`.
"""
function random_information_set(C::AbstractLinearCode; rng::AbstractRNG = Random.seed!())
    # Safely use the getter to ensure G is populated
    G = generator_matrix(C)
    
    shuffle_perm_julia = shuffle(rng, collect(1:C.n)) 
    shuffle_perm = perm(shuffle_perm_julia)
    
    # apply shuffle permutation
    permuted_mat = G[:, invperm(shuffle_perm_julia)] 
    
    # transform to rref. The first k columns of the rref matrix are pivot columns
    _, _, rref_perm = CodingTheory._rref_col_swap_perm(permuted_mat)
    inv_perm = inv(shuffle_perm * rref_perm)
    
    # apply the inverse permutation to the pivots
    return on_sets(collect(1:C.k), inv_perm) 
end

# ==============================================================================
# 9. STANDARD FORM UTILITY
# ==============================================================================

# Internal function, no docstring needed
function _standard_form(G::CTMatrixTypes)
    rnk, G_stand, P = _rref_col_swap(G, 1:nrows(G), 1:ncols(G))
    F = base_ring(G_stand)
    
    nrows(G_stand) > rnk && (G_stand = _remove_empty(G_stand, :rows))
    A_tr = transpose(view(G_stand, :, (nrows(G_stand) + 1):ncols(G_stand)))
    
    H_stand = hcat(order(F) == 2 ? A_tr : -A_tr, identity_matrix(F, ncols(G_stand) - nrows(G_stand)))
    return G_stand, H_stand, P, rnk
end

# ==============================================================================
# 10. REPL DISPLAY
# ==============================================================================

function Base.show(io::IO, C::AbstractLinearCode)
    print(io, "[$(C.n), $(C.k)")
    !ismissing(C.d) && print(io, ", $(C.d)")
    typeof(C) <: AbstractBCHCode && print(io, "; $(C.b)")
    print(io, "]_$(order(C.F)) ")
    
    if isa(C, ReedSolomonCode)
        println(io, "Reed-Solomon code")
    elseif isa(C, BCHCode)
        println(io, "BCH code")
    elseif isa(C, CyclicCode)
        println(io, "cyclic code")
    elseif isa(C, ReedMullerCode)
        println(io, "Reed-Muller code RM($(C.r), $(C.m))")
    elseif isa(C, QuasiCyclicCode)
        println(io, "quasi-cyclic code of index $(C.l)")
    elseif isa(C, GeneralizedReedSolomonCode)
        println(io, "generalized Reed-Solomon code")
    elseif isa(C, GoppaCode)
        println(io, "Goppa code")
    elseif isa(C, TwistedReedSolomonCode)
        println(io, "twisted Reed-Solomon code")
    else
        println(io, "linear code")
    end

    if get(io, :compact, true)
        if typeof(C) <: AbstractCyclicCode
            println(io, "$(order(C.F))-Cyclotomic cosets: ")
            len = length(qcosets_reps(C))
            if len == 1
                println("\tC_$(qcosets_reps(C)[1])")
            else
                for (i, x) in enumerate(qcosets_reps(C))
                    if i == 1
                        print(io, "\tC_$x ∪ ")
                    elseif i == 1 && i == len
                        println(io, "\tC_$x")
                    elseif i != len
                        print(io, "C_$x ∪ ")
                    else
                        println(io, "C_$x")
                    end
                end
            end
            println(io, "Generator polynomial:")
            println(io, "\t", generator_polynomial(C))
        elseif isa(C, GoppaCode)
            len = length(C.L)
            if len ≤ 20
                println(io, "Evaluated at:")
                print(io, "\t[")
                for i in 1:len
                    if i ≠ len
                        print(C.L[i], ", ")
                    else
                        println(C.L[i], "]")
                    end
                end
            end
            println(io, "Goppa polynomial:")
            println(io, "\t", C.g)
        elseif isa(C, TwistedReedSolomonCode)
            println(io, "Number of twists: $(C.l)")
            if C.l ≤ 20
                println(io, "Twist vector:")
                print(io, "\t[")
                for i in 1:C.l
                    if i ≠ C.l
                        print(C.t[i], ", ")
                    else
                        println(C.t[i], "]")
                    end
                end
                println(io, "Hook vector:")
                print(io, "\t[")
                for i in 1:C.l
                    if i ≠ C.l
                        print(C.h[i], ", ")
                    else
                        println(C.h[i], "]")
                    end
                end
                println(io, "Coefficient vector:")
                print(io, "\t[")
                for i in 1:C.l
                    if i ≠ C.l
                        print(C.η[i], ", ")
                    else
                        println(C.η[i], "]")
                    end
                end
            end
            if C.n ≤ 20
                println(io, "Evaluated at:")
                print(io, "\t[")
                for i in 1:C.n
                    if i ≠ C.n
                        print(C.α[i], ", ")
                    else
                        println(C.α[i], "]")
                    end
                end
            end
        end

        if C.n ≤ 30 && C.k ≠ 0
            if isa(C, QuasiCyclicCode)
                if C.A_type == :G
                    M = generator_matrix(C)
                    nr, nc = size(M)
                    println(io, "Generator matrix: $(nr) × $(nc)")
                else
                    M = parity_check_matrix(C)
                    nr, nc = size(M)
                    println(io, "Parity-check matrix: $(nr) × $(nc)")
                end
                for i in 1:nr
                    print(io, "\t")
                    for j in 1:nc
                        if j != nc
                            print(io, "$(M[i, j]) ")
                        elseif j == nc && i != nr
                            println(io, "$(M[i, j])")
                        else
                            print(io, "$(M[i, j])")
                        end
                    end
                end
            else
                G = generator_matrix(C)
                nr, nc = size(G)
                println(io, "Generator matrix: $nr × $nc")
                println(io, "\t" * replace(replace(replace(repr(MIME("text/plain"), G), 
                r"\n" => "\n\t"), r"\[" => ""), r"\]" => ""))
            end
        end
    end
end

# ==============================================================================
# 11. SINGLETON BOUND
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the Singleton bound ``d \\leq n - k + 1`` or ``k \\leq n - d + 1`` depending on the interpretation of `a`.
"""
Singleton_bound(n::Int, a::Int) = 0 <= a <= n ? (n - a + 1) : error("Invalid parameters for the Singleton bound. Received n = $n, k/d = $a")

"""
$(TYPEDSIGNATURES)

Return the Singleton bound on the minimum distance of the code (``d \\leq n - k + 1``).
"""
Singleton_bound(C::AbstractLinearCode) = Singleton_bound(C.n, C.k)

# ==============================================================================
# 12. ENCODING & SYNDROME
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the encoding of `v` into `C`.
"""
function encode(C::AbstractLinearCode, v::Union{CTMatrixTypes, Vector{Int}})
    w = isa(v, Vector{Int}) ? matrix(C.F, 1, length(v), v) : v
    G = generator_matrix(C) # Safe getter
    nr = nrows(G)
    
    (size(w) != (1, nr) && size(w) != (nr, 1)) &&
        throw(ArgumentError("Vector has incorrect dimension; expected length $nr, received: $(size(v))."))
    base_ring(w) == C.F || throw(ArgumentError("Vector must have the same base ring as the generator matrix."))
    
    nrows(w) != 1 || return w * G
    return transpose(w) * G
end

"""
$(TYPEDSIGNATURES)

Return the syndrome of `v` with respect to `C`.
"""
function syndrome(C::AbstractLinearCode, v::Union{CTMatrixTypes, Vector{Int}, Vector{fpFieldElem}, Vector{FpFieldElem}})
    w = isa(v, Union{Vector{Int}, Vector{fpFieldElem}, Vector{FpFieldElem}}) ? matrix(C.F, length(v), 1, v) : v
    H = parity_check_matrix(C) # Safe getter
    nc = ncols(H)
    
    (size(w) != (nc, 1) && size(w) != (1, nc)) &&
        throw(ArgumentError("Vector has incorrect dimension; expected length $nc, received: $(size(v))."))
        
    if base_ring(w) != C.F
        if order(base_ring(w)) == order(C.F)
            @warn "Fields are of different types, but have the same order."
        else
            throw(ArgumentError("Vector must have the same base ring as the parity-check matrix."))
        end
    end
    
    return nrows(w) == 1 ? H * transpose(w) : H * w
end

# ==============================================================================
# 13. CODE MEMBERSHIP & SUBCODES
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return whether or not `v` is a codeword of `C`.
"""
function Base.in(v::Union{CTMatrixTypes, Vector{Int}, Vector{fpFieldElem}, Vector{FpFieldElem}}, C::AbstractLinearCode)
    return iszero(syndrome(C, v))
end
# # Also cover the unicode symbol
# ∈(v::Union{CTMatrixTypes, Vector{Int}, Vector{fpFieldElem}, Vector{FpFieldElem}}, C::AbstractLinearCode) = in(v, C)

"""
$(TYPEDSIGNATURES)

Return whether or not `C1` is a subcode of `C2`.
"""
function ⊆(C1::AbstractLinearCode, C2::AbstractLinearCode)
    C1.F == C2.F || (order(C1.F) == order(C2.F) ? (@warn "Fields are of different types, but have the same order.") : (return false;))
    (C1.n == C2.n && C1.k <= C2.k) || return false

    G1 = generator_matrix(C1) # Safe getter
    return all(view(G1, r:r, 1:C1.n) ∈ C2 for r in axes(G1, 1))
end
⊂(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊆ C2
is_subcode(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊆ C2

# ==============================================================================
# 14. DUAL CODES (Cache-Aware & Lazy)
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the (Euclidean) dual of the code `C`. Evaluates lazily by swapping the 
cached generator and parity-check matrices.
"""
function dual(C::AbstractLinearCode)
    if typeof(C) <: AbstractCyclicCode
        return CyclicCode(Int(order(C.F)), C.n, dual_qcosets(Int(order(C.F)), C.n, C.qcosets))
        
    elseif isa(C, GeneralizedReedSolomonCode)
        d = C.k + 1
        new_cache = Dict{Symbol, Any}()
        old_cache = getfield(C, :cache)
        haskey(old_cache, :G) && (new_cache[:H] = old_cache[:G])
        haskey(old_cache, :H) && (new_cache[:G] = old_cache[:H])
        
        return GeneralizedReedSolomonCode(C.F, C.n, C.n - C.k, d, d, d,
            deepcopy(C.dual_scalars), deepcopy(C.scalars), deepcopy(C.eval_pts), new_cache)
            
    elseif isa(C, MatrixProductCode)
        nr, nc = size(C.A)
        nr == nc || return LinearCode.dual(C) 
        
        D = Vector{LinearCode}()
        for i in 1:length(C.Cvec)
            push!(D, dual(C.Cvec[i]))
        end
        
        try
            A_inv = inv(C.A)
            return MatrixProductCode(D, transpose(A_inv))
        catch
            return LinearCode.dual(C) 
        end
        
    elseif isa(C, ReedMullerCode)
        d = 2^(C.r + 1)
        new_cache = Dict{Symbol, Any}()
        old_cache = getfield(C, :cache)
        haskey(old_cache, :G) && (new_cache[:H] = old_cache[:G])
        haskey(old_cache, :H) && (new_cache[:G] = old_cache[:H])
        
        return ReedMullerCode(C.F, C.n, C.n - C.k, d, d, d, C.m - C.r - 1, C.m, new_cache)
        
    else
        # Standard LinearCode Fallback - The Lazy Swap
        old_cache = getfield(C, :cache)
        new_cache = Dict{Symbol, Any}()
        
        haskey(old_cache, :G) && (new_cache[:H] = old_cache[:G])
        haskey(old_cache, :H) && (new_cache[:G] = old_cache[:H])
        
        d_new, l_new, u_new = missing, 1, C.n
        
        # Handle Weight Enumerator MacWilliams Identity securely without polynomials
        if haskey(old_cache, :weight_enum)
            dual_wt_enum = MacWilliams_transform(old_cache[:weight_enum], C.k, Int(order(C.F)))
            new_cache[:weight_enum] = dual_wt_enum
            
            non_zero_wts = filter(>(0), collect(keys(dual_wt_enum.counts)))
            if !isempty(non_zero_wts)
                d_val = minimum(non_zero_wts)
                d_new, l_new, u_new = d_val, d_val, d_val
            end
        else
            ub1 = haskey(new_cache, :G) ? _min_wt_row(new_cache[:G])[1] : C.n
            ub2 = haskey(new_cache, :H) ? _min_wt_row(new_cache[:H])[1] : C.n
            u_new = min(ub1, ub2)
        end

        return LinearCode(C.F, C.n, C.n - C.k, d_new, l_new, u_new, new_cache)
    end
end
Euclidean_dual(C::AbstractLinearCode) = dual(C)

# ==============================================================================
# 15. HERMITIAN DUAL
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the Hermitian dual of a code defined over a quadratic extension.
"""
function Hermitian_dual(C::AbstractLinearCode)
    if isa(C, MatrixProductCode)
        # the inner functions here should complain if not quadratic, so don't have to check here
        nr, nc = size(C.A)
        # Fallback to standard construction if non-square
        nr == nc || return dual(LinearCode(Hermitian_conjugate_matrix(generator_matrix(C))))
        
        D = Vector{LinearCode}()
        for i in 1:length(C.Cvec)
            push!(D, Hermitian_dual(C.Cvec[i]))
        end

        try
            A_inv = inv(Hermitian_conjugate_matrix(C.A))
            return MatrixProductCode(D, transpose(A_inv))
        catch
            # Fallback if matrix inversion fails
            return dual(LinearCode(Hermitian_conjugate_matrix(generator_matrix(C))))
        end
    else
        return dual(LinearCode(Hermitian_conjugate_matrix(generator_matrix(C))))
    end
end

# ==============================================================================
# 16. EQUIVALENCE CHECKERS
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return `true` if `C1 ⊆ C2` and `C2 ⊆ C1`.
"""
are_equivalent(C1::AbstractLinearCode, C2::AbstractLinearCode) = (C1 ⊆ C2) && (C2 ⊆ C1)

"""
$(TYPEDSIGNATURES)

Return `true` if `are_equivalent(C, dual(C))`.
"""
is_self_dual(C::AbstractLinearCode) = are_equivalent(C, dual(C))

"""
$(TYPEDSIGNATURES)

Return `true` if `C ⊆ dual(C)`.
"""
is_self_orthogonal(C::AbstractLinearCode) = C ⊆ dual(C)

"""
$(TYPEDSIGNATURES)

Return `true` if `C ⊆ dual(C)`.
"""
is_weakly_self_dual(C::AbstractLinearCode) = is_self_orthogonal(C)

"""
$(TYPEDSIGNATURES)

Return `true` if `C ⊆ Euclidean_dual(C)`.
"""
is_Euclidean_self_orthogonal(C::AbstractLinearCode) = is_self_orthogonal(C)

"""
$(TYPEDSIGNATURES)

Return `true` if `dual(C) ⊆ C`.
"""
is_dual_containing(C::AbstractLinearCode) = dual(C) ⊆ C

"""
$(TYPEDSIGNATURES)

Return `true` if `Euclidean_dual(C) ⊆ C`.
"""
is_Euclidean_dual_containing(C::AbstractLinearCode) = is_dual_containing(C)

"""
$(TYPEDSIGNATURES)

Return `true` if `are_equivalent(C, Hermitian_dual(C))`.
"""
is_Hermitian_self_dual(C::AbstractLinearCode) = are_equivalent(C, Hermitian_dual(C))

"""
$(TYPEDSIGNATURES)

Return `true` if `C ⊆ Hermitian_dual(C)`.
"""
is_Hermitian_self_orthogonal(C::AbstractLinearCode) = C ⊆ Hermitian_dual(C)

"""
$(TYPEDSIGNATURES)

Return `true` if `C ⊆ Hermitian_dual(C)`.
"""
is_Hermitian_weakly_self_dual(C::AbstractLinearCode) = is_Hermitian_self_orthogonal(C)

"""
$(TYPEDSIGNATURES)

Return `true` if `Hermitian_dual(C) ⊆ C`.
"""
is_Hermitian_dual_containing(C::AbstractLinearCode) = Hermitian_dual(C) ⊆ C

# ==============================================================================
# 15.5 l-GALOIS DUAL & PROPERTIES
# ==============================================================================

"""
Internal helper function to apply the `l`-th power of the Frobenius automorphism
to every element of a matrix.
"""
function _l_Galois_conjugate_matrix(M::CTMatrixTypes, l::Int)
    F = base_ring(M)
    p = characteristic(F)
    power = p^l
    
    nr, nc = nrows(M), ncols(M)
    M_conj = zero_matrix(F, nr, nc)
    
    for r in 1:nr, c in 1:nc
        if !iszero(M[r, c])
            M_conj[r, c] = M[r, c]^power
        end
    end
    return M_conj
end

"""
$(TYPEDSIGNATURES)

Return the `l`-Galois dual of a code defined over a finite field.
"""
function l_Galois_dual(C::AbstractLinearCode, l::Int)
    if isa(C, MatrixProductCode)
        nr, nc = size(C.A)
        # Fallback to standard construction if non-square
        nr == nc || return dual(LinearCode(_l_Galois_conjugate_matrix(generator_matrix(C), l)))
        
        D = Vector{LinearCode}()
        for i in 1:length(C.Cvec)
            push!(D, l_Galois_dual(C.Cvec[i], l))
        end

        try
            A_inv = inv(_l_Galois_conjugate_matrix(C.A, l))
            return MatrixProductCode(D, transpose(A_inv))
        catch
            # Fallback if matrix inversion fails
            return dual(LinearCode(_l_Galois_conjugate_matrix(generator_matrix(C), l)))
        end
    else
        return dual(LinearCode(_l_Galois_conjugate_matrix(generator_matrix(C), l)))
    end
end

"""
$(TYPEDSIGNATURES)

Return `true` if `are_equivalent(C, l_Galois_dual(C, l))`.
"""
is_l_Galois_self_dual(C::AbstractLinearCode, l::Int) = are_equivalent(C, l_Galois_dual(C, l))

"""
$(TYPEDSIGNATURES)

Return `true` if `C ⊆ l_Galois_dual(C, l)`.
"""
is_l_Galois_self_orthogonal(C::AbstractLinearCode, l::Int) = C ⊆ l_Galois_dual(C, l)

"""
$(TYPEDSIGNATURES)

Return `true` if `C ⊆ l_Galois_dual(C, l)`.
"""
is_l_Galois_weakly_self_dual(C::AbstractLinearCode, l::Int) = is_l_Galois_self_orthogonal(C, l)

"""
$(TYPEDSIGNATURES)

Return `true` if `l_Galois_dual(C, l) ⊆ C`.
"""
is_l_Galois_dual_containing(C::AbstractLinearCode, l::Int) = l_Galois_dual(C, l) ⊆ C

"""
$(TYPEDSIGNATURES)

Return the `l`-Galois hull of `C` and its dimension.
"""
function l_Galois_hull(C::AbstractLinearCode, l::Int)
    D = l_Galois_dual(C, l)
    G = generator_matrix(C)
    H = parity_check_matrix(D)
    F = field(C)
    VS = vector_space(F, C.n)
    
    U, U_to_VS = sub(VS, [VS(G[i, :]) for i in 1:nrows(G)])
    W, _ = sub(VS, [VS(H[i, :]) for i in 1:nrows(H)])
    I, I_to_W = intersect(U, W)
    
    if !iszero(AbstractAlgebra.dim(I))
        I_basis = [U_to_VS(I_to_W(g)) for g in gens(I)]
        G_I = reduce(vcat, I_basis)
        F_basis = [[F(G_I[j][i]) for i in 1:C.n] for j in 1:AbstractAlgebra.dim(I)]
        G_hull = matrix(F, length(F_basis), length(F_basis[1]), reduce(vcat, F_basis))
        return LinearCode(G_hull), AbstractAlgebra.dim(I)
    else
        return missing, 0 
    end
end

"""
$(TYPEDSIGNATURES)

Return `true` if `C` is linear complementary `l`-Galois dual 
(i.e., the dimension of `l_Galois_hull(C, l)` is zero).
"""
function is_l_Galois_LCD(C::AbstractLinearCode, l::Int)
    _, dim_hull_C = l_Galois_hull(C, l)
    return dim_hull_C == 0
end

"""
$(TYPEDSIGNATURES)

Return the characteristic polynomial of `C`.
"""
function characteristic_polynomial(C::AbstractLinearCode)
    _, x = polynomial_ring(Nemo.QQ, :x)
    D = dual(C)
    supD = support(D)
    q = Int(order(C.F))
    # BUG FIX: n and k were undefined
    return q^(C.n - C.k) * prod(1 - x / j for j in supD if j > 0)
end

"""
$(TYPEDSIGNATURES)

Return the code `C` as a vector space object.
"""
function vector_space(C::AbstractLinearCode)
    V = AbstractAlgebra.vector_space(C.F, C.n)
    G = generator_matrix(C) # Safe getter
    return sub(V, [V(view(G, i:i, :)) for i in 1:nrows(G)])
end

"""
$(TYPEDSIGNATURES)

Return `true` if `C` is even.
"""
function is_even(C::AbstractLinearCode)
    Int(order(C.F)) == 2 || throw(ArgumentError("Even-ness is only defined for binary codes."))
    G = generator_matrix(C) # Safe getter
    return all(wt(view(G, r:r, :)) % 2 == 0 for r in 1:nrows(G))
end

"""
$(TYPEDSIGNATURES)

Return `true` if `C` is doubly-even.
"""
function is_doubly_even(C::AbstractLinearCode)
    Int(order(C.F)) == 2 || throw(ArgumentError("Even-ness is only defined for binary codes."))
    G = generator_matrix(C) # Safe getter
    nr = nrows(G)
    all(wt(view(G, r:r, :)) % 4 == 0 for r in 1:nr) || (return false)
    all(wt(view(G, r1:r1, :) + view(G, r2:r2, :)) % 4 == 0 for r1 in 1:nr, r2 in 1:nr) || (return false)
    return true
end

"""
$(TYPEDSIGNATURES)

Return `true` if `C` is triply-even.
"""
function is_triply_even(C::AbstractLinearCode)
    Int(order(C.F)) == 2 || throw(ArgumentError("Even-ness is only defined for binary codes."))
    G = _Flint_matrix_to_Julia_int_matrix(generator_matrix(C)) # Safe getter
    nr, _ = size(G)
    all(wt(view(G, r:r, :)) % 8 == 0 for r in 1:nr) || (return false)
    all(wt(view(G, r1:r1, :) .* view(G, r2:r2, :)) % 4 == 0 for r1 in 1:nr, r2 in 1:nr) || (return false)
    all(wt(view(G, r1:r1, :) .* view(G, r2:r2, :) .* view(G, r3:r3, :)) % 2 == 0 for r1 in 1:nr, r2 in 1:nr, r3 in 1:nr) || (return false)
    return true
end

"""
$(TYPEDSIGNATURES)

Return the elements of `C`. If `only_print` is `true`, the elements are only 
printed to the console and not returned.
"""
function words(C::AbstractLinearCode, only_print::Bool = false)
    G_actual = generator_matrix(C) # Force evaluation if missing
    words_vec = only_print ? nothing : Vector{typeof(G_actual)}()
    
    P_stand = standard_form_permutation(C) # Safe getter
    G = ismissing(P_stand) ? generator_matrix(C, true) : generator_matrix(C, true) * P_stand
    E = base_ring(G)

    if iszero(G)
        row = zero_matrix(E, 1, C.n)
        only_print ? println(row) : push!(words_vec, row)
        return words_vec
    end

    for iter in Nemo.AbstractAlgebra.ProductIterator([E for _ in 1:nrows(G)], inplace = true)
        row = iter[1] * view(G, 1:1, :)
        for r in 2:nrows(G)
            if !iszero(iter[r])
                row += iter[r] * view(G, r:r, :)
            end
        end
        only_print ? println(row) : push!(words_vec, row)
    end
    return words_vec
end
codewords(C::AbstractLinearCode, only_print::Bool = false) = words(C, only_print)
elements(C::AbstractLinearCode, only_print::Bool = false) = words(C, only_print)

"""
$(TYPEDSIGNATURES)

Return the (Euclidean) hull of `C` and its dimension.
"""
function hull(C::AbstractLinearCode)
    G = generator_matrix(C)
    H = parity_check_matrix(C)
    F = field(C)
    VS = AbstractAlgebra.vector_space(F, C.n)
    U, U_to_VS = sub(VS, [VS(view(G, i:i, :)) for i in 1:nrows(G)])
    W, _ = sub(VS, [VS(view(H, i:i, :)) for i in 1:nrows(H)])
    I, I_to_W = intersect(U, W)
    
    if !iszero(AbstractAlgebra.dim(I))
        I_basis = [U_to_VS(I_to_W(g)) for g in gens(I)]
        G_I = reduce(vcat, I_basis)
        F_basis = [[F(G_I[j][i]) for i in 1:C.n] for j in 1:AbstractAlgebra.dim(I)]
        G_hull = matrix(F, length(F_basis), length(F_basis[1]), reduce(vcat, F_basis))
        return LinearCode(G_hull), AbstractAlgebra.dim(I)
    else
        return missing, 0 
    end
end
Euclidean_hull(C::AbstractLinearCode) = hull(C)

"""
$(TYPEDSIGNATURES)

Return the Hermitian hull of `C` and its dimension.
The Hermitian hull of a code is the intersection of it and its Hermitian dual.
"""
function Hermitian_hull(C::AbstractLinearCode)
    D = Hermitian_dual(C)
    G = generator_matrix(C) # Safe getter
    H = parity_check_matrix(D) # Safe getter
    F = field(C)
    VS = AbstractAlgebra.vector_space(F, C.n)
    
    # Use view() to prevent row-copy allocations
    U, U_to_VS = sub(VS, [VS(view(G, i:i, :)) for i in 1:nrows(G)])
    W, _ = sub(VS, [VS(view(H, i:i, :)) for i in 1:nrows(H)])
    I, I_to_W = intersect(U, W)
    
    if !iszero(AbstractAlgebra.dim(I))
        I_basis = [U_to_VS(I_to_W(g)) for g in gens(I)]
        G_I = reduce(vcat, I_basis)
        F_basis = [[F(G_I[j][i]) for i in 1:C.n] for j in 1:AbstractAlgebra.dim(I)]
        G_hull = matrix(F, length(F_basis), length(F_basis[1]), reduce(vcat, F_basis))
        return LinearCode(G_hull), AbstractAlgebra.dim(I)
    else
        return missing, 0 
    end
end

"""
$(TYPEDSIGNATURES)

Return `true` if `C` is linear complementary dual.
"""
function is_LCD(C::AbstractLinearCode)
    _, dim_hull_C = hull(C)
    return dim_hull_C == 0
end

"""
$(TYPEDSIGNATURES)

Return `true` if `C` is linear complementary Hermitian dual.
A code is linear complementary Hermitian dual if the dimension of `Hermitian_hull(C)` is zero.
"""
function is_Hermitian_LCD(C::AbstractLinearCode)
    _, dim_hull_C = Hermitian_hull(C)
    return dim_hull_C == 0
end

"""
$(TYPEDSIGNATURES)

Return `true` if `C` contains a self-dual subcode.
"""
function contains_self_dual_subcode(C::AbstractLinearCode)
    q = Int(order(C.F))
    factors = Nemo.factor(q)
    (p, t) = first(factors) # BUG FIX: Safer extraction from the factor dictionary

    if p == 2
        return iseven(C.n) && is_self_orthogonal(dual(C))
    elseif q % 4 == 1 
        return iseven(C.n) && is_self_orthogonal(dual(C))
    elseif q % 4 == 3 
        return C.n % 4 == 0 && is_self_orthogonal(dual(C))
    else
        error("Unknown case for finite field of order $(p)^$(t)")
    end
end
