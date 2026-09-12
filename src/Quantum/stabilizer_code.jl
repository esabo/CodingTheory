# Copyright (c) 2021 - 2026 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

function _css_symplectic_matrix(
    X::CTMatrixTypes, Z::CTMatrixTypes, is_sparse::Bool
)
    if X isa SparseMatrixCSC && Z isa SparseMatrixCSC &&
       eltype(X) <: Integer && eltype(Z) <: Integer
        T = promote_type(eltype(X), eltype(Z))
        return vcat(
            hcat(X, spzeros(T, nrows(X), ncols(Z))),
            hcat(spzeros(T, nrows(Z), ncols(X)), Z)
        )
    elseif X isa SparseMatrixCSC && Z isa SparseMatrixCSC
        F = _code_matrix_base_ring(
            nrows(X) > 0 && !iszero(X) ? X : Z)
        return _sparse_code_matrix(direct_sum(
            _dense_code_matrix(X, F), _dense_code_matrix(Z, F)))
    end
    return direct_sum(X, Z)
end

"""
    StabilizerCodeCSS(X_matrix::CTMatrixTypes, Z_matrix::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)
    CSSCode(X_matrix::CTMatrixTypes, Z_matrix::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing}= missing, logs_alg::Symbol = :stnd_frm)

Return a CSS code whose `X`-stabilizers are given by `X_matrix`, `Z`-stabilizers by `Z_matrix`.
"""
function StabilizerCodeCSS(X_matrix::T, Z_matrix::T; char_vec::Union{Vector{zzModRingElem}, Missing} = 
    missing, logs_alg::Symbol = :stnd_frm) where T <: CTMatrixTypes

    is_sparse = _is_sparse_code_matrix(X_matrix)
    X_matrix = _normalize_quantum_matrix(X_matrix)
    Z_matrix = _normalize_quantum_matrix(Z_matrix)
    logs_alg ∈ (:stnd_frm, :sys_eqs) || throw(ArgumentError("Unrecognized logicals algorithm. Use :stnd_frm or :sys_eqs."))
    n = ncols(X_matrix)
    n > 0 || throw(ArgumentError("The stabilizer matrices must have positive length."))
    n == ncols(Z_matrix) || throw(ArgumentError("Both matrices must have the same length in the CSS construction."))
    F = _code_matrix_base_ring(X_matrix)
    F == _code_matrix_base_ring(Z_matrix) || throw(ArgumentError("Both matrices must be over the same base field."))
    
    X_clean = _remove_empty(deepcopy(X_matrix), :rows)
    Z_clean = _remove_empty(deepcopy(Z_matrix), :rows)
    if is_sparse
        X_clean = _sparse_code_matrix(X_clean)
        Z_clean = _sparse_code_matrix(Z_clean)
    end
    stabs = _css_symplectic_matrix(X_clean, Z_clean, is_sparse)
    are_symplectic_orthogonal(stabs, stabs) ||
        throw(ArgumentError("The given matrices are not symplectic orthogonal."))

    p = Int(characteristic(F))
    clean_char_vec = ismissing(char_vec) ? zzModRingElem[] : _process_char_vec(char_vec, p, 2 * n)

    X_rank = _additive_rank(X_clean, F)
    Z_rank = _additive_rank(Z_clean, F)
    rnk = X_rank + Z_rank

    dim_code = _quantum_dimension(F, n, rnk)
    
    over_comp = (nrows(X_clean) > X_rank) || (nrows(Z_clean) > Z_rank)
    
    X_final = X_clean
    Z_final = Z_clean

    cache = Dict{Symbol, Any}(
        :stabs => stabs,
        :overcomplete => over_comp,
        :logs_alg => degree(F) == 1 ? logs_alg : :sys_eqs
    )

    result = StabilizerCodeCSS(
        F, n, dim_code, X_final, Z_final, clean_char_vec, cache)
    return _seed_quantum_singleton_bound!(result)
end

"""
$(TYPEDSIGNATURES)

Return the CSS stabilizer code whose trimmed `X`- and `Z`-stabilizer matrices
are `X_matrix` and `Z_matrix`. This is an alias for `StabilizerCodeCSS`.
"""
CSSCode(X_matrix::T, Z_matrix::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm) where T <: CTMatrixTypes = StabilizerCodeCSS(X_matrix, Z_matrix,
    char_vec = char_vec, logs_alg = logs_alg)

"""
$(TYPEDSIGNATURES)

Return the stabilizer code whose stabilizers is determined by `stabs`.
"""
function StabilizerCode(stabs::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm)

    is_sparse = _is_sparse_code_matrix(stabs)
    stabs = _normalize_quantum_matrix(stabs)
    logs_alg ∈ (:stnd_frm, :sys_eqs) || throw(ArgumentError("Unrecognized logicals algorithm. Use :stnd_frm or :sys_eqs."))
    F = _code_matrix_base_ring(stabs)
    p = Int(characteristic(F))
    iseven(ncols(stabs)) ||
        throw(ArgumentError("A symplectic stabilizer matrix must have an even number of columns."))
    n = div(ncols(stabs), 2)
    n > 0 || throw(ArgumentError("The stabilizer matrix must have positive length."))
    clean_char_vec = ismissing(char_vec) ? zzModRingElem[] : _process_char_vec(char_vec, p, 2 * n)
    
    stabs_final = _remove_empty(deepcopy(stabs), :rows)
    is_sparse && (stabs_final = _sparse_code_matrix(stabs_final))
    are_symplectic_orthogonal(stabs_final, stabs_final) || throw(ArgumentError("The given stabilizers are not symplectic orthogonal."))
    
    rnk = _additive_rank(stabs_final, F)
    dim_code = _quantum_dimension(F, n, rnk)
    over_comp = nrows(stabs_final) > rnk
    
    is_css_S, X_stabs_dense, Z_stabs_dense = if rnk == 0
        (true, zero_matrix(F, 0, n), zero_matrix(F, 0, n))
    else
        robust_CSS_split(_dense_code_matrix(stabs_final, F))
    end

    cache = Dict{Symbol, Any}(
        :stabs => stabs_final,
        :overcomplete => over_comp,
        :logs_alg => degree(F) == 1 ? logs_alg : :sys_eqs
    )
    
    if is_css_S
        X_stabs = is_sparse ? _sparse_code_matrix(X_stabs_dense) : X_stabs_dense
        Z_stabs = is_sparse ? _sparse_code_matrix(Z_stabs_dense) : Z_stabs_dense
        result = StabilizerCodeCSS(
            F, n, dim_code, X_stabs, Z_stabs, clean_char_vec, cache)
        return _seed_quantum_singleton_bound!(result)
    else
        result = StabilizerCode(
            F, n, dim_code, stabs_final, clean_char_vec, cache)
        return _seed_quantum_singleton_bound!(result)
    end
end

"""
    StabilizerCodeCSS(C1::AbstractLinearCode, C2::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)
    CSSCode(C1::AbstractLinearCode, C2::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)
"""
function StabilizerCodeCSS(C1::AbstractLinearCode, C2::AbstractLinearCode;
    char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

    C2 ⊆ C1 || throw(ArgumentError("The second argument must be a subset of the first in the CSS construction."))
    D2 = dual(C2)
    
    S = StabilizerCodeCSS(D2.H, C1.H, char_vec=char_vec, logs_alg=logs_alg)
    
    S.cache[:X_orig_code] = D2
    S.cache[:Z_orig_code] = C1
    
    S.cache[:l_bound_dx] = D2.l_bound
    S.cache[:l_bound_dz] = C1.l_bound
    S.cache[:l_bound] = min(C1.l_bound, D2.l_bound)
    return S
end
CSSCode(C1::AbstractLinearCode, C2::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) = StabilizerCodeCSS(C1, C2, char_vec = char_vec, logs_alg = logs_alg)

"""
    StabilizerCodeCSS(C::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)
    CSSCode(C::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)
"""
function StabilizerCodeCSS(C::LinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm)

    D = dual(C)
    C ⊆ D || throw(ArgumentError("The single code CSS construction requires C ⊆ C^⟂."))
    
    S = StabilizerCodeCSS(D.H, D.H, char_vec=char_vec, logs_alg=logs_alg)
    S.cache[:X_orig_code] = D
    S.cache[:Z_orig_code] = D
    
    S.cache[:l_bound_dx] = D.l_bound
    S.cache[:l_bound_dz] = D.l_bound
    S.cache[:l_bound] = D.l_bound
    return S
end
CSSCode(C::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) = StabilizerCodeCSS(C, char_vec = char_vec, logs_alg = logs_alg)

function _quadratic_code_to_symplectic(
    C::AbstractLinearCode, K::CTFieldTypes,
    basis::Union{Missing, Vector{<:CTFieldElem}}=missing
)
    E = C.F
    order(E) == order(K)^2 ||
        throw(ArgumentError("The classical code must be over a quadratic extension of the requested symplectic field."))
    β = ismissing(basis) ? first(primitive_basis(E, K)) : basis
    length(β) == 2 ||
        throw(ArgumentError("A quadratic extension basis must contain two elements."))
    is_basis(E, K, β)[1] ||
        throw(ArgumentError("The supplied elements are not a basis of the quadratic extension."))

    G = generator_matrix(C)
    additive_generators = vcat(β[1] * G, β[2] * G)
    expanded = expand_matrix(additive_generators, K, β)
    n = C.n
    symplectic = hcat(expanded[:, 1:2:2n], expanded[:, 2:2:2n])
    prime_field = _prime_subfield(K)
    prime_basis = degree(K) == 1 ?
        [one(K)] : first(primitive_basis(K, prime_field))
    return reduce(vcat, [α * symplectic for α in prime_basis])
end

"""
$(TYPEDSIGNATURES)

Return the symplectic stabilizer code over `F` associated with a Hermitian
self-orthogonal linear code over the quadratic extension of `F`. The optional
`basis` is an ordered extension basis; a primitive basis is used by default.
"""
function StabilizerCode(
    C::AbstractLinearCode, F::CTFieldTypes;
    basis::Union{Missing, Vector{<:CTFieldElem}}=missing,
    char_vec::Union{Vector{zzModRingElem}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm
)
    is_Hermitian_self_orthogonal(C) ||
        throw(ArgumentError("The classical code must be Hermitian self-orthogonal."))
    stabs = _quadratic_code_to_symplectic(C, F, basis)
    S = StabilizerCode(stabs; char_vec=char_vec, logs_alg=logs_alg)
    S.cache[:quadratic_code] = C
    S.cache[:quadratic_basis] =
        ismissing(basis) ? first(primitive_basis(C.F, F)) : basis
    return S
end

"""
$(TYPEDSIGNATURES)

"""
function StabilizerCodeCSS(S_Pauli::Vector{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm) where T <: Union{String, Vector{Char}}

    stabs = _Pauli_string_to_symplectic(_process_strings(S_Pauli))
    iszero(stabs) && throw(ArgumentError("The processed Pauli strings returned a set of empty stabilizer generators."))
    stabs = _remove_empty(stabs, :rows)
    are_symplectic_orthogonal(stabs, stabs) || throw(ArgumentError("The given stabilizers are not symplectic orthogonal."))
    
    return StabilizerCode(stabs, char_vec = char_vec, logs_alg = logs_alg)
end
CSSCode(S_Pauli::Vector{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{String, Vector{Char}} = StabilizerCodeCSS(S_Pauli, char_vec = char_vec, logs_alg = logs_alg)

"""
$(TYPEDSIGNATURES)

"""
function StabilizerCode(S_Pauli::Vector{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm) where T <: Union{String, Vector{Char}}

    S_Pauli_stripped = _process_strings(S_Pauli)
    stabs = _Pauli_string_to_symplectic(S_Pauli_stripped)
    iszero(stabs) && throw(ArgumentError("The processed Pauli strings returned a set of empty stabilizer generators."))
    return StabilizerCode(stabs, char_vec = char_vec, logs_alg = logs_alg)
end

"""
$(TYPEDSIGNATURES)

"""
function StabilizerCodeCSS(S::AbstractStabilizerCode; logs_alg::Symbol = :stnd_frm)
    Z = deepcopy(stabilizers(S))
    Z[:, 1:S.n] = stabilizers(S)[:, S.n + 1:2 * S.n]
    Z[:, S.n + 1:2 * S.n] = -stabilizers(S)[:, 1:S.n]
    return CSSCode(stabilizers(S), Z, logs_alg = logs_alg)
end
CSSCode(S::AbstractStabilizerCode; logs_alg::Symbol = :stnd_frm) = StabilizerCodeCSS(S; logs_alg = logs_alg)

"""
$(TYPEDSIGNATURES)

"""
StabilizerCode(S::T; logs_alg::Symbol = :stnd_frm) where {T <: AbstractStabilizerCode} = StabilizerCode(CSSTrait(T), S, logs_alg)
function StabilizerCode(::IsCSS, S::AbstractStabilizerCode, logs_alg::Symbol)
    iseven(S.n) || throw(ArgumentError("Only valid for codes of even length"))
    nrows(S.X_stabs) == nrows(S.Z_stabs) || throw(ArgumentError("Requires an equal number of X and Z stabilizers"))
    n = div(S.n, 2)
    return StabilizerCode(hcat(S.X_stabs[:, 1:n], S.Z_stabs[:, 1:n]), logs_alg = logs_alg)
end
StabilizerCode(::IsNotCSS, S::AbstractStabilizerCode, logs_alg::Symbol) = throw(ArgumentError("Only valid for CSS codes of even length"))

"""
$(TYPEDSIGNATURES)

"""
function StabilizerCode(S::AbstractSubsystemCode)
    typeof(S) <: AbstractStabilizerCode && return S
    # For a subsystem code, promoting all gauges to logicals yields the underlying stabilizer code
    return promote_gauges_to_logical(S, collect(1:length(gauges(S))))
end

#############################
      # getter functions
#############################

"""
$(TYPEDSIGNATURES)

Return the currently stored lower bound on the minimum distance.
"""
minimum_distance_lower_bound(S::AbstractStabilizerCode) = get(S.cache, :l_bound, missing)

"""
$(TYPEDSIGNATURES)

Return the currently stored upper bound on the minimum distance.
"""
minimum_distance_upper_bound(S::AbstractStabilizerCode) = get(S.cache, :u_bound, missing)

"""
$(TYPEDSIGNATURES)

Return the currently stored lower bound on the minimum `X`-distance.
"""
X_minimum_distance_lower_bound(S::T) where T <: AbstractStabilizerCode = X_minimum_distance_lower_bound(CSSTrait(T), S)
X_minimum_distance_lower_bound(::IsCSS, S::AbstractStabilizerCode) = get(S.cache, :l_bound_dx, missing)
X_minimum_distance_lower_bound(::IsNotCSS, S::AbstractStabilizerCode) = error("Only valid for CSS codes")

"""
$(TYPEDSIGNATURES)

Return the currently stored upper bound on the minimum `X`-distance.
"""
X_minimum_distance_upper_bound(S::T) where T <: AbstractStabilizerCode = X_minimum_distance_upper_bound(CSSTrait(T), S)
X_minimum_distance_upper_bound(::IsCSS, S::AbstractStabilizerCode) = get(S.cache, :u_bound_dx, missing)
X_minimum_distance_upper_bound(::IsNotCSS, S::AbstractStabilizerCode) = error("Only valid for CSS codes")

"""
$(TYPEDSIGNATURES)

Return the currently stored lower bound on the minimum `Z`-distance.
"""
Z_minimum_distance_lower_bound(S::T) where T <: AbstractStabilizerCode = Z_minimum_distance_lower_bound(CSSTrait(T), S)
Z_minimum_distance_lower_bound(::IsCSS, S::AbstractStabilizerCode) = get(S.cache, :l_bound_dz, missing)
Z_minimum_distance_lower_bound(::IsNotCSS, S::AbstractStabilizerCode) = error("Only valid for CSS codes")

"""
$(TYPEDSIGNATURES)

Return the currently stored upper bound on the minimum `Z`-distance.
"""
Z_minimum_distance_upper_bound(S::T) where T <: AbstractStabilizerCode = Z_minimum_distance_upper_bound(CSSTrait(T), S)
Z_minimum_distance_upper_bound(::IsCSS, S::AbstractStabilizerCode) = get(S.cache, :u_bound_dz, missing)
Z_minimum_distance_upper_bound(::IsNotCSS, S::AbstractStabilizerCode) = error("Only valid for CSS codes")


#############################
      # setter functions
#############################

"""
$(TYPEDSIGNATURES)

Set the minimum distance of the code to `d`.
"""
function set_minimum_distance!(S::AbstractStabilizerCode, d::Int)
    0 < d <= S.n || throw(DomainError("The minimum distance of a code must be ≥ 1; received: d = $d."))
    if dimension(S) > 0
        singleton = quantum_Singleton_bound(S)
        d <= singleton ||
            throw(DomainError(d,
                "The distance exceeds the quantum Singleton bound $singleton."))
    end
    
    curr_u_bound = get(S.cache, :u_bound, S.n)
    curr_l_bound = get(S.cache, :l_bound, 1)
    
    curr_u_bound < d && (@warn "The distance set is greater than the current upper bound of $curr_u_bound")
    d < curr_l_bound && (@warn "The distance set is less than the current lower bound of $curr_l_bound")
    
    S.cache[:d] = d
    S.cache[:l_bound] = d
    S.cache[:u_bound] = d

    if CSSTrait(typeof(S)) == IsCSS()
        !haskey(S.cache, :dx) && (S.cache[:l_bound_dx] = d)
        !haskey(S.cache, :dz) && (S.cache[:l_bound_dz] = d)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Set the minimum `X`-distance of the code to `d`.
"""
set_X_minimum_distance!(S::T, d::Int) where T <: AbstractStabilizerCode = set_X_minimum_distance!(CSSTrait(T), S, d)
function set_X_minimum_distance!(::IsCSS, S::AbstractStabilizerCode, d::Int)
    0 < d <= S.n || throw(DomainError("The minimum distance of a code must be ≥ 1; received: d = $d."))
    
    curr_u_bound_dx = get(S.cache, :u_bound_dx, S.n)
    curr_l_bound_dx = get(S.cache, :l_bound_dx, 1)
    
    curr_u_bound_dx < d && (@warn "The distance set is greater than the current upper bound of $curr_u_bound_dx")
    d < curr_l_bound_dx && (@warn "The distance set is less than the current lower bound of $curr_l_bound_dx")
    
    S.cache[:dx] = d
    S.cache[:l_bound_dx] = d
    S.cache[:u_bound_dx] = d

    if haskey(S.cache, :dz)
        S.cache[:d] = min(d, S.cache[:dz])
        S.cache[:u_bound] = S.cache[:d]
        S.cache[:l_bound] = S.cache[:d]
    else
        curr_u_bound = get(S.cache, :u_bound, S.n)
        d < curr_u_bound && (S.cache[:u_bound] = d)
    end
    return nothing
end
set_X_minimum_distance!(::IsNotCSS, S::AbstractStabilizerCode, d::Int) =
    error("Only valid for CSS codes")

"""
$(TYPEDSIGNATURES)

Set the minimum `Z`-distance of the code to `d`.
"""
set_Z_minimum_distance!(S::T, d::Int) where T <: AbstractStabilizerCode = set_Z_minimum_distance!(CSSTrait(T), S, d)
function set_Z_minimum_distance!(::IsCSS, S::AbstractStabilizerCode, d::Int)
    0 < d <= S.n || throw(DomainError("The minimum distance of a code must be ≥ 1; received: d = $d."))
    
    curr_u_bound_dz = get(S.cache, :u_bound_dz, S.n)
    curr_l_bound_dz = get(S.cache, :l_bound_dz, 1)
    
    curr_u_bound_dz < d && (@warn "The distance set is greater than the current upper bound of $curr_u_bound_dz")
    d < curr_l_bound_dz && (@warn "The distance set is less than the current lower bound of $curr_l_bound_dz")
    
    S.cache[:dz] = d
    S.cache[:l_bound_dz] = d
    S.cache[:u_bound_dz] = d

    if haskey(S.cache, :dx)
        S.cache[:d] = min(S.cache[:dx], d)
        S.cache[:u_bound] = S.cache[:d]
        S.cache[:l_bound] = S.cache[:d]
    else
        curr_u_bound = get(S.cache, :u_bound, S.n)
        d < curr_u_bound && (S.cache[:u_bound] = d)
    end
    return nothing
end
set_Z_minimum_distance!(::IsNotCSS, S::AbstractStabilizerCode, d::Int) =
    error("Only valid for CSS codes")


#############################
     # general functions
#############################

function _logicals(
    stabs::CTMatrixTypes, dual_gens::CTMatrixTypes,
    logs_alg::Symbol = :sys_eqs
)
    logs_alg ∈ (:sys_eqs,) || throw(ArgumentError("Unrecognized logicals algorithm. Use :sys_eqs."))

    F = _code_matrix_base_ring(dual_gens)
    L = _additive_quotient_space(stabs, dual_gens, F)
    logs = _pair_operators(L, false)
    # verify
    n = div(ncols(L), 2)
    logs_mat = _pairs_to_matrix(
        logs, F, n, _is_sparse_code_matrix(L))
    are_symplectic_orthogonal(stabs, logs_mat) || error("Computed logicals do not commute with the codespace.")
    prod = _trace_symplectic_product_matrix(
        logs_mat, logs_mat, F)
    expected = zero_matrix(base_ring(prod), nrows(logs_mat), nrows(logs_mat))
    for i in 1:2:nrows(logs_mat)
        expected[i, i + 1] = one(base_ring(prod))
        expected[i + 1, i] = -one(base_ring(prod))
    end
    prod == expected ||
        error("Computed logicals do not have canonical commutation relations.")
    return logs, logs_mat
end

"""
$(TYPEDSIGNATURES)

Return a random CSS code with an equal number of `X` and `Z` stabilizers.
"""
function random_CSS_code(n::Int, k::Int)
    d = _rand_single_sector_boundary(n, k)
    dt = transpose(d)
    Xrref = rref(dt)
    Xind = [findfirst(!iszero, Xrref[2][i, :])[2] for i in 1:Xrref[1]]
    Zrref = rref(d)
    Zind = [findfirst(!iszero, Zrref[2][i, :])[2] for i in 1:Zrref[1]]
    return CSSCode(d[Xind, :], dt[Zind, :])
end

function _random_symplectic_pairs(rng::AbstractRNG, F::CTFieldTypes, n::Int)
    n > 0 || throw(DomainError(n, "The code length must be positive."))
    prime_field = _prime_subfield(F)
    additive_dimension = 2 * degree(F) * n
    coefficients = matrix(
        prime_field,
        rand(rng, prime_field, additive_dimension, additive_dimension)
    )
    while rank(coefficients) != additive_dimension
        coefficients = matrix(
            prime_field,
            rand(rng, prime_field, additive_dimension, additive_dimension)
        )
    end
    M = _lift_prime_matrix(coefficients, F) *
        _additive_ambient_basis(F, 2n)
    return _make_pairs(M)
end

"""
$(TYPEDSIGNATURES)

Return a random, not necessarily uniformly sampled, ``[[n, k]]`` stabilizer
code over `F`.
"""
function random_stabilizer_code(
    rng::AbstractRNG, F::CTFieldTypes, n::Int, k::Union{Int, Rational};
    char_vec::Union{Vector{zzModRingElem}, Missing}=missing
)
    0 <= k <= n || throw(DomainError(k, "Expected 0 ≤ k ≤ n."))
    num_stabs_rat = degree(F) * (n - k)
    denominator(Rational{BigInt}(num_stabs_rat)) == 1 ||
        throw(DomainError(k, "The requested additive dimension is incompatible with the field."))
    num_stabs = Int(num_stabs_rat)
    pairs = _random_symplectic_pairs(rng, F, n)
    stabs = num_stabs == 0 ? zero_matrix(F, 0, 2n) :
        reduce(vcat, [pairs[i][1] for i in 1:num_stabs])
    return StabilizerCode(stabs; char_vec=char_vec)
end
random_stabilizer_code(F::CTFieldTypes, n::Int, k::Union{Int, Rational}; kwargs...) =
    random_stabilizer_code(Random.default_rng(), F, n, k; kwargs...)
random_stabilizer_code(rng::AbstractRNG, n::Int, k::Union{Int, Rational}; kwargs...) =
    random_stabilizer_code(rng, Oscar.Nemo.Native.GF(2), n, k; kwargs...)
random_stabilizer_code(n::Int, k::Union{Int, Rational}; kwargs...) =
    random_stabilizer_code(Random.default_rng(), n, k; kwargs...)

"""
$(TYPEDSIGNATURES)

Return `true` if `S` is a CSS-T code.
"""
is_CSS_T_code(S::AbstractStabilizerCode; verbose::Bool = false) = is_CSS_T_code(CSSTrait(
    typeof(S)), S; verbose = verbose)
function is_CSS_T_code(::IsCSS, S::AbstractStabilizerCode; verbose::Bool = false)
    C_X = LinearCode(X_stabilizers(S)) # C_2
    C_Z = LinearCode(Z_stabilizers(S)) # dual(C_1)
    
    num_thrds = Threads.nthreads()
    verbose && println("Detected $num_thrds threads.")
    # should probably floor some logs here
    power = 0
    for i in 1:20
        if 2^i > num_thrds
            power = i - 1
            break
        end
    end

    p = Int(characteristic(C_X.F))
    n = S.n
    k = C_X.k
    z = zeros(Int, n)
    G = _Flint_matrix_to_Julia_int_matrix(transpose(C_X.G))
    for r in 1:k
        verbose && println("r: $r")
        flag = Threads.Atomic{Bool}(true)
        Threads.@threads for m in 1:num_thrds
            c = deepcopy(z)
            u = zeros(Int, k)
            for support in Combinatorics.combinations(1:k, r)
                if flag[]
                    fill!(u, 0)
                    @inbounds for i in support
                        u[i] = 1
                    end
                    LinearAlgebra.mul!(c, G, u)
                    wt_c = 0
                    supp_c = Int[]
                    @inbounds for j in 1:n
                        c[j] % p != 0 && (wt_c += 1; append!(supp_c, j);)
                    end

                    if !iseven(wt_c)
                        Threads.atomic_cas!(flag, true, false)
                        verbose && println("not even weight")
                        break
                    end
                    verbose && println("shortening on ", setdiff(1:n, supp_c))
                    
                    temp = shorten(C_Z, setdiff(1:n, supp_c))
                    if !contains_self_dual_subcode(temp)
                        Threads.atomic_cas!(flag, true, false)
                        verbose && println("m: $m, no subcode")
                        break
                    end
                else
                    break
                end
            end
            if !flag[]
                break
            end
        end
        if !flag[]
            return false
        end
    end
    return true
end
is_CSS_T_code(::IsNotCSS, S::AbstractStabilizerCode; verbose::Bool = false) = false

"""
$(TYPEDSIGNATURES)

Return `true` if the CSS code `S` is triorthogonal.
Cached to avoid recomputing on subsequent checks.
"""
function is_triorthogonal(S::AbstractStabilizerCodeCSS; verbose::Bool=false)
    haskey(S.cache, :is_triorthogonal) && return S.cache[:is_triorthogonal]
    
    result = is_triorthogonal(X_stabilizers(S), verbose)
    S.cache[:is_triorthogonal] = result
    return result
end
