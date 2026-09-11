# Copyright (c) 2023 - 2026 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    SubsystemCode(G::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing} = missing)

Return the subsystem code whose gauge group is determined by `G`.
"""
function SubsystemCode(G::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing} = missing)

    is_sparse = _is_sparse_code_matrix(G)
    G = _normalize_quantum_matrix(G)
    iszero(G) && throw(ArgumentError("The gauge matrix is empty."))
    G = _remove_empty(G, :rows)

    F = _code_matrix_base_ring(G)
    p = Int(characteristic(F))
    n = div(ncols(G), 2)
    clean_char_vec = ismissing(char_vec) ? zzModRingElem[] : _process_char_vec(char_vec, p, 2 * n)

    G_dense = _dense_code_matrix(G, F)

    # Stabilizers are the trace-symplectic center of the additive gauge group.
    stabs = _additive_center(G, F)
    stabs_dense = _dense_code_matrix(stabs, F)

    rnk_stabs = _additive_rank(stabs, F)
    rnk_gauge = _additive_rank(G, F)
    if rnk_stabs == rnk_gauge
        println("Stabilizer code detected.")
        return StabilizerCode(is_sparse ? _sparse_code_matrix(stabs_dense) : stabs_dense, char_vec = char_vec)
    end

    k, r = _subsystem_dimensions(F, n, rnk_stabs, rnk_gauge)
    is_css_S, X_stabs_dense, Z_stabs_dense = if rnk_stabs == 0
        (true, zero_matrix(F, 0, n), zero_matrix(F, 0, n))
    else
        robust_CSS_split(stabs_dense)
    end
    is_css_G, _, _ = robust_CSS_split(G_dense)

    # Bare logicals live in C(G)/S; dressed gauges live in G/S.
    centralizer = _additive_centralizer(G, F)
    bare_logs = _pair_operators(
        _additive_quotient_space(stabs, centralizer, F), is_css_S)
    gauge_ops = _pair_operators(
        _additive_quotient_space(stabs, G, F), is_css_S)
    length(gauge_ops) == r || error("Failed to extract $r symplectic gauge pairs.")
    expected_logical_pairs = _logical_pair_count(F, k)
    length(bare_logs) == expected_logical_pairs ||
        error("Failed to extract $expected_logical_pairs symplectic logical pairs.")

    if is_sparse
        stabs = _sparse_code_matrix(stabs_dense)
        gauge_ops = [(_sparse_code_matrix(pair[1]),
                      _sparse_code_matrix(pair[2])) for pair in gauge_ops]
        bare_logs = [(_sparse_code_matrix(pair[1]),
                      _sparse_code_matrix(pair[2])) for pair in bare_logs]
    end

    cache = Dict{Symbol, Any}(
        :stabs => stabs,
        :gauge_ops => gauge_ops,
        :g_ops_mat => _pairs_to_matrix(gauge_ops, F, n, is_sparse),
        :overcomplete => nrows(stabs_dense) > rnk_stabs,
        :logs_alg => :sys_eqs
    )
    if k > 0
        cache[:logicals] = bare_logs
        cache[:logs_mat] = _pairs_to_matrix(bare_logs, F, n, is_sparse)
    end

    if is_css_S && is_css_G
        X_stabs = is_sparse ? _sparse_code_matrix(X_stabs_dense) : X_stabs_dense
        Z_stabs = is_sparse ? _sparse_code_matrix(Z_stabs_dense) : Z_stabs_dense
        result = SubsystemCodeCSS(
            F, n, k, r, X_stabs, Z_stabs, gauge_ops, clean_char_vec, cache)
        return _seed_quantum_singleton_bound!(result)
    else
        result = SubsystemCode(
            F, n, k, r, stabs, gauge_ops, clean_char_vec, cache)
        return _seed_quantum_singleton_bound!(result)
    end
end

"""
    SubsystemCode(C::AbstractLinearCode, F; basis=missing, ...)

Construct the subsystem code over `F` whose gauge space is the additive
symplectic image of a linear code over the quadratic extension of `F`.
"""
function SubsystemCode(
    C::AbstractLinearCode, F::CTFieldTypes;
    basis::Union{Missing, Vector{<:CTFieldElem}}=missing,
    char_vec::Union{Vector{zzModRingElem}, Missing}=missing
)
    gauge_group = _quadratic_code_to_symplectic(C, F, basis)
    S = SubsystemCode(gauge_group; char_vec=char_vec)
    S.cache[:quadratic_code] = C
    S.cache[:quadratic_basis] =
        ismissing(basis) ? first(primitive_basis(C.F, F)) : basis
    return S
end

"""
    SubsystemCodeCSS(X_gauges, Z_gauges; char_vec=missing)
    CSSSubsystemCode(X_gauges, Z_gauges; char_vec=missing)

Construct a CSS subsystem code from trimmed `X`- and `Z`-type gauge
generators. The constructor may return a stabilizer code when the gauge group
is abelian.
"""
function SubsystemCodeCSS(
    X_gauges::T, Z_gauges::T;
    char_vec::Union{Vector{zzModRingElem}, Missing}=missing
) where {T <: CTMatrixTypes}
    is_sparse = _is_sparse_code_matrix(X_gauges)
    X_gauges = _normalize_quantum_matrix(X_gauges)
    Z_gauges = _normalize_quantum_matrix(Z_gauges)
    ncols(X_gauges) == ncols(Z_gauges) ||
        throw(ArgumentError("The X and Z gauge matrices must have the same length."))
    _code_matrix_base_ring(X_gauges) == _code_matrix_base_ring(Z_gauges) ||
        throw(ArgumentError("The X and Z gauge matrices must use the same field."))
    gauges = _css_symplectic_matrix(X_gauges, Z_gauges, false)
    is_sparse && (gauges = _sparse_code_matrix(gauges))
    return SubsystemCode(gauges; char_vec=char_vec)
end

function SubsystemCodeCSS(
    C_X::AbstractLinearCode, C_Z::AbstractLinearCode;
    char_vec::Union{Vector{zzModRingElem}, Missing}=missing
)
    C_X.F == C_Z.F ||
        throw(ArgumentError("The X and Z gauge codes must use the same field."))
    C_X.n == C_Z.n ||
        throw(ArgumentError("The X and Z gauge codes must have the same length."))
    return SubsystemCodeCSS(
        generator_matrix(C_X), generator_matrix(C_Z); char_vec=char_vec)
end

CSSSubsystemCode(args...; kwargs...) = SubsystemCodeCSS(args...; kwargs...)

"""
    random_subsystem_code([rng], F, n, k, r; char_vec=missing)

Construct a random (not guaranteed uniformly sampled) `[[n,k,r]]` subsystem
code over `F`.
"""
function random_subsystem_code(
    rng::AbstractRNG, F::CTFieldTypes, n::Int,
    k::Union{Int, Rational}, r::Int;
    char_vec::Union{Vector{zzModRingElem}, Missing}=missing
)
    0 <= k <= n || throw(DomainError(k, "Expected 0 ≤ k ≤ n."))
    num_stabs_rat = degree(F) * (n - k) - r
    denominator(Rational{BigInt}(num_stabs_rat)) == 1 ||
        throw(DomainError(k, "The requested additive dimension is incompatible with the field."))
    num_stabs = Int(num_stabs_rat)
    0 <= r && num_stabs >= 0 ||
        throw(DomainError(r, "The gauge dimension is incompatible with n and k."))
    pairs = _random_symplectic_pairs(rng, F, n)
    rows = CTMatrixTypes[]
    append!(rows, [pairs[i][1] for i in 1:num_stabs])
    for i in (num_stabs + 1):(num_stabs + r)
        push!(rows, pairs[i][1], pairs[i][2])
    end
    gauge_group = isempty(rows) ? zero_matrix(F, 0, 2n) : reduce(vcat, rows)
    isempty(rows) &&
        return StabilizerCode(gauge_group; char_vec=char_vec)
    return SubsystemCode(gauge_group; char_vec=char_vec)
end
random_subsystem_code(
    F::CTFieldTypes, n::Int, k::Union{Int, Rational}, r::Int; kwargs...
) =
    random_subsystem_code(Random.default_rng(), F, n, k, r; kwargs...)
random_subsystem_code(
    rng::AbstractRNG, n::Int, k::Union{Int, Rational}, r::Int; kwargs...
) = random_subsystem_code(
    rng, Oscar.Nemo.Native.GF(2), n, k, r; kwargs...)
random_subsystem_code(n::Int, k::Union{Int, Rational}, r::Int; kwargs...) =
    random_subsystem_code(Random.default_rng(), n, k, r; kwargs...)

"""
    SubsystemCode(G_Pauli::Vector{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing) where T <: Union{String, Vector{Char}}

Return the subsystem code whose gauge group is determined by the vector of Pauli strings `G_Pauli`.
"""
function SubsystemCode(G_Pauli::Vector{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing) where T <: Union{String, Vector{Char}}

    G_Pauli_stripped = _process_strings(G_Pauli)
    G = _Pauli_string_to_symplectic(G_Pauli_stripped)
    iszero(G) && error("The processed Pauli strings returned a set of empty gauge group generators.")
    return SubsystemCode(G, char_vec = char_vec)
end

"""
    SubsystemCode(S::CTMatrixTypes, L::CTMatrixTypes, G::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing} = missing)

Return the subsystem code whose stabilizers are given by `S`, (bare) logical operators
by `L`, gauge operators (not including stabilizers) by `G`.
"""
function SubsystemCode(S::CTMatrixTypes, L::CTMatrixTypes, G::CTMatrixTypes;
    char_vec::Union{Vector{zzModRingElem}, Missing} = missing)

    is_sparse = _is_sparse_code_matrix(S)
    S = _normalize_quantum_matrix(S)
    L = _normalize_quantum_matrix(L)
    G = _normalize_quantum_matrix(G)
    iszero(S) && error("The stabilizer matrix is empty.")
    S = _remove_empty(S, :rows)
    n = div(ncols(S), 2)
    F = _code_matrix_base_ring(S)
    p = Int(characteristic(F))

    S_dense = _dense_code_matrix(S, F)
    L_dense = _dense_code_matrix(L, F)
    G_dense = _dense_code_matrix(G, F)

    are_symplectic_orthogonal(S_dense, S_dense) || error("The given stabilizers are not symplectic orthogonal.")

    is_css_S, X_stabs_dense, Z_stabs_dense = robust_CSS_split(S_dense)

    # logicals validation
    if iszero(L_dense) || nrows(L_dense) == 0
        log_pairs_dense = Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}()
        logs_mat_dense = zero_matrix(F, 0, 2n)
    else
        L_dense = _remove_empty(L_dense, :rows)
        are_symplectic_orthogonal(S_dense, L_dense) || error("Logicals do not commute with the code.")
        prod = _trace_symplectic_product_matrix(L_dense, L_dense, F)
        iszero(prod) && error("Logicals should not be symplectic self-orthogonal.")
        
        is_css_L, L_X, L_Z = robust_CSS_split(L_dense)
        if is_css_S && is_css_L
            log_pairs_dense = _make_CSS_pairs(L_X, L_Z)
        else
            log_pairs_dense = _make_pairs(L_dense)
        end
        logs_mat_dense = reduce(vcat, [reduce(vcat, log_pairs_dense[i]) for i in 1:length(log_pairs_dense)])
    end

    # gauge validation
    if iszero(G_dense) || nrows(G_dense) == 0
        g_ops_pairs_dense = Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}()
        g_ops_mat_dense = zero_matrix(F, 0, 2n)
    else
        G_dense = _remove_empty(G_dense, :rows)
        are_symplectic_orthogonal(S_dense, G_dense) || error("Gauges do not commute with the code.")
        if !iszero(logs_mat_dense)
            are_symplectic_orthogonal(logs_mat_dense, G_dense) || error("Gauges do not commute with the logicals.")
        end
        prod = _trace_symplectic_product_matrix(G_dense, G_dense, F)
        iszero(prod) && error("Gauges should not be symplectic self-orthogonal.")
        
        is_css_G, G_X, G_Z = robust_CSS_split(G_dense)
        if is_css_S && is_css_G
            g_ops_pairs_dense = _make_CSS_pairs(G_X, G_Z)
        else
            g_ops_pairs_dense = _make_pairs(G_dense)
        end
        g_ops_mat_dense = reduce(vcat, [reduce(vcat, g_ops_pairs_dense[i]) for i in 1:length(g_ops_pairs_dense)])
    end

    clean_char_vec = ismissing(char_vec) ? zzModRingElem[] : _process_char_vec(char_vec, p, 2 * n)
    rnk_S = _additive_rank(S_dense, F)
    rnk_G_only = (iszero(G_dense) || nrows(G_dense) == 0) ?
        0 : _additive_rank(G_dense, F)
    k, r = _subsystem_dimensions(F, n, rnk_S, rnk_S + rnk_G_only)
    length(g_ops_pairs_dense) == r ||
        error("Expected $r independent symplectic gauge pairs.")
    expected_logical_pairs = _logical_pair_count(F, k)
    length(log_pairs_dense) == expected_logical_pairs ||
        error("Expected $expected_logical_pairs independent symplectic logical pairs.")

    S_final = is_sparse ? _sparse_code_matrix(S_dense) : S_dense
    g_ops_mat_final = is_sparse ? _sparse_code_matrix(g_ops_mat_dense) : g_ops_mat_dense
    logs_mat_final = is_sparse ? _sparse_code_matrix(logs_mat_dense) : logs_mat_dense

    log_pairs = Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}()
    for i in 1:length(log_pairs_dense)
        push!(log_pairs, (is_sparse ? _sparse_code_matrix(log_pairs_dense[i][1]) : log_pairs_dense[i][1],
                          is_sparse ? _sparse_code_matrix(log_pairs_dense[i][2]) : log_pairs_dense[i][2]))
    end

    g_ops_pairs = Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}()
    for i in 1:length(g_ops_pairs_dense)
        push!(g_ops_pairs, (is_sparse ? _sparse_code_matrix(g_ops_pairs_dense[i][1]) : g_ops_pairs_dense[i][1],
                            is_sparse ? _sparse_code_matrix(g_ops_pairs_dense[i][2]) : g_ops_pairs_dense[i][2]))
    end

    cache = Dict{Symbol, Any}(
        :stabs => S_final,
        :logicals => log_pairs,
        :logs_mat => logs_mat_final,
        :gauge_ops => g_ops_pairs,
        :g_ops_mat => g_ops_mat_final,
        :overcomplete => nrows(S_dense) > rnk_S,
        :logs_alg => :provided
    )

    if is_css_S
        X_stabs = is_sparse ? _sparse_code_matrix(X_stabs_dense) : X_stabs_dense
        Z_stabs = is_sparse ? _sparse_code_matrix(Z_stabs_dense) : Z_stabs_dense
        return SubsystemCodeCSS(F, n, k, r, X_stabs, Z_stabs, g_ops_pairs, clean_char_vec, cache)
    else
        return SubsystemCode(F, n, k, r, S_final, g_ops_pairs, clean_char_vec, cache)
    end
end

"""
    SubsystemCode(S_Pauli::Vector{T}, L_Pauli::Vector{T}, G_Pauli::Vector{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing) where T <: Union{String, Vector{Char}}
"""
function SubsystemCode(S_Pauli::Vector{T}, L_Pauli::Vector{T}, G_Pauli::Vector{T};
    char_vec::Union{Vector{zzModRingElem}, Missing} = missing) where T <: Union{String, Vector{Char}}

    S = _Pauli_string_to_symplectic(_process_strings(S_Pauli))
    iszero(S) && error("The processed Pauli strings returned a set of empty stabilizer generators.")
    L = _Pauli_string_to_symplectic(_process_strings(L_Pauli))
    iszero(L) && error("The processed Pauli strings returned a set of empty logical generators.")
    G = _Pauli_string_to_symplectic(_process_strings(G_Pauli))
    iszero(G) && error("The processed Pauli strings returned a set of empty gauge group generators.")
    return SubsystemCode(S, L, G, char_vec = char_vec)
end

#############################
      # getter functions
#############################

"""
    field(S::AbstractSubsystemCode)

Return the base ring of the code.
"""
field(S::AbstractSubsystemCode) = S.F

"""
    length(S::AbstractSubsystemCode)
    num_qubits(S::AbstractSubsystemCode)

Return the length of the code.
"""
length(S::AbstractSubsystemCode) = S.n
num_qubits(S::AbstractSubsystemCode) = S.n

"""
    dimension(S::AbstractSubsystemCode)

Return the dimension of the code.
"""
dimension(S::AbstractSubsystemCode) = S.k

"""
    cardinality(S::AbstractSubsystemCode)

Return the cardinality of the stabilizer group of the code.
"""
cardinality(S::AbstractSubsystemCode) =
    BigInt(characteristic(S.F))^_additive_rank(stabilizers(S), S.F)

"""
    rate(S::AbstractSubsystemCode)

Return the rate, `R = k/n`, of the code.
"""
rate(S::AbstractSubsystemCode) = S.k / S.n

"""
    is_overcomplete(S::AbstractSubsystemCode)

Return `true` if `S` has an overcomplete set of stabilizers.
"""
is_overcomplete(S::AbstractSubsystemCode) = get(S.cache, :overcomplete, false)

"""
    is_CSS(S::AbstractSubsystemCode)

Return `true` if `S` is CSS.
"""
is_CSS(S::T) where {T <: AbstractSubsystemCode} = is_CSS(CSSTrait(T), S)
is_CSS(::IsCSS, S::AbstractSubsystemCode) = true
is_CSS(::IsNotCSS, S::AbstractSubsystemCode) = false

"""
    signs(S::AbstractSubsystemCode)

Return the signs of the stabilizers of the code. Lazily computes them if a non-trivial character vector exists.
"""
function signs(S::AbstractSubsystemCode)
    haskey(S.cache, :signs) && return S.cache[:signs]
    
    stabs = stabilizers(S)
    nr = nrows(stabs)
    
    if isempty(S.char_vec)
        R, _ = residue_ring(Nemo.ZZ, Int(characteristic(S.F)) == 2 ? 4 : Int(characteristic(S.F)))
        S.cache[:signs] = [R(0) for _ in 1:nr]
    else
        S.cache[:signs] = _get_signs(stabs, S.char_vec)
    end
    return S.cache[:signs]
end

"""
    X_signs(S::AbstractSubsystemCode)

Return the signs of the `X` stabilizers of the CSS code.
"""
X_signs(S::T) where {T <: AbstractSubsystemCode} = X_signs(CSSTrait(T), S)
X_signs(::IsCSS, S::AbstractSubsystemCode) = get!(S.cache, :X_signs) do
    signs(S)[1:nrows(S.X_stabs)]
end
X_signs(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    Z_signs(S::AbstractSubsystemCode)

Return the signs of the `Z` stabilizers of the CSS code.
"""
Z_signs(S::T) where {T <: AbstractSubsystemCode} = Z_signs(CSSTrait(T), S)
Z_signs(::IsCSS, S::AbstractSubsystemCode) = get!(S.cache, :Z_signs) do
    signs(S)[nrows(S.X_stabs)+1:end]
end
Z_signs(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    stabilizers(S::AbstractSubsystemCode; standform::Bool = false)

Return the stabilizer matrix of the code. Computes the unified matrix for CSS codes if missing.

# Notes
- If the optional parameter `standform` is set to `true`, the standard form of the
  stabilizer matrix is returned instead (computed densely, stored sparsely if applicable).
"""
function stabilizers(S::AbstractSubsystemCode; standform::Bool = false)
    if standform
        degree(S.F) == 1 ||
            error("Stabilizer standard form currently requires a prime-field symplectic representation.")
        if !haskey(S.cache, :stabs_stand)
            # Standard form algorithms require dense matrices
            stabs_base = stabilizers(S)
            is_sparse = _is_sparse_code_matrix(stabs_base)
            stabs_dense = _dense_code_matrix(stabs_base, S.F)
            
            stabs_stand, P_stand, stand_r, stand_k, _ = _standard_form_stabilizer(stabs_dense)
            
            S.cache[:stabs_stand] =
                is_sparse ? _sparse_code_matrix(stabs_stand) : stabs_stand
            S.cache[:P_stand] = ismissing(P_stand) ? missing :
                (is_sparse ? _sparse_code_matrix(P_stand) : P_stand)
            S.cache[:stand_r] = stand_r
            S.cache[:stand_k] = stand_k
        end
        return S.cache[:stabs_stand]
    end
    
    # Base stabs fetch logic
    haskey(S.cache, :stabs) && return S.cache[:stabs]
    hasproperty(S, :stabs) && return S.stabs
    
    # Lazily build unified stabs for CSS codes
    if CSSTrait(typeof(S)) == IsCSS()
        X_stabs = _dense_code_matrix(S.X_stabs, S.F)
        Z_stabs = _dense_code_matrix(S.Z_stabs, S.F)
        unified_stabs = direct_sum(X_stabs, Z_stabs)
        S.cache[:stabs] = unified_stabs
        return unified_stabs
    end
    error("Stabilizers not found in cache or struct.")
end

"""
    standard_form_permutation(S::AbstractSubsystemCode)

Return the permutation matrix required to permute the columns of the code matrices to have the same
row space as the matrices in standard form. Returns `missing` if no such permutation is required.
"""
function standard_form_permutation(S::AbstractSubsystemCode)
    haskey(S.cache, :P_stand) && return S.cache[:P_stand]
    stabilizers(S, standform=true) # Triggers the standard form computation
    return get(S.cache, :P_stand, missing)
end

"""
    X_stabilizers(S::AbstractSubsystemCode)

Return the `X`-stabilizer matrix of the CSS code.
"""
X_stabilizers(S::T) where {T <: AbstractSubsystemCode} = X_stabilizers(CSSTrait(T), S)
X_stabilizers(::IsCSS, S::AbstractSubsystemCode) = S.X_stabs
X_stabilizers(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    Z_stabilizers(S::AbstractSubsystemCode)

Return the `Z`-stabilizer matrix of the CSS code.
"""
Z_stabilizers(S::T) where {T <: AbstractSubsystemCode} = Z_stabilizers(CSSTrait(T), S)
Z_stabilizers(::IsCSS, S::AbstractSubsystemCode) = S.Z_stabs
Z_stabilizers(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    standard_form_A(S::AbstractSubsystemCode)

Return the named matrix `A` from the standard form of the stabilizer matrix.
"""
function standard_form_A(S::AbstractSubsystemCode)
    stabs_stand = stabilizers(S; standform = true)
    r = S.cache[:stand_r]
    return stabs_stand[1:r, r + 1:S.n]
end

"""
    standard_form_A1(S::AbstractSubsystemCode)
    
Return the named matrix `A1` from the standard form of the stabilizer matrix.
"""
function standard_form_A1(S::AbstractSubsystemCode)
    stabs_stand = stabilizers(S; standform = true)
    r = S.cache[:stand_r]
    k = S.cache[:stand_k]
    return stabs_stand[1:r, r + 1:S.n - k]
end

"""
    standard_form_A2(S::AbstractSubsystemCode)
    
Return the named matrix `A2` from the standard form of the stabilizer matrix.
"""
function standard_form_A2(S::AbstractSubsystemCode)
    stabs_stand = stabilizers(S; standform = true)
    r = S.cache[:stand_r]
    k = S.cache[:stand_k]
    return stabs_stand[1:r, S.n - k + 1:S.n]
end

"""
    standard_form_B(S::AbstractSubsystemCode)
    
Return the named matrix `B` from the standard form of the stabilizer matrix.
"""
function standard_form_B(S::AbstractSubsystemCode)
    stabs_stand = stabilizers(S; standform = true)
    r = S.cache[:stand_r]
    return stabs_stand[1:r, S.n + 1:S.n + r]
end

"""
    standard_form_C1(S::AbstractSubsystemCode)
    
Return the named matrix `C1` from the standard form of the stabilizer matrix.
"""
function standard_form_C1(S::AbstractSubsystemCode)
    stabs_stand = stabilizers(S; standform = true)
    r = S.cache[:stand_r]
    k = S.cache[:stand_k]
    return stabs_stand[1:r, S.n + r + 1:2 * S.n - k]
end

"""
    standard_form_C2(S::AbstractSubsystemCode)
    
Return the named matrix `C2` from the standard form of the stabilizer matrix.
"""
function standard_form_C2(S::AbstractSubsystemCode)
    stabs_stand = stabilizers(S; standform = true)
    r = S.cache[:stand_r]
    k = S.cache[:stand_k]
    return stabs_stand[1:r, 2 * S.n - k + 1:2 * S.n]
end

"""
    standard_form_D(S::AbstractSubsystemCode)
    
Return the named matrix `D` from the standard form of the stabilizer matrix.
"""
function standard_form_D(S::AbstractSubsystemCode)
    stabs_stand = stabilizers(S; standform = true)
    r = S.cache[:stand_r]
    k = S.cache[:stand_k]
    return stabs_stand[r + 1:S.n - k, S.n + 1:S.n + r]
end

"""
    standard_form_E(S::AbstractSubsystemCode)
    
Return the named matrix `E` from the standard form of the stabilizer matrix.
"""
function standard_form_E(S::AbstractSubsystemCode)
    stabs_stand = stabilizers(S; standform = true)
    r = S.cache[:stand_r]
    k = S.cache[:stand_k]
    return stabs_stand[r + 1:S.n - k, 2 * S.n - k + 1:2 * S.n]
end

"""
    metacheck(S::AbstractSubsystemCode)

Return the metacheck matrix of the code, if it has been set; otherwise returns missing.
"""
metacheck(S::T) where {T <: AbstractSubsystemCode} = metacheck(CSSTrait(T), S)
metacheck(::IsCSS, S::AbstractSubsystemCode) = error("Use `X_metacheck` or `Z_metacheck` for CSS codes.")
metacheck(::IsNotCSS, S::AbstractSubsystemCode) = get(S.cache, :metacheck, missing)

"""
    X_metacheck(S::AbstractSubsystemCode)

Return the `X`-metacheck matrix of the CSS code, if it has been set; otherwise returns missing.
"""
X_metacheck(S::T) where {T <: AbstractSubsystemCode} = X_metacheck(CSSTrait(T), S)
X_metacheck(::IsCSS, S::AbstractSubsystemCode) = get(S.cache, :X_metacheck, missing)
X_metacheck(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    Z_metacheck(S::AbstractSubsystemCode)

Return the `Z`-metacheck matrix of the CSS code, if it has been set; otherwise returns missing.
"""
Z_metacheck(S::T) where {T <: AbstractSubsystemCode} = Z_metacheck(CSSTrait(T), S)
Z_metacheck(::IsCSS, S::AbstractSubsystemCode) = get(S.cache, :Z_metacheck, missing)
Z_metacheck(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    logicals(S::AbstractSubsystemCode)
    logical_operators(S::AbstractSubsystemCode)
    bare_logicals(S::AbstractSubsystemCode)
    bare(S::AbstractSubsystemCode)

Return a vector of logical operator generator pairs for `S`.
"""
logicals(S::T) where {T <: AbstractSubsystemCode} = logicals(LogicalTrait(T), S)
function logicals(::HasLogicals, S::AbstractSubsystemCode)
    S.k == 0 && return Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}()
    haskey(S.cache, :logicals) && return S.cache[:logicals]

    GaugeTrait(typeof(S)) == HasNoGauges() ||
        error("Logicals not found in cache. For subsystem codes, logicals should be seeded during initialization.")

    logs_alg = get(S.cache, :logs_alg, :stnd_frm)
    stabs = stabilizers(S)
    is_sparse = _is_sparse_code_matrix(stabs)
    stabs_dense = _dense_code_matrix(stabs, S.F)

    if logs_alg == :sys_eqs
        dual_gens = _additive_centralizer(stabs, S.F)
        logs_dense, logs_mat_dense = _logicals(stabs, dual_gens, :sys_eqs)
        logs = Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}()
        for i in eachindex(logs_dense)
            push!(logs, (is_sparse ? _sparse_code_matrix(logs_dense[i][1]) : logs_dense[i][1],
                         is_sparse ? _sparse_code_matrix(logs_dense[i][2]) : logs_dense[i][2]))
        end
        S.cache[:logicals] = logs
        S.cache[:logs_mat] =
            is_sparse ? _sparse_code_matrix(logs_mat_dense) : logs_mat_dense
        return logs
    end

    stabs_stand = stabilizers(S; standform = true)
    r = S.cache[:stand_r]
    k = S.cache[:stand_k]
    P = S.cache[:P_stand]
    stabs_stand_dense = _dense_code_matrix(stabs_stand, S.F)
    P_dense = (P === missing) ? missing : _dense_code_matrix(P, S.F)
    logs_dense = _make_pairs(_logicals_standard_form(stabs_stand_dense, S.n, k, r, P_dense))
    logs = Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}()
    for i in eachindex(logs_dense)
        push!(logs, (is_sparse ? _sparse_code_matrix(logs_dense[i][1]) : logs_dense[i][1],
                     is_sparse ? _sparse_code_matrix(logs_dense[i][2]) : logs_dense[i][2]))
    end
    S.cache[:logicals] = logs
    S.cache[:logs_mat] = reduce(vcat, [reduce(vcat, logs[i]) for i in eachindex(logs)])
    return logs
end
logicals(::HasNoLogicals, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no logicals.")
logical_operators(S::AbstractSubsystemCode) = logicals(S)
bare_logicals(S::AbstractSubsystemCode) = logicals(S)
bare(S::AbstractSubsystemCode) = logicals(S)

"""
    logicals_matrix(S::AbstractSubsystemCode)

Returns the result of `logicals(S)` as a vertically concatenated matrix.
"""
logicals_matrix(S::T) where {T <: AbstractSubsystemCode} = logicals_matrix(LogicalTrait(T), S)
function logicals_matrix(::HasLogicals, S::AbstractSubsystemCode)
    S.k == 0 && return zero_matrix(S.F, 0, 2 * S.n)
    haskey(S.cache, :logs_mat) && return S.cache[:logs_mat]
    logicals(S) # Trigger computation
    return S.cache[:logs_mat]
end
logicals_matrix(::HasNoLogicals, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no logicals.")

"""
    logicals_standard_form(S::AbstractSubsystemCode)

Return a matrix of logical operators as determined by the stabilizers in standard form.
"""
logicals_standard_form(S::T) where {T <: AbstractSubsystemCode} = logicals_standard_form(LogicalTrait(T), S)
function logicals_standard_form(::HasLogicals, S::AbstractSubsystemCode)
    degree(S.F) == 1 ||
        error("Standard-form logical extraction requires a prime-field symplectic representation; use logicals(S) for additive extension-field codes.")
    stabs_stand = stabilizers(S; standform = true)
    
    stabs_dense = _dense_code_matrix(stabs_stand, S.F)
    P = S.cache[:P_stand]
    P_dense = P === missing ? missing : _dense_code_matrix(P, S.F)
    
    logs_stand = _logicals_standard_form(stabs_dense, S.n, S.cache[:stand_k], S.cache[:stand_r], P_dense)
    return _is_sparse_code_matrix(stabs_stand) ?
        _sparse_code_matrix(logs_stand) : logs_stand
end
logicals_standard_form(::HasNoLogicals, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no logicals.")

"""
    gauges(S::AbstractSubsystemCode)
    gauge_operators(S::AbstractSubsystemCode)

Return a vector of gauge operator generator pairs for `S`.
"""
gauges(S::T) where {T <: AbstractSubsystemCode} = gauges(GaugeTrait(T), S)
function gauges(::HasGauges, S::AbstractSubsystemCode)
    # Automatically seeded in constructor for subsystem codes
    haskey(S.cache, :gauge_ops) && return S.cache[:gauge_ops]
    hasproperty(S, :gauge_ops) && return S.gauge_ops
    error("Gauge operators not found in cache.")
end
gauges(::HasNoGauges, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no gauges.")
gauge_operators(S::AbstractSubsystemCode) = gauges(S)

"""
    gauges_matrix(S::AbstractSubsystemCode)
    gauge_operators_matrix(S::AbstractSubsystemCode)

Return the result of `gauges(S)` as a vertically concatenated matrix.
"""
gauges_matrix(S::T) where {T <: AbstractSubsystemCode} = gauges_matrix(GaugeTrait(T), S)
function gauges_matrix(::HasGauges, S::AbstractSubsystemCode)
    haskey(S.cache, :g_ops_mat) && return S.cache[:g_ops_mat]
    hasproperty(S, :g_ops_mat) && return S.g_ops_mat
    error("Gauge matrix not found in cache.")
end
gauges_matrix(::HasNoGauges, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no gauges.")
gauge_operators_matrix(S::AbstractSubsystemCode) = gauges_matrix(S)

"""
    gauge_group(S::AbstractSubsystemCode)
    gauge_group_matrix(S::AbstractSubsystemCode)
    gauge_generators_matrix(S::AbstractSubsystemCode)
    gauge_group_generators_matrix(S::AbstractSubsystemCode)

Return a matrix giving a (maybe overcomplete) basis for the gauge group.
"""
gauge_group(S::T) where {T <: AbstractSubsystemCode} = gauge_group(GaugeTrait(T), S)
gauge_group(::HasGauges, S::AbstractSubsystemCode) =
    _vcat_code_matrices(S.F, stabilizers(S), gauges_matrix(S))
gauge_group(::HasNoGauges, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no gauges.")
gauge_group_matrix(S::AbstractSubsystemCode) = gauge_group(S)
gauge_generators_matrix(S::AbstractSubsystemCode) = gauge_group(S)
gauge_group_generators_matrix(S::AbstractSubsystemCode) = gauge_group(S)

# # -----------------------------------------------------------------------------
# # Distance Catchers
# # -----------------------------------------------------------------------------

# bare_minimum_distance_lower_bound(S::T) where T <: AbstractSubsystemCode = bare_minimum_distance_lower_bound(GaugeTrait(T), S)
# bare_minimum_distance_lower_bound(::HasGauges, S::AbstractSubsystemCode) = get(S.cache, :l_bound_bare, missing)
# bare_minimum_distance_lower_bound(::HasNoGauges, S::AbstractSubsystemCode) = error("Only valid for subsystem codes; use `minimum_distance_lower_bound` for stabilizer codes.")

# bare_minimum_distance_upper_bound(S::T) where T <: AbstractSubsystemCode = bare_minimum_distance_upper_bound(GaugeTrait(T), S)
# bare_minimum_distance_upper_bound(::HasGauges, S::AbstractSubsystemCode) = get(S.cache, :u_bound_bare, missing)
# bare_minimum_distance_upper_bound(::HasNoGauges, S::AbstractSubsystemCode) = error("Only valid for subsystem codes; use `minimum_distance_lower_bound` for stabilizer codes.")

# dressed_minimum_distance_lower_bound(S::T) where T <: AbstractSubsystemCode = dressed_minimum_distance_lower_bound(GaugeTrait(T), S)
# dressed_minimum_distance_lower_bound(::HasGauges, S::AbstractSubsystemCode) = get(S.cache, :l_bound_dressed, missing)
# dressed_minimum_distance_lower_bound(::HasNoGauges, S::AbstractSubsystemCode) = error("Only valid for subsystem codes; use `minimum_distance_lower_bound` for stabilizer codes.")

# dressed_minimum_distance_upper_bound(S::T) where T <: AbstractSubsystemCode = dressed_minimum_distance_upper_bound(GaugeTrait(T), S)
# dressed_minimum_distance_upper_bound(::HasGauges, S::AbstractSubsystemCode) = get(S.cache, :u_bound_dressed, missing)
# dressed_minimum_distance_upper_bound(::HasNoGauges, S::AbstractSubsystemCode) = error("Only valid for subsystem codes; use `minimum_distance_lower_bound` for stabilizer codes.")

# bare_X_minimum_distance_lower_bound(S::T) where T <: AbstractSubsystemCode = bare_X_minimum_distance_lower_bound(GaugeTrait(T), CSSTrait(T), S)
# bare_X_minimum_distance_lower_bound(::HasGauges, ::IsCSS, S::AbstractSubsystemCode) = get(S.cache, :l_bound_dx_bare, missing)
# bare_X_minimum_distance_lower_bound(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes")
# bare_X_minimum_distance_lower_bound(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for subsystem codes")
# bare_X_minimum_distance_lower_bound(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode) = error("Only valid for subsystem codes; use `X_minimum_distance_lower_bound` for stabilizer codes.")

# bare_X_minimum_distance_upper_bound(S::T) where T <: AbstractSubsystemCode = bare_X_minimum_distance_upper_bound(GaugeTrait(T), CSSTrait(T), S)
# bare_X_minimum_distance_upper_bound(::HasGauges, ::IsCSS, S::AbstractSubsystemCode) = get(S.cache, :u_bound_dx_bare, missing)
# bare_X_minimum_distance_upper_bound(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes")
# bare_X_minimum_distance_upper_bound(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for subsystem codes")
# bare_X_minimum_distance_upper_bound(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode) = error("Only valid for subsystem codes; use `X_minimum_distance_upper_bound` for stabilizer codes.")

# dressed_X_minimum_distance_lower_bound(S::T) where T <: AbstractSubsystemCode = dressed_X_minimum_distance_lower_bound(GaugeTrait(T), CSSTrait(T), S)
# dressed_X_minimum_distance_lower_bound(::HasGauges, ::IsCSS, S::AbstractSubsystemCode) = get(S.cache, :l_bound_dx_dressed, missing)
# dressed_X_minimum_distance_lower_bound(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes")
# dressed_X_minimum_distance_lower_bound(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for subsystem codes")
# dressed_X_minimum_distance_lower_bound(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode) = error("Only valid for subsystem codes; use `X_minimum_distance_lower_bound` for stabilizer codes.")

# dressed_X_minimum_distance_upper_bound(S::T) where T <: AbstractSubsystemCode = dressed_X_minimum_distance_upper_bound(GaugeTrait(T), CSSTrait(T), S)
# dressed_X_minimum_distance_upper_bound(::HasGauges, ::IsCSS, S::AbstractSubsystemCode) = get(S.cache, :u_bound_dx_dressed, missing)
# dressed_X_minimum_distance_upper_bound(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes")
# dressed_X_minimum_distance_upper_bound(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for subsystem codes")
# dressed_X_minimum_distance_upper_bound(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode) = error("Only valid for subsystem codes; use `X_minimum_distance_upper_bound` for stabilizer codes.")

# bare_Z_minimum_distance_lower_bound(S::T) where T <: AbstractSubsystemCode = bare_Z_minimum_distance_lower_bound(GaugeTrait(T), CSSTrait(T), S)
# bare_Z_minimum_distance_lower_bound(::HasGauges, ::IsCSS, S::AbstractSubsystemCode) = get(S.cache, :l_bound_dz_bare, missing)
# bare_Z_minimum_distance_lower_bound(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes")
# bare_Z_minimum_distance_lower_bound(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for subsystem codes")
# bare_Z_minimum_distance_lower_bound(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode) = error("Only valid for subsystem codes; use `Z_minimum_distance_lower_bound` for stabilizer codes.")

# bare_Z_minimum_distance_upper_bound(S::T) where T <: AbstractSubsystemCode = bare_Z_minimum_distance_upper_bound(GaugeTrait(T), CSSTrait(T), S)
# bare_Z_minimum_distance_upper_bound(::HasGauges, ::IsCSS, S::AbstractSubsystemCode) = get(S.cache, :u_bound_dz_bare, missing)
# bare_Z_minimum_distance_upper_bound(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes")
# bare_Z_minimum_distance_upper_bound(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for subsystem codes")
# bare_Z_minimum_distance_upper_bound(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode) = error("Only valid for subsystem codes; use `Z_minimum_distance_upper_bound` for stabilizer codes.")

# dressed_Z_minimum_distance_lower_bound(S::T) where T <: AbstractSubsystemCode = dressed_Z_minimum_distance_lower_bound(GaugeTrait(T), CSSTrait(T), S)
# dressed_Z_minimum_distance_lower_bound(::HasGauges, ::IsCSS, S::AbstractSubsystemCode) = get(S.cache, :l_bound_dz_dressed, missing)
# dressed_Z_minimum_distance_lower_bound(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes")
# dressed_Z_minimum_distance_lower_bound(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for subsystem codes")
# dressed_Z_minimum_distance_lower_bound(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode) = error("Only valid for subsystem codes; use `Z_minimum_distance_lower_bound` for stabilizer codes.")

# dressed_Z_minimum_distance_upper_bound(S::T) where T <: AbstractSubsystemCode = dressed_Z_minimum_distance_upper_bound(GaugeTrait(T), CSSTrait(T), S)
# dressed_Z_minimum_distance_upper_bound(::HasGauges, ::IsCSS, S::AbstractSubsystemCode) = get(S.cache, :u_bound_dz_dressed, missing)
# dressed_Z_minimum_distance_upper_bound(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes")
# dressed_Z_minimum_distance_upper_bound(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for subsystem codes")
# dressed_Z_minimum_distance_upper_bound(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode) = error("Only valid for subsystem codes; use `Z_minimum_distance_upper_bound` for stabilizer codes.")

"""
    num_X_stabs(S::AbstractSubsystemCode)

Return the number of `X` stabilizers of the CSS code.
"""
num_X_stabs(S::T) where {T <: AbstractSubsystemCode} = num_X_stabs(CSSTrait(T), S)
num_X_stabs(::IsCSS, S::AbstractSubsystemCode) = nrows(X_stabilizers(S))
num_X_stabs(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    num_Z_stabs(S::AbstractSubsystemCode)

Return the number of `Z` stabilizers of the CSS code.
"""
num_Z_stabs(S::T) where {T <: AbstractSubsystemCode} = num_Z_stabs(CSSTrait(T), S)
num_Z_stabs(::IsCSS, S::AbstractSubsystemCode) = nrows(Z_stabilizers(S))
num_Z_stabs(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    character_vector(S::AbstractSubsystemCode)

Return the character vector of the code.
"""
function character_vector(S::AbstractSubsystemCode)
    hasfield(typeof(S), :char_vec) && return getfield(S, :char_vec)
    return get(S.cache, :char_vec, zzModRingElem[])
end

# TODO: quantum Singletonbound k <= n - 2d + 2
# MDS/optimal for subsystem codes: k + r <= n - 2d + 2

"""
    relative_distance(S::AbstractSubsystemCode)

Return the relative minimum distance, `δ = d / n` of the code if `d` is known,
otherwise errors.
"""
function relative_distance(S::AbstractSubsystemCode)
    d = if GaugeTrait(typeof(S)) == HasGauges()
        get(S.cache, :d_dressed, get(S.cache, :d, missing))
    else
        get(S.cache, :d, missing)
    end
    ismissing(d) && error("Missing minimum distance for this code.")
    return d / S.n
end

"""
    dressed(S::AbstractSubsystemCode)
    dressed_operators(S::AbstractSubsystemCode
    dressed_logicals(S::AbstractSubsystemCode)

Return a vector of pairs generators for the dressed operators of `S`.

# Notes
- Here, the dressed operators are the logicals and the gauge operators.
"""
function dressed(S::T) where {T <: AbstractSubsystemCode}
    if LogicalTrait(T) == HasNoLogicals()
        error("Type $T has no logicals.")
    elseif GaugeTrait(T) == HasNoGauges()
        error("Type $T has no gauges.")
    end
    return logicals(S) ∪ gauges(S)
end
dressed_operators(S::AbstractSubsystemCode) = dressed(S)
dressed_logicals(S::AbstractSubsystemCode) = dressed(S)

"""
    bare_minimum_distance_lower_bound(S::AbstractSubsystemCode)
Return the currently stored lower bound on the bare minimum distance.
"""
bare_minimum_distance_lower_bound(S::T) where T <: AbstractSubsystemCode =
    bare_minimum_distance_lower_bound(GaugeTrait(T), S)
bare_minimum_distance_lower_bound(::HasGauges, S::AbstractSubsystemCode) = get(S.cache, :l_bound_bare, missing)
bare_minimum_distance_lower_bound(::HasNoGauges, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes; use `minimum_distance_lower_bound` for stabilizer codes.")

"""
    bare_minimum_distance_upper_bound(S::AbstractSubsystemCode)
Return the currently stored upper bound on the bare minimum distance.
"""
bare_minimum_distance_upper_bound(S::T) where T <: AbstractSubsystemCode =
    bare_minimum_distance_upper_bound(GaugeTrait(T), S)
bare_minimum_distance_upper_bound(::HasGauges, S::AbstractSubsystemCode) = get(S.cache, :u_bound_bare, missing)
bare_minimum_distance_upper_bound(::HasNoGauges, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes; use `minimum_distance_lower_bound` for stabilizer codes.")

"""
    dressed_minimum_distance_lower_bound(S::AbstractSubsystemCode)
Return the currently stored lower bound on the dressed minimum distance.
"""
dressed_minimum_distance_lower_bound(S::T) where T <: AbstractSubsystemCode =
    dressed_minimum_distance_lower_bound(GaugeTrait(T), S)
dressed_minimum_distance_lower_bound(::HasGauges, S::AbstractSubsystemCode) = get(S.cache, :l_bound_dressed, missing)
dressed_minimum_distance_lower_bound(::HasNoGauges, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes; use `minimum_distance_lower_bound` for stabilizer codes.")

"""
    dressed_minimum_distance_upper_bound(S::AbstractSubsystemCode)
Return the currently stored upper bound on the dressed minimum distance.
"""
dressed_minimum_distance_upper_bound(S::T) where T <: AbstractSubsystemCode =
    dressed_minimum_distance_upper_bound(GaugeTrait(T), S)
dressed_minimum_distance_upper_bound(::HasGauges, S::AbstractSubsystemCode) = get(S.cache, :u_bound_dressed, missing)
dressed_minimum_distance_upper_bound(::HasNoGauges, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes; use `minimum_distance_lower_bound` for stabilizer codes.")

"""
    bare_X_minimum_distance_lower_bound(S::AbstractSubsystemCode)
Return the currently stored lower bound on the bare `X`-minimum distance.
"""
bare_X_minimum_distance_lower_bound(S::T) where T <: AbstractSubsystemCode = bare_X_minimum_distance_lower_bound(GaugeTrait(T), CSSTrait(T), S)
bare_X_minimum_distance_lower_bound(::HasGauges, ::IsCSS, S::AbstractSubsystemCode) =
    get(S.cache, :l_bound_dx_bare, missing)
bare_X_minimum_distance_lower_bound(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode) =
    error("Only valid for CSS codes")
bare_X_minimum_distance_lower_bound(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes")
bare_X_minimum_distance_lower_bound(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes; use `X_minimum_distance_lower_bound` for stabilizer codes.")

"""
    bare_X_minimum_distance_upper_bound(S::AbstractSubsystemCode)
Return the currently stored upper bound on the bare `X`-minimum distance.
"""
bare_X_minimum_distance_upper_bound(S::T) where T <: AbstractSubsystemCode = bare_X_minimum_distance_upper_bound(GaugeTrait(T), CSSTrait(T), S)
bare_X_minimum_distance_upper_bound(::HasGauges, ::IsCSS, S::AbstractSubsystemCode) =
    get(S.cache, :u_bound_dx_bare, missing)
bare_X_minimum_distance_upper_bound(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode) =
    error("Only valid for CSS codes")
bare_X_minimum_distance_upper_bound(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes")
bare_X_minimum_distance_upper_bound(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes; use `X_minimum_distance_upper_bound` for stabilizer codes.")

"""
    dressed_X_minimum_distance_lower_bound(S::AbstractSubsystemCode)
Return the currently stored lower bound on the dressed `X`-minimum distance.
"""
dressed_X_minimum_distance_lower_bound(S::T) where T <: AbstractSubsystemCode = dressed_X_minimum_distance_lower_bound(GaugeTrait(T), CSSTrait(T), S)
dressed_X_minimum_distance_lower_bound(::HasGauges, ::IsCSS, S::AbstractSubsystemCode) =
    get(S.cache, :l_bound_dx_dressed, missing)
dressed_X_minimum_distance_lower_bound(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode) =
    error("Only valid for CSS codes")
dressed_X_minimum_distance_lower_bound(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes")
dressed_X_minimum_distance_lower_bound(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes; use `X_minimum_distance_lower_bound` for stabilizer codes.")

"""
    dressed_X_minimum_distance_upper_bound(S::AbstractSubsystemCode)
Return the currently stored upper bound on the dressed `X`-minimum distance.
"""
dressed_X_minimum_distance_upper_bound(S::T) where T <: AbstractSubsystemCode = dressed_X_minimum_distance_upper_bound(GaugeTrait(T), CSSTrait(T), S)
dressed_X_minimum_distance_upper_bound(::HasGauges, ::IsCSS, S::AbstractSubsystemCode) =
    get(S.cache, :u_bound_dx_dressed, missing)
dressed_X_minimum_distance_upper_bound(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode) =
    error("Only valid for CSS codes")
dressed_X_minimum_distance_upper_bound(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes")
dressed_X_minimum_distance_upper_bound(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes; use `X_minimum_distance_upper_bound` for stabilizer codes.")

"""
    bare_Z_minimum_distance_lower_bound(S::AbstractSubsystemCode)
Return the currently stored lower bound on the bare `Z`-minimum distance.
"""
bare_Z_minimum_distance_lower_bound(S::T) where T <: AbstractSubsystemCode = bare_Z_minimum_distance_lower_bound(GaugeTrait(T), CSSTrait(T), S)
bare_Z_minimum_distance_lower_bound(::HasGauges, ::IsCSS, S::AbstractSubsystemCode) =
    get(S.cache, :l_bound_dz_bare, missing)
bare_Z_minimum_distance_lower_bound(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode) =
    error("Only valid for CSS codes")
bare_Z_minimum_distance_lower_bound(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes")
bare_Z_minimum_distance_lower_bound(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes; use `Z_minimum_distance_lower_bound` for stabilizer codes.")

"""
    bare_Z_minimum_distance_upper_bound(S::AbstractSubsystemCode)
Return the currently stored upper bound on the bare `Z`-minimum distance.
"""
bare_Z_minimum_distance_upper_bound(S::T) where T <: AbstractSubsystemCode = bare_Z_minimum_distance_upper_bound(GaugeTrait(T), CSSTrait(T), S)
bare_Z_minimum_distance_upper_bound(::HasGauges, ::IsCSS, S::AbstractSubsystemCode) =
    get(S.cache, :u_bound_dz_bare, missing)
bare_Z_minimum_distance_upper_bound(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode) =
    error("Only valid for CSS codes")
bare_Z_minimum_distance_upper_bound(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes")
bare_Z_minimum_distance_upper_bound(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes; use `Z_minimum_distance_upper_bound` for stabilizer codes.")

"""
    dressed_Z_minimum_distance_lower_bound(S::AbstractSubsystemCode)
Return the currently stored lower bound on the dressed `Z`-minimum distance.
"""
dressed_Z_minimum_distance_lower_bound(S::T) where T <: AbstractSubsystemCode = dressed_Z_minimum_distance_lower_bound(GaugeTrait(T), CSSTrait(T), S)
dressed_Z_minimum_distance_lower_bound(::HasGauges, ::IsCSS, S::AbstractSubsystemCode) =
    get(S.cache, :l_bound_dz_dressed, missing)
dressed_Z_minimum_distance_lower_bound(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode) =
    error("Only valid for CSS codes")
dressed_Z_minimum_distance_lower_bound(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes")
dressed_Z_minimum_distance_lower_bound(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes; use `Z_minimum_distance_lower_bound` for stabilizer codes.")

"""
    dressed_Z_minimum_distance_upper_bound(S::AbstractSubsystemCode)
Return the currently stored upper bound on the dressed `Z`-minimum distance.
"""
dressed_Z_minimum_distance_upper_bound(S::T) where T <: AbstractSubsystemCode = dressed_Z_minimum_distance_upper_bound(GaugeTrait(T), CSSTrait(T), S)
dressed_Z_minimum_distance_upper_bound(::HasGauges, ::IsCSS, S::AbstractSubsystemCode) =
    get(S.cache, :u_bound_dz_dressed, missing)
dressed_Z_minimum_distance_upper_bound(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode) =
    error("Only valid for CSS codes")
dressed_Z_minimum_distance_upper_bound(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes")
dressed_Z_minimum_distance_upper_bound(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes; use `Z_minimum_distance_upper_bound` for stabilizer codes.")

#############################
      # setter functions
#############################

"""
    set_signs(S::AbstractSubsystemCode, char_vec::Vector{zzModRingElem})
    set_signs!(S::AbstractSubsystemCode, char_vec::Vector{zzModRingElem})

Set the character vector of `S` to `char_vec` and update the signs.
"""
function set_signs!(S::AbstractSubsystemCode, char_vec::Vector{zzModRingElem})
    if !isempty(char_vec)
        length(char_vec) == 2 * S.n || throw(ArgumentError("Characteristic vector is of improper length for the code."))
        R = base_ring(char_vec[1])
        R_expect = residue_ring(Nemo.ZZ, Int(characteristic(S.F)) == 2 ? 4 : Int(characteristic(S.F)))[1]
        modulus(R) == modulus(R_expect) || throw(ArgumentError("Phases are not in the correct ring."))
    end

    S.char_vec = char_vec
    
    # Fully Lazy Cache Invalidation
    delete!(S.cache, :signs)
    delete!(S.cache, :X_signs)
    delete!(S.cache, :Z_signs)
    return nothing
end
set_signs(S::AbstractSubsystemCode, char_vec::Vector{zzModRingElem}) = (S_new = deepcopy(S); set_signs!(S_new, char_vec); return S_new)

"""
    set_stabilizers(S::AbstractSubsystemCode, stabs::CTMatrixTypes)
    set_stabilizers!(S::AbstractSubsystemCode, stabs::CTMatrixTypes)

Set the stabilizers of `S` to `stabs`.
"""
function set_stabilizers!(S::AbstractSubsystemCode, stabs::CTMatrixTypes)
    is_sparse = _is_sparse_code_matrix(stabs)
    stabs = _normalize_quantum_matrix(stabs)
    iszero(stabs) && throw(ArgumentError("The stabilizers cannot be zero."))
    order(S.F) == order(_code_matrix_base_ring(stabs)) || throw(ArgumentError("The stabilizers must be over the same field as the code."))

    stabs = _remove_empty(stabs, :rows)
    _code_matrix_base_ring(stabs) == S.F ||
        (stabs = change_base_ring(S.F, _dense_code_matrix(
            stabs, _code_matrix_base_ring(stabs))))
    is_sparse && (stabs = _sparse_code_matrix(stabs))
    if _additive_row_spaces_equal(stabilizers(S), stabs, S.F)
        if is_CSS(S)
            S.cache[:stabs] = stabs
        else
            S.stabs = stabs
        end
        
        expected_rank_rat = degree(S.F) * (S.n - S.k) -
            (GaugeTrait(typeof(S)) == HasGauges() ? S.r : 0)
        denominator(Rational{BigInt}(expected_rank_rat)) == 1 ||
            error("The code parameters imply a nonintegral additive stabilizer rank.")
        expected_rank = Int(expected_rank_rat)
        S.cache[:overcomplete] = nrows(stabs) != expected_rank
    else
        error("The current stabilizers are not equivalent to the input.")
    end
    
    # Invalidate cached signs
    set_signs!(S, S.char_vec)

    if CSSTrait(typeof(S)) == IsCSS()
        stabs_dense = _dense_code_matrix(stabs, S.F)
        flag, X_stabs, Z_stabs = robust_CSS_split(stabs_dense)
        flag || error("Detected equivalent stabilizers but is no longer CSS.")
        css_stabs = _css_symplectic_matrix(X_stabs, Z_stabs, false)
        _, X_signs, Z_signs = _determine_signs_CSS(
            css_stabs, S.char_vec, nrows(X_stabs), nrows(Z_stabs))
        S.X_stabs = is_sparse ? _sparse_code_matrix(X_stabs) : X_stabs
        S.Z_stabs = is_sparse ? _sparse_code_matrix(Z_stabs) : Z_stabs
        S.cache[:X_signs] = X_signs
        S.cache[:Z_signs] = Z_signs
    end
    
    # Clear standard form caches
    _invalidate_stabilizer_form_cache!(S)
    return nothing
end
set_stabilizers(S::AbstractSubsystemCode, stabs::CTMatrixTypes) = (S_new = deepcopy(S); set_stabilizers!(S_new, stabs); return S_new)

"""
    set_X_stabilizers!(S::AbstractSubsystemCode, X_stabs::CTMatrixTypes; trimmed::Bool = true)
"""
set_X_stabilizers!(S::T, X_stabs::CTMatrixTypes; trimmed::Bool = true) where {T <: AbstractSubsystemCode} = set_X_stabilizers!(CSSTrait(T), S, X_stabs, trimmed)
function set_X_stabilizers!(::IsCSS, S::AbstractSubsystemCode, X_stabs::CTMatrixTypes, trimmed::Bool)
    is_sparse = _is_sparse_code_matrix(X_stabs)
    X_stabs = _normalize_quantum_matrix(X_stabs)
    iszero(X_stabs) && throw(ArgumentError("The stabilizers cannot be zero."))
    order(S.F) == order(_code_matrix_base_ring(X_stabs)) || throw(ArgumentError("The stabilizers must be over the same field as the code."))
    
    if trimmed
        ncols(X_stabs) == S.n || throw(ArgumentError("Trimmed set and input of wrong size"))
        X_trimmed = X_stabs
    else
        ncols(X_stabs) == 2 * S.n || throw(ArgumentError("Trimmed not set and input of wrong size"))
        iszero(X_stabs[:, S.n + 1:end]) || throw(ArgumentError("Input is not in CSS form"))
        X_trimmed = X_stabs[:, 1:S.n]
    end

    X_trimmed = _remove_empty(X_trimmed, :rows)
    _code_matrix_base_ring(X_trimmed) == S.F ||
        (X_trimmed = change_base_ring(S.F, X_trimmed))
    is_sparse && (X_trimmed = _sparse_code_matrix(X_trimmed))

    if _additive_row_spaces_equal(S.X_stabs, X_trimmed, S.F)
        S.X_stabs = X_trimmed
        nrows(X_trimmed) != _additive_rank(X_trimmed, S.F) &&
            (S.cache[:overcomplete] = true)
    else
        error("The current stabilizers are not equivalent to the input.")
    end

    # Invalidate caches
    delete!(S.cache, :stabs)
    delete!(S.cache, :stabs_stand)
    delete!(S.cache, :P_stand)
    delete!(S.cache, :stand_r)
    delete!(S.cache, :stand_k)
    
    set_signs!(S, S.char_vec)
    return nothing
end
set_X_stabilizers!(::IsNotCSS, S::AbstractSubsystemCode, X_stabs::CTMatrixTypes, trimmed::Bool) = error("X stabilizers are only defined for CSS codes")
set_X_stabilizers(S::T, X_stabs::CTMatrixTypes; trimmed::Bool = true) where {T <: AbstractSubsystemCode} = (S_new = deepcopy(S); set_X_stabilizers!(S_new, X_stabs, trimmed=trimmed); return S_new)

"""
    set_Z_stabilizers!(S::AbstractSubsystemCode, Z_stabs::CTMatrixTypes; trimmed::Bool = true)
"""
set_Z_stabilizers!(S::T, Z_stabs::CTMatrixTypes; trimmed::Bool = true) where {T <: AbstractSubsystemCode} = set_Z_stabilizers!(CSSTrait(T), S, Z_stabs, trimmed)
function set_Z_stabilizers!(::IsCSS, S::AbstractSubsystemCode, Z_stabs::CTMatrixTypes, trimmed::Bool)
    is_sparse = _is_sparse_code_matrix(Z_stabs)
    Z_stabs = _normalize_quantum_matrix(Z_stabs)
    iszero(Z_stabs) && throw(ArgumentError("The stabilizers cannot be zero."))
    order(S.F) == order(_code_matrix_base_ring(Z_stabs)) || throw(ArgumentError("The stabilizers must be over the same field as the code."))
    
    if trimmed
        ncols(Z_stabs) == S.n || throw(ArgumentError("Trimmed set and input of wrong size"))
        Z_trimmed = Z_stabs
    else
        ncols(Z_stabs) == 2 * S.n || throw(ArgumentError("Trimmed not set and input of wrong size"))
        iszero(Z_stabs[:, 1:S.n]) || throw(ArgumentError("Input is not in CSS form"))
        Z_trimmed = Z_stabs[:, S.n + 1:end]
    end

    Z_trimmed = _remove_empty(Z_trimmed, :rows)
    _code_matrix_base_ring(Z_trimmed) == S.F ||
        (Z_trimmed = change_base_ring(S.F, Z_trimmed))
    is_sparse && (Z_trimmed = _sparse_code_matrix(Z_trimmed))
    if _additive_row_spaces_equal(S.Z_stabs, Z_trimmed, S.F)
        S.Z_stabs = Z_trimmed
        nrows(Z_trimmed) != _additive_rank(Z_trimmed, S.F) &&
            (S.cache[:overcomplete] = true)
    else
        error("The current stabilizers are not equivalent to the input.")
    end

    # Invalidate caches
    delete!(S.cache, :stabs)
    delete!(S.cache, :stabs_stand)
    delete!(S.cache, :P_stand)
    delete!(S.cache, :stand_r)
    delete!(S.cache, :stand_k)
    
    set_signs!(S, S.char_vec)
    return nothing
end
set_Z_stabilizers!(::IsNotCSS, S::AbstractSubsystemCode, Z_stabs::CTMatrixTypes, trimmed::Bool) = error("Z stabilizers are only defined for CSS codes")
set_Z_stabilizers(S::T, Z_stabs::CTMatrixTypes; trimmed::Bool = true) where {T <: AbstractSubsystemCode} = (S_new = deepcopy(S); set_Z_stabilizers!(S_new, Z_stabs, trimmed=trimmed); return S_new)

"""
    set_logicals!(S::AbstractSubsystemCode, L::CTMatrixTypes)
"""
set_logicals!(S::T, L::W) where {T <: AbstractSubsystemCode, W <: CTMatrixTypes} = set_logicals!(LogicalTrait(T), S, L)
function set_logicals!(::HasLogicals, S::AbstractSubsystemCode, L::W) where {W <: CTMatrixTypes}
    num_pairs = _logical_pair_count(S.F, S.k)
    size(L) == (2 * num_pairs, 2 * S.n) ||
        throw(ArgumentError("Provided matrix is of incorrect size for the logical space."))
    iseven(ncols(L)) || throw(ArgumentError("Expected a symplectic input but the input matrix has an odd number of columns."))
    S.F == _code_matrix_base_ring(L) || throw(ArgumentError("The logicals must be over the same field as the code."))

    current_space =
        _vcat_code_matrices(S.F, logicals_matrix(S), stabilizers(S))
    proposed_space = _vcat_code_matrices(S.F, L, stabilizers(S))
    (_additive_row_space_contains(current_space, proposed_space, S.F) &&
     _additive_row_space_contains(proposed_space, current_space, S.F)) ||
        error("The current logicals are not additively equivalent to the input.")

    is_sparse = _is_sparse_code_matrix(L)
    L_dense = _dense_code_matrix(L, S.F)
    dense_pairs = _make_pairs(deepcopy(L_dense))
    logs = Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}()
    for pair in dense_pairs
        push!(logs, is_sparse ?
            (_sparse_code_matrix(pair[1]), _sparse_code_matrix(pair[2])) : pair)
    end
    S.cache[:logicals] = logs
    S.cache[:logs_mat] = _pairs_to_matrix(logs, S.F, S.n, is_sparse)
    S.cache[:logs_alg] = :provided
    return nothing
end
set_logicals!(::HasNoLogicals, S::AbstractSubsystemCode, L::CTMatrixTypes) = error("Type $(typeof(S)) has no logicals.")
set_logicals(S::T, L::CTMatrixTypes) where {T <: AbstractSubsystemCode} = (S_new = deepcopy(S); set_logicals!(S_new, L); return S_new)

"""
    set_metacheck!(S::AbstractSubsystemCode, M::CTMatrixTypes)
"""
set_metacheck!(S::T, M::U) where {T <: AbstractSubsystemCode, U <: CTMatrixTypes} = set_metacheck!(CSSTrait(T), S, M)
set_metacheck!(::IsCSS, S::AbstractSubsystemCode, M::CTMatrixTypes) = error("Use `set_X_metacheck` and `set_Z_metacheck` for CSS codes.")
function set_metacheck!(::IsNotCSS, S::AbstractSubsystemCode, M::CTMatrixTypes)
    iszero(M * stabilizers(S)) ? (S.cache[:metacheck] = M) : error("Invalid metacheck for code")
    return nothing
end
set_metacheck(S::T, M::U) where {T <: AbstractSubsystemCode, U <: CTMatrixTypes} = (S_new = deepcopy(S); set_metacheck!(S_new, M); return S_new)

"""
    set_X_metacheck!(S::AbstractSubsystemCode, M::CTMatrixTypes)
"""
set_X_metacheck!(S::T, M::U) where {T <: AbstractSubsystemCode, U <: CTMatrixTypes} = set_X_metacheck!(CSSTrait(T), S, M)
function set_X_metacheck!(::IsCSS, S::AbstractSubsystemCode, M::CTMatrixTypes)
    iszero(M * S.X_stabs) ? (S.cache[:X_metacheck] = M) : error("Invalid metacheck for code")
    return nothing
end
set_X_metacheck!(::IsNotCSS, S::AbstractSubsystemCode, M::CTMatrixTypes) = error("Only valid for CSS codes.")
set_X_metacheck(S::T, M::U) where {T <: AbstractSubsystemCode, U <: CTMatrixTypes} = (S_new = deepcopy(S); set_X_metacheck!(S_new, M); return S_new)

"""
    set_Z_metacheck!(S::AbstractSubsystemCode, M::CTMatrixTypes)
"""
set_Z_metacheck!(S::T, M::U) where {T <: AbstractSubsystemCode, U <: CTMatrixTypes} = set_Z_metacheck!(CSSTrait(T), S, M)
function set_Z_metacheck!(::IsCSS, S::AbstractSubsystemCode, M::CTMatrixTypes)
    iszero(M * S.Z_stabs) ? (S.cache[:Z_metacheck] = M) : error("Invalid metacheck for code")
    return nothing
end
set_Z_metacheck!(::IsNotCSS, S::AbstractSubsystemCode, M::CTMatrixTypes) = error("Only valid for CSS codes.")
set_Z_metacheck(S::T, M::U) where {T <: AbstractSubsystemCode, U <: CTMatrixTypes} = (S_new = deepcopy(S); set_Z_metacheck!(S_new, M); return S_new)

# -----------------------------------------------------------------------------
# Distance Setters
# -----------------------------------------------------------------------------

"""
    set_bare_minimum_distance!(S::AbstractSubsystemCode, d::Int)
"""
set_bare_minimum_distance!(S::T, d::Int) where T <: AbstractSubsystemCode = set_bare_minimum_distance!(GaugeTrait(T), S, d)
function set_bare_minimum_distance!(::HasGauges, S::AbstractSubsystemCode, d::Int)
    0 < d <= S.n || throw(DomainError("The minimum distance of a code must be ≥ 1; received: d = $d."))
    
    u_bound_bare = get(S.cache, :u_bound_bare, S.n)
    l_bound_bare = get(S.cache, :l_bound_bare, 1)
    
    u_bound_bare < d && (@warn "The distance set is greater than the current upper bound of $u_bound_bare")
    d < l_bound_bare && (@warn "The distance set is less than the current lower bound of $l_bound_bare")
    
    S.cache[:d_bare] = d
    S.cache[:l_bound_bare] = d
    S.cache[:u_bound_bare] = d

    if !haskey(S.cache, :d_dressed)
        u_bound_dressed = get(S.cache, :u_bound_dressed, S.n)
        d < u_bound_dressed && (S.cache[:u_bound_dressed] = d)
    else
        S.cache[:d_dressed] ≤ d || (@warn "The bare distance is a bound on the dressed distance, but this is false for the new set parameters.")
    end

    if CSSTrait(typeof(S)) == IsCSS()
        !haskey(S.cache, :dx_bare) && (S.cache[:l_bound_dx_bare] = d)
        !haskey(S.cache, :dz_bare) && (S.cache[:l_bound_dz_bare] = d) 
    end
    return nothing
end
set_bare_minimum_distance!(::HasNoGauges, S::AbstractSubsystemCode, d::Int) = error("Only valid for subsytem codes; use `set_minimum_distance!` for stabilizer codes.")

"""
    set_bare_X_minimum_distance!(S::AbstractStabilizerCode, d::Int)
"""
set_bare_X_minimum_distance!(S::T, d::Int) where T <: AbstractSubsystemCode = set_bare_X_minimum_distance!(GaugeTrait(T), CSSTrait(T), S, d)
function set_bare_X_minimum_distance!(::HasGauges, ::IsCSS, S::AbstractSubsystemCode, d::Int)
    0 < d <= S.n || throw(DomainError("The minimum distance of a code must be ≥ 1; received: d = $d."))
    
    u_bound_dx_bare = get(S.cache, :u_bound_dx_bare, S.n)
    l_bound_dx_bare = get(S.cache, :l_bound_dx_bare, 1)
    
    u_bound_dx_bare < d && (@warn "The distance set is greater than the current upper bound of $u_bound_dx_bare")
    d < l_bound_dx_bare && (@warn "The distance set is less than the current lower bound of $l_bound_dx_bare")
    
    S.cache[:dx_bare] = d
    S.cache[:l_bound_dx_bare] = d
    S.cache[:u_bound_dx_bare] = d

    if !haskey(S.cache, :dx_dressed)
        u_bound_dx_dressed = get(S.cache, :u_bound_dx_dressed, S.n)
        d < u_bound_dx_dressed && (S.cache[:u_bound_dx_dressed] = d)
    else
        S.cache[:dx_dressed] ≤ d || (@warn "The bare distance is a bound on the dressed distance, but this is false for the new set parameters.")
    end

    if haskey(S.cache, :dz_bare)
        S.cache[:d_bare] = min(d, S.cache[:dz_bare])
        S.cache[:u_bound_bare] = S.cache[:d_bare]
        S.cache[:l_bound_bare] = S.cache[:d_bare]
    else
        u_bound_bare = get(S.cache, :u_bound_bare, S.n)
        d < u_bound_bare && (S.cache[:u_bound_bare] = d)
    end
    return nothing
end
set_bare_X_minimum_distance!(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode, d::Int) = error("Only valid for CSS codes")
set_bare_X_minimum_distance!(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode, d::Int) = error("Not valid for this code")
set_bare_X_minimum_distance!(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode, d::Int) = error("Only valid for subsytem codes; use `set_X_minimum_distance!` for stabilizer codes.")

"""
    set_bare_Z_minimum_distance!(S::AbstractStabilizerCode, d::Int)
"""
set_bare_Z_minimum_distance!(S::T, d::Int) where T <: AbstractSubsystemCode = set_bare_Z_minimum_distance!(GaugeTrait(T), CSSTrait(T), S, d)
function set_bare_Z_minimum_distance!(::HasGauges, ::IsCSS, S::AbstractSubsystemCode, d::Int)
    0 < d <= S.n || throw(DomainError("The minimum distance of a code must be ≥ 1; received: d = $d."))
    
    u_bound_dz_bare = get(S.cache, :u_bound_dz_bare, S.n)
    l_bound_dz_bare = get(S.cache, :l_bound_dz_bare, 1)
    
    u_bound_dz_bare < d && (@warn "The distance set is greater than the current upper bound of $u_bound_dz_bare")
    d < l_bound_dz_bare && (@warn "The distance set is less than the current lower bound of $l_bound_dz_bare")
    
    S.cache[:dz_bare] = d
    S.cache[:l_bound_dz_bare] = d
    S.cache[:u_bound_dz_bare] = d

    if !haskey(S.cache, :dz_dressed)
        u_bound_dz_dressed = get(S.cache, :u_bound_dz_dressed, S.n)
        d < u_bound_dz_dressed && (S.cache[:u_bound_dz_dressed] = d)
    else
        S.cache[:dz_dressed] ≤ d || (@warn "The bare distance is a bound on the dressed distance, but this is false for the new set parameters.")
    end

    if haskey(S.cache, :dx_bare)
        S.cache[:d_bare] = min(S.cache[:dx_bare], d)
        S.cache[:u_bound_bare] = S.cache[:d_bare]
        S.cache[:l_bound_bare] = S.cache[:d_bare]
    else
        u_bound_bare = get(S.cache, :u_bound_bare, S.n)
        d < u_bound_bare && (S.cache[:u_bound_bare] = d)
    end
    return nothing
end
set_bare_Z_minimum_distance!(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode, d::Int) = error("Only valid for CSS codes")
set_bare_Z_minimum_distance!(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode, d::Int) = error("Not valid for this code")
set_bare_Z_minimum_distance!(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode, d::Int) = error("Only valid for subsytem codes; use `set_Z_minimum_distance!` for stabilizer codes.")

"""
    set_dressed_minimum_distance!(S::AbstractSubsystemCode, d::Int)
"""
set_dressed_minimum_distance!(S::T, d::Int) where T <: AbstractSubsystemCode = set_dressed_minimum_distance!(GaugeTrait(T), S, d)
function set_dressed_minimum_distance!(::HasGauges, S::AbstractSubsystemCode, d::Int)
    0 < d <= S.n || throw(DomainError("The minimum distance of a code must be ≥ 1; received: d = $d."))
    if dimension(S) > 0 && _subsystem_singleton_is_proven(S)
        singleton = quantum_Singleton_bound(S)
        d <= singleton ||
            throw(DomainError(d,
                "The dressed distance exceeds the subsystem Singleton bound $singleton."))
    end
    
    u_bound_dressed = get(S.cache, :u_bound_dressed, S.n)
    l_bound_dressed = get(S.cache, :l_bound_dressed, 1)
    
    u_bound_dressed < d && (@warn "The distance set is greater than the current upper bound of $u_bound_dressed")
    d < l_bound_dressed && (@warn "The distance set is less than the current lower bound of $l_bound_dressed")
    
    S.cache[:d_dressed] = d
    S.cache[:l_bound_dressed] = d
    S.cache[:u_bound_dressed] = d

    if haskey(S.cache, :d_dressed) && haskey(S.cache, :d_bare)
        S.cache[:d_dressed] ≤ S.cache[:d_bare] || (@warn "The bare distance is a bound on the dressed distance, but this is false for the new set parameters.")
    end

    if CSSTrait(typeof(S)) == IsCSS()
        !haskey(S.cache, :dx_dressed) && (S.cache[:l_bound_dx_dressed] = d)
        !haskey(S.cache, :dz_dressed) && (S.cache[:l_bound_dz_dressed] = d) 
    end
    return nothing
end
set_dressed_minimum_distance!(::HasNoGauges, S::AbstractSubsystemCode, d::Int) = error("Only valid for subsytem codes; use `set_minimum_distance!` for stabilizer codes.")


"""
    set_dressed_X_minimum_distance!(S::AbstractStabilizerCode, d::Int)
"""
set_dressed_X_minimum_distance!(S::T, d::Int) where T <: AbstractSubsystemCode = set_dressed_X_minimum_distance!(GaugeTrait(T), CSSTrait(T), S, d)
function set_dressed_X_minimum_distance!(::HasGauges, ::IsCSS, S::AbstractSubsystemCode, d::Int)
    0 < d <= S.n || throw(DomainError("The minimum distance of a code must be ≥ 1; received: d = $d."))
    
    u_bound_dx_dressed = get(S.cache, :u_bound_dx_dressed, S.n)
    l_bound_dx_dressed = get(S.cache, :l_bound_dx_dressed, 1)
    
    u_bound_dx_dressed < d && (@warn "The distance set is greater than the current upper bound of $u_bound_dx_dressed")
    d < l_bound_dx_dressed && (@warn "The distance set is less than the current lower bound of $l_bound_dx_dressed")
    
    S.cache[:dx_dressed] = d
    S.cache[:l_bound_dx_dressed] = d
    S.cache[:u_bound_dx_dressed] = d

    if haskey(S.cache, :dx_bare)
        d ≤ S.cache[:dx_bare] || (@warn "The bare distance is a bound on the dressed distance, but this is false for the new set parameters.")
    end

    if haskey(S.cache, :dz_dressed)
        S.cache[:d_dressed] = min(d, S.cache[:dz_dressed])
        S.cache[:u_bound_dressed] = S.cache[:d_dressed]
        S.cache[:l_bound_dressed] = S.cache[:d_dressed]
    else
        u_bound_dressed = get(S.cache, :u_bound_dressed, S.n)
        d < u_bound_dressed && (S.cache[:u_bound_dressed] = d)
    end
    return nothing
end
set_dressed_X_minimum_distance!(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode, d::Int) = error("Only valid for CSS codes")
set_dressed_X_minimum_distance!(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode, d::Int) = error("Not valid for this code")
set_dressed_X_minimum_distance!(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode, d::Int) = error("Only valid for subsytem codes; use `set_X_minimum_distance!` for stabilizer codes.")

"""
    set_dressed_Z_minimum_distance!(S::AbstractStabilizerCode, d::Int)
"""
set_dressed_Z_minimum_distance!(S::T, d::Int) where T <: AbstractSubsystemCode = set_dressed_Z_minimum_distance!(GaugeTrait(T), CSSTrait(T), S, d)
function set_dressed_Z_minimum_distance!(::HasGauges, ::IsCSS, S::AbstractSubsystemCode, d::Int)
    0 < d <= S.n || throw(DomainError("The minimum distance of a code must be ≥ 1; received: d = $d."))
    
    u_bound_dz_dressed = get(S.cache, :u_bound_dz_dressed, S.n)
    l_bound_dz_dressed = get(S.cache, :l_bound_dz_dressed, 1)
    
    u_bound_dz_dressed < d && (@warn "The distance set is greater than the current upper bound of $u_bound_dz_dressed")
    d < l_bound_dz_dressed && (@warn "The distance set is less than the current lower bound of $l_bound_dz_dressed")
    
    S.cache[:dz_dressed] = d
    S.cache[:l_bound_dz_dressed] = d
    S.cache[:u_bound_dz_dressed] = d

    if haskey(S.cache, :dz_bare)
        d ≤ S.cache[:dz_bare] || (@warn "The bare distance is a bound on the dressed distance, but this is false for the new set parameters.")
    end

    if haskey(S.cache, :dx_dressed)
        S.cache[:d_dressed] = min(S.cache[:dx_dressed], d)
        S.cache[:u_bound_dressed] = S.cache[:d_dressed]
        S.cache[:l_bound_dressed] = S.cache[:d_dressed]
    else
        u_bound_dressed = get(S.cache, :u_bound_dressed, S.n)
        d < u_bound_dressed && (S.cache[:u_bound_dressed] = d)
    end
    return nothing
end
set_dressed_Z_minimum_distance!(::HasGauges, ::IsNotCSS, S::AbstractSubsystemCode, d::Int) = error("Only valid for CSS codes")
set_dressed_Z_minimum_distance!(::HasNoGauges, ::IsNotCSS, S::AbstractSubsystemCode, d::Int) = error("Not valid for this code")
set_dressed_Z_minimum_distance!(::HasNoGauges, ::IsCSS, S::AbstractSubsystemCode, d::Int) = error("Only valid for subsytem codes; use `set_Z_minimum_distance!` for stabilizer codes.")

#############################
     # general functions
#############################

function _quantum_dimension(F::CTFieldTypes, n::Int, additive_rank::Int)
    k = Rational{BigInt}(n) -
        Rational{BigInt}(additive_rank, degree(F))
    return denominator(k) == 1 ? Int(numerator(k)) : k
end

function _logical_pair_count(F::CTFieldTypes, k)
    count = Rational{BigInt}(k) * degree(F)
    denominator(count) == 1 ||
        error("The logical dimension is incompatible with the base field.")
    return Int(numerator(count))
end

function _lift_prime_matrix(A::CTMatrixTypes, F::CTFieldTypes)
    K = base_ring(A)
    K == F && return A
    inclusion = embed(K, F)
    return matrix(F, nrows(A), ncols(A),
        [inclusion(A[r, c]) for r in 1:nrows(A) for c in 1:ncols(A)])
end

function _lift_prime_element(x::CTFieldElem, F::CTFieldTypes)
    parent(x) == F && return x
    return embed(parent(x), F)(x)
end

function _scale_code_matrix(A::SparseMatrixCSC, scalar)
    scaled = copy(A)
    values = SparseArrays.nonzeros(scaled)
    map!(x -> scalar * x, values, values)
    return scaled
end
_scale_code_matrix(A::CTMatrixTypes, scalar) = scalar * A

function _trace_symplectic_product_matrix(
    A::CTMatrixTypes, B::CTMatrixTypes, F::CTFieldTypes
)
    ncols(A) == ncols(B) && iseven(ncols(A)) ||
        throw(ArgumentError("Symplectic matrices must have the same even length."))
    A_dense = _dense_code_matrix(A, F)
    B_dense = _dense_code_matrix(B, F)
    n = div(ncols(A_dense), 2)
    products = hcat(A_dense[:, n + 1:end], -A_dense[:, 1:n]) *
        transpose(B_dense)
    degree(F) == 1 && return products

    prime_field = _prime_subfield(F)
    return matrix(prime_field, nrows(products), ncols(products),
        [_subfield_preimage(
            prime_field, F, CodingTheory.tr(products[r, c], prime_field))
         for r in 1:nrows(products) for c in 1:ncols(products)])
end

function _additive_row_space_contains(
    big::CTMatrixTypes, small::CTMatrixTypes, F::CTFieldTypes
)
    big_expanded = _additive_expansion(big, F)
    small_expanded = _additive_expansion(small, F)
    return rank(vcat(big_expanded, small_expanded)) == rank(big_expanded)
end

function _additive_row_spaces_equal(
    A::CTMatrixTypes, B::CTMatrixTypes, F::CTFieldTypes
)
    return _additive_row_space_contains(A, B, F) &&
        _additive_row_space_contains(B, A, F)
end

function _additive_quotient_space(
    small::CTMatrixTypes, big::CTMatrixTypes, F::CTFieldTypes
)
    _additive_row_space_contains(big, small, F) ||
        throw(ArgumentError("The first additive row space is not contained in the second."))
    current = _additive_expansion(small, F)
    current_rank = rank(current)
    selected = Int[]
    big_expanded = _additive_expansion(big, F)
    for r in 1:nrows(big)
        candidate = vcat(current, big_expanded[r:r, :])
        candidate_rank = rank(candidate)
        if candidate_rank > current_rank
            push!(selected, r)
            current = candidate
            current_rank = candidate_rank
        end
    end
    source = big isa SMat ? _dense_code_matrix(big, F) : big
    quotient = isempty(selected) ? source[1:0, :] : source[selected, :]
    return big isa SMat ? _sparse_code_matrix(quotient) : quotient
end

function _additive_ambient_basis(F::CTFieldTypes, num_coordinates::Int)
    m = degree(F)
    prime_field = _prime_subfield(F)
    basis = m == 1 ? [one(F)] : first(primitive_basis(F, prime_field))
    ambient = zero_matrix(F, m * num_coordinates, num_coordinates)
    for c in 1:num_coordinates, j in 1:m
        ambient[(c - 1) * m + j, c] = basis[j]
    end
    return ambient
end

function _additive_centralizer(G::CTMatrixTypes, F::CTFieldTypes)
    ambient = _additive_ambient_basis(F, ncols(G))
    commutation = _trace_symplectic_product_matrix(G, ambient, F)
    coefficients = _rowspace_kernel(commutation)
    centralizer = _lift_prime_matrix(coefficients, F) * ambient
    return _is_sparse_code_matrix(G) ?
        _sparse_code_matrix(centralizer) : centralizer
end

function _additive_center(G::CTMatrixTypes, F::CTFieldTypes)
    commutation = _trace_symplectic_product_matrix(G, G, F)
    coefficients = _rowspace_kernel(transpose(commutation))
    center = _lift_prime_matrix(coefficients, F) *
        _dense_code_matrix(G, F)
    center = _remove_empty(center, :rows)
    return _is_sparse_code_matrix(G) ? _sparse_code_matrix(center) : center
end

function _rowspace_kernel(A::CTMatrixTypes)
    K_cols = kernel(A, side = :right)
    rnk = rank(K_cols)
    if ncols(K_cols) == rnk
        return transpose(K_cols)
    end
    F = base_ring(A)
    nr = nrows(K_cols)
    K = zero_matrix(F, rnk, nr)
    for r in 1:nr, c in 1:rnk
        !iszero(K_cols[r, c]) && (K[c, r] = K_cols[r, c])
    end
    return K
end

function _symplectic_orthogonal_complement(G::CTMatrixTypes)
    return _additive_centralizer(G, _code_matrix_base_ring(G))
end

function _subsystem_dimensions(F, n::Int, rnk_stabs::Int, rnk_gauge::Int)
    rnk_gauge >= rnk_stabs || error("Gauge group rank cannot be smaller than the stabilizer rank.")
    iseven(rnk_gauge - rnk_stabs) || error("Gauge operators outside the stabilizer group must come in symplectic pairs.")
    r = div(rnk_gauge - rnk_stabs, 2)
    k = _quantum_dimension(F, n, rnk_stabs + r)
    k >= 0 || error("The supplied generators define a negative logical dimension.")
    return k, r
end

function _pair_operators(basis::CTMatrixTypes, prefer_css::Bool)
    (iszero(basis) || nrows(basis) == 0) &&
        return Vector{Tuple{typeof(basis), typeof(basis)}}()
    is_sparse = _is_sparse_code_matrix(basis)
    work = is_sparse ?
        _dense_code_matrix(basis, _code_matrix_base_ring(basis)) : basis
    pairs = nothing
    if prefer_css
        is_css, L_X, L_Z = robust_CSS_split(work)
        if is_css
            try
                pairs = _make_CSS_pairs(L_X, L_Z)
            catch
            end
        end
    end
    isnothing(pairs) && (pairs = _make_pairs(work))
    is_sparse || return pairs
    return [(_sparse_code_matrix(pair[1]), _sparse_code_matrix(pair[2]))
            for pair in pairs]
end

function _vcat_code_matrices(F::CTFieldTypes, matrices::CTMatrixTypes...)
    any(_is_sparse_code_matrix, matrices) || return reduce(vcat, matrices)
    dense = [_dense_code_matrix(M, F) for M in matrices]
    return _sparse_code_matrix(reduce(vcat, dense))
end

function _pairs_to_matrix(pairs, F, n::Int, is_sparse::Bool)
    if isempty(pairs)
        mat = zero_matrix(F, 0, 2 * n)
        return is_sparse ? _sparse_code_matrix(mat) : mat
    end
    if is_sparse
        dense_parts = [_dense_code_matrix(part, F)
                       for pair in pairs for part in pair]
        return _sparse_code_matrix(reduce(vcat, dense_parts))
    end
    return reduce(vcat, [reduce(vcat, pairs[i]) for i in eachindex(pairs)])
end

const _STABILIZER_FORM_CACHE_KEYS = (
    :stabs_stand, :P_stand, :stand_r, :stand_k, :signs, :X_signs, :Z_signs
)

function _invalidate_stabilizer_form_cache!(S::AbstractSubsystemCode)
    for key in _STABILIZER_FORM_CACHE_KEYS
        delete!(S.cache, key)
    end
    return nothing
end

function _process_char_vec(char_vec::Union{Vector{zzModRingElem}, Missing}, p::Int, n::Int)
    if !ismissing(char_vec)
        n == length(char_vec) || throw(ArgumentError("The characteristic value is of incorrect length."))
        R, _ = residue_ring(Nemo.ZZ, p == 2 ? 4 : p)
        for s in char_vec
            modulus(s) == modulus(R) || throw(ArgumentError("Phases are not in the correct ring."))
        end
        return char_vec
    else
        return zzModRingElem[]
    end
end

function _determine_signs(S::CTMatrixTypes, char_vec::Vector{zzModRingElem})
    if isempty(char_vec)
        F = _code_matrix_base_ring(S)
        R, _ = residue_ring(Nemo.ZZ, Int(characteristic(F)) == 2 ? 4 : Int(characteristic(F)))
        return [R(0) for _ in 1:nrows(S)]
    else
        return _get_signs(S, char_vec)
    end
end

function _determine_signs_CSS(S::CTMatrixTypes, char_vec::Vector{zzModRingElem}, X_size::Int, Z_size::Int)
    if isempty(char_vec)
        F = _code_matrix_base_ring(S)
        R, _ = residue_ring(Nemo.ZZ, Int(characteristic(F)) == 2 ? 4 : Int(characteristic(F)))
        signs = [R(0) for _ in 1:nrows(S)]
        X_signs = [R(0) for _ in 1:X_size]
        Z_signs = [R(0) for _ in 1:Z_size]
    else
        signs = _get_signs(S, char_vec)
        X_signs = signs[1:X_size]
        Z_signs = signs[X_size + 1:end]
    end
    return signs, X_signs, Z_signs
end

function _get_signs(A::CTMatrixTypes, char_vec::Vector{zzModRingElem})
    R = base_ring(char_vec[1])
    nc = ncols(A)
    length(char_vec) == nc || throw(ArgumentError("Input to _get_signs is expected to be in symplectic form and of the same length as the characteristic vector."))
    
    isempty(char_vec) && return [R(0) for _ in 1:div(nc, 2)]
    signs = Vector{zzModRingElem}()
    for r in 1:nrows(A)
        parity = R(0)
        for c = 1:nc
            !iszero(A[r, c]) && (parity += char_vec[c])
        end
        push!(signs, parity)
    end
    return signs
end

function _symplectic_row(S::AbstractSubsystemCode, v::CTMatrixTypes)
    _code_matrix_base_ring(v) == S.F ||
        throw(ArgumentError("The operator must use the same field as the code."))
    size(v) == (1, 2S.n) && return v
    size(v) == (2S.n, 1) && return transpose(v)
    throw(ArgumentError(
        "Expected a 1 × $(2S.n) row or $(2S.n) × 1 column."))
end

function _in_row_space(M::CTMatrixTypes, v::CTMatrixTypes)
    iszero(v) && return true
    nrows(M) == 0 && return iszero(v)
    F = _code_matrix_base_ring(v)
    return _additive_row_space_contains(M, v, F)
end

"""
    symplectic_weight(v)

Return the Pauli weight of a symplectic vector, counting a nonzero `X` or `Z`
component on a coordinate once.
"""
function symplectic_weight(v::CTMatrixTypes)
    (nrows(v) == 1 || ncols(v) == 1) ||
        throw(ArgumentError("Expected a symplectic vector."))
    iseven(length(v)) ||
        throw(ArgumentError("A symplectic vector must have even length."))
    row = nrows(v) == 1 ? v : transpose(v)
    n = div(ncols(row), 2)
    return count(
        q -> !iszero(row[1, q]) || !iszero(row[1, n + q]), 1:n)
end

"""
    normalizer_matrix(S)
    stabilizer_centralizer_matrix(S)

Return a row basis for the symplectic centralizer (Pauli normalizer) of the
stabilizer group.
"""
normalizer_matrix(S::AbstractSubsystemCode) =
    _additive_centralizer(stabilizers(S), S.F)
stabilizer_centralizer_matrix(S::AbstractSubsystemCode) = normalizer_matrix(S)

"""
    gauge_centralizer_matrix(S)
    bare_normalizer_matrix(S)

Return a row basis for the symplectic centralizer of the gauge group.
"""
gauge_centralizer_matrix(S::AbstractSubsystemCode) =
    gauge_centralizer_matrix(GaugeTrait(typeof(S)), S)
gauge_centralizer_matrix(::HasGauges, S::AbstractSubsystemCode) =
    _additive_centralizer(gauge_group(S), S.F)
gauge_centralizer_matrix(::HasNoGauges, S::AbstractSubsystemCode) =
    normalizer_matrix(S)
bare_normalizer_matrix(S::AbstractSubsystemCode) = gauge_centralizer_matrix(S)

"""
    is_stabilizer(S, v)

Return whether `v` belongs to the row space of the stabilizer generators.
"""
is_stabilizer(S::AbstractSubsystemCode, v::CTMatrixTypes) =
    _in_row_space(stabilizers(S), _symplectic_row(S, v))

"""
    is_normalizer(S, v)

Return whether `v` commutes with every stabilizer.
"""
is_normalizer(S::AbstractSubsystemCode, v::CTMatrixTypes) =
    are_symplectic_orthogonal(stabilizers(S), _symplectic_row(S, v))

"""
    is_bare_normalizer(S, v)

Return whether `v` commutes with the full gauge group.
"""
is_bare_normalizer(S::AbstractSubsystemCode, v::CTMatrixTypes) =
    are_symplectic_orthogonal(
        GaugeTrait(typeof(S)) == HasGauges() ? gauge_group(S) : stabilizers(S),
        _symplectic_row(S, v))

"""
    robust_CSS_split(stabs::CTMatrixTypes)

Splits a stabilizer matrix into pure-X and pure-Z generators using left nullspaces.
Returns `(is_css, pure_X, pure_Z)`.
"""
function robust_CSS_split(stabs::CTMatrixTypes)
    n = div(ncols(stabs), 2)
    F = _code_matrix_base_ring(stabs)
    nrows(stabs) == 0 &&
        return true, zero_matrix(F, 0, n), zero_matrix(F, 0, n)
    S_X = stabs[:, 1:n]
    S_Z = stabs[:, n+1:end]

    X_expanded = _additive_expansion(S_X, F)
    Z_expanded = _additive_expansion(S_Z, F)
    K_Z = _rowspace_kernel(transpose(Z_expanded))
    K_X = _rowspace_kernel(transpose(X_expanded))
    pure_X = _remove_empty(
        _lift_prime_matrix(K_Z, F) * _dense_code_matrix(S_X, F), :rows)
    pure_Z = _remove_empty(
        _lift_prime_matrix(K_X, F) * _dense_code_matrix(S_Z, F), :rows)

    is_css = _additive_rank(pure_X, F) + _additive_rank(pure_Z, F) ==
        _additive_rank(stabs, F)
    if _is_sparse_code_matrix(stabs)
        pure_X = _sparse_code_matrix(pure_X)
        pure_Z = _sparse_code_matrix(pure_Z)
    end
    return is_css, pure_X, pure_Z
end

"""
    is_CSS_symplectic(stabs::CTMatrixTypes)

Returns `true` if the given symplectic matrix spans a CSS code.
"""
function is_CSS_symplectic(stabs::CTMatrixTypes)
    is_css, _, _ = robust_CSS_split(stabs)
    return is_css
end

"""
    _make_CSS_pairs(L_X::CTMatrixTypes, L_Z::CTMatrixTypes)

Creates transversal logical pairs for a CSS code.
Forces the X and Z logical bases to satisfy `L_X * transpose(L_Z) = I`.
"""
function _make_CSS_pairs(L_X::CTMatrixTypes, L_Z::CTMatrixTypes)
    F = _code_matrix_base_ring(L_X)
    k = nrows(L_X)
    n = ncols(L_X)
    
    # Commutation matrix
    X_symplectic = hcat(L_X, zero_matrix(F, nrows(L_X), n))
    Z_symplectic = hcat(zero_matrix(F, nrows(L_Z), n), L_Z)
    C = _trace_symplectic_product_matrix(X_symplectic, Z_symplectic, F)
    
    # For a valid CSS code, the bare X and Z logicals must form a non-degenerate pairing
    flag, C_inv = is_invertible_with_inverse(C)
    flag || error("Provided L_X and L_Z do not form a full, non-degenerate dual basis.")
    
    # With the package convention <(x|z),(x'|z')> = zx' - xz',
    # the minus sign makes <L_X_paired[i], L_Z[j]> = δ_ij.
    L_X_paired = _lift_prime_matrix(C_inv, F) * L_X
    
    pairs = Vector{Tuple{typeof(L_X), typeof(L_X)}}()
    z_pad = zero_matrix(F, 1, n)
    for i in 1:k
        lx = hcat(L_X_paired[i:i, :], z_pad)
        lz = hcat(z_pad, L_Z[i:i, :])
        push!(pairs, (lx, lz))
    end
    
    return pairs
end

"""
    _make_pairs(L::CTMatrixTypes)

Pairs symplectic logical operators using Symplectic Gram-Schmidt.
"""
function _make_pairs(L::CTMatrixTypes)
    iseven(nrows(L)) ||
        error("Cannot make a symplectic basis from an odd-dimensional space.")
    logs = Vector{Tuple{typeof(L), typeof(L)}}()

    while nrows(L) >= 2
        n = div(ncols(L), 2)
        F = _code_matrix_base_ring(L)
        prod = _trace_symplectic_product_matrix(L, L, F)
        partner = findfirst(c -> !iszero(prod[1, c]), 2:nrows(L))
        isnothing(partner) &&
            error("Cannot make symplectic basis; input space is degenerate.")
        partner += 1

        a = deepcopy(L[1:1, :])
        scale = inv(prod[1, partner])
        b = _scale_code_matrix(
            L[partner:partner, :], _lift_prime_element(scale, F))
        remaining = setdiff(1:nrows(L), [1, partner])
        reduced = L[remaining, :]
        for (new_row, old_row) in enumerate(remaining)
            reduced[new_row:new_row, :] =
                L[old_row:old_row, :] -
                _scale_code_matrix(
                    a,
                    _lift_prime_element(
                        prod[old_row, partner] * scale, F)
                ) +
                _scale_code_matrix(
                    b, _lift_prime_element(prod[old_row, 1], F))
        end

        push!(logs, (a, b))
        L = reduced
    end

    return logs
end

_test_logicals_relationships(S::T) where {T <: AbstractSubsystemCode} = _test_logicals_relationships(LogicalTrait(T), S)
function _test_logicals_relationships(::HasLogicals, S::AbstractSubsystemCode)
    L = logicals_matrix(S)
    prod = hcat(L[:, S.n + 1:end], -L[:, 1:S.n]) * transpose(L)
    display(prod)
    return nothing
end
_test_logicals_relationships(::HasNoLogicals, S) = error("Type $(typeof(S)) has no logicals.")

"""
    is_logical(S::AbstractSubsystemCode, v::CTMatrixTypes)
"""
is_logical(S::T, v::CTMatrixTypes) where {T <: AbstractSubsystemCode} = is_logical(LogicalTrait(T), S, v)
function is_logical(::HasLogicals, S::AbstractSubsystemCode, v::CTMatrixTypes)
    row = _symplectic_row(S, v)
    is_normalizer(S, row) || return false
    trivial_group = GaugeTrait(typeof(S)) == HasGauges() ?
        gauge_group(S) : stabilizers(S)
    return !_in_row_space(trivial_group, row)
end
is_logical(::HasNoLogicals, S::AbstractSubsystemCode, v::CTMatrixTypes) = error("Type $(typeof(S)) has no logicals.")

"""
    is_bare_logical(S, v)

Return whether `v` is a nontrivial bare logical: it centralizes the gauge
group but is not a stabilizer.
"""
function is_bare_logical(S::AbstractSubsystemCode, v::CTMatrixTypes)
    row = _symplectic_row(S, v)
    return is_bare_normalizer(S, row) && !is_stabilizer(S, row)
end

"""
    is_gauge(S::AbstractSubsystemCode, v::CTMatrixTypes)
"""
is_gauge(S::T, v::CTMatrixTypes) where {T <: AbstractSubsystemCode} = is_gauge(GaugeTrait(T), S, v)
function is_gauge(::HasGauges, S::AbstractSubsystemCode, v::CTMatrixTypes)
    return _in_row_space(gauge_group(S), _symplectic_row(S, v))
end
is_gauge(::HasNoGauges, S::AbstractSubsystemCode, v::CTMatrixTypes) = error("Type $(typeof(S)) has no gauges.")

function _minimum_nonzero_symplectic_weight(M::CTMatrixTypes, n::Int)
    F = _code_matrix_base_ring(M)
    independent = _additive_quotient_space(M[1:0, :], M, F)
    rnk = nrows(independent)
    rnk == 0 && return n + 1
    basis = _dense_code_matrix(independent, F)
    prime_field = _prime_subfield(F)
    p = Int(characteristic(F))
    cardinality = BigInt(p)^rnk
    cardinality <= typemax(Int) ||
        error("The generator group is too large for exhaustive symplectic enumeration.")
    field_elements = collect(prime_field)
    best = n + 1
    for index in 1:(Int(cardinality) - 1)
        coefficients = digits(index, base=p, pad=rnk)
        word = zero_matrix(F, 1, 2n)
        for i in 1:rnk
            coefficients[i] == 0 && continue
            scalar = _lift_prime_element(
                field_elements[coefficients[i] + 1], F)
            word += scalar * basis[i:i, :]
        end
        best = min(best, symplectic_weight(word))
        best == 1 && break
    end
    return best
end

function _minimum_CSS_group_weight(M::CTMatrixTypes, n::Int; alg::Symbol=:auto)
    is_css, X, Z = robust_CSS_split(M)
    is_css || return _minimum_nonzero_symplectic_weight(M, n)
    F = _code_matrix_base_ring(M)
    if degree(F) > 1
        distances = Int[]
        nrows(X) > 0 &&
            push!(distances, _minimum_nonzero_symplectic_weight(
                hcat(X, zero_matrix(F, nrows(X), n)), n))
        nrows(Z) > 0 &&
            push!(distances, _minimum_nonzero_symplectic_weight(
                hcat(zero_matrix(F, nrows(Z), n), Z), n))
        return isempty(distances) ? n + 1 : minimum(distances)
    end
    distances = Int[]
    if rank(X) > 0
        d_X, _ = minimum_distance(LinearCode(X); alg=alg)
        push!(distances, d_X)
    end
    if rank(Z) > 0
        d_Z, _ = minimum_distance(LinearCode(Z); alg=alg)
        push!(distances, d_Z)
    end
    return isempty(distances) ? n + 1 : minimum(distances)
end

"""
    minimum_stabilizer_weight(S; alg=:auto)

Return the minimum Pauli weight of a nonidentity stabilizer. For CSS codes this
uses the classical minimum-distance machinery on each sector. General
symplectic groups are enumerated exactly.
"""
function minimum_stabilizer_weight(
    S::AbstractSubsystemCode; alg::Symbol=:auto
)
    M = stabilizers(S)
    return CSSTrait(typeof(S)) == IsCSS() ?
        _minimum_CSS_group_weight(M, S.n; alg=alg) :
        _minimum_nonzero_symplectic_weight(M, S.n)
end

"""
    minimum_gauge_weight(S; alg=:auto)

Return the minimum Pauli weight of a nonidentity element of the gauge group.
"""
function minimum_gauge_weight(S::AbstractSubsystemCode; alg::Symbol=:auto)
    M = GaugeTrait(typeof(S)) == HasGauges() ? gauge_group(S) : stabilizers(S)
    return CSSTrait(typeof(S)) == IsCSS() ?
        _minimum_CSS_group_weight(M, S.n; alg=alg) :
        _minimum_nonzero_symplectic_weight(M, S.n)
end

function _purity_distance(S::AbstractSubsystemCode, distance)
    !ismissing(distance) && return distance
    key = GaugeTrait(typeof(S)) == HasGauges() ? :d_dressed : :d
    d = get(S.cache, key, missing)
    ismissing(d) &&
        error("Purity requires an exact minimum distance; pass `distance` or compute and cache it first.")
    return d
end

"""
    is_pure(S; distance=missing, alg=:auto)
    is_degenerate(S; distance=missing, alg=:auto)

Determine purity from the exact code distance and the minimum weight of the
stabilizer group (stabilizer codes) or gauge group (subsystem codes). This does
not require a full weight enumerator, but the minimum-group-weight computation
can still be exponential.
"""
function is_pure(
    S::AbstractSubsystemCode; distance::Union{Int, Missing}=missing,
    alg::Symbol=:auto
)
    d = _purity_distance(S, distance)
    group_distance = GaugeTrait(typeof(S)) == HasGauges() ?
        minimum_gauge_weight(S; alg=alg) :
        minimum_stabilizer_weight(S; alg=alg)
    return group_distance >= d
end
is_degenerate(S::AbstractSubsystemCode; kwargs...) = !is_pure(S; kwargs...)

"""
    syndrome(S::AbstractSubsystemCode, v::CTMatrixTypes)
"""
function syndrome(S::AbstractSubsystemCode, v::CTMatrixTypes)
    (size(v) != (2 * S.n, 1) && size(v) != (1, 2 * S.n)) &&
        throw(ArgumentError("Vector to be tested is of incorrect dimension; expected length $(2 * S.n), received: $(size(v))."))
    stabs = stabilizers(S)
    nrows(v) != 1 || return stabs * transpose(v)
    return stabs * v
end

"""
    X_syndrome(S::AbstractSubsystemCode, v::CTMatrixTypes)
"""
X_syndrome(S::T, v::CTMatrixTypes) where {T <: AbstractSubsystemCode} = X_syndrome(CSSTrait(T), S, v)
function X_syndrome(::IsCSS, S::AbstractSubsystemCode, v::CTMatrixTypes)
    length(v) == 2 * S.n && (v = v[S.n + 1:end])
    (size(v) != (S.n, 1) && size(v) != (1, S.n)) &&
        error("Vector to be tested is of incorrect dimension; expected length $(S.n), received: $(size(v)).")
    base_ring(v) == S.F || error("Vector must have the same base ring as the stabilizers.")
    nrows(v) != 1 || return S.X_stabs * transpose(v)
    return S.X_stabs * v
end
X_syndrome(::IsNotCSS, S::AbstractSubsystemCode, v::CTMatrixTypes) = error("Only valid for CSS codes.")

"""
    Z_syndrome(S::AbstractSubsystemCode, v::CTMatrixTypes)
"""
Z_syndrome(S::T, v::CTMatrixTypes) where {T <: AbstractSubsystemCode} = Z_syndrome(CSSTrait(T), S, v)
function Z_syndrome(::IsCSS, S::AbstractSubsystemCode, v::CTMatrixTypes)
    length(v) == 2 * S.n && (v = v[1:S.n])
    (size(v) != (S.n, 1) && size(v) != (1, S.n)) &&
        error("Vector to be tested is of incorrect dimension; expected length $(S.n), received: $(size(v)).")
    base_ring(v) == S.F || error("Vector must have the same base ring as the stabilizers.")
    nrows(v) != 1 || return S.Z_stabs * transpose(v)
    return S.Z_stabs * v
end
Z_syndrome(::IsNotCSS, S::AbstractSubsystemCode, v::CTMatrixTypes) = error("Only valid for CSS codes.")

"""
    promote_logicals_to_gauge!(S::AbstractSubsystemCode, pairs::Vector{Int})
"""
promote_logicals_to_gauge!(S::T, pairs::Vector{Int}) where {T <: AbstractSubsystemCode} = promote_logicals_to_gauge!(LogicalTrait(T), S, pairs)
function promote_logicals_to_gauge!(::HasLogicals, S::AbstractSubsystemCode, pairs::Vector{Int})
    pairs = sort!(unique!(pairs))
    logs = logicals(S)
    g_ops = gauges(S)
    
    append!(g_ops, logs[pairs])
    S.cache[:gauge_ops] = g_ops
    S.cache[:g_ops_mat] = reduce(vcat, [reduce(vcat, g_ops[i]) for i in 1:length(g_ops)])
    
    deleteat!(logs, pairs)
    S.cache[:logicals] = logs
    S.cache[:logs_mat] = isempty(logs) ? zero_matrix(S.F, 0, 2 * S.n) : reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
    
    S.r = S.r + length(pairs)
    p = Int(characteristic(S.F))
    
    if isinteger(S.k)
        S.k = S.k - length(pairs)
    else
        S.k = S.k / BigInt(p)^(length(pairs))
    end
    return nothing
end
promote_logicals_to_gauge!(::HasNoLogicals, S::AbstractSubsystemCode, pairs::Vector{Int}) = error("Type $(typeof(S)) has no logicals.")
promote_logicals_to_gauge(S::T, pairs::Vector{Int}) where {T <: AbstractSubsystemCode} = (S_new = deepcopy(S); promote_logicals_to_gauge!(S_new, pairs); return S_new)

"""
    promote_gauges_to_logical!(S::AbstractSubsystemCode, pairs::Vector{Int})
"""
promote_gauges_to_logical!(S::T, pairs::Vector{Int}) where {T <: AbstractSubsystemCode} = promote_gauges_to_logical!(LogicalTrait(T), S, pairs)
function promote_gauges_to_logical!(::HasLogicals, S::AbstractSubsystemCode, pairs::Vector{Int})
    pairs = sort!(unique!(pairs))
    logs = logicals(S)
    g_ops = gauges(S)
    
    append!(logs, g_ops[pairs])
    S.cache[:logicals] = logs
    S.cache[:logs_mat] = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
    
    deleteat!(g_ops, pairs)
    S.cache[:gauge_ops] = g_ops
    S.cache[:g_ops_mat] = isempty(g_ops) ? zero_matrix(S.F, 0, 2 * S.n) : reduce(vcat, [reduce(vcat, g_ops[i]) for i in 1:length(g_ops)])
    
    S.r = S.r - length(pairs)
    p = Int(characteristic(S.F))
    
    if isinteger(S.k)
        S.k = S.k + length(pairs)
    else
        S.k = S.k * BigInt(p)^(length(pairs))
    end
    return nothing
end
promote_gauges_to_logical!(::HasNoLogicals, S::AbstractSubsystemCode, pairs::Vector{Int}) = error("Type $(typeof(S)) has no logicals.")
promote_gauges_to_logical(S::T, pairs::Vector{Int}) where {T <: AbstractSubsystemCode} = (S_new = deepcopy(S); promote_gauges_to_logical!(S_new, pairs); return S_new)

"""
    swap_X_Z_logicals!(S::AbstractSubsystemCode, pairs::Vector{Int})
"""
swap_X_Z_logicals!(S::T, pairs::Vector{Int}) where {T <: AbstractSubsystemCode} = swap_X_Z_logicals!(LogicalTrait(T), S, pairs)
function swap_X_Z_logicals!(::HasLogicals, S::AbstractSubsystemCode, pairs::Vector{Int})
    pairs = sort!(unique!(pairs))
    logs = logicals(S)
    for i in pairs
        logs[i] = (logs[i][2], logs[i][1])
    end
    S.cache[:logicals] = logs
    S.cache[:logs_mat] = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
    return nothing
end
swap_X_Z_logicals!(::HasNoLogicals, S::AbstractSubsystemCode, pairs::Vector{Int}) = error("Type $(typeof(S)) has no logicals.")

"""
    swap_X_Z_gauge_operators!(S::AbstractSubsystemCode, pairs::Vector{Int})
"""
swap_X_Z_gauge_operators!(S::T, pairs::Vector{Int}) where {T <: AbstractSubsystemCode} = swap_X_Z_gauge_operators!(GaugeTrait(T), S, pairs)
function swap_X_Z_gauge_operators!(::HasGauges, S::AbstractSubsystemCode, pairs::Vector{Int})
    pairs = sort!(unique!(pairs))
    g_ops = gauges(S)
    for i in pairs
        g_ops[i] = (g_ops[i][2], g_ops[i][1])
    end
    S.cache[:gauge_ops] = g_ops
    S.cache[:g_ops_mat] = reduce(vcat, [reduce(vcat, g_ops[i]) for i in 1:length(g_ops)])
    return nothing
end
swap_X_Z_gauge_operators!(::HasNoGauges, S::AbstractSubsystemCode, pairs::Vector{Int}) = error("Type $(typeof(S)) has no gauges.")

"""
    are_equivalent(S1::T, S2::T) where T <: AbstractSubsystemCode
"""
function are_equivalent(S1::T, S2::T) where {T <: AbstractSubsystemCode}
    (S1.n == S2.n && S1.k == S2.k) || return false
    Int(order(S1.F)) == Int(order(S2.F)) || return false
    if GaugeTrait(T) == HasGauges()
        S1.r == S2.r || return false
    end

    _has_equivalent_row_spaces(stabilizers(S1), stabilizers(S2)) || return false

    if LogicalTrait(T) == HasLogicals()
        _has_equivalent_row_spaces(vcat(logicals_matrix(S1), stabilizers(S1)), vcat(logicals_matrix(S2), stabilizers(S2))) || return false
    end

    if GaugeTrait(T) == HasGauges()
        return _has_equivalent_row_spaces(vcat(gauges_matrix(S1), stabilizers(S1)), vcat(gauges_matrix(S2), stabilizers(S2)))
    else
        return true
    end
end

"""
    fix_gauge(::HasGauges, S::AbstractSubsystemCode, pair::Int, which::Symbol)
"""
fix_gauge(S::T, pair::Int, which::Symbol) where {T <: AbstractSubsystemCode} = fix_gauge(GaugeTrait(T), S, pair, which)
function fix_gauge(::HasGauges, S::AbstractSubsystemCode, pair::Int, which::Symbol)
    g_ops = gauges(S)
    if which == :X
        return augment(S, g_ops[pair][1], verbose = false)
    elseif which == :Z
        return augment(S, g_ops[pair][2], verbose = false)
    else
        throw(ArgumentError("Unknown type $which"))
    end
end
fix_gauge(::HasNoGauges, S::AbstractSubsystemCode, pair::Int, which::Symbol) = error("Type $(typeof(S)) has no gauges.")

"""
    fix_all_gauges(S::AbstractSubsystemCode; choice::Symbol = :X)

Return the stabilizer code obtained by promoting one commuting half of every
gauge pair to stabilizers.
"""
function fix_all_gauges(S::AbstractSubsystemCode; choice::Symbol = :X)
    choice ∈ (:X, :Z) || throw(ArgumentError("Choice must be :X or :Z"))
    GaugeTrait(typeof(S)) == HasGauges() ||
        throw(ArgumentError("Gauge fixing requires a subsystem code."))

    pair_index = choice == :X ? 1 : 2
    selected = [pair[pair_index] for pair in gauges(S)]
    fixed_gauges = isempty(selected) ? zero_matrix(S.F, 0, 2S.n) :
        reduce(vcat, selected)
    new_stabs = vcat(stabilizers(S), fixed_gauges)
    are_symplectic_orthogonal(new_stabs, new_stabs) ||
        error("The selected gauge generators do not form a commuting subgroup.")

    char_vec = isempty(S.char_vec) ? missing : S.char_vec
    fixed = StabilizerCode(new_stabs; char_vec=char_vec)
    if S.k > 0
        fixed.cache[:logicals] = deepcopy(logicals(S))
        fixed.cache[:logs_mat] = deepcopy(logicals_matrix(S))
        fixed.cache[:logs_alg] = :provided
    end
    fixed.cache[:gauge_fixed_from] = S
    fixed.cache[:gauge_fixing_choice] = choice

    dressed_lower = get(S.cache, :l_bound_dressed, missing)
    !ismissing(dressed_lower) && (fixed.cache[:l_bound] = dressed_lower)
    bare_upper = get(S.cache, :u_bound_bare, missing)
    !ismissing(bare_upper) && (fixed.cache[:u_bound] = bare_upper)
    return fixed
end

function show(io::IO, S::AbstractSubsystemCode)
    if isa(S.k, Integer)
        print(io, "[[$(S.n), $(S.k)")
    else
        print(io, "(($(S.n), $(S.k)")
    end
    
    if typeof(S) <: AbstractStabilizerCode
        d = get(S.cache, :d, missing)
        !ismissing(d) && print(io, ", $(d)")
    else
        print(io, ", $(S.r)")
        d_dressed = get(S.cache, :d_dressed, missing)
        !ismissing(d_dressed) && print(io, ", $(d_dressed)")
    end
    
    if isa(S.k, Integer)
        print(io, "]]_$(order(S.F))")
    else
        print(io, "))_$(order(S.F))")
    end
    
    if iszero(S.k)
        if is_CSS(S)
            println(io, typeof(S) <: AbstractStabilizerCode ? " CSS graph state" : " CSS subsystem graph state")
        else
            println(io, typeof(S) <: AbstractStabilizerCode ? " graph state" : " subsystem graph state")
        end
    else
        if isa(S, StabilizerCodeCSS)
            println(io, " CSS stabilizer code")
        elseif isa(S, StabilizerCode)
            println(io, " stabilizer code")
        elseif isa(S, SubsystemCodeCSS)
            println(io, " CSS subsystem code")
        else
            println(io, " subsystem code")
        end
    end
    
    if get(io, :compact, true) && S.n <= 30
        if isa(S, SubsystemCodeCSS) || isa(S, StabilizerCodeCSS) || isa(S, GraphStateStabilizerCSS) || isa(S, GraphStateSubsystemCSS)
            num_X = nrows(S.X_stabs)
            overcomp = get(S.cache, :overcomplete, false)
            X_sgn = X_signs(S)
            
            println(io, overcomp ? "X-stabilizer matrix (overcomplete): $num_X × $(S.n)" : "X-stabilizer matrix: $num_X × $(S.n)")
            for r in 1:num_X
                print(io, "\t chi($(X_sgn[r])) ")
                for c in 1:S.n
                    c != S.n ? print(io, "$(S.X_stabs[r, c]) ") : println(io, "$(S.X_stabs[r, c])")
                end
            end
            println(io, " ")

            num_Z = nrows(S.Z_stabs)
            Z_sgn = Z_signs(S)
            println(io, overcomp ? "Z-stabilizer matrix (overcomplete): $num_Z × $(S.n)" : "Z-stabilizer matrix: $num_Z × $(S.n)")
            for r in 1:num_Z
                print(io, "\t chi($(Z_sgn[r])) ")
                for c in 1:S.n
                    c != S.n ? print(io, "$(S.Z_stabs[r, c]) ") : println(io, "$(S.Z_stabs[r, c])")
                end
            end
        else
            stabs = stabilizers(S)
            num_stabs = nrows(stabs)
            overcomp = get(S.cache, :overcomplete, false)
            sgn = signs(S)
            
            println(io, overcomp ? "Stabilizer matrix (overcomplete): $num_stabs × $(2 * S.n)" : "Stabilizer matrix: $num_stabs × $(2 * S.n)")
            for r in 1:num_stabs
                print(io, "\t chi($(sgn[r])) ")
                for c in 1:2 * S.n
                    c != 2 * S.n ? print(io, "$(stabs[r, c]) ") : println(io, "$(stabs[r, c])")
                    c == S.n && print(io, "| ")
                end
            end
        end
    end
end

"""
    permute_code(S::AbstractSubsystemCode, σ::Union{PermGroupElem, Perm{Int}, Vector{Int}})
    permute_code!(S::AbstractSubsystemCode, σ::Union{PermGroupElem, Perm{Int}, Vector{Int}})
"""
function permute_code!(S::AbstractSubsystemCode, σ::Union{PermGroupElem, Perm{Int}, Vector{Int}})
    perm1 = transpose(permutation_matrix(S.F, typeof(σ) <: Perm ? σ.d : σ))
    perm = perm1 ⊕ perm1
    
    if is_CSS(S)
        S.X_stabs = S.X_stabs * perm1
        S.Z_stabs = S.Z_stabs * perm1
        if haskey(S.cache, :stabs)
            S.cache[:stabs] = direct_sum(S.X_stabs, S.Z_stabs)
        end
    else
        S.stabs = S.stabs * perm
    end
    
    !isempty(S.char_vec) && (S.char_vec .= data.(Array(transpose(perm))) * S.char_vec)
    
    # Invalidate Standard Form
    delete!(S.cache, :stabs_stand)
    delete!(S.cache, :P_stand)

    W = typeof(S)
    if LogicalTrait(W) == HasLogicals() && haskey(S.cache, :logicals)
        logs = S.cache[:logicals]
        for i in 1:length(logs)
            logs[i] = (logs[i][1] * perm, logs[i][2] * perm)
        end
        S.cache[:logs_mat] = S.cache[:logs_mat] * perm
    end

    if GaugeTrait(W) == HasGauges() && haskey(S.cache, :gauge_ops)
        g_ops = S.cache[:gauge_ops]
        for i in 1:length(g_ops)
            g_ops[i] = (g_ops[i][1] * perm, g_ops[i][2] * perm)
        end
        S.cache[:g_ops_mat] = S.cache[:g_ops_mat] * perm
    end
    return nothing
end
permute_code(S::AbstractSubsystemCode, σ::Union{PermGroupElem, Perm{Int}, Vector{Int}}) = (S_new = deepcopy(S); permute_code!(S_new, σ); return S_new)

"""
    augment(S::AbstractSubsystemCode, row::CTMatrixTypes; verbose::Bool = true)
"""
function _legacy_augment(S::AbstractSubsystemCode, row::CTMatrixTypes; verbose::Bool = true)
    stabs = stabilizers(S)
    typeof(stabs) == typeof(row) || throw(ArgumentError("Vector of different type than stabilizers"))
    iszero(row) && return S
    nrows(row) == 1 || throw(ArgumentError("Only one stabilizer may be passed in at a time."))

    prod = hcat(stabs[:, S.n + 1:end], -stabs[:, 1:S.n]) * transpose(row)
    if iszero(prod)
        verbose && println("Vector is already in the stabilizer group. Nothing to update.")    
        S_new = deepcopy(S)
        if is_CSS(S)
            S_new.cache[:stabs] = vcat(stabs, row)
        else
            S_new.stabs = vcat(stabs, row)
        end
        S_new.cache[:overcomplete] = true
        return S_new
    else
        stabs_to_keep = Vector{Int}()
        for i in 1:nrows(stabs)
            iszero(prod[i]) && push!(stabs_to_keep, i)
        end
        if isempty(stabs_to_keep)
            verbose && println("The vector anticommutes with all stabilizers. The new stabilizer group is just the vector.")
            new_stabs = row
        else
            update = setdiff(1:nrows(stabs), stabs_to_keep)
            if verbose
                isempty(update) ? println("No stabilizers requiring updating") : (println("Stabilizers requiring updating:"); display(update))
            end
            new_stabs = isempty(update) ? stabs : stabs[stabs_to_keep, :]
        end
    end

    if LogicalTrait(typeof(S)) == HasLogicals()
        L_mat = logicals_matrix(S)
        prod = hcat(L_mat[:, S.n + 1:end], -L_mat[:, 1:S.n]) * transpose(row)
        logs_to_keep = Vector{Int}()
        log_pairs_to_keep = Vector{Int}()
        pair = 1
        for i in 1:2:nrows(L_mat)
            if iszero(prod[i]) && iszero(prod[i + 1])
                append!(logs_to_keep, i, i + 1)
                push!(log_pairs_to_keep, pair)
            end
            pair += 1
        end

        if isempty(logs_to_keep)
            verbose && println("The vector anticommutes with all logical pairs.")
            logs_mat_new = zero_matrix(S.F, 1, 2 * S.n)
        else
            update = setdiff(1:length(logicals(S)), log_pairs_to_keep)
            if verbose
                isempty(update) ? println("No logical pairs requiring updating") : (println("Logical pairs requiring updating:"); display(update))
            end
            logs_mat_new = isempty(update) ? L_mat : L_mat[logs_to_keep, :]
        end
    else
        logs_mat_new = zero_matrix(S.F, 1, 2 * S.n)
    end

    if GaugeTrait(typeof(S)) == HasGauges()
        G_mat = gauges_matrix(S)
        prod = hcat(G_mat[:, S.n + 1:end], -G_mat[:, 1:S.n]) * transpose(row)
        g_ops_to_keep = Vector{Int}()
        g_op_pairs_to_keep = Vector{Int}()
        pair = 1
        for i in 1:2:nrows(G_mat)
            if iszero(prod[i]) && iszero(prod[i + 1])
                append!(g_ops_to_keep, i, i + 1)
                push!(g_op_pairs_to_keep, pair)
            end
            pair += 1
        end

        if isempty(g_ops_to_keep)
            verbose && println("The vector anticommutes with all gauge operator pairs.")
            gauge_ops_new = zero_matrix(S.F, 1, 2 * S.n)
        else
            update = setdiff(1:length(gauges(S)), g_op_pairs_to_keep)
            if verbose
                isempty(update) ? println("No gauge operator pairs requiring updating") : (println("Gauge operator pairs requiring updating:"); display(update))
            end
            gauge_ops_new = isempty(update) ? G_mat : G_mat[g_ops_to_keep, :]
        end
    else
        gauge_ops_new = zero_matrix(S.F, 1, 2 * S.n)
    end

    temp = _remove_empty(vcat(new_stabs, logs_mat_new, gauge_ops_new), :rows)
    temp = kernel(hcat(temp[:, S.n + 1:end], -temp[:, 1:S.n]), side = :right)
    rnk_temp = rank(temp)
    if ncols(temp) == rnk_temp
        temp = transpose(temp)
    else
        nr = nrows(temp)
        temp_tr = zero_matrix(base_ring(temp), rnk_temp, nr)
        for r in 1:nr
            for c in 1:rnk_temp
                !iszero(temp[r, c]) && (temp_tr[c, r] = temp[r, c])
            end
        end
        temp = temp_tr
    end
    
    temp = _quotient_space(new_stabs, temp, :sys_eqs)
    new_logs = _make_pairs(temp)
    return SubsystemCode(new_stabs, vcat(logs_mat_new, reduce(vcat, [reduce(vcat, new_logs[i]) for i in 1:length(new_logs)])), gauge_ops_new, char_vec = isempty(S.char_vec) ? missing : S.char_vec)
end

"""
    expurgate(S::AbstractStabilizerCode, rows::Vector{Int}; verbose::Bool = true)
"""
function _legacy_expurgate(S::AbstractSubsystemCode, rows::Vector{Int}; verbose::Bool = true)
    stabs = stabilizers(S)
    num_stabs = nrows(stabs)
    rows ⊆ 1:num_stabs || throw(ArgumentError("Argument `rows` not a subset of the number of stabilizers."))
    rows == 1:num_stabs && throw(ArgumentError("Cannot remove all stabilizers"))

    verbose && println("Removing stabilizers: $rows")
    new_stabs = stabs[setdiff(1:num_stabs, rows), :]
    temp = new_stabs
    
    if LogicalTrait(typeof(S)) == HasLogicals()
        temp = vcat(temp, logicals_matrix(S))
    end
    if GaugeTrait(typeof(S)) == HasGauges()
        temp = vcat(temp, gauges_matrix(S))
    end
    
    H = kernel(hcat(temp[:, S.n + 1:end], -temp[:, 1:S.n]), side = :right)
    rnk_H = rank(H)
    if ncols(H) == rnk_H
        H_tr = transpose(H)
    else
        nr = nrows(H)
        H_tr = zero_matrix(base_ring(H), rnk_H, nr)
        for r in 1:nr
            for c in 1:rnk_H
                !iszero(H[r, c]) && (H_tr[c, r] = H[r, c])
            end
        end
    end

    new_logs = _quotient_space(new_stabs, H_tr, :sys_eqs)
    if iszero(new_logs)
        verbose && println("No new logicals need to be added")
        S_new = deepcopy(S)
        if is_CSS(S)
            S_new.cache[:stabs] = new_stabs
        else
            S_new.stabs = new_stabs
        end
        rank(new_stabs) == nrows(new_stabs) ? (S_new.cache[:overcomplete] = false) : (S_new.cache[:overcomplete] = true)
        return S_new
    else
        new_log_pairs = _make_pairs(new_logs)
        verbose && (println("New logical pairs:"); display(new_log_pairs))
        
        char_vec = isempty(S.char_vec) ? missing : S.char_vec
        new_logs_mat = reduce(vcat, [reduce(vcat, new_log_pairs[i]) for i in 1:length(new_log_pairs)])
        L_mat_combined = haskey(S.cache, :logs_mat) ? vcat(S.cache[:logs_mat], new_logs_mat) : new_logs_mat
        
        if GaugeTrait(typeof(S)) == HasGauges()
            return SubsystemCode(new_stabs, L_mat_combined, gauges_matrix(S), char_vec = char_vec)
        else
            S_new = StabilizerCode(new_stabs, char_vec = char_vec)
            set_logicals!(S_new, L_mat_combined)
            return S_new
        end
    end
end

function _standard_form_stabilizer(M::CTMatrixTypes)
    stabs = deepcopy(M)
    _rref_no_col_swap!(stabs, 1:size(stabs, 1), 1:size(stabs, 2))
    nr = size(stabs, 1)
    for i in size(stabs, 1):-1:1
        nr = i
        iszero(stabs[i, :]) || break
    end
    if nr != size(stabs, 1)
        stabs = stabs[1:nr, :]
    end

    n = div(size(stabs, 2), 2)
    k = n - nr
    r, P1 = _rref_symp_col_swap!(stabs, 1:nr, 1:n)
    _, P2 = _rref_symp_col_swap!(stabs, (r + 1):nr, (n + r + 1):2n)

    P = if ismissing(P1) && ismissing(P2)
        missing
    elseif ismissing(P1)
        P2
    elseif ismissing(P2)
        P1
    else
        P2 * P1
    end
    return stabs, P, r, k, n - k
end

function _logicals_standard_form(stabs::CTMatrixTypes, n::Int, k::Int, r::Int, P::Union{Missing, CTMatrixTypes})
    F = base_ring(stabs)
    F_one = F(1)
    logs = zero_matrix(F, 2k, 2n)

    E = stabs[r + 1:n - k, 2n - k + 1:2n]
    C1 = stabs[1:r, n + r + 1:2n - k]
    C1_E = C1 * E
    for i in 1:k
        logs[i, n - k + i] = F_one
        logs[k + i, 2n - k + i] = F_one

        for j in 1:(n - k - r)
            logs[i, j + r] = stabs[r + j, 2n - k + i]
        end

        for j in 1:r
            logs[i, n + j] = C1_E[j, i] + stabs[j, 2n - k + i]
            logs[k + i, n + j] = stabs[j, n - k + i]
        end
    end
    return ismissing(P) ? logs : logs * P
end
