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

    iszero(G) && throw(ArgumentError("The gauge matrix is empty."))
    G = _remove_empty(G, :rows)

    F = base_ring(G)
    p = Int(characteristic(F))
    n = div(ncols(G), 2)
    clean_char_vec = ismissing(char_vec) ? zzModRingElem[] : _process_char_vec(char_vec, p, 2 * n)

    is_sparse = G isa SparseMatrixCSC
    G_dense = is_sparse ? matrix(F, G) : G

    # 1. Compute Stabilizers
    P = hcat(G_dense[:, n + 1:end], -G_dense[:, 1:n]) * transpose(G_dense)
    K_cols = kernel(transpose(P), side = :right)
    rnk_K = rank(K_cols)
    if ncols(K_cols) == rnk_K
        K = transpose(K_cols)
    else
        nr = nrows(K_cols)
        K = zero_matrix(F, rnk_K, nr)
        for r in 1:nr, c in 1:rnk_K
            !iszero(K_cols[r, c]) && (K[c, r] = K_cols[r, c])
        end
    end
    stabs_dense = _remove_empty(K * G_dense, :rows)
    iszero(stabs_dense) && error("Error computing the stabilizer group of the subsystem code; ker G ∩ G has dimension zero.")

    if rank(stabs_dense) == rank(G_dense)
        println("Stabilizer code detected.")
        return StabilizerCode(is_sparse ? sparse(stabs_dense) : stabs_dense, char_vec = char_vec)
    end

    # 2. Robust CSS Split of Stabilizers
    is_css_S, X_stabs_dense, Z_stabs_dense = robust_CSS_split(stabs_dense)

    # 3. Compute Bare Logicals
    ker_G_cols = kernel(hcat(G_dense[:, n + 1:end], -G_dense[:, 1:n]), side = :right)
    rnk_ker_G = rank(ker_G_cols)
    if ncols(ker_G_cols) == rnk_ker_G
        ker_G = transpose(ker_G_cols)
    else
        nr = nrows(ker_G_cols)
        ker_G = zero_matrix(F, rnk_ker_G, nr)
        for r in 1:nr, c in 1:rnk_ker_G
            !iszero(ker_G_cols[r, c]) && (ker_G[c, r] = ker_G_cols[r, c])
        end
    end

    BL = _quotient_space(ker_G, stabs_dense, :sys_eqs)
    if iszero(BL) || nrows(BL) == 0
        bare_logs_dense = Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}()
        logs_mat_dense = zero_matrix(F, 0, 2n)
    else
        is_css_BL, L_X, L_Z = robust_CSS_split(BL)
        if is_css_S && is_css_BL
            bare_logs_dense = _make_CSS_pairs(L_X, L_Z)
        else
            bare_logs_dense = _make_pairs(BL)
        end
        logs_mat_dense = reduce(vcat, [reduce(vcat, bare_logs_dense[i]) for i in 1:length(bare_logs_dense)])
    end

    # 4. Compute Gauge Operators
    GO = _quotient_space(G_dense, stabs_dense, :sys_eqs)
    if iszero(GO) || nrows(GO) == 0
        gauge_ops_dense = Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}()
        g_ops_mat_dense = zero_matrix(F, 0, 2n)
    else
        is_css_GO, G_X, G_Z = robust_CSS_split(GO)
        if is_css_S && is_css_GO
            gauge_ops_dense = _make_CSS_pairs(G_X, G_Z)
        else
            gauge_ops_dense = _make_pairs(GO)
        end
        g_ops_mat_dense = reduce(vcat, [reduce(vcat, gauge_ops_dense[i]) for i in 1:length(gauge_ops_dense)])
    end

    # 5. Fast Dimensions
    rnk_stabs = rank(stabs_dense)
    r = length(gauge_ops_dense)
    top = BigInt(order(F))^n
    k = top // BigInt(p)^(rnk_stabs + r)
    isinteger(k) && (k = round(Int, log(BigInt(p), k));)

    # 6. Sparse Reseeding & Cache Construction
    stabs = is_sparse ? sparse(stabs_dense) : stabs_dense
    g_ops_mat = is_sparse ? sparse(g_ops_mat_dense) : g_ops_mat_dense
    
    gauge_ops = Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}()
    for i in 1:length(gauge_ops_dense)
        push!(gauge_ops, (is_sparse ? sparse(gauge_ops_dense[i][1]) : gauge_ops_dense[i][1], 
                          is_sparse ? sparse(gauge_ops_dense[i][2]) : gauge_ops_dense[i][2]))
    end

    cache = Dict{Symbol, Any}(
        :stabs => stabs,
        :g_ops_mat => g_ops_mat
    )

    if k > 0
        logs_mat = is_sparse ? sparse(logs_mat_dense) : logs_mat_dense
        bare_logs = Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}()
        for i in 1:length(bare_logs_dense)
            push!(bare_logs, (is_sparse ? sparse(bare_logs_dense[i][1]) : bare_logs_dense[i][1], 
                              is_sparse ? sparse(bare_logs_dense[i][2]) : bare_logs_dense[i][2]))
        end
        cache[:logicals] = bare_logs
        cache[:logs_mat] = logs_mat
    end

    if is_css_S
        X_stabs = is_sparse ? sparse(X_stabs_dense) : X_stabs_dense
        Z_stabs = is_sparse ? sparse(Z_stabs_dense) : Z_stabs_dense
        return SubsystemCodeCSS(F, n, k, r, X_stabs, Z_stabs, gauge_ops, clean_char_vec, cache)
    else
        return SubsystemCode(F, n, k, r, stabs, gauge_ops, clean_char_vec, cache)
    end
end

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

    iszero(S) && error("The stabilizer matrix is empty.")
    S = _remove_empty(S, :rows)
    n = div(ncols(S), 2)
    F = base_ring(S)
    p = Int(characteristic(F))

    is_sparse = S isa SparseMatrixCSC
    S_dense = is_sparse ? matrix(F, S) : S
    L_dense = is_sparse ? matrix(F, L) : L
    G_dense = is_sparse ? matrix(F, G) : G

    are_symplectic_orthogonal(S_dense, S_dense) || error("The given stabilizers are not symplectic orthogonal.")

    is_css_S, X_stabs_dense, Z_stabs_dense = robust_CSS_split(S_dense)

    # logicals validation
    if iszero(L_dense) || nrows(L_dense) == 0
        log_pairs_dense = Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}()
        logs_mat_dense = zero_matrix(F, 0, 2n)
    else
        L_dense = _remove_empty(L_dense, :rows)
        are_symplectic_orthogonal(S_dense, L_dense) || error("Logicals do not commute with the code.")
        prod = hcat(L_dense[:, n + 1:end], -L_dense[:, 1:n]) * transpose(L_dense)
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
        prod = hcat(G_dense[:, n + 1:end], -G_dense[:, 1:n]) * transpose(G_dense)
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
    rnk_G = rank(G_dense)
    top = BigInt(order(F))^n
    r = div(rnk_G, 2)
    rnk_S = rank(S_dense)
    k = top // BigInt(p)^(rnk_S + r)
    isinteger(k) && (k = round(Int, log(BigInt(p), k));)
    
    S_final = is_sparse ? sparse(S_dense) : S_dense
    g_ops_mat_final = is_sparse ? sparse(g_ops_mat_dense) : g_ops_mat_dense
    logs_mat_final = is_sparse ? sparse(logs_mat_dense) : logs_mat_dense

    log_pairs = Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}()
    for i in 1:length(log_pairs_dense)
        push!(log_pairs, (is_sparse ? sparse(log_pairs_dense[i][1]) : log_pairs_dense[i][1], 
                          is_sparse ? sparse(log_pairs_dense[i][2]) : log_pairs_dense[i][2]))
    end

    g_ops_pairs = Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}()
    for i in 1:length(g_ops_pairs_dense)
        push!(g_ops_pairs, (is_sparse ? sparse(g_ops_pairs_dense[i][1]) : g_ops_pairs_dense[i][1], 
                            is_sparse ? sparse(g_ops_pairs_dense[i][2]) : g_ops_pairs_dense[i][2]))
    end

    cache = Dict{Symbol, Any}(
        :stabs => S_final,
        :logicals => log_pairs,
        :logs_mat => logs_mat_final,
        :g_ops_mat => g_ops_mat_final,
        :overcomplete => nrows(S_dense) > rnk_S
    )
    
    if is_css_S
        X_stabs = is_sparse ? sparse(X_stabs_dense) : X_stabs_dense
        Z_stabs = is_sparse ? sparse(Z_stabs_dense) : Z_stabs_dense
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

# CSS construction, Euclidean and Hermitian
# min dist is min dressed logical operator weight

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
cardinality(S::AbstractSubsystemCode) = BigInt(characteristic(S.F))^(S.n - S.k - get(S.cache, :r, S.r))

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
        if !haskey(S.cache, :stabs_stand)
            # Standard form algorithms require dense matrices
            stabs_base = stabilizers(S)
            is_sparse = stabs_base isa SparseMatrixCSC
            stabs_dense = _dense_code_matrix(stabs_base, S.F)
            
            stabs_stand, P_stand, stand_r, stand_k, _ = _standard_form_stabilizer(stabs_dense)
            
            S.cache[:stabs_stand] = is_sparse ? sparse(stabs_stand) : stabs_stand
            S.cache[:P_stand] = is_sparse ? sparse(P_stand) : P_stand
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
    
    # If not in cache, we lazily compute for stabilizer codes using standard form
    if GaugeTrait(typeof(S)) == HasNoGauges()
        stabs_stand = stabilizers(S; standform = true)
        r = S.cache[:stand_r]
        k = S.cache[:stand_k]
        P = S.cache[:P_stand]
        
        # Dense Compute / Sparse Store
        stabs_dense = stabs_stand isa SparseMatrixCSC ? matrix(S.F, stabs_stand) : stabs_stand
        P_dense = (P === missing) ? missing : (P isa SparseMatrixCSC ? matrix(S.F, P) : P)
        
        logs_dense = _make_pairs(_logicals_standard_form(stabs_dense, S.n, k, r, P_dense))
        
        is_sparse = stabs_stand isa SparseMatrixCSC
        logs = Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}()
        for i in 1:length(logs_dense)
            push!(logs, (is_sparse ? sparse(logs_dense[i][1]) : logs_dense[i][1], 
                         is_sparse ? sparse(logs_dense[i][2]) : logs_dense[i][2]))
        end
        
        S.cache[:logicals] = logs
        S.cache[:logs_mat] = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
        return logs
    else
        error("Logicals not found in cache. For subsystem codes, logicals should be seeded during initialization.")
    end
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
    stabs_stand = stabilizers(S; standform = true)
    
    stabs_dense = stabs_stand isa SparseMatrixCSC ? matrix(S.F, stabs_stand) : stabs_stand
    P = S.cache[:P_stand]
    P_dense = (P === missing) ? missing : (P isa SparseMatrixCSC ? matrix(S.F, P) : P)
    
    logs_stand = _logicals_standard_form(stabs_dense, S.n, S.cache[:stand_k], S.cache[:stand_r], P_dense)
    return stabs_stand isa SparseMatrixCSC ? sparse(logs_stand) : logs_stand
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
gauge_group(::HasGauges, S::AbstractSubsystemCode) = vcat(stabilizers(S), gauges_matrix(S))
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
num_X_stabs(::IsCSS, S::AbstractSubsystemCode) = nrows(S.X_stabs)
num_X_stabs(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    num_Z_stabs(S::AbstractSubsystemCode)

Return the number of `Z` stabilizers of the CSS code.
"""
num_Z_stabs(S::T) where {T <: AbstractSubsystemCode} = num_Z_stabs(CSSTrait(T), S)
num_Z_stabs(::IsCSS, S::AbstractSubsystemCode) = nrows(S.Z_stabs)
num_Z_stabs(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    character_vector(S::AbstractSubsystemCode)

Return the character vector of the code.
"""
character_vector(S::AbstractSubsystemCode) = S.char_vec

# TODO: quantum Singletonbound k <= n - 2d + 2
# MDS/optimal for subsystem codes: k + r <= n - 2d + 2

"""
    relative_distance(S::AbstractSubsystemCode)

Return the relative minimum distance, `δ = d / n` of the code if `d` is known,
otherwise errors.
"""
function relative_distance(S::AbstractSubsystemCode)
    !ismissing(S.d) || error("Missing minimum distance for this code.")
    return S.d / S.n
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
    if LogicalTrait(T) == HasNoLogicals
        error("Type $T has no logicals.")
    elseif GaugeTrait(T) == HasNoGauges
        error("Type $T has no gauges.")
    end
    return S.logicals ∪ S.gauge_ops
end
dressed_operators(S::AbstractSubsystemCode) = dressed(S)
dressed_logicals(S::AbstractSubsystemCode) = dressed(S)

"""
    bare_minimum_distance_lower_bound(S::AbstractSubsystemCode)
Return the currently stored lower bound on the bare minimum distance.
"""
bare_minimum_distance_lower_bound(S::T) where T <: AbstractSubsystemCode =
    bare_minimum_distance_lower_bound(GaugeTrait(T), S)
bare_minimum_distance_lower_bound(::HasGauges, S::AbstractSubsystemCode) = S.l_bound_bare
bare_minimum_distance_lower_bound(::HasNoGauges, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes; use `minimum_distance_lower_bound` for stabilizer codes.")

"""
    bare_minimum_distance_upper_bound(S::AbstractSubsystemCode)
Return the currently stored upper bound on the bare minimum distance.
"""
bare_minimum_distance_upper_bound(S::T) where T <: AbstractSubsystemCode =
    bare_minimum_distance_upper_bound(GaugeTrait(T), S)
bare_minimum_distance_upper_bound(::HasGauges, S::AbstractSubsystemCode) = S.u_bound_bare
bare_minimum_distance_upper_bound(::HasNoGauges, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes; use `minimum_distance_lower_bound` for stabilizer codes.")

"""
    dressed_minimum_distance_lower_bound(S::AbstractSubsystemCode)
Return the currently stored lower bound on the dressed minimum distance.
"""
dressed_minimum_distance_lower_bound(S::T) where T <: AbstractSubsystemCode =
    dressed_minimum_distance_lower_bound(GaugeTrait(T), S)
dressed_minimum_distance_lower_bound(::HasGauges, S::AbstractSubsystemCode) = S.l_bound_dressed
dressed_minimum_distance_lower_bound(::HasNoGauges, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes; use `minimum_distance_lower_bound` for stabilizer codes.")

"""
    dressed_minimum_distance_upper_bound(S::AbstractSubsystemCode)
Return the currently stored upper bound on the dressed minimum distance.
"""
dressed_minimum_distance_upper_bound(S::T) where T <: AbstractSubsystemCode =
    dressed_minimum_distance_upper_bound(GaugeTrait(T), S)
dressed_minimum_distance_upper_bound(::HasGauges, S::AbstractSubsystemCode) = S.u_bound_dressed
dressed_minimum_distance_upper_bound(::HasNoGauges, S::AbstractSubsystemCode) =
    error("Only valid for subsystem codes; use `minimum_distance_lower_bound` for stabilizer codes.")

"""
    bare_X_minimum_distance_lower_bound(S::AbstractSubsystemCode)
Return the currently stored lower bound on the bare `X`-minimum distance.
"""
bare_X_minimum_distance_lower_bound(S::T) where T <: AbstractSubsystemCode = bare_X_minimum_distance_lower_bound(GaugeTrait(T), CSSTrait(T), S)
bare_X_minimum_distance_lower_bound(::HasGauges, ::IsCSS, S::AbstractSubsystemCode) =
    S.l_bound_dx_bare
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
    S.u_bound_dx_bare
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
    S.l_bound_dx_dressed
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
    S.u_bound_dx_dressed
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
    S.l_bound_dZ_bare
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
    S.u_bound_dz_bare
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
    S.l_bound_dz_dressed
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
    S.u_bound_dz_dressed
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
    iszero(stabs) && throw(ArgumentError("The stabilizers cannot be zero."))
    order(S.F) == order(base_ring(stabs)) || throw(ArgumentError("The stabilizers must be over the same field as the code."))

    stabs = _remove_empty(stabs, :rows)
    stabs = change_base_ring(S.F, stabs)
    if _has_equivalent_row_spaces(stabilizers(S), stabs)
        if is_CSS(S)
            S.cache[:stabs] = stabs
        else
            S.stabs = stabs
        end
        
        expected_rank = GaugeTrait(typeof(S)) == HasNoGauges() ? (S.n - S.k) : (S.n - S.k - S.r)
        nrows(stabs) != expected_rank ? (S.cache[:overcomplete] = true) : (S.cache[:overcomplete] = false)
    else
        error("The current stabilizers are not equivalent to the input.")
    end
    
    # Invalidate cached signs
    set_signs!(S, S.char_vec)

    if CSSTrait(typeof(S)) == IsCSS()
        flag, X_stabs, X_signs, Z_stabs, Z_signs = _is_CSS_symplectic(stabs, signs(S), true)
        flag || error("Detected equivalent stabilizers but is no longer CSS.")
        S.X_stabs = X_stabs
        S.Z_stabs = Z_stabs
        S.cache[:X_signs] = X_signs
        S.cache[:Z_signs] = Z_signs
    end
    
    # Clear standard form caches
    delete!(S.cache, :stabs_stand)
    delete!(S.cache, :P_stand)
    delete!(S.cache, :stand_r)
    delete!(S.cache, :stand_k)
    return nothing
end
set_stabilizers(S::AbstractSubsystemCode, stabs::CTMatrixTypes) = (S_new = deepcopy(S); set_stabilizers!(S_new, stabs); return S_new)

"""
    set_X_stabilizers!(S::AbstractSubsystemCode, X_stabs::CTMatrixTypes; trimmed::Bool = true)
"""
set_X_stabilizers!(S::T, X_stabs::CTMatrixTypes; trimmed::Bool = true) where {T <: AbstractSubsystemCode} = set_X_stabilizers!(CSSTrait(T), S, X_stabs, trimmed)
function set_X_stabilizers!(::IsCSS, S::AbstractSubsystemCode, X_stabs::CTMatrixTypes, trimmed::Bool)
    iszero(X_stabs) && throw(ArgumentError("The stabilizers cannot be zero."))
    order(S.F) == order(base_ring(X_stabs)) || throw(ArgumentError("The stabilizers must be over the same field as the code."))
    
    if trimmed
        ncols(X_stabs) == S.n || throw(ArgumentError("Trimmed set and input of wrong size"))
        X_trimmed = X_stabs
    else
        ncols(X_stabs) == 2 * S.n || throw(ArgumentError("Trimmed not set and input of wrong size"))
        iszero(X_stabs[:, S.n + 1:end]) || throw(ArgumentError("Input is not in CSS form"))
        X_trimmed = X_stabs[:, 1:S.n]
    end

    X_trimmed = _remove_empty(X_trimmed, :rows)
    X_trimmed = change_base_ring(S.F, X_trimmed)
    
    if _has_equivalent_row_spaces(S.X_stabs, X_trimmed)
        S.X_stabs = X_trimmed
        nrows(X_trimmed) != rank(X_trimmed) && (S.cache[:overcomplete] = true)
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
    iszero(Z_stabs) && throw(ArgumentError("The stabilizers cannot be zero."))
    order(S.F) == order(base_ring(Z_stabs)) || throw(ArgumentError("The stabilizers must be over the same field as the code."))
    
    if trimmed
        ncols(Z_stabs) == S.n || throw(ArgumentError("Trimmed set and input of wrong size"))
        Z_trimmed = Z_stabs
    else
        ncols(Z_stabs) == 2 * S.n || throw(ArgumentError("Trimmed not set and input of wrong size"))
        iszero(Z_stabs[:, 1:S.n]) || throw(ArgumentError("Input is not in CSS form"))
        Z_trimmed = Z_stabs[:, S.n + 1:end]
    end

    Z_trimmed = _remove_empty(Z_trimmed, :rows)
    Z_trimmed = change_base_ring(S.F, Z_trimmed)
    if _has_equivalent_row_spaces(S.Z_stabs, Z_trimmed)
        S.Z_stabs = Z_trimmed
        nrows(Z_trimmed) != rank(Z_trimmed) && (S.cache[:overcomplete] = true)
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
    size(L) == (2 * S.k, 2 * S.n) || throw(ArgumentError("Provided matrix is of incorrect size for the logical space."))
    iseven(ncols(L)) || throw(ArgumentError("Expected a symplectic input but the input matrix has an odd number of columns."))
    S.F == base_ring(L) || throw(ArgumentError("The logicals must be over the same field as the code."))
    
    _has_equivalent_row_spaces(vcat(logicals_matrix(S), stabilizers(S)), vcat(L, stabilizers(S))) || error("The current logicals are not equivalent to the input.")

    prod = hcat(L[:, S.n + 1:end], -L[:, 1:S.n]) * transpose(L)
    iszero(prod) && throw(ArgumentError("Provided logicals should not be symplectic self-orthogonal."))
    nc_pr = ncols(prod)
    prod_Jul = _Flint_matrix_to_Julia_int_matrix(prod)
    cols = [sum(prod_Jul[:, i]) for i in 1:nc_pr]
    sum(cols) == nc_pr || throw(ArgumentError("Incorrect commutation relationships between provided logicals."))

    F = base_ring(L)
    F_one = F(1)
    logs = Vector{Tuple{W, W}}()
    if Int(order(F)) != 2
        while nrows(L) >= 2
            y = findfirst(x -> x > 0, prod_Jul[:, 1])
            y = [F(prod[y, 1]), y]
            if y[1] != F_one
                push!(logs, (L[1:1, :], y[1]^-1 * L[y[2]:y[2], :]))
            else
                push!(logs, (L[1:1, :], L[y[2]:y[2], :]))
            end
            L = L[setdiff(1:size(L, 1), [1, y[2]]), :]
        end
    else
        while nrows(L) >= 2
            y = findfirst(x -> x > 0, prod_Jul[:, 1])
            y = [F(prod[y, 1]), y]
            push!(logs, (L[1:1, :], L[y[2]:y[2], :]))
            L = L[setdiff(1:size(L, 1), [1, y[2]]), :]
        end
    end
    
    S.cache[:logicals] = logs
    S.cache[:logs_mat] = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
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
        R, _ = residue_ring(Nemo.ZZ, Int(characteristic(base_ring(S))) == 2 ? 4 : Int(characteristic(base_ring(S))))
        return [R(0) for _ in 1:nrows(S)]
    else
        return _get_signs(S, char_vec)
    end
end

function _determine_signs_CSS(S::CTMatrixTypes, char_vec::Vector{zzModRingElem}, X_size::Int, Z_size::Int)
    if isempty(char_vec)
        R, _ = residue_ring(Nemo.ZZ, Int(characteristic(base_ring(S))) == 2 ? 4 : Int(characteristic(base_ring(S))))
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

"""
    robust_CSS_split(stabs::CTMatrixTypes)

Splits a stabilizer matrix into pure-X and pure-Z generators using left nullspaces.
Returns `(is_css, pure_X, pure_Z)`.
"""
function robust_CSS_split(stabs::CTMatrixTypes)
    n = div(ncols(stabs), 2)
    F = base_ring(stabs)
    S_X = stabs[:, 1:n]
    S_Z = stabs[:, n+1:end]
    
    # 1. Pure X stabilizers: Combinations of rows where the Z part vanishes
    K_Z_cols = kernel(transpose(S_Z), side=:right)
    rnk_K_Z = rank(K_Z_cols)
    if ncols(K_Z_cols) == rnk_K_Z
        K_Z = transpose(K_Z_cols)
    else
        # Flint bug workaround
        nr = nrows(K_Z_cols)
        K_Z = zero_matrix(F, rnk_K_Z, nr)
        for r in 1:nr, c in 1:rnk_K_Z
            !iszero(K_Z_cols[r, c]) && (K_Z[c, r] = K_Z_cols[r, c])
        end
    end
    pure_X = _remove_empty(K_Z * S_X, :rows)
    
    # 2. Pure Z stabilizers: Combinations of rows where the X part vanishes
    K_X_cols = kernel(transpose(S_X), side=:right)
    rnk_K_X = rank(K_X_cols)
    if ncols(K_X_cols) == rnk_K_X
        K_X = transpose(K_X_cols)
    else
        nr = nrows(K_X_cols)
        K_X = zero_matrix(F, rnk_K_X, nr)
        for r in 1:nr, c in 1:rnk_K_X
            !iszero(K_X_cols[r, c]) && (K_X[c, r] = K_X_cols[r, c])
        end
    end
    pure_Z = _remove_empty(K_X * S_Z, :rows)
    
    # 3. Validation
    is_css = (rank(pure_X) + rank(pure_Z) == rank(stabs))
    
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
    F = base_ring(L_X)
    k = nrows(L_X)
    n = ncols(L_X)
    
    # Commutation matrix
    C = L_X * transpose(L_Z)
    
    # For a valid CSS code, the bare X and Z logicals must form a non-degenerate pairing
    flag, C_inv = is_invertible_with_inverse(C)
    flag || error("Provided L_X and L_Z do not form a full, non-degenerate dual basis.")
    
    # Align L_X to exactly match L_Z
    L_X_paired = C_inv * L_X
    
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
    F = base_ring(L)
    n = div(ncols(L), 2)
    logs = Vector{Tuple{typeof(L), typeof(L)}}()
    
    while nrows(L) >= 2
        prod = hcat(L[:, n + 1:end], -L[:, 1:n]) * transpose(L)
        num_prod = ncols(prod)
        first = 0
        
        for c in 1:num_prod
            if !iszero(prod[1, c])
                if iszero(first)
                    first = c
                    if !isone(prod[1, c])
                        L[first:first, :] *= inv(prod[1, c])
                    end
                else
                    L[c:c, :] += inv(prod[1, c]) * L[first:first, :]
                end
            end
        end
        
        iszero(first) && error("Cannot make symplectic basis; input logicals are degenerate.")
        
        for c in 2:num_prod
            if !iszero(prod[first, c])
                L[c:c, :] += inv(prod[first, c]) * L[1:1, :]
            end
        end
        
        push!(logs, (L[1:1, :], L[first:first, :]))
        L = L[setdiff(1:nrows(L), [1, first]), :]
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
    L = logicals_matrix(S)
    nc = ncols(L)
    are_symplectic_orthogonal(stabilizers(S), v) || return false
    size(v) == (1, nc) && (return !are_symplectic_orthogonal(L, v))
    size(v) == (nc, 1) && (return !are_symplectic_orthogonal(L, transpose(v)))
    throw(ArgumentError("Vector to be tested is of incorrect dimension."))
end
is_logical(::HasNoLogicals, S::AbstractSubsystemCode, v::CTMatrixTypes) = error("Type $(typeof(S)) has no logicals.")

"""
    is_gauge(S::AbstractSubsystemCode, v::CTMatrixTypes)
"""
is_gauge(S::T, v::CTMatrixTypes) where {T <: AbstractSubsystemCode} = is_gauge(GaugeTrait(T), S, v)
function is_gauge(::HasGauges, S::AbstractSubsystemCode, v::CTMatrixTypes)
    G_mat = gauges_matrix(S)
    nc = ncols(G_mat)
    are_symplectic_orthogonal(stabilizers(S), v) || return false
    size(v) == (1, nc) && (return !iszero(G_mat * transpose(v)))
    size(v) == (nc, 1) && (return !iszero(G_mat * v))
    throw(ArgumentError("Vector to be tested is of incorrect dimension."))
end
is_gauge(::HasNoGauges, S::AbstractSubsystemCode, v::CTMatrixTypes) = error("Type $(typeof(S)) has no gauges.")

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

Returns a `GaugeFixedCode` by promoting a maximal independent commuting subset 
of the gauge operators to stabilizers.
"""
function fix_all_gauges(S::AbstractSubsystemCode; choice::Symbol = :X)
    choice ∈ (:X, :Z) || throw(ArgumentError("Choice must be :X or :Z"))
    
    # Extract properties safely (using hasproperty for custom/lazy subsystem structs)
    F = hasproperty(S, :cache) && haskey(S.cache, :F) ? S.cache[:F] : Oscar.Nemo.Native.GF(2)
    
    n_new = S.n
    k_new = S.k
    
    # Distance bounds: Gauge fixing cannot decrease distance.
    # Therefore, the subsystem l_bound directly becomes the new l_bound.
    l_bound = hasproperty(S, :l_bound) ? S.l_bound : 1
    u_bound = hasproperty(S, :u_bound) ? S.u_bound : S.n
    
    # Because gauge fixing might INCREASE the distance (depending on the choice), 
    # the exact distance d is lost and drops back to a lower bound.
    d_exact = missing 
    if hasproperty(S, :d) && !ismissing(S.d)
        l_bound = max(l_bound, S.d)
    end
    
    cache = Dict{Symbol, Any}(:F => F)
    
    return GaugeFixedCode(
        S, choice,
        n_new, k_new,
        d_exact, l_bound, u_bound,
        cache
    )
end

function stabilizers(S::GaugeFixedCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    sub_S = S.subsystem_code
    pair_idx = S.choice == :X ? 1 : 2
    
    # 1. Get original stabilizers
    base_stabs = stabilizers(sub_S)
    
    # 2. Extract the gauge pairs and select the requested commuting half
    if hasproperty(sub_S, :r) && sub_S.r > 0
        g_pairs = gauge_operators(sub_S)
        fixed_gauges = reduce(vcat, [g[pair_idx] for g in g_pairs])
        new_stabs = vcat(base_stabs, fixed_gauges)
    else
        # If r == 0, it was already a stabilizer code, so we do nothing
        new_stabs = base_stabs
    end
    
    S.cache[:stabilizers] = new_stabs
    return new_stabs
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
            println(" ")

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

# function _all_stabilizers(S::AbstractStabilizerCode, only_print::Bool = false)
#     E = quadraticfield(S)
#     stabs = stabilizers(S)
#     all = Vector{typeof(stabs)}()
    
#     for iter in Base.Iterators.product([0:(Int64(characteristic(S.F)) - 1) for _ in 1:nrows(stabs)]...)
#         stab = E(iter[1]) * stabs[1, :]
#         for r in 2:nrows(stabs)
#             !iszero(iter[r]) && (stab += E(iter[r]) * stabs[r, :])
#         end
#         if only_print
#             println(stab)
#         else
#             push!(all, stab)
#         end
#     end
#     only_print ? return : return all
# end

# all_stabilizers(S::AbstractSubsystemCode) = _all_stabilizers(S, false)
# elements(S::AbstractSubsystemCode) = all_stabilizers(S)
# print_all_stabilizers(S::AbstractSubsystemCode) = _all_stabilizers(S, true)
# print_all_elements(S::AbstractSubsystemCode) = print_all_stabilizers(S)

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
function augment(S::AbstractSubsystemCode, row::CTMatrixTypes; verbose::Bool = true)
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
            iszero(prod[i]) && append!(stabs_to_keep, i, i + 1)
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
    
    temp = _quotient_space(temp, new_stabs, :sys_eqs)
    new_logs = _make_pairs(temp)
    return SubsystemCode(new_stabs, vcat(logs_mat_new, reduce(vcat, [reduce(vcat, new_logs[i]) for i in 1:length(new_logs)])), gauge_ops_new, char_vec = isempty(S.char_vec) ? missing : S.char_vec)
end

"""
    expurgate(S::AbstractStabilizerCode, rows::Vector{Int}; verbose::Bool = true)
"""
function expurgate(S::AbstractSubsystemCode, rows::Vector{Int}; verbose::Bool = true)
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

    new_logs = _quotient_space(H_tr, new_stabs, :sys_eqs)
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
