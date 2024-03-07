# Copyright (c) 2023 - 2024 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    SubsystemCode(G::CTMatrixTypes; char_vec::Union{Vector{nmod}, Missing} = missing, logs_alg::Symbol = :sys_eqs)

Return the subsystem code whose gauge group is determined by `G`.
"""
function SubsystemCode(G::CTMatrixTypes; char_vec::Union{Vector{nmod}, Missing} = missing,
    logs_alg::Symbol = :sys_eqs)

    logs_alg ∈ (:sys_eqs, :VS) || throw(ArgumentError("Unrecognized logicals algorithm"))
    iszero(G) && throw(ArgumentError("The gauge matrix is empty."))
    G = _remove_empty(G, :rows)

    F = base_ring(G)
    p = Int(characteristic(F))
    n = div(ncols(G), 2)
    char_vec = _process_char_vec(char_vec, p, 2 * n)

    # G = < S, gauge ops >
    # C(G) = < S, logs >
    # S = C(G) \cap G
    # C(S) = < S, logs, gauges >
    # logs = C(G) \ S (bare operators)
    # gauge ops = G - S
    # dressed operators = < logs, gauge ops >

    # stabilizer group: ker G ∩ G
    rnk_ker_G, ker_G = right_kernel(hcat(G[:, n + 1:end], -G[:, 1:n]))
    if ncols(ker_G) == rnk_ker_G
        ker_G = transpose(ker_G)
    else
        # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
        nr = nrows(ker_G)
        ker_G_tr = zero_matrix(base_ring(ker_G), rnk_ker_G, nr)
        for r in 1:nr
            for c in 1:rnk_ker_G
                !iszero(ker_G[r, c]) && (ker_G_tr[c, r] = ker_G[r, c];)
            end
        end
        ker_G = ker_G_tr
    end
    V = vector_space(base_ring(G), ncols(ker_G))
    ker_G_VS, ker_G_to_V = sub(V, [V(ker_G[i, :]) for i in 1:nrows(ker_G)])
    G_VS, _ = sub(V, [V(G[i, :]) for i in 1:nrows(G)])
    I, I_to_ker_G = intersect(ker_G_VS, G_VS)
    if !iszero(AbstractAlgebra.dim(I))
        I_basis = [ker_G_to_V(I_to_ker_G(g)) for g in gens(I)]
        F_basis = [[F(I_basis[j][i]) for i in 1:AbstractAlgebra.dim(parent(I_basis[1]))] for j in 1:length(I_basis)]
        stabs = matrix(F, length(F_basis), 2 * n, reduce(vcat, F_basis))
    else
        error("Error computing the stabilizer group of the subsystem code; ker G ∩ G has dimension zero.")
    end

    # check if this subsystem code is really a stabilizer code
    # this should be fine because if they are the same dimension then they are the same
    # don't need to check row space equivalence
    if is_isomorphic(I, G_VS)
        println("Stabilizer code detected.")
        # this is duplicating some of the initial work done in this function,
        # but it seems appropriate not to paste the rest of that constructor here
        return StabilizerCode(stabs, char_vec = char_vec)
    end
    stabs_stand, P_stand, stand_r, stand_k, _ = _standard_form_stabilizer(stabs)

    # bare logicals (reps): ker G / S
    BL = _quotient_space(ker_G, stabs, logs_alg)
    graph_state = false
    if !iszero(BL)
        bare_logs = _make_pairs(BL)
        # verify
        bare_logs_mat = reduce(vcat, [reduce(vcat, bare_logs[i]) for i in 1:length(bare_logs)])
        are_symplectic_orthogonal(stabs, bare_logs_mat) || error("Computed logicals do not commute with the codespace.")
        prod = hcat(bare_logs_mat[:, n + 1:end], -bare_logs_mat[:, 1:n]) * transpose(bare_logs_mat)
        sum(FpmattoJulia(prod), dims=1) == ones(Int, 1, size(prod, 1)) || error("Computed logicals do not have the right commutation relations.")
    else
        graph_state = true
    end
    
    # gauge operators (reps): G / S
    # this can't be zero if S != G
    GO = _quotient_space(G, stabs, logs_alg)
    gauge_ops = _make_pairs(GO)
    # verify
    gauge_ops_mat = reduce(vcat, [reduce(vcat, gauge_ops[i]) for i in 1:length(gauge_ops)])
    are_symplectic_orthogonal(stabs, gauge_ops_mat) || error("Computed gauge operators do not commute with the codespace.")
    if !graph_state
        are_symplectic_orthogonal(bare_logs_mat, gauge_ops_mat) || error("Computed gauge operators do not commute with the computed logicals.")
    end
    prod = hcat(gauge_ops_mat[:, n + 1:end], -gauge_ops_mat[:, 1:n]) * transpose(gauge_ops_mat)
    sum(FpmattoJulia(prod), dims=1) == ones(Int, 1, size(prod, 1)) || error("Computed gauge operators do not have the right commutation relations.")

    # since S is computed, it automatically has full rank and is not overcomplete
    # q^n / p^k but rows is n - k
    top = BigInt(order(F))^n
    r = length(gauge_ops)
    k =  top // BigInt(p)^(nrows(stabs) + r)
    isinteger(k) && (k = round(Int, log(BigInt(p), k));)
    # r = top // BigInt(p)^length(gauge_ops)
    # isinteger(r) && (r = round(Int, log(BigInt(p), r));)
    signs = _determine_signs(stabs, char_vec)

    args = _is_CSS_symplectic(stabs, signs, true)
    if args[1]
        if graph_state
            return GraphStateSubsystemCSS(F, n, 0, r, missing, missing, missing, stabs, args[2],
                args[4], missing, missing, signs, args[3], args[5], char_vec, missing, false,
                gauge_ops, gauge_ops_mat, stabs_stand, stand_r, stand_k, P_stand, missing, missing)
        end
        return SubsystemCodeCSS(F, n, k, r, missing, stabs, args[2], args[4], missing, missing,
            signs, args[3], args[5], bare_logs, bare_logs_mat, char_vec, gauge_ops, gauge_ops_mat,
            false, stabs_stand, stand_r, stand_k, P_stand, missing, missing)
    else
        if graph_state
            return GraphStateSubsystem(F, n, 0, r, missing, stabs, char_vec, signs, missing, false,
                gauge_ops, gauge_ops_mat, stabs_stand, stand_r, stand_k, P_stand, missing)
        end
        return SubsystemCode(F, n, k, r, missing, stabs, signs, bare_logs, bare_logs_mat, char_vec,
            gauge_ops, gauge_ops_mat, false, stabs_stand, stand_r, stand_k, P_stand, missing)
    end
end

"""
    SubsystemCode(G_Pauli::Vector{T}; char_vec::Union{Vector{nmod}, Missing} = missing, logs_alg::Symbol = :sys_eqs) where T <: Union{String, Vector{Char}}

Return the subsystem code whose gauge group is determined by the vector of Pauli strings `G_Pauli`.

# Notes
* Any +/- 1 characters in front of each stabilizer are stripped. No check is done
  to make sure these signs agree with the ones computed using the character vector.
"""
function SubsystemCode(G_Pauli::Vector{T}; char_vec::Union{Vector{nmod}, Missing} = missing,
    logs_alg::Symbol = :sys_eqs) where T <: Union{String, Vector{Char}}

    logs_alg ∈ (:sys_eqs, :VS) || throw(ArgumentError("Unrecognized logicals algorithm"))
    G_Pauli_stripped = _process_strings(G_Pauli)
    G = _Pauli_string_to_symplectic(G_Pauli_stripped)
    iszero(G) && error("The processed Pauli strings returned a set of empty gauge group generators.")
    return SubsystemCode(G, char_vec = char_vec, logs_alg = logs_alg)
end

"""
    SubsystemCode(S::CTMatrixTypes, L::CTMatrixTypes, G::CTMatrixTypes; char_vec::Union{Vector{nmod}, Missing} = missing)

Return the subsystem code whose stabilizers are given by `S`, (bare) logical operators
by `L`, gauge operators (not including stabilizers) by `G`.
"""
function SubsystemCode(S::CTMatrixTypes, L::CTMatrixTypes, G::CTMatrixTypes;
    char_vec::Union{Vector{nmod}, Missing} = missing)

    iszero(S) && error("The stabilizer matrix is empty.")
    S = _remove_empty(S, :rows)
    n = div(ncols(S), 2)
    # check S commutes with itself
    are_symplectic_orthogonal(S, S) || error("The given stabilizers are not symplectic orthogonal.")
    stabs_stand, P_stand, stand_r, stand_k, rnk = _standard_form_stabilizer(S)

    # logicals
    iszero(L) && error("The logicals are empty.")
    L = _remove_empty(L, :rows)
    # check S commutes with L
    are_symplectic_orthogonal(S, L) || error("Logicals do not commute with the code.")
    # check L doesn't commute with itself
    prod = hcat(L[:, n + 1:end], -L[:, 1:n]) * transpose(L)
    iszero(prod) && error("Logicals should not be symplectic self-orthogonal.")
    # check L is properly in pairs
    nc_pr = ncols(prod)
    # need an integer sum and this is cheaper than Nemo.ZZ
    prod_Jul = FpmattoJulia(prod)
    cols = [sum(prod_Jul[:, i]) for i in 1:nc_pr]
    sum(cols) == nc_pr || println("Detected logicals not in anticommuting pairs.")
    log_pairs = _make_pairs(L)
    logs_mat = reduce(vcat, [reduce(vcat, log_pairs[i]) for i in 1:length(log_pairs)])

    # gauge operators
    iszero(G) && error("The gauges are empty.")
    G = _remove_empty(G, :rows)
    # check S commutes with G
    are_symplectic_orthogonal(S, G) || error("Gauges do not commute with the code.")
    # check L commutes with G
    are_symplectic_orthogonal(logs_mat, G) || error("Gauges do not commute with the logicals.")
    # check G doesn't commute with itself
    prod = hcat(G[:, n + 1:end], -G[:, 1:n]) * transpose(G)
    iszero(prod) && error("Gauges should not be symplectic self-orthogonal.")
    # check G is properly in pairs
    nc_pr = ncols(prod)
    # need an integer sum and this is cheaper than Nemo.ZZ
    prod_Jul = FpmattoJulia(prod)
    cols = [sum(prod_Jul[:, i]) for i in 1:nc_pr]
    sum(cols) == nc_pr || println("Detected gauges not in anticommuting pairs.")
    # display(prod)
    g_ops_pairs = _make_pairs(G)
    g_ops_mat = reduce(vcat, [reduce(vcat, g_ops_pairs[i]) for i in 1:length(g_ops_pairs)])

    F = base_ring(g_ops_mat)
    p = Int(characteristic(F))
    n = div(ncols(G), 2)
    char_vec = _process_char_vec(char_vec, p, 2 * n)
    nrows(S) > rnk ? (overcomp = true) : (overcomp = false)
    rnk_G = rank(G)

    top = BigInt(order(F))^n
    # TODO: handle overcompleteness in inputs
    r = div(rnk_G, 2)
    k =  top // BigInt(p)^(rnk + r)
    isinteger(k) && (k = round(Int, log(BigInt(p), k));)
    signs = _determine_signs(S, char_vec)

    args = _is_CSS_symplectic(S, signs, true)
    if args[1]
        return SubsystemCodeCSS(F, n, k, r, missing, S, args[2], args[4], missing, missing, signs,
            args[3], args[5], log_pairs, logs_mat, char_vec, g_ops_pairs, g_ops_mat, false,
            stabs_stand, stand_r, stand_k, P_stand, missing, missing)
    else
        return SubsystemCode(F, n, k, r, missing, S, signs, log_pairs, logs_mat, char_vec,
            g_ops_pairs, g_ops_mat, false, stabs_stand, stand_r, stand_k, P_stand, missing)
    end
end

# if people want to make a graph code go through the other constructor
"""
    SubsystemCode(S_Pauli::Vector{T}, L_Pauli::Vector{T}, G_Pauli::Vector{T}; char_vec::Union{Vector{nmod}, Missing} = missing) where T <: Union{String, Vector{Char}}

Return the subsystem code whose stabilizers are given by the vectors of Pauli strings `S_Pauli`,
(bare) logical operators by `L_Pauli`, gauge operators (not including stabilizers) by `G_Pauli`.    
"""
function SubsystemCode(S_Pauli::Vector{T}, L_Pauli::Vector{T}, G_Pauli::Vector{T};
    char_vec::Union{Vector{nmod}, Missing} = missing) where T <: Union{String, Vector{Char}}

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
cardinality(S::AbstractSubsystemCode) = BigInt(characteristic(S.F))^(S.n - S.k)

"""
    rate(S::AbstractSubsystemCode)

Return the rate, `R = k/n`, of the code.
"""
rate(S::AbstractSubsystemCode) = S.k / S.n

"""
    signs(S::AbstractSubsystemCode)

Return the signs of the stabilizers of the code.
"""
signs(S::AbstractSubsystemCode) = S.signs

"""
    X_signs(S::AbstractSubsystemCode)

Return the signs of the `X` stabilizers of the CSS code.
"""
X_signs(S::T) where {T <: AbstractSubsystemCode} = X_signs(CSSTrait(T), S)
X_signs(::IsCSS, S::AbstractSubsystemCode) = S.X_signs
X_signs(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    Z_signs(S::AbstractSubsystemCode)

Return the signs of the `Z` stabilizers of the CSS code.
"""
Z_signs(S::T) where {T <: AbstractSubsystemCode} = Z_signs(CSSTrait(T), S)
Z_signs(::IsCSS, S::AbstractSubsystemCode) = S.Z_signs
Z_signs(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    stabilizers(S::AbstractSubsystemCode; standform::Bool = false)

Return the stabilizer matrix of the code.

# Notes
- If the optional parameter `standform` is set to `true`, the standard form of the
  stabilizer matrix is returned instead.
"""
stabilizers(S::AbstractSubsystemCode; standform::Bool = false) = standform ?
    (return S.stabs_stand) : (return S.stabs)

"""
    standard_form_permutation(S::AbstractSubsystemCode)

Return the permutation matrix required to permute the columns of the code matrices to have the same
row space as the matrices in standard form. Returns `missing` is no such permutation is required.
"""
standard_form_permutation(S::AbstractSubsystemCode) = S.P_stand

"""
    standard_form_A(S::AbstractSubsystemCode)

Return the named matrix `A` from the standard form of the stabilizer matrix.
"""
standard_form_A(S::AbstractSubsystemCode) = S.stabs_stand[1:S.stand_r, S.stand_r + 1:S.n]

"""
    standard_form_A1(S::AbstractSubsystemCode)
    
Return the named matrix `A1` from the standard form of the stabilizer matrix.
"""
standard_form_A1(S::AbstractSubsystemCode) = S.stabs_stand[1:S.stand_r, S.stand_r + 1:S.n -
    S.stand_k]

"""
    standard_form_A2(S::AbstractSubsystemCode)
    
Return the named matrix `A2` from the standard form of the stabilizer matrix.
"""
standard_form_A2(S::AbstractSubsystemCode) = S.stabs_stand[1:S.stand_r, S.n - S.stand_k + 1:S.n]

"""
    standard_form_B(S::AbstractSubsystemCode)
    
Return the named matrix `B` from the standard form of the stabilizer matrix.
"""
standard_form_B(S::AbstractSubsystemCode) = S.stabs_stand[1:S.stand_r, S.n + 1:S.n + S.stand_r]

"""
    standard_form_C1(S::AbstractSubsystemCode)
    
Return the named matrix `C1` from the standard form of the stabilizer matrix.
"""
standard_form_C1(S::AbstractSubsystemCode) = S.stabs_stand[1:S.stand_r, S.n + S.stand_r + 1:2 *
    S.n - S.stand_k]

"""
    standard_form_C2(S::AbstractSubsystemCode)
    
Return the named matrix `C2` from the standard form of the stabilizer matrix.
"""
standard_form_C2(S::AbstractSubsystemCode) = S.stabs_stand[1:S.stand_r, 2 * S.n - S.stand_k + 1:2 *
    S.n]

"""
    standard_form_D(S::AbstractSubsystemCode)
    
Return the named matrix `D` from the standard form of the stabilizer matrix.
"""
standard_form_D(S::AbstractSubsystemCode) = S.stabs_stand[S.stand_r + 1:S.n - S.stand_k, S.n +
    1:S.n + S.stand_r]

"""
    standard_form_E(S::AbstractSubsystemCode)
    
Return the named matrix `E` from the standard form of the stabilizer matrix.
"""
standard_form_E(S::AbstractSubsystemCode) = S.stabs_stand[S.stand_r + 1:S.n - S.stand_k, 2 * S.n -
    S.stand_k + 1:2 * S.n]

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
    metacheck(S::AbstractSubsystemCode)

Return the metacheck matrix of the code, if it has been set; otherwise returns missing.
"""
metacheck(S::T) where {T <: AbstractSubsystemCode} = metacheck(CSSTrait(T), S)
metacheck(::IsCSS, S::AbstractSubsystemCode) =
    error("Use `X_metacheck` or `Z_metacheck` for CSS codes.")
metacheck(::IsNotCSS, S::AbstractSubsystemCode) = S.metacheck

"""
    X_metacheck(S::AbstractSubsystemCode)

Return the `X`-metacheck matrix of the CSS code, if it has been set; otherwise returns missing.
"""
X_metacheck(S::T) where {T <: AbstractSubsystemCode} = X_metacheck(CSSTrait(T), S)
X_metacheck(::IsCSS, S::AbstractSubsystemCode) = S.X_metacheck
X_metacheck(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    Z_metacheck(S::AbstractSubsystemCode)

Return the `Z`-metacheck matrix of the CSS code, if it has been set; otherwise returns missing.
"""
Z_metacheck(S::T) where {T <: AbstractSubsystemCode} = Z_metacheck(CSSTrait(T), S)
Z_metacheck(::IsCSS, S::AbstractSubsystemCode) = S.Z_metacheck
Z_metacheck(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

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

"""
    is_overcomplete(S::AbstractSubsystemCode)

Return `true` if `S` has an overcomplete set of stabilizers.
"""
is_overcomplete(S::AbstractSubsystemCode) = S.overcomplete

# TODO: redo here
"""
    is_CSS(S::AbstractSubsystemCode)

Return `true` if `S` is CSS.
"""
is_CSS(S::T) where {T <: AbstractSubsystemCode} = is_CSS(CSSTrait(T), S)
is_CSS(::IsCSS, S::AbstractSubsystemCode) = true
is_CSS(::IsNotCSS, S::AbstractSubsystemCode) = false

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
    logicals(S::AbstractSubsystemCode)
    logical_operators(S::AbstractSubsystemCode)
    bare_logicals(S::AbstractSubsystemCode)
    bare(S::AbstractSubsystemCode)

Return a vector of logical operator generator pairs for `S`.
"""
logicals(S::T) where {T <: AbstractSubsystemCode} = logicals(LogicalTrait(T), S)
logicals(::HasLogicals, S::AbstractSubsystemCode) = S.logicals
logicals(::HasNoLogicals, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no logicals.")
logical_operators(S::AbstractSubsystemCode) = logicals(S)
bare_logicals(S::AbstractSubsystemCode) = logicals(S)
bare(S::AbstractSubsystemCode) = logicals(S)

"""
    logicals_matrix(S::AbstractSubsystemCode)

Returns the result of `logicals(S)` as a vertically concatenated matrix.
"""
logicals_matrix(S::T) where {T <: AbstractSubsystemCode} = logicals_matrix(LogicalTrait(T), S)
logicals_matrix(::HasLogicals, S::AbstractSubsystemCode) = S.logs_mat
logicals_matrix(::HasNoLogicals, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no logicals.")

"""
    logicals_standard_form(S::AbstractSubsystemCode)

Return the a matrix of logical operators as determined by the stabilizers in standard form.
"""
logicals_standard_form(S::T) where {T <: AbstractSubsystemCode} = logicals_standard_form(LogicalsTrait(T), S)
logicals_standard_form(::HasLogicals, S::AbstractSubsystemCode) = _logicals_standard_form(S.stabs_stand, S.n, S.stand_k,S.stand_r, S.P_stand)
logicals_standard_form(::HasNoLogicals, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no logicals.")

"""
    gauges(S::AbstractSubsystemCode)
    gauge_operators(S::AbstractSubsystemCode

Return a vector of gauge operator generator pairs for `S`.

# Notes
- Here, gauge operators refers to the gauge group minus the stabilizers.
"""
gauges(S::T) where {T <: AbstractSubsystemCode} = gauges(GaugeTrait(T), S)
gauges(::HasGauges, S::AbstractSubsystemCode) = S.gauge_ops
gauges(::HasNoGauges, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no gauges.")
gauge_operators(S::AbstractSubsystemCode) = gauges(S)

"""
    gauges_matrix(S::AbstractSubsystemCode)
    gauge_operators_matrix(S::AbstractSubsystemCode)

Return the result of `gauges(S)` as a vertically concatenated matrix.
"""
gauges_matrix(S::T) where {T <: AbstractSubsystemCode} = gauges_matrix(GaugeTrait(T), S)
gauges_matrix(::HasGauges, S::AbstractSubsystemCode) = S.g_ops_mat
gauges_matrix(::HasNoGauges, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no gauges.")
gauge_operators_matrix(S::AbstractSubsystemCode) = gauges_matrix(S)

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
    elseif GaugeTrait(S) == HasNoGauges
        error("Type $T has no gauges.")
    end
    return S.logicals ∪ S.gauge_ops
end
dressed_operators(S::AbstractSubsystemCode) = dressed(S)
dressed_logicals(S::AbstractSubsystemCode) = dressed(S)

"""
    gauge_group(S::AbstractSubsystemCode)
    gauge_group_matrix(S::AbstractSubsystemCode)
    gauge_generators_matrix(S::AbstractSubsystemCode)
    gauge_group_generators_matrix(S::AbstractSubsystemCode)

Return a matrix giving a (maybe overcomplete) basis for the gauge group.

# Notes
* Here, this is the stabilizers and the gauge operators.
"""
gauge_group(S::T) where {T <: AbstractSubsystemCode} = gauge_group(GaugeTrait(T), S)
gauge_group(::HasGauges, S::AbstractSubsystemCode) = vcat(S.stabs, S.g_ops_mat)
gauge_group(::HasNoGauges, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no gauges.")
gauge_group_matrix(S::AbstractSubsystemCode) = gauge_group(S)
gauge_generators_matrix(S::AbstractSubsystemCode) = gauge_group(S)
gauge_group_generators_matrix(S::AbstractSubsystemCode) = gauge_group(S)

#############################
      # setter functions
#############################

"""
    set_signs(S::AbstractSubsystemCode, char_vec::Vector{nmod})
    set_signs!(S::AbstractSubsystemCode, char_vec::Vector{nmod})

Set the character vector of `S` to `char_vec` and update the signs.
"""
function set_signs!(S::AbstractSubsystemCode, char_vec::Vector{nmod})
    R = base_ring(character_vector(S))
    length(char_vec) == 2 * S.n || throw(ArgumentError("Characteristic vector is of improper length for the code."))
    for s in char_vec
        modulus(s) == modulus(R) || throw(ArgumentError("Phases are not in the correct ring."))
    end

    S.signs = _get_signs(S.stabilizers, char_vec)
    S.char_vec = char_vec
end
set_signs(S::AbstractSubsystemCode, char_vec::Vector{nmod}) = (S_new = deepcopy(S);
    return setsign!s(S_new, char_vec))

"""
    set_stabilizers(S::AbstractSubsystemCode, stabs::CTMatrixTypes)
    set_stabilizers!(S::AbstractSubsystemCode, stabs::CTMatrixTypes)

Set the stabilizers of `S` to `stabs`.

# Notes
- A check is done to make sure `stabs` are equivalent to the current set of stabilizers.
"""
function set_stabilizers!(S::AbstractSubsystemCode, stabs::CTMatrixTypes)
    iszero(stabs) && throw(ArgumentError("The stabilizers cannot be zero."))
    order(S.F) == order(base_ring(stabs)) || throw(ArgumentError("The stabilizers must be over the same field as the code."))

    stabs = _remove_empty(stabs, :rows)
    stabs = change_base_ring(S.F, stabs)
    if _has_equivalent_row_spaces(S.stabs, stabs)
        S.stabs = stabs
        nrows(stabs) != S.k && (S.overcomplete = true;)
    else
        error("The current stabilizers are not equivalent to the input.")
    end
    S.signs = _determine_signs(stabs, S.char_vec)

    if CSSTrait(typeof(S)) == IsCSS()
        flag, X_stabs, X_signs, Z_stabs, Z_signs = _is_CSS_symplectic(stabs, S.signs, true)
        flag || error("Detected equivalent stabilizers but is no longer CSS.")
        S.X_stabs = X_stabs
        S.X_signs = X_signs
        S.Z_stabs = Z_stabs
        S.Z_signs = Z_signs
    end
    return nothing
end
set_stabilizers(S::AbstractSubsystemCode, stabs::CTMatrixTypes) = (S_new = deepcopy(S);
    set_stabilizers!(S_new, stabs); return S_new)

"""
    set_X_stabilizers(S::AbstractSubsystemCode, X_stabs::CTMatrixTypes; trimmed::Bool = true)
    set_X_stabilizers!(S::AbstractSubsystemCode, X_stabs::CTMatrixTypes; trimmed::Bool = true)

Set the `X` stabilizers of `S` to `X_stabs`. If `trimmed` is `true`, `X_stabs` are assumed
to have `length(S)` columns; otherwise, they are assumed to be in symplectic form with
`2 * length(S)` columns.

# Notes
* A check is done to make sure `stabs` are equivalent to the current set of stabilizers.
"""
set_X_stabilizers!(S::T, X_stabs::CTMatrixTypes; trimmed::Bool = true) where {T <:
    AbstractSubsystemCode} = set_X_stabilizers!(CSSTrait(T), S, X_stabs, trimmed)
function set_X_stabilizers!(::IsCSS, S::AbstractSubsystemCode, X_stabs::CTMatrixTypes,
    trimmed::Bool)

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
        nrows(X_trimmed) != rank(X_trimmed) && (S.overcomplete = true;)
    else
        error("The current stabilizers are not equivalent to the input.")
    end
    S.X_signs = _determine_signs(X_trimmed, S.char_vec)

    S.stabs = vcat(hcat(X_trimmed, zero_matrix(S.F, nrows(X_trimmed), S.n)),
        hcat(zero_matrix(S.F, nrows(S.Z_stabs), S.n), S.Z_stabs))
    # stabs_stand, P_stand, stand_r, stand_k, _ = _standard_form_stabilizer(S.stabs)
    # S.stabs_stand = stabs_stand
    # S.stand_r = stand_r
    # S.stand_k = stand_k
    # S.P_stand = P_stand
    S.signs = S.X_signs ∪ S.Z_signs
    return nothing
end
set_X_stabilizers!(::IsNotCSS, S::AbstractSubsystemCode, X_stabs::CTMatrixTypes, trimmed::Bool) =
    error("X stabilizers are only defined for CSS codes")
set_X_stabilizers(S::T, X_stabs::CTMatrixTypes; trimmed::Bool = true) where {T <:
    AbstractSubsystemCode} = set_X_stabilizers(CSSTrait(T), S, X_stabs, trimmed)
set_X_stabilizers(::IsCSS, S::AbstractSubsystemCode, X_stabs::CTMatrixTypes, trimmed::Bool) =
    (S_new = deepcopy(S); set_X_stabilizers!(IsCSS(), S_new, X_stabs,  trimmed); return S_new)
set_X_stabilizers(::IsNotCSS, S::AbstractSubsystemCode, X_stabs::CTMatrixTypes, trimmed::Bool) =
    error("X stabilizers are only defined for CSS codes")

"""
    set_Z_stabilizers(S::AbstractSubsystemCode, Z_stabs::CTMatrixTypes; trimmed::Bool = true)
    set_Z_stabilizers!(S::AbstractSubsystemCode, Z_stabs::CTMatrixTypes; trimmed::Bool = true)

Set the `Z` stabilizers of `S` to `Z_stabs`. If `trimmed` is `true`, `Z_stabs` are assumed
to have `length(S)` columns; otherwise, they are assumed to be in symplectic form with
`2 * length(S)` columns.

# Notes
* A check is done to make sure `stabs` are equivalent to the current set of stabilizers.
"""
set_Z_stabilizers!(S::T, Z_stabs::CTMatrixTypes; trimmed::Bool = true) where {T <:
    AbstractSubsystemCode} = set_Z_stabilizers!(CSSTrait(T), S, Z_stabs, trimmed)
function set_Z_stabilizers!(::IsCSS, S::AbstractSubsystemCode, Z_stabs::CTMatrixTypes,
    trimmed::Bool)

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
        nrows(Z_trimmed) != rank(Z_trimmed) && (S.overcomplete = true;)
    else
        error("The current stabilizers are not equivalent to the input.")
    end
    S.Z_signs = _determine_signs(Z_trimmed, S.char_vec)

    S.stabs = vcat(hcat(S.X_stabs, zero_matrix(S.F, nrows(S.X_stabs), S.n)),
        hcat(zero_matrix(S.F, nrows(Z_trimmed), S.n), Z_trimmed))
    S.signs = S.X_signs ∪ S.Z_signs
    return nothing
end
set_Z_stabilizers!(::IsNotCSS, S::AbstractSubsystemCode, Z_stabs::CTMatrixTypes, trimmed::Bool) =
    error("Z stabilizers are only defined for CSS codes")
set_Z_stabilizers(S::T, Z_stabs::CTMatrixTypes; trimmed::Bool = true) where {T <:
    AbstractSubsystemCode} = set_Z_stabilizers(CSSTrait(T), S, Z_stabs, trimmed)
set_Z_stabilizers(::IsCSS, S::AbstractSubsystemCode, Z_stabs::CTMatrixTypes, trimmed::Bool) = (S_new = deepcopy(S); set_Z_stabilizers!(IsCSS(), S_new, Z_stabs, trimmed); return S_new)
set_Z_stabilizers(::IsNotCSS, S::AbstractSubsystemCode, Z_stabs::CTMatrixTypes, trimmed::Bool) =
    error("Z stabilizers are only defined for CSS codes")

"""
    set_logicals(S::AbstractSubsystemCode, L::CTMatrixTypes)
    set_logicals!(S::AbstractSubsystemCode, L::CTMatrixTypes)

Set the logical operators of `S` to `L`.

# Notes
* A check is done to make sure `L` are eqivalent to the current set of logicals (up to stabilizers).
"""
set_logicals!(S::T, L::W) where {T <: AbstractSubsystemCode, W <: CTMatrixTypes} =
    set_logicals!(LogicalTrait(T), S, L)
function set_logicals!(::HasLogicals, S::AbstractSubsystemCode, L::W) where {W <: CTMatrixTypes}
    size(L) == (2 * S.k, 2 * S.n) || throw(ArgumentError("Provided matrix is of incorrect size for the logical space."))
    iseven(ncols(L)) || throw(ArgumentError("Expected a symplectic input but the input matrix has an odd number of columns."))
    S.F == base_ring(L) || throw(ArgumentError("The logicals must be over the same field as the code."))
    
    _has_equivalent_row_spaces(vcat(S.logs_mat, S.stabs), vcat(L, S.stabs)) || error("The current logicals are not equivalent to the input.")
    # are_symplectic_orthogonal(symplecticstabilizers(S), Lsym, true) ||
    #     error("Provided logicals do not commute with the code.")

    # the columns in prod give the commutation relationships between the provided
    # logical operators; they ideally should only consist of {X_1, Z_i} pairs
    # so there should only be one nonzero element in each column
    prod = hcat(L[:, S.n + 1:end], -L[:, 1:S.n]) * transpose(L)
    iszero(prod) && throw(ArgumentError("Provided logicals should not be symplectic self-orthogonal."))
    nc_pr = ncols(prod)
    # need an integer sum and this is cheaper than Nemo.ZZ
    prod_Jul = FpmattoJulia(prod)
    cols = [sum(prod_Jul[:, i]) for i in 1:nc_pr]
    sum(cols) == nc_pr || throw(ArgumentError("Incorrect commutation relationships between provided logicals."))

    # pairs are row i, and then whatever column is nonzero, and then shift such that it is one
    F = base_ring(L)
    F_one = F(1)
    logs = Vector{Tuple{W, W}}()
    if Int(order(F)) != 2
        # this does indeed grow smaller each iteration
        while nrows(L) >= 2
            y = findfirst(x -> x > 0, prod_Jul[:, 1])
            y = [F(prod[y, 1]), y]
            if y[1] != F_one
                push!(logs, (L[1, :], y[1]^-1 * L[y[2], :]))
            else
                push!(logs, (L[1, :], L[y[2], :]))
            end
            L = L[setdiff(1:size(L, 1), [1, y[2]]), :]
        end
    else
        # this does indeed grow smaller each iteration
        while nrows(L) >= 2
            y = findfirst(x -> x > 0, prod_Jul[:, 1])
            y = [F(prod[y, 1]), y]
            push!(logs, (L[1, :], L[y[2], :]))
            L = L[setdiff(1:size(L, 1), [1, y[2]]), :]
        end
    end
    S.logicals = logs
    S.logs_mat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
end
set_logicals!(::HasNoLogicals, S::AbstractSubsystemCode, L::CTMatrixTypes) =
    error("Type $(typeof(S)) has no logicals.")
set_logicals(S::T, L::CTMatrixTypes) where {T <: AbstractSubsystemCode} =
    set_logicals(LogicalTrait(T), S, L)
set_logicals(::HasLogicals, S::AbstractSubsystemCode, L::CTMatrixTypes) =
    (S_new = deepcopy(S); set_logicals!(S_new, L); return S_new)
set_logicals(::HasNoLogicals, S::AbstractSubsystemCode, L::CTMatrixTypes) =
    error("Type $(typeof(S)) has no logicals.")

"""
    set_metacheck(S::AbstractSubsystemCode, M::CTMatrixTypes)
    set_metacheck!(S::AbstractSubsystemCode, M::CTMatrixTypes)

Set the metacheck matrix of `S` to `M`.
"""
set_metacheck!(S::T, M::U) where {T <: AbstractSubsystemCode, U <: CTMatrixTypes} =
    set_metacheck!(CSSTrait(T), S, M)
set_metacheck!(::IsCSS, S::AbstractSubsystemCode, M::CTMatrixTypes) =
    error("Use `set_X_metacheck` and `set_Z_metacheck` for CSS codes.")
function set_metacheck!(::IsNotCSS, S::AbstractSubsystemCode, M::CTMatrixTypes)
    iszero(M * S.stabs) ? (S.metacheck = M;) : error("Invalid metacheck for code")
    return nothing
end
set_metacheck(S::T, M::U) where {T <: AbstractSubsystemCode, U <: CTMatrixTypes} =
    set_metacheck(CSSTrait(T), S, M)
set_metacheck(::IsCSS, S::AbstractSubsystemCode, M::CTMatrixTypes) =
    error("Use `set_X_metacheck` and `set_Z_metacheck` for CSS codes.")
set_metacheck(::IsNotCSS, S::AbstractSubsystemCode, M::CTMatrixTypes) =
    (S_new = deepcopy(S); set_metacheck!(IsNotCSS(), S_new, M); return S_new)

"""
    set_X_metacheck(S::AbstractSubsystemCode, M::CTMatrixTypes)
    set_X_metacheck!(S::AbstractSubsystemCode, M::CTMatrixTypes)

Set the `X`-metacheck matrix of `S` to `M`.
"""
set_X_metacheck!(S::T, M::U) where {T <: AbstractSubsystemCode, U <: CTMatrixTypes} =
    set_X_metacheck!(CSSTrait(T), S, M)
function set_X_metacheck!(::IsCSS, S::AbstractSubsystemCode, M::CTMatrixTypes)
    iszero(M * S.X_stabs) ? (S.X_metacheck = M;) : error("Invalid metacheck for code")
    return nothing
end
set_X_metacheck!(::IsNotCSS, S::AbstractSubsystemCode, M::CTMatrixTypes) = 
    error("Only valid for CSS codes.")
set_X_metacheck(S::T, M::U) where {T <: AbstractSubsystemCode, U <: CTMatrixTypes} =
    set_X_metacheck(CSSTrait(T), S, M)
set_X_metacheck(::IsCSS, S::AbstractSubsystemCode, M::CTMatrixTypes) =
    (S_new = deepcopy(S); set_X_metacheck!(IsNotCSS(), S_new, M); return S_new)
set_X_metacheck(::IsNotCSS, S::AbstractSubsystemCode, M::CTMatrixTypes) =
    error("Only valid for CSS codes.")

"""
    set_Z_metacheck(S::AbstractSubsystemCode, M::CTMatrixTypes)
    set_Z_metacheck!(S::AbstractSubsystemCode, M::CTMatrixTypes)

Set the `Z`-metacheck matrix of `S` to `M`.
"""
set_Z_metacheck!(S::T, M::U) where {T <: AbstractSubsystemCode, U <: CTMatrixTypes} =
    set_Z_metacheck!(CSSTrait(T), S, M)
function set_Z_metacheck!(::IsCSS, S::AbstractSubsystemCode, M::CTMatrixTypes)
    iszero(M * S.Z_stabs) ? (S.Z_metacheck = M;) : error("Invalid metacheck for code")
    return nothing
end
set_Z_metacheck!(::IsNotCSS, S::AbstractSubsystemCode, M::CTMatrixTypes) = 
    error("Only valid for CSS codes.")
set_Z_metacheck(S::T, M::U) where {T <: AbstractSubsystemCode, U <: CTMatrixTypes} =
    set_Z_metacheck(CSSTrait(T), S, M)
set_Z_metacheck(::IsCSS, S::AbstractSubsystemCode, M::CTMatrixTypes) =
    (S_new = deepcopy(S); set_Z_metacheck!(IsNotCSS(), S_new, M); return S_new)
set_Z_metacheck(::IsNotCSS, S::AbstractSubsystemCode, M::CTMatrixTypes) =
    error("Only valid for CSS codes.")

"""
    set_minimum_distance!(S::AbstractSubsystemCode, d::Int)
Set the minimum distance of the code to `d`.

# Notes
* The only check done on the value of `d` is that `1 ≤ d ≤ n`.
"""
function set_minimum_distance!(S::AbstractSubsystemCode, d::Int)
    # TODO: should check bounds like Singleton for possibilities
    d > 0 && d <= S.n || throw(DomainError("The minimum distance of a code must be ≥ 1; received: d = $d."))
    S.d = d
end

#############################
     # general functions
#############################

function _process_char_vec(char_vec::Union{Vector{nmod}, Missing}, p::Int, n::Int)
    if !ismissing(char_vec)
        n == length(char_vec) || throw(ArgumentError("The characteristic value is of incorrect length."))
        if p == 2
            R = residue_ring(Nemo.ZZ, 4)
        else
            R = residue_ring(Nemo.ZZ, p)
        end
        for s in char_vec
            modulus(s) == modulus(R) || throw(ArgumentError("Phases are not in the correct ring."))
        end
    else
        if p == 2
            R = residue_ring(Nemo.ZZ, 4)
        else
            R = residue_ring(Nemo.ZZ, p)
        end
        char_vec = [R(0) for _ in 1:n]
    end
    return char_vec
end

function _determine_signs(S::CTMatrixTypes, char_vec::Vector{nmod})
    if iszero(char_vec)
        R = parent(char_vec[1])
        signs = [R(0) for _ in 1:nrows(S)]
    else
        signs = _get_signs(S, char_vec)
    end
    return signs
end

function _determine_signs_CSS(S::CTMatrixTypes, char_vec::Vector{nmod}, X_size::Int, Z_size::Int)
    if iszero(char_vec)
        R = parent(char_vec[1])
        signs = [R(0) for _ in 1:nrows(S)]
        X_signs = [R(0) for _ in 1:X_size]
        Z_signs = [R(0) for _ in 1:Z_size]
    else
        signs = _get_signs(S, char_vec)
        X_signs = signs[1:X_size, :]
        Z_signs = signs[X_size + 1:end, :]
    end
    return signs, X_signs, Z_signs
end

function _get_signs(A::CTMatrixTypes, char_vec::Vector{nmod})
    R = base_ring(char_vec[1])
    nc = ncols(A)
    length(char_vec) == nc || throw(ArgumentError("Input to _get_signs is expected to be in symplectic form and of the same length as the characteristic vector."))
    
    iszero(char_vec) && return [R(0) for _ in 1:div(nc, 2)]
    signs = Vector{Int64}()
    for r in 1:nrows(A)
        parity = R(0)
        for c = 1:nc
            !iszero(A[r, c]) && (parity += char_vec[c];)
        end
        append!(signs, parity)
    end
    return signs
end

function _split_vectors_CSS(S::T, signs::Vector{nmod}) where {T <: CTMatrixTypes}
    X_stabs = Vector{T}()
    X_signs = Vector{nmod}()
    Z_stabs = Vector{T}()
    Z_signs = Vector{nmod}()
    mixed_stabs = Vector{T}()
    mixed_signs = Vector{nmod}()

    half = div(ncols(S), 2)
    for r in 1:nrows(S)
        # use views?
        s = S[r, :]
        if iszero(s)
            continue
        else
            s_x = iszero(s[1, 1:half])
            s_z = iszero(s[1, half + 1:end])
            if (s_x && !s_z)
                push!(Z_stabs, s)
                push!(Z_signs, signs[r])
            elseif !s_x && s_z
                push!(X_stabs, s)
                push!(X_signs, signs[r])
            elseif !s_x && !s_z
                push!(mixed_stabs, s)
                push!(mixed_signs, signs[r])
            end
        end
    end

    if !isempty(X_stabs)
        X_stabs = reduce(vcat, X_stabs)
    end
    if !isempty(Z_stabs)
        Z_stabs = reduce(vcat, Z_stabs)
    end
    if !isempty(mixed_stabs)
        mixed_stabs = reduce(vcat, mixed_stabs)
    end
    return X_stabs, X_signs, Z_stabs, Z_signs, mixed_stabs, mixed_signs
end

# """
#     splitstabilizers(S::AbstractSubsystemCode)

# Return the set of `X`-only stabilizers and their signs, the set of `Z`-only
# stabilizers and their signs, and the remaining stabilizers and their signs.

# # Notes
# * This function returns six objects of alternating types `CTMatrixTypes` and
# `Vector{Int}` for the three sets of stabilizers and signs, respectively.
# An empty set of stabilizers is returned as type `Vector{CTMatrixTypes}`.
# """
# splitstabilizers(S::AbstractSubsystemCode) = _split_vectors_CSS(S.stabs, S.signs)

# TODO: rethink how I'm returning all of this and the bottom trim stuff
# probably need to simply redo the above to simply start with zero matrix
# and then either set first or add to it (return matrix[2:end, L])
# TODO: need more robust CSS detection, what if I add and X and Z stabilizer and use it implace of the Z
# TODO: what's going on here
function _is_CSS_symplectic(stabs::T, signs::Vector{nmod}, trim::Bool=true) where T <: CTMatrixTypes
    X_stabs, X_signs, Z_stabs, Z_signs, mixed_stabs, mixed_signs = _split_vectors_CSS(stabs, signs)
    if typeof(mixed_stabs) <: Vector{T}
        if trim
            half = div(ncols(stabs), 2)
            return true, X_stabs[:, 1:half], X_signs, Z_stabs[:, half + 1:end], Z_signs
        else
            return true, X_stabs, X_signs, Z_stabs, Z_signs
        end
    else
        if trim
            half = div(ncols(stabs), 2)
            !(typeof(X_stabs) <: Vector{T}) && (X_stabs = X_stabs[:, 1:half];)
            !(typeof(Z_stabs) <: Vector{T}) && (Z_stabs = Z_stabs[:, half + 1:end];)
            return false, X_stabs, X_signs, Z_stabs, Z_signs, mixed_stabs, mixed_signs
        else
            return false, X_stabs, X_signs, Z_stabs, Z_signs, mixed_stabs, mixed_signs
        end
    end
end

# using this function for logical and gauge operators
function _make_pairs(L::T) where T <: CTMatrixTypes
    F = base_ring(L)
    n = div(ncols(L), 2)
    logs = Vector{Tuple{T, T}}()
    if Int(order(F)) != 2
        # this does indeed grow smaller each iteration
        while nrows(L) >= 2
            # the columns in prod give the commutation relationships between the provided
            # logical operators; they ideally should only consist of {X_1, Z_i} pairs
            # so there should only be one nonzero element in each column
            prod = hcat(L[:, n + 1:end], -L[:, 1:n]) * transpose(L)
            # println("before")
            # display(prod)
            num_prod = ncols(prod)
            first = 0
            for c in 1:num_prod
                if !iszero(prod[1, c])
                    if iszero(first)
                        first = c
                        if !isone(prod[1, c])
                            L[first, :] *= F(prod[1, c]^-1)
                        end
                    else
                        L[c, :] += F(prod[1, c]^-1) * L[first, :]
                    end
                end
            end
            iszero(first) && error("Cannot make symplectic basis. Often this is due to the fact that the stabilizers are not maximal and therefore the centralizer still containing part of the isotropic subspace.")
            for c in 2:num_prod
                if !iszero(prod[first, c])
                    L[c, :] += F(prod[first, c]^-1) * L[1, :]
                end
            end
            prod = hcat(L[:, n + 1:end], -L[:, 1:n]) * transpose(L)
            # println("after")
            # display(prod)
            push!(logs, (L[1, :], L[first, :]))
            L = L[setdiff(1:nrows(L), [1, first]), :]
        end
        # display(logs)
    else
        # this does indeed grow smaller each iteration
        while nrows(L) >= 2
            # the columns in prod give the commutation relationships between the provided
            # logical operators; they ideally should only consist of {X_1, Z_i} pairs
            # so there should only be one nonzero element in each column
            prod = hcat(L[:, n + 1:end], -L[:, 1:n]) * transpose(L)
            # println("before")
            # display(prod)
            num_prod = ncols(prod)
            first = 0
            for c in 1:num_prod
                if !iszero(prod[1, c])
                    if iszero(first)
                        first = c
                    else
                        L[c, :] += L[first, :]
                    end
                end
            end
            iszero(first) && error("Cannot make symplectic basis. Often this is due to the fact that the stabilizers are not maximal and therefore the centralizer still containing part of the isotropic subspace.")
            for c in 2:num_prod
                if !iszero(prod[first, c])
                    L[c, :] += L[1, :]
                end
            end
            prod = hcat(L[:, n + 1:end], -L[:, 1:n]) * transpose(L)
            # println("after")
            # display(prod)
            push!(logs, (L[1, :], L[first, :]))
            L = L[setdiff(1:nrows(L), [1, first]), :]
        end
        # display(logs)
    end
    return logs
end

_test_logicals_relationships(S::T) where {T <: AbstractSubsystemCode} =
    _test_logicals_relationships(LogicalTrait(T), S)
function _test_logicals_relationships(::HasLogicals, S::AbstractSubsystemCode)
    prod = hcat(S.logs_mat[:, S.n + 1:end], -S.logs_mat[:, 1:S.n]) * transpose(S.logs_mat)
    display(prod)
    return nothing
end
_test_logicals_relationships(::HasNoLogicals, S) = error("Type $(typeof(S)) has no logicals.")

"""
    is_logical(S::AbstractSubsystemCode, v::CTMatrixTypes)

Return `true` if the vector `v` is a logical operator for `S`.
"""
is_logical(S::T, v::CTMatrixTypes) where {T <: AbstractSubsystemCode} = is_logical(LogicalTrait(T),
    S, v)
function is_logical(::HasLogicals, S::AbstractSubsystemCode, v::CTMatrixTypes)
    nc = ncols(S.logs_mat)
    are_symplectic_orthogonal(S.stabs, v) || return false
    size(v) == (1, nc) && (return !are_symplectic_orthogonal(S.logs_mat, v);)
    size(v) == (nc, 1) && (return !are_symplectic_orthogonal(S.logs_mat, transpose(v));)
    throw(ArgumentError("Vector to be tested is of incorrect dimension."))
end
is_logical(::HasNoLogicals, S::AbstractSubsystemCode, v::CTMatrixTypes) =
    error("Type $(typeof(S)) has no logicals.")

"""
    is_gauge(S::AbstractSubsystemCode, v::CTMatrixTypes)

Return `true` if the vector `v` is a gauge operator for `S`.
"""
is_gauge(S::T, v::CTMatrixTypes) where {T <: AbstractSubsystemCode} = is_gauge(GaugeTrait(T), S, v)
function is_gauge(::HasGauges, S::AbstractSubsystemCode, v::CTMatrixTypes)
    nc = ncols(S.g_ops_mat)
    are_symplectic_orthogonal(S.stabs, v) || return false
    size(v) == (1, nc) && (return !iszero(S.g_ops_mat * transpose(v));)
    size(v) == (nc, 1) && (return !iszero(S.g_ops_mat * v);)
    throw(ArgumentError("Vector to be tested is of incorrect dimension."))
end
is_gauge(::HasNoGauges, S::AbstractSubsystemCode, v::CTMatrixTypes) = error("Type $(typeof(S)) has no gauges.")

# TODO: make uniform with approach in function above and linear code
"""
    syndrome(S::AbstractSubsystemCode, v::CTMatrixTypes)

Return the syndrome of the vector `v` with respect to the stabilizers of `S`.
"""
function syndrome(S::AbstractSubsystemCode, v::CTMatrixTypes)
    (size(v) != (2 * S.n, 1) && size(v) != (1, 2 * S.n)) &&
        throw(ArgumentError("Vector to be tested is of incorrect dimension; expected length $(2 * n), received: $(size(v))."))
    # base_ring(v) == field(S) || error("Vector must have the same base ring as the stabilizers.")

    nrows(v) != 1 || return S.stabs * transpose(v)
    return S.stabs * v
end

"""
    X_syndrome(S::AbstractSubsystemCode, v::CTMatrixTypes)

Return the syndrome of the vector `v` with respect to the `X` stabilizers of the
CSS code.
"""
X_syndrome(S::T, v::CTMatrixTypes) where {T <: AbstractSubsystemCode} = X_syndrome(CSSTrait(T),
    S, v)
function X_syndrome(::IsCSS, S::AbstractSubsystemCode, v::CTMatrixTypes)
    length(v) == 2 * S.n && (v = v[S.n + 1:end];)
    (size(v) != (S.n, 1) && size(v) != (1, S.n)) &&
        error("Vector to be tested is of incorrect dimension; expected length $n, received: $(size(v)).")
    base_ring(v) == S.F || error("Vector must have the same base ring as the stabilizers.")

    nrows(v) != 1 || return S.X_stabs * transpose(v)
    return S.X_stabs * v
end
X_syndrome(::IsNotCSS, S::AbstractSubsystemCode, v::CTMatrixTypes) =
    error("Only valid for CSS codes.")

"""
    Z_syndrome(S::AbstractSubsystemCode, v::CTMatrixTypes)

Return the syndrome of the vector `v` with respect to the `Z` stabilizers of the
CSS code.
"""
Z_syndrome(S::T, v::CTMatrixTypes) where {T <: AbstractSubsystemCode} = Z_syndrome(CSSTrait(T),
    S, v)
function Z_syndrome(::IsCSS, S::AbstractSubsystemCode, v::CTMatrixTypes)
    length(v) == 2 * S.n && (v = v[1:S.n];)
    (size(v) != (n, 1) && size(v) != (1, n)) &&
        error("Vector to be tested is of incorrect dimension; expected length $n, received: $(size(v)).")
    base_ring(v) == S.F || error("Vector must have the same base ring as the stabilizers.")

    nrows(v) != 1 || return S.Z_stabs * transpose(v)
    return S.Z_stabs * v
end
Z_syndrome(::IsNotCSS, S::AbstractSubsystemCode, v::CTMatrixTypes) =
    error("Only valid for CSS codes.")

"""
    promote_logicals_to_gauge(S::AbstractSubsystemCode, pairs::Vector{Int})
    promote_logicals_to_gauge!(S::AbstractSubsystemCode, pairs::Vector{Int})

Add the logical pairs in `pairs` to the gauge operators.
"""
promote_logicals_to_gauge!(S::T, pairs::Vector{Int}) where {T <: AbstractSubsystemCode} =
    promote_logicals_to_gauge!(LogicalTrait(T), S, pairs)
function promote_logicals_to_gauge!(::HasLogicals, S::AbstractSubsystemCode, pairs::Vector{Int})
    pairs = sort!(unique!(pairs))
    # will let this error naturally if pairs contains invalid elements
    S.gauge_ops = S.gauge_ops ∪ S.logicals[pairs]
    S.g_ops_mat = reduce(vcat, [reduce(vcat, S.gauge_ops[i]) for i in 1:length(S.gauge_ops)])
    S.logicals = S.logicals[setdiff!(append!(collect(1:length(S.logicals)), pairs))]
    S.logs_mat = reduce(vcat, [reduce(vcat, S.logicals[i]) for i in 1:length(S.logicals)])
    S.r = S.r + length(pairs)

    if isinteger(S.k)
        S.k = S.k - length(pairs)
    else
        k = BigInt(order(S.F))^S.n // BigInt(p)^(S.stand_k + length(pairs))
        isinteger(k) ? (S.k = round(Int, log(BigInt(p), k));) : (S.k = k;)
    end
    return nothing
end
promote_logicals_to_gauge!(::HasNoLogicals, S::AbstractSubsystemCode, pairs::Vector{Int}) =
    error("Type $(typeof(S)) has no logicals.")
# TODO: this can drop to a graph state, which is why I did what I did before

promote_logicals_to_gauge(S::T, pairs::Vector{Int}) where {T <: AbstractSubsystemCode} =
    promote_logicals_to_gauge(LogicalTrait(T), S, pairs)
promote_logicals_to_gauge(::HasLogicals, S::AbstractSubsystemCode, pairs::Vector{Int}) =
    (S_new = deepcopy(S); promote_logicals_to_gauge!(HasLogicals(), S_new, pairs); return S_new)
promote_logicals_to_gauge(::HasNoLogicals, S::AbstractSubsystemCode, pairs::Vector{Int}) =
    error("Type $(typeof(S)) has no logicals.")

"""
    promote_gauges_to_logical(S::AbstractSubsystemCode, pairs::Vector{Int})
    promote_gauges_to_logical(S::AbstractSubsystemCode, pairs::Vector{Int})

Add the gauge pairs in `pairs` to the logical operators.
"""
promote_gauges_to_logical!(S::T, pairs::Vector{Int}) where {T <: AbstractSubsystemCode} =
    promote_gauges_to_logical!(LogicalTrait(T), S, pairs)
function promote_gauges_to_logical!(::HasLogicals, S::AbstractSubsystemCode, pairs::Vector{Int})
    pairs = sort!(unique!(pairs))
    # will let this error naturally if pairs contains invalid elements
    S.logicals = S.logicals ∪ S.gauge_ops[pairs]
    S.logs_mat = reduce(vcat, [reduce(vcat, S.logicals[i]) for i in 1:length(S.logicals)])
    S.gauge_ops = S.gauge_ops[setdiff!(append!(collect(1:length(S.gauge_ops)), pairs))]
    S.g_ops_mat = reduce(vcat, [reduce(vcat, S.gauge_ops[i]) for i in 1:length(S.gauge_ops)])
    S.r = S.r - length(pairs)

    if isinteger(S.k)
        S.k = S.k + length(pairs)
    else
        k =  BigInt(order(S.F))^S.n // BigInt(p)^(S.stand_k + S.r)
        isinteger(k) ? (S.k = round(Int, log(BigInt(p), k));) : (S.k = k;)
    end
    return nothing
end
promote_gauges_to_logical!(::HasNoLogicals, S::AbstractSubsystemCode, pairs::Vector{Int}) =
    error("Type $(typeof(S)) has no logicals.")
# TODO: this can drop to a graph state, which is why I did what I did before

promote_gauges_to_logical(S::T, pairs::Vector{Int}) where {T <: AbstractSubsystemCode} =
    promote_gauges_to_logical(LogicalTrait(T), S, pairs)
promote_gauges_to_logical(::HasLogicals, S::AbstractSubsystemCode, pairs::Vector{Int}) =
    (S_new = deepcopy(S); promote_gauges_to_logical!(HasLogicals(), S_new, pairs); return S_new)
promote_gauges_to_logical(::HasNoLogicals, S::AbstractSubsystemCode, pairs::Vector{Int}) =
    error("Type $(typeof(S)) has no logicals.")

"""
    swap_X_Z_logicals!(S::AbstractSubsystemCode, pairs::Vector{Int})

Swap the `X` and `Z` logicals specified by `pairs`.    
"""
swap_X_Z_logicals!(S::T, pairs::Vector{Int}) where {T <: AbstractSubsystemCode} =
    swap_X_Z_logicals!(LogicalTrait(T), S, pairs)
function swap_X_Z_logicals!(::HasLogicals, S::AbstractSubsystemCode, pairs::Vector{Int})
    # let indexing check for inbounds naturally
    pairs = sort!(unique!(pairs))
    temp::typeof(S.logicals[1][1])
    for i in pairs
        temp = S.logicals[i][1]
        S.logicals[i][1] = S.logicals[i][2]
        S.logicals[i][2] = temp
    end
    S.logs_mat = reduce(vcat, [reduce(vcat, S.logicals[i]) for i in 1:length(S.logicals)])
    return nothing
end
swap_X_Z_logicals!(::HasNoLogicals, S::AbstractSubsystemCode, pairs::Vector{Int}) =
    error("Type $(typeof(S)) has no logicals.")

"""
    swap_X_Z_gauge_operators!(S::AbstractSubsystemCode, pairs::Vector{Int})

Swap the `X` and `Z` gauge operators specified by `pairs`.
"""
swap_X_Z_gauge_operators!(S::T, pairs::Vector{Int}) where {T <: AbstractSubsystemCode} =
    swap_X_Z_gauge_operators!(GaugeTrait(T), S, pairs)
function swap_X_Z_gauge_operators!(::HasGauges, S::AbstractSubsystemCode, pairs::Vector{Int})
    # let indexing check for inbounds naturally
    pairs = sort!(unique!(pairs))
    temp::typeof(S.logicals[1][1])
    for i in pairs
        temp = S.gauge_ops[i][1]
        S.gauge_ops[i][1] = S.gauge_ops[i][2]
        S.gauge_ops[i][2] = temp
    end
    S.g_ops_mat = reduce(vcat, [reduce(vcat, S.gauge_ops[i]) for i in 1:length(S.gauge_ops)])
    return nothing
end
swap_X_Z_gauge_operators!(::HasNoGauges, S::AbstractSubsystemCode, pairs::Vector{Int}) =
    error("Type $(typeof(S)) has no gauges.")

"""
    are_equivalent(S1::T, S2::T) where T <: AbstractSubsystemCode

Return `true` if the codes are equivalent as symplectic vector spaces.

# Note
* This is not intended to detect if `S1` and `S2` are permutation equivalent.
"""
function are_equivalent(S1::T, S2::T) where {T <: AbstractSubsystemCode}
    (S1.n == S2.n && S1.k == S2.k) || return false
    # can't compare fields directly because they are compared based on ptr addresses
    Int(order(S1.F)) == Int(order(S2.F)) || return false
    if GaugeTrait(T) == HasGauges
        S1.r == S2.r || return false
    end

    # # test stabilizers
    _has_equivalent_row_spaces(S1.stabs, S2.stabs) || return false

    if LogicalTrait(T) == HasLogicals()
        # test logicals
        _has_equivalent_row_spaces(vcat(S1.logs_mat, S1.stabs), vcat(S2.logs_mat, S2.stabs)) || return false
    end

    if GaugeTrait(T) == HasGauges()
        # test gauge operators
        return _has_equivalent_row_spaces(vcat(S1.g_ops_mat, S1.stabs), vcat(S2.g_ops_mat, S2.stabs))
    else
        return true
    end
end

"""
    fix_gauge(::HasGauges, S::AbstractSubsystemCode, pair::Int, which::Symbol)

Return the code induced by adding either the first of `gauge_operators(S)[pair]` to the stabilizers
if `which = :X` or the second if `which = :Z`.
"""
fix_gauge(S::T, pair::Int, which::Symbol) where {T <: AbstractSubsystemCode} =
    fix_gauge(GaugeTrait(T), S, pair, which)
function fix_gauge(::HasGauges, S::AbstractSubsystemCode, pair::Int, which::Symbol)
    if which == :X
        return augment(S, S.gauge_ops[pair][1], verbose = false)
    elseif which == :Z
        return augment(S, S.gauge_ops[pair][2], verbose = false)
    else
        throw(ArgumentError("Unknown type $which"))
    end
end
fix_gauge(::HasNoGauges, S::AbstractSubsystemCode, pair::Int, which::Symbol) =
    error("Type $(typeof(S)) has no gauges.")

function show(io::IO, S::AbstractSubsystemCode)
    if isa(S.k, Integer)
        print(io, "[[$(S.n), $(S.k)")
    else
        print(io, "(($(S.n), $(S.k)")
    end
    !(typeof(S) <: AbstractStabilizerCode) && print(io, ", $(S.r)")
    !ismissing(S.d) && print(io, ", $(S.d)")
    if isa(S.k, Integer)
        print(io, "]]_$(order(S.F))")
    else
        print(io, "))_$(order(S.F))")
    end
    if iszero(S.k)
        if isa(S, GraphStateSubsystemCSS)
            println(io, " CSS subsystem graph state")
        elseif isa(S, GraphStateStabilizerCSS)
            println(io, " CSS graph state")
        else
            println(io, " graph state")
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
        if isa(S, SubsystemCodeCSS) || isa(S, StabilizerCodeCSS) || isa(S, GraphStateStabilizerCSS)
            || isa(S, GraphStateSubsystemCSS)
            
            num_X = num_X_stabs(S)
            if S.overcomplete
                println(io, "X-stabilizer matrix (overcomplete): $num_X × $(S.n)")
            else
                println(io, "X-stabilizer matrix: $num_X × $(S.n)")
            end
            for r in 1:num_X
                print(io, "\t chi($(S.X_signs[r])) ")
                for c in 1:S.n
                    if c != S.n
                        print(io, "$(S.X_stabs[r, c]) ")
                    elseif c == S.n && r != num_X
                        println(io, "$(S.X_stabs[r, c])")
                    else
                        print(io, "$(S.X_stabs[r, c])")
                    end
                end
            end
            println(" ")

            num_Z = num_Z_stabs(S)
            if S.overcomplete
                println(io, "Z-stabilizer matrix (overcomplete): $num_Z × $(S.n)")
            else
                println(io, "Z-stabilizer matrix: $num_Z × $(S.n)")
            end
            for r in 1:num_Z
                print(io, "\t chi($(S.Z_signs[r])) ")
                for c in 1:S.n
                    if c != S.n
                        print(io, "$(S.Z_stabs[r, c]) ")
                    elseif c == S.n && r != num_Z
                        println(io, "$(S.Z_stabs[r, c])")
                    else
                        print(io, "$(S.Z_stabs[r, c])")
                    end
                end
            end
        else
            num_stabs = nrows(S.stabs)
            if S.overcomplete
                println(io, "Stabilizer matrix (overcomplete): $num_stabs × $(2 * S.n)")
            else
                println(io, "Stabilizer matrix: $num_stabs × $(2 * S.n)")
            end
            for r in 1:num_stabs
                print(io, "\t chi($(S.signs[r])) ")
                for c in 1:2 * S.n
                    if c != 2 * S.n
                        print(io, "$(S.stabs[r, c]) ")
                    elseif c == 2 * S.n && r != num_stabs
                        println(io, "$(S.stabs[r, c])")
                    else
                        print(io, "$(S.stabs[r, c])")
                    end
                    if c == S.n
                        print(io, "| ")
                    end
                end
            end
        end
    end
end

# TODO: should probably split these in two for type stability
# TODO: iterate over this faster
# TODO: remove quadratic
function _all_stabilizers(S::AbstractStabilizerCode, only_print::Bool = false)
    E = quadraticfield(S)
    all = Vector{typeof(S.stabs)}()
    stabs = S.stabs
    for iter in Base.Iterators.product([0:(Int64(characteristic(S.F)) - 1)
        for _ in 1:nrows(stabs)]...)

        stab = E(iter[1]) * stabs[1, :]
        for r in 2:nrows(stabs)
            if !iszero(iter[r])
                stab += E(iter[r]) * stabs[r, :]
            end
        end
        if only_print
            println(stab)
        else
            push!(all, stab)
        end
    end
    if only_print
        return
    else
        return all
    end
end

"""
    all_stabilizers(S::AbstractSubsystemCode)
    elements(S::AbstractSubsystemCode)

Return a vector of all the elements of the stabilizer group of `S`.
"""
all_stabilizers(S::AbstractSubsystemCode) = _all_stabilizers(S, false)
elements(S::AbstractSubsystemCode) = all_stabilizers(S)

"""
    print_all_stabilizers(S::AbstractSubsystemCode)
    print_all_elements(S::AbstractSubsystemCode)

Print all the elements of the stabilizer group of `S`.
"""
print_all_stabilizers(S::AbstractSubsystemCode) = _all_stabilizers(S, true)
print_all_elements(S::AbstractSubsystemCode) = print_all_stabilizers(S)


# function Singletonbound()
#     n - k ≧ 2(d - 1) for stabilizer
#     k + r ≤ n - 2d + 2 for subsystem
# they should be the same but with r = 0 for stabilizer
# end

"""
    permute_code(S::AbstractSubsystemCode, σ::Union{PermGroupElem, Perm{Int}, Vector{Int}})
    permute_code!(S::AbstractSubsystemCode, σ::Union{PermGroupElem, Perm{Int}, Vector{Int}})

Return the code permuted by `σ`.

# Notes
* If `σ` is a vector, it is interpreted as the desired column order for the stabilizers of `S`.
"""
function permute_code!(S::AbstractSubsystemCode, σ::Union{PermGroupElem, Perm{Int}, Vector{Int}})
    perm1 = transpose(permutation_matrix(S.F, typeof(σ) <: Perm ? σ.d : σ))
    perm = perm1 ⊕ perm1
    S.stabs = S.stabs * perm
    S.char_vec .= data.(Array(transpose(perm))) * S.char_vec
    S.stabs_stand = S.stabs_stand * perm
    S.P_stand = ismissing(S.P_stand) ? perm : S.P_stand * perm

    W = typeof(S)
    if CSSTrait(W) == IsCSS()
        S.X_stabs = S.X_stabs * perm1
        S.Z_stabs = S.Z_stabs * perm1
    end

    if LogicalTrait(W) == HasLogicals()
        S.logs_mat = S.logs_mat * perm
        for (i, pair) in enumerate(S.logicals)
            S.logicals[i] = (pair[1] * perm, pair[2] * perm)
        end
    end

    if GaugeTrait(W) == HasGauges()
        S.g_ops_mat = S.g_ops_mat * perm
        for (i, pair) in enumerate(S.gauge_ops)
            S.gauge_ops[i] = (pair[1] * perm, pair[2] * perm)
        end
    end
    return nothing
end
permute_code(S::AbstractSubsystemCode, σ::Union{PermGroupElem, Perm{Int}, Vector{Int}}) =
    (S_new = copy(S); permute_code!(S_new, σ); return S_new)

"""
    augment(S::AbstractSubsystemCode, row::CTMatrixTypes; verbose::Bool = true, logs_alg::Symbol = :sys_eqs)

Return the code created by added `row` to the stabilizers of `S`.

# Notes
* The goal of this function is to track how the logical operators update given the new stabilizer.
  The unaffected logical operators are kept during the update and only those which don't commute
  with the new stabilizer are recomputed.
"""
function augment(S::AbstractSubsystemCode, row::CTMatrixTypes; verbose::Bool = true, logs_alg::Symbol = :sys_eqs)

    logs_alg ∈ (:sys_eqs, :VS) || throw(ArgumentError("Unrecognized logicals algorithm"))
    typeof(S.stabs) == typeof(row) || throw(ArgumentError("Vector of different (Julia) type than stabilizers"))
    iszero(row) && return S
    nrows(row) == 1 || throw(ArgumentError("Only one stabilizer may be passed in at a time."))

    # stabilizers
    prod = hcat(S.stabs[:, S.n + 1:end], -S.stabs[:, 1:S.n]) * transpose(row)
    if iszero(prod)
        verbose && println("Vector is already in the stabilizer group. Nothing to update.")    
        S_new = deepcopy(S)
        S_new.stabs = vcat(S.stabs, row)
        S_new.overcomplete = true
        return S_new
    else
        stabs_to_keep = Vector{Int}()
        for i in 1:nrows(S.stabs)
            iszero(prod[i]) && append!(stabs_to_keep, i, i + 1)
        end

        if isempty(stabs_to_keep)
            verbose && println("The vector anticommutes with all stabilizers. The new stabilizer group is just the vector.")
            stabs = row
        else
            update = setdiff(1:nrows(S.stabs), stabs_to_keep)
            if verbose
                if isempty(update)
                    println("No stabilizers requiring updating")
                else
                    println("Stabilizers requiring updating:")
                    display(update)
                end
            end
            isempty(update) ? (stabs = S.stabs;) : (stabs = S.stabs[stabs_to_keep, :];)
        end
    end

    # logicals
    if LogicalTrait(typeof(S)) == HasLogicals()
        prod = hcat(S.logs_mat[:, S.n + 1:end], -S.logs_mat[:, 1:S.n]) * transpose(row)
        logs_to_keep = Vector{Int}()
        log_pairs_to_keep = Vector{Int}()
        pair = 1
        for i in 1:2:nrows(S.logs_mat)
            # this is predicated on the idea that the pairs are stacked together in this matrix
            # TODO: write a unit test checking this never changes
            if iszero(prod[i]) && iszero(prod[i + 1])
                append!(logs_to_keep, i, i + 1)
                append!(log_pairs_to_keep, pair)
            end
            pair += 1
        end

        if isempty(logs_to_keep)
            verbose && println("The vector anticommutes with all logical pairs.")
            logs = zero_matrix(S.F, 1, 2 * S.n)
        else
            update = setdiff(1:length(S.logicals), log_pairs_to_keep)
            if verbose
                if isempty(update)
                    println("No logical pairs requiring updating")
                else
                    println("Logical pairs requiring updating:")
                    display(update)
                end
            end
            isempty(update) ? (logs = S.logs_mat;) : (logs = S.logs_mat[logs_to_keep, :];)
        end
    else
        logs = zero_matrix(S.F, 1, 2 * S.n)
    end

    # gauges
    if GaugeTrait(typeof(S)) == HasGauges()
        prod = hcat(S.g_ops_mat[:, S.n + 1:end], -S.g_ops_mat[:, 1:S.n]) * transpose(row)
        g_ops_to_keep = Vector{Int}()
        g_op_pairs_to_keep = Vector{Int}()
        pair = 1
        for i in 1:2:nrows(S.g_ops_mat)
            # this is predicated on the idea that the pairs are stacked together in this matrix
            # TODO: write a unit test checking this never changes
            if iszero(prod[i]) && iszero(prod[i + 1])
                append!(g_ops_to_keep, i, i + 1)
                append!(g_op_pairs_to_keep, pair)
            end
            pair += 1
        end

        if isempty(g_ops_to_keep)
            verbose && println("The vector anticommutes with all gauge operator pairs.")
            gauge_ops = zero_matrix(S.F, 1, 2 * S.n)
        else
            update = setdiff(1:length(S.gauge_ops), g_op_pairs_to_keep)
            if verbose
                if isempty(update)
                    println("No gauge operator pairs requiring updating")
                else
                    println("Gauge operator pairs requiring updating:")
                    display(update)
                end
            end
            isempty(update) ? (gauge_ops = S.gauge_ops;) : (gauge_ops = S.gauge_ops[g_ops_to_keep, :];)
        end
    else
        gauge_ops = zero_matrix(S.F, 1, 2 * S.n)
    end

    # compute newly opened degrees of freedom
    temp = _remove_empty(vcat(stabs, logs, gauge_ops), :rows)
    rnk_temp, temp = right_kernel(hcat(temp[:, S.n + 1:end], -temp[:, 1:S.n]))
    if ncols(temp) == rnk_temp
        temp = transpose(temp)
    else
        # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
        nr = nrows(temp)
        temp_tr = zero_matrix(base_ring(temp), rnk_temp, nr)
        for r in 1:nr
            for c in 1:rnk_temp
                !iszero(temp[r, c]) && (temp_tr[c, r] = temp[r, c];)
            end
        end
        temp = temp_tr
    end
    # BUG: did I rename new_sym_stabs to temp_tr?
    temp = _quotient_space(temp, new_sym_stabs, logs_alg)
    new_logs = _make_pairs(temp)
    return SubsystemCode(stabs, vcat(logs, new_logs), gauge_ops, S.char_vec)
end

"""
    expurgate(S::AbstractStabilizerCode, rows::Vector{Int}; verbose::Bool = true, logs_alg::Symbol = :sys_eqs)

Return the code created by removing the stabilizers indexed by `rows`.

# Notes
* The goal of this function is to track how the logical operators update through this process.
  Here, the original logical pairs are kept and an appropriate number of new pairs are added.
"""
function expurgate(S::AbstractSubsystemCode, rows::Vector{Int}; verbose::Bool = true, logs_alg::Symbol = :sys_eqs)

    logs_alg ∈ (:sys_eqs, :VS) || throw(ArgumentError("Unrecognized logicals algorithm"))
    num_stabs = nrows(S.stabs)
    rows ⊆ 1:num_stabs || throw(ArgumentError("Argument `rows` not a subset of the number of stabilizers."))
    rows == 1:num_stabs && throw(ArgumentError("Cannot remove all stabilizers"))

    verbose && println("Removing stabilizers: $rows")
    new_stabs = S.stabs[setdiff(1:num_stabs, rows), :]
    temp = new_stabs
    if LogicalTrait(typeof(S)) == HasLogicals()
        temp = vcat(temp, S.logs_mat)
    end
    if GaugeTrait(typeof(S)) == HasGauges()
        temp = vcat(temp, S.g_ops_mat)
    end
    rnk_H, H = right_kernel(hcat(temp[:, S.n + 1:end], -temp[:, 1:S.n]))
    if ncols(H) == rnk_H
        H_tr = transpose(H)
    else
        # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
        nr = nrows(H)
        H_tr = zero_matrix(base_ring(H), rnk_H, nr)
        for r in 1:nr
            for c in 1:rnk_H
                !iszero(H[r, c]) && (H_tr[c, r] = H[r, c];)
            end
        end
    end

    new_logs = _quotient_space(H_tr, new_stabs, logs_alg)
    if iszero(new_logs)
        verbose && println("No new logicals need to be add")
        S_new = deepcopy(S)
        S_new.stabs = new_stabs
        rank(new_stabs) == nrows(new_stabs) ? (S.overcomplete = false;) : (S.overcomplete = true;)
        return S_new
    else
        new_log_pairs = _make_pairs(new_logs)
        verbose && println("New logical pairs:")
        verbose && display(new_log_pairs)
        if GaugeTrait(typeof(S)) == HasGauges()
            return SubsystemCode(new_stabs, vcat(S.logs_mat, reduce(vcat, reduce(vcat, new_log_pairs))), S.g_ops_mat, S.char_vec)
        else
            S_new = StabilizerCode(new_stabs, char_vec = S.char_vec)
            set_logicals!(S_new, vcat(S.logs_mat, reduce(vcat, reduce(vcat, new_log_pairs))))
            return S_new
        end
    end
end

function _standard_form_stabilizer(M::CTMatrixTypes)
    stabs = deepcopy(M)
    # if the stabilizer is overdetermined, remove unnecessary rows
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
    # put S in standard form
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

function _logicals_standard_form(stabs::CTMatrixTypes, n::Int, k::Int, r::Int, P::Union{Missing,
    CTMatrixTypes})

    F = base_ring(stabs)
    F_one = F(1)
    logs = zero_matrix(F, 2k, 2n)

    E = stabs[r + 1:n - k, 2n - k + 1:2n]
    C1 = stabs[1:r, n + r + 1:2n - k]
    C1_E = C1 * E
    for i in 1:k
        # put I in a couple of places
        logs[i, n - k + i] = F_one
        logs[k + i, 2n - k + i] = F_one

        # put in E^T
        for j in 1:(n - k - r)
            logs[i, j + r] = stabs[r + j, 2n - k + i]
        end

        for j in 1:r
            # put in E^T * C1^T + C2^T
            logs[i, n + j] = C1_E[j, i] + stabs[j, 2n - k + i]
            # put in A2^T
            logs[k + i, n + j] = stabs[j, n - k + i]
        end
    end
    return ismissing(P) ? logs : logs * P
end
