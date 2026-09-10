# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

function _quantum_qudit_indices(S::AbstractSubsystemCode, qudits)
    indices = qudits isa Integer ?
        [Int(qudits)] : sort!(unique!(collect(Int, qudits)))
    all(q -> 1 <= q <= S.n, indices) ||
        throw(ArgumentError("Qudit indices must lie in 1:$(S.n)."))
    return indices
end

function _compatible_quantum_fields(
    F::CTFieldTypes, K::CTFieldTypes
)
    return characteristic(F) == characteristic(K) &&
        degree(F) == degree(K) && order(F) == order(K)
end

function _coerce_quantum_matrix(
    M::CTMatrixTypes, F::CTFieldTypes
)
    source = _code_matrix_base_ring(M)
    source == F && return _dense_code_matrix(M, F)
    _compatible_quantum_fields(source, F) ||
        throw(ArgumentError("The operator uses an incompatible field."))
    degree(F) == 1 ||
        throw(ArgumentError(
            "Equivalent extension fields require an explicit basis map."))
    return matrix(F, nrows(M), ncols(M), [
        F(M[r, c] isa Integer ?
            Int(M[r, c]) : Int(lift(Nemo.ZZ, M[r, c])))
        for r in 1:nrows(M) for c in 1:ncols(M)
    ])
end

function _coerce_symplectic_rows(
    S::AbstractSubsystemCode, M::CTMatrixTypes
)
    normalized = _normalize_quantum_matrix(M)
    if size(normalized) == (2S.n, 1)
        normalized = transpose(normalized)
    end
    ncols(normalized) == 2S.n ||
        throw(ArgumentError(
            "Expected symplectic rows of length $(2S.n)."))
    return _coerce_quantum_matrix(normalized, S.F)
end

function _require_phase_free_surgery(S::AbstractSubsystemCode)
    isempty(character_vector(S)) ||
        throw(ArgumentError(
            "This operation currently requires an empty character vector; " *
            "use phase-free generators or save the phase convention separately."))
end

function _symplectic_direct_sum(
    A::CTMatrixTypes, n_A::Int,
    B::CTMatrixTypes, n_B::Int,
    F::CTFieldTypes,
)
    A_dense = _dense_code_matrix(A, F)
    B_dense = _dense_code_matrix(B, F)
    Z_A_B = zero_matrix(F, nrows(A_dense), n_B)
    Z_B_A = zero_matrix(F, nrows(B_dense), n_A)
    result = vcat(
        hcat(A_dense[:, 1:n_A], Z_A_B,
             A_dense[:, n_A + 1:end], Z_A_B),
        hcat(Z_B_A, B_dense[:, 1:n_B],
             Z_B_A, B_dense[:, n_B + 1:end]),
    )
    return (_is_sparse_code_matrix(A) || _is_sparse_code_matrix(B)) ?
        _sparse_code_matrix(result) : result
end

function _direct_sum_character_vector(
    A::AbstractSubsystemCode, B::AbstractSubsystemCode
)
    a = character_vector(A)
    b = character_vector(B)
    isempty(a) && isempty(b) && return missing
    p = Int(characteristic(A.F))
    R, _ = residue_ring(Nemo.ZZ, p == 2 ? 4 : p)
    a_full = isempty(a) ? fill(R(0), 2A.n) : a
    b_full = isempty(b) ? fill(R(0), 2B.n) : b
    return vcat(
        a_full[1:A.n], b_full[1:B.n],
        a_full[A.n + 1:end], b_full[B.n + 1:end],
    )
end

"""
    quantum_direct_sum(A, B)
    A ⊕ B

Return the independent direct sum of two stabilizer/subsystem codes. Additive
generators and sparse storage are preserved; parameters are derived by the
normal constructors rather than assumed to be integral.
"""
function quantum_direct_sum(
    A::AbstractSubsystemCode, B::AbstractSubsystemCode
)
    _compatible_quantum_fields(A.F, B.F) ||
        throw(ArgumentError("Quantum direct sums require the same base field."))
    has_gauges = GaugeTrait(typeof(A)) == HasGauges() ||
        GaugeTrait(typeof(B)) == HasGauges()
    G_A = GaugeTrait(typeof(A)) == HasGauges() ?
        gauge_group(A) : stabilizers(A)
    G_B = GaugeTrait(typeof(B)) == HasGauges() ?
        gauge_group(B) : stabilizers(B)
    sparse_A = _is_sparse_code_matrix(G_A)
    sparse_B = _is_sparse_code_matrix(G_B)
    G_A = _coerce_quantum_matrix(G_A, A.F)
    G_B = _coerce_quantum_matrix(G_B, A.F)
    sparse_A && (G_A = _sparse_code_matrix(G_A))
    sparse_B && (G_B = _sparse_code_matrix(G_B))
    generators = _symplectic_direct_sum(G_A, A.n, G_B, B.n, A.F)
    char_vec = _direct_sum_character_vector(A, B)
    if has_gauges
        iszero(generators) &&
            return StabilizerCode(generators; char_vec=char_vec)
        return SubsystemCode(generators; char_vec=char_vec)
    end
    return StabilizerCode(generators; char_vec=char_vec, logs_alg=:sys_eqs)
end

⊕(A::AbstractSubsystemCode, B::AbstractSubsystemCode) =
    quantum_direct_sum(A, B)
direct_sum(A::AbstractSubsystemCode, B::AbstractSubsystemCode) =
    quantum_direct_sum(A, B)

function _project_symplectic_coordinates(
    M::CTMatrixTypes, n::Int, keep::Vector{Int}, F::CTFieldTypes
)
    dense = _dense_code_matrix(M, F)
    projected = dense[:, vcat(keep, n .+ keep)]
    return _is_sparse_code_matrix(M) ?
        _sparse_code_matrix(projected) : projected
end

function _construct_projected_quantum_code(
    S::AbstractSubsystemCode, generators::CTMatrixTypes
)
    new_n = div(ncols(generators), 2)
    (nrows(generators) == 0 || iszero(generators)) &&
        return StabilizerCode(
            zero_matrix(S.F, 0, 2new_n); logs_alg=:sys_eqs)
    if GaugeTrait(typeof(S)) == HasGauges()
        return SubsystemCode(generators)
    end
    if are_symplectic_orthogonal(generators, generators)
        return StabilizerCode(generators; logs_alg=:sys_eqs)
    end
    # Projection can turn commuting stabilizers into a nonabelian gauge group.
    return SubsystemCode(generators)
end

"""
    puncture(S, qudits)

Delete physical coordinates from every generator. A punctured stabilizer code
may become a subsystem code when the projected generators no longer commute.
Character-vector phases are intentionally rejected until phase transport is
implemented.
"""
function puncture(S::AbstractSubsystemCode, qudits)
    _require_phase_free_surgery(S)
    removed = _quantum_qudit_indices(S, qudits)
    isempty(removed) && return deepcopy(S)
    length(removed) < S.n ||
        throw(ArgumentError("Cannot puncture every physical qudit."))
    keep = setdiff(1:S.n, removed)
    generators = GaugeTrait(typeof(S)) == HasGauges() ?
        gauge_group(S) : stabilizers(S)
    projected =
        _project_symplectic_coordinates(generators, S.n, keep, S.F)
    return _construct_projected_quantum_code(S, projected)
end

function _additively_shorten_generators(
    M::CTMatrixTypes, n::Int, removed::Vector{Int}, F::CTFieldTypes
)
    dense = _dense_code_matrix(M, F)
    constrained = dense[:, vcat(removed, n .+ removed)]
    expanded = _additive_expansion(constrained, F)
    coefficients = _rowspace_kernel(transpose(expanded))
    shortened = _lift_prime_matrix(coefficients, F) * dense
    shortened = _remove_empty(shortened, :rows)
    keep = setdiff(1:n, removed)
    projected = shortened[:, vcat(keep, n .+ keep)]
    return _is_sparse_code_matrix(M) ?
        _sparse_code_matrix(projected) : projected
end

"""
    shorten(S, qudits)

Restrict to additive generator combinations acting trivially on `qudits`,
then delete those coordinates. Unlike `puncture`, shortening a stabilizer code
always remains a stabilizer code.
"""
function shorten(S::AbstractSubsystemCode, qudits)
    _require_phase_free_surgery(S)
    removed = _quantum_qudit_indices(S, qudits)
    isempty(removed) && return deepcopy(S)
    length(removed) < S.n ||
        throw(ArgumentError("Cannot shorten on every physical qudit."))
    generators = GaugeTrait(typeof(S)) == HasGauges() ?
        gauge_group(S) : stabilizers(S)
    shortened =
        _additively_shorten_generators(generators, S.n, removed, S.F)
    return _construct_projected_quantum_code(S, shortened)
end

const _LOCAL_CLIFFORD_DISTANCE_KEYS = (
    :d, :d_bare, :d_dressed, :l_bound, :u_bound,
    :l_bound_bare, :u_bound_bare, :l_bound_dressed, :u_bound_dressed,
)

function _copy_local_clifford_distance_cache!(
    target::AbstractSubsystemCode, source::AbstractSubsystemCode
)
    for key in _LOCAL_CLIFFORD_DISTANCE_KEYS
        haskey(source.cache, key) &&
            (target.cache[key] = deepcopy(source.cache[key]))
    end
    return target
end

"""
    local_fourier(S, qudits)
    swap_X_Z(S, qudits)

Apply the single-qudit Fourier Clifford `(x,z) -> (-z,x)` on the selected
coordinates. For binary codes this is the Hadamard `X`/`Z` swap. Pauli weight
and cached exact distances are preserved.
"""
function local_fourier(S::AbstractSubsystemCode, qudits)
    _require_phase_free_surgery(S)
    selected = _quantum_qudit_indices(S, qudits)
    isempty(selected) && return deepcopy(S)
    generators = GaugeTrait(typeof(S)) == HasGauges() ?
        gauge_group(S) : stabilizers(S)
    transformed = deepcopy(_dense_code_matrix(generators, S.F))
    for q in selected
        x = deepcopy(transformed[:, q:q])
        z = deepcopy(transformed[:, S.n + q:S.n + q])
        transformed[:, q:q] = -z
        transformed[:, S.n + q:S.n + q] = x
    end
    _is_sparse_code_matrix(generators) &&
        (transformed = _sparse_code_matrix(transformed))
    result = _construct_projected_quantum_code(S, transformed)
    return _copy_local_clifford_distance_cache!(result, S)
end

swap_X_Z(S::AbstractSubsystemCode, qudits) =
    local_fourier(S, qudits)

function _centralizing_additive_subgroup(
    generators::CTMatrixTypes, row::CTMatrixTypes, F::CTFieldTypes
)
    dense = _dense_code_matrix(generators, F)
    commutation =
        _trace_symplectic_product_matrix(dense, row, F)
    coefficients = _rowspace_kernel(transpose(commutation))
    centralizing = _lift_prime_matrix(coefficients, F) * dense
    centralizing = _remove_empty(centralizing, :rows)
    if !_additive_row_space_contains(centralizing, row, F)
        centralizing = vcat(centralizing, row)
    end
    return _is_sparse_code_matrix(generators) ?
        _sparse_code_matrix(centralizing) : centralizing
end

"""
    augment(S, row; verbose=true)

Impose a Pauli generator as a new stabilizer constraint. Existing generator
combinations that anticommute with `row` are removed by an additive kernel
calculation before `row` is added. This implements stabilizer measurement and
subsystem gauge fixing without assuming adjacent row pairs.
"""
function augment(
    S::AbstractSubsystemCode, row::CTMatrixTypes; verbose::Bool=true
)
    (size(row) == (1, 2S.n) || size(row) == (2S.n, 1)) ||
        throw(ArgumentError("Expected one symplectic row of length $(2S.n)."))
    row_dense = _coerce_symplectic_rows(S, row)
    iszero(row_dense) && return deepcopy(S)
    generators = GaugeTrait(typeof(S)) == HasGauges() ?
        gauge_group(S) : stabilizers(S)
    if is_stabilizer(S, row_dense)
        verbose && println("Generator is already a stabilizer; nothing changed.")
        return deepcopy(S)
    end
    updated = _centralizing_additive_subgroup(
        generators, row_dense, S.F)
    char_vec = isempty(character_vector(S)) ?
        missing : character_vector(S)
    result = if GaugeTrait(typeof(S)) == HasGauges()
        SubsystemCode(updated; char_vec=char_vec)
    else
        StabilizerCode(updated; char_vec=char_vec, logs_alg=:sys_eqs)
    end
    verbose && println("Added one stabilizer constraint.")
    return result
end

"""
    expurgate(S, rows; verbose=true)

Remove the selected stabilizer presentation rows and reconstruct the code from
the remaining stabilizers and existing gauge pairs. Additive dimensions,
logicals, and caches are recomputed.
"""
function expurgate(
    S::AbstractSubsystemCode, rows::Vector{<:Integer};
    verbose::Bool=true
)
    num_stabs = nrows(stabilizers(S))
    selected = sort!(unique!(Int.(rows)))
    all(r -> 1 <= r <= num_stabs, selected) ||
        throw(ArgumentError("Stabilizer row indices are out of range."))
    isempty(selected) && return deepcopy(S)
    stabs = stabilizers(S)
    stabs_dense = _dense_code_matrix(stabs, S.F)
    kept = stabs_dense[setdiff(1:num_stabs, selected), :]
    _is_sparse_code_matrix(stabs) &&
        (kept = _sparse_code_matrix(kept))
    generators = if GaugeTrait(typeof(S)) == HasGauges()
        _vcat_code_matrices(S.F, kept, gauges_matrix(S))
    else
        kept
    end
    char_vec = isempty(character_vector(S)) ?
        missing : character_vector(S)
    result = if iszero(generators) || nrows(generators) == 0
        StabilizerCode(
            zero_matrix(S.F, 0, 2S.n);
            char_vec=char_vec, logs_alg=:sys_eqs)
    elseif GaugeTrait(typeof(S)) == HasGauges()
        SubsystemCode(generators; char_vec=char_vec)
    else
        StabilizerCode(generators; char_vec=char_vec, logs_alg=:sys_eqs)
    end
    verbose && println("Removed $(length(selected)) stabilizer rows.")
    return result
end

expurgate(
    S::AbstractSubsystemCode, row::Integer; kwargs...
) = expurgate(S, [row]; kwargs...)

"""
    gauge_code(S, additional_generators)

Promote additional Pauli generators into the gauge group of `S`. The
subsystem constructor recomputes the center and all protected/gauge
dimensions additively.
"""
function gauge_code(
    S::AbstractStabilizerCode, additional_generators::CTMatrixTypes
)
    _require_phase_free_surgery(S)
    additional = _coerce_symplectic_rows(S, additional_generators)
    generators = _vcat_code_matrices(
        S.F, stabilizers(S), additional)
    return SubsystemCode(generators)
end
