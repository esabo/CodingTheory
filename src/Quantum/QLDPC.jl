# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

function _quantum_degree_distribution(M::SparseMatrixCSC, n::Int)
    nc = ncols(M)
    nc ∈ (n, 2n) ||
        throw(ArgumentError("A quantum generator matrix must have n or 2n columns."))
    qubit_degrees = zeros(Int, n)
    generator_weights = zeros(Int, nrows(M))
    seen = falses(nrows(M))
    touched = Int[]
    num_sectors = nc == n ? 1 : 2

    for q in 1:n
        empty!(touched)
        for sector in 0:(num_sectors - 1)
            c = q + sector * n
            for ptr in SparseArrays.nzrange(M, c)
                iszero(SparseArrays.nonzeros(M)[ptr]) && continue
                r = SparseArrays.rowvals(M)[ptr]
                if !seen[r]
                    seen[r] = true
                    push!(touched, r)
                    qubit_degrees[q] += 1
                    generator_weights[r] += 1
                end
            end
        end
        for r in touched
            seen[r] = false
        end
    end
    return qubit_degrees, generator_weights
end

function _quantum_degree_distribution(M::SMat, n::Int)
    nc = ncols(M)
    nc ∈ (n, 2n) ||
        throw(ArgumentError("A quantum generator matrix must have n or 2n columns."))
    qubit_degrees = zeros(Int, n)
    generator_weights = zeros(Int, nrows(M))
    for (r, row) in enumerate(M)
        support = Set{Int}()
        for (c, value) in row
            !iszero(value) && push!(support, mod1(c, n))
        end
        generator_weights[r] = length(support)
        for q in support
            qubit_degrees[q] += 1
        end
    end
    return qubit_degrees, generator_weights
end

function _quantum_degree_distribution(M::CTMatrixTypes, n::Int)
    nc = ncols(M)
    nc ∈ (n, 2n) ||
        throw(ArgumentError("A quantum generator matrix must have n or 2n columns."))
    qubit_degrees = zeros(Int, n)
    generator_weights = zeros(Int, nrows(M))
    num_sectors = nc == n ? 1 : 2
    for r in 1:nrows(M), q in 1:n
        acts = false
        for sector in 0:(num_sectors - 1)
            if !iszero(M[r, q + sector * n])
                acts = true
                break
            end
        end
        if acts
            qubit_degrees[q] += 1
            generator_weights[r] += 1
        end
    end
    return qubit_degrees, generator_weights
end

function _quantum_generators(S::AbstractSubsystemCode, generators::Symbol)
    if generators == :stabilizers
        return stabilizers(S)
    elseif generators == :gauges
        GaugeTrait(typeof(S)) == HasGauges() ||
            throw(ArgumentError("A stabilizer code has no gauge-operator generators."))
        return gauges_matrix(S)
    elseif generators == :gauge_group
        GaugeTrait(typeof(S)) == HasGauges() ||
            throw(ArgumentError("A stabilizer code has no nontrivial gauge group."))
        return gauge_group(S)
    end
    throw(ArgumentError(
        "Expected `generators` to be `:stabilizers`, `:gauges`, or `:gauge_group`."))
end

"""
$(TYPEDSIGNATURES)

Return the number of supplied generators acting nontrivially on each physical
qubit. A `Y`-type action contributes one, not two.
"""
function qubit_degrees(
    S::AbstractSubsystemCode; generators::Symbol=:stabilizers
)
    M = _quantum_generators(S, generators)
    return first(_quantum_degree_distribution(M, S.n))
end

"""
$(TYPEDSIGNATURES)

Return the Pauli weight of each supplied stabilizer, gauge, or gauge-group
generator.
"""
function generator_weights(
    S::AbstractSubsystemCode; generators::Symbol=:stabilizers
)
    M = _quantum_generators(S, generators)
    return last(_quantum_degree_distribution(M, S.n))
end

"""
$(TYPEDSIGNATURES)

Return a vector containing the Pauli weight of each stabilizer generator of `S`.
"""
stabilizer_weights(S::AbstractSubsystemCode) =
    generator_weights(S; generators=:stabilizers)

"""
$(TYPEDSIGNATURES)

Return a vector containing the Pauli weight of each gauge generator of `S`.
"""
gauge_weights(S::AbstractSubsystemCode) =
    generator_weights(S; generators=:gauges)

"""
$(TYPEDSIGNATURES)

Return a vector containing the Pauli weight of each generator of the gauge group
of `S`.
"""
gauge_group_weights(S::AbstractSubsystemCode) =
    generator_weights(S; generators=:gauge_group)

"""
$(TYPEDSIGNATURES)

Return the number of `X`-stabilizer generators acting on each physical qubit of
`S`.
"""
X_qubit_degrees(S::AbstractSubsystemCode) =
    first(_quantum_degree_distribution(X_stabilizers(S), S.n))

"""
$(TYPEDSIGNATURES)

Return the number of `Z`-stabilizer generators acting on each physical qubit of
`S`.
"""
Z_qubit_degrees(S::AbstractSubsystemCode) =
    first(_quantum_degree_distribution(Z_stabilizers(S), S.n))

"""
$(TYPEDSIGNATURES)

Return a vector containing the Pauli weight of each `X`-stabilizer generator of
`S`.
"""
X_stabilizer_weights(S::AbstractSubsystemCode) =
    last(_quantum_degree_distribution(X_stabilizers(S), S.n))

"""
$(TYPEDSIGNATURES)

Return a vector containing the Pauli weight of each `Z`-stabilizer generator of
`S`.
"""
Z_stabilizer_weights(S::AbstractSubsystemCode) =
    last(_quantum_degree_distribution(Z_stabilizers(S), S.n))

"""
$(TYPEDSIGNATURES)

Return the maximum generator weight and maximum qubit degree of the supplied
presentation. These are presentation-dependent quantities; redundant
generators are retained.
"""
function quantum_LDPC_parameters(
    S::AbstractSubsystemCode; generators::Symbol=:stabilizers
)
    M = _quantum_generators(S, generators)
    qubits, checks = _quantum_degree_distribution(M, S.n)
    return (
        max_generator_weight=isempty(checks) ? 0 : maximum(checks),
        max_qubit_degree=isempty(qubits) ? 0 : maximum(qubits),
    )
end

"""
$(TYPEDSIGNATURES)

                    generators=:stabilizers)

Return whether this presentation obeys the supplied LDPC degree bounds.
LDPC is an asymptotic family property, so both finite-size thresholds are
required rather than chosen by the library.
"""
function is_quantum_LDPC(
    S::AbstractSubsystemCode;
    max_generator_weight::Int,
    max_qubit_degree::Int,
    generators::Symbol=:stabilizers
)
    max_generator_weight >= 0 ||
        throw(DomainError(max_generator_weight, "Degree bounds must be nonnegative."))
    max_qubit_degree >= 0 ||
        throw(DomainError(max_qubit_degree, "Degree bounds must be nonnegative."))
    parameters = quantum_LDPC_parameters(S; generators=generators)
    return parameters.max_generator_weight <= max_generator_weight &&
        parameters.max_qubit_degree <= max_qubit_degree
end

# Match the established classical LDPC naming for the stabilizer presentation.
variable_degree_distribution(S::AbstractSubsystemCode) = qubit_degrees(S)
check_degree_distribution(S::AbstractSubsystemCode) = stabilizer_weights(S)
degree_distributions(S::AbstractSubsystemCode) =
    (qubit_degrees(S), stabilizer_weights(S))

"""
$(TYPEDSIGNATURES)

Return the stabilizer-generator degree of each physical qubit of `S`.
"""
qubit_degree_distribution(S::AbstractSubsystemCode) = qubit_degrees(S)

"""
$(TYPEDSIGNATURES)

Return the Pauli weight of each stabilizer generator of `S`.
"""
stabilizer_weight_distribution(S::AbstractSubsystemCode) =
    stabilizer_weights(S)

"""
$(TYPEDSIGNATURES)

Return the minimum stabilizer-generator degree among the physical qubits of `S`,
or zero when `S` has no physical qubits.
"""
minimum_qubit_degree(S::AbstractSubsystemCode) =
    isempty(qubit_degrees(S)) ? 0 : minimum(qubit_degrees(S))

"""
$(TYPEDSIGNATURES)

Return the maximum stabilizer-generator degree among the physical qubits of `S`,
or zero when `S` has no physical qubits.
"""
maximum_qubit_degree(S::AbstractSubsystemCode) =
    isempty(qubit_degrees(S)) ? 0 : maximum(qubit_degrees(S))

"""
$(TYPEDSIGNATURES)

Return the maximum Pauli weight among the stabilizer generators of `S`, or zero
when the presentation has no stabilizer generators.
"""
maximum_stabilizer_weight(S::AbstractSubsystemCode) =
    isempty(stabilizer_weights(S)) ? 0 : maximum(stabilizer_weights(S))

variable_degree_distribution(M::CTMatrixTypes) =
    first(_degree_distribution(M))
check_degree_distribution(M::CTMatrixTypes) =
    last(_degree_distribution(M))
degree_distributions(M::CTMatrixTypes) = _degree_distribution(M)
column_bound(M::CTMatrixTypes) =
    isempty(variable_degree_distribution(M)) ? 0 :
    maximum(variable_degree_distribution(M))
row_bound(M::CTMatrixTypes) =
    isempty(check_degree_distribution(M)) ? 0 :
    maximum(check_degree_distribution(M))
column_row_bounds(M::CTMatrixTypes) = (column_bound(M), row_bound(M))
limited(M::CTMatrixTypes) = max(column_bound(M), row_bound(M))
"""
$(TYPEDSIGNATURES)

Return the number of nonzero entries in `M`, equivalently the number of edges
in its Tanner graph.
"""
num_edges(M::CTMatrixTypes) = sum(check_degree_distribution(M))
density(M::CTMatrixTypes) =
    isempty(M) ? 0.0 : num_edges(M) / (nrows(M) * ncols(M))
function is_regular(M::CTMatrixTypes)
    columns, rows = degree_distributions(M)
    return (isempty(columns) || all(==(first(columns)), columns)) &&
        (isempty(rows) || all(==(first(rows)), rows))
end

"""
$(TYPEDSIGNATURES)

Return the `X`-stabilizer degree of each physical qubit of `S`.
"""
X_variable_degree_distribution(S::AbstractSubsystemCode) =
    X_qubit_degrees(S)

"""
$(TYPEDSIGNATURES)

Return the `Z`-stabilizer degree of each physical qubit of `S`.
"""
Z_variable_degree_distribution(S::AbstractSubsystemCode) =
    Z_qubit_degrees(S)

"""
$(TYPEDSIGNATURES)

Return the Pauli weight of each `X`-stabilizer generator of `S`.
"""
X_check_degree_distribution(S::AbstractSubsystemCode) =
    X_stabilizer_weights(S)

"""
$(TYPEDSIGNATURES)

Return the Pauli weight of each `Z`-stabilizer generator of `S`.
"""
Z_check_degree_distribution(S::AbstractSubsystemCode) =
    Z_stabilizer_weights(S)

"""
$(TYPEDSIGNATURES)

Return `(qubit_degrees, stabilizer_weights)` for the `X` sector of `S`.
"""
X_degree_distributions(S::AbstractSubsystemCode) =
    (X_qubit_degrees(S), X_stabilizer_weights(S))

"""
$(TYPEDSIGNATURES)

Return `(qubit_degrees, stabilizer_weights)` for the `Z` sector of `S`.
"""
Z_degree_distributions(S::AbstractSubsystemCode) =
    (Z_qubit_degrees(S), Z_stabilizer_weights(S))

"""
$(TYPEDSIGNATURES)

Return the maximum number of `X`-stabilizer generators acting on a physical
qubit of `S`, or zero when `S` has no physical qubits.
"""
X_column_bound(S::AbstractSubsystemCode) =
    isempty(X_qubit_degrees(S)) ? 0 : maximum(X_qubit_degrees(S))

"""
$(TYPEDSIGNATURES)

Return the maximum number of `Z`-stabilizer generators acting on a physical
qubit of `S`, or zero when `S` has no physical qubits.
"""
Z_column_bound(S::AbstractSubsystemCode) =
    isempty(Z_qubit_degrees(S)) ? 0 : maximum(Z_qubit_degrees(S))

"""
$(TYPEDSIGNATURES)

Return the maximum Pauli weight of an `X`-stabilizer generator of `S`, or zero
when the `X` sector has no generators.
"""
X_row_bound(S::AbstractSubsystemCode) =
    isempty(X_stabilizer_weights(S)) ? 0 : maximum(X_stabilizer_weights(S))

"""
$(TYPEDSIGNATURES)

Return the maximum Pauli weight of a `Z`-stabilizer generator of `S`, or zero
when the `Z` sector has no generators.
"""
Z_row_bound(S::AbstractSubsystemCode) =
    isempty(Z_stabilizer_weights(S)) ? 0 : maximum(Z_stabilizer_weights(S))

"""
$(TYPEDSIGNATURES)

Return `(column_bound, row_bound)` for the `X` sector of `S`.
"""
X_column_row_bounds(S::AbstractSubsystemCode) =
    (X_column_bound(S), X_row_bound(S))

"""
$(TYPEDSIGNATURES)

Return `(column_bound, row_bound)` for the `Z` sector of `S`.
"""
Z_column_row_bounds(S::AbstractSubsystemCode) =
    (Z_column_bound(S), Z_row_bound(S))

"""
$(TYPEDSIGNATURES)

Return the larger of the column and row bounds for the `X` sector of `S`.
"""
X_limited(S::AbstractSubsystemCode) = max(X_column_bound(S), X_row_bound(S))

"""
$(TYPEDSIGNATURES)

Return the larger of the column and row bounds for the `Z` sector of `S`.
"""
Z_limited(S::AbstractSubsystemCode) = max(Z_column_bound(S), Z_row_bound(S))

"""
$(TYPEDSIGNATURES)

Return the fraction of entries in the `X`-stabilizer matrix of `S` that are
nonzero.
"""
X_density(S::AbstractSubsystemCode) = density(X_stabilizers(S))

"""
$(TYPEDSIGNATURES)

Return the fraction of entries in the `Z`-stabilizer matrix of `S` that are
nonzero.
"""
Z_density(S::AbstractSubsystemCode) = density(Z_stabilizers(S))

"""
$(TYPEDSIGNATURES)

Return whether all columns of the `X`-stabilizer matrix have one common degree
and all rows have one common degree.
"""
X_is_regular(S::AbstractSubsystemCode) = is_regular(X_stabilizers(S))

"""
$(TYPEDSIGNATURES)

Return whether all columns of the `Z`-stabilizer matrix have one common degree
and all rows have one common degree.
"""
Z_is_regular(S::AbstractSubsystemCode) = is_regular(Z_stabilizers(S))

column_bound(S::AbstractSubsystemCode) = maximum_qubit_degree(S)
row_bound(S::AbstractSubsystemCode) = maximum_stabilizer_weight(S)
column_row_bounds(S::AbstractSubsystemCode) =
    (column_bound(S), row_bound(S))
limited(S::AbstractSubsystemCode) = max(column_bound(S), row_bound(S))

"""
$(TYPEDSIGNATURES)

Return the number of edges in the stabilizer Tanner graph of `S`. A nontrivial
Pauli action, including `Y`, contributes one edge.
"""
num_edges(S::AbstractSubsystemCode) = sum(stabilizer_weights(S))
function density(S::AbstractSubsystemCode)
    num_checks = nrows(stabilizers(S))
    return num_checks == 0 ? 0.0 : num_edges(S) / (num_checks * S.n)
end
function is_regular(S::AbstractSubsystemCode)
    if CSSTrait(typeof(S)) == IsCSS()
        return X_is_regular(S) && Z_is_regular(S)
    end
    qubits, checks = degree_distributions(S)
    return (isempty(qubits) || all(==(first(qubits)), qubits)) &&
        (isempty(checks) || all(==(first(checks)), checks))
end

"""
$(TYPEDSIGNATURES)

Return `(X_row_bound, X_column_bound, Z_row_bound, Z_column_bound)` for a CSS
code `S`, or `(row_bound, column_bound)` for a non-CSS code.
"""
function check_weights(S::AbstractSubsystemCode)
    if CSSTrait(typeof(S)) == IsCSS()
        return (
            X_row_bound(S), X_column_bound(S),
            Z_row_bound(S), Z_column_bound(S),
        )
    end
    return (row_bound(S), column_bound(S))
end

"""
$(TYPEDSIGNATURES)

Return `(row_bound, column_bound)` for the matrix `M`.
"""
check_weights(M::CTMatrixTypes) = (row_bound(M), column_bound(M))

"""
$(TYPEDSIGNATURES)

Return whether the `X` sector of `S` obeys the supplied check-weight and
qubit-degree bounds.
"""
function is_X_LDPC(
    S::AbstractSubsystemCode; check_bound::Int, column_bound::Int
)
    return X_row_bound(S) <= check_bound &&
        X_column_bound(S) <= column_bound
end

"""
$(TYPEDSIGNATURES)

Return whether the `Z` sector of `S` obeys the supplied check-weight and
qubit-degree bounds.
"""
function is_Z_LDPC(
    S::AbstractSubsystemCode; check_bound::Int, column_bound::Int
)
    return Z_row_bound(S) <= check_bound &&
        Z_column_bound(S) <= column_bound
end

"""
$(TYPEDSIGNATURES)

Return whether the stabilizer presentation of `S` obeys the supplied
check-weight and qubit-degree bounds. For a CSS code, both sectors must obey
the bounds.
"""
function is_LDPC(
    S::AbstractSubsystemCode; check_bound::Int, column_bound::Int
)
    check_bound >= 0 && column_bound >= 0 ||
        throw(DomainError((check_bound, column_bound),
            "Degree bounds must be nonnegative."))
    if CSSTrait(typeof(S)) == IsCSS()
        return is_X_LDPC(
            S; check_bound=check_bound, column_bound=column_bound) &&
            is_Z_LDPC(
                S; check_bound=check_bound, column_bound=column_bound)
    end
    return row_bound(S) <= check_bound &&
        CodingTheory.column_bound(S) <= column_bound
end

"""
$(TYPEDSIGNATURES)

Return whether the rows and columns of `M` obey the supplied degree bounds.
"""
is_LDPC(M::CTMatrixTypes; check_bound::Int, column_bound::Int) =
    row_bound(M) <= check_bound && CodingTheory.column_bound(M) <= column_bound

"""
$(TYPEDSIGNATURES)

Return the classical `LDPCCode` whose parity-check matrix is the
`X`-stabilizer matrix of `S`.
"""
X_LDPC_code(S::AbstractSubsystemCode) = LDPCCode(X_stabilizers(S))

"""
$(TYPEDSIGNATURES)

Return the classical `LDPCCode` whose parity-check matrix is the
`Z`-stabilizer matrix of `S`.
"""
Z_LDPC_code(S::AbstractSubsystemCode) = LDPCCode(Z_stabilizers(S))

"""
$(TYPEDSIGNATURES)

Return `(X_LDPC_code, Z_LDPC_code)`, the classical LDPC codes formed from the
two stabilizer sectors of `S`.
"""
LDPC_codes(S::AbstractSubsystemCode) = (X_LDPC_code(S), Z_LDPC_code(S))
