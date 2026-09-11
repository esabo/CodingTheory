# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

# This file contains finite, verifiable expansion computations.  All
# matrix-based routines below are binary: replacing nonzero field elements by
# bits is not valid over a larger field.

# Re-export the Graphs implementation rather than creating a competing generic
# with the same name.
const laplacian_matrix = Grphs.laplacian_matrix

function _binary_support_matrix(H)
    nr, nc = size(H)
    B = falses(nr, nc)
    if eltype(H) <: Integer || eltype(H) <: Bool
        for c in 1:nc, r in 1:nr
            value = H[r, c]
            (iszero(value) || isone(value)) ||
                throw(ArgumentError("A binary matrix may contain only zeros and ones."))
            B[r, c] = isone(value)
        end
        return B
    end

    F = _code_matrix_base_ring(H)
    Int(order(F)) == 2 ||
        throw(ArgumentError("Expansion computations currently require a matrix over GF(2)."))
    for c in 1:nc, r in 1:nr
        B[r, c] = !iszero(H[r, c])
    end
    return B
end

function _binary_support_matrix(H::SparseMatrixCSC)
    nr, nc = size(H)
    values = SparseArrays.nonzeros(H)
    if eltype(H) <: Integer || eltype(H) <: Bool
        all(value -> iszero(value) || isone(value), values) ||
            throw(ArgumentError("A binary matrix may contain only zeros and ones."))
    else
        F = _code_matrix_base_ring(H)
        Int(order(F)) == 2 ||
            throw(ArgumentError(
                "Expansion computations currently require a matrix over GF(2)."))
    end
    B = falses(nr, nc)
    rows = SparseArrays.rowvals(H)
    for column in 1:nc
        for index in SparseArrays.nzrange(H, column)
            iszero(values[index]) || (B[rows[index], column] = true)
        end
    end
    return B
end

function _quantum_check_support(S::AbstractSubsystemCode, check_type::Symbol)
    if CSSTrait(typeof(S)) == IsCSS()
        check_type in (:X, :Z, :both) ||
            throw(ArgumentError("check_type must be :X, :Z, or :both."))
        if check_type == :X
            return _binary_support_matrix(X_stabilizers(S))
        elseif check_type == :Z
            return _binary_support_matrix(Z_stabilizers(S))
        end
        return vcat(
            _binary_support_matrix(X_stabilizers(S)),
            _binary_support_matrix(Z_stabilizers(S)),
        )
    end

    check_type == :both ||
        throw(ArgumentError("A non-CSS code supports only check_type=:both."))
    stabs = _binary_support_matrix(stabilizers(S))
    n = length(S)
    return stabs[:, 1:n] .| stabs[:, n + 1:2n]
end

"""
$(TYPEDSIGNATURES)

Return the symmetric normalized Laplacian. Isolated vertices have zero
diagonal, following the convention used by `Graphs.jl`.
"""
function normalized_laplacian_matrix(G::Grphs.AbstractGraph)
    A = Grphs.adjacency_matrix(G)
    degrees = Grphs.degree(G)
    inv_sqrt = [iszero(d) ? 0.0 : inv(sqrt(Float64(d))) for d in degrees]
    diagonal = Float64[iszero(d) ? 0.0 : 1.0 for d in degrees]
    D_inv_sqrt = spdiagm(0 => inv_sqrt)
    return spdiagm(0 => diagonal) - D_inv_sqrt * A * D_inv_sqrt
end

"""
$(TYPEDSIGNATURES)

Return the second-smallest eigenvalue of the combinatorial (or normalized)
Laplacian. It is zero for a disconnected graph and for graphs with fewer than
two vertices.
"""
function algebraic_connectivity(
    G::Grphs.AbstractGraph; normalized::Bool=false,
)
    Grphs.nv(G) >= 2 || return 0.0
    L = normalized ? normalized_laplacian_matrix(G) : Grphs.laplacian_matrix(G)
    values = eigvals(Symmetric(Matrix{Float64}(L)))
    return max(0.0, Float64(values[2]))
end

"""
$(TYPEDSIGNATURES)

Return the second-smallest eigenvalue of the normalized graph Laplacian.
"""
normalized_spectral_gap(G::Grphs.AbstractGraph) =
    algebraic_connectivity(G; normalized=true)

"""
$(TYPEDSIGNATURES)

Return the largest absolute adjacency eigenvalue after removing the trivial
spectral-radius eigenvalues. For a connected `d`-regular bipartite graph this
removes both `d` and `-d`.
"""
function nontrivial_adjacency_spectral_radius(G::Grphs.AbstractGraph)
    Grphs.nv(G) >= 2 || return 0.0
    values = eigvals(Symmetric(Matrix{Float64}(Grphs.adjacency_matrix(G))))
    radius = maximum(abs, values)
    tolerance = 100 * eps(Float64) * max(1.0, radius)
    nontrivial = filter(value -> abs(abs(value) - radius) > tolerance, values)
    return isempty(nontrivial) ? 0.0 : maximum(abs, nontrivial)
end

"""
$(TYPEDSIGNATURES)

Return a unit Fiedler vector of the combinatorial Laplacian. For a graph with
fewer than two vertices, return a zero vector.
"""
function fiedler_vector(G::Grphs.AbstractGraph)
    n = Grphs.nv(G)
    n >= 2 || return zeros(Float64, n)
    decomposition =
        eigen(Symmetric(Matrix{Float64}(Grphs.laplacian_matrix(G))))
    return decomposition.vectors[:, 2]
end

function _quantum_tanner_graph(S::AbstractSubsystemCode, check_type::Symbol)
    if CSSTrait(typeof(S)) == IsCSS()
        check_type in (:X, :Z, :both) ||
            throw(ArgumentError("check_type must be :X, :Z, or :both."))
        return check_type == :X ? Tanner_graph_X(S)[1] :
            check_type == :Z ? Tanner_graph_Z(S)[1] : Tanner_graph(S)[1]
    end
    check_type == :both ||
        throw(ArgumentError("A non-CSS code supports only check_type=:both."))
    return Tanner_graph(S)[1]
end

function fiedler_vector(
    S::AbstractSubsystemCode, check_type::Symbol=:both,
)
    key = Symbol(:fiedler_vector_, check_type)
    return get!(S.cache, key) do
        fiedler_vector(_quantum_tanner_graph(S, check_type))
    end
end

"""
$(TYPEDSIGNATURES)

Return whether the selected Tanner graph is connected.
"""
is_topologically_connected(
    S::AbstractSubsystemCode, check_type::Symbol=:both,
) = Grphs.is_connected(_quantum_tanner_graph(S, check_type))

"""
$(TYPEDSIGNATURES)

Return the smallest edge-boundary ratio found by sweeping prefixes of a
Fiedler ordering. This is an upper bound on the graph's edge expansion.
"""
function estimated_edge_expansion(G::Grphs.AbstractGraph)
    n = Grphs.nv(G)
    n >= 2 || return 0.0
    order = sortperm(fiedler_vector(G))
    inside = falses(n)
    cut = 0
    best = Inf
    for set_size in 1:fld(n, 2)
        vertex = order[set_size]
        inside[vertex] = true
        for neighbor in Grphs.neighbors(G, vertex)
            cut += inside[neighbor] ? -1 : 1
        end
        best = min(best, cut / set_size)
    end
    return isfinite(best) ? best : 0.0
end

"""
$(TYPEDSIGNATURES)

Return the smallest external-vertex-boundary ratio found by a Fiedler sweep.
This is an upper bound on vertex expansion.
"""
function estimated_vertex_expansion(G::Grphs.AbstractGraph)
    n = Grphs.nv(G)
    n >= 2 || return 0.0
    order = sortperm(fiedler_vector(G))
    inside = falses(n)
    outside_neighbor_count = zeros(Int, n)
    boundary = 0
    best = Inf
    for set_size in 1:fld(n, 2)
        vertex = order[set_size]
        if outside_neighbor_count[vertex] > 0
            boundary -= 1
        end
        inside[vertex] = true
        for neighbor in Grphs.neighbors(G, vertex)
            inside[neighbor] && continue
            iszero(outside_neighbor_count[neighbor]) && (boundary += 1)
            outside_neighbor_count[neighbor] += 1
        end
        best = min(best, boundary / set_size)
    end
    return isfinite(best) ? best : 0.0
end

"""
$(TYPEDSIGNATURES)

Return Cheeger lower and upper bounds for edge expansion:
`λ₂/2 ≤ h(G) ≤ sqrt(2Δλ₂)`.
"""
function edge_expansion_bounds(G::Grphs.AbstractGraph)
    Grphs.nv(G) >= 2 || return (0.0, 0.0)
    λ₂ = algebraic_connectivity(G)
    Δ = maximum(Grphs.degree(G); init=0)
    return (λ₂ / 2, sqrt(2Δ * λ₂))
end

function _bitpack_columns(H)
    B = _binary_support_matrix(H)
    nr, nc = size(B)
    words = cld(nr, 64)
    packed = zeros(UInt64, nc, words)
    for c in 1:nc, r in 1:nr
        B[r, c] || continue
        packed[c, (r - 1) ÷ 64 + 1] |= UInt64(1) << ((r - 1) % 64)
    end
    return packed, words, nr, nc
end

_packed_weight(words) = sum(count_ones, words; init=0)

function _validate_expansion_parameters(n::Int, γ::Real, A::Real)
    0 <= γ <= 1 || throw(DomainError(γ, "γ must lie in [0, 1]."))
    A >= 0 || throw(DomainError(A, "A must be nonnegative."))
    return min(n, floor(Int, γ * n))
end

function _first_expansion_violation(H, max_size::Int, A::Real)
    packed, words, _, nc = _bitpack_columns(H)
    0 <= max_size <= nc ||
        throw(DomainError(max_size, "max_size must lie between zero and the number of columns."))
    max_size == 0 && return nothing
    neighborhood = zeros(UInt64, max_size + 1, words)
    support = Vector{Int}(undef, max_size)

    function visit(start::Int, depth::Int)
        depth == max_size && return nothing
        for column in start:nc
            next_depth = depth + 1
            for word in 1:words
                neighborhood[next_depth + 1, word] =
                    neighborhood[depth + 1, word] | packed[column, word]
            end
            support[next_depth] = column
            neighbor_count = _packed_weight(@view neighborhood[next_depth + 1, :])
            neighbor_count < A * next_depth &&
                return copy(@view support[1:next_depth])

            # Neighborhoods only grow. Once the largest requested right-hand
            # side is met, every descendant is certified.
            neighbor_count >= A * max_size && continue
            witness = visit(column + 1, next_depth)
            isnothing(witness) || return witness
        end
        return nothing
    end
    return visit(1, 0)
end

"""
$(TYPEDSIGNATURES)

Return a subset of columns violating `|N(S)| ≥ A|S|` for
`1 ≤ |S| ≤ floor(γ*ncols(H))`, or `nothing` if every subset passes.
This is an exact, exponential-time computation over the support graph.
"""
function expansion_witness(H, γ::Real, A::Real)
    max_size = _validate_expansion_parameters(size(H, 2), γ, A)
    return _first_expansion_violation(H, max_size, A)
end

"""
$(TYPEDSIGNATURES)

Return whether one-sided vertex expansion from columns to rows holds exactly.
"""
is_expander(H, γ::Real, A::Real) = isnothing(expansion_witness(H, γ, A))

function is_expander(
    S::AbstractSubsystemCode, γ::Real, A::Real,
    check_type::Symbol=:both,
)
    return is_expander(_quantum_check_support(S, check_type), γ, A)
end

"""
$(TYPEDSIGNATURES)

Return the exact minimum value of `|N(S)|/|S|` separately for every nonempty
subset size through `max_subset_size`.
"""
function bipartite_expansion_profile(
    H; max_subset_size::Integer=fld(size(H, 2), 2),
)
    packed, words, _, nc = _bitpack_columns(H)
    max_size = Int(max_subset_size)
    0 <= max_size <= nc ||
        throw(DomainError(max_subset_size,
            "max_subset_size must lie between zero and the number of columns."))
    minima = fill(Inf, max_size)
    neighborhood = zeros(UInt64, max_size + 1, words)
    function visit(start::Int, depth::Int)
        depth == max_size && return
        for column in start:nc
            next_depth = depth + 1
            for word in 1:words
                neighborhood[next_depth + 1, word] =
                    neighborhood[depth + 1, word] | packed[column, word]
            end
            minima[next_depth] = min(
                minima[next_depth],
                _packed_weight(@view neighborhood[next_depth + 1, :]) / next_depth,
            )
            visit(column + 1, next_depth)
        end
    end
    visit(1, 0)
    return minima
end

function _matrix_from_bipartition(
    G::Grphs.AbstractGraph, left::AbstractVector{<:Integer},
    right::AbstractVector{<:Integer},
)
    isempty(intersect(left, right)) ||
        throw(ArgumentError("The left and right vertex sets must be disjoint."))
    right_index = Dict(vertex => row for (row, vertex) in enumerate(right))
    H = falses(length(right), length(left))
    for (column, vertex) in enumerate(left)
        for neighbor in Grphs.neighbors(G, vertex)
            row = get(right_index, neighbor, 0)
            iszero(row) || (H[row, column] = true)
        end
    end
    return H
end

"""
$(TYPEDSIGNATURES)

Return whether the selected left-to-right bipartite adjacency satisfies the
specified one-sided expansion inequality.
"""
function is_bipartite_expander(
    G::Grphs.AbstractGraph, left::AbstractVector{<:Integer},
    right::AbstractVector{<:Integer}, γ::Real, A::Real,
)
    return is_expander(_matrix_from_bipartition(G, left, right), γ, A)
end

"""
$(TYPEDSIGNATURES)

Return whether both orientations of a binary bipartite adjacency matrix
satisfy their respective one-sided expansion inequalities.
"""
function is_left_right_expander(
    H, γ_left::Real, A_left::Real, γ_right::Real, A_right::Real,
)
    is_expander(H, γ_left, A_left) || return false
    return is_expander(transpose(_binary_support_matrix(H)), γ_right, A_right)
end

function is_left_right_expander(
    G::Grphs.AbstractGraph, left::AbstractVector{<:Integer},
    right::AbstractVector{<:Integer},
    γ_left::Real, A_left::Real, γ_right::Real, A_right::Real,
)
    H = _matrix_from_bipartition(G, left, right)
    return is_left_right_expander(H, γ_left, A_left, γ_right, A_right)
end

is_left_right_expander(
    C::AbstractLinearCode, γ_left::Real, A_left::Real,
    γ_right::Real, A_right::Real,
) = is_left_right_expander(
    parity_check_matrix(C), γ_left, A_left, γ_right, A_right)

function _sweep_bipartite_expansion(H, max_size::Int)
    B = _binary_support_matrix(H)
    nr, nc = size(B)
    max_size == 0 && return 0.0
    G = Grphs.SimpleGraph(nr + nc)
    for c in 1:nc, r in 1:nr
        B[r, c] && Grphs.add_edge!(G, c, nc + r)
    end
    vector = fiedler_vector(G)
    orders = (sortperm(@view vector[1:nc]), sortperm(@view vector[1:nc]; rev=true))
    best = Inf
    for order in orders
        neighborhood = falses(nr)
        for set_size in 1:max_size
            neighborhood .|= @view B[:, order[set_size]]
            best = min(best, count(neighborhood) / set_size)
        end
    end
    return best
end

"""
$(TYPEDSIGNATURES)

Return a one-sided column-to-row vertex-expansion estimate from forward and reverse
Fiedler sweeps. The result is an upper bound on the exact minimum.
"""
function estimated_bipartite_vertex_expansion(
    H; max_subset_size::Integer=fld(size(H, 2), 2),
)
    max_size = Int(max_subset_size)
    0 <= max_size <= size(H, 2) ||
        throw(DomainError(max_subset_size,
            "max_subset_size must lie between zero and the number of columns."))
    return _sweep_bipartite_expansion(H, max_size)
end

function estimated_bipartite_vertex_expansion(
    S::AbstractSubsystemCode, check_type::Symbol=:both;
    max_subset_fraction::Real=0.5,
)
    0 <= max_subset_fraction <= 1 ||
        throw(DomainError(max_subset_fraction,
            "max_subset_fraction must lie in [0, 1]."))
    H = _quantum_check_support(S, check_type)
    max_size = floor(Int, max_subset_fraction * size(H, 2))
    return estimated_bipartite_vertex_expansion(H; max_subset_size=max_size)
end

for function_name in (:estimated_edge_expansion, :estimated_vertex_expansion,
                      :edge_expansion_bounds)
    @eval function $function_name(
        S::AbstractSubsystemCode, check_type::Symbol=:both,
    )
        key = Symbol($(QuoteNode(function_name)), :_, check_type)
        return get!(S.cache, key) do
            $function_name(_quantum_tanner_graph(S, check_type))
        end
    end
end

function _syndrome_coset_distances(H, max_error_weight::Int)
    packed, words, _, nc = _bitpack_columns(H)
    0 <= max_error_weight <= nc ||
        throw(DomainError(max_error_weight,
            "max_error_weight must lie between zero and the number of columns."))
    zero_syndrome = ntuple(_ -> UInt64(0), words)
    distances = Dict{Tuple, Int}(zero_syndrome => 0)
    frontier = Tuple[zero_syndrome]
    for distance in 1:max_error_weight
        next_frontier = Tuple[]
        for syndrome in frontier, column in 1:nc
            next_syndrome =
                ntuple(word -> syndrome[word] ⊻ packed[column, word], words)
            haskey(distances, next_syndrome) && continue
            distances[next_syndrome] = distance
            push!(next_frontier, next_syndrome)
        end
        frontier = next_frontier
        isempty(frontier) && break
    end
    return distances
end

"""
$(TYPEDSIGNATURES)

Return, for every observed syndrome weight, the largest minimum error weight
among syndrome cosets reached through `max_error_weight`. This computes
reduced error weight, not the weight of an arbitrary representative.
"""
function confinement_profile(H; max_error_weight::Integer=4)
    distances = _syndrome_coset_distances(H, Int(max_error_weight))
    profile = Dict{Int, Int}()
    for (syndrome, distance) in distances
        syndrome_weight = _packed_weight(syndrome)
        profile[syndrome_weight] =
            max(get(profile, syndrome_weight, 0), distance)
    end
    return profile
end

function confinement_profile(
    S::AbstractSubsystemCode, check_type::Symbol=:X;
    max_error_weight::Integer=4,
)
    return confinement_profile(_quantum_check_support(S, check_type);
        max_error_weight=max_error_weight)
end

"""
$(TYPEDSIGNATURES)

Return the exact minimum of `|syndrome(e)|/dist(e, ker(H))` over nontrivial
syndrome cosets whose leader has weight at most `max_error_weight`.
"""
function deterministic_QLTC_soundness(H; max_error_weight::Integer=4)
    distances = _syndrome_coset_distances(H, Int(max_error_weight))
    best = Inf
    for (syndrome, distance) in distances
        iszero(distance) && continue
        best = min(best, _packed_weight(syndrome) / distance)
    end
    return best
end

function deterministic_QLTC_soundness(
    S::AbstractSubsystemCode, check_type::Symbol=:X;
    max_error_weight::Integer=4,
)
    return deterministic_QLTC_soundness(_quantum_check_support(S, check_type);
        max_error_weight=max_error_weight)
end

"""
$(TYPEDSIGNATURES)

Return whether the exact restricted deterministic QLTC soundness is at least
`threshold`.
"""
function verify_QLTC_soundness(
    H, threshold::Real; max_error_weight::Integer=4,
)
    threshold >= 0 ||
        throw(DomainError(threshold, "The threshold must be nonnegative."))
    return deterministic_QLTC_soundness(
        H; max_error_weight=max_error_weight) >= threshold
end

function verify_QLTC_soundness(
    S::AbstractSubsystemCode, threshold::Real,
    check_type::Symbol=:X; max_error_weight::Integer=4,
)
    return verify_QLTC_soundness(
        _quantum_check_support(S, check_type), threshold;
        max_error_weight=max_error_weight)
end

"""
$(TYPEDSIGNATURES)

Return whether every nontrivial syndrome coset with leader weight at most
`max_error_weight` has syndrome weight strictly above the target.
"""
function verify_confinement(
    H, target_syndrome_weight::Integer; max_error_weight::Integer=4,
)
    target_syndrome_weight >= 0 ||
        throw(DomainError(target_syndrome_weight,
            "The target syndrome weight must be nonnegative."))
    distances = _syndrome_coset_distances(H, Int(max_error_weight))
    return all(iszero(distance) ||
        _packed_weight(syndrome) > target_syndrome_weight
        for (syndrome, distance) in distances)
end

function verify_confinement(
    S::AbstractSubsystemCode, target_syndrome_weight::Integer,
    check_type::Symbol=:X; max_error_weight::Integer=4,
)
    return verify_confinement(
        _quantum_check_support(S, check_type), target_syndrome_weight;
        max_error_weight=max_error_weight)
end

"""
$(TYPEDSIGNATURES)

Return the confinement profile restricted to syndrome weights below `cutoff`.
The profile is exact for syndrome cosets with leaders through
`max_error_weight`.
"""
function evaluate_single_shot_soundness(
    H, cutoff::Integer; max_error_weight::Integer=5,
)
    cutoff >= 1 || throw(DomainError(cutoff, "The cutoff must be positive."))
    profile = confinement_profile(H; max_error_weight=max_error_weight)
    return Dict(weight => distance for (weight, distance) in profile
        if weight < cutoff)
end

function evaluate_single_shot_soundness(
    S::AbstractSubsystemCode, cutoff::Integer,
    check_type::Symbol=:both; max_error_weight::Integer=5,
)
    return evaluate_single_shot_soundness(
        _quantum_check_support(S, check_type), cutoff;
        max_error_weight=max_error_weight)
end

"""
$(TYPEDSIGNATURES)

Return the confinement profile through `max_error_weight`.
"""
evaluate_confinement(H, max_error_weight::Integer) =
    confinement_profile(H; max_error_weight=max_error_weight)

function evaluate_confinement(
    S::AbstractSubsystemCode, max_error_weight::Integer,
    check_type::Symbol=:both,
)
    return confinement_profile(S, check_type;
        max_error_weight=max_error_weight)
end

"""
$(TYPEDSIGNATURES)

Return the standard consequences of a *certified* `(γ,A)` one-sided
expansion bound for a binary left-regular parity-check matrix. By default the
claim is checked exactly; use `verify=false` only when it was certified
elsewhere.
"""
function sipser_spielman_guarantees(
    H, γ::Real, A::Real; verify::Bool=true,
)
    max_set_size = _validate_expansion_parameters(size(H, 2), γ, A)
    B = _binary_support_matrix(H)
    column_weights = vec(sum(B; dims=1))
    isempty(column_weights) &&
        throw(ArgumentError("The parity-check matrix must have at least one column."))
    all(==(column_weights[1]), column_weights) ||
        throw(ArgumentError("Sipser--Spielman guarantees require left regularity."))
    degree = column_weights[1]
    degree > 0 ||
        throw(ArgumentError("The left degree must be positive."))
    verify && !is_expander(B, γ, A) &&
        throw(ArgumentError("The matrix is not a ($γ, $A)-expander."))
    relative_expansion = A / degree
    linear_distance = relative_expansion > 1 // 2
    guaranteed_decoding = relative_expansion > 3 // 4
    return (
        epsilon=relative_expansion,
        degree=degree,
        has_linear_distance=linear_distance,
        guaranteed_distance=linear_distance ? max_set_size + 1 : 0,
        has_guaranteed_decoding=guaranteed_decoding,
        bit_flip_radius=guaranteed_decoding ? fld(max_set_size, 2) : 0,
    )
end

function sipser_spielman_guarantees(
    S::AbstractSubsystemCode, γ::Real, A::Real,
    check_type::Symbol=:both; verify::Bool=true,
)
    return sipser_spielman_guarantees(
        _quantum_check_support(S, check_type), γ, A; verify=verify)
end

function _packed_column_tuple(H)
    packed, words, _, nc = _bitpack_columns(H)
    return [ntuple(word -> packed[column, word], words) for column in 1:nc]
end

function _binary_image(incidence; max_rank::Int)
    columns = _packed_column_tuple(incidence)
    words = isempty(columns) ? cld(size(incidence, 1), 64) : length(columns[1])
    zero_vector = ntuple(_ -> UInt64(0), words)
    image = Tuple[zero_vector]
    seen = Set{Tuple}((zero_vector,))
    rank = 0
    for column in columns
        candidates = Tuple[
            ntuple(word -> vector[word] ⊻ column[word], words)
            for vector in image
        ]
        all(candidate -> candidate in seen, candidates) && continue
        rank += 1
        rank <= max_rank ||
            throw(ArgumentError(
                "The incoming boundary rank exceeds max_boundary_rank=$max_rank."))
        for candidate in candidates
            candidate in seen && continue
            push!(seen, candidate)
            push!(image, candidate)
        end
    end
    return image
end

"""
$(TYPEDSIGNATURES)

Return the exact restricted ratio
`|boundary*x| / dist(x, image(incoming_boundary))` over binary vectors `x`
of Hamming weight at most `max_weight`. The matrices represent
`Cₖ → Cₖ₋₁` and `Cₖ₊₁ → Cₖ`, respectively. A nontrivial cocycle therefore
correctly gives ratio zero.
"""
function cosystolic_expansion(
    boundary, incoming_boundary;
    max_weight::Integer=4, max_boundary_rank::Integer=20,
)
    B = _binary_support_matrix(boundary)
    incoming = _binary_support_matrix(incoming_boundary)
    size(B, 2) == size(incoming, 1) ||
        throw(ArgumentError("The two boundary maps have incompatible dimensions."))
    iszero(mod.(Int.(B) * Int.(incoming), 2)) ||
        throw(ArgumentError("The supplied maps do not form a chain complex."))
    n = size(B, 2)
    max_wt = Int(max_weight)
    0 <= max_wt <= n ||
        throw(DomainError(max_weight, "max_weight must lie in 0:n."))
    max_rank = Int(max_boundary_rank)
    max_rank >= 0 ||
        throw(DomainError(max_boundary_rank,
            "max_boundary_rank must be nonnegative."))

    syndrome_columns = _packed_column_tuple(B)
    vector_words = cld(n, 64)
    image = _binary_image(incoming; max_rank=max_rank)
    best = Inf
    for weight in 1:max_wt
        for support in Combinatorics.combinations(1:n, weight)
            vector = zeros(UInt64, vector_words)
            syndrome = zeros(UInt64, isempty(syndrome_columns) ?
                cld(size(B, 1), 64) : length(syndrome_columns[1]))
            for coordinate in support
                vector[(coordinate - 1) ÷ 64 + 1] |=
                    UInt64(1) << ((coordinate - 1) % 64)
                for word in eachindex(syndrome)
                    syndrome[word] ⊻= syndrome_columns[coordinate][word]
                end
            end
            distance = minimum(
                sum(count_ones(vector[word] ⊻ representative[word])
                    for word in eachindex(vector); init=0)
                for representative in image
            )
            iszero(distance) && continue
            best = min(best, _packed_weight(syndrome) / distance)
        end
    end
    return best
end
