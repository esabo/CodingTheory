# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

function _css_distance_cache_keys(which::Symbol)
    if which == :X
        return (
            exact=:dx, lower=:l_bound_dx, upper=:u_bound_dx,
            witness=:X_minimum_distance_witness,
            upper_witness=:X_minimum_distance_upper_bound_witness)
    elseif which == :Z
        return (
            exact=:dz, lower=:l_bound_dz, upper=:u_bound_dz,
            witness=:Z_minimum_distance_witness,
            upper_witness=:Z_minimum_distance_upper_bound_witness)
    elseif which == :full
        return (
            exact=:d, lower=:l_bound, upper=:u_bound,
            witness=:minimum_distance_witness,
            upper_witness=:minimum_distance_upper_bound_witness)
    end
    throw(ArgumentError("Expected `which` to be `:full`, `:X`, or `:Z`."))
end

function _validate_css_bound_witness(
    S::AbstractStabilizerCodeCSS, which::Symbol, d::Int, witness::CTMatrixTypes
)
    size(witness) == (1, 2 * S.n) ||
        throw(ArgumentError("A CSS distance witness must be a 1 × $(2 * S.n) symplectic row."))
    base_ring(witness) == S.F ||
        throw(ArgumentError("The witness must have the same base field as the code."))
    witness_weight = which == :full ?
        count(i -> !iszero(witness[1, i]) ||
            !iszero(witness[1, S.n + i]), 1:S.n) :
        wt(witness)
    witness_weight == d ||
        throw(ArgumentError("The witness has weight $witness_weight, not the supplied bound $d."))
    is_logical(S, witness) ||
        throw(ArgumentError("The supplied upper-bound witness is not a logical operator."))
    which == :X && !iszero(witness[:, S.n + 1:2 * S.n]) &&
        throw(ArgumentError("An X-distance witness must have zero Z support."))
    which == :Z && !iszero(witness[:, 1:S.n]) &&
        throw(ArgumentError("A Z-distance witness must have zero X support."))
    return nothing
end

function _refresh_css_full_bounds!(S::AbstractStabilizerCodeCSS)
    lower_x = get(S.cache, :l_bound_dx, 1)
    lower_z = get(S.cache, :l_bound_dz, 1)
    upper_x = get(S.cache, :u_bound_dx, S.n)
    upper_z = get(S.cache, :u_bound_dz, S.n)
    lower = min(lower_x, lower_z)
    upper = min(upper_x, upper_z, get(S.cache, :u_bound, S.n))
    lower <= upper || error("Inconsistent CSS distance bounds: lower $lower exceeds upper $upper.")
    S.cache[:l_bound] = lower
    S.cache[:u_bound] = upper

    if lower == upper
        S.cache[:d] = lower
        if !haskey(S.cache, :minimum_distance_upper_bound_witness)
            source = upper_x <= upper_z ? :X : :Z
            source_keys = _css_distance_cache_keys(source)
            if haskey(S.cache, source_keys.upper_witness)
                S.cache[:minimum_distance_witness] =
                    S.cache[source_keys.upper_witness]
            elseif haskey(S.cache, source_keys.witness)
                S.cache[:minimum_distance_witness] =
                    S.cache[source_keys.witness]
            end
        else
            S.cache[:minimum_distance_witness] =
                S.cache[:minimum_distance_upper_bound_witness]
        end
    end
    return lower, upper
end

function _close_css_sector_if_proven!(
    S::AbstractStabilizerCodeCSS, which::Symbol
)
    keys = _css_distance_cache_keys(which)
    lower = get(S.cache, keys.lower, 1)
    upper = get(S.cache, keys.upper, S.n)
    if lower == upper
        S.cache[keys.exact] = lower
        if haskey(S.cache, keys.upper_witness)
            S.cache[keys.witness] = S.cache[keys.upper_witness]
        end
    end
    _refresh_css_full_bounds!(S)
    return lower, upper
end

"""
$(TYPEDSIGNATURES)

Return `nothing` after tightening a certified lower bound. A full-distance lower bound also applies to
both CSS sectors because `d = min(d_X, d_Z)`.
"""
function set_minimum_distance_lower_bound!(
    S::AbstractStabilizerCodeCSS, lower::Int; which::Symbol=:full
)
    keys = _css_distance_cache_keys(which)
    upper = get(S.cache, keys.upper, S.n)
    1 <= lower <= upper ||
        throw(DomainError(lower,
            "The lower bound must lie between 1 and the current upper bound $upper."))

    if which == :full
        set_minimum_distance_lower_bound!(S, lower; which=:X)
        set_minimum_distance_lower_bound!(S, lower; which=:Z)
    else
        S.cache[keys.lower] = max(get(S.cache, keys.lower, 1), lower)
        _close_css_sector_if_proven!(S, which)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return `nothing` after tightening a CSS distance upper bound using a validated logical witness.
"""
function set_minimum_distance_upper_bound!(
    S::AbstractStabilizerCodeCSS, upper::Int, witness::CTMatrixTypes;
    which::Symbol=:full
)
    keys = _css_distance_cache_keys(which)
    lower = get(S.cache, keys.lower, 1)
    lower <= upper <= S.n ||
        throw(DomainError(upper,
            "The upper bound must lie between the current lower bound $lower and $(S.n)."))
    _validate_css_bound_witness(S, which, upper, witness)

    if upper <= get(S.cache, keys.upper, S.n)
        S.cache[keys.upper] = upper
        S.cache[keys.upper_witness] = witness
    end

    if which == :full
        zero_half = zero_matrix(S.F, 1, S.n)
        X_part = hcat(witness[:, 1:S.n], zero_half)
        Z_part = hcat(zero_half, witness[:, S.n + 1:2 * S.n])
        !iszero(X_part) && is_logical(S, X_part) &&
            set_minimum_distance_upper_bound!(
                S, wt(X_part), X_part; which=:X)
        !iszero(Z_part) && is_logical(S, Z_part) &&
            set_minimum_distance_upper_bound!(
                S, wt(Z_part), Z_part; which=:Z)
        _refresh_css_full_bounds!(S)
    else
        _close_css_sector_if_proven!(S, which)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return `nothing` after tightening the certified lower bound on the
`X`-distance. The bound records that no nontrivial `X` logical has smaller
weight.
"""
set_X_minimum_distance_lower_bound!(S::AbstractStabilizerCodeCSS, lower::Int) =
    set_minimum_distance_lower_bound!(S, lower; which=:X)

"""
$(TYPEDSIGNATURES)

Return `nothing` after tightening the certified lower bound on the
`Z`-distance. The bound records that no nontrivial `Z` logical has smaller
weight.
"""
set_Z_minimum_distance_lower_bound!(S::AbstractStabilizerCodeCSS, lower::Int) =
    set_minimum_distance_lower_bound!(S, lower; which=:Z)

"""
$(TYPEDSIGNATURES)

Return `nothing` after tightening the `X`-distance upper bound using
`witness`, which is validated as an `X` logical operator of weight `upper`.
"""
set_X_minimum_distance_upper_bound!(
    S::AbstractStabilizerCodeCSS, upper::Int, witness::CTMatrixTypes
) = set_minimum_distance_upper_bound!(S, upper, witness; which=:X)

"""
$(TYPEDSIGNATURES)

Return `nothing` after tightening the `Z`-distance upper bound using
`witness`, which is validated as a `Z` logical operator of weight `upper`.
"""
set_Z_minimum_distance_upper_bound!(
    S::AbstractStabilizerCodeCSS, upper::Int, witness::CTMatrixTypes
) = set_minimum_distance_upper_bound!(S, upper, witness; which=:Z)

"""
$(TYPEDSIGNATURES)

Return user-supplied physical-qubit permutation generators available to
distance preprocessors. Quantum code families do not generate these
automatically yet.
"""
distance_automorphisms(S::AbstractStabilizerCodeCSS) =
    get(S.cache, :distance_automorphisms, Vector{Vector{Int}}())

"""
$(TYPEDSIGNATURES)

Return `nothing` after registering physical-qubit permutations for distance searches. Validation checks
that each permutation preserves both CSS stabilizer row spaces.
"""
function set_distance_automorphisms!(
    S::AbstractStabilizerCodeCSS,
    permutations::Vector{<:AbstractVector{<:Integer}};
    validate::Bool=true
)
    normalized = Vector{Vector{Int}}()
    expected = collect(1:S.n)
    H_X, H_Z = X_stabilizers(S), Z_stabilizers(S)
    for permutation in permutations
        σ = Int.(permutation)
        length(σ) == S.n && sort(σ) == expected ||
            throw(ArgumentError("Each quantum-code automorphism must permute 1:$(S.n)."))
        if validate
            _has_equivalent_row_spaces(H_X, H_X[:, σ]) ||
                throw(ArgumentError("A supplied permutation does not preserve the X-stabilizer row space."))
            _has_equivalent_row_spaces(H_Z, H_Z[:, σ]) ||
                throw(ArgumentError("A supplied permutation does not preserve the Z-stabilizer row space."))
        end
        push!(normalized, σ)
    end
    S.cache[:distance_automorphisms] = unique(normalized)
    return nothing
end
