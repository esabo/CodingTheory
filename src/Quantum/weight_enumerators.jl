# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
    ShorLaflammeWeightEnumerator

Hamming-weight Shor--Laflamme enumerators. `A` is the stabilizer distribution
and `B` is its trace-symplectic dual (the stabilizer normalizer
distribution). Only coefficient dictionaries are stored; operator lists are
never retained.
"""
struct ShorLaflammeWeightEnumerator
    n::Int
    A::HammingWeightEnumerator
    B::HammingWeightEnumerator
end

function Base.show(io::IO, enumerator::ShorLaflammeWeightEnumerator)
    print(io, "Shor–Laflamme Hamming enumerator\nA = ")
    show(io, enumerator.A)
    print(io, "\nB = ")
    show(io, enumerator.B)
end

function _independent_additive_generators(
    M::CTMatrixTypes, F::CTFieldTypes
)
    empty = zero_matrix(F, 0, ncols(M))
    return _additive_quotient_space(empty, M, F)
end

function _hamming_distribution_additive(
    M::CTMatrixTypes, n::Int, F::CTFieldTypes;
    max_terms::Integer=10_000_000,
)
    basis = _dense_code_matrix(
        _independent_additive_generators(M, F), F)
    r = nrows(basis)
    p = Int(characteristic(F))
    num_terms = BigInt(p)^r
    num_terms <= max_terms ||
        throw(ArgumentError(
            "Exact additive enumeration requires $num_terms terms, " *
            "exceeding max_terms=$max_terms."))

    counts = Dict{Int, BigInt}(0 => BigInt(1))
    r == 0 && return counts
    prime_field = _prime_subfield(F)
    digits = zeros(Int, r)
    word = zero_matrix(F, 1, ncols(basis))

    for _ in 2:Int(num_terms)
        position = 1
        while true
            old_digit = digits[position]
            new_digit = mod(old_digit + 1, p)
            digits[position] = new_digit
            delta = mod(new_digit - old_digit, p)
            if delta != 0
                scalar = _lift_prime_element(
                    prime_field(delta), F)
                word += scalar * basis[position:position, :]
            end
            new_digit != 0 && break
            position += 1
        end
        weight = symplectic_weight(word)
        counts[weight] = get(counts, weight, BigInt(0)) + 1
    end
    return counts
end

function _quantum_Krawtchouk(
    j::Int, i::Int, n::Int, alphabet_size::Int
)
    total = BigInt(0)
    lower = max(0, j - (n - i))
    upper = min(j, i)
    for ell in lower:upper
        total += (-BigInt(1))^ell *
            BigInt(alphabet_size - 1)^(j - ell) *
            binomial(BigInt(i), ell) *
            binomial(BigInt(n - i), j - ell)
    end
    return total
end

function _trace_symplectic_MacWilliams(
    A::HammingWeightEnumerator, stabilizer_cardinality::BigInt,
    q::Int,
)
    counts = Dict{Int, BigInt}()
    alphabet_size = q^2
    for j in 0:A.n
        coefficient = sum(
            get(A.counts, i, BigInt(0)) *
                _quantum_Krawtchouk(j, i, A.n, alphabet_size)
            for i in 0:A.n
        )
        coefficient % stabilizer_cardinality == 0 ||
            error("The trace-symplectic MacWilliams transform was not integral.")
        value = div(coefficient, stabilizer_cardinality)
        !iszero(value) && (counts[j] = value)
    end
    return HammingWeightEnumerator(A.n, counts)
end

function _validate_SL_enumerator(
    enumerator::ShorLaflammeWeightEnumerator,
    stabilizer_cardinality::BigInt,
    q::Int,
)
    get(enumerator.A.counts, 0, BigInt(0)) == 1 ||
        error("A valid stabilizer enumerator must have A₀ = 1.")
    get(enumerator.B.counts, 0, BigInt(0)) == 1 ||
        error("A valid normalizer enumerator must have B₀ = 1.")
    all(0 <= w <= enumerator.n && count >= 0
        for (w, count) in enumerator.A.counts) ||
        error("The A enumerator contains an invalid coefficient.")
    all(0 <= w <= enumerator.n && count >= 0
        for (w, count) in enumerator.B.counts) ||
        error("The B enumerator contains an invalid coefficient.")
    all(get(enumerator.A.counts, w, BigInt(0)) <=
        get(enumerator.B.counts, w, BigInt(0)) for w in 0:enumerator.n) ||
        error("The Shor–Laflamme inequality Aᵢ ≤ Bᵢ is violated.")
    sum(values(enumerator.A.counts)) == stabilizer_cardinality ||
        error("The A enumerator cardinality does not match the stabilizer.")
    ambient_cardinality = BigInt(q)^(2 * enumerator.n)
    ambient_cardinality % stabilizer_cardinality == 0 ||
        error("The stabilizer cardinality does not divide the Pauli space.")
    sum(values(enumerator.B.counts)) ==
        div(ambient_cardinality, stabilizer_cardinality) ||
        error("The B enumerator cardinality does not match the normalizer.")
    return enumerator
end

function _cache_SL_distance!(
    S::AbstractStabilizerCode,
    enumerator::ShorLaflammeWeightEnumerator,
)
    S.k == 0 && return
    candidates = [
        w for w in 1:S.n
        if get(enumerator.B.counts, w, BigInt(0)) >
           get(enumerator.A.counts, w, BigInt(0))
    ]
    isempty(candidates) || set_minimum_distance!(S, minimum(candidates))
    return
end

function _set_SL_weight_enumerator!(
    S::AbstractSubsystemCode, A_counts::Dict{Int, BigInt},
)
    A = HammingWeightEnumerator(S.n, copy(A_counts))
    B = _trace_symplectic_MacWilliams(
        A, cardinality(S), Int(order(S.F)))
    enumerator = _validate_SL_enumerator(
        ShorLaflammeWeightEnumerator(S.n, A, B),
        cardinality(S), Int(order(S.F)))
    S.cache[:weight_enum_A] = A
    S.cache[:weight_enum_B] = B
    S.cache[:weight_dist_A] = A.counts
    S.cache[:weight_dist_B] = B.counts
    S.cache[:SL_weight_enum] = enumerator
    S isa AbstractStabilizerCode &&
        _cache_SL_distance!(S, enumerator)
    return S
end

"""
    Shor_Laflamme_weight_enumerator(S; max_terms=10_000_000)

Compute or retrieve the Hamming-weight Shor--Laflamme pair `(A,B)`. `A` is
enumerated from the additive stabilizer group and `B` is obtained by the
trace-symplectic MacWilliams transform. The cache stores only
`HammingWeightEnumerator` coefficient data.

For subsystem codes this describes the full stabilized codespace; it is not a
bare/dressed protected-subsystem enumerator.
"""
function Shor_Laflamme_weight_enumerator(
    S::AbstractSubsystemCode; max_terms::Integer=10_000_000,
)
    return get!(S.cache, :SL_weight_enum) do
        A_counts = _hamming_distribution_additive(
            stabilizers(S), S.n, S.F; max_terms=max_terms)
        A = HammingWeightEnumerator(S.n, A_counts)
        B = _trace_symplectic_MacWilliams(
            A, cardinality(S), Int(order(S.F)))
        S.cache[:weight_enum_A] = A
        S.cache[:weight_enum_B] = B
        S.cache[:weight_dist_A] = A.counts
        S.cache[:weight_dist_B] = B.counts
        result = _validate_SL_enumerator(
            ShorLaflammeWeightEnumerator(S.n, A, B),
            cardinality(S), Int(order(S.F)))
        S isa AbstractStabilizerCode &&
            _cache_SL_distance!(S, result)
        result
    end
end

SL_weight_enumerator(S::AbstractSubsystemCode; kwargs...) =
    Shor_Laflamme_weight_enumerator(S; kwargs...)
shor_laflamme_weight_enumerator(S::AbstractSubsystemCode; kwargs...) =
    Shor_Laflamme_weight_enumerator(S; kwargs...)

function _SL_quotient_enumerator(
    enumerator::ShorLaflammeWeightEnumerator
)
    counts = Dict{Int, BigInt}()
    for w in 0:enumerator.n
        value = get(enumerator.B.counts, w, BigInt(0)) -
            get(enumerator.A.counts, w, BigInt(0))
        !iszero(value) && (counts[w] = value)
    end
    return HammingWeightEnumerator(enumerator.n, counts)
end

"""
    weight_enumerator(S; set=:all, max_terms=10_000_000)

Return Hamming-weight quantum enumerators. `set=:all` returns the
`ShorLaflammeWeightEnumerator`; `:stabilizers`/`:A`, `:normalizer`/`:B`, and
`:quotient` select one `HammingWeightEnumerator`.
"""
function weight_enumerator(
    S::AbstractSubsystemCode;
    set::Symbol=:all,
    type::Symbol=:Hamming,
    alg::Symbol=:auto,
    max_terms::Integer=10_000_000,
    verbose::Bool=false,
)
    type in (:Hamming, :hamming) ||
        throw(ArgumentError(
            "Only Hamming quantum enumerators are supported; received `$type`."))
    alg in (:auto, :bruteforce) ||
        throw(ArgumentError(
            "Quantum enumerator algorithm `$alg` is not implemented."))
    set in (:all, :stabilizers, :A, :normalizer, :B, :quotient) ||
        throw(ArgumentError("Unknown quantum enumerator set `$set`."))
    verbose && println("Computing Shor–Laflamme Hamming enumerators.")
    enumerator =
        Shor_Laflamme_weight_enumerator(S; max_terms=max_terms)
    set == :all && return enumerator
    set in (:stabilizers, :A) && return enumerator.A
    set in (:normalizer, :B) && return enumerator.B
    return _SL_quotient_enumerator(enumerator)
end

function weight_distribution(
    S::AbstractSubsystemCode;
    set::Symbol=:all,
    alg::Symbol=:auto,
    max_terms::Integer=10_000_000,
    verbose::Bool=false,
)
    enumerator = weight_enumerator(
        S; set=set, alg=alg, max_terms=max_terms, verbose=verbose)
    if enumerator isa ShorLaflammeWeightEnumerator
        return (A=copy(enumerator.A.counts), B=copy(enumerator.B.counts))
    end
    return copy(enumerator.counts)
end

function weight_distribution_array(
    S::AbstractSubsystemCode;
    set::Symbol=:stabilizers,
    max_terms::Integer=10_000_000,
)
    counts =
        weight_distribution(S; set=set, max_terms=max_terms)
    counts isa NamedTuple &&
        throw(ArgumentError("Select set=:A, :B, or :quotient for an array."))
    result = zeros(BigInt, S.n + 1)
    for (weight, count) in counts
        result[weight + 1] = count
    end
    return result
end
