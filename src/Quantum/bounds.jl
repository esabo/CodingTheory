# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

const _quantum_Krawtchouk_matrix_cache =
    Dict{Tuple{Int, Int, Bool}, Matrix{BigInt}}()
const _quantum_Krawtchouk_matrix_cache_lock = ReentrantLock()

function _quantum_bound_dimensions(
    n::Integer, k::Union{Integer, Rational},
    r::Union{Integer, Rational}=0,
    field_degree::Integer=1,
)
    n >= 1 || throw(DomainError(n, "The code length must be positive."))
    field_degree >= 1 ||
        throw(DomainError(field_degree, "The field degree must be positive."))
    k_rat = Rational{BigInt}(k)
    r_rat = Rational{BigInt}(r)
    k_rat > 0 ||
        throw(DomainError(k, "Quantum distance bounds require a positive protected dimension."))
    r_rat >= 0 || throw(DomainError(r, "The gauge dimension cannot be negative."))
    effective_r = r_rat / field_degree
    k_rat + effective_r <= n ||
        throw(DomainError(k + effective_r,
            "The protected and gauge dimensions cannot exceed n."))
    return Int(n), k_rat, effective_r
end

function _quantum_bound_field(S::AbstractSubsystemCode)
    hasfield(typeof(S), :F) && return getfield(S, :F)
    haskey(S.cache, :F) && return S.cache[:F]
    hasfield(typeof(S), :outer_code) &&
        return _quantum_bound_field(getfield(S, :outer_code))
    return _code_matrix_base_ring(stabilizers(S))
end

function _quantum_gauge_dimension(S::AbstractSubsystemCode)
    GaugeTrait(typeof(S)) == HasNoGauges() && return 0
    r = S.r
    ismissing(r) &&
        throw(ArgumentError("The gauge dimension must be known to apply this bound."))
    return r
end

function _prime_power_data(q::Integer)
    q >= 2 || throw(DomainError(q, "The alphabet size must be at least two."))
    factors = Nemo.factor(q)
    length(factors) == 1 ||
        throw(DomainError(q, "The alphabet size must be a prime power."))
    p, m = first(factors)
    return Int(p), Int(m)
end

function _protected_space_size(
    q::Integer, k::Union{Integer, Rational},
    r::Union{Integer, Rational}=0,
)
    p, m = _prime_power_data(q)
    # `k` is log_q(K), while library `r` counts prime-field gauge pairs:
    # K = p^(m*k), R = p^r.
    exponent = m * Rational{BigInt}(k) + Rational{BigInt}(r)
    denominator(exponent) == 1 ||
        throw(DomainError(k + r,
            "q^(k + r) must be an integer for a q-ary quantum code."))
    return BigInt(p)^Int(numerator(exponent))
end

"""
    quantum_Singleton_bound(n, k; r=0)
    quantum_Singleton_bound(S)
    Singleton_bound(S)

Return the quantum Singleton upper bound on distance. For a stabilizer code,
`k ≤ n - 2d + 2`; for a subsystem code,
`k + r/degree(F) ≤ n - 2d + 2`, because library `r` counts prime-field
gauge pairs.
The code must encode a positive-dimensional protected subsystem. For
subsystem code objects the bound is applied automatically only over prime
fields, or when Fq-linearity or purity is certified; use
`assume_applicable=true` to evaluate the formula outside those cases.
"""
function quantum_Singleton_bound(
    n::Integer, k::Union{Integer, Rational};
    r::Union{Integer, Rational}=0,
    field_degree::Integer=1,
)
    n_int, k_rat, effective_r =
        _quantum_bound_dimensions(n, k, r, field_degree)
    return min(n_int, floor(Int, (n_int - k_rat - effective_r + 2) / 2))
end

function _subsystem_singleton_is_proven(S::AbstractSubsystemCode)
    GaugeTrait(typeof(S)) == HasNoGauges() && return true
    degree(_quantum_bound_field(S)) == 1 && return true
    return get(S.cache, :Fq_linear, false) || get(S.cache, :pure, false)
end

function quantum_Singleton_bound(
    S::AbstractSubsystemCode; assume_applicable::Bool=false,
)
    assume_applicable || _subsystem_singleton_is_proven(S) ||
        throw(ArgumentError(
            "The subsystem Singleton bound is proved for prime-field, " *
            "Fq-linear, or pure subsystem codes; pass " *
            "assume_applicable=true to evaluate the formula explicitly."))
    return quantum_Singleton_bound(S.n, dimension(S);
        r=_quantum_gauge_dimension(S),
        field_degree=degree(_quantum_bound_field(S)))
end
Singleton_bound(S::AbstractSubsystemCode) = quantum_Singleton_bound(S)

"""
    is_quantum_MDS(S)
    is_MDS(S)

Return whether the known exact (dressed, for subsystem codes) distance meets
the quantum Singleton bound. Return `missing` when the exact distance is not
stored.
"""
function is_quantum_MDS(S::AbstractSubsystemCode)
    d = if GaugeTrait(typeof(S)) == HasGauges()
        get(S.cache, :d_dressed,
            haskey(S.cache, :dx_dressed) && haskey(S.cache, :dz_dressed) ?
                min(S.cache[:dx_dressed], S.cache[:dz_dressed]) : missing)
    else
        get(S.cache, :d,
            haskey(S.cache, :dx) && haskey(S.cache, :dz) ?
                min(S.cache[:dx], S.cache[:dz]) :
                (hasfield(typeof(S), :d) ? getfield(S, :d) : missing))
    end
    ismissing(d) && return missing
    return d == quantum_Singleton_bound(S)
end
is_MDS(S::AbstractSubsystemCode) = is_quantum_MDS(S)

"""
    quantum_Hamming_volume(n, t, q)

Return the number of q-ary Pauli errors of weight at most `t`.
"""
function quantum_Hamming_volume(n::Integer, t::Integer, q::Integer)
    n >= 1 || throw(DomainError(n, "The code length must be positive."))
    0 <= t <= n ||
        throw(DomainError(t, "The error radius must lie between zero and n."))
    _prime_power_data(q)
    q2_minus_one = BigInt(q)^2 - 1
    return sum(
        binomial(BigInt(n), j) * q2_minus_one^j for j in 0:Int(t);
        init=BigInt(0))
end

function _quantum_hamming_capacity(
    n::Integer, k::Union{Integer, Rational}, q::Integer,
    r::Union{Integer, Rational},
)
    _, m = _prime_power_data(q)
    n_int, k_rat, _ = _quantum_bound_dimensions(n, k, r, m)
    protected_and_gauge = _protected_space_size(q, k_rat, r)
    total = BigInt(q)^n_int
    total % protected_and_gauge == 0 ||
        throw(DomainError(k + r, "q^(n - k - r) must be an integer."))
    return n_int, div(total, protected_and_gauge)
end

"""
    satisfies_quantum_Hamming_bound(n, k, d, q; r=0)

Test the sphere-packing inequality for a pure q-ary `[[n,k,d]]` stabilizer
code or pure `[[n,k,r,d]]` subsystem code. This is not valid for an impure
code.
"""
function satisfies_quantum_Hamming_bound(
    n::Integer, k::Union{Integer, Rational}, d::Integer, q::Integer;
    r::Union{Integer, Rational}=0,
)
    n_int, capacity = _quantum_hamming_capacity(n, k, q, r)
    1 <= d <= n_int ||
        throw(DomainError(d, "The distance must lie between one and n."))
    t = (Int(d) - 1) ÷ 2
    return quantum_Hamming_volume(n_int, t, q) <= capacity
end

"""
    quantum_Hamming_bound(n, k, q; r=0)
    quantum_Hamming_bound(S; assume_pure=false)

Return the largest distance not excluded by the pure quantum Hamming
(sphere-packing) bound. The code method requires cached purity or an explicit
`assume_pure=true`, because impure codes need not obey this bound.
"""
function quantum_Hamming_bound(
    n::Integer, k::Union{Integer, Rational}, q::Integer;
    r::Union{Integer, Rational}=0,
)
    n_int, capacity = _quantum_hamming_capacity(n, k, q, r)
    q2_minus_one = BigInt(q)^2 - 1
    term = BigInt(1)
    volume = BigInt(1)
    radius = 0
    for t in 1:n_int
        term = div(term * (n_int - t + 1) * q2_minus_one, t)
        volume += term
        volume <= capacity || break
        radius = t
    end
    return min(n_int, 2radius + 2)
end

function quantum_Hamming_bound(
    S::AbstractSubsystemCode; assume_pure::Bool=false,
)
    cached_purity = get(S.cache, :pure, missing)
    assume_pure || cached_purity === true ||
        throw(ArgumentError(
            "The quantum Hamming bound applies only to pure codes; " *
            "cache purity first or pass assume_pure=true."))
    q = Int(order(_quantum_bound_field(S)))
    return quantum_Hamming_bound(
        S.n, dimension(S), q; r=_quantum_gauge_dimension(S))
end

"""
    quantum_Gilbert_Varshamov_exists(n, k, d, q; r=0, variant=:additive)
    quantum_Gilbert_Varshamov_bound(n, k, q; r=0, variant=:additive)
    quantum_Gilbert_Varshamov_bound(S)

Test, or return the largest `d` guaranteed by, a finite quantum
Gilbert--Varshamov existence bound. `variant=:additive` supports additive
stabilizer and subsystem parameters. `:linear` is the Fq²-linear Ketkar bound,
and `:pure_feng_ma` is the pure linear Feng--Ma bound. These are parameter
benchmarks, not certified lower bounds for a concrete code, and therefore
never change its distance cache.
"""
function quantum_Gilbert_Varshamov_exists(
    n::Integer, k::Union{Integer, Rational}, d::Integer, q::Integer;
    r::Union{Integer, Rational}=0, variant::Symbol=:additive,
)
    p, m = _prime_power_data(q)
    n_int, k_rat, _ = _quantum_bound_dimensions(n, k, r, m)
    1 <= d <= n_int ||
        throw(DomainError(d, "The target distance must lie between one and n."))
    d == 1 && return true
    variant in (:additive, :linear, :pure_feng_ma) ||
        throw(ArgumentError(
            "Expected variant=:additive, :linear, or :pure_feng_ma."))

    if variant != :additive
        iszero(r) ||
            throw(ArgumentError("$variant is a stabilizer, not subsystem, bound."))
        denominator(k_rat) == 1 ||
            throw(ArgumentError("$variant requires integral k."))
        k_int = Int(k_rat)
        k_int >= (variant == :pure_feng_ma ? 2 : 1) ||
            throw(DomainError(k,
                "$variant requires k ≥ $(variant == :pure_feng_ma ? 2 : 1)."))
        iseven(n_int - k_int) ||
            throw(ArgumentError("$variant requires n ≡ k (mod 2)."))
        q2_minus_one = BigInt(q)^2 - 1
        reduced_volume = sum(
            binomial(BigInt(n_int), j) * q2_minus_one^(j - 1)
            for j in 1:Int(d - 1); init=BigInt(0))
        if variant == :linear
            multiplier =
                BigInt(q)^(n_int + k_int) - BigInt(q)^(n_int - k_int)
            return multiplier * reduced_volume < BigInt(q)^(2n_int) - 1
        end
        numerator_rhs = BigInt(q)^(n_int - k_int + 2) - 1
        numerator_rhs % q2_minus_one == 0 ||
            error("The Feng--Ma right-hand side was not integral.")
        return reduced_volume < div(numerator_rhs, q2_minus_one)
    end

    K = _protected_space_size(q, k_rat)
    R = _protected_space_size(q, 0, r)
    qn = BigInt(q)^n_int
    qn * R % K == 0 ||
        throw(DomainError(k, "q^n R/K must be an integer."))

    multiplier = qn * K * R - div(qn * R, K)
    right = (BigInt(q)^(2n_int) - 1) * (p - 1)
    volume = sum(
        binomial(BigInt(n_int), j) * (BigInt(q)^2 - 1)^j
        for j in 1:Int(d - 1); init=BigInt(0))
    return multiplier * volume < right
end

function quantum_Gilbert_Varshamov_bound(
    n::Integer, k::Union{Integer, Rational}, q::Integer;
    r::Union{Integer, Rational}=0, variant::Symbol=:additive,
)
    p, m = _prime_power_data(q)
    n_int, _, _ = _quantum_bound_dimensions(n, k, r, m)
    best = 1
    for d in 2:n_int
        quantum_Gilbert_Varshamov_exists(
            n, k, d, q; r=r, variant=variant) || break
        best = d
    end
    return best
end

function quantum_Gilbert_Varshamov_bound(S::AbstractSubsystemCode)
    q = Int(order(_quantum_bound_field(S)))
    return quantum_Gilbert_Varshamov_bound(
        S.n, dimension(S), q; r=_quantum_gauge_dimension(S))
end

"""
    quantum_stabilizer_generator_weight_lower_bound(n, k; min_distance=2)
    quantum_stabilizer_generator_weight_lower_bound(S; min_distance=2)

Return `ceil(2n/(n-k))`, the finite lower bound of Wei et al. on
the optimal maximum generator weight of a binary stabilizer code, strengthened
to six by their Proposition 20 at `(n,k,d)=(12,7,2)`. The theorem requires
`k ≥ 1` and distance at least two.
"""
function quantum_stabilizer_generator_weight_lower_bound(
    n::Integer, k::Integer; min_distance::Integer=2,
)
    n >= 1 || throw(DomainError(n, "The code length must be positive."))
    1 <= k < n ||
        throw(DomainError(k, "This bound requires 1 ≤ k < n."))
    min_distance >= 2 ||
        throw(DomainError(min_distance, "This bound requires distance at least two."))
    n == 12 && k == 7 && min_distance == 2 && return BigInt(6)
    return cld(2 * BigInt(n), BigInt(n - k))
end

function quantum_stabilizer_generator_weight_lower_bound(
    S::AbstractStabilizerCode; min_distance::Integer=2,
)
    Int(order(_quantum_bound_field(S))) == 2 ||
        throw(ArgumentError("The low-check-weight theorem is binary."))
    isinteger(dimension(S)) ||
        throw(ArgumentError("The low-check-weight theorem requires integral k."))
    return quantum_stabilizer_generator_weight_lower_bound(
        S.n, Int(dimension(S)); min_distance=min_distance)
end

"""
    quantum_check_weight_dimension_bound(n, check_weight; min_distance=2)

Return the corresponding upper bound
`k ≤ n - ceil(2n/check_weight)` for a binary stabilizer code of distance at
least two.
"""
function quantum_check_weight_dimension_bound(
    n::Integer, check_weight::Integer; min_distance::Integer=2,
)
    n >= 1 || throw(DomainError(n, "The code length must be positive."))
    check_weight >= 1 ||
        throw(DomainError(check_weight, "The check weight must be positive."))
    min_distance >= 2 ||
        throw(DomainError(min_distance, "This bound requires distance at least two."))
    return max(BigInt(0), BigInt(n) - cld(2 * BigInt(n), BigInt(check_weight)))
end

"""
    quantum_low_weight_stabilizer_distance_bound(n, k, check_weight)
    quantum_low_weight_stabilizer_distance_bound(S)

Return the finite distance upper bound for a binary stabilizer presentation
whose generators all have weight at most three. A weight-at-most-two
presentation has `d ≤ 1`; at weight three, `d ≤ 2`, strengthened to `d ≤ 1`
when `k/n > 1/4`. Return `missing` above weight three.
"""
function quantum_low_weight_stabilizer_distance_bound(
    n::Integer, k::Integer, check_weight::Integer,
)
    n >= 1 || throw(DomainError(n, "The code length must be positive."))
    1 <= k <= n ||
        throw(DomainError(k, "This bound requires 1 ≤ k ≤ n."))
    check_weight >= 0 ||
        throw(DomainError(check_weight, "The check weight cannot be negative."))
    check_weight <= 2 && return 1
    check_weight == 3 && return 4 * BigInt(k) > n ? 1 : 2
    return missing
end

function quantum_low_weight_stabilizer_distance_bound(
    S::AbstractStabilizerCode,
)
    Int(order(_quantum_bound_field(S))) == 2 ||
        throw(ArgumentError("The low-check-weight theorem is binary."))
    isinteger(dimension(S)) ||
        throw(ArgumentError("The low-check-weight theorem requires integral k."))
    return quantum_low_weight_stabilizer_distance_bound(
        S.n, Int(dimension(S)), maximum_stabilizer_weight(S))
end

"""
    quantum_CSS_subsystem_weight_two_distance_bound(n, k)
    satisfies_quantum_CSS_subsystem_weight_two_bounds(n, k, d_X, d_Z)

Return `min(floor(sqrt(n)), floor(n/k))`, the dressed-distance upper bound
for a binary CSS subsystem code presented by gauge checks of weight at most
two. This parameter API does not inspect a code object because CodingTheory
does not yet retain the original gauge-check presentation required by the
theorem. The predicate preserves the stronger asymmetric statements
`d_X*d_Z ≤ n`, `k*d_X ≤ n`, and `k*d_Z ≤ n`.
"""
function quantum_CSS_subsystem_weight_two_distance_bound(
    n::Integer, k::Integer,
)
    n >= 1 || throw(DomainError(n, "The code length must be positive."))
    1 <= k <= n ||
        throw(DomainError(k, "This bound requires 1 ≤ k ≤ n."))
    return min(isqrt(BigInt(n)), fld(BigInt(n), BigInt(k)))
end

function satisfies_quantum_CSS_subsystem_weight_two_bounds(
    n::Integer, k::Integer, d_X::Integer, d_Z::Integer,
)
    n >= 1 || throw(DomainError(n, "The code length must be positive."))
    1 <= k <= n ||
        throw(DomainError(k, "This bound requires 1 ≤ k ≤ n."))
    d_X >= 1 && d_Z >= 1 ||
        throw(DomainError((d_X, d_Z), "The X and Z distances must be positive."))
    return BigInt(d_X) * d_Z <= n &&
        BigInt(k) * d_X <= n &&
        BigInt(k) * d_Z <= n
end

"""
    quantum_stabilizer_check_weight_existence_bound(n, k, d)

Return a constructive upper bound on the optimal maximum generator weight:
three when `d=2` and `n ≥ 4k`, or four when `d≥3` and `n ≥ k*d^2`.
Return `missing` when these constructions do not apply.
"""
function quantum_stabilizer_check_weight_existence_bound(
    n::Integer, k::Integer, d::Integer,
)
    n >= 1 || throw(DomainError(n, "The code length must be positive."))
    1 <= k <= n ||
        throw(DomainError(k, "This bound requires 1 ≤ k ≤ n."))
    d >= 2 || throw(DomainError(d, "This bound requires d ≥ 2."))
    d == 2 && return n >= 4 * BigInt(k) ? BigInt(3) : missing
    return BigInt(n) >= BigInt(k) * BigInt(d)^2 ? BigInt(4) : missing
end

"""
    quantum_stabilizer_group_average_weight(n, A₁=0)
    quantum_stabilizer_group_total_weight(n, k, A₁=0)

Evaluate Wei et al.'s exact stabilizer-group identity for distance at least
two:
`sum(j*A_j)/sum(A_j) = (3n-A₁)/4`. The total-weight form uses
`sum(A_j)=2^(n-k)` and returns exact `Rational{BigInt}` arithmetic.
"""
function quantum_stabilizer_group_average_weight(
    n::Integer, A₁::Integer=0,
)
    n >= 1 || throw(DomainError(n, "The code length must be positive."))
    0 <= A₁ <= n ||
        throw(DomainError(A₁, "A₁ must lie between zero and n."))
    return (3 * BigInt(n) - A₁) // 4
end

function quantum_stabilizer_group_total_weight(
    n::Integer, k::Integer, A₁::Integer=0,
)
    0 <= k <= n || throw(DomainError(k, "The dimension must lie in 0:n."))
    return BigInt(2)^(n - k) *
        quantum_stabilizer_group_average_weight(n, A₁)
end

"""
    quantum_Krawtchouk_matrix(n; alphabet_size=4, signed_columns=false)

Return an exact Krawtchouk matrix
`M[i + 1, j + 1] = P_i(j; n)` as `Matrix{BigInt}`. With
`signed_columns=true`, column `j` is multiplied by `(-1)^j`, giving the
signed matrix used in the Shor--Laflamme LP constraints. Use
`alphabet_size=2` for binary CSS constituent codes and `4` for binary Pauli
weight enumerators.
"""
function quantum_Krawtchouk_matrix(
    n::Integer; alphabet_size::Integer=4, signed_columns::Bool=false,
)
    n >= 0 || throw(DomainError(n, "The code length cannot be negative."))
    alphabet_size >= 2 ||
        throw(DomainError(alphabet_size, "The alphabet size must be at least two."))
    n_int = Int(n)
    q = Int(alphabet_size)
    key = (n_int, q, signed_columns)
    cached = lock(_quantum_Krawtchouk_matrix_cache_lock) do
        get(_quantum_Krawtchouk_matrix_cache, key, nothing)
    end
    isnothing(cached) || return copy(cached)

    M = Matrix{BigInt}(undef, n_int + 1, n_int + 1)
    for j in 0:n_int
        M[1, j + 1] = 1
        n_int == 0 && continue
        M[2, j + 1] = BigInt(q - 1) * n_int - q * j
        for ell in 1:(n_int - 1)
            numerator =
                ((q - 1) * (n_int - ell) + ell - q * j) *
                    M[ell + 1, j + 1] -
                (q - 1) * (n_int - ell + 1) * M[ell, j + 1]
            value, remainder = divrem(numerator, ell + 1)
            iszero(remainder) ||
                error("The Krawtchouk recurrence produced a nonintegral value.")
            M[ell + 2, j + 1] = value
        end
        signed_columns && isodd(j) && (M[:, j + 1] .*= -1)
    end
    stored = lock(_quantum_Krawtchouk_matrix_cache_lock) do
        get!(_quantum_Krawtchouk_matrix_cache, key, M)
    end
    return copy(stored)
end

function _quantum_CSS_weight_enumerator_LP(
    args...; kwargs...
)
    throw(ArgumentError(
        "Load both JuMP and Tulip to use quantum weight-enumerator LP bounds."))
end

"""
    quantum_CSS_weight_enumerator_LP(n, k_X, k_Z, d; kwargs...)

Solve the exact-coefficient CSS split-enumerator LP of Wang et al. Here
`C_X` and `C_Z` have dimensions `k_X` and `k_Z`, so the quantum dimension is
`k_X + k_Z - n`. Pass `check_weight=w` to add the paper's cumulative
low-check-weight constraints, and `exclude_weight_one=true` for its second
post-processing branch. `model_hook(model, A)` supports script-level
extensions; entries `A[1:n+1]` and `A[n+2:2n+2]` are the X and Z dual-code
enumerators, respectively.
"""
function quantum_CSS_weight_enumerator_LP(
    n::Integer, k_X::Integer, k_Z::Integer, d::Integer; kwargs...
)
    n >= 1 || throw(DomainError(n, "The code length must be positive."))
    0 <= k_X <= n && 0 <= k_Z <= n ||
        throw(DomainError((k_X, k_Z), "CSS constituent dimensions must lie in 0:n."))
    k_X + k_Z > n ||
        throw(DomainError((k_X, k_Z), "The CSS code must encode at least one qubit."))
    1 <= d <= n || throw(DomainError(d, "The distance must lie in 1:n."))
    return _quantum_CSS_weight_enumerator_LP(
        Int(n), Int(k_X), Int(k_Z), Int(d); kwargs...)
end

function _quantum_stabilizer_dimension_LP_bound(
    args...; kwargs...
)
    throw(ArgumentError(
        "Load both JuMP and Tulip to use quantum weight-enumerator LP bounds."))
end

"""
    quantum_stabilizer_dimension_LP_bound(n, d, check_weight; kwargs...)

Search the Wang et al. general-stabilizer LP for the largest feasible quantum
dimension `k`. The result retains every trial and distinguishes numerical
infeasibility from an exact certificate.
"""
function quantum_stabilizer_dimension_LP_bound(
    n::Integer, d::Integer, check_weight::Integer; kwargs...
)
    n >= 1 || throw(DomainError(n, "The code length must be positive."))
    1 <= d <= n || throw(DomainError(d, "The distance must lie in 1:n."))
    1 <= check_weight <= n ||
        throw(DomainError(check_weight, "The check weight must lie in 1:n."))
    return _quantum_stabilizer_dimension_LP_bound(
        Int(n), Int(d), Int(check_weight); kwargs...)
end

function _quantum_CSS_dimension_LP_bound(
    args...; kwargs...
)
    throw(ArgumentError(
        "Load both JuMP and Tulip to use quantum weight-enumerator LP bounds."))
end

"""
    quantum_CSS_dimension_LP_bound(n, d, check_weight; kwargs...)

Search all CSS constituent-dimension splits in the Wang et al. LP and return
the largest feasible quantum dimension. `exclude_weight_one` selects the
second branch used in the paper's monotonic post-processing.
"""
function quantum_CSS_dimension_LP_bound(
    n::Integer, d::Integer, check_weight::Integer; kwargs...
)
    n >= 1 || throw(DomainError(n, "The code length must be positive."))
    1 <= d <= n || throw(DomainError(d, "The distance must lie in 1:n."))
    1 <= check_weight <= n ||
        throw(DomainError(check_weight, "The check weight must lie in 1:n."))
    return _quantum_CSS_dimension_LP_bound(
        Int(n), Int(d), Int(check_weight); kwargs...)
end

"""
    quantum_check_weight_LP_postprocess(with_weight_one, without_weight_one)

Apply Wang et al. Equation (18) to finite raw LP tables. Dictionary keys are
`(n,d,w)` and values are maximum feasible `k`. The first table permits
weight-one checks; the second imposes `A₁=0` (on both CSS enumerators when
applicable). The available keys define the finite `n′` and `w′` horizons.
Only points for which both intermediate closures are defined are returned.
"""
function quantum_check_weight_LP_postprocess(
    with_weight_one::AbstractDict{<:NTuple{3, Int}, <:Integer},
    without_weight_one::AbstractDict{<:NTuple{3, Int}, <:Integer},
)
    points = union(keys(with_weight_one), keys(without_weight_one))
    k1 = Dict{NTuple{3, Int}, Int}()
    k2 = Dict{NTuple{3, Int}, Int}()
    for point in points
        n, d, w = point
        upper_n = Int[
            value for ((n_prime, d_prime, w_prime), value) in with_weight_one
            if n_prime >= n && d_prime == d && w_prime == w
        ]
        lower_n = Int[
            value for ((n_prime, d_prime, w_prime), value) in without_weight_one
            if n_prime <= n && d_prime == d && w_prime == w
        ]
        isempty(upper_n) || (k1[point] = minimum(upper_n))
        isempty(lower_n) || (k2[point] = maximum(lower_n))
    end

    result = Dict{NTuple{3, Int}, Int}()
    candidates = intersect(keys(k1), keys(k2))
    for point in candidates
        n, d, w = point
        closure = Int[
            min(k1[other], k2[other])
            for other in candidates
            if other[1] == n && other[2] <= d && other[3] >= w
        ]
        isempty(closure) || (result[point] = minimum(closure))
    end
    return result
end

"""
    QuantumLPResult

Result of an arbitrary-precision quantum weight-enumerator LP. `status` is
`:feasible`, `:infeasible_numerical`, or `:unknown`; numerical infeasibility
is deliberately not presented as an exact certificate. Coefficients are
constructed exactly as `BigInt`, independently normalized by constraint row,
and converted to `BigFloat` only at the optimizer boundary.
"""
struct QuantumLPResult
    status::Symbol
    formulation::Symbol
    precision::Int
    termination_status::Symbol
    enumerator::Union{Nothing, Vector{BigFloat}}
    max_normalized_violation::BigFloat
end

function _quantum_weight_enumerator_LP(
    args...; kwargs...
)
    throw(ArgumentError(
        "Load both JuMP and Tulip to use quantum weight-enumerator LP bounds."))
end

"""
    quantum_weight_enumerator_LP(n, k, d; formulation=:standard, kwargs...)

Solve the binary Shor--Laflamme/Rains weight-enumerator feasibility LP.
`formulation=:standard` uses the MacWilliams and shadow constraints.
`:coarse` additionally imposes the cumulative check-growth inequalities for
`check_weight`. `:refined` adds the low-weight-generator constraints of Wei
et al. and requires `check_weight` and `num_max_weight_generators`.
Set `include_shadow=false` to reproduce Wang et al.'s general-stabilizer LP
without the additional Rains shadow inequalities.

The optimizer runs in arbitrary precision through Tulip. Increase
`precisions`, for example to `(256, 512, 1024)`, when the returned status is
`:unknown` or when numerical infeasibility needs stronger evidence.
For family-specific experiments, `model_hook(model, A)` may add JuMP
variables and constraints after the exact standard rows are installed. Such
custom rows are enforced by the optimizer but are not included in
`max_normalized_violation`.
`cumulative_lower_bounds` accepts `cutoff => count` pairs imposing
`sum(A[0:cutoff]) ≥ count`; this is the architecture/support-union interface
used in Wei et al.'s geometry-aware LP.
"""
function quantum_weight_enumerator_LP(
    n::Integer, k::Integer, d::Integer; kwargs...
)
    n >= 1 || throw(DomainError(n, "The code length must be positive."))
    1 <= k < n || throw(DomainError(k, "This LP requires 1 ≤ k < n."))
    1 <= d <= n || throw(DomainError(d, "The distance must lie in 1:n."))
    return _quantum_weight_enumerator_LP(
        Int(n), Int(k), Int(d); kwargs...)
end

function _quantum_stabilizer_generator_weight_LP_bound(
    args...; kwargs...
)
    throw(ArgumentError(
        "Load both JuMP and Tulip to use quantum weight-enumerator LP bounds."))
end

"""
    quantum_stabilizer_generator_weight_LP_bound(n, k, d; kwargs...)

Search the refined Wei et al. LPs for a lower bound on the maximum generator
weight of a binary `[[n,k,d]]` stabilizer code. The returned named tuple
records whether excluded weights are numerical LP conclusions or exact
certificates. `connected=true` adds the paper's connected-overlap inequality;
it asserts that tensor-factor decompositions have already been excluded.
"""
function quantum_stabilizer_generator_weight_LP_bound(
    n::Integer, k::Integer, d::Integer; kwargs...
)
    n >= 1 || throw(DomainError(n, "The code length must be positive."))
    1 <= k < n || throw(DomainError(k, "This LP requires 1 ≤ k < n."))
    2 <= d <= n || throw(DomainError(d, "This LP requires 2 ≤ d ≤ n."))
    return _quantum_stabilizer_generator_weight_LP_bound(
        Int(n), Int(k), Int(d); kwargs...)
end

function _seed_quantum_singleton_bound!(S::AbstractSubsystemCode)
    dimension(S) > 0 || return S
    _subsystem_singleton_is_proven(S) || return S
    bound = quantum_Singleton_bound(S)
    S.cache[:u_bound] = min(get(S.cache, :u_bound, S.n), bound)
    if GaugeTrait(typeof(S)) == HasGauges()
        S.cache[:u_bound_dressed] =
            min(get(S.cache, :u_bound_dressed, S.n), bound)
    end
    if S isa AbstractStabilizerCode &&
       Int(order(_quantum_bound_field(S))) == 2 &&
       isinteger(dimension(S)) && dimension(S) < S.n
        low_weight_bound = quantum_low_weight_stabilizer_distance_bound(S)
        ismissing(low_weight_bound) ||
            (S.cache[:u_bound] = min(S.cache[:u_bound], low_weight_bound))
    end
    return S
end
