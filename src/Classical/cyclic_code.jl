# Copyright (c) 2021 - 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################
"""
$(TYPEDSIGNATURES)

Return the generator matrix of the cyclic code.
Evaluates lazily by extracting the coefficients of the generator polynomial 
`C.g` and mapping them into a strictly typed dense matrix (`fpMatrix` or `FqMatrix`).
"""
function generator_matrix(C::AbstractCyclicCode, stand_form::Bool = false)
    cache = getfield(C, :cache)
    
    if !haskey(cache, :G)
        n, k = C.n, C.k
        deg_g = n - k
        
        # Enforce modern native matrix types (fpMatrix for prime, FqMatrix for extension)
        p_char = Int(characteristic(C.F))
        q_order = Int(order(C.F))
        # modern_F = (q_order == p_char) ? Oscar.Nemo.Native.GF(p_char) : GF(p_char, degree(C.F), :α)
        
        G_mat = zero_matrix(C.F, k, n)
        for i in 1:k
            for j in 0:deg_g
                c_val = coeff(C.g, j)
                # Safe casting from polynomial coefficients to native field elements
                G_mat[i, i + j] = C.F(c_val)
            end
        end
        
        cache[:G] = G_mat
    end
    
    if stand_form
        if !haskey(cache, :G_stand)
            G_stand, H_stand, P_stand, _ = _standard_form(cache[:G])
            cache[:G_stand] = G_stand
            cache[:H_stand] = H_stand
            cache[:P_stand] = P_stand
        end
        return cache[:G_stand]
    end
    
    return cache[:G]
end

"""
$(TYPEDSIGNATURES)

Return the parity-check matrix of the cyclic code.
Evaluates lazily by extracting the reversed coefficients of the parity polynomial `C.h`.
"""
function parity_check_matrix(C::AbstractCyclicCode, stand_form::Bool = false)
    cache = getfield(C, :cache)
    
    if !haskey(cache, :H)
        n, k = C.n, C.k
        deg_h = k
        
        # Enforce modern native matrix types
        p_char = Int(characteristic(C.F))
        q_order = Int(order(C.F))
        # modern_F = (q_order == p_char) ? Oscar.Nemo.Native.GF(p_char) : GF(p_char, degree(C.F), :α)
        
        H_mat = zero_matrix(C.F, n - k, n)
        for i in 1:(n - k)
            for j in 0:deg_h
                # The parity check rows are formed by the reversed coefficients of h(x)
                c_val = coeff(C.h, deg_h - j)
                H_mat[i, i + j] = C.F(c_val)
            end
        end
        
        cache[:H] = H_mat
    end
    
    if stand_form
        if !haskey(cache, :H_stand)
            # Forcing the standard form generator implicitly populates H_stand in the cache
            generator_matrix(C, true)
        end
        return cache[:H_stand]
    end
    
    return cache[:H]
end

function _cyclic_algebra_from_def_set(q::Int, n::Int, def_set::Vector{Int}, cosets::Vector{Vector{Int}})
    factors = Nemo.factor(q)
    length(factors) == 1 || throw(DomainError(q, "There is no finite field of order $q."))
    (p, t), = factors

    # F = Oscar.Nemo.Native.GF(p, t, :α)
    F = GF(p, t, :α)
    deg = ord(n, q)
    # E = Oscar.Nemo.Native.GF(p, t * deg, :α)
    E = GF(p, t * deg, :α)
    α = (t * deg == 1) ? E(2) : gen(E)
    
    R, x = polynomial_ring(E, :x)
    β = α^(div(BigInt(q)^deg - 1, n))

    k = n - length(def_set)
    com_cosets = complement_qcosets(q, n, cosets)
    
    g = _generator_polynomial(R, β, def_set)
    h = _generator_polynomial(R, β, reduce(vcat, com_cosets))
    e = _idempotent(g, h, n)
    
    # Compute bounds eagerly only because we need them to auto-detect BCH/RS
    δ, b, runs = _BCH_bound_math(def_set, n)
    
    coeffs = collect(coefficients(g))
    upper = count(!iszero, coeffs)
    
    return F, E, R, β, k, g, h, e, δ, b, upper, runs
end

"""
$(TYPEDSIGNATURES)

Return the CyclicCode of length `n` over `GF(q)` with `q`-cyclotomic cosets `cosets`.
Evaluates lazily and auto-detects BCH and Reed-Solomon parameters.
"""
function CyclicCode(q::Int, n::Int, cosets::Vector{Vector{Int}})
    (q <= 1 || n <= 1) && throw(DomainError((q, n), "Invalid parameters passed to CyclicCode constructor."))
    
    def_set = isempty(cosets) ? Int[] : sort!(reduce(vcat, cosets))
    qcosets_reps = sort!([arr[1] for arr in cosets])
    
    F, E, R, β, k, g, h, e, δ, b, upper, runs = _cyclic_algebra_from_def_set(q, n, def_set, cosets)

    cache = Dict{Symbol, Any}(
        :BCH_bound => δ,
        :bch_offset => b,
        :bch_runs => runs,
        :d => missing
    )
    
    if δ >= 2 && def_set == defining_set(collect(b:(b + δ - 2)), q, n, true)
        deg = ord(n, q)
        if deg == 1 && n == q - 1
            d = n - k + 1
            cache[:HT_bound] = d
            cache[:d] = d
            cache[:designed_b] = b
            cache[:designed_delta] = δ
            
            counts = Dict{Int, BigInt}(0 => 1)
            for w in d:n
                sum_val = BigInt(0)
                for j in 0:(w - d)
                    term = ((-1)^j) * binomial(w - 1, j) * (BigInt(q)^(w - d - j))
                    sum_val += term
                end
                counts[w] = binomial(n, w) * (q - 1) * sum_val
            end
            cache[:weight_enum] = HammingWeightEnumerator(n, counts)
            cache[:d] = d
            
            # Struct strictly expects 15 args: ..., n, k, b, d, ...
            return ReedSolomonCode(F, E, R, β, n, k, b, d, d, d, cosets, qcosets_reps, def_set, g, h, e, cache)
        end

        cache[:designed_b] = b
        cache[:designed_delta] = δ
        # Struct strictly expects 15 args: ..., n, k, b, δ, ...
        return BCHCode(F, E, R, β, n, k, b, δ, 1, n - k + 1, cosets, qcosets_reps, def_set, g, h, e, cache)
    end

    # Standard CyclicCode strictly expects 15 args: ..., n, k, l_bound, u_bound, ...
    return CyclicCode(F, E, R, β, n, k, 1, upper, cosets, qcosets_reps, def_set, g, h, e, cache)
end

"""
$(TYPEDSIGNATURES)

Return the length `n` cyclic code generated by the polynomial `g`.
"""
function CyclicCode(n::Int, g::Union{fpPolyRingElem, FqPolyRingElem, fqPolyRepPolyRingElem})
    n > 0 || throw(DomainError(n, "Invalid parameters passed to CyclicCode constructor."))
    R = parent(g)
    flag, _ = divides(gen(R)^n - 1, g)
    flag || throw(ArgumentError("Given polynomial does not divide x^$n - 1."))

    F = base_ring(R)
    q = Int(order(F))
    deg = ord(n, q)
    p = Int(characteristic(F))
    t = Int(degree(F))
    
    E = GF(p, t * deg, :α)
    α = (t * deg == 1) ? E(2) : gen(E)
    β = α^(div(q^deg - 1, n))
    
    R_E, _ = polynomial_ring(E, :y)
    g_E = (t == 1 && typeof(g) == fpPolyRingElem) ? R_E(E.(lift.(Ref(ZZ), collect(coefficients(g))))) : R_E([E(i) for i in collect(coefficients(g))])

    # Build a rapid lookup table for the powers of β
    beta_powers = Dict{typeof(β), Int}()
    curr_power = one(E)
    for i in 0:(n - 1)
        beta_powers[curr_power] = i
        curr_power *= β
    end

    # Look up each root's exponent directly
    rt_indices = [beta_powers[rt] for rt in roots(g_E)]
    cosets = defining_set(sort!(rt_indices), q, n, false)
    
    return CyclicCode(q, n, cosets)
end

"""
$(TYPEDSIGNATURES)

Return the BCHCode of length `n` over `GF(q)` with design distance `δ` and designed offset `b`.
"""
function BCHCode(q::Int, n::Int, δ::Int, b::Int = 0)
    δ >= 2 || throw(DomainError(δ, "BCH codes require δ ≥ 2."))

    if n == q - 1
        return ReedSolomonCode(q, δ, b)
    end
    
    # 1. Define the sets exactly as designed
    cosets = defining_set(collect(b:(b + δ - 2)), q, n, false)
    def_set = isempty(cosets) ? Int[] : sort!(reduce(vcat, cosets))
    qcosets_reps = sort!([arr[1] for arr in cosets])
    
    # 2. Drop to the internal engine for heavy algebra
    F, E, R, β, k, g, h, e, true_delta, true_b, upper, runs = _cyclic_algebra_from_def_set(q, n, def_set, cosets)

    # 3. Prime the cache with the TRUE mathematical bounds
    cache = Dict{Symbol, Any}(
        :BCH_bound => true_delta,
        :bch_offset => true_b,
        :bch_runs => runs,
        :d => missing
    )

    # 4. Construct the struct explicitly, preserving the DESIGNED b and δ
    return BCHCode(F, E, R, β, n, k, b, δ, δ, n - k + 1, cosets, qcosets_reps, def_set, g, h, e, cache)
end

"""
$(TYPEDSIGNATURES)

Return the ReedSolomonCode over `GF(q)` with distance `d` and designed offset `b`.
"""
function ReedSolomonCode(q::Int, d::Int, b::Int = 0)
    d >= 2 || throw(DomainError(d, "Reed Solomon codes require d ≥ 2."))
    q > 4 || throw(DomainError(q, "Invalid or too small parameters passed to ReedSolomonCode constructor."))
    
    n = q - 1
    
    # 1. Define the sets exactly as designed
    cosets = defining_set(collect(b:(b + d - 2)), q, n, false)
    def_set = isempty(cosets) ? Int[] : sort!(reduce(vcat, cosets))
    qcosets_reps = sort!([arr[1] for arr in cosets])
    
    # 2. Drop to the internal engine
    F, E, R, β, k, g, h, e, true_delta, true_b, upper, runs = _cyclic_algebra_from_def_set(q, n, def_set, cosets)

    # 3. Prime the cache (true_delta should equal d here, but we cache it for safety)
    cache = Dict{Symbol, Any}(
        :BCH_bound => true_delta,
        :bch_offset => true_b,
        :bch_runs => runs,
        :HT_bound => d,
        :d => d
    )
    
    # 4. Inject the mathematically known MDS Weight Enumerator
    counts = Dict{Int, BigInt}(0 => 1)
    for w in d:n
        sum_val = BigInt(0)
        for j in 0:(w - d)
            term = ((-1)^j) * binomial(w - 1, j) * (BigInt(q)^(w - d - j))
            sum_val += term
        end
        counts[w] = binomial(n, w) * (q - 1) * sum_val
    end
    cache[:weight_enum] = HammingWeightEnumerator(n, counts)

    # 5. Construct the struct explicitly, preserving the DESIGNED b and d
    return ReedSolomonCode(F, E, R, β, n, k, b, d, d, d, cosets, qcosets_reps, def_set, g, h, e, cache)
end

"""
$(TYPEDSIGNATURES)

Return the Cyclic code of length `n` over `GF(q)` defined by the roots given as 
integer exponents in `elements` (where the roots are `β^i` for `i ∈ elements`).
If `type = :nonzeros` is passed, `elements` are treated as the non-root exponents.
"""
function CyclicCode(q::Int, n::Int, elements::Vector{Int}; type::Symbol = :zeros)
    if type == :zeros
        zeros_idx = elements
    elseif type == :nonzeros
        zeros_idx = setdiff(0:(n-1), elements)
    else
        throw(ArgumentError("The `type` keyword argument must be either `:zeros` or `:nonzeros`."))
    end
    
    return CyclicCode(q, n, defining_set(zeros_idx, q, n, false))
end

"""
$(TYPEDSIGNATURES)

Return the Cyclic code of length `n` over `GF(q)` defined by the exact field elements 
passed in `elements`.
If `type = :nonzeros`, they are treated as the non-roots.
"""
function CyclicCode(q::Int, n::Int, elements::Vector{<:CTFieldElem}; type::Symbol = :zeros)
    isempty(elements) && return CyclicCode(q, n, Int[]; type=type)
    
    factors = Nemo.factor(q)
    (p, t), = factors
    deg = ord(n, q)
    E = GF(p, t * deg, :α)
    α = (t * deg == 1) ? E(2) : gen(E)
    β = α^(div(BigInt(q)^deg - 1, n))
    
    # Build a rapid lookup table for the powers of β
    beta_powers = Dict{typeof(β), Int}()
    curr_power = one(E)
    for i in 0:(n - 1)
        beta_powers[curr_power] = i
        curr_power *= β
    end
    
    indices = Int[]
    for elem in elements
        parent(elem) == E || throw(ArgumentError("Field elements must belong to the splitting field."))
        haskey(beta_powers, elem) || throw(ArgumentError("Element is not a valid power of β within range n."))
        push!(indices, beta_powers[elem])
    end
    
    return CyclicCode(q, n, indices; type=type)
end

# covered nicely in van Lint and Betten et al
"""
$(TYPEDSIGNATURES)

Return the cyclic code whose roots are the quadratic residues of `q`, `n`.
"""
QuadraticResidueCode(q::Int, n::Int) = CyclicCode(q, n, defining_set(quadratic_residues(q, n)[1], q, n, false))

"""
$(TYPEDSIGNATURES)

Return the fire code with generator polynomial `(x^(2l - 1) + 1) * p`.
"""
function FireCode(p::Union{fpPolyRingElem, FqPolyRingElem}, l::Int)
    Oscar.is_irreducible(p) || throw(ArgumentError("The polynomial `p` must be irreducible."))
    x = gen(parent(p))
    1 ≤ l ≤ degree(p) || throw(DomainError(l, "This construction requires 1 ≤ l ≤ degree(p)."))
    isone(gcd(p, x^(2l - 1) + 1)) || throw(ArgumentError("This construction requires gcd(p, x^(2l - 1) + 1) = 1."))
    
    g = (x^(2l - 1) + 1) * p
    F = base_ring(p)
    q = Int(order(F))
    m = degree(p)
    
    # Fast period calculation for irreducible p(x).
    # The period e is the smallest integer such that p(x) | x^e - 1.
    e = 1
    curr_pow = x
    max_period = q^m - 1
    for i in 1:max_period
        if isone(curr_pow)
            e = i
            break
        end
        curr_pow = mod(curr_pow * x, p)
    end
    
    # Because gcd(p, x^(2l-1)+1) = 1, the period of g(x) is the LCM of their individual periods.
    n = lcm(2l - 1, e)
    
    return CyclicCode(n, g)
end

#############################
      # getter functions
#############################

"""
$(TYPEDSIGNATURES)

Return the splitting field of the generator polynomial.
"""
splitting_field(C::AbstractCyclicCode) = C.E

"""
$(TYPEDSIGNATURES)

Return the polynomial ring of the generator polynomial.
"""
polynomial_ring(C::AbstractCyclicCode) = C.R

"""
$(TYPEDSIGNATURES)

Return the primitive root of the splitting field.
"""
primitive_root(C::AbstractCyclicCode) = C.β

"""
$(TYPEDSIGNATURES)

Return the offset of the BCH code.
"""
offset(C::AbstractBCHCode) = C.b

"""
$(TYPEDSIGNATURES)

Return the design distance of the BCH code.
"""
design_distance(C::AbstractBCHCode) = C.δ

"""
$(TYPEDSIGNATURES)

Return the q-cyclotomic cosets of the cyclic code.
"""
qcosets(C::AbstractCyclicCode) = C.qcosets

"""
$(TYPEDSIGNATURES)

Return the set of representatives for the q-cyclotomic cosets of the cyclic code.
"""
qcosets_reps(C::AbstractCyclicCode) = C.qcosets_reps

"""
$(TYPEDSIGNATURES)

Return the defining set (as integer exponents) of the cyclic code.
"""
defining_set(C::AbstractCyclicCode) = C.def_set

"""
$(TYPEDSIGNATURES)

Return the zeros (as field elements) of `C`.
"""
zeros(C::AbstractCyclicCode) = [C.β^i for i in C.def_set]

"""
$(TYPEDSIGNATURES)

Return the nonzeros (as field elements) of `C`.
"""
nonzeros(C::AbstractCyclicCode) = [C.β^i for i in setdiff(0:C.n - 1, C.def_set)]

"""
$(TYPEDSIGNATURES)

Return the generator polynomial of the cyclic code.
"""
generator_polynomial(C::AbstractCyclicCode) = C.g

"""
$(TYPEDSIGNATURES)

Return the parity-check polynomial of the cyclic code.
"""
parity_check_polynomial(C::AbstractCyclicCode) = C.h

"""
$(TYPEDSIGNATURES)

Return the idempotent (polynomial) of the cyclic code.
"""
idempotent(C::AbstractCyclicCode) = C.e

# """
# $(TYPEDSIGNATURES)

# Return the BCH bound for `C`.
# """
# BCH_bound(C::AbstractCyclicCode) = C.δ

# """
#     HT_bound(C::AbstractCyclicCode)

# Return the Hartmann-Tzeng refinement to the BCH bound for `C`.

# This is a lower bound on the minimum distance of `C`.
# """
# HT_bound(C::AbstractCyclicCode) = C.HT

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

function _generator_polynomial(R::CTPolyRing, β::CTFieldElem, Z::Vector{Int})
    g = one(R)
    x = gen(R)
    for i in Z
        g *= (x - β^i)
    end
    return g
end
_generator_polynomial(R::CTPolyRing, β::CTFieldElem, qcosets::Vector{Vector{Int}}) = _generator_polynomial(R, β, reduce(vcat, qcosets))

function _idempotent(g::CTPolyRingElem, h::CTPolyRingElem, n::Int)
    # Solve 1 = a(x)g(x) + b(x)h(x) for a(x), then e(x) = a(x)g(x) mod x^n - 1
    d, a, b = gcdx(g, h)
    return mod(g * a, gen(parent(g))^n - 1)
end

function _classify_factors(poly::fpPolyRingElem)
    # the problem with reverse(poly) == poly is that coefficents(poly) is not a fixed size
    n = degree(poly)
    F = base_ring(parent(poly))
    F0 = F(0)
    facs = [x[1] for x in collect(factor(poly))]
    self_reciprocal_facs = Vector{fpPolyRingElem}()
    pairs = Vector{Tuple{fpPolyRingElem, fpPolyRingElem}}()
    temp = zeros(F, 1, n + 1)
    temp2 = zeros(F, 1, n + 1)
    for f in facs
        f_coeffs = collect(coefficients(f))
        temp[1, end - length(f_coeffs) + 1:end] .= reverse(f_coeffs)
        temp2[1, end - length(f_coeffs) + 1:end] .= f_coeffs
        if temp == temp2
            push!(self_reciprocal_facs, f)
        else
            for f2 in facs
                if f != f2
                    f_coeffs = collect(coefficients(f2))
                    temp2[1, 1:end - length(f_coeffs)] .= F0
                    temp2[1, end - length(f_coeffs) + 1:end] .= f_coeffs
                    if temp == temp2
                        push!(pairs, (f, f2))
                        break
                    end
                end
            end
        end
        temp[1, :] .= F0
        temp2[1, :] .= F0
    end

    return self_reciprocal_facs, pairs
end

# TODO: make flat optional throughout
"""
$(TYPEDSIGNATURES)

Return the set of `q`-cyclotomic cosets of the numbers in `nums` modulo `n`.
If `flat` is true, returns a single sorted array of the defining set.
"""
function defining_set(nums::Vector{Int}, q::Int, n::Int, flat::Bool = true)
    arr = Vector{Vector{Int}}()
    seen = Set{Int}()
    
    for x in nums
        if !(x in seen)
            Cx = cyclotomic_coset(x, q, n)
            push!(arr, Cx)
            union!(seen, Cx)
        end
    end

    flat && return sort!(reduce(vcat, arr))
    return arr
end

"""
$(TYPEDSIGNATURES)

Return the Mattson--Solomon polynomial of a vector `v` over `F`.

# Notes
* The Mattson-Solomon transform is the finite field equivalent of the 
  Discrete Fourier Transform (DFT).
"""
function MattsonSolomon_transform(v::Vector{<:CTFieldElem}, α::CTFieldElem)
    n = length(v)
    E = parent(α)
    R, z = polynomial_ring(E, "z")
    
    MS = zero(R)
    for j in 1:n
        # A_j = sum_{i=1}^n v_i * α^(i * j)
        A_j = sum(v[i] * α^(i * j) for i in 1:n)
        MS += A_j * z^(n - j)
    end
    
    return MS
end

"""
$(TYPEDSIGNATURES)

Return the vector recovered by applying the inverse Mattson--Solomon transform
to the Mattson--Solomon polynomial `MS`.
"""
function inverse_MattsonSolomon_transform(MS::CTPolyRingElem, n::Int, α::CTFieldElem)
    E = parent(α)
    coeffs = [coeff(MS, n - j) for j in 1:n]
    
    v = elem_type(E)[]
    n_inv = inv(E(n))
    
    for i in 1:n
        v_i = sum(coeffs[j] * α^(-i * j) for j in 1:n)
        push!(v, n_inv * v_i)
    end
    
    return v
end

"""
$(TYPEDSIGNATURES)

Compute the BCH design distance `δ`, the offset `b`, and the consecutive runs 
for a cyclic code of length `n` with defining set `def_set`.
"""
function _BCH_bound_math(def_set::Vector{Int}, n::Int)
    isempty(def_set) && return 1, 0, Vector{Int}[]
    
    extended_set = vcat(def_set, def_set .+ n)
    
    runs = Vector{Vector{Int}}()
    current_run = [extended_set[1]]
    
    for i in 2:length(extended_set)
        if extended_set[i] == extended_set[i-1] + 1
            push!(current_run, extended_set[i])
        elseif extended_set[i] != extended_set[i-1]
            push!(runs, current_run)
            current_run = [extended_set[i]]
        end
    end
    push!(runs, current_run)
    
    run_lens = length.(runs)
    max_len, ind = findmax(run_lens)
    
    δ = min(n, max_len + 1)
    offset = mod(runs[ind][1], n)
    
    return δ, offset, runs
end

"""
$(TYPEDSIGNATURES)

Compute the Hartmann-Tzeng lower bound on the minimum distance.

# Notes
* Searches for arithmetic progressions of shifted runs. 
* Evaluated against all consecutive root runs, not just the maximal BCH run.
"""
function _HT_bound_math(def_set::Vector{Int}, n::Int, runs::Vector{Vector{Int}})
    HT = 1
    def_set_fast = Set(def_set)
    
    for A in runs
        δ_local = length(A) + 1
        base_run = [mod(x, n) for x in A]
        
        for c in 1:(n - 1)
            # The shift multiplier must be coprime-ish to the local design distance
            if gcd(c, n) < δ_local
                s = 1
                while true
                    shifted_run = [mod(x + s * c, n) for x in base_run]
                    if all(x -> x in def_set_fast, shifted_run)
                        s += 1
                    else
                        break
                    end
                end
                
                if HT < δ_local + (s - 1)
                    HT = δ_local + (s - 1)
                end
            end
        end
    end
    
    return min(HT, n)
end

"""
$(TYPEDSIGNATURES)

Compute the Roos lower bound on the minimum distance.

# Notes
* The Roos bound generalizes the HT bound by allowing the shifting set `B` 
  to be arbitrary, provided its size `s+1` satisfies `max(B) - min(B) <= δ + s - 2`.
"""
function _Roos_bound_math(def_set::Vector{Int}, n::Int, runs::Vector{Vector{Int}})
    Roos = 1
    def_set_fast = Set(def_set)
    
    for A in runs
        δ_local = length(A) + 1
        base_run = [mod(x, n) for x in A]
        
        # 1. Find all valid shifts c in 0:n-1 such that base_run + c is in the defining set
        valid_c = Int[]
        for c in 0:(n - 1)
            shifted_run = [mod(x + c, n) for x in base_run]
            if all(x -> x in def_set_fast, shifted_run)
                push!(valid_c, c)
            end
        end
        
        # 2. For each valid starting shift, normalize the set B to start at 0
        for c_start in valid_c
            normalized_c = sort!([mod(c - c_start, n) for c in valid_c])
            
            # 3. Greedily find the maximum size s+1 of a subset B satisfying the span bound
            for k_idx in 1:length(normalized_c)
                v = normalized_c[k_idx]
                s = k_idx - 1 # size of subset B is s + 1
                
                # Roos Condition: max(B) <= δ_local + |B| - 3
                if v <= δ_local + s - 2
                    if Roos < δ_local + s
                        Roos = δ_local + s
                    end
                end
            end
        end
    end
    
    return min(Roos, n)
end

"""
$(TYPEDSIGNATURES)

Return the BCH bound for `C`. Computes lazily and caches the result.
"""
function BCH_bound(C::AbstractCyclicCode)
    cache = getfield(C, :cache)
    if !haskey(cache, :BCH_bound)
        δ, offset, runs = _BCH_bound_math(C.def_set, C.n)
        cache[:BCH_bound] = δ
        cache[:bch_offset] = offset
        cache[:bch_runs] = runs # Stored to accelerate HT_bound
    end
    return cache[:BCH_bound]
end

"""
$(TYPEDSIGNATURES)

Return the offset of the BCH bound for `C`. Computes lazily.
"""
function BCH_offset(C::AbstractCyclicCode)
    cache = getfield(C, :cache)
    if !haskey(cache, :bch_offset)
        BCH_bound(C) # Guarantees offset is populated
    end
    return cache[:bch_offset]
end

"""
$(TYPEDSIGNATURES)

Return the Hartmann-Tzeng bound for `C`. Computes lazily and caches the result.
"""
function HT_bound(C::AbstractCyclicCode)
    cache = getfield(C, :cache)
    if !haskey(cache, :HT_bound)
        BCH_bound(C) # Guarantees base runs are populated in cache
        runs = cache[:bch_runs]
        
        cache[:HT_bound] = _HT_bound_math(C.def_set, C.n, runs)
    end
    return cache[:HT_bound]
end

"""
$(TYPEDSIGNATURES)

Return the Roos bound for `C`. Computes lazily and caches the result.
"""
function Roos_bound(C::AbstractCyclicCode)
    cache = getfield(C, :cache)
    if !haskey(cache, :Roos_bound)
        BCH_bound(C) # Guarantees base runs are populated
        runs = cache[:bch_runs]
        
        cache[:Roos_bound] = _Roos_bound_math(C.def_set, C.n, runs)
    end
    return cache[:Roos_bound]
end

"""
$(TYPEDSIGNATURES)

Return the defining set of the dual code of length `n` and defining set `def_set`.

# Notes
* Mathematically, if C has roots β^i for i ∈ Z, the dual has roots β^(-i) for i ∉ Z.
"""
dual_defining_set(def_set::Vector{Int}, n::Int) = sort!([mod(n - i, n) for i in setdiff(0:n - 1, def_set)])

"""
$(TYPEDSIGNATURES)

Return `true` and the equivalent cyclic code object if `C` is a cyclic code;
otherwise, return `false, missing`.
"""
function is_cyclic(C::AbstractLinearCode)::Tuple{Bool, Union{Missing, AbstractCyclicCode}}
    typeof(C) <: AbstractCyclicCode && return true, C
    
    ord_F = Int(order(C.F))
    gcd(C.n, ord_F) == 1 || return false, missing
    
    G = generator_matrix(C)
    H_trans = transpose(parity_check_matrix(C))
    nr, nc = nrows(G), ncols(G)
    
    # Extract the exact concrete element type to prevent Vector{Any} fallbacks
    T = typeof(G[1, 1])
    
    # Strictly typed preallocation for row shifts
    shifted_coeffs = Vector{T}(undef, nc)
    for r in 1:nr
        shifted_coeffs[1] = G[r, nc]
        for c in 2:nc
            shifted_coeffs[c] = G[r, c - 1]
        end
        
        shifted_row = matrix(C.F, 1, nc, shifted_coeffs)
        iszero(shifted_row * H_trans) || return false, missing
    end
    
    R, _ = polynomial_ring(C.F, :x)
    
    # Strictly typed preallocation for polynomial coefficients
    g_coeffs = Vector{T}(undef, nc)
    for c in 1:nc
        g_coeffs[c] = G[1, c]
    end
    g = R(g_coeffs)
    
    for r in 2:nr
        for c in 1:nc
            g_coeffs[c] = G[r, c]
        end
        g = gcd(g, R(g_coeffs))
    end
    
    g = divexact(g, leading_coefficient(g))
    
    return true, CyclicCode(C.n, g)
end

"""
$(TYPEDSIGNATURES)

Return the cyclic code whose cyclotomic cosets are the complement of `C`'s.
"""
function complement(C::AbstractCyclicCode)
    q = Int(order(C.F))
    comp_cosets = complement_qcosets(q, C.n, C.qcosets)
    return CyclicCode(q, C.n, comp_cosets)
end

"""
$(TYPEDSIGNATURES)

Return whether or not `C1` is a subcode of `C2`.
A cyclic code is a subcode of another if and only if its defining set is a superset of the other's.
"""
⊆(C1::AbstractCyclicCode, C2::AbstractCyclicCode) = C2.def_set ⊆ C1.def_set
"""
$(TYPEDSIGNATURES)

Return whether `C1` is a subcode of `C2`, not necessarily properly. This is an
alias for `⊆`, matching the behavior for general linear codes; use `⊊` to test
proper containment.
"""
⊂(C1::AbstractCyclicCode, C2::AbstractCyclicCode) = C1 ⊆ C2

"""
$(TYPEDSIGNATURES)

Return `true` if `C1` is a subcode of `C2`, and `false` otherwise.
For cyclic codes this holds exactly when the defining set of `C1` contains
the defining set of `C2`.
"""
is_subcode(C1::AbstractCyclicCode, C2::AbstractCyclicCode) = C1 ⊆ C2

"""
$(TYPEDSIGNATURES)

Return `true` if `C1` and `C2` are mathematically equal (same field, length, defining sets, and primitive root).
"""
==(C1::AbstractCyclicCode, C2::AbstractCyclicCode) = C1.F == C2.F && C1.n == C2.n && C1.def_set == C2.def_set && C1.β == C2.β

# function μa(C::CyclicCode)
#     # check gcd(a, n) = 1
#     # technically changes g(x) and e(x) but the q-cosets are the same?
# end

"""
$(TYPEDSIGNATURES)

Return the dual of the cyclic code `C`. 
This operation evaluates in O(1) time using the dual defining set properties.
"""
function dual(C::AbstractCyclicCode)
    q = Int(order(C.F))
    dual_def = dual_defining_set(C.def_set, C.n)
    return CyclicCode(q, C.n, dual_def; type=:zeros)
end

"""
$(TYPEDSIGNATURES)

Return the intersection code of `C1` and `C2`.
"""
function ∩(C1::AbstractCyclicCode, C2::AbstractCyclicCode)
    C1.F == C2.F && C1.n == C2.n && C1.β == C2.β || throw(ArgumentError("Cannot intersect codes over different fields, lengths, or primitive roots."))
    
    # The intersection has generator polynomial lcm(g1, g2), which corresponds to the union of defining sets.
    return CyclicCode(Int(order(C1.F)), C1.n, union(C1.def_set, C2.def_set); type=:zeros)
end

"""
$(TYPEDSIGNATURES)

Return the addition code of `C1` and `C2`.
"""
function +(C1::AbstractCyclicCode, C2::AbstractCyclicCode)
    C1.n == C2.n || throw(ArgumentError("Codes must have the same length."))
    C1.F == C2.F || throw(ArgumentError("Codes must be over the same base field."))
    
    # The defining set of the sum of two cyclic codes is the intersection of their defining sets.
    new_def_set = intersect(C1.def_set, C2.def_set)
    
    # If the intersection is empty, this simply yields the ambient space (generator polynomial = 1).
    # We pass it directly to the constructor which already handles empty defining sets correctly.
    q = Int(order(C1.F))
    return CyclicCode(q, C1.n, new_def_set)
end

"""
$(TYPEDSIGNATURES)

Return whether or not `C == dual(C)`.
"""
is_self_dual(C::AbstractCyclicCode) = C == dual(C)

"""
$(TYPEDSIGNATURES)

Return `true` if the BCH code is narrow-sense (offset `b` is 1 or 0 depending on convention).
"""
is_narrowsense(C::AbstractBCHCode) = iszero(C.b) || isone(C.b)

"""
$(TYPEDSIGNATURES)

Return `true` if the cyclic code is reversible.
"""
is_reversible(C::AbstractCyclicCode) = [mod(C.n - i, C.n) for i in C.def_set] ⊆ C.def_set

"""
$(TYPEDSIGNATURES)

Return `true` if the cyclic code is degenerate.
A cyclic code is degenerate if the parity-check polynomial divides `x^r - 1` for some `r < n`.
"""
function is_degenerate(C::AbstractCyclicCode)
    x = gen(C.R)
    for r in 1:C.n - 1
        flag, _ = divides(x^r - 1, C.h)
        flag && return true
    end
    return false
end

"""
$(TYPEDSIGNATURES)

Return `true` if the BCH code is primitive.
"""
is_primitive(C::AbstractBCHCode) = C.n == Int(order(C.F)) - 1

"""
$(TYPEDSIGNATURES)

Return `true` if the BCH code is antiprimitive.
"""
is_antiprimitive(C::AbstractBCHCode) = C.n == Int(order(C.F)) + 1

"""
$(TYPEDSIGNATURES)

Return the entrywise (Schur / Hadamard) product of cyclic codes `C1` and `C2`.

# Notes
* By the Mattson-Solomon Transform, the non-zeros of the Schur product 
  are the Minkowski sum of the non-zeros of the constituent codes.
"""
function entrywise_product_code(C1::AbstractCyclicCode, C2::AbstractCyclicCode)
    C1.F == C2.F && C1.n == C2.n && C1.β == C2.β || throw(ArgumentError("Cannot compute Schur product of codes over different fields, lengths, or primitive roots."))
    
    # 1. Extract non-zeros (as integer exponents)
    N1 = setdiff(0:(C1.n - 1), C1.def_set)
    N2 = setdiff(0:(C2.n - 1), C2.def_set)
    
    # 2. Compute the Minkowski sum of the non-zeros modulo n
    N_sum = Set{Int}()
    for x in N1
        for y in N2
            push!(N_sum, mod(x + y, C1.n))
        end
    end
    
    # 3. The defining set (zeros) is the complement of the new non-zeros
    def_set_new = setdiff(0:(C1.n - 1), N_sum)
    
    # 4. Route back into the O(1) lazy constructor
    q = Int(order(C1.F))
    return CyclicCode(q, C1.n, def_set_new; type=:zeros)
end

"""
$(TYPEDSIGNATURES)

Return the entrywise (Schur / Hadamard) product of `C` with itself.
"""
entrywise_product_code(C::AbstractCyclicCode) = entrywise_product_code(C, C)

# Aliases
*(C1::AbstractCyclicCode, C2::AbstractCyclicCode) = entrywise_product_code(C1, C2)
*(C::AbstractCyclicCode) = entrywise_product_code(C)

"""
$(TYPEDSIGNATURES)

Return the entrywise product code of `C1` and `C2`, or the entrywise square
of `C`. This is an alias for `entrywise_product_code`.
"""
Schur_product_code(C1::AbstractCyclicCode, C2::AbstractCyclicCode) = entrywise_product_code(C1, C2)
Schur_product_code(C::AbstractCyclicCode) = entrywise_product_code(C)

"""
$(TYPEDSIGNATURES)

Return the entrywise product code of `C1` and `C2`, or the entrywise square
of `C`. This is an alias for `Schur_product_code`.
"""
Hadamard_product_code(C1::AbstractCyclicCode, C2::AbstractCyclicCode) = entrywise_product_code(C1, C2)
Hadamard_product_code(C::AbstractCyclicCode) = entrywise_product_code(C)

"""
$(TYPEDSIGNATURES)

Return the entrywise product code of `C1` and `C2`, or the entrywise square
of `C`. This is an alias for `Schur_product_code`.
"""
componentwise_product_code(C1::AbstractCyclicCode, C2::AbstractCyclicCode) = entrywise_product_code(C1, C2)
componentwise_product_code(C::AbstractCyclicCode) = entrywise_product_code(C)

"""
$(TYPEDSIGNATURES)

Return the Trace representation parameters of the cyclic code `C`.

# Notes
* Returns a vector of the primitive root powers that define the independent 
  trace components of the codewords.
* A codeword `c` can be generated by `c_i = sum_j Tr_{E/F}( A_j * β^(i * j) )` 
  where `A_j` are arbitrary elements in the splitting field.
"""
function trace_representation(C::AbstractCyclicCode)
    # The trace representation is defined exactly by the non-zeros of the code
    # reduced to their cyclotomic coset representatives.
    non_zeros = setdiff(0:(C.n - 1), C.def_set)
    q = Int(order(C.F))
    
    # Extract the unique coset representatives for the non-zeros
    trace_reps = Int[]
    seen = Set{Int}()
    
    for nz in non_zeros
        if !(nz in seen)
            Cx = cyclotomic_coset(nz, q, C.n)
            push!(trace_reps, Cx[1])
            union!(seen, Cx)
        end
    end
    
    return sort!(trace_reps)
end

# ==============================================================================
# MULTIPLIERS & EQUIVALENCE
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Apply the multiplier `a` to the vector `v` of length `n`.

# Notes
* The multiplier maps the coordinate index `i` to `a * i (mod n)`.
* Requires `gcd(a, n) == 1` to ensure the mapping is a valid permutation.
"""
function apply_multiplier(v::Vector{T}, a::Int) where T
    n = length(v)
    gcd(a, n) == 1 || throw(ArgumentError("Multiplier 'a' must be coprime to the vector length 'n'."))
    
    v_new = similar(v)
    for i in 0:(n - 1)
        # Using 0-based index math, then adjusting to 1-based for Julia arrays
        new_idx = mod(a * i, n) + 1
        v_new[new_idx] = v[i + 1]
    end
    
    return v_new
end

"""
$(TYPEDSIGNATURES)

Return a new cyclic code by applying the multiplier `a` to the cyclic code `C`.

# Notes
* Algebraically, this multiplies the defining set of the code by `a (mod n)`.
* This implies the codewords of the new code are a permutation of the original codewords.
"""
function apply_multiplier(C::AbstractCyclicCode, a::Int)
    gcd(a, C.n) == 1 || throw(ArgumentError("Multiplier 'a' must be coprime to the code length 'n'."))
    
    # Multiply the defining set by 'a' modulo 'n'
    new_def_set = sort!(unique([mod(a * x, C.n) for x in C.def_set]))
    
    q = Int(order(C.F))
    return CyclicCode(q, C.n, new_def_set; type=:zeros)
end

"""
$(TYPEDSIGNATURES)

Return `true` and the multiplier `a` if the cyclic codes `C1` and `C2` are multiplier equivalent.
Otherwise, return `false` and `missing`.

# Notes
* Two cyclic codes are multiplier equivalent if there exists some multiplier `a` 
  coprime to `n` such that `μ_a(C1) == C2`.
* Multiplier equivalence implies the codes are permutation equivalent and share 
  the same weight enumerator and minimum distance.
"""
function is_multiplier_equivalent(C1::AbstractCyclicCode, C2::AbstractCyclicCode)
    C1.n == C2.n || return false, missing
    C1.F == C2.F || return false, missing
    C1.k == C2.k || return false, missing # Must have the same dimension
    
    # Brute force search over the multiplicative group of integers modulo n
    for a in 1:(C1.n - 1)
        if gcd(a, C1.n) == 1
            # Fast defining set check
            test_set = sort!(unique([mod(a * x, C1.n) for x in C1.def_set]))
            if test_set == C2.def_set
                return true, a
            end
        end
    end
    
    return false, missing
end

"""
$(TYPEDSIGNATURES)

Return the set of all multipliers `a` that map the cyclic code `C` strictly to itself.

# Notes
* This set forms a subgroup of the multiplicative group `(Z/nZ)*`.
* These multipliers correspond to the automorphisms of the cyclic code.
"""
function multiplier_group(C::AbstractCyclicCode)
    group = Int[]
    for a in 1:(C.n - 1)
        if gcd(a, C.n) == 1
            test_set = sort!(unique([mod(a * x, C.n) for x in C.def_set]))
            if test_set == C.def_set
                push!(group, a)
            end
        end
    end
    return group
end

"""
$(TYPEDSIGNATURES)

Return the multiplier group of `C` as a formal subgroup of the symmetric group `S_n`.

# Notes
* Returns `(H, f)`, where `H` is the subgroup and `f` is the inclusion morphism `H -> S_n`.
* The multiplier action `i -> a * i (mod n)` is internally shifted to 1-based indexing 
  to match Oscar's permutation group standards.
"""
function multiplier_subgroup_Sn(C::AbstractCyclicCode)
    M = multiplier_group(C)
    Sn = symmetric_group(C.n)
    
    # Map each multiplier to a formal permutation in S_n
    perms = elem_type(Sn)[]
    for a in M
        # 0-based coordinate math shifted to 1-based permutation array
        p_array = [mod(a * (i - 1), C.n) + 1 for i in 1:C.n]
        push!(perms, Sn(p_array))
    end
    
    # Return the formal GAP subgroup and its injection
    return sub(Sn, perms)
end

"""
$(TYPEDSIGNATURES)

Return the multiplier group of `C` as a formal subgroup of the unit group `Z_n^x`.

# Notes
* Returns `(H, inc)`, where `H` is the abstract abelian subgroup and `inc` is the injection.
* To map an element `h ∈ H` back to an integer, use the unit group isomorphism:
  `R, _ = residue_ring(ZZ, C.n); U, f = unit_group(R); int_val = lift(f(inc(h)))`
"""
function multiplier_subgroup_Zn(C::AbstractCyclicCode)
    M = multiplier_group(C)
    
    # Construct the residue ring Z/nZ and its abstract unit group
    R, _ = residue_ring(ZZ, C.n)
    U, f = unit_group(R)
    
    # Find the abstract group elements corresponding to our integer multipliers
    preimages = [f \ R(a) for a in M]
    
    # Return the formal subgroup
    return sub(U, preimages)
end

"""
$(TYPEDSIGNATURES)

Search for a valid `m`-adic splitting of the non-zero `q`-cyclotomic cosets modulo `n`.
Return `(true, a, S)` where `a` is the cycling multiplier and `S` is an array of `m` defining sets.
If no such splitting exists, returns `(false, missing, missing)`.
"""
function _find_polyadic_splittings(q::Int, n::Int, m::Int)
    gcd(q, n) == 1 || throw(ArgumentError("gcd(q, n) must be 1 for cyclotomic cosets."))
    
    # 1. Isolate the non-zero cosets
    all_cosets = all_cyclotomic_cosets(q, n, to_sort=true, verbose=false)
    nz_cosets = filter(c -> c != [0], all_cosets)
    
    # Fast failure: The number of non-zero cosets must be divisible by m
    if length(nz_cosets) % m != 0
        return false, missing, missing
    end

    # 2. Search for a multiplier 'a' that generates orbits of exactly length 'm'
    for a in 1:(n - 1)
        if gcd(a, n) == 1
            orbits = Vector{Vector{Vector{Int}}}()
            unvisited = copy(nz_cosets)
            valid_multiplier = true

            while !isempty(unvisited)
                start_coset = unvisited[1]
                curr_orbit = [start_coset]
                curr_coset = start_coset

                # Apply the multiplier 'a' repeatedly
                while true
                    next_rep = mod(a * curr_coset[1], n)
                    idx = findfirst(c -> next_rep in c, nz_cosets)
                    isnothing(idx) && error("Mathematical mapping failure.")
                    next_coset = nz_cosets[idx]

                    if next_coset == start_coset
                        break # Closed the orbit
                    end
                    if next_coset in curr_orbit
                        valid_multiplier = false # Orbit looped improperly
                        break
                    end
                    push!(curr_orbit, next_coset)
                    curr_coset = next_coset # FIX: Advance the orbit tracker!
                end

                # For an m-adic splitting, EVERY orbit must have a length of exactly m
                if !valid_multiplier || length(curr_orbit) != m
                    valid_multiplier = false
                    break
                end

                push!(orbits, curr_orbit)
                setdiff!(unvisited, curr_orbit)
            end

            # 3. If valid, distribute the orbits across m defining sets
            if valid_multiplier
                S = [Int[] for _ in 1:m]
                for orb in orbits
                    for i in 1:m
                        append!(S[i], orb[i])
                    end
                end
                
                # Sort the generated defining sets
                for i in 1:m
                    sort!(S[i])
                end
                
                return true, a, S
            end
        end
    end
    
    return false, missing, missing
end

"""
$(TYPEDSIGNATURES)

Return the family of `m` Polyadic codes of length `n` over `GF(q)`.

# Notes
* Searches for a multiplier `μ_a` that splits the non-zero roots into `m` cycling sets.
* Returns a NamedTuple containing the array of `m` codes and the multiplier `a` used.
* If `include_zero = true`, the root `0` is added to all defining sets (yielding the even-like subcodes).
"""
function PolyadicCodes(q::Int, n::Int, m::Int; include_zero::Bool = false)
    found, a, S = _find_polyadic_splittings(q, n, m)
    
    if !found
        error("No polyadic splitting of order $m exists for q=$q and n=$n.")
    end
    
    if include_zero
        for i in 1:m
            S[i] = sort!(vcat(S[i], [0]))
        end
    end
    
    # Leverage our O(1) lazy constructor
    codes = [CyclicCode(q, n, S[i]; type=:zeros) for i in 1:m]
    
    return (codes = codes, multiplier = a)
end

"""
$(TYPEDSIGNATURES)

Return the pair of Duadic codes of length `n` over `GF(q)`.

# Notes
* Duadic codes exist if and only if there is a multiplier of order 2 that splits the roots.
* Quadratic Residue (QR) codes are a special, highly symmetric case of Duadic codes.
"""
DuadicCodes(q::Int, n::Int; include_zero::Bool = false) = PolyadicCodes(q, n, 2; include_zero=include_zero)

"""
$(TYPEDSIGNATURES)

Return the triplet of Triadic codes of length `n` over `GF(q)`.
"""
TriadicCodes(q::Int, n::Int; include_zero::Bool = false) = PolyadicCodes(q, n, 3; include_zero=include_zero)

"""
$(TYPEDSIGNATURES)

Return the quad of Tetradic codes of length `n` over `GF(q)`.
"""
TetradicCodes(q::Int, n::Int; include_zero::Bool = false) = PolyadicCodes(q, n, 4; include_zero=include_zero)

# ==============================================================================
# CONSTITUENTS AND IRREDUCIBILITY
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the irreducible cyclic constituents of the cyclic code `C`.

# Notes
* By the Chinese Remainder Theorem, every cyclic code is a direct sum of 
  irreducible cyclic codes.
* These constituents correspond to the individual `q`-cyclotomic cosets 
  that make up the non-zeros (the trace representation) of `C`.
"""
function constituents(C::AbstractCyclicCode)
    # The constituents are defined by the non-zeros of the code
    non_zeros = setdiff(0:(C.n - 1), C.def_set)
    isempty(non_zeros) && return typeof(C)[] # The zero code has no constituents
    
    q = Int(order(C.F))
    seen = Set{Int}()
    constituent_codes = typeof(C)[]
    
    for nz in non_zeros
        if !(nz in seen)
            # Find the cyclotomic coset for this specific non-zero
            coset = cyclotomic_coset(nz, q, C.n)
            union!(seen, coset)
            
            # An irreducible constituent has exactly this coset as its ONLY non-zeros.
            # Therefore, its defining set (zeros) is the complement of this coset.
            constituent_def_set = setdiff(0:(C.n - 1), coset)
            
            # Instantiated instantly via the O(1) lazy constructor
            push!(constituent_codes, CyclicCode(q, C.n, constituent_def_set; type=:zeros))
        end
    end
    
    return constituent_codes
end

"""
$(TYPEDSIGNATURES)

Return `true` if the cyclic code `C` is irreducible.

# Notes
* An irreducible cyclic code has no non-trivial cyclic subcodes.
* Algebraically, this occurs if and only if its non-zeros form exactly 
  one `q`-cyclotomic coset (meaning it has exactly one constituent).
"""
function is_irreducible(C::AbstractCyclicCode)
    # A code is irreducible if it is composed of exactly one constituent
    return length(constituents(C)) == 1
end

"""
$(TYPEDSIGNATURES)

Return all irreducible cyclic codes of length `n` over `GF(q)`.

# Notes
* This decomposes the entire ambient space `F_q[x]/<x^n - 1>` into its 
  minimal ideals.
"""
function ambient_constituents(q::Int, n::Int)
    gcd(q, n) == 1 || throw(ArgumentError("Code length `n` must be coprime to the field size `q`."))
    
    all_cosets = all_cyclotomic_cosets(q, n, to_sort=true, verbose=false)
    
    # Typed correctly to allow auto-detected BCHCode structs
    constituent_codes = AbstractCyclicCode[] 
    for coset in all_cosets
        constituent_def_set = setdiff(0:(n - 1), coset)
        push!(constituent_codes, CyclicCode(q, n, constituent_def_set; type=:zeros))
    end
    
    return constituent_codes
end
