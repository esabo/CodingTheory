# Copyright (c) 2022 - 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

function generator_matrix(C::GeneralizedReedSolomonCode, stand_form::Bool = false)
    cache = getfield(C, :cache)
    if !haskey(cache, :G)
        G = zero_matrix(C.F, C.k, C.n)
        for c in 1:C.n
            for r in 1:C.k
                G[r, c] = C.scalars[c] * C.eval_pts[c]^(r - 1)
            end
        end
        cache[:G] = G
    end
    if stand_form
        if !haskey(cache, :G_stand)
            G_stand, H_stand, P, _ = _standard_form(cache[:G])
            cache[:G_stand] = G_stand
            cache[:H_stand] = H_stand
            cache[:P_stand] = P
        end
        return cache[:G_stand]
    end
    return cache[:G]
end

function parity_check_matrix(C::GeneralizedReedSolomonCode, stand_form::Bool = false)
    cache = getfield(C, :cache)
    if !haskey(cache, :H)
        H = zero_matrix(C.F, C.n - C.k, C.n)
        for c in 1:C.n
            for r in 1:(C.n - C.k)
                H[r, c] = C.dual_scalars[c] * C.eval_pts[c]^(r - 1)
            end
        end
        cache[:H] = H
    end
    if stand_form
        generator_matrix(C, true)
        return cache[:H_stand]
    end
    return cache[:H]
end

"""
    GeneralizedReedSolomonCode(k::Int, v::Vector{FqFieldElem}, γ::Vector{FqFieldElem})

Return the dimension `k` Generalized Reed-Solomon code with scalars `v` and
evaluation points `γ`.

# Notes
* The vectors `v` and `γ` must have the same length and every element must be over the same field.
* The elements of `v` need not be distinct but must be nonzero.
* The elements of `γ` must be distinct.
"""
function GeneralizedReedSolomonCode(k::Int, v::Vector{<:CTFieldElem}, γ::Vector{<:CTFieldElem})
    n = length(v)
    1 <= k <= n || throw(DomainError(k, "The dimension of the code must be between `1` and `n`."))
    n == length(γ) || throw(DomainError((n, length(γ)), "Lengths of scalars and evaluation points must be equal."))
    
    F = parent(v[1])
    1 <= n <= Int(order(F)) || throw(DomainError(n, "The length of the code must be between `1` and the order of the field."))
    
    for (i, x) in enumerate(v)
        iszero(x) && throw(ArgumentError("The elements of `v` must be nonzero."))
        parent(x) == F || throw(ArgumentError("The elements of `v` must be over the same field."))
        parent(γ[i]) == F || throw(ArgumentError("The elements of `γ` must be over the same field as `v`."))
    end
    length(unique(γ)) == n || throw(ArgumentError("The elements of `γ` must be distinct."))

    # Compute dual scalars via Lagrange interpolation
    w = elem_type(F)[]
    for i in 1:n
        push!(w, (v[i] * prod(γ[j] - γ[i] for j in 1:n if j ≠ i))^-1)
    end

    d = n - k + 1
    cache = Dict{Symbol, Any}()
    return GeneralizedReedSolomonCode(F, n, k, d, d, d, v, w, γ, cache)
end

# using notation of MacWilliams & Sloane, p. 340
"""
    GeneralizedReedSolomonCode(C::AbstractGoppaCode)

Return the generalized Reed-Solomon code associated with the Goppa code `C`.
"""
GeneralizedReedSolomonCode(C::AbstractGoppaCode) = GeneralizedReedSolomonCode(C.n - degree(C.g), [C.g(C.L[i]) * prod(C.L[j] - C.L[i] for j in 1:C.n if i ≠ j)^(-1) for i in 1:C.n], C.L)

"""
    GeneralizedReedSolomonCode(C::AbstractAlternateCode)

Return the generalized Reed-Solomon code associated with the alternate code `C`.
"""
GeneralizedReedSolomonCode(C::AbstractAlternateCode) = dual(GeneralizedReedSolomonCode(C.k, C.scalars, C.eval_pts))

"""
$(TYPEDSIGNATURES)

Return a random Generalized Reed-Solomon code of length `n` and dimension `k` over `F`.
Generates random distinct evaluation points and random non-zero scalars.
"""
function RandomGeneralizedReedSolomonCode(F::CTFieldTypes, n::Int, k::Int)
    1 <= k <= n || throw(DomainError((k, n), "Require 1 <= k <= n."))
    n <= Int(order(F)) || throw(DomainError(n, "Length cannot exceed the order of the field."))
    
    # 1. Random distinct evaluation points
    eval_pts = elem_type(F)[]
    while length(eval_pts) < n
        pt = rand(F)
        !(pt in eval_pts) && push!(eval_pts, pt)
    end
    
    # 2. Random non-zero scalars
    scalars = elem_type(F)[]
    while length(scalars) < n
        s = rand(F)
        !iszero(s) && push!(scalars, s)
    end
    
    return GeneralizedReedSolomonCode(k, scalars, eval_pts)
end

"""
$(TYPEDSIGNATURES)

Convert a cyclic `ReedSolomonCode` into its equivalent `GeneralizedReedSolomonCode` representation.
"""
function GeneralizedReedSolomonCode(C::ReedSolomonCode)
    n = C.n
    k = C.k
    β = primitive_root(C)
    b = offset(C)
    
    # 1. The GRS evaluation points for a cyclic RS code are consecutive powers of the primitive root
    γ = [β^(i - 1) for i in 1:n]
    
    # 2. The cyclic parity-check matrix implies dual GRS scalars of w_i = γ_i^b
    w = [γ[i]^b for i in 1:n]
    
    # 3. Compute the primal GRS scalars via Lagrange interpolation
    v = elem_type(C.F)[]
    for i in 1:n
        push!(v, (w[i] * prod(γ[j] - γ[i] for j in 1:n if j != i))^-1)
    end
    
    return GeneralizedReedSolomonCode(k, v, γ)
end

"""
$(TYPEDSIGNATURES)

Return the generalized Srivastava code over `F`. Evaluates lazily.
"""
function GeneralizedSrivastavaCode(F::CTFieldTypes, a::Vector{T}, w::Vector{T}, z::Vector{T}, t::Int) where T <: CTFieldElem
    isempty(a) && throw(ArgumentError("The input vector `a` cannot be empty."))
    isempty(w) && throw(ArgumentError("The input vector `w` cannot be empty."))
    isempty(z) && throw(ArgumentError("The input vector `z` cannot be empty."))
    t > 0 || throw(DomainError(t, "The parameter `t` must be positive"))
    
    n = length(a)
    n == length(z) || throw(ArgumentError("Vectors `a` and `z` must be the same length"))
    s = length(w)
    length(unique([a; w])) == n + s || throw(ArgumentError("Elements of `a` and `w` must be distinct"))
    any(iszero, z) && throw(DomainError(z, "Elements of `z` must be nonzero"))
    
    E = parent(a[1])
    all(parent(pt) == E for pt in a) || throw(ArgumentError("All elements of the input vector `a` must be over the same base ring."))
    all(parent(pt) == E for pt in w) || throw(ArgumentError("All elements of the input vector `w` must be over the same base ring as `a`."))
    all(parent(pt) == E for pt in z) || throw(ArgumentError("All elements of the input vector `z` must be over the same base ring as `a`."))
    
    flag, _ = is_subfield(F, E)
    flag || throw(ArgumentError("Input field is not a subfield of the base ring of the input vectors"))

    H = zero_matrix(E, s * t, n)
    for l in 1:s
        count = 1
        for r in (l - 1) * s + 1:(l - 1) * s + t
            for c in 1:n
                H[r, c] = z[c] * (a[c] - w[l])^(-count)
            end
            count += 1
        end
    end

    basis, _ = primitive_basis(E, F)
    if typeof(E) === typeof(F)
        H_exp = transpose(expand_matrix(transpose(H), F, basis))
    else
        H_exp = change_base_ring(F, transpose(expand_matrix(transpose(H), Oscar.Nemo.Native.GF(Int(order(F))), basis)))
    end
    
    # We must eagerly evaluate the subcode dimension because of the expansion
    rnk_H = rank(H_exp)
    k = n - rnk_H
    
    cache = Dict{Symbol, Any}(:H => H_exp)
    return GeneralizedSrivastavaCode(F, E, n, k, missing, s * t + 1, n, a, w, z, t, cache)
end

"""
    SrivastavaCode(F::CTFieldTypes, a::Vector{T}, w::Vector{T}, z::Vector{T}, t::Int) where T <: CTFieldElem

Return the Srivastava code over `F` given `a`, `w`, and `z`.

# Notes
- These inputs are defined on page 357 of MacWilliams & Sloane
"""
SrivastavaCode(F::CTFieldTypes, a::Vector{T}, w::Vector{T}, z::Vector{T}) where T <: CTFieldElem = GeneralizedSrivastavaCode(F, a, w, z, 1)

"""
$(TYPEDSIGNATURES)

Return the Generalized BCH code over `F` with evaluation points `γ`, design distance `δ`, 
and offset `b`.

# Notes
* A Generalized BCH code is a subfield subcode of a specific Generalized Reed-Solomon code,
  making it a special case of an Alternant code where the dual scalars are `w_i = γ_i^b`.
"""
function GeneralizedBCHCode(F::CTFieldTypes, γ::Vector{<:CTFieldElem}, δ::Int, b::Int=1)
    δ >= 2 || throw(DomainError(δ, "Design distance must be >= 2."))
    E = parent(γ[1])
    all(parent(pt) == E for pt in γ) || throw(ArgumentError("All evaluation points must be over the same field."))
    length(unique(γ)) == length(γ) || throw(ArgumentError("Evaluation points must be distinct."))
    
    flag, _ = is_subfield(F, E)
    flag || throw(ArgumentError("The base field F must be a subfield of the field of evaluation points."))
    
    n = length(γ)
    
    # 1. The defining parity-check rows for GBCH imply dual scalars w_i = γ_i^b
    w = [γ[i]^b for i in 1:n]
    
    # 2. Compute the primal GRS scalars via Lagrange interpolation
    v = elem_type(E)[]
    for i in 1:n
        push!(v, (w[i] * prod(γ[j] - γ[i] for j in 1:n if j != i))^-1)
    end
    
    # 3. The parent GRS code guarantees distance δ, so k_GRS = n - δ + 1
    k_GRS = n - δ + 1
    k_GRS > 0 || throw(DomainError(δ, "Design distance is too large for the number of evaluation points."))
    
    # 4. A GBCH code is simply the Alternant code over F derived from this GRS code
    # This automatically invokes the O(1) lazy subfield subcode logic.
    return AlternateCode(F, k_GRS, v, γ)
end

#############################
      # getter functions
#############################

"""
    scalars(C::GeneralizedReedSolomonCode)

Return the scalars `v` of the Generalized Reed-Solomon code `C`.
"""
scalars(C::GeneralizedReedSolomonCode) = C.scalars

"""
    dual_scalars(C::GeneralizedReedSolomonCode)

Return the scalars of the dual of the Generalized Reed-Solomon code `C`.
"""
dual_scalars(C::GeneralizedReedSolomonCode) = C.dual_scalars

"""
    evaluation_points(C::GeneralizedReedSolomonCode)

Return the evaluation points `γ` of the Generalized Reed-Solomon code `C`.
"""
evaluation_points(C::GeneralizedReedSolomonCode) = C.eval_pts

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

"""
    is_primitive(C::AbstractGeneralizedSrivastavaCode)

Return `true` if `C` is primitive
"""
is_primitive(C::AbstractGeneralizedSrivastavaCode) = C.n == Int(order(C.E)) - length(C.w)

"""
$(TYPEDSIGNATURES)

Compute the sequence of syndromes for the received vector `y` with respect to 
the Generalized Reed-Solomon code `C`.
"""
function syndromes(C::GeneralizedReedSolomonCode, y::Union{Vector{Int}, Vector{<:CTFieldElem}})
    F = C.F
    length(y) == C.n || throw(ArgumentError("Received vector length does not match code length."))
    
    # S_j = sum(y_i * dual_scalar_i * eval_pt_i^(j-1))
    n_k = C.n - C.k
    S = elem_type(F)[]
    
    for j in 1:n_k
        val = zero(F)
        for i in 1:C.n
            val += F(y[i]) * C.dual_scalars[i] * C.eval_pts[i]^(j - 1)
        end
        push!(S, val)
    end
    return S
end

"""
$(TYPEDSIGNATURES)

Compute the syndrome polynomial `S(z)` for the received vector `y` with respect to 
the Generalized Reed-Solomon code `C`. This polynomial is the standard input 
for the Berlekamp-Massey or Sugiyama Euclidean decoding algorithms.
"""
function syndrome_polynomial(C::GeneralizedReedSolomonCode, y::Union{Vector{Int}, Vector{<:CTFieldElem}})
    F = C.F
    R, z = polynomial_ring(F, "z")
    S_seq = syndromes(C, y)
    
    # S(z) = S_1 + S_2 z + S_3 z^2 + ... + S_{n-k} z^{n-k-1}
    poly = zero(R)
    for j in 1:length(S_seq)
        poly += S_seq[j] * z^(j - 1)
    end
    
    return poly
end

# We can safely map Alternate codes to their parent GRS codes to calculate syndromes over the extension field E.
syndromes(C::AbstractAlternateCode, y) = syndromes(GeneralizedReedSolomonCode(C), y)
syndrome_polynomial(C::AbstractAlternateCode, y) = syndrome_polynomial(GeneralizedReedSolomonCode(C), y)
