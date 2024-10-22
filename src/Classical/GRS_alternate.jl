# Copyright (c) 2022 - 2024 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    GeneralizedReedSolomonCode(k::Int, v::Vector{FqFieldElem}, γ::Vector{FqFieldElem})

Return the dimension `k` Generalized Reed-Solomon code with scalars `v` and
evaluation points `γ`.

# Notes
* The vectors `v` and `γ` must have the same length and every element must be over the same field.
* The elements of `v` need not be distinct but must be nonzero.
* The elements of `γ` must be distinct.
"""
function GeneralizedReedSolomonCode(k::Int, v::Vector{FqFieldElem}, γ::Vector{FqFieldElem})
    n = length(v)
    1 <= k <= n || throw(DomainError("The dimension of the code must be between `1` and `n`."))
    n == length(γ) || throw(DomainError("Lengths of scalars and evaluation points must be equal."))
    F = parent(v[1])
    1 <= n <= Int(order(F)) || throw(DomainError("The length of the code must be between `1` and the order of the field."))
    for (i, x) in enumerate(v)
        iszero(x) && throw(ArgumentError("The elements of `v` must be nonzero."))
        parent(x) == F || throw(ArgumentError("The elements of `v` must be over the same field."))
        parent(γ[i]) == F || throw(ArgumentError("The elements of `γ` must be over the same field as `v`."))
    end
    length(unique(γ)) == n || throw(ArgumentError("The elements of `γ` must be distinct."))

    G = zero_matrix(F, k, n)
    for c in 1:n
        for r in 1:k
            G[r, c] = v[c] * γ[c]^(r - 1)
        end
    end

    # evaluation points are the same for the parity check, but scalars are now
    # w_i = (v_i * \prod(γ_j - γ_i))^-1
    # which follows from Lagrange interpolation
    w = [F(0) for _ in 1:n]
    for i in 1:n
        w[i] = (v[i] * prod(γ[j] - γ[i] for j in 1:n if j ≠ i))^-1
    end

    H = zero_matrix(F, n - k, n)
    for c in 1:n
        for r in 1:n - k
            H[r, c] = w[c] * γ[c]^(r - 1)
        end
    end

    iszero(G * transpose(H)) || error("Calculation of dual scalars failed in constructor.")
    G_stand, H_stand, P, rnk = _standard_form(G)
    rnk == k || error("Computed rank not equal to desired input rank")
    d = n - k + 1
    return GeneralizedReedSolomonCode(F, n, k, d, d, d, v, w, γ, G, H,
        G_stand, H_stand, P, missing)
end

# using notation of MacWilliams & Sloane, p. 340
"""
    GeneralizedReedSolomonCode(C::AbstractGoppaCode)

Return the generalized Reed-Solomon code associated with the Goppa code `C`.
"""
GeneralizedReedSolomonCode(C::AbstractGoppaCode) = GeneralizedReedSolomonCode(C.n - degree(C.g), [C.g(C.L[i]) * prod(C.L[j] - C.L[i] for j in 1:C.n if i ≠ j)^(-1) for i in 1:C.n], C.L)

# doesn't have it's own struct
"""
    AlternateCode(F::CTFieldTypes, k::Int, v::Vector{FqFieldElem}, γ::Vector{FqFieldElem})

Return the alternate code as a subfield subcode of `GRS_k(v, γ)^⟂` over `F`.
"""
function AlternateCode(F::CTFieldTypes, k::Int, v::Vector{FqFieldElem},
    γ::Vector{FqFieldElem})

    E = parent(v[1])
    flag, _ = is_subfield(F, E)
    flag || throw(ArgumentError("Given field is not a subfield of the base ring of the vectors"))

    # let this do the rest of the type checking
    GRS = GeneralizedReedSolomonCode(k, v, γ)
    GRS = dual(GRS)
    basis, _ = primitive_basis(E, F)
    C = subfield_subcode(GRS, F, basis)
    return AlternateCode(F, E, C.n, C.k, C.d, C.l_bound, C.u_bound, v, γ, C.G, C.H, C.G_stand,
        C.H_stand, C.P_stand, C.weight_enum)
end

"""
   AlternateCode(F::CTFieldTypes, C::GeneralizedReedSolomonCode)

Return the subfield subcode of `C` over `F`.
"""
function AlternateCode(F::CTFieldTypes, C::GeneralizedReedSolomonCode)
    flag, _ = is_subfield(F, C.F)
    flag || throw(ArgumentError("Given field is not a subfield of the base ring of the code"))

    basis, _ = primitive_basis(C.F, F)
    return subfield_subcode(C, F, basis)
end

"""
    GeneralizedReedSolomonCode(C::AbstractAlternateCode)

Return the generalized Reed-Solomon code associated with the alternate code `C`.
"""
GeneralizedReedSolomonCode(C::AbstractAlternateCode) = dual(GeneralizedReedSolomonCode(C.k, C.scalars, C.eval_pts))

"""
    GeneralizedSrivastavaCode(F::CTFieldTypes, a::Vector{T}, w::Vector{T}, z::Vector{T}, t::Int) where T <: CTFieldElem

Return the generalized Srivastava code over `F` given `a`, `w`, `z`, and `t`.

# Notes
- These inputs are defined on page 357 of MacWilliams & Sloane
"""
function GeneralizedSrivastavaCode(F::CTFieldTypes, a::Vector{T}, w::Vector{T}, z::Vector{T},
    t::Int) where T <: CTFieldElem

    is_empty(a) && throw(ArgumentError("The input vector `a` cannot be empty."))
    is_empty(w) && throw(ArgumentError("The input vector `w` cannot be empty."))
    is_empty(z) && throw(ArgumentError("The input vector `z` cannot be empty."))
    is_positive(t) || throw(DomainError(t, "The parameter `t` must be positive"))
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
    flag || throw(ArgumentError("Input field is not a subfield of the base ring of the input veectors"))

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
        H_exp = change_base_ring(F, transpose(expand_matrix(transpose(H), GF(Int(order(F))), basis)))
    end
    C = LinearCode(H_exp, true)
    C2 = GeneralizedSrivastavaCode(F, E, C.n, C.k, C.d, C.l_bound, C.u_bound, a, w, z, t, C.G, C.H,
        C.G_stand, C.H_stand, C.P_stand, C.weight_enum)
    ismissing(C2.d) && set_distance_lower_bound!(C2, s * t + 1)

    return C2
end

"""
    SrivastavaCode(F::CTFieldTypes, a::Vector{T}, w::Vector{T}, z::Vector{T}, t::Int) where T <: CTFieldElem

Return the Srivastava code over `F` given `a`, `w`, and `z`.

# Notes
- These inputs are defined on page 357 of MacWilliams & Sloane
"""
SrivastavaCode(F::CTFieldTypes, a::Vector{T}, w::Vector{T}, z::Vector{T}) where T <: CTFieldElem =
    GeneralizedSrivastavaCode(F, a, w, z, 1)

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

# TODO
# write conversion function from Reed-Solomon code to GRS
# generalized BCH codes
