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
    1 <= k <= n || throw(DomainError("The dimension of the code must be between 1 and n."))
    n == length(γ) || throw(DomainError("Lengths of scalars and evaluation points must be equal."))
    F = parent(v[1])
    1 <= n <= Int(order(F)) || throw(DomainError("The length of the code must be between 1 and the order of the field."))
    for (i, x) in enumerate(v)
        iszero(x) && throw(ArgumentError("The elements of v must be nonzero."))
        parent(x) == F || throw(ArgumentError("The elements of v must be over the same field."))
        parent(γ[i]) == F || throw(ArgumentError("The elements of γ must be over the same field as v."))
    end
    length(unique(γ)) == n || throw(ArgumentError("The elements of γ must be distinct."))

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

function GeneralizedSrivastavaCode(F, a, w, z, t::Int)
    is_positive(t) || throw(DomainError(t, "The parameter `t` must be positive"))
    n = length(a)
    n == length(z) || throw(ArgumentError("Vectors `a` and `z` must be the same length"))
    s = length(w)
    length(unique([a; w])) == n + s || throw(ArgumentError("Elements of `a` and `w` must be distinct"))
    any(iszero, z) && throw(DomainError(z, "Elements of `z` must be nonzero"))
    E = parent(a[1])
    flag, _ = is_subfield(F, E)
    flag || throw(ArgumentError("Input field is not a subfield of the base ring of the input veectors"))
    # TODO check all elements in the vectors have the same base ring
    # TODO repeat this for the constructors above

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
    return LinearCode(H_exp, true)

    # compute G
    G = kernel(H_exp, side = :right)
    k = rank(G)
    if ncols(G) == k
        G_tr = transpose(G)
    else
        # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
        nr = nrows(G)
        G_tr = zero_matrix(base_ring(G), k, nr)
        for r in 1:nr
            for c in 1:rnk_G
                !iszero(G[r, c]) && (G_tr[c, r] = G[r, c];)
            end
        end
    end
    G = G_tr
    n = ncols(G)

    G_stand, H_stand, P, _ = _standard_form(G)
    ub1, _ = _min_wt_row(G)
    ub2, _ = _min_wt_row(G_stand)
    ub = min(ub1, ub2)

    C = GeneralizedSrivastavaCode(F, E, n, k, missing, 1, ub, G, H_exp, G_stand, H_stand, P,
        missing, L, g)
    if BigInt(order(base_ring(G)))^min(k, n - k) <= 1.5e5
        C.weight_enum = if 2k <= n
            _weight_enumerator_BF(C.G_stand)
        else
            MacWilliams_identity(dual(C), _weight_enumerator_BF(C.H_stand))
        end
        d = minimum(filter(is_positive, first.(exponent_vectors(CWE_to_HWE(C.weight_enum).polynomial))))
        set_minimum_distance!(C, d)
    else
        l_bound = s * t + 1
        
        set_distance_lower_bound!(C, l_bound)
    end

    return C
end

SrivastavaCode(F, a, w, z) = GeneralizedSrivastavaCode(F, a, w, [z], 1)

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

# TODO
# write conversion function from Reed-Solomon code to GRS
# generalized BCH codes
