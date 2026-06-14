# Copyright (c) 2024 - 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
$(TYPEDSIGNATURES)

Return the Goppa code `Γ(L, g)` over `F`.
"""
function GoppaCode(F::CTFieldTypes, L::Vector{<:CTFieldElem}, g::CTPolyRingElem)
    isempty(L) && throw(ArgumentError("The input vector `L` cannot be empty."))
    E = parent(L[1])
    all(parent(pt) == E for pt in L) || throw(ArgumentError("All elements of the input vector `L` must be over the same base ring."))
    E == base_ring(g) || throw(ArgumentError("Input vector must be over the same base ring as the Goppa polynomial."))
    
    rts = roots(g)
    isempty(L ∩ rts) || throw(ArgumentError("The input vector must not contain any roots of the Goppa polynomial."))
    is_subfield(F, E)[1] || throw(ArgumentError("The input field is not a subfield of the base ring of the polynomial."))
    
    n = length(L)
    t = degree(g)
    H = zero_matrix(E, t, n)
    for c in 1:n
        for r in 1:t
            H[r, c] = L[c]^(r - 1) * g(L[c])^(-1)
        end
    end

    basis, _ = primitive_basis(E, F)
    if typeof(E) === typeof(F)
        H_exp = transpose(expand_matrix(transpose(H), F, basis))
    else
        H_exp = change_base_ring(F, transpose(expand_matrix(transpose(H), Oscar.Nemo.Native.GF(Int(order(F))), basis)))
    end

    # The exact dimension requires knowing the rank of the expanded parity-check matrix
    rnk_H = rank(H_exp)
    k = n - rnk_H

    l_bound = t + 1
    if Int(order(F)) == 2
        facs = factor(g)
        deg_g2 = 0
        for i in collect(values(facs.fac))
            deg_g2 += iseven(i) ? i : i + 1
        end
        l_bound = deg_g2 + 1
    end

    cache = Dict{Symbol, Any}(:H => H_exp)
    return GoppaCode(F, E, n, k, missing, l_bound, n, L, g, cache)
end

# ==============================================================================
# CRYPTOGRAPHIC GENERATORS
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return a random Goppa code of length `n` and Goppa polynomial degree `t`
over the base field `F`, using the extension field `E`. 

# Notes
* This is heavily used to generate public/private keypairs for McEliece cryptosystems.
"""
function RandomGoppaCode(F::CTFieldTypes, E::CTFieldTypes, n::Int, t::Int)
    flag, m = is_extension(E, F)
    flag || throw(ArgumentError("The field E must be an extension of F."))
    n <= BigInt(order(F))^m || throw(DomainError(n, "Length n cannot exceed the size of the extension field E."))
    t >= 1 || throw(DomainError(t, "Goppa polynomial degree t must be at least 1."))

    # 1. Generate a random Support L of unique elements from E
    L = elem_type(E)[]
    while length(L) < n
        pt = rand(E)
        if !(pt in L)
            push!(L, pt)
        end
    end

    # 2. Generate a random Irreducible Polynomial g(x) over E
    Rx, x = polynomial_ring(E, "x")
    local g
    while true
        # Generate random coefficients for degree 0 to t-1
        coeffs = [rand(E) for _ in 1:t]
        push!(coeffs, one(E)) # Ensure it is monic
        g = Rx(coeffs)
        
        if Oscar.is_irreducible(g)
            # Guarantee no roots exist in the support L
            has_root = false
            for pt in L
                if iszero(g(pt))
                    has_root = true
                    break
                end
            end
            if !has_root
                break
            end
        end
    end

    # 3. Construct the code using our optimized O(1) lazy architecture
    return GoppaCode(F, L, g)
end

#############################
      # getter functions
#############################

"""
    Goppa_polynomial(C::AbstractGoppaCode)

Return the Goppa polynomial of `C`.
"""
Goppa_polynomial(C::AbstractGoppaCode) = C.g

"""
    extension(C::AbstractGoppaCode)

Return the field over which the Goppa polynomial is defined.
"""
extension_field(C::AbstractGoppaCode) = C.E

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

"""
    is_seperable(C::AbstractGoppaCode)

Return true if the Goppa polynomail is seperable.
"""
is_irreducible(C::AbstractGoppaCode) = Oscar.is_irreducible(C.g)

"""
    is_seperable(C::AbstractGoppaCode)

Return true if the Goppa polynomail is seperable.
"""
is_seperable(C::AbstractGoppaCode) = is_squarefree(C.g)

"""
    nonzeros(C::AbstractGoppaCode)

Return the set `L` of the `Γ(L, g)` Goppa code.
"""
nonzeros(C::AbstractGoppaCode) = C.L

"""
    is_cumulative(C::AbstractGoppaCode)

Return true if the Goppa polynomail is of the form `g(z) = (z - β)^r`.
"""
function is_cumulative(C::AbstractGoppaCode)
    fac = factor(C.g)
    return length(fac.facs) == 1
end
