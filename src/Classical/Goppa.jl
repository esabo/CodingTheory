# Copyright (c) 2024 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
# constructors
#############################

function GoppaCode(F::CTFieldTypes, L::Vector{FqFieldElem}, g::FqPolyRingElem)
    is_empty(L) && throw(ArgumentError("The input vector `L` cannot be empty."))
    E = parent(L[1])
    all(parent(pt) == E for pt in L) || throw(
        ArgumentError(
            "All elements of the input vector `L` must be over the same base ring.",
        ),
    )
    E == base_ring(g) || throw(
        ArgumentError(
            "Input vector must be over the same base ring as the Goppa polynomial.",
        ),
    )
    rts = roots(g)
    isempty(L ∩ rts) || throw(
        ArgumentError(
            "The input vector must not contain any roots of the Goppa polynomial.",
        ),
    )
    is_subfield(F, E)[1] || throw(
        ArgumentError(
            "The input field is not a subfield of the base ring of the polynomial.",
        ),
    )

    n = length(L)
    t = degree(g)
    H = zero_matrix(E, t, n)
    for c = 1:n
        for r = 1:t
            H[r, c] = L[c]^(r - 1) * g(L[c])^(-1)
        end
    end

    basis, _ = primitive_basis(E, F)
    if typeof(E) === typeof(F)
        H_exp = transpose(expand_matrix(transpose(H), F, basis))
    else
        H_exp = change_base_ring(
            F,
            transpose(expand_matrix(transpose(H), GF(Int(order(F))), basis)),
        )
    end

    # just call LinearCode constructor here?

    # compute G
    G = kernel(H_exp, side = :right)
    k = rank(G)
    if ncols(G) == k
        G_tr = transpose(G)
    else
        # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
        nr = nrows(G)
        G_tr = zero_matrix(base_ring(G), k, nr)
        for r = 1:nr
            for c = 1:rnk_G
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

    C = GoppaCode(F, E, n, k, missing, 1, ub, G, H_exp, G_stand, H_stand, P, missing, L, g)
    if BigInt(order(base_ring(G)))^min(k, n - k) <= 1.5e5
        C.weight_enum = if 2k <= n
            _weight_enumerator_BF(C.G_stand)
        else
            MacWilliams_identity(dual(C), _weight_enumerator_BF(C.H_stand))
        end
        d = minimum(
            filter(
                is_positive,
                first.(exponent_vectors(CWE_to_HWE(C.weight_enum).polynomial)),
            ),
        )
        set_minimum_distance!(C, d)
    else
        l_bound = t + 1
        if Int(order(F)) == 2
            facs = factor(g)
            deg_g2 = 0
            for i in collect(values(facs.fac))
                iseven(i) ? (deg_g2 += i;) : (deg_g2 += i + 1;)
            end
            l_bound = deg_g2 + 1
        end
        set_distance_lower_bound!(C, l_bound)
    end

    return C
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
    length(fac.facs) == 1 ? (return true;) : (return false;)
end
