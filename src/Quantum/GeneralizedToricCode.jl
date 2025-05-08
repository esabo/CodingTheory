# Copyright (c) 2025 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

function GeneralizedToricCode(f::CTLRPolyElem, g::CTLRPolyElem)
    LR = parent(f)
    LR == parent(g) || throw(ArgumentError("The polynomials must be over the same ring."))
    length(symbols(LR)) == 2 || throw(ArgumentError("The polynomials must be over a Laurent polynomial ring in two variables."))
    return GeneralizedToricCode(LR, base_ring(LR), f, g)
end

function FiniteGeneralizedToricCode(f::CTLRPolyElem, g::CTLRPolyElem, a1::Tuple{Int, Int},
    a2::Tuple{Int, Int})

    LR = parent(f)
    LR == parent(g) || throw(ArgumentError("The polynomials must be over the same ring."))
    length(symbols(LR)) == 2 || throw(ArgumentError("The polynomials must be over a Laurent polynomial ring in two variables."))
    return FiniteGeneralizedToricCode(LR, base_ring(LR), f, g, a1, a2)
end

#############################
      # getter functions
#############################

Laurent_polynomial_ring(S::AbstractGeneralizedToricCode) = S.LR

field(S::AbstractGeneralizedToricCode) = S.F

defining_polynomials(S::AbstractGeneralizedToricCode) = S.f, S.g

twist_vectors(S::FiniteGeneralizedToricCode) = S.a1, S.a2

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

function maximum_dimension(S::AbstractGeneralizedToricCode)
    R2 = Oscar._polyringquo(S.LR)
    R = codomain(R2)
    (x, y) = gens(S.LR)
    if isa(S, FiniteGeneralizedToricCode)
        I = ideal(S.LR, [S.f, S.g, x^S.a1[1] * y^S.a1[2] - 1, x^S.a2[1] * y^S.a2[2] - 1])
    else
        I = ideal(S.LR, [S.f, S.g])
    end
    II = ideal(R, R2.(gens(I)))
    Q, ϕ = quo(R, II)
    return 2 * vector_space_dimension(Q)
end

function CSSCode(S::FiniteGeneralizedToricCode)
    R2 = Oscar._polyringquo(S.LR)
    R = codomain(R2)
    (x, y) = gens(S.LR)
    swap = hom(S.LR, S.LR, [x^-1, y^-1])
    f_anti = swap(S.f)
    g_anti = swap(S.g)
    f_R2 = R2(S.f)
    g_R2 = R2(S.g)
    f_anti_R2 = R2(f_anti)
    g_anti_R2 = R2(g_anti)
    I = ideal(S.LR, [x^S.a1[1] * y^S.a1[2] - 1, x^S.a2[1] * y^S.a2[2] - 1])
    II = ideal(R, R2.(gens(I)))
    Q, ϕ = quo(R, II)

    mono = monomial_basis(Q)
    len_mon = length(mono)
    n = 2 * len_mon
    LR_edge_index = Dict{fpMPolyRingElem, Int}(mono[i] => i for i in 1:len_mon)
    TB_edge_index = Dict{fpMPolyRingElem, Int}(mono[i] => i + len_mon for i in 1:len_mon)

    Fone = S.F(1)
    row = 1
    X_stabs = zero_matrix(S.F, len_mon, n)
    Z_stabs = zero_matrix(S.F, len_mon, n)
    for edge in mono
        f_shift = simplify(ϕ(edge * f_R2))
        for term in terms(f_shift.f)
            # X_12
            X_stabs[row, TB_edge_index[term]] = Fone
        end
        g_shift = simplify(ϕ(edge * g_R2))
        for term in terms(g_shift.f)
            # X_14
            X_stabs[row, LR_edge_index[term]] = Fone
        end

        g_shift = simplify(ϕ(edge * g_anti_R2))
        for term in terms(g_shift.f)
            # Z_12
            Z_stabs[row, TB_edge_index[term]] = Fone
        end
        f_shift = simplify(ϕ(edge * f_anti_R2))
        for term in terms(f_shift.f)
            # Z_14
            Z_stabs[row, LR_edge_index[term]] = Fone
        end
        row += 1
    end
    # println("X stabs")
    # display(X_stabs)
    # println("Z stabs")
    # display(Z_stabs)

    return CSSCode(X_stabs, Z_stabs)
end

function show(io::IO, S::AbstractGeneralizedToricCode)
    if isa(S, FiniteGeneralizedToricCode)
        println(io, "Finite Generalized Toric Code:")
        println(io, "\tf: $(S.f)")
        println(io, "\tg: $(S.g)")
        println(io, "\ta1: $(S.a1)")
        println(io, "\ta2: $(S.a2)")
    else
        println(io, "Generalized Toric Code:")
        println(io, "\tf: $(S.f)")
        println(io, "\tg: $(S.g)")
    end
end
