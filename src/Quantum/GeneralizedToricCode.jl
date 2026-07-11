# Copyright (c) 2025 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

function _BB_ansatz_check(f::T) where T <: Union{MPolyQuoRingElem{FqMPolyRingElem}, MPolyQuoRingElem{fpMPolyRingElem}}

    count_x = 0
    count_y = 0
    for e in exponents(f.f)
        if e[1] != 0 && e[2] != 0
            return false
        end
        if e[1] != 0
            count_x += 1
        end
        if e[2] != 0
            count_y += 1
        end
    end
    if count_x > 1 && count_y > 1
        return false
    end
    return true
end

"""
    BivariateBicycleCode(a::MPolyQuoRingElem{FqMPolyRingElem}, b::MPolyQuoRingElem{FqMPolyRingElem})

Return the bivariate bicycle code defined by the residue ring elements `a` and `b`.
"""
function BivariateBicycleCode(a::T, b::T) where T <: Union{MPolyQuoRingElem{FqMPolyRingElem},
    MPolyQuoRingElem{fpMPolyRingElem}}

    R = parent(a)
    R == parent(b) || throw(DomainError("Polynomials must have the same parent."))
    F = base_ring(base_ring(a))
    order(F) == 2 || throw(DomainError("This code family is currently only defined over binary fields."))
    length(symbols(parent(a))) == 2 || throw(DomainError("Polynomials must be over two variables."))
    g = gens(modulus(R))
    length(g) == 2 || throw(DomainError("Residue rings must have only two generators."))

    m = -1
    l = -1
    for g1 in g
        exps = collect(exponents(g1))
        length(exps) == 2 || throw(ArgumentError("Moduli of the incorrect form."))
        iszero(exps[2]) || throw(ArgumentError("Moduli of the incorrect form."))
        !iszero(exps[1][1]) && !iszero(exps[1][2]) && throw(ArgumentError("Moduli of the incorrect form."))
        if iszero(exps[1][1])
            m = exps[1][2]
        else
            l = exps[1][1]
        end
    end

    if !_BB_ansatz_check(a) || !_BB_ansatz_check(b)
        throw(ArgumentError("Polynomials do not satisfy the bivariate bicycle code ansatz."))
    end

    # already has the modulus built into R
    I = ideal(R, [a, b])
    Q, _ = quo(R, I)
    k_dim = 2 * vector_space_dimension(Q)

    return BivariateBicycleCode(R, F, 2 * l * m, k_dim, a, b, l, m)
end

"""
    CoprimeBivariateBicycleCode(a::MPolyQuoRingElem{FqMPolyRingElem}, b::MPolyQuoRingElem{FqMPolyRingElem})

Return the coprime bivariate bicycle code defined by the residue ring elements `a` and `b`.

# Note

- This is defined in https://arxiv.org/pdf/2408.10001.
"""
function CoprimeBivariateBicycleCode(a::ResElem, b::ResElem)
    R = parent(a)
    S = base_ring(a)
    R == parent(b) || throw(DomainError("Polynomials must have the same parent."))
    F = base_ring(S)
    order(F) == 2 || throw(DomainError("This code family is currently only defined over binary fields."))
    length(gens(S)) == 1 || throw(DomainError("Polynomials must be over one variable."))
    f = modulus(R)
    deg_P = degree(f)
    f == gen(S)^deg_P - 1 || throw(ArgumentError("Residue ring not of the form π^(l * m) - 1."))

    # BUG this is not particularly true since l or m could be factored itself and yet still be coprime
    facs = Nemo.factor(deg_P)
    length(facs) == 2 || throw(ArgumentError("Residue ring not of the form π^(l * m) - 1."))
    k = collect(keys(facs.fac))
    v = collect(values(facs.fac))
    l = k[1]^v[1]
    m = k[2]^v[2]

    k_dim = 2 * degree(gcd(a, b, f))
    return CoprimeBivariateBicycleCode(R, F, 2 * l * m, k_dim, a, b, l, m)
end

"""
    GeneralizedToricCode(f::CTLRPolyElem, g::CTLRPolyElem)

Return the generalized toric code defined by the polynomials `f` and `g`.
"""
function GeneralizedToricCode(f::CTLRPolyElem, g::CTLRPolyElem)
    LR = parent(f)
    LR == parent(g) || throw(ArgumentError("The polynomials must be over the same ring."))
    length(symbols(LR)) == 2 || throw(ArgumentError("The polynomials must be over a Laurent polynomial ring in two variables."))
    # TODO check polynomial ansatz here

    R2 = Oscar._polyringquo(LR)
    R = codomain(R2)
    (x, y) = gens(LR)
    I = ideal(LR, [f, g])
    II = ideal(R, R2.(gens(I)))
    Q, _ = quo(R, II)
    n = 2 * length(monomial_basis(Q))

    return GeneralizedToricCode(LR, base_ring(LR), n, f, g)
end

"""
    FiniteGeneralizedToricCode(f::CTLRPolyElem, g::CTLRPolyElem, a1::Tuple{Int, Int}, a2::Tuple{Int, Int})

Return the generalized toric code defined by the polynomials `f` and `g` and twist vectors `a1`, `a2`.
"""
function FiniteGeneralizedToricCode(f::CTLRPolyElem, g::CTLRPolyElem, a1::Tuple{Int, Int},
    a2::Tuple{Int, Int})

    LR = parent(f)
    LR == parent(g) || throw(ArgumentError("The polynomials must be over the same ring."))
    length(symbols(LR)) == 2 || throw(ArgumentError("The polynomials must be over a Laurent polynomial ring in two variables."))
    # TODO check polynomial ansatz here

    R2 = Oscar._polyringquo(LR)
    R = codomain(R2)
    (x, y) = gens(LR)
    I = ideal(LR, [f, g, x^a1[1] * y^a1[2] - 1, x^a2[1] * y^a2[2] - 1])
    II = ideal(R, R2.(gens(I)))
    Q, _ = quo(R, II)
    k_dim = 2 * vector_space_dimension(Q)
    n = 2 * length(monomial_basis(Q))

    return FiniteGeneralizedToricCode(LR, base_ring(LR), n, k_dim, f, g, a1, a2)
end

"""
    MonomialCode(f::CTLRPolyElem, g::CTLRPolyElem)

Return the monomial code defined by the polynomials `f` and `g`.
"""
function MonomialCode(f::CTLRPolyElem, g::CTLRPolyElem)
    LR = parent(f)
    LR == parent(g) || throw(ArgumentError("The polynomials must be over the same ring."))
    length(symbols(LR)) == 2 || throw(ArgumentError("The polynomials must be over a Laurent polynomial ring in two variables."))

    R2 = Oscar._polyringquo(LR)
    R = codomain(R2)
    (x, y) = gens(LR)
    I = ideal(LR, [f, g])
    II = ideal(R, R2.(gens(I)))
    Q, _ = quo(R, II)
    n = 2 * length(monomial_basis(Q))

    return MonomialCode(LR, base_ring(LR), n, f, g)
end

"""
    FiniteMonomialCode(f::CTLRPolyElem, g::CTLRPolyElem, a1::Tuple{Int, Int}, a2::Tuple{Int, Int})

Return the monomial code defined by the polynomials `f` and `g` and twist vectors `a1`, `a2`.
"""
function FiniteMonomialCode(f::CTLRPolyElem, g::CTLRPolyElem, a1::Tuple{Int, Int},
    a2::Tuple{Int, Int})

    LR = parent(f)
    LR == parent(g) || throw(ArgumentError("The polynomials must be over the same ring."))
    length(symbols(LR)) == 2 || throw(ArgumentError("The polynomials must be over a Laurent polynomial ring in two variables."))

    R2 = Oscar._polyringquo(LR)
    R = codomain(R2)
    (x, y) = gens(LR)
    I = ideal(LR, [f, g, x^a1[1] * y^a1[2] - 1, x^a2[1] * y^a2[2] - 1])
    II = ideal(R, R2.(gens(I)))
    Q, _ = quo(R, II)
    k_dim = 2 * vector_space_dimension(Q)
    n = 2 * length(monomial_basis(Q))

    return FiniteMonomialCode(LR, base_ring(LR), n, k_dim, f, g, a1, a2)
end

#############################
      # getter functions
#############################

polynomial_ring(S::AbstractBivariateBicycleCode) = S.R

Laurent_polynomial_ring(S::Union{AbstractMonomialCode, AbstractGeneralizedToricCode}) = S.LR

field(S::AbstractMonomialCode) = S.F

"""
    defining_polynomials(S::AbstractMonomialCode) -> Tuple{CTLRPolyElem, CTLRPolyElem}

Return the polynomials defining the monomial code `S`.
"""
defining_polynomials(S::AbstractMonomialCode) = S.f, S.g

"""
    twist_vectors(S::AbstractMonomialCode) -> Tuple{Tuple{Int, Int}, Tuple{Int, Int}}

Return the twist vectors of the monomial code `S` if they are defined.
"""
function twist_vectors(S::AbstractMonomialCode)
    if hasproperty(S, :a1) && hasproperty(S, :a2)
        return S.a1, S.a2
    elseif hasproperty(S, :l) && hasproperty(S, :m)
        return (S.l, 0), (0, S.m)
    else
        throw(ArgumentError("The twist vectors are not defined for this code.."))
    end
end

length(S::AbstractMonomialCode) = S.n

"""
    dimension(S::AbstractMonomialCode) -> Int

Return the dimension of the monomial code `S` if it is defined.
"""
function dimension(S::AbstractMonomialCode)
    if hasproperty(S, :k)
        return S.k
    else
        throw(ArgumentError("The dimension is not defined for this code. Use `maximum_dimension` instead."))
    end
end


#############################
      # setter functions
#############################

#############################
     # general functions
#############################

"""
    maximum_dimension(S::Union{MonomialCode, GeneralizedToricCode}) -> Int

Return the maximum dimension of the monomial code `S`.
"""
function maximum_dimension(S::Union{MonomialCode, GeneralizedToricCode})
    R2 = Oscar._polyringquo(S.LR)
    R = codomain(R2)
    (x, y) = gens(S.LR)
    I = ideal(S.LR, [S.f, S.g])
    II = ideal(R, R2.(gens(I)))
    Q, ϕ = quo(R, II)
    return 2 * vector_space_dimension(Q)
end

"""
    CSSCode(S::BivariateBicycleCode)

Return the CSS code defined by the bivariate bicycle code `S`.
"""
function CSSCode(S::BivariateBicycleCode)
    x = matrix(S.F, [mod1(i + 1, S.l) == j ? 1 : 0 for i in 1:S.l, j in 1:S.l]) ⊗ identity_matrix(S.F, S.m)
    y = identity_matrix(S.F, S.l) ⊗ matrix(S.F, [mod1(i + 1, S.m) == j ? 1 : 0 for i in 1:S.m, j in 1:S.m])

    A = zero_matrix(S.F, S.l * S.m, S.l * S.m)
    for ex in exponents(lift(S.f))
        # iszero(ex[1]) || iszero(ex[2]) || throw(ArgumentError("Polynomial `a` must not have any `xy` terms"))
        power, which = findmax(ex)
        if which == 1
            A += x^power
        elseif which == 2
            A += y^power
        end
    end

    B = zero_matrix(S.F, S.l * S.m, S.l * S.m)
    for ex in exponents(lift(S.g))
        # iszero(ex[1]) || iszero(ex[2]) || throw(ArgumentError("Polynomial `b` must not have any `xy` terms"))
        power, which = findmax(ex)
        if which == 1
            B += x^power
        elseif which == 2
            B += y^power
        end
    end


    return CSSCode(hcat(A, B), hcat(transpose(B), transpose(A)))
end

"""
    CSSCode(S::CoprimeBivariateBicycleCode)

Return the coprime bivariate bicycle code defined by the residue ring elements `a` and `b`.

# Note

- This is defined in https://arxiv.org/pdf/2408.10001.
"""
function CSSCode(S::CoprimeBivariateBicycleCode)
    deg_P = degree(modulus(S.R))
    x = matrix(S.F, [mod1(i + 1, S.l) == j ? 1 : 0 for i in 1:S.l, j in 1:S.l]) ⊗ identity_matrix(S.F, S.m)
    y = identity_matrix(S.F, S.l) ⊗ matrix(S.F, [mod1(i + 1, S.m) == j ? 1 : 0 for i in 1:S.m, j in 1:S.m])

    P = x * y
    A = zero_matrix(S.F, deg_P, deg_P)
    exps = findall(i -> !is_zero(i), collect(coefficients(lift(S.a)))) .- 1
    for ex in exps
        A += P^ex
    end

    B = zero_matrix(S.F, deg_P, deg_P)
    exps = findall(i -> !is_zero(i), collect(coefficients(lift(S.b)))) .- 1
    for ex in exps
        B += P^ex
    end

    return CSSCode(hcat(A, B), hcat(transpose(B), transpose(A)))
end

"""
    CSSCode(S::FiniteGeneralizedToricCode)

Return the CSS code defined by the finite generalized toric code `S`.
"""
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

    return CSSCode(X_stabs, Z_stabs)
end

# TODO add show methods for other codes
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

"""
    is_pure(S::AbstractMonomialCode) -> Bool

Return `true` if the logicals of the monomial code `S` are pure.
"""
function is_pure(S::AbstractMonomialCode)
    # TODO return true if l and m are odd?
    I_f = ideal(S.R, [S.f])
    I_g = ideal(S.R, [S.g])
    return I_f ∩ I_g == I_f * I_g
end

# """
#     is_principal(S::AbstractMonomialCode) -> Bool

# Return `true` if the logicals of the monomial code `S` are principal.
# """
# function is_principal(S::AbstractMonomialCode)
#     TODO return true if l and m are odd?
#     TODO: check if pure
#     I_f = ideal(S.R, [S.f])
#     I_g = ideal(S.R, [S.g])
#     # TODO bug here with types and annihilator
#     return is_principal(annihilator(I_f)) && is_principal(annihilator(I_g))
# end

#############################
         # 3D Trial
#############################

# this one is really the same as above
function Generalized3DToricCode(f::CTLRPolyElem, g::CTLRPolyElem)
    LR = parent(f)
    LR == parent(g) || throw(ArgumentError("The polynomials must be over the same ring."))
    length(symbols(LR)) == 3 || throw(ArgumentError("The polynomials must be over a Laurent polynomial ring in three variables."))
    return Generalized3DToricCode(LR, base_ring(LR), f, g)
end

function FiniteGeneralized3DToricCode(f::CTLRPolyElem, g::CTLRPolyElem, a1::Tuple{Int, Int},
    a2::Tuple{Int, Int}, l::Int)

    l > 0 || throw(DomainError(l, "Parameter l must be positive"))
    LR = parent(f)
    LR == parent(g) || throw(ArgumentError("The polynomials must be over the same ring."))
    length(symbols(LR)) == 3 || throw(ArgumentError("The polynomials must be over a Laurent polynomial ring in three variables."))
    return FiniteGeneralized3DToricCode(LR, base_ring(LR), f, g, a1, a2, l)
end

function maximum_dimension(S::AbstractGeneralized3DToricCode)
    R2 = Oscar._polyringquo(S.LR)
    R = codomain(R2)
    (x, y, z) = gens(S.LR)
    if isa(S, FiniteGeneralized3DToricCode)
        I = ideal(S.LR, [S.f, S.g, x^S.a1[1] * y^S.a1[2] - 1, x^S.a2[1] * y^S.a2[2] - 1, z^S.l - 1])
    else
        I = ideal(S.LR, [S.f, S.g])
    end
    II = ideal(R, R2.(gens(I)))
    Q, ϕ = quo(R, II)
    return 2 * vector_space_dimension(Q)
end

function CSSCode(S::FiniteGeneralized3DToricCode)
    R2 = Oscar._polyringquo(S.LR)
    R = codomain(R2)
    (x, y, z) = gens(S.LR)
    swap = hom(S.LR, S.LR, [x^-1, y^-1, z^-1])
    f_anti = swap(S.f)
    g_anti = swap(S.g)
    f_R2 = R2(S.f)
    g_R2 = R2(S.g)
    f_anti_R2 = R2(f_anti)
    g_anti_R2 = R2(g_anti)
    I = ideal(S.LR, [x^S.a1[1] * y^S.a1[2] - 1, x^S.a2[1] * y^S.a2[2] - 1, z^S.l - 1])
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

    return CSSCode(X_stabs, Z_stabs)
end
