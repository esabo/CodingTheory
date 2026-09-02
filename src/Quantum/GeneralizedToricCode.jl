# Copyright (c) 2025 - 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    BivariateBicycleCode(a::CTLRPolyElem, b::CTLRPolyElem, a1::Tuple{Int, Int}, a2::Tuple{Int, Int})

Return a lazy, twisted Bivariate Bicycle Code defined by the polynomials `a` and `b` and twist vectors `a1`, `a2`.
"""
function BivariateBicycleCode(a::CTLRPolyElem, b::CTLRPolyElem, a1::Tuple{Int, Int}, a2::Tuple{Int, Int})
    LR = parent(a)
    LR == parent(b) || throw(ArgumentError("The polynomials must be over the same ring."))
    length(symbols(LR)) == 2 || throw(ArgumentError("The polynomials must be over a Laurent polynomial ring in two variables."))
    
    F = base_ring(LR)
    
    # Algebraically compute n and k in O(1) time without building matrices
    R2 = Oscar._polyringquo(LR)
    R = codomain(R2)
    (x, y) = gens(LR)
    
    I = ideal(LR, [a, b, x^a1[1] * y^a1[2] - 1, x^a2[1] * y^a2[2] - 1])
    II = ideal(R, R2.(gens(I)))
    Q, _ = quo(R, II)
    
    k_dim = 2 * vector_space_dimension(Q)
    n_dim = 2 * length(monomial_basis(Q))
    
    cache = Dict{Symbol, Any}()
    
    return BivariateBicycleCode(
        LR, F, a, b, a1, a2, 
        n_dim, k_dim, missing, 1, n_dim, 
        cache
    )
end

"""
    BivariateBicycleCode(a::CTLRPolyElem, b::CTLRPolyElem, l::Int, m::Int)

Return a lazy, standard (untwisted) Bivariate Bicycle Code where `a1 = (l, 0)` and `a2 = (0, m)`.
"""
function BivariateBicycleCode(a::CTLRPolyElem, b::CTLRPolyElem, l::Int, m::Int)
    return BivariateBicycleCode(a, b, (l, 0), (0, m))
end

"""
    GeneralizedToricCode(a::CTLRPolyElem, b::CTLRPolyElem, a1::Tuple{Int, Int}, a2::Tuple{Int, Int})
    GeneralizedToricCode(a::CTLRPolyElem, b::CTLRPolyElem, l::Int, m::Int)

Return a lazy `BivariateBicycleCode`. The Generalized Toric Code family is modeled 
as a strict subset of Bivariate Bicycle codes with defined twist vectors.
"""
function GeneralizedToricCode(a::CTLRPolyElem, b::CTLRPolyElem, a1::Tuple{Int, Int}, a2::Tuple{Int, Int})
    return BivariateBicycleCode(a, b, a1, a2)
end

function GeneralizedToricCode(a::CTLRPolyElem, b::CTLRPolyElem, l::Int, m::Int)
    return BivariateBicycleCode(a, b, l, m)
end

function _evaluate_BB_matrices!(S::BivariateBicycleCode)
    haskey(S.cache, :H_X) && return
    
    R2 = Oscar._polyringquo(S.LR)
    R = codomain(R2)
    (x, y) = gens(S.LR)
    
    swap = hom(S.LR, S.LR, [x^-1, y^-1])
    a_anti = swap(S.a)
    b_anti = swap(S.b)
    
    a_R2 = R2(S.a)
    b_R2 = R2(S.b)
    a_anti_R2 = R2(a_anti)
    b_anti_R2 = R2(b_anti)
    
    I = ideal(S.LR, [x^S.a1[1] * y^S.a1[2] - 1, x^S.a2[1] * y^S.a2[2] - 1])
    II = ideal(R, R2.(gens(I)))
    Q, ϕ = quo(R, II)

    mono = monomial_basis(Q)
    len_mon = length(mono)
    n = S.n
    
    LR_edge_index = Dict(mono[i] => i for i in 1:len_mon)
    TB_edge_index = Dict(mono[i] => i + len_mon for i in 1:len_mon)

    Fone = S.F(1)
    row = 1
    
    # Preallocate matrices
    X_stabs = zero_matrix(S.F, len_mon, n)
    Z_stabs = zero_matrix(S.F, len_mon, n)
    
    for edge in mono
        a_shift = simplify(ϕ(edge * a_R2))
        for term in terms(a_shift.f)
            X_stabs[row, TB_edge_index[term]] = Fone
        end
        
        b_shift = simplify(ϕ(edge * b_R2))
        for term in terms(b_shift.f)
            X_stabs[row, LR_edge_index[term]] = Fone
        end

        b_shift_anti = simplify(ϕ(edge * b_anti_R2))
        for term in terms(b_shift_anti.f)
            Z_stabs[row, TB_edge_index[term]] = Fone
        end
        
        a_shift_anti = simplify(ϕ(edge * a_anti_R2))
        for term in terms(a_shift_anti.f)
            Z_stabs[row, LR_edge_index[term]] = Fone
        end
        row += 1
    end
    
    S.cache[:H_X] = X_stabs
    S.cache[:H_Z] = Z_stabs
end

# Accessors remain unchanged but call the updated evaluator:
function X_stabilizers(S::BivariateBicycleCode)
    _evaluate_BB_matrices!(S)
    return S.cache[:H_X]
end

function Z_stabilizers(S::BivariateBicycleCode)
    _evaluate_BB_matrices!(S)
    return S.cache[:H_Z]
end

function stabilizers(S::BivariateBicycleCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    H_X = X_stabilizers(S)
    H_Z = Z_stabilizers(S)
    F = S.F
    
    stabs = vcat(hcat(H_X, zero_matrix(F, size(H_X, 1), S.n)),
                 hcat(zero_matrix(F, size(H_Z, 1), S.n), H_Z))
                 
    S.cache[:stabilizers] = stabs
    return stabs
end

#############################
      # getter functions
#############################

Laurent_polynomial_ring(S::BivariateBicycleCode) = S.LR

field(S::BivariateBicycleCode) = S.F

"""
    defining_polynomials(S::BivariateBicycleCode) -> Tuple{CTLRPolyElem, CTLRPolyElem}

Return the polynomials defining the monomial code `S`.
"""
defining_polynomials(S::BivariateBicycleCode) = S.a, S.b

"""
    twist_vectors(S::BivariateBicycleCode) -> Tuple{Tuple{Int, Int}, Tuple{Int, Int}}

Return the twist vectors of the monomial code `S` if they are defined.
"""
twist_vectors(S::BivariateBicycleCode) = S.a1, S.a2

length(S::BivariateBicycleCode) = S.n

"""
    dimension(S::BivariateBicycleCode) -> Int

Return the dimension of the monomial code `S` if it is defined.
"""
dimension(S::BivariateBicycleCode) = S.k


#############################
      # setter functions
#############################

#############################
     # general functions
#############################

"""
    maximum_dimension(S::BivariateBicycleCode) -> Int

Return the maximum dimension of the monomial code `S`.
"""
function maximum_dimension(S::BivariateBicycleCode)
    R2 = Oscar._polyringquo(S.LR)
    R = codomain(R2)
    (x, y) = gens(S.LR)
    I = ideal(S.LR, [S.a, S.b])
    II = ideal(R, R2.(gens(I)))
    Q, ϕ = quo(R, II)
    return 2 * vector_space_dimension(Q)
end

function show(io::IO, S::BivariateBicycleCode)
    println(io, "Bivariate Bicycle Code:")
    println(io, "\tl: $(S.l)")
    println(io, "\tm: $(S.m)")
    println(io, "\ta: $(S.a)")
    println(io, "\tb: $(S.b)")
    println(io, "\ta1: $(S.a1)")
    println(io, "\ta2: $(S.a2)")
end

"""
    is_pure(S::AbstractMonomialCode) -> Bool

Return `true` if the logicals of the monomial code `S` are pure.
"""
function is_pure(S::BivariateBicycleCode)
    R2 = Oscar._polyringquo(S.LR)
    R = codomain(R2)
    (x, y) = gens(S.LR)
    
    # Quotient out the twists to get the proper base ring R
    I_twist = ideal(S.LR, [x^S.a1[1] * y^S.a1[2] - 1, x^S.a2[1] * y^S.a2[2] - 1])
    II_twist = ideal(R, R2.(gens(I_twist)))
    Q, ϕ = quo(R, II_twist)
    
    I_a = ideal(Q, [ϕ(R2(S.a))])
    I_b = ideal(Q, [ϕ(R2(S.b))])
    
    return I_a ∩ I_b == I_a * I_b
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

"""
    Generalized3DToricCode(a::CTLRPolyElem, b::CTLRPolyElem, a1::Tuple{Int, Int}, a2::Tuple{Int, Int}, l_z::Int)

Return a lazy, twisted 3D Generalized Toric Code defined by the polynomials `a` and `b`.
"""
function Generalized3DToricCode(a::CTLRPolyElem, b::CTLRPolyElem, a1::Tuple{Int, Int}, a2::Tuple{Int, Int}, l_z::Int)
    l_z > 0 || throw(DomainError(l_z, "Parameter l_z must be positive"))
    LR = parent(a)
    LR == parent(b) || throw(ArgumentError("The polynomials must be over the same ring."))
    length(symbols(LR)) == 3 || throw(ArgumentError("The polynomials must be over a Laurent polynomial ring in three variables."))
    
    F = base_ring(LR)
    
    R2 = Oscar._polyringquo(LR)
    R = codomain(R2)
    (x, y, z) = gens(LR)
    
    I = ideal(LR, [a, b, x^a1[1] * y^a1[2] - 1, x^a2[1] * y^a2[2] - 1, z^l_z - 1])
    II = ideal(R, R2.(gens(I)))
    Q, _ = quo(R, II)
    
    k_dim = 2 * vector_space_dimension(Q)
    n_dim = 2 * length(monomial_basis(Q))
    
    cache = Dict{Symbol, Any}()
    
    return Generalized3DToricCode(
        LR, F, a, b, a1, a2, l_z,
        n_dim, k_dim, missing, 1, n_dim, 
        cache
    )
end

"""
    Generalized3DToricCode(a::CTLRPolyElem, b::CTLRPolyElem, l_x::Int, l_y::Int, l_z::Int)

Return a lazy, standard (untwisted) 3D Generalized Toric Code.
"""
function Generalized3DToricCode(a::CTLRPolyElem, b::CTLRPolyElem, l_x::Int, l_y::Int, l_z::Int)
    return Generalized3DToricCode(a, b, (l_x, 0), (0, l_y), l_z)
end

function _evaluate_3D_matrices!(S::Generalized3DToricCode)
    haskey(S.cache, :H_X) && return
    
    R2 = Oscar._polyringquo(S.LR)
    R = codomain(R2)
    (x, y, z) = gens(S.LR)
    
    swap = hom(S.LR, S.LR, [x^-1, y^-1, z^-1])
    a_anti = swap(S.a)
    b_anti = swap(S.b)
    
    a_R2 = R2(S.a)
    b_R2 = R2(S.b)
    a_anti_R2 = R2(a_anti)
    b_anti_R2 = R2(b_anti)
    
    I = ideal(S.LR, [x^S.a1[1] * y^S.a1[2] - 1, x^S.a2[1] * y^S.a2[2] - 1, z^S.l_z - 1])
    II = ideal(R, R2.(gens(I)))
    Q, ϕ = quo(R, II)

    mono = monomial_basis(Q)
    len_mon = length(mono)
    n = S.n
    
    LR_edge_index = Dict(mono[i] => i for i in 1:len_mon)
    TB_edge_index = Dict(mono[i] => i + len_mon for i in 1:len_mon)

    Fone = S.F(1)
    row = 1
    
    X_stabs = zero_matrix(S.F, len_mon, n)
    Z_stabs = zero_matrix(S.F, len_mon, n)
    
    for edge in mono
        a_shift = simplify(ϕ(edge * a_R2))
        for term in terms(a_shift.f)
            X_stabs[row, TB_edge_index[term]] = Fone
        end
        
        b_shift = simplify(ϕ(edge * b_R2))
        for term in terms(b_shift.f)
            X_stabs[row, LR_edge_index[term]] = Fone
        end

        b_shift_anti = simplify(ϕ(edge * b_anti_R2))
        for term in terms(b_shift_anti.f)
            Z_stabs[row, TB_edge_index[term]] = Fone
        end
        
        a_shift_anti = simplify(ϕ(edge * a_anti_R2))
        for term in terms(a_shift_anti.f)
            Z_stabs[row, LR_edge_index[term]] = Fone
        end
        row += 1
    end
    
    S.cache[:H_X] = X_stabs
    S.cache[:H_Z] = Z_stabs
end

function X_stabilizers(S::Generalized3DToricCode)
    _evaluate_3D_matrices!(S)
    return S.cache[:H_X]
end

function Z_stabilizers(S::Generalized3DToricCode)
    _evaluate_3D_matrices!(S)
    return S.cache[:H_Z]
end

function stabilizers(S::Generalized3DToricCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    H_X = X_stabilizers(S)
    H_Z = Z_stabilizers(S)
    F = S.F
    
    stabs = vcat(hcat(H_X, zero_matrix(F, size(H_X, 1), S.n)),
                 hcat(zero_matrix(F, size(H_Z, 1), S.n), H_Z))
                 
    S.cache[:stabilizers] = stabs
    return stabs
end

function maximum_dimension(S::AbstractGeneralized3DToricCode)
    R2 = Oscar._polyringquo(S.LR)
    R = codomain(R2)
    (x, y, z) = gens(S.LR)
    if isa(S, FiniteGeneralized3DToricCode)
        I = ideal(S.LR, [S.a, S.b, x^S.a1[1] * y^S.a1[2] - 1, x^S.a2[1] * y^S.a2[2] - 1, z^S.l_z - 1])
    else
        I = ideal(S.LR, [S.a, S.b])
    end
    II = ideal(R, R2.(gens(I)))
    Q, ϕ = quo(R, II)
    return 2 * vector_space_dimension(Q)
end

"""
    CoprimeBivariateBicycleCode(a::CTPolyElem, b::CTPolyElem, N::Int)

Return a lazy Coprime Bivariate Bicycle Code (Generalized Bicycle Code) 
constructed directly from univariate polynomials `a` and `b` modulo `z^N - 1`.

When grid dimensions `l` and `m` are coprime, the bivariate group algebra is 
isomorphic to the univariate group algebra over Z_N (where N = lm). This allows 
the code to be constructed directly using a single shift parameter N, avoiding 
computationally expensive 2D algebraic geometries.
"""
function CoprimeBivariateBicycleCode(a::CTPolyRingElem, b::CTPolyRingElem, N::Int)
    R = parent(a)
    R == parent(b) || throw(ArgumentError("The polynomials must be over the same ring."))
    
    F = base_ring(R)
    
    # Algebraically compute k in O(1) time
    # The rank of circulant matrix [A | B] is N - deg(gcd(a, b, z^N - 1))
    # Therefore, k = 2N - 2 * rank = 2 * deg(gcd(a, b, z^N - 1))
    z = gen(R)
    mod_poly = z^N - 1
    g = gcd(a, gcd(b, mod_poly))
    
    k_dim = 2 * degree(g)
    n_dim = 2 * N
    
    cache = Dict{Symbol, Any}()
    
    return CoprimeBivariateBicycleCode(
        R, F, a, b, N,
        n_dim, k_dim, missing, 1, n_dim, 
        cache
    )
end

function _evaluate_CoprimeBB_matrices!(S::CoprimeBivariateBicycleCode)
    haskey(S.cache, :H_X) && return
    
    N = S.N
    F = S.F
    R = S.R
    z = gen(R)
    mod_poly = z^N - 1
    
    # Preallocate circulant matrices
    A = zero_matrix(F, N, N)
    B = zero_matrix(F, N, N)
    
    for i in 1:N
        # Row i corresponds to shifting the polynomial by z^(i-1)
        a_shift = (z^(i - 1) * S.a) % mod_poly
        for d in 0:(N - 1)
            A[i, d + 1] = coeff(a_shift, d)
        end
        
        b_shift = (z^(i - 1) * S.b) % mod_poly
        for d in 0:(N - 1)
            B[i, d + 1] = coeff(b_shift, d)
        end
    end
    
    # Standard Generalized Bicycle formulation:
    # H_X = [A | B]
    H_X = hcat(A, B)
    
    # H_Z = [B^T | A^T]
    H_Z = hcat(transpose(B), transpose(A))
    
    S.cache[:H_X] = H_X
    S.cache[:H_Z] = H_Z
end

function X_stabilizers(S::CoprimeBivariateBicycleCode)
    _evaluate_CoprimeBB_matrices!(S)
    return S.cache[:H_X]
end

function Z_stabilizers(S::CoprimeBivariateBicycleCode)
    _evaluate_CoprimeBB_matrices!(S)
    return S.cache[:H_Z]
end

function stabilizers(S::CoprimeBivariateBicycleCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    H_X = X_stabilizers(S)
    H_Z = Z_stabilizers(S)
    F = S.F
    
    stabs = vcat(hcat(H_X, zero_matrix(F, size(H_X, 1), S.n)),
                 hcat(zero_matrix(F, size(H_Z, 1), S.n), H_Z))
                 
    S.cache[:stabilizers] = stabs
    return stabs
end

polynomial_ring(S::CoprimeBivariateBicycleCode) = S.R
field(S::CoprimeBivariateBicycleCode) = S.F
defining_polynomials(S::CoprimeBivariateBicycleCode) = S.a, S.b
Base.length(S::CoprimeBivariateBicycleCode) = S.n
dimension(S::CoprimeBivariateBicycleCode) = S.k

function Base.show(io::IO, S::CoprimeBivariateBicycleCode)
    println(io, "Coprime Bivariate Bicycle Code (Generalized Bicycle Code):")
    println(io, "\tN: $(S.N)")
    println(io, "\ta: $(S.a)")
    println(io, "\tb: $(S.b)")
    println(io, "\t[n, k]: [$(S.n), $(S.k)]")
end
