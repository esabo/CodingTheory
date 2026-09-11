# Copyright (c) 2025 - 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
$(TYPEDSIGNATURES)

Return the algebraic three-dimensional generalized toric-code datum defined by
`a` and `b`. This object is not a finite stabilizer code.
"""
function Generalized3DToricCode(a::CTLRPolyElem, b::CTLRPolyElem)
    _check_BB_pair(a, b)
    LR = parent(a)
    length(symbols(LR)) == 3 ||
        throw(ArgumentError("The polynomials must be in three Laurent variables."))
    return Generalized3DToricCode(LR, base_ring(base_ring(LR)), a, b)
end

"""
$(TYPEDSIGNATURES)

Return a finite, twisted three-dimensional generalized toric code.
"""
function FiniteGeneralized3DToricCode(
    a::CTLRPolyElem, b::CTLRPolyElem,
    a1::Tuple{Int, Int}, a2::Tuple{Int, Int}, l_z::Int
)
    _check_BB_pair(a, b)
    l_z > 0 || throw(DomainError(l_z, "The z period must be positive."))
    det = a1[1] * a2[2] - a1[2] * a2[1]
    iszero(det) &&
        throw(ArgumentError("The xy lattice vectors must be linearly independent."))
    LR = parent(a)
    length(symbols(LR)) == 3 ||
        throw(ArgumentError("The polynomials must be in three Laurent variables."))
    F = base_ring(base_ring(LR))
    order(F) == 2 ||
        throw(ArgumentError("Generalized 3D toric codes are currently implemented over binary fields."))

    R2 = Oscar._polyringquo(LR)
    R = codomain(R2)
    x, y, z = gens(LR)
    lattice = [x^a1[1] * y^a1[2] - 1,
               x^a2[1] * y^a2[2] - 1,
               z^l_z - 1]
    Q_lattice, _ = quo(R, ideal(R, R2.(lattice)))
    n = 2 * length(monomial_basis(Q_lattice))
    Q_code, _ = quo(R, ideal(R, R2.([a, b, lattice...])))
    k = 2 * vector_space_dimension(Q_code)
    twisted = !iszero(a1[2]) || !iszero(a2[1])
    result = FiniteGeneralized3DToricCode(
        LR, F, a, b, a1, a2, l_z, n, k, twisted, _BB_cache(F, n))
    return _seed_quantum_singleton_bound!(result)
end

"""
$(TYPEDSIGNATURES)

Return a finite, untwisted three-dimensional generalized toric code.
"""
function FiniteGeneralized3DToricCode(
    a::CTLRPolyElem, b::CTLRPolyElem, l_x::Int, l_y::Int, l_z::Int
)
    l_x > 0 || throw(DomainError(l_x, "The x period must be positive."))
    l_y > 0 || throw(DomainError(l_y, "The y period must be positive."))
    S = FiniteGeneralized3DToricCode(
        a, b, (l_x, 0), (0, l_y), l_z)
    S.twisted = false
    return S
end

"""
$(TYPEDSIGNATURES)

Return the corresponding member of the generalized three-dimensional
toric-code family.
Two polynomial arguments construct the algebraic object; lattice arguments
construct a `FiniteGeneralized3DToricCode`.
"""
BBCode3D(a::CTLRPolyElem, b::CTLRPolyElem) =
    Generalized3DToricCode(a, b)
BBCode3D(a::CTLRPolyElem, b::CTLRPolyElem,
         a1::Tuple{Int, Int}, a2::Tuple{Int, Int}, l_z::Int) =
    FiniteGeneralized3DToricCode(a, b, a1, a2, l_z)
BBCode3D(a::CTLRPolyElem, b::CTLRPolyElem,
         l_x::Int, l_y::Int, l_z::Int) =
    FiniteGeneralized3DToricCode(a, b, l_x, l_y, l_z)

function _evaluate_3D_matrices!(S::FiniteGeneralized3DToricCode)
    haskey(S.cache, :X_stabs) && return

    R2 = Oscar._polyringquo(S.LR)
    R = codomain(R2)
    x, y, z = gens(S.LR)
    anti = hom(S.LR, S.LR, [x^-1, y^-1, z^-1])
    a_R2, b_R2 = R2(S.a), R2(S.b)
    a_anti_R2, b_anti_R2 = R2(anti(S.a)), R2(anti(S.b))
    lattice = [x^S.a1[1] * y^S.a1[2] - 1,
               x^S.a2[1] * y^S.a2[2] - 1,
               z^S.l_z - 1]
    Q, phi = quo(R, ideal(R, R2.(lattice)))
    mono = monomial_basis(Q)
    q = length(mono)
    left = Dict(mono[i] => i for i in 1:q)
    right = Dict(mono[i] => i + q for i in 1:q)
    H_X = zero_matrix(S.F, q, 2q)
    H_Z = zero_matrix(S.F, q, 2q)
    one_F = S.F(1)

    for (row, edge) in enumerate(mono)
        for term in terms(simplify(phi(edge * a_R2)).f)
            H_X[row, right[term]] = one_F
        end
        for term in terms(simplify(phi(edge * b_R2)).f)
            H_X[row, left[term]] = one_F
        end
        for term in terms(simplify(phi(edge * b_anti_R2)).f)
            H_Z[row, right[term]] = one_F
        end
        for term in terms(simplify(phi(edge * a_anti_R2)).f)
            H_Z[row, left[term]] = one_F
        end
    end
    S.cache[:H_X] = S.cache[:X_stabs] = H_X
    S.cache[:H_Z] = S.cache[:Z_stabs] = H_Z
end

function X_stabilizers(S::FiniteGeneralized3DToricCode)
    _evaluate_3D_matrices!(S)
    return S.cache[:X_stabs]
end

function Z_stabilizers(S::FiniteGeneralized3DToricCode)
    _evaluate_3D_matrices!(S)
    return S.cache[:Z_stabs]
end

function stabilizers(S::FiniteGeneralized3DToricCode; standform::Bool = false)
    _evaluate_3D_matrices!(S)
    return invoke(stabilizers, Tuple{AbstractSubsystemCode}, S;
                  standform = standform)
end

character_vector(S::FiniteGeneralized3DToricCode) = S.cache[:char_vec]
defining_polynomials(
    S::Union{Generalized3DToricCode, FiniteGeneralized3DToricCode}
) = S.a, S.b
Laurent_polynomial_ring(
    S::Union{Generalized3DToricCode, FiniteGeneralized3DToricCode}
) = S.LR
twist_vectors(S::FiniteGeneralized3DToricCode) = S.a1, S.a2

function maximum_dimension(S::Generalized3DToricCode)
    R2 = Oscar._polyringquo(S.LR)
    R = codomain(R2)
    Q, _ = quo(R, ideal(R, R2.([S.a, S.b])))
    return 2 * vector_space_dimension(Q)
end
maximum_dimension(S::FiniteGeneralized3DToricCode) = S.k

function Base.show(io::IO, S::Generalized3DToricCode)
    print(io, "Generalized 3D toric-code datum over $(S.LR) with a = $(S.a), b = $(S.b)")
end

function Base.show(io::IO, S::FiniteGeneralized3DToricCode)
    kind = S.twisted ? "twisted " : ""
    print(io, "$(kind)finite generalized 3D toric code [[$(S.n), $(S.k)]]_$(order(S.F))")
end
