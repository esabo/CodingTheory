# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

using AbstractAlgebra
using Nemo

import LinearAlgebra: tr

function ⊕(A::T, B::T) where T <: Union{fq_nmod_mat, gfp_mat}
    base_ring(A) == base_ring(B) || error("Matrices must be over the same base ring in directsum.")

    return vcat(hcat(A, zero_matrix(base_ring(B), size(A, 1), size(B, 2))),
        hcat(zero_matrix(base_ring(A), size(B, 1), size(A, 2)), B))
end

directsum(A::T, B::T) where T <: Union{fq_nmod_mat, gfp_mat} = A ⊕ B
⊗(A::T, B::T) where T <: Union{fq_nmod_mat, gfp_mat} = kronecker_product(A, B)
kron(A::T, B::T) where T <: Union{fq_nmod_mat, gfp_mat} = kronecker_product(A, B)
tensorproduct(A::T, B::T) where T <: Union{fq_nmod_mat, gfp_mat} = kronecker_product(A, B)
kroneckerproduct(A::T, B::T) where T <: Union{fq_nmod_mat, gfp_mat} = kronecker_product(A, B)
# nrows(A::T) where T = size(A, 1)
# ncols(A::T) where T = size(A, 2)

# I think we should avoid length checking here and return it for entire matrix if given
Hammingweight(v::T) where T <: Union{fq_nmod_mat, gfp_mat, Vector{S}} where S <: Integer = count(i->(i != 0), v)
weight(v::T) where T <: Union{fq_nmod_mat, gfp_mat, Vector{S}} where S <: Integer = Hammingweight(v)
wt(v::T) where T <: Union{fq_nmod_mat, gfp_mat, Vector{S}} where S <: Integer = Hammingweight(v)
Hammingdistance(u::T, v::T) where T <: Union{fq_nmod_mat, gfp_mat, Vector{S}} where S <: Integer = Hammingweight(u .- v)
distance(u::T, v::T) where T <: Union{fq_nmod_mat, gfp_mat, Vector{S}} where S <: Integer = Hammingweight(u .- v)
dist(u::T, v::T) where T <: Union{fq_nmod_mat, gfp_mat, Vector{S}} where S <: Integer = Hammingweight(u .- v)

function tr(x::fq_nmod, K::FqNmodFiniteField, verify::Bool=false)
    L = parent(x)
    if verify
        # shouldn't need Int casting here but just in case...
        Int64(characteristic(L)) == Int64(characteristic(K)) || error("The given field is not a subfield of the base ring of the element.")
        degree(L) % degree(K) == 0 || error("The given field is not a subfield of the base ring of the element.")
    end
    n = div(degree(L), degree(K))
    return sum([x^(q^i) for i in 0:(n - 1)])
end

function _expandelement(x::fq_nmod, K::FqNmodFiniteField, basis::Vector{fq_nmod}, verify::Bool=false)
    return [tr(x * i, K, verify) for i in basis]
end

function _expandrow(row::Vector{fq_nmod}, K::FqNmodFiniteField, basis::Vector{fq_nmod}, verify::Bool=false)
    return hcat([_expandelement(row[i], K, basis, verify)]...)
end

function expandmatrix(M::fq_nmod_mat, K::FqNmodFiniteField, basis::Vector{fq_nmod})
    L = base_ring(M)
    L != K || return M
    # shouldn't need Int casting here but just in case...
    Int64(characteristic(L)) == Int64(characteristic(K)) || error("The given field is not a subfield of the base ring of the element.")
    degree(L) % degree(K) == 0 || error("The given field is not a subfield of the base ring of the element.")
    n = div(degree(L), degree(K))
    n == length(basis) || error("Provided basis is of incorrect size for the given field and subfield.")
    # should really check if it is a basis
    return vcat([_expandrow(M[r, :], K, basis) for r in 1:size(M, 1)])
end

function symplecticinnerproduct(u::fq_nmod_mat, v::fq_nmod_mat)
    (size(u, 1) == 1 || size(u, 2) == 1) || error("First argument of symplectic inner product is not a vector: dims = $(size(u, 1)).")
    (size(v, 1) == 1 || size(v, 2) == 1) || error("Second argument of symplectic inner product is not a vector: dims = $(size(v, 1)).")
    length(u) == length(v) || error("Vectors must be the same length in symplectic inner product.")
    iseven(length(u)) || error("Vectors must have even length in symplectic inner product.")
    base_ring(u) == base_ring(v) || error("Vectors must be over the same field in symplectic inner product.")
    ncols = div(length(u), 2)
    return sum([u[i + ncols] * v[i] - v[i + ncols] * u[i] for i in 1:ncols])
end
# SIP(u::fq_nmod_mat, v::fq_nmod_mat) = symplecticinnerproduct(u, v)

# since changed method here, need to add error checks on dimensions
function aresymplecticorthogonal(A::fq_nmod_mat, B::fq_nmod_mat)
    AEuc = hcat(A[:, div(size(A, 2), 2) + 1:end], -A[:, 1:div(size(A, 2), 2)])
    iszero(AEuc * transpose(B)) || return false
    return true
end

# function traceinnerproduct(u::fq_nmod_mat, v::fq_nmod_mat)
#
# end
# TIP(u::fq_nmod_mat, v::fq_nmod_mat) = traceinnerproduct(u, v)

function Hermitianinnerproduct(u::fq_nmod_mat, v::fq_nmod_mat)
    (size(u, 1) == 1 || size(u, 2) == 1) || error("First argument of Hermitian inner product is not a vector: dims = $(size(u, 1)).")
    (size(v, 1) == 1 || size(v, 2) == 1) || error("Second argument of Hermitian inner product is not a vector: dims = $(size(v, 1)).")
    length(u) == length(v) || error("Vectors must be the same length in Hermitian inner product.")
    base_ring(u) == base_ring(v) || error("Vectors must be over the same field in Hermitian inner product.")
    q2 = order(base_ring(u))
    issquare(q2) || error("The Hermitian inner product is only defined over quadratic field extensions.")
    q = Int64(sqrt(q2))
    return sum([u[i] * v[i]^q for i in 1:length(u)])
end
# HIP(u::fq_nmod_mat, v::fq_nmod_mat) = Hermitianinnerproduct(u, v)

function Hermitianconjugatematrix(A::fq_nmod_mat)
    B = copy(A)
    q2 = order(base_ring(A))
    issquare(q2) || error("The Hermitian conjugate is only defined over quadratic field extensions.")
    q = Int64(sqrt(q2))
    return B .^ q
end

function printstringarray(A::Vector{String}, withoutIs=false)
    for a in A
        if !withoutIs
            println(a)
        else
            for i in a
                if i == 'I'
                    print(' ')
                else
                    print(i)
                end
            end
            print('\n')
        end
    end
end
printchararray(A::Vector{Vector{Char}}, withoutIs=false) = printstringarray(setchartostringarray(A), withoutIs)
printsymplecticarray(A::Vector{Vector{T}}, withoutIs=false) where T <: Integer = printstringarray(setsymplectictostringarray(A), withoutIs)

function entropy(x::Real)
    x != 0 || return 0
    (0 < x <= 1 - 1 / q) || error("Number should be in the range [0, 1 - 1/order(field)].")
    F = parent(x)
    q = order(F)
    return x * (log(q, q - 1) - log(q, x)) - (1 - x) * log(q, 1 - x)
end

function FpmattoJulia(M::fq_nmod_mat)
    degree(base_ring(M)) == 1 || error("Cannot promote higher order elements to the integers.")
    # Fp = [i for i in 0:Int64(characteristic(base_ring(M)))]
    A = zeros(Int64, size(M))
    for r in 1:size(M, 1)
        for c in 1:size(M, 2)
            # A[r, c] = Fp[findfirst(x->x==M[r, c], Fp)]
            A[r, c] = coeff(M[r, c], 0)
        end
    end
    return A
end

function pseudoinverse(M::fq_nmod_mat, verify::Bool=true)
    # let this fail elsewhere if not actually over a quadratic extension
    if degree(base_ring(M)) != 1
        M = deepcopy(quadratictosymplectic(M)')
    else
        M = deepcopy(M')
    end

    nr, nc = size(M)
    MS = MatrixSpace(base_ring(M), nr, nr)
    _, E = rref(hcat(M, MS(1)))
    E = E[:, (nc + 1):end]
    pinv = E[1:nc, :]
    dual = E[nc + 1:nr, :]

    if verify
        _, Mrref = rref(M)
        MScols = MatrixSpace(base_ring(M), nc, nc)
        E * M == Mrref || error("Pseudoinverse calculation failed (transformation incorrect).")
        Mrref[1:nc, 1:nc] == MScols(1) || error("Pseudoinverse calculation failed (failed to get I).")
        iszero(Mrref[nc + 1:nr, :]) || error("Pseudoinverse calculation failed (failed to get zero).")
        pinv * M == MScols(1) || error("Pseudoinverse calculation failed (eq 1).")
        M' * pinv' == MScols(1) || error("Pseudoinverse calculation failed (eq 2).")
        iszero(M' * dual') || error("Failed to correctly compute dual (rhs).")
        iszero(dual * M) || error("Failed to correctly compute dual (lhs).")
    end
    return pinv
end
