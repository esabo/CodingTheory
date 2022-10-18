# Copyright (c) 2022 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

mutable struct MatrixProductCode <: AbstractMatrixProductCode
    F::Union{FqNmodFiniteField}
    n::Integer
    k::Integer
    d::Union{Integer, Missing}
    G::fq_nmod_mat
    Gorig::Union{fq_nmod_mat, Missing}
    H::fq_nmod_mat
    Horig::Union{fq_nmod_mat, Missing}
    Gstand::fq_nmod_mat
    Hstand::fq_nmod_mat
    weightenum::Union{WeightEnumerator, Missing}
    C::Vector{AbstractLinearCode}
    A::fq_nmod_mat
end

"""
    MatrixProductCode(C::Vector{AbstractLinearCode}, A::fq_nmod_mat)

Return the matrix product code defined by the vector of linear codes `C` and matrix `A`.
"""
# TODO: need to make sure codes are not overcomplete
function MatrixProductCode(C::Vector{AbstractLinearCode}, A::fq_nmod_mat)
    isempty(C) && throw(ArgumentError("Vector of linear codes cannot be empty."))
    iszero(A) && throw(ArgumentError("Matrix A cannot be zero."))
    s, l = size(A)
    s == length(C) || throw(ArgumentError("Number of rows of A must be equal to the number of codes."))
    F = C[1].F
    n = C[1].n
    for i in 2:s
        F == C[i].F || throw(ArgumentError("All codes must have the same base ring."))
        n == C[i].n || throw(ArgumentError("All codes must have the same length."))
    end
    F == base_ring(A) || throw(ArgumentError("Codes and matrix must have the same base ring."))

    # H formula holds in special case only
    G = zero_matrix(F, sum([C[i].k for i in 1:s], l * n))
    # H = zero_matrix(F, sum([nrows(C[i].H) for i in 1:s], l * n))
    curr = 1
    # currH = 1
    # need to do in this row/column order
    for r in 1:s
        for c in 1:l
            G[curr:curr + C[r].k, 1 + (r - 1) * n:r * n] = A[r, c] * generatormatrix(C[r])
            # H[currH:currH + nrows(C[r].H), 1 + (r - 1) * n:r * n] = A[r, c] * paritycheckmatrix(C[r])
        end
        curr += C[r].k
        # currH += nrows(C[r].H)
    end
    
    _, H = right_kernel(G)
    Gstand, Hstand = _standardform(G)
    return MatrixProductCode(F, n, nrows(G), missing, G, missing, transpose(H), missing, Gstand,
        Hstand, missing, C, A)
end

"""
    dual(C::MatrixProductCode)

Return the dual of `C`. If the dual is also a matrix product code, a matrix product code
will be returned.
"""
function dual(C::MatrixProductCode)
    nr, nc = size(C.A)
    # probably not going to work
    nr == nc || return LinearCode.dual(C)
    D = Vector{LinearCode}()
    for i in 1:length(C.C)
        push!(D, dual(C.C[i]))
    end
    
    try
        Ainv = inv(C.A)
    catch
        return LinearCode.dual(C)
    end
    return MatrixProductCode(D, transpose(Ainv))
end

"""
Hermitiandual(C::MatrixProductCode)

Return the Hermitian dual of `C`. If the dual is also a matrix product code, a matrix product code
will be returned.
"""
function Hermitiandual(C::MatrixProductCode)
    # the inner functions here should complain if not quadratic, so don't have to check here
    nr, nc = size(C.A)
    # probably not going to work
    nr == nc || return LinearCode.Hermitiandual(C)
    D = Vector{LinearCode}()
    for i in 1:length(C.C)
        push!(D, Hermitiandual(C.C[i]))
    end

    try
        Ainv = inv(Hermitianconjugatematrix(C.A))
    catch
        return LinearCode.dual(C)
    end
    return MatrixProductCode(D, transpose(Ainv))
end
