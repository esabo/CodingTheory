# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

using AbstractAlgebra
using Nemo

function ⊕(A::T, B::T) where T <: Union{fq_nmod_mat, gfp_mat}
    if base_ring(A) != base_ring(B)
        error("Matrices must be over the same base ring in directsum.")
    end

    return vcat(hcat(A, zero_matrix(base_ring(B), size(A, 1), size(B, 2))),
        hcat(zero_matrix(base_ring(A), size(B, 1), size(A, 2)), B))
end

function directsum(A::T, B::T) where T <: Union{fq_nmod_mat, gfp_mat}
    return A ⊕ B
end

function ⊗(A::T, B::T) where T <: Union{fq_nmod_mat, gfp_mat}
    if base_ring(A) != base_ring(B)
        error("Matrices in tensor product must be over the same ring.")
    end

    M = MatrixSpace(base_ring(A), size(A, 1) * size(B, 1), size(A, 2) * size(B, 2))
    C = M(0)
    if iszero(A) || iszero(B)
        return M(0)
    else
        (rB, cB) = size(B)
        for ar in 1:size(A, 1)
            for ac in 1:size(A, 2)
                C[(1 + (ar - 1) * rB):(ar * rB), (1 + (ac - 1) * cB):(ac * cB)] = A[ar, ac] * B
            end
        end
    end

    return C
end

function kron(A::T, B::T) where T <: Union{fq_nmod_mat, gfp_mat}
    return A ⊗ B
end

function tensorproduct(A::T, B::T) where T <: Union{fq_nmod_mat, gfp_mat}
    return A ⊗ B
end

# # to get primitive element out of field in Nemo
# function primitive_element(F::T; n_quo::Int = -1) where T <: Union{FqFiniteField, FqNmodFiniteField, Nemo.GaloisField, Nemo.GaloisFmpzField}
#     n = order(F)-1
#     k = fmpz(1)
#     if n_quo != -1
#         if !divisible(n, n_quo)
#             return F(1)
#         end
#         n, k = ppio(n, fmpz(n_quo))
#     end
#     primefactors_n = collect(keys(factor(n).fac))
#
#     x = rand(F)^Int(k)
#     while iszero(x)
#         x = rand(F)^Int(k)
#     end
#     while true
#         found = true
#         for l in primefactors_n
#             if isone(x^Int(div(n, l)))
#                 found = false
#                 break
#             end
#         end
#         if found
#             break
#         end
#         x = rand(F)^Int(k)
#         while iszero(x)
#             x = rand(F)^Int(k)
#         end
#     end
#     return x
# end
