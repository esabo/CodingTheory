# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

using AbstractAlgebra
using Nemo

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

# function hamming_distance(w₁::T, w₂::T) where T <: AbstractWord
#     if ! isequal(length(w₁), length(w₂))
#         throw(error("Cannot compute Hamming Distance on strings of unequal length."))
#     end
#
#     distance = 0
#
# 	for (s₁, s₂) in zip(w₁, w₂)
# 		if s₁ ≠ s₂
# 			distance += 1
# 		end
# 	end
#
#     return distance
# end
#
# """
# ```julia
# sizeof_perfect_code(q::Number, n::Number, d::Number) -> Number
# ```
# Calculates the number of gigabytes required to store a perfect code of parameters q, n, and d.
# """
# function sizeof_perfect_code(q::Int, n::Int, d::Int)
# 	return (sizeof(ntuple(_ -> gensym(), n)) * hamming_bound(big.([q, n, d])...)) / (2^30)
# end
# sizeof_perfect_code(q::Number, n::Number, d::Number) = sizeof_perfect_code(round.(BigInt, [q, n, d])...)
#
# """
# ```julia
# sizeof_all_words(q::Number, n::Number) -> Number
# ```
# Calculates the number of gigabytes required to store all unique words of length n from an alphabet of size q.
# """
# function sizeof_all_words(q::Int, n::Int)
# 	return (sizeof(ntuple(_ -> gensym(), n)) * big(q)^n) / (2^30)
# end
# sizeof_all_words(q::Number, n::Number) = sizeof_all_words(round.(BigInt, (q, n))...)
#
# """
# ```julia
# get_codewords(G::AbstractArray, m::Int) -> Codewords{M}
# ```
# Get codewords of a code from the _generating matrix_ under a finite field of modulo `m`.  Precisely, computes all linear combinations of the rows of the generating matrix.
#
# Parameters:
#   - `G::AbstractArray`: A matrix of Ints which generates the code.
#   - `m::Int`: The bounds of the finite field (i.e., the molulus you wish to work in).
#
# Returns:
#   - `Codewords{M}`: An array of codewords, each of length `M`.  Each codewords is a tuple, and each character in said word is a symbol.
# """
#
# function get_codewords(G::AbstractArray, m::Int)
# 	codewords = Vector()
# 	rows = Vector(undef, size(G, 2))
#
# 	for i in 1:size(G, 1)
# 		rows[i] = [G[i, j] for j in 1:size(G, 2)]
# 	end
#
# 	for c in Base.Iterators.product([0:m-1 for _ in 1:size(G, 1)]...)
# 		word = Ref(c[1]) .* rows[1]
#
# 		for i in 2:size(G, 1)
# 			word = mod.(word .+ (Ref(c[i]) .* rows[i]), m)
# 		end
#
# 		push!(codewords, word)
# 	end
#
# 	return codewords
# end

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
