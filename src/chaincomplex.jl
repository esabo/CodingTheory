# Copyright (c) 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    ChainComplex(chain::Vector{T}) where T <: CTMatrixTypes

Return a chain complex based on the boundary maps in `chain`.
"""
function ChainComplex(chain::Vector{T}) where T <: CTMatrixTypes
    F = base_ring(chain[1])
    all(x -> base_ring(x) == F, chain) || throw(ArgumentError("All inputs must be over the same base ring"))
    println([size(mat) for mat in chain])

    # check to make sure matrices define a valid chain complex
    len = length(chain)
    # this is one place people may actually not index linearly
    for i in 1:len - 1
        iszero(chain[i] * chain[i + 1]) || throw(ArgumentError("Boundary maps must satisfy `iszero(chain[i] * chain[i + 1])`"))
    end
    return ChainComplex{T}(F, len, chain)
end

"""
    ChainComplex(F::CTFieldTypes, len::Integer, boundaries::Vector{T}) where T <: CTMatrixTypes

Return the chain complex based on the input data.
"""
ChainComplex(F::CTFieldTypes, len::Integer, boundaries::Vector{T}) where T <: CTMatrixTypes = ChainComplex{T}(F, len, boundaries)

"""
    ChainComplex(S::AbstractStabilizerCode)

Return the 3-term chain complex associated with the CSS code `S`.
"""
ChainComplex(S::AbstractStabilizerCode) = ChainComplex(CSSTrait(typeof(S)), S)
ChainComplex(::IsCSS, S::AbstractStabilizerCode) = ChainComplex(S.F, 2, [S.X_stabs, transpose(S.Z_stabs)])
ChainComplex(::IsNotCSS, S::AbstractStabilizerCode) = throw(ArgumentError("This is only defined for CSS codes"))

"""
    ChainComplex(C::AbstractLinearCode)

Return the 2-term chain complex associated with the parity-check matrix of `C`.
"""
ChainComplex(C::AbstractLinearCode) = ChainComplex([parity_check_matrix(C)])

#############################
      # getter functions
#############################

"""
    field(chain::ChainComplex)

Return the base ring of the boundary maps of the chain complex.
"""
field(chain::ChainComplex) = chain.F

"""
    length(chain::ChainComplex)

Return the number of spaces in the chain complex.

* Note
- This is the number of boundary maps plus one.
"""
length(chain::ChainComplex) = chain.length

"""
    boundaries(chain::ChainComplex)

Return the boundary maps of the chain complex.

* Note
- boundaries(chain)[i] is the map $∂_i: C_i \to C_{i - 1}$.
"""
boundaries(chain::ChainComplex) = chain.boundaries

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

"""
    cochain(chain::ChainComplex)

Return the dual of the chain complex.
"""
cochain(chain::ChainComplex) = ChainComplex([transpose(∂) for ∂ in reverse(chain.boundaries)])

"""
    ⊗(chain_A::ChainComplex{T}, chain_B::ChainComplex{T}) where T <: CTMatrixTypes
    tensor_product(chain_A::ChainComplex{T}, chain_B::ChainComplex{T}) where T <: CTMatrixTypes

Return the total complex of the tensor product of the two chains.
"""
function ⊗(chain_A::ChainComplex{T}, chain_B::ChainComplex{T}) where T <: CTMatrixTypes
    chain_A.F == chain_B.F || throw(ArgumentError("Both chains must be over the same base ring"))

    identity_maps_A = [[identity_matrix(chain_A.F, nrows(mat)) for mat in chain_A.boundaries]; identity_matrix(chain_A.F, ncols(chain_A.boundaries[end]))]
    identity_maps_B = [[identity_matrix(chain_A.F, nrows(mat)) for mat in chain_B.boundaries]; identity_matrix(chain_A.F, ncols(chain_B.boundaries[end]))]
    boundaries = Vector{T}()
    n_max = chain_A.length + chain_B.length

    ns = [[(i, n - 1 - i) for i in 0:chain_A.length if n - 1 - i in 0:chain_B.length] for n in 1:chain_A.length + chain_B.length]
    for (n, indices) in reverse(collect(enumerate(ns)))
        num_cols = n == n_max ? ncols(chain_A.boundaries[end]) * ncols(chain_B.boundaries[end]) : nrows(boundaries[1])
        ∂ = zero_matrix(chain_A.F, 0, num_cols)
        for (i, j) in indices
            temp = if i == chain_A.length
                temp2 = (-1)^i * identity_maps_A[i + 1] ⊗ chain_B.boundaries[j + 1]
                hcat(zero_matrix(chain_A.F, nrows(temp2), num_cols - ncols(temp2)), temp2)
            elseif j == chain_B.length
                temp2 = chain_A.boundaries[i + 1] ⊗ identity_maps_B[j + 1]
                hcat(temp2, zero_matrix(chain_A.F, nrows(temp2), num_cols - ncols(temp2)))
            else
                temp2 = chain_A.boundaries[i + 1] ⊗ identity_maps_B[j + 1]
                temp3 = (-1)^i * identity_maps_A[i + 1] ⊗ chain_B.boundaries[j + 1]
                num_cols1 = 0
                for l in j + 2:chain_B.length
                    k = n - l + 1
                    if k in eachindex(chain_A.boundaries)
                        num_cols1 += nrows(chain_A.boundaries[k]) * ncols(chain_B.boundaries[l])
                    end
                end
                num_cols2 = num_cols - num_cols1 - ncols(temp2) - ncols(temp3)
                hcat(zero_matrix(chain_A.F, nrows(temp2), num_cols1), temp3, temp2, zero_matrix(chain_A.F, nrows(temp2), num_cols2))
            end
            ∂ = vcat(∂, temp)
        end
        pushfirst!(boundaries, ∂)
    end

    return ChainComplex(boundaries)
end
tensor_product(chain_A::ChainComplex, chain_B::ChainComplex) = ⊗(chain_A, chain_B)
