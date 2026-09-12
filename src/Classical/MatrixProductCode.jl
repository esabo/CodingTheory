# Copyright (c) 2022 - 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
      # lazy getters
#############################

function generator_matrix(C::MatrixProductCode, stand_form::Bool = false)
    cache = getfield(C, :cache)
    if !haskey(cache, :G)
        s, l = size(C.A)
        n_sub = C.C[1].n
        G = zero_matrix(C.F, C.k, l * n_sub)
        
        curr = 1
        for r in 1:s
            G_sub = generator_matrix(C.C[r])
            
            for c in 1:l
                G[curr:(curr + C.C[r].k - 1), (1 + (c - 1) * n_sub):(c * n_sub)] = C.A[r, c] * G_sub
            end
            curr += C.C[r].k
        end
        cache[:G] = G
    end
    
    if stand_form
        if !haskey(cache, :G_stand)
            G_stand, H_stand, P, _ = _standard_form(cache[:G])
            cache[:G_stand] = G_stand
            cache[:H_stand] = H_stand
            cache[:P_stand] = P
        end
        return cache[:G_stand]
    end
    return cache[:G]
end

#############################
        # constructors
#############################

"""
$(TYPEDSIGNATURES)

Return the matrix product code defined by the vector of linear codes `C` and matrix `A`.
Evaluates lazily.
"""
function MatrixProductCode(C::Vector{<:AbstractLinearCode}, A::CTMatrixTypes)
    isempty(C) && throw(ArgumentError("Vector of linear codes cannot be empty."))
    iszero(A) && throw(ArgumentError("Matrix A cannot be zero."))
    
    s, l = size(A)
    s == length(C) || throw(ArgumentError("Number of rows of A must be equal to the number of codes."))
    
    # Prevent overcomplete/degenerate codes that collapse the dimension
    rank(A) == s || throw(ArgumentError("Matrix A must have full row rank to prevent overcomplete codes."))
    
    F = C[1].F
    n_sub = C[1].n
    for i in 2:s
        F == C[i].F || throw(ArgumentError("All codes must have the same base ring."))
        n_sub == C[i].n || throw(ArgumentError("All codes must have the same length."))
    end
    F == base_ring(A) || throw(ArgumentError("Codes and matrix must have the same base ring."))

    k_new = sum(c.k for c in C)
    n_new = l * n_sub
    
    cache = Dict{Symbol, Any}()
    return MatrixProductCode(C, A, F, n_new, k_new, missing, 1, n_new, cache)
end

# ==============================================================================
# CRYPTOGRAPHIC GENERATORS
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return a random Matrix Product Code defined by the vector of linear codes `C` 
and a randomly generated full-rank `s × l` defining matrix `A`.
"""
function RandomMatrixProductCode(C::Vector{<:AbstractLinearCode}, l::Int)
    isempty(C) && throw(ArgumentError("Vector of linear codes cannot be empty."))
    s = length(C)
    l >= s || throw(DomainError(l, "The number of columns l must be >= s to ensure a full-rank defining matrix."))
    
    F = C[1].F
    
    # 1. Generate a random full-rank defining matrix
    A = zero_matrix(F, s, l)
    while true
        for i in 1:s
            for j in 1:l
                A[i, j] = rand(F)
            end
        end
        if rank(A) == s
            break
        end
    end
    
    # 2. Construct the code
    return MatrixProductCode(C, A)
end

#############################
      # getter functions
#############################

"""
$(TYPEDSIGNATURES)

Return the constituent linear codes used to construct the matrix product code.
"""
constituent_codes(C::MatrixProductCode) = C.C

"""
$(TYPEDSIGNATURES)

Return the defining matrix `A` of the matrix product code.
"""
defining_matrix(C::MatrixProductCode) = C.A

#############################
      # setter functions
#############################

#############################
     # general functions
#############################