# Copyright (c) 2021 - 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

function _Reed_Muller_generator_matrix(r::Int, m::Int, alt::Bool=false)
    F = Oscar.Nemo.Native.GF(2)
    if r == 1 && m == 1 && !alt
        return matrix(F, 2, 2, [1, 1, 0, 1])
    elseif r == m
        return identity_matrix(F, 2^m)
    elseif r == 0
        return matrix(F, ones(Int, 1, 2^m))
    else
        Grm1 = _Reed_Muller_generator_matrix(r, m - 1, alt)
        Gr1m1 = _Reed_Muller_generator_matrix(r - 1, m - 1, alt)
        return vcat(hcat(Grm1, Grm1), hcat(zero_matrix(F, nrows(Gr1m1), ncols(Gr1m1)), Gr1m1))
    end
end

function generator_matrix(C::ReedMullerCode, stand_form::Bool = false)
    cache = getfield(C, :cache)
    if !haskey(cache, :G)
        cache[:G] = _Reed_Muller_generator_matrix(C.r, C.m, cache[:alt])
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

function parity_check_matrix(C::ReedMullerCode, stand_form::Bool = false)
    cache = getfield(C, :cache)
    if !haskey(cache, :H)
        H_mat = _Reed_Muller_generator_matrix(C.m - C.r - 1, C.m, cache[:alt])
        size(H_mat) == (C.n - C.k, C.k) && (H_mat = transpose(H_mat))
        cache[:H] = H_mat
    end
    if stand_form
        generator_matrix(C, true)
        return cache[:H_stand]
    end
    return cache[:H]
end

"""
$(TYPEDSIGNATURES)

Return the ``\\mathcal{RM}(r, m)`` Reed-Muller code.

# Notes
* If `alt` is `true`, the identity is used for the generator matrix for ``\\mathcal{RM}(1, 1)``, as in common in some sources.
  Otherwise, `[1 1; 0 1]` is used, as is common in other sources.
"""
function ReedMullerCode(r::Int, m::Int, alt::Bool=false)
    0 ≤ r < m || throw(DomainError((r, m), "Reed-Muller codes require 0 ≤ r < m."))
    m < 64 || throw(DomainError(m, "This Reed-Muller code requires the implementation of BigInts. Change if necessary."))

    F = Oscar.Nemo.Native.GF(2)
    n = 2^m
    k = sum(binomial(m, i) for i in 0:r)
    d = 2^(m - r)

    # Store the `alt` parameter in the cache so the lazy getter can access it
    cache = Dict{Symbol, Any}(:alt => alt)
    
    # We know the weight enumerator mathematically, inject it into the cache directly
    if r == 1
        # RM(1, m) weight enumerator: 1 weight 0, (2^{m+1} - 2) weight 2^{m-1}, 1 weight 2^m
        counts = Dict{Int, BigInt}(0 => 1, 2^(m - 1) => BigInt(2^(m + 1) - 2), n => 1)
        cache[:weight_enum] = HammingWeightEnumerator(n, counts)
    end

    return ReedMullerCode(F, n, k, d, d, d, r, m, cache)
end

"""
$(TYPEDSIGNATURES)

Return a random permuted Reed-Muller code `RM(r, m)`. 

# Notes
* Applies a random column permutation to the standard `RM(r, m)` generator matrix.
* Used in cryptographic settings to hide the affine geometric structure of the code.
"""
function RandomPermutedReedMullerCode(r::Int, m::Int)
    C_RM = ReedMullerCode(r, m)
    
    # Generate a random permutation of length 2^m
    n = C_RM.n
    rand_perm = randperm(n) # Requires `using Random`
    
    # We use our lazy permute_code function to safely shuffle the columns
    # without destroying the minimum distance bounds.
    return permute_code(C_RM, rand_perm)
end

# ==============================================================================
# BOOLEAN FUNCTIONS & COSETS
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the truth table (as a vector) of a random Boolean function in `m` variables 
of algebraic degree at most `r`.

# Notes
* This is mathematically equivalent to returning a random codeword from `RM(r, m)`.
"""
function RandomBooleanFunction(r::Int, m::Int)
    C = ReedMullerCode(r, m)
    G = generator_matrix(C)
    F = C.F
    
    # A random Boolean function of degree <= r is just a random linear combination
    # of the basis vectors of RM(r, m).
    coeffs = matrix(F, 1, C.k, [rand(F) for _ in 1:C.k])
    
    return coeffs * G
end

"""
$(TYPEDSIGNATURES)

Return a random coset representative of the `RM(r, m)` code.

# Notes
* When `r = 1`, the weight of the lowest-weight element in this coset 
  defines the nonlinearity of the representative Boolean function.
"""
function RandomRMCoset(r::Int, m::Int)
    C = ReedMullerCode(r, m)
    F = C.F
    
    # Generate a completely random Boolean function of unrestricted degree (up to m)
    # This acts as our error vector / coset leader.
    v = matrix(F, 1, C.n, [rand(F) for _ in 1:C.n])
    
    # In a full computational algebra system, you would return a formal Coset object here.
    # For now, we return the representative vector.
    return v
end

#############################
      # getter functions
#############################

"""
    order(C::ReedMullerCode)
    RM_r(C::ReedMullerCode)

Return the order, `r`, of the ``\\mathcal{RM}(r, m)`` Reed-Muller code.
"""
order(C::ReedMullerCode) = C.r
RM_r(C::ReedMullerCode) = order(C)

"""
    number_of_variables(C::ReedMullerCode)
    RM_m(C::ReedMullerCode)

Return the number of variables, `m`, of the ``\\mathcal{RM}(r, m)`` Reed-Muller code.
"""
number_of_variables(C::ReedMullerCode) = C.m
RM_m(C::ReedMullerCode) = C.m

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

