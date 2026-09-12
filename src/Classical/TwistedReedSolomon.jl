# Copyright (c) 2024 - 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
      # lazy getters
#############################

function generator_matrix(C::TwistedReedSolomonCode, stand_form::Bool = false)
    cache = getfield(C, :cache)
    if !haskey(cache, :G)
        _, x = polynomial_ring(C.F, :x)
        G = zero_matrix(C.F, C.k, C.n)
        
        for i in 0:(C.k - 1)
            g_i = x^i
            for j in 1:C.l
                C.h[j] == i && (g_i += C.η[j] * x^(C.k - 1 + C.t[j]))
            end
            
            for c in 1:C.n
                G[i + 1, c] = g_i(C.α[c])
            end
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

function parity_check_matrix(C::TwistedReedSolomonCode, stand_form::Bool = false)
    cache = getfield(C, :cache)
    if !haskey(cache, :H)
        _, x = polynomial_ring(C.F, :x)
        
        t_dual = C.k .- C.h
        h_dual = (C.n - C.k) .- C.t
        H = zero_matrix(C.F, C.n - C.k, C.n)
        
        for i in 0:(C.n - C.k - 1)
            g_i = x^i
            for j in 1:C.l
                h_dual[j] == i && (g_i += -C.η[j] * x^(C.n - C.k - 1 + t_dual[j]))
            end
            
            for c in 1:C.n
                H[i + 1, c] = g_i(C.α[c])
            end
        end
        cache[:H] = H
    end
    
    if stand_form
        generator_matrix(C, true)
        return cache[:H_stand]
    end
    return cache[:H]
end

#############################
        # constructors
#############################

"""
$(TYPEDSIGNATURES)

Return the twisted Reed-Solomon code defined in `https://arxiv.org/abs/2107.06945`.
Evaluates lazily.
"""
function TwistedReedSolomonCode(k::Int, α::Vector{T}, t::Vector{Int}, h::Vector{Int}, η::Vector{T}) where T <: CTFieldElem
    l = length(t)
    l == length(h) || throw(ArgumentError("Input vectors `t`, `h`, and `η` must have the same length"))
    l == length(η) || throw(ArgumentError("Input vectors `t`, `h`, and `η` must have the same length"))
    length(unique(collect(zip(h, t)))) == l || throw(ArgumentError("The tuples `(h[i], t[i])` must be distinct"))
    
    n = length(α)
    length(unique(α)) == n || throw(ArgumentError("The elements of `α` must be distinct"))
    1 ≤ k ≤ n || throw(DomainError(k, "The dimension of the code must satisfy `1 ≤ k ≤ length(α)`"))
    all(1 ≤ x ≤ n - k for x in t) || throw(DomainError(t, "Elements of `t` must satisfy `1 ≤ t[i] ≤ n - k`"))
    all(0 ≤ x ≤ k - 1 for x in h) || throw(DomainError(h, "Elements of `h` must satisfy `0 ≤ h[i] ≤ k - 1`"))
    
    F = parent(α[1])
    all(parent(x) == F for x in α) || throw(DomainError(α, "All elements of `α` must be over the same base ring"))
    all(parent(x) == F for x in η) || throw(DomainError(η, "All elements of `η` must be over the same base ring as the elements of `α`"))

    cache = Dict{Symbol, Any}()
    return TwistedReedSolomonCode(F, n, k, missing, 1, n - k + 1, α, t, h, η, l, cache)
end

"""
$(TYPEDSIGNATURES)

Return a random Twisted Reed-Solomon code of length `n`, dimension `k`, and `l` twists over `F`.
"""
function RandomTwistedReedSolomonCode(F::CTFieldTypes, n::Int, k::Int, l::Int)
    1 <= k <= n || throw(DomainError((k, n), "Require 1 <= k <= n."))
    n <= Int(order(F)) || throw(DomainError(n, "Length cannot exceed the order of the field."))
    
    # 1. Distinct evaluation points
    α = elem_type(F)[]
    while length(α) < n
        pt = rand(F)
        !(pt in α) && push!(α, pt)
    end
    
    # 2. Generate valid unique (h, t) pairs
    # Constraints: 1 <= t <= n-k, 0 <= h <= k-1
    possible_pairs = [(H, T) for H in 0:(k-1) for T in 1:(n-k)]
    l <= length(possible_pairs) || throw(DomainError(l, "Too many twists requested for the given dimension constraints."))
    
    # Sample without replacement
    chosen_pairs = Tuple{Int, Int}[]
    while length(chosen_pairs) < l
        pair = rand(possible_pairs)
        !(pair in chosen_pairs) && push!(chosen_pairs, pair)
    end
    
    h = [p[1] for p in chosen_pairs]
    t = [p[2] for p in chosen_pairs]
    
    # 3. Non-zero twist coefficients
    η = elem_type(F)[]
    while length(η) < l
        s = rand(F)
        !iszero(s) && push!(η, s)
    end
    
    return TwistedReedSolomonCode(k, α, t, h, η)
end

#############################
      # getter functions
#############################

"""
$(TYPEDSIGNATURES)

Return the twist vector of `C`.
"""
twist_vector(C::AbstractTwistedReedSolomonCode) = C.t

"""
$(TYPEDSIGNATURES)

Return the hook vector of `C`.
"""
hook_vector(C::AbstractTwistedReedSolomonCode) = C.h

"""
$(TYPEDSIGNATURES)

Return the coefficient vector of `C`.
"""
coefficient_vector(C::AbstractTwistedReedSolomonCode) = C.η

"""
$(TYPEDSIGNATURES)

Return the number of twists of `C`.
"""
number_of_twists(C::AbstractTwistedReedSolomonCode) = C.l

#############################
     # general functions
#############################

"""
$(TYPEDSIGNATURES)

Return the dual of the Twisted Reed-Solomon code `C`. 
This operation is `O(1)` and perfectly tracks the dual's twists, hooks, and coefficients.
"""
function dual(C::AbstractTwistedReedSolomonCode)
    t_dual = C.k .- C.h
    h_dual = (C.n - C.k) .- C.t
    η_dual = -C.η
    k_dual = C.n - C.k
    
    # Optional: If matrices are cached, we can safely swap them
    cache = Dict{Symbol, Any}()
    old_cache = getfield(C, :cache)
    haskey(old_cache, :H) && (cache[:G] = old_cache[:H])
    haskey(old_cache, :G) && (cache[:H] = old_cache[:G])
    
    return TwistedReedSolomonCode(C.F, C.n, k_dual, missing, 1, C.k + 1, C.α, t_dual, h_dual, η_dual, C.l, cache)
end
