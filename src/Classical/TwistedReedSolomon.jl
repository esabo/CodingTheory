# Copyright (c) 2024 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
# constructors
#############################

"""
    TwistedReedSolomonCode(k::Int, α::Vector{T}, t::Vector{Int}, h::Vector{Int}, η::Vector{T}) where T <: CTFieldElem

Return the twisted Reed-Solomon code defined in `https://arxiv.org/abs/2107.06945`.
"""
function TwistedReedSolomonCode(
    k::Int,
    α::Vector{T},
    t::Vector{Int},
    h::Vector{Int},
    η::Vector{T},
) where {T<:CTFieldElem}

    l = length(t)
    l == length(h) ||
        throw(ArgumentError("Input vectors `t`, `h`, and `η` must have the same length"))
    l == length(η) ||
        throw(ArgumentError("Input vectors `t`, `h`, and `η` must have the same length"))
    length(unique(collect(zip(h, t)))) == l ||
        throw(ArgumentError("The tuples `(h[i], t[i])` must be distinct"))
    n = length(α)
    length(unique(α)) == n || throw(ArgumentError("The elements of `α` must be distinct"))
    1 ≤ k ≤ n ||
        throw(DomainError(k, "The dimension of the code must satisfy `1 ≤ k ≤ length(α)`"))
    all(1 ≤ x ≤ n - k for x in t) ||
        throw(DomainError(t, "Elements of `t` must satsify `1 ≤ t[i] ≤ n - k`"))
    all(0 ≤ x ≤ k - 1 for x in h) ||
        throw(DomainError(h, "Elements of `h` must satsify `0 ≤ h[i] ≤ k - 1`"))
    F = parent(α[1])
    all(parent(x) == F for x in α) ||
        throw(DomainError(α, "All elements of `α` must be over the same base ring"))
    all(parent(x) == F for x in η) || throw(
        DomainError(
            η,
            "All elements of `η` must be over the same base ring as the elements of `α`",
        ),
    )

    _, x = polynomial_ring(F, :x)
    G = zero_matrix(F, k, n)
    for i = 0:(k-1)
        g_i = x^i
        for j = 1:l
            h[j] == i && (g_i += η[j] * x^(k - 1 + t[j]);)
        end

        for c = 1:n
            G[i+1, c] = g_i(α[c])
        end
    end
    # display(G)

    t_dual = k .- h
    h_dual = (n - k) .- t
    H = zero_matrix(F, n - k, n)
    for i = 0:(n-k-1)
        g_i = x^i
        for j = 1:l
            h_dual[j] == i && (g_i += -η[j] * x^(n - k - 1 + t_dual[j]);)
        end

        for c = 1:n
            H[i+1, c] = g_i(α[c])
        end
    end
    # println(" ")
    # display(H)

    # println(" ")
    # display(G * transpose(H))

    C = LinearCode(G, H, false)
    return TwistedReedSolomonCode(
        C.F,
        C.n,
        C.k,
        C.d,
        C.l_bound,
        C.u_bound,
        C.G,
        C.H,
        C.G_stand,
        C.H_stand,
        C.P_stand,
        C.weight_enum,
        α,
        t,
        h,
        η,
        l,
    )
end

#############################
# getter functions
#############################

"""
    twist_vector(C::AbstractTwistedReedSolomonCode)

Return the twist vector of `C`.
"""
twist_vector(C::AbstractTwistedReedSolomonCode) = C.t

"""
    hook_vector(C::AbstractTwistedReedSolomonCode)

Return the hook vector of `C`.
"""
hook_vector(C::AbstractTwistedReedSolomonCode) = C.h

"""
    coefficient_vector(C::AbstractTwistedReedSolomonCode)

Return the coefficient vector of `C`.
"""
coefficient_vector(C::AbstractTwistedReedSolomonCode) = C.η

"""
    number_of_twists(C::AbstractTwistedReedSolomonCode)

Return the number of twists of `C`.
"""
number_of_twists(C::AbstractTwistedReedSolomonCode) = C.l

#############################
# setter functions
#############################

#############################
# general functions
#############################

# TODO check if all of these are in linear_code and branches on type
# TODO can do better wrt d, bounds, and weight enumerator
function dual(C::AbstractTwistedReedSolomonCode)
    return TwistedReedSolomonCode(
        C.F,
        C.n,
        C.n - C.k,
        missing,
        1,
        C.n,
        C.H,
        C.G,
        C.H_stand,
        C.G_stand,
        C.P_stand,
        missing,
        C.α,
        C.k .- C.h,
        (C.n - C.k) .- C.t,
        -C.η,
        C.l,
    )
end
