# Copyright (c) 2022-2024 Benjamin Ide, David Marquis 
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

struct SubsetGrayCode
    # details about the iterator using this struct are in the comments in the iterate() function below
    n::Int # length of codewords
    k::Int # weight of codewords
    len::BigInt # length of the iterator
    init_rank::BigInt # iteration will start with the subset of this rank
end

function SubsetGrayCode(n::Int, k::Int)
    init_rank = big(1)
    num_threads = 1
    return _subset_gray_code_from_num_threads(n, k, init_rank, num_threads)
end

function _subset_gray_code_from_num_threads(n::Int, k::Int, init_rank::BigInt, num_threads::Int)
    bin = extended_binomial(n, k)
    if bin % num_threads != 0 
        throw(ArgumentError("num_threads=$num_threads, must divide binomial($n $k)"))
    end
    len = fld(bin, num_threads)
    return SubsetGrayCode(n, k, len, init_rank)
end

function _subset_gray_codes_from_num_threads(n::Int, k::Int, num_threads::Int)
    # This function splits a single iterator into several pieces. 
    # The intended usage is to do the iteration with multiple threads
    bin = extended_binomial(n, k)
    i = 0
    len = fld(bin, num_threads)

    itrs = fill(SubsetGrayCode(1,1,1,1), num_threads)
    for i in 0:num_threads - 1
        init_rank = 1 + i * len
        itrs[i+1] = SubsetGrayCode(n, k, len, init_rank)
    end
    return itrs
end

Base.IteratorEltype(::SubsetGrayCode) = Base.HasEltype()
Base.eltype(::SubsetGrayCode) = Array{Int, 1}
Base.IteratorSize(::SubsetGrayCode) = Base.HasLength()
@inline function Base.length(G::SubsetGrayCode)
    return G.len
end
Base.in(v::Vector{Int}, G::SubsetGrayCode) = length(v) == G.k 

@inline function Base.isdone(G::SubsetGrayCode, state) 
    (_, rank, _) = state
    return rank == G.len
end

@inline function Base.iterate(G::SubsetGrayCode)
    #=
    The iterator's state is: (v, rank, inds)

    v: a subset of k elements which we represent as an ordered vector. 
    The intended usage for this iterator is to consider v as the nonzero indexs of a weight k, 
    length n, vector w over F_2 (Note that w is not part of the state and we do not ever need to 
    explicitly compute w).

    rank: an int in [1, choose(G.n G.k)]. 
    rank and v have are in one-to-one correspondance and we can go between them with the functions 
    _subset_unrank and _subset_rank.

    inds: an int array of length 3. It represents the indexs of w that were changed when updating 
    the vector v. 
    ind[i] can be either and index of w (and so in [1, G.n]) or be -1, meaning no index of w is 
    stored at the position.

    Note that inds is part of the iterator's state inds only to prevent reallocation in each 
    iteration.
    =#
    # subset_vec = zeros(Int, G.k)
    # CodingTheory._subset_unrank!(G.init_rank, UInt(G.n), subset_vec)

    # a = @allocated begyyin
    subset_vec = zeros(Int, G.k)
    if G.init_rank == 1
        subset_vec = Int.(collect(1 : G.k))
        # first_vec = vcat(ones(Int, r), zeros(Int, k - r))
        # for i in 1:G.k
        #     subset_vec[i] = 1
        # end
    else
        CodingTheory._subset_unrank!(G.init_rank, UInt(G.n), subset_vec)
    end
    # end; a > 0 && @show a

    inds = fill(-1, 3) 
    (inds, (subset_vec, 1, inds))
end

@inline function Base.iterate(G::SubsetGrayCode, state)
    # Based on Algorithm 2.13 in kreher1999combinatorial
    v, rank, inds = state

    if rank == G.len
        return nothing 
    end
    rank += 1

    @inbounds begin
        inds[1] = -1
        inds[2] = -1
        inds[3] = -1

        i = 1
        while i <= G.k && v[i] == i
            i = i + 1
        end
        if Base.rem(G.k - i, 2) != 0
            if i == 1
                v[1] != (v[1] - 1) && _update_indices!(inds, v[1], v[1] - 1) 
                v[1] = (v[1] - 1)
            else 
                v[i - 1] != i && _update_indices!(inds,v[i - 1], i) 
                v[i - 1] = i
                if i > 2
                    v[i - 2] != (i - 1) && _update_indices!(inds,v[i - 2], i - 1) 
                    v[i - 2] = (i - 1)
                end
            end
        else 
            if i == G.k
                if G.n != v[i]
                    if (i > 1)
                      v[i - 1] != v[i] && _update_indices!(inds,v[i - 1], v[i]) 
                      v[i - 1] = v[i]
                    end
                    v[i] != (v[i] + 1) && _update_indices!(inds,v[i], v[i] + 1) 
                    v[i] = (v[i] + 1)
                else
                    v[i] != i && _update_indices!(inds, v[i], i) 
                    v[i] = i
                end
            else
                if v[i + 1] != (v[i] + 1)
                    if i > 1
                        v[i - 1] != v[i] && _update_indices!(inds, v[i - 1], v[i]) 
                        v[i - 1] = v[i]
                    end
                    v[i] != (v[i] + 1) && _update_indices!(inds, v[i], v[i] + 1) 
                    v[i] = (v[i] + 1)
                else
                    v[i + 1] != v[i] && _update_indices!(inds, v[i + 1], v[i]) 
                    v[i + 1] = v[i]
                    v[i] != i && _update_indices!(inds, v[i], i) 
                    v[i] = i
                end
            end
        end
    end
    return (inds, (v, rank, inds))
end

function _update_indices!(indices::Vector{Int}, x::Int, y::Int)
    _update_indices!(indices, x)
    _update_indices!(indices, y)
    return nothing
end

function _update_indices!(indices::Vector{Int}, x::Int)
    if x == indices[1]
        indices[1] = -1
    elseif x == indices[2]
        indices[2] = -1
    elseif x == indices[3]
        indices[3] = -1
    elseif indices[1] == -1
        indices[1] = x
    elseif indices[2] == -1
        indices[2] = x
    elseif indices[3] == -1
        indices[3] = x
    else 
        error("No index positions remaining")
    end
    return nothing
end

function _subset_rank(v::Vector{UInt}, k::UInt)
    #=
    Based on Algorithm 2.11 in kreher1999combinatorial
    Results are undefined if the entries of v arent in {1,..,n} for n>=k
    =#
    r = BigInt(0)
    s = BigInt(1)
    for i in k:-1:1
        r = r + extended_binomial(v[i], i) * s
        s = -s
    end
    if (k % 2) == 1
        r = r - 1
    end
    r = r + 1 # we use 1-indexed ranks to follow the usual Julia convention
    return r
end

function _subset_unrank_to_vec!(r::BigInt, k::UInt, vec::Vector{Int})
    # returns a {0,1} vector with nonzero entries corresponding to the subset returned by _subset_unrank!
    n = length(vec)
    subset_vec = Int.(zeros(k))
    CodingTheory._subset_unrank!(r, UInt(n), subset_vec)
    for i in subset_vec
        vec[i] = 1
    end
end

function _subset_unrank!(r::BigInt, n::UInt, T::Vector{Int})
    # Based on Algorithm 2.12 in kreher1999combinatorial
    r = r - 1 

    k = length(T)
    subset_size_str::String = "subset size k = $k must be smaller than the set size n = $n"
    k > n && throw(ArgumentError(subset_size_str))
    # count_bins = 0
    # bnd = extended_binomial(n, k)
    # # count_bins += 1
    # rank_size_str::String = "rank must be in [0, choose(n, k) - 1] = $bnd but is $(r)"
    # # r > bnd && throw(ArgumentError(rank_size_str))
    # if r > bnd 
    #     println(rank_size_str)
    #     println("bug exit _subset_unrank")
    #     return
    # end
  
    x = 0
    i = 0
    y = 0 
    x = n
    for i::UInt in k:-1:1
        y = extended_binomial(x, i)
        # println((Int16.(x),Int16.(i)))
        # count_bins += 1
        while y > r
           x = x - 1
           y = extended_binomial(x, i)
        #    println((x,i))
        #    count_bins += 1
        end 
        T[i] = x + 1
        r = extended_binomial(x + 1, i) - r - 1
        # println((x+1,i))
        # count_bins += 1
    end 
    # println("num calls to binomial(a,b) in 1 unrank = $(count_bins)")
end

struct GrayCode
    n::Int # length of codewords
    k::Int # weight of codewords
    ns::Int # n for the subcode
    ks::Int # k for the subcode
    prefix::Vector{Int}
    prefix_length::Int
    mutate::Bool
end

GrayCode(n::Int, k::Int; mutate::Bool = false) = GrayCode(n, k, Int[]; mutate = mutate)

function GrayCode(n::Int, k::Int, prefix::Vector{Int}; mutate::Bool = false)
    GrayCode(n, k, n - length(prefix), k - count(prefix .!= 0), prefix,
        length(prefix), mutate)
end

Base.IteratorEltype(::GrayCode) = Base.HasEltype()
Base.eltype(::GrayCode) = Array{Int, 1}
Base.IteratorSize(::GrayCode) = Base.HasLength()
@inline function Base.length(G::GrayCode)
    if 0 <= G.ks <= G.ns
        factorial(big(G.ns)) ÷ (factorial(big(G.ks)) * factorial(big(G.ns - G.ks)))
    else
        0
    end
end
Base.in(v::Vector{Int}, G::GrayCode) = length(v) == G.n && count(v .!= 0) == G.k && (view(v, 1:G.prefix_length) == G.prefix)

@inline function Base.iterate(G::GrayCode)
    0 <= G.ks <= G.ns || return nothing

    g = [i <= G.ks ? 1 : 0 for i in 1:G.ns + 1]
    τ = collect(2:G.ns + 2)
    τ[1] = G.ks + 1
    # to force stopping with returning the only valid vector when ks == 0 and ns > 0
    iszero(G.ks) && (τ[1] = G.ns + 1;)
    v = [G.prefix; g[end - 1:-1:1]]
    ((G.mutate ? v : copy(v);), (g, τ, G.ks, v))
end

@inline function Base.iterate(G::GrayCode, state)
    #=
    based on Algorithm 2 from bitner1976efficient and modified to support iterating from a binary 
    prefix vector 
    =#
    g, τ, t, v = state
    @inbounds begin
        i = τ[1]
        i < G.ns + 1 || return nothing

        τ[1] = τ[i]
        τ[i] = i + 1
        if g[i] == 1
            if t != 0
                g[t] = g[t] == 0 ? 1 : 0
                if t < G.ns + 1
                    v[G.prefix_length + G.ns + 1 - t] = g[t]
                end
            else
                g[i - 1] = g[i - 1] == 0 ? 1 : 0
                if i - 1 < G.ns + 1
                    v[G.prefix_length + G.ns + 2 - i] = g[i - 1]
                end
            end
            t = t + 1
        else
            if t != 1
                g[t - 1] = g[t - 1] == 0 ? 1 : 0
                if t - 1 < G.ns + 1
                    v[G.prefix_length + G.ns + 2 - t] = g[t - 1]
                end
            else
                g[i - 1] = g[i - 1] == 0 ? 1 : 0
                if i - 1 < G.ns + 1
                    v[G.prefix_length + G.ns + 2 - i] = g[i - 1]
                end
            end
            t = t - 1
        end

        g[i] = g[i] == 0 ? 1 : 0
        if i < G.ns + 1
            v[G.prefix_length + G.ns + 1 - i] = g[i]
        end

        if t == i - 1 || t == 0
            t = t + 1
        else
            t = t - g[i - 1]
            τ[i - 1] = τ[1]
            if t == 0
                τ[1] = i - 1
            else
                τ[1] = t + 1
            end
        end
    end

    ((G.mutate ? v : copy(v);), (g, τ, t, v))
end
