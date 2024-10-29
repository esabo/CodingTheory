# Copyright (c) 2024 David Marquis 
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
based on Algo 2.13 from K & S 
"""
struct SubsetGrayCode
    n::Int # length of codewords
    k::Int # weight of codewords
    len::Int # length of the iterator
end

function SubsetGrayCode(n::Int, k::Int)
    if 0 <= k <= n
        len = factorial(big(n)) ÷ (factorial(big(k)) * factorial(big(n - k)))
    else
        len = 0
    end
    return SubsetGrayCode(n, k, len)
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
    # The iterator's state is: (v, rank, inds)

    # v: an ordered collection of length G.k. It represents the nonzero indexs of a weight k length n vector w over F_2 
    # (This w is not part of this iterators state and we do not ever need to explicitly compute w)

    # rank: an int in [1, choose(G.n G.k)]. 
    # rank and v have a 1-1 correspondance and we can go between them with the functions _subset_unrank
    # and _subset_rank.

    # inds: an int array of length 3. It represents the indexs of w that were changed when updating the vector v. 
    # ind[i] can be either -1, meaning no index of w is stored at the position, or in [1, G.n].

    # Note that inds is part of the iterator's state inds only to prevent reallocation in each iteration.
    rank = Int(1)
    v = collect(1:G.k) 
    inds = Vector([-1, -1, -1])
    (inds, (v, rank, inds))
end

@inline function Base.iterate(G::SubsetGrayCode, state)
  # Based on Algorithm 2.13 in S & K
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
        while (i <= G.k) && (v[i] == i) 
          i = i+1
        end
        if Base.rem(G.k-i, 2) != 0 
            if (i == 1)
              v[1] != (v[1]-1) && _update_indexs!(inds, v[1], v[1]-1) 
              v[1] = (v[1]-1)
            else 
              v[i-1] != i && _update_indexs!(inds,v[i-1], i) 
              v[i-1] = i
              if i > 2
                v[i-2] != (i-1) && _update_indexs!(inds,v[i-2], i-1) 
                v[i-2] = (i-1)
              end
            end
        else 
          if i == G.k
            if ( G.n != v[i] )
              if (i > 1)
                v[i-1] != v[i] && _update_indexs!(inds,v[i-1], v[i]) 
                v[i-1] = v[i]
              end
              v[i] != (v[i]+1) && _update_indexs!(inds,v[i], v[i]+1) 
              v[i] = (v[i]+1)
            else
              v[i] != i && _update_indexs!(inds, v[i], i) 
              v[i] = i
            end
          else
            if ( v[i+1] != (v[i] + 1) )
              if i > 1
                v[i-1] != v[i] && _update_indexs!(inds, v[i-1], v[i]) 
                v[i-1] = v[i]
              end
              v[i] != (v[i]+1) && _update_indexs!(inds, v[i], v[i]+1) 
              v[i] = (v[i]+1)
            else
              v[i+1] != v[i] && _update_indexs!(inds, v[i+1], v[i]) 
              v[i+1] = v[i]
              v[i] != i && _update_indexs!(inds, v[i], i) 
              v[i] = i
            end
          end
        end
    end
    return (inds, (v, rank, inds))
end

@inline function rest(G::SubsetGrayCode, rank)
  #TODO there will be a future PR for integrating this (its used for multithreaded weight calculation)
  #=
  _subset_unrank(rank, G.n, vec)
  inds = Vector{Int}([-1,-1,-1]) 
  state = (vec, rank, inds)
  return Base.rest(G, state)
  =#
end 

function _update_indexs!(indexs::Vector{Int}, x::Int, y::Int)
  _update_indexs!(indexs, x)
  _update_indexs!(indexs, y)
  return nothing
end

function _update_indexs!(indexs::Vector{Int}, x::Int)
  if x == indexs[1]
    indexs[1] = -1
  elseif x == indexs[2]
    indexs[2] = -1
  elseif x == indexs[3]
    indexs[3] = -1
  elseif indexs[1] == -1
    indexs[1] = x
  elseif indexs[2] == -1
    indexs[2] = x
  elseif indexs[3] == -1
    indexs[3] = x
  else 
    throw("No index positions remaining")
  end
  return nothing
end

function _subset_rank(v::Vector{UInt}, k::UInt)
  # Based on algorithm 2.11 in S & K
  # Results are undefined if the entries of v arent in {1,..,n} for n>=k
  r = BigInt(0)
  s = BigInt(1)
  for i in k:-1:1
    r = r + extended_binomial(v[i], i) * s;
    s = -s;
  end
  if (k % 2) == 1
    r = r-1
  end
  return r
end

function _subset_unrank(r::BigInt, n::UInt, T::Vector{UInt})
  # Based on algorithm 2.12 in S & K
  k = length(T)
  subset_size_str="subset size k=($k) must be smaller than the set size n=($n)"
  k > n && throw(ArgumentError(subset_size_str))
  bnd = binomial(n, k)
  rank_size_str="rank must be in [0, choose(n, k)-1]=$(bnd)"
  r > bnd && throw(ArgumentError(rank_size_str))
  
  x = 0
  i = 0
  y = 0
  x = n
  for i::UInt in k:-1:1
    y = extended_binomial(x,i)
    while y > r
      x = x-1
      y = extended_binomial(x,i)
    end 
    T[i] = x+1
    r = extended_binomial(x+1,i) - r - 1
  end 
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
Base.in(v::Vector{Int}, G::GrayCode) = length(v) == G.n && count(v .!= 0) == G.k && view(v, 1:G.prefix_length) == G.prefix

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

"""
initally based on Algo 2 from Bitner, Ehrlich, and Reingold 1976. Modified to support prefix 
"""
@inline function Base.iterate(G::GrayCode, state)
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
