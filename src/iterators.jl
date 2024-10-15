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

@inline function Base.iterate(G::SubsetGrayCode)
  # The invariants of the iterator's state are:
  # - v is an ordered collection of length G.k,
  # - rank is an int in [1, choose(G.n G.k)],
  # - inds is an array of length 3.
  # We only pass inds so that it is not allocated again in each iteration.
    rank = Int(1)
    v = collect(1:G.k) 
    inds = Vector([-1,-1,-1])
    (inds, (v, rank, inds))
end

@inline function Base.iterate(G::SubsetGrayCode, state)
    v, rank, inds = state

    if rank == G.len 
      return nothing 
    end
    rank += 1

    @inbounds begin
        # inds will store the set of indexs in the binary vector that were changed when updating the vector v
        inds[1] = -1
        inds[2] = -1
        inds[3] = -1

        i = 1
        while (i <= G.k) && (v[i] == i) 
          i = i+1
        end
        if Base.rem(G.k-i, 2) != 0 
            if (i == 1)
              v[1] != (v[1]-1) && update_tup!(inds, v[1], v[1]-1) 
              v[1] = (v[1]-1)
            else 
              v[i-1] != i && update_tup!(inds,v[i-1], i) 
              v[i-1] = i
              if i > 2
                v[i-2] != (i-1) && update_tup!(inds,v[i-2], i-1) 
                v[i-2] = (i-1)
              end
            end
        else 
          if i == G.k
            if ( G.n != v[i] )
              if (i > 1)
                v[i-1] != v[i] && update_tup!(inds,v[i-1], v[i]) 
                v[i-1] = v[i]
              end
              v[i] != (v[i]+1) && update_tup!(inds,v[i], v[i]+1) 
              v[i] = (v[i]+1)
            else
              v[i] != i && update_tup!(inds, v[i], i) 
              v[i] = i
            end
          else
            if ( v[i+1] != (v[i] + 1) )
              if i > 1
                v[i-1] != v[i] && update_tup!(inds, v[i-1], v[i]) 
                v[i-1] = v[i]
              end
              v[i] != (v[i]+1) && update_tup!(inds, v[i], v[i]+1) 
              v[i] = (v[i]+1)
            else
              v[i+1] != v[i] && update_tup!(inds, v[i+1], v[i]) 
              v[i+1] = v[i]
              v[i] != i && update_tup!(inds, v[i], i) 
              v[i] = i
            end
          end
        end
    end
    return (inds, (v, rank, inds))
end

@inline function rest(G::SubsetGrayCode, rank)
  _kSubsetRevDoorUnrank(rank, G.n, vec)
  inds = Vector{Int}([-1,-1,-1]) 
  state = (vec, rank, inds)
  return Base.rest(G, state)
end 

function update_tup!(ind_pairs::Vector{Int}, x::Int, y::Int)
  update_tup!(ind_pairs, x)
  update_tup!(ind_pairs, y)
  return nothing
end

function update_tup!(ind_pairs::Vector{Int}, x::Int)
  if x == ind_pairs[1]
    ind_pairs[1] = -1
  elseif x == ind_pairs[2]
    ind_pairs[2] = -1
  elseif x == ind_pairs[3]
    ind_pairs[3] = -1
  elseif ind_pairs[1] == -1
    ind_pairs[1] = x
  elseif ind_pairs[2] == -1
    ind_pairs[2] = x
  elseif ind_pairs[3] == -1
    ind_pairs[3] = x
  end
  return nothing
end

function _kSubsetRevDoorRank(v::Vector{UInt}, k::UInt)
  # usage: this can be an expensive function to call and should only be used for testing
  # Based on Algorithm 2.11 in S & K
  # Results are undefined if the entries of v arent in {1,..,k}
  r = UInt(0)
  s = 1 
  for i in k:-1:1
    r = r + extended_binomial(v[i], i) * s;
    s = -s;
  end
  if (k % 2) == 1
    r = r-1
  end
  return r
end

function _kSubsetRevDoorUnrank(r::UInt, n::UInt, T::Vector{UInt})
  k = length(T)
  subset_size_str="subset size k=($k) must be smaller than the set size n=($n)"
  k > n && throw(ArgumentError(subset_size_str))
  bnd = binomial(n, 2)
  rank_size_str="rank must be in [0, choose(n, r)-1]=$(bnd)"
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