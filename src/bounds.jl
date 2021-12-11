"""
```julia
rate(q::Integer, M::Integer, n::Integer) -> Real
```

Calculate the rate of a code.  That is, how efficient the code is.

Parameters:
  - `q::Integer`: the number of symbols in the code.
  - `M::Integer`: the size/number of elements in the code.
  - `n::Integer`: The word length.

Returns:
  - `Real`: Rate of the code.

---

### Examples

```julia
julia> rate(3, 5, 4) # the rate of the code which has 3 symbols, 5 words in the code, and word length of 4 (e.g., Σ = {A, B, C}, C = {ABBA,CABA,BBBB,CAAB,ACBB})
0.3662433801794817
```
"""
rate(q::T, M::T, n::T) where {T <: Integer} = log(q, M) / n

__spheres(q::T, n::T, r::T) where {T <: Integer} = sum(Integer[((big(q) - 1)^i) * binomial(big(n), big(i)) for i in 0:r])
__sphere_bound(round_func::Function, q::T, n::T, d::T) where {T <: Integer} = round_func((big(q)^n) / __spheres(q, n, d))

"""
```julia
sphere_covering_bound(q::Integer, n::Integer, d::Integer) -> Integer
```

Computes the sphere covering bound of a ``[n, d]_q``-code.

Parameters:
  - `q::Integer`: the number of symbols in the code.
  - `n::Integer`: the word length.
  - `d::Integer`: the distance of the code.

Returns:
  - `Integer`: the sphere covering bound.

---

### Examples

```julia
julia> sphere_covering_bound(5,7,3)
215
```
"""
sphere_covering_bound(q::T, n::T, d::T) where {T <: Integer} = __sphere_bound(ceil, q, n, d - 1)

"""
```julia
sphere_packing_bound(q::Integer, n::Integer, d::Integer) -> Integer
sphere_packing_bound(q::Integer, n::Integer, d::Integer, ::Rounding) -> Real
```

Computes the sphere packing bound of a ``[n, d]_q``-code.  The sphere packing bound is also known as the hamming bound.  You can use `hamming_bound` to compute the same thing.

Parameters:
  - `q::Integer`: the number of symbols in the code.
  - `n::Integer`: the word length.
  - `d::Integer`: the distance of the code.
  - `::Rounding`: use the argument `no_round` in this position to preserve the rounding of the code &mdash; which usually by default rounds down.

Returns:
  - `Integer`: the sphere packing bound.

---

### Examples

```julia
julia> sphere_packing_bound(5,7,3)
2693
```
"""
sphere_packing_bound(q::T, n::T, d::T) where T <: Integer =
	__sphere_bound(a -> floor(T, a), q, n, floor(T, (d - 1) / 2))
sphere_packing_bound(q::T, n::T, d::T, ::Rounding) where T <: Integer =
	__sphere_bound(identity, q, n, floor(T, (d - 1) / 2))
hamming_bound(q::T, n::T, d::T) where T <: Integer =
	sphere_packing_bound(q, n, d)
hamming_bound(q::T, n::T, d::T, ::Rounding) where T <: Integer =
	sphere_packing_bound(q, n, d, no_round)

"""
```julia
sphere_packing_bound(q::Integer, n::Integer, d::Integer) -> Real
```

Computes the Singleton bound of a ``[n, d]_q``-code.

Parameters:
  - `q::Integer`: the number of symbols in the code.
  - `n::Integer`: the word length.
  - `d::Integer`: the distance of the code.

Returns:
  - `Real`: the Singleton bound.  Can round down, as it is an equivalent to the Hamming bound in that it is an upper bound.
"""
# promote()
# _T = typeof(T)

singleton_bound(q::T, n::T, d::T) where T <: Integer =
	floor(T, float(big(q))^(big(n) - big(d) + 1))
singleton_bound(q::T, n::T, d::T, ::Rounding) where T <: Integer =
	float(big(q))^(big(n) - big(d) + 1)


gilbert_varshamov_bound(q::T, n::T, d::T) where T <: Integer =
	__sphere_bound(a -> floor(T, a), q, n, d - 1)
gilbert_varshamov_bound(q::T, n::T, d::T, ::Rounding) where T <: Integer =
	__sphere_bound(identity, q, n, d - 1)

function __plotkin_bound_core(round_func::Function, q::T, n::T, d::T) where T <: Integer
	if ! isequal(q, 2)
		throw(error("The Plotkin bound only works for the binary code."))
	end

	if iseven(d) && 2d > n
		return round_func((d) / (2d + 1 - n))
	elseif isodd(d) && 2d + 1 > n
		return round_func((d + 1) / (2d + 1 - n))
	elseif iseven(d)
		return T(4d)
		# return A_2(2d, d) ≤ 4d
	elseif isodd(d)
		return T(4d + 4)
		# return A_2(2d + 1, d) ≤ 4d + 4
	end
end

plotkin_bound(q::T, n::T, d::T) where T <: Integer =
	__plotkin_bound_core(a -> floor(T, a), q, n, d)
plotkin_bound(q::T, n::T, d::T, ::Rounding) where T <: Integer =
	__plotkin_bound_core(identity, q, n, d, no_round)

elias_bassalygo_bound(q::T, n::T, d::T) where T <: Integer =
elias_bassalygo_bound(q::T, n::T, d::T, ::Rounding) where T <: Integer =

function __johnson_bound_core(round_func::Function, q::T, n::T, d::T) where T <: Integer
	if isinteger((d - 1) / 2) # is odd
		t = T((d - 1) / 2)
		__sphere_bound(round_func, q, n, t) # if d = 2t + 1
	elseif isinteger(d / 2)
		t = T(d / 2)
		__sphere_bound(round_func, q, n, t)
	end
end

@doc raw"""
```julia
construct_ham_matrix(r::Int, q::Int) -> Matrix
```

Construct a Hamming parity-check matrix.

Parameters:
  - `r::Int`: number of rows of a parity check matrix.
  - `q:::Int`: The size of the alphabet of the code.

Returns:
  - `Matrix`: The Hamming matrix, denoted as ``\text{Ham}(r, q)``

---

### Examples

```julia
julia> construct_ham_matrix(3,2)
3×7 Array{Int64,2}:
 0  0  0  1  1  1  1
 0  1  1  0  0  1  1
 1  0  1  0  1  0  1
```
"""
function construct_ham_matrix(r::Int, q::Int)
    ncols = Int(floor((q^r - 1) / (q - 1)))
    M = Matrix{Int}(undef, r, ncols)

    for i in 1:ncols
        M[:, i] = reverse(digits(parse(Int, string(i, base = q)), pad = r), dims = 1)
    end

    return M
end

"""
```julia
isperfect(n::Int, k::Int, d::Int, q::Int) -> Bool
```

Checks if a code is perfect.  That is, checks if the number of words in the code is exactly the "Hamming bound", or the "Sphere Packing Bound".

Parameters:
  - `q:::Int`: The size of the alphabet of the code.
  - `n::Int`: The length of the words in the code (block length).
  - `d::Int`: The distance of the code.
  - `k::Int`: The dimension of the code.

Returns:
  - `Bool`: true or false

---

### Examples

```julia
julia> isperfect(11, 6, 5, 3)
true
```
"""
function isperfect(n::T, k::T, d::T, q::T) where T <: Int
	isprimepower(q) || throw(error("Cannot check if the code is perfect with q not a prime power."))
    M = q^k

    isequal(sphere_packing_bound(q, n, d), M) && return true
	return false
end

"""
```julia
ishammingbound(r::Int, q::Int) -> Bool
```

Checks if the code is a perfect code that is of the form of a generalised Hamming code.

Parameters:
  - `r::Int`: number of rows of a parity check matrix.
  - `q::Int`: The size of the alphabet of the code.

Returns:
  - `Bool`: true or false
"""
function ishammingperfect(r::Int, q::Int)
    n = 2^r - 1
    k = n - r
    M = q^k
    d = size(construct_ham_matrix(r, q), 1) # the number of rows of the hamming matrix (which is, by design, linearly independent)
    d = r
    # r is dim of dueal code; dim of code itself is block length minus r
    # println(n)
    # println((q^r - 1) / (q - 1))

    isequal(n, (q^r - 1) / (q - 1)) && \
		isequal(d, 3) && \
		isequal(M, q^(((q^r - 1) / (q - 1)) - r)) && \
        return true
    return false
end

"""
```julia
ishammingperfect(n::Int, k::Int, d::Int, q::Int) -> Bool
ishammingperfect(q::Int, n::Int, d::Int) -> Bool
```

Checks if the code is a perfect code that is of the form of a generalised Hamming code.

Parameters:
  - `q:::Int`: The size of the alphabet of the code.
  - `n::Int`: The length of the words in the code (block length).
  - `d::Int`: The distance of the code.
  - `k::Int`: The dimension of the code.

Returns:
  - `Bool`: true or false

---

### Examples

```julia
julia> isgolayperfect(11, 6, 5, 3) # this is one of golay's perfect codes
true
```
"""
function ishammingperfect(n::T, k::T, d::T, q::T) where T <: Int
    isprimepower(q) || return false

    M = q^k
    r = log(ℯ, ((n * log(ℯ, 1)) / (log(ℯ, 2))) + 1) / log(ℯ, 2)

    if isequal(n, (q^(r - 1)) / (q - 1)) && isequal(d, 3) && isequal(M, q^(((q^r - 1) / (q - 1)) - r))
        return true
    end

    return false
end
function ishammingperfect(q::Int, n::Int, d::Int)
	isprimepower(q) || return false # we are working in finite fields, so q must be a prime power
	d ≠ 3 && return false

	r = 1
	while ((q^r - 1) / (q - 1)) < n
		r = r + 1
	end

	return ifelse(isequal(((q^r - 1) / (q - 1)), n), true, false)
end

"""
```julia
isgolayperfect(n::Int, k::Int, d::Int, q::Int) -> Bool
```

Golay found two perfect codes.  `isgolayperfect` checks if a code of block length n, distance d, alphabet size q, and dimension k, is a perfect code as described by Golay.

Parameters:
  - `n::Int`: The block length of words in the code (e.g., word "abc" has block length 3).
  - `k::Int`: The dimension of the code.
  - `d::Int`: The distance of the code (i.e., the minimum distance between codewords in the code).
  - `q::Int`: An Int that is a prime power.  The modulus of the finite field.

Returns:
  - `Bool`: true or false.

---

### Examples

```julia
julia> isgolayperfect(11, 6, 5, 3) # this is one of golay's perfect codes
true
```
"""
function isgolayperfect(n::T, k::T, d::T, q::T) where T <: Int
	isprimepower(q) ||  false # we are working in finite fields, so q must be a prime power
    M = q^k
    (isequal(q, 2) && isequal(n, 23) && isequal(d, 7) && isequal(M, 2^12)) && return true
    (isequal(q, 3) && isequal(n, 11) && isequal(d, 5) && isequal(M, 3^6)) && return true
    return false
end
