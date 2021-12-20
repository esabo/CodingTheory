# CodingTheory
 A basic coding theory library for Julia.

TODO:
iterator for linear code
function something(C::T) where T <: AbstractLinearCode
    # object has matrix C.G
    # need to return all linear combinations of sums of vectors
    # can easily be extremely large to return these in a Vector so I'm thinking just an interator to go through them
    # then one can store them outside the function if necessary or merely generate them to process one-by-one
    # other julia library actually does this given a Vector of Vectors (but has no concept of a matrix or a code itself)
    # maybe we just take the below implementation?
	# don't need to mod anything since C.G is over a finite field and that will automatically be taken care of
end

function words(C::T) where T <: AbstractLinearCode
	# collect the above
end
codewords(C:T) where T <: AbstractLinearCode = words(C)

"""
```julia
get_codewords(G::AbstractArray, m::Int) -> Codewords{M}
```
Get codewords of a code from the _generating matrix_ under a finite field of modulo `m`.  Precisely, computes all linear combinations of the rows of the generating matrix.

Parameters:
  - `G::AbstractArray`: A matrix of Ints which generates the code.
  - `m::Int`: The bounds of the finite field (i.e., the molulus you wish to work in).

Returns:
  - `Codewords{M}`: An array of codewords, each of length `M`.  Each codewords is a tuple, and each character in said word is a symbol.
"""

function get_codewords(G::AbstractArray, m::Int)
	codewords = Vector()
	rows = Vector(undef, size(G, 2))

	for i in 1:size(G, 1)
		rows[i] = [G[i, j] for j in 1:size(G, 2)]
	end

	for c in Base.Iterators.product([0:m-1 for _ in 1:size(G, 1)]...)
		word = Ref(c[1]) .* rows[1]

		for i in 2:size(G, 1)
			word = mod.(word .+ (Ref(c[i]) .* rows[i]), m)
		end

		push!(codewords, word)
	end

	return codewords
end





I will consider this a brute force solution to an NP-Hard problem and do a more sophisticated one later
I will also have to call cardinality(C) and put a max limit on this number or else it errors
function weightdistribution(C::T) where T <: AbstractLinearCode
	# make Vector of zeros from length 1:length(C)
	# iterate through words above and +1 in spot i if vector has Hamming weight = i
end



in Hamming.jl
I do need to replace all my Int64's with general integers...
function Hammingparitycheckmatrix(m::Int64, q::Int64)
	# return a matrix whose columns are all nonzero m-tuples from GF(q) with first nonzero entry == 1
	# in binary this is merely all binary strings as columns
	# below is other library's code for this where his r = my m (m is the standard letter in my textbooks)
	# however, our implementation needs to start with
	F, _ = FiniteField(q, 1, "α")
	# and then end with
	return matrix(F, (other matrix he creates below))
	# unless it's easier to make a zero one and then chance column slices then
	M = MatrixSpace(F, m, n)
	# where n is as below

end

function construct_ham_matrix(r::Int, q::Int)
    ncols = Int(floor((q^r - 1) / (q - 1)))
    M = Matrix{Int}(undef, r, ncols)

    for i in 1:ncols
        M[:, i] = reverse(digits(parse(Int, string(i, base = q)), pad = r), dims = 1)
    end

    return M
end
