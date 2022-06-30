include("linearcode.jl")

"""
    repetitioncode(q::Integer, n::Integer)

Return the `[n, 1, n]` repetition code over `GF(q)`.
"""
function repetitioncode(q::Integer, n::Integer)
    F, _ = FiniteField(q, 1, "α")
    G = matrix(F, 1, n, [1 for i in 1:n])
    M2 = MatrixSpace(F, n - 1, 1)
    M3 = MatrixSpace(F, n - 1, n - 1)
    H = hcat(M2([1 for i in 1:(n - 1)]), M3(1))
    return LinearCode(F, n, 1, n, G, G, H, H, G, H, missing)
end

# function Hammingcode(p::Integer, r::Integer)
#
# end

# function tetracode()
#     return Hammingcode(3, 2)
# end
