#############################
        # Repetition
#############################

"""
    RepetitionCode(q::Integer, n::Integer)

Return the `[n, 1, n]` repetition code over `GF(q)`.
"""
function RepetitionCode(q::Integer, n::Integer)
    F, _ = FiniteField(q, 1, "α")
    G = matrix(F, 1, n, [1 for i in 1:n])
    M2 = MatrixSpace(F, n - 1, 1)
    M3 = MatrixSpace(F, n - 1, n - 1)
    H = hcat(M2([1 for i in 1:(n - 1)]), M3(1))
    return LinearCode(F, n, 1, n, G, G, H, H, G, H, missing)
end

#############################
          # Tetra
#############################

# function tetracode()
#     return Hammingcode(3, 2)
# end

#############################
         # Hamming
#############################

# function Hammingcode(p::Integer, r::Integer)
#
# end

# function construct_ham_matrix(r::Int, q::Int)
#     ncols = Int(floor((q^r - 1) / (q - 1)))
#     M = Matrix{Int}(undef, r, ncols)
#
#     for i in 1:ncols
#         M[:, i] = reverse(digits(parse(Int, string(i, base = q)), pad = r), dims = 1)
#     end
#
#     return M
# end





#############################
         # Simplex
#############################

# if kept as separate type, then can detect this and return simplexcode constrcutor
# instead of regular dual
"""
    SimplexCode(q::Int, r::Int)

Return the `[(q^r - 1)/(q - 1), r]` simplex code over `GF(q)`.

Generator matrices for the binary codes are constructed using the standard
recursive definition. The higher fields return `dual(HammingCode(q, r))`.

# Notes
* This is currently only implemented for binary codes.
"""
function SimplexCode(q::Int, r::Int)
    2 ≤ r || error("Simplex codes require 2 ≤ r; received r = $r.")
    r < 64 || error("The weight enumerator for the simplex codes for r > 64 require BigInts. Implement if necessary.")
    q == 2 || error("Nonbinary simplex codes have not yet been implemented.")

    # actually really need to check here that q^r is not over sizeof(Int)
    if !isprime(q)
        factors = factor(q)
        if length(factors) != 1
            error("There is no finite field of order $(prod(factors)).")
        end
    end

    # the known weight distribution is Hamming and not complete
    q > 2 && return dual(HammingCode(q, r))

    # binary simplex codes
    F, _ = FiniteField(2, 1, "α");
    G2 = matrix(F, [0 1 1; 1 0 1]);
    if r == 2
        C = LinearCode(G2)
        # all nonzero codewords have weights q^{r - 1}
        # should have q^r - 1 nonzero codewords
        # coeff, 0's, 1's
        C.weightenum = WeightEnumerator([[1, 3, 0], [3, 1, 2]], "complete")
        setminimumdistance!(C, 3)
        return C
    else
        Grm1 = G2
        for i in 3:r
            zs = matrix(F, size(Grm1, 1), 1, zeros(Int, size(Grm1, 1), 1))
            bot = hcat(Grm1, zs, Grm1)
            zs = matrix(F, 1, size(Grm1, 2), zeros(Int, 1, size(Grm1, 2)))
            os = matrix(F, 1, size(Grm1, 2) + 1, ones(Int, 1, size(Grm1, 2) + 1))
            top = hcat(zs, os)
            Grm1 = vcat(top, bom)
        end
        C = LinearCode(Grm1)
        # all nonzero codewords have weights q^{r - 1}
        # should have q^r - 1 nonzero codewords
        # coeff, 0's, 1's
        C.weightenum = WeightEnumerator([[1, 2^r, 0], [2^r - 1, 2^r - 2^(r - 1) - 1,
            2^(r - 1)]], "complete")
        setminimumdistance!(C, 2^(r - 1))
        return C
    end
end

#############################
          # Golay
#############################

"""
    ExtendedGolayCode(p::Integer)

Return the [24, 12, 8] extended binary Golay code if `p == 2` or the [12, 6, 6]
extended ternary Golay code if `p == 3`.

The [23, 12, 7] binary Golay code may be obtained by puncturing this code in
any coordinate.

# Notes
* This is constructed using an explicit matrix. The Golay codes are constructed
 by calling `puncture(ExtendedGolayCode(p), 1)`. All single punctures are
 equivalent.
"""
function ExtendedGolayCode(p::Integer)
    if p == 2
        F, _ = FiniteField(2, 1, "α")
        M = MatrixSpace(F, 12 , 12)
        A = M([0 1 1 1 1 1 1 1 1 1 1 1;
             1 1 1 0 1 1 1 0 0 0 1 0;
             1 1 0 1 1 1 0 0 0 1 0 1;
             1 0 1 1 1 0 0 0 1 0 1 1;
             1 1 1 1 0 0 0 1 0 1 1 0;
             1 1 1 0 0 0 1 0 1 1 0 1;
             1 1 0 0 0 1 0 1 1 0 1 1;
             1 0 0 0 1 0 1 1 0 1 1 1;
             1 0 0 1 0 1 1 0 1 1 1 0;
             1 0 1 0 1 1 0 1 1 1 0 0;
             1 1 0 1 1 0 1 1 1 0 0 0;
             1 0 1 1 0 1 1 1 0 0 0 1])
        G = hcat(M(1), A)
        H = hcat(-transpose(A), M(1))
        return LinearCode(F, 24, 12, 8, G, G, H, H, G, H, missing)
    elseif p == 3
        F, _ = FiniteField(3, 1, "α")
        M = MatrixSpace(F, 6 , 6)
        A = M([0 1 1 1 1 1;
               1 0 1 -1 -1 1;
               1 1 0 1 -1 -1;
               1 -1 1 0 1 -1;
               1 -1 -1 1 0 1;
               1 1 -1 -1 1 0])
        G = hcat(M(1), A)
        H = hcat(-transpose(A), M(1))
        return LinearCode(F, 12, 6, 6, G, G, H, H, G, H, missing)
    else
        error("Golay code not implemented for q = $q.")
    end
end

"""
    GolayCode(p::Integer)

Return the `[23, 12, 7]`` binary Golay code if `p == 2` or the `[11, 6, 5]`
ternary Golay code if `p == 3`.

The `[24, 12, 8]` extended binary Golay code may be obtained by extending this
code. Adding an overall parity check to the ternary Golay code will give either
a `[12, 6, 6]` if punctured and extended in the first coordinate or a
`[12, 6, 5]` code otherwise.

# Notes
* These codes are constructed by calling `puncture(ExtendedGolayGode(p), 1)`.
 All single punctures are equivalent.
"""
function GolayCode(p::Integer)
    return puncture(ExtendedGolayCode(p), 1)
end
