# Copyright (c) 2021, 2022, 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
# Misc
#############################

function ZeroCode(F::CTFieldTypes, n::Integer)
    n > 0 || throw(ArgumentError("Code length must be positive (received n = $n)"))
    _, vars = polynomial_ring(Nemo.ZZ, Int(order(F)))
    return LinearCode(
        F,
        n,
        0,
        0,
        0,
        0,
        zero_matrix(F, 1, n),
        identity_matrix(F, n),
        zero_matrix(F, 0, n),
        identity_matrix(F, n),
        missing,
        WeightEnumerator(vars[1]^n, :complete),
    )
end

function ZeroCode(q::Integer, n::Integer)
    F = if is_prime(q)
        Oscar.Nemo.Native.GF(q)
    else
        factors = Nemo.factor(q)
        length(factors) == 1 ||
            throw(DomainError("There is no finite field of order $q"))
        p, t = first(factors)
        GF(p, t, :α)
    end
    return ZeroCode(F, n)
end

ZeroCode(n::Integer) = ZeroCode(GF(2), n)
IdentityCode(F::CTFieldTypes, n::Integer) = dual(ZeroCode(F, n))
IdentityCode(q::Integer, n::Integer) = dual(ZeroCode(q, n))
IdentityCode(n::Integer) = dual(ZeroCode(n))

"""
    RepetitionCode(q::Int, n::Int)

Return the `[n, 1, n]` repetition code over `GF(q)`.
"""
function RepetitionCode(q::Int, n::Int)
    F = if is_prime(q)
        Oscar.Nemo.Native.GF(q)
    else
        factors = Nemo.factor(q)
        length(factors) == 1 ||
            throw(DomainError("There is no finite field of order $q"))
        (p, t), = factors
        GF(p, t, :α)
    end
    G = matrix(F, ones(Int, 1, n))
    H = hcat(matrix(F, ones(Int, n - 1, 1)), identity_matrix(F, n - 1))
    G_stand, H_stand, P, _ = _standard_form(G)
    _, vars = polynomial_ring(Nemo.ZZ, Int(order(F)))
    W = WeightEnumerator(sum(v^n for v in vars), :complete)
    return LinearCode(F, n, 1, n, n, n, G, H, G_stand, H_stand, P, W)
end

"""
    SingleParityCheckCode(q::Int, n::Int)

Return the `[n, n-1, 2]` single parity check code over `GF(q)`.
"""
function SingleParityCheckCode(q::Int, n::Int)
    iseven(q) && (return dual(RepetitionCode(q, n));)
    F = if is_prime(q)
        Oscar.Nemo.Native.GF(q)
    else
        factors = Nemo.factor(q)
        length(factors) == 1 ||
            throw(DomainError("There is no finite field of order $q"))
        (p, t), = factors
        GF(p, t, :α)
    end
    G = hcat(matrix(F, .-ones(Int, n - 1, 1)), identity_matrix(F, n - 1))
    return LinearCode(G)
end
SPCCode(q::Int, n::Int) = SingleParityCheckCode(q, n)

"""
    Hexacode()

Return the `[6, 3, 4]` hexacode over `GF(4)`.
"""
function Hexacode()
    # TODO: is Hexacode a Hamming code?
    F = GF(2, 2, :ω)
    ω = gen(F)
    G = matrix(F, [1 0 0 1 ω ω; 0 1 0 ω 1 ω; 0 0 1 ω ω 1])
    H = matrix(F, [1 ω ω 1 0 0; ω 1 ω 0 1 0; ω ω 1 0 0 1])
    G_stand, H_stand, P, rnk = _standard_form(G)
    return LinearCode(
        F,
        6,
        3,
        4,
        4,
        4,
        G,
        H,
        G_stand,
        H_stand,
        P,
        _weight_enumerator_BF(G_stand),
    )
end

#############################
# Hamming
#############################
# unclear if this should be promoted to its own type so r can be extracted

"""
    HammingCode(q::Int, r::Int)

Return the `[(q^r - 1)/(q - 1), (q^r - 1)/(q - 1) - r, 3]` Hamming code over `GF(q)`.

# Notes
* This is currently only implemented for binary codes.
"""
function HammingCode(q::Int, r::Int)
    2 ≤ r || throw(DomainError("Hamming codes require r ≥ 2; received r = $r."))
    r < 64 || throw(
        DomainError(
            "This Hamming code requires the implmentation of BigInts. Change if necessary.",
        ),
    )
    q == 2 || throw(DomainError("Nonbinary Hamming codes have not yet been implemented."))
    factors = Nemo.factor(q)
    length(factors) == 1 || throw(ArgumentError("There is no finite field of order $q."))

    if q == 2
        F = Oscar.Nemo.Native.GF(2)
        # there are faster ways to do this using trees, but the complexity and
        # overhead is not worth it for the sizes required here
        H = matrix(
            F,
            reduce(hcat, [reverse(digits(i, base = 2, pad = r)) for i = 1:(2^r-1)]),
        )
        C = LinearCode(H, true, false)
        set_minimum_distance!(C, 3)
        R, vars = polynomial_ring(Nemo.ZZ, 2)
        C.weight_enum = WeightEnumerator(
            divexact(
                (vars[2] + vars[1])^C.n +
                C.n *
                (vars[2] + vars[1])^div(C.n - 1, 2) *
                (vars[1] - vars[2])^div(C.n + 1, 2),
                C.n + 1,
            ),
            :complete,
        )
        return C
    end

    # will have to check this example much further in the book because maybe I need
    # to staircase this 1 up from the bottom to the top

    # for non-binary, let's try a ProductIterator over all tuples of size two smaller
    # store them and then construct a matrix where we add 0 then 1 to the top of
    # each computed tuple
    # then compute all tuples the size larger and add a 1 to the top
end

"""
    TetraCode()

Return the `[4, 2, 3]` tetra code over `GF(3)`.

# Notes
* This is equiavlent to the `Hamming(3, 2, 3)` code, but the construction here is
  based on the commonly presented generator and parity-check matrices.
"""
function TetraCode()
    F = Oscar.Nemo.Native.GF(3)
    G = matrix(F, [1 0 1 1; 0 1 1 -1])
    H = matrix(F, [-1 -1 1 0; -1 1 0 1])
    G_stand, H_stand, P, rnk = _standard_form(G)
    R, vars = polynomial_ring(Nemo.ZZ, 3)
    CWE = WeightEnumerator(
        vars[1]^4 +
        vars[1]*vars[2]^3 +
        3*vars[1]*vars[2]^2*vars[3] +
        3*vars[1]*vars[2]*vars[3]^2 +
        vars[1]*vars[3]^3,
        :complete,
    )
    return LinearCode(F, 4, 2, 3, 3, 3, G, H, G_stand, H_stand, P, CWE)
end

#############################
# Simplex
#############################

# if kept as separate type, then can detect this and return simplexcode constrcutor
# instead of regular dual
"""
    SimplexCode(q::Int, r::Int)

Return the `[(q^r - 1)/(q - 1), r]` simplex code over `GF(q)`.

# Notes
* Generator matrices for the binary codes are constructed using the standard
  recursive definition. The higher fields return `dual(HammingCode(q, r))`.
* This is currently only implemented for binary codes.
"""
function SimplexCode(q::Int, r::Int)
    2 ≤ r || throw(DomainError("Simplex codes require 2 ≤ r; received r = $r."))
    r < 64 || throw(
        DomainError(
            "The weight enumerator for the simplex codes for r > 64 require BigInts. Implement if necessary.",
        ),
    )
    q == 2 || throw(DomainError("Nonbinary simplex codes have not yet been implemented."))

    # actually really need to check here that q^r is not over sizeof(Int)
    factors = Nemo.factor(q)
    length(factors) == 1 || throw(ArgumentError("There is no finite field of order $q."))

    # the known weight distribution is Hamming and not complete
    q > 2 && return dual(HammingCode(q, r))

    # binary simplex codes
    F = Oscar.Nemo.Native.GF(2)
    G2 = matrix(F, [0 1 1; 1 0 1]);
    if r == 2
        C = LinearCode(G2, false, false)
    else
        Grm1 = G2
        for i = 3:r
            zs = matrix(F, nrows(Grm1), 1, zeros(Int, nrows(Grm1), 1))
            bot = hcat(Grm1, zs, Grm1)
            zs = matrix(F, 1, ncols(Grm1), zeros(Int, 1, ncols(Grm1)))
            os = matrix(F, 1, ncols(Grm1) + 1, ones(Int, 1, ncols(Grm1) + 1))
            top = hcat(zs, os)
            Grm1 = vcat(top, bot)
        end
        C = LinearCode(Grm1, false, false)
    end
    # all nonzero codewords have weights q^{r - 1}
    # should have q^r - 1 nonzero codewords
    R, vars = polynomial_ring(Nemo.ZZ, 2)
    C.weight_enum = WeightEnumerator(
        vars[1]^(2^r - 1) + (2^r - 1) * vars[1]^(2^r - 2^(r - 1) - 1) * vars[2]^(2^(r - 1)),
        :complete,
    )
    set_minimum_distance!(C, 2^(r - 1))
    return C
end

#############################
# Golay
#############################

"""
    ExtendedGolayCode(p::Int)

Return the [24, 12, 8] extended binary Golay code if `p == 2` or the [12, 6, 6]
extended ternary Golay code if `p == 3`.
"""
function ExtendedGolayCode(p::Int)
    if p == 2
        F = Oscar.Nemo.Native.GF(2)
        A = matrix(
            F,
            [
                0 1 1 1 1 1 1 1 1 1 1 1;
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
                1 0 1 1 0 1 1 1 0 0 0 1
            ],
        )
        G = hcat(identity_matrix(F, 12), A)
        H = hcat(-transpose(A), identity_matrix(F, 12))
        R, vars = polynomial_ring(Nemo.ZZ, 2)
        wt_enum = WeightEnumerator(
            vars[1]^24 +
            759*vars[2]^8*vars[1]^16 +
            2576 * vars[2]^12 * vars[1]^12 +
            759*vars[1]^8*vars[2]^16 +
            vars[2]^24,
            :complete,
        )
        return LinearCode(F, 24, 12, 8, 8, 8, G, H, G, H, missing, wt_enum)
    elseif p == 3
        F = Oscar.Nemo.Native.GF(3)
        A = matrix(
            F,
            [
                0 1 1 1 1 1;
                1 0 1 -1 -1 1;
                1 1 0 1 -1 -1;
                1 -1 1 0 1 -1;
                1 -1 -1 1 0 1;
                1 1 -1 -1 1 0
            ],
        )
        G = hcat(identity_matrix(F, 6), A)
        H = hcat(-transpose(A), identity_matrix(F, 6))
        R, vars = polynomial_ring(Nemo.ZZ, 2)
        # this looks like Hamming and not complete
        # wt_enum = WeightEnumerator(vars[1]^12 + 264*vars[2]^6*vars[1]^6 +
        #     440*vars[2]^9*vars[1]^3 + 24*vars[2]^12, :complete)
        return LinearCode(F, 12, 6, 6, 6, 6, G, H, G, H, missing, _weight_enumerator_BF(G))
    else
        throw(ArgumentError("Golay code not implemented for q = $q."))
    end
end

# The `[24, 12, 8]` extended binary Golay code may be obtained by extending this
# code. Adding an overall parity check to the ternary Golay code will give either
# a `[12, 6, 6]` if punctured and extended in the first coordinate or a
# `[12, 6, 5]` code otherwise.
"""
    GolayCode(p::Int)

Return the `[23, 12, 7]`` binary Golay code if `p == 2` or the `[11, 6, 5]`
ternary Golay code if `p == 3`.
"""
function GolayCode(p::Int)
    C = puncture(ExtendedGolayCode(p), [1])
    p == 2 ? (C.d = 7;) : (C.d = 5;)
    return C
end

#############################
# Hadamard
#############################

# """
#     HadamardCode(m)
#     WalshHadamardCode(m)
#     WalshCode(m)
#
# Return the `[2^m, m, 2^{m - 1}]` binary Hadamard code.
#
# Hadamard codes are generally constructed using Hadamard matrices `H_n`. When `n`
# is a power of two, the codes are linear. Constructing `H_{2^m}` then mapping `+/- 1`
# to `{0, 1}` gives a generator matrix which, up to permutation, is equivalent to
# using all `2^m` binary strings as column vectors. The construction here uses this
# later definition. These are of course equivalent to the first-order Reed-Muller
# codes `RM(1, m)`, but we do not use that form of the generator matrix here.
#
# Note that some engineering fields define the Hadamard code to be the `[2^m, m +
# 1, 2^{m - 1}]` augmented Hadamard code.
# """
# function HadamardCode(m)
#     m < 64 || error("This Hadamard code requires the implmentation of BigInts. Change if necessary.")
#
#     F = Oscar.Nemo.Native.GF(2)
#     G ...
#     C = LinearCode(G)
#     R, vars = polynomial_ring(Nemo.ZZ, 2)
#     # each non-zero codeword has a Hamming weight of exactly 2^{k-1}
#     C.weight_enum = WeightEnumerator(vars[1]^(2^m) + (2^m - 1) * vars[1]^2 *
#         vars[2]^(2^m - 1), :complete)
#     return C
# end
# WalshHadamardCode(m) = HadamardCode(m)
# WalshCode(m) = HadamardCode(m)

#############################
# Best Known Linear Codes
#############################

function best_known_linear_code(n::Int, k::Int)
    C_GAP = GAP.Globals.BestKnownLinearCode(n, k, GAP.Globals.GF(2))
    G = GAP.Globals.GeneratorMat(C_GAP)
    dims = GAP.Globals.DimensionsMat(G)
    g = matrix(GF(2), [GAP.Globals.Int(G[i, j]) for i = 1:dims[1], j = 1:dims[2]])
    H = GAP.Globals.CheckMat(C_GAP)
    dims = GAP.Globals.DimensionsMat(H)
    h = matrix(GF(2), [GAP.Globals.Int(H[i, j]) for i = 1:dims[1], j = 1:dims[2]])
    return LinearCode(g, h)
end
