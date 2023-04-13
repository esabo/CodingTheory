# Copyright (c) 2021, 2023 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
    order(C::ReedMullerCode)
    RMr(C::ReedMullerCode)

Return the order, `r`, of the `RM(r, m)` Reed-Muller code.
"""
order(C::ReedMullerCode) = C.r
RMr(C::ReedMullerCode) = order(C)

"""
    numberofvariables(C::ReedMullerCode)
    RMm(C::ReedMullerCode)

Return the number of variables, `m`, of the `RM(r, m)` Reed-Muller code.
"""
numberofvariables(C::ReedMullerCode) = C.m
RMm(C::ReedMullerCode) = C.m

function show(io::IO, C::AbstractReedMullerCode)
    if ismissing(C.d)
        println(io, "[$(C.n), $(C.k)]_$(order(C.F)) Reed-Muller code RM($(C.r), $(C.m)).")
    else
        println(io, "[$(C.n), $(C.k), $(C.d)]_$(order(C.F)) Reed-Muller code RM($(C.r), $(C.m)).")
    end
    if get(io, :compact, false) && C.n <= 30
        G = generatormatrix(C)
        nr, nc = size(G)
        println(io, "Generator matrix: $nr × $nc")
        for i in 1:nr
            print(io, "\t")
            for j in 1:nc
                if j != nc
                    print(io, "$(G[i, j]) ")
                elseif j == nc && i != nr
                    println(io, "$(G[i, j])")
                else
                    print(io, "$(G[i, j])")
                end
            end
        end
        # if !ismissing(C.weightenum)
        #     println(io, "\nComplete weight enumerator:")
        #     print(io, "\t", C.weightenum.polynomial)
        # end
    end
end

"""
    ReedMullergeneratormatrix(q::Integer, r::Integer, m::Integer)

Return the recursive form of the generator matrix for the `RM(r, m)` Reed-Muller
code over `GF(q)`.
"""
function ReedMullergeneratormatrix(q::Integer, r::Integer, m::Integer)
    (0 ≤ r ≤ m) || throw(DomainError("Reed-Muller codes require 0 ≤ r ≤ m, received r = $r and m = $m."))

    if q == 2
        F, _ = FiniteField(2, 1, "α")
        if r == m
            M = MatrixSpace(F, 2^m, 2^m)
            return M(1)
        elseif r == 0
            M = MatrixSpace(F, 1, 2^m)
            return M([1 for i in 1:2^m])
        else
            Grm1 = ReedMullergeneratormatrix(q, r, m - 1)
            Gr1m1 = ReedMullergeneratormatrix(q, r - 1, m - 1)
            M = MatrixSpace(F, nrows(Gr1m1), ncols(Gr1m1))
            return vcat(hcat(Grm1, Grm1), hcat(M(0), Gr1m1))
        end
    else
        throw(ArgumentError("Nonbinary Reed-Muller codes have not yet been implemented."))
    end
end

"""
    ReedMullerCode(q::Integer, r::Integer, m::Integer)

Return the `RM(r, m)` Reed-Muller code over `GF(q)`.
"""
function ReedMullerCode(q::Integer, r::Integer, m::Integer)
    (0 ≤ r < m) || throw(DomainError("Reed-Muller codes require 0 ≤ r < m, received r = $r and m = $m."))
    m < 64 || throw(DomainError("This Reed-Muller code requires the implmentation of BigInts. Change if necessary."))
    q == 2 || throw(DomainError("Nonbinary Reed-Muller codes have not yet been implemented."))
    factors = factor(q)
    length(factors) == 1 || throw(ArgumentError("There is no finite field of order $q."))
    # (p, t), = factors

    G = ReedMullergeneratormatrix(q, r, m)
    H = ReedMullergeneratormatrix(q, m - r - 1, m)
    Gstand, Hstand, P, rnk = _standardform(G)

    # verify
    ncols(G) == 2^m || error("Generator matrix computed in ReedMuller has the wrong number of columns; received: $(ncols(G)), expected: $(BigInt(2)^m).")
    k = sum([binomial(m, i) for i in 0:r])
    nrows(G) == k || error("Generator matrix computed in ReedMuller has the wrong number of rows; received: $(nrows(G)), expected: $k.")
    size(H) == (2^m - k, k) && (H = transpose(H);)
    d = 2^(m - r)

    if q == 2 && r == 1
        R, vars = PolynomialRing(Nemo.ZZ, 2)
        poly = vars[1]^(2^m) + (2^(m + 1) - 2) * vars[1]^(2^m - 2^(m - 1))*vars[2]^(2^(m - 1)) + vars[2]^(2^m)
        return ReedMullerCode(base_ring(G), ncols(G), nrows(G), d, d, d, r, m, G, missing,
            H, missing, Gstand, Hstand, P, WeightEnumerator(poly, "complete"))
    end

    return ReedMullerCode(base_ring(G), ncols(G), nrows(G), d, d, d, r, m, G, missing,
        H, missing, Gstand, Hstand, P, missing)
end

"""
    dual(C::ReedMullerCode)

Return the dual of the Reed-Muller code `C`.
"""
function dual(C::ReedMullerCode)
    # really only put this here to remind me later that I hardcoded binary
    if Int(characteristic(C.F)) == 2
        d = 2^(C.r + 1)
        return ReedMullerCode(C.F, C.n, C.n - C.k, d, d, d,
            C.m - C.r - 1, C.m, C.H, missing, C.G, missing, C.Hstand, C.Gstand,
            C.P, missing)
    end

end

# "On products and powers of linear codes under componentwise multiplication"
# Hugues Randriambololona
"""
    entrywiseproductcode(C::ReedMullerCode, D::ReedMullerCode)
    *(C::ReedMullerCode, D::ReedMullerCode)
    Schurproductcode(C::ReedMullerCode, D::ReedMullerCode)
    Hadamardproductcode(C::ReedMullerCode, D::ReedMullerCode)
    componentwiseproductcode(C::ReedMullerCode, D::ReedMullerCode)

Return the entrywise product of `C` and `D`.

Note that this is known to often be the full ambient space.
"""
function entrywiseproductcode(C::ReedMullerCode, D::ReedMullerCode)
    C.F == D.F || throw(ArgumentError("Codes must be over the same field in the Schur product."))
    C.n == D.n || throw(ArgumentError("Codes must have the same length in the Schur product."))

    r = C.r + D.r
    if r <= C.n
        return ReedMullerCode(Int(characteristic(C.F)), r, C.m)
    else
        return ReedMullerCode(Int(characteristic(C.F)), C.m, C.m)
    end
end
*(C::ReedMullerCode, D::ReedMullerCode) = entrywiseproductcode(C, D)
Schurproductcode(C::ReedMullerCode, D::ReedMullerCode) = entrywiseproductcode(C, D)
Hadamardproductcode(C::ReedMullerCode, D::ReedMullerCode) = entrywiseproductcode(C, D)
componentwiseproductcode(C::ReedMullerCode, D::ReedMullerCode) = entrywiseproductcode(C, D)
