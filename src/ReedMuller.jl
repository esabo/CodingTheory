# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

mutable struct ReedMullerCode <: AbstractReedMullerCode
    F::Union{FqNmodFiniteField, Nemo.GaloisField, AbstractAlgebra.GFField{Int64}}
    n::Integer
    k::Integer
    d::Union{Integer, Missing}
    r::Integer
    m::Integer
    G::Union{gfp_mat, fq_nmod_mat}
    Gorig::Union{gfp_mat, fq_nmod_mat, Missing}
    H::Union{gfp_mat, fq_nmod_mat}
    Horig::Union{gfp_mat, fq_nmod_mat, Missing}
    Gstand::Union{gfp_mat, fq_nmod_mat}
    Hstand::Union{gfp_mat, fq_nmod_mat}
    weightenum::Union{WeightEnumerator, Missing}
end

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
    if get(io, :compact, false)
        println(io, "[$(C.n), $(C.k), $(minimumdistance(C))]_$(order(C.F)) Reed-Muller code RM($(C.r), $(C.m)).")
    else
        println(io, "[$(C.n), $(C.k), $(minimumdistance(C))]_$(order(C.F)) Reed-Muller code RM($(C.r), $(C.m)).")
        println(io, "Generator matrix: $(C.k) × $(C.n)")
        for i in 1:C.k
            print(io, "\t")
            for j in 1:C.n
                if j != C.n
                    print(io, "$(C.G[i, j]) ")
                elseif j == C.n && i != C.k
                    println(io, "$(C.G[i, j])")
                else
                    print(io, "$(C.G[i, j])")
                end
            end
        end
        if !ismissing(C.weightenum)
            println(io, "\nComplete weight enumerator:")
            print(io, "\t", C.weightenum.polynomial)
        end
    end
end

"""
    ReedMullergeneratormatrix(q::Integer, r::Integer, m::Integer)

Return the recursive form of the generator matrix for the `RM(r, m)` Reed-Muller
code over `GF(q)`.
"""
function ReedMullergeneratormatrix(q::Integer, r::Integer, m::Integer)
    (0 ≤ r ≤ m) || error("Reed-Muller codes require 0 ≤ r ≤ m, received r = $r and m = $m.")

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
        error("Nonbinary Reed-Muller codes have not yet been implemented.")
    end
end

"""
    ReedMullerCode(q::Integer, r::Integer, m::Integer)

Return the `RM(r, m)` Reed-Muller code over `GF(q)`.
"""
function ReedMullerCode(q::Integer, r::Integer, m::Integer)
    (0 ≤ r < m) || error("Reed-Muller codes require 0 ≤ r < m, received r = $r and m = $m.")
    m < 64 || error("This Reed-Muller code requires the implmentation of BigInts. Change if necessary.")
    q == 2 || error("Nonbinary Reed-Muller codes have not yet been implemented.")
    factors = factor(q)
    length(factors) == 1 || error("There is no finite field of order $q.")
    (p, t), = factors

    G = ReedMullergeneratormatrix(q, r, m)
    H = ReedMullergeneratormatrix(q, m - r - 1, m)
    Gstand, Hstand = _standardform(G)

    # verify
    ncols(G) == 2^m || error("Generator matrix computed in ReedMuller has the wrong number of columns; received: $(ncols(G)), expected: $(BigInt(2)^m).")
    k = sum([binomial(m, i) for i in 0:r])
    nrows(G) == k || error("Generator matrix computed in ReedMuller has the wrong number of rows; received: $(nrows(G)), expected: $k.")
    if size(H) == (2^m - k, k)
        H = transpose(H)
    end
    for r in 1:nrows(Gstand)
        iszero(Gstand[r, :] * transpose(H)) || error("Column swap appeared in _standardform.")
    end

    if q == 2 && r == 1
        R, vars = PolynomialRing(Nemo.ZZ, 2)
        poly = vars[1]^(2^m) + (2^(m + 1) - 2) * vars[1]^(2^m - 2^(m - 1))*vars[2]^(2^(m - 1)) + vars[2]^(2^m)
        return ReedMullerCode(base_ring(G), ncols(G), nrows(G), 2^(m - r), r, m, G, missing,
            H, missing, Gstand, Hstand, WeightEnumerator(poly, "complete"))
    end

    return ReedMullerCode(base_ring(G), ncols(G), nrows(G), 2^(m - r), r, m, G, missing,
        H, missing, Gstand, Hstand, missing)
end

"""
    dual(C::ReedMullerCode)

Return the dual of the Reed-Muller code `C`.
"""
function dual(C::ReedMullerCode)
    # really only put this here to remind me later that I hardcoded binary
    if Int(characteristic(C.F)) == 2
        return ReedMullerCode(C.F, C.n, C.n - C.k, 2^(C.r + 1),
            C.m - C.r - 1, C.m, C.H, missing, C.G, missing, C.Hstand, C.Gstand,
            missing)
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
    C.F == D.F || error("Codes must be over the same field in the Schur product.")
    C.n == D.n || error("Codes must have the same length in the Schur product.")

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
