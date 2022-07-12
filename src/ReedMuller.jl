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
        println(io, "[$(length(C)), $(dimension(C)), $(minimumdistance(C))]_$(order(field(C))) Reed-Muller code RM($(r(C)), $(m(C))).")
    else
        println(io, "[$(length(C)), $(dimension(C)), $(minimumdistance(C))]_$(order(field(C))) Reed-Muller code RM($(r(C)), $(m(C))).")
        println(io, "Generator matrix: $(dimension(C)) × $(length(C))")
        for i in 1:dimension(C)
            print(io, "\t")
            for j in 1:length(C)
                if j != length(C)
                    print(io, "$(C.G[i, j]) ")
                elseif j == length(C) && i != dimension(C)
                    println(io, "$(C.G[i, j])")
                else
                    print(io, "$(C.G[i, j])")
                end
            end
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
            M = MatrixSpace(F, size(Gr1m1, 1), size(Gr1m1, 2))
            return vcat(hcat(Grm1, Grm1), hcat(M(0), Gr1m1))
        end
    else
        error("Nonbinary Reed-Muller codes have not yet been implemented.")
    end
end

"""
    ReedMullerCode(q::Integer, r::Integer, m::Integer, verify::Bool=true)

Return the `RM(r, m)` Reed-Muller code over `GF(q)`.

If the optional parameter `verify` is set to `true`, basic checks are done to
ensure correctness.
"""
function ReedMullerCode(q::Integer, r::Integer, m::Integer, verify::Bool=true)
    (0 ≤ r < m) || error("Reed-Muller codes require 0 ≤ r < m, received r = $r and m = $m.")
    m < 64 || error("This Reed-Muller code requires the implmentation of BigInts. Change if necessary.")
    q == 2 || error("Nonbinary Reed-Muller codes have not yet been implemented.")

    if !isprime(q)
        factors = factor(q)
        if length(factors) != 1
            error("There is no finite field of order $(prod(factors)).")
        end
    end

    G = ReedMullergeneratormatrix(q, r, m)
    H = ReedMullergeneratormatrix(q, m - r - 1, m)
    Gstand, Hstand = _standardform(G)

    if verify
        size(G, 2) == 2^m || error("Generator matrix computed in ReedMuller has the wrong number of columns; received: $(size(G, 2)), expected: $(BigInt(2)^m).")
        k = sum([binomial(m, i) for i in 0:r])
        size(G, 1) == k || error("Generator matrix computed in ReedMuller has the wrong number of rows; received: $(size(G, 1)), expected: $k.")
        if size(H) == (2^m - k, k)
            H = transpose(H)
        end
        for r in 1:size(Gstand, 1)
            iszero(Gstand[r, :] * transpose(H)) || error("Column swap appeared in _standardform.")
        end
    end

    return ReedMullerCode(base_ring(G), size(G, 2), size(G, 1), 2^(m - r), r, m, G, missing,
        H, missing, Gstand, Hstand, missing)
end

"""
    dual(C::ReedMullerCode)

Return the dual of the Reed-Muller code `C`.
"""
function dual(C::ReedMullerCode)
    # really only put this here to remind me later that I hardcoded binary
    if Int(characteristic(field(C))) == 2
        return ReedMullerCode(field(C), length(C), length(C) - dimension(C), 2^(order(C) + 1),
            RMm(C) - order(C) - 1, RMm(C), paritycheckmatrix(C), originalparitycheckmatrix(C), generatormatrix(C),
            originalgeneratormatrix(C), paritycheckmatrix(C, true), generatormatrix(C, true), missing)
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
    field(C) == field(D) || error("Codes must be over the same field in the Schur product.")
    length(C) == length(D) || error("Codes must have the same length in the Schur product.")

    r = order(C) + order(D)
    if r <= length(C)
        return ReedMullerCode(Int(characteristic(field(C))), r, RMm(C))
    else
        return ReedMullerCode(Int(characteristic(field(C))), RMm(C), RMm(C))
    end
end
*(C::ReedMullerCode, D::ReedMullerCode) = entrywiseproductcode(C, D)
Schurproductcode(C::ReedMullerCode, D::ReedMullerCode) = entrywiseproductcode(C, D)
Hadamardproductcode(C::ReedMullerCode, D::ReedMullerCode) = entrywiseproductcode(C, D)
componentwiseproductcode(C::ReedMullerCode, D::ReedMullerCode) = entrywiseproductcode(C, D)
