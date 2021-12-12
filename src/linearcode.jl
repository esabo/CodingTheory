# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

ENV["NEMO_PRINT_BANNER"] = "false"

using AbstractAlgebra
using Nemo

import Base: length, in, ⊆, /
import AbstractAlgebra: quo

include("utils.jl")

abstract type AbstractLinearCode end

mutable struct LinearCode <: AbstractLinearCode
    F::Union{FqNmodFiniteField, Nemo.GaloisField, AbstractAlgebra.GFField{Int64}}
    n::Integer
    k::Integer
    d::Union{Integer, Missing}
    G::Union{gfp_mat, fq_nmod_mat}
    Gorig::Union{gfp_mat, fq_nmod_mat, Missing}
    H::Union{gfp_mat, fq_nmod_mat}
    Horig::Union{gfp_mat, fq_nmod_mat, Missing}
    Gstand::Union{gfp_mat, fq_nmod_mat}
    Hstand::Union{gfp_mat, fq_nmod_mat}
end

# need to make sure this isn't doing column swaps
function _standardform(G::Union{gfp_mat, fq_nmod_mat})
    # hnf(C.G) shoud also work in most situations
    _, Gstand = rref(G)
    A = G[:, (size(G, 1) + 1):size(G, 2)]
    T = MatrixSpace(base_ring(G), size(G, 2) - size(G, 1), size(G, 2) - size(G, 1))
    Hstand = hcat(-transpose(A), T(1))
    return Gstand, Hstand
end

function LinearCode(G::Union{gfp_mat, fq_nmod_mat}, parity::Bool=false, verify::Bool=true)
    # if passed a VS, just ask for the basis
    !iszero(G) || error("Zero matrix passed into LinearCode constructor.")
    if typeof(base_ring(G)) != FqNmodFiniteField && typeof(base_ring(G)) != Nemo.GaloisField && typeof(base_ring(G)) != AbstractAlgebra.GFField{Int64}
        error("Codes over general rings are currently not supported, only fields.")
    end

    Gorig = G
    if rank(G) < size(G, 1)
        _, G = rref(G)
    end

    del = Vector{Int64}()
    for r in 1:size(G, 1)
        if iszero(G[r, :])
            append!(del, r)
        end
    end

    # not necessary but faster to skip
    if !isempty(del)
        G = G[setdiff(1:size(G, 1), del), :]
    end

    _, H = right_kernel(G)
    # note the H here is transpose of the standard definition
    !(!iszero(G * H) || !iszero(H' * G')) || error("Generator and parity check matrices are not transpose orthogonal.")
    Gstand, Hstand = _standardform(G)

    for r in 1:size(Gstand, 1)
        iszero(Gstand[r, :] * H) || error("Column swap appeared in _standardform.")
    end

    if !parity
        return LinearCode(base_ring(G), size(G, 2), size(G, 1), missing, G, Gorig,
            deepcopy(H'), missing, Gstand, Hstand)
    else
        return LinearCode(base_ring(G), size(G, 2), size(G, 1), missing, deepcopy(H'), missing,
            G, Gorig, Hstand, Gstand)
    end
end

function Base.show(io::IO, C::AbstractLinearCode)
    if get(io, :compact, false)
        println(io, "[$(length(C)), $(dimension(C))]_$(order(field(C))) linear code.")
    else
        println(io, "[$(length(C)), $(dimension(C))]_$(order(field(C))) linear code.")
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

field(C::AbstractLinearCode) = C.F
length(C::AbstractLinearCode) = C.n
n(C::AbstractLinearCode) = C.n
dimension(C::AbstractLinearCode) = C.k
k(C::AbstractLinearCode) = C.k
cardinality(C::AbstractLinearCode) = BigInt(order(field(C)))^dimension(C)
rate(C::AbstractLinearCode) = dimension(C) / length(C)

function minimumdistance(C::AbstractLinearCode)
    !ismissing(C.d) || error("Unknown minimum distance for this code.")
    return C.d
end
d(C::AbstractLinearCode) = minimumdistance(C)


function setminimumdistance!(C::T, d::Integer) where T <: AbstractLinearCode
    !(d <= 0) || error("The minimum distance of a code must be ≥ 1; received: d = $d.")
    C.d = d
end

function relativedistance(C::AbstractLinearCode)
    !ismissing(C.d) || error("Unknown minimum distance for this code.")
    return C.d / C.n
end

function generatormatrix(C::AbstractLinearCode, standform::Bool=false)
    !standform || return C.Gstand
    return C.G
end

function originalgeneratormatrix(C::AbstractLinearCode)
    !ismissing(C.Gorig) || return C.G # is ths the behavior I want?
    return C.Gorig
end

function paritycheckmatrix(C::AbstractLinearCode, standform::Bool=false)
    !standform || return C.Hstand
    return C.H
end

function originalparitycheckmatrix(C::AbstractLinearCode)
    !ismissing(C.Horig) || return C.H # is ths the behavior I want?
    return C.Horig
end

function genus(C::AbstractLinearCode)
    !ismissing(C.d) || error("Unknown minimum distance for this code.")
    return length(C) + 1 - dimension(C) - minimumdistance(C)
end

function Singletonbound(n::Integer, a::Integer)
    # d ≤ n - k + 1 or k ≤ n - d + 1
    if n >= 0 && a >= 0 && n >= a
        return n - a + 1
    else
        error("Invalid parameters for the Singleton bound. Received n = $n, k/d = $a")
    end
end

function Singletonbound(C::AbstractLinearCode)
    return Singletonbound(length(C), dimension(C))
end

function isMDS(C::AbstractLinearCode)
    minimumdistance(C) != Singletonbound(C) || return true
    return false
end

function numbercorrectableerrors(C::AbstractLinearCode)
    !ismissing(minimumdistance(C)) || return -1 # what do I want to return here?
    return Int64(floor((minimumdistance(C) - 1) / 2))
end

function encode(v::Union{gfp_mat, fq_nmod_mat}, C::AbstractLinearCode)
    !(size(v) != (1, dimension(C)) && size(v) != (dimension(C), 1)) ||
        error("Vector to be encoded is of incorrect dimension; expected length $(C.k), received: $(size(v)).")
    base_ring(v) == field(C) || error("Vector must have the same base ring as the generator matrix.")

    size(v, 1) != 1 || return v * generatormatrix(C)
    return v' * generatormatrix(C)
end

function encode(v::Vector{Int64}, C::AbstractLinearCode)
    length(v) == dimension(C) ||
        error("Vector to be encoded is of incorrect length; expected length $(C.k), received: $(size(v)).")
    return encode(matrix(field(C), v'), C)
end

function syndrome(v::Union{gfp_mat, fq_nmod_mat}, C::AbstractLinearCode)
    !(size(v) != (length(C), 1) && size(v) != (1, length(C))) ||
        error("Vector to be tested is of incorrect dimension; expected length $(C.n), received: $(size(v)).")
    base_ring(v) == field(C) || error("Vector must have the same base ring as the generator matrix.")

    size(v, 1) != 1 || return paritycheckmatrix(C) * v'
    return paritycheckmatrix(C) * v
end

function syndrome(v::Vector{Int64}, C::AbstractLinearCode)
    length(v) == length(C) ||
        error("Vector to be tested is of incorrect dimension; expected length $(C.n), received: $(size(v)).")
    return syndrome(matrix(field(C), v'), C)
end

# need to remove Base.:(in)(
in(v::Union{gfp_mat, fq_nmod_mat}, C::AbstractLinearCode) = iszero(syndrome(v, C))

function Base.:(⊆)(C1::AbstractLinearCode, C2::AbstractLinearCode)
    if field(C1) != field(C2) || length(C1) != length(C2) || dimension(C1) > dimension(C2)
        return false
    end

    # eachrow doesn't work on these objects
    for r in 1:dimension(C1)
        if C1.G[r, :] ∉ C2
            return false
        end
    end

    return true
end

issubcode(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊆ C2

# Credit to Tommy Hofmann of the AbstractAlgebra/Nemo/Hecke package for debugging
# and providing an elegant implementation here
function codecomplement(C1::AbstractLinearCode, C2::AbstractLinearCode)
    C1 ⊆ C2 || error("First code must be a subset of the second code.")
    F = field(C1)
    V = VectorSpace(F, length(C1))
    U, UtoV = sub(V, [V(C1.G[i, :]) for i in 1:nrows(C1.G)])
    W, WtoV = sub(V, [V(C2.G[i, :]) for i in 1:nrows(C2.G)])
    gensofUinW = [preimage(WtoV, UtoV(g)) for g in gens(U)]
    UinW, UinWtoW = sub(W, gensofUinW)
    Q, WtoQ = quo(W, UinW)
    C2modC1basis = [WtoV(x) for x in [preimage(WtoQ, g) for g in gens(Q)]]

    Fbasis = [[F(C2modC1basis[j][i]) for i in 1:dim(parent(C2modC1basis[1]))] for j in 1:length(C2modC1basis)]
    G = matrix(F, length(Fbasis), length(Fbasis[1]), vcat(Fbasis...))
    for r in 1:length(Fbasis)
        v = G[r, :]
        (v ∈ C2 && v ∉ C1) || error("Error in creation of basis for C2 / C1.")
    end
    return LinearCode(G)
end
quo(C1::AbstractLinearCode, C2::AbstractLinearCode) = codecomplement(C1, C2)
quotient(C1::AbstractLinearCode, C2::AbstractLinearCode) = codecomplement(C1, C2)
/(C2::S, C1::T) where S <: AbstractLinearCode where T <: AbstractLinearCode = codecomplement(C1, C2)
#Base.:(/)

function dual(C::AbstractLinearCode)
    return LinearCode(field(C), length(C), length(C) - dimension(C), missing,
        paritycheckmatrix(C), originalparitycheckmatrix(C), generatormatrix(C),
        originalgeneratormatrix(C), paritycheckmatrix(C, true), generatormatrix(C, true))
end

isequivalent(C1::AbstractLinearCode, C2::AbstractLinearCode) = (C1 ⊆ C2 && C2 ⊆ C1)
isselfdual(C::AbstractLinearCode) = isequivalent(C, dual(C))
isselforthogonal(C::AbstractLinearCode) = C ⊆ dual(C)

function ⊕(C1::AbstractLinearCode, C2::AbstractLinearCode)
    field(C1) == field(C2) || error("Codes must be over the same field in the direct sum.")

    G = directsum(generatormatrix(C1), generatormatrix(C2))
    Gstand, Hstand = _standardform(G)
    return LinearCode(field(C1), length(C1), dimension(C1), missing, G, missing,
        directsum(C1.H, C2.H), missing, Gstand, Hstand)
end

directsum(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊕ C2

function ⊗(C1::AbstractLinearCode, C2::AbstractLinearCode)
    field(C1) == field(C2) || error("Codes must be over the same field in the tensor product.")

    G = generatormatrix(C1) ⊗ generatormatrix(C2)
    Gstand, Hstand = _standardform(G)
    return LinearCode(field(C1), length(C1) * length(C2), dimension(C1) * dimension(C2),
        missing, G, missing, paritycheckmatrix(C1) ⊗ paritycheckmatrix(C2), missing,
        Gstand, Hstand)
end

kron(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊗ C2
tensorproduct(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊗ C2
directproduct(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊗ C2
productcode(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊗ C2

# # support is currently not implemented
# function support(C::AbstractLinearCode)
#     # returns the nonzero coefficient locations of the weight distribution
# end
#
# # support is currently not implemented
# function characteristicpolynomial(C::T) where T <: Union{LinearCode, CyclicCode}
#     # Return the characteristic polynomial of a linear code, as defined in [Lin1999]_.
#     R, x = PolynomialRing(Nemo.QQ, "x")
#     D = dual(C)
#     supD = support(D)
#     q = order(base_ring(C.G))
#     return q^(n - k) * prod([1 - x / j for j in supD if j > 0])
# end

function extend(C::AbstractLinearCode)
    # add an extra column to G such that the sum of the coordinates of each
    # row is 0
    # this is the most common form of extending a code

    col = [sum(C.G[i, :]) for i in 1:dimension(C)]
    M = MatrixSpace(field(C), dimension(C), 1)
    G = hcat(generatormatrix(C), M([characteristic(field(C)) - i for i in col]))
    M = MatrixSpace(field(C), 1, length(C) + 1)
    toprow = M([1 for _ in 1:(length(C) + 1)])
    M = MatrixSpace(field(C), length(C) - dimension(C), 1)
    rightcol = M(0)
    H = hcat(paritycheckmatrix(C), rightcol)
    H = vcat(toprow, H)
    Gstand, Hstand = _standardform(G)
    return LinearCode(field(C), size(G, 2), size(G, 1), missing,
        G, generatormatrix(C), H, paritycheckmatrix(C), Gstand, Hstand)
end

function puncture(C::AbstractLinearCode, cols::Vector{Int64})
    # delete columns from generator matrix and then remove any potentially resulting zero rows
    # to decode a punctured code it ~inserts zeros into punctured locations then uses original
    # decoder + erasures

    !isempty(cols) || return C

    G = C.G[:, setdiff(1:length(C), cols)]
    if rank(G) < size(G, 1)
        _, G = rref(G)
    end

    del = Vector{Int64}()
    for r in 1:size(G, 1)
        if iszero(G[r, :])
            append!(del, r)
        end
    end

    # not necessary but faster to skip
    if !isempty(del)
        G = G[setdiff(1:dimension(C), del), :]
    end
    _, H = right_kernel(G)
    # note the H here is transpose of the standard definition

    !(!iszero(G * H) || !iszero(H' * G')) ||
        error("Generator and parity check matrices are not transpose orthogonal.")

    Gstand, Hstand = _standardform(G)
    return LinearCode(field(C), size(G, 2), size(G, 1), missing, G, generatormatrix(C),
        H', paritycheckmatrix(C), Gstand, Hstand)
end

function expurgate(C::AbstractLinearCode, rows::Vector{Int64})
    # delete rows from generator matrix and then remove any potentially resulting zero columns

    !isempty(rows) || return C

    G = C.G[setdiff(1:dimension(C), rows), :]
    # instead of calling LinearCode(G) here I want to keep C.G has Gorig so I am
    # repeating this constructor here such that I have more control over that
    if rank(G) < size(G, 1)
        _, G = rref(G)
    end

    del = Vector{Int64}()
    for r in 1:size(G, 1)
        if iszero(G[r, :])
            append!(del, r)
        end
    end

    # not necessary but faster to skip
    if !isempty(del)
        G = G[setdiff(1:dimension(C), del), :]
    end
    _, H = right_kernel(G)
    # note the H here is transpose of the standard definition

    !(!iszero(G * H) || !iszero(H' * G')) ||
        error("Generator and parity check matrices are not transpose orthogonal.")

    Gstand, Hstand = _standardform(G)
    return LinearCode(field(C), size(G, 2), size(G, 1), missing, G, generatormatrix(C),
        H', paritycheckmatrix(C), Gstand, Hstand)
end

function augment(C::AbstractLinearCode, M::Union{gfp_mat, fq_nmod_mat})
    # adds rows (matrix M) to bottom of generator matrix of C

    length(C) == size(M, 2) || error("Rows to add must have the same number of columns as the generator matrix.")
    field(C) == base_ring(M) || error("Rows to augment must have the same base field as the code.")

    G = vcat(generatormatrix(C), M)
    # instead of calling LinearCode(G) here I want to keep C.G has Gorig so I am
    # repeating this constructor here such that I have more control over that
    if rank(G) < size(G, 1)
        _, G = rref(G)
    end

    del = Vector{Int64}()
    for r in 1:size(G, 1)
        if iszero(G[r, :])
            append!(del, r)
        end
    end

    # not necessary but faster to skip
    if !isempty(del)
        G = G[setdiff(1:dimension(C), del), :]
    end
    _, H = right_kernel(G)
    # note the H here is transpose of the standard definition

    !(!iszero(G * H) || !iszero(H' * G')) ||
        error("Generator and parity check matrices are not transpose orthogonal.")

    Gstand, Hstand = _standardform(G)
    return LinearCode(field(C), size(G, 2), size(G, 1), missing, G, generatormatrix(C),
        H', paritycheckmatrix(C), Gstand, Hstand)
end

function shorten(C::AbstractLinearCode, L::Vector{Int64})
    # expurgate followed by puncture
    # Using the fact that the shortened code C_L = ((C^⟂)^L)^⟂ where the exponent is
    # code punctured on L

    !isempty(L) || return C
    return dual(puncture(dual(C), L))
end

function lengthen(C::AbstractLinearCode)
    row = matrix(field(C), [1 for _ in 1:length(C)]')
    newcode = augment(C, row)
    return extend(newcode)
end

function uuplusv(C1::AbstractLinearCode, C2::AbstractLinearCode, verify::Bool=true)
    # Returns the Plotkin `(u|u+v)`-construction with u = C1 and v = C2
    field(C1) == field(C2) || error("Base field must be the same in the Plotkin (u|u + v)-construction.")
    length(C1) == length(C2) || error("Both codes must have the same length in the Plotkin (u|u + v)-construction.")

    G = vcat(hcat(generatormatrix(C1), generatormatrix(C1)), hcat(parent(generatormatrix(C2))(0), generatormatrix(C2)))
    H = vcat(hcat(paritycheckmatrix(C1), parent(paritycheckmatrix(C1))(0)), hcat(-paritycheckmatrix(C2), paritycheckmatrix(C2)))
    Gstand, Hstand = _standardform(G)

    if ismissing(C1.d) || ismissing(C2.d)
        d = missing
    else
        d = minimum([2 * minimumdistance(C1), minimumdistance(C2)])
    end

    if verify
        size(G, 2) == 2 * length(C1) ||
            error("Something went wrong in the Plotkin (u|u + v)-construction;
                length ($(size(G, 2))) and expected length ($(2 * length(C1))) are not equal.")
        size(G, 1) == dimension(C1) + dimension(C2) ||
            error("Something went wrong in the Plotkin (u|u + v)-construction;
                dimension ($(size(G, 1))) and expected dimension ($(dimension(C1) + dimension(C2))) are not equal.")
    end

    return LinearCode(field(C1), size(G, 2), size(G, 1), d, G, missing, H, missing, Gstand, Hstand)
end
Plotkinconstruction(C1::AbstractLinearCode, C2::AbstractLinearCode) = uuplusv(C1, C2)

# function tracecode(C::AbstractLinearCode)
#
# end

# function subfieldsubcode(C::AbstractLinearCode)
#     # return dual(tracecode(dual(C)))
# end
