# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

using AbstractAlgebra
using Nemo

import Base: length#, ⊆

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

function standardform(G::Union{gfp_mat, fq_nmod_mat})
    # hnf(C.G) shoud also work in most situations
    _, Gstand = rref(G)
    A = G[:, (size(G, 1) + 1):size(G, 2)]
    T = MatrixSpace(base_ring(G), size(G, 2) - size(G, 1), size(G, 2) - size(G, 1))
    Hstand = hcat(-transpose(A), T(1))
    return Gstand, Hstand
end

function LinearCode(G::Union{gfp_mat, fq_nmod_mat}, parity::Bool=false, verify::Bool=true)
    # if passed a VS, just ask for the basis
    if iszero(G)
        error("Zero matrix passed into LinearCode constructor.")
    end

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

    if !iszero(G * H) || !iszero(H' * G')
        error("Generator and parity check matrices are not transpose orthogonal.")
    end

    Gstand, Hstand = standardform(G)

    if !parity
        return LinearCode(base_ring(G), size(G, 2), size(G, 1), missing, G, Gorig,
            H', missing, Gstand, Hstand)
    else
        return LinearCode(base_ring(G), size(G, 2), size(G, 1), missing, H', missing,
            G, Gorig, Hstand, Gstand)
    end
end

function Base.show(io::IO, C::T) where T <: AbstractLinearCode
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

function field(C::T) where T <: AbstractLinearCode
    return C.F
end

function length(C::T) where T <: AbstractLinearCode
    return C.n
end

function dimension(C::T) where T <: AbstractLinearCode
    return C.k
end

function minimumdistance(C::T) where T <: AbstractLinearCode
    if ismissing(C.d)
        error("Unknown minimum distance for this code.")
    else
        return C.d
    end
end

function cardinality(C::T) where T <: AbstractLinearCode
    return BigInt(order(field(C)))^dimension(C)
end

function rate(C::T) where T <: AbstractLinearCode
    return dimension(C) / length(C)
end

function setminimumdistance!(C::T, d::Integer) where T <: AbstractLinearCode
    if d <= 0
        error("The minimum distance of a code must be ≥ 1; received: d = $d.")
    end
    C.d = d
end

function relativedistance(C::T) where T <: AbstractLinearCode
    if ismissing(C.d)
        error("Unknown minimum distance for this code.")
    else
        return C.d / C.n
    end
end

function generatormatrix(C::T, standform::Bool=false) where T <: AbstractLinearCode
    if standform
        return C.Gstand
    else
        return C.G
    end
end

function originalgeneratormatrix(C::T) where T <: AbstractLinearCode
    if ismissing(C.Gorig)
        return C.G # is ths the behavior I want?
    else
        return C.Gorig
    end
end

function paritycheckmatrix(C::T, standform::Bool=false) where T <: AbstractLinearCode
    if standform
        return C.Hstand
    else
        return C.H
    end
end

function originalparitycheckmatrix(C::T) where T <: AbstractLinearCode
    if ismissing(C.Horig)
        return C.H # is ths the behavior I want?
    else
        return C.Horig
    end
end

function genus(C::T) where T <: AbstractLinearCode
    if ismissing(C.d)
        error("Unknown minimum distance for this code.")
    else
        return length(C) + 1 - dimension(C) - minimumdistance(C)
    end
end

function Singletonbound(n::Integer, a::Integer)
    # d ≤ n - k + 1 or k ≤ n - d + 1
    if n >= 0 && a >= 0 && n >= a
        return n - a + 1
    else
        error("Invalid parameters for the Singleton bound. Received n = $n, k/d = $a")
    end
end

function Singletonbound(C::T) where T <: AbstractLinearCode
    return Singletonbound(length(C), dimension(C))
end

function isMDS(C::T) where T <: AbstractLinearCode
    if minimumdistance(C) == Singletonbound(C)
        return true
    else
        return false
    end
end

function numbercorrectableerrors(C::T) where T <: AbstractLinearCode
    if ismissing(minimumdistance(C))
        return -1 # what do I want to return here?
    else
        return Int64(floor((minimumdistance(C) - 1) / 2))
    end
end

function encode(v::Union{gfp_mat, fq_nmod_mat}, C::T) where T <: AbstractLinearCode
    if size(v) != (1, dimension(C)) && size(v) != (dimension(C), 1)
        error("Vector to be encoded is of incorrect dimension; expected length $(C.k), received: $(size(v)).")
    end

    if base_ring(v) != field(C)
        error("Vector must have the same base ring as the generator matrix.")
    end

    if size(v, 1) == 1
        return v * generatormatrix(C)
    else
        return v' * generatormatrix(C)
    end
end

function encode(v::Vector{Int64}, C::T) where T <: AbstractLinearCode
    if length(v) != dimension(C)
        error("Vector to be encoded is of incorrect length; expected length $(C.k), received: $(size(v)).")
    end

    return encode(matrix(field(C), v'), C)
end

function syndrome(v::Union{gfp_mat, fq_nmod_mat}, C::T) where T <: AbstractLinearCode
    if size(v) != (length(C), 1) && size(v) != (1, length(C))
        error("Vector to be tested is of incorrect dimension; expected length $(C.n), received: $(size(v)).")
    end

    if base_ring(v) != field(C)
        error("Vector must have the same base ring as the generator matrix.")
    end

    if size(v, 1) == 1
        return paritycheckmatrix(C) * v'
    else
        return paritycheckmatrix(C) * v
    end
end

function syndrome(v::Vector{Int64}, C::T) where T <: AbstractLinearCode
    if length(v) != length(C)
        error("Vector to be tested is of incorrect dimension; expected length $(C.n), received: $(size(v)).")
    end

    return syndrome(matrix(field(C), v'), C)
end

function Base.:(in)(v::Union{gfp_mat, fq_nmod_mat}, C::T) where T <: AbstractLinearCode
    # should cast this zero up but probably faster to let it cast down
    return syndrome(v, C) == 0
end

function Base.:(⊆)(C1::S, C2::T)  where S <: AbstractLinearCode where T <: AbstractLinearCode
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

function issubcode(C1::S, C2::T)  where S <: AbstractLinearCode where T <: AbstractLinearCode
    return C1 ⊆ C2
end

function dual(C::T) where T <: AbstractLinearCode
    return LinearCode(field(C), length(C), length(C) - dimension(C), missing,
        paritycheckmatrix(C), originalparitycheckmatrix(C), generatormatrix(C),
        originalgeneratormatrix(C), paritycheckmatrix(C, true), generatormatrix(C, true))
end

function isequivalent(C1::S, C2::T)  where S <: AbstractLinearCode where T <: AbstractLinearCode
    return C1 ⊆ C2 && C2 ⊆ C1
end

function isselfdual(C::T) where T <: AbstractLinearCode
    return isequivalent(C, dual(C))
end

function isselforthogonal(C::T) where T <: AbstractLinearCode
    # A code is self-orthogonal if it is a subcode of its dual.
    return C ⊆ dual(C)
end

function ⊕(C1::S, C2::T)  where S <: AbstractLinearCode where T <: AbstractLinearCode
    if field(C1) != field(C2)
        error("Codes must be over the same field in the direct sum.")
    end

    G = directsum(generatormatrix(C1), generatormatrix(C2))
    Gstand, Hstand = standardform(G)
    return LinearCode(field(C1), length(C1), dimension(C1), missing, G, missing,
        directsum(C1.H, C2.H), missing, Gstand, Hstand)
end

function directsum(C1::S, C2::T)  where S <: AbstractLinearCode where T <: AbstractLinearCode
    return C1 ⊕ C2
end

function ⊗(C1::S, C2::T)  where S <: AbstractLinearCode where T <: AbstractLinearCode
    if field(C1) != field(C2)
        error("Codes must be over the same field in the tensor product.")
    end

    G = generatormatrix(C1) ⊗ generatormatrix(C2)
    Gstand, Hstand = standardform(G)
    return LinearCode(field(C1), length(C1) * length(C2), dimension(C1) * dimension(C2),
        missing, G, missing, paritycheckmatrix(C1) ⊗ paritycheckmatrix(C2), missing,
        Gstand, Hstand)
end

function kron(C1::S, C2::T)  where S <: AbstractLinearCode where T <: AbstractLinearCode
    return C1 ⊗ C2
end

function tensorproduct(C1::S, C2::T)  where S <: AbstractLinearCode where T <: AbstractLinearCode
    return C1 ⊗ C2
end

function directproduct(C1::S, C2::T)  where S <: AbstractLinearCode where T <: AbstractLinearCode
    return C1 ⊗ C2
end

function productcode(C1::S, C2::T)  where S <: AbstractLinearCode where T <: AbstractLinearCode
    return C1 ⊗ C2
end


# # support is currently not implemented
# function support(C::T) where T <: AbstractLinearCode
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

function extend(C::T) where T <: AbstractLinearCode
    # add an extra column to G such that the sum of the coordinates of each
    # row is 0
    # this is the most common form of extending a code

    col = [sum(C.G[i, :]) for i in 1:dimension(C)]
    M = MatrixSpace(field(C), dimension(C), 1)
    G = hcat(generatormatrix(C), M([characteristic(field(C)) - i for i in col]))
    M = MatrixSpace(field(C), 1, length(C) + 1)
    toprow = M([1 for i in 1:(length(C) + 1)])
    M = MatrixSpace(field(C), length(C) - dimension(C), 1)
    rightcol = M(0)
    H = hcat(paritycheckmatrix(C), rightcol)
    H = vcat(toprow, H)
    Gstand, Hstand = standardform(G)
    return LinearCode(field(C), size(G, 2), size(G, 1), missing,
        G, generatormatrix(C), H, paritycheckmatrix(C), Gstand, Hstand)
end

function puncture(C::T, cols::Vector{Int64}) where T <: AbstractLinearCode
    # delete columns from generator matrix and then remove any potentially resulting zero rows
    # to decode a punctured code it ~inserts zeros into punctured locations then uses original
    # decoder + erasures

    if isempty(cols)
        return C
    end

    G = C.G[:, setdiff(1:length(C), cols)]
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

    if !iszero(G * H) || !iszero(H' * G')
        error("Generator and parity check matrices are not transpose orthogonal.")
    end

    Gstand, Hstand = standardform(G)
    return LinearCode(field(C), size(G, 2), size(G, 1), missing, G, generatormatrix(C),
        H', paritycheckmatrix(C), Gstand, Hstand)
end

function expurgate(C::T, rows::Vector{Int64}) where T <: AbstractLinearCode
    # delete rows from generator matrix and then remove any potentially resulting zero columns

    if isempty(rows)
        return C
    end

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

    if !iszero(G * H) || !iszero(H' * G')
        error("Generator and parity check matrices are not transpose orthogonal.")
    end

    Gstand, Hstand = standardform(G)
    return LinearCode(field(C), size(G, 2), size(G, 1), missing, G, generatormatrix(C),
        H', paritycheckmatrix(C), Gstand, Hstand)
end

function augment(C::T, M::Union{gfp_mat, fq_nmod_mat}) where T <: AbstractLinearCode
    # adds rows (matrix M) to bottom of generator matrix of C

    if length(C) != length(M)
        error("Rows to add must have the same number of columns as the generator matrix.")
    end

    if field(C) != base_ring(M)
        error("Rows to augment must have the same base field as the code.")
    end

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

    if !iszero(G * H) || !iszero(H' * G')
        error("Generator and parity check matrices are not transpose orthogonal.")
    end

    Gstand, Hstand = standardform(G)
    return LinearCode(field(C), size(G, 2), size(G, 1), missing, G, generatormatrix(C),
        H', paritycheckmatrix(C), Gstand, Hstand)
end

function shorten(C::T, L::Vector{Int64}) where T <: AbstractLinearCode
    # expurgate followed by puncture
    # Using the fact that the shortened code C_L = ((C^⟂)^L)^⟂ where the exponent is
    # code punctured on L

    if isempty(L)
        return C
    end

    return dual(puncture(dual(C), L))
end

# should save these objects and then create a new LinearCode carefully which C.G as G original
function lengthen(C::T) where T <: AbstractLinearCode
    row = matrix(field(C), [1 for i in 1:length(C)]')
    newcode = augment(C, row)
    return extend(newcode)
end

function uuplusv(C1::S, C2::T, verify::Bool=true)  where S <: AbstractLinearCode where T <: AbstractLinearCode
    # Returns the Plotkin `(u|u+v)`-construction with u = C1 and v = C2
    if field(C1) != field(C2)
        error("Base field must be the same in the Plotkin (u|u + v)-construction.")
    end

    if length(C1) != length(C2)
        error("Both codes must have the same length in the Plotkin (u|u + v)-construction.")
    end

    G = vcat(hcat(generatormatrix(C1), generatormatrix(C1)), hcat(parent(generatormatrix(C2))(0), generatormatrix(C2)))
    H = vcat(hcat(paritycheckmatrix(C1), parent(paritycheckmatrix(C1))(0)), hcat(-paritycheckmatrix(C2), paritycheckmatrix(C2)))
    Gstand, Hstand = standardform(G)

    if ismissing(C1.d) || ismissing(C2.d)
        d = missing
    else
        d = minimum([2 * minimumdistance(C1), minimumdistance(C2)])
    end

    if verify
        if size(G, 2) != 2 * length(C1)
            error("Something went wrong in the Plotkin (u|u + v)-construction;
                length ($(size(G, 2))) and expected length ($(2 * length(C1))) are not equal.")
        end
        if size(G, 1) != dimension(C1) + dimension(C2)
            error("Something went wrong in the Plotkin (u|u + v)-construction;
                dimension ($(size(G, 1))) and expected dimension ($(dimension(C1) + dimension(C2))) are not equal.")
        end
    end

    return LinearCode(field(C1), size(G, 2), size(G, 1), d, G, missing, H, missing, Gstand, Hstand)
end

function Plotkinconstruction(C1::S, C2::T) where S <: AbstractLinearCode where T <: AbstractLinearCode
    # Returns the Plotkin `(u|u+v)`-construction with u = C1 and v = C2)
    return uuplusv(C1, C2)
end

# function tracecode(C::T) where T <: AbstractLinearCode
#
# end

# function subfieldsubcode(C::T) where T <: AbstractLinearCode
#     # return dual(tracecode(dual(C)))
# end
