# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.


# can't get around including this here and not in weight_dist without modules
struct WeightEnumerator
    polynomial::Vector{Vector{Int64}}
    Type::String
end

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
    weightenum::Union{WeightEnumerator, Missing}
end

# need to make sure this isn't doing column swaps
function _standardform(G::Union{gfp_mat, fq_nmod_mat})
    # hnf(C.G) shoud also work in most situations
    _, Gstand = rref(G) # breaks on this line
    A = Gstand[:, (size(Gstand, 1) + 1):size(Gstand, 2)]
    T = MatrixSpace(base_ring(Gstand), size(Gstand, 2) - size(Gstand, 1),
        size(Gstand, 2) - size(Gstand, 1))
    Hstand = hcat(-transpose(A), T(1))
    return Gstand, Hstand
end

"""
    LinearCode(G::fq_nmod_mat, parity::Bool=false, verify::Bool=true)

Return the linear code constructed with generator matrix `G`.

If `G` is not full rank, a row-reduced form is computed for the generator matrix.
The dimension of the code is the number of rows of the full-rank matrix, and the
length the number of columns. If the optional paramater `parity` is set to `true`,
a LinearCode is built with `G` as the parity check matrix. If the optional parameter
`verify` is set to `true`, basic checks, such as checking `GH == 0` and checking
for column swaps in the standard form, are done to ensure correctness.

# Arguments
* `G`: a matrix over a finite field of type `FqNmodFiniteField`
* `parity`: a boolean
* `verify`: a boolean

# Notes
* At the moment, no convention is used for G = 0 and an error is thrown.
* Zero columns are not removed.
* Row reduction is based on Nemo's rref function. It should be monitored to make
  sure that this does not introduce column swapping in a future version.

# Examples
```julia
julia> F, _ = FiniteField(2, 1, "α");
julia> G = matrix(F, 4, 7, [1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1]);
julia> C = LinearCode(G)
[7, 4]_2 linear code.
Generator matrix: 4 × 7
        1 0 0 0 0 1 1
        0 1 0 0 1 0 1
        0 0 1 0 1 1 0
        0 0 0 1 1 1 1
```
"""
function LinearCode(G::Union{gfp_mat, fq_nmod_mat}, parity::Bool=false, verify::Bool=true)
    # if passed a VS, just ask for the basis
    !iszero(G) || error("Zero matrix passed into LinearCode constructor.")
    if typeof(base_ring(G)) != FqNmodFiniteField && typeof(base_ring(G)) != Nemo.GaloisField && typeof(base_ring(G)) != AbstractAlgebra.GFField{Int64}
        error("Codes over general rings are currently not supported, only fields.")
    end

    Gorig = G
    rk = rank(G)
    !iszero(rk) || error("Rank zero matrix passed into LinearCode constructor.")
    if rk < size(G, 1)
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

    # note the H here is transpose of the standard definition
    _, H = right_kernel(G)

    Gstand, Hstand = _standardform(G)
    if verify
        size(transpose(H), 1) == size(G, 2) - size(G, 1) || error("Parity check matrix is not size n - k.")
        !(!iszero(G * H) || !iszero(transpose(H) * transpose(G))) || error("Generator and parity check matrices are not transpose orthogonal.")
        for r in 1:size(Gstand, 1)
            iszero(Gstand[r, :] * H) || error("Column swap appeared in _standardform.")
        end
    end

    if !parity
        return LinearCode(base_ring(G), size(G, 2), size(G, 1), missing, G, Gorig,
            deepcopy(transpose(H)), missing, Gstand, Hstand, missing)
    else
        return LinearCode(base_ring(G), size(G, 2), size(G, 1), missing, deepcopy(transpose(H)), missing,
            G, Gorig, Hstand, Gstand, missing)
    end
end

# need to unify notation here to call minimumdistance but since this was removed
# from this file and put into weight_dist, keep this until exporting
function show(io::IO, C::AbstractLinearCode)
    if get(io, :compact, false)
        if ismissing(C.d)
            println(io, "[$(length(C)), $(dimension(C))]_$(order(field(C))) linear code.")
        else
            println(io, "[$(length(C)), $(dimension(C)), $(C.d)]_$(order(field(C))) linear code.")
        end
    else
        if ismissing(C.d)
            println(io, "[$(length(C)), $(dimension(C))]_$(order(field(C))) linear code.")
        else
            println(io, "[$(length(C)), $(dimension(C)), $(C.d)]_$(order(field(C))) linear code.")
        end
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
    field(C::AbstractLinearCode)

Return the base ring of the generator matrix.
"""
field(C::AbstractCode) = C.F

"""
    length(C::AbstractLinearCode)

Return the length of the code.
"""
length(C::AbstractCode) = C.n

"""
    dimension(C::AbstractLinearCode)

Return the dimension of the code.
"""
dimension(C::AbstractLinearCode) = C.k

"""
    cardinality(C::AbstractLinearCode)

Return the cardinality of the code.

No size checking is done on the parameters of the code, returns a BitInt by default.
"""
cardinality(C::AbstractLinearCode) = BigInt(order(field(C)))^dimension(C)

"""
    rate(C::AbstractLinearCode)

Return the rate, `R = k/n', of the code.
"""
rate(C::AbstractLinearCode) = dimension(C) / length(C)

# """
#     minimumdistance(C::AbstractLinearCode)
#
# Return the minimum distance of the linear code if known, otherwise returns missing.
# """
# minimumdistance(C::AbstractLinearCode) = C.d

"""
    setminimumdistance(C::AbstractLinearCode, d::Integer)

Set the minimum distance of the code to `d`.

The only check done on the value of `d` is that `1 ≤ d ≤ n`.
"""
function setminimumdistance!(C::AbstractLinearCode, d::Integer)
    d > 0 && d <= length(C) || error("The minimum distance of a code must be ≥ 1; received: d = $d.")
    C.d = d
end

"""
    relativedistance(C::AbstractLinearCode)

Return the relative minimum distance, `δ = d / n` of the code if `d` is known,
otherwise errors.
"""
function relativedistance(C::AbstractLinearCode)
    !ismissing(C.d) || error("Unknown minimum distance for this code.")
    return C.d / C.n
end

"""
    generatormatrix(C::AbstractLinearCode, standform::Bool=false)

Return the generator matrix of the code.

If the optional parameter `standform` is set to `true`, the standard form of the
generator matrix is returned instead.
"""
function generatormatrix(C::AbstractLinearCode, standform::Bool=false)
    !standform || return C.Gstand
    return C.G
end

"""
    originalgeneratormatrix(C::AbstractLinearCode)

Return the original generator matrix used to create the code.

This should not be interpreted as any kind of generator to the code. For example,
if `G` is not a full-rank matrix, then this function acting on `C = LinearCode(G)`
will return `G` instead of the generator matrix for `C`. This may be used to diagnose
or understand the resulting code after its creation. An important use of this
function is in connection with standard code modification procedures such as
extending, puncturing, shortening, etc where this returns the matrix on which these
procedures have been applied.
"""
function originalgeneratormatrix(C::AbstractLinearCode)
    !ismissing(C.Gorig) || return C.G # is ths the behavior I want?
    return C.Gorig
end

"""
    paritycheckmatrix(C::AbstractLinearCode, standform::Bool=false)

Return the parity check matrix of the code.

If the optional parameter `standform` is set to `true`, the standard form of the
parity check matrix is returned instead.
"""
function paritycheckmatrix(C::AbstractLinearCode, standform::Bool=false)
    !standform || return C.Hstand
    return C.H
end

"""
    originalparitycheckmatrix(C::AbstractLinearCode)

Return the original parity check matrix used to create the code.

See `originalgeneratormatrix`.
"""
function originalparitycheckmatrix(C::AbstractLinearCode)
    !ismissing(C.Horig) || return C.H # is ths the behavior I want?
    return C.Horig
end

"""
    genus(C::AbstractLinearCode)

Return the genus, `n + 1 - k - d`, of the code.
"""
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

"""
    isMDS(C::AbstractLinearCode)

Return `true` if code is maximum distance separable (MDS).

A linear code is MDS if it saturates the Singleton bound, `d = n - k + 1`.
"""
function isMDS(C::AbstractLinearCode)
    minimumdistance(C) != Singletonbound(C) || return true
    return false
end

"""
    numbercorrectableerrors(C::AbstractLinearCode)

Return the number of correctable errors for the code.

The number of correctable errors is `t = floor((d - 1) / 2)`.
"""
function numbercorrectableerrors(C::AbstractLinearCode)
    !ismissing(minimumdistance(C)) || return -1 # what do I want to return here?
    return Int64(floor((minimumdistance(C) - 1) / 2))
end

"""
    encode(v::Union{fq_nmod_mat, Vector{Int64}}, C::AbstractLinearCode)

Return `v * G`, where `G` is the generator matrix of `C`.

# Arguments
* `v`: Either a `1 × k` or a `k × 1` vector.
"""
function encode(v::Union{gfp_mat, fq_nmod_mat}, C::AbstractLinearCode)
    !(size(v) != (1, dimension(C)) && size(v) != (dimension(C), 1)) ||
        error("Vector to be encoded is of incorrect dimension; expected length $(C.k), received: $(size(v)).")
    base_ring(v) == field(C) || error("Vector must have the same base ring as the generator matrix.")

    size(v, 1) != 1 || return v * generatormatrix(C)
    return transpose(v) * generatormatrix(C)
end

function encode(v::Vector{Int64}, C::AbstractLinearCode)
    length(v) == dimension(C) ||
        error("Vector to be encoded is of incorrect length; expected length $(C.k), received: $(size(v)).")
    return encode(matrix(field(C), transpose(v)), C)
end

"""
    syndrome(v::Union{fq_nmod_mat, Vector{Int64}}, C::AbstractLinearCode)

Return `H * v`, where `H` is the parity check matrix of `C`.

# Arguments
* `v`: Either a `1 × k` or a `k × 1` vector.
"""
function syndrome(v::Union{gfp_mat, fq_nmod_mat}, C::AbstractLinearCode)
    !(size(v) != (length(C), 1) && size(v) != (1, length(C))) ||
        error("Vector to be tested is of incorrect dimension; expected length $(C.n), received: $(size(v)).")
    base_ring(v) == field(C) || error("Vector must have the same base ring as the generator matrix.")

    size(v, 1) != 1 || return paritycheckmatrix(C) * transpose(v)
    return paritycheckmatrix(C) * v
end

function syndrome(v::Vector{Int64}, C::AbstractLinearCode)
    length(v) == length(C) ||
        error("Vector to be tested is of incorrect dimension; expected length $(C.n), received: $(size(v)).")
    return syndrome(matrix(field(C), transpose(v)), C)
end

"""
    in(v::Union{fq_nmod_mat, Vector{Int64}}, C::AbstractLinearCode)

Return whether or not `v` is a codeword of `C`.

The vector `v` is a valid codeword of `C` if and only if the syndrome of `v` is zero.
"""
in(v::Union{gfp_mat, fq_nmod_mat}, C::AbstractLinearCode) = iszero(syndrome(v, C))
in(v::Vector{Int64}, C::AbstractLinearCode) = iszero(syndrome(v, C))

"""
    ⊆(C1::AbstractLinearCode, C2::AbstractLinearCode)
    ⊂(C1::AbstractLinearCode, C2::AbstractLinearCode)
    issubcode(C1::AbstractLinearCode, C2::AbstractLinearCode)

Return whether or not `C1` is a subcode of `C2`.

A code `C1` is a subcode of another code `C2` if each row of
the generator matrix of `C1` is a valid codeword of `C2`.
"""
function ⊆(C1::AbstractLinearCode, C2::AbstractLinearCode)
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
⊂(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊆ C2
issubcode(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊆ C2

"""
    codecomplement(C1::AbstractLinearCode, C2::AbstractLinearCode)
    quo(C1::AbstractLinearCode, C2::AbstractLinearCode)
    quotient(C1::AbstractLinearCode, C2::AbstractLinearCode)
    /(C2::AbstractLinearCode, C1::AbstractLinearCode)

Return the code `C2 / C1` given `C1 ⊆ C2`.

Credit to Tommy Hofmann of the AbstractAlgebra/Nemo/Hecke packages for help with
debugging and providing the most elegant implementation now used here.
"""
function codecomplement(C1::AbstractLinearCode, C2::AbstractLinearCode)
    C1 ⊆ C2 || error("Need C1 ⊆ C2 in codecomplement(C1, C2) or C2 / C1.")
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
/(C2::AbstractLinearCode, C1::AbstractLinearCode) = codecomplement(C1, C2)

"""
    dual(C::AbstractLinearCode)

Return the dual of the code `C`.

All of the necessary information for the dual is already stored in a LinearCode
object, so this implementation merely swaps data fields, e.g., `G <-> H`, without
doing any new computation.
"""
function dual(C::AbstractLinearCode)
    return LinearCode(field(C), length(C), length(C) - dimension(C), missing,
        deepcopy(paritycheckmatrix(C)), deepcopy(originalparitycheckmatrix(C)),
        deepcopy(generatormatrix(C)), deepcopy(originalgeneratormatrix(C)),
        deepcopy(paritycheckmatrix(C, true)), deepcopy(generatormatrix(C, true)),
        missing)
end

"""
    Hermitiandual(C::AbstractLinearCode)

Return the Hermitian dual of a code defined over a quadratic extension.
"""
Hermitiandual(C::AbstractLinearCode) = dual(LinearCode(Hermitianconjugatematrix(generatormatrix(C))))
# also valid to do LinearCode(Hermitianconjugatematrix(generatormatrix(dual(C))))

"""
    isequivalent(C1::AbstractLinearCode, C2::AbstractLinearCode)

Return `true` if `C1 ⊆ C2` and `C2 ⊆ C1`.
"""
isequivalent(C1::AbstractLinearCode, C2::AbstractLinearCode) = (C1 ⊆ C2 && C2 ⊆ C1)

"""
    isselfdual(C::AbstractLinearCode)

Return `true` if `isequivalent(C, dual(C))`.
"""
isselfdual(C::AbstractLinearCode) = isequivalent(C, dual(C))

"""
    isselforthogonal(C::AbstractLinearCode)
    isweaklyselfdual(C::AbstractLinearCode)

Return `true` if `C ⊆ dual(C)`.
"""
isselforthogonal(C::AbstractLinearCode) = C ⊆ dual(C)
isweaklyselfdual(C::AbstractLinearCode) = isselforthogonal(C)

"""
    ⊕(C1::AbstractLinearCode, C2::AbstractLinearCode)
    directsum(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊕ C2

Return the direct sum code of `C1` and `C2`.

The direct sum code has generator matrix `G1 ⊕ G2` and parity check matrix `H1 ⊕ H2`.
"""
function ⊕(C1::AbstractLinearCode, C2::AbstractLinearCode)
    field(C1) == field(C2) || error("Codes must be over the same field in the direct sum.")

    G = directsum(generatormatrix(C1), generatormatrix(C2))
    # this should just be the direct sum of these variables as well
    Gstand, Hstand = _standardform(G)
    return LinearCode(field(C1), length(C1), dimension(C1), missing, G, missing,
        directsum(C1.H, C2.H), missing, Gstand, Hstand, missing)
end
directsum(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊕ C2

"""
    ⊗(C1::AbstractLinearCode, C2::AbstractLinearCode)
    kron(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊗ C2
    tensorproduct(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊗ C2
    directproduct(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊗ C2
    productcode(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊗ C2

Return the (direct/tensor) product code of `C1` and `C2`.

The product code has generator matrix `G1 ⊗ G2` and parity check matrix `H1 ⊗ H2`.

# Notes
* The resulting product is not checked for any zero columns.
"""
function ⊗(C1::AbstractLinearCode, C2::AbstractLinearCode)
    field(C1) == field(C2) || error("Codes must be over the same field in the tensor product.")

    G = generatormatrix(C1) ⊗ generatormatrix(C2)
    Gstand, Hstand = _standardform(G)
    return LinearCode(field(C1), length(C1) * length(C2), dimension(C1) * dimension(C2),
        missing, G, missing, paritycheckmatrix(C1) ⊗ paritycheckmatrix(C2), missing,
        Gstand, Hstand, missing)
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

"""
    extend(C::AbstractLinearCode)

Return the extended code of `C`.

This implementation chooses the most common form of extending a code, which is to
add an extra column to the generator matrix such that the sum of the coordinates
of each row is 0.
"""
function extend(C::AbstractLinearCode)
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
        G, generatormatrix(C), H, paritycheckmatrix(C), Gstand, Hstand, missing)
end

"""
    puncture(C::AbstractLinearCode, cols::Vector{Int64})

Return the code of `C` punctured at the columns in `cols`.

Deletes the columns from the generator matrix and then removes any potentially
resulting zero rows.
"""
function puncture(C::AbstractLinearCode, cols::Vector{Int64})
    # to decode a punctured code ~kinda insert zeros into punctured locations
    # then uses original decoder + erasures

    !isempty(cols) || return C
    cols ⊆ 1:length(C) || error("Columns to puncture are not a subset of the index set.")
    length(cols) != length(C) || error("Cannot puncture all columns of a generator matrix.")

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

    !(!iszero(G * H) || !iszero(transpose(H) * transpose(G))) ||
        error("Generator and parity check matrices are not transpose orthogonal.")

    Gstand, Hstand = _standardform(G)
    return LinearCode(field(C), size(G, 2), size(G, 1), missing, G, generatormatrix(C),
        transpose(H), paritycheckmatrix(C), Gstand, Hstand, missing)
end

"""
    expurgate(C::AbstractLinearCode, rows::Vector{Int64})

Return the code of `C` expuragated at the rows in `rows`.

Deletes the rows from the generator matrix and then removes any potentially
resulting zero columns.
"""
function expurgate(C::AbstractLinearCode, rows::Vector{Int64})
    !isempty(rows) || return C
    rows ⊆ 1:dimension(C) || error("Rows to expurgate are not a subset of the index set.")
    length(rows) != dimension(C) || error("Cannot expurgate all rows of a generator matrix.")

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

    !(!iszero(G * H) || !iszero(transpose(H) * transpose(G))) ||
        error("Generator and parity check matrices are not transpose orthogonal.")

    Gstand, Hstand = _standardform(G)
    return LinearCode(field(C), size(G, 2), size(G, 1), missing, G, generatormatrix(C),
        transpose(H), paritycheckmatrix(C), Gstand, Hstand, missing)
end

"""
    augment(C::AbstractLinearCode, M::fq_nmod_mat)

Return the code of `C` whose generator matrix is augmented with `M`.

Vertically joins the matrix `M` to the bottom of the generator matrix of `C`.
"""
function augment(C::AbstractLinearCode, M::Union{gfp_mat, fq_nmod_mat})
    length(C) == size(M, 2) || error("Rows to augment must have the same number of columns as the generator matrix.")
    field(C) == base_ring(M) || error("Rows to augment must have the same base field as the code.")

    G = vcat(generatormatrix(C), M)
    # instead of calling LinearCode(G) here I want to keep C.G has Gorig so I am
    # repeating this constructor here such that I have more control over that
    rk = rank(G)
    !iszero(rk) || error("Rank zero matrix passed into LinearCode constructor.")
    if rk < size(G, 1)
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

    !(!iszero(G * H) || !iszero(transpose(H) * transpose(G))) ||
        error("Generator and parity check matrices are not transpose orthogonal.")

    Gstand, Hstand = _standardform(G)
    for r in 1:size(Gstand, 1)
        iszero(Gstand[r, :] * H) || error("Column swap appeared in _standardform.")
    end

    return LinearCode(field(C), size(G, 2), size(G, 1), missing, G, generatormatrix(C),
        transpose(H), paritycheckmatrix(C), Gstand, Hstand, missing)
end

"""
    shorten(C::AbstractLinearCode, L::Vector{Int64})

Return the code of `C` shortened on the indices `L`.

Shortening is expurgating followed by puncturing. This implementation uses the
theorem that the dual of code shortened on `L` is equal to the puncture of the
dual code on `L`, i.e., `dual(puncture(dual(C), L))`.
"""
function shorten(C::AbstractLinearCode, L::Vector{Int64})
    # most parameter checks here done in puncture
    !isempty(L) || return C
    return dual(puncture(dual(C), L))
end

"""
    lengthen(C::AbstractLinearCode)

Return the lengthened code of `C`.

This augments the all 1's row and then extends.
"""
function lengthen(C::AbstractLinearCode)
    row = matrix(field(C), transpose([1 for _ in 1:length(C)]))
    newcode = augment(C, row)
    return extend(newcode)
end

"""
    uuplusv(C1::AbstractLinearCode, C2::AbstractLinearCode, verify::Bool=true)
    Plotkinconstruction(C1::AbstractLinearCode, C2::AbstractLinearCode, verify::Bool=true)

Return the Plotkin- or so-called (u | u + v)-construction with `u ∈ C1` and `v ∈ C2`.

# Arguments
* `verify`: Runs a verification step on the sizes of the constructed code versus the
theoretically guarenteed values.
"""
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

    return LinearCode(field(C1), size(G, 2), size(G, 1), d, G, missing, H, missing,
        Gstand, Hstand, missing)
end
Plotkinconstruction(C1::AbstractLinearCode, C2::AbstractLinearCode, verify::Bool=true) = uuplusv(C1, C2, verify)

"""
    subcode(C::AbstractLinearCode, k::Integer)

Return a `k`-dimensional subcode of `C`.
"""
function subcode(C::AbstractLinearCode, k::Integer)
    k >= 1 && k < dimension(C) || error("Cannot construct a $k-dimensional subcode of an $(dimension(C))-dimensional code.")
    k != dimension(C) || return C
    return LinearCode(generatormatrix(C)[1:k, :])
end

"""
    subcode(C::AbstractLinearCode, rows::Vector{Int64})

Return a subcode of `C` using the rows of the generator matrix of `C` listed in
`rows`.
"""
function subcode(C::AbstractLinearCode, rows::Vector{Int64})
    !isempty(rows) || error("Row index set empty in subcode.")
    rows ⊆ 1:dimension(C) || error("Rows are not a subset of the index set.")
    length(rows) != dimension(C) || return C
    return LinearCode(generatormatri(C)[setdiff(1:dimension(C), rows), :])
end

"""
    juxtaposition(C1::AbstractLinearCode, C2::AbstractLinearCode)

Return the code generated by the horizontal concatenation of the generator
matrices of `C1` then `C2`.
"""
function juxtaposition(C1::AbstractLinearCode, C2::AbstractLinearCode)
    field(C1) == field(C2) || error("Cannot juxtapose two codes over different fields.")
    dimension(C1) == dimension(C2) || error("Cannot juxtapose two codes of different dimensions.")
    return LinearCode(hcat(generatormatrix(C1), generatormatrix(C2)))
end

"""
    constructionX(C1::AbstractLinearCode, C2::AbstractLinearCode, C3::AbstractLinearCode)

Return the code generated by the construction X procedure.

Let `C1` be an [n, k, d], `C2` be an [n, k - l, d + e], and `C3` be an [m, l, e] linear code
with `C2 ⊂ C1` be proper. Construction X creates a [n + m, k, d + e] linear code.
"""
function constructionX(C1::AbstractLinearCode, C2::AbstractLinearCode, C3::AbstractLinearCode)
    C1 ⊆ C2 || error("The first code must be a subcode of the second in construction X.")
    !isequivalent(C1, C2) || error("The first code must be a proper subcode of the second in construction X.")
    # the above line checks field(C1) == field(C2)
    field(C1) == field(C3) || error("All codes must be over the same base ring in construction X.")
    dimension(C2) == dimension(C1) + dimension(C3) ||
        error("The dimension of the second code must be the sum of the dimensions of the first and third codes.")
    # could do some verification steps on parameters
    return LinearCode(vcat(hcat(generatormatrix(C1 / C2), generatormatrix(C3)),
        hcat(generatormatrix(C1), zero_matrix(field(C1), dimension(C1), length(C3)))))
end

"""
    constructionX3(C1::AbstractLinearCode, C2::AbstractLinearCode, C3::AbstractLinearCode,
        C4::AbstractLinearCode, C5::AbstractLinearCode))

Return the code generated by the construction X3 procedure.

Let C1 = [n, k1, d1], C2 = [n, k2, d2], C3 = [n, k3, d3], C4 = [n4, k2 - k1, d4], and
C5 = [n5, k3 - k2, d5] with `C1 ⊂ C2 ⊂ C3`. Construction X3 creates an [n + n4 + n5, k3, d]
linear code with d ≥ min{d1, d2 + d4, d3 + d5}.
"""
function constructionX3(C1::AbstractLinearCode, C2::AbstractLinearCode, C3::AbstractLinearCode,
    C4::AbstractLinearCode, C5::AbstractLinearCode)

    C1 ⊆ C2 || error("The first code must be a subcode of the second in construction X3.")
    !isequivalent(C1, C2) || error("The first code must be a proper subcode of the second in construction X3.")
    C2 ⊆ C3 || error("The second code must be a subcode of the third in construction X3.")
    !isequivalent(C2, C3) || error("The second code must be a proper subcode of the third in construction X3.")
    # the above lines check field(C1) == field(C2) == field(C3)
    field(C1) == field(C4)  == field(C5) || error("All codes must be over the same base ring in construction X3.")
    dimension(C3) == dimension(C2) + dimension(C4) ||
        error("The dimension of the third code must be the sum of the dimensions of the second and fourth codes in construction X3.")
    dimension(C2) == dimension(C1) + dimension(C5) ||
        error("The dimension of the second code must be the sum of the dimensions of the first and fifth codes in construction X3.")

    C2modC1 = C2 / C1
    C3modC2 = C3 / C2
    F = field(C1)

    # could do some verification steps on parameters
    G = vcat(hcat(generatormatrix(C1), zero_matrix(F, dimension(C1), length(C4)), zero_matrix(F, dimension(C1), length(C5))),
    hcat(generatormatrix(C2modC1), generatormatrix(C4), zero_matrix(F, dimension(C4), length(C5))),
    hcat(generatormatrix(C3modC2), zero_matrix(F, dimension(C5), length(C4)), generatormatrix(C5)))
    return LinearCode(G)
end

"""
    upluswvpluswuplusvplusw(C1::AbstractLinearCode, C2::AbstractLinearCode)

Return the code generated by the (u + w | v + w | u + v + w)-construction.

Let C1 = [n, k1, d1] and C2 = [n, k2, d2]. This construction produces an [3n, 2k1 + k2]
linear code. For binary codes, wt(u + w | v + w | u + v + w) = 2 wt(u ⊻ v) - wt(w) + 4s,
where s = |{i | u_i = v_i = 0, w_i = 1}|.
"""
function upluswvpluswuplusvplusw(C1::AbstractLinearCode, C2::AbstractLinearCode)
    field(C1) == field(C2) || error("All codes must be over the same base ring in the (u + w | v + w | u + v + w)-construction.")
    length(C1) == length(C2) || error("All codes must be the same length in the (u + w | v + w | u + v + w)-construction.")

    # could do some verification steps on parameters
    G = vcat(hcat(generatormatrix(C1), zero_matrix(field(C1), dimension(C1), length(C1)), generatormatrix(C1)),
    hcat(generatormatrix(C2), generatormatrix(C2), generatormatrix(C2)),
    hcat(zero_matrix(field(C1), dimension(C1), length(C1)), generatormatrix(C1), generatormatrix(C1)))
    return LinearCode(G)
end

"""
    expandedcode(C::AbstractLinearCode, K::FqNmodFiniteField, basis::Vector{fq_nmod})

Return the expanded code of `C` constructed by exapnding the generator matrix
to the subfield `K` using the provided `basis` for the base ring of `M` over `K`.

No check is done to ensure that `basis` is indeed a basis for the extension.
"""
function expandedcode(C::AbstractLinearCode, K::FqNmodFiniteField, basis::Vector{fq_nmod})
    return LinearCode(expandmatrix(generatormatrix(C), K, basis))
end

# R.Pellikaan, On decoding by error location and dependent sets of error
# positions, Discrete Mathematics, 106107 (1992), 369-381.
# the Schur product of vector spaces is highly basis dependent and is often the
# full vector space (an [n, n, 1] code)
"""
    entrywiseproductcode(C::AbstractLinearCode, D::AbstractLinearCode)
    *(C::AbstractLinearCode, D::AbstractLinearCode)
    Schurproductcode(C::AbstractLinearCode, D::AbstractLinearCode)
    Hadamardproductcode(C::AbstractLinearCode, D::AbstractLinearCode)
    componentwiseproductcode(C::AbstractLinearCode, D::AbstractLinearCode)

Return the entrywise product of `C` and `D`.

Note that this is known to often be the full ambient space.
"""
function entrywiseproductcode(C::AbstractLinearCode, D::AbstractLinearCode)
    field(C) == field(D) || error("Codes must be over the same field in the Schur product.")
    length(C) == length(D) || error("Codes must have the same length in the Schur product.")

    GC = generatormatrix(C)
    GD = generatormatrix(D)
    indices = Vector{Tuple{Int, Int}}()
    for i in 1:dimension(C)
        for j in 1:dimension(D)
            i <= j && push!(indices, (i, j))
        end
    end
    return LinearCode(matrix(field(C), vcat([GC[i, :] .* GD[j, :] for (i, j) in indices]...)))

    # verify C ⊂ it
end
*(C::AbstractLinearCode, D::AbstractLinearCode) = entrywiseproductcode(C, D)
Schurproductcode(C::AbstractLinearCode, D::AbstractLinearCode) = entrywiseproductcode(C, D)
Hadamardproductcode(C::AbstractLinearCode, D::AbstractLinearCode) = entrywiseproductcode(C, D)
componentwiseproductcode(C::AbstractLinearCode, D::AbstractLinearCode) = entrywiseproductcode(C, D)

"""
    VectorSpace(C::AbstractLinearCode)

Return the code `C` as a vector space object.
"""
function VectorSpace(C::AbstractLinearCode)
    V = VectorSpace(field(C), length(C))
    G = generatormatrix(C)
    return sub(V, [V(G[i, :]) for i in 1:nrows(G)])
end

# needs significant testing, works in a very special case
function evensubcode(C::AbstractLinearCode)
    F = field(C)
    # add checks here that this is only implemented for binary codes

    VC, ψ = VectorSpace(C)
    GFVS = VectorSpace(F, 1)
    homo1quad = ModuleHomomorphism(VC, GFVS, matrix(F, dim(VC), 1,
        vcat([F(weight(ψ(g).v) % 2) for g in gens(VC)]...)))
    evensub, ϕ1 = kernel(homo1quad)
    !iszero(dim(evensub)) && return LinearCode(vcat([ψ(ϕ1(g)).v for g in
        gens(evensub)]...))
    return missing
end

# needs significant testing, works in a very special case
function doublyevensubcode(C::AbstractLinearCode)
    F = field(C)
    # add checks here that this is only implemented for binary codes

    VC, ψ = VectorSpace(C)
    GFVS = VectorSpace(F, 1)

    # first get the even subspace
    homo1quad = ModuleHomomorphism(VC, GFVS, matrix(F, dim(VC), 1,
        vcat([F(weight(ψ(g).v) % 2) for g in gens(VC)]...)))
    evensub, ϕ1 = kernel(homo1quad)

    if !iszero(dim(evensub))
        # now control the overlap (Ward's divisibility theorem)
        homo2bi = ModuleHomomorphism(evensub, evensub, matrix(F, dim(evensub),
            dim(evensub),
            vcat([F(weight(matrix(F, 1, length(C), ψ(ϕ1(gens(evensub)[i])).v .*
                ψ(ϕ1(gens(evensub)[j])).v)) % 2)
            for i in 1:dim(evensub), j in 1:dim(evensub)]...)))
        evensubwoverlap, μ1 = kernel(homo2bi)

        if !iszero(dim(evensubwoverlap))
            # now apply the weight four condition
            homo2quad = ModuleHomomorphism(evensubwoverlap, GFVS, matrix(F,
                dim(evensubwoverlap), 1, vcat([F(div(weight(ψ(ϕ1(μ1(g))).v),
                2) % 2) for g in gens(evensubwoverlap)]...)))
            foursub, ϕ2 = kernel(homo2quad)

            if !iszero(dim(foursub))
                return LinearCode(vcat([ψ(ϕ1(μ1(ϕ2(g)))).v for g in gens(foursub)]...))
            end
        end
    end
    return missing
end

# currently incomplete
# function triplyevensubcode(C::AbstractLinearCode)
#     F = field(C)
#     VC, ψ = VectorSpace(C)
#     GFVS = VectorSpace(F, 1)
#
#     # first get the even subspace
#     homo1quad = ModuleHomomorphism(VC, GFVS, matrix(F, dim(VC), 1, vcat([F(weight(ψ(g).v) % 2) for g in gens(VC)]...)))
#     evensub, ϕ1 = kernel(homo1quad)
#
#     # now control the overlap (Ward's divisibility theorem)
#     homo2bi = ModuleHomomorphism(evensub, evensub, matrix(F, dim(evensub), dim(evensub),
#         vcat([F(weight(matrix(F, 1, length(C), ψ(ϕ1(gens(evensub)[i])).v .* ψ(ϕ1(gens(evensub)[j])).v)) % 2)
#         for i in 1:dim(evensub), j in 1:dim(evensub)]...)))
#     evensubwoverlap, μ1 = kernel(homo2bi)
#
#     # now apply the weight four condition
#     homo2quad = ModuleHomomorphism(evensubwoverlap, GFVS, matrix(F, dim(evensubwoverlap), 1,
#         vcat([F(div(weight(ψ(ϕ1(μ1(g))).v), 2) % 2) for g in gens(evensubwoverlap)]...)))
#     foursub, ϕ2 = kernel(homo2quad)
#
#     # now control the triple overlap
#     ###########
#     # how do we actually represent this?
#
#
#     # now apply the weight eight condition
#     homo3 = ModuleHomomorphism(foursub, GFVS, matrix(F, dim(foursub), 1, vcat([F(div(weight(ψ(ϕ1(ϕ2(g))).v), 4) % 2) for g in gens(foursub)]...)))
#     eightsub, ϕ3 = kernel(homo3)
#     # # println(dim(eightsub))
#     # # println(vcat([ψ(ϕ1(ϕ2(ϕ3(g)))).v for g in gens(eightsub)]...))
#     #
#     # q, γ = quo(foursub, eightsub)
#     # println(dim(q))
#     # # println([preimage(γ, g).v for g in gens(q)])
#     # return LinearCode(vcat([ψ(ϕ1(ϕ2(preimage(γ, g)))).v for g in gens(q)]...))
#     # # return LinearCode(vcat([ψ(ϕ1(ϕ2(ϕ3(g)))).v for g in gens(eightsub)]...))
#
#     return LinearCode(vcat([ψ(ϕ1(μ1(ϕ2(g)))).v for g in gens(foursub)]...))
# end
