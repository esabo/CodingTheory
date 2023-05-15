# Copyright (c) 2021, 2022, 2023 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    LinearCode(G::fq_nmod_mat, parity::Bool=false)

Return the linear code constructed with generator matrix `G`.

# Notes
* If `G` is not full rank, a row-reduced form is computed for the generator matrix.
  The dimension of the code is the number of rows of the full-rank matrix, and the
  length the number of columns. If the optional paramater `parity` is set to `true`,
  a linear code is built with `G` as the parity-check matrix.
* At the moment, no convention is used for G = 0 and an error is thrown.
* Zero columns are not removed.

# Examples
```julia
julia> F = GF(2);
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
function LinearCode(G::CTMatrixTypes, parity::Bool=false)
    iszero(G) && throw(ArgumentError("Zero matrix passed into LinearCode constructor."))

    Gnew = deepcopy(G)
    Gnew = _removeempty(Gnew, :rows)
    Gstand, Hstand, P, k = _standardform(Gnew)
    if ismissing(P)
        _, H = right_kernel(Gnew)
        # note the H here is transpose of the standard definition
        # remove empty for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
        H = _removeempty(transpose(H), :rows)
    else
        H = Hstand * transpose(P)
    end

    if parity
        ub1, _ = _minwtrow(H)
        ub2, _ = _minwtrow(Hstand)
        ub = minimum([ub1, ub2])
        # treat G as the parity-check matrix H
        return LinearCode(base_ring(Gnew), ncols(H), nrows(Hstand), missing, 1, ub, H, Gnew, Hstand, Gstand, transpose(P), missing)
    else
        ub1, _ = _minwtrow(Gnew)
        ub2, _ = _minwtrow(Gstand)
        ub = minimum([ub1, ub2])
        return LinearCode(base_ring(Gnew), ncols(Gnew), k, missing, 1, ub, Gnew, H, Gstand, Hstand, P, missing)
    end
end

"""
LinearCode(V::AbstractAlgebra.Generic.FreeModule{fq_nmod}, parity::Bool=false)

Return the linear code constructed from the basis of the vector space `V`. If the optional paramater
`parity` is set to `true`, a linear code is built with `G` as the parity-check matrix.
"""
LinearCode(V::AbstractAlgebra.Generic.FreeModule{fq_nmod}, parity::Bool=false) = LinearCode(basis(V), parity)

#############################
      # getter functions
#############################

"""
    field(C::AbstractLinearCode)

Return the base ring of the generator matrix.
"""
field(C::AbstractLinearCode) = C.F

"""
    length(C::AbstractLinearCode)

Return the length of the code.
"""
length(C::AbstractLinearCode) = C.n

"""
    dimension(C::AbstractLinearCode)

Return the dimension of the code.
"""
dimension(C::AbstractLinearCode) = C.k

"""
    cardinality(C::AbstractLinearCode)

Return the cardinality of the code.
"""
cardinality(C::AbstractLinearCode) = BigInt(order(C.F))^C.k

"""
    rate(C::AbstractLinearCode)

Return the rate, `R = k/n`, of the code.
"""
rate(C::AbstractLinearCode) = C.k / C.n

"""
    generatormatrix(C::AbstractLinearCode, standform::Bool=false)

Return the generator matrix of the code.

# Notes
* If the optional parameter `standform` is set to `true`, the standard form of the
generator matrix is returned instead.
"""
generatormatrix(C::AbstractLinearCode, standform::Bool=false) = standform ? (return C.Gstand;) : (return C.G;)

"""
    paritycheckmatrix(C::AbstractLinearCode, standform::Bool=false)

Return the parity-check matrix of the code.

# Notes
* If the optional parameter `standform` is set to `true`, the standard form of the
  parity-check matrix is returned instead.
"""
paritycheckmatrix(C::AbstractLinearCode, standform::Bool=false) = standform ? (return C.Hstand;) : (return C.H;)

"""
    standardformpermutation(C::AbstractLinearCode)

Return the permutation matrix required to permute the columns of the code matrices to have the same
row space as the matrices in standard form. Returns `missing` is no such permutation is required.
"""
# TODO: perhaps should just return the identity?
standardformpermutation(C::AbstractLinearCode) = C.Pstand

"""
    relativedistance(C::AbstractLinearCode)

Return the relative minimum distance, `δ = d / n` of the code if `d` is known,
otherwise return `missing`.
"""
relativedistance(C::AbstractLinearCode) = ismissing(C.d) ? (return missing) : (return C.d // C.n)

"""
    genus(C::AbstractLinearCode)

Return the genus, `n + 1 - k - d`, of the code.
"""
genus(C::AbstractLinearCode) = ismissing(C.d) ? (return missing) : (return C.n + 1 - C.k - minimumdistance(C))

"""
    minimumdistancelowerbound(C::AbstractLinearCode)

Return the current lower bound on the minimum distance of `C`.
"""
minimumdistancelowerbound(C::AbstractLinearCode) = C.lbound

"""
    minimumdistanceupperbound(C::AbstractLinearCode)

Return the current upper bound on the minimum distance of `C`.
"""
minimumdistanceupperbound(C::AbstractLinearCode) = C.ubound

"""
    isMDS(C::AbstractLinearCode)

Return `true` if code is maximum distance separable (MDS).
"""
isMDS(C::AbstractLinearCode) = ismissing(C.d) ? (return missing) :
    (minimumdistance(C) != Singletonbound(C.n, C.k) ? (return true) : (return false))

"""
    numbercorrectableerrors(C::AbstractLinearCode)

Return the number of correctable errors for the code.

# Notes
* The number of correctable errors is `t = floor((d - 1) / 2)`.
"""
numbercorrectableerrors(C::AbstractLinearCode) = ismissing(C.d) ? (return missing) : (return Int(floor((minimumdistance(C) - 1) / 2)))

#############################
      # setter functions
#############################

"""
    setdistancelowerbound!(C::AbstractLinearCode, l::Int)

Set the lower bound on the minimum distance of `C`, if `l` is better than the current bound.
"""
function setdistancelowerbound!(C::AbstractLinearCode, l::Int)
    1 <= l <= C.ubound || throw(DomainError("The lower bound must be between 1 and the upper bound."))
    C.lbound < l && (C.lbound = l;)
    if C.lbound == C.ubound
        @warn "The new lower bound is equal to the upper bound; setting the minimum distance."
        C.d = C.lbound
    end
end

"""
    setdistanceupperbound!(C::AbstractLinearCode, u::Int)

Set the upper bound on the minimum distance of `C`, if `u` is better than the current bound.
"""
function setdistanceupperbound!(C::AbstractLinearCode, u::Int)
    C.lbound <= u <= C.n || throw(DomainError("The upper bound must be between the lower bound and the code length."))
    u < C.ubound && (C.ubound = u;)
    if C.lbound == C.ubound
        @warn "The new upper bound is equal to the lower bound; setting the minimum distance."
        C.d = C.lbound
    end
end

"""
    setminimumdistance(C::AbstractLinearCode, d::Int)

Set the minimum distance of the code to `d`.

# Notes
* The only check done on the value of `d` is that `1 ≤ d ≤ n`.
"""
function setminimumdistance!(C::AbstractLinearCode, d::Int)
    d > 0 && d <= C.n || throw(DomainError("The minimum distance of a code must be ≥ 1; received: d = $d."))
    C.d = d
end

#############################
     # general functions
#############################

function _standardform(G::CTMatrixTypes)
    rnk, Gstand, P = _rref_col_swap(G, 1:nrows(G), 1:ncols(G))
    nrows(Gstand) > rnk && (Gstand = _removeempty(Gstand, :rows);)
    A = Gstand[:, (nrows(Gstand) + 1):ncols(Gstand)]
    T = MatrixSpace(base_ring(Gstand), ncols(Gstand) - nrows(Gstand),
        ncols(Gstand) - nrows(Gstand))
    Hstand = hcat(-transpose(A), T(1))
    return Gstand, Hstand, P, rnk
end

# a bit odd to handle all subtypes in the supertype but saves hundreds
# of repeated lines and keeps uniform
function show(io::IO, C::AbstractLinearCode)
    print(io, "[$(C.n), $(C.k)")
    !ismissing(C.d) && print(io, ", $(C.d)")
    typeof(C) <: AbstractBCHCode && print(io, "; $(C.b)")
    print(io, "]_$(order(C.F)) ")
    if isa(C, ReedSolomonCode)
        println(io, "Reed-Solomon code")
    elseif isa(C, BCHCode)
        println(io, "BCH code")
    elseif isa(C, CyclicCode)
        println(io, "cyclic code")
    elseif isa(C, ReedMullerCode)
        println(io, "Reed-Muller code RM($(C.r), $(C.m))")
    elseif isa(C, QuasiCyclicCode)
        println(io, "quasi-cyclic code of index $(C.l)")
    elseif isa(C, GeneralizedReedSolomonCode)
        println(io, "generalized Reed-Solomon code")
    else
        println(io, "linear code")
    end

    if get(io, :compact, true)
        if typeof(C) <: AbstractCyclicCode
            println(io, "$(order(C.F))-Cyclotomic cosets: ")
            len = length(qcosetsreps(C))
            if len == 1
                println("\tC_$(qcosetsreps(C)[1])")
            else
                for (i, x) in enumerate(qcosetsreps(C))
                    if i == 1
                        print(io, "\tC_$x ∪ ")
                    elseif i == 1 && i == len
                        println(io, "\tC_$x")
                    elseif i != len
                        print(io, "C_$x ∪ ")
                    else
                        println(io, "C_$x")
                    end
                end
            end
            println(io, "Generator polynomial:")
            println(io, "\t", generatorpolynomial(C))
        end

        if C.n <= 30
            if isa(C, QuasiCyclicCode)
                if C.Atype == 'G'
                    M = generatormatrix(C)
                    nr, nc = size(M)
                    println(io, "Generator matrix: $(nr) × $(nc)")
                else
                    M = paritycheckmatrix(C)
                    nr, nc = size(M)
                    println(io, "Parity-check matrix: $(nr) × $(nc)")
                end
                for i in 1:nr
                    print(io, "\t")
                    for j in 1:nc
                        if j != nc
                            print(io, "$(M[i, j]) ")
                        elseif j == nc && i != nr
                            println(io, "$(M[i, j])")
                        else
                            print(io, "$(M[i, j])")
                        end
                    end
                end
            else
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
            end
        end
        # if !ismissing(C.weightenum)
        #     println(io, "\nComplete weight enumerator:")
        #     print(io, "\t", polynomial(C.weightenum))
        # end
    end
end

# TODO: symbol :k or :d
function Singletonbound(n::Int, a::Int)
    # d ≤ n - k + 1 or k ≤ n - d + 1
    if n >= 0 && a >= 0 && n >= a
        return n - a + 1
    else
        error("Invalid parameters for the Singleton bound. Received n = $n, k/d = $a")
    end
end
Singletonbound(C::AbstractLinearCode) = Singletonbound(C.n, C.k)

"""
    encode(v::Union{fq_nmod_mat, Vector{Int}}, C::AbstractLinearCode)

Return `v * G`, where `G` is the generator matrix of `C`.
"""
# TODO: check quantum functions and make uniform - prefer C then v on all such functions
# might be a breaking fix though, also check runtests.jl
function encode(v::CTMatrixTypes, C::AbstractLinearCode)
    G = generatormatrix(C)
    nr = nrows(G)
    (size(v) != (1, nr) && size(v) != (nr, 1)) &&
        throw(ArgumentError("Vector has incorrect dimension; expected length $nr, received: $(size(v))."))
    base_ring(v) == C.F || throw(ArgumentError("Vector must have the same base ring as the generator matrix."))
    nrows(v) != 1 || return v * G
    return transpose(v) * G
end
# TODO: combine these two functions
function encode(v::Vector{Int}, C::AbstractLinearCode)
    length(v) == C.k ||
        throw(ArgumentError("Vector has incorrect length; expected length $(C.k), received: $(size(v))."))
    return encode(matrix(C.F, transpose(v)), C)
end

"""
    syndrome(C::AbstractLinearCode, v::Union{fq_nmod_mat, Vector{Int}})

Return `Hv`, where `H` is the parity-check matrix of `C`.
"""
function syndrome(C::AbstractLinearCode, v::CTMatrixTypes)
    H = paritycheckmatrix(C)
    nc = ncols(H)
    (size(v) != (nc, 1) && size(v) != (1, nc)) &&
        throw(ArgumentError("Vector has incorrect dimension; expected length $nc, received: $(size(v))."))
    if base_ring(v) != C.F
        if order(base_ring(v)) == order(C.F)
            @warn "Fields are of different types, but have the same order."
        else
            throw(ArgumentError("Vector must have the same base ring as the parity-check matrix."))
        end
    end
    nrows(v) != 1 || return H * transpose(v)
    return H * v
end
# TODO: combine these two functions
function syndrome(C::AbstractLinearCode, v::Vector{Int})
    length(v) == C.n ||
        throw(ArgumentError(("Vector to be tested is of incorrect dimension; expected length $(C.n), received: $(size(v)).")))
    return syndrome(C, matrix(C.F, transpose(v)))
end

"""
    in(v::Union{fq_nmod_mat, Vector{Int}}, C::AbstractLinearCode)

Return whether or not `v` is a codeword of `C`.
"""
in(v::CTMatrixTypes, C::AbstractLinearCode) = iszero(syndrome(C, v))
in(v::Vector{Int}, C::AbstractLinearCode) = iszero(syndrome(C, v))

"""
    ⊆(C1::AbstractLinearCode, C2::AbstractLinearCode)
    ⊂(C1::AbstractLinearCode, C2::AbstractLinearCode)
    issubcode(C1::AbstractLinearCode, C2::AbstractLinearCode)

Return whether or not `C1` is a subcode of `C2`.
"""
function ⊆(C1::AbstractLinearCode, C2::AbstractLinearCode)
    if C1.F != C2.F || C1.n != C2.n || C1.k > C2.k
        if order(C1.F) == order(C2.F)
            @warn "Fields are of different types, but have the same order."
        else
            return false
        end
    end

    # eachrow doesn't work on these objects
    G1 = generatormatrix(C1)
    nr = nrows(G1)
    for r in 1:nr
        if G1[r, :] ∉ C2
            return false
        end
    end
    return true
end
⊂(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊆ C2
issubcode(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊆ C2

"""
    dual(C::AbstractLinearCode)
    Euclideandual(C::AbstractLinearCode)

Return the (Euclidean) dual of the code `C`.
"""
function dual(C::AbstractLinearCode)
    G = generatormatrix(C)
    Gstand = generatormatrix(C, true)
    H = paritycheckmatrix(C)
    Hstand = paritycheckmatrix(C, true)
    ub1, _ = _minwtrow(H)
    ub2, _ = _minwtrow(Hstand)
    ub = minimum([ub1, ub2])

    if !ismissing(C.weightenum)
        dualwtenum = MacWilliamsIdentity(C, C.weightenum)
        HWE = CWEtoHWE(dualwtenum)
        d = minimum([collect(exponent_vectors(polynomial(HWE)))[i][1]
            for i in 1:length(polynomial(HWE))])
        return LinearCode(C.F, C.n, C.n - C.k, d, 1, ub, deepcopy(H), deepcopy(G),
            deepcopy(Hstand), deepcopy(Gstand), deepcopy(C.P), dualwtenum)
    else
        return LinearCode(C.F, C.n, C.n - C.k, missing, 1, ub, deepcopy(H),
            deepcopy(G), deepcopy(Hstand), deepcopy(Gstand), deepcopy(C.P), missing)
    end
end
Euclideandual(C::AbstractLinearCode) = dual(C)

"""
    Hermitiandual(C::AbstractLinearCode)

Return the Hermitian dual of a code defined over a quadratic extension.
"""
Hermitiandual(C::AbstractLinearCode) = dual(LinearCode(Hermitianconjugatematrix(C.G)))
# also valid to do LinearCode(Hermitianconjugatematrix(generatormatrix(dual(C))))

"""
    areequivalent(C1::AbstractLinearCode, C2::AbstractLinearCode)

Return `true` if `C1 ⊆ C2` and `C2 ⊆ C1`.
"""
areequivalent(C1::AbstractLinearCode, C2::AbstractLinearCode) = (C1 ⊆ C2 && C2 ⊆ C1)

"""
    isselfdual(C::AbstractLinearCode)

Return `true` if `areequivalent(C, dual(C))`.
"""
isselfdual(C::AbstractLinearCode) = areequivalent(C, dual(C))

"""
    isselforthogonal(C::AbstractLinearCode)
    isweaklyselfdual(C::AbstractLinearCode)
    isEuclideanselforthogonal(C::AbstractLinearCode)

Return `true` if `C ⊆ dual(C)`.
"""
isselforthogonal(C::AbstractLinearCode) = C ⊆ dual(C)
isweaklyselfdual(C::AbstractLinearCode) = isselforthogonal(C)
isEuclideanselforthogonal(C::AbstractLinearCode) = isselforthogonal(C)

"""
    isdualcontaining(C::AbstractLinearCode)
    isEuclideandualcontaining(C::AbstractLinearCode)

Return `true` if `dual(C) ⊆ C`.
"""
isdualcontaining(C::AbstractLinearCode) = dual(C) ⊆ C
isEuclideandualcontaining(C::AbstractLinearCode) = isdualcontaining(C)

"""
    isHermitianselfdual(C::AbstractLinearCode)

Return `true` if `areequivalent(C, Hermitiandual(C))`.
"""
isHermitianselfdual(C::AbstractLinearCode) = areequivalent(C, Hermitiandual(C))

"""
    isHermitianselforthogonal(C::AbstractLinearCode)
    isHermitianweaklyselfdual(C::AbstractLinearCode)

Return `true` if `C ⊆ Hermitiandual(C)`.
"""
isHermitianselforthogonal(C::AbstractLinearCode) = C ⊆ Hermitiandual(C)
isHermitianweaklyselfdual(C::AbstractLinearCode) = isHermitianselforthogonal(C)

"""
    isHermitiandualcontaining(C::AbstractLinearCode)

Return `true` if `Hermitiandual(C) ⊆ C`.
"""
isHermitiandualcontaining(C::AbstractLinearCode) = Hermitiandual(C) ⊆ C

# TODO: add l-Galois dual and self-dual/orthogonal based functions

"""
    characteristicpolynomial(C::AbstractLinearCode)

Return the characteristic polynomial of `C`.

# Notes
* The characteristic polynomial is defined in [Lin1999]_
"""
function characteristicpolynomial(C::AbstractLinearCode)
    _, x = PolynomialRing(Nemo.QQ, "x")
    D = dual(C)
    supD = support(D)
    q = Int(order(C.F))
    return q^(n - k) * prod([1 - x / j for j in supD if j > 0])
end

"""
    VectorSpace(C::AbstractLinearCode)

Return the code `C` as a vector space object.
"""
function VectorSpace(C::AbstractLinearCode)
    V = VectorSpace(C.F, C.n)
    G = generatormatrix(C)
    return sub(V, [V(G[i, :]) for i in 1:nrows(G)])
end

"""
    iseven(C::AbstractLinearCode)

Return `true` if `C` is even.
"""
function iseven(C::AbstractLinearCode)
    Int(order(C.F)) == 2 || throw(ArgumentError("Even-ness is only defined for binary codes."))
    
    # A binary code generated by G is even if and only if each row of G has
    # even weight.
    G = generatormatrix(C)
    for r in 1:nrows(G)
        wt(G[r, :]) % 2 == 0 || (return false;)
    end
    return true
end

"""
    isdoublyeven(C::AbstractLinearCode)

Return `true` if `C` is doubly-even.
"""
function isdoublyeven(C::AbstractLinearCode)
    Int(order(C.F)) == 2 || throw(ArgumentError("Even-ness is only defined for binary codes."))

    # A binary code generated by G is doubly-even if and only if each row of G
    # has weight divisible by 4 and the sum of any two rows of G has weight
    # divisible by 4.
    G = generatormatrix(C)
    nr = nrows(G)
    for r in 1:nr
        # TODO: allocates a lot, should redo calculation here
        wt(G[r, :]) % 4 == 0 || (return false;)
    end
    for r1 in 1:nr
        for r2 in 1:nr
            # or by Ward's thm can do * is % 2 == 0
            wt(G[r1, :] + G[r2, :]) % 4 == 0 || (return false;)
        end
    end
    return true
end

"""
    istriplyeven(C::AbstractLinearCode)

Return `true` if `C` is triply-even.
"""
function istriplyeven(C::AbstractLinearCode)
    Int(order(C.F)) == 2 || throw(ArgumentError("Even-ness is only defined for binary codes."))

    # following Ward's divisibility theorem
    G = FpmattoJulia(generatormatrix(C))
    nr, _ = size(G)
    for r in 1:nr
        wt(G[r, :]) % 8 == 0 || (return false;)
    end
    for r1 in 1:nr
        for r2 in 1:nr
            wt(G[r1, :] .* G[r2, :]) % 4 == 0 || (return false;)
        end
    end
    for r1 in 1:nr
        for r2 in 1:nr
            for r3 in 1:nr
                wt(G[r1, :] .* G[r2, :] .* G[r3, :])% 2 == 0 || (return false;)
            end
        end
    end
    return true
end

"""
    permutecode(C::AbstractLinearCode, σ::Union{Perm{T}, Vector{T}}) where T <: Int
    permutecode!(C::AbstractLinearCode, σ::Union{Perm{T}, Vector{T}}) where T <: Int

Return the code whose generator matrix is `C`'s with the columns permuted by `σ`.

# Notes
* If `σ` is a vector, it is interpreted as the desired column order for the generator matrix of `C`.
"""
# TODO: write one for QCC codes
function permutecode!(C::AbstractLinearCode, σ::Union{Perm{T}, Vector{T}}) where T <: Int
    G = generatormatrix(C)
    if typeof(σ) <: Perm
        # a straight-forward multiplication messes up Gstand, Hstand
        return LinearCode(G * matrix(C.F, Array(matrix_repr(σ))))
    else
        length(unique(σ)) == C.n || error("Incorrect number of digits in permutation.")
        (1 == minimum(σ) && C.n == maximum(σ)) || error("Digits are not in the range `1:$(C.n)`.")

        return LinearCode(G[:, σ])
    end
end
permutecode(C::AbstractLinearCode, σ::Union{Perm{T}, Vector{T}}) where T <: Int = (Cnew = deepcopy(C); return permutecode!(Cnew);)

"""
    words(C::AbstractLinearCode, onlyprint::Bool=false)
    codewords(C::AbstractLinearCode, onlyprint::Bool=false)
    elements(C::AbstractLinearCode, onlyprint::Bool=false)

Return the elements of `C`.

# Notes
* If `onlyprint` is `true`, the elements are only printed to the console and not
  returned.
"""
function words(C::AbstractLinearCode, onlyprint::Bool=false)
    words = Vector{fq_nmod_mat}()
    G = generatormatrix(C)
    E = base_ring(G)
    # Nemo.AbstractAlgebra.ProductIterator
    for iter in Base.Iterators.product([0:(Int(characteristic(E)) - 1) for _ in 1:nrows(G)]...)
        row = E(iter[1]) * G[1, :]
        for r in 2:nrows(G)
            if !iszero(iter[r])
                row += E(iter[r]) * G[r, :]
            end
        end
        if !onlyprint
            push!(words, row)
        else
            println(row)
        end
    end
    if !onlyprint
        return words
    else
        return
    end
end
codewords(C::AbstractLinearCode, onlyprint::Bool=false) = words(C, onlyprint)
elements(C::AbstractLinearCode, onlyprint::Bool=false) = words(C, onlyprint)

"""
    hull(C::AbstractLinearCode)
    Euclideanhull(C::AbstractLinearCode)

Return the (Euclidean) hull of `C` and its dimension.

# Notes
* The hull of a code is the intersection of it and its dual.
"""
function hull(C::AbstractLinearCode)
    G = generatormatrix(C)
    H = paritycheckmatrix(C)
    VS = VectorSpace(F, C.n)
    U, _ = sub(VS, [VS(G[i, :]) for i in 1:nrows(G)])
    W, WtoVS = sub(VS, [VS(H[i, :]) for i in 1:nrows(H)])
    I, _ = intersect(U, W)
    if !iszero(AbstractAlgebra.dim(I))
        Ibasis = [WtoVS(g) for g in gens(I)]
        GI = reduce(vcat, Ibasis)
        return LinearCode(GI), AbstractAlgebra.dim(I)
    else
        return missing, 0 # is this the behavior I want?
    end
end
Euclideanhull(C::AbstractLinearCode) = hull(C)

"""
    Hermitianhull::AbstractLinearCode)
    
Return the Hermitian hull of `C` and its dimension.

# Notes
* The Hermitian hull of a code is the intersection of it and its Hermitian dual.
"""
function Hermitianhull(C::AbstractLinearCode)
    D = Hermitiandual(C)
    G = generatormatrix(C)
    H = generatormatrix(D)
    VS = VectorSpace(F, C.n)
    U, _ = sub(VS, [VS(G[i, :]) for i in 1:nrows(G)])
    W, WtoVS = sub(VS, [VS(H[i, :]) for i in 1:nrows(H)])
    I, _ = intersect(U, W)
    if !iszero(AbstractAlgebra.dim(I))
        Ibasis = [WtoVS(g) for g in gens(I)]
        GI = reduce(vcat, Ibasis)
        return LinearCode(GI), AbstractAlgebra.dim(I)
    else
        return missing, 0 # is this the behavior I want?
    end
end

"""
    isLCD(C::AbstractLinearCode)

Return `true` if `C` is linear complementary dual.

# Notes
* A code is linear complementary dual if the dimension of `hull(C)` is zero.
"""
function isLCD(C::AbstractLinearCode)
    # return true if dim(hull) = 0
    _, dimhullC = hull(C)
    dimhullC == 0 ? (return true;) : (return false;)
end

"""
    isHermitianLCD(C::AbstractLinearCode)

Return `true` if `C` is linear complementary Hermitian dual.

# Notes
* A code is linear complementary Hermitian dual if the dimension of `Hermitianhull(C)` is zero.
"""
function isHermitianLCD(C::AbstractLinearCode)
    # return true if dim(hull) = 0
    _, dimhullC = Hermitianhull(C)
    dimhullC == 0 ? (return true;) : (return false;)
end

# TODO: add l-Galois hull and isLCD functions for that dual
