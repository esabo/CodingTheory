# Copyright (c) 2021, 2022, 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    LinearCode(G::CTMatrixTypes, parity::Bool=false, bruteforceWE::Bool=true)

Return the linear code constructed with generator matrix `G`. If the optional paramater `parity` is
set to `true`, a linear code is built with `G` as the parity-check matrix. If the optional parameter
`bruteforceWE` is `true`, the weight enumerator and (and therefore the distance) is calculated when
there are fewer than 1.5e5 codewords.
"""
function LinearCode(G::CTMatrixTypes, parity::Bool=false, bruteforceWE::Bool=true)
    iszero(G) && return parity ? IdentityCode(base_ring(G), ncols(G)) : ZeroCode(base_ring(G), ncols(G))

    Gnew = deepcopy(G)
    Gnew = _removeempty(Gnew, :rows)

    C = if parity
        rnkH, H = right_kernel(Gnew)
        if ncols(H) == rnkH
            Htr = transpose(H)
        else
            # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
            nr = nrows(H)
            Htr = zero_matrix(base_ring(H), rnkH, nr)
            for r in 1:nr
                for c in 1:rnkH
                    !iszero(H[r, c]) && (Htr[c, r] = H[r, c];)
                end
            end
        end
        H = Htr
        Hstand, Gstand, P, k = _standardform(H)
        k == ncols(H) && return IdentityCode(base_ring(Gnew), ncols(H))
        ub1, _ = _minwtrow(H)
        ub2, _ = _minwtrow(Hstand)
        ub = min(ub1, ub2)
        # treat G as the parity-check matrix H
        LinearCode(base_ring(Gnew), ncols(H), nrows(Hstand), missing, 1, ub, H, Gnew, Hstand, Gstand, P, missing)
    else
        Gstand, Hstand, P, k = _standardform(Gnew)
        if k == ncols(G)
            C = IdentityCode(base_ring(Gnew), ncols(G))
            C.G = Gnew
            return C
        end
        H = ismissing(P) ? Hstand : Hstand * P
        ub1, _ = _minwtrow(Gnew)
        ub2, _ = _minwtrow(Gstand)
        ub = min(ub1, ub2)
        LinearCode(base_ring(Gnew), ncols(Gnew), k, missing, 1, ub, Gnew, H, Gstand, Hstand, P, missing)
    end

    if bruteforceWE && BigInt(order(base_ring(G)))^min(k, ncols(G) - k) <= 1.5e5
        C.weightenum = if 2k <= ncols(G)
            _weightenumeratorBF(C.Gstand)
        else
            MacWilliamsIdentity(dual(C), _weightenumeratorBF(C.Hstand))
        end
        d = minimum(filter(ispositive, first.(exponent_vectors(CWEtoHWE(C.weightenum).polynomial))))
        setminimumdistance!(C, d)
    end

    return C
end

# TODO: add doc strings
function LinearCode(G::T, H::T, bruteforceWE::Bool=true) where T <: CTMatrixTypes
    ncols(G) == ncols(H) ||
        throw(ArgumentError("The number of columns of G and H should be the same (received ncols(G) = $(ncols(G)), ncols(H) = $(ncols(H)))"))
    base_ring(G) == base_ring(H) || throw(ArgumentError("G and H are not over the same field"))
    Gnew = _removeempty(G, :rows)
    Hnew = _removeempty(H, :rows)
    iszero(Gnew * transpose(Hnew)) || throw(ArgumentError("H isn't orthogonal to G"))
    Gstand, Hstand, P, k = _standardform(Gnew)
    rank(H) == ncols(G) - k || throw(ArgumentError("The given matrix H is not a parity check matrix for G"))

    ub1, _ = _minwtrow(Gnew)
    ub2, _ = _minwtrow(Gstand)
    ub = min(ub1, ub2)
    C = LinearCode(base_ring(Gnew), ncols(Gnew), k, missing, 1, ub, Gnew, Hnew, Gstand, Hstand, P, missing)

    if bruteforceWE && BigInt(order(base_ring(G)))^min(k, ncols(G) - k) <= 1.5e5
        C.weightenum = if 2k <= ncols(G)
            _weightenumeratorBF(C.Gstand)
        else
            MacWilliamsIdentity(dual(C), _weightenumeratorBF(C.Hstand))
        end
        d = minimum(filter(ispositive, first.(exponent_vectors(CWEtoHWE(C.weightenum).polynomial))))
        setminimumdistance!(C, d)
    end

    return C
end

function LinearCode(G::Matrix{Int}, q::Int, parity::Bool=false)
    factors = Nemo.factor(q)
    (length(factors) == 1 && q > 1) || throw(ArgumentError("There is no finite field of order $q."))

    p, m = first(factors)
    F = m == 1 ? GF(p) : GF(p, m, :ω)
    G1 = matrix(F, G)
    rref!(G1)
    return LinearCode(G1, parity)
end

function LinearCode(Gs::Vector{<:CTMatrixTypes})
    s = size(Gs[1])
    all(s == size(Gs[i]) for i = 2:length(Gs)) || throw(ArgumentError("Not all vectors in `Gs` were the same size."))

    G = reduce(vcat, Gs)
    rref!(G)
    return LinearCode(G)
end

function LinearCode(Gs::Vector{Vector{Int}}, q::Int, parity::Bool=false)
    s = size(Gs[1])
    all(s == size(Gs[i]) for i = 2:length(Gs)) || throw(ArgumentError("Not all vectors in `Gs` were the same size."))

    return LinearCode(reduce(vcat, Gs), q, parity)
end

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

Return the rate, ``R = k/n``, of the code.
"""
rate(C::AbstractLinearCode) = C.k / C.n

"""
    generatormatrix(C::AbstractLinearCode, standform::Bool=false)

Return the generator matrix of the code. If the optional parameter `standform`
is set to `true`, the standard form of the generator matrix is returned instead.
"""
function generatormatrix(C::AbstractLinearCode, standform::Bool=false)
    if isa(C, QuasiCyclicCode)
        standform && !ismissing(C.Gstand) && (return C.Gstand;)
        if ismissing(C.G)
            if C.Atype == :G
                G = lift(C.A)
            else
                _, G = right_kernel(lift(C.A))
                G = transpose(G)
            end
            C.G = G
        end
        if standform
            C.Gstand, C.Hstand, C.Pstand, _ = _standardform(C.G)
            return C.Gstand
        end
        return C.G
    else
        standform ? (return C.Gstand;) : (return C.G;)
    end
end

"""
    paritycheckmatrix(C::AbstractLinearCode, standform::Bool=false)

Return the parity-check matrix of the code. If the optional parameter
`standform` is set to `true`, the standard form of the parity-check matrix
is returned instead.
"""
function paritycheckmatrix(C::AbstractLinearCode, standform::Bool=false)
    if isa(C, QuasiCyclicCode)
        if standform
            ismissing(C.Hstand) || (return C.Hstand;)
            if ismissing(C.G)
                if C.Atype == :G
                    G = lift(C.A)
                else
                    _, G = right_kernel(lift(C.A))
                    G = transpose(G)
                end
                C.G = G
            end
            C.Gstand, C.Hstand, C.Pstand, _ = _standardform(C.G)
            return C.Hstand
        elseif ismissing(C.H)
            C.H = C.Atype == :H ? lift(C.A) : transpose(right_kernel(lift(C.A))[2])
        end
        return C.H
    elseif isa(C, LDPCCode)
        standform ? (return C.C.Hstand;) : (return C.C.H;)
    else
        standform ? (return C.Hstand;) : (return C.H;)
    end
end

"""
    standardformpermutation(C::AbstractLinearCode)

Return the permutation matrix required to permute the columns of the code matrices to have the same
row space as the matrices in standard form. Returns `missing` is no such permutation is required.
"""
standardformpermutation(C::AbstractLinearCode) = C.Pstand

"""
    relativedistance(C::AbstractLinearCode)

Return the relative minimum distance, ``\\delta = d / n`` of the code if ``d`` is known;
otherwise return `missing`.
"""
relativedistance(C::AbstractLinearCode) = C.d / C.n

"""
    genus(C::AbstractLinearCode)

Return the genus, ``n + 1 - k - d``, of the code.
"""
genus(C::AbstractLinearCode) = C.n + 1 - C.k - C.d

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
isMDS(C::AbstractLinearCode) = ismissing(C.d) ? missing : C.d != Singletonbound(C.n, C.k)

"""
    numbercorrectableerrors(C::AbstractLinearCode)

Return the number of correctable errors for the code.

# Notes
* The number of correctable errors is ``t = \\floor{(d - 1) / 2}``.
"""
numbercorrectableerrors(C::AbstractLinearCode) = ismissing(C.d) ? missing : Int(fld(C.d - 1, 2))

"""
    isovercomplete(C::AbstractLinearCode, which::Symbol=:G)

Return `true` if the generator matrix is over complete, or if the optional parameter is
set to :H and the parity-check matrix is over complete.
"""
function isovercomplete(C::AbstractLinearCode, which::Symbol=:G)
    if which == :G
        return nrows(C.G) > C.k
    elseif which == :H
        return nrows(C.H) > C.n - C.k
    end
    throw(ArgumentError("Received symbol $which, expected :G or :H"))
end

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
        @info "The new lower bound is equal to the upper bound; setting the minimum distance."
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
        @info "The new upper bound is equal to the lower bound; setting the minimum distance."
        C.d = C.lbound
    end
end

"""
    setminimumdistance(C::AbstractLinearCode, d::Int)

Set the minimum distance of the code to `d`.

# Notes
* The only check done on the value of `d` is that ``1 \\leq d \\leq n``.
"""
function setminimumdistance!(C::AbstractLinearCode, d::Int)
    0 < d <= C.n || throw(DomainError("The minimum distance of a code must be ≥ 1; received: d = $d."))
    C.d = d
    C.lbound = d
    C.ubound = d
end

function changefield!(C::T, F::CTFieldTypes) where T <: AbstractLinearCode
    # TODO: this is doesn't work for cyclic codes yet
    T <: AbstractCyclicCode && @error "Not implemented for cyclic codes yet"
    C.G = change_base_ring(F, C.G)
    C.H = change_base_ring(F, C.H)
    C.Gstand = change_base_ring(F, C.Gstand)
    C.Hstand = change_base_ring(F, C.Hstand)
    ismissing(C.Pstand) || (C.Pstand = change_base_ring(F, C.Pstand);)

    if order(F) != order(C.F)
        # TODO: should be able to compute the new weightenum easily as well
        C.weightenum = missing
    end

    C.F = F
    
    return nothing
end

function changefield(C::AbstractLinearCode, F::CTFieldTypes)
    C2 = deepcopy(C)
    changefield!(C2, F)
    return C2
end

#############################
     # general functions
#############################

function _standardform(G::CTMatrixTypes)
    rnk, Gstand, P = _rref_col_swap(G, 1:nrows(G), 1:ncols(G))
    F = base_ring(Gstand)
    nrows(Gstand) > rnk && (Gstand = _removeempty(Gstand, :rows);)
    Atr = transpose(view(Gstand, :, (nrows(Gstand) + 1):ncols(Gstand)))
    Hstand = hcat(order(F) == 2 ? Atr : -Atr, identity_matrix(F, ncols(Gstand) - nrows(Gstand)))
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
                if C.Atype == :G
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

"""
    Singletonbound(n::Int, a::Int)

Return the Singleton bound ``d \\leq n - k + 1`` or ``k \\leq n - d + 1`` depending on the interpretation of `a`.
"""
Singletonbound(n::Int, a::Int) = 0 <= a <= n ? (return n - a + 1) : 
    error("Invalid parameters for the Singleton bound. Received n = $n, k/d = $a")

"""
    Singletonbound(C::AbstractLinearCode)

Return the Singleton bound on the minimum distance of the code (``d \\leq n - k + 1``).
"""
Singletonbound(C::AbstractLinearCode) = Singletonbound(C.n, C.k)

"""
    encode(C::AbstractLinearCode, v::Union{CTMatrixTypes, Vector{Int}})

Return the encoding of `v` into `C`
"""
function encode(C::AbstractLinearCode, v::Union{CTMatrixTypes, Vector{Int}})
    w = isa(v, Vector{Int}) ? matrix(C.F, 1, length(v), v) : v
    G = generatormatrix(C)
    nr = nrows(G)
    (size(w) != (1, nr) && size(w) != (nr, 1)) &&
        throw(ArgumentError("Vector has incorrect dimension; expected length $nr, received: $(size(v))."))
    base_ring(w) == C.F || throw(ArgumentError("Vector must have the same base ring as the generator matrix."))
    nrows(w) != 1 || return w * G
    return transpose(w) * G
end

"""
    syndrome(C::AbstractLinearCode, v::Union{CTMatrixTypes, Vector{Int}})

Return the syndrome of `v` with respect to `C`.
"""
function syndrome(C::AbstractLinearCode, v::Union{CTMatrixTypes, Vector{Int}})
    w = isa(v, Vector{Int}) ? matrix(C.F, length(v), 1, v) : v
    H = paritycheckmatrix(C)
    nc = ncols(H)
    (size(w) != (nc, 1) && size(w) != (1, nc)) &&
        throw(ArgumentError("Vector has incorrect dimension; expected length $nc, received: $(size(v))."))
    if base_ring(w) != C.F
        if order(base_ring(w)) == order(C.F)
            @warn "Fields are of different types, but have the same order."
        else
            throw(ArgumentError("Vector must have the same base ring as the parity-check matrix."))
        end
    end
    return nrows(w) == 1 ? H * transpose(w) : H * w
end

"""
    in(v::Union{CTMatrixTypes, Vector{Int}}, C::AbstractLinearCode)

Return whether or not `v` is a codeword of `C`.
"""
in(v::Union{CTMatrixTypes, Vector{Int}}, C::AbstractLinearCode) = iszero(syndrome(C, v))

"""
    ⊆(C1::AbstractLinearCode, C2::AbstractLinearCode)
    ⊂(C1::AbstractLinearCode, C2::AbstractLinearCode)
    issubcode(C1::AbstractLinearCode, C2::AbstractLinearCode)

Return whether or not `C1` is a subcode of `C2`.
"""
function ⊆(C1::AbstractLinearCode, C2::AbstractLinearCode)
    C1.F == C2.F || (order(C1.F) == order(C2.F) ? (@warn "Fields are of different types, but have the same order.") : (return false;))
    (C1.n == C2.n && C1.k <= C2.k) || return false

    G1 = generatormatrix(C1)
    return all(view(G1, r:r, 1:C1.n) ∈ C2 for r in axes(G1, 1))
end
⊂(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊆ C2
issubcode(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊆ C2

"""
    dual(C::AbstractLinearCode)
    Euclideandual(C::AbstractLinearCode)

Return the (Euclidean) dual of the code `C`.
"""
function dual(C::AbstractLinearCode)
    if typeof(C) <: AbstractCyclicCode
        return CyclicCode(Int(order(C.F)), C.n, dualqcosets(Int(order(C.F)), C.n, C.qcosets))
    elseif isa(C, GeneralizedReedSolomonCode)
        d = C.k + 1
        return GeneralizedReedSolomonCode(C.F, C.n, C.n - C.k, d, d, d,
            deepcopy(C.dualscalars), deepcopy(C.scalars), deepcopy(C.evaluationpoints),
            deepcopy(C.H), deepcopy(C.G), deepcopy(C.Hstand),
            deepcopy(C.Gstand), deepcopy(C.Pstand), missing)
    elseif isa(C, MatrixProductCode)
        nr, nc = size(C.A)
        # probably not going to work
        nr == nc || return LinearCode.dual(C)
        D = Vector{LinearCode}()
        for i in 1:length(C.C)
            push!(D, dual(C.C[i]))
        end
        
        try
            Ainv = inv(C.A)
        catch
            return LinearCode.dual(C)
        end
        return MatrixProductCode(D, transpose(Ainv))
    elseif isa(C, ReedMullerCode)
        d = 2^(C.r + 1)
        return ReedMullerCode(C.F, C.n, C.n - C.k, d, d, d, C.m - C.r - 1, C.m, C.H, C.G,
            C.Hstand, C.Gstand, C.Pstand, missing)
    else
        G = deepcopy(generatormatrix(C))
        H = deepcopy(paritycheckmatrix(C))

        Hstand, Gstand, P, _ = _standardform(C.H)
        Pstandold = ismissing(C.Pstand) ? identity_matrix(C.F, C.n) : C.Pstand
        Pstand = ismissing(P) ? Pstandold : Pstandold * P

        # Possible alternative to above, might be faster:
        # Gstandold = generatormatrix(C, true)
        # Gstand = hcat(view(Gstandold, :, C.k + 1:C.n), view(Gstandold, :, 1:C.k))
        # Hstandold = paritycheckmatrix(C, true)
        # Hstand = hcat(view(Hstandold, :, C.k + 1:C.n), view(Hstandold, :, 1:C.k))
        # Pstandold = ismissing(C.Pstand) ? identity_matrix(C.F, C.n) : C.Pstand
        # Pstand = vcat(view(C.Pstand, C.k + 1:C.n, :), view(C.Pstand, 1:C.k, :))

        if !ismissing(C.weightenum)
            dualwtenum = MacWilliamsIdentity(C, C.weightenum)
            dualHWEpoly = CWEtoHWE(dualwtenum).polynomial
            d = minimum(filter(ispositive, first.(exponent_vectors(dualHWEpoly))))
            return LinearCode(C.F, C.n, C.n - C.k, d, d, d, H, G,
                Hstand, Gstand, Pstand, dualwtenum)
        else
            ub1, _ = _minwtrow(H)
            ub2, _ = _minwtrow(Hstand)
            ub = min(ub1, ub2)
            return LinearCode(C.F, C.n, C.n - C.k, missing, 1, ub, H,
                            G, Hstand, Gstand, Pstand, missing)
        end
    end
end
Euclideandual(C::AbstractLinearCode) = dual(C)

"""
    Hermitiandual(C::AbstractLinearCode)

Return the Hermitian dual of a code defined over a quadratic extension.
"""
function Hermitiandual(C::AbstractLinearCode)
    if isa(C, MatrixProductCode)
        # the inner functions here should complain if not quadratic, so don't have to check here
        nr, nc = size(C.A)
        # probably not going to work
        nr == nc || return LinearCode.Hermitiandual(C)
        D = Vector{LinearCode}()
        for i in 1:length(C.C)
            push!(D, Hermitiandual(C.C[i]))
        end

        try
            Ainv = inv(Hermitianconjugatematrix(C.A))
        catch
            return LinearCode.dual(C)
        end
        return MatrixProductCode(D, transpose(Ainv))
    else
        return dual(LinearCode(Hermitianconjugatematrix(C.G)))
    end
end

"""
    areequivalent(C1::AbstractLinearCode, C2::AbstractLinearCode)

Return `true` if `C1 ⊆ C2` and `C2 ⊆ C1`.
"""
areequivalent(C1::AbstractLinearCode, C2::AbstractLinearCode) = (C1 ⊆ C2) && (C2 ⊆ C1)

"""
    arepermutationequivalent(C1::AbstractLinearCode, C2::AbstractLinearCode)

Returns a `Tuple{Bool, Union{Missing, Perm}}`. The first entry is
`true` if `C1` and `C2` are permutation equivalent codes. The second
entry gives the permutation that, when applied to the columns of the
generator matrix of C1, would give C2.
"""
function arepermutationequivalent(C1::AbstractLinearCode, C2::AbstractLinearCode)
    G1 = generatormatrix(C1, true)
    G2 = generatormatrix(C2, true)
    if G1 == G2
        if ismissing(C1.Pstand) && ismissing(C2.Pstand)
            (true, Perm(1:C1.n))
        elseif ismissing(C1.Pstand)
            (true, Perm(Int.(data.(transpose(C2.Pstand))) * collect(1:C1.n)))
        elseif ismissing(C2.Pstand)
            (true, Perm(Int.(data.(C1.Pstand)) * collect(1:C1.n)))
        else
            (true, Perm(Int.(data.(transpose(C2.Pstand) * C1.Pstand)) * collect(1:C1.n)))
        end
    else
        (false, missing)
    end
end

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
# LGaloisdual(C::AbstractLinearCode) = dual(LinearCode(LGaloisconjugatematrix(C.G)))
# need LGaloisconjugatematrix(C.G)

"""
    characteristicpolynomial(C::AbstractLinearCode)

Return the characteristic polynomial of `C`.

# Notes
* The characteristic polynomial is defined in [Lin1999]_
"""
function characteristicpolynomial(C::AbstractLinearCode)
    _, x = PolynomialRing(Nemo.QQ, :x)
    D = dual(C)
    supD = support(D)
    q = Int(order(C.F))
    return q^(n - k) * prod(1 - x / j for j in supD if j > 0)
end

"""
    VectorSpace(C::AbstractLinearCode)

Return the code `C` as a vector space object.
"""
function VectorSpace(C::AbstractLinearCode)
    V = VectorSpace(C.F, C.n)
    G = generatormatrix(C)
    return sub(V, [V(view(G, i:i, :)) for i in 1:nrows(G)])
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
    return all(wt(view(G, r:r, :)) % 2 == 0 for r in 1:nrows(G))
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
    all(wt(view(G, r:r, :)) % 4 == 0 for r in 1:nr) || (return false;)
    all(wt(view(G, r1:r1, :) + view(G, r2:r2, :)) % 4 == 0
        for r1 in r:nr, r2 in 1:nr) || (return false;)
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
    all(wt(view(G, r:r, :)) % 8 == 0
        for r in 1:nr) || (return false;)
    all(wt(view(G, r1:r1, :) .* view(G, r2:r2, :)) % 4 == 0
        for r1 in 1:nr, r2 in 1:nr) || (return false;)
    all(wt(view(G, r1:r1, :) .* view(G, r2:r2, :) .* view(G, r3:r3, :)) % 2 == 0
        for r1 in 1:nr, r2 in 1:nr, r3 in 1:nr) || (return false;)
    return true
end

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
    words = onlyprint ? nothing : Vector{typeof(C.G)}()
    G = ismissing(C.Pstand) ? generatormatrix(C, true) : generatormatrix(C, true) * C.Pstand
    E = base_ring(G)

    if iszero(G)
        row = zero_matrix(E, 1, C.n)
        onlyprint ? println(row) : push!(words, row)
        return words
    end

    # for iter in Iterators.product(Iterators.repeated(E, nrows(G))...)
    for iter in Nemo.AbstractAlgebra.ProductIterator([E for _ in 1:nrows(G)], inplace = true)
        row = iter[1] * view(G, 1:1, :)
        for r in 2:nrows(G)
            if !iszero(iter[r])
                row += iter[r] * view(G, r:r, :)
            end
        end
        onlyprint ? println(row) : push!(words, row)
    end
    return words
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
