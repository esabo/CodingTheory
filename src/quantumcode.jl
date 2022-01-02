# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

include("linearcode.jl")
include("trellis.jl")

# using additive over stabilizer
abstract type AbstractAdditiveCode <: AbstractCode end
abstract type AbstractQuantumCode <: AbstractAdditiveCode end
abstract type AbstractCSSCode <: AbstractQuantumCode end
# can also build a AbstractQuantumLinearCode if the additive code is also linear
# will probably have to tweak all of these later

mutable struct CSSCode <: AbstractCSSCode
    F::FqNmodFiniteField # base field (symplectic)
    E::FqNmodFiniteField # additive field
    n::Integer
    k::Integer
    d::Union{Integer, Missing}
    dx::Union{Integer, Missing}
    dz::Union{Integer, Missing}
    stabs::fq_nmod_mat
    Xstabs::fq_nmod_mat
    Zstabs::fq_nmod_mat
    Xorigcode::Union{LinearCode, Missing}
    ZorigCode::Union{LinearCode, Missing}
    signs::Vector{Int64} # make bools, uint8's?
    Xsigns::Vector{Int64} # make bools, uint8's?
    Zsigns::Vector{Int64} # make bools, uint8's?
    dualgens::fq_nmod_mat
    logspace::Union{fq_nmod_mat, Missing}
    charvec::Vector{Int64} # make bools, uint8's?
    Hwtenum::Union{WeightEnumerator, Missing}
    Pwtenum::Union{WeightEnumerator, Missing}
    trellis::Union{Trellis, Missing}
end

mutable struct QuantumCode <: AbstractQuantumCode
    F::FqNmodFiniteField # base field (symplectic)
    E::FqNmodFiniteField # additive field
    n::Integer
    k::Union{Integer, Rational{Int64}}
    d::Union{Integer, Missing}
    stabs::fq_nmod_mat
    dualgens::fq_nmod_mat
    logspace::Union{fq_nmod_mat, Missing}
    charvec::Vector{Int64} # make bools, uint8's?
    signs::Vector{Int64}
    Hwtenum::Union{WeightEnumerator, Missing}
    Pwtenum::Union{WeightEnumerator, Missing}
    Xwtenum::Union{WeightEnumerator, Missing}
    Zwtenum::Union{WeightEnumerator, Missing}
    trellis::Union{Trellis, Missing}
end

field(S::AbstractQuantumCode) = S.F
quadraticfield(S::AbstractQuantumCode) = S.E
length(S::AbstractQuantumCode) = S.n
numqubits(S::AbstractQuantumCode) = S.n
dimension(S::AbstractQuantumCode) = S.k
signs(S::AbstractQuantumCode) = S.signs
Xsigns(S::CSSCode) = S.Xsigns
Zsigns(S::CSSCode) = S.Zsigns
stabilizers(S::AbstractQuantumCode) = S.stabs
Xstabilizers(S::CSSCode) = S.Xstabs
Zstabilizers(S::CSSCode) = S.Zstabs
numXstabs(S::CSSCode) = size(S.Xstabs, 1)
numZstabs(S::CSSCode) = size(S.Zstabs, 1)
normalizermatrix(S::AbstractQuantumCode) = S.dualgens
charactervector(S::AbstractQuantumCode) = S.charvec
minimumdistanceX(S::CSSCode) = S.dx
minimumdistanceZ(S::CSSCode) = S.dz
minimumdistance(S::CSSCode) = S.d

function quadratictosymplectic(M::fq_nmod_mat)
    E = base_ring(M)
    iseven(degree(E)) || error("The base ring of the given matrix is not a quadratic extension.")
    F, _ = FiniteField(Int64(characteristic(E)), div(degree(E), 2), "ω")
    nrows = size(M, 1)
    ncols = size(M, 2)
    Msym = zero_matrix(F, nrows, 2 * ncols)
    for r in 1:nrows
        for c in 1:ncols
            if !iszero(M[r, c])
                Msym[r, c] = F(coeff(M[r, c], 0))
                Msym[r, c + ncols] = F(coeff(M[r, c], 1))
            end
        end
    end
    return Msym
end

function symplectictoquadratic(M::fq_nmod_mat)
    iseven(size(M, 2)) || error("Input to symplectictoquadratic is not of even length.")
    nrows = size(M, 1)
    ncols = div(size(M, 2), 2)
    F = base_ring(M)
    E, ω = FiniteField(Int64(characteristic(F)), 2 * degree(F), "ω")
    ϕ = embed(F, E)
    Mquad = zero_matrix(E, nrows, ncols)
    for r in 1:nrows
        for c in 1:ncols
            Mquad[r, c] = ϕ(M[r, c]) + ϕ(M[r, c + ncols]) * ω
        end
    end
    return Mquad
end

function _Paulistringtosymplectic(str::T) where T <: Union{String, Vector{Char}}
    n = length(str)
    F, _ = FiniteField(2, 1, "ω")
    sym = zero_matrix(F, 1, 2 * n)
    for (i, c) in enumerate(str)
        if c == 'X'
            sym[1, i] = F(1)
        elseif c == 'Z'
            sym[1, i + n] = F(1)
        elseif c == 'Y'
            sym[1, i] = 1
            sym[1, i + n] = F(1)
        elseif c != 'I'
            error("Encountered non-{I, X, Y, Z} character in Pauli string. This function is only defined for binary strings.")
        end
    end
    return sym
end
_Paulistringtosymplectic(A::Vector{T}) where T <: Union{String, Vector{Char}} = vcat([_Paulistringtosymplectic(s) for s in A]...)
_Paulistringstofield(str::T) where T <: Union{String, Vector{Char}} = symplectictoquadratic(_Paulistringtosymplectic(str))
_Paulistringstofield(A::Vector{T}) where T <: Union{String, Vector{Char}} = vcat([_Paulistringstofield(s) for s in A]...)

# need symplectictoPaulistring
# quadratictoPaulistring

symplecticstabilizers(S::AbstractQuantumCode) = quadratictosymplectic(S)

function _getsigns(A::fq_nmod_mat, charvec::Vector{Int64})
    if length(charvec) == 2 * size(A, 2)
        A = quadratictosymplectic(A)
    end
    length(charvec) == size(A, 2) || error("Input to _getsigns is expected to be in symplectic form and of the same length as the characteristic vector.")
    if charvec == ones(size(A, 2))
        return ones(Int64, div(size(A, 2), 2))
    end

    signs = Vector{Int64}()
    for r in 1:size(A, 1)
        parity = 1
        for c = 1:size(A, 2)
            if !iszero(A[r, c]) && charvec[c] == -1
                parity *= -1
            end
        end
        append!(signs, parity)
    end
    return signs
end

function _processstrings(SPauli::Vector{T}, charvec::Union{Vector{Int64}, Vector{Any}}=[]) where T <: Union{String, Vector{Char}}
    Paulisigns = Vector{Int64}()
    SPaulistripped = Vector{String}()
    for (i, s) in enumerate(SPauli)
        if s[1] ∈ ['I', 'X', 'Y', 'Z']
            append!(Paulisigns, 1)
            push!(SPaulistripped, s)
        elseif s[1] == '+'
            append!(Paulisigns, 1)
            push!(SPaulistripped, s[2:end])
        elseif s[1] == '-'
            append!(Paulisigns, -1)
            push!(SPaulistripped, s[2:end])
        else
            error("The first element of Pauli string $i is neither a Pauli character or +/-: $s.")
        end
    end

    for s in SPaulistripped
        for i in s
            if i ∉ ['I', 'X', 'Y', 'Z']
                error("Element of provided Pauli string is not a Pauli character: $s.")
            end
        end
    end

    n = length(SPaulistripped[1])
    for s in SPaulistripped
        length(s) == n || error("Provided Pauli strings are not all of the same length.")
    end

    if isempty(charvec)
        charvec = ones(Int64, 2 * n)
    else
        if length(charvec) != 2 * n
            error("Unexpected length of character vector. Expected: $(2 * n), received: $(length(charvec)).")
        end

        for s in charvec
            (s == 1 || s == -1) || error("Qubit phases must be +/- 1; received: $s")
        end
    end
    return SPaulistripped, charvec
end

function splitsymplecticstabilizers(S::fq_nmod_mat, signs::Vector{Int64}=[])
    Xstabs = []
    Xsigns = Vector{Int64}()
    Zstabs = []
    Zsigns = Vector{Int64}()
    mixedstabs = []
    mixedsigns = Vector{Int64}()

    half = div(size(S, 2), 2)
    for r in 1:size(S, 1)
        # use views?
        s = S[r, :]
        if iszero(s)
            continue
        else
            sx = iszero(s[1, 1:half])
            sz = iszero(s[1, half + 1:end])
            if (sx && !sz)
                push!(Zstabs, s)
                append!(Zsigns, signs[r])
            elseif !sx && sz
                push!(Xstabs, s)
                append!(Xsigns, signs[r])
            elseif !sx && !sz
                push!(mixedstabs, s)
                append!(mixedsigns, signs[r])
            end
        end
    end

    if !isempty(Xstabs)
        Xstabs = vcat(Xstabs...)
    end
    if !isempty(Zstabs)
        Zstabs = vcat(Zstabs...)
    end
    if !isempty(mixedstabs)
        mixedstabs = vcat(mixedstabs...)
    end
    return Xstabs, Xsigns, Zstabs, Zsigns, mixedstabs, mixedsigns
end

function isCSSsymplectic(S::fq_nmod_mat, signs::Vector{Int64}=[], trim::Bool=true)
    Xstabs, Xsigns, Zstabs, Zsigns, mixedstabs, mixedsigns = splitsymplecticstabilizers(S, signs)
    if isempty(mixedstabs)
        if trim
            half = div(size(Xstabs, 2), 2)
            return true, Xstabs[:, 1:half], Xsigns, Zstabs[:, half + 1:end], Zsigns
        else
            return true, Xstabs, Xsigns, Zstabs, Zsigns
        end
    else
        if trim
            half = div(size(Xstabs, 2), 2)
            return false, Xstabs[:, 1:half], Xsigns, Zstabs[:, half + 1:end], Zsigns, mixedstabs, mixedsigns
        else
            return false, Xstabs, Xsigns, Zstabs, Zsigns, mixedstabs, mixedsigns
        end
    end
end

function CSSCode(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{Int64}, Vector{Any}}=[])
    length(C1) ==  length(C2) || error("Both codes must have the same length in the CSS construction.")
    field(C1) == field(C2) || error("Both codes must be over the same base field in the CSS construction.")
    C2 ⊆ C1 || error("The second argument must be a subset of the first in the CSS construction.")
    if !isempty(charvec)
        for s in charvec
            (s == 1 || s == -1) || error("X-stabilizer signs must be +/- 1; received: $s")
        end
    end

    # C2 ⊆ C1
    # k = k1 - k2
    # d >= minimum(d1, d2^⟂)
    # X - H(C2^⟂), Z - H(C1)
    D2 = dual(C2)
    S = directsum(paritycheckmatrix(D2), paritycheckmatrix(C1))
    Sq2 = symplectictoquadratic(S)
    if isempty(charvec)
        charvec = ones(Int64, 2 * length(C1))
        signs = ones(Int64, size(S, 1))
        Xsigns = ones(Int64, size(paritycheckmatrix(D2), 1))
        Zsigns = ones(Int64, size(paritycheckmatrix(C1), 1))
    else
        signs = _getsigns(S, charvec)
        Xsigns = signs[1:size(paritycheckmatrix(D2), 1), :]
        Zsigns = signs[size(paritycheckmatrix(D2), 1) + 1:end, :]
    end
    dualgens = directsum(generatormatrix(C1), generatormatrix(D2))

    return CSSCode(field(C1), base_ring(Sq2), length(C1), dimension(C1) - dimension(C2), missing,
        minimumdistance(D2), minimumdistance(C1), Sq2, paritycheckmatrix(D2), paritycheckmatrix(C1),
        C2, C1, signs, Xsigns, Zsigns, dualgens, missing, charvec, missing, missing, missing)
end

function CSSCode(C::LinearCode, charvec::Union{Vector{Int64}, Vector{Any}}=[])
    # this should have Xstabs = Zstabs
    D = dual(C)
    C ⊆ D || error("The single code CSS construction requires C ⊆ C^⟂.")
    if !isempty(charvec)
        for s in charvec
            (s == 1 || s == -1) || error("X-stabilizer signs must be +/- 1; received: $s")
        end
    end

    # C2 ⊆ C1
    # k = k1 - k2
    # d >= minimum(d1, d2^⟂)
    # X - H(C2^⟂), Z - H(C1)
    S = directsum(paritycheckmatrix(D), paritycheckmatrix(D))
    Sq2 = symplectictoquadratic(S)
    k = dimension(D) - dimension(C)
    if isempty(charvec)
        charvec = ones(Int64, 2 * length(C))
        signs = ones(Int64, size(S, 1))
        Xsigns = ones(Int64, size(paritycheckmatrix(D), 1))
        Zsigns = ones(Int64, size(paritycheckmatrix(D), 1))
    else
        signs = _getsigns(S, charvec)
        Xsigns = signs[1:size(paritycheckmatrix(D), 1), :]
        Zsigns = signs[size(paritycheckmatrix(D), 1) + 1:end, :]
    end
    dualgens = directsum(generatormatrix(D), generatormatrix(D))

    return CSSCode(field(C), base_ring(Sq2), length(C), k, missing, minimumdistance(D),
        minimumdistance(D), Sq2, paritycheckmatrix(D), paritycheckmatrix(D), C, D,
        signs, Xsigns, Zsigns, dualgens, missing, charvec, missing, missing, missing)
end

function CSSCode(Xmatrix::fq_nmod_mat, Zmatrix::fq_nmod_mat, charvec::Union{Vector{Int64}, Vector{Any}}=[])
    # this should have Zstabs * Xstabs' = 0
    # set dual containing if Xstabs == Zstabs, else not
    size(Xmatrix, 2) ==  size(Zmatrix, 2) || error("Both matrices must have the same length in the CSS construction.")
    base_ring(Xmatrix) == base_ring(Zmatrix) || error("Both matrices must be over the same base field in the CSS construction.")
    # TODO can remove rank check now that we have full character vector
    rank(Xmatrix) == size(Xmatrix, 1) || error("Provided X-stabilizer matrix is not full rank.")
    rank(Zmatrix) == size(Zmatrix, 1) || error("Provided Z-stabilizer matrix is not full rank.")
    if !isempty(charvec)
        for s in charvec
            (s == 1 || s == -1) || error("X-stabilizer signs must be +/- 1; received: $s")
        end
    end

    S = directsum(Xmatrix, Zmatrix)
    aresymplecticorthogonal(S, S) || error("The given matrices are not symplectic orthogonal.")
    Sq2 = symplectictoquadratic(S)
    if isempty(charvec)
        charvec = ones(Int64, size(S, 2))
        signs = ones(Int64, size(S, 1))
        Xsigns = ones(Int64, size(Xmatrix, 1))
        Zsigns = ones(Int64, size(Zmatrix, 1))
    else
        signs = _getsigns(S, charvec)
        Xsigns = signs[1:size(Xmatrix, 1), :]
        Zsigns = signs[size(Xmatrix, 1) + 1:end, :]
    end

    G = hcat(S[:, size(Xmatrix, 2) + 1:end], -S[:, 1:size(Xmatrix, 2)])
    _, H = right_kernel(G)
    size(H', 1) == 2 * size(Sq2, 2) - size(Sq2, 1) || error("Normalizer matrix is not size n + k.")
    # note the H here is transpose of the standard definition
    !(!iszero(G * H) || !iszero(H' * G')) || error("Normalizer matrix is not transpose, symplectic orthogonal.")
    dualgens = symplectictoquadratic(H')

    return CSSCode(base_ring(Xmatrix), base_ring(Sq2), size(Xmatrix, 2),
        size(Xmatrix, 2) - size(Xmatrix, 1) - size(Zmatrix, 1),
        missing, missing, missing, Sq2, missing, missing, missing, missing, signs,
        Xsigns, Zsigns, dualgens, missing, charvec, missing, missing, missing)
end

function CSSCode(SPauli::Vector{T}, charvec::Union{Vector{Int64}, Vector{Any}}=[]) where T <: Union{String, Vector{Char}}
    SPaulistripped, charvec = _processstrings(SPauli, charvec)
    S = _Paulistringtosymplectic(SPaulistripped)
    aresymplecticorthogonal(S, S) || error("The given stabilizers are not symplectic orthogonal.")
    signs = _getsigns(S, charvec)
    Sq2 = symplectictoquadratic(S)

    del = Vector{Int64}()
    for r in 1:size(Sq2, 1)
        if iszero(Sq2[r, :])
            append!(del, r)
        end
    end

    flag1 = false
    # not necessary but faster to skip
    if !isempty(del)
        Sq2 = Sq2[setdiff(1:size(Sq2, 1), del), :]
        signs = signs[setdiff(1:size(Sq2, 1), del)]
        flag1 = true
    end

    ###############
    # need to figure out how to track signs as it drops rank
    # would be easier to reassign given character vector instead of stabilizer signs
    ###############

    # if rank(Sq2) != size(Sq2, 1)
    # rk = rank(G)
    # !iszero(rk) || error("Rank zero matrix passed into LinearCode constructor.")
    # if rk < size(G, 1)
    #     _, G = rref(G)
    # end

    # find dual
    if flag1
        S = quadratictosymplectic(Sq2)
    end
    G = hcat(S[:, size(Sq2, 2) + 1:end], -S[:, 1:size(Sq2, 2)])
    _, H = right_kernel(G)
    size(H', 1) == 2 * size(Sq2, 2) - size(Sq2, 1) || error("Normalizer matrix is not size n + k.")
    # note the H here is transpose of the standard definition
    !(!iszero(G * H) || !iszero(H' * G')) || error("Normalizer matrix is not transpose, symplectic orthogonal.")
    dualgens = symplectictoquadratic(H')

    args = isCSSsymplectic(S, signs, true)
    if args[1]
        return CSSCode(base_ring(S), base_ring(Sq2), size(Sq2, 2), size(Sq2, 2) - size(Sq2, 1),
            missing, missing, missing, Sq2, args[2], args[4], missing, missing, signs,
            args[3], args[5], dualgens, missing, charvec, missing, missing, missing)
    else
        error("Provided Pauli strings are not CSS.")
    end
end

# entanglement-assisted is not symplectic orthogonal
function QuantumCode(SPauli::Vector{T}, charvec::Union{Vector{Int64}, Vector{Any}}=[]) where T <: Union{String, Vector{Char}}
    SPaulistripped, charvec = _processstrings(SPauli, charvec)
    S = _Paulistringtosymplectic(SPaulistripped)
    aresymplecticorthogonal(S, S) || error("The given stabilizers are not symplectic orthogonal.")
    signs = _getsigns(S, charvec)
    Sq2 = symplectictoquadratic(S)

    # make del zero rows, cols helper functions
    # which will return a flag if this was ever done
    del = Vector{Int64}()
    for r in 1:size(Sq2, 1)
        if iszero(Sq2[r, :])
            append!(del, r)
        end
    end

    flag1 = false
    # not necessary but faster to skip
    if !isempty(del)
        Sq2 = Sq2[setdiff(1:size(Sq2, 1), del), :]
        signs = signs[setdiff(1:size(Sq2, 1), del)]
        flag1 = true
    end

    # need to code an additive only rref but then apply _getsigns
    # if rank(Sq2) != size(Sq2, 1)
    # rk = rank(G)
    # !iszero(rk) || error("Rank zero matrix passed into LinearCode constructor.")
    # if rk < size(G, 1)
    #     _, G = rref(G)
    # end

    # find dual
    if flag1
        S = quadratictosymplectic(Sq2)
    end
    G = hcat(S[:, size(Sq2, 2) + 1:end], -S[:, 1:size(Sq2, 2)])
    _, H = right_kernel(G)
    size(H', 1) == 2 * size(Sq2, 2) - size(Sq2, 1) || error("Normalizer matrix is not size n + k.")
    # note the H here is transpose of the standard definition
    !(!iszero(G * H) || !iszero(H' * G')) || error("Normalizer matrix is not transpose, symplectic orthogonal.")
    dualgens = symplectictoquadratic(H')

    # q^n / p^k but rows is n - k
    dimcode = Int64(order(base_ring(S)))^size(Sq2, 2) // Int64(characteristic(base_ring(S)))^(size(Sq2, 1))
    if isinteger(dimcode)
        dimcode = Int64(log(Int64(characteristic(base_ring(S))), Int64(dimcode)))
    end

    args = isCSSsymplectic(S, signs, true)
    if args[1]
        return CSSCode(base_ring(S), base_ring(Sq2), size(Sq2, 2), dimcode,
            missing, missing, missing, Sq2, args[2], args[4], missing, missing, signs,
            args[3], args[5], dualgens, missing, charvec, missing, missing, missing)
    else
        return QuantumCode(base_ring(S), base_ring(Sq2), size(Sq2, 2), dimcode,
            missing, Sq2, dualgens, missing, charvec, signs, missing, missing, missing, missing, missing)
    end
end

# function QuantumCode()
#
#     # determine K
#     K = q^n // characteristic()^size(Ssomething, 1)
#     if isinteger(K)
#         K = Int(log(Int64(K), q))
#     end
# end

# this returns a basis for the space S^⟂ / S
# work would need to be done to extract appropriate X_i, Z_i pairs from this
# can implement Gottesman algorithm later
function logicalspace(S::AbstractQuantumCode, alg::String="quotient")
    if !ismissing(Q.logicalspace)
        return Q.logicalspace
    end

    alg ∈ ["quotient"] || error("Algorithm '$alg' not specified in logicalspace.")
    if alg == "quotient"
        F = quadraticfield(S)
        G = stabilizers(S)
        Gdual = normalizermatrix(S)
        V = VectorSpace(F, length(S))
        U, UtoV = sub(V, [V(G[i, :]) for i in 1:size(G, 1)])
        W, WtoV = sub(V, [V(Gdual[i, :]) for i in 1:size(Gdual, 1)])
        gensofUinW = [preimage(WtoV, UtoV(g)) for g in gens(U)]
        UinW, UinWtoW = sub(W, gensofUinW)
        Q, WtoQ = quo(W, UinW)
        C2modC1basis = [WtoV(x) for x in [preimage(WtoQ, g) for g in gens(Q)]]
        Fbasis = [[F(C2modC1basis[j][i]) for i in 1:dim(parent(C2modC1basis[1]))] for j in 1:length(C2modC1basis)]
        G2 = matrix(F, length(Fbasis), length(Fbasis[1]), vcat(Fbasis...))

        # returns linear object, make additive if necessary
        if size(G2, 1) == dimension(S)
            Q.logicalspace = vcat(G2, gen(F) * G2)
            return Q.logicalspace
        elseif size(G2, 1) == 2 * dimension(S)
            Q.logicalspace = G2
            return Q.logicalspace
        else
            error("Logical space produced of incorrect dimension; expected: ", 2 * dimension(S), ", received: ", size(G2, 1))
        end
    end
end

function setlogicals!(Q::AbstractQuantumCode, L::fq_nmod_mat)
    size(L) == (length(Q) + dimension(Q), length(Q)) || error("Provided matrix is of incorrect size for the logical space.")
    rank(L) == size(L, 1) || error("Provided matrix is not of full rank.")
    aresymplecticorthogonal(stabilizers(Q), L) || error("Provided matrix does not commute with the code.")
    !aresymplecticorthogonal(L, L) || error("Provided matrix should not be symplectic self-orthogonal.")
    Q.logspace = L
end

# update for roots of unity
function changesigns!(Q::AbstractQuantumCode, charvec::Vector{Int64})
    length(charvec) == 2 * length(Q) || error("Characteristic vector is of improper length for the code.")
    for s in charvec
        (s == 1 || s == -1) || error("Qubit phases must be +/- 1; received: $s")
    end
    Q.signs = _getsigns(stabilizers(Q), charvec)
    Q.charvec = charvec
end


# add for minimum distance check
function show(io::IO, S::AbstractQuantumCode)
    if get(io, :compact, false)
        if typeof(S) <: CSSCode
            if typeof(dimension(S)) <: Integer
                println(io, "[[$(length(S)), $(dimension(S))]]_$(order(field(S))) CSS code.")
            else # don't think this can ever be reached
                println(io, "(($(length(S)), $(dimension(S))))_$(order(field(S))) CSS code.")
            end
        else
            if typeof(dimension(S)) <: Integer
                println(io, "[[$(length(S)), $(dimension(S))]]_$(order(field(S))) quantum (additive) code.")
            else
                println(io, "(($(length(S)), $(dimension(S))))_$(order(field(S))) quantum (additive) code.")
            end
        end
    else
        if typeof(S) <: CSSCode
            if typeof(dimension(S)) <: Integer
                println(io, "[[$(length(S)), $(dimension(S))]]_$(order(field(S))) CSS code.")
            else # don't think this can ever be reached
                println(io, "(($(length(S)), $(dimension(S))))_$(order(field(S))) CSS code.")
            end
            println(io, "X-stabilizer matrix: $(numXstabs(S)) × $(length(S))")
            for i in 1:numXstabs(S)
                print(io, "\t")
                for j in 1:length(S)
                    if j != length(S)
                        print(io, "$(S.Xstabs[i, j]) ")
                    elseif j == length(S) && i != length(S) - dimension(S)
                        println(io, "$(S.Xstabs[i, j])")
                    else
                        print(io, "$(S.Xstabs[i, j])")
                    end
                end
            end
            println(io, "Z-stabilizer matrix: $(numZstabs(S)) × $(length(S))")
            for i in 1:numZstabs(S)
                print(io, "\t")
                for j in 1:length(S)
                    if j != length(S)
                        print(io, "$(S.Zstabs[i, j]) ")
                    elseif j == length(S) && i != length(S) - dimension(S)
                        println(io, "$(S.Zstabs[i, j])")
                    else
                        print(io, "$(S.Zstabs[i, j])")
                    end
                end
            end
        else
            if typeof(dimension(S)) <: Integer
                println(io, "[[$(length(S)), $(dimension(S))]]_$(order(field(S))) quantum (additive) code.")
            else
                println(io, "(($(length(S)), $(dimension(S))))_$(order(field(S))) quantum (additive) code.")
            end
            println(io, "Stabilizer matrix: $(size(S.stabs, 1)) × $(length(S))")
            for i in 1:size(S.stabs, 1)
                print(io, "\t")
                for j in 1:length(S)
                    if j != length(S)
                        print(io, "$(S.stabs[i, j]) ")
                    elseif j == length(S) && i != length(S) - dimension(S)
                        println(io, "$(S.stabs[i, j])")
                    else
                        print(io, "$(S.stabs[i, j])")
                    end
                end
            end
        end
    end
end

# # these only CSS otherwise it's just numlogicals
# function numXlogicals(S::CSSCode)
#     return
# end
#
# function numZlogicals(S::CSSCode)
#     return
# end
#
# # just use 2k
# numlogicals(S::AbstractQuantumCode) = numXlogicals(S) + numZlogicals(S)

function Xsyndrome(S::CSSCode, v::fq_nmod_mat)
    n = length(S)
    if length(v) == 2 * n
        v = v[n + 1:end]
    end
    !(size(v) != (n, 1) && size(v) != (1, n)) ||
        error("Vector to be tested is of incorrect dimension; expected length $n, received: $(size(v)).")
    base_ring(v) == field(S) || error("Vector must have the same base ring as the stabilizers.")

    size(v, 1) != 1 || return Xstabilizers(S) * v'
    return Xstabilizers(S) * v
end

function Zsyndrome(S::CSSCode, v::fq_nmod_mat)
    n = length(S)
    if length(v) == 2 * n
        v = v[1:n]
    end
    !(size(v) != (n, 1) && size(v) != (1, n)) ||
        error("Vector to be tested is of incorrect dimension; expected length $n, received: $(size(v)).")
    base_ring(v) == field(S) || error("Vector must have the same base ring as the stabilizers.")

    size(v, 1) != 1 || return Zstabilizers(S) * v'
    return Zstabilizers(S) * v
end

function syndrome(S::AbstractQuantumCode, v::fq_nmod_mat)
    n = length(S)
    !(size(v) != (2 * n, 1) && size(v) != (1, 2 * n)) ||
        error("Vector to be tested is of incorrect dimension; expected length $(2 * n), received: $(size(v)).")
    # base_ring(v) == field(S) || error("Vector must have the same base ring as the stabilizers.")

    size(v, 1) != 1 || return symplecticstabilizers(S) * v'
    return symplecticstabilizers(S) * v
end

#############################
        # Known codes
#############################

function fivequbitcode()
    # is a perfect code
    return QuantumCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"])
end

# should also do a test for other CSS construction via Hamming code and actually make that one default
function Steanecode()
    return CSSCode(["XXXXIII", "XXIIXXI", "XIXIXIX", "ZZZZIII", "ZZIIZZI", "ZIZIZIZ"])
    # also ZZIZZII, ZIZZIZI, IZZZIIZ, XXIXXII, XIXXIXI, IXXXIIX
end

function Shorcode()
    return CSSCode(["ZZIIIIIII", "IZZIIIIII", "IIIZZIIII", "IIIIZZIII", "IIIIIIZZI", "IIIIIIIZZ", "XXXXXXIII", "IIIXXXXXX"])
end

function Singletonbound()
    n - k ≧ 2(d - 1)
end

#############################
        # Trellises
#############################

# think of more scenarios
# could allow general trellises given partial stabilizers for use in trellis product
function trellisprofiles(Q::AbstractQuantumCode, type::String="weight")
    type ∈ ["weight", "decoding"] || error("Unknown type parameter in trellisprofiles.")
    if type == "weight"
        return _trellisprofiles(Q.dualgens, stabilizers(Q))
    else
        return _trellisprofiles(stabilizers(Q), Q.dualgens)
    end
end

# need to fix all of these stabilizer grabs
# I need to grab the stabs and make them quadratic before passing them
function trellisprofilesCSS(Q::CSSCode, type::String="weight", Pauli::Char=' ')
    type ∈ ["weight", "decoding"] || error("Unknown type parameter in trellisprofilesCSS.")
    Pauli ∈ ['X', 'Z'] || error("Pauli parameter in trellisprofilesCSS must be either 'X' or 'Z'.")
    if type == "weight"
        if Pauli == 'X'
            # this should be wrtV = SperpZ, wrtE = SX
            return _trellisprofiles(Q.dualgens, stabilizers(Q))
        else
            # this should be wrtV = SperpX, wrtE = SZ
            return _trellisprofiles(Q.dualgens, stabilizers(Q))
        end
    else
        if Pauli == 'X'
            # this should be wrtV = SX, wrtE = SperpZ
            return _trellisprofiles(stabilizers(Q), Q.dualgens)
        else
            # this should be wrtV = SZ, wrtE = SperpX
            return _trellisprofiles(stabilizers(Q), Q.dualgens)
        end
    end
end


############# need TOF here
function syndrometrellis(Q::AbstractQuantumCode, type::String="weight")
    type ∈ ["weight", "decoding"] || error("Unknown type parameter in syndrometrellis.")
    profiles = trellisprofiles(Q, type)
    if type == "weight"
        return _syndrometrellis(Q.dualgens, stabilizers(Q), profiles)
    else
        return _syndrometrellis(stabilizers(Q), Q.dualgens, profiles)
    end
end

# need to fix all of these stabilizer grabs
# I need to grab the stabs and make them quadratic before passing them
# how to pass CSS into this function? should separate?
function syndrometrellisCSSsymplectic(Q::CSSCode, type::String="weight", Pauli::Char=' ')
    type ∈ ["weight", "decoding"] || error("Unknown type parameter in syndrometrellisCSSsymplectic.")
    Pauli ∈ ['X', 'Z'] || error("Pauli parameter in syndrometrellisCSSsymplectic must be either 'X' or 'Z'.")
    profiles = trellisprofilesCSS(Q, type, Pauli)

    if type == "weight"
        if Pauli == 'X'
            # this should be wrtV = SperpZ, wrtE = SX
            return _syndrometrellis(Q.dualgens, stabilizers(Q), profiles)
        else
            # this should be wrtV = SperpX, wrtE = SZ
            return _syndrometrellis(Q.dualgens, stabilizers(Q), profiles)
        end
    else
        if Pauli == 'X'
            # this should be wrtV = SX, wrtE = SperpZ
            return _syndrometrellis(stabilizers(Q), Q.dualgens, profiles)
        else
            # this should be wrtV = SZ, wrtE = SperpX
            return _syndrometrellis(stabilizers(Q), Q.dualgens, profiles)
        end
    end
end

function _reducepoly(poly::Vector{Vector{Int64}})
    reducedpoly = Vector{Vector{Int64}}()
    processed = trues(length(poly))
    # println("received: $poly")
    # println("processed: $processed")
    for (i, term) in enumerate(poly)
        if processed[i]
            for j in (i + 1):length(poly)
                if processed[j]
                    # println("comparing $term and $(poly[j])")
                    if term[2] == poly[j][2] && term[3] == poly[j][3] && term[4] == poly[j][4]
                        term[1] += poly[j][1]
                        processed[j] = false
                    end
                end
            end
            push!(reducedpoly, term)
            processed[i] = false
        end
    end
    return reducedpoly
end

function islessPoly(a::Vector{Int64}, b::Vector{Int64})
    if isless(a[2:4], b[2:4])
        return true
    elseif a[2:4] == b[2:4]
        if isless(a[1], b[1])
            return true
        else
            return false
        end
    else
        return false
    end
end

function _HammingweightenumeratorQ(T::Trellis, cleanV::Bool=true)
    V = vertices(T)
    E = edges(T)
    for i in 2:length(V)
        for (j, v) in enumerate(V[i])
            outer = Vector{Vector{Int64}}()
            for e in E[i - 1][j]
                inner = deepcopy(V[i - 1][e.outvertex].polynomial)
                # println("edge: $(e.label), inner: $inner")
                # if e.label == "I"
                #     continue
                if e.label != "I"
                    for term in inner
                        term[2] += 1
                    end
                else
                    for term in inner
                        term[3] += 1
                    end
                end
                append!(outer, inner)
            end
            # println(outer)
            v.polynomial = _reducepoly(outer)
            # println(v.polynomial)
        end
    end
    W = WeightEnumerator(sort!(V[end][1].polynomial, lt=islessPoly))

    if cleanV
        for i in 2:length(V)
            for v in V[i]
                v.polynomial = [[0, 0, 0, 0]]
            end
        end
    end
    return W
end

function _PWEtoHWE(S::AbstractQuantumCode)
    poly = deepcopy(S.Pwtenum.polynomial)
    for term in poly
        tot = abs(term[2] + term[3] + term[4])
        term[2] = tot
        term[3] = length(S) - tot
        term[4] = 0
    end
    W = WeightEnumerator(_reducepoly(poly))
    if ismissing(S.Hwtenum)
        S.Hwtenum = W
    end
    return W
end

function _PWEtoXWE(S::AbstractQuantumCode)
    poly = deepcopy(S.Pwtenum.polynomial)
    for term in poly
        term[3] = length(S) - abs(term[2])
        term[4] = 0
    end
    W = WeightEnumerator(_reducepoly(poly))
    if ismissing(S.Xwtenum)
        S.Xwtenum = W
    end
    return W
end

function _PWEtoZWE(S::AbstractQuantumCode)
    poly = deepcopy(S.Pwtenum.polynomial)
    for term in poly
        term[2] = length(S) - abs(term[4])
        term[3] = 0
    end
    W = WeightEnumerator(_reducepoly(poly))
    if ismissing(S.Zwtenum)
        S.Zwtenum = W
    end
    return W
end

function Hammingweightenumerator(S::AbstractQuantumCode, keeptrellis::Bool=true, cleanV::Bool=true)
    if ismissing(S.Hwtenum)
        if !ismissing(S.Pwtenum)
            return _PWEtoHWE(S)
        elseif ismissing(trellis(S))
            STOF = trellisorientedform(stabilizers(S))
            nTOF = trellisorientedform(normalizermatrix(S))
            profiles = _trellisprofiles(nTOF, STOF)
            T = _syndrometrellis(profiles, nTOF, STOF, charactervector(S), ' ', false)
            if keeptrellis
                S.trellis = T
            end
        end
        S.Hwtenum = _HammingweightenumeratorQ(T)
    end
    return S.Hwtenum
end

function _Pauliweightenumerator(T::Trellis, cleanV::Bool=true)
    V = vertices(T)
    E = edges(T)
    for i in 2:length(V)
        for (j, v) in enumerate(V[i])
            outer = Vector{Vector{Int64}}()
            for e in E[i - 1][j]
                inner = deepcopy(V[i - 1][e.outvertex].polynomial)
                if e.label == "X"
                    for term in inner
                        if term[2] < 0 || term[3] < 0 || term[4] < 0
                            term[2] -= 1
                        else
                            term[2] += 1
                        end
                    end
                elseif e.label == "Y"
                    for term in inner
                        if term[2] < 0 || term[3] < 0 || term[4] < 0
                            term[3] -= 1
                        else
                            term[3] += 1
                        end
                    end
                elseif e.label == "Z"
                    for term in inner
                        if term[2] < 0 || term[3] < 0 || term[4] < 0
                            term[4] -= 1
                        else
                            term[4] += 1
                        end
                    end
                elseif e.label == "-X"
                    for term in inner
                        term[2] *= -1
                        term[3] *= -1
                        term[4] *= -1
                        if (term[2] < 0 || term[3] < 0 || term[4] < 0) || (term[2] == 0 && term[3] == 0 && term[4] == 0)
                            term[2] -= 1
                        else
                            term[2] += 1
                        end
                    end
                elseif e.label == "-Y"
                    for term in inner
                        term[2] *= -1
                        term[3] *= -1
                        term[4] *= -1
                        if (term[2] < 0 || term[3] < 0 || term[4] < 0) || (term[2] == 0 && term[3] == 0 && term[4] == 0)
                            term[3] -= 1
                        else
                            term[3] += 1
                        end
                    end
                elseif e.label == "-Z"
                    for term in inner
                        term[2] *= -1
                        term[3] *= -1
                        term[4] *= -1
                        if (term[2] < 0 || term[3] < 0 || term[4] < 0) || (term[2] == 0 && term[3] == 0 && term[4] == 0)
                            term[4] -= 1
                        else
                            term[4] += 1
                        end
                    end
                end
                append!(outer, inner)
            end
            v.polynomial = _reducepoly(outer)
        end
    end
    # sort is somehow broken
    # println(V[end][1].polynomial)
    # W = WeightEnumerator(sort!(V[end][1].polynomial, lt=islessPoly))
    W = WeightEnumerator(V[end][1].polynomial)
    # println(W)

    if cleanV
        for i in 2:length(V)
            for v in V[i]
                v.polynomial = [[0, 0, 0, 0]]
            end
        end
    end
    return W
end

# do I want to compute this here? probably, for consistency
function trellis(S::AbstractQuantumCode)
    if ismissing(S.trellis)
        return missing
    else
        return S.trellis
    end
end

function Pauliweightenumerator(S::AbstractQuantumCode, keeptrellis::Bool=true, cleanV::Bool=true)
    if ismissing(S.Pwtenum)
        if ismissing(trellis(S))
            STOF = trellisorientedform(stabilizers(S))
            nTOF = trellisorientedform(normalizermatrix(S))
            profiles = _trellisprofiles(nTOF, STOF)
            T = _syndrometrellis(profiles, nTOF, STOF, charactervector(S), ' ', false)
            if keeptrellis
                S.trellis = T
            end
        end
        S.Pwtenum = _Pauliweightenumerator(T)
    end
    _PWEtoHWE(S)
    _PWEtoXWE(S)
    _PWEtoZWE(S)
    return S.Pwtenum
end
