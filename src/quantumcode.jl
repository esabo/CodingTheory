# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

include("linearcode.jl")

# using additive over stabilizer
abstract type AbstractAdditiveCode <: AbstractCode end
abstract type AbstractQuantumCode <: AbstractAdditiveCode end
abstract type AbstractCSSCode <: AbstractQuantumCode end
# can also build a AbstractQuantumLinearCode if the additive code is also linear
# will probably have to tweak all of these later

# do I want to store properties like non-orthogonal or just have them as helper functions?
struct CSSCode <: AbstractCSSCode
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
end

struct QuantumCode <: AbstractQuantumCode
    F::FqNmodFiniteField # base field (symplectic)
    E::FqNmodFiniteField # additive field
    n::Integer
    k::Union{Integer, Rational{Int64}} # use isinteger in show
    d::Union{Integer, Missing}
    stabs::fq_nmod_mat
    dualgens::fq_nmod_mat
    signs::Vector{Int64} # make bools, uint8's?
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

function minimumdistanceX(S::CSSCode)
    !ismissing(S.dx) || error("Unknown minimum distance for this code.")
    return S.dx
end

function minimumdistanceZ(S::CSSCode)
    !ismissing(S.dz) || error("Unknown minimum distance for this code.")
    return S.dz
end

function minimumdistance(S::AbstractQuantumCode)
    !ismissing(S.d) || error("Unknown minimum distance for this code.")
    return S.d
end

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

function _processsigns(SPauli::Vector{T}, signs::Vector{Int64}=[]) where T <: Union{String, Vector{Char}}
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

    if isempty(signs)
        signs = Paulisigns
    else
        if length(signs) != length(SPaulistripped)
            error("Unexpected number of signs in CSSCode constructor. Expected: $(length(SPaulistripped)), received: $(length(signs)).")
        end

        for s in signs
            (s == 1 || s == -1) || error("Stabilizer signs must be +/- 1; received: $s")
        end

        temp = [1 for i in 1:length(signs)]
        if Paulisigns != temp
            signs == Paulisigns || error("Received contradictory sign information from the Pauli strings and the signs parameter.")
        end
    end
    return SPaulistripped, signs
end

function splitstabilizers(S::fq_nmod_mat, signs::Vector{Int64}=[])
    # iseven(size(S, 2)) || ()
    # if not, check quadratic degree and send to symplectic and then do


    ### need to handle signs first

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

function isCSS(S::fq_nmod_mat, signs::Vector{Int64}=[], trim::Bool=true)
    Xstabs, Xsigns, Zstabs, Zsigns, mixedstabs, mixedsigns = splitstabilizers(S, signs)
    if isempty(mixedstabs)
        if trim
            half = div(size(Xstabs, 2), 2)
            return true, Xstabs[:, 1:half], Xsigns, Zstabs[:, half + 1:end], Zsigns
        else
            return true, Xstabs, Xsigns, Zstabs, Zsigns
        end
    else
        if trim
            return false, Xstabs[:, 1:half], Xsigns, Zstabs[:, half + 1:end], Zsigns, mixedstabs, mixedsigns
        else
            return false, Xstabs, Xsigns, Zstabs, Zsigns, mixedstabs, mixedsigns
        end
    end
end

function CSSCode(C1::AbstractLinearCode, C2::AbstractLinearCode)
    length(C1) ==  length(C2) || error("Both codes must have the same length in the CSS construction.")
    field(C1) == field(C2) || error("Both codes must be over the same base field in the CSS construction.")
    C2 ⊆ C1 || error("The second argument must be a subset of the first in the CSS construction.")

    # C2 ⊆ C1
    # k = k1 - k2
    # d >= minimum(d1, d2^⟂)
    # X - H(C2^⟂), Z - H(C1)
    D2 = dual(C2)
    Sq2 = symplectictoquadratic(directsum(paritycheckmatrix(D2), paritycheckmatrix(C1)))
    return CSSCode(field(C1), base_ring(Sq2), length(C1), dimension(C1) - dimension(C2), missing,
        minimumdistance(D2), minimumdistance(C1), Sq2, paritycheckmatrix(D2), paritycheckmatrix(C1), C2, C1)
end

function CSSCode(C::LinearCode, signs::Vector{Int64}=[])
    # this should have Xstabs = Zstabs
    D = dual(C)
    C ⊆ D || error("The single code CSS construxtion requires C ⊆ C^⟂.")
    if !isempty(signs)
        for s in signs
            (s == 1 || s == -1) || error("Stabilizer signs must be +/- 1; received: $s")
        end
    end


    # C2 ⊆ C1
    # k = k1 - k2
    # d >= minimum(d1, d2^⟂)
    # X - H(C2^⟂), Z - H(C1)

    # # should do some checks for equality of the two methods
    # return CSSCode(C, D)
    Sq2 = symplectictoquadratic(directsum(paritycheckmatrix(D), paritycheckmatrix(D)))
    k = dimension(D) - dimension(C)
    if isempty(signs)
        signs = [1 for i in 1:k]
    end
    return CSSCode(field(C), base_ring(Sq2), length(C), k, missing, minimumdistance(D),
        minimumdistance(C), Sq2, paritycheckmatrix(D), paritycheckmatrix(D), C, D,
        signs)
end

function CSSCode(Xmatrix::fq_nmod_mat, Zmatrix::fq_nmod_mat, Xsigns::Vector{Int64}=[],
    Zsigns::Vector{Int64}=[])
    # this should have Zstabs * Xstabs' = 0
    # set dual containing if Xstabs == Zstabs, else not
    size(Xmatrix, 2) ==  size(Zmatrix, 2) || error("Both matrices must have the same length in the CSS construction.")
    base_ring(Xmatrix) == base_ring(Zmatrix) || error("Both matrices must be over the same base field in the CSS construction.")
    rank(Xmatrix) == size(Xmatrix, 1) || error("Provided X-stabilizer matrix is not full rank.")
    rank(Zmatrix) == size(Zmatrix, 1) || error("Provided Z-stabilizer matrix is not full rank.")
    if !isempty(Xsigns)
        for s in Xsigns
            (s == 1 || s == -1) || error("X-stabilizer signs must be +/- 1; received: $s")
        end
    else
        Xsigns = [1 for i in 1:size(Xmatrix, 1)]
    end
    if !isempty(Zsigns)
        for s in Zsigns
            (s == 1 || s == -1) || error("Z-stabilizer signs must be +/- 1; received: $s")
        end
    else
        Zsigns = [1 for i in 1:size(Zmatrix, 1)]
    end

    S = directsum(Xmatrix, Zmatrix)
    for r1 in 1:size(S, 1)
        for r2 in 1:size(S, 1)
            # views?
            iszero(symplecticinnerproduct(S[r1, :], S[r2, :])) || error("The given matrices are not symplectic orthogonal.")
        end
    end
    Sq2 = symplectictoquadratic(S)
    return CSSCode(base_ring(Xmatrix), base_ring(Sq2), size(Xmatrix, 2), size(Xmatrix, 1) + size(Zmatrix, 1),
        missing, missing, missing, Sq2, missing, missing, Xsigns, Zsigns)
end

function CSSCode(SPauli::Vector{T}, signs::Vector{Int64}=[]) where T <: Union{String, Vector{Char}}
    SPaulistripped, signs = _processsigns(SPauli, signs)
    S = _Paulistringtosymplectic(SPaulistripped)
    for r1 in 1:size(S, 1)
        for r2 in 1:size(S, 1)
            iszero(symplecticinnerproduct(S[r1, :], S[r2, :])) || error("The given stabilizers are not symplectic orthogonal.")
        end
    end
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
    # note the H here is transpose of the standard definition
    !(!iszero(G * H) || !iszero(H' * G')) || error("Dual stabilizer matrix is not transpose, symplectic orthogonal.")

    args = isCSS(S, signs, true)
    if args[1]
        return CSSCode(base_ring(S), base_ring(Sq2), size(Sq2, 2), size(Sq2, 1), missing, missing,
            missing, Sq2, args[2], args[4], missing, missing, signs, args[3], args[5])
    else
        error("Provided Pauli strings are not CSS.")
    end
end

# entanglement-assisted is not symplectic orthogonal
function QuantumCode(SPauli::Vector{T}, signs::Vector{Int64}=[]) where T <: Union{String, Vector{Char}}
    SPaulistripped, signs = _processsigns(SPauli, signs)
    S = _Paulistringtosymplectic(SPaulistripped)
    for r1 in 1:size(S, 1)
        for r2 in 1:size(S, 1)
            iszero(symplecticinnerproduct(S[r1, :], S[r2, :])) || error("The given stabilizers are not symplectic orthogonal.")
        end
    end
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
    # note the H here is transpose of the standard definition
    !(!iszero(G * H) || !iszero(H' * G')) || error("Dual stabilizer matrix is not transpose, symplectic orthogonal.")
    Htemp = hcat(H[:, size(Sq2, 2) + 1:end], -H[:, 1:size(Sq2, 2)])
    dualgens = symplectictoquadratic(Htemp')
    dualgens2 = symplectictoquadratic(H')
    # println(dualgens2)
    # println(dualgens)

    # println(size(S))
    # println(size(Htemp))
    # println(size(H))
    #
    # println("testing H'")
    # for r1 in 1:size(S, 1)
    #     for r2 in 1:size(H', 1)
    #         # println("size S:", size(S[r1, :]))
    #         # println("size dualgens:", size(dualgens[r2, :]))
    #         println(iszero(symplecticinnerproduct(S[r1, :], H'[r2, :])))
    #     end
    # end

    # println("testing Htemp")
    # for r1 in 1:size(S, 1)
    #     for r2 in 1:size(Htemp', 1)
    #         # println("size S:", size(S[r1, :]))
    #         # println("size Htemp':", size(Htemp'[r2, :]))
    #         println(iszero(symplecticinnerproduct(S[r1, :], Htemp'[r2, :])))
    #     end
    # end

    ###########

    args = isCSS(S, signs, true)
    if args[1]
        return CSSCode(base_ring(S), base_ring(Sq2), size(Sq2, 2), size(Sq2, 1), missing, missing,
            missing, Sq2, args[2], args[4], missing, missing, signs, args[3], args[5])
    else
        return QuantumCode(base_ring(S), base_ring(Sq2), size(Sq2, 2), size(Sq2, 1), missing,
            Sq2, dualgens, signs)
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
                    elseif j == length(S) && i != dimension(S)
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
                    elseif j == length(S) && i != dimension(S)
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
                    elseif j == length(S) && i != dimension(S)
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
