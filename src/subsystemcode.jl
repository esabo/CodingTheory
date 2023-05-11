# Copyright (c) 2023 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    SubsystemCode(G::fq_nmod_mat, charvec::Union{Vector{nmod}, Missing}=missing)

Return the subsystem code whose gauge group is determined by `G` and signs by `charvec`.
"""
function SubsystemCode(G::fq_nmod_mat, charvec::Union{Vector{nmod}, Missing}=missing)
    iszero(G) && throw(ArgumentError("The gauge matrix is empty."))
    G = _removeempty(G, :rows)

    F = base_ring(G)
    p = Int(characteristic(F))
    n = div(ncols(G), 2)
    charvec = _processcharvec(charvec, p, n)

    # G = < S, gauge ops >
    # C(G) = < S, logs >
    # S = C(G) \cap G
    # C(S) = < S, logs, gauges >
    # logs = C(G) \ S (bare operators)
    # gauge ops = G - S
    # dressed operators = < logs, gauge ops >

    # stabilizer group: ker G ∩ G
    _, kerG = right_kernel(hcat(G[:, n + 1:end], -G[:, 1:n]))
    kerG = transpose(kerG)
    V = VectorSpace(base_ring(G), ncols(kerG))
    kerGVS, kerGtoV = sub(V, [V(kerG[i, :]) for i in 1:nrows(kerG)])
    GVS, _ = sub(V, [V(G[i, :]) for i in 1:nrows(G)])
    I, ItokerG = intersect(kerGVS, GVS)
    if !iszero(AbstractAlgebra.dim(I))
        Ibasis = [kerGtoV(ItokerG(g)) for g in gens(I)]
        Fbasis = [[F(Ibasis[j][i]) for i in 1:AbstractAlgebra.dim(parent(Ibasis[1]))] for j in 1:length(Ibasis)]
        S = matrix(F, length(Fbasis), length(Fbasis[1]), reduce(vcat, Fbasis))
    else
        error("Error computing the stabilizer group of the subsystem code; ker G ∩ G has dimension zero.")
    end

    # check if this subsystem code is really a stabilizer code
    # this should be fine because if they are the same dimension then they are the same
    # don't need to check row space equivalence
    if is_isomorphic(I, GVS)
        println("Stabilizer code detected.")
        # this is duplicating some of the initial work done in this function,
        # but it seems appropriate not to paste the rest of that constructor here
        return StabilizerCode(S, charvec)
    end

    # bare logicals (reps): ker G / S
    BL = _quotientspace(kerG, S)
    graphstate = false
    if !iszero(BL)
        barelogs = _makepairs(BL)
        # verify
        barelogsmat = reduce(vcat, [reduce(vcat, barelogs[i]) for i in 1:length(barelogs)])
        aresymplecticorthogonal(S, barelogsmat) || error("Computed logicals do not commute with the codespace.")
        prod = hcat(barelogsmat[:, n + 1:end], -barelogsmat[:, 1:n]) * transpose(barelogsmat)
        sum(FpmattoJulia(prod), dims=1) == ones(Int, 1, size(prod, 1)) || error("Computed logicals do not have the right commutation relations.")
    else
        graphstate = true
    end
    
    # gauge operators (reps): G / S
    # this can't be zero if S != G
    GO = _quotientspace(G, S)
    gaugeops = _makepairs(GO)
    # verify
    gaugeopsmat = reduce(vcat, [reduce(vcat, gaugeops[i]) for i in 1:length(gaugeops)])
    aresymplecticorthogonal(S, gaugeopsmat) || error("Computed gauge operators do not commute with the codespace.")
    aresymplecticorthogonal(barelogsmat, gaugeopsmat) || error("Computed gauge operators do not commute with the computed logicals.")
    prod = hcat(gaugeopsmat[:, n + 1:end], -gaugeopsmat[:, 1:n]) * transpose(gaugeopsmat)
    sum(FpmattoJulia(prod), dims=1) == ones(Int, 1, size(prod, 1)) || error("Computed gauge operators do not have the right commutation relations.")

    # since S is computed, it automatically has full rank and is not overcomplete
    # q^n / p^k but rows is n - k
    top = BigInt(order(F))^n
    r = length(gaugeops)
    k =  top // BigInt(p)^(nrows(S) + r)
    isinteger(k) && (k = round(Int, log(BigInt(p), k));)
    # r = top // BigInt(p)^length(gaugeops)
    # isinteger(r) && (r = round(Int, log(BigInt(p), r));)

    # determine signs
    signs = _determinesigns(S, charvec)

    args = _isCSSsymplectic(S, signs, true)
    if args[1]
        if graphstate
            return GraphStateSubsystemCSS(F, n, 0, r, missing, missing, missing, S, args[2], args[4],
                missing, missing, signs, args[3], args[4], charvec, missing, false, gaugeops, gaugeopsmat)
        end
        return SubsystemCodeCSS(F, n, k, r, missing, S, args[2], args[4], missing, missing, signs,
            args[3], args[5], barelogs, barelogsmat, charvec, gaugeops, gaugeopsmat, false)
    else
        if graphstate
            return GraphStateSubsystem(F, n, 0, r, missing, S, charvec, signs, missing, false, gaugeops,
                gaugeopsmat)
        end
        return SubsystemCode(F, n, k, r, missing, S, signs, barelogs, barelogsmat, charvec, gaugeops,
            gaugeopsmat, false)
    end
end

"""
    SubsystemCode(GPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}

Return the subsystem code whose gauge group is determined by the vector of Pauli strings `GPauli`
and signs by `charvec`.

# Notes
* Any +/- 1 characters in front of each stabilizer are stripped. No check is done
  to make sure these signs agree with the ones computed using the character vector.
"""
function SubsystemCode(GPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}
    GPaulistripped = _processstrings(GPauli)
    G = _Paulistringtosymplectic(GPaulistripped)
    iszero(G) && error("The processed Pauli strings returned a set of empty gauge group generators.")
    return SubsystemCode(G, charvec)
end

"""
    SubsystemCode(S::fq_nmod_mat, L::CTMatrixTypes, G::CTMatrixTypes, charvec::Union{Vector{nmod}, Missing}=missing)

Return the subsystem code whose stabilizers are given by `S`, (bare) logical operators
by `L`, gauge operators (not including stabilizers) by `G`, and signs by `charvec`.
"""
function SubsystemCode(S::fq_nmod_mat, L::CTMatrixTypes, G::CTMatrixTypes,
    charvec::Union{Vector{nmod}, Missing}=missing)

    # all quantities, including catting all vectors if need be
    iszero(S) && error("The stabilizer matrix is empty.")
    S = _removeempty(S, :rows)
    n = div(ncols(S), 2)
    # check S commutes with itself
    aresymplecticorthogonal(S, S) || error("The given stabilizers are not symplectic orthogonal.")

    # logicals
    iszero(L) && error("The logicals are empty.")
    L = _removeempty(L, :rows)
    # check S commutes with L
    aresymplecticorthogonal(S, L) || error("Logicals do not commute with the code.")
    # check L doesn't commute with itself
    prod = hcat(L[:, n + 1:end], -L[:, 1:n]) * transpose(L)
    iszero(prod) && error("Logicals should not be symplectic self-orthogonal.")
    # check L is properly in pairs
    ncpr = ncols(prod)
    # need an integer sum and this is cheaper than Nemo.ZZ
    prodJul = FpmattoJulia(prod)
    cols = [sum(prodJul[:, i]) for i in 1:ncpr]
    sum(cols) == ncpr || println("Detected logicals not in anticommuting pairs.")
    logpairs = _makepairs(L)
    logsmat = reduce(vcat, [reduce(vcat, logpairs[i]) for i in 1:length(logpairs)])

    # gauge operators
    iszero(G) && error("The gauges are empty.")
    G = _removeempty(G, :rows)
    # check S commutes with G
    aresymplecticorthogonal(S, G) || error("Gauges do not commute with the code.")
    # check L commutes with G
    aresymplecticorthogonal(logsmatrix, G) || error("Gauges do not commute with the logicals.")
    # check G doesn't commute with itself
    prod = hcat(G[:, n + 1:end], -G[:, 1:n]) * transpose(G)
    iszero(prod) && error("Gauges should not be symplectic self-orthogonal.")
    # check G is properly in pairs
    ncpr = ncols(prod)
    # need an integer sum and this is cheaper than Nemo.ZZ
    prodJul = FpmattoJulia(prod)
    cols = [sum(prodJul[:, i]) for i in 1:ncpr]
    sum(cols) == ncpr || println("Detected gauges not in anticommuting pairs.")
    # display(prod)
    gopspairs = _makepairs(G)
    gopsmat = reduce(vcat, [reduce(vcat, gopspairs[i]) for i in 1:length(gopspairs)])

    # display(G)
    F = base_ring(gopsmatrix)
    p = Int(characteristic(F))
    n = div(ncols(G), 2)
    charvec = _processcharvec(charvec, p, n)

    # determine if the provided set of stabilizers are redundant
    rkS = rank(S)
    if nrows(S) > rkS
        overcomp = true
    else
        overcomp = false
    end
    rkG = rank(G)

    top = BigInt(order(F))^n
    # TODO: handle overcompleteness in inputs
    r = div(rkG, 2)
    k =  top // BigInt(p)^(rkS + r)
    isinteger(k) && (k = round(Int, log(BigInt(p), k));)

    # determine signs
    signs = _determinesigns(S, charvec)

    args = _isCSSsymplectic(S, signs, true)
    if args[1]
        return SubsystemCodeCSS(F, n, k, r, missing, S, args[2], args[4], missing, missing, signs,
            args[3], args[5], logpairs, logsmat, charvec, gopspairs, gopsmat, false)
    else
        return SubsystemCode(F, n, k, r, missing, S, signs, logpairs, logsmat, charvec, gopspairs, gopsmat, false)
    end

end

"""
    SubsystemCode(SPauli::Vector{T}, LPauli::Vector{T}, GPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}

Return the subsystem code whose stabilizers are given by the vectors of Pauli strings `SPauli`, (bare)
logical operators by `LPauli`, gauge operators (not including stabilizers) by `GPauli`, and signs by `charvec`.    
"""
# if people want to make a graph code go through the other constructor
function SubsystemCode(SPauli::Vector{T}, LPauli::Vector{T}, GPauli::Vector{T},
    charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}

    S = _Paulistringtosymplectic(_processstrings(SPauli))
    iszero(S) && error("The processed Pauli strings returned a set of empty stabilizer generators.")
    L = _Paulistringtosymplectic(_processstrings(LPauli))
    iszero(L) && error("The processed Pauli strings returned a set of empty logical generators.")
    G = _Paulistringtosymplectic(_processstrings(GPauli))
    iszero(G) && error("The processed Pauli strings returned a set of empty gauge group generators.")
    return SubsystemCode(S, L, G, charvec)
end

# CSS construction, Euclidean and Hermitian
# min dist is min dressed logical operator weight

#############################
      # getter functions
#############################

"""
    field(S::AbstractSubsystemCode)

Return the base ring of the code.
"""
field(S::AbstractSubsystemCode) = S.F

"""
    length(S::AbstractSubsystemCode)
    numqubits(S::AbstractSubsystemCode)

Return the length of the code.
"""
length(S::AbstractSubsystemCode) = S.n
numqubits(S::AbstractSubsystemCode) = S.n

"""
    dimension(S::AbstractSubsystemCode)

Return the dimension of the code.
"""
dimension(S::AbstractSubsystemCode) = S.k

"""
    cardinality(S::AbstractSubsystemCode)

Return the cardinality of the stabilizer group of the code.
"""
cardinality(S::AbstractSubsystemCode) = BigInt(characteristic(S.F))^(S.n - S.k)

"""
    rate(S::AbstractSubsystemCode)

Return the rate, `R = k/n`, of the code.
"""
rate(S::AbstractSubsystemCode) = S.k / S.n

"""
    signs(S::AbstractSubsystemCode)

Return the signs of the stabilizers of the code.
"""
signs(S::AbstractSubsystemCode) = S.signs

"""
    Xsigns(S::AbstractSubsystemCode)

Return the signs of the `X` stabilizers of the CSS code.
"""
Xsigns(S::T) where {T <: AbstractSubsystemCode} = Xsigns(CSSTrait(T), S)
Xsigns(::IsCSS, S::AbstractSubsystemCode) = S.Xsigns
Xsigns(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    Zsigns(S::AbstractSubsystemCode)

Return the signs of the `Z` stabilizers of the CSS code.
"""
Zsigns(S::T) where {T <: AbstractSubsystemCode} = Zsigns(CSSTrait(T), S)
Zsigns(::IsCSS, S::AbstractSubsystemCode) = S.Zsigns
Zsigns(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    stabilizers(S::AbstractSubsystemCode)

Return the stabilizer matrix of the code.
"""
stabilizers(S::AbstractSubsystemCode) = S.stabs

"""
    Xstabilizers(S::AbstractSubsystemCode)

Return the `X`-stabilizer matrix of the CSS code.
"""
Xstabilizers(S::T) where {T <: AbstractSubsystemCode} = Xstabilizers(CSSTrait(T), S)
Xstabilizers(::IsCSS, S::AbstractSubsystemCode) = S.Xstabs
Xstabilizers(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    Zstabilizers(S::AbstractSubsystemCode)

Return the `Z`-stabilizer matrix of the CSS code.
"""
Zstabilizers(S::T) where {T <: AbstractSubsystemCode} = Zstabilizers(CSSTrait(T), S)
Zstabilizers(::IsCSS, S::AbstractSubsystemCode) = S.Zstabs
Zstabilizers(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    numXstabs(S::AbstractSubsystemCode)

Return the number of `X` stabilizers of the CSS code.
"""
numXstabs(S::T) where {T <: AbstractSubsystemCode} = numXstabs(CSSTrait(T), S)
numXstabs(::IsCSS, S::AbstractSubsystemCode) = nrows(S.Xstabs)
numXstabs(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    numZstabs(S::AbstractSubsystemCode)

Return the number of `Z` stabilizers of the CSS code.
"""
numZstabs(S::T) where {T <: AbstractSubsystemCode} = numZstabs(CSSTrait(T), S)
numZstabs(::IsCSS, S::AbstractSubsystemCode) = nrows(S.Zstabs)
numZstabs(::IsNotCSS, S::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    charactervector(S::AbstractSubsystemCode)

Return the character vector of the code.
"""
charactervector(S::AbstractSubsystemCode) = S.charvec

"""
    isovercomplete(S::AbstractSubsystemCode)

Return `true` if `S` has an overcomplete set of stabilizers.
"""
isovercomplete(S::AbstractSubsystemCode) = S.overcomplete

"""
    isCSS(S::AbstractSubsystemCode)

Return `true` is `S` is CSS.
"""
isCSS(S::T) where {T <: AbstractSubsystemCode} = isCSS(CSSTrait(T), S)
isCSS(::IsCSS, S::AbstractSubsystemCode) = true
isCSS(::IsNotCSS, S::AbstractSubsystemCode) = false

"""
    relativedistance(S::AbstractSubsystemCode)

Return the relative minimum distance, `δ = d / n` of the code if `d` is known,
otherwise errors.
"""
function relativedistance(S::AbstractSubsystemCode)
    !ismissing(S.d) || error("Missing minimum distance for this code.")
    return S.d / S.n
end

"""
    logicals(S::AbstractSubsystemCode)
    logicaloperators(S::AbstractSubsystemCode)
    barelogicals(S::AbstractSubsystemCode)
    bare(S::AbstractSubsystemCode)

Return a vector of logical operator generator pairs for `S`.
"""
logicals(S::T) where {T <: AbstractSubsystemCode} = logicals(LogicalTrait(T), S)
logicals(::HasLogicals, S::AbstractSubsystemCode) = S.logicals
logicals(::HasNoLogicals, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no logicals.")
logicaloperators(S::AbstractSubsystemCode) = logicals(S)
barelogicals(S::AbstractSubsystemCode) = logicals(S)
bare(S::AbstractSubsystemCode) = logicals(S)

"""
    logicalsmatrix(S::AbstractSubsystemCode)

Returns the result of `logicals(S)` as a vertically concatenated matrix.
"""
logicalsmatrix(S::T) where {T <: AbstractSubsystemCode} = logicalsmatrix(LogicalTrait(T), S)
logicalsmatrix(::HasLogicals, S::AbstractSubsystemCode) = S.logsmat
logicalsmatrix(::HasNoLogicals, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no logicals.")    

"""
    gauges(S::AbstractSubsystemCode)
    gaugeoperators(S::AbstractSubsystemCode

Return a vector of gauge operator generator pairs for `S`.

# Notes
* Here, gauge operators refers to the gauge group minus the stabilizers.
"""
gauges(S::T) where {T <: AbstractSubsystemCode} = gauges(GaugeTrait(T), S)
gauges(::HasGauges, S::AbstractSubsystemCode) = S.gaugeops
gauges(::HasNoGauges, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no gauges.")
gaugeoperators(S::AbstractSubsystemCode) = gauges(S)

"""
    gaugesmatrix(S::AbstractSubsystemCode)
    gaugeoperatorsmatrix(S::AbstractSubsystemCode)

Return the result of `gauges(S)` as a vertically concatenated matrix.
"""
gaugesmatrix(S::T) where {T <: AbstractSubsystemCode} = gaugesmatrix(GaugeTrait(T), S)
gaugesmatrix(::HasGauges, S::AbstractSubsystemCode) = S.gopsmat
gaugesmatrix(::HasNoGauges, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no gauges.")
gaugeoperatorsmatrix(S::AbstractSubsystemCode) = gaugesmatrix(S)

"""
    dressed(S::AbstractSubsystemCode)
    dressedoperators(S::AbstractSubsystemCode
    dressedlogicals(S::AbstractSubsystemCode)

Return a vector of pairs generators for the dressed operators of `S`.

# Notes
* Here, the dressed operators are the logicals and the gauge operators.
"""
function dresssed(S::T) where {T <: AbstractSubsystemCode}
    if LogicalTrait(T) == HasNoLogicals
        error("Type $T has no logicals.")
    elseif GaugeTrait(S) == HasNoGauges
        error("Type $T has no gauges.")
    end
    return S.logicals ∪ S.gaugeops
end
dressedoperators(S::AbstractSubsystemCode) = dressed(S)
dressedlogicals(S::AbstractSubsystemCode) = dressed(S)

"""
    gaugegroup(S::AbstractSubsystemCode)
    gaugegroupmatrix(S::AbstractSubsystemCode)
    gaugegeneratorsmatrix(S::AbstractSubsystemCode)
    gaugegroupgeneratorsmatrix(S::AbstractSubsystemCode)

Return a matrix giving a (maybe overcomplete) basis for the gauge group.

# Notes
* Here, this is the stabilizers and the gauge operators.
"""
gaugegroup(S::T) where {T <: AbstractSubsystemCode} = gaugegroup(GaugeTrait(T), S)
gaugegroup(::HasGauges, S::AbstractSubsystemCode) = vcat(S.stabs, S.gopsmat)
gaugegroup(::HasNoGauges, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no gauges.")
gaugegroupmatrix(S::AbstractSubsystemCode) = gaugegroup(S)
gaugegeneratorsmatrix(S::AbstractSubsystemCode) = gaugegroup(S)
gaugegroupgeneratorsmatrix(S::AbstractSubsystemCode) = gaugegroup(S)

#############################
      # setter functions
#############################

"""
    changesigns(S::AbstractSubsystemCode, charvec::Vector{nmod})
    changesigns!(S::AbstractSubsystemCode, charvec::Vector{nmod})

Set the character vector of `S` to `charvec` and update the signs.
"""
function changesigns!(S::AbstractSubsystemCode, charvec::Vector{nmod})
    R = base_ring(charactervector(S))
    length(charvec) == 2 * S.n || throw(ArgumentError("Characteristic vector is of improper length for the code."))
    for s in charvec
        modulus(s) == modulus(R) || throw(ArgumentError("Phases are not in the correct ring."))
    end

    S.signs = _getsigns(S.stabilizers, charvec)
    S.charvec = charvec
end

changesigns(S::AbstractSubsystemCode, charvec::Vector{nmod}) = (Snew = deepcopy(S); return changesigns!(Snew, charvec))

"""
    setstabilizers(S::AbstractSubsystemCode, stabs::fq_nmod_mat)
    setstabilizers!(S::AbstractSubsystemCode, stabs::fq_nmod_mat)

Set the stabilizers of `S` to `stabs`.

# Notes
* A check is done to make sure `stabs` are equivalent to the current set of stabilizers.
"""
function setstabilizers!(S::AbstractSubsystemCode, stabs::fq_nmod_mat)
    iszero(stabs) && throw(ArgumentError("The stabilizers cannot be zero."))
    order(S.F) == order(base_ring(stabs)) || throw(ArgumentError("The stabilizers must be over the same field as the code."))

    stabs = _removeempty(stabs, :rows)
    stabs = change_base_ring(S.F, stabs)
    if _hasequivalentrowspaces(S.stabs, stabs)
        S.stabs = stabs
        nrows(stabs) != S.k && (S.overcomplete = true;)
    else
        error("The current stabilizers are not equivalent to the input.")
    end
    S.signs = _determinesigns(stabs, charvec)

    if CSSTrait(typeof(S)) == IsCSS()
        flag, Xstabs, Xsigns, Zstabs, Zsigns = _isCSSsymplectic(stabs, S.signs, true)
        flag || error("Detected equivalent stabilizers but is no longer CSS.")
        S.Xstabs = Xstabs
        S.Xsigns = Xsigns
        S.Zstabs = Zstabs
        S.Zsigns = Zsigns
    end
    return nothing
end

setstabilizers(S::AbstractSubsystemCode, stabs::fq_nmod_mat) = (Snew = deepcopy(S); return setstabilizers!(Snew, stabs))

"""
    setXstabilizers(S::AbstractSubsystemCode, Xstabs::fq_nmod_mat)
    setXstabilizers!(S::AbstractSubsystemCode, Xstabs::fq_nmod_mat)

Set the `X` stabilizers of `S` to `Xstabs`.

# Notes
* A check is done to make sure `stabs` are equivalent to the current set of stabilizers.
"""
setXstabilizers!(S::T, Xstabs::fq_nmod_mat) where {T <: AbstractSubsystemCode} = setXstabilizers!(CSSTrait(T), S, Xstabs)
function setXstabilizers!(::IsCSS, S::AbstractSubsystemCode, Xstabs::fq_nmod_mat)
    iszero(Xstabs) && throw(ArgumentError("The stabilizers cannot be zero."))
    order(S.F) == order(base_ring(Xstabs)) || throw(ArgumentError("The stabilizers must be over the same field as the code."))

    Xstabs = _removeempty(Xstabs, :rows)
    Xstabs = change_base_ring(S.F, Xstabs)
    if _hasequivalentrowspaces(S.Xstabs, Xstabs)
        S.Xstabs = Xstabs
        nrows(Xstabs) != rank(Xstabs) && (S.overcomplete = true;)
    else
        error("The current stabilizers are not equivalent to the input.")
    end
    S.Xsigns = _determinesigns(Xstabs, charvec)
    return nothing
end
setXstabilizers!(::IsNotCSS, S::AbstractSubsystemCode, Xstabs::fq_nmod_mat) = error("X stabilizers are only defined for CSS codes")

setXstabilizers(S::T, Xstabs::fq_nmod_mat) where {T <: AbstractSubsystemCode} = setXstabilizers!(CSSTrait(T), S, Xstabs)
setXstabilizers(::IsCSS, S::AbstractSubsystemCode, Xstabs::fq_nmod_mat) = (Snew = deepcopy(S); return setXstabilizers!(Snew, Xstabs))
setXstabilizers(::IsNotCSS, S::AbstractSubsystemCode, Xstabs::fq_nmod_mat) = error("X stabilizers are only defined for CSS codes")

"""
    setZstabilizers(S::AbstractSubsystemCode, Zstabs::fq_nmod_mat)
    setZstabilizers!(S::AbstractSubsystemCode, Zstabs::fq_nmod_mat)

Set the `Z` stabilizers of `S` to `Zstabs`.

# Notes
* A check is done to make sure `stabs` are equivalent to the current set of stabilizers.
"""
setZstabilizers!(S::T, Zstabs::fq_nmod_mat) where {T <: AbstractSubsystemCode} = setZstabilizers!(CSSTrait(T), S, Zstabs)
function setZstabilizers!(::IsCSS, S::AbstractSubsystemCode, Zstabs::fq_nmod_mat)
    iszero(Zstabs) && throw(ArgumentError("The stabilizers cannot be zero."))
    order(S.F) == order(base_ring(Zstabs)) || throw(ArgumentError("The stabilizers must be over the same field as the code."))

    Zstabs = _removeempty(Zstabs, :rows)
    Zstabs = change_base_ring(S.F, Zstabs)
    if _hasequivalentrowspaces(S.Zstabs, Zstabs)
        S.Zstabs = Zstabs
        nrows(Zstabs) != rank(Zstabs) && (S.overcomplete = true;)
    else
        error("The current stabilizers are not equivalent to the input.")
    end
    S.Zsigns = _determinesigns(Zstabs, charvec)
    return nothing
end
setZstabilizers!(::IsNotCSS, S::AbstractSubsystemCode, Zstabs::fq_nmod_mat) = error("Z stabilizers are only defined for CSS codes")

setZstabilizers(S::T, Zstabs::fq_nmod_mat) where {T <: AbstractSubsystemCode} = setZstabilizers!(CSSTrait(T), S, Zstabs)
setZstabilizers(::IsCSS, S::AbstractSubsystemCode, Zstabs::fq_nmod_mat) = (Snew = deepcopy(S); return setZstabilizers!(Snew, Zstabs))
setZstabilizers(::IsNotCSS, S::AbstractSubsystemCode, Zstabs::fq_nmod_mat) = error("Z stabilizers are only defined for CSS codes")

"""
    setlogicals(S::AbstractSubsystemCode, L::fq_nmod_mat)
    setlogicals!(S::AbstractSubsystemCode, L::fq_nmod_mat)

Set the logical operators of `S` to `L`.

# Notes
* A check is done to make sure `L` are eqivalent to the current set of logicals (up to stabilizers).
"""
setlogicals!(S::T, L::fq_nmod_mat) where {T <: AbstractSubsystemCode} = setlogicals!(LogicalTrait(T), S, L)
function setlogicals!(::HasLogicals, S::AbstractSubsystemCode, L::fq_nmod_mat)
    size(L) == (2 * S.k, 2 * S.n) || throw(ArgumentError("Provided matrix is of incorrect size for the logical space."))
    iseven(ncols(L)) || throw(ArgumentError("Expected a symplectic input but the input matrix has an odd number of columns."))
    S.F == base_ring(L) || throw(ArgumentError("The logicals must be over the same field as the code."))
    
    _hasequivalentrowspaces(vcat(S.logsmat, S.stabs), vcat(L, S.stabs)) || error("The current logicals are not equivalent to the input.")
    # aresymplecticorthogonal(symplecticstabilizers(S), Lsym, true) ||
    #     error("Provided logicals do not commute with the code.")

    # the columns in prod give the commutation relationships between the provided
    # logical operators; they ideally should only consist of {X_1, Z_i} pairs
    # so there should only be one nonzero element in each column
    prod = hcat(L[:, S.n + 1:end], -L[:, 1:S.n]) * transpose(L)
    iszero(prod) && throw(ArgumentError("Provided logicals should not be symplectic self-orthogonal."))
    ncpr = ncols(prod)
    # need an integer sum and this is cheaper than Nemo.ZZ
    prodJul = FpmattoJulia(prod)
    cols = [sum(prodJul[:, i]) for i in 1:ncpr]
    sum(cols) == ncpr || throw(ArgumentError("Incorrect commutation relationships between provided logicals."))

    # pairs are row i, and then whatever column is nonzero, and then shift such that it is one
    F = base_ring(L)
    logs = Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}()
    # this does indeed grow smaller each iteration
    while nrows(L) >= 2
        y = findfirst(x->x>0, prodJul[:, 1])
        y = [F(prod[y, 1]), y]
        if y[1] != F(1)
            push!(logs, (L[1, :], y[1]^-1 * L[y[2], :]))
        else
            push!(logs, (L[1, :], L[y[2], :]))
        end
        L = L[setdiff(1:size(L, 1), [1, y[2]]), :]
    end
    S.logicals = logs
    S.logsmat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
end
setlogicals!(::HasNoLogicals, S::AbstractSubsystemCode, L::fq_nmod_mat) = error("Type $(typeof(S)) has no logicals.")

setlogicals(S::T, L::fq_nmod_mat) where {T <: AbstractSubsystemCode} = setlogicals(LogicalTrait(T), S, L)
setlogicals(::HasLogicals, S::AbstractSubsystemCode, L::fq_nmod_mat) = (Snew = deepcopy(S); return setlogicals!(Snew, L))
setlogicals(::HasNoLogicals, S::AbstractSubsystemCode, L::fq_nmod_mat) = error("Type $(typeof(S)) has no logicals.")

"""
    setminimumdistance(S::AbstractSubsystemCode, d::Int)
Set the minimum distance of the code to `d`.

# Notes
* The only check done on the value of `d` is that `1 ≤ d ≤ n`.
"""
function setminimumdistance!(S::AbstractSubsystemCode, d::Int)
    # TODO: should check bounds like Singleton for possibilities
    d > 0 && d <= S.n || throw(ArgumentError("The minimum distance of a code must be ≥ 1; received: d = $d."))
    S.d = d
end

#############################
     # general functions
#############################

function _processcharvec(charvec::Union{Vector{nmod}, Missing}, p::Int, n::Int)
    if !ismissing(charvec)
        n == length(charvec) || throw(ArgumentError("The characteristic value is of incorrect length."))
        if p == 2
            R = residue_ring(Nemo.ZZ, 4)
        else
            R = residue_ring(Nemo.ZZ, p)
        end
        for s in charvec
            modulus(s) == modulus(R) || throw(ArgumentError("Phases are not in the correct ring."))
        end
    else
        if p == 2
            R = residue_ring(Nemo.ZZ, 4)
        else
            R = residue_ring(Nemo.ZZ, p)
        end
        charvec = [R(0) for _ in 1:n]
    end
    return charvec
end

function _determinesigns(S::fq_nmod_mat, charvec::Vector{nmod})
    if iszero(charvec)
        R = parent(charvec[1])
        signs = [R(0) for _ in 1:nrows(S)]
    else
        signs = _getsigns(S, charvec)
    end
    return signs
end

function _determinesignsCSS(S::fq_nmod_mat, charvec::Vector{nmod}, Xsize::Int, Zsize::Int)
    if iszero(charvec)
        R = parent(charvec[1])
        signs = [R(0) for _ in 1:nrows(S)]
        Xsigns = [R(0) for _ in 1:Xsize]
        Zsigns = [R(0) for _ in 1:Zsize]
    else
        signs = _getsigns(S, charvec)
        Xsigns = signs[1:Xsize, :]
        Zsigns = signs[Xsize + 1:end, :]
    end
    return signs, Xsigns, Zsigns
end

function _getsigns(A::fq_nmod_mat, charvec::Vector{nmod})
    R = base_ring(charvec[1])
    nc = ncols(A)
    length(charvec) == nc || throw(ArgumentError("Input to _getsigns is expected to be in symplectic form and of the same length as the characteristic vector."))
    
    iszero(charvec) && return [R(0) for _ in 1:div(nc, 2)]
    signs = Vector{Int64}()
    for r in 1:nrows(A)
        parity = R(0)
        for c = 1:nc
            !iszero(A[r, c]) && (parity += charvec[c];)
        end
        append!(signs, parity)
    end
    return signs
end

function _splitsymplecticstabilizers(S::fq_nmod_mat, signs::Vector{nmod})
    Xstabs = Vector{fq_nmod_mat}()
    Xsigns = Vector{nmod}()
    Zstabs = Vector{fq_nmod_mat}()
    Zsigns = Vector{nmod}()
    mixedstabs = Vector{fq_nmod_mat}()
    mixedsigns = Vector{nmod}()

    half = div(ncols(S), 2)
    for r in 1:nrows(S)
        # use views?
        s = S[r, :]
        if iszero(s)
            continue
        else
            sx = iszero(s[1, 1:half])
            sz = iszero(s[1, half + 1:end])
            if (sx && !sz)
                push!(Zstabs, s)
                push!(Zsigns, signs[r])
            elseif !sx && sz
                push!(Xstabs, s)
                push!(Xsigns, signs[r])
            elseif !sx && !sz
                push!(mixedstabs, s)
                push!(mixedsigns, signs[r])
            end
        end
    end

    if !isempty(Xstabs)
        Xstabs = reduce(vcat, Xstabs)
    end
    if !isempty(Zstabs)
        Zstabs = reduce(vcat, Zstabs)
    end
    if !isempty(mixedstabs)
        mixedstabs = reduce(vcat, mixedstabs)
    end
    return Xstabs, Xsigns, Zstabs, Zsigns, mixedstabs, mixedsigns
end

"""
    splitstabilizers(S::AbstractSubsystemCode)

Return the set of `X`-only stabilizers and their signs, the set of `Z`-only
stabilizers and their signs, and the remaining stabilizers and their signs.

# Notes
* This function returns six objects of alternating types `fq_nmod_mat` and
`Vector{Int}` for the three sets of stabilizers and signs, respectively.
An empty set of stabilizers is returned as type `Vector{fq_nmod_mat}`.
"""
splitstabilizers(S::AbstractSubsystemCode) = _splitsymplecticstabilizers(symplecticstabilizerS(S), S.signs)

# TODO: rethink how I'm returning all of this and the bottom trim stuff
# probably need to simply redo the above to simply start with zero matrix
# and then either set first or add to it (return matrix[2:end, L])
function _isCSSsymplectic(S::fq_nmod_mat, signs::Vector{nmod}, trim::Bool=true)
    Xstabs, Xsigns, Zstabs, Zsigns, mixedstabs, mixedsigns = _splitsymplecticstabilizers(S, signs)
    if typeof(mixedstabs) <: Vector{fq_nmod_mat}
        if trim
            half = div(ncols(S), 2)
            return true, Xstabs[:, 1:half], Xsigns, Zstabs[:, half + 1:end], Zsigns
        else
            return true, Xstabs, Xsigns, Zstabs, Zsigns
        end
    else
        if trim
            half = div(ncols(S), 2)
            !(typeof(Xstabs) <: Vector{fq_nmod_mat}) && (Xstabs = Xstabs[:, 1:half];)
            !(typeof(Zstabs) <: Vector{fq_nmod_mat}) && (Zstabs = Zstabs[:, half + 1:end];)
            return false, Xstabs, Xsigns, Zstabs[:, half + 1:end], Zsigns, mixedstabs, mixedsigns
        else
            return false, Xstabs, Xsigns, Zstabs, Zsigns, mixedstabs, mixedsigns
        end
    end
end

# using this function for logical and gauge operators
function _makepairs(L::fq_nmod_mat)
    F = base_ring(L)
    n = div(ncols(L), 2)
    logs = Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}()
    # this does indeed grow smaller each iteration
    while nrows(L) >= 2
        # the columns in prod give the commutation relationships between the provided
        # logical operators; they ideally should only consist of {X_1, Z_i} pairs
        # so there should only be one nonzero element in each column
        prod = hcat(L[:, n + 1:end], -L[:, 1:n]) * transpose(L)
        # println("before")
        # display(prod)
        nprod = ncols(prod)
        first = 0
        for c in 1:nprod
            if !iszero(prod[1, c])
                if iszero(first)
                    first = c
                    if !isone(prod[1, c])
                        L[first, :] *= F(prod[1, c]^-1)
                    end
                else
                    L[c, :] += F(prod[1, c]^-1) * L[first, :]
                end
            end
        end
        iszero(first) && error("Cannot make symplectic basis. Often this is due to the fact that the stabilizers are not maximal and therefore the centralizer still containing part of the isotropic subspace.")
        for c in 2:nprod
            if !iszero(prod[first, c])
                L[c, :] += F(prod[first, c]^-1) * L[1, :]
            end
        end
        prod = hcat(L[:, n + 1:end], -L[:, 1:n]) * transpose(L)
        # println("after")
        # display(prod)
        push!(logs, (L[1, :], L[first, :]))
        L = L[setdiff(1:nrows(L), [1, first]), :]
    end
    # display(logs)
    return logs
end

_testlogicalsrelationships(S::T) where {T <: AbstractSubsystemCode} = _testlogicalsrelationships(LogicalTrait(T), S)
function _testlogicalsrelationships(::HasLogicals, S::AbstractSubsystemCode)
    prod = hcat(S.logsmat[:, S.n + 1:end], -S.logsmat[:, 1:S.n]) * transpose(S.logsmat)
    display(prod)
    return nothing
end
_testlogicalsrelationships(::HasNoLogicals, S) = error("Type $(typeof(S)) has no logicals.")

"""
    islogical(S::AbstractSubsystemCode, v::fq_nmod_mat)

Return `true` if the vector `v` is a logical operator for `S`.
"""
islogical(S::T, v::fq_nmod_mat) where {T <: AbstractSubsystemCode} = islogical(LogicalTrait(T), S, v)
function islogical(::HasLogicals, S::AbstractSubsystemCode, v::fq_nmod_mat)
    nc = ncols(S.logsmat)
    aresymplecticorthogonal(S.stabs, v) || return false
    size(v) == (1, nc) && (return !iszero(S.logsmat * transpose(v));)
    size(v) == (nc, 1) && (return !iszero(S.logsmat * v);)
    throw(ArgumentError("Vector to be tested is of incorrect dimension."))
end
islogical(::HasNoLogicals, S::AbstractSubsystemCode, v::fq_nmod_mat) = error("Type $(typeof(S)) has no logicals.")

"""
    syndrome(S::AbstractSubsystemCode, v::fq_nmod_mat)

Return the syndrome of the vector `v` with respect to the stabilizers of `S`.
"""
# TODO: make uniform with approach in function above
function syndrome(S::AbstractSubsystemCode, v::fq_nmod_mat)
    (size(v) != (2 * S.n, 1) && size(v) != (1, 2 * S.n)) &&
        throw(ArgumentError("Vector to be tested is of incorrect dimension; expected length $(2 * n), received: $(size(v))."))
    # base_ring(v) == field(S) || error("Vector must have the same base ring as the stabilizers.")

    nrows(v) != 1 || return S.stabs * transpose(v)
    return S.stabs * v
end

"""
    Xsyndrome(S::AbstractSubsystemCode, v::fq_nmod_mat)

Return the syndrome of the vector `v` with respect to the `X` stabilizers of the
CSS code.
"""
Xsyndrome(S::T, v::fq_nmod_mat) where {T <: AbstractSubsystemCode} = Xsyndrome(CSSTrait(T), S, v)
function Xsyndrome(::IsCSS, S::AbstractSubsystemCode, v::fq_nmod_mat)
    length(v) == 2 * S.n && (v = v[S.n + 1:end];)
    (size(v) != (S.n, 1) && size(v) != (1, S.n)) &&
        error("Vector to be tested is of incorrect dimension; expected length $n, received: $(size(v)).")
    base_ring(v) == S.F || error("Vector must have the same base ring as the stabilizers.")

    nrows(v) != 1 || return S.Xstabs * transpose(v)
    return S.Xstabs * v
end
Xsyndrome(::IsNotCSS, S::AbstractSubsystemCode, v::fq_nmod_mat) = error("Only valid for CSS codes.")

"""
    Zsyndrome(S::AbstractSubsystemCode, v::fq_nmod_mat)

Return the syndrome of the vector `v` with respect to the `Z` stabilizers of the
CSS code.
"""
Zsyndrome(S::T, v::fq_nmod_mat) where {T <: AbstractSubsystemCode} = Zsyndrome(CSSTrait(T), S, v)
function Zsyndrome(::IsCSS, S::AbstractSubsystemCode, v::fq_nmod_mat)
    length(v) == 2 * S.n && (v = v[1:S.n];)
    (size(v) != (n, 1) && size(v) != (1, n)) &&
        error("Vector to be tested is of incorrect dimension; expected length $n, received: $(size(v)).")
    base_ring(v) == S.F || error("Vector must have the same base ring as the stabilizers.")

    nrows(v) != 1 || return S.Zstabs * transpose(v)
    return S.Zstabs * v
end
Zsyndrome(::IsNotCSS, S::AbstractSubsystemCode, v::fq_nmod_mat) = error("Only valid for CSS codes.")

"""
    promotelogicalstogauge(S::AbstractSubsystemCode, pairs::Vector{Int})

Return a new `AbstractSubsystemCode` where the logical pairs in `pairs` are
now considered as gauge operators.
"""
promotelogicalstogauge(S::T, pairs::Vector{Int}) where {T <: AbstractSubsystemCode} = logicalpromotelogicalstogauges(LogicalTrait(T), S, pairs)
function promotelogicalstogauge(::HasLogicals, S::AbstractSubsystemCode, pairs::Vector{Int})
    pairs = sort!(unique!(pairs))
    stabs = S.stabs
    logs = S.logicals
    # will let this error naturally if pairs contains invalid elements
    gaugeops = S.gaugeops ∪ logs[pairs]
    gopsmat = reduce(vcat, [reduce(vcat, gaugeops[i]) for i in 1:length(gaugeops)])
    logs = logs[setdiff![1:S.k, pairs]]
    logsmat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
    # recompute k
    top = BigInt(order(F))^n
    if S.overcomplete
        overcomp = true
        rkS = rank(stabs)
        k =  top // BigInt(p)^(rkS + length(pairs))
    else
        overcomp = false
        k =  top // BigInt(p)^(nrows(stabs) + length(pairs))
    end
    isinteger(k) && (k = round(Int, log(BigInt(p), k));)
    # compute r
    r = top // BigInt(p)^length(pairs)
    isinteger(r) && (r = round(Int, log(BigInt(p), r));)
    return SubsystemCode(S.F, S.n, k, r, S.d, stabs, logs, logsmat, S.charvec, S.signs, gaugeops,
        gopsmat, overcomplete)
end
promotelogicalstogauge(::HasNoLogicals, S::AbstractSubsystemCode, pairs::Vector{Int}) = error("Type $(typeof(S)) has no logicals.")

"""
    swapXZlogicals!(S::AbstractSubsystemCode, pairs::Vector{Int})

Swap the `X` and `Z` logicals specified by `pairs`.    
"""
swapXZlogicals!(S::T, pairs::Vector{Int}) where {T <: AbstractSubsystemCode} = swapXZlogicals!(LogicalTrait(T), S, pairs)
function swapXZlogicals!(::HasLogicals, S::AbstractSubsystemCode, pairs::Vector{Int})
    # let indexing check for inbounds naturally
    pairs = sort!(unique!(pairs))
    temp::fq_nmod_mat
    for i in pairs
        temp = S.logicals[i][1]
        S.logicals[i][1] = S.logicals[i][2]
        S.logicals[i][2] = temp
    end
    S.logsmat = reduce(vcat, [reduce(vcat, S.logicals[i]) for i in 1:length(S.logicals)])
    return nothing
end
swapXZlogicals!(::HasNoLogicals, S::AbstractSubsystemCode, pairs::Vector{Int}) = error("Type $(typeof(S)) has no logicals.")

"""
    swapXZgaugeoperators!(S::AbstractSubsystemCode, pairs::Vector{Int})

Swap the `X` and `Z` gauge operators specified by `pairs`.
"""
swapXZgaugeoperators!(S::T, pairs::Vector{Int}) where {T <: AbstractSubsystemCode} = swapXZgaugeoperators!(GaugeTrait(T), S, pairs)
function swapXZgaugeoperators!(::HasGauges, S::AbstractSubsystemCode, pairs::Vector{Int})
    # let indexing check for inbounds naturally
    pairs = sort!(unique!(pairs))
    temp::fq_nmod_mat
    for i in pairs
        temp = S.gaugeops[i][1]
        S.gaugeops[i][1] = S.gaugeops[i][2]
        S.gaugeops[i][2] = temp
    end
    S.gopsmat = reduce(vcat, [reduce(vcat, S.gaugeops[i]) for i in 1:length(S.gaugeops)])
    return nothing
end
swapXZgaugeoperators!(::HasNoGauges, S::AbstractSubsystemCode, pairs::Vector{Int}) = error("Type $(typeof(S)) has no gauges.")

"""
    areequivalent(S1::T, S2::T) where T <: AbstractSubsystemCode

Return `true` if the codes are equivalent as symplectic vector spaces.

# Note
* This is not intended to detect if `S1` and `S2` are permutation equivalent.
"""
function areequivalent(S1::T, S2::T) where T <: AbstractSubsystemCode
    (S1.n == S2.n && S1.k == S2.k) || return false
    # can't compare fields directly because they are compared based on ptr addresses
    Int(order(S1.F)) == Int(order(S2.F)) || return false
    if GaugeTrait(T) == HasGauges
        S1.r == S2.r || return false
    end

    # # test stabilizers
    _hasequivalentrowspaces(S1.stabs, S2.stabs) || return false

    if LogicalTrait(T) == HasLogicals()
        # test logicals
        _hasequivalentrowspaces(vcat(S1.logsmat, S1.stabs), vcat(S2.logsmat, S2.stabs)) || return false
    end

    if GaugeTrait(T) == HasGauges()
        # test gauge operators
        return _hasequivalentrowspaces(vcat(S1.gopsmat, S1.stabs), vcat(S2.gopsmat, S2.stabs))
    else
        return true
    end
end

"""
    fixgauge(::HasGauges, S::AbstractSubsystemCode, pair::Int, which::Symbol)

Return the code induced by adding either the first of `gaugeoperators(S)[pair]` to the stabilizers
if `which = :X` or the second if `which = :Z`.
"""
fixgauge(S::T, pair::Int, which::Symbol) where {T <: AbstractSubsystemCode} = fixgauge(GaugeTrait(T), S, pair, which)
function fixgauge(::HasGauges, S::AbstractSubsystemCode, pair::Int, which::Symbol)
    if which == :X
        return augment(S, S.gaugeops[pair][1], false)
    elseif which == :Z
        return augment(S, S.gaugeops[pair][2], false)
    else
        throw(ArgumentError("Unknown type $which"))
    end
end
fixgauge(::HasNoGauges, S::AbstractSubsystemCode, pair::Int, which::Symbol) = error("Type $(typeof(S)) has no gauges.")

function show(io::IO, S::AbstractSubsystemCode)
    if isa(S.k, Integer)
        print(io, "[[$(S.n), $(S.k)")
    else
        print(io, "(($(S.n), $(S.k)")
    end
    !(typeof(S) <: AbstractStabilizerCode) && print(io, ", $(S.r)")
    !ismissing(S.d) && print(io, ", $(S.d)")
    if isa(S.k, Integer)
        print(io, "]]_$(order(S.F))")
    else
        print(io, "))_$(order(S.F))")
    end
    if iszero(S.k)
        if isa(S, GraphStateStabilizerCSS)
            println(io, " CSS graph state.")
        else
            println(io, " graph state.")
        end
    else
        if isa(S, StabilizerCodeCSS)
            println(io, " CSS stabilizer code.")
        elseif isa(S, StabilizerCode)
            println(io, " stabilizer code.")
        elseif isa(S, SubsystemCodeCSS)
            println(io, " CSS subsystem code.")
        else
            println(io, " subsystem code.")
        end
    end
    
    if get(io, :compact, true) && S.n <= 30
        if isa(S, SubsystemCodeCSS) || isa(S, StabilizerCodeCSS) || isa(S, GraphStateStabilizerCSS)
            numX = numXstabs(S)
            if S.overcomplete
                println(io, "X-stabilizer matrix (overcomplete): $numX × $(S.n)")
            else
                println(io, "X-stabilizer matrix: $numX × $(S.n)")
            end
            for r in 1:numX
                print(io, "\t chi($(S.Xsigns[r])) ")
                for c in 1:S.n
                    if c != S.n
                        print(io, "$(S.Xstabs[r, c]) ")
                    elseif c == S.n && r != numX
                        println(io, "$(S.Xstabs[r, c])")
                    else
                        print(io, "$(S.Xstabs[r, c])")
                    end
                end
            end
            println(" ")

            numZ = numZstabs(S)
            if S.overcomplete
                println(io, "Z-stabilizer matrix (overcomplete): $numZ × $(S.n)")
            else
                println(io, "Z-stabilizer matrix: $numZ × $(S.n)")
            end
            for r in 1:numZ
                print(io, "\t chi($(S.Zsigns[r])) ")
                for c in 1:S.n
                    if c != S.n
                        print(io, "$(S.Zstabs[r, c]) ")
                    elseif c == S.n && r != numZ
                        println(io, "$(S.Zstabs[r, c])")
                    else
                        print(io, "$(S.Zstabs[r, c])")
                    end
                end
            end
        else
            numstabs = nrows(S.stabs)
            if S.overcomplete
                println(io, "Stabilizer matrix (overcomplete): $numstabs × $(S.n)")
            else
                println(io, "Stabilizer matrix: $numstabs × $(S.n)")
            end
            for r in 1:numstabs
                print(io, "\t chi($(S.signs[r])) ")
                for c in 1:S.n
                    if c != S.n
                        print(io, "$(S.stabs[r, c]) ")
                    elseif c == S.n && r != numstabs
                        println(io, "$(S.stabs[r, c])")
                    else
                        print(io, "$(S.stabs[r, c])")
                    end
                end
            end
        end
    end
end

# TODO: should probably split these in two for type stability
# TODO: iterate over this faster
# TODO: remove quadratic
function _allstabilizers(S::AbstractStabilizerCode, onlyprint::Bool=false)
    E = quadraticfield(S)
    all = Vector{fq_nmod_mat}()
    stabs = S.stabs
    for iter in Base.Iterators.product([0:(Int64(characteristic(S.F)) - 1) for _ in 1:nrows(stabs)]...)
        stab = E(iter[1]) * stabs[1, :]
        for r in 2:nrows(stabs)
            if !iszero(iter[r])
                stab += E(iter[r]) * stabs[r, :]
            end
        end
        if onlyprint
            println(stab)
        else
            push!(all, stab)
        end
    end
    if onlyprint
        return
    else
        return all
    end
end

"""
    allstabilizers(S::AbstractSubsystemCode)
    elements(S::AbstractSubsystemCode)

Return a vector of all the elements of the stabilizer group of `S`.
"""
allstabilizers(S::AbstractSubsystemCode) = _allstabilizers(S, false)
elements(S::AbstractSubsystemCode) = allstabilizers(S)

"""
    printallstabilizers(S::AbstractSubsystemCode)
    printallelements(S::AbstractSubsystemCode)

Print all the elements of the stabilizer group of `S`.
"""
printallstabilizers(S::AbstractSubsystemCode) = _allstabilizers(S, true)
printallelements(S::AbstractSubsystemCode) = printallstabilizers(S)


# function Singletonbound()
#     n - k ≧ 2(d - 1) for stabilizer
#     k + r ≤ n - 2d + 2 for subsystem
# they should be the same but with r = 0 for stabilizer
# end

"""
    permutecode(S::AbstractSubsystemCode, σ::Union{Perm{T}, Vector{T}}) where T <: Int
    permutecode!(S::AbstractSubsystemCode, σ::Union{Perm{T}, Vector{T}}) where T <: Int

Return the code permuted by `σ`.

# Notes
* If `σ` is a vector, it is interpreted as the desired column order for the generator matrix of `C`.
"""
function permutecode!(S::AbstractSubsystemCode, σ::Union{Perm{T}, Vector{T}}) where T <: Int
    if typeof(σ) <: Perm
        perm1 = matrix(S.F, Array(matrix_repr(σ)))
        perm = perm1 ⊕ perm1
        S.stabs = S.stabs * perm
        S.charvec = S.charvec * perm

        W = typeof(S)
        if CSSTrait(W) == IsCSS()
            S.Xstabs = S.Xstabs * perm1
            S.Zstabs = S.Zstabs * perm1
        end

        if LogicalTrait(W) == HasLogicals()
            S.logsmat = S.logsmat * perm
            for (i, pair) in enumerate(S.logicals)
                S.logicals[i] = (pair[1] * perm, pair[2] * perm)
            end
        end

        if GaugeTrait(W) == HasGauges()
            S.gopsmat = S.gopsmat * perm
            for (i, pair) in enumerate(S.gaugeops)
                S.gaugeops[i] = (pair[1] * perm, pair[2] * perm)
            end
        end
    else
        length(unique(σ)) == S.n || throw(ArgumentError("Incorrect number of digits in permutation."))
        (1 == minimum(σ) && S.n == maximum(σ)) || throw(ArgumentError("Digits are not in the range `1:$(S.n)`."))

        σ = σ ∪ (σ .+ S.n)
        S.stabs = S.stabs[:, σ]
        S.charvec = S.charvec[:, σ]
        W = typeof(S)
        if CSSTrait(W) == IsCSS()
            S.Xstabs = S.Xstabs[:, σ]
            S.Zstabs = S.Zstabs[:, σ]
        end

        if LogicalTrait(W) == HasLogicals()
            S.logsmat = S.logsmat[:, σ]
            for (i, pair) in enumerate(S.logicals)
                S.logicals[i] = (pair[1][:, σ], pair[2][:, σ])
            end
        end

        if GaugeTrait(W) == HasGauges()
            S.gopsmat = S.gopsmat[:, σ]
            for (i, pair) in enumerate(S.gaugeops)
                S.gaugeops[i] = (pair[1][:, σ], pair[2][:, σ])
            end
        end
    end
    return S
end

permutecode(S::AbstractSubsystemCode, σ::Union{Perm{T}, Vector{T}}) where T <: Int = (Snew = deepcopy(S); return permutecode!(Snew, σ);)

"""
    augment(S::AbstractSubsystemCode, row::fq_nmod_mat, verbose::Bool=true)

Return the code created by added `row` to the stabilizers of `S`.

# Notes
* The goal of this function is to track how the logical operators update given the new stabilizer.
  The unaffected logical operators are kept during the update and only those which don't commute
  with the new stabilizer are recomputed. Use `verbose` to better 
"""
function augment(S::AbstractSubsystemCode, row::fq_nmod_mat, verbose::Bool=true)
    iszero(row) && return S
    nrows(row) == 1 || throw(ArgumentError("Only one stabilizer may be passed in at a time."))

    # stabilizers
    prod = hcat(S.stabs[:, S.n + 1:end], -S.stabs[:, 1:S.n]) * transpose(row)
    if iszero(prod)
        verbose && println("Vector is already in the stabilizer group. Nothing to update.")    
        Snew = deepcopy(S)
        Snew.stabs = vcat(S.stabs, row)
        Snew.overcomplete = true
        return Snew
    else
        stabstokeep = Vector{Int}()
        for i in 1:nrows(S.stabs)
            iszero(prod[i]) && append!(stabstokeep, i, i + 1)
        end

        if isempty(stabstokeep)
            verbose && println("The vector anticommutes with all stabilizers. The new stabilizer group is just the vector.")
            stabs = row
        else
            update = setdiff(1:nrows(S.stabs), stabstokeep)
            if verbose
                if isempty(update)
                    println("No stabilizers requiring updating")
                else
                    println("Stabilizers requiring updating:")
                    display(update)
                end
            end
            isempty(update) ? (stabs = S.stabs;) : (stabs = S.stabs[stabstokeep, :];)
        end
    end

    # logicals
    if LogicalTrait(typeof(S)) == HasLogicals()
        prod = hcat(S.logsmat[:, S.n + 1:end], -S.logsmat[:, 1:S.n]) * transpose(row)
        logstokeep = Vector{Int}()
        logpairstokeep = Vector{Int}()
        pair = 1
        for i in 1:2:nrows(S.logsmat)
            # this is predicated on the idea that the pairs are stacked together in this matrix
            # TODO: write a unit test checking this never changes
            if iszero(prod[i]) && iszero(prod[i + 1])
                append!(logstokeep, i, i + 1)
                append!(logpairstokeep, pair)
            end
            pair += 1
        end

        if isempty(logstokeep)
            verbose && println("The vector anticommutes with all logical pairs.")
            logs = zero_matrix(S.F, 1, 2 * S.n)
        else
            update = setdiff(1:length(S.logicals), logpairstokeep)
            if verbose
                if isempty(update)
                    println("No logical pairs requiring updating")
                else
                    println("Logical pairs requiring updating:")
                    display(update)
                end
            end
            isempty(update) ? (logs = S.logsmat;) : (logs = S.logsmat[logstokeep, :];)
        end
    else
        logs = zero_matrix(S.F, 1, 2 * S.n)
    end

    # gauges
    if GaugeTrait(typeof(S)) == HasGauges()
        prod = hcat(S.gopsmat[:, S.n + 1:end], -S.gopsmat[:, 1:S.n]) * transpose(row)
        gopstokeep = Vector{Int}()
        goppairstokeep = Vector{Int}()
        pair = 1
        for i in 1:2:nrows(S.gopsmat)
            # this is predicated on the idea that the pairs are stacked together in this matrix
            # TODO: write a unit test checking this never changes
            if iszero(prod[i]) && iszero(prod[i + 1])
                append!(gopstokeep, i, i + 1)
                append!(goppairstokeep, pair)
            end
            pair += 1
        end

        if isempty(gopstokeep)
            verbose && println("The vector anticommutes with all gauge operator pairs.")
            gaugeops = zero_matrix(S.F, 1, 2 * S.n)
        else
            update = setdiff(1:length(S.gaugeops), goppairstokeep)
            if verbose
                if isempty(update)
                    println("No gauge operator pairs requiring updating")
                else
                    println("Gauge operator pairs requiring updating:")
                    display(update)
                end
            end
            isempty(update) ? (gaugeops = S.gaugeops;) : (gaugeops = S.gaugeops[gopstokeep, :];)
        end
    else
        gaugeops = zero_matrix(S.F, 1, 2 * S.n)
    end

    # compute newly opened degrees of freedom
    temp = _removeempty(vcat(stabs, logs, gaugeops), :rows)
    _, temp = right_kernel(hcat(temp[:, S.n + 1:end], -temp[:, 1:S.n]))
    temp = _quotientspace(transpose(temp), newsymstabs)
    newlogs = _makepairs(temp)
    return SubsystemCode(stabs, vcat(logs, newlogs), gaugeops, S.charvec)
end

"""
    expurgate(S::AbstractStabilizerCode, rows::Vector{Int}, verbose::Bool=true)

Return the code created by removing the stabilizers indexed by `rows`.

# Notes
* The goal of this function is to track how the logical operators update through this process.
  Here, the original logical pairs are kept and an appropriate number of new pairs are added.
"""
function expurgate(S::AbstractSubsystemCode, rows::Vector{Int}, verbose::Bool=true)
    numstabs = nrows(S.stabs)
    rows ⊆ 1:numstabs || throw(ArgumentError("Argument `rows` not a subset of the number of stabilizers."))
    rows == 1:numstabs && throw(ArgumentError("Cannot remove all stabilizers"))

    verbose && println("Removing stabilizers: $rows")
    newstabs = S.stabs[setdiff(1:numstabs, rows), :]
    temp = newstabs
    if LogicalTrait(typeof(S)) == HasLogicals()
        temp = vcat(temp, S.logsmat)
    end
    if GaugeTrait(typeof(S)) == HasGauges()
        temp = vcat(temp, S.gopsmat)
    end
    _, H = right_kernel(hcat(temp[:, S.n + 1:end], -temp[:, 1:S.n]))
    H = transpose(H)
    newlogs = _quotientspace(H, newstabs)
    if iszero(newlogs)
        verbose && println("No new logicals need to be add")
        Snew = deepcopy(S)
        Snew.stabs = newstabs
        rank(newstabs) == nrows(newstabs) ? (S.overcomplete = false;) : (S.overcomplete = true;)
        return Snew
    else
        newlogpairs = _makepairs(newlogs)
        verbose && println("New logical pairs:")
        verbose && display(newlogpairs)
        if GaugeTrait(typeof(S)) == HasGauges()
            return SubsystemCode(newstabs, vcat(S.logsmat, reduce(vcat, reduce(vcat, newlogpairs))), S.gopsmat, S.charvec)
        else
            Snew = StabilizerCode(newstabs, S.charvec)
            setlogicals!(Snew, vcat(S.logsmat, reduce(vcat, reduce(vcat, newlogpairs))))
            return Snew
        end
    end
end
