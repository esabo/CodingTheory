# Copyright (c) 2023 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

function SubsystemCode(Gq2::fq_nmod_mat, symp::Bool=false, charvec::Union{Vector{nmod}, Missing}=missing)
    iszero(Gq2) && error("The gauge matrix is empty.")
    Gq2 = _removeempty(Gq2, "rows")
    if symp
        G = Gq2
        # this will error properly if not correct
        Gq2 = symplectictoquadratic(Gq2)
    else
        iseven(degree(base_ring(Gq2))) || error("The base ring of the given matrix is not a quadratic extension.")
        G = quadratictosymplectic(Gq2)
    end

    F = base_ring(G)
    p = Int(characteristic(F))
    n = ncols(Gq2)
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
        S = matrix(F, length(Fbasis), length(Fbasis[1]), vcat(Fbasis...))
        Sq2 = symplectictoquadratic(S)
        # display(Sq2)
    else
        error("Error computing the stabilizer group of the subsystem code; ker G ∩ G has dimension zero.")
    end

    # check if this subsystem code is really a stabilizer code
    if is_isomorphic(I, GVS)
        println("Stabilizer code detected.")
        # this is duplicating some of the initial work done in this function,
        # but it seems appropriate not to paste the rest of that constructor here
        return StabilizerCode(S, true, charvec)
    end

    # bare logicals (reps): ker G / S
    # _, kerS = right_kernel(hcat(S[:, n + 1:end], -S[:, 1:n]))
    # kerS = transpose(kerS)
    # display(kerS)
    # println(" ")
    # display(G)
    BL = _quotientspace(kerG, S)
    graphstate = false
    if !iszero(BL)
        BLq2 = symplectictoquadratic(BL)
        # display(BLq2)
        barelogs = _makepairs(BLq2)
        # verify
        barelogsmat = vcat([vcat(barelogs[i]...) for i in 1:length(barelogs)]...)
        aresymplecticorthogonal(Sq2, barelogsmat) || error("Computed logicals do not commute with the codespace.")
        BLsymmat = quadratictosymplectic(barelogsmat)
        prod = hcat(BLsymmat[:, n + 1:end], -BLsymmat[:, 1:n]) * transpose(BLsymmat)
        sum(FpmattoJulia(prod), dims=1) == ones(Int, 1, size(prod, 1)) || error("Computed logicals do not have the right commutation relations.")
    else
        graphstate = true
    end
    
    # gauge operators (reps): G / S
    # this can't be zero if S != G
    GO = _quotientspace(G, S)
    GOq2 = symplectictoquadratic(GO)
    gaugeops = _makepairs(GOq2)
    # verify
    gaugeopsmat = vcat([vcat(gaugeops[i]...) for i in 1:length(gaugeops)]...)
    aresymplecticorthogonal(Sq2, gaugeopsmat) || error("Computed gauge operators do not commute with the codespace.")
    aresymplecticorthogonal(barelogsmat, gaugeopsmat) || error("Computed gauge operators do not commute with the computed logicals.")
    GOsymmat = quadratictosymplectic(gaugeopsmat)
    prod = hcat(GOsymmat[:, n + 1:end], -GOsymmat[:, 1:n]) * transpose(GOsymmat)
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
            return GraphStateSubsystemCSS(F, base_ring(Sq2), n, 0, r, missing, missing, missing, Sq2,
                args[2], args[4], missing, missing, signs, args[3], args[4], charvec, missing, false, gaugeops)
        end
        return SubsystemCodeCSS(F, base_ring(Sq2), n, k, r, missing, Sq2, args[2], args[4],
            missing, missing, signs, args[3], args[5], barelogs, charvec, gaugeops, false)
    else
        if graphstate
            return GraphStateSubsystem(F, base_ring(Sq2), n, 0, r, missing, Sq2, charvec, signs,
                missing, false, gaugeops)
        end
        return SubsystemCode(F, base_ring(Sq2), n, k, r, missing, Sq2, signs, barelogs, charvec, gaugeops, false)
    end
end

function SubsystemCode(GPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}
    GPaulistripped, charvec = _processstrings(GPauli, charvec)
    G = _Paulistringtosymplectic(GPaulistripped)
    iszero(G) && error("The processed Pauli strings returned a set of empty gauge group generators.")
    return SubsystemCode(G, true, charvec)
end

function SubsystemCode(Sq2::fq_nmod_mat, Lq2::Union{fq_nmod_mat, Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}},
    Gq2::Union{fq_nmod_mat, Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}}, symp::Bool=false, charvec::Union{Vector{nmod}, Missing}=missing)

    # start by getting both symplectic and quadratic representations of
    # all quantities, including catting all vectors if need be
    iszero(Sq2) && error("The stabilizer matrix is empty.")
    Sq2 = _removeempty(Sq2, "rows")
    n = ncols(Sq2)
    if symp
        S = Sq2
        n = div(n, 2)
        # this will error properly if not correct
        Sq2 = symplectictoquadratic(Sq2)
    else
        iseven(degree(base_ring(Sq2))) || error("The base ring of the given matrix is not a quadratic extension.")
        S = quadratictosymplectic(Sq2)
    end
    # check S commutes with itself
    aresymplecticorthogonal(S, S, true) || error("The given stabilizers are not symplectic orthogonal.")

    # logicals
    if typeof(Lq2) <: Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}
        if symp
            # in this case I have the log pairings in symplectic form
            logsmatrix = vcat([vcat(Lq2[i]...) for i in 1:length(Lq2)]...)
            iszero(logsmatrix) && error("The logicals are empty.")
            logsmatrix = _removeempty(logsmatrix, "rows")
            # check S commutes with L
            aresymplecticorthogonal(S, logsmatrix, true) || error("Logicals do not commute with the code.")
            # check L doesn't commute with itself
            prod = hcat(logsmatrix[:, n + 1:end], -logsmatrix[:, 1:n]) * transpose(logsmatrix)
            iszero(prod) && error("Logicals should not be symplectic self-orthogonal.")
            # check L is properly in pairs
            ncpr = ncols(prod)
            # need an integer sum and this is cheaper than Nemo.ZZ
            prodJul = FpmattoJulia(prod)
            cols = [sum(prodJul[:, i]) for i in 1:ncpr]
            if sum(cols) != ncpr
                println("Detected logicals not in anticommuting pairs.")
                logpairsq2 = _makepairs(symplectictoquadratic(logsmatrix))
            else
                # already in proper pairs, just get quadratic versions
                logpairsq2::Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}
                for pair in Lq2
                    push!(logpairsq2, (symplectictoquadratic(pair[1]), symplectictoquadratic(pair[2])))
                end
            end
        else
            # in this case I have the log pairings in quadratic form
            logsmatrixq2 = vcat([vcat(Lq2[i]...) for i in 1:length(Lq2)]...)
            iszero(logsmatrixq2) && error("The logicals are empty.")
            logsmatrixq2 = _removeempty(logsmatrixq2, "rows")
            logsmatrix = symplectictoquadratic(logsmatrixq2)
            # check S commutes with L
            aresymplecticorthogonal(S, logsmatrix, true) || error("Logicals do not commute with the code.")
            # check L doesn't commute with itself
            prod = hcat(logsmatrix[:, n + 1:end], -logsmatrix[:, 1:n]) * transpose(logsmatrix)
            iszero(prod) && error("Logicals should not be symplectic self-orthogonal.")
            # check L is properly in pairs
            ncpr = ncols(prod)
            # need an integer sum and this is cheaper than Nemo.ZZ
            prodJul = FpmattoJulia(prod)
            cols = [sum(prodJul[:, i]) for i in 1:ncpr]
            if sum(cols) != ncpr
                println("Detected logicals not in anticommuting pairs.")
                logpairsq2 = _makepairs(logsmatrixq2)
            else
                # the pairs that were passed in were good
                logpairsq2 = Lq2
            end
        end
    else
        iszero(Lq2) && error("The logicals are empty.")
        Lq2 = _removeempty(Lq2, "rows")
        # in this case logs were passed in as a matrix
        if symp
            logsmatrix = Lq2
            # check S commutes with L
            aresymplecticorthogonal(S, logsmatrix, true) || error("Logicals do not commute with the code.")
            # check L doesn't commute with itself
            prod = hcat(logsmatrix[:, n + 1:end], -logsmatrix[:, 1:n]) * transpose(logsmatrix)
            iszero(prod) && error("Logicals should not be symplectic self-orthogonal.")
            # check L is properly in pairs
            ncpr = ncols(prod)
            # need an integer sum and this is cheaper than Nemo.ZZ
            prodJul = FpmattoJulia(prod)
            cols = [sum(prodJul[:, i]) for i in 1:ncpr]
            sum(cols) == ncpr || println("Detected logicals not in anticommuting pairs.")
            logpairsq2 = _makepairs(symplectictoquadratic(logsmatrix))
        else
            logsmatrix = quadratictosymplectic(Lq2)
            # check S commutes with L
            aresymplecticorthogonal(S, logsmatrix, true) || error("Logicals do not commute with the code.")
            # check L doesn't commute with itself
            prod = hcat(logsmatrix[:, n + 1:end], -logsmatrix[:, 1:n]) * transpose(logsmatrix)
            iszero(prod) && error("Logicals should not be symplectic self-orthogonal.")
            # check L is properly in pairs
            ncpr = ncols(prod)
            # need an integer sum and this is cheaper than Nemo.ZZ
            prodJul = FpmattoJulia(prod)
            cols = [sum(prodJul[:, i]) for i in 1:ncpr]
            sum(cols) == ncpr || println("Detected logicals not in anticommuting pairs.")
            logpairsq2 = _makepairs(symplectictoquadratic(logsmatrix))
        end
    end

    # gauge operators
    if typeof(Gq2) <: Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}
        if symp
            # in this case I have the gauge pairings in symplectic form
            gopsmatrix = vcat([vcat(Gq2[i]...) for i in 1:length(Gq2)]...)
            iszero(gopsmatrix) && error("The gauges are empty.")
            gopsmatrix = _removeempty(gopsmatrix, "rows")
            # check S commutes with G
            aresymplecticorthogonal(S, gopsmatrix, true) || error("Gauges do not commute with the code.")
            # check L commutes with G
            aresymplecticorthogonal(logsmatrix, gopsmatrix, true) || error("Gauges do not commute with the logicals.")
            # check G doesn't commute with itself
            prod = hcat(gopsmatrix[:, n + 1:end], -gopsmatrix[:, 1:n]) * transpose(gopsmatrix)
            iszero(prod) && error("Gauges should not be symplectic self-orthogonal.")
            # check G is properly in pairs
            ncpr = ncols(prod)
            # need an integer sum and this is cheaper than Nemo.ZZ
            prodJul = FpmattoJulia(prod)
            cols = [sum(prodJul[:, i]) for i in 1:ncpr]
            if sum(cols) != ncpr
                println("Detected gauges not in anticommuting pairs.")
                gopspairsq2 = _makepairs(symplectictoquadratic(gopsmatrix))
            else
                # already in proper pairs, just get quadratic versions
                gopspairsq2::Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}
                for pair in Gq2
                    push!(gopspairsq2, (symplectictoquadratic(pair[1]), symplectictoquadratic(pair[2])))
                end
            end
        else
            # in this case I have the log pairings in quadratic form
            gopsmatrixq2 = vcat([vcat(Fq2[i]...) for i in 1:length(Gq2)]...)
            iszero(gopsmatrixq2) && error("The gauges are empty.")
            gopsmatrixq2 = _removeempty(gopsmatrixq2, "rows")
            gopsmatrix = symplectictoquadratic(gopsmatrixq2)
            # check S commutes with G
            aresymplecticorthogonal(S, gopsmatrix, true) || error("Gauges do not commute with the code.")
            # check L commutes with G
            aresymplecticorthogonal(logsmatrix, gopsmatrix, true) || error("Gauges do not commute with the logicals.")
            # check G doesn't commute with itself
            prod = hcat(gopsmatrix[:, n + 1:end], -gopsmatrix[:, 1:n]) * transpose(gopsmatrix)
            iszero(prod) && error("Gauges should not be symplectic self-orthogonal.")
            # check G is properly in pairs
            ncpr = ncols(prod)
            # need an integer sum and this is cheaper than Nemo.ZZ
            prodJul = FpmattoJulia(prod)
            cols = [sum(prodJul[:, i]) for i in 1:ncpr]
            if sum(cols) != ncpr
                println("Detected gauges not in anticommuting pairs.")
                gopspairsq2 = _makepairs(symplectictoquadratic(gopsmatrix))
            else
                # the pairs that were passed in were good
                gopspairsq2 = Gq2
            end
        end
    else
        iszero(Gq2) && error("The gauges are empty.")
        Gq2 = _removeempty(Gq2, "rows")
        # in this case gauges were passed in as a matrix
        if symp
            gopsmatrix = Gq2
            # check S commutes with G
            aresymplecticorthogonal(S, gopsmatrix, true) || error("Gauges do not commute with the code.")
            # check L commutes with G
            aresymplecticorthogonal(logsmatrix, gopsmatrix, true) || error("Gauges do not commute with the logicals.")
            # check G doesn't commute with itself
            prod = hcat(gopsmatrix[:, n + 1:end], -gopsmatrix[:, 1:n]) * transpose(gopsmatrix)
            iszero(prod) && error("Gauges should not be symplectic self-orthogonal.")
            # check G is properly in pairs
            ncpr = ncols(prod)
            # need an integer sum and this is cheaper than Nemo.ZZ
            prodJul = FpmattoJulia(prod)
            cols = [sum(prodJul[:, i]) for i in 1:ncpr]
            sum(cols) == ncpr || println("Detected gauges not in anticommuting pairs.")
            # display(prod)
            gopspairsq2 = _makepairs(symplectictoquadratic(gopsmatrix))
        else
            gopsmatrix = quadratictosymplectic(Gq2)
            # check S commutes with G
            aresymplecticorthogonal(S, gopsmatrix, true) || error("Gauges do not commute with the code.")
            # check L commutes with G
            aresymplecticorthogonal(logsmatrix, gopsmatrix, true) || error("Gauges do not commute with the logicals.")
            # check G doesn't commute with itself
            prod = hcat(gopsmatrix[:, n + 1:end], -gopsmatrix[:, 1:n]) * transpose(gopsmatrix)
            iszero(prod) && error("Gauges should not be symplectic self-orthogonal.")
            # check G is properly in pairs
            ncpr = ncols(prod)
            # need an integer sum and this is cheaper than Nemo.ZZ
            prodJul = FpmattoJulia(prod)
            cols = [sum(prodJul[:, i]) for i in 1:ncpr]
            sum(cols) == ncpr || println("Detected gauges not in anticommuting pairs.")
            gopspairsq2 = _makepairs(symplectictoquadratic(gopsmatrix))
        end
    end

    # display(G)
    F = base_ring(gopsmatrix)
    p = Int(characteristic(F))
    # display(Gq2)
    symp ? (n = div(ncols(Gq2), 2);) : (n = ncols(Gq2);)
    # println(n)
    charvec = _processcharvec(charvec, p, n)

    # determine if the provided set of stabilizers are redundant
    rkS = rank(S)
    if nrows(S) > rkS
        overcomp = true
    else
        overcomp = false
    end
    rkG = rank(gopsmatrix)

    top = BigInt(order(F))^n
    # TODO: handle overcompleteness in inputs
    r = div(rkG, 2)
    k =  top // BigInt(p)^(rkS + r)
    isinteger(k) && (k = round(Int, log(BigInt(p), k));)

    # determine signs
    signs = _determinesigns(S, charvec)

    args = _isCSSsymplectic(S, signs, true)
    if args[1]
        return SubsystemCodeCSS(F, base_ring(Sq2), n, k, r, missing, Sq2, args[2], args[4],
            missing, missing, signs, args[3], args[5], logpairsq2, charvec, gopspairsq2, false)
    else
        return SubsystemCode(F, base_ring(Sq2), n, k, r, missing, Sq2, signs, logpairsq2, charvec, gopspairsq2, false)
    end

end

# if people want to make a graph code go through the other constructor
function SubsystemCode(SPauli::Vector{T}, LPauli::Vector{T}, GPauli::Vector{T},
    charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}

    # TODO: look into stripping the charvec processing from this so not to repeat
    SPaulistripped, charvec = _processstrings(SPauli, charvec)
    LPaulistripped, _ = _processstrings(LPauli, charvec)
    GPaulistripped, _ = _processstrings(GPauli, charvec)

    S = _Paulistringtosymplectic(SPaulistripped)
    iszero(S) && error("The processed Pauli strings returned a set of empty stabilizer generators.")
    L = _Paulistringtosymplectic(LPaulistripped)
    iszero(L) && error("The processed Pauli strings returned a set of empty logical generators.")
    G = _Paulistringtosymplectic(GPaulistripped)
    iszero(G) && error("The processed Pauli strings returned a set of empty gauge group generators.")
    return SubsystemCode(S, L, G, true, charvec)
end

# CSS construction, Euclidean and Hermitian
# min dist is min dressed logical operator weight

# common codes to implement: GeneralizedBaconShor, CompassCode

#############################
      # getter functions
#############################

"""
    field(S::AbstractSubsystemCode)

Return the base ring of the code as a Nemo object.
"""
field(S::AbstractSubsystemCode) = S.F

"""
    quadraticfield(S::AbstractSubsystemCode)

Return the quadratic field of the code as a Nemo object.
"""
quadraticfield(S::AbstractSubsystemCode) = S.E

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

No size checking is done on the parameters of the code, returns a BitInt by default.
"""
cardinality(S::AbstractSubsystemCode) = BigInt(characteristic(S.F))^(S.n - S.k)

"""
    rate(S::AbstractSubsystemCode)

Return the rate, `R = k/n', of the code.
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
    symplecticstabilizers(S::AbstractSubsystemCode)

Return the stabilizer matrix of the code in symplectic form.
"""
symplecticstabilizers(S::AbstractSubsystemCode) = quadratictosymplectic(S.stabs)

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

# Notes
* This is intended to be a simple function wrapper for `typeof(S)` since the
 constructor for `SubsystemCode` automatically returns a `SubsystemCodeCSS` if possible.
 Manually changing the elements of the struct `S` without using the helper
 functions provided here is therefore not recommended.
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
logicalsmatrix(::HasLogicals, S::AbstractSubsystemCode) = vcat([vcat(S.logicals[i]...) for i in 1:S.k]...)
logicalsmatrix(::HasNoLogicals, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no logicals.")    

"""
    gauges(S::AbstractSubsystemCode)
    gaugeoperators(S::AbstractSubsystemCode

Return a vector of gauge operator generator pairs for `S`.

Here, gauge operators refers to the gauge group minus the stabilizers.
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
gaugesmatrix(::HasGauges, S::AbstractSubsystemCode) =
    vcat([vcat(S.gaugeops[i]...) for i in 1:length(S.gaugeops)]...)
gaugesmatrix(::HasNoGauges, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no gauges.")
gaugeoperatorsmatrix(S::AbstractSubsystemCode) = gaugesmatrix(S)

"""
    dressed(S::AbstractSubsystemCode)
    dressedoperators(S::AbstractSubsystemCode
    dressedlogicals(S::AbstractSubsystemCode)

Return a vector of pairs generators for the dressed operators of `S`.

Here, the dressed operators are the logicals and the gauge operators.
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

Here, this is the stabilizers and the gauge operators.
"""
gaugegroup(S::T) where {T <: AbstractSubsystemCode} = gaugegroup(GaugeTrait(T), S)
gaugegroup(::HasGauges, S::AbstractSubsystemCode) =
    vcat(S.stabs, vcat([vcat(S.gaugeops[i]...) for i in 1:length(S.gaugeops)]...))
gaugegroup(::HasNoGauges, S::AbstractSubsystemCode) = error("Type $(typeof(S)) has no gauges.")
gaugegroupmatrix(S::AbstractSubsystemCode) = gaugegroup(S)
gaugegeneratorsmatrix(S::AbstractSubsystemCode) = gaugegroup(S)
gaugegroupgeneratorsmatrix(S::AbstractSubsystemCode) = gaugegroup(S)

#############################
      # setter functions
#############################

"""
    changesigns!(S::AbstractSubsystemCode, charvec::Vector{nmod})

Set the character vector of `S` to `charvec` and update the signs.
"""
function changesigns!(S::AbstractSubsystemCode, charvec::Vector{nmod})
    R = base_ring(charactervector(S))
    length(charvec) == 2 * S.n || error("Characteristic vector is of improper length for the code.")
    for s in charvec
        modulus(s) == modulus(R) || error("Phases are not in the correct ring.")
    end
    S.signs = _getsigns(S.stabilizers, charvec)
    S.charvec = charvec
end

"""
    setstabilizers!(S::AbstractSubsystemCode, stabs::fq_nmod_mat, symp::Bool=true)

Set the stabilizers of `S` to `stabs`.

If the optional parameter `symp` is set to `true`, `stabs` is assumed to be in
symplectic form over the base field of `S`. A check is done to make sure `stabs`
is isomorphic to the current set of stabilizers.
"""
function setstabilizers!(S::AbstractSubsystemCode, stabs::fq_nmod_mat, symp::Bool=true)
    iszero(stabs) && throw(ArgumentError("The stabilizers cannot be zero."))
    stabs = _removeempty(stabs, "rows")
    if symp
        order(S.F) == order(base_ring(stabs)) || throw(ArgumentError("The stabilizers must be over the same field as the code."))
        stabs = change_base_ring(S.F, stabs)
        if _isisomorphic(symplecticstabilizers(S), stabs)
            S.stabs = symplectictoquadratic(stabs)
            nrows(stabs) != S.k && (S.overcomplete = true;)
        else
            error("The current stabilizers are not isomorphic to the input.")
        end
    else
        order(S.E) == order(base_ring(stabs)) || throw(ArgumentError("The stabilizers must be over the same field as the code."))
        stabs = change_base_ring(S.E, stabs)
        if _isisomorphic(symplecticstabilizers(S), quadratictosymplectic(stabs))
            S.stabs = stabs
            nrows(stabs) != S.k && (S.overcomplete = true;)
        else
            error("The current stabilizers are not isomorphic to the input.")
        end
    end

    # TODO: update signs
    return nothing
end

"""
    setlogicals!(S::AbstractSubsystemCode, L::fq_nmod_mat, symp::Bool=false)

Set the logical operators of `S` to `L`.

If the optional parameter `symp` is set to `true`, `L` is assumed to be in
symplectic form over the base field of `S`. A check is done to make sure `L`
is isomorphic to the current set of logicals.
"""
setlogicals!(S::T, L::fq_nmod_mat, symp::Bool=false) where {T <: AbstractSubsystemCode} = setlogicals!(LogicalTrait(T), S, L, symp)
# TODO: also take in vector format, this probably should just be copy and paste from above
function setlogicals!(::HasLogicals, S::AbstractSubsystemCode, L::fq_nmod_mat, symp::Bool=false)
    if symp
        Lsym = L
        size(L) == (2 * S.k, 2 * S.n) || error("Provided matrix is of incorrect size for the logical space.")
        iseven(ncols(L)) || error("Expected a symplectic input but the input matrix has an odd number of columns.")
        S.F == base_ring(Lsym) || error("The logicals must be over the same field as the code.")
        L = symplectictoquadratic(L)
    else
        E = base_ring(L)
        size(L) == (2 * S.k, S.n) || error("Provided matrix is of incorrect size for the logical space.")
        S.E == E || error("The logicals must be over the same field as the code.")
        iseven(degree(E)) || error("The base ring of the given matrix is not a quadratic extension.")
        Lsym = quadratictosymplectic(L)
    end
    _isisomorphic(Lsym, quadratictosymplectic(logicalsmatrix(S))) || error("The current logicals are not isomorphic to the input.")
    # aresymplecticorthogonal(symplecticstabilizers(S), Lsym, true) ||
    #     error("Provided logicals do not commute with the code.")

    # the columns in prod give the commutation relationships between the provided
    # logical operators; they ideally should only consist of {X_1, Z_i} pairs
    # so there should only be one nonzero element in each column
    prod = hcat(Lsym[:, S.n + 1:end], -Lsym[:, 1:S.n]) * transpose(Lsym)
    iszero(prod) && error("Provided logicals should not be symplectic self-orthogonal.")
    ncpr = ncols(prod)
    # need an integer sum and this is cheaper than Nemo.ZZ
    prodJul = FpmattoJulia(prod)
    cols = [sum(prodJul[:, i]) for i in 1:ncpr]
    sum(cols) == ncpr || error("Incorrect commutation relationships between provided logicals.")

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
end
setlogicals!(::HasNoLogicals, S::AbstractSubsystemCode, L::fq_nmod_mat, symp::Bool=false) = error("Type $(typeof(S)) has no logicals.")

"""
    setminimumdistance(S::AbstractSubsystemCode, d::Int)

Set the minimum distance of the code to `d`.

The only check done on the value of `d` is that `1 ≤ d ≤ n`.
"""
function setminimumdistance!(S::AbstractSubsystemCode, d::Int)
    # TODO: should check bounds like Singleton for possibilities
    d > 0 && d <= S.n || error("The minimum distance of a code must be ≥ 1; received: d = $d.")
    S.d = d
end

#############################
     # general functions
#############################

function _processcharvec(charvec::Union{Vector{nmod}, Missing}, p::Int, n::Int)
    if !ismissing(charvec)
        2 * n == length(charvec) || error("The characteristic value is of incorrect length.")
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        for s in charvec
            modulus(s) == modulus(R) || error("Phases are not in the correct ring.")
        end
    else
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        charvec = [R(0) for _ in 1:2 * n]
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
    if length(charvec) == 2 * nc
        A = quadratictosymplectic(A)
    end
    length(charvec) == nc || error("Input to _getsigns is expected to be in symplectic form and of the same length as the characteristic vector.")
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

"""
    splitstabilizers(S::AbstractSubsystemCode)

Return the set of `X`-only stabilizers and their signs, the set of `Z`-only
stabilizers and their signs, and the remaining stabilizers and their signs.

This function returns six objects of alternating types `fq_nmod_mat` and
`Vector{Int}` for the three sets of stabilizers and signs, respectively.
An empty set of stabilizers is returned as type `Vector{fq_nmod_mat}`.
"""
# TODO: type stability is probably shit here
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
    E = base_ring(L)
    n = ncols(L)
    logs = Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}()
    # this does indeed grow smaller each iteration
    while nrows(L) >= 2
        Lsym = quadratictosymplectic(L)
        # the columns in prod give the commutation relationships between the provided
        # logical operators; they ideally should only consist of {X_1, Z_i} pairs
        # so there should only be one nonzero element in each column
        prod = hcat(Lsym[:, n + 1:end], -Lsym[:, 1:n]) * transpose(Lsym)
        nprod = ncols(prod)
        first = 0
        for c in 1:nprod
            if !iszero(prod[1, c])
                if iszero(first)
                    first = c
                    if !isone(prod[1, c])
                        L[first, :] *= E(prod[1, c]^-1)
                    end
                else
                    L[c, :] += E(prod[1, c]^-1) * L[first, :]
                end
            end
        end
        # TODO: need to check when iszero(first)
        for c in 2:nprod
            if !iszero(prod[first, c])
                L[c, :] += E(prod[first, c]^-1) * L[1, :]
            end
        end
        push!(logs, (L[1, :], L[first, :]))
        L = L[setdiff(1:nrows(L), [1, first]), :]
    end
    return logs
end

_testlogicalsrelationships(S::T) where {T <: AbstractSubsystemCode} = _testlogicalsrelationships(LogicalTrait(T), S)
function _testlogicalsrelationships(::HasLogicals, S::AbstractSubsystemCode)
    L = vcat([vcat(S.logicals[i]...) for i in 1:S.k]...)
    Lsym = quadratictosymplectic(L);
    prod = hcat(Lsym[:, S.n + 1:end], -Lsym[:, 1:S.n]) * transpose(Lsym)
    display(prod)
    return nothing
end
_testlogicalsrelationships(::HasNoLogicals, S) = error("Type $(typeof(S)) has no logicals.")

"""
    islogical(S::AbstractSubsystemCode, v::fq_nmod_mat)

Return `true` if the vector `v` anticommutes with any of the logical operators.
"""
# TODO: check for type stability here
islogical(S::T, v::fq_nmod_mat) where {T <: AbstractSubsystemCode} = islogical(LogicalTrait(T), S, v)
function islogical(::HasLogicals, S::AbstractSubsystemCode, v::fq_nmod_mat)
    logsmat = quadratictosymplectic(logicalsmatrix(S))
    nc = ncols(logsmat)
    size(v) == (1, nc) && (return !iszero(logsmat * transpose(v));)
    size(v) == (nc, 1) && (return !iszero(logsmat * v);)
    error("Vector to be tested is of incorrect dimension.")
end
islogical(::HasNoLogicals, S::AbstractSubsystemCode, v::fq_nmod_mat) = error("Type $(typeof(S)) has no logicals.")

"""
    syndrome(S::AbstractSubsystemCode, v::fq_nmod_mat)

Return the syndrome of the vector `v` with respect to the stabilizers of `S`.
"""
# TODO: make uniform with approach in function above
function syndrome(S::AbstractSubsystemCode, v::fq_nmod_mat)
    (size(v) != (2 * S.n, 1) && size(v) != (1, 2 * S.n)) &&
        error("Vector to be tested is of incorrect dimension; expected length $(2 * n), received: $(size(v)).")
    # base_ring(v) == field(S) || error("Vector must have the same base ring as the stabilizers.")

    nrows(v) != 1 || return symplecticstabilizers(S) * transpose(v)
    return symplecticstabilizers(S) * v
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
# TODO: graph states, CSS codes
# can fix both by simply recalling the constructor
promotelogicalstogauge(S::T, pairs::Vector{Int}) where {T <: AbstractSubsystemCode} = logicalpromotelogicalstogauges(LogicalTrait(T), S, pairs)
function promotelogicalstogauge(::HasLogicals, S::AbstractSubsystemCode, pairs::Vector{Int})
    pairs = sort!(unique!(pairs))
    stabs = S.stabs
    logs = S.logicals
    # will let this error naturally if pairs contains invalid elements
    gaugeops = logs[pairs]
    logs = logs[setdiff![1:S.k, pairs]]
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
    return SubsystemCode(S.F, S.E, S.n, k, r, S.d, stabs, logs, S.charvec, S.signs, gaugeops, overcomplete)
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
    return nothing
end
swapXZgaugeoperators!(::HasNoGauges, S::AbstractSubsystemCode, pairs::Vector{Int}) = error("Type $(typeof(S)) has no gauges.")

"""
    isisomorphic(S1::T, S2::T) where T <: AbstractSubsystemCode

Return `true` if the subsystem or stabilizer codes are isomorphic as
symplectic vector spaces.

* Note
- This is not intended to detect if `S1` and `S2` are permutation equivalent.
"""
# TODO: redo hierarchy based on conceptualparent(s) or traits
function isisomorphic(S1::T, S2::T) where T <: AbstractSubsystemCode
    (S1.n == S2.n && S1.k == S2.k) || return false
    # can't compare fields directly because they are compared based on ptr addresses
    (Int(order(S1.F)) == Int(order(S2.F)) && Int(order(S1.E)) == Int(order(S2.E))) || return false
    if GaugeTrait(T) == HasGauges
        S1.r == S2.r || return false
    end


    # V = VectorSpace(S1.F, 2 * S1.n)
    # # test stabilizers
    # S1symstabs = symplecticstabilizers(S1)
    # S2symstabs = symplecticstabilizers(S2)
    # S1VS, _ = sub(V, [V(S1symstabs[i, :]) for i in 1:nrows(S1symstabs)])
    # S2VS, _ = sub(V, [V(S2symstabs[i, :]) for i in 1:nrows(S2symstabs)])
    # is_isomorphic(S1VS, S2VS) || return false
    _isisomorphic(symplecticstabilizers(S1), symplecticstabilizers(S2)) || return false

    # graph states
    if LogicalTrait(T) == HasLogicals()
        # test logicals
        # S1symlogs = quadratictosymplectic(logicalsmatrix(S1))
        # S2symlogs = quadratictosymplectic(logicalsmatrix(S2))
        # S1logsVS, _ = sub(V, [V(S1symlogs[i, :]) for i in 1:nrows(S1symlogs)])
        # S2logsVS, _ = sub(V, [V(S2symlogs[i, :]) for i in 1:nrows(S2symlogs)])
        # is_isomorphic(S1logsVS, S2logsVS) || return false
        _isisomorphic(quadratictosymplectic(logicalsmatrix(S1)), quadratictosymplectic(logicalsmatrix(S2))) || return false
    else
        return true
    end

    if GaugeTrait(T) == HasGauges()
        # test gauge operators
        # S1symgops = quadratictosymplectic(gaugesmatrix(S1))
        # S2symgops = quadratictosymplectic(gaugesmatrix(S2))
        # S1gopsVS, _ = sub(V, [V(S1symgops[i, :]) for i in 1:nrows(S1symgops)])
        # S2gopsVS, _ = sub(V, [V(S2symgops[i, :]) for i in 1:nrows(S2symgops)])
        return is_isomorphic(quadratictosymplectic(gaugesmatrix(S1)), quadratictosymplectic(gaugesmatrix(S2)))
    else
        return true
    end
end

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
