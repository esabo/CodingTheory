# Copyright (c) 2021, 2022, 2023 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    StabilizerCodeCSS(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)
    CSSCode(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)

Return the CSS code given by the CSS construction on two linear codes `C1 = [n, k1, d1]`
and `C2 = [n, k2, d2]` with `C2 ⊆ C1` and whose signs by `charvec`.

# Notes
* The resulting code has dimension `k = k1 - k2` and minimum distance
  `d >= min(d1, d2^⟂)`. The `X` stabilizers are given by the parity-check matrix
  of `C2^⟂`, `H(C2^⟂)`, and the `Z` stabilizers by `H(C1)`.
"""
function StabilizerCodeCSS(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)
    C2 ⊆ C1 || throw(ArgumentError("The second argument must be a subset of the first in the CSS construction."))
    p = Int(characteristic(C1.F))
    charvec = _processcharvec(charvec, p, C1.n)

    # C2 ⊆ C1
    # k = k1 - k2
    # d >= minimum(d1, d2^⟂)
    # X - H(C2^⟂), Z - H(C1)
    D2 = dual(C2)
    S = directsum(D2.H, C1.H)
    logs, logsmat = _logicals(S, directsum(C1.G, D2.G))

    # determine signs
    signs, Xsigns, Zsigns = _determinesignsCSS(S, charvec, nrows(D2.H), nrows(C1.H))

    # q^n / p^k but rows is n - k
    rkS = rank(S)
    if rkS != C1.n
        dimcode = BigInt(order(C1.F))^ncols(C1.n) // BigInt(p)^rkS
        isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

        return StabilizerCodeCSS(C1.F, C1.n, dimcode, missing, missing, missing, S, D2.H, C1.H,
            C2, C1, signs, Xsigns, Zsigns, logs, logsmat, charvec, missing, missing, missing,
            false, missing)
    else
        return GraphStateStabilizerCSS(C1.F, C1.n, 0, missing, D2.d, C1.d, S, D2.H, C1.H, C2, C1,
            signs, Xsigns, Zsigns, charvec, missing, false)
    end
end
CSSCode(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing) =
    StabilizerCodeCSS(C1, C2, charvec)

"""
    StabilizerCodeCSS(C::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)
    CSSCode(C::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)

Return the CSS code given by the CSS construction on a self-orthogonal linear code
`C`, i.e., `C ⊆ C^⟂`, and whose signs by `charvec`.

# Notes
* Setting `C1 = C^⟂` and `C2 = C`, the resulting code has dimension `k = k1 - k2`
  and minimum distance `d >= min(d1, d2^⟂)`. The `X` stabilizers are given by the
  parity-check matrix of `C2^⟂`, `H(C2^⟂)`, and the `Z` stabilizers by `H(C1)`.
"""
function StabilizerCodeCSS(C::LinearCode, charvec::Union{Vector{nmod}, Missing}=missing)
    # this should have Xstabs = Zstabs
    D = dual(C)
    C ⊆ D || throw(ArgumentError("The single code CSS construction requires C ⊆ C^⟂."))
    p = Int(characteristic(D.F))
    charvec = _processcharvec(charvec, p, D.n)

    # C2 ⊆ C1
    # k = k1 - k2
    # d >= minimum(d1, d2^⟂)
    # X - H(C2^⟂), Z - H(C1)
    S = directsum(D.H, D.H)
    logs, logsmat = _logicals(S, directsum(D.G, D.G))

    # determine signs
    nr = nrows(D.H)
    signs, Xsigns, Zsigns = _determinesignsCSS(S, charvec, nr, nr)

    # q^n / p^k but rows is n - k
    rkS = rank(S)
    if rkS != D.n
        dimcode = BigInt(order(D.F))^D.n // BigInt(p)^rkS
        isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

        return StabilizerCodeCSS(D.F, D.n, dimcode, missing, missing, missing, S, D.H, D.H, C,
            D, signs, Xsigns, Zsigns, logs, logsmat, charvec, missing, missing, missing, false,
            missing)
    else
        return GraphStateStabilizerCSS(D.F, D.n, 0, missing, D.d, D.d, S, D.H, D.H, C, D, signs,
            Xsigns, Zsigns, charvec, missing, false)
    end
end
CSSCode(C::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing) = StabilizerCodeCSS(C, charvec)

"""
    StabilizerCodeCSS(Xmatrix::fq_nmod_mat, Zmatrix::fq_nmod_mat, charvec::Union{Vector{nmod}, Missing}=missing)
    CSSCode(Xmatrix::fq_nmod_mat, Zmatrix::fq_nmod_mat, charvec::Union{Vector{nmod}, Missing}=missing)

Return a CSS code whose `X`-stabilizers are given by `Xmatrix`, `Z`-stabilizers by `Zmatrix`, and signs by `charvec`.
"""
function StabilizerCodeCSS(Xmatrix::fq_nmod_mat, Zmatrix::fq_nmod_mat, charvec::Union{Vector{nmod}, Missing}=missing)
    iszero(Xmatrix) && throw(ArgumentError("The `X` stabilizer matrix is empty."))
    iszero(Zmatrix) && throw(ArgumentError("The `Z` stabilizer matrix is empty."))
    n = ncols(Xmatrix)
    n ==  ncols(Zmatrix) || throw(ArgumentError("Both matrices must have the same length in the CSS construction."))
    F = base_ring(Xmatrix)
    F == base_ring(Zmatrix) || throw(ArgumentError("Both matrices must be over the same base field in the CSS construction."))
    iszero(Zmatrix * transpose(Xmatrix)) || throw(ArgumentError("The given matrices are not symplectic orthogonal."))
    p = Int(characteristic(F))
    charvec = _processcharvec(charvec, p, n)

    Xmatrix = _removeempty(Xmatrix, :rows)
    Zmatrix = _removeempty(Zmatrix, :rows)

    # determine if the provided set of stabilizers are redundant
    Xrank = rank(Xmatrix)
    Zrank = rank(Zmatrix)
    if nrows(Xmatrix) > Xrank || nrows(Zmatrix) > Zrank
        overcomp = true
    else
        overcomp = false
    end

    S = directsum(Xmatrix, Zmatrix)
    signs, Xsigns, Zsigns = _determinesignsCSS(S, charvec, nrows(Xmatrix), nrows(Zmatrix))

    # find generators for S^⟂
    # note the H here is transpose of the standard definition
    _, H = right_kernel(hcat(S[:, n + 1:end], -S[:, 1:n]))
    # n + (n - Srank)
    ncols(H) == 2 * n - Xrank - Zrank || error("Normalizer matrix is not size n + k.")
    logs, logsmat = _logicals(S, transpose(H))

    # q^n / p^k but rows is n - k
    rkS = Xrank + Zrank
    if rkS != n
        dimcode = BigInt(order(F))^n // BigInt(p)^rkS
        isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

        return StabilizerCodeCSS(F, n, dimcode, missing, missing, missing, S, Xmatrix, Zmatrix,
            missing, missing, signs, Xsigns, Zsigns, logs, logsmat, charvec, missing, missing,
            missing, overcomp, missing)
    else
        return GraphStateStabilizerCSS(F, n, 0, missing, missing, missing, S, Xmatrix, Zmatrix,
            missing, missing, signs, Xsigns, Zsigns, charvec, missing, overcomp)
    end
end
CSSCode(Xmatrix::fq_nmod_mat, Zmatrix::fq_nmod_mat, charvec::Union{Vector{nmod}, Missing}=missing) =
    StabilizerCodeCSS(Xmatrix, Zmatrix, charvec)

"""
    StabilizerCodeCSS(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}
    CSSCode(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}

Return the CSS code whose stabilizers are determined by the vector of Pauli strings `SPauli` and signs by `charvec`.

# Notes
* Any +/- 1 characters in front of each stabilizer are stripped. No check is done
  to make sure these signs agree with the ones computed using the character vector.
"""
function StabilizerCodeCSS(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}
    S = _Paulistringtosymplectic(_processstrings(SPauli))
    iszero(S) && throw(ArgumentError("The processed Pauli strings returned a set of empty stabilizer generators."))
    S = _removeempty(S, :rows)
    # the reason we repeat here and not call another constructor is the else
    # statement at the bottom of this function
    # would also need to compute down to signs to call _isCSSsymplectic
    # which would allow us to call the other constructor
    aresymplecticorthogonal(S, S) || throw(ArgumentError("The given stabilizers are not symplectic orthogonal."))
    n = div(ncols(S), 2)

    F = base_ring(S)
    p = Int(characteristic(F))
    charvec = _processcharvec(charvec, p, 2 * n)
    signs = _determinesigns(S, charvec)
    
    # determine if the provided set of stabilizers are redundant
    rkS = rank(S)
    if nrows(S) > rkS
        overcomp = true
    else
        overcomp = false
    end

    # find generators for S^⟂
    # note the H here is transpose of the standard definition
    _, H = right_kernel(hcat(S[:, n + 1:end], -S[:, 1:n]))
    # n + (n - rkS)
    ncols(H) == 2 * n - rkS || error("Normalizer matrix is not size n + k.")
    logs, logsmat = _logicals(S, transpose(H))

    # q^n / p^k but rows is n - k
    args = _isCSSsymplectic(S, signs, true)
    if args[1]
        if rkS != n
            dimcode = BigInt(order(F))^n // BigInt(p)^rkS
            isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

            return StabilizerCodeCSS(F, n, dimcode, missing, missing, missing, S, args[2],
                args[4], missing, missing, signs, args[3], args[5], logs, logsmat, charvec,
                missing, missing, missing, overcomp, missing)
        else
            return GraphStateStabilizerCSS(F, n, 0, missing, missing, missing, S, args[2],
                args[4], missing, missing, signs, args[3], args[5], charvec, missing, overcomp)
        end
    else
        error("Provided Pauli strings are not CSS.")
    end
end
CSSCode(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String,
    Vector{Char}} = StabilizerCodeCSS(SPauli, charvec)

"""
    StabilizerCodeCSS(S::StabilizerCode)
    CSSCode(S::StabilizerCode)

Return the `[[2n, 2k, S.d <= d <= 2 S.d]]` CSS code derived by splitting the stabilizers of `S`.
"""
function StabilizerCodeCSS(S::StabilizerCode)
	X = S.stabs[:, 1:S.n]
	Z = S.stabs[:, S.n + 1:end]
	return StabilizerCodeCSS(X, Z, S.charvec)
end
CSSCode(S::StabilizerCode) = StabilizerCodeCSS(S)

# entanglement-assisted is not symplectic orthogonal
"""
    StabilizerCode(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}

Return the stabilizer code whose stabilizers are determined by the vector of Pauli strings `SPauli` and signs by `charvec`.

# Notes
* Any +/- 1 characters in front of each stabilizer are stripped. No check is done
  to make sure these signs agree with the ones computed using the character vector.
"""
function StabilizerCode(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}
    SPaulistripped = _processstrings(SPauli)
    S = _Paulistringtosymplectic(SPaulistripped)
    iszero(S) && throw(ArgumentError("The processed Pauli strings returned a set of empty stabilizer generators."))
    return StabilizerCode(S, charvec)
end

"""
    StabilizerCode(S::fq_nmod_mat, charvec::Union{Vector{nmod}, Missing}=missing)

Return the stabilizer code whose stabilizers is determined by `S` and signs by `charvec`.
"""
function StabilizerCode(S::fq_nmod_mat, charvec::Union{Vector{nmod}, Missing}=missing)
    iszero(S) && throw(ArgumentError("The stabilizer matrix is empty."))
    S = _removeempty(S, :rows)
    aresymplecticorthogonal(S, S) || throw(ArgumentError("The given stabilizers are not symplectic orthogonal."))

    F = base_ring(S)
    p = Int(characteristic(F))
    n = div(ncols(S), 2)
    charvec = _processcharvec(charvec, p, n)
    signs = _determinesigns(S, charvec)

    # determine if the provided set of stabilizers are redundant
    rkS = rank(S)
    if nrows(S) > rkS
        overcomp = true
    else
        overcomp = false
    end

    # find generators for S^⟂
    # note the H here is transpose of the standard definition
    _, H = right_kernel(hcat(S[:, n + 1:end], -S[:, 1:n]))
    # n + (n - rkS)
    ncols(H) == 2 * n - rkS || error("Normalizer matrix is not size n + k.")
    logs, logsmat = _logicals(S, transpose(H))

    # q^n / p^k but rows is n - k
    dimcode = BigInt(order(F))^n // BigInt(p)^rkS
    isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

    args = _isCSSsymplectic(S, signs, true)
    if args[1]
        if rkS != n
            return StabilizerCodeCSS(F, n, dimcode, missing, missing, missing, S, args[2],
                args[4], missing, missing, signs, args[3], args[5], logs, logsmat, charvec,
                missing, missing, missing, overcomp, missing)
        else
            return GraphStateStabilizerCSS(F, n, 0, missing, missing, missing, S, args[2],
                args[4], missing, missing, signs, args[3], args[5], charvec, missing, overcomp)
        end
    else
        if rkS != n
            return StabilizerCode(F, n, dimcode, missing, S, logs, logsmat, charvec, signs, missing,
                missing, missing, overcomp, missing)
        else
            return GraphState(F, n, 0, missing, S, charvec, signs, missing, overcomp)
        end
    end
end

function _logicals(stabs::fq_nmod_mat, dualgens::fq_nmod_mat)
    L = _quotientspace(dualgens, stabs)
    logs = _makepairs(L)
    # verify
    n = div(ncols(L), 2)
    logsmat = vcat([vcat(logs[i]...) for i in 1:length(logs)]...)
    aresymplecticorthogonal(stabs, logsmat) || error("Computed logicals do not commute with the codespace.")
    prod = hcat(logsmat[:, n + 1:end], -logsmat[:, 1:n]) * transpose(logsmat)
    sum(FpmattoJulia(prod), dims=1) == ones(Int, 1, size(prod, 1)) ||
        error("Computed logicals do not have the right commutation relations.")
    return logs, logsmat
end

"""
    augment(S::AbstractStabilizerCode, row::fq_nmod_mat, verbose::Bool=true)

Return the code created by added `row` to the stabilizers of `S`.

# Notes
* The goal of this function is to track how the logical operators update given the new stabilizer.
  The unaffected logical operators are kept during the update and only those which don't commute
  with the new stabilizer are recomputed. Use `verbose` to better 
"""
# TODO: move to subsystem
function augment(S::AbstractStabilizerCode, row::fq_nmod_mat, verbose::Bool=true)
    # typeof(S) ∈ [GraphState, GraphStateStabilizerCSS] && return S
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
    # this redoes extra work but let's start by seeing if it works
    return SubsystemCode(stabs, vcat(logs, newlogs), gaugeops, S.charvec)

    # if LogicalTrait(typeof(S)) == HasLogicals() && GaugeTrait(typeof(S)) == HasGauges()

    # elseif LogicalTrait(typeof(S)) == HasNoLogicals() && GaugeTrait(typeof(S)) == HasGauges()
    #     # do the above but don't include logs
    # elseif LogicalTrait(typeof(S)) == HasLogicals() && GaugeTrait(typeof(S)) == HasNoGauges()
    #     # do the above but don't include gauges
    # elseif LogicalTrait(typeof(S)) == HasNoLogicals() && GaugeTrait(typeof(S)) == HasNoGauges()
    #     # being explicit here for the reader instead of using else but this is a compile time
    #     # choice so won't slow anything down

    #     # pass stabs to StabilizerCode constructor
    # end



    # rankS = rank(S.stabs)
    # newsymstabs = vcat(S.stabs, row)
    # ranknewS = rank(newsymstabs)
    # if rankS == ranknewS
    #     verbose && println("Row is already in the stabilizer group. Nothing to update.")
    #     Snew = deepcopy(S)
    #     Snew.stabs = newsymstabs
    #     Snew.overcomp = true
    #     return Snew
    # elseif S.k == 1
    #     verbose && println("Row is not in the stabilizer group; the result is a graph state.")
    #     ranknewS == nrows(newsymstabs) ? (overcomp = false;) : (overcomp = true;)
    #     if GaugeTrait(typeof(S)) == HasGauges()
    #         return GraphStateSubsystem(S.F, S.n, 0, S.r, missing, newsymstabs, S.charvec, _determinesigns(newsymstabs,
    #             charvec), missing, overcomp, S.gaugeops, S.gopsmat)
    #     else
    #         return GraphState(S.F, S.n, 0, missing, newsymstabs, S.charvec, _determinesigns(newsymstabs,
    #             charvec), missing, overcomp)
    #     end
    # end


    #     if isempty(logstokeep)
    #         verbose && println("Row does not commute with any logicals. The entire code needs to be recomputed from scratch.")
    #         if GaugeTrait(typeof(S)) == HasGauges()
    #             return SubsystemCode(vcat(gaugegroupmatrix(S), row))
    #         else
    #             return StabilizerCode(newsymstabs, S.charvec)
    #         end
    #     end

        
    #     # kernel should contain newstabs and new logs but not old logs
        
    #     # verify
    #     fulllogs = [S.logicals[logpairstokeep]; newlogs]
    #     logsmat = vcat([vcat(fulllogs[i]...) for i in 1:length(fulllogs)]...)
    #     aresymplecticorthogonal(newsymstabs, logsmat) || error("Computed logicals do not commute with the codespace.")
    #     prod = hcat(logsmat[:, S.n + 1:end], -logsmat[:, 1:S.n]) * transpose(logsmat)
    #     sum(FpmattoJulia(prod), dims=1) == ones(Int, 1, size(prod, 1)) || error("Computed logicals do not have the right commutation relations.")
    #     # set and return if good
    #     verbose && println("New logicals:")
    #     verbose && display(newlogs)
    #     Snew = StabilizerCode(newsymstabs, S.charvec)
    #     setlogicals!(Snew, fulllogs)
    # end
    # return Snew
end

"""
    expurgate(S::AbstractStabilizerCode, rows::Vector{Int}, verbose::Bool=true)

Return the code created by removing the stabilizers indexed by `rows`.

# Notes
* The goal of this function is to track how the logical operators update through this process.
  Here, the original logical pairs are kept and an appropriate number of new pairs are added.
"""
# TODO: move to subsystem, unify prints with above
function expurgate(S::AbstractSubsystemCode, rows::Vector{Int}, verbose::Bool=true)
    numstabs = nrows(S.stabs)
    rows ⊆ 1:numstabs || throw(ArgumentError("Argument `rows` not a subset of the number of stabilizers."))
    verbose && println("Removing stabilizers: $rows")
    newstabs = S.stabs[setdiff(1:numstabs, rows), :]
    typeof(S) <: AbstractStabilizerCode ? (Snew = StabilizerCode(newstabs, S.charvec);) :
        (Snew = SubsystemCode(newstabs, S.charvec);)
    
    # never going to be a graph state because it just gained a logical
    # never going to be zero because StabilizerCode would have errored
    # this is going to check if it was previously a graph state
    if LogicalTrait(typeof(S)) == HasLogicals()
        logs = logicals(S)
        logsmatrix = logicalsmatrix(S)
        small = vcat(newstabs, logsmatrix)
        # println(size(small))
        _, H = right_kernel(hcat(small[:, S.n + 1:end], -small[:, 1:S.n]))
        H = transpose(H)
        # println(size(H))
        # dualgenssym = hcat(H[:, S.n + 1:end], -H[:, 1:S.n])
        # println(size(dualgenssym))
        temp = _quotientspace(H, newstabs)
        # temp = hcat(temp[:, S.n + 1:end], -temp[:, 1:S.n])
        # using dualgenssym then switching temp here just switches {X, Z} to {Z, X}
        # but the vectors remain the same for some reason
        newlogs = _makepairs(temp)
        # verify
        fulllogs = [logs; newlogs]
        logsmatrix = vcat([vcat(fulllogs[i]...) for i in 1:length(fulllogs)]...)
        aresymplecticorthogonal(newstabs, logsmatrix) || error("Computed logicals do not commute with the codespace.")
        prod = hcat(logsmatrix[:, S.n + 1:end], -logsmatrix[:, 1:S.n]) * transpose(logsmatrix)
        sum(FpmattoJulia(prod), dims=1) == ones(Int, 1, size(prod, 1)) || error("Computed logicals do not have the right commutation relations.")
        # set and return if good
        verbose && println("New logicals:")
        verbose && display(newlogs)
        setlogicals!(Snew, logsmatrix)
        return Snew
    else
        verbose && println("Started with all graph state. New logicals:")
        verbose && display(logicals(Snew))
    end
end

#############################
#   Generator Coefficients  #
#############################


# # to run this, need an error model,
# # need wrtV,
# # need clarification on general stabilizer codes
# # need to process Pauli here
# # sum of syndrome and logical?
# function generatorcoefficients(Q::AbstractStabilizerCode, θ::Union{Float64, Missing},
#     Paui::Char=' ')
#
#     synlen = size(stabilizers(Q), 1)
#     numsyns = BigInt(characteristic(field(Q)))^synlen
#     # do some checks here for feasibility
#     W, T = Pauliweightenumerator(Q, Pauli, true, true)
#     prevsyn = zeros(Int64, synlen)
#     logs = logicalspace(Q)
#
#     for i in 1:size(logs, 1)
#         γ = logs[i, :]
#         for j in 1:numsyns
#             μ = digits(j, base=2, pad=synlen)
#             shifttrellis!(T, μ .+ γ, wrtV, err_model, chractervector(Q), Pauli)
#             push!(W, _Pauliweightenumerator(T))
#             _reducepoly(W)
#         end
#     end
#
#     if ismissing(θ)
#         # @vars θ
#         n = length(Q)
#         sum = Complex(0.0) # typeof(sum) = Sym
#         for term in W
#             if term[2] < 0 || term[3] < 0 || term[4] < 0
#                 sum -= term[1] * cos(θ / 2)^abs(term[2] + term[3]) * (1im * sin(θ / 2))^abs(term[3] + term[4])
#             else
#                 sum += term[1] * cos(θ / 2)^abs(term[2] + term[3]) * (1im * sin(θ / 2))^abs(term[3] + term[4])
#             end
#         end
#         return sum, W
#     else
#         n = length(Q)
#         sum = Complex(0.0)
#         for term in W
#             if term[2] < 0 || term[3] < 0 || term[4] < 0
#                 sum -= term[1] * cos(θ / 2)^(n - abs(term[2] + term[3])) * (1im * sin(θ / 2))^abs(term[3] + term[4])
#             else
#                 sum += term[1] * cos(θ / 2)^(n - abs(term[2] + term[3])) * (1im * sin(θ / 2))^abs(term[3] + term[4])
#             end
#         end
#         return sum, W
#     end
# end

# function plotgeneratorcoefficients(W, θmin, Qmax)
#
# end
#
# function plotgeneratorcoefficients(Q, θmin, Qmax)
#
# end
