# Copyright (c) 2021, 2022 Eric Sabo
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

Return a CSS code using the CSS construction on two linear codes `C1 = [n, k1, d1]`
and `C2 = [n, k2, d2]` with `C2 ⊆ C1`.

The resulting code has dimension `k = k1 - k2` and minimum distance
`d >= min(d1, d2^⟂)`. The `X` stabilizers are given by the parity-check matrix
of `C2^⟂`, `H(C2^⟂)`, and the `Z` stabilizers by `H(C1)`.

# Arguments
* `C1`: a linear code
* `C2`: a subcode of `C1`
* `charvec`: a length `2n` vector with elements in the `Z/(2p)` if
  `chracteristic(field(C1))` is 2 and `Z/(p)` otherwise. The first `n` elements
  specify the exponents of the `X` phases and second `n` the exponents of the
  `Z` phases; a missing argument will be set to the all-zero vector

# Notes
* A `+1` phase should be entered as `0` since the character vector stores the
  exponents.
* Stabilizer signs are automatically computed given the character vector.
"""
function StabilizerCodeCSS(C1::AbstractLinearCode, C2::AbstractLinearCode,
    charvec::Union{Vector{nmod}, Missing}=missing)

    C2 ⊆ C1 || error("The second argument must be a subset of the first in the CSS construction.")
    p = Int(characteristic(C1.F))
    charvec = _processcharvec(charvec, p, C1.n)

    # C2 ⊆ C1
    # k = k1 - k2
    # d >= minimum(d1, d2^⟂)
    # X - H(C2^⟂), Z - H(C1)
    D2 = dual(C2)
    S = directsum(D2.H, C1.H)
    Sq2 = symplectictoquadratic(S)
    logs = _logicals(S, directsum(C1.G, D2.G))

    # determine signs
    signs, Xsigns, Zsigns = _determinesignsCSS(S, charvec, nrows(D2.H), nrows(C1.H))

    # q^n / p^k but rows is n - k
    rkS = rank(S)
    if rkS != C1.n
        dimcode = BigInt(order(C1.F))^ncols(Sq2) // BigInt(p)^rkS
        isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

        return StabilizerCodeCSS(C1.F, base_ring(Sq2), C1.n, dimcode, missing, missing, missing, Sq2, D2.H, C1.H,
            C2, C1, signs, Xsigns, Zsigns, logs, charvec, missing,
            missing, missing, false, missing)
    else
        return GraphStateStabilizerCSS(C1.F, base_ring(Sq2), C1.n, 0, missing, D2.d, C1.d, Sq2, D2.H, C1.H,
            C2, C1, signs, Xsigns, Zsigns, charvec, missing, false)
    end
end
CSSCode(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing) =
    StabilizerCodeCSS(C1, C2, charvec)

"""
    StabilizerCodeCSS(C::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)
    CSSCode(C::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)

Return a CSS code using the CSS construction on a self-orthogonal linear code
`C`, i.e., `C ⊆ C^⟂`.

Setting `C1 = C^⟂` and `C2 = C`, the resulting code has dimension `k = k1 - k2`
and minimum distance `d >= min(d1, d2^⟂)`. The `X` stabilizers are given by the
parity-check matrix of `C2^⟂`, `H(C2^⟂)`, and the `Z` stabilizers by `H(C1)`.

It is often desirable in quantum error correction to work with a set of
overcomplete stabilizers. Therefore this constructor does not simplify any
provided set of stabilizers. The dimension of the code is computed based on the
rank, and the user should not use the matrix dimension of the stabilizers to
determine such quantities. Use `isovercomplete` to determine if an
`AbstractStabilizerCode` is overcomplete.

# Arguments
* `C`: a self-orthogonal linear code
* `charvec`: a length `2n` vector with elements in the `Z/(2p)` if
  `chracteristic(field(C))` is 2 and `Z/(p)` otherwise. The first `n` elements
  specify the exponents of the `X` phases and second `n` the exponents of the
  `Z` phases; a missing argument will be set to the all-zero vector

# Notes
* A `+1` phase should be entered as `0` since the character vector stores the
  exponents.
* Stabilizer signs are automatically computed given the character vector.
"""
function StabilizerCodeCSS(C::LinearCode, charvec::Union{Vector{nmod}, Missing}=missing)
    # this should have Xstabs = Zstabs
    D = dual(C)
    C ⊆ D || error("The single code CSS construction requires C ⊆ C^⟂.")
    p = Int(characteristic(D.F))
    charvec = _processcharvec(charvec, p, D.n)

    # C2 ⊆ C1
    # k = k1 - k2
    # d >= minimum(d1, d2^⟂)
    # X - H(C2^⟂), Z - H(C1)
    S = directsum(D.H, D.H)
    Sq2 = symplectictoquadratic(S)
    logs = _logicals(S, directsum(D.G, D.G))

    # determine signs
    nr = nrows(D.H)
    signs, Xsigns, Zsigns = _determinesignsCSS(S, charvec, nr, nr)

    # q^n / p^k but rows is n - k
    rkS = rank(S)
    if rkS != D.n
        dimcode = BigInt(order(D.F))^D.n // BigInt(p)^rkS
        isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

        return StabilizerCodeCSS(D.F, base_ring(Sq2), D.n, dimcode, missing, missing, missing, Sq2,
            D.H, D.H, C, D, signs, Xsigns, Zsigns, logs, charvec,
            missing, missing, missing, false, missing)
    else
        return GraphStateStabilizerCSS(D.F, base_ring(Sq2), D.n, 0, missing, D.d, D.d, Sq2,
            D.H, D.H, C, D, signs, Xsigns, Zsigns, charvec, missing, false)
    end
end
CSSCode(C::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing) =
    StabilizerCodeCSS(C, charvec)

"""
    StabilizerCodeCSS(Xmatrix::fq_nmod_mat, Zmatrix::fq_nmod_mat, charvec::Union{Vector{nmod}, Missing}=missing)
    CSSCode(Xmatrix::fq_nmod_mat, Zmatrix::fq_nmod_mat, charvec::Union{Vector{nmod}, Missing}=missing)

Return a CSS code using the matrix `Xmatrix` as the `X` stabilizer matrix
and `Zmatrix` as the `Z` stabilizer matrix.

It is often desirable in quantum error correction to work with a set of
overcomplete stabilizers. Therefore this constructor does not simplify any
provided set of stabilizers. The dimension of the code is computed based on the
rank, and the user should not use the matrix dimension of the stabilizers to
determine such quantities. Use `isovercomplete` to determine if an
`AbstractStabilizerCode` is overcomplete.

# Arguments
* `Xmatrix`: a matrix over a finite field of type `FqNmodFiniteField`
* `Zmatrix`: a matrix over a finite field of type `FqNmodFiniteField`
* `charvec`: a length `2n` vector with elements in the `Z/(2p)` if
  `chracteristic(field(C1))` is 2 and `Z/(p)` otherwise. The first `n` elements
  specify the exponents of the `X` phases and second `n` the exponents of the
  `Z` phases; a missing argument will be set to the all-zero vector

# Notes
* A `+1` phase should be entered as `0` since the character vector stores the
  exponents.
* Stabilizer signs are automatically computed given the character vector.
* The orthogonality of the `X` and `Z` stabilizers are automatically checked and
  will error upon failure.
"""
function StabilizerCodeCSS(Xmatrix::fq_nmod_mat, Zmatrix::fq_nmod_mat,
    charvec::Union{Vector{nmod}, Missing}=missing)

    iszero(Xmatrix) && error("The `X` stabilizer matrix is empty.")
    iszero(Zmatrix) && error("The `Z` stabilizer matrix is empty.")
    n = ncols(Xmatrix)
    n ==  ncols(Zmatrix) || error("Both matrices must have the same length in the CSS construction.")
    F = base_ring(Xmatrix)
    F == base_ring(Zmatrix) || error("Both matrices must be over the same base field in the CSS construction.")
    iszero(Zmatrix * transpose(Xmatrix)) || error("The given matrices are not symplectic orthogonal.")
    p = Int(characteristic(F))
    charvec = _processcharvec(charvec, p, n)

    Xmatrix = _removeempty(Xmatrix, "rows")
    Zmatrix = _removeempty(Zmatrix, "rows")

    # determine if the provided set of stabilizers are redundant
    Xrank = rank(Xmatrix)
    Zrank = rank(Zmatrix)
    if nrows(Xmatrix) > Xrank || nrows(Zmatrix) > Zrank
        overcomp = true
    else
        overcomp = false
    end

    S = directsum(Xmatrix, Zmatrix)
    Sq2 = symplectictoquadratic(S)
    signs, Xsigns, Zsigns = _determinesignsCSS(S, charvec, nrows(Xmatrix), nrows(Zmatrix))

    # find generators for S^⟂
    # note the H here is transpose of the standard definition
    _, H = right_kernel(hcat(S[:, n + 1:end], -S[:, 1:n]))
    # n + (n - Srank)
    ncols(H) == 2 * n - Xrank - Zrank || error("Normalizer matrix is not size n + k.")
    logs = _logicals(S, transpose(H))

    # q^n / p^k but rows is n - k
    rkS = Xrank + Zrank
    if rkS != n
        dimcode = BigInt(order(F))^n // BigInt(p)^rkS
        isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

        return StabilizerCodeCSS(F, base_ring(Sq2), n, dimcode, missing, missing, missing, Sq2,
            Xmatrix, Zmatrix, missing, missing, signs, Xsigns, Zsigns, logs,
            charvec, missing, missing, missing, overcomp, missing)
    else
        return GraphStateStabilizerCSS(F, base_ring(Sq2), n, 0, missing, missing, missing, Sq2,
            Xmatrix, Zmatrix, missing, missing, signs, Xsigns, Zsigns, charvec, missing,
            overcomp)
    end
end
CSSCode(Xmatrix::fq_nmod_mat, Zmatrix::fq_nmod_mat, charvec::Union{Vector{nmod}, Missing}=missing) =
    StabilizerCodeCSS(Xmatrix, Zmatrix, charvec)

"""
    StabilizerCodeCSS(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}
    CSSCode(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}

Return a CSS code using the vector of Pauli strings `SPauli` as stabilizers.

It is often desirable in quantum error correction to work with a set of
overcomplete stabilizers. Therefore this constructor does not simplify any
provided set of stabilizers. The dimension of the code is computed based on the
rank, and the user should not use the matrix dimension of the stabilizers to
determine such quantities. Use `isovercomplete` to determine if an
`AbstractStabilizerCode` is overcomplete.

# Arguments
* `SPauli`: a vector of Strings or Char vectors containing the letters {I, X, Y, Z}
* `charvec`: a length `2n` vector with elements in the `Z/(2p)` if
  `chracteristic(field(C1))` is 2 and `Z/(p)` otherwise. The first `n` elements
  specify the exponents of the `X` phases and second `n` the exponents of the
  `Z` phases; a missing argument will be set to the all-zero vector

# Notes
* A `+1` phase should be entered as `0` since the character vector stores the
  exponents.
* Stabilizer signs are automatically computed given the character vector.
* The orthogonality of the `X` and `Z` stabilizers are automatically checked and
  will error upon failure.
* Any +/- 1 characters in front of each stabilizer are stripped. No check is done
  to make sure these signs agree with the ones computed using the character vector.
* Will error when the provided strings are not CSS.
"""
function StabilizerCodeCSS(SPauli::Vector{T}, charvec::Union{Vector{nmod},
    Missing}=missing) where T <: Union{String, Vector{Char}}

    S = _Paulistringtosymplectic(_processstrings(SPauli))
    iszero(S) && error("The processed Pauli strings returned a set of empty stabilizer generators.")
    S = _removeempty(S, "rows")
    # the reason we repeat here and not call another constructor is the else
    # statement at the bottom of this function
    # would also need to compute down to signs to call _isCSSsymplectic
    # which would allow us to call the other constructor
    aresymplecticorthogonal(S, S, true) || error("The given stabilizers are not symplectic orthogonal.")
    Sq2 = symplectictoquadratic(S)
    n = ncols(Sq2)

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
    logs = _logicals(S, transpose(H))

    # q^n / p^k but rows is n - k
    args = _isCSSsymplectic(S, signs, true)
    if args[1]
        if rkS != n
            dimcode = BigInt(order(F))^n // BigInt(p)^rkS
            isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

            return StabilizerCodeCSS(F, base_ring(Sq2), n, dimcode, missing, missing, missing,
                Sq2, args[2], args[4], missing, missing, signs, args[3], args[5],
                logs, charvec, missing, missing, missing, overcomp, missing)
        else
            return GraphStateStabilizerCSS(F, base_ring(Sq2), n, 0, missing, missing, missing,
                Sq2, args[2], args[4], missing, missing, signs, args[3], args[5],
                charvec, missing, overcomp)
        end
    else
        error("Provided Pauli strings are not CSS.")
    end
end
CSSCode(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}} =
 StabilizerCodeCSS(SPauli, charvec)

"""
    StabilizerCodeCSS(S::StabilizerCode)
    CSSCode(S::StabilizerCode)

Return the `[[2n, 2k, S.d <= d <= 2 S.d]]` CSS code derived by splitting the stabilizers of `S`.
"""
# TODO: wait, gotta be a typo here since it only has S.n columns
function StabilizerCodeCSS(S::StabilizerCode)
	X = S.stabs[:, 1:S.n]
	Z = S.stabs[:, S.n + 1:end]
	return StabilizerCodeCSS(X, Z, S.charvec)
end
CSSCode(S::StabilizerCode) = StabilizerCodeCSS(S)

# entanglement-assisted is not symplectic orthogonal
"""
    StabilizerCode(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}

Return a stabilizer code using the vector of Pauli strings `SPauli` as stabilizers.

It is often desirable in quantum error correction to work with a set of
overcomplete stabilizers. Therefore this constructor does not simplify any
provided set of stabilizers. The dimension of the code is computed based on the
rank, and the user should not use the matrix dimension of the stabilizers to
determine such quantities. Use `isovercomplete` to determine if an
`AbstractStabilizerCode` is overcomplete.

# Arguments
* `SPauli`: a vector of Strings or Char vectors containing the letters {I, X, Y, Z}
* `charvec`: a length `2n` vector with elements in the `Z/(2p)` if
  `chracteristic(field(C1))` is 2 and `Z/(p)` otherwise. The first `n` elements
  specify the exponents of the `X` phases and second `n` the exponents of the
  `Z` phases; a missing argument will be set to the all-zero vector

# Notes
* A `+1` phase should be entered as `0` since the character vector stores the
  exponents.
* Stabilizer signs are automatically computed given the character vector.
* The orthogonality of the stabilizers are automatically checked and will error
  upon failure.
* Any +/- 1 characters in front of each stabilizer are stripped. No check is done
  to make sure these signs agree with the ones computed using the character vector.
"""
function StabilizerCode(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}
    SPaulistripped = _processstrings(SPauli)
    S = _Paulistringtosymplectic(SPaulistripped)
    iszero(S) && error("The processed Pauli strings returned a set of empty stabilizer generators.")
    return StabilizerCode(S, true, charvec)
end

"""
    StabilizerCode(Sq2::fq_nmod_mat, symp::Bool=false, charvec::Union{Vector{nmod}, Missing}=missing)

Return a stabilizer code using the matrix `Sq2` as the stabilizer matrix.

The matrix `Sq2` is assumed to be an `n` column matrix over the quadratic extension.
If the optional parameter `symp` is set to `true`, `Sq2` is assumed to be in
symplectic form over the base field.

It is often desirable in quantum error correction to work with a set of
overcomplete stabilizers. Therefore this constructor does not simplify any
provided set of stabilizers. The dimension of the code is computed based on the
rank, and the user should not use the matrix dimension of the stabilizers to
determine such quantities. Use `isovercomplete` to determine if an
`AbstractStabilizerCode` is overcomplete.

# Arguments
* `Sq2`: a matrix over a finite field of type `FqNmodFiniteField`
* `symp`: a boolean
* `charvec`: a length `2n` vector with elements in the `Z/(2p)` if
  `chracteristic(field(C1))` is 2 and `Z/(p)` otherwise. The first `n` elements
  specify the exponents of the `X` phases and second `n` the exponents of the
  `Z` phases; a missing argument will be set to the all-zero vector

# Notes
* A `+1` phase should be entered as `0` since the character vector stores the
  exponents.
* Stabilizer signs are automatically computed given the character vector.
* The orthogonality of the stabilizers are automatically checked and will error
  upon failure.
"""
function StabilizerCode(Sq2::fq_nmod_mat, symp::Bool=false,
    charvec::Union{Vector{nmod}, Missing}=missing)

    iszero(Sq2) && error("The stabilizer matrix is empty.")
    Sq2 = _removeempty(Sq2, "rows")
    if symp
        S = Sq2
        # this will error properly if not correct
        Sq2 = symplectictoquadratic(Sq2)
    else
        iseven(degree(base_ring(Sq2))) || error("The base ring of the given matrix is not a quadratic extension.")
        S = quadratictosymplectic(Sq2)
    end
    aresymplecticorthogonal(S, S, true) || error("The given stabilizers are not symplectic orthogonal.")

    F = base_ring(S)
    p = Int(characteristic(F))
    n = ncols(Sq2)
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
    logs = _logicals(S, transpose(H))

    # q^n / p^k but rows is n - k
    dimcode = BigInt(order(F))^n // BigInt(p)^rkS
    isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

    args = _isCSSsymplectic(S, signs, true)
    if args[1]
        if rkS != n
            return StabilizerCodeCSS(F, base_ring(Sq2), n, dimcode, missing, missing, missing,
                Sq2, args[2], args[4], missing, missing, signs, args[3], args[5],
                logs, charvec, missing, missing, missing, overcomp, missing)
        else
            return GraphStateStabilizerCSS(F, base_ring(Sq2), n, 0, missing, missing, missing,
                Sq2, args[2], args[4], missing, missing, signs, args[3], args[5],
                charvec, missing, overcomp)
        end
    else
        if rkS != n
            return StabilizerCode(F, base_ring(Sq2), n, dimcode, missing, Sq2, logs,
                charvec, signs, missing, missing, missing, overcomp, missing)
        else
            return GraphState(F, base_ring(Sq2), n, 0, missing, Sq2, charvec, signs,
                missing, overcomp)
        end
    end
end

# slow? but works without permutations
function _logicals(stabs::fq_nmod_mat, dualgens::fq_nmod_mat)
    L = symplectictoquadratic(_quotientspace(dualgens, stabs))
    logs = _makepairs(L)
    # verify
    n = ncols(L)
    logsmat = vcat([vcat(logs[i]...) for i in 1:length(logs)]...)
    F = base_ring(stabs)
    Lsym = map_entries(x -> F(coeff(x, 0)), quadratictosymplectic(logsmat))
    aresymplecticorthogonal(stabs, Lsym, true) || error("Computed logicals do not commute with the codespace.")
    prod = hcat(Lsym[:, n + 1:end], -Lsym[:, 1:n]) * transpose(Lsym)
    sum(FpmattoJulia(prod), dims=1) == ones(Int, 1, size(prod, 1)) ||
        error("Computed logicals do not have the right commutation relations.")
    return logs
end

"""
    augment(S::AbstractStabilizerCode, row::fq_nmod_mat, symp::Bool=false, verbose::Bool=true)

Return the code created by added `row` to the stabilizers of `S`.

* Notes:
- The goal of this function is to track how the logical operators update given the new stabilizer.
  The unaffected logical operators are kept during the update and only those which don't commute
  with the new stabilizer are recomputed. Use `verbose` to better 
"""
# TODO: redo for subsystem codes and traits
function augment(S::AbstractStabilizerCode, row::fq_nmod_mat, symp::Bool=false, verbose::Bool=true)
    # typeof(S) ∈ [GraphState, GraphStateStabilizerCSS] && return S
    iszero(row) && return S
    nrows(row) == 1 || throw(ArgumentError("Only one stabilizer may be passed in at a time."))

    # it might be extremely problematic how I create fields throughout since they are pointer equality
    if symp
        ncols(row) == 2 * S.n || throw(ArgumentError("Symplectic flag set but row has incorrect number of columns."))
        base_ring(row) == S.F || throw(ArgumentError("Row must be over the same ring as the code."))
        rowq2 = symplectictoquadratic(row)
    else
        ncols(row) = S.n || throw(ArgumentError("Row has incorrect number of columns."))
        base_ring(row) == S.E || throw(ArgumentError("Row must be over the same ring as the code."))
        rowq2 = row
        row = quadratictosymplectic(rowq2)
    end

    # this is a more theoretically pure way to do this, probably slower than the iszero(prod) below
    # but this avoids possibly computing the logicals if they aren't already known and aren't needed
    symstabs = quadratictosymplectic(S.stabs)
    rankS = rank(symstabs)
    newstabs = vcat(S.stabs, rowq2)
    newsymstabs = vcat(symstabs, row)
    ranknewS = rank(newsymstabs)
    if rankS == ranknewS
        verbose && println("Row is already in the stabilizer group. Nothing to update.")    
        return S
    elseif S.k == 1
        verbose && println("Row is not in the stabilizer group; the result is a graph state.")
        return StabilizerCode(newstabs, false, S.charvec)
    end

    # not a stabilizer and multiple logical pairs
    logs = logicals(S)
    logsmat = logicalsmatrix(S)
    Lsym = quadratictosymplectic(logsmat)
    LsymEuc = hcat(Lsym[:, S.n + 1:end], -Lsym[:, 1:S.n])
    prod = LsymEuc * transpose(row)
    # if iszero(prod)
    #     verbose && println("Row is already in the stabilizer group. Nothing to update.")    
    #     return S
    # end
    numlogs = nrows(logsmat)
    logstokeep = Vector{Int}()
    logpairstokeep = Vector{Int}()
    # safer than 1:2k
    pair = 1
    for i in 1:2:numlogs
        if iszero(prod[i]) && iszero(prod[i + 1])
            append!(logstokeep, i, i + 1)
            append!(logpairstokeep, pair)
        end
        pair += 1
    end
    if isempty(logstokeep)
        verbose && println("Row does not commute with any logicals. The entire code needs to be recomputed from scratch.")
        Snew = StabilizerCode(newstabs, false, S.charvec)
        _ = logicals(Snew)
        return Snew
    end

    verbose && println("Logical pairs not requiring updating:")
    verbose && display(logs[logpairstokeep]) # how do I want to best display this?
    temp = quadratictosymplectic(vcat(newstabs, logsmat[logstokeep, :]))
    # kernel should contain newstabs and new logs but not old logs
    _, H = right_kernel(hcat(temp[:, S.n + 1:end], -temp[:, 1:S.n]))
    H = transpose(H)
    # unclear why this shouldn't be made not symplectic again
    temp = _quotientspace(H, quadratictosymplectic(newstabs))
    # temp = hcat(temp[:, S.n + 1:end], -temp[:, 1:S.n])
    newlogs = _makepairs(symplectictoquadratic(temp))
    # verify
    fulllogs = [logs[logpairstokeep]; newlogs]
    logsmat = vcat([vcat(fulllogs[i]...) for i in 1:length(fulllogs)]...)
    aresymplecticorthogonal(newstabs, logsmat) || error("Computed logicals do not commute with the codespace.")
    Lsym = quadratictosymplectic(logsmat);
    prod = hcat(Lsym[:, S.n + 1:end], -Lsym[:, 1:S.n]) * transpose(Lsym)
    sum(FpmattoJulia(prod), dims=1) == ones(Int, 1, size(prod, 1)) || error("Computed logicals do not have the right commutation relations.")
    # set and return if good
    verbose && println("New logicals:")
    verbose && display(newlogs)
    Snew = StabilizerCode(newstabs, false, S.charvec)
    Snew.logicals = fulllogs
    return Snew
end

"""
    expurgate(S::AbstractStabilizerCode, rows::Vector{Int}, verbose::Bool=true)

Return the code created by removing the stabilizers indexed by `rows`.

* Notes:
- The goal of this function is to track how the logical operators update through this process.
  Here, the original logical pairs are kept and an appropriate number of new pairs are added.
"""
# TODO: redo for subsystem codes and traits
function expurgate(S::AbstractStabilizerCode, rows::Vector{Int}, verbose::Bool=true)
    numstabs = nrows(S.stabs)
    rows ⊆ 1:numstabs || throw(ArgumentError("Argument rows not a subset of the number of stabilizers."))
    verbose && println("Removing stabilizers: $rows")
    newstabs = S.stabs[setdiff(1:numstabs, rows), :]
    Snew = StabilizerCode(newstabs, false, S.charvec)
    # if typeof(S) ∉ [GraphState, GraphStateStabilizerCSS]
    if true
        logs = logicals(S)
        logsmatrix = logicalsmatrix(S)
        small = quadratictosymplectic(vcat(newstabs, logsmatrix))
        # println(size(small))
        _, H = right_kernel(hcat(small[:, S.n + 1:end], -small[:, 1:S.n]))
        H = transpose(H)
        # println(size(H))
        # dualgenssym = hcat(H[:, S.n + 1:end], -H[:, 1:S.n])
        # println(size(dualgenssym))
        temp = _quotientspace(H, quadratictosymplectic(newstabs))
        # temp = hcat(temp[:, S.n + 1:end], -temp[:, 1:S.n])
        # using dualgenssym then switching temp here just switches {X, Z} to {Z, X}
        # but the vectors remain the same for some reason
        newlogs = _makepairs(symplectictoquadratic(temp))
        # verify
        fulllogs = [logs; newlogs]
        logsmatrix = vcat([vcat(fulllogs[i]...) for i in 1:length(fulllogs)]...)
        aresymplecticorthogonal(newstabs, logsmatrix) || error("Computed logicals do not commute with the codespace.")
        Lsym = quadratictosymplectic(logsmatrix);
        prod = hcat(Lsym[:, S.n + 1:end], -Lsym[:, 1:S.n]) * transpose(Lsym)
        sum(FpmattoJulia(prod), dims=1) == ones(Int, 1, size(prod, 1)) || error("Computed logicals do not have the right commutation relations.")
        # set and return if good
        verbose && println("New logicals:")
        verbose && display(newlogs)
        Snew.logicals = [fulllogs; newlogs]
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
