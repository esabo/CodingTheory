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

    # new stuff from Ben
    # TODO: replace logicals above with getting logicals from here
    standardform, perms, standr = _standardformstabilizer(S)

    # determine signs
    signs, Xsigns, Zsigns = _determinesignsCSS(S, charvec, nrows(D2.H), nrows(C1.H))

    # q^n / p^k but rows is n - k
    rkS = rank(S)
    if rkS != C1.n
        dimcode = BigInt(order(C1.F))^C1.n // BigInt(p)^rkS
        isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

        return StabilizerCodeCSS(C1.F, C1.n, dimcode, missing, missing, missing, S, D2.H, C1.H,
            C2, C1, signs, Xsigns, Zsigns, logs, logsmat, charvec, missing, missing, missing,
            false, missing, standardform, standr, perms)
    else
        return GraphStateStabilizerCSS(C1.F, C1.n, 0, missing, D2.d, C1.d, S, D2.H, C1.H, C2, C1,
            signs, Xsigns, Zsigns, charvec, missing, false, standardform, standr, perms)
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

    # new stuff from Ben
    # TODO: replace logicals above with getting logicals from here
    standardform, perms, standr = _standardformstabilizer(S)

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
            missing, standardform, standr, perms)
    else
        return GraphStateStabilizerCSS(D.F, D.n, 0, missing, D.d, D.d, S, D.H, D.H, C, D, signs,
            Xsigns, Zsigns, charvec, missing, false, standardform, standr, perms)
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

    # new stuff from Ben
    # TODO: replace logicals above with getting logicals from here
    standardform, perms, standr = _standardformstabilizer(S)

    # q^n / p^k but rows is n - k
    rkS = Xrank + Zrank
    if rkS != n
        dimcode = BigInt(order(F))^n // BigInt(p)^rkS
        isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

        return StabilizerCodeCSS(F, n, dimcode, missing, missing, missing, S, Xmatrix, Zmatrix,
            missing, missing, signs, Xsigns, Zsigns, logs, logsmat, charvec, missing, missing,
            missing, overcomp, missing, standardform, standr, perms)
    else
        return GraphStateStabilizerCSS(F, n, 0, missing, missing, missing, S, Xmatrix, Zmatrix,
            missing, missing, signs, Xsigns, Zsigns, charvec, missing, overcomp, standardform, standr, perms)
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

    # new stuff from Ben
    # TODO: replace logicals above with getting logicals from here
    standardform, perms, standr = _standardformstabilizer(S)

    # q^n / p^k but rows is n - k
    args = _isCSSsymplectic(S, signs, true)
    if args[1]
        if rkS != n
            dimcode = BigInt(order(F))^n // BigInt(p)^rkS
            isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

            return StabilizerCodeCSS(F, n, dimcode, missing, missing, missing, S, args[2],
                args[4], missing, missing, signs, args[3], args[5], logs, logsmat, charvec,
                missing, missing, missing, overcomp, missing, standardform, standr, perms)
        else
            return GraphStateStabilizerCSS(F, n, 0, missing, missing, missing, S, args[2],
                args[4], missing, missing, signs, args[3], args[5], charvec, missing, overcomp, standardform, standr, perms)
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

    # new stuff from Ben
    # TODO: replace logicals above with getting logicals from here
    standardform, perms, standr = _standardformstabilizer(S)

    # q^n / p^k but rows is n - k
    dimcode = BigInt(order(F))^n // BigInt(p)^rkS
    isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

    args = _isCSSsymplectic(S, signs, true)
    if args[1]
        if rkS != n
            return StabilizerCodeCSS(F, n, dimcode, missing, missing, missing, S, args[2],
                args[4], missing, missing, signs, args[3], args[5], logs, logsmat, charvec,
                missing, missing, missing, overcomp, missing, standardform, standr, perms)
        else
            return GraphStateStabilizerCSS(F, n, 0, missing, missing, missing, S, args[2],
                args[4], missing, missing, signs, args[3], args[5], charvec, missing, overcomp, standardform, standr, perms)
        end
    else
        if rkS != n
            return StabilizerCode(F, n, dimcode, missing, S, logs, logsmat, charvec, signs, missing,
                missing, missing, overcomp, missing, standardform, standr, perms)
        else
            return GraphState(F, n, 0, missing, S, charvec, signs, missing, overcomp, standardform, standr, perms)
        end
    end
end

function _logicals(stabs::fq_nmod_mat, dualgens::fq_nmod_mat)
    L = _quotientspace(dualgens, stabs)
    logs = _makepairs(L)
    # verify
    n = div(ncols(L), 2)
    logsmat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
    aresymplecticorthogonal(stabs, logsmat) || error("Computed logicals do not commute with the codespace.")
    prod = hcat(logsmat[:, n + 1:end], -logsmat[:, 1:n]) * transpose(logsmat)
    sum(FpmattoJulia(prod), dims=1) == ones(Int, 1, size(prod, 1)) ||
        error("Computed logicals do not have the right commutation relations.")
    return logs, logsmat
end

function _standardformstabilizer(M::CTMatrixTypes)
    @assert iseven(size(M, 2))

    S = deepcopy(M)

    # If the stabilizer is overdetermined, remove unnecessary rows
    _rref_no_col_swap!(S, 1:size(S, 1), 1:size(S, 2))
    nr = size(S, 1)
    for i in size(S, 1):-1:1
        nr = i
        iszero(S[i, :]) || break
    end
    if nr != size(S, 1)
        S = S[1:nr, :]
    end

    n = div(size(S, 2), 2)
    k = n - nr

    # put S in standard form
    r, P1 = _rref_col_swap!(S, 1:nr, 1:n)
    _, P2 = _rref_col_swap!(S, (r + 1):nr, (n + r + 1):2n)

    P = if ismissing(P1) && ismissing(P2)
        missing
    elseif ismissing(P1)
        P2
    elseif ismissing(P2)
        P1
    else
        P1 * P2
    end

    return S, P, r
end

function _logicalsstandardform(C::AbstractSubsystemCode)
    n = C.n
    k = C.k
    r = C.standr
    S = C.standardform
    R = base_ring(S)
    logs = zero_matrix(R, 2k, 2n)
    E = S[(r + 1):size(S,1), (2n - k + 1):2n]
    C1 = S[1:r, (n + r + 1):(2n - k)]
    C1E = C1 * E
    for i in 1:k
        for j in 1:(n - k - r)
            # E^T
            logs[i, j + r] = S[r + j, 2n - k + i]
        end
        for j in 1:k
            # I in a couple of places
            logs[i, n - k - r + j] = one(R)
            logs[k + i, 2n - k + j] = one(R)
        end
        for j in 1:r
            # E^T * C1^T + C2^T
            logs[i, n + j] = C1E[j, i] + S[j, 2n - k + i]

            # A2^T
            logs[k + i, n + j] = S[j, n - k + i]
        end
    end

    return logs
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
