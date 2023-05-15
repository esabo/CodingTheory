# Copyright (c) 2021, 2022, 2023 Eric Sabo, Benjamin Ide
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
function StabilizerCodeCSS(C1::AbstractLinearCode, C2::AbstractLinearCode,
    charvec::Union{Vector{nmod}, Missing}=missing, logsalg::Symbol=:stndfrm)

    logsalg ∈ [:stndfrm, :VS, :syseqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    C2 ⊆ C1 || throw(ArgumentError("The second argument must be a subset of the first in the CSS construction."))
    p = Int(characteristic(C1.F))
    charvec = _processcharvec(charvec, p, C1.n)

    # C2 ⊆ C1
    # k = k1 - k2
    # d >= minimum(d1, d2^⟂)
    # X - H(C2^⟂), Z - H(C1)
    D2 = dual(C2)
    stabs = directsum(D2.H, C1.H)
    stabsstand, Pstand, standr, standk, rnk = _standardformstabilizer(stabs)
    if !iszero(standk)
        if logsalg == :stndfrm
            logs = _makepairs(_logicalsstandardform(stabsstand, C1.n, standk, standr, Pstand))
            logsmat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
        else
            logs, logsmat = _logicals(stabs, directsum(C1.G, D2.G), logsalg)
        end
    end
    signs, Xsigns, Zsigns = _determinesignsCSS(stabs, charvec, nrows(D2.H), nrows(C1.H))

    # q^n / p^k but rows is n - k
    if !iszero(standk)
        dimcode = BigInt(order(C1.F))^C1.n // BigInt(p)^rnk
        isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

        return StabilizerCodeCSS(C1.F, C1.n, dimcode, missing, missing, missing, stabs, D2.H, C1.H,
            C2, C1, signs, Xsigns, Zsigns, logs, logsmat, charvec, missing, missing, missing,
            false, missing, stabsstand, standr, standk, Pstand)
    else
        return GraphStateStabilizerCSS(C1.F, C1.n, 0, missing, D2.d, C1.d, stabs, D2.H, C1.H, C2, C1,
            signs, Xsigns, Zsigns, charvec, missing, false, stabsstand, standr, standk, Pstand)
    end
end
CSSCode(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing,
    logsalg::Symbol=:stndfrm) = StabilizerCodeCSS(C1, C2, charvec, logsalg)

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
function StabilizerCodeCSS(C::LinearCode, charvec::Union{Vector{nmod}, Missing}=missing,
    logsalg::Symbol=:stndfrm)

    logsalg ∈ [:stndfrm, :VS, :syseqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    # this should have Xstabs = Zstabs
    D = dual(C)
    C ⊆ D || throw(ArgumentError("The single code CSS construction requires C ⊆ C^⟂."))
    p = Int(characteristic(D.F))
    charvec = _processcharvec(charvec, p, D.n)

    # C2 ⊆ C1
    # k = k1 - k2
    # d >= minimum(d1, d2^⟂)
    # X - H(C2^⟂), Z - H(C1)
    stabs = directsum(D.H, D.H)
    stabsstand, Pstand, standr, standk, rnk = _standardformstabilizer(stabs)
    if !iszero(standk)
        if logsalg == :stndfrm
            logs = _makepairs(_logicalsstandardform(stabsstand, C.n, standk, standr, Pstand))
            logsmat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
        else
            logs, logsmat = _logicals(stabs, directsum(D.G, D.G), logsalg)
        end
    end
    nr = nrows(D.H)
    signs, Xsigns, Zsigns = _determinesignsCSS(stabs, charvec, nr, nr)

    # q^n / p^k but rows is n - k
    if !iszero(standk)
        dimcode = BigInt(order(D.F))^D.n // BigInt(p)^rnk
        isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

        return StabilizerCodeCSS(D.F, D.n, dimcode, missing, missing, missing, stabs, D.H, D.H, C,
            D, signs, Xsigns, Zsigns, logs, logsmat, charvec, missing, missing, missing, false,
            missing, stabsstand, standr, standk, Pstand)
    else
        return GraphStateStabilizerCSS(D.F, D.n, 0, missing, D.d, D.d, stabs, D.H, D.H, C, D, signs,
            Xsigns, Zsigns, charvec, missing, false, stabsstand, standr, standk, Pstand)
    end
end
CSSCode(C::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing,
    logsalg::Symbol=:stndfrm) = StabilizerCodeCSS(C, charvec, logsalg)

"""
    StabilizerCodeCSS(Xmatrix::fq_nmod_mat, Zmatrix::fq_nmod_mat, charvec::Union{Vector{nmod}, Missing}=missing)
    CSSCode(Xmatrix::fq_nmod_mat, Zmatrix::fq_nmod_mat, charvec::Union{Vector{nmod}, Missing}=missing)

Return a CSS code whose `X`-stabilizers are given by `Xmatrix`, `Z`-stabilizers by `Zmatrix`, and signs by `charvec`.
"""
function StabilizerCodeCSS(Xmatrix::T, Zmatrix::T, charvec::Union{Vector{nmod}, Missing}=missing,
    logsalg::Symbol=:stndfrm) where T <: CTMatrixTypes

    logsalg ∈ [:stndfrm, :VS, :syseqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
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
    nrows(Xmatrix) > Xrank || nrows(Zmatrix) > Zrank ? (overcomp = true) : (overcomp = false)
    stabs = directsum(Xmatrix, Zmatrix)

    stabsstand, Pstand, standr, standk, rnk = _standardformstabilizer(stabs)
    if !iszero(standk)
        if logsalg == :stndfrm
            logs = _makepairs(_logicalsstandardform(stabsstand, n, standk, standr, Pstand))
            logsmat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
        else
            # find generators for S^⟂
            # note the H here is transpose of the standard definition
            _, H = right_kernel(hcat(stabs[:, n + 1:end], -stabs[:, 1:n]))
            # remove empty for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
            H = _removeempty(transpose(H), :rows)
            # n + (n - Srank)
            nrows(H) == 2 * n - Xrank - Zrank || error("Normalizer matrix is not size n + k.")
            logs, logsmat = _logicals(stabs, H)
        end
    end
    signs, Xsigns, Zsigns = _determinesignsCSS(stabs, charvec, nrows(Xmatrix), nrows(Zmatrix))

    # q^n / p^k but rows is n - k
    if !iszero(standk)
        dimcode = BigInt(order(F))^n // BigInt(p)^rnk
        isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

        return StabilizerCodeCSS(F, n, dimcode, missing, missing, missing, stabs, Xmatrix, Zmatrix,
            missing, missing, signs, Xsigns, Zsigns, logs, logsmat, charvec, missing, missing,
            missing, overcomp, missing, stabsstand, standr, standk, Pstand)
    else
        return GraphStateStabilizerCSS(F, n, 0, missing, missing, missing, stabs, Xmatrix, Zmatrix,
            missing, missing, signs, Xsigns, Zsigns, charvec, missing, overcomp, stabsstand, standr,
            standk, Pstand)
    end
end
CSSCode(Xmatrix::T, Zmatrix::T, charvec::Union{Vector{nmod}, Missing}=missing,
    logsalg::Symbol=:stndfrm) where T <: CTMatrixTypes = StabilizerCodeCSS(Xmatrix, Zmatrix, charvec, logsalg)

"""
    StabilizerCodeCSS(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}
    CSSCode(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}

Return the CSS code whose stabilizers are determined by the vector of Pauli strings `SPauli` and signs by `charvec`.

# Notes
* Any +/- 1 characters in front of each stabilizer are stripped. No check is done
  to make sure these signs agree with the ones computed using the character vector.
"""
function StabilizerCodeCSS(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing,
    logsalg::Symbol=:stndfrm) where T <: Union{String, Vector{Char}}

    logsalg ∈ [:stndfrm, :VS, :syseqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    stabs = _Paulistringtosymplectic(_processstrings(SPauli))
    iszero(stabs) && throw(ArgumentError("The processed Pauli strings returned a set of empty stabilizer generators."))
    stabs = _removeempty(stabs, :rows)
    # the reason we repeat here and not call another constructor is the else
    # statement at the bottom of this function
    # would also need to compute down to signs to call _isCSSsymplectic
    # which would allow us to call the other constructor
    aresymplecticorthogonal(stabs, stabs) || throw(ArgumentError("The given stabilizers are not symplectic orthogonal."))
    
    n = div(ncols(stabs), 2)
    stabsstand, Pstand, standr, standk, rnk = _standardformstabilizer(stabs)
    if !iszero(standk)
        if logsalg == :stndfrm
            logs = _makepairs(_logicalsstandardform(stabsstand, n, standk, standr, Pstand))
            logsmat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
        else
            # find generators for S^⟂
            # note the H here is transpose of the standard definition
            _, H = right_kernel(hcat(stabs[:, n + 1:end], -stabs[:, 1:n]))
            # remove empty for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
            H = _removeempty(transpose(H), :rows)
            # n + (n - Srank)
            nrows(H) == 2 * n - Xrank - Zrank || error("Normalizer matrix is not size n + k.")
            logs, logsmat = _logicals(stabs, H)
        end
    end

    F = base_ring(stabs)
    p = Int(characteristic(F))
    charvec = _processcharvec(charvec, p, 2 * n)
    signs = _determinesigns(stabs, charvec)
    nrows(stabs) > rnk ? (overcomp = true) : (overcomp = false)

    # q^n / p^k but rows is n - k
    args = _isCSSsymplectic(stabs, signs, true)
    if args[1]
        if !iszero(standk)
            dimcode = BigInt(order(F))^n // BigInt(p)^rkS
            isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

            return StabilizerCodeCSS(F, n, dimcode, missing, missing, missing, stabs, args[2],
                args[4], missing, missing, signs, args[3], args[5], logs, logsmat, charvec,
                missing, missing, missing, overcomp, missing, stabsstand, standr, standk, Pstand)
        else
            return GraphStateStabilizerCSS(F, n, 0, missing, missing, missing, stabs, args[2],
                args[4], missing, missing, signs, args[3], args[5], charvec, missing, overcomp, stabsstand,
                standr, standk, Pstand)
        end
    else
        error("Provided Pauli strings are not CSS.")
    end
end
CSSCode(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing,
    logsalg::Symbol=:stndfrm) where T <: Union{String, Vector{Char}} = StabilizerCodeCSS(SPauli, charvec, logsalg)

"""
    StabilizerCodeCSS(S::StabilizerCode)
    CSSCode(S::StabilizerCode)

Return the `[[2n, 2k, S.d <= d <= 2 S.d]]` CSS code derived by splitting the stabilizers of `S`.
"""
function StabilizerCodeCSS(S::StabilizerCode, logsalg::Symbol=:stndfrm)
	X = S.stabs[:, 1:S.n]
	Z = S.stabs[:, S.n + 1:end]
	return StabilizerCodeCSS(X, Z, S.charvec, logsalg)
end
CSSCode(S::StabilizerCode, logsalg::Symbol=:stndfrm) = StabilizerCodeCSS(S, logsalg)

# entanglement-assisted is not symplectic orthogonal
"""
    StabilizerCode(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}

Return the stabilizer code whose stabilizers are determined by the vector of Pauli strings `SPauli` and signs by `charvec`.

# Notes
* Any +/- 1 characters in front of each stabilizer are stripped. No check is done
  to make sure these signs agree with the ones computed using the character vector.
"""
function StabilizerCode(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing,
    logsalg::Symbol=:stndfrm) where T <: Union{String, Vector{Char}}

    SPaulistripped = _processstrings(SPauli)
    stabs = _Paulistringtosymplectic(SPaulistripped)
    iszero(stabs) && throw(ArgumentError("The processed Pauli strings returned a set of empty stabilizer generators."))
    return StabilizerCode(stabs, charvec, logsalg)
end

"""
    StabilizerCode(S::fq_nmod_mat, charvec::Union{Vector{nmod}, Missing}=missing)

Return the stabilizer code whose stabilizers is determined by `S` and signs by `charvec`.
"""
function StabilizerCode(stabs::CTMatrixTypes, charvec::Union{Vector{nmod}, Missing}=missing,
    logsalg::Symbol=:stndfrm)

    logsalg ∈ [:stndfrm, :VS, :syseqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    iszero(stabs) && throw(ArgumentError("The stabilizer matrix is empty."))
    stabs = _removeempty(stabs, :rows)
    aresymplecticorthogonal(stabs, stabs) || throw(ArgumentError("The given stabilizers are not symplectic orthogonal."))
    
    n = div(ncols(stabs), 2)
    stabsstand, Pstand, standr, standk, rnk = _standardformstabilizer(stabs)
    if !iszero(standk)
        if logsalg == :stndfrm
            logs = _makepairs(_logicalsstandardform(stabsstand, n, standk, standr, Pstand))
            logsmat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
        else
            _, H = right_kernel(hcat(stabs[:, n + 1:end], -stabs[:, 1:n]))
            # remove empty for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
            H = _removeempty(transpose(H), :rows)
            # n + (n - Srank)
            nrows(H) == 2 * n - rnk || error("Normalizer matrix is not size n + k.")
            logs, logsmat = _logicals(stabs, H)
        end
    end

    F = base_ring(stabs)
    p = Int(characteristic(F))
    charvec = _processcharvec(charvec, p, n)
    signs = _determinesigns(stabs, charvec)
    nrows(stabs) > rnk ? (overcomp = true) : (overcomp = false)

    # q^n / p^k but rows is n - k
    dimcode = BigInt(order(F))^n // BigInt(p)^rnk
    isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

    args = _isCSSsymplectic(stabs, signs, true)
    if args[1]
        if !iszero(standk)
            return StabilizerCodeCSS(F, n, dimcode, missing, missing, missing, stabs, args[2],
                args[4], missing, missing, signs, args[3], args[5], logs, logsmat, charvec,
                missing, missing, missing, overcomp, missing, stabsstand, standr, standk, Pstand)
        else
            return GraphStateStabilizerCSS(F, n, 0, missing, missing, missing, stabs, args[2],
                args[4], missing, missing, signs, args[3], args[5], charvec, missing, overcomp, stabsstand,
                standr, standk, Pstand)
        end
    else
        if !iszero(standk)
            return StabilizerCode(F, n, dimcode, missing, stabs, logs, logsmat, charvec, signs, missing,
                missing, missing, overcomp, missing, stabsstand, standr, standk, Pstand)
        else
            return GraphStateStabilizer(F, n, 0, missing, stabs, charvec, signs, missing, overcomp, stabsstand,
                standr, standk, Pstand)
        end
    end
end

function _logicals(stabs::T, dualgens::T, logsalgs::Symbol=:syseqs) where T <: CTMatrixTypes
    logsalg ∈ [:syseqs, :VS] || throw(ArgumentError("Unrecognized logicals algorithm"))

    L = _quotientspace(dualgens, stabs, logsalgs)
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
