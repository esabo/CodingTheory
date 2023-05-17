# Copyright (c) 2022, 2023 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    HypergraphProductCode(A::CTMatrixTypes, B::CTMatrixTypes, charvec::Union{Vector{nmod}, Missing}=missing, logsalg::Symbol=:stndfrm)

Return the hypergraph product code of matrices `A` and `B` whose signs are determined by `charvec`.
"""
function HypergraphProductCode(A::CTMatrixTypes, B::CTMatrixTypes, charvec::Union{Vector{nmod},
    Missing}=missing, logsalg::Symbol=:stndfrm)

    logsalg ∈ [:stndfrm, :VS, :syseqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    F = base_ring(A)
    F == base_ring(B) || throw(ArgumentError("Matrices need to be over the same base ring"))

    # note that orthogonality of C1 and C2 here is not necessary because
    # HX * transpose(HZ) = C1.H \otimes transpose(C2.H) + C1.H \otimes transpose(C2.H) = 0
    # in characteristic 2
    Atr = transpose(A)
    Btr = transpose(B)
    # branch for speedup
    if Int(order(F)) == 2
        HX = hcat(A ⊗ identity_matrix(F, ncols(B)), identity_matrix(F, ncols(Atr)) ⊗ Btr)
    else
        HX = hcat(A ⊗ identity_matrix(F, ncols(B)), -identity_matrix(F, ncols(Atr)) ⊗ Btr)
    end
    HZ = hcat(identity_matrix(F, ncols(A)) ⊗ B, Atr ⊗ identity_matrix(F, ncols(Btr)))
    n = ncols(HX)
    stabs = directsum(HX, HZ)
    stabsstand, Pstand, standr, standk, rnk = _standardformstabilizer(stabs)
    if !iszero(standk)
        if logsalg == :stndfrm
            logs = _makepairs(_logicalsstandardform(stabsstand, n, standk, standr, Pstand))
            logsmat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
        else
            rnkH, H = right_kernel(hcat(stabs[:, n + 1:end], -stabs[:, 1:n]))
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
            # n + (n - Srank)
            nrows(Htr) == 2 * n - rnk || error("Normalizer matrix is not size n + k.")
            logs, logsmat = _logicals(stabs, Htr, logsalg)
        end
    end

    F = base_ring(stabs)
    p = Int(characteristic(F))
    charvec = _processcharvec(charvec, p, 2 * n)
    signs, Xsigns, Zsigns = _determinesignsCSS(stabs, charvec, nrows(HX), nrows(HZ))
    overcomp = nrows(stabs) > rnk

    # q^n / p^k but rows is n - k
    k = BigInt(order(F))^n // BigInt(p)^rnk
    isinteger(k) && (k = round(Int, log(BigInt(p), k));)
    # (ismissing(C1.d) || ismissing(C2.d)) ? d = missing : d = minimum([C1.d, C2.d])

    return HypergraphProductCode(F, n, k, missing, missing, missing, stabs, HX, HZ, missing,
        missing, signs, Xsigns, Zsigns, logs, logsmat, charvec, overcomp, stabsstand, standr,
        standk, Pstand, missing, missing)
end

"""
    HypergraphProductCode(C::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)

Return the (symmetric) hypergraph product code of `C` whose signs are determined by `charvec`.
"""
function HypergraphProductCode(C::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing,
    logsalg::Symbol=:stndfrm)

    S = HypergraphProductCode(paritycheckmatrix(C), paritycheckmatrix(C), charvec, logsalg)
    S.C1 = C
    S.C2 = C
    S.d = C.d
    return S
    # Htr = transpose(C.H)
    # eye = identity_matrix(C.F, C.n)
    # eyetr = identity_matrix(C.F, ncols(Htr))
    # # branch for speed
    # Int(order(C.F)) == 2 ? (HX = hcat(C.H ⊗ eye, eyetr ⊗ Htr)) : (HX = hcat(C.H ⊗ eye, -eyetr ⊗ Htr);)
    # HZ = hcat(eye ⊗ C.H, Htr ⊗ eyetr)
    # n = ncols(HX)
    # stabs = directsum(HX, HZ)
    # stabsstand, Pstand, standr, standk, rnk = _standardformstabilizer(stabs)
    # if !iszero(standk)
    #     if logsalg == :stndfrm
    #         logs = _makepairs(_logicalsstandardform(stabsstand, n, standk, standr, Pstand))
    #         logsmat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
    #     else
    #         rnkH, H = right_kernel(hcat(stabs[:, n + 1:end], -stabs[:, 1:n]))
    #         if ncols(H) == rnkH
    #             Htr = transpose(H)
    #         else
    #             # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
    #             nr = nrows(H)
    #             Htr = zero_matrix(base_ring(H), rnkH, nr)
    #             for r in 1:nr
    #                 for c in 1:rnkH
    #                     !iszero(H[r, c]) && (Htr[c, r] = H[r, c];)
    #                 end
    #             end
    #         end
    #         # n + (n - Srank)
    #         nrows(Htr) == 2 * n - rnk || error("Normalizer matrix is not size n + k.")
    #         logs, logsmat = _logicals(stabs, Htr, logsalg)
    #     end
    # end

    # F = base_ring(stabs)
    # p = Int(characteristic(F))
    # charvec = _processcharvec(charvec, p, 2 * n)
    # signs, Xsigns, Zsigns = _determinesignsCSS(stabs, charvec, nrows(HX), nrows(HZ))
    # overcomp = nrows(stabs) > rnk

    # # q^n / p^k but rows is n - k
    # dimcode = BigInt(order(F))^n // BigInt(p)^rnk
    # isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

    # return HypergraphProductCode(C.F, n, dimcode, C.d, missing, missing, stabs, HX, HZ, C, C,
    #     signs, Xsigns, Zsigns, logs, logsmat, charvec, overcomp, stabsstand, standr, standk, Pstand,
    #     missing, missing)
end

"""
    HypergraphProductCode(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)

Return the hypergraph product code of `C1` and `C2` whose signs are determined by `charvec`.
"""
function HypergraphProductCode(C1::AbstractLinearCode, C2::AbstractLinearCode,
    charvec::Union{Vector{nmod}, Missing}=missing, logsalg::Symbol=:stndfrm)

    S = HypergraphProductCode(paritycheckmatrix(C1), paritycheckmatrix(C2), charvec, logsalg)
    S.C1 = C1
    S.C2 = C1
    (ismissing(C1.d) || ismissing(C2.d)) ? (S.d = missing;) : (S.d = minimum([C1.d, C2.d]);)
    return S

    # logsalg ∈ [:stndfrm, :VS, :syseqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    # Int(order(C1.F)) == Int(order(C2.F)) || throw(ArgumentError("Codes need to be over the same base ring"))

    # # note that orthogonality of C1 and C2 here is not necessary because
    # # HX * transpose(HZ) = C1.H \otimes transpose(C2.H) + C1.H \otimes transpose(C2.H) = 0
    # # in characteristic 2
    # H1tr = transpose(C1.H)
    # H2tr = transpose(C2.H)
    # HX = hcat(H1 ⊗ identity_matrix(C2.F, C2.n), -identity_matrix(C1.F, ncols(H1tr)) ⊗ H2tr)
    # HZ = hcat(identity_matrix(C1.F, C1.n) ⊗ H2, H1tr ⊗ identity_matrix(C2.F, ncols(H2tr)))
    # n = ncols(HX)
    # stabs = directsum(HX, HZ)
    # stabsstand, Pstand, standr, standk, rnk = _standardformstabilizer(stabs)
    # if !iszero(standk)
    #     if logsalg == :stndfrm
    #         logs = _makepairs(_logicalsstandardform(stabsstand, n, standk, standr, Pstand))
    #         logsmat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
    #     else
    #         rnkH, H = right_kernel(hcat(stabs[:, n + 1:end], -stabs[:, 1:n]))
    #         if ncols(H) == rnkH
    #             Htr = transpose(H)
    #         else
    #             # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
    #             nr = nrows(H)
    #             Htr = zero_matrix(base_ring(H), rnkH, nr)
    #             for r in 1:nr
    #                 for c in 1:rnkH
    #                     !iszero(H[r, c]) && (Htr[c, r] = H[r, c];)
    #                 end
    #             end
    #         end
    #         # n + (n - Srank)
    #         nrows(Htr) == 2 * n - rnk || error("Normalizer matrix is not size n + k.")
    #         logs, logsmat = _logicals(stabs, Htr, logsalg)
    #     end
    # end

    # F = base_ring(stabs)
    # p = Int(characteristic(F))
    # charvec = _processcharvec(charvec, p, 2 * n)
    # signs, Xsigns, Zsigns = _determinesignsCSS(stabs, charvec, nrows(HX), nrows(HZ))
    # overcomp = nrows(stabs) > rnk

    # # q^n / p^k but rows is n - k
    # dimcode = BigInt(order(F))^n // BigInt(p)^rnk
    # isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)
    # (ismissing(C1.d) || ismissing(C2.d)) ? d = missing : d = minimum([C1.d, C2.d])

    # return HypergraphProductCode(C.F, n, dimcode, d, missing, missing, stabs, HX, HZ, C1, C2,
    #     signs, Xsigns, Zsigns, logs, logsmat, charvec, overcomp, stabsstand, standr, standk, Pstand,
    #     missing, missing)
end

# unable to yield quantum LDPC code families with non constant minimum distance
"""
    GeneralizedShorCode(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing, logsalg::Symbol=:stndfrm)
    BaconCasaccinoConstruction(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing, logsalg::Symbol=:stndfrm)

Return the generalized Shor code of `C1` and `C2` with `C1⟂ ⊆ C2` whose signs are determined by `charvec`.
"""
function GeneralizedShorCode(C1::AbstractLinearCode, C2::AbstractLinearCode,
    charvec::Union{Vector{nmod}, Missing}=missing, logsalg::Symbol=:stndfrm)

    logsalg ∈ [:stndfrm, :VS, :syseqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    Int(order(C1.F)) == 2 || error("Generalized Shor codes are only defined for binary codes.")
    Int(order(C2.F)) == 2 || error("Generalized Shor codes are only defined for binary codes.")
    dual(C1) ⊆ C2 || error("Generalized Shor codes require the dual of the first code is a subset of the second.")

    # HX stays sparse if H1 is but HZ does not if C1.d is large
    HX = paritycheckmatrix(C1) ⊗ identity_matrix(C2.F, C2.n)
    HZ = generatormatrix(C1) ⊗ paritycheckmatrix(C2)
    S = CSSCode(HX, HZ, charvec, logsalg)
    (ismissing(C1.d) || ismissing(C2.d)) ? (S.d = missing;) : (S.d = minimum([C1.d, C2.d]);)
    return S
end
BaconCasaccinoConstruction(C1::AbstractLinearCode, C2::AbstractLinearCode,
    charvec::Union{Vector{nmod}, Missing}=missing, logsalg::Symbol=:stndfrm) =
    GeneralizedShorCode(C1, C2, charvec, logsalg)

"""
    HyperBicycleCodeCSS(a::Vector{fq_nmod_mat}, b::Vector{fq_nmod_mat}, χ::Int, charvec::Union{Vector{nmod}, Missing}=missing, logsalg::Symbol=:stndfrm)

Return the hyperbicycle CSS code of `a` and `b` given `χ` whose signs are determined by `charvec`.

# Arguments
* a: A vector of length `c` of binary matrices of the same dimensions.
* b: A vector of length `c` of binary matrices of the same dimensions,
  potentially different from those of `a`.
* χ: A strictly positive integer coprime with `c`.
"""
function HyperBicycleCodeCSS(a::Vector{T}, b::Vector{T}, χ::Int,
    charvec::Union{Vector{nmod}, Missing}=missing, logsalg::Symbol=:stndfrm) where T <: CTMatrixTypes

    logsalg ∈ [:stndfrm, :VS, :syseqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    χ > 0 || throw(ArgumentError("Required χ > 0."))
    c = length(a)
    gcd(c, χ) == 1 || throw(ArgumentError("The length of the input vectors must be coprime with χ."))
    c == length(b) || throw(ArgumentError("Input vectors must have same length."))
    k1, n1 = size(a[1])
    k2, n2 = size(b[1])
    F = base_ring(a[1])
    Int(order(F)) == 2 || throw(ArgumentError("Hyperbicycle codes require binary inputs."))
    for i in 1:c
        F == base_ring(a[i]) || throw(ArgumentError("Inputs must share the same base ring."))
        (k1, n1) == size(a[i]) || throw(ArgumentError("First set of matrices must all have the same dimensions."))
        F == base_ring(b[i]) || throw(ArgumentError("Inputs must share the same base ring."))
        (k2, n2) == size(b[i]) || throw(ArgumentError("Second set of matrices must all have the same dimensions."))
    end

    # Julia creates a new scope for every iteration of the for loop,
    # so the else doesn't work without declaring these variables outside here
    H1::T
    H2::T
    HT1::T
    HT2::T
    for i in 1:c
        Sχi = zero_matrix(F, c, c)
        Ii = zero_matrix(F, c, c)
        for c1 in 1:c
            for r in 1:c
                if (c1 - r) % c == 1
                    Ii[r, c1] = F(1)
                end
                if (c1 - r) % c == (r - 1) * (χ - 1) % c
                    Sχi[r, c1] = F(1)
                end
            end
        end
        Iχi = Sχi * Ii
        ITχi = transpose(Sχi) * transpose(Ii)
        if i == 1
            H1 = Iχi ⊗ a[i]
            H2 = b[i] ⊗ Iχi
            HT1 = ITχi ⊗ transpose(a[i])
            HT2 = transpose(b[i]) ⊗ ITχi
        else
            H1 += Iχi ⊗ a[i]
            H2 += b[i] ⊗ Iχi
            HT1 += ITχi ⊗ transpose(a[i])
            HT2 += transpose(b[i]) ⊗ ITχi
        end
    end

    Ek1 = identity_matrix(F, k1)
    Ek2 = identity_matrix(F, k2)
    En1 = identity_matrix(F, n1)
    En2 = identity_matrix(F, n2)

    GX = hcat(Ek2 ⊗ H1, H2 ⊗ Ek1)
    GZ = hcat(HT2 ⊗ En1, En2 ⊗ HT1)
    return CSSCode(GX, GZ, charvec, logsalg)
    # equations 41 and 42 of the paper give matrices from which the logicals may be chosen
    # not really worth it, just use standard technique from StabilizerCode.jl
end

"""
    HyperBicycleCode(a::Vector{fq_nmod_mat}, b::Vector{fq_nmod_mat}, χ::Int, charvec::Union{Vector{nmod}, Missing}=missing, logsalg::Symbol=:stndfrm)

Return the hyperbicycle CSS code of `a` and `b` given `χ` whose signs are determined by `charvec`.

# Arguments
* a: A vector of length `c` of binary matrices of the same dimensions.
* b: A vector of length `c` of binary matrices of the same dimensions,
  potentially different from those of `a`.
* χ: A strictly positive integer coprime with `c`.
"""
function HyperBicycleCode(a::Vector{T}, b::Vector{T}, χ::Int,
    charvec::Union{Vector{nmod}, Missing}=missing, logsalg::Symbol=:stndfrm) where T <: CTMatrixTypes

    logsalg ∈ (:stndfrm, :VS, :syseqs) || throw(ArgumentError("Unrecognized logicals algorithm"))
    χ > 0 || throw(ArgumentError("Required χ > 0."))
    c = length(a)
    gcd(c, χ) == 1 || throw(ArgumentError("The length of the input vectors must be coprime with χ."))
    c == length(b) || throw(ArgumentError("Input vectors must have same length."))
    k1, n1 = size(a[1])
    k2, n2 = size(b[1])
    F = base_ring(a[1])
    Int(order(F)) == 2 || throw(ArgumentError("Hyperbicycle codes require binary inputs."))
    for i in 1:c
        F == base_ring(a[i]) || throw(ArgumentError("Inputs must share the same base ring."))
        (k1, n1) == size(a[i]) || throw(ArgumentError("First set of matrices must all have the same dimensions."))
        F == base_ring(b[i]) || throw(ArgumentError("Inputs must share the same base ring."))
        (k2, n2) == size(b[i]) || throw(ArgumentError("Second set of matrices must all have the same dimensions."))
    end

    for i in 1:c
        Sχi = zero_matrix(F, c, c)
        Ii = zero_matrix(F, c, c)
        for c1 in 1:c
            for r in 1:c
                if (c1 - r) % c == 1
                    Ii[r, c1] = F(1)
                end
                if (c1 - r) % c == (r - 1) * (χ - 1) % c
                    Sχi[r, c1] = F(1)
                end
            end
        end
        Iχi = Sχi * Ii
        if i == 1
            H1 = Iχi ⊗ a[i]
            H2 = b[i] ⊗ Iχi
        else
            H1 += Iχi ⊗ a[i]
            H2 += b[i] ⊗ Iχi
        end
    end

    Ek1 = identity_matrix(F, k1)
    Ek2 = identity_matrix(F, k2)
    stabs = hcat(Ek2 ⊗ H1, H2 ⊗ Ek1)
    return StabilizerCode(stabs, charvec, logsalg)
    # equations 41 and 42 of the paper give matrices from which the logicals may be chosen
    # not really worth it, just use standard technique from StabilizerCode.jl
end

"""
    GeneralizedBicycleCode(A::fq_nmod_mat, B::fq_nmod_mat, charvec::Union{Vector{nmod}, Missing}=missing, logsalg::Symbol=:stndfrm)

Return the generealized bicycle code given by `A` and `B` whose signs are determined by `charvec`.
"""
function GeneralizedBicycleCode(A::T, B::T,
    charvec::Union{Vector{nmod}, Missing}=missing, logsalg::Symbol=:stndfrm) where T <: CTMatrixTypes

    logsalg ∈ [:stndfrm, :VS, :syseqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    F = base_ring(A)
    F == base_ring(B) || throw(ArgumentError("Arguments must be over the same base ring."))
    (iszero(A) || iszero(B)) && throw(ArgumentError("Arguments should not be zero."))
    # this will take care of the sizes being square
    iszero(A * B - B * A) || throw(ArgumentError("Arguments must commute."))

    HX = hcat(A, B)
    # branch for speedup
    HZ = Int(order(F)) == 2 ? hcat(transpose(B), transpose(A)) : hcat(transpose(B), -transpose(A))
    return CSSCode(HX, HZ, charvec, logsalg)
end

"""
    GeneralizedBicycleCode(a::T, b::T, charvec::Union{Vector{nmod}, Missing}=missing, logsalg::Symbol=:stndfrm) where T <: ResElem

Return the generealized bicycle code determined by `a` and `b` whose signs are determined by `charvec`.

# Notes
* `l x l` circulant matrices are constructed using the coefficients of the polynomials
  `a` and `b` in `F_q[x]/(x^l - 1)` (`gcd(q, l) = 1`) as the first column
"""
function GeneralizedBicycleCode(a::T, b::T, charvec::Union{Vector{nmod}, Missing}=missing,
    logsalg::Symbol=:stndfrm) where T <: ResElem

    logsalg ∈ [:stndfrm, :VS, :syseqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    parent(a) == parent(b) || throw(ArgumentError("Both objects must be defined over the same residue ring."))

    return GeneralizedBicycleCode(polytocircmatrix(a), polytocircmatrix(b), charvec, logsalg)
end

# function BicycleCode(A::fq_nmot_mat)
#     m, n = size(A)
#     m == n || throw(ArgumentError("Input matrix must be square."))
#     # should probably check for F_2

#     H = hcat(A, transpose(A))
#     return CSSCode(H, H)
# end
    
"""
    GeneralizedHypergraphProductCodeMatrices(A::MatElem{T}, b::T) where T <: ResElem

Return the pre-lifted matrices `HX` and `HZ` of the generalized hypergraph product code of `A` and `b`.

# Arguments
* `A` - an `m x n` matrix with coefficents in `F_2[x]/(x^m - 1)`
* `b` - a polynomial over the same residue ring

# Notes
* Use `LiftedGeneralizedHypergraphProductCode` to return a quantum code over the base ring directly.
"""
function GeneralizedHypergraphProductCodeMatrices(A::MatElem{T}, b::T) where T <: ResElem

    @warn "Commutativity of A and b required but not yet enforced."
    S = base_ring(b)
    F = base_ring(S)
    Int(order(F)) == 2 || throw(ArgumentError("The generalized hypergraph product is only defined over GF(2)."))
    R = parent(A[1, 1])
    R == parent(b) || throw(ArgumentError("Both objects must be defined over the same residue ring."))
    m, n = size(A)
    (m != 1 && n != 1) || throw(ArgumentError("First input matrix must not be a vector."))
    f = modulus(R)
    l = degree(f)
    f == gen(S)^l - 1 || throw(ArgumentError("Residue ring not of the form x^l - 1."))
    # gcd(l, Int(characteristic(F))) == 1 || throw(ArgumentError("Residue ring over F_q[x] must be defined by x^l - 1 with gcd(l, q) = 1."))
    
    Atr = transpose(A)
    for c in 1:m
        for r in 1:n
            hcoeffs = collect(coefficients(Nemo.lift(Atr[r, c])))
            for _ in 1:l - length(hcoeffs)
                push!(hcoeffs, F(0))
            end
            hcoeffs[2:end] = reverse(hcoeffs[2:end])
            Atr[r, c] = R(S(hcoeffs))
        end
    end
    bcoeffs = collect(coefficients(Nemo.lift(b)))
    for _ in 1:l - length(bcoeffs)
        push!(bcoeffs, F(0))
    end
    bcoeffs[2:end] = reverse(bcoeffs[2:end])
    btr = R(S(bcoeffs))
    Mn = MatrixSpace(R, n, n)
    HZ = hcat(Mn(btr), Atr)
    Mm = MatrixSpace(R, m, m)
    HX = hcat(A, Mm(b))
    return HX, HZ
end

"""
    LiftedGeneralizedHypergraphProductCode(A::MatElem{T}, b::T, charvec::Union{Vector{nmod}, Missing}=missing, logsalg::Symbol=:stndfrm) where T <: ResElem

Return the lifted generalized hypergraph product code of `A` and `b`.

# Arguments
* `A` - an `m x n` matrix with coefficents in `F_2[x]/(x^m - 1)`
* `b` - a polynomial over the same residue ring
"""
function LiftedGeneralizedHypergraphProductCode(A::MatElem{T}, b::T, charvec::Union{Vector{nmod}, Missing}=missing,
    logsalg::Symbol=:stndfrm) where T <: ResElem

    HX, HZ = GeneralizedHypergraphProductCodeMatrices(A, b)
    return CSSCode(lift(HX), lift(HZ), charvec, logsalg)
end

"""
    QuasiCyclicLiftedProductCodeMatrices(A::MatElem{T}, B::MatElem{T}) where T <: ResElem

Return the pre-lifted matrices `HX` and `HZ` for the lifted quasi-cyclic lifted product code.

# Arguments
* `A` - an `m x n1` matrix with coefficents in `F_2[x]/(x^m - 1)`
* `B` - an `m x n2` matrix with coefficents in the same residue ring

# Notes
* Use `QuasiCyclicLiftedProductCode` to return a quantum code over the base ring directly.
"""
function QuasiCyclicLiftedProductCodeMatrices(A::MatElem{T}, B::MatElem{T}) where T <: ResElem
    @warn "Commutativity of A and b required but not yet enforced."
    S = base_ring(A[1, 1])
    F = base_ring(S)
    Int(order(F)) == 2 || throw(ArgumentError("The quasi-cyclic lifted product is only defined over GF(2)."))
    R = parent(A[1, 1])
    R == parent(B[1, 1]) || throw(ArgumentError("Both objects must be defined over the same residue ring."))
    f = modulus(R)
    l = degree(f)
    f == gen(S)^l - 1 || throw(ArgumentError("Residue ring not of the form x^l - 1."))
    
    k1, n1 = size(A)
    Atr = transpose(A)
    for c in 1:k1
        for r in 1:n1
            hcoeffs = collect(coefficients(Nemo.lift(Atr[r, c])))
            for _ in 1:l - length(hcoeffs)
                push!(hcoeffs, F(0))
            end
            hcoeffs[2:end] = reverse(hcoeffs[2:end])
            Atr[r, c] = R(S(hcoeffs))
        end
    end

    k2, n2 = size(B)
    Btr = transpose(B)
    for c in 1:k2
        for r in 1:n2
            hcoeffs = collect(coefficients(Nemo.lift(Btr[r, c])))
            for _ in 1:l - length(hcoeffs)
                push!(hcoeffs, F(0))
            end
            hcoeffs[2:end] = reverse(hcoeffs[2:end])
            Btr[r, c] = R(S(hcoeffs))
        end
    end

    Ek1 = identity_matrix(R, k1)
    Ek2 = identity_matrix(R, k2)
    En1 = identity_matrix(R, n1)
    En2 = identity_matrix(R, n2)

    HX = hcat(kronecker_product(A, Ek2), kronecker_product(Ek1, B))
    HZ = hcat(kronecker_product(En1, Btr), kronecker_product(Atr, En2))
    return HX, HZ
end

"""
    QuasiCyclicLiftedProductCode(A::MatElem{T}, B::MatElem{T}, charvec::Union{Vector{nmod}, Missing}=missing, logsalg::Symbol=:stndfrm) where T <: ResElem

Return the quasi-cyclic lifted product code given by the matrices `A` and `B` and whose signs are determined by `charvec`.
"""
function QuasiCyclicLiftedProductCode(A::MatElem{T}, B::MatElem{T}, charvec::Union{Vector{nmod}, Missing}=missing,
    logsalg::Symbol=:stndfrm) where T <: ResElem

    HX, HZ = QuasiCyclicLiftedProductCodeMatrices(A, B)
    return CSSCode(lift(HX), lift(HZ), charvec, logsalg)
end

"""
    BiasTailoredQuasiCyclicLiftedProductCodeMatrices(A::MatElem{T}, B::MatElem{T}) where T <: ResElem

Return the pre-lifted stabilizer matrix for bias-tailored lifted product code of `A` and `B`.

# Arguments
* `A` - an `m x n1` matrix with coefficents in `F_2[x]/(x^m - 1)`
* `B` - an `m x n2` matrix with coefficents in the same residue ring

# Notes
* Use `BiasTailoredQuasiCyclicLiftedProductCode` to return a quantum code over the base ring directly.
"""
function BiasTailoredQuasiCyclicLiftedProductCodeMatrices(A::MatElem{T}, B::MatElem{T}) where T <: ResElem
    @warn "Commutativity of A and b required but not yet enforced."
    S = base_ring(A[1, 1])
    F = base_ring(S)
    Int(order(F)) == 2 || throw(ArgumentError("The quasi-cyclic lifted product is only defined over GF(2)."))
    R = parent(A[1, 1])
    R == parent(B[1, 1]) || throw(ArgumentError("Both objects must be defined over the same residue ring."))
    f = modulus(R)
    l = degree(f)
    f == gen(S)^l - 1 || throw(ArgumentError("Residue ring not of the form x^l - 1."))
    
    k1, n1 = size(A)
    Atr = transpose(A)
    for c in 1:k1
        for r in 1:n1
            hcoeffs = collect(coefficients(Nemo.lift(Atr[r, c])))
            for _ in 1:l - length(hcoeffs)
                push!(hcoeffs, F(0))
            end
            hcoeffs[2:end] = reverse(hcoeffs[2:end])
            Atr[r, c] = R(S(hcoeffs))
        end
    end

    k2, n2 = size(B)
    Btr = transpose(B)
    for c in 1:k2
        for r in 1:n2
            hcoeffs = collect(coefficients(Nemo.lift(Btr[r, c])))
            for _ in 1:l - length(hcoeffs)
                push!(hcoeffs, F(0))
            end
            hcoeffs[2:end] = reverse(hcoeffs[2:end])
            Btr[r, c] = R(S(hcoeffs))
        end
    end

    Ek1 = identity_matrix(R, k1)
    Ek2 = identity_matrix(R, k2)
    En1 = identity_matrix(R, n1)
    En2 = identity_matrix(R, n2)

    A12 = kronecker_product(Atr, Ek2)
    A13 = kronecker_product(En1, B)
    A21 = kronecker_product(A, En2)
    A24 = kronecker_product(Ek1, Btr)
    return vcat(hcat(zeros(A21), A12, A13, zeros(A24)), hcat(A21, zeros(A12), zeros(A13), A24))
end

"""
    BiasTailoredQuasiCyclicLiftedProductCodeMatrices(A::MatElem{T}, B::MatElem{T}, charvec::Union{Vector{nmod}, Missing}=missing, logsalg::Symbol=:stndfrm) where T <: ResElem

Return the bias-tailored lifted product code of `A` and `B` whose signs are given by `charvec`.

# Arguments
* `A` - an `m x n1` matrix with coefficents in `F_2[x]/(x^m - 1)`
* `B` - an `m x n2` matrix with coefficents in the same residue ring
"""
function BiasTailoredQuasiCyclicLiftedProductCode(A::MatElem{T}, B::MatElem{T}, charvec::Union{Vector{nmod}, Missing}=missing,
    logsalg::Symbol=:stndfrm) where T <: ResElem

    stabs = BiasTailoredQuasiCyclicLiftedProductCodeMatrices(A, B)
    return StabilizerCode(lift(stabs), charvec, logsalg)
end

#############################
      # getter functions
#############################

#############################
      # setter functions
#############################

#############################
     # general functions
#############################
