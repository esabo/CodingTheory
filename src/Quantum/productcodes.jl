# Copyright (c) 2022, 2023 Eric Sabo, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    HypergraphProductCode(A::CTMatrixTypes, B::CTMatrixTypes, char_vec::Union{Vector{nmod}, Missing}=missing, logs_alg::Symbol=:stnd_frm)

Return the hypergraph product code of matrices `A` and `B` whose signs are determined by `char_vec`.
"""
function HypergraphProductCode(A::CTMatrixTypes, B::CTMatrixTypes, char_vec::Union{Vector{nmod},
        Missing}=missing, logs_alg::Symbol=:stnd_frm)

    logs_alg ∈ [:stnd_frm, :VS, :sys_eqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    F = base_ring(A)
    F == base_ring(B) || throw(ArgumentError("Matrices need to be over the same base ring"))

    # note that orthogonality of C1 and C2 here is not necessary because
    # H_X * transpose(H_Z) = C1.H \otimes transpose(C2.H) + C1.H \otimes transpose(C2.H) = 0
    # in characteristic 2
    A_tr = transpose(A)
    B_tr = transpose(B)
    # branch for speedup
    if Int(order(F)) == 2
        H_X = hcat(A ⊗ identity_matrix(F, ncols(B)), identity_matrix(F, ncols(A_tr)) ⊗ B_tr)
    else
        H_X = hcat(A ⊗ identity_matrix(F, ncols(B)), -identity_matrix(F, ncols(A_tr)) ⊗ B_tr)
    end
    H_Z = hcat(identity_matrix(F, ncols(A)) ⊗ B, A_tr ⊗ identity_matrix(F, ncols(B_tr)))
    n = ncols(H_X)
    stabs = direct_sum(H_X, H_Z)
    stabs_stand, P_stand, stand_r, stand_k, rnk = _standard_form_stabilizer(stabs)
    if !iszero(stand_k)
        if logs_alg == :stnd_frm
            logs = _make_pairs(_logicals_standard_form(stabs_stand, n, stand_k, stand_r, P_stand))
            logs_mat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
        else
            rnk_H, H = right_kernel(hcat(stabs[:, n + 1:end], -stabs[:, 1:n]))
            if ncols(H) == rnk_H
                H_tr = transpose(H)
            else
                # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
                nr = nrows(H)
                H_tr = zero_matrix(base_ring(H), rnk_H, nr)
                for r in 1:nr
                    for c in 1:rnk_H
                        !iszero(H[r, c]) && (H_tr[c, r] = H[r, c];)
                    end
                end
            end
            # n + (n - Srank)
            nrows(H_tr) == 2 * n - rnk || error("Normalizer matrix is not size n + k.")
            logs, logs_mat = _logicals(stabs, H_tr, logs_alg)
        end
    end

    F = base_ring(stabs)
    p = Int(characteristic(F))
    char_vec = _process_char_vec(char_vec, p, 2 * n)
    signs, X_signs, Z_signs = _determine_signs_CSS(stabs, char_vec, nrows(H_X), nrows(H_Z))
    over_comp = nrows(stabs) > rnk

    # q^n / p^k but rows is n - k
    k = BigInt(order(F))^n // BigInt(p)^rnk
    isinteger(k) && (k = round(Int, log(BigInt(p), k));)
    # (ismissing(C1.d) || ismissing(C2.d)) ? d = missing : d = minimum([C1.d, C2.d])

    return HypergraphProductCode(F, n, k, missing, missing, missing, stabs, H_X, H_Z, missing,
        missing, signs, X_signs, Z_signs, logs, logs_mat, char_vec, over_comp, stabs_stand, stand_r,
        stand_k, P_stand, missing, missing)
end

"""
    HypergraphProductCode(C::AbstractLinearCode, char_vec::Union{Vector{nmod}, Missing}=missing)

Return the (symmetric) hypergraph product code of `C` whose signs are determined by `char_vec`.
"""
function HypergraphProductCode(C::AbstractLinearCode, char_vec::Union{Vector{nmod}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm)

    S = HypergraphProductCode(parity_check_matrix(C), parity_check_matrix(C), char_vec, logs_alg)
    S.C1 = C
    S.C2 = C
    S.d = C.d
    return S
    # H_tr = transpose(C.H)
    # eye = identity_matrix(C.F, C.n)
    # eyetr = identity_matrix(C.F, ncols(H_tr))
    # # branch for speed
    # Int(order(C.F)) == 2 ? (H_X = hcat(C.H ⊗ eye, eyetr ⊗ H_tr)) : (H_X = hcat(C.H ⊗ eye, -eyetr ⊗ H_tr);)
    # H_Z = hcat(eye ⊗ C.H, H_tr ⊗ eyetr)
    # n = ncols(H_X)
    # stabs = direct_sum(H_X, H_Z)
    # stabs_stand, P_stand, stand_r, stand_k, rnk = _standard_form_stabilizer(stabs)
    # if !iszero(stand_k)
    #     if logs_alg == :stnd_frm
    #         logs = _make_pairs(_logicals_standard_form(stabs_stand, n, stand_k, stand_r, P_stand))
    #         logs_mat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
    #     else
    #         rnk_H, H = right_kernel(hcat(stabs[:, n + 1:end], -stabs[:, 1:n]))
    #         if ncols(H) == rnk_H
    #             H_tr = transpose(H)
    #         else
    #             # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
    #             nr = nrows(H)
    #             H_tr = zero_matrix(base_ring(H), rnk_H, nr)
    #             for r in 1:nr
    #                 for c in 1:rnk_H
    #                     !iszero(H[r, c]) && (H_tr[c, r] = H[r, c];)
    #                 end
    #             end
    #         end
    #         # n + (n - Srank)
    #         nrows(H_tr) == 2 * n - rnk || error("Normalizer matrix is not size n + k.")
    #         logs, logs_mat = _logicals(stabs, H_tr, logs_alg)
    #     end
    # end

    # F = base_ring(stabs)
    # p = Int(characteristic(F))
    # char_vec = _process_char_vec(char_vec, p, 2 * n)
    # signs, X_signs, Z_signs = _determine_signs_CSS(stabs, char_vec, nrows(H_X), nrows(H_Z))
    # over_comp = nrows(stabs) > rnk

    # # q^n / p^k but rows is n - k
    # dimcode = BigInt(order(F))^n // BigInt(p)^rnk
    # isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)

    # return HypergraphProductCode(C.F, n, dimcode, C.d, missing, missing, stabs, H_X, H_Z, C, C,
    #     signs, X_signs, Z_signs, logs, logs_mat, char_vec, over_comp, stabs_stand, stand_r, stand_k, P_stand,
    #     missing, missing)
end

"""
    HypergraphProductCode(C1::AbstractLinearCode, C2::AbstractLinearCode, char_vec::Union{Vector{nmod}, Missing}=missing)

Return the hypergraph product code of `C1` and `C2` whose signs are determined by `char_vec`.
"""
function HypergraphProductCode(C1::AbstractLinearCode, C2::AbstractLinearCode,
    char_vec::Union{Vector{nmod}, Missing}=missing, logs_alg::Symbol=:stnd_frm)

    S = HypergraphProductCode(parity_check_matrix(C1), parity_check_matrix(C2), char_vec, logs_alg)
    S.C1 = C1
    S.C2 = C2
    (ismissing(C1.d) || ismissing(C2.d)) ? (S.d = missing;) : (S.d = minimum([C1.d, C2.d]);)
    return S

    # logs_alg ∈ [:stnd_frm, :VS, :sys_eqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    # Int(order(C1.F)) == Int(order(C2.F)) || throw(ArgumentError("Codes need to be over the same base ring"))

    # # note that orthogonality of C1 and C2 here is not necessary because
    # # H_X * transpose(H_Z) = C1.H \otimes transpose(C2.H) + C1.H \otimes transpose(C2.H) = 0
    # # in characteristic 2
    # H1tr = transpose(C1.H)
    # H2tr = transpose(C2.H)
    # H_X = hcat(H1 ⊗ identity_matrix(C2.F, C2.n), -identity_matrix(C1.F, ncols(H1tr)) ⊗ H2tr)
    # H_Z = hcat(identity_matrix(C1.F, C1.n) ⊗ H2, H1tr ⊗ identity_matrix(C2.F, ncols(H2tr)))
    # n = ncols(H_X)
    # stabs = direct_sum(H_X, H_Z)
    # stabs_stand, P_stand, stand_r, stand_k, rnk = _standard_form_stabilizer(stabs)
    # if !iszero(stand_k)
    #     if logs_alg == :stnd_frm
    #         logs = _make_pairs(_logicals_standard_form(stabs_stand, n, stand_k, stand_r, P_stand))
    #         logs_mat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
    #     else
    #         rnk_H, H = right_kernel(hcat(stabs[:, n + 1:end], -stabs[:, 1:n]))
    #         if ncols(H) == rnk_H
    #             H_tr = transpose(H)
    #         else
    #             # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
    #             nr = nrows(H)
    #             H_tr = zero_matrix(base_ring(H), rnk_H, nr)
    #             for r in 1:nr
    #                 for c in 1:rnk_H
    #                     !iszero(H[r, c]) && (H_tr[c, r] = H[r, c];)
    #                 end
    #             end
    #         end
    #         # n + (n - Srank)
    #         nrows(H_tr) == 2 * n - rnk || error("Normalizer matrix is not size n + k.")
    #         logs, logs_mat = _logicals(stabs, H_tr, logs_alg)
    #     end
    # end

    # F = base_ring(stabs)
    # p = Int(characteristic(F))
    # char_vec = _process_char_vec(char_vec, p, 2 * n)
    # signs, X_signs, Z_signs = _determine_signs_CSS(stabs, char_vec, nrows(H_X), nrows(H_Z))
    # over_comp = nrows(stabs) > rnk

    # # q^n / p^k but rows is n - k
    # dimcode = BigInt(order(F))^n // BigInt(p)^rnk
    # isinteger(dimcode) && (dimcode = round(Int, log(BigInt(p), dimcode));)
    # (ismissing(C1.d) || ismissing(C2.d)) ? d = missing : d = minimum([C1.d, C2.d])

    # return HypergraphProductCode(C.F, n, dimcode, d, missing, missing, stabs, H_X, H_Z, C1, C2,
    #     signs, X_signs, Z_signs, logs, logs_mat, char_vec, over_comp, stabs_stand, stand_r, stand_k, P_stand,
    #     missing, missing)
end

# unable to yield quantum LDPC code families with non constant minimum distance
"""
    GeneralizedShorCode(C1::AbstractLinearCode, C2::AbstractLinearCode, char_vec::Union{Vector{nmod}, Missing}=missing, logs_alg::Symbol=:stnd_frm)
    BaconCasaccinoConstruction(C1::AbstractLinearCode, C2::AbstractLinearCode, char_vec::Union{Vector{nmod}, Missing}=missing, logs_alg::Symbol=:stnd_frm)

Return the generalized Shor code of `C1` and `C2` with `C1⟂ ⊆ C2` whose signs are determined by `char_vec`.
"""
function GeneralizedShorCode(C1::AbstractLinearCode, C2::AbstractLinearCode,
    char_vec::Union{Vector{nmod}, Missing}=missing, logs_alg::Symbol=:stnd_frm)

    logs_alg ∈ [:stnd_frm, :VS, :sys_eqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    Int(order(C1.F)) == 2 || error("Generalized Shor codes are only defined for binary codes.")
    Int(order(C2.F)) == 2 || error("Generalized Shor codes are only defined for binary codes.")
    dual(C1) ⊆ C2 || error("Generalized Shor codes require the dual of the first code is a subset of the second.")

    # H_X stays sparse if H1 is but H_Z does not if C1.d is large
    H_X = parity_check_matrix(C1) ⊗ identity_matrix(C2.F, C2.n)
    H_Z = generator_matrix(C1) ⊗ parity_check_matrix(C2)
    S = CSSCode(H_X, H_Z, char_vec, logs_alg)
    (ismissing(C1.d) || ismissing(C2.d)) ? (S.d = missing;) : (S.d = minimum([C1.d, C2.d]);)
    return S
end
BaconCasaccinoConstruction(C1::AbstractLinearCode, C2::AbstractLinearCode,
    char_vec::Union{Vector{nmod}, Missing}=missing, logs_alg::Symbol=:stnd_frm) =
    GeneralizedShorCode(C1, C2, char_vec, logs_alg)

"""
    HyperBicycleCodeCSS(a::Vector{fq_nmod_mat}, b::Vector{fq_nmod_mat}, χ::Int, char_vec::Union{Vector{nmod}, Missing}=missing, logs_alg::Symbol=:stnd_frm)

Return the hyperbicycle CSS code of `a` and `b` given `χ` whose signs are determined by `char_vec`.

# Arguments
* a: A vector of length `c` of binary matrices of the same dimensions.
* b: A vector of length `c` of binary matrices of the same dimensions,
  potentially different from those of `a`.
* χ: A strictly positive integer coprime with `c`.
"""
function HyperBicycleCodeCSS(a::Vector{T}, b::Vector{T}, χ::Int,
    char_vec::Union{Vector{nmod}, Missing}=missing, logs_alg::Symbol=:stnd_frm) where T <: CTMatrixTypes

    logs_alg ∈ [:stnd_frm, :VS, :sys_eqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
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
    return CSSCode(GX, GZ, char_vec, logs_alg)
    # equations 41 and 42 of the paper give matrices from which the logicals may be chosen
    # not really worth it, just use standard technique from StabilizerCode.jl
end

"""
    HyperBicycleCode(a::Vector{fq_nmod_mat}, b::Vector{fq_nmod_mat}, χ::Int, char_vec::Union{Vector{nmod}, Missing}=missing, logs_alg::Symbol=:stnd_frm)

Return the hyperbicycle CSS code of `a` and `b` given `χ` whose signs are determined by `char_vec`.

# Arguments
* a: A vector of length `c` of binary matrices of the same dimensions.
* b: A vector of length `c` of binary matrices of the same dimensions,
  potentially different from those of `a`.
* χ: A strictly positive integer coprime with `c`.
"""
function HyperBicycleCode(a::Vector{T}, b::Vector{T}, χ::Int,
    char_vec::Union{Vector{nmod}, Missing}=missing, logs_alg::Symbol=:stnd_frm) where T <: CTMatrixTypes

    logs_alg ∈ (:stnd_frm, :VS, :sys_eqs) || throw(ArgumentError("Unrecognized logicals algorithm"))
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
    return StabilizerCode(stabs, char_vec, logs_alg)
    # equations 41 and 42 of the paper give matrices from which the logicals may be chosen
    # not really worth it, just use standard technique from StabilizerCode.jl
end

"""
    GeneralizedBicycleCode(A::fq_nmod_mat, B::fq_nmod_mat, char_vec::Union{Vector{nmod}, Missing}=missing, logs_alg::Symbol=:stnd_frm)

Return the generealized bicycle code given by `A` and `B` whose signs are determined by `char_vec`.
"""
function GeneralizedBicycleCode(A::T, B::T,
    char_vec::Union{Vector{nmod}, Missing}=missing, logs_alg::Symbol=:stnd_frm) where T <: CTMatrixTypes

    logs_alg ∈ [:stnd_frm, :VS, :sys_eqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    F = base_ring(A)
    F == base_ring(B) || throw(ArgumentError("Arguments must be over the same base ring."))
    (iszero(A) || iszero(B)) && throw(ArgumentError("Arguments should not be zero."))
    # this will take care of the sizes being square
    iszero(A * B - B * A) || throw(ArgumentError("Arguments must commute."))

    H_X = hcat(A, B)
    # branch for speedup
    H_Z = Int(order(F)) == 2 ? hcat(transpose(B), transpose(A)) : hcat(transpose(B), -transpose(A))
    return CSSCode(H_X, H_Z, char_vec, logs_alg)
end

"""
    GeneralizedBicycleCode(a::T, b::T, char_vec::Union{Vector{nmod}, Missing}=missing, logs_alg::Symbol=:stnd_frm) where T <: ResElem

Return the generealized bicycle code determined by `a` and `b` whose signs are determined by `char_vec`.

# Notes
* `l x l` circulant matrices are constructed using the coefficients of the polynomials
  `a` and `b` in `F_q[x]/(x^l - 1)` (`gcd(q, l) = 1`) as the first column
"""
function GeneralizedBicycleCode(a::T, b::T, char_vec::Union{Vector{nmod}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm) where T <: ResElem

    logs_alg ∈ [:stnd_frm, :VS, :sys_eqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    parent(a) == parent(b) || throw(ArgumentError("Both objects must be defined over the same residue ring."))

    return GeneralizedBicycleCode(poly_to_circ_matrix(a), poly_to_circ_matrix(b), char_vec, logs_alg)
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

Return the pre-lifted matrices `H_X` and `H_Z` of the generalized hypergraph product code of `A` and `b`.

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
    
    A_tr = transpose(A)
    for c in 1:m
        for r in 1:n
            h_coeffs = collect(coefficients(Nemo.lift(A_tr[r, c])))
            for _ in 1:l - length(h_coeffs)
                push!(h_coeffs, F(0))
            end
            h_coeffs[2:end] = reverse(h_coeffs[2:end])
            A_tr[r, c] = R(S(h_coeffs))
        end
    end
    b_coeffs = collect(coefficients(Nemo.lift(b)))
    for _ in 1:l - length(b_coeffs)
        push!(b_coeffs, F(0))
    end
    b_coeffs[2:end] = reverse(b_coeffs[2:end])
    B_tr = R(S(b_coeffs))
    Mn = MatrixSpace(R, n, n)
    H_Z = hcat(Mn(B_tr), A_tr)
    Mm = MatrixSpace(R, m, m)
    H_X = hcat(A, Mm(b))
    return H_X, H_Z
end

"""
    LiftedGeneralizedHypergraphProductCode(A::MatElem{T}, b::T, char_vec::Union{Vector{nmod}, Missing}=missing, logs_alg::Symbol=:stnd_frm) where T <: ResElem

Return the lifted generalized hypergraph product code of `A` and `b`.

# Arguments
* `A` - an `m x n` matrix with coefficents in `F_2[x]/(x^m - 1)`
* `b` - a polynomial over the same residue ring
"""
function LiftedGeneralizedHypergraphProductCode(A::MatElem{T}, b::T, char_vec::Union{Vector{nmod}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm) where T <: ResElem

    H_X, H_Z = GeneralizedHypergraphProductCodeMatrices(A, b)
    return CSSCode(lift(H_X), lift(H_Z), char_vec, logs_alg)
end

"""
    QuasiCyclicLiftedProductCodeMatrices(A::MatElem{T}, B::MatElem{T}) where T <: ResElem

Return the pre-lifted matrices `H_X` and `H_Z` for the lifted quasi-cyclic lifted product code.

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
    A_tr = transpose(A)
    for c in 1:k1
        for r in 1:n1
            h_coeffs = collect(coefficients(Nemo.lift(A_tr[r, c])))
            for _ in 1:l - length(h_coeffs)
                push!(h_coeffs, F(0))
            end
            h_coeffs[2:end] = reverse(h_coeffs[2:end])
            A_tr[r, c] = R(S(h_coeffs))
        end
    end

    k2, n2 = size(B)
    B_tr = transpose(B)
    for c in 1:k2
        for r in 1:n2
            h_coeffs = collect(coefficients(Nemo.lift(B_tr[r, c])))
            for _ in 1:l - length(h_coeffs)
                push!(h_coeffs, F(0))
            end
            h_coeffs[2:end] = reverse(h_coeffs[2:end])
            B_tr[r, c] = R(S(h_coeffs))
        end
    end

    Ek1 = identity_matrix(R, k1)
    Ek2 = identity_matrix(R, k2)
    En1 = identity_matrix(R, n1)
    En2 = identity_matrix(R, n2)

    H_X = hcat(kronecker_product(A, Ek2), kronecker_product(Ek1, B))
    H_Z = hcat(kronecker_product(En1, B_tr), kronecker_product(A_tr, En2))
    return H_X, H_Z
end

"""
    QuasiCyclicLiftedProductCode(A::MatElem{T}, B::MatElem{T}, char_vec::Union{Vector{nmod}, Missing}=missing, logs_alg::Symbol=:stnd_frm) where T <: ResElem

Return the quasi-cyclic lifted product code given by the matrices `A` and `B` and whose signs are determined by `char_vec`.
"""
function QuasiCyclicLiftedProductCode(A::MatElem{T}, B::MatElem{T}, char_vec::Union{Vector{nmod}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm) where T <: ResElem

    H_X, H_Z = QuasiCyclicLiftedProductCodeMatrices(A, B)
    return CSSCode(lift(H_X), lift(H_Z), char_vec, logs_alg)
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
    A_tr = transpose(A)
    for c in 1:k1
        for r in 1:n1
            h_coeffs = collect(coefficients(Nemo.lift(A_tr[r, c])))
            for _ in 1:l - length(h_coeffs)
                push!(h_coeffs, F(0))
            end
            h_coeffs[2:end] = reverse(h_coeffs[2:end])
            A_tr[r, c] = R(S(h_coeffs))
        end
    end

    k2, n2 = size(B)
    B_tr = transpose(B)
    for c in 1:k2
        for r in 1:n2
            h_coeffs = collect(coefficients(Nemo.lift(B_tr[r, c])))
            for _ in 1:l - length(h_coeffs)
                push!(h_coeffs, F(0))
            end
            h_coeffs[2:end] = reverse(h_coeffs[2:end])
            B_tr[r, c] = R(S(h_coeffs))
        end
    end

    Ek1 = identity_matrix(R, k1)
    Ek2 = identity_matrix(R, k2)
    En1 = identity_matrix(R, n1)
    En2 = identity_matrix(R, n2)

    A12 = kronecker_product(A_tr, Ek2)
    A13 = kronecker_product(En1, B)
    A21 = kronecker_product(A, En2)
    A24 = kronecker_product(Ek1, B_tr)
    return vcat(hcat(zeros(A21), A12, A13, zeros(A24)), hcat(A21, zeros(A12), zeros(A13), A24))
end

"""
    BiasTailoredQuasiCyclicLiftedProductCodeMatrices(A::MatElem{T}, B::MatElem{T}, char_vec::Union{Vector{nmod}, Missing}=missing, logs_alg::Symbol=:stnd_frm) where T <: ResElem

Return the bias-tailored lifted product code of `A` and `B` whose signs are given by `char_vec`.

# Arguments
* `A` - an `m x n1` matrix with coefficents in `F_2[x]/(x^m - 1)`
* `B` - an `m x n2` matrix with coefficents in the same residue ring
"""
function BiasTailoredQuasiCyclicLiftedProductCode(A::MatElem{T}, B::MatElem{T}, char_vec::Union{Vector{nmod}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm) where T <: ResElem

    stabs = BiasTailoredQuasiCyclicLiftedProductCodeMatrices(A, B)
    return StabilizerCode(lift(stabs), char_vec, logs_alg)
end

"""
    SPCDFoldProductCode(D::Int, s::Int = 1)
    SingleParityCheckDFoldProductCode(D::Int, s::Int = 1) = SPCDFoldProductCode(D, s)

Return the single-parity-check `D`-fold product code.

* Note
- This is defined in https://arxiv.org/abs/2209.13474
"""
function SPCDFoldProductCode(D::Int, s::Int = 1)
    vec_S = Vector{AbstractStabilizerCode}()
    for i in 1:D
        for l in 1:D
            if l == (i - 1) * D + i
                # RepetitionCode(2, x) is dual(SPCCode(2, x))
                push!(vec_S, CSSCode(RepetitionCode(2, 2 * s)))
            else
                push!(vec_S, CSSCode(RepetitionCode(2, 2)))
            end
        end
    end
    
    S = symmetric_product(vec_S)
    set_minimum_distance!(S, 2^D)
    S.d_x = S.d
    S.d_x = S.d
    S.is_pure = true
    # metacheck distance = 3
    return S
end
SingleParityCheckDFoldProductCode(D::Int, s::Int = 1) = SPCDFoldProductCode(D, s)

#############################
      # getter functions
#############################

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

"""
    Quintavalle_basis(C::HypergraphProductCode)

Return a symplectic canonical basis for the logical operators of `C`.
    
* Note
- This implements https://doi.org/10.48550/arXiv.2204.10812.
"""
function Quintavalle_basis(C::HypergraphProductCode)
    H1 = parity_check_matrix(C.C1)
    H2 = parity_check_matrix(C.C2)

    # c - complement
    ker_H1, im_H1_tr_c = strongly_lower_triangular_reduction(H1)
    ker_H1_tr, im_H1_c = strongly_lower_triangular_reduction(transpose(H1))
    ker_H2, im_H2_tr_c = strongly_lower_triangular_reduction(H2)
    ker_H2_tr, im_H2_c = strongly_lower_triangular_reduction(transpose(H2))
    F = C.F
    lx = zero_matrix(F, C.k, C.n)
    lz = deepcopy(lx)

    l = 1
    temp = zero_matrix(F, 1, nrows(ker_H1_tr) * nrows(ker_H2_tr))
    tr_im_H1_tr_c = transpose(im_H1_tr_c)
    tr_ker_H2 = transpose(ker_H2)
    tr_ker_H1 = transpose(ker_H1)
    tr_im_H2_tr_c = transpose(im_H2_tr_c)
    for i in 1:nrows(tr_ker_H1)
        for h in 1:nrows(tr_ker_H2)
            lx[l, :] = hcat(tr_im_H1_tr_c[i, :] ⊗ tr_ker_H2[h, :], temp)
            lz[l, :] = hcat(tr_ker_H1[i, :] ⊗ tr_im_H2_tr_c[h, :], temp)
            l += 1
        end
    end

    temp = zero_matrix(F, 1, nrows(ker_H1) * nrows(ker_H2))
    tr_ker_H1_tr = transpose(ker_H1_tr)
    tr_im_H2_c = transpose(im_H2_c)
    tr_im_H1_c = transpose(im_H1_c)
    tr_ker_H2_tr = transpose(ker_H2_tr)
    for i in 1:nrows(tr_ker_H1_tr)
        for h in 1:nrows(tr_ker_H2_tr)
            lx[l, :] = hcat(temp, tr_ker_H1_tr[i, :] ⊗ tr_im_H2_c[h, :])
            lz[l, :] = hcat(temp, tr_im_H1_c[i, :] ⊗ tr_ker_H2_tr[h, :])
            l += 1
        end
    end
    return lx, lz
end

# TODO: present the stabilizers in the docs and mention how to switch X and Z by
# using the switch on the inputs beforehand
"""
    asymmetric_product(S1::T, S2::T) where {T <: AbstractSubsystemCode}

Return the asymmetric 2-fold product quantum CSS code of the CSS codes `S1` and `S2`.

* Note
- This is defined in https://arxiv.org/abs/2209.13474
"""
asymmetric_product(S1::T, S2::T) where {T <: AbstractSubsystemCode} = asymmetric_product(CSSTrait(T), S1, S2)
function asymmetric_product(::IsCSS, S1::AbstractSubsystemCode, S2::AbstractSubsystemCode)
    # TODO: check fields are the same or convertable
    H_X = vcat(S1.X_stabs ⊗ identity_matrix(S1.F, S2.n), identity_matrix(S1.F, S1.n) ⊗ S2.X_stabs)
    H_Z = S1.Z_stabs ⊗ S2.Z_stabs
    return CSSCode(H_X, H_Z)
end
asymmetric_product(::IsNotCSS, S1::AbstractSubsystemCode, S2::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    symmetric_product(vec_S::Vector{T}) where {T <: AbstractSubsystemCode}

Return the symmetric `D`-fold product quantum CSS code, where `D` is
the square-root of the length of the vector of CSS codes `vec_S`.

* Note
- This is defined in https://arxiv.org/abs/2209.13474
"""
function symmetric_product(vec_S::Vector{T}) where {T <: AbstractSubsystemCode}
    isempty(vec_S) && throw(ArgumentError("Input vector of CSS codes cannot be empty"))
    for S in vec_S
        if CSSTrait(typeof(S)) == IsNotCSS()
            return symmetric_product(IsNotCSS(), vec_S)
        end
    end
    return symmetric_product(IsCSS(), vec_S)
end
function symmetric_product(::IsCSS, vec_S::Vector{T}) where {T <: AbstractSubsystemCode}
    # TODO: check fields are the same or convertable
    length(vec_S) >= 4 || throw(DomainError("The length of the input vector must be at least 4"))
    D = sqrt(length(vec_S))
    isinteger(D) ? (D = Int(D);) : throw(ArgumentError("The number of CSS codes must be D^2"))

    # since the sizes of the identities are different for every entry,
    # there's not too much of a better way to do this
    F = vec_S[1].F
    # X stabilizers
    H_X = nothing
    for j in 0:D - 1
        temp_row = nothing
        for l in 1:D^2
            if j * D + 1 <= l <= (j + 1) * D
                if l == 1
                    temp_row = vec_S[l].X_stabs
                else
                    temp_row = temp_row ⊗ vec_S[l].X_stabs
                end
            else
                if l == 1
                    temp_row = identity_matrix(F, vec_S[l].n)
                else
                    temp_row = temp_row ⊗ identity_matrix(F, vec_S[l].n)
                end
            end
        end
        if j == 0
            H_X = temp_row
        else
            H_X = vcat(H_X, temp_row)
        end
    end

    # Z stabilizers
    H_Z = nothing
    for j in 0:D - 1
        temp_row = nothing
        for l in 1:D^2
            if (l - 1) % D == j
                if l == 1
                    temp_row = vec_S[l].Z_stabs
                else
                    temp_row = temp_row ⊗ vec_S[l].Z_stabs
                end
            else
                if l == 1
                    temp_row = identity_matrix(F, vec_S[l].n)
                else
                    temp_row = temp_row ⊗ identity_matrix(F, vec_S[l].n)
                end
            end
        end
        if j == 0
            H_Z = temp_row
        else
            H_Z = vcat(H_Z, temp_row)
        end
    end
    
    return CSSCode(H_X, H_Z)
end
symmetric_product(::IsNotCSS, vec_S::Vector{T}) where {T <: AbstractSubsystemCode} = error("Only valid for CSS codes.")

# extend to subsystem codes?
⊠(S::AbstractStabilizerCode, U::Union{Missing, T} = missing) where T <: CTMatrixTypes = ⊠(CSSTrait(typeof(S)), S, U)
function ⊠(::IsCSS, S::AbstractStabilizerCode, U::Union{Missing, T} = missing) where T <: CTMatrixTypes
    num_stabs = num_X_stabs(S)
    num_stabs == num_Z_stabs(S) || throw(ArgumentError("The code must have the same number of X and Z stabilizers"))
    if !ismissing(U)
        T == typeof(S.X_stabs) || throw(ArgumentError("Input matrix must have the same type as the stabilizers"))
        base_ring(U) == S.F || throw(ArgumentError("Input matrix must have the same base ring as the stabilizers"))
    end
    
    δ = zero_matrix(S.F, S.n, S.n)
    for i in 1:num_stabs
        for j in 1:num_stabs
            if ismissing(U)
                δ += transpose(S.Z_stabs[i, :]) * S.X_stabs[j, :]
            else
                δ += U[i, j] * transpose(S.Z_stabs[i, :]) * S.X_stabs[j, :]
            end
        end
    end
    chain = ChainComplex([δ])
    prod = chain ⊗ chain
    # extract code from this
    return prod
    # kernel of boundary is
    # return CSSCode(..., ...)
end
⊠(::IsNotCSS, S::AbstractStabilizerCode, U::Union{Missing, T} = missing) where T <: CTMatrixTypes = throw(ArgumentError("This is only defined for CSS codes"))
homological_product(S::AbstractStabilizerCode, U::Union{Missing, T} = missing) where T <: CTMatrixTypes = ⊠(S, U)
