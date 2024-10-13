# Copyright (c) 2022, 2023 Eric Sabo, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    HypergraphProductCode(A::CTMatrixTypes, B::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing}= missing, logs_alg::Symbol = :stnd_frm)

Return the hypergraph product code of matrices `A` and `B`.

# Example

[[1922, 50, 16]] Hypergraph Product Code from Appendix B, Example C1 of [panteleev2021degenerate](@cite).

```jldoctest
julia> using CodingTheory, Oscar;

julia> F, x = polynomial_ring(Oscar.Nemo.Native.GF(2), :x);

julia> l = 31;

julia> R, = residue_ring(F, x^l -1);

julia> h = R(1 + x^2 + x^5);

julia> A = residue_polynomial_to_circulant_matrix(h);

julia> code = HypergraphProductCode(A, A);

julia> length(code), dimension(code)
(1922, 50)
```
"""
function HypergraphProductCode(A::CTMatrixTypes, B::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem},
    Missing} = missing, logs_alg::Symbol = :stnd_frm)

    logs_alg ∈ (:stnd_frm, :VS, :sys_eqs) || throw(ArgumentError("Unrecognized logicals algorithm"))
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
    # TODO is this distance formula not correct?
    # (ismissing(C1.d) || ismissing(C2.d)) ? d = missing : d = minimum([C1.d, C2.d])
    X_logs = reduce(vcat, [log[1][:, 1:n] for log in logs])
    Z_logs = reduce(vcat, [log[2][:, n + 1:end] for log in logs])
    _, mat = rref(vcat(H_X, X_logs))
    anti = _remove_empty(mat, :rows) * transpose(Z_logs)
    u_bound_dx, _ = _min_wt_row(mat[findall(!iszero(anti[i:i, :]) for i in axes(anti, 1)), :])
    _, mat = rref(vcat(H_Z, Z_logs))
    anti = _remove_empty(mat, :rows) * transpose(X_logs)
    u_bound_dz, _ = _min_wt_row(mat[findall(!iszero(anti[i:i, :]) for i in axes(anti, 1)), :])
    return HypergraphProductCode(F, n, k, missing, missing, missing, 1, min(u_bound_dx,
        u_bound_dz), 1, u_bound_dx, 1, u_bound_dz, stabs, H_X, H_Z, missing, missing, signs,
        X_signs, Z_signs, logs, logs_mat, char_vec, over_comp, stabs_stand, stand_r, stand_k,
        P_stand, missing, missing, missing, missing)
end

"""
    HypergraphProductCode(C::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

Return the (symmetric) hypergraph product code of `C`.
"""
function HypergraphProductCode(C::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} =
    missing, logs_alg::Symbol = :stnd_frm)

    S = HypergraphProductCode(parity_check_matrix(C), parity_check_matrix(C), char_vec = char_vec, logs_alg = logs_alg)
    S.C1 = C
    S.C2 = C
    ismissing(C.d) || set_minimum_distance!(S, C.d)
    return S
end

"""
    HypergraphProductCode(C1::AbstractLinearCode, C2::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

Return the hypergraph product code of `C1` and `C2`.
"""
function HypergraphProductCode(C1::AbstractLinearCode, C2::AbstractLinearCode;
    char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

    S = HypergraphProductCode(parity_check_matrix(C1), parity_check_matrix(C2), char_vec =
        char_vec, logs_alg = logs_alg)
    S.C1 = C1
    S.C2 = C2
    (ismissing(C1.d) || ismissing(C2.d)) ? (S.d = missing;) :
        (set_minimum_distance!(S, min(C1.d, C2.d));)
    return S
end

# unable to yield quantum LDPC code families with non constant minimum distance
"""
    GeneralizedShorCode(C1::AbstractLinearCode, C2::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)
    BaconCasaccinoConstruction(C1::AbstractLinearCode, C2::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

Return the generalized Shor code of `C1` and `C2` with `C1⟂ ⊆ C2`.
"""
function GeneralizedShorCode(C1::AbstractLinearCode, C2::AbstractLinearCode;
    char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

    logs_alg ∈ (:stnd_frm, :VS, :sys_eqs) || throw(ArgumentError("Unrecognized logicals algorithm"))
    Int(order(C1.F)) == 2 || error("Generalized Shor codes are only defined for binary codes.")
    Int(order(C2.F)) == 2 || error("Generalized Shor codes are only defined for binary codes.")
    dual(C1) ⊆ C2 || error("Generalized Shor codes require the dual of the first code is a subset of the second.")

    # H_X stays sparse if H1 is but H_Z does not if C1.d is large
    H_X = parity_check_matrix(C1) ⊗ identity_matrix(C2.F, C2.n)
    H_Z = generator_matrix(C1) ⊗ parity_check_matrix(C2)
    S = CSSCode(H_X, H_Z, char_vec = char_vec, logs_alg = logs_alg)
    (ismissing(C1.d) || ismissing(C2.d)) ? (S.d = missing;) :
        (set_minimum_distance!(S, min(C1.d, C2.d));)
    return S
end
BaconCasaccinoConstruction(C1::AbstractLinearCode, C2::AbstractLinearCode,
    char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) =
    GeneralizedShorCode(C1, C2, char_vec = char_vec, logs_alg = logs_alg)

"""
    HyperBicycleCodeCSS(a::Vector{CTMatrixTypes}, b::Vector{CTMatrixTypes}, χ::Int; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

Return the hyperbicycle CSS code of `a` and `b` given `χ`.

# Arguments
- a: A vector of length `c` of binary matrices of the same dimensions.
- b: A vector of length `c` of binary matrices of the same dimensions,
  potentially different from those of `a`.
- χ: A strictly positive integer coprime with `c`.

# Example

[[900, 50, 14]] CSS Hyperbicycle Code from Example 6 of [Kovalev_2013](@cite).

```jldoctest
julia> S, x = polynomial_ring(Oscar.Nemo.Native.GF(2), :x);

julia> l = 30; χ = 1;

julia> R, = residue_ring(S, x^l - 1);

julia> h = R(1 + x + x^3 + x^5);

julia> A = residue_polynomial_to_circulant_matrix(h);

julia> a1 = A[1:15, 1:15];

julia> a2 = A[1:15, 16:30];

julia> code = HyperBicycleCodeCSS([a1, a2], [a1, a2], χ);

julia> length(code), dimension(code)
(900, 50)
```
"""
function HyperBicycleCodeCSS(a::Vector{T}, b::Vector{T}, χ::Int; char_vec::Union{Vector{zzModRingElem},
    Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: CTMatrixTypes

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

    H1 = zero_matrix(F, c * k1, c * n1)
    H2 = zero_matrix(F, c * k2, c * n2)
    HT1 = zero_matrix(F, c * n1, c * k1)
    HT2 = zero_matrix(F, c * n2, c * k2)
    Ic = identity_matrix(F, c)
    Sχ = Ic[mod1.(1:χ:c * χ, c), :]
    for i in 1:c
        Ii = vcat(Ic[i:end, :], Ic[1:i - 1, :])
        Iχi = Sχ * Ii
        ITχi = transpose(Sχ) * transpose(Ii)
        H1 += Iχi ⊗ a[i]
        H2 += b[i] ⊗ Iχi
        HT1 += ITχi ⊗ transpose(a[i])
        HT2 += transpose(b[i]) ⊗ ITχi
    end

    Ek1 = identity_matrix(F, k1)
    Ek2 = identity_matrix(F, k2)
    En1 = identity_matrix(F, n1)
    En2 = identity_matrix(F, n2)

    GX = hcat(Ek2 ⊗ H1, H2 ⊗ Ek1)
    GZ = hcat(HT2 ⊗ En1, En2 ⊗ HT1)
    return CSSCode(GX, GZ, char_vec = char_vec, logs_alg = logs_alg)
end

"""
    HyperBicycleCode(a::Vector{CTMatrixTypes}, b::Vector{CTMatrixTypes}, χ::Int; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

Return the hyperbicycle non-CSS code of `a` and `b` given `χ`.

# Arguments
- a: A vector of length `c` of binary matrices of the same dimensions.
- b: A vector of length `c` of binary matrices of the same dimensions,
  potentially different from those of `a`.
- χ: A strictly positive integer coprime with `c`.

# Example

[[289, 81, 5]] non-CSS Hyperbicycle Code from Example 13 of [Kovalev_2013](@cite).

```jldoctest
julia> using CodingTheory, Oscar;

julia> S, x = polynomial_ring(Oscar.Nemo.Native.GF(2), :x);

julia> l = 17; χ = 1;

julia> R, = residue_ring(S, x^l - 1);

julia> h = R(x^4 * (1 + x + x^3 + x^6 + x^8 + x^9));

julia> A = residue_polynomial_to_circulant_matrix(h);

julia> code = HyperBicycleCode([A], [A], χ);

julia> length(code), dimension(code)
(289, 81)
```
"""
function HyperBicycleCode(a::Vector{T}, b::Vector{T}, χ::Int, char_vec::Union{Vector{zzModRingElem},
    Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: CTMatrixTypes

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

    H1 = zero_matrix(F, c * k1, c * n1)
    H2 = zero_matrix(F, c * k2, c * n2)
    HT1 = zero_matrix(F, c * n1, c * k1)
    HT2 = zero_matrix(F, c * n2, c * k2)
    Ic = identity_matrix(F, c)
    Sχ = Ic[mod1.(1:χ:c * χ, c), :]
    for i in 1:c
        Ii = vcat(Ic[i:end, :], Ic[1:i - 1, :])
        Iχi = Sχ * Ii
        ITχi = transpose(Sχ) * transpose(Ii)
        H1 += Iχi ⊗ a[i]
        H2 += b[i] ⊗ Iχi
        HT1 += ITχi ⊗ transpose(a[i])
        HT2 += transpose(b[i]) ⊗ ITχi
    end

    (H1 == HT1 && H2 == HT2) || throw(ArgumentError("H_i must equal H̃_i for i = 1, 2."))

    Ek1 = identity_matrix(F, k1)
    Ek2 = identity_matrix(F, k2)

    stabs = hcat(Ek2 ⊗ H1, H2 ⊗ Ek1)
    return StabilizerCode(stabs, char_vec = char_vec, logs_alg = logs_alg)
end

"""
    GeneralizedBicycleCode(A::CTMatrixTypes, B::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

Return the generealized bicycle code given by `A` and `B`.

# Example

[[254, 28, 14 ≤ d ≤ 20]] Generalized Bicycle Code from Appendix B, Example A1 of [panteleev2021degenerate](@cite).

```jldoctest
julia> using CodingTheory, Oscar;

julia> F = Oscar.Nemo.Native.GF(2);

julia> S, x = polynomial_ring(F, :x);

julia> l = 127;

julia> R, _ = residue_ring(S, x^l - 1);

julia> a = 1 + x^15 + x^20 + x^28 + x^66;

julia> b = 1 + x^58 + x^59 + x^100 + x^121;

julia> code = GeneralizedBicycleCode(R(a), R(b));

julia> length(code), dimension(code)
(254, 28)
```
"""
function GeneralizedBicycleCode(A::T, B::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm) where T <: CTMatrixTypes

    logs_alg ∈ (:stnd_frm, :VS, :sys_eqs) || throw(ArgumentError("Unrecognized logicals algorithm"))
    F = base_ring(A)
    F == base_ring(B) || throw(ArgumentError("Arguments must be over the same base ring."))
    (iszero(A) || iszero(B)) && throw(ArgumentError("Arguments should not be zero."))
    # this will take care of the sizes being square
    iszero(A * B - B * A) || throw(ArgumentError("Arguments must commute."))

    H_X = hcat(A, B)
    # branch for speedup
    H_Z = Int(order(F)) == 2 ? hcat(transpose(B), transpose(A)) : hcat(transpose(B), -transpose(A))
    return CSSCode(H_X, H_Z, char_vec = char_vec, logs_alg = logs_alg)
end

"""
    GeneralizedBicycleCode(a::T, b::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: ResElem

Return the generealized bicycle code determined by `a` and `b`.

# Notes
- `l x l` circulant matrices are constructed using the coefficients of the polynomials
  `a` and `b` in `F_q[x]/(x^l - 1)` (`gcd(q, l) = 1`) as the first column
"""
function GeneralizedBicycleCode(a::T, b::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm) where T <: ResElem

    logs_alg ∈ (:stnd_frm, :VS, :sys_eqs) || throw(ArgumentError("Unrecognized logicals algorithm"))
    parent(a) == parent(b) || throw(ArgumentError("Both objects must be defined over the same residue ring."))

    return GeneralizedBicycleCode(residue_polynomial_to_circulant_matrix(a),
        residue_polynomial_to_circulant_matrix(b), char_vec = char_vec, logs_alg = logs_alg)
end

"""
    GeneralizedBicycleCode(a::T, b::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: CTGroupAlgebra

Return the generealized bicycle code determined by `a` and `b`.

# Notes
- `|G| x |G|` circulant matrices are constructed using the coefficients of the elements in the group algebra `FG` as` the first column
"""
function GeneralizedBicycleCode(a::T, b::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm) where T <: CTGroupAlgebra

    logs_alg ∈ (:stnd_frm, :VS, :sys_eqs) || throw(ArgumentError("Unrecognized logicals algorithm"))
    parent(a) == parent(b) || throw(ArgumentError("Both objects must be defined over the same residue ring."))

    return GeneralizedBicycleCode(group_algebra_element_to_circulant_matrix(a),
        group_algebra_element_to_circulant_matrix(b), char_vec = char_vec, logs_alg = logs_alg)
end

# function BicycleCode(A::fq_nmot_mat)
#     m, n = size(A)
#     m == n || throw(ArgumentError("Input matrix must be square."))
#     # should probably check for F_2

#     H = hcat(A, transpose(A))
#     return CSSCode(H, H)
# end
    
"""
    generalized_hypergraph_product_matrices(A::MatElem{T}, b::T) where T <: ResElem
    GHGP_matrices(A::MatElem{T}, b::T) where T <: ResElem
    lifted_product_matrices(A::MatElem{T}, b::T) where T <: ResElem

Return the pre-lifted matrices `H_X` and `H_Z` of the generalized hypergraph product code of `A` and `b`.

# Arguments
- `A` - an `m x n` matrix with elements in `F_2[x]/(x^m - 1)`
- `b` - a polynomial over the same residue ring

# Notes
- Use `LiftedProductCode` to return a quantum code over the base ring directly.
"""
function generalized_hypergraph_product_matrices(A::MatElem{T}, b::T) where T <: ResElem

    # @warn "Commutativity of A and b required but not yet enforced."
    S = base_ring(b)
    F = base_ring(S)
    # Int(order(F)) == 2 || throw(ArgumentError("The generalized hypergraph product is only defined over GF(2)."))
    R = parent(A[1, 1])
    R == parent(b) || throw(ArgumentError("Both objects must be defined over the same residue ring."))
    m, n = size(A)
    (m != 1 && n != 1) || throw(ArgumentError("First input matrix must not be a vector."))
    f = modulus(R)
    l = degree(f)
    f == gen(S)^l - 1 || throw(ArgumentError("Residue ring not of the form x^l - 1."))
    # gcd(l, Int(characteristic(F))) == 1 || throw(ArgumentError("Residue ring over F_q[x] must be defined by x^l - 1 with gcd(l, q) = 1."))
    
    A_tr = _CT_adjoint(A)
    # b_coeffs = collect(coefficients(Nemo.lift(b)))
    # for _ in 1:l - length(b_coeffs)
    #     push!(b_coeffs, F(0))
    # end
    # b_coeffs[2:end] = reverse(b_coeffs[2:end])
    # B_tr = R(S(b_coeffs))
    B_tr = _CT_adjoint(matrix(R, 1, 1, [b]))[1, 1]
    Mn = matrix_space(R, n, n)
    H_Z = hcat(Mn(B_tr), A_tr)
    Mm = matrix_space(R, m, m)
    # TODO: check extending this function past F_2 makes sense
    # branch for speed
    if Int(order(F)) == 2
        H_X = hcat(A, Mm(b))
    else
        H_X = hcat(A, -Mm(b))
    end
    return H_X, H_Z
end
GHGP_matrices(A::MatElem{T}, b::T) where T <: ResElem =
    generalized_hypergraph_product_matrices(A, b)
lifted_product_matrices(A::MatElem{T}, b::T) where T <: ResElem =
    generalized_hypergraph_product_matrices(A, b)

"""
    generalized_hypergraph_product_matrices(A::MatElem{T}, b::T) where T <: CTGroupAlgebra
    GHGP_matrices(A::MatElem{T}, b::T) where T <: CTGroupAlgebra
    lifted_product_matrices(A::MatElem{T}, b::T) where T <: CTGroupAlgebra

Return the pre-lifted matrices `H_X` and `H_Z` of the generalized hypergraph product code of `A` and `b`.

# Arguments
- `A` - a matrix with elements in a group algebra
- `b` - a polynomial over the same group algebra

# Notes
- Use `LiftedProductCode` to return a quantum code over the base ring directly.
"""
function generalized_hypergraph_product_matrices(A::MatElem{T}, b::T) where T <: CTGroupAlgebra
    FG = parent(A[1, 1])
    parent(b) == FG || throw(ArgumentError("Inputs must be over the same group algebra"))

    A_tr = _CT_adjoint(A)
    B_tr = _CT_adjoint(matrix(FG, 1, 1, [b]))[1, 1]
    m, n = size(A)
    Mn = matrix_space(FG, n, n)
    H_Z = hcat(Mn(B_tr), A_tr)
    Mm = matrix_space(FG, m, m)
    F = base_ring(FG)
    # TODO: check extending this function past F_2 makes sense
    # branch for speed
    if Int(order(F)) == 2
        H_X = hcat(A, Mm(b))
    else
        H_X = hcat(A, -Mm(b))
    end
    return H_X, H_Z
end
GHGP_matrices(A::MatElem{T}, b::T) where T <: CTGroupAlgebra =
    generalized_hypergraph_product_matrices(A, b)
lifted_product_matrices(A::MatElem{T}, b::T) where T <: CTGroupAlgebra =
    generalized_hypergraph_product_matrices(A, b)

"""
    GeneralizedHypergraphProductCode(A::MatElem{T}, b::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}
    LiftedProductCode(A::MatElem{T}, b::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
        logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

Return the lifted (generalized hypergraph) product code of `A` and `b`.

# Arguments
- `A` - either an `m x n` matrix with elements in `F_2[x]/(x^m - 1)` or a group algebra
- `b` - a polynomial over the same residue ring or group algebra
"""
function GeneralizedHypergraphProductCode(A::MatElem{T}, b::T; char_vec::Union{Vector{zzModRingElem},
    Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

    H_X, H_Z = generalized_hypergraph_product_matrices(A, b)
    return CSSCode(lift(H_X), lift(H_Z), char_vec = char_vec, logs_alg = logs_alg)
end
LiftedProductCode(A::MatElem{T}, b::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm) where T <: ResElem =
    GeneralizedHypergraphProductCode(A, b, char_vec = char_vec, logs_alg = logs_alg)

"""
    lifted_product_matrices(A::MatElem{T}, B::MatElem{T}) where T <: ResElem

Return the pre-lifted matrices `H_X` and `H_Z` for the lifted quasi-cyclic lifted product code.

# Arguments
- `A` - an `m x n1` matrix with elements in `F_2[x]/(x^m - 1)`
- `B` - an `m x n2` matrix with elements in the same residue ring

# Notes
- Use `LiftedProductCode` to return a quantum code over the base ring directly.
"""
function lifted_product_matrices(A::MatElem{T}, B::MatElem{T}) where T <: ResElem
    # @warn "Commutativity of A and b required but not yet enforced."
    S = base_ring(A[1, 1])
    F = base_ring(S)
    Int(order(F)) == 2 || throw(ArgumentError("The quasi-cyclic lifted product is only defined over GF(2)."))
    R = parent(A[1, 1])
    R == parent(B[1, 1]) || throw(ArgumentError("Both objects must be defined over the same residue ring."))
    f = modulus(R)
    l = degree(f)
    f == gen(S)^l - 1 || throw(ArgumentError("Residue ring not of the form x^l - 1."))
    
    A_tr = _CT_adjoint(A)
    B_tr = _CT_adjoint(B)

    k1, n1 = size(A)
    k2, n2 = size(B)
    Ek1 = identity_matrix(R, k1)
    Ek2 = identity_matrix(R, k2)
    En1 = identity_matrix(R, n1)
    En2 = identity_matrix(R, n2)

    H_X = hcat(A ⊗ Ek2, Ek1 ⊗ B)
    H_Z = hcat(En1 ⊗ B_tr, A_tr ⊗ En2)
    return H_X, H_Z
end

"""
    lifted_product_matrices(A::MatElem{T}, B::MatElem{T}) where T <: ResElem

Return the pre-lifted matrices `H_X` and `H_Z` for the lifted quasi-cyclic lifted product code.

# Arguments
- `A` - a matrix with elements in a group algebra
- `B` - a matrix with coefficents in the same group algebra

# Notes
- Use `LiftedProductCode` to return a quantum code over the base ring directly.
"""
function lifted_product_matrices(A::MatElem{T}, B::MatElem{T}) where T <: CTGroupAlgebra
    FG = parent(A[1, 1])
    parent(B[1, 1]) == FG || throw(ArgumentError("Inputs must be over the same group algebra"))

    A_tr = _CT_adjoint(A)
    B_tr = _CT_adjoint(B)

    k1, n1 = size(A)
    k2, n2 = size(B)
    Ek1 = identity_matrix(FG, k1)
    Ek2 = identity_matrix(FG, k2)
    En1 = identity_matrix(FG, n1)
    En2 = identity_matrix(FG, n2)

    H_X = hcat(A ⊗ Ek2, Ek1 ⊗ B)
    H_Z = hcat(En1 ⊗ B_tr, A_tr ⊗ En2)
    return H_X, H_Z
end

"""
    LiftedProductCode(A::MatElem{T}, B::MatElem{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

Return the lifted product code given by the matrices `A` and `B`.

# Example

[[882, 24, 18 ≤ d ≤ 24]] Lifted Product Code from Appendix B, Example B1 of [panteleev2021degenerate](@cite).

```jldoctest
julia> using CodingTheory, Oscar;

julia> F = Oscar.Nemo.Native.GF(2);

julia> S, x = polynomial_ring(F, :x);

julia> l = 63;

julia> R, _ = residue_ring(S, x^l - 1);

julia> A = matrix(R, 7, 7,
           [x^27, 0   , 0   , 0   , 0   , 1   , x^54,
            x^54, x^27, 0   , 0   , 0   , 0   , 1   ,
            1   , x^54, x^27, 0   , 0   , 0   , 0   ,
            0   , 1   , x^54, x^27, 0   , 0   , 0   ,
            0   , 0   , 1   , x^54, x^27, 0   , 0   ,
            0   , 0   , 0   , 1   , x^54, x^27, 0   ,
            0   , 0   , 0   , 0   , 1   , x^54, x^27]);

julia> b = R(1 + x + x^6);

julia> code = LiftedProductCode(A, b);

julia> length(code), dimension(code)
(882, 24)
```
"""
function LiftedProductCode(A::MatElem{T}, B::MatElem{T}; char_vec::Union{Vector{zzModRingElem}, Missing} =
    missing, logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

    H_X, H_Z = lifted_product_matrices(A, B)
    return CSSCode(lift(H_X), lift(H_Z), char_vec = char_vec, logs_alg = logs_alg)
end

"""
    bias_tailored_lifted_product_matrices(A::MatElem{T}, B::MatElem{T}) where T <: ResElem

Return the pre-lifted stabilizer matrix for bias-tailored lifted product code of `A` and `B`.

# Arguments
- `A` - an `m x n1` matrix with elements in `F_2[x]/(x^m - 1)`
- `B` - an `m x n2` matrix with elements in the same residue ring

# Notes
- Use `BiasTailoredLiftedProductCode` to return a quantum code over the base ring directly.
"""
function bias_tailored_lifted_product_matrices(A::MatElem{T}, B::MatElem{T}) where T <: ResElem
    # @warn "Commutativity of A and b required but not yet enforced."
    S = base_ring(A[1, 1])
    F = base_ring(S)
    Int(order(F)) == 2 || throw(ArgumentError("The quasi-cyclic lifted product is only defined over GF(2)."))
    R = parent(A[1, 1])
    R == parent(B[1, 1]) || throw(ArgumentError("Both objects must be defined over the same residue ring."))
    f = modulus(R)
    l = degree(f)
    f == gen(S)^l - 1 || throw(ArgumentError("Residue ring not of the form x^l - 1."))
    
    A_tr = _CT_adjoint(A)
    B_tr = _CT_adjoint(B)

    k1, n1 = size(A)
    k2, n2 = size(B)
    Ek1 = identity_matrix(R, k1)
    Ek2 = identity_matrix(R, k2)
    En1 = identity_matrix(R, n1)
    En2 = identity_matrix(R, n2)

    A12 = A_tr ⊗ Ek2
    A13 = En1 ⊗ B
    A21 = A ⊗ En2
    A24 = Ek1 ⊗ B_tr
    return vcat(hcat(zero(A21), A12, A13, zero(A24)), hcat(A21, zero(A12), zero(A13), A24))
end

"""
    bias_tailored_lifted_product_matrices(A::MatElem{T}, B::MatElem{T}) where T <: CTGroupAlgebra

Return the pre-lifted stabilizer matrix for bias-tailored lifted product code of `A` and `B`.

# Arguments
- `A` - a matrix with elements in a group algebra
- `B` - a matrix with elements in the same group algebra

# Notes
- Use `BiasTailoredLiftedProductCode` to return a quantum code over the base ring directly.

# Example

[[882, 24, d ≤ 24]] BiasTailored Lifted Product Code from Appendix B of [roffe2023bias](@cite).

```jldoctest
julia> using CodingTheory, Oscar;

julia> F = Oscar.Nemo.Native.GF(2);

julia> S, x = polynomial_ring(F, :x);

julia> l = 63;

julia> R, _ = residue_ring(S, x^l - 1);

julia> A1 = matrix(R, 1, 1, [1 + x^1 + x^6]);

julia> A2 = matrix(R, 7, 7,
           [x^36, 0   , 0   , 0   , 0   , 1   , x^9 ,
            x^9 , x^36, 0   , 0   , 0   , 0   , 1   ,
            1   , x^9 , x^36, 0   , 0   , 0   , 0   ,
            0   , 1   , x^9 , x^36, 0   , 0   , 0   ,
            0   , 0   , 1   , x^9 , x^36, 0   , 0   ,
            0   , 0   , 0   , 1   , x^9 , x^36, 0   ,
            0   , 0   , 0   , 0   , 1   , x^9 , x^36]);

julia> code = BiasTailoredLiftedProductCode(A1, A2);

julia> length(code), dimension(code)
(882, 24)
```
"""
function bias_tailored_lifted_product_matrices(A::MatElem{T}, B::MatElem{T}) where T <: CTGroupAlgebra

    FG = parent(A[1, 1])
    parent(B[1, 1]) == FG || throw(ArgumentError("Inputs must be over the same group algebra"))

    A_tr = _CT_adjoint(A)
    B_tr = _CT_adjoint(B)

    k1, n1 = size(A)
    k2, n2 = size(B)
    Ek1 = identity_matrix(R, k1)
    Ek2 = identity_matrix(R, k2)
    En1 = identity_matrix(R, n1)
    En2 = identity_matrix(R, n2)

    A12 = A_tr ⊗ Ek2
    A13 = En1 ⊗ B
    A21 = A ⊗ En2
    A24 = Ek1 ⊗ B_tr
    return vcat(hcat(zero(A21), A12, A13, zero(A24)), hcat(A21, zero(A12), zero(A13), A24))
end

"""
    BiasTailoredLiftedProductCode(A::MatElem{T}, B::MatElem{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

Return the bias-tailored lifted product code of `A` and `B`.

# Arguments
- `A` - either an `m x n` matrix with elements in `F_2[x]/(x^m - 1)` or a group algebra
- `B` - an `m x n2` matrix with elements in the same parent as `A`
"""
function BiasTailoredLiftedProductCode(A::MatElem{T}, B::MatElem{T}; char_vec::Union{Vector{zzModRingElem},
    Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

    stabs = bias_tailored_lifted_product_matrices(A, B)
    return StabilizerCode(lift(stabs), char_vec = char_vec, logs_alg = logs_alg)
end

"""
    SPCDFoldProductCode(D::Int, s::Int = 1)
    SingleParityCheckDFoldProductCode(D::Int, s::Int = 1) = SPCDFoldProductCode(D, s)

Return the single-parity-check `D`-fold product code.

# Note
- This is defined in https://arxiv.org/abs/2209.13474

# Example

[512, 174, 8]] Symmetric 2-fold product CSS code from [ostrev2024classical](@cite)

```jldoctest
julia> using CodingTheory, Oscar;

julia> F = Oscar.Nemo.Native.GF(2);

julia> h = matrix(F, [1 1]);

julia> id = identity_matrix(F, 2);

julia> H_X = vcat(
             h ⊗ h ⊗ h ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id,
             id ⊗ id ⊗ id ⊗ h ⊗ h ⊗ h ⊗ id ⊗ id ⊗ id,
             id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ h ⊗ h ⊗ h);

julia> H_Z = vcat(
             h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id,
             id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id,
             id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h);

julia> code = SPCDFoldProductCode(3);

julia> length(code), dimension(code)
(512, 174)
```
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
    set_X_minimum_distance!(S, 2^D)
    set_Z_minimum_distance!(S, 2^D)
    S.pure = true
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
    
# Note
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
            lx[l:l, :] = hcat(tr_im_H1_tr_c[i:i, :] ⊗ tr_ker_H2[h:h, :], temp)
            lz[l:l, :] = hcat(tr_ker_H1[i:i, :] ⊗ tr_im_H2_tr_c[h:h, :], temp)
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
            lx[l:l, :] = hcat(temp, tr_ker_H1_tr[i:i, :] ⊗ tr_im_H2_c[h:h, :])
            lz[l:l, :] = hcat(temp, tr_im_H1_c[i:i, :] ⊗ tr_ker_H2_tr[h:h, :])
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

# Note
- This is defined in https://arxiv.org/abs/2209.13474
"""
asymmetric_product(S1::T, S2::T) where {T <: AbstractSubsystemCode} = asymmetric_product(
    CSSTrait(T), S1, S2)
function asymmetric_product(::IsCSS, S1::AbstractSubsystemCode, S2::AbstractSubsystemCode)
    # TODO: check fields are the same or convertable
    H_X = vcat(S1.X_stabs ⊗ identity_matrix(S1.F, S2.n), identity_matrix(S1.F, S1.n) ⊗ S2.X_stabs)
    H_Z = S1.Z_stabs ⊗ S2.Z_stabs
    return CSSCode(H_X, H_Z)
end
asymmetric_product(::IsNotCSS, S1::AbstractSubsystemCode, S2::AbstractSubsystemCode) =
    error("Only valid for CSS codes.")

"""
    symmetric_product(vec_S::Vector{T}) where {T <: AbstractSubsystemCode}

Return the symmetric `D`-fold product quantum CSS code, where `D` is
the square-root of the length of the vector of CSS codes `vec_S`.

# Note
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
symmetric_product(::IsNotCSS, vec_S::Vector{T}) where {T <: AbstractSubsystemCode} =
    error("Only valid for CSS codes.")

# has this been extended to subsystem codes?
"""
    homological_product(S1::AbstractStabilizerCode, S2::AbstractStabilizerCode, U::CTMatrixTypes = identity_matrix(S1.F, S1.n), V::CTMatrixTypes = identity_matrix(S2.F, S2.n))
    ⊠(S1::AbstractStabilizerCode, S2::AbstractStabilizerCode) = homological_product(S1, S2)

# Note
- This is the single-sector homological product. Use ⊗ for the more general product.
"""
function homological_product(S1::AbstractStabilizerCode, S2::AbstractStabilizerCode,
    U::CTMatrixTypes = identity_matrix(S1.F, S1.n), V::CTMatrixTypes = identity_matrix(S2.F, S2.n))

    return homological_product(CSSTrait(typeof(S1)), CSSTrait(typeof(S2)), S1, S2, U, V)
end

function homological_product(::IsCSS, ::IsCSS, S1::AbstractStabilizerCode, S2::AbstractStabilizerCode, U::CTMatrixTypes, V::CTMatrixTypes)

    num_stabs1 = num_X_stabs(S1)
    num_stabs1 == num_Z_stabs(S1) || throw(ArgumentError("The first code didn't have the same number of X and Z stabilizers"))
    num_stabs2 = num_X_stabs(S2)
    num_stabs2 == num_Z_stabs(S2) || throw(ArgumentError("The second code didn't have the same number of X and Z stabilizers"))
    nrows(U) == ncols(U) == S1.n || throw(ArgumentError("U is the wrong size for the code S1"))
    nrows(V) == ncols(V) == S2.n || throw(ArgumentError("V is the wrong size for the code S2"))
    isinvertible(U) || throw(ArgumentError("U must be invertible"))
    isinvertible(V) || throw(ArgumentError("V must be invertible"))
    F = S1.F
    F == S2.F == base_ring(U) == base_ring(V) || throw(ArgumentError("S1, S2, U, and V should all have the same base ring"))

    δ1 = zero_matrix(F, S1.n, S1.n)
    for i in 1:num_stabs1
        for j in 1:num_stabs1
            # TODO are X and Z in the correct order here?
            δ1 += U[i, j] * transpose(S1.Z_stabs[i, :]) * S1.X_stabs[j, :]
        end
    end

    δ2 = zero_matrix(F, S2.n, S2.n)
    for i in 1:num_stabs2
        for j in 1:num_stabs2
            # TODO are X and Z in the correct order here?
            δ2 += V[i, j] * transpose(S2.Z_stabs[i, :]) * S2.X_stabs[j, :]
        end
    end

    i1 = identity_matrix(F, S1.n)
    i2 = identity_matrix(F, S2.n)
    ∂ = δ1 ⊗ i2 - i1 ⊗ δ2

    return CSSCode(∂, transpose(∂))
end
homological_product(::IsNotCSS, ::IsNotCSS, S1::AbstractStabilizerCode, S2::AbstractStabilizerCode,
    U, V) = throw(ArgumentError("This is only defined for CSS codes"))
homological_product(::IsNotCSS, ::IsCSS, S1::AbstractStabilizerCode, S2::AbstractStabilizerCode, U,
    V) = throw(ArgumentError("This is only defined for CSS codes"))
homological_product(::IsCSS, ::IsNotCSS, S1::AbstractStabilizerCode, S2::AbstractStabilizerCode, U,
    V) = throw(ArgumentError("This is only defined for CSS codes"))
@doc (@doc homological_product)
⊠(S1::AbstractStabilizerCode, S2::AbstractStabilizerCode) = homological_product(S1, S2)

function _rand_single_sector_boundary(n::Int, k::Int)
    num_stabs = divexact(n - k, 2)
    U = _rand_invertible_matrix(GF(2), n)
    d0 = zero_matrix(GF(2), n, n)
    d0[k + 1:k + num_stabs, k + num_stabs + 1:end] = identity_matrix(GF(2), num_stabs)
    return U * d0 * inv(U)
end

"""
   random_homological_product_code(n1::Int, k1::Int, n2::Int, k2::Int)

Return a random homological product code.

# Note
- This implements the construction in https://arxiv.org/abs/1311.0885.
"""
function random_homological_product_code(n1::Int, k1::Int, n2::Int, k2::Int)
    d1 = _rand_single_sector_boundary(n1, k1)
    d2 = _rand_single_sector_boundary(n2, k2)
    i1 = identity_matrix(GF(2), nrows(d1))
    i2 = identity_matrix(GF(2), nrows(d2))
    d = d1 ⊗ i2 - i1 ⊗ d2
    return CSSCode(d, transpose(d))
end

"""
    BivariateBicycleCode(a::MPolyQuoRingElem{FqMPolyRingElem}, b::MPolyQuoRingElem{FqMPolyRingElem})

Return the bivariate bicycle code defined by the residue ring elements `a` and `b`.

# Note
- This is defined in https://arxiv.org/pdf/2308.07915

# Example

[[360, 12, ≤24]] Bivariate Bicycle Code from Table 3 of [bravyi2024high](@cite).

```jldoctest
julia> using CodingTheory, Oscar;

julia> S, (x, y) = polynomial_ring(Oscar.Nemo.Native.GF(2), [:x, :y]);

julia> l = 30; m = 6;

julia> R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]));

julia> a = R(x^9 + y + y^2);

julia> b = R(y^3 + x^25 + x^26);

julia> code = BivariateBicycleCode(a, b);

julia> length(code), dimension(code)
(360, 12)
```
"""
function BivariateBicycleCode(a::T, b::T) where T <: Union{MPolyQuoRingElem{FqMPolyRingElem}, MPolyQuoRingElem{fpMPolyRingElem}}
    R = parent(a)
    R == parent(b) || throw(DomainError("Polynomials must have the same parent."))
    F = base_ring(base_ring(a))
    order(F) == 2 || throw(DomainError("This code family is currently only defined over binary fields."))
    length(symbols(parent(a))) == 2 || throw(DomainError("Polynomials must be over two variables."))
    g = gens(modulus(R))
    length(g) == 2 || throw(DomainError("Residue rings must have only two generators."))

    m = -1
    l = -1
    for g1 in g
        exps = collect(exponents(g1))
        length(exps) == 2 || throw(ArgumentError("Moduli of the incorrect form."))
        iszero(exps[2]) || throw(ArgumentError("Moduli of the incorrect form."))
        !iszero(exps[1][1]) && !iszero(exps[1][2]) && throw(ArgumentError("Moduli of the incorrect form."))
        if iszero(exps[1][1])
            m = exps[1][2]
        else
            l = exps[1][1]
        end
    end

    x = matrix(F, [mod1(i + 1, l) == j ? 1 : 0 for i in 1:l, j in 1:l]) ⊗ identity_matrix(F, m)
    y = identity_matrix(F, l) ⊗ matrix(F, [mod1(i + 1, m) == j ? 1 : 0 for i in 1:m, j in 1:m])

    A = zero_matrix(F, l * m, l * m)
    for ex in exponents(lift(a))
        iszero(ex[1]) || iszero(ex[2]) || throw(ArgumentError("Polynomial `a` must not have any `xy` terms"))
        power, which = findmax(ex)
        if which == 1
            A += x^power
        elseif which == 2
            A += y^power
        end
    end

    B = zero_matrix(F, l * m, l * m)
    for ex in exponents(lift(b))
        iszero(ex[1]) || iszero(ex[2]) || throw(ArgumentError("Polynomial `b` must not have any `xy` terms"))
        power, which = findmax(ex)
        if which == 1
            B += x^power
        elseif which == 2
            B += y^power
        end
    end


    return CSSCode(hcat(A, B), hcat(transpose(B), transpose(A)))
end

"""
    CoprimeBivariateBicycleCode(a::MPolyQuoRingElem{FqMPolyRingElem}, b::MPolyQuoRingElem{FqMPolyRingElem})

Return the coprime bivariate bicycle code defined by the residue ring elements `a` and `b`.

# Note

- This is defined in https://arxiv.org/pdf/2408.10001.
"""
function CoprimeBivariateBicycleCode(a::ResElem, b::ResElem)
    R = parent(a)
    S = base_ring(a)
    R == parent(b) || throw(DomainError("Polynomials must have the same parent."))
    F = base_ring(S)
    order(F) == 2 || throw(DomainError("This code family is currently only defined over binary fields."))
    length(gens(S)) == 1 || throw(DomainError("Polynomials must be over one variable."))
    f = modulus(R)
    deg_P = degree(f)
    f == gen(S)^deg_P - 1 || throw(ArgumentError("Residue ring not of the form π^(l * m) - 1."))

    facs = Nemo.factor(deg_P)
    length(facs) == 2 || throw(ArgumentError("Residue ring not of the form π^(l * m) - 1."))
    k = collect(keys(facs.fac))
    v = collect(values(facs.fac))
    l = k[1]^v[1]
    m = k[2]^v[2]
    # l = 3
    # m = 7
    # println("l = $l, m = $m")
    # actually this check is guaranteed to be true since factor breaks it into its prime
    # factorization
    # gcd(l, m) == 1 || throw(ArgumentError("l and m must be coprime"))

    x = matrix(F, [mod1(i + 1, l) == j ? 1 : 0 for i in 1:l, j in 1:l]) ⊗ identity_matrix(F, m)
    y = identity_matrix(F, l) ⊗ matrix(F, [mod1(i + 1, m) == j ? 1 : 0 for i in 1:m, j in 1:m])

    P = x * y
    A = zero_matrix(F, deg_P, deg_P)
    exps = findall(i -> !is_zero(i), collect(coefficients(lift(a)))) .- 1
    for ex in exps
        A += P^ex
    end

    B = zero_matrix(F, deg_P, deg_P)
    exps = findall(i -> !is_zero(i), collect(coefficients(lift(b)))) .- 1
    for ex in exps
        B += P^ex
    end

    return CSSCode(hcat(A, B), hcat(transpose(B), transpose(A)))
end
