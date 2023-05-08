# Copyright (c) 2022, 2023 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

# TODO: _processcharvec, _determinesigns, remove dualgens etc from constructor
"""
    HypergraphProductCode(C::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)

Return the (symmetric) hypergraph product code of `C` whose signs are determined by `charvec`.

# Notes
* The hypergraph product is defined in "J. Tillich, G. Zémor. Quantum LDPC codes
  with positive rate and minimum distance proportional to n^(1/2). (2013)
  arXiv:0903.0566v2"
"""
function HypergraphProductCode(C::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)

    # Int(order(C.F)) == 2 || error("Hypergraph product codes are only defined for binary codes.")

    # perform product with H as is, but need to actually compute the parameters
    # of the transpose code (LinearCode(H^T, true)) in order to use formulas
    Htr = transpose(C.H)
    M = MatrixSpace(C.F, C.n, C.n)
    eye = M(1)
    Mtr = MatrixSpace(C.F, ncols(Htr), ncols(Htr))
    eyetr = Mtr(1)
    HX = hcat(C.H ⊗ eye, -eyetr ⊗ Htr)
    HZ = hcat(eye ⊗ C.H, Htr ⊗ eyetr)
    S = directsum(HX, HZ)
    Sq2 = symplectictoquadratic(S)
    # use the below as a unit test
    # return CSSCode(HX, HZ, charvec)

    p = Int(characteristic(C.F))
    nsym = ncols(S)
    if !ismissing(charvec)
        nsyms == length(charvec) || throw(ArgumentError("The characteristic value is of incorrect length."))
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        for s in charvec
            modulus(s) == modulus(R) || throw(ArgumentError("Phases are not in the correct ring."))
        end
    else
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        charvec = [R(0) for _ in 1:nsym]
    end

    # determine signs
    nrS = nrows(S)
    nr = nrows(HX)
    if iszero(charvec)
        signs = [R(0) for _ in 1:nrS]
        Xsigns = [R(0) for _ in 1:nr]
        Zsigns = [R(0) for _ in 1:nrows(HZ)]
    else
        signs = _getsigns(S, charvec)
        Xsigns = signs[1:nr, :]
        Zsigns = signs[nr + 1:end, :]
    end

    # find generators for S^⟂
    # note the H here is transpose of the standard definition
    # this line redefines the variable H, but that's fine, keep for consistency
    n = ncols(Sq2)
    _, H = right_kernel(hcat(S[:, n + 1:end], -S[:, 1:n]))
    # n + (n - k)
    ncols(H) == nsym - nrS || error("Normalizer matrix is not size n + k.")
    dualgens = symplectictoquadratic(transpose(hcat(H[:, n + 1:end], -H[:, 1:n])))

    # q^n / p^k but rows is n - k
    dimcode = BigInt(order(C.F))^n // BigInt(p)^nrS
    isinteger(dimcode) && (dimcode = Int(log(BigInt(p), dimcode));)

    # TODO: should I store this? I'm thinking not
    # this is an (l, q)-QLDPC code
    # l, q = columnrowweights(Sq2)

    return HypergraphProductCode(C.F, base_ring(Sq2), n, rank(S), C.d, missing,
        missing, Sq2, HX, HZ, C, missing, signs, Xsigns, Zsigns, dualgens, missing,
        missing, charvec, missing, missing, false, missing)

    # length C.n^2 + (n^T)^2
    # dimension C.k^2 + (k^T)^2,
    # minimum distance min(d, d^T)
    # stabilizers have row weights of the form i + j, where i and j are the
    # row and column weights of the H, respecitvely

end

"""
    HypergraphProductCode(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)

Return the hypergraph product code of `C1` and `C2` whose signs are determined by `charvec`.
"""
function HypergraphProductCode(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)
    # Int(order(C1.F)) == 2 || error("Hypergraph product codes are only defined for binary codes.")
    # Int(order(C2.F)) == 2 || error("Hypergraph product codes are only defined for binary codes.")

    # note that orthogonality of C1 and C2 here is not necessary because
    # HX * transpose(HZ) = C1.H \otimes transpose(C2.H) + C1.H \otimes transpose(C2.H) = 0
    # in characteristic 2
    H1tr = transpose(C1.H)
    H2tr = transpose(C2.H)
    Mn1 = MatrixSpace(C1.F, C1.n, C1.n)
    Mnk1 = MatrixSpace(C1.F, ncols(H1tr), ncols(H1tr))
    Mn2 = MatrixSpace(C2.F, C2.n, C2.n)
    Mnk2 = MatrixSpace(C2.F, ncols(H2tr), ncols(H2tr))
    HX = hcat(H1 ⊗ Mn2(1), -Mnk1(1) ⊗ H2tr)
    HZ = hcat(Mn1(1) ⊗ H2, H1tr ⊗ Mnk2(1))
    S = directsum(HX, HZ)
    Sq2 = symplectictoquadratic(S)

    p = Int(characteristic(C1.F))
    nsym = ncols(S)
    if !ismissing(charvec)
        nsym == length(charvec) || throw(ArgumentError("The characteristic value is of incorrect length."))
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        for s in charvec
            modulus(s) == modulus(R) || throw(ArgumentError("Phases are not in the correct ring."))
        end
    else
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        charvec = [R(0) for _ in 1:nsym]
    end

    # determine signs
    nrS = nrows(S)
    nr = nrows(HX)
    if iszero(charvec)
        signs = [R(0) for _ in 1:nrS]
        Xsigns = [R(0) for _ in 1:nr]
        Zsigns = [R(0) for _ in 1:nrows(HZ)]
    else
        signs = _getsigns(S, charvec)
        Xsigns = signs[1:nr, :]
        Zsigns = signs[nr + 1:end, :]
    end

    # find generators for S^⟂
    # note the H here is transpose of the standard definition
    n = ncols(Sq2)
    _, H = right_kernel(hcat(S[:, n + 1:end], -S[:, 1:n]))
    # n + (n - k)
    ncols(H) == nsym - nrS || error("Normalizer matrix is not size n + k.")
    dualgens = symplectictoquadratic(transpose(hcat(H[:, n + 1:end], -H[:, 1:n])))

    # q^n / p^k but rows is n - k
    dimcode = BigInt(order(C.F))^ncols(Sq2) // BigInt(p)^nrS
    isinteger(dimcode) && (dimcode = Int(log(BigInt(p), dimcode));)
    (ismissing(C1.d) || ismissing(C2.d)) ? d = missing : d = minimum([C1.d, C2.d])

    # TODO: should I store this? I'm thinking not
    # this is an (l, q)-QLDPC code
    # l, q = columnrowweights(Sq2)

    return HypergraphProductCode(C1.F, base_ring(Sq2), n, rank(S), d, missing,
        missing, Sq2, HX, HZ, C1, C2, signs, Xsigns, Zsigns, dualgens, missing,
        missing, charvec, missing, missing, false, missing)

    # [[(C1.n)^2 + (C2^T.n)^2, (C1.k)^2 + (C2^T.k)^2, min(d, d^T)]]
    # dX = min(d^T_1, d_2), dZ = min(d1, d^T_2), d = min(dX, dZ)
end

# TODO
# function HypergraphProductCode(A::fq_nmod_mat, B::fq_nmod_mat)

# end

# unable to yield quantum LDPC code families with non constant minimum distance
"""
    GeneralizedShorCode(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)
    BaconCasaccinoConstruction(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)

Return the generalized Shor code of `C1` and `C2` with `C1⟂ ⊆ C2` whose signs are determined by `charvec`.

# Notes
* The generalized Shor code is defined in "D. Bacon and A. Casaccino. Quantum
  error correcting subsystem codes from two classical linear codes. (2006)
  http://arxiv.org/abs/quant-ph/0610088"
"""
function GeneralizedShorCode(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)

    Int(order(C1.F)) == 2 || error("Generalized Shor codes are only defined for binary codes.")
    Int(order(C2.F)) == 2 || error("Generalized Shor codes are only defined for binary codes.")
    dual(C1) ⊆ C2 || error("Generalized Shor codes require the dual of the first code is a subset of the second.")

    # HX stays sparse if H1 is but HZ does not if C1.d is large
    Mn2 = MatrixSpace(C2.F, C2.n, C2.n)
    HX = C1.H ⊗ Mn2(1)
    HZ = C1.G ⊗ C2.H
    (ismissing(C1.d) || ismissing(C2.d)) ? d = missing : d = minimum([C1.d, C2.d])
    return CSSCode(HX, HZ, charvec)

    # [[n1 * n2, k1 * k2, min(d1, d2)]]
end
BaconCasaccinoConstruction(C1::AbstractLinearCode, C2::AbstractLinearCode,
    charvec::Union{Vector{nmod}, Missing}=missing) = GeneralizedShorCode(C1, C2, charvec)

"""
    HyperBicycleCodeCSS(a::Vector{fq_nmod_mat}, b::Vector{fq_nmod_mat}, χ::Int)

Return the hyperbicycle CSS code of `a` and `b` given `χ`.

# Arguments
* a: A vector of length `c` of binary matrices of the same dimensions.
* b: A vector of length `c` of binary matrices of the same dimensions,
  potentially different from those of `a`.
* χ: A strictly positive integer coprime with `c`.

# Notes
* Hyperbicycle codes are found in "Quantum ``hyperbicycle'' low-density parity check
  codes with finite rate" and "Quantum Kronecker sum-product low-density parity-check
  codes with finite rate".
"""
# TODO: charvec
function HyperBicycleCodeCSS(a::Vector{fq_nmod_mat}, b::Vector{fq_nmod_mat}, χ::Int)
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

    # julia creates a new scope for every iteration of the for loop,
    # so the else doesn't work without declaring these variables outside here
    H1 = missing
    H2 = missing
    HT1 = missing
    HT2 = missing
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

    Mk1 = MatrixSpace(F, k1, k1)
    Ek1 = Mk1(1)
    Mk2 = MatrixSpace(F, k2, k2)
    Ek2 = Mk2(1)
    Mn1 = MatrixSpace(F, n1, n1)
    En1 = Mn1(1)
    Mn2 = MatrixSpace(F, n2, n2)
    En2 = Mn2(1)

    GX = hcat(Ek2 ⊗ H1, H2 ⊗ Ek1)
    GZ = hcat(HT2 ⊗ En1, En2 ⊗ HT1)
    return CSSCode(GX, GZ)

    # equations 41 and 42 of the paper give matrices from which the logicals may be chosen
    # not really worth it, just use standard technique from StabilizerCode.jl
end

"""
    HyperBicycleCode(a::Vector{fq_nmod_mat}, b::Vector{fq_nmod_mat}, χ::Int)

Return the hyperbicycle CSS code of `a` and `b` given `χ`.

# Arguments
* a: A vector of length `c` of binary matrices of the same dimensions.
* b: A vector of length `c` of binary matrices of the same dimensions,
  potentially different from those of `a`.
* χ: A strictly positive integer coprime with `c`.

# Notes
* Hyperbicycle codes are found in "Quantum ``hyperbicycle'' low-density parity check
  codes with finite rate" and "Quantum Kronecker sum-product low-density parity-check
  codes with finite rate".
"""
# TODO: charvec
function HyperBicycleCode(a::Vector{fq_nmod_mat}, b::Vector{fq_nmod_mat}, χ::Int)
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

    Mk1 = MatrixSpace(F, k1, k1)
    Ek1 = Mk1(1)
    Mk2 = MatrixSpace(F, k2, k2)
    Ek2 = Mk2(1)

    G = hcat(Ek2 ⊗ H1, H2 ⊗ Ek1)
    return StabilizerCode(G, true)

    # equations 41 and 42 of the paper give matrices from which the logicals may be chosen
    # not really worth it, just use standard technique from StabilizerCode.jl
end

"""
    GeneralizedBicycleCode(A::fq_nmod_mat, B::fq_nmod_mat)

Return the generealized bicycle code given by `A` and `B`.

# Notes
* Generealized bicycle codes are are found in "Quantum ``hyperbicycle'' low-density parity check
  codes with finite rate", "Quantum kronecker sum-product low-density parity- check codes with
  finite rate", and "Degenerate Quantum LDPC Codes With Good Finite Length Performance".
"""
# TODO: charvec
function GeneralizedBicycleCode(A::fq_nmod_mat, B::fq_nmod_mat)
    base_ring(A) == base_ring(B) || throw(ArgumentError("Arguments must be over the same base ring."))
    (iszero(A) || iszero(B)) && throw(ArgumentError("Arguments should not be zero."))
    # this will take care of the sizes being square
    iszero(A * B - B * A) || throw(ArgumentError("Arguments must commute."))
    HX = hcat(A, B)
    HZ = hcat(transpose(B), -transpose(A))
    return CSSCode(HX, HZ)
end

"""
    GeneralizedBicycleCode(a::AbstractAlgebra.Generic.Res{fq_nmod_poly}, b::AbstractAlgebra.Generic.Res{fq_nmod_poly})

Return the generealized bicycle code determined by `a` and `b`.

# Notes
* `l x l` circulant matrices are constructed using the coefficients of the polynomials
  `a` and `b` in `F_q[x]/(x^l - 1)` (`gcd(q, l) = 1`) as the first column
* Generealized bicycle codes are are found in "Quantum ``hyperbicycle'' low-density parity check
  codes with finite rate", "Quantum kronecker sum-product low-density parity- check codes with
  finite rate", and "Degenerate Quantum LDPC Codes With Good Finite Length Performance".
"""
function GeneralizedBicycleCode(a::AbstractAlgebra.Generic.Res{fq_nmod_poly}, b::AbstractAlgebra.Generic.Res{fq_nmod_poly})
    parent(a) == parent(b) || throw(ArgumentError("Both objects must be defined over the same residue ring."))
    return GeneralizedBicycleCode(polytocircmatrix(a), polytocircmatrix(b))
end

# function BicycleCode(A::fq_nmot_mat)
#     m, n = size(A)
#     m == n || throw(ArgumentError("Input matrix must be square."))
#     # should probably check for F_2

#     H = hcat(A, transpose(A))
#     return CSSCode(H, H)
# end
    
"""
    GeneralizedHypergraphProductCode(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}}, b::AbstractAlgebra.Generic.Res{fq_nmod_poly})

Return the two matrices `HX` and `HZ` of generalized hypergraph product of `A` and `b`
over the reside ring.

# Arguments
* `A` - an `m x n` matrix with coefficents in a residue ring over `GF(2)`.
* `b` - a polynomial over the same residue ring

# Notes
* The generalized hypergraph product is defined in "Degenerate Quantum LDPC Codes With Good Finite Length Performance".
* To return a quantum code directly, use `LiftedGeneralizedHypergraphProductCode`.
"""
# TODO: charvec
function GeneralizedHypergraphProductCode(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}}, b::AbstractAlgebra.Generic.Res{fq_nmod_poly})
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
    LiftedGeneralizedHypergraphProductCode(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}}, b::AbstractAlgebra.Generic.Res{fq_nmod_poly})

Return the CSS code produced by lifting the generalized hypergraph product of `A` and `b`
over the underlying base field.
"""
# TODO: charvec
function LiftedGeneralizedHypergraphProductCode(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}}, b::AbstractAlgebra.Generic.Res{fq_nmod_poly})
    HX, HZ = GeneralizedHypergraphProductCode(A, b)
    return CSSCode(lift(HX), lift(HZ))
end

# TODO: charvec
function QuasiCyclicLiftedProductCode(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}},
    B::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}})

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

    Mk1 = MatrixSpace(R, k1, k1)
    Ek1 = Mk1(1)
    Mk2 = MatrixSpace(R, k2, k2)
    Ek2 = Mk2(1)
    Mn1 = MatrixSpace(R, n1, n1)
    En1 = Mn1(1)
    Mn2 = MatrixSpace(R, n2, n2)
    En2 = Mn2(1)

    HX = hcat(kronecker_product(A, Ek2), kronecker_product(Ek1, B))
    HZ = hcat(kronecker_product(En1, Btr), kronecker_product(Atr, En2))
    return HX, HZ
end

# TODO: charvec
function LiftedQuasiCyclicLiftedProductCode(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}},
    B::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}})

    HX, HZ = QuasiCyclicLiftedProductCode(A, B)
    return CSSCode(lift(HX), lift(HZ))
end

"""
    BiasTailoredQuasiCyclicLiftedProductCode(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}},
B::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}})

Return the bias-tailored lifted product code of `A` and `B` with entries over the residue ring
`F_2[x]/(x^m - 1)`.

# Notes
* The bias-tailored lifted product is defined in `Bias-tailored quantum LDPC codes`.
"""
# TODO: charvec
function BiasTailoredQuasiCyclicLiftedProductCode(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}},
    B::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}})

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

    Mk1 = MatrixSpace(R, k1, k1)
    Ek1 = Mk1(1)
    Mk2 = MatrixSpace(R, k2, k2)
    Ek2 = Mk2(1)
    Mn1 = MatrixSpace(R, n1, n1)
    En1 = Mn1(1)
    Mn2 = MatrixSpace(R, n2, n2)
    En2 = Mn2(1)

    A12 = kronecker_product(Atr, Ek2)
    A13 = kronecker_product(En1, B)
    A21 = kronecker_product(A, En2)
    A24 = kronecker_product(Ek1, Btr)
    return vcat(hcat(zeros(A21), A12, A13, zeros(A24)), hcat(A21, zeros(A12), zeros(A13), A24))
end

# TODO: charvec
function LiftedBiasTailoredQuasiCyclicLiftedProductCode(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}},
    B::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}})

    S = BiasTailoredQuasiCyclicLiftedProductCode(A, B)
    return StabilizerCode(lift(S), true)
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
