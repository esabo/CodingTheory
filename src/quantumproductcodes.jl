# need to make an AbstractQuantumLDPCCode
# then AbstractQuantumLDPCCSSCode
# then this goes below that
# maybe make an AbastractLDPCCode <: AbstractLinearCode
# then this can take in an LDPC code

# J. Tillich, G. Zémor. "Quantum LDPC codes with positive rate and minimum distance
# proportional to n^(1/2)". (2013) arXiv:0903.0566v2
mutable struct HypergraphProductCode <: AbstractHypergraphProductCode
    F::FqNmodFiniteField # base field (symplectic)
    E::FqNmodFiniteField # additive field
    n::Integer
    k::Union{Integer, Rational{BigInt}}
    d::Union{Integer, Missing}
    dx::Union{Integer, Missing}
    dz::Union{Integer, Missing}
    stabs::fq_nmod_mat
    Xstabs::fq_nmod_mat
    Zstabs::fq_nmod_mat
    C1::Union{LinearCode, Missing}
    C2::Union{LinearCode, Missing}
    signs::Vector{nmod}
    Xsigns::Vector{nmod}
    Zsigns::Vector{nmod}
    dualgens::fq_nmod_mat
    logspace::Union{fq_nmod_mat, Missing}
    logicals::Union{Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}, Missing}
    charvec::Vector{nmod}
    sCWEstabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    sCWEdual::Union{WeightEnumerator, Missing} # S^⟂
    overcomplete::Bool
    Lsigns::Union{Vector{nmod}, Missing}
end

"""
    HypergraphProductCode(C::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)

Return the (symmetric) hypergraph product code of `C`.

The hypergraph product is defined in "J. Tillich, G. Zémor. Quantum LDPC codes
with positive rate and minimum distance proportional to n^(1/2). (2013)
arXiv:0903.0566v2"

# Arguments
* `C`: a binary linear code
* `charvec`: a length `2n` vector with elements in the `Z/(4)`. The first `n`
  elements specify the exponents of the `X` phases and second `n` the exponents
  of the `Z` phases; a missing argument will be set to the all-zero vector

# Notes
* A `+1` phase should be entered as `0` since the character vector stores the
  exponents.
* Stabilizer signs are automatically computed given the character vector.
"""
function HypergraphProductCode(C::AbstractLinearCode, charvec::Union{Vector{nmod},
        Missing}=missing)

    Int(order(C.F)) == 2 || error("Hypergraph product codes are only defined for binary codes.")

    # perform product with H as is, but need to actually compute the parameters
    # of the transpose code (LinearCode(H^T, true)) in order to use formulas
    Htr = transpose(C.H)
    M = MatrixSpace(C.F, C.n, C.n)
    eye = M(1)
    Mtr = MatrixSpace(C.F, ncols(Htr), ncols(Htr))
    eyetr = Mtr(1)
    HX = hcat(C.H ⊗ eye, eyetr ⊗ Htr)
    HZ = hcat(eye ⊗ C.H, Htr ⊗ eyetr)
    S = directsum(HX, HZ)
    Sq2 = symplectictoquadratic(S)
    # use the below as a unit test
    # return CSSCode(HX, HZ, charvec)

    p = Int(characteristic(C.F))
    nsym = ncols(S)
    if !ismissing(charvec)
        nsyms == length(charvec) || error("The characteristic value is of incorrect length.")
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

Return the hypergraph product code of `C1` and `C2`.

The hypergraph product is defined in "J. Tillich, G. Zémor. Quantum LDPC codes
with positive rate and minimum distance proportional to n^(1/2). (2013)
arXiv:0903.0566v2"

# Arguments
* `C1`: a binary linear code
* `C2`: a binary linear code
* `charvec`: a length `2n` vector with elements in the `Z/(4)`. The first `n`
  elements specify the exponents of the `X` phases and second `n` the exponents
  of the `Z` phases; a missing argument will be set to the all-zero vector

# Notes
* A `+1` phase should be entered as `0` since the character vector stores the
  exponents.
* Stabilizer signs are automatically computed given the character vector.
"""
function HypergraphProductCode(C1::AbstractLinearCode, C2::AbstractLinearCode,
        charvec::Union{Vector{nmod}, Missing}=missing)

    Int(order(C1.F)) == 2 || error("Hypergraph product codes are only defined for binary codes.")
    Int(order(C2.F)) == 2 || error("Hypergraph product codes are only defined for binary codes.")

    # note that orthogonality of C1 and C2 here is not necessary because
    # HX * transpose(HZ) = C1.H \otimes transpose(C2.H) + C1.H \otimes transpose(C2.H) = 0
    # in characteristic 2
    H1tr = transpose(C1.H)
    H2tr = transpose(C2.H)
    Mn1 = MatrixSpace(C1.F, C1.n, C1.n)
    Mnk1 = MatrixSpace(C1.F, ncols(H1tr), ncols(H1tr))
    Mn2 = MatrixSpace(C2.F, C2.n, C2.n)
    Mnk2 = MatrixSpace(C2.F, ncols(H2tr), ncols(H2tr))
    HX = hcat(H1 ⊗ Mn2(1), Mnk1(1) ⊗ H2tr)
    HZ = hcat(Mn1(1) ⊗ H2, H1tr ⊗ Mnk2(1))
    S = directsum(HX, HZ)
    Sq2 = symplectictoquadratic(S)

    p = Int(characteristic(C1.F))
    nsym = ncols(S)
    if !ismissing(charvec)
        nsym == length(charvec) || error("The characteristic value is of incorrect length.")
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

# unable to yield quantum LDPC code families with non constant minimum distance
"""
    GeneralizedShorCode(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)
    BaconCasaccinoConstruction(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{nmod}, Missing}=missing)

Return the generalized Shor code of `C1` and `C2` with `C1⟂ ⊆ C2`.

The generalized Shor code is defined in "D. Bacon and A. Casaccino. Quantum
error correcting subsystem codes from two classical linear codes. (2006)
http://arxiv.org/abs/quant-ph/0610088"

# Arguments
* `C1`: a binary linear code
* `C2`: a binary linear code
* `charvec`: a length `2n` vector with elements in the `Z/(4)`. The first `n`
  elements specify the exponents of the `X` phases and second `n` the exponents
  of the `Z` phases; a missing argument will be set to the all-zero vector

# Notes
* A `+1` phase should be entered as `0` since the character vector stores the
  exponents.
* Stabilizer signs are automatically computed given the character vector.
"""
function GeneralizedShorCode(C1::AbstractLinearCode, C2::AbstractLinearCode,
        charvec::Union{Vector{nmod}, Missing}=missing)

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
    charvec::Union{Vector{nmod}, Missing}=missing) = GeneralizedShorCode(C1, C2,
    charvec)

"""
    HyperBicycleCodeCSS(a::Vector{fq_nmod_mat}, b::Vector{fq_nmod_mat}, χ::Int)

Return the CSS hyperbicycle construction of "Quantum Kronecker sum-product low-density
parity-check codes with finite rate".

Inputs:
- a: A vector of length `c` of binary `fq_nmod_mat` matrices of the same dimensions.
- b: A vector of length `c` of binary `fq_nmod_mat` matrices of the same dimensions,
  potentially different from those of `a`.
- χ: A strictly positive integer coprime with `c`.
"""
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
        k1, n1 == size(a[i]) || throw(ArgumentError("First set of imput matrices must all have the same dimensions."))
        F == base_ring(b[i]) || throw(ArgumentError("Inputs must share the same base ring."))
        k2, n2 == size(b[i]) || throw(ArgumentError("Second set of imput matrices must all have the same dimensions."))
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
    # not really worth it, just use standard technique from quantumcode.jl
end

"""
    HyperBicycleCode(a::Vector{fq_nmod_mat}, b::Vector{fq_nmod_mat}, χ::Int)

Return the hyperbicycle construction of "Quantum Kronecker sum-product low-density
parity-check codes with finite rate".

Inputs:
- a: A vector of length `c` of binary `fq_nmod_mat` matrices of the same dimensions.
- b: A vector of length `c` of binary `fq_nmod_mat` matrices of the same dimensions,
  potentially different from those of `a`.
- χ: A strictly positive integer coprime with `c`.
"""
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
        k1, n1 == size(a[i]) || throw(ArgumentError("First set of imput matrices must all have the same dimensions."))
        F == base_ring(b[i]) || throw(ArgumentError("Inputs must share the same base ring."))
        k2, n2 == size(b[i]) || throw(ArgumentError("Second set of imput matrices must all have the same dimensions."))
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
    return QuantumCode(G, true)

    # equations 41 and 42 of the paper give matrices from which the logicals may be chosen
    # not really worth it, just use standard technique from quantumcode.jl
end
