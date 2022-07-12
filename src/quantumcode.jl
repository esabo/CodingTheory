# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

mutable struct CSSCode <: AbstractCSSCode
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
    Xorigcode::Union{LinearCode, Missing}
    ZorigCode::Union{LinearCode, Missing}
    signs::Vector{nmod} # make bools, uint8's?
    Xsigns::Vector{nmod} # make bools, uint8's?
    Zsigns::Vector{nmod} # make bools, uint8's?
    dualgens::fq_nmod_mat
    logspace::Union{fq_nmod_mat, Missing}
    logicals::Union{Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}, Missing}
    charvec::Vector{nmod} # make bools, uint8's?
    sCWEstabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    sCWEdual::Union{WeightEnumerator, Missing} # S^⟂
    overcomplete::Bool
end

mutable struct QuantumCode <: AbstractStabilizerCode
    F::FqNmodFiniteField # base field (symplectic)
    E::FqNmodFiniteField # additive field
    n::Integer
    k::Union{Integer, Rational{BigInt}}
    d::Union{Integer, Missing}
    stabs::fq_nmod_mat
    dualgens::fq_nmod_mat
    logspace::Union{fq_nmod_mat, Missing}
    logicals::Union{Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}, Missing}
    charvec::Vector{nmod} # make bools, uint8's?
    signs::Vector{nmod}
    sCWEstabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    sCWEdual::Union{WeightEnumerator, Missing} # S^⟂
    overcomplete::Bool
end

"""
    field(S::AbstractStabilizerCode)

Return the base ring of the code as a Nemo object.
"""
field(S::AbstractStabilizerCode) = S.F

"""
    quadraticfield(S::AbstractStabilizerCode)

Return the quadratic field of the code as a Nemo object.
"""
quadraticfield(S::AbstractStabilizerCode) = S.E

"""
    length(S::AbstractStabilizerCode)
    numqubits(S::AbstractStabilizerCode)

Return the length of the code.
"""
length(S::AbstractStabilizerCode) = S.n
numqubits(S::AbstractStabilizerCode) = S.n

"""
    dimension(S::AbstractStabilizerCode)

Return the dimension of the code.
"""
dimension(S::AbstractStabilizerCode) = S.k

"""
    cardinality(S::AbstractStabilizerCode)

Return the cardinality of the stabilizer group of the code.

No size checking is done on the parameters of the code, returns a BitInt by default.
"""
cardinality(S::AbstractStabilizerCode) = BigInt(order(field(C)))^dimension(C)

"""
    rate(S::AbstractStabilizerCode)

Return the rate, `R = k/n', of the code.
"""
rate(S::AbstractStabilizerCode) = dimension(S) / length(S)

"""
    signs(S::AbstractStabilizerCode)

Return the signs of the stabilizers of the code.
"""
signs(S::AbstractStabilizerCode) = S.signs

"""
    Xsigns(S::CSSCode)

Return the signs of the `X` stabilizers of the CSS code.
"""
Xsigns(S::CSSCode) = S.Xsigns

"""
    Zsigns(S::CSSCode)

Return the signs of the `Z` stabilizers of the CSS code.
"""
Zsigns(S::CSSCode) = S.Zsigns

"""
    stabilizers(S::AbstractStabilizerCode)

Return the stabilizer matrix of the code.
"""
stabilizers(S::AbstractStabilizerCode) = S.stabs

"""
    symplecticstabilizers(S::AbstractStabilizerCode)

Return the stabilizer matrix of the code in symplectic form.
"""
symplecticstabilizers(S::AbstractStabilizerCode) = quadratictosymplectic(stabilizers(S))

"""
    Xstabilizers(S::CSSCode)

Return the `X` stabilizer matrix of the CSS code.
"""
Xstabilizers(S::CSSCode) = S.Xstabs

"""
    Zstabilizers(S::CSSCode)

Return the `Z` stabilizer matrix of the CSS code.
"""
Zstabilizers(S::CSSCode) = S.Zstabs

numXstabs(S::CSSCode) = nrows(S.Xstabs)
numZstabs(S::CSSCode) = nrows(S.Zstabs)

"""
    normalizermatrix(S::AbstractStabilizerCode)

Return the normalizer matrix of the code.
"""
normalizermatrix(S::AbstractStabilizerCode) = S.dualgens

"""
    charactervector(S::AbstractStabilizerCode)

Return the character vector of the code.
"""
charactervector(S::AbstractStabilizerCode) = S.charvec

"""
    setminimumdistance(S::AbstractStabilizerCode, d::Integer)

Set the minimum distance of the code to `d`.

The only check done on the value of `d` is that `1 ≤ d ≤ n`.
"""
function setminimumdistance!(S::AbstractStabilizerCode, d::Integer)
    # TODO: should check bounds like Singleton for possibilities
    d > 0 && d <= length(S) || error("The minimum distance of a code must be ≥ 1; received: d = $d.")
    S.d = d
end

"""
    relativedistance(S::AbstractStabilizerCode)

Return the relative minimum distance, `δ = d / n` of the code if `d` is known,
otherwise errors.
"""
function relativedistance(S::AbstractStabilizerCode)
    !ismissing(S.d) || error("Missing minimum distance for this code.")
    return S.d / S.n
end

"""
    isovercomplete(S::AbstractStabilizerCode)

Return `true` if `S` has an overcomplete set of stabilizers.
"""
isovercomplete(S::AbstractStabilizerCode) = return S.overcomplete

function _getsigns(A::fq_nmod_mat, charvec::Vector{nmod})
    R = base_ring(charvec[1])
    if length(charvec) == 2 * ncols(A)
        A = quadratictosymplectic(A)
    end
    length(charvec) == ncols(A) || error("Input to _getsigns is expected to be in symplectic form and of the same length as the characteristic vector.")
    if iszero(charvec)
        return [R(0) for _ in 1:div(ncols(A), 2)]
    end

    signs = Vector{Int64}()
    for r in 1:nrows(A)
        parity = R(0)
        for c = 1:ncols(A)
            parity += charvec[c]
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
    splitstabilizers(S::AbstractStabilizerCode)

Return the set of `X`-only stabilizers and their signs, the set of `Z`-only
stabilizers and their signs, and the remaining stabilizers and their signs.

This function returns six objevts of alternating types `fq_nmod_mat` and
`Vector{Int}` for the three sets of stabilizers and signs, respectively.
An empty set of stabilizers is returned as type `Vector{fq_nmod_mat}`.
"""
function splitstabilizers(S::AbstractStabilizerCode)
    return _splitsymplecticstabilizers(symplecticstabilizers(S), signs(S))
end

function _isCSSsymplectic(S::fq_nmod_mat, signs::Vector{nmod}, trim::Bool=true)
    Xstabs, Xsigns, Zstabs, Zsigns, mixedstabs, mixedsigns = _splitsymplecticstabilizers(S, signs)
    if typeof(mixedstabs) <: Vector{fq_nmod_mat}
        if trim
            half = div(ncols(Xstabs), 2)
            return true, Xstabs[:, 1:half], Xsigns, Zstabs[:, half + 1:end], Zsigns
        else
            return true, Xstabs, Xsigns, Zstabs, Zsigns
        end
    else
        if trim
            half = div(ncols(Xstabs), 2)
            return false, Xstabs[:, 1:half], Xsigns, Zstabs[:, half + 1:end], Zsigns, mixedstabs, mixedsigns
        else
            return false, Xstabs, Xsigns, Zstabs, Zsigns, mixedstabs, mixedsigns
        end
    end
end

"""
    isCSS(S::AbstractStabilizerCode)

Return `true` is `S` is CSS.

# Notes
* This is intended to be a simple function wrapper for `typeof(S)` since the
 constructor for `QuantumCode` automatically returns a `CSSCode` if possible.
 Manually changing the elements of the struct `S` without using the helper
 functions provided here is not recommended.
"""
function isCSS(S::AbstractStabilizerCode)
    typeof(S) <: AbstractCSSCode && return true
    return false
end
"""
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
function CSSCode(C1::AbstractLinearCode, C2::AbstractLinearCode,
    charvec::Union{Vector{nmod}, Missing}=missing)

    C2 ⊆ C1 || error("The second argument must be a subset of the first in the CSS construction.")
    if !ismissing(charvec)
        2 * length(C1) == length(charvec) || error("The characteristic value is of incorrect length.")
        p = Int(characteristic(field(C1)))
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        for s in charvec
            modulus(s) == modulus(R) || error("Phases are not in the correct ring.")
        end
    else
        p = Int(characteristic(field(C1)))
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        charvec = [R(0) for _ in 1:2 * length(C1)]
    end

    # C2 ⊆ C1
    # k = k1 - k2
    # d >= minimum(d1, d2^⟂)
    # X - H(C2^⟂), Z - H(C1)
    D2 = dual(C2)
    S = directsum(paritycheckmatrix(D2), paritycheckmatrix(C1))
    Sq2 = symplectictoquadratic(S)
    dualgens = directsum(generatormatrix(C1), generatormatrix(D2))

    # determine signs
    if iszero(charvec)
        signs = [R(0) for _ in 1:nrows(S)]
        Xsigns = [R(0) for _ in 1:nrows(paritycheckmatrix(D2))]
        Zsigns = [R(0) for _ in 1:nrows(paritycheckmatrix(C1))]
    else
        signs = _getsigns(S, charvec)
        Xsigns = signs[1:nrows(paritycheckmatrix(D2)), :]
        Zsigns = signs[nrows(paritycheckmatrix(D2)) + 1:end, :]
    end

    # q^n / p^k but rows is n - k
    dimcode = BigInt(order(base_ring(S)))^ncols(Sq2) // BigInt(characteristic(base_ring(S)))^rank(S)
    if isinteger(dimcode)
        dimcode = Int64(log(BigInt(characteristic(base_ring(S))), dimcode))
    end

    return CSSCode(field(C1), base_ring(Sq2), length(C1), dimcode, missing,
        minimumdistance(D2), minimumdistance(C1), Sq2, paritycheckmatrix(D2),
        paritycheckmatrix(C1), C2, C1, signs, Xsigns, Zsigns, dualgens, missing,
        missing, charvec, missing, missing, false)
end

"""
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
  `chracteristic(field(C1))` is 2 and `Z/(p)` otherwise. The first `n` elements
  specify the exponents of the `X` phases and second `n` the exponents of the
  `Z` phases; a missing argument will be set to the all-zero vector

# Notes
* A `+1` phase should be entered as `0` since the character vector stores the
  exponents.
* Stabilizer signs are automatically computed given the character vector.
"""
function CSSCode(C::LinearCode, charvec::Union{Vector{nmod}, Missing}=missing)
    # this should have Xstabs = Zstabs
    D = dual(C)
    C ⊆ D || error("The single code CSS construction requires C ⊆ C^⟂.")
    if !ismissing(charvec)
        2 * length(C) == length(charvec) || error("The characteristic value is of incorrect length.")
        p = Int(characteristic(field(C)))
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        for s in charvec
            modulus(s) == modulus(R) || error("Phases are not in the correct ring.")
        end
    else
        p = Int(characteristic(field(C1)))
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        charvec = [R(0) for _ in 1:2 * length(C)]
    end

    # C2 ⊆ C1
    # k = k1 - k2
    # d >= minimum(d1, d2^⟂)
    # X - H(C2^⟂), Z - H(C1)
    S = directsum(paritycheckmatrix(D), paritycheckmatrix(D))
    Sq2 = symplectictoquadratic(S)
    dualgens = directsum(generatormatrix(D), generatormatrix(D))

    # determine signs
    if iszero(charvec)
        signs = [R(0) for _ in 1:nrows(S)]
        Xsigns = [R(0) for _ in 1:nrows(paritycheckmatrix(D))]
        Zsigns = [R(0) for _ in 1:nrows(paritycheckmatrix(D))]
    else
        signs = _getsigns(S, charvec)
        Xsigns = signs[1:nrows(paritycheckmatrix(D)), :]
        Zsigns = signs[nrows(paritycheckmatrix(D)) + 1:end, :]
    end

    # q^n / p^k but rows is n - k
    dimcode = BigInt(order(base_ring(S)))^ncols(Sq2) // BigInt(characteristic(base_ring(S)))^rank(S)
    if isinteger(dimcode)
        dimcode = Int64(log(BigInt(characteristic(base_ring(S))), dimcode))
    end

    return CSSCode(field(C), base_ring(Sq2), length(C), dimcode, missing,
        minimumdistance(D), minimumdistance(D), Sq2, paritycheckmatrix(D),
        paritycheckmatrix(D), C, D, signs, Xsigns, Zsigns, dualgens, missing,
        missing, charvec, missing, missing, overcomp)
end

"""
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
function CSSCode(Xmatrix::fq_nmod_mat, Zmatrix::fq_nmod_mat,
    charvec::Union{Vector{nmod}, Missing}=missing)

    iszero(Xmatrix) && error("The `X` stabilizer matrix is empty.")
    iszero(Zmatrix) && error("The `Z` stabilizer matrix is empty.")
    n = ncols(Xmatrix)
    n ==  ncols(Zmatrix) || error("Both matrices must have the same length in the CSS construction.")
    base_ring(Xmatrix) == base_ring(Zmatrix) || error("Both matrices must be over the same base field in the CSS construction.")
    iszero(Zmatrix * transpose(Xmatrix)) || error("The given matrices are not symplectic orthogonal.")
    if !ismissing(charvec)
        2 * n == length(charvec) || error("The characteristic value is of incorrect length.")
        p = Int(characteristic(base_ring(Xmatrix)))
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        for s in charvec
            modulus(s) == modulus(R) || error("Phases are not in the correct ring.")
        end
    else
        p = Int(characteristic(base_ring(Xmatrix)))
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        charvec = [R(0) for _ in 1:2 * n]
    end

    # remove any empty rows
    _removeempty(Xmatrix, "rows")
    _removeempty(Zmatrix, "rows")

    # determine if the provided set of stabilizers are redundant
    Xrank = rank(Xmatrix)
    Zrank = rank(Zmatrix)
    if nrows(Xmatrix) > Xrank || nrows(Zmatrix) > Zrank
        overcomp = true
    else
        overcomp = false
    end

    S = directsum(Xmatrix, Zmatrix)
    # aresymplecticorthogonal(S, S, true) || error("The given matrices are not symplectic orthogonal.")
    Sq2 = symplectictoquadratic(S)

    # determine signs
    if iszero(charvec)
        signs = [R(0) for _ in 1:nrows(S)]
        Xsigns = [R(0) for _ in 1:nrows(Xmatrix)]
        Zsigns = [R(0) for _ in 1:nrows(Zmatrix)]
    else
        signs = _getsigns(S, charvec)
        Xsigns = signs[1:nrows(Xmatrix), :]
        Zsigns = signs[nrows(Xmatrix) + 1:end, :]
    end

    # find generators for S^⟂
    G = hcat(S[:, n + 1:end], -S[:, 1:n])
    _, H = right_kernel(G)
    ncols(H) == n + k || error("Normalizer matrix is not size n + k.")
    # note the H here is transpose of the standard definition
    iszero(G * H) || error("Normalizer matrix is not symplectic orthogonal.")
    dualgens = symplectictoquadratic(transpose(H))

    # q^n / p^k but rows is n - k
    dimcode = BigInt(order(base_ring(S)))^ncols(Sq2) // BigInt(characteristic(base_ring(S)))^(Xrank + Zrank)
    if isinteger(dimcode)
        dimcode = Int64(log(BigInt(characteristic(base_ring(S))), dimcode))
    end

    return CSSCode(base_ring(Xmatrix), base_ring(Sq2), n, dimcode, missing,
        missing, missing, Sq2, Xmatrix, Zmatrix, missing, missing, signs, Xsigns,
        Zsigns, dualgens, missing, missing, charvec, missing, missing, overcomp)
end

"""
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
function CSSCode(SPauli::Vector{T}, charvec::Union{Vector{nmod},
    Missing}=missing) where T <: Union{String, Vector{Char}}

    SPaulistripped, charvec = _processstrings(SPauli, charvec)
    S = _Paulistringtosymplectic(SPaulistripped)
    iszero(S) && error("The processed Pauli strings returned a set of empty stabilizers.")
    aresymplecticorthogonal(S, S, true) || error("The given stabilizers are not symplectic orthogonal.")
    if !ismissing(charvec)
        ncols(S) == length(charvec) || error("The characteristic value is of incorrect length.")
        p = Int(characteristic(base_ring(S)))
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        for s in charvec
            modulus(s) == modulus(R) || error("Phases are not in the correct ring.")
        end
    else
        p = Int(characteristic(base_ring(S)))
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        charvec = [R(0) for _ in 1:ncols(S)]
    end

    # remove any empty rows
    _removeempty!(S, "rows")
    Sq2 = symplectictoquadratic(S)

    # determine if the provided set of stabilizers are redundant
    Srank = rank(S)
    if nrows(S) > Srank
        overcomp = true
    else
        overcomp = false
    end

    # determine signs
    if iszero(charvec)
        signs = [R(0) for _ in 1:ncols(S)]
    else
        signs = _getsigns(S, charvec)
    end

    # find generators for S^⟂
    G = hcat(S[:, ncols(Sq2) + 1:end], -S[:, 1:ncols(Sq2)])
    _, H = right_kernel(G)
    # why'd I do 2n - k instead of n + k here?
    ncols(H) == 2 * ncols(Sq2) - nrows(Sq2) || error("Normalizer matrix is not size n + k.")
    # note the H here is transpose of the standard definition
    iszero(G * H) || error("Normalizer matrix is not symplectic orthogonal.")
    dualgens = symplectictoquadratic(transpose(H))

    # q^n / p^k but rows is n - k
    dimcode = BigInt(order(base_ring(S)))^ncols(Sq2) // BigInt(characteristic(base_ring(S)))^Srank
    if isinteger(dimcode)
        dimcode = Int64(log(BigInt(characteristic(base_ring(S))), dimcode))
    end

    args = _isCSSsymplectic(S, signs, true)
    if args[1]
        return CSSCode(base_ring(S), base_ring(Sq2), ncols(Sq2), dimcode, missing,
            missing, missing, Sq2, args[2], args[4], missing, missing, signs,
            args[3], args[5], dualgens, missing, missing, charvec, missing,
            missing, overcomp)
    else
        error("Provided Pauli strings are not CSS.")
    end
end

# entanglement-assisted is not symplectic orthogonal
"""
    QuantumCode(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}

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
function QuantumCode(SPauli::Vector{T}, charvec::Union{Vector{nmod}, Vector{Any}}=[]) where T <: Union{String, Vector{Char}}

    SPaulistripped, charvec = _processstrings(SPauli, charvec)
    S = _Paulistringtosymplectic(SPaulistripped)
    iszero(S) && error("The processed Pauli strings returned a set of empty stabilizers.")
    aresymplecticorthogonal(S, S, true) || error("The given stabilizers are not symplectic orthogonal.")
    if !ismissing(charvec)
        ncols(S) == length(charvec) || error("The characteristic value is of incorrect length.")
        p = Int(characteristic(base_ring(S)))
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        for s in charvec
            modulus(s) == modulus(R) || error("Phases are not in the correct ring.")
        end
    else
        p = Int(characteristic(base_ring(S)))
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        charvec = [R(0) for _ in 1:ncols(S)]
    end

    # remove any empty rows
    _removeempty!(S, "rows")
    Sq2 = symplectictoquadratic(S)

    # determine if the provided set of stabilizers are redundant
    Srank = rank(S)
    if size(S, 1) > Srank
        overcomp = true
    else
        overcomp = false
    end

    # determine signs
    if iszero(charvec)
        signs = [R(0) for _ in 1:nrows(S)]
    else
        signs = _getsigns(S, charvec)
    end

    # find generators for S^⟂
    G = hcat(S[:, ncols(Sq2) + 1:end], -S[:, 1:ncols(Sq2)])
    _, H = right_kernel(G)
    # why'd I do 2n - k instead of n + k here?
    ncols(H) == 2 * ncols(Sq2) - nrows(Sq2) || error("Normalizer matrix is not size n + k.")
    # note the H here is transpose of the standard definition
    iszero(G * H) || error("Normalizer matrix is not symplectic orthogonal.")
    dualgens = symplectictoquadratic(transpose(H))

    # q^n / p^k but rows is n - k
    dimcode = BigInt(order(base_ring(S)))^ncols(Sq2) // BigInt(characteristic(base_ring(S)))^Srank
    if isinteger(dimcode)
        dimcode = Int64(log(BigInt(characteristic(base_ring(S))), dimcode))
    end

    args = _isCSSsymplectic(S, signs, true)
    if args[1]
        return CSSCode(base_ring(S), base_ring(Sq2), ncols(Sq2), dimcode,
            missing, missing, missing, Sq2, args[2], args[4], missing, missing, signs,
            args[3], args[5], dualgens, missing, missing, charvec, missing, missing,
            overcomp)
    else
        return QuantumCode(base_ring(S), base_ring(Sq2), ncols(Sq2), dimcode,
            missing, Sq2, dualgens, missing, missing, charvec, signs, missing,
            missing, overcomp)
    end
end

"""
    QuantumCode(Sq2::fq_nmod_mat, symp::Bool=false, charvec::Union{Vector{nmod}, Missing}=missing)

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
function QuantumCode(Sq2::fq_nmod_mat, symp::Bool=false,
    charvec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}

    iszero(Sq2) && error("The stabilizer matrix is empty.")
    # remove any empty rows
    _removeempty!(Sq2, "rows")

    if symp
        S = Sq2
        # this will error properly if not correct
        Sq2 = symplectictoquadratic(Sq2)
    else
        E = base_ring(Sq2)
        iseven(degree(E)) || error("The base ring of the given matrix is not a quadratic extension.")
        S = quadratictosymplectic(Sq2)
    end
    aresymplecticorthogonal(S, S, symp) || error("The given stabilizers are not symplectic orthogonal.")

    if !ismissing(charvec)
        ncols(S) == length(charvec) || error("The characteristic value is of incorrect length.")
        p = Int(characteristic(base_ring(Sq2)))
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        for s in charvec
            modulus(s) == modulus(R) || error("Phases are not in the correct ring.")
        end
    else
        p = Int(characteristic(base_ring(Sq2)))
        if p == 2
            R = ResidueRing(Nemo.ZZ, 4)
        else
            R = ResidueRing(Nemo.ZZ, p)
        end
        charvec = [R(0) for _ in 1:ncols(S)]
    end

    # determine if the provided set of stabilizers are redundant
    Srank = rank(S)
    if ncols(S) > Srank
        overcomp = true
    else
        overcomp = false
    end

    # determine signs
    if iszero(charvec)
        signs = [R(0) for _ in 1:nrows(S)]
    else
        signs = _getsigns(S, charvec)
    end

    # find generators for S^⟂
    G = hcat(S[:, ncols(Sq2) + 1:end], -S[:, 1:ncols(Sq2)])
    _, H = right_kernel(G)
    # why'd I do 2n - k instead of n + k here?
    ncols(H) == 2 * ncols(Sq2) - nrows(Sq2) || error("Normalizer matrix is not size n + k.")
    # note the H here is transpose of the standard definition
    iszero(G * H) || error("Normalizer matrix is not symplectic orthogonal.")
    dualgens = symplectictoquadratic(transpose(H))

    # q^n / p^k but rows is n - k
    dimcode = BigInt(order(base_ring(S)))^ncols(Sq2) // BigInt(characteristic(base_ring(S)))^Srank
    if isinteger(dimcode)
        dimcode = Int64(log(BigInt(characteristic(base_ring(S))), dimcode))
    end

    args = _isCSSsymplectic(S, signs, true)
    if args[1]
        return CSSCode(base_ring(S), base_ring(Sq2), ncols(Sq2), dimcode,
            missing, missing, missing, Sq2, args[2], args[4], missing, missing, signs,
            args[3], args[5], dualgens, missing, missing, charvec, missing, missing,
            overcomp)
    else
        return QuantumCode(base_ring(S), base_ring(Sq2), ncols(Sq2), dimcode,
            missing, Sq2, dualgens, missing, missing, charvec, signs, missing,
            missing, overcomp)
    end
end

"""
    logicalspace(S::AbstractStabilizerCode)

Return a basis for the space `S^⟂ / S`.

Note that this is not necessarily a symplectic basis with `{X_L, Z_L}` pairs.
"""
function logicalspace(S::AbstractStabilizerCode)
    if !ismissing(S.logspace)
        return S.logspace
    end

    F = quadraticfield(S)
    G = stabilizers(S)
    Gdual = normalizermatrix(S)
    V = VectorSpace(F, length(S))
    U, UtoV = sub(V, [V(G[i, :]) for i in 1:nrows(G)])
    W, WtoV = sub(V, [V(Gdual[i, :]) for i in 1:nrows(Gdual)])
    gensofUinW = [preimage(WtoV, UtoV(g)) for g in gens(U)]
    UinW, UinWtoW = sub(W, gensofUinW)
    Q, WtoQ = quo(W, UinW)
    C2modC1basis = [WtoV(x) for x in [preimage(WtoQ, g) for g in gens(Q)]]
    Fbasis = [[F(C2modC1basis[j][i]) for i in 1:AbstractAlgebra.dim(parent(C2modC1basis[1]))] for j in 1:length(C2modC1basis)]
    G2 = matrix(F, length(Fbasis), length(Fbasis[1]), vcat(Fbasis...))

    # returns linear object, make additive if necessary
    if nrows(G2) == dimension(S)
        S.logspace = vcat(G2, gen(F) * G2)
        return S.logspace
    elseif nrows(G2) == 2 * dimension(S)
        S.logspace = G2
        return S.logspace
    else
        error("Logical space produced of incorrect dimension; expected: ",
            2 * dimension(S), ", received: ", nrows(G2))
    end
end

"""
    setlogicals!(S::AbstractStabilizerCode, L::fq_nmod_mat, symp::Bool=false)

Set the logical operators of `S` to `L`.

If the optional parameter `symp` is set to `true`, `L` is assumed to be in
symplectic form over the base field of `S`.
"""
function setlogicals!(S::AbstractStabilizerCode, L::fq_nmod_mat, symp::Bool=false)
    F = base_ring(L)
    if symp
        size(L) == (2 * dimension(S), 2 * length(S)) || error("Provided matrix is of incorrect size for the logical space.")
        iseven(ncols(L)) || error("Expected a symplectic input but the input matrix has an odd number of columns.")
        field(S) == F || error("The logicals must be over the same field as the code.")
    else
        size(L) == (2 * dimension(S), length(S)) || error("Provided matrix is of incorrect size for the logical space.")
        quadraticfield(Sq2) == F || error("The logicals must be over the same field as the code.")
        iseven(F) || error("The base ring of the given matrix is not a quadratic extension.")
        L = quadratictosymplectic(L)
    end
    aresymplecticorthogonal(symplecticstabilizers(S), L) || error("Provided logicals do not commute with the code.")

    # the columns in prod give the commutation relationships between the provided
    # logical operators; they ideally should only consist of {X_1, Z_i} pairs
    # so there should only be one nonzero element in each column
    prod = hcat(A[:, length(S) + 1:end], -A[:, 1:length(S)]) * transpose(L)
    iszero(prod) || error("Provided logical should not be symplectic self-orthogonal.")
    cols = sum(prod, 2)
    sum(cols) == length(cols) || error("Incorrect commutation relationships between provided logicals.")

    # pairs are row i, and then whatever column is nonzero, and then shift such that it is one
    logs = Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}()
    while nrows(L) >= 2
        y = findfirst(x->x>0, prod[:, 1])
        if y[1] != F(1)
            push!(logs, (L[1, :], y[1]^-1 * L[y[2], :]))
        else
            push!(logs, (L[1, :], L[y[2], :]))
        end
        L = L[setdiff(1:size(L, 1), [1, y[2]]), :]
    end

    # # can actually do this via the matrix multiplication above
    # logs = Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}()
    # while !iszero(size(L, 1))
    #     rows = []
    #     X1 = L[1, :]
    #     for i in 2:size(L, 1)
    #         SIP = sympecticinnerproduct(Xi, L[i, :])
    #         if !iszero(SIP)
    #             if SIP != typeof(SIP)(1)
    #                 L[i, :] *= SIP^-1
    #                 push!(rows, i)
    #                 length(rows) == 1 || error("Incorrect commutation relationships between provided logicals.")
    #             end
    #         end
    #     end
    #     push!(logs, (X1, L[rows[1], :]))
    #     L = L[setdiff(1:size(L, 1), [1, rows[1]]), :]
    # end

    S.logicals = logs
end

"""
    changesigns!(S::AbstractStabilizerCode, charvec::Vector{nmod})

Set the character vector of `S` to `charvec` and update the signs.
"""
function changesigns!(S::AbstractStabilizerCode, charvec::Vector{nmod})
    R = base_ring(charactervector(S))
    length(charvec) == 2 * length(S) || error("Characteristic vector is of improper length for the code.")
    for s in charvec
        modulus(s) == modulus(R) || error("Phases are not in the correct ring.")
    end
    S.signs = _getsigns(stabilizers(S), charvec)
    S.charvec = charvec
end

function show(io::IO, S::AbstractStabilizerCode)
    if get(io, :compact, false)
        if typeof(S) <: CSSCode
            if typeof(dimension(S)) <: Integer
                if ismissing(minimumdistance(S))
                    println(io, "[[$(length(S)), $(dimension(S))]]_$(order(field(S))) CSS code.")
                else
                    println(io, "[[$(length(S)), $(dimension(S)), $(minimumdistance(S))]]_$(order(field(S))) CSS code.")
                end
            else # don't think this can ever be reached
                if ismissing(minimumdistance(S))
                    println(io, "(($(length(S)), $(dimension(S))))_$(order(field(S))) CSS code.")
                else
                    println(io, "(($(length(S)), $(dimension(S)), $(minimumdistance(S))))_$(order(field(S))) CSS code.")
                end
            end
        else
            if typeof(dimension(S)) <: Integer
                if ismissing(minimumdistance(S))
                    println(io, "[[$(length(S)), $(dimension(S))]]_$(order(field(S))) stabilizer code.")
                else
                    println(io, "[[$(length(S)), $(dimension(S)), $(minimumdistance(S))]]_$(order(field(S))) stabilizer code.")
                end
            else
                if ismissing(minimumdistance(S))
                    println(io, "(($(length(S)), $(dimension(S))))_$(order(field(S))) stabilizer code.")
                else
                    println(io, "(($(length(S)), $(dimension(S)), $(minimumdistance(S))))_$(order(field(S))) stabilizer code.")
                end
            end
        end
    else
        if typeof(S) <: CSSCode
            if typeof(dimension(S)) <: Integer
                if ismissing(minimumdistance(S))
                    println(io, "[[$(length(S)), $(dimension(S))]]_$(order(field(S))) CSS code.")
                else
                    println(io, "[[$(length(S)), $(dimension(S)), $(minimumdistance(S))]]_$(order(field(S))) CSS code.")
                end
            else # don't think this can ever be reached
                if ismissing(minimumdistance(S))
                    println(io, "(($(length(S)), $(dimension(S))))_$(order(field(S))) CSS code.")
                else
                    println(io, "(($(length(S)), $(dimension(S)), $(minimumdistance(S))))_$(order(field(S))) CSS code.")
                end
            end
            if isovercomplete(S)
                println(io, "X-stabilizer matrix (overcomplete): $(numXstabs(S)) × $(length(S))")
            else
                println(io, "X-stabilizer matrix: $(numXstabs(S)) × $(length(S))")
            end
            for i in 1:numXstabs(S)
                print(io, "\t chi($(Xsigns[i])) ")
                for j in 1:length(S)
                    if j != length(S)
                        print(io, "$(S.Xstabs[i, j]) ")
                    elseif j == length(S) && i != length(S) - dimension(S)
                        println(io, "$(S.Xstabs[i, j])")
                    else
                        print(io, "$(S.Xstabs[i, j])")
                    end
                end
            end
            if isovercomplete(S)
                println(io, "Z-stabilizer matrix (overcomplete): $(numZstabs(S)) × $(length(S))")
            else
                println(io, "Z-stabilizer matrix: $(numZstabs(S)) × $(length(S))")
            end
            for i in 1:numZstabs(S)
                print(io, "\t chi($(Zsigns[i])) ")
                for j in 1:length(S)
                    if j != length(S)
                        print(io, "$(S.Zstabs[i, j]) ")
                    elseif j == length(S) && i != length(S) - dimension(S)
                        println(io, "$(S.Zstabs[i, j])")
                    else
                        print(io, "$(S.Zstabs[i, j])")
                    end
                end
            end
        else
            if typeof(dimension(S)) <: Integer
                if ismissing(minimumdistance(S))
                    println(io, "[[$(length(S)), $(dimension(S))]]_$(order(field(S))) stabilizer code.")
                else
                    println(io, "[[$(length(S)), $(dimension(S)), $(minimumdistance(S))]]_$(order(field(S))) stabilizer code.")
                end
            else
                if ismissing(minimumdistance(S))
                    println(io, "(($(length(S)), $(dimension(S))))_$(order(field(S))) stabilizer code.")
                else
                    println(io, "(($(length(S)), $(dimension(S)), $(minimumdistance(S))))_$(order(field(S))) stabilizer code.")
                end
            end
            if isovercomplete(S)
                println(io, "Stabilizer matrix (overcomplete): $(nrows(S.stabs)) × $(length(S))")
            else
                println(io, "Stabilizer matrix: $(nrows(S.stabs)) × $(length(S))")
            end
            for i in 1:nrows(S.stabs)
                print(io, "\t chi($(signs[i])) ")
                for j in 1:length(S)
                    if j != length(S)
                        print(io, "$(S.stabs[i, j]) ")
                    elseif j == length(S) && i != length(S) - dimension(S)
                        println(io, "$(S.stabs[i, j])")
                    else
                        print(io, "$(S.stabs[i, j])")
                    end
                end
            end
            if !ismissing(S.sCWEstabs)
                println(io, "Signed complete weight enumerator for the stabilizer:")
                println(io, S.sCWEstabs)
            end
            if !ismissing(S.sCWEdual)
                println(io, "Signed complete weight enumerator for the normalizer:")
                println(io, S.sCWEdual)
            end
        end
    end
end

"""
    Xsyndrome(S::CSSCode, v::fq_nmod_mat)

Return the syndrome of the vector `v` with respect to the `X` stabilizers of the
CSS code.
"""
function Xsyndrome(S::CSSCode, v::fq_nmod_mat)
    n = length(S)
    if length(v) == 2 * n
        v = v[n + 1:end]
    end
    !(size(v) != (n, 1) && size(v) != (1, n)) ||
        error("Vector to be tested is of incorrect dimension; expected length $n, received: $(size(v)).")
    base_ring(v) == field(S) || error("Vector must have the same base ring as the stabilizers.")

    nrows(v) != 1 || return Xstabilizers(S) * transpose(v)
    return Xstabilizers(S) * v
end

"""
    Zsyndrome(S::CSSCode, v::fq_nmod_mat)

Return the syndrome of the vector `v` with respect to the `Z` stabilizers of the
CSS code.
"""
function Zsyndrome(S::CSSCode, v::fq_nmod_mat)
    n = length(S)
    if length(v) == 2 * n
        v = v[1:n]
    end
    !(size(v) != (n, 1) && size(v) != (1, n)) ||
        error("Vector to be tested is of incorrect dimension; expected length $n, received: $(size(v)).")
    base_ring(v) == field(S) || error("Vector must have the same base ring as the stabilizers.")

    nrows(v) != 1 || return Zstabilizers(S) * transpose(v)
    return Zstabilizers(S) * v
end

"""
    syndrome(S::AbstractStabilizerCode, v::fq_nmod_mat)

Return the syndrome of the vector `v` with respect to the stabilizers of `S`.
"""
function syndrome(S::AbstractStabilizerCode, v::fq_nmod_mat)
    n = length(S)
    !(size(v) != (2 * n, 1) && size(v) != (1, 2 * n)) ||
        error("Vector to be tested is of incorrect dimension; expected length $(2 * n), received: $(size(v)).")
    # base_ring(v) == field(S) || error("Vector must have the same base ring as the stabilizers.")

    nrows(v) != 1 || return symplecticstabilizers(S) * transpose(v)
    return symplecticstabilizers(S) * v
end

# function Singletonbound()
#     n - k ≧ 2(d - 1)
# end

# will have to get signs on them separately
# could make this considerably faster by using combinatorics package and skipping 0's
# using Combinatorics
# iter = combinations(1:size(stabs, 2))
# but this only in base 2
"""
    allstabilizers(Q::AbstractStabilizerCode, onlyprint::Bool=false)
    elements(Q::AbstractStabilizerCode, onlyprint::Bool=false)

Return the elements of `S`.

If `onlyprint` is `true`, the elements are only printed to the console and not returned.
"""
function allstabilizers(Q::AbstractStabilizerCode, onlyprint::Bool=false)
    E = quadraticfield(Q)
    all = Vector{fq_nmod_mat}()
    stabs = stabilizers(Q)
    for iter in Base.Iterators.product([0:(Int64(characteristic(field(Q))) - 1) for _ in 1:nrows(stabs)]...)
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
elements(Q::AbstractStabilizerCode, onlyprint::Bool=false) = allstabilizers(Q, onlyprint)

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
