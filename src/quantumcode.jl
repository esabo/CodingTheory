# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.


mutable struct CSSCode <: AbstractCSSCode
    F::FqNmodFiniteField # base field (symplectic)
    E::FqNmodFiniteField # additive field
    n::Integer
    k::Integer
    d::Union{Integer, Missing}
    dx::Union{Integer, Missing}
    dz::Union{Integer, Missing}
    stabs::fq_nmod_mat
    Xstabs::fq_nmod_mat
    Zstabs::fq_nmod_mat
    Xorigcode::Union{LinearCode, Missing}
    ZorigCode::Union{LinearCode, Missing}
    signs::Vector{Int64} # make bools, uint8's?
    Xsigns::Vector{Int64} # make bools, uint8's?
    Zsigns::Vector{Int64} # make bools, uint8's?
    dualgens::fq_nmod_mat
    logspace::Union{fq_nmod_mat, Missing}
    logicals::Union{Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}, Missing}
    charvec::Vector{Int64} # make bools, uint8's?
    dualweightenum::Union{WeightEnumerator, Missing}
    logsweightenum::Union{WeightEnumerator, Missing}
    Pauliweightenum::Union{WeightEnumerator, Missing}
end

mutable struct QuantumCode <: AbstractStabilizerCode
    F::FqNmodFiniteField # base field (symplectic)
    E::FqNmodFiniteField # additive field
    n::Integer
    k::Union{Integer, Rational{Int64}}
    d::Union{Integer, Missing}
    stabs::fq_nmod_mat
    dualgens::fq_nmod_mat
    logspace::Union{fq_nmod_mat, Missing}
    logicals::Union{Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}, Missing}
    charvec::Vector{Int64} # make bools, uint8's?
    signs::Vector{Int64}
    dualweightenum::Union{WeightEnumerator, Missing}
    logsweightenum::Union{WeightEnumerator, Missing}
    Pauliweightenum::Union{WeightEnumerator, Missing}
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

numXstabs(S::CSSCode) = size(S.Xstabs, 1)
numZstabs(S::CSSCode) = size(S.Zstabs, 1)

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
    d > 0 && d <= length(S) || error("The minimum distance of a code must be ≥ 1; received: d = $d.")
    S.d = d
end

"""
    relativedistance(S::AbstractStabilizerCode)

Return the relative minimum distance, `δ = d / n` of the code if `d` is known,
otherwise errors.
"""
function relativedistance(S::AbstractStabilizerCode)
    !ismissing(S.d) || error("Unknown minimum distance for this code.")
    return S.d / S.n
end

function _getsigns(A::fq_nmod_mat, charvec::Vector{Int64})
    if length(charvec) == 2 * size(A, 2)
        A = quadratictosymplectic(A)
    end
    length(charvec) == size(A, 2) || error("Input to _getsigns is expected to be in symplectic form and of the same length as the characteristic vector.")
    if charvec == ones(size(A, 2))
        return ones(Int64, div(size(A, 2), 2))
    end

    signs = Vector{Int64}()
    for r in 1:size(A, 1)
        parity = 1
        for c = 1:size(A, 2)
            if !iszero(A[r, c]) && charvec[c] == -1
                parity *= -1
            end
        end
        append!(signs, parity)
    end
    return signs
end

function _splitsymplecticstabilizers(S::fq_nmod_mat, signs::Vector{Int64})
    Xstabs = []
    Xsigns = Vector{Int64}()
    Zstabs = []
    Zsigns = Vector{Int64}()
    mixedstabs = []
    mixedsigns = Vector{Int64}()

    half = div(size(S, 2), 2)
    for r in 1:size(S, 1)
        # use views?
        s = S[r, :]
        if iszero(s)
            continue
        else
            sx = iszero(s[1, 1:half])
            sz = iszero(s[1, half + 1:end])
            if (sx && !sz)
                push!(Zstabs, s)
                append!(Zsigns, signs[r])
            elseif !sx && sz
                push!(Xstabs, s)
                append!(Xsigns, signs[r])
            elseif !sx && !sz
                push!(mixedstabs, s)
                append!(mixedsigns, signs[r])
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
"""
function splitstabilizers(S::AbstractStabilizerCode)
    return _splitsymplecticstabilizers(symplecticstabilizers(S), signs(S))
end

# TODO: remove this function and make this a parameter in the structs
function _isCSSsymplectic(S::fq_nmod_mat, signs::Vector{Int64}=[], trim::Bool=true)
    Xstabs, Xsigns, Zstabs, Zsigns, mixedstabs, mixedsigns = _splitsymplecticstabilizers(S, signs)
    if isempty(mixedstabs)
        if trim
            half = div(size(Xstabs, 2), 2)
            return true, Xstabs[:, 1:half], Xsigns, Zstabs[:, half + 1:end], Zsigns
        else
            return true, Xstabs, Xsigns, Zstabs, Zsigns
        end
    else
        if trim
            half = div(size(Xstabs, 2), 2)
            return false, Xstabs[:, 1:half], Xsigns, Zstabs[:, half + 1:end], Zsigns, mixedstabs, mixedsigns
        else
            return false, Xstabs, Xsigns, Zstabs, Zsigns, mixedstabs, mixedsigns
        end
    end
end

"""
    splitstabilizers(S::AbstractStabilizerCode)

Return the set of `X`-only stabilizers and their signs, the set of `Z`-only
stabilizers and their signs, and the remaining stabilizers and their signs.
"""
function isCSS(S::AbstractStabilizerCode)
    return
end

"""
    CSSCode(C1::AbstractLinearCode, C2::AbstractLinearCode, charvec::Union{Vector{Int64}, Vector{Any}}=[])

Return a CSS code using the CSS construction on two linear codes `C1 = [n, k1, d1]`
and `C2 = [n, k2, d2]` with `C2 ⊆ C1`.

The resulting code has dimension `k = k1 - k2` and minimum distance
`d >= min(d1, d2^⟂)`. The `X` stabilizers are given by the parity-check matrix
of `C2^⟂`, `H(C2^⟂)`, and the `Z` stabilizers by `H(C1)`.

# Arguments
* `C1`: a linear code
* `C2`: a subcode of `C1`
* `charvec`: a length `2n` vector with elements in {+/-1} whose first `n` elements
  specify the `X` phases and second `n` the `Z` phases; a missing or empty argument
  will be set to the all-ones vector

# Notes
* Stabilizer signs are automatically computed given the character vector.
"""
function CSSCode(C1::AbstractLinearCode, C2::AbstractLinearCode,
    charvec::Union{Vector{Int64}, Vector{Any}}=[])

    length(C1) ==  length(C2) || error("Both codes must have the same length in the CSS construction.")
    field(C1) == field(C2) || error("Both codes must be over the same base field in the CSS construction.")
    C2 ⊆ C1 || error("The second argument must be a subset of the first in the CSS construction.")
    if !isempty(charvec)
        for s in charvec
            (s == 1 || s == -1) || error("X-stabilizer signs must be +/- 1; received: $s")
        end
    end

    # C2 ⊆ C1
    # k = k1 - k2
    # d >= minimum(d1, d2^⟂)
    # X - H(C2^⟂), Z - H(C1)
    D2 = dual(C2)
    S = directsum(paritycheckmatrix(D2), paritycheckmatrix(C1))
    Sq2 = symplectictoquadratic(S)
    if isempty(charvec)
        charvec = ones(Int64, 2 * length(C1))
        signs = ones(Int64, size(S, 1))
        Xsigns = ones(Int64, size(paritycheckmatrix(D2), 1))
        Zsigns = ones(Int64, size(paritycheckmatrix(C1), 1))
    else
        signs = _getsigns(S, charvec)
        Xsigns = signs[1:size(paritycheckmatrix(D2), 1), :]
        Zsigns = signs[size(paritycheckmatrix(D2), 1) + 1:end, :]
    end
    dualgens = directsum(generatormatrix(C1), generatormatrix(D2))

    return CSSCode(field(C1), base_ring(Sq2), length(C1), dimension(C1) - dimension(C2),
        missing, minimumdistance(D2), minimumdistance(C1), Sq2, paritycheckmatrix(D2),
        paritycheckmatrix(C1), C2, C1, signs, Xsigns, Zsigns, dualgens, missing,
        missing, charvec, missing, missing, missing)
end

"""
    CSSCode(C::AbstractLinearCode, charvec::Union{Vector{Int64}, Vector{Any}}=[])

Return a CSS code using the CSS construction on a self-orthogonal linear code
`C`, i.e., `C ⊆ C^⟂`.

Setting `C1 = C^⟂` and `C2 = C`, the resulting code has dimension `k = k1 - k2`
and minimum distance `d >= min(d1, d2^⟂)`. The `X` stabilizers are given by the
parity-check matrix of `C2^⟂`, `H(C2^⟂)`, and the `Z` stabilizers by `H(C1)`.

# Arguments
* `C`: a self-orthogonal linear code
* `charvec`: a length `2n` vector with elements in {+/-1} whose first `n` elements
  specify the `X` phases and second `n` the `Z` phases; a missing or empty argument
  will be set to the all-ones vector

# Notes
* Stabilizer signs are automatically computed given the character vector.
"""
function CSSCode(C::LinearCode, charvec::Union{Vector{Int64}, Vector{Any}}=[])
    # this should have Xstabs = Zstabs
    D = dual(C)
    C ⊆ D || error("The single code CSS construction requires C ⊆ C^⟂.")
    if !isempty(charvec)
        for s in charvec
            (s == 1 || s == -1) || error("X-stabilizer signs must be +/- 1; received: $s")
        end
    end

    # C2 ⊆ C1
    # k = k1 - k2
    # d >= minimum(d1, d2^⟂)
    # X - H(C2^⟂), Z - H(C1)
    S = directsum(paritycheckmatrix(D), paritycheckmatrix(D))
    Sq2 = symplectictoquadratic(S)
    k = dimension(D) - dimension(C)
    if isempty(charvec)
        charvec = ones(Int64, 2 * length(C))
        signs = ones(Int64, size(S, 1))
        Xsigns = ones(Int64, size(paritycheckmatrix(D), 1))
        Zsigns = ones(Int64, size(paritycheckmatrix(D), 1))
    else
        signs = _getsigns(S, charvec)
        Xsigns = signs[1:size(paritycheckmatrix(D), 1), :]
        Zsigns = signs[size(paritycheckmatrix(D), 1) + 1:end, :]
    end
    dualgens = directsum(generatormatrix(D), generatormatrix(D))

    return CSSCode(field(C), base_ring(Sq2), length(C), k, missing, minimumdistance(D),
        minimumdistance(D), Sq2, paritycheckmatrix(D), paritycheckmatrix(D), C, D,
        signs, Xsigns, Zsigns, dualgens, missing, missing, charvec, missing,
        missing, missing)
end

"""
    CSSCode(Xmatrix::fq_nmod_mat, Zmatrix::fq_nmod_mat, charvec::Union{Vector{Int64}, Vector{Any}}=[])

Return a CSS code using the matrix `Xmatrix` as the `X` stabilizer matrix
and `Zmatrix` as the `Z` stabilizer matrix.

# Arguments
* `Xmatrix`: a matrix over a finite field of type `FqNmodFiniteField`
* `Zmatrix`: a matrix over a finite field of type `FqNmodFiniteField`
* `charvec`: a length `2n` vector with elements in {+/-1} whose first `n` elements
  specify the `X` phases and second `n` the `Z` phases; a missing or empty argument
  will be set to the all-ones vector

# Notes
* Stabilizer signs are automatically computed given the character vector.
* The orthogonality of the `X` and `Z` stabilizers are automatically checked and
  will error upon failure.
"""
function CSSCode(Xmatrix::fq_nmod_mat, Zmatrix::fq_nmod_mat,
    charvec::Union{Vector{Int64}, Vector{Any}}=[])

    # this should have Zstabs * Xstabs' = 0
    # set dual containing if Xstabs == Zstabs, else not
    size(Xmatrix, 2) ==  size(Zmatrix, 2) || error("Both matrices must have the same length in the CSS construction.")
    base_ring(Xmatrix) == base_ring(Zmatrix) || error("Both matrices must be over the same base field in the CSS construction.")
    # TODO can remove rank check now that we have full character vector
    rank(Xmatrix) == size(Xmatrix, 1) || error("Provided X-stabilizer matrix is not full rank.")
    rank(Zmatrix) == size(Zmatrix, 1) || error("Provided Z-stabilizer matrix is not full rank.")
    if !isempty(charvec)
        for s in charvec
            (s == 1 || s == -1) || error("X-stabilizer signs must be +/- 1; received: $s")
        end
    end

    S = directsum(Xmatrix, Zmatrix)
    aresymplecticorthogonal(S, S) || error("The given matrices are not symplectic orthogonal.")
    Sq2 = symplectictoquadratic(S)
    if isempty(charvec)
        charvec = ones(Int64, size(S, 2))
        signs = ones(Int64, size(S, 1))
        Xsigns = ones(Int64, size(Xmatrix, 1))
        Zsigns = ones(Int64, size(Zmatrix, 1))
    else
        signs = _getsigns(S, charvec)
        Xsigns = signs[1:size(Xmatrix, 1), :]
        Zsigns = signs[size(Xmatrix, 1) + 1:end, :]
    end

    G = hcat(S[:, size(Xmatrix, 2) + 1:end], -S[:, 1:size(Xmatrix, 2)])
    _, H = right_kernel(G)
    size(H', 1) == 2 * size(Sq2, 2) - size(Sq2, 1) || error("Normalizer matrix is not size n + k.")
    # note the H here is transpose of the standard definition
    !(!iszero(G * H) || !iszero(H' * G')) || error("Normalizer matrix is not transpose, symplectic orthogonal.")
    dualgens = symplectictoquadratic(H')

    return CSSCode(base_ring(Xmatrix), base_ring(Sq2), size(Xmatrix, 2),
        size(Xmatrix, 2) - size(Xmatrix, 1) - size(Zmatrix, 1),
        missing, missing, missing, Sq2, Xmatrix, Zmatrix, missing, missing, signs,
        Xsigns, Zsigns, dualgens, missing, missing, charvec, missing, missing,
        missing)
end

"""
    CSSCode(SPauli::Vector{T}, charvec::Union{Vector{Int64}, Vector{Any}}=[]) where T <: Union{String, Vector{Char}}

Return a CSS code using the vector of Pauli strings `SPauli` as stabilizers.

# Arguments
* `SPauli`: a vector of Strings or Char vectors containing the letters {I, X, Y, Z}
* `charvec`: a length `2n` vector with elements in {+/-1} whose first `n` elements
  specify the `X` phases and second `n` the `Z` phases; a missing or empty argument
  will be set to the all-ones vector

# Notes
* Stabilizer signs are automatically computed given the character vector.
* The orthogonality of the `X` and `Z` stabilizers are automatically checked and
  will error upon failure.
* Any +/- 1 characters in front of each stabilizer are stripped. No check is done
  to make sure these signs agree with the ones computed using the character vector.
* Will error when the provided strings are not CSS.
"""
function CSSCode(SPauli::Vector{T}, charvec::Union{Vector{Int64}, Vector{Any}}=[])
    where T <: Union{String, Vector{Char}}

    SPaulistripped, charvec = _processstrings(SPauli, charvec)
    S = _Paulistringtosymplectic(SPaulistripped)
    aresymplecticorthogonal(S, S) || error("The given stabilizers are not symplectic orthogonal.")
    if isempty(charvec)
        charvec = ones(Int64, size(S, 2))
    end
    signs = _getsigns(S, charvec)
    Sq2 = symplectictoquadratic(S)

    del = Vector{Int64}()
    for r in 1:size(Sq2, 1)
        if iszero(Sq2[r, :])
            append!(del, r)
        end
    end

    # not necessary but faster to skip
    if !isempty(del)
        Sq2 = Sq2[setdiff(1:size(Sq2, 1), del), :]
        S = S[setdiff(1:size(Sq2, 1), del), :]
        signs = signs[setdiff(1:size(Sq2, 1), del)]
    end

    # if rank(Sq2) != size(Sq2, 1)
    # rk = rank(G)
    # !iszero(rk) || error("Rank zero matrix passed into LinearCode constructor.")
    # if rk < size(G, 1)
    #     _, G = rref(G)
    # end

    # find dual
    G = hcat(S[:, size(Sq2, 2) + 1:end], -S[:, 1:size(Sq2, 2)])
    _, H = right_kernel(G)
    size(H', 1) == 2 * size(Sq2, 2) - size(Sq2, 1) || error("Normalizer matrix is not size n + k.")
    # note the H here is transpose of the standard definition
    !(!iszero(G * H) || !iszero(H' * G')) || error("Normalizer matrix is not transpose, symplectic orthogonal.")
    dualgens = symplectictoquadratic(H')

    args = _isCSSsymplectic(S, signs, true)
    if args[1]
        return CSSCode(base_ring(S), base_ring(Sq2), size(Sq2, 2), size(Sq2, 2) - size(Sq2, 1),
            missing, missing, missing, Sq2, args[2], args[4], missing, missing, signs,
            args[3], args[5], dualgens, missing, missing, charvec, missing, missing,
            missing)
    else
        error("Provided Pauli strings are not CSS.")
    end
end

# entanglement-assisted is not symplectic orthogonal
"""
    QuantumCode(SPauli::Vector{T}, charvec::Union{Vector{Int64}, Vector{Any}}=[]) where T <: Union{String, Vector{Char}}

Return a stabilizer code using the vector of Pauli strings `SPauli` as stabilizers.

# Arguments
* `SPauli`: a vector of Strings or Char vectors containing the letters {I, X, Y, Z}
* `charvec`: a length `2n` vector with elements in {+/-1} whose first `n` elements
  specify the `X` phases and second `n` the `Z` phases; a missing or empty argument
  will be set to the all-ones vector

# Notes
* Stabilizer signs are automatically computed given the character vector.
* The orthogonality of the stabilizers are automatically checked and will error
  upon failure.
* Any +/- 1 characters in front of each stabilizer are stripped. No check is done
  to make sure these signs agree with the ones computed using the character vector.
"""
function QuantumCode(SPauli::Vector{T}, charvec::Union{Vector{Int64}, Vector{Any}}=[])
    where T <: Union{String, Vector{Char}}

    SPaulistripped, charvec = _processstrings(SPauli, charvec)
    S = _Paulistringtosymplectic(SPaulistripped)
    aresymplecticorthogonal(S, S) || error("The given stabilizers are not symplectic orthogonal.")
    if isempty(charvec)
        charvec = ones(Int64, size(S, 2))
    end
    signs = _getsigns(S, charvec)
    Sq2 = symplectictoquadratic(S)

    # make del zero rows, cols helper functions
    # which will return a flag if this was ever done
    del = Vector{Int64}()
    for r in 1:size(Sq2, 1)
        if iszero(Sq2[r, :])
            append!(del, r)
        end
    end

    # not necessary but faster to skip
    if !isempty(del)
        Sq2 = Sq2[setdiff(1:size(Sq2, 1), del), :]
        S = S[setdiff(1:size(Sq2, 1), del), :]
        signs = signs[setdiff(1:size(Sq2, 1), del)]
    end

    # need to code an additive only rref but then apply _getsigns
    # if rank(Sq2) != size(Sq2, 1)
    # rk = rank(G)
    # !iszero(rk) || error("Rank zero matrix passed into LinearCode constructor.")
    # if rk < size(G, 1)
    #     _, G = rref(G)
    # end

    # find dual
    G = hcat(S[:, size(Sq2, 2) + 1:end], -S[:, 1:size(Sq2, 2)])
    _, H = right_kernel(G)
    size(H', 1) == 2 * size(Sq2, 2) - size(Sq2, 1) || error("Normalizer matrix is not size n + k.")
    # note the H here is transpose of the standard definition
    !(!iszero(G * H) || !iszero(H' * G')) || error("Normalizer matrix is not transpose, symplectic orthogonal.")
    dualgens = symplectictoquadratic(H')

    # q^n / p^k but rows is n - k
    dimcode = Int64(order(base_ring(S)))^size(Sq2, 2) // Int64(characteristic(base_ring(S)))^(size(Sq2, 1))
    if isinteger(dimcode)
        dimcode = Int64(log(Int64(characteristic(base_ring(S))), Int64(dimcode)))
    end

    args = _isCSSsymplectic(S, signs, true)
    if args[1]
        return CSSCode(base_ring(S), base_ring(Sq2), size(Sq2, 2), dimcode,
            missing, missing, missing, Sq2, args[2], args[4], missing, missing, signs,
            args[3], args[5], dualgens, missing, missing, charvec, missing, missing,
            missing)
    else
        return QuantumCode(base_ring(S), base_ring(Sq2), size(Sq2, 2), dimcode,
            missing, Sq2, dualgens, missing, missing, charvec, signs, missing,
            missing, missing)
    end
end

"""
    QuantumCode(Sq2::fq_nmod_mat, symp::Bool=false, charvec::Union{Vector{Int64}, Vector{Any}}=[])

Return a stabilizer code using the matrix `Sq2` as the stabilizer matrix.

The matrix `Sq2` is assumed to be an `n` column matrix over the quadratic extension.
If the optional parameter `symp` is set to `true`, `Sq2` is assumed to be a `2n`
column matrix over the base field.

# Arguments
* `Sq2`: a matrix over a finite field of type `FqNmodFiniteField`
* `symp`: a boolean
* `charvec`: a length `2n` vector with elements in {+/-1} whose first `n` elements
  specify the `X` phases and second `n` the `Z` phases; a missing or empty argument
  will be set to the all-ones vector

# Notes
* Stabilizer signs are automatically computed given the character vector.
* The orthogonality of the stabilizers are automatically checked and will error
  upon failure.
"""
function QuantumCode(Sq2::fq_nmod_mat, symp::Bool=false, charvec::Union{Vector{Int64}, Vector{Any}}=[]) where T <: Union{String, Vector{Char}}
    if symp
        iseven(size(Sq2, 2)) || error("Expected a symplectic input but the input matrix has an odd number of columns.")
        if !isempty(charvec)
            size(Sq2, 2) == length(charvec) || error("The characteristic value is of incorrect length.")
        end
        S = Sq2
        Sq2 = symplectictoquadratic(Sq2)
    else
        E = base_ring(Sq2)
        iseven(degree(E)) || error("The base ring of the given matrix is not a quadratic extension.")
        S = quadratictosymplectic(Sq2)
    end
    if isempty(charvec)
        charvec = ones(Int64, size(S, 2))
    end
    aresymplecticorthogonal(S, S) || error("The given stabilizers are not symplectic orthogonal.")
    signs = _getsigns(S, charvec)
    # make del zero rows, cols helper functions
    # which will return a flag if this was ever done
    del = Vector{Int64}()
    for r in 1:size(Sq2, 1)
        if iszero(Sq2[r, :])
            append!(del, r)
        end
    end

    # not necessary but faster to skip
    if !isempty(del)
        Sq2 = Sq2[setdiff(1:size(Sq2, 1), del), :]
        S = S[setdiff(1:size(Sq2, 1), del), :]
        signs = signs[setdiff(1:size(Sq2, 1), del)]
    end

    # need to code an additive only rref but then apply _getsigns
    # if rank(Sq2) != size(Sq2, 1)
    # rk = rank(G)
    # !iszero(rk) || error("Rank zero matrix passed into LinearCode constructor.")
    # if rk < size(G, 1)
    #     _, G = rref(G)
    # end

    # find dual
    G = hcat(S[:, size(Sq2, 2) + 1:end], -S[:, 1:size(Sq2, 2)])
    _, H = right_kernel(G)
    size(H', 1) == 2 * size(Sq2, 2) - size(Sq2, 1) || error("Normalizer matrix is not size n + k.")
    # note the H here is transpose of the standard definition
    !(!iszero(G * H) || !iszero(H' * G')) || error("Normalizer matrix is not transpose, symplectic orthogonal.")
    dualgens = symplectictoquadratic(H')

    try
        temp = div(Int64(order(base_ring(Sq2))), Int64(characteristic(base_ring(Sq2))))
        dimcode = temp^(size(Sq2, 2) - size(Sq2, 1))
    catch
        # over flows hard here - make rational{BigInt}?
        dimcode = Int64(order(base_ring(Sq2)))^size(Sq2, 2) // Int64(characteristic(base_ring(Sq2)))^(size(Sq2, 1))
    end
    if isinteger(dimcode)
        dimcode = Int64(log(Int64(characteristic(base_ring(S))), Int64(dimcode)))
    end

    args = _isCSSsymplectic(S, signs, true)
    if args[1]
        return CSSCode(base_ring(S), base_ring(Sq2), size(Sq2, 2), dimcode,
            missing, missing, missing, Sq2, args[2], args[4], missing, missing, signs,
            args[3], args[5], dualgens, missing, missing, charvec, missing, missing,
            missing)
    else
        return QuantumCode(base_ring(S), base_ring(Sq2), size(Sq2, 2), dimcode,
            missing, Sq2, dualgens, missing, missing, charvec, signs, missing,
            missing, missing)
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
    U, UtoV = sub(V, [V(G[i, :]) for i in 1:size(G, 1)])
    W, WtoV = sub(V, [V(Gdual[i, :]) for i in 1:size(Gdual, 1)])
    gensofUinW = [preimage(WtoV, UtoV(g)) for g in gens(U)]
    UinW, UinWtoW = sub(W, gensofUinW)
    Q, WtoQ = quo(W, UinW)
    C2modC1basis = [WtoV(x) for x in [preimage(WtoQ, g) for g in gens(Q)]]
    Fbasis = [[F(C2modC1basis[j][i]) for i in 1:AbstractAlgebra.dim(parent(C2modC1basis[1]))] for j in 1:length(C2modC1basis)]
    G2 = matrix(F, length(Fbasis), length(Fbasis[1]), vcat(Fbasis...))

    # returns linear object, make additive if necessary
    if size(G2, 1) == dimension(S)
        S.logspace = vcat(G2, gen(F) * G2)
        return S.logspace
    elseif size(G2, 1) == 2 * dimension(S)
        S.logspace = G2
        return S.logspace
    else
        error("Logical space produced of incorrect dimension; expected: ",
            2 * dimension(S), ", received: ", size(G2, 1))
    end
end

# # need a canonical phase to set Weyl pairs to
# # does not check which elements in Weyl pair are X or Z
# function logicaloperators(Q::AbstractStabilizerCode)
#     logspace = logicalspace(Q)
#     Weylpairs = Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}()
#     todo = [quadratictosymplectic(logspace)[i, :] for i in 1:size(logspace, 1)]
#
#     count = 1
#     while !isempty(todo)
#         curr1 = popfirst!(todo)
#         commutes = Vector{fq_nmod_mat}()
#         doesnt = Vector{Tuple{fq_nmod_mat, fq_nmod}}()
#         for i in 1:length(todo)
#             SIP = symplecticinnerproduct(curr1, todo[i])
#             if iszero(SIP)
#                 push!(commutes, todo[i])
#             else
#                 push!(doesnt, (todo[i], SIP))
#             end
#         end
#
#         if !isempty(doesnt)
#             curr2 = popfirst!(doesnt)
#             curr2, curr2SIP = curr2[1], curr2[2]
#         else
#             error("Logical basis element anti-commutes with no elements.")
#         end
#
#         todo = Vector{fq_nmod_mat}()
#         for i in 1:length(doesnt)
#             elm, SIP = doesnt[i]
#             elm .-= SIP * inv(curr2SIP) * curr2
#             push!(todo, elm)
#         end
#
#         push!(Weylpairs, [symplectictoquadratic(curr1), symplectictoquadratic(curr2)])
#         count < 50 ||error("Logical operator loop has gone through 50 iterations. Increase if more than this many logical qudits.")
#         count += 1
#     end
#
#     length(Weylpairs) == dimension(Q) || error("Completed symplectic Gram-Schmidt without generating all logicals.")
#     # isempty(todo) || error("Completed symplectic Gram-Schmidt without using all basis elements.")
#
#     Q.logicals = Weylpairs
#     return Weylpairs
# end

# can't check rank here because a self-dual code will have half rank since function is not additive
# need to think about what this means for an irrational dimension
# should check commutation relations
# TODO process this heavily and don't put this in logspace
"""
    setlogicals!(S::AbstractStabilizerCode, L::fq_nmod_mat)

Set the logical operators of `S` to `L`.
"""
function setlogicals!(S::AbstractStabilizerCode, L::fq_nmod_mat)
    size(L) == (2 * dimension(S), length(S)) || error("Provided matrix is of incorrect size for the logical space.")
    # rank(L) == size(L, 1) || error("Provided matrix is not of full rank.")
    # aresymplecticorthogonal(stabilizers(S), L) || error("Provided matrix does not commute with the code.")
    # !aresymplecticorthogonal(L, L) || error("Provided matrix should not be symplectic self-orthogonal.")
    S.logspace = L
end

# update for roots of unity
"""
    changesigns!(S::AbstractStabilizerCode, charvec::Vector{Int64})

Set the character vector of `S` to `charvec` and update the signs.
"""
function changesigns!(S::AbstractStabilizerCode, charvec::Vector{Int64})
    length(charvec) == 2 * length(S) || error("Characteristic vector is of improper length for the code.")
    for s in charvec
        (s == 1 || s == -1) || error("Qubit phases must be +/- 1; received: $s")
    end
    S.signs = _getsigns(stabilizers(S), charvec)
    S.charvec = charvec
end

# add for minimum distance check
function show(io::IO, S::AbstractStabilizerCode)
    if get(io, :compact, false)
        if typeof(S) <: CSSCode
            if typeof(dimension(S)) <: Integer
                println(io, "[[$(length(S)), $(dimension(S))]]_$(order(field(S))) CSS code.")
            else # don't think this can ever be reached
                println(io, "(($(length(S)), $(dimension(S))))_$(order(field(S))) CSS code.")
            end
        else
            if typeof(dimension(S)) <: Integer
                println(io, "[[$(length(S)), $(dimension(S))]]_$(order(field(S))) quantum (additive) code.")
            else
                println(io, "(($(length(S)), $(dimension(S))))_$(order(field(S))) quantum (additive) code.")
            end
        end
    else
        if typeof(S) <: CSSCode
            if typeof(dimension(S)) <: Integer
                println(io, "[[$(length(S)), $(dimension(S))]]_$(order(field(S))) CSS code.")
            else # don't think this can ever be reached
                println(io, "(($(length(S)), $(dimension(S))))_$(order(field(S))) CSS code.")
            end
            println(io, "X-stabilizer matrix: $(numXstabs(S)) × $(length(S))")
            for i in 1:numXstabs(S)
                print(io, "\t")
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
            println(io, "Z-stabilizer matrix: $(numZstabs(S)) × $(length(S))")
            for i in 1:numZstabs(S)
                print(io, "\t")
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
                println(io, "[[$(length(S)), $(dimension(S))]]_$(order(field(S))) quantum (additive) code.")
            else
                println(io, "(($(length(S)), $(dimension(S))))_$(order(field(S))) quantum (additive) code.")
            end
            println(io, "Stabilizer matrix: $(size(S.stabs, 1)) × $(length(S))")
            for i in 1:size(S.stabs, 1)
                print(io, "\t")
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

    size(v, 1) != 1 || return Xstabilizers(S) * v'
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

    size(v, 1) != 1 || return Zstabilizers(S) * v'
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

    size(v, 1) != 1 || return symplecticstabilizers(S) * v'
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
    allstabilizers(S::AbstractStabilizerCode)

Return the set of all stabilizers using brute-force.
"""
function allstabilizers(Q::AbstractStabilizerCode)
    E = quadraticfield(Q)
    all = Vector{fq_nmod_mat}()
    stabs = stabilizers(Q)
    for iter in Base.Iterators.product([0:(Int64(characteristic(field(Q))) - 1) for _ in 1:size(stabs, 1)]...)
        stab = E(iter[1]) * stabs[1, :]
        for r in 2:size(stabs, 1)
            if !iszero(iter[r])
                stab += E(iter[r]) * stabs[r, :]
            end
        end
        push!(all, stab)
    end
    return all
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
