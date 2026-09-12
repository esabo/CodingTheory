# Copyright (c) 2022 - 2026 Eric Sabo, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
$(TYPEDSIGNATURES)

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
function HyperBicycleCodeCSS(a::Vector{T}, b::Vector{T}, χ::Int;
    char_vec::Union{Vector{zzModRingElem}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm) where T <: CTMatrixTypes
    χ > 0 || throw(ArgumentError("Required χ > 0."))
    c = length(a)
    c > 0 || throw(ArgumentError("Input vectors must be nonempty."))
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

    # Calculate physical qubits O(1)
    n_new = c * (k2 * n1 + n2 * k1)
    
    H1, H2, HT1, HT2 =
        _build_hyperbicycle_blocks(a, b, χ, F, c, k1, n1, k2, n2)
    H_X = hcat(kronecker_product(identity_matrix(F, k2), H1),
        kronecker_product(H2, identity_matrix(F, k1)))
    H_Z = hcat(kronecker_product(HT2, identity_matrix(F, n1)),
        kronecker_product(identity_matrix(F, n2), HT1))
    k = n_new - rank(H_X) - rank(H_Z)
    cache = _family_cache(n_new; logs_alg=logs_alg)
    cache[:H_X], cache[:H_Z] = H_X, H_Z
    result = HyperBicycleCodeCSS(F, a, b, χ, n_new, k,
        _family_char_vec(char_vec, F, n_new), cache)
    return _seed_quantum_singleton_bound!(result)
end

"""
$(TYPEDSIGNATURES)

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
function HyperBicycleCode(a::Vector{T}, b::Vector{T}, χ::Int;
    char_vec::Union{Vector{zzModRingElem}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm) where T <: CTMatrixTypes
    χ > 0 || throw(ArgumentError("Required χ > 0."))
    c = length(a)
    c > 0 || throw(ArgumentError("Input vectors must be nonempty."))
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

    H1, H2, HT1, HT2 =
        _build_hyperbicycle_blocks(a, b, χ, F, c, k1, n1, k2, n2)
    (H1 == HT1 && H2 == HT2) ||
        throw(ArgumentError("H_i must equal H̃_i for i = 1, 2 for Non-CSS HyperBicycle."))
    stabs = hcat(kronecker_product(identity_matrix(F, k2), H1),
        kronecker_product(H2, identity_matrix(F, k1)))
    iseven(ncols(stabs)) ||
        throw(ArgumentError("The assembled symplectic matrix has odd width."))
    n_new = div(ncols(stabs), 2)
    k = n_new - rank(stabs)
    cache = _family_cache(n_new; logs_alg=logs_alg)
    cache[:stabilizers] = stabs
    result = HyperBicycleCode(F, a, b, χ, n_new, k,
        _family_char_vec(char_vec, F, n_new), cache)
    return _seed_quantum_singleton_bound!(result)
end

# Helper function that executes the internal Hyperbicycle sum math
function _build_hyperbicycle_blocks(a, b, χ, F, c, k1, n1, k2, n2)
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
        H1 += kronecker_product(Iχi, a[i])
        H2 += kronecker_product(b[i], Iχi)
        HT1 += kronecker_product(ITχi, transpose(a[i]))
        HT2 += kronecker_product(transpose(b[i]), ITχi)
    end
    return H1, H2, HT1, HT2
end

function X_stabilizers(S::HyperBicycleCodeCSS)
    haskey(S.cache, :H_X) && return S.cache[:H_X]
    
    a, b = S.a, S.b
    k1, n1 = size(a[1])
    k2, n2 = size(b[1])
    F = S.F
    
    H1, H2, _, _ = _build_hyperbicycle_blocks(a, b, S.χ, F, length(a), k1, n1, k2, n2)
    
    Ek1 = identity_matrix(F, k1)
    Ek2 = identity_matrix(F, k2)
    
    H_X = hcat(kronecker_product(Ek2, H1), kronecker_product(H2, Ek1))
    S.cache[:H_X] = H_X
    return H_X
end

function Z_stabilizers(S::HyperBicycleCodeCSS)
    haskey(S.cache, :H_Z) && return S.cache[:H_Z]
    
    a, b = S.a, S.b
    k1, n1 = size(a[1])
    k2, n2 = size(b[1])
    F = S.F
    
    _, _, HT1, HT2 = _build_hyperbicycle_blocks(a, b, S.χ, F, length(a), k1, n1, k2, n2)
    
    En1 = identity_matrix(F, n1)
    En2 = identity_matrix(F, n2)
    
    H_Z = hcat(kronecker_product(HT2, En1), kronecker_product(En2, HT1))
    S.cache[:H_Z] = H_Z
    return H_Z
end

function stabilizers(S::HyperBicycleCodeCSS)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    H_X = X_stabilizers(S)
    H_Z = Z_stabilizers(S)
    n = S.n
    F = S.F
    
    stabs = vcat(hcat(H_X, zero_matrix(F, size(H_X, 1), n)),
                 hcat(zero_matrix(F, size(H_Z, 1), n), H_Z))
                 
    S.cache[:stabilizers] = stabs
    return stabs
end

function stabilizers(S::HyperBicycleCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    a, b = S.a, S.b
    k1, n1 = size(a[1])
    k2, n2 = size(b[1])
    F = S.F
    
    H1, H2, HT1, HT2 = _build_hyperbicycle_blocks(a, b, S.χ, F, length(a), k1, n1, k2, n2)
    (H1 == HT1 && H2 == HT2) || throw(ArgumentError("H_i must equal H̃_i for i = 1, 2 for Non-CSS HyperBicycle."))
    
    Ek1 = identity_matrix(F, k1)
    Ek2 = identity_matrix(F, k2)

    stabs = hcat(kronecker_product(Ek2, H1), kronecker_product(H2, Ek1))
    S.cache[:stabilizers] = stabs
    return stabs
end
