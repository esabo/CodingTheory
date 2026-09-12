# Copyright (c) 2022 - 2026 Eric Sabo, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
$(TYPEDSIGNATURES)

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
function GeneralizedBicycleCode(A::T, B::T;
    char_vec::Union{Vector{zzModRingElem}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm) where T <: CTMatrixTypes
    F = base_ring(A)
    F == base_ring(B) || throw(ArgumentError("Arguments must be over the same base ring."))
    (iszero(A) || iszero(B)) && throw(ArgumentError("Arguments should not be zero."))
    
    m, n1 = size(A)
    m2, n2 = size(B)
    m == n1 == m2 == n2 || throw(ArgumentError("A and B must be square matrices of the same dimensions."))
    iszero(A * B - B * A) || throw(ArgumentError("Arguments must commute."))

    # Physical qubits n = 2m since H_X = [A | B]
    n_new = 2 * m
    
    H_X = hcat(A, B)
    H_Z = Int(characteristic(F)) == 2 ?
        hcat(transpose(B), transpose(A)) :
        hcat(transpose(B), -transpose(A))
    k = n_new - rank(H_X) - rank(H_Z)
    cache = _family_cache(n_new; logs_alg=logs_alg)
    cache[:H_X], cache[:H_Z] = H_X, H_Z
    result = GeneralizedBicycleCode(F, A, B, n_new, k,
        _family_char_vec(char_vec, F, n_new), cache)
    return _seed_quantum_singleton_bound!(result)
end

"""
$(TYPEDSIGNATURES)

Return the generealized bicycle code determined by `a` and `b`.

# Notes
- `l x l` circulant matrices are constructed using the coefficients of the polynomials
  `a` and `b` in `F_q[x]/(x^l - 1)` (`gcd(q, l) = 1`) as the first column
"""
function GeneralizedBicycleCode(a::T, b::T; kwargs...) where T <: ResElem
    parent(a) == parent(b) || throw(ArgumentError("Both objects must be defined over the same residue ring."))

    S = GeneralizedBicycleCode(residue_polynomial_to_circulant_matrix(a),
        residue_polynomial_to_circulant_matrix(b); kwargs...)
    S.cache[:polynomials] = (a, b)
    return S
end

"""
$(TYPEDSIGNATURES)

Return the generealized bicycle code determined by `a` and `b`.

# Notes
- `|G| x |G|` circulant matrices are constructed using the coefficients of the elements in the group algebra `FG` as` the first column
"""
function GeneralizedBicycleCode(a::T, b::T; kwargs...) where T <: CTGroupAlgebra
    parent(a) == parent(b) || throw(ArgumentError("Both objects must be defined over the same residue ring."))

    S = GeneralizedBicycleCode(group_algebra_element_to_circulant_matrix(a),
        group_algebra_element_to_circulant_matrix(b); kwargs...)
    S.cache[:polynomials] = (a, b)
    return S
end

function X_stabilizers(S::GeneralizedBicycleCode)
    haskey(S.cache, :H_X) && return S.cache[:H_X]
    
    # H_X = [A | B]
    H_X = hcat(S.A, S.B)
    S.cache[:H_X] = H_X
    return H_X
end

function Z_stabilizers(S::GeneralizedBicycleCode)
    haskey(S.cache, :H_Z) && return S.cache[:H_Z]
    
    F = S.F
    
    # H_Z = [B^T | ±A^T]
    if Int(order(F)) == 2
        H_Z = hcat(transpose(S.B), transpose(S.A))
    else
        H_Z = hcat(transpose(S.B), -transpose(S.A))
    end
    
    S.cache[:H_Z] = H_Z
    return H_Z
end

function stabilizers(S::GeneralizedBicycleCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    H_X = X_stabilizers(S)
    H_Z = Z_stabilizers(S)
    n = S.n
    F = S.F
    
    is_sparse = H_X isa SparseMatrixCSC || H_Z isa SparseMatrixCSC
    
    if is_sparse
        stabs = vcat(hcat(H_X, spzeros(F, size(H_X, 1), n)),
                     hcat(spzeros(F, size(H_Z, 1), n), H_Z))
    else
        stabs = vcat(hcat(H_X, zero_matrix(F, size(H_X, 1), n)),
                     hcat(zero_matrix(F, size(H_Z, 1), n), H_Z))
    end
                 
    S.cache[:stabilizers] = stabs
    return stabs
end

"""
$(TYPEDSIGNATURES)

Return a lazy Bicycle code given by the square matrix `A`. 
This is equivalent to a Generalized Bicycle Code where `B = A^T`.
"""
function BicycleCode(A::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm)
    
    m, n = size(A)
    m == n || throw(ArgumentError("Input matrix must be square."))
    
    F = base_ring(A)
    Int(order(F)) == 2 || throw(ArgumentError("Bicycle codes are strictly defined over GF(2)."))
    
    # A Bicycle code is exactly a Generalized Bicycle Code with B = A^T
    return GeneralizedBicycleCode(A, transpose(A); char_vec = char_vec, logs_alg = logs_alg)
end

"""
$(TYPEDSIGNATURES)

Return the lazy Bicycle code determined by the residue ring element `a`.
"""
function BicycleCode(a::ResElem; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm)
    
    A = residue_polynomial_to_circulant_matrix(a)
    
    S = BicycleCode(A; char_vec = char_vec, logs_alg = logs_alg)
    # Save the generating polynomial for researchers
    S.cache[:polynomials] = (a, )
    return S
end

"""
$(TYPEDSIGNATURES)

Return the lazy Bicycle code determined by the group algebra element `a`.
"""
function BicycleCode(a::CTGroupAlgebra; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm)
    
    A = group_algebra_element_to_circulant_matrix(a)
    
    S = BicycleCode(A; char_vec = char_vec, logs_alg = logs_alg)
    # Save the generating group algebra element for researchers
    S.cache[:polynomials] = (a, )
    return S
end
