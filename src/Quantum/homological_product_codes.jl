# Copyright (c) 2022 - 2026 Eric Sabo, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

# has this been extended to subsystem codes?
"""
    homological_product(S1::AbstractStabilizerCode, S2::AbstractStabilizerCode, U::CTMatrixTypes = identity_matrix(S1.F, S1.n), V::CTMatrixTypes = identity_matrix(S2.F, S2.n))
    ⊠(S1::AbstractStabilizerCode, S2::AbstractStabilizerCode) = homological_product(S1, S2)

Return the single-sector homological product code of `S1` and `S2`.

# Note
- This is the single-sector homological product. Use ⊗ for the more general product.
"""
function homological_product(S1::AbstractStabilizerCode, S2::AbstractStabilizerCode,
    U::CTMatrixTypes = identity_matrix(field(S1), num_X_stabs(S1)),
    V::CTMatrixTypes = identity_matrix(field(S2), num_X_stabs(S2));
    char_vec::Union{Vector{zzModRingElem}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm)

    return homological_product(CSSTrait(typeof(S1)), CSSTrait(typeof(S2)),
        S1, S2, U, V; char_vec=char_vec, logs_alg=logs_alg)
end

function homological_product(::IsCSS, ::IsCSS, S1::AbstractStabilizerCode,
    S2::AbstractStabilizerCode, U::CTMatrixTypes, V::CTMatrixTypes;
    char_vec::Union{Vector{zzModRingElem}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm)

    num_stabs1 = num_X_stabs(S1)
    num_stabs1 == num_Z_stabs(S1) || throw(ArgumentError("The first code didn't have the same number of X and Z stabilizers"))
    num_stabs2 = num_X_stabs(S2)
    num_stabs2 == num_Z_stabs(S2) || throw(ArgumentError("The second code didn't have the same number of X and Z stabilizers"))
    nrows(U) == ncols(U) == num_stabs1 ||
        throw(ArgumentError("U must be square with one row per X stabilizer of S1."))
    nrows(V) == ncols(V) == num_stabs2 ||
        throw(ArgumentError("V must be square with one row per X stabilizer of S2."))
    
    F = S1.F
    F == S2.F == base_ring(U) == base_ring(V) || throw(ArgumentError("S1, S2, U, and V should all have the same base ring"))

    # Physical qubits scale multiplicatively
    n_new = S1.n * S2.n
    
    cache = _family_cache(n_new; logs_alg=logs_alg)
    result = HomologicalProductCode(F, S1, S2, U, V, n_new, 0,
        _family_char_vec(char_vec, F, n_new), cache)
    result.k = n_new - rank(X_stabilizers(result)) - rank(Z_stabilizers(result))
    return _seed_quantum_singleton_bound!(result)
end

homological_product(::IsNotCSS, ::IsNotCSS, S1::AbstractStabilizerCode, S2::AbstractStabilizerCode, U, V; kwargs...) = throw(ArgumentError("This is only defined for CSS codes"))
homological_product(::IsNotCSS, ::IsCSS, S1::AbstractStabilizerCode, S2::AbstractStabilizerCode, U, V; kwargs...) = throw(ArgumentError("This is only defined for CSS codes"))
homological_product(::IsCSS, ::IsNotCSS, S1::AbstractStabilizerCode, S2::AbstractStabilizerCode, U, V; kwargs...) = throw(ArgumentError("This is only defined for CSS codes"))

@doc (@doc homological_product)
⊠(S1::AbstractStabilizerCode, S2::AbstractStabilizerCode) = homological_product(S1, S2)

function X_stabilizers(S::HomologicalProductCode)
    haskey(S.cache, :H_X) && return S.cache[:H_X]
    
    S1, S2 = S.S1, S.S2
    U, V = S.U, S.V
    F = S.F
    
    num_stabs1 = num_X_stabs(S1)
    num_stabs2 = num_X_stabs(S2)
    
    X1, Z1 = X_stabilizers(S1), Z_stabilizers(S1)
    δ1 = zero_matrix(F, S1.n, S1.n)
    for i in 1:num_stabs1
        for j in 1:num_stabs1
            δ1 += U[i, j] * transpose(Z1[i:i, :]) * X1[j:j, :]
        end
    end

    X2, Z2 = X_stabilizers(S2), Z_stabilizers(S2)
    δ2 = zero_matrix(F, S2.n, S2.n)
    for i in 1:num_stabs2
        for j in 1:num_stabs2
            δ2 += V[i, j] * transpose(Z2[i:i, :]) * X2[j:j, :]
        end
    end

    i1 = identity_matrix(F, S1.n)
    i2 = identity_matrix(F, S2.n)
    
    # ∂ = δ1 ⊗ i2 - i1 ⊗ δ2
    ∂ = kronecker_product(δ1, i2) - kronecker_product(i1, δ2)
    
    S.cache[:H_X] = ∂
    return ∂
end

function Z_stabilizers(S::HomologicalProductCode)
    haskey(S.cache, :H_Z) && return S.cache[:H_Z]
    
    # For this construction, H_Z is just transpose(H_X)
    H_X = X_stabilizers(S)
    H_Z = transpose(H_X)
    
    S.cache[:H_Z] = H_Z
    return H_Z
end

function stabilizers(S::HomologicalProductCode)
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
