# Copyright (c) 2022 - 2026 Eric Sabo, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

function LiftedProductCode(A::MatElem{T}, B::MatElem{T}; char_vec::Union{Vector{zzModRingElem}, Missing} =
    missing, logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

    R = parent(A[1, 1])
    R == parent(B[1, 1]) ||
        throw(ArgumentError("A and B must be defined over the same ring."))
    test_lift = lift(matrix(R, 1, 1, [A[1, 1]]))
    l = size(test_lift, 1)
    F = base_ring(test_lift)
    
    k1, n1 = size(A)
    k2, n2 = size(B)
    
    # Calculate physical qubits O(1)
    n_pre = k2 * n1 + n2 * k1
    n_new = l * n_pre
    
    cache = _family_cache(n_new; logs_alg=logs_alg)
    result = LiftedProductCode(F, A, B, n_new, 0,
        _family_char_vec(char_vec, F, n_new), cache)
    result.k = n_new - rank(X_stabilizers(result)) - rank(Z_stabilizers(result))
    return _seed_quantum_singleton_bound!(result)
end

function LiftedProductCode(A::MatElem{T}, b::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

    # A Generalized Hypergraph Product is exactly a Lifted Product where B is a 1x1 matrix.
    R = parent(A[1, 1])
    R == parent(b) || throw(ArgumentError("A and b must be defined over the same ring."))
    B_mat = matrix(R, 1, 1, [b])
    
    return LiftedProductCode(A, B_mat; char_vec = char_vec, logs_alg = logs_alg)
end
GeneralizedHypergraphProductCode(A, b; kwargs...) = LiftedProductCode(A, b; kwargs...)

function BiasTailoredLiftedProductCode(A::MatElem{T}, B::MatElem{T}; char_vec::Union{Vector{zzModRingElem},
    Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

    R = parent(A[1, 1])
    R == parent(B[1, 1]) ||
        throw(ArgumentError("A and B must be defined over the same ring."))
    test_lift = lift(matrix(R, 1, 1, [A[1, 1]]))
    l = size(test_lift, 1)
    F = base_ring(test_lift)
    
    k1, n1 = size(A)
    k2, n2 = size(B)
    
    # Calculate physical qubits O(1)
    n_pre = n1 * n2 + k1 * k2
    n_new = l * n_pre
    
    cache = _family_cache(n_new; logs_alg=logs_alg)
    result = BiasTailoredLiftedProductCode(F, A, B, n_new, 0,
        _family_char_vec(char_vec, F, n_new), cache)
    result.k = n_new - rank(stabilizers(result))
    return _seed_quantum_singleton_bound!(result)
end

function X_stabilizers(S::LiftedProductCode)
    haskey(S.cache, :H_X) && return S.cache[:H_X]
    
    A, B = S.A, S.B
    R = parent(A[1, 1])
    k1, n1 = size(A)
    k2, n2 = size(B)
    F = S.F
    
    Ek1 = identity_matrix(R, k1)
    Ek2 = identity_matrix(R, k2)
    
    if Int(characteristic(F)) == 2
        H_X_pre = hcat(kronecker_product(A, Ek2), kronecker_product(Ek1, B))
    else
        H_X_pre = hcat(kronecker_product(A, Ek2), kronecker_product(-Ek1, B))
    end
    
    H_X = lift(H_X_pre)
    S.cache[:H_X] = H_X
    return H_X
end

function Z_stabilizers(S::LiftedProductCode)
    haskey(S.cache, :H_Z) && return S.cache[:H_Z]
    
    A, B = S.A, S.B
    R = parent(A[1, 1])
    k1, n1 = size(A)
    k2, n2 = size(B)
    
    A_tr = _CT_adjoint(A)
    B_tr = _CT_adjoint(B)
    
    En1 = identity_matrix(R, n1)
    En2 = identity_matrix(R, n2)
    
    H_Z_pre = hcat(kronecker_product(En1, B_tr), kronecker_product(A_tr, En2))
    
    H_Z = lift(H_Z_pre)
    S.cache[:H_Z] = H_Z
    return H_Z
end

function stabilizers(S::LiftedProductCode)
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

function stabilizers(S::BiasTailoredLiftedProductCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    A, B = S.A, S.B
    R = parent(A[1, 1])
    k1, n1 = size(A)
    k2, n2 = size(B)
    
    A_tr = _CT_adjoint(A)
    B_tr = _CT_adjoint(B)
    
    Ek1 = identity_matrix(R, k1)
    Ek2 = identity_matrix(R, k2)
    En1 = identity_matrix(R, n1)
    En2 = identity_matrix(R, n2)
    
    A12 = kronecker_product(A_tr, Ek2) 
    A13 = kronecker_product(En1, B)    
    A21 = kronecker_product(A, En2)    
    A24 = kronecker_product(Ek1, B_tr) 
    
    # Safe zero blocks avoiding strict identity mirroring
    Z1 = zero_matrix(R, nrows(A12), ncols(A21))
    Z2 = zero_matrix(R, nrows(A13), ncols(A24))
    Z3 = zero_matrix(R, nrows(A21), ncols(A12))
    Z4 = zero_matrix(R, nrows(A24), ncols(A13))
    
    H_pre = vcat(hcat(Z1, A12, A13, Z2), 
                 hcat(A21, Z3, Z4, A24))
    
    stabs = lift(H_pre)
    S.cache[:stabilizers] = stabs
    return stabs
end
