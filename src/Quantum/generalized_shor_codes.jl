# Copyright (c) 2022 - 2026 Eric Sabo, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

# unable to yield quantum LDPC code families with non constant minimum distance
"""
$(TYPEDSIGNATURES)

Return the generalized Shor code of `C1` and `C2` with
``C1^\\perp \\subseteq C2``.
"""
function _generalized_shor_gauge_matrix(C1::AbstractLinearCode, C2::AbstractLinearCode)
    C1.F == C2.F || throw(ArgumentError("Codes must be over the same field."))
    F = C1.F
    H1 = parity_check_matrix(C1)
    H2 = parity_check_matrix(C2)
    
    m1, n1 = size(H1)
    m2, n2 = size(H2)
    n = C1.n * C2.n
    
    is_sparse = H1 isa SparseMatrixCSC || H2 isa SparseMatrixCSC
    F_one = F(1)
    
    if is_sparse
        H1_sp = H1 isa SparseMatrixCSC ? H1 : sparse(H1)
        H2_sp = H2 isa SparseMatrixCSC ? H2 : sparse(H2)
        
        I_n1 = spdiagm(0 => fill(F_one, n1))
        I_n2 = spdiagm(0 => fill(F_one, n2))
        
        G_X_mat = kron(H1_sp, I_n2)
        G_Z_mat = kron(I_n1, H2_sp)
        
        G_X_sym = hcat(G_X_mat, spzeros(F, size(G_X_mat, 1), n))
        G_Z_sym = hcat(spzeros(F, size(G_Z_mat, 1), n), G_Z_mat)
    else
        I_n1 = identity_matrix(F, n1)
        I_n2 = identity_matrix(F, n2)
        
        G_X_mat = kronecker_product(H1, I_n2)
        G_Z_mat = kronecker_product(I_n1, H2)
        
        G_X_sym = hcat(G_X_mat, zero_matrix(F, size(G_X_mat, 1), n))
        G_Z_sym = hcat(zero_matrix(F, size(G_Z_mat, 1), n), G_Z_mat)
    end
    
    return vcat(G_X_sym, G_Z_sym)
end

function GeneralizedShorCode(C1::AbstractLinearCode, C2::AbstractLinearCode;
    char_vec::Union{Vector{zzModRingElem}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm)
    dual(C1) ⊆ C2 ||
        throw(ArgumentError("The inputs must satisfy dual(C1) ⊆ C2."))
    G = _generalized_shor_gauge_matrix(C1, C2)
    temp = SubsystemCode(G; char_vec=char_vec)
    n_new = C1.n * C2.n
    expected_k = C1.k * C2.k
    temp.k == expected_k ||
        throw(ArgumentError("The input codes do not satisfy the generalized Shor dimension condition."))
    d = (!ismissing(C1.d) && !ismissing(C2.d)) ? min(C1.d, C2.d) : missing
    cache = _family_cache(n_new; d=d,
        l_bound=min(C1.l_bound, C2.l_bound),
        u_bound=min(C1.u_bound, C2.u_bound), logs_alg=logs_alg)
    cache[:stabs] = stabilizers(temp)
    cache[:X_stabs] = X_stabilizers(temp)
    cache[:Z_stabs] = Z_stabilizers(temp)
    cache[:gauge_ops] = gauges(temp)
    cache[:g_ops_mat] = gauges_matrix(temp)
    result = GeneralizedShorCode(C1.F, C1, C2, n_new, temp.k, temp.r,
        character_vector(temp), cache)
    return _seed_quantum_singleton_bound!(result)
end
"""
$(TYPEDSIGNATURES)

Return the Bacon--Casaccino subsystem code obtained from the classical codes
`C1` and `C2`, requiring ``C1^\\perp \\subseteq C2``. Its ``X`` gauge
generators replicate checks of `C1` across columns, and its ``Z`` gauge
generators replicate checks of `C2` across rows.
"""
BaconCasaccinoConstruction(C1::AbstractLinearCode, C2::AbstractLinearCode; kwargs...) =
    GeneralizedShorCode(C1, C2; kwargs...)

stabilizers(S::GeneralizedShorCode; standform::Bool=false) =
    invoke(stabilizers, Tuple{AbstractSubsystemCode}, S; standform=standform)
X_stabilizers(S::GeneralizedShorCode) = S.cache[:X_stabs]
Z_stabilizers(S::GeneralizedShorCode) = S.cache[:Z_stabs]
