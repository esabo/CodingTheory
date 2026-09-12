# Copyright (c) 2022 - 2026 Eric Sabo, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
$(TYPEDSIGNATURES)

Return a lazy `HypergraphProductCode`. Computes parameters `n`, `k`, and bounds 
instantly without generating the quantum parity check matrices.
"""
function HypergraphProductCode(C1::AbstractLinearCode, C2::AbstractLinearCode;
    char_vec::Union{Vector{zzModRingElem}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm)
    C1.F == C2.F || throw(ArgumentError("Codes must be over the same field."))
    
    # Parity check matrices for classical codes are generally already cached or cheap.
    # We only need them to grab the exact dimensions m1 and m2.
    H1 = parity_check_matrix(C1)
    H2 = parity_check_matrix(C2)
    m1, n1 = size(H1)
    m2, n2 = size(H2)
    
    # 1. Calculate Physical Qubits (n)
    n_new = n1 * n2 + m1 * m2
    
    # 2. Calculate Logical Qubits (k)
    C1T, C2T = transpose(C1), transpose(C2)
    k_new = C1.k * C2.k + C1T.k * C2T.k
    
    primal_sector = C1.k > 0 && C2.k > 0
    transposed_sector = C1T.k > 0 && C2T.k > 0
    dx_codes = AbstractLinearCode[]
    dz_codes = AbstractLinearCode[]
    primal_sector && (push!(dx_codes, C2); push!(dz_codes, C1))
    transposed_sector && (push!(dx_codes, C1T); push!(dz_codes, C2T))
    dx_lower = isempty(dx_codes) ? 1 : minimum(C.l_bound for C in dx_codes)
    dx_upper = isempty(dx_codes) ? n_new : minimum(C.u_bound for C in dx_codes)
    dz_lower = isempty(dz_codes) ? 1 : minimum(C.l_bound for C in dz_codes)
    dz_upper = isempty(dz_codes) ? n_new : minimum(C.u_bound for C in dz_codes)
    active_codes = (dx_codes..., dz_codes...)
    exact = map(C -> C.d, active_codes)
    d = !isempty(exact) && all(x -> !ismissing(x), exact) ? minimum(exact) : missing
    cache = _family_cache(n_new; d=d,
        l_bound=min(dx_lower, dz_lower), u_bound=min(dx_upper, dz_upper),
        logs_alg=logs_alg)
    cache[:l_bound_dx], cache[:u_bound_dx] = dx_lower, dx_upper
    cache[:l_bound_dz], cache[:u_bound_dz] = dz_lower, dz_upper

    result = HypergraphProductCode(C1.F, C1, C2, C1T, C2T, n_new, k_new,
        _family_char_vec(char_vec, C1.F, n_new), cache)
    return _seed_quantum_singleton_bound!(result)
end
HypergraphProductCode(C::AbstractLinearCode; kwargs...) =
    HypergraphProductCode(C, C; kwargs...)
HypergraphProductCode(H1::CTMatrixTypes, H2::CTMatrixTypes; kwargs...) =
    HypergraphProductCode(LinearCode(H1, true), LinearCode(H2, true); kwargs...)
HypergraphProductCode(H::CTMatrixTypes; kwargs...) =
    HypergraphProductCode(H, H; kwargs...)

function X_stabilizers(S::HypergraphProductCode)
    haskey(S.cache, :H_X) && return S.cache[:H_X]
    
    F = S.F
    H1 = parity_check_matrix(S.C1)
    H2 = parity_check_matrix(S.C2)
    m1, n1 = size(H1)
    m2, n2 = size(H2)
    
    is_sparse = H1 isa SparseMatrixCSC || H2 isa SparseMatrixCSC
    F_one = F(1)
    
    if is_sparse
        H1_sp = H1 isa SparseMatrixCSC ? H1 : sparse(H1)
        H2_sp = H2 isa SparseMatrixCSC ? H2 : sparse(H2)
        I_n2 = spdiagm(0 => fill(F_one, n2))
        I_m1 = spdiagm(0 => fill(F_one, m1))
        
        H_X = hcat(kron(H1_sp, I_n2), kron(I_m1, transpose(H2_sp)))
    else
        I_n2 = identity_matrix(F, n2)
        I_m1 = identity_matrix(F, m1)
        
        H_X = hcat(kronecker_product(H1, I_n2), kronecker_product(I_m1, transpose(H2)))
    end
    
    S.cache[:H_X] = H_X
    return H_X
end

function Z_stabilizers(S::HypergraphProductCode)
    haskey(S.cache, :H_Z) && return S.cache[:H_Z]
    
    F = S.F
    H1 = parity_check_matrix(S.C1)
    H2 = parity_check_matrix(S.C2)
    m1, n1 = size(H1)
    m2, n2 = size(H2)
    
    is_sparse = H1 isa SparseMatrixCSC || H2 isa SparseMatrixCSC
    F_one = F(1)
    
    if is_sparse
        H1_sp = H1 isa SparseMatrixCSC ? H1 : sparse(H1)
        H2_sp = H2 isa SparseMatrixCSC ? H2 : sparse(H2)
        I_n1 = spdiagm(0 => fill(F_one, n1))
        I_m2 = spdiagm(0 => fill(F_one, m2))
        
        if Int(characteristic(F)) == 2
            H_Z = hcat(kron(I_n1, H2_sp), kron(transpose(H1_sp), I_m2))
        else
            H_Z = hcat(kron(I_n1, H2_sp), kron(-transpose(H1_sp), I_m2))
        end
    else
        I_n1 = identity_matrix(F, n1)
        I_m2 = identity_matrix(F, m2)
        
        if Int(characteristic(F)) == 2
            H_Z = hcat(kronecker_product(I_n1, H2), kronecker_product(transpose(H1), I_m2))
        else
            H_Z = hcat(kronecker_product(I_n1, H2), kronecker_product(-transpose(H1), I_m2))
        end
    end
    
    S.cache[:H_Z] = H_Z
    return H_Z
end

function stabilizers(S::HypergraphProductCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    H_X = X_stabilizers(S)
    H_Z = Z_stabilizers(S)
    n = S.n
    F = S.F
    
    is_sparse = H_X isa SparseMatrixCSC || H_Z isa SparseMatrixCSC
    
    # Assemble standard symplectic parity check matrix
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
