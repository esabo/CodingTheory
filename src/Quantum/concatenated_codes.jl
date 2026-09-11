# Copyright (c) 2022 - 2026 Eric Sabo, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
    concatenate(S_outer::AbstractStabilizerCode, S_inner::AbstractStabilizerCode)

Return a lazy `QuantumConcatenatedCode`. Computes parameters instantly 
without generating the massive concatenated stabilizer matrices.
"""
function _family_char_vec(char_vec, F, n)
    return ismissing(char_vec) ? zzModRingElem[] :
        _process_char_vec(char_vec, Int(characteristic(F)), 2n)
end

function _family_cache(n; d=missing, l_bound=1, u_bound=n, logs_alg=:stnd_frm)
    return Dict{Symbol, Any}(
        :d => d, :l_bound => l_bound, :u_bound => u_bound,
        :logs_alg => logs_alg)
end

function concatenate(S_outer::AbstractStabilizerCode, S_inner::AbstractStabilizerCode;
    char_vec::Union{Vector{zzModRingElem}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm)
    S_inner.k == 1 || throw(ArgumentError("The inner code must encode exactly 1 logical qubit (k=1)."))
    F = field(S_outer)
    F == field(S_inner) || throw(ArgumentError("The codes must be over the same field."))
    
    # 1. Physical and Logical Qubits
    n_new = S_outer.n * S_inner.n
    k_new = S_outer.k 
    
    # 2. Distance Bounds (Multiplicative)
    l_bound = S_outer.l_bound * S_inner.l_bound
    u_bound = S_outer.u_bound * S_inner.u_bound
    
    d_exact = missing
    if !ismissing(S_outer.d) && !ismissing(S_inner.d)
        d_exact = S_outer.d * S_inner.d
    end
    
    cache = _family_cache(n_new; d=d_exact, l_bound=l_bound,
        u_bound=u_bound, logs_alg=logs_alg)
    result = QuantumConcatenatedCode(F, S_outer, S_inner, n_new, k_new,
        _family_char_vec(char_vec, F, n_new), cache)
    return _seed_quantum_singleton_bound!(result)
end

function stabilizers(S::QuantumConcatenatedCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    S_outer = S.outer_code
    S_inner = S.inner_code
    F = S.F
    
    n_out = S_outer.n
    n_in = S_inner.n
    n_new = S.n
    
    stabs_out = stabilizers(S_outer)
    stabs_in = stabilizers(S_inner)
    
    L_in = logicals_matrix(S_inner)
    L_X_in = L_in[1:1, :]
    L_Z_in = L_in[2:2, :]
    
    num_stabs_in = nrows(stabs_in)
    num_stabs_out = nrows(stabs_out)
    
    is_sparse = stabs_out isa SparseMatrixCSC || stabs_in isa SparseMatrixCSC
    
    if is_sparse
        stabs_new_inner = spzeros(F, n_out * num_stabs_in, 2 * n_new)
        stabs_new_outer = spzeros(F, num_stabs_out, 2 * n_new)
    else
        stabs_new_inner = zero_matrix(F, n_out * num_stabs_in, 2 * n_new)
        stabs_new_outer = zero_matrix(F, num_stabs_out, 2 * n_new)
    end
    
    # 1. Inner Code Stabilizers (Applied to each block)
    curr_row = 1
    for i in 1:n_out
        offset_X = (i - 1) * n_in
        offset_Z = n_new + (i - 1) * n_in
        for r in 1:num_stabs_in
            for c in 1:n_in
                stabs_new_inner[curr_row, offset_X + c] = stabs_in[r, c]
                stabs_new_inner[curr_row, offset_Z + c] = stabs_in[r, n_in + c]
            end
            curr_row += 1
        end
    end
    
    # 2. Outer Code Stabilizers (Mapped through inner logicals)
    for r in 1:num_stabs_out
        for i in 1:n_out
            x_val = stabs_out[r, i]
            z_val = stabs_out[r, n_out + i]
            
            offset_X = (i - 1) * n_in
            offset_Z = n_new + (i - 1) * n_in
            for c in 1:n_in
                # Symplectic dot product applying the outer check to the inner logical block
                stabs_new_outer[r, offset_X + c] += x_val * L_X_in[1, c] + z_val * L_Z_in[1, c]
                stabs_new_outer[r, offset_Z + c] += x_val * L_X_in[1, n_in + c] + z_val * L_Z_in[1, n_in + c]
            end
        end
    end
    
    stabs_final = vcat(stabs_new_inner, stabs_new_outer)
    
    S.cache[:stabilizers] = stabs_final
    return stabs_final
end
