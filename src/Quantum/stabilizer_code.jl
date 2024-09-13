# Copyright (c) 2021 - 2024 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    StabilizerCodeCSS(C1::AbstractLinearCode, C2::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)
    CSSCode(C1::AbstractLinearCode, C2::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

Return the CSS code given by the CSS construction on two linear codes `C1`
and `C2` with `C2 ⊆ C1`.

# Notes
* If `C1 = [n, k1, d1]` and `C2 = [n, k2, d2]`, the resulting code has dimension `k = k1 - k2`
  and minimum distance `d >= min(d1, d2^⟂)`. The `X` stabilizers are given by the parity-check matrix
  of `C2^⟂`, `H(C2^⟂)`, and the `Z` stabilizers by `H(C1)`.
"""
function StabilizerCodeCSS(C1::AbstractLinearCode, C2::AbstractLinearCode;
    char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

    logs_alg ∈ [:stnd_frm, :VS, :sys_eqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    C2 ⊆ C1 || throw(ArgumentError("The second argument must be a subset of the first in the CSS construction."))
    p = Int(characteristic(C1.F))
    char_vec = _process_char_vec(char_vec, p, 2 * C1.n)

    # C2 ⊆ C1
    # k = k1 - k2
    # d >= minimum(d1, d2^⟂)
    # X - H(C2^⟂), Z - H(C1)
    D2 = dual(C2)
    stabs = direct_sum(D2.H, C1.H)
    stabs_stand, P_stand, stand_r, stand_k, rnk = _standard_form_stabilizer(stabs)
    if !iszero(stand_k)
        if logs_alg == :stnd_frm
            logs = _make_pairs(_logicals_standard_form(stabs_stand, C1.n, stand_k, stand_r, P_stand))
            logs_mat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
        else
            logs, logs_mat = _logicals(stabs, direct_sum(C1.G, D2.G), logs_alg)
        end
    end
    signs, X_signs, Z_signs = _determine_signs_CSS(stabs, char_vec, nrows(D2.H), nrows(C1.H))

    # q^n / p^k but rows is n - k
    if !iszero(stand_k)
        dim_code = BigInt(order(C1.F))^C1.n // BigInt(p)^rnk
        isinteger(dim_code) && (dim_code = round(Int, log(BigInt(p), dim_code));)

        _, mat = rref(vcat(D2.H, C1.G))
        anti = _remove_empty(mat, :rows) * transpose(C1.G)
        u_bound_dx, _ = _min_wt_row(mat[findall(!iszero(anti[i:i, :]) for i in axes(anti, 1)), :])
        _, mat = rref(vcat(C1.H, D2.G))
        anti = _remove_empty(mat, :rows) * transpose(D2.G)
        u_bound_dz, _ = _min_wt_row(mat[findall(!iszero(anti[i:i, :]) for i in axes(anti, 1)), :])
        u_bound = min(u_bound_dx, u_bound_dz)
        if !ismissing(C1.d) && is !ismissing(D2.d)
            u_bound_classical = min(C1.d, D2.d)
            u_bound_classical < u_bound && (u_bound = u_bound_classical;)
        end
        return StabilizerCodeCSS(C1.F, C1.n, dim_code, missing, missing, missing, 1,
            u_bound, 1, u_bound_dx, 1, u_bound_dz, stabs, D2.H, C1.H, C2, C1, signs, X_signs,
            Z_signs, logs, logs_mat, char_vec, missing, missing, missing, false, missing,
            stabs_stand, stand_r, stand_k, P_stand, missing, missing)
    else
        # graph state distance is defined to be smallest weight stabilizer
        # TODO check if d's are known already from classical
        _, mat = rref(D2.H)
        mat = _remove_empty(mat, :rows)
        u_bound_dx, _ = _min_wt_row(mat)
        _, mat = rref(C1.H)
        mat = _remove_empty(mat, :rows)
        u_bound_dz, _ = _min_wt_row(Z_mat)
        u_bound = min(u_bound_dx, u_bound_dz)
        if !ismissing(C1.d) && is !ismissing(D2.d)
            u_bound_classical = min(C1.d, D2.d)
            u_bound_classical < u_bound && (u_bound = u_bound_classical;)
        end
        # TODO set upper bounds to distance if known
        return GraphStateStabilizerCSS(C1.F, C1.n, 0, missing, D2.d, C1.d, 1, u_bound, 1, u_bound_dx, 1, u_bound_dz, stabs, D2.H, C1.H, C2,
            C1, signs, X_signs, Z_signs, char_vec, missing, false, stabs_stand, stand_r, stand_k,
            P_stand, missing, missing)
    end
end
CSSCode(C1::AbstractLinearCode, C2::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing}
    = missing, logs_alg::Symbol = :stnd_frm) = StabilizerCodeCSS(C1, C2, char_vec = char_vec,
    logs_alg = logs_alg)

"""
    StabilizerCodeCSS(C::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)
    CSSCode(C::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

Return the CSS code given by the CSS construction on a self-orthogonal linear code
`C`, i.e., `C ⊆ C^⟂`.

# Notes
* Setting `C1 = C^⟂` and `C2 = C`, the resulting code has dimension `k = k1 - k2`
  and minimum distance `d >= min(d1, d2^⟂)`. The `X` stabilizers are given by the
  parity-check matrix of `C2^⟂`, `H(C2^⟂)`, and the `Z` stabilizers by `H(C1)`.
"""
function StabilizerCodeCSS(C::LinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm)

    logs_alg ∈ [:stnd_frm, :VS, :sys_eqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    # this should have Xstabs = Zstabs
    D = dual(C)
    C ⊆ D || throw(ArgumentError("The single code CSS construction requires C ⊆ C^⟂."))
    p = Int(characteristic(D.F))
    char_vec = _process_char_vec(char_vec, p, 2 * D.n)

    # C2 ⊆ C1
    # k = k1 - k2
    # d >= minimum(d1, d2^⟂)
    # X - H(C2^⟂), Z - H(C1)
    stabs = direct_sum(D.H, D.H)
    stabs_stand, P_stand, stand_r, stand_k, rnk = _standard_form_stabilizer(stabs)
    if !iszero(stand_k)
        if logs_alg == :stnd_frm
            logs = _make_pairs(_logicals_standard_form(stabs_stand, C.n, stand_k, stand_r, P_stand))
            logs_mat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
        else
            logs, logs_mat = _logicals(stabs, direct_sum(D.G, D.G), logs_alg)
        end
    end
    nr = nrows(D.H)
    signs, X_signs, Z_signs = _determine_signs_CSS(stabs, char_vec, nr, nr)

    # q^n / p^k but rows is n - k
    if !iszero(stand_k)
        dim_code = BigInt(order(D.F))^D.n // BigInt(p)^rnk
        isinteger(dim_code) && (dim_code = round(Int, log(BigInt(p), dim_code));)

        X_logs = reduce(vcat, [log[1][:, 1:n] for log in logs])
        Z_logs = reduce(vcat, [log[2][:, n + 1:end] for log in logs])
        _, mat = rref(vcat(D.H, X_logs))
        anti = _remove_empty(mat, :rows) * transpose(Z_logs)
        u_bound_dx, _ = _min_wt_row(mat[findall(!iszero(anti[i:i, :]) for i in axes(anti, 1)), :])
        _, mat = rref(vcat(D.H, Z_logs))
        anti = _remove_empty(mat, :rows) * transpose(X_logs)
        u_bound_dz, _ = _min_wt_row(mat[findall(!iszero(anti[i:i, :]) for i in axes(anti, 1)), :])
        return StabilizerCodeCSS(D.F, D.n, dim_code, missing, missing, missing, 1,
            min(u_bound_dx, u_bound_dz), 1, u_bound_dx, 1, u_bound_dz, stabs, D.H, D.H, C,
            D, signs, X_signs, Z_signs, logs, logs_mat, char_vec, missing, missing, missing, false,
            missing, stabs_stand, stand_r, stand_k, P_stand, missing, missing)
    else
        # graph state distance is defined to be smallest weight stabilizer
        if ismissing(D.d)
            _, mat = rref(D.H)
            mat = _remove_empty(mat, :rows)
            u_bound, _ = _min_wt_row(mat)
            return GraphStateStabilizerCSS(D.F, D.n, 0, missing, missing, missing, 1, u_bound, 1, 
                u_bound, 1, u_bound, stabs, D.H, D.H, C, D, signs, X_signs, Z_signs, char_vec,
                missing, false, stabs_stand, stand_r, stand_k, P_stand, missing, missing)
        else
            return GraphStateStabilizerCSS(D.F, D.n, 0, D.d, D.d, D.d, D.d, D.d, D.d, D.d, D.d,
                D.d, stabs, D.H, D.H, C, D, signs, X_signs, Z_signs, char_vec, missing, false,
                stabs_stand, stand_r, stand_k, P_stand, missing, missing)
        end
    end
end
CSSCode(C::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm) = StabilizerCodeCSS(C, char_vec = char_vec, logs_alg = logs_alg)

"""
    StabilizerCodeCSS(X_matrix::fqPolyRepMatrix, Z_matrix::fqPolyRepMatrix; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)
    CSSCode(X_matrix::fqPolyRepMatrix, Z_matrix::fqPolyRepMatrix; char_vec::Union{Vector{zzModRingElem}, Missing}= missing, logs_alg::Symbol = :stnd_frm)

Return a CSS code whose `X`-stabilizers are given by `X_matrix`, `Z`-stabilizers by `Z_matrix`.
"""
function StabilizerCodeCSS(X_matrix::T, Z_matrix::T; char_vec::Union{Vector{zzModRingElem}, Missing} = 
    missing, logs_alg::Symbol = :stnd_frm) where T <: CTMatrixTypes

    logs_alg ∈ [:stnd_frm, :VS, :sys_eqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    iszero(X_matrix) && throw(ArgumentError("The `X` stabilizer matrix is empty."))
    iszero(Z_matrix) && throw(ArgumentError("The `Z` stabilizer matrix is empty."))
    n = ncols(X_matrix)
    n == ncols(Z_matrix) || throw(ArgumentError("Both matrices must have the same length in the CSS construction."))
    F = base_ring(X_matrix)
    F == base_ring(Z_matrix) || throw(ArgumentError("Both matrices must be over the same base field in the CSS construction."))
    iszero(Z_matrix * transpose(X_matrix)) || throw(ArgumentError("The given matrices are not symplectic orthogonal."))
    p = Int(characteristic(F))
    char_vec = _process_char_vec(char_vec, p, 2 * n)

    X_matrix = _remove_empty(X_matrix, :rows)
    Z_matrix = _remove_empty(Z_matrix, :rows)

    # determine if the provided set of stabilizers are redundant
    X_rank = rank(X_matrix)
    Z_rank = rank(Z_matrix)
    nrows(X_matrix) > X_rank || nrows(Z_matrix) > Z_rank ? (over_comp = true) : (over_comp = false)
    stabs = direct_sum(X_matrix, Z_matrix)

    stabs_stand, P_stand, stand_r, stand_k, rnk = _standard_form_stabilizer(stabs)
    if !iszero(stand_k)
        if logs_alg == :stnd_frm
            logs = _make_pairs(_logicals_standard_form(stabs_stand, n, stand_k, stand_r, P_stand))
            logs_mat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
        else
            H = kernel(hcat(stabs[:, n + 1:end], -stabs[:, 1:n]), side = :right)
            rnk_H = rank(H)
            if ncols(H) == rnk_H
                H_tr = transpose(H)
            else
                # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
                nr = nrows(H)
                H_tr = zero_matrix(base_ring(H), rnk_H, nr)
                for r in 1:nr
                    for c in 1:rnk_H
                        !iszero(H[r, c]) && (H_tr[c, r] = H[r, c];)
                    end
                end
            end
            # n + (n - Srank)
            nrows(H) == 2 * n - X_rank - Z_rank || error("Normalizer matrix is not size n + k.")
            logs, logs_mat = _logicals(stabs, H_tr, logs_alg)
        end
    end
    signs, X_signs, Z_signs = _determine_signs_CSS(stabs, char_vec, nrows(X_matrix), nrows(Z_matrix))

    # q^n / p^k but rows is n - k
    if !iszero(stand_k)
        dim_code = BigInt(order(F))^n // BigInt(p)^rnk
        isinteger(dim_code) && (dim_code = round(Int, log(BigInt(p), dim_code));)

        X_logs = reduce(vcat, [log[1][:, 1:n] for log in logs])
        Z_logs = reduce(vcat, [log[2][:, n + 1:end] for log in logs])
        _, mat = rref(vcat(X_matrix, X_logs))
        anti = _remove_empty(mat, :rows) * transpose(Z_logs)
        u_bound_dx, _ = _min_wt_row(mat[findall(!iszero(anti[i:i, :]) for i in axes(anti, 1)), :])
        _, mat = rref(vcat(Z_matrix, Z_logs))
        anti = _remove_empty(mat, :rows) * transpose(X_logs)
        u_bound_dz, _ = _min_wt_row(mat[findall(!iszero(anti[i:i, :]) for i in axes(anti, 1)), :])
        return StabilizerCodeCSS(F, n, dim_code, missing, missing, missing, 1, min(u_bound_dx,
            u_bound_dz), 1, u_bound_dx, 1, u_bound_dz, stabs, X_matrix, Z_matrix,
            missing, missing, signs, X_signs, Z_signs, logs, logs_mat, char_vec, missing, missing,
            missing, over_comp, missing, stabs_stand, stand_r, stand_k, P_stand, missing, missing)
    else
        # graph state distance is defined to be smallest weight stabilizer
        _, mat = rref(X_matrix)
        mat = _remove_empty(mat, :rows)
        u_bound_dx, _ = _min_wt_row(mat)
        _, mat = rref(Z_matrix)
        mat = _remove_empty(mat, :rows)
        u_bound_dz, _ = _min_wt_row(mat)
        return GraphStateStabilizerCSS(F, n, 0, missing, missing, missing, 1, min(u_bound_dx,
            u_bound_dz), 1, u_bound_dx, 1, u_bound_dz, stabs, X_matrix, Z_matrix, missing,
            missing, signs, X_signs, Z_signs, char_vec, missing, over_comp, stabs_stand, stand_r,
            stand_k, P_stand, missing, missing)
    end
end
CSSCode(X_matrix::T, Z_matrix::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm) where T <: CTMatrixTypes = StabilizerCodeCSS(X_matrix, Z_matrix,
    char_vec = char_vec, logs_alg = logs_alg)

"""
    StabilizerCodeCSS(S_Pauli::Vector{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{String, Vector{Char}}
    CSSCode(S_Pauli::Vector{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{String, Vector{Char}}

Return the CSS code whose stabilizers are determined by the vector of Pauli strings `S_Pauli`.

# Notes
* Any +/- 1 characters in front of each stabilizer are stripped. No check is done
  to make sure these signs agree with the ones computed using the character vector.
"""
function StabilizerCodeCSS(S_Pauli::Vector{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm) where T <: Union{String, Vector{Char}}

    logs_alg ∈ [:stnd_frm, :VS, :sys_eqs] || throw(ArgumentError("Unrecognized logicals algorithm"))
    stabs = _Pauli_string_to_symplectic(_process_strings(S_Pauli))
    iszero(stabs) && throw(ArgumentError("The processed Pauli strings returned a set of empty stabilizer generators."))
    stabs = _remove_empty(stabs, :rows)
    # the reason we repeat here and not call another constructor is the else
    # statement at the bottom of this function
    # would also need to compute down to signs to call _is_CSS_symplectic
    # which would allow us to call the other constructor
    are_symplectic_orthogonal(stabs, stabs) || throw(ArgumentError("The given stabilizers are not symplectic orthogonal."))
    
    n = div(ncols(stabs), 2)
    stabs_stand, P_stand, stand_r, stand_k, rnk = _standard_form_stabilizer(stabs)
    if !iszero(stand_k)
        if logs_alg == :stnd_frm
            logs = _make_pairs(_logicals_standard_form(stabs_stand, n, stand_k, stand_r, P_stand))
            logs_mat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
        else
            H = kernel(hcat(stabs[:, n + 1:end], -stabs[:, 1:n]), side = :right)
            rnk_H = rank(H)
            if ncols(H) == rnk_H
                H_tr = transpose(H)
            else
                # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
                nr = nrows(H)
                H_tr = zero_matrix(base_ring(H), rnk_H, nr)
                for r in 1:nr
                    for c in 1:rnk_H
                        !iszero(H[r, c]) && (H_tr[c, r] = H[r, c];)
                    end
                end
            end
            # n + (n - Srank)
            nrows(H_tr) == 2 * n - rnk || error("Normalizer matrix is not size n + k.")
            logs, logs_mat = _logicals(stabs, H_tr, logs_alg)
        end
    end

    F = base_ring(stabs)
    p = Int(characteristic(F))
    char_vec = _process_char_vec(char_vec, p, 2 * n)
    signs = _determine_signs(stabs, char_vec)
    over_comp = nrows(stabs) > rnk

    # q^n / p^k but rows is n - k
    args = _is_CSS_symplectic(stabs, signs, true)
    if args[1]
        if !iszero(stand_k)
            dim_code = BigInt(order(F))^n // BigInt(p)^rnk
            isinteger(dim_code) && (dim_code = round(Int, log(BigInt(p), dim_code));)

            X_logs = reduce(vcat, [log[1][:, 1:n] for log in logs])
            Z_logs = reduce(vcat, [log[2][:, n + 1:end] for log in logs])
            _, mat = rref(vcat(args[2], X_logs))
            anti = _remove_empty(mat, :rows) * transpose(Z_logs)
            u_bound_dx, _ = _min_wt_row(mat[findall(!iszero(anti[i:i, :]) for i in axes(anti, 1)), :])
            _, mat = rref(vcat(args[4], Z_logs))
            anti = _remove_empty(mat, :rows) * transpose(X_logs)
            u_bound_dz, _ = _min_wt_row(mat[findall(!iszero(anti[i:i, :]) for i in axes(anti, 1)), :])
            return StabilizerCodeCSS(F, n, dim_code, missing, missing, missing, 1, min(u_bound_dx,
                u_bound_dz), 1, u_bound_dx, 1, u_bound_dz, stabs, args[2], args[4], missing,
                missing, signs, args[3], args[5], logs, logs_mat, char_vec, missing, missing,
                missing, over_comp, missing, stabs_stand, stand_r, stand_k, P_stand, missing,
                missing)
        else
            # graph state distance is defined to be smallest weight stabilizer
            _, mat = rref(args[2])
            mat = _remove_empty(mat, :rows)
            u_bound_dx, _ = _min_wt_row(mat)
            _, mat = rref(args[4])
            mat = _remove_empty(mat, :rows)
            u_bound_dz, _ = _min_wt_row(mat)
            return GraphStateStabilizerCSS(F, n, 0, missing, missing, missing, 1, min(u_bound_dx,
                u_bound_dz), 1, u_bound_dx, 1, u_bound_dz, stabs, args[2], args[4], missing,
                missing, signs, args[3], args[5], char_vec, missing, over_comp, stabs_stand,
                stand_r, stand_k, P_stand, missing, missing)
        end
    else
        error("Provided Pauli strings are not CSS.")
    end
end
CSSCode(S_Pauli::Vector{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm) where T <: Union{String, Vector{Char}} =
    StabilizerCodeCSS(S_Pauli, char_vec = char_vec, logs_alg = logs_alg)

"""
    StabilizerCodeCSS(S::AbstractStabilizerCode; logs_alg::Symbol = :stnd_frm)
    CSSCode(S::AbstractStabilizerCode; logs_alg::Symbol = :stnd_frm)

Return the `[[2n, 2k, d <= d' <= 2d]]` CSS code derived from the `[n, k, d]` stabilizer
code `S` using the technique of arXiv:2311.18003.
"""
function StabilizerCodeCSS(S::AbstractStabilizerCode; logs_alg::Symbol = :stnd_frm)
    # correct size, fastest object creation, reset all entries anyway
    Z = deepcopy(S.stabs)
    Z[:, 1:S.n] = S.stabs[:, S.n + 1:2 * S.n]
    Z[:, S.n + 1:2 * S.n] = -S.stabs[:, 1:S.n]
    return CSSCode(S.stabs, Z, logs_alg = logs_alg)
end
CSSCode(S::AbstractStabilizerCode; logs_alg::Symbol = :stnd_frm) = StabilizerCodeCSS(S; logs_alg =
    logs_alg)

"""
    StabilizerCode(S::AbstractStabilizerCode; logs_alg::Symbol = :stnd_frm)

Return the `[[n, k, d' <= d]]` stabilizer code from the `[[2n, 2k, d]]` CSS code using the
inverse of the technique of arXiv:2311.18003. 
"""
StabilizerCode(S::T; logs_alg::Symbol = :stnd_frm) where {T <: AbstractStabilizerCode} =
    StabilizerCode(CSSTrait(T), S, logs_alg)
function StabilizerCode(::IsCSS, S::AbstractStabilizerCode, logs_alg::Symbol)
    iseven(S.n) || throw(ArgumentError("Only valid for codes of even length"))
    nrows(S.X_stabs) == nrows(S.Z_stabs) || throw(ArgumentError("Requires an equal number of X and Z stabilizers"))

    n = div(S.n, 2)
    return StabilizerCode(hcat(S.X_stabs[:, 1:n], S.Z_stabs[:, 1:n]), logs_alg = logs_alg)
end
StabilizerCode(::IsNotCSS, S::AbstractStabilizerCode, logs_alg::Symbol) = throw(ArgumentError("Only valid for CSS codes of even length"))

# entanglement-assisted is not symplectic orthogonal
"""
    StabilizerCode(S_Pauli::Vector{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{String, Vector{Char}}

Return the stabilizer code whose stabilizers are determined by the vector of Pauli strings `S_Pauli`.

# Notes
* Any +/- 1 characters in front of each stabilizer are stripped. No check is done
  to make sure these signs agree with the ones computed using the character vector.
"""
function StabilizerCode(S_Pauli::Vector{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm) where T <: Union{String, Vector{Char}}

    S_Pauli_stripped = _process_strings(S_Pauli)
    stabs = _Pauli_string_to_symplectic(S_Pauli_stripped)
    iszero(stabs) && throw(ArgumentError("The processed Pauli strings returned a set of empty stabilizer generators."))
    return StabilizerCode(stabs, char_vec = char_vec, logs_alg = logs_alg)
end

"""
    StabilizerCode(stabs::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

Return the stabilizer code whose stabilizers is determined by `stabs`.
"""
function StabilizerCode(stabs::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm)

    logs_alg ∈ (:stnd_frm, :VS, :sys_eqs) || throw(ArgumentError("Unrecognized logicals algorithm"))
    iszero(stabs) && throw(ArgumentError("The stabilizer matrix is empty."))
    stabs = _remove_empty(stabs, :rows)
    are_symplectic_orthogonal(stabs, stabs) || throw(ArgumentError("The given stabilizers are not symplectic orthogonal."))
    
    n = div(ncols(stabs), 2)
    stabs_stand, P_stand, stand_r, stand_k, rnk = _standard_form_stabilizer(stabs)
    if !iszero(stand_k)
        if logs_alg == :stnd_frm
            logs = _make_pairs(_logicals_standard_form(stabs_stand, n, stand_k, stand_r, P_stand))
            logs_mat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
        else
            H = kernel(hcat(stabs[:, n + 1:end], -stabs[:, 1:n]), side = :right)
            rnk_H = rank(H)
            if ncols(H) == rnk_H
                H_tr = transpose(H)
            else
                # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
                nr = nrows(H)
                H_tr = zero_matrix(base_ring(H), rnk_H, nr)
                for r in 1:nr
                    for c in 1:rnk_H
                        !iszero(H[r, c]) && (H_tr[c, r] = H[r, c];)
                    end
                end
            end
            # n + (n - Srank)
            nrows(H_tr) == 2 * n - rnk || error("Normalizer matrix is not size n + k.")
            logs, logs_mat = _logicals(stabs, H_tr, logs_alg)
        end
    end

    F = base_ring(stabs)
    p = Int(characteristic(F))
    char_vec = _process_char_vec(char_vec, p, 2 * n)
    signs = _determine_signs(stabs, char_vec)
    nrows(stabs) > rnk ? (over_comp = true) : (over_comp = false)

    # q^n / p^k but rows is n - k
    dim_code = BigInt(order(F))^n // BigInt(p)^rnk
    isinteger(dim_code) && (dim_code = round(Int, log(BigInt(p), dim_code));)

    args = _is_CSS_symplectic(stabs, signs, true)
    if args[1]
        if !iszero(stand_k)
            X_logs = reduce(vcat, [log[1][:, 1:n] for log in logs])
            Z_logs = reduce(vcat, [log[2][:, n + 1:end] for log in logs])
            _, mat = rref(vcat(args[2], X_logs))
            anti = _remove_empty(mat, :rows) * transpose(Z_logs)
            u_bound_dx, _ = _min_wt_row(mat[findall(!iszero(anti[i:i, :]) for i in axes(anti, 1)), :])
            _, mat = rref(vcat(args[4], Z_logs))
            anti = _remove_empty(mat, :rows) * transpose(X_logs)
            u_bound_dz, _ = _min_wt_row(mat[findall(!iszero(anti[i:i, :]) for i in axes(anti, 1)), :])
            return StabilizerCodeCSS(F, n, dim_code, missing, missing, missing, 1, min(u_bound_dx,
                u_bound_dz), 1, u_bound_dx, 1, u_bound_dz, stabs, args[2],  args[4], missing,
                missing, signs, args[3], args[5], logs, logs_mat, char_vec, missing, missing,
                missing, over_comp, missing, stabs_stand, stand_r, stand_k, P_stand, missing,
                missing)
        else
            # graph state distance is defined to be smallest weight stabilizer
            _, mat = rref(args[2])
            mat = _remove_empty(mat, :rows)
            u_bound_dx, _ = _min_wt_row(mat)
            _, mat = rref(args[4])
            mat = _remove_empty(mat, :rows)
            u_bound_dz, _ = _min_wt_row(mat)
            return GraphStateStabilizerCSS(F, n, 0, missing, missing, missing, 1, min(u_bound_dx,
                u_bound_dz), 1, u_bound_dx, 1, u_bound_dz, stabs, args[2], args[4], missing,
                missing, signs, args[3], args[5], char_vec, missing, over_comp, stabs_stand,
                stand_r, stand_k, P_stand, missing, missing)
        end
    else
        if !iszero(stand_k)
            # TODO symplectic multiply here
            _, mat = rref(vcat(stabs, logs_mat))
            mat = _remove_empty(mat, :rows)
            u_bound, _ = _min_wt_row(mat)
            return StabilizerCode(F, n, dim_code, missing, 1, u_bound, stabs, logs, logs_mat, 
                char_vec, signs, missing, missing, missing, over_comp, missing, stabs_stand,
                stand_r, stand_k, P_stand, missing)
        else
            # graph state distance is defined to be smallest weight stabilizer
            _, mat = rref(stabs)
            mat = _remove_empty(mat, :rows)
            u_bound, _ = _min_wt_row(mat)
            return GraphStateStabilizer(F, n, 0, missing, 1, u_bound, stabs, char_vec, signs,
                missing, over_comp, stabs_stand, stand_r, stand_k, P_stand, missing)
        end
    end
end

"""
    StabilizerCode(S::AbstractSubsystemCode)

Return the stabilizer code from the subsystem code by promoting all of the gauge operators to logical operators.
"""
function StabilizerCode(S::AbstractSubsystemCode)
    typeof(S) <: AbstractStabilizerCode && return S
    return promote_gauges_to_logical(S, collect(1:length(S.gauge_ops)))
end
# function StabilizerCodeCSS(S::AbstractSubsystemCodeCSS)
#     # promote all gauges to logicals
# end

#############################
      # getter functions
#############################

"""
    minimum_distance_lower_bound(S::AbstractStabilizerCode)
Return the currently stored lower bound on the minimum distance.
"""
minimum_distance_lower_bound(S::AbstractStabilizerCode) = S.l_bound

"""
    minimum_distance_upper_bound(S::AbstractStabilizerCode)
Return the currently stored upper bound on the minimum distance.
"""
minimum_distance_upper_bound(S::AbstractStabilizerCode) = S.u_bound

"""
    X_minimum_distance_lower_bound(S::AbstractStabilizerCode)
Return the currently stored lower bound on the minimum `X`-distance.
"""
X_minimum_distance_lower_bound(S::T) where T <: AbstractStabilizerCode = X_minimum_distance_lower_bound(CSSTrait(T), S)
X_minimum_distance_lower_bound(::IsCSS, S::AbstractStabilizerCode) = S.l_bound_dx
X_minimum_distance_lower_bound(::IsNotCSS, S::AbstractStabilizerCode) = error("Only valid for CSS codes")

"""
    X_minimum_distance_upper_bound(S::AbstractStabilizerCode)
Return the currently stored upper bound on the minimum `X`-distance.
"""
X_minimum_distance_upper_bound(S::T) where T <: AbstractStabilizerCode = X_minimum_distance_upper_bound(CSSTrait(T), S)
X_minimum_distance_upper_bound(::IsCSS, S::AbstractStabilizerCode) = S.u_bound_dx
X_minimum_distance_upper_bound(::IsNotCSS, S::AbstractStabilizerCode) = error("Only valid for CSS codes")

"""
    Z_minimum_distance_lower_bound(S::AbstractStabilizerCode)
Return the currently stored lower bound on the minimum `Z`-distance.
"""
Z_minimum_distance_lower_bound(S::T) where T <: AbstractStabilizerCode = Z_minimum_distance_lower_bound(CSSTrait(T), S)
Z_minimum_distance_lower_bound(::IsCSS, S::AbstractStabilizerCode) = S.l_bound_dz
Z_minimum_distance_lower_bound(::IsNotCSS, S::AbstractStabilizerCode) = error("Only valid for CSS codes")

"""
    Z_minimum_distance_upper_bound(S::AbstractStabilizerCode)
Return the currently stored upper bound on the minimum `Z`-distance.
"""
Z_minimum_distance_upper_bound(S::T) where T <: AbstractStabilizerCode = Z_minimum_distance_upper_bound(CSSTrait(T), S)
Z_minimum_distance_upper_bound(::IsCSS, S::AbstractStabilizerCode) = S.u_bound_dz
Z_minimum_distance_upper_bound(::IsNotCSS, S::AbstractStabilizerCode) = error("Only valid for CSS codes")

#############################
      # setter functions
#############################

"""
    set_minimum_distance!(S::AbstractStabilizerCode, d::Int)

Set the minimum distance of the code to `d`.
"""
function set_minimum_distance!(S::AbstractStabilizerCode, d::Int)
    0 < d <= S.n || throw(DomainError("The minimum distance of a code must be ≥ 1; received: d = $d."))
    S.u_bound < d && (@warn "The distance set is greater than the current upper bound of $(S.u_bound)")
    d < S.l_bound && (@warn "The distance set is less than the current lower bound of $(S.l_bound)")
    S.d = d
    S.l_bound = d
    S.u_bound = d

    if CSSTrait(typeof(S)) == IsCSS()
        ismissing(S.dx) && (S.l_bound_dx = d;)
        ismissing(S.dz) && (S.l_bound_dz = d;)
    end
    return nothing
end

"""
    set_X_minimum_distance!(S::AbstractStabilizerCode, d::Int)

Set the minimum `X`-distance of the code to `d`.
"""
set_X_minimum_distance!(S::T, d::Int) where T <: AbstractStabilizerCode = set_X_minimum_distance!(CSSTrait(T), S, d)
function set_X_minimum_distance!(::IsCSS, S::AbstractStabilizerCode, d::Int)
    0 < d <= S.n || throw(DomainError("The minimum distance of a code must be ≥ 1; received: d = $d."))
    S.u_bound_dx < d && (@warn "The distance set is greater than the current upper bound of $(S.u_bound_dx)")
    d < S.l_bound_dx && (@warn "The distance set is less than the current lower bound of $(S.l_bound_dx)")
    S.dx = d
    S.l_bound_dx = d
    S.u_bound_dx = d

    if !ismissing(S.dz)
        S.d = min(S.dx, S.dz)
        S.u_bound = S.d
        S.l_bound = S.d
    else
        S.dx < S.u_bound && (S.u_bound = S.dx;)
    end
    return nothing
end
set_X_minimum_distance!(::IsNotCSS, S::AbstractStabilizerCode, d::Int) =
    error("Only valid for CSS codes")

"""
    set_Z_minimum_distance!(S::AbstractStabilizerCode, d::Int)

Set the minimum `Z`-distance of the code to `d`.
"""
set_Z_minimum_distance!(S::T, d::Int) where T <: AbstractStabilizerCode = set_Z_minimum_distance!(CSSTrait(T), S, d)
function set_Z_minimum_distance!(::IsCSS, S::AbstractStabilizerCode, d::Int)
    0 < d <= S.n || throw(DomainError("The minimum distance of a code must be ≥ 1; received: d = $d."))
    S.u_bound_dz < d && (@warn "The distance set is greater than the current upper bound of $(S.u_bound_dz)")
    d < S.l_bound_dz && (@warn "The distance set is less than the current lower bound of $(S.l_bound_dz)")
    S.dz = d
    S.l_bound_dz = d
    S.u_bound_dz = d

    if !ismissing(S.dx)
        S.d = min(S.dx, S.dz)
        S.u_bound = S.d
        S.l_bound = S.d
    else
        S.dz < S.u_bound && (S.u_bound = S.dz;)
    end
    return nothing
end
set_Z_minimum_distance!(::IsNotCSS, S::AbstractStabilizerCode, d::Int) =
    error("Only valid for CSS codes")

#############################
     # general functions
#############################

function _logicals(stabs::T, dual_gens::T, logs_alg::Symbol = :sys_eqs) where {T <: CTMatrixTypes}
    logs_alg ∈ [:sys_eqs, :VS] || throw(ArgumentError("Unrecognized logicals algorithm"))

    L = _quotient_space(dual_gens, stabs, logs_alg)
    logs = _make_pairs(L)
    # verify
    n = div(ncols(L), 2)
    logs_mat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
    are_symplectic_orthogonal(stabs, logs_mat) || error("Computed logicals do not commute with the codespace.")
    prod = hcat(logs_mat[:, n + 1:end], -logs_mat[:, 1:n]) * transpose(logs_mat)
    sum(_Flint_matrix_to_Julia_int_matrix(prod), dims=1) == ones(Int, 1, size(prod, 1)) ||
        error("Computed logicals do not have the right commutation relations.")
    return logs, logs_mat
end

"""
   random_CSS_code(n::Int, k::Int)

Return a random CSS code with an equal number of X and Z stabilizers.
"""
function random_CSS_code(n::Int, k::Int)
    d = _rand_single_sector_boundary(n, k)
    dt = transpose(d)
    Xrref = rref(dt)
    Xind = [findfirst(!iszero, Xrref[2][i, :])[2] for i in 1:Xrref[1]]
    Zrref = rref(d)
    Zind = [findfirst(!iszero, Zrref[2][i, :])[2] for i in 1:Zrref[1]]
    return CSSCode(d[Xind, :], dt[Zind, :])
end

"""
    is_CSS_T_code(S::AbstractStabilizerCode)

Return `true` if `S` is a CSS-T code.
"""
is_CSS_T_code(S::AbstractStabilizerCode; verbose::Bool = false) = is_CSS_T_code(CSSTrait(
    typeof(S)), S; verbose = verbose)
function is_CSS_T_code(::IsCSS, S::AbstractStabilizerCode; verbose::Bool = false)
    C_X = LinearCode(S.X_stabs) # C_2
    C_Z = LinearCode(S.Z_stabs) # dual(C_1)
    
    num_thrds = Threads.nthreads()
    verbose && println("Detected $num_thrds threads.")
    # should probably floor some logs here
    power = 0
    for i in 1:20
        if 2^i > num_thrds
            power = i - 1
            break
        end
    end

    p = Int(characteristic(C_X.F))
    n = S.n
    k = C_X.k
    z = zeros(Int, n)
    G = _Flint_matrix_to_Julia_int_matrix(transpose(C_X.G))
    for r in 1:k
        verbose && println("r: $r")
        flag = Threads.Atomic{Bool}(true)
        Threads.@threads for m in 1:num_thrds # 
            # println("here")
            c = deepcopy(z)
            # prefix = digits(m - 1, base = 2, pad = power)
            # println(k, ", ", r, ", ", prefix)
            # TODO: fix this for SteaneCode - prefix is larger than k so iterator is empty
            # for u in GrayCode(k, r, prefix, mutate = true)
            for u in GrayCode(k, r, Int[], mutate = true)
                # println("u:", u)
                if flag[]
                    LinearAlgebra.mul!(c, G, u)
                    wt_c = 0
                    supp_c = Int[]
                    @inbounds for j in 1:n
                        c[j] % p != 0 && (wt_c += 1; append!(supp_c, j);)
                    end

                     # for binary codes we can do this with a homomorphism, but this is non-linear for q-ary
                    # since we have to enumerate the codes anyway, we might as well just do this
                    # note, this is not even-like for q-ary
                    if !iseven(wt_c)
                        Threads.atomic_cas!(flag, true, false)
                        verbose && println("not even weight")
                        break
                    end
                    verbose && println("shortening on ", setdiff(1:n, supp_c))
                    # shorten C_Z on the complement of supp_x
                    temp = shorten(C_Z, setdiff(1:n, supp_c))
                    # println("thread: $m ", temp)
                    # display(temp)
                    if !contains_self_dual_subcode(temp)
                        Threads.atomic_cas!(flag, true, false)
                        verbose && println("m: $m, no subcode")
                        break
                    end
                else
                    break
                end
            end
            if !flag[]
                break
            end
        end
        if !flag[]
            return false
        end
    end
    return true
end
is_CSS_T_code(::IsNotCSS, S::AbstractStabilizerCode; verbose::Bool = false) = false

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
