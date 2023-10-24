# Copyright (c) 2021, 2022, 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    StabilizerCodeCSS(C1::AbstractLinearCode, C2::AbstractLinearCode, char_vec::Union{Vector{nmod}, Missing}=missing)
    CSSCode(C1::AbstractLinearCode, C2::AbstractLinearCode, char_vec::Union{Vector{nmod}, Missing}=missing)

Return the CSS code given by the CSS construction on two linear codes `C1`
and `C2` with `C2 ⊆ C1` and whose signs by `char_vec`.

# Notes
* If `C1 = [n, k1, d1]` and `C2 = [n, k2, d2]`, the resulting code has dimension `k = k1 - k2`
  and minimum distance `d >= min(d1, d2^⟂)`. The `X` stabilizers are given by the parity-check matrix
  of `C2^⟂`, `H(C2^⟂)`, and the `Z` stabilizers by `H(C1)`.
"""
function StabilizerCodeCSS(C1::AbstractLinearCode, C2::AbstractLinearCode,
    char_vec::Union{Vector{nmod}, Missing}=missing, logs_alg::Symbol=:stnd_frm)

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

        return StabilizerCodeCSS(C1.F, C1.n, dim_code, missing, missing, missing, stabs, D2.H, C1.H,
            C2, C1, signs, X_signs, Z_signs, logs, logs_mat, char_vec, missing, missing, missing,
            false, missing, stabs_stand, stand_r, stand_k, P_stand)
    else
        return GraphStateStabilizerCSS(C1.F, C1.n, 0, missing, D2.d, C1.d, stabs, D2.H, C1.H, C2, C1,
            signs, X_signs, Z_signs, char_vec, missing, false, stabs_stand, stand_r, stand_k, P_stand)
    end
end
CSSCode(C1::AbstractLinearCode, C2::AbstractLinearCode, char_vec::Union{Vector{nmod}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm) = StabilizerCodeCSS(C1, C2, char_vec, logs_alg)

"""
    StabilizerCodeCSS(C::AbstractLinearCode, char_vec::Union{Vector{nmod}, Missing}=missing)
    CSSCode(C::AbstractLinearCode, char_vec::Union{Vector{nmod}, Missing}=missing)

Return the CSS code given by the CSS construction on a self-orthogonal linear code
`C`, i.e., `C ⊆ C^⟂`, and whose signs by `char_vec`.

# Notes
* Setting `C1 = C^⟂` and `C2 = C`, the resulting code has dimension `k = k1 - k2`
  and minimum distance `d >= min(d1, d2^⟂)`. The `X` stabilizers are given by the
  parity-check matrix of `C2^⟂`, `H(C2^⟂)`, and the `Z` stabilizers by `H(C1)`.
"""
function StabilizerCodeCSS(C::LinearCode, char_vec::Union{Vector{nmod}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm)

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

        return StabilizerCodeCSS(D.F, D.n, dim_code, missing, missing, missing, stabs, D.H, D.H, C,
            D, signs, X_signs, Z_signs, logs, logs_mat, char_vec, missing, missing, missing, false,
            missing, stabs_stand, stand_r, stand_k, P_stand)
    else
        return GraphStateStabilizerCSS(D.F, D.n, 0, missing, D.d, D.d, stabs, D.H, D.H, C, D, signs,
            X_signs, Z_signs, char_vec, missing, false, stabs_stand, stand_r, stand_k, P_stand)
    end
end
CSSCode(C::AbstractLinearCode, char_vec::Union{Vector{nmod}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm) = StabilizerCodeCSS(C, char_vec, logs_alg)

"""
    StabilizerCodeCSS(X_matrix::fq_nmod_mat, Z_matrix::fq_nmod_mat, char_vec::Union{Vector{nmod}, Missing}=missing)
    CSSCode(X_matrix::fq_nmod_mat, Z_matrix::fq_nmod_mat, char_vec::Union{Vector{nmod}, Missing}=missing)

Return a CSS code whose `X`-stabilizers are given by `X_matrix`, `Z`-stabilizers by `Z_matrix`, and signs by `char_vec`.
"""
function StabilizerCodeCSS(X_matrix::T, Z_matrix::T, char_vec::Union{Vector{nmod}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm) where T <: CTMatrixTypes

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
            rnk_H, H = right_kernel(hcat(stabs[:, n + 1:end], -stabs[:, 1:n]))
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

        return StabilizerCodeCSS(F, n, dim_code, missing, missing, missing, stabs, X_matrix, Z_matrix,
            missing, missing, signs, X_signs, Z_signs, logs, logs_mat, char_vec, missing, missing,
            missing, over_comp, missing, stabs_stand, stand_r, stand_k, P_stand)
    else
        return GraphStateStabilizerCSS(F, n, 0, missing, missing, missing, stabs, X_matrix, Z_matrix,
            missing, missing, signs, X_signs, Z_signs, char_vec, missing, over_comp, stabs_stand, stand_r,
            stand_k, P_stand)
    end
end
CSSCode(X_matrix::T, Z_matrix::T, char_vec::Union{Vector{nmod}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm) where T <: CTMatrixTypes = StabilizerCodeCSS(X_matrix, Z_matrix, char_vec, logs_alg)

"""
    StabilizerCodeCSS(S_Pauli::Vector{T}, char_vec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}
    CSSCode(S_Pauli::Vector{T}, char_vec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}

Return the CSS code whose stabilizers are determined by the vector of Pauli strings `S_Pauli` and signs by `char_vec`.

# Notes
* Any +/- 1 characters in front of each stabilizer are stripped. No check is done
  to make sure these signs agree with the ones computed using the character vector.
"""
function StabilizerCodeCSS(S_Pauli::Vector{T}, char_vec::Union{Vector{nmod}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm) where T <: Union{String, Vector{Char}}

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
            rnk_H, H = right_kernel(hcat(stabs[:, n + 1:end], -stabs[:, 1:n]))
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

            return StabilizerCodeCSS(F, n, dim_code, missing, missing, missing, stabs, args[2],
                args[4], missing, missing, signs, args[3], args[5], logs, logs_mat, char_vec,
                missing, missing, missing, over_comp, missing, stabs_stand, stand_r, stand_k, P_stand)
        else
            return GraphStateStabilizerCSS(F, n, 0, missing, missing, missing, stabs, args[2],
                args[4], missing, missing, signs, args[3], args[5], char_vec, missing, over_comp, stabs_stand,
                stand_r, stand_k, P_stand)
        end
    else
        error("Provided Pauli strings are not CSS.")
    end
end
CSSCode(S_Pauli::Vector{T}, char_vec::Union{Vector{nmod}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm) where T <: Union{String, Vector{Char}} = StabilizerCodeCSS(S_Pauli, char_vec, logs_alg)

"""
    StabilizerCodeCSS(S::StabilizerCode)
    CSSCode(S::StabilizerCode)

Return the `[[2n, 2k, S.d <= d <= 2 S.d]]` CSS code derived by splitting the stabilizers of `S`.
"""
StabilizerCodeCSS(S::StabilizerCode, logs_alg::Symbol=:stnd_frm) = return StabilizerCodeCSS(S.stabs[:, 1:S.n], S.stabs[:, S.n+1:end], S.char_vec, logs_alg)
CSSCode(S::StabilizerCode, logs_alg::Symbol=:stnd_frm) = StabilizerCodeCSS(S, logs_alg)

# entanglement-assisted is not symplectic orthogonal
"""
    StabilizerCode(S_Pauli::Vector{T}, char_vec::Union{Vector{nmod}, Missing}=missing) where T <: Union{String, Vector{Char}}

Return the stabilizer code whose stabilizers are determined by the vector of Pauli strings `S_Pauli` and signs by `char_vec`.

# Notes
* Any +/- 1 characters in front of each stabilizer are stripped. No check is done
  to make sure these signs agree with the ones computed using the character vector.
"""
function StabilizerCode(S_Pauli::Vector{T}, char_vec::Union{Vector{nmod}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm) where T <: Union{String, Vector{Char}}

    S_Pauli_stripped = _process_strings(S_Pauli)
    stabs = _Pauli_string_to_symplectic(S_Pauli_stripped)
    iszero(stabs) && throw(ArgumentError("The processed Pauli strings returned a set of empty stabilizer generators."))
    return StabilizerCode(stabs, char_vec, logs_alg)
end

"""
    StabilizerCode(S::fq_nmod_mat, char_vec::Union{Vector{nmod}, Missing}=missing)

Return the stabilizer code whose stabilizers is determined by `S` and signs by `char_vec`.
"""
function StabilizerCode(stabs::CTMatrixTypes, char_vec::Union{Vector{nmod}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm)

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
            rnk_H, H = right_kernel(hcat(stabs[:, n + 1:end], -stabs[:, 1:n]))
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
            return StabilizerCodeCSS(F, n, dim_code, missing, missing, missing, stabs, args[2],
                args[4], missing, missing, signs, args[3], args[5], logs, logs_mat, char_vec,
                missing, missing, missing, over_comp, missing, stabs_stand, stand_r, stand_k, P_stand)
        else
            return GraphStateStabilizerCSS(F, n, 0, missing, missing, missing, stabs, args[2],
                args[4], missing, missing, signs, args[3], args[5], char_vec, missing, over_comp, stabs_stand,
                stand_r, stand_k, P_stand)
        end
    else
        if !iszero(stand_k)
            return StabilizerCode(F, n, dim_code, missing, stabs, logs, logs_mat, char_vec, signs, missing,
                missing, missing, over_comp, missing, stabs_stand, stand_r, stand_k, P_stand)
        else
            return GraphStateStabilizer(F, n, 0, missing, stabs, char_vec, signs, missing, over_comp, stabs_stand,
                stand_r, stand_k, P_stand)
        end
    end
end

"""
    StabilizerCode(S::AbstractSubsystemCode)

Return the stabilizer code from the subsystem code by promoting all of the gauge operators to logical operators.
"""
function StabilizerCode(S::AbstractSubsystemCode)
    typeof(S) <: AbstractStabilizerCode && return S
    return promotegaugestological(S, collect(1:length(S.gaugeops)))
end
# function StabilizerCodeCSS(S::AbstractSubsystemCodeCSS)
#     # promote all gauges to logicals
# end

#############################
      # getter functions
#############################

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

function _logicals(stabs::T, dual_gens::T, logs_alg::Symbol=:sys_eqs) where {T<:CTMatrixTypes}
    logs_alg ∈ [:sys_eqs, :VS] || throw(ArgumentError("Unrecognized logicals algorithm"))

    L = _quotient_space(dual_gens, stabs, logs_alg)
    logs = _make_pairs(L)
    # verify
    n = div(ncols(L), 2)
    logs_mat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
    are_symplectic_orthogonal(stabs, logs_mat) || error("Computed logicals do not commute with the codespace.")
    prod = hcat(logs_mat[:, n + 1:end], -logs_mat[:, 1:n]) * transpose(logs_mat)
    sum(FpmattoJulia(prod), dims=1) == ones(Int, 1, size(prod, 1)) ||
        error("Computed logicals do not have the right commutation relations.")
    return logs, logs_mat
end

# TODO: switch to n, k inputs
"""
   random_CSS_code(num_logicals::Int, num_stabs::Int)

Return a random CSS code with an equal number of X and Z stabilizers.
"""
function random_CSS_code(num_logicals::Int, num_stabs::Int)
    d = _rand_single_sector_boundary(num_logicals, num_stabs)
    return CSSCode(d, transpose(d))
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
