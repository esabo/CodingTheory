# Copyright (c) 2021 - 2026 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    StabilizerCodeCSS(X_matrix::CTMatrixTypes, Z_matrix::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)
    CSSCode(X_matrix::CTMatrixTypes, Z_matrix::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing}= missing, logs_alg::Symbol = :stnd_frm)

Return a CSS code whose `X`-stabilizers are given by `X_matrix`, `Z`-stabilizers by `Z_matrix`.
"""
function StabilizerCodeCSS(X_matrix::T, Z_matrix::T; char_vec::Union{Vector{zzModRingElem}, Missing} = 
    missing, logs_alg::Symbol = :stnd_frm) where T <: CTMatrixTypes

    logs_alg ∈ (:stnd_frm, :sys_eqs) || throw(ArgumentError("Unrecognized logicals algorithm. Use :stnd_frm or :sys_eqs."))
    iszero(X_matrix) && throw(ArgumentError("The `X` stabilizer matrix is empty."))
    iszero(Z_matrix) && throw(ArgumentError("The `Z` stabilizer matrix is empty."))
    
    n = ncols(X_matrix)
    n == ncols(Z_matrix) || throw(ArgumentError("Both matrices must have the same length in the CSS construction."))
    F = _code_matrix_base_ring(X_matrix)
    F == _code_matrix_base_ring(Z_matrix) || throw(ArgumentError("Both matrices must be over the same base field."))
    
    is_sparse = X_matrix isa SparseMatrixCSC
    X_dense = _dense_code_matrix(X_matrix, F)
    Z_dense = _dense_code_matrix(Z_matrix, F)
    
    iszero(Z_dense * transpose(X_dense)) || throw(ArgumentError("The given matrices are not symplectic orthogonal."))
    
    p = Int(characteristic(F))
    clean_char_vec = ismissing(char_vec) ? zzModRingElem[] : _process_char_vec(char_vec, p, 2 * n)

    X_dense = _remove_empty(X_dense, :rows)
    Z_dense = _remove_empty(Z_dense, :rows)

    X_rank = rank(X_dense)
    Z_rank = rank(Z_dense)
    rnk = X_rank + Z_rank
    
    dim_code = BigInt(order(F))^n // BigInt(p)^rnk
    isinteger(dim_code) && (dim_code = round(Int, log(BigInt(p), dim_code));)
    
    over_comp = (nrows(X_dense) > X_rank) || (nrows(Z_dense) > Z_rank)
    
    X_final = is_sparse ? _remove_empty(deepcopy(X_matrix), :rows) : X_dense
    Z_final = is_sparse ? _remove_empty(deepcopy(Z_matrix), :rows) : Z_dense

    cache = Dict{Symbol, Any}(
        :overcomplete => over_comp,
        :logs_alg => logs_alg
    )

    return StabilizerCodeCSS(F, n, dim_code, X_final, Z_final, clean_char_vec, cache)
end
CSSCode(X_matrix::T, Z_matrix::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm) where T <: CTMatrixTypes = StabilizerCodeCSS(X_matrix, Z_matrix,
    char_vec = char_vec, logs_alg = logs_alg)

"""
    StabilizerCode(stabs::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

Return the stabilizer code whose stabilizers is determined by `stabs`.
"""
function StabilizerCode(stabs::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm)

    logs_alg ∈ (:stnd_frm, :sys_eqs) || throw(ArgumentError("Unrecognized logicals algorithm. Use :stnd_frm or :sys_eqs."))
    iszero(stabs) && throw(ArgumentError("The stabilizer matrix is empty."))
    
    F = base_ring(stabs)
    p = Int(characteristic(F))
    n = div(ncols(stabs), 2)
    clean_char_vec = ismissing(char_vec) ? zzModRingElem[] : _process_char_vec(char_vec, p, 2 * n)
    
    is_sparse = stabs isa SparseMatrixCSC
    stabs_dense = is_sparse ? matrix(F, stabs) : stabs
    stabs_dense = _remove_empty(stabs_dense, :rows)
    are_symplectic_orthogonal(stabs_dense, stabs_dense) || throw(ArgumentError("The given stabilizers are not symplectic orthogonal."))
    
    rnk = rank(stabs_dense)
    dim_code = BigInt(order(F))^n // BigInt(p)^rnk
    isinteger(dim_code) && (dim_code = round(Int, log(BigInt(p), dim_code));)
    over_comp = nrows(stabs_dense) > rnk
    
    is_css_S, X_stabs_dense, Z_stabs_dense = robust_CSS_split(stabs_dense)

    stabs_final = is_sparse ? sparse(stabs_dense) : stabs_dense

    cache = Dict{Symbol, Any}(
        :stabs => stabs_final,
        :overcomplete => over_comp,
        :logs_alg => logs_alg
    )
    
    if is_css_S
        X_stabs = is_sparse ? sparse(X_stabs_dense) : X_stabs_dense
        Z_stabs = is_sparse ? sparse(Z_stabs_dense) : Z_stabs_dense
        return StabilizerCodeCSS(F, n, dim_code, X_stabs, Z_stabs, clean_char_vec, cache)
    else
        return StabilizerCode(F, n, dim_code, stabs_final, clean_char_vec, cache)
    end
end

"""
    StabilizerCodeCSS(C1::AbstractLinearCode, C2::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)
    CSSCode(C1::AbstractLinearCode, C2::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)
"""
function StabilizerCodeCSS(C1::AbstractLinearCode, C2::AbstractLinearCode;
    char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

    C2 ⊆ C1 || throw(ArgumentError("The second argument must be a subset of the first in the CSS construction."))
    D2 = dual(C2)
    
    S = StabilizerCodeCSS(D2.H, C1.H, char_vec=char_vec, logs_alg=logs_alg)
    
    S.cache[:X_orig_code] = D2
    S.cache[:Z_orig_code] = C1
    
    if !ismissing(C1.d) && !ismissing(D2.d)
        d = min(C1.d, D2.d)
        S.cache[:l_bound_dx] = D2.d
        S.cache[:l_bound_dz] = C1.d
        S.cache[:l_bound] = d
    end
    return S
end
CSSCode(C1::AbstractLinearCode, C2::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) = StabilizerCodeCSS(C1, C2, char_vec = char_vec, logs_alg = logs_alg)

"""
    StabilizerCodeCSS(C::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)
    CSSCode(C::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)
"""
function StabilizerCodeCSS(C::LinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm)

    D = dual(C)
    C ⊆ D || throw(ArgumentError("The single code CSS construction requires C ⊆ C^⟂."))
    
    S = StabilizerCodeCSS(D.H, D.H, char_vec=char_vec, logs_alg=logs_alg)
    S.cache[:X_orig_code] = D
    S.cache[:Z_orig_code] = D
    
    if !ismissing(D.d)
        S.cache[:l_bound_dx] = D.d
        S.cache[:l_bound_dz] = D.d
        S.cache[:l_bound] = D.d
    end
    return S
end
CSSCode(C::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) = StabilizerCodeCSS(C, char_vec = char_vec, logs_alg = logs_alg)

"""
    StabilizerCodeCSS(S_Pauli::Vector{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{String, Vector{Char}}
"""
function StabilizerCodeCSS(S_Pauli::Vector{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm) where T <: Union{String, Vector{Char}}

    stabs = _Pauli_string_to_symplectic(_process_strings(S_Pauli))
    iszero(stabs) && throw(ArgumentError("The processed Pauli strings returned a set of empty stabilizer generators."))
    stabs = _remove_empty(stabs, :rows)
    are_symplectic_orthogonal(stabs, stabs) || throw(ArgumentError("The given stabilizers are not symplectic orthogonal."))
    
    return StabilizerCode(stabs, char_vec = char_vec, logs_alg = logs_alg)
end
CSSCode(S_Pauli::Vector{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{String, Vector{Char}} = StabilizerCodeCSS(S_Pauli, char_vec = char_vec, logs_alg = logs_alg)

"""
    StabilizerCode(S_Pauli::Vector{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{String, Vector{Char}}
"""
function StabilizerCode(S_Pauli::Vector{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm) where T <: Union{String, Vector{Char}}

    S_Pauli_stripped = _process_strings(S_Pauli)
    stabs = _Pauli_string_to_symplectic(S_Pauli_stripped)
    iszero(stabs) && throw(ArgumentError("The processed Pauli strings returned a set of empty stabilizer generators."))
    return StabilizerCode(stabs, char_vec = char_vec, logs_alg = logs_alg)
end

"""
    StabilizerCodeCSS(S::AbstractStabilizerCode; logs_alg::Symbol = :stnd_frm)
"""
function StabilizerCodeCSS(S::AbstractStabilizerCode; logs_alg::Symbol = :stnd_frm)
    Z = deepcopy(stabilizers(S))
    Z[:, 1:S.n] = stabilizers(S)[:, S.n + 1:2 * S.n]
    Z[:, S.n + 1:2 * S.n] = -stabilizers(S)[:, 1:S.n]
    return CSSCode(stabilizers(S), Z, logs_alg = logs_alg)
end
CSSCode(S::AbstractStabilizerCode; logs_alg::Symbol = :stnd_frm) = StabilizerCodeCSS(S; logs_alg = logs_alg)

"""
    StabilizerCode(S::AbstractStabilizerCode; logs_alg::Symbol = :stnd_frm)
"""
StabilizerCode(S::T; logs_alg::Symbol = :stnd_frm) where {T <: AbstractStabilizerCode} = StabilizerCode(CSSTrait(T), S, logs_alg)
function StabilizerCode(::IsCSS, S::AbstractStabilizerCode, logs_alg::Symbol)
    iseven(S.n) || throw(ArgumentError("Only valid for codes of even length"))
    nrows(S.X_stabs) == nrows(S.Z_stabs) || throw(ArgumentError("Requires an equal number of X and Z stabilizers"))
    n = div(S.n, 2)
    return StabilizerCode(hcat(S.X_stabs[:, 1:n], S.Z_stabs[:, 1:n]), logs_alg = logs_alg)
end
StabilizerCode(::IsNotCSS, S::AbstractStabilizerCode, logs_alg::Symbol) = throw(ArgumentError("Only valid for CSS codes of even length"))

"""
    StabilizerCode(S::AbstractSubsystemCode)
"""
function StabilizerCode(S::AbstractSubsystemCode)
    typeof(S) <: AbstractStabilizerCode && return S
    # For a subsystem code, promoting all gauges to logicals yields the underlying stabilizer code
    return promote_gauges_to_logical(S, collect(1:length(gauges(S))))
end

#############################
      # getter functions
#############################

"""
    minimum_distance_lower_bound(S::AbstractStabilizerCode)
Return the currently stored lower bound on the minimum distance.
"""
minimum_distance_lower_bound(S::AbstractStabilizerCode) = get(S.cache, :l_bound, missing)

"""
    minimum_distance_upper_bound(S::AbstractStabilizerCode)
Return the currently stored upper bound on the minimum distance.
"""
minimum_distance_upper_bound(S::AbstractStabilizerCode) = get(S.cache, :u_bound, missing)

"""
    X_minimum_distance_lower_bound(S::AbstractStabilizerCode)
Return the currently stored lower bound on the minimum `X`-distance.
"""
X_minimum_distance_lower_bound(S::T) where T <: AbstractStabilizerCode = X_minimum_distance_lower_bound(CSSTrait(T), S)
X_minimum_distance_lower_bound(::IsCSS, S::AbstractStabilizerCode) = get(S.cache, :l_bound_dx, missing)
X_minimum_distance_lower_bound(::IsNotCSS, S::AbstractStabilizerCode) = error("Only valid for CSS codes")

"""
    X_minimum_distance_upper_bound(S::AbstractStabilizerCode)
Return the currently stored upper bound on the minimum `X`-distance.
"""
X_minimum_distance_upper_bound(S::T) where T <: AbstractStabilizerCode = X_minimum_distance_upper_bound(CSSTrait(T), S)
X_minimum_distance_upper_bound(::IsCSS, S::AbstractStabilizerCode) = get(S.cache, :u_bound_dx, missing)
X_minimum_distance_upper_bound(::IsNotCSS, S::AbstractStabilizerCode) = error("Only valid for CSS codes")

"""
    Z_minimum_distance_lower_bound(S::AbstractStabilizerCode)
Return the currently stored lower bound on the minimum `Z`-distance.
"""
Z_minimum_distance_lower_bound(S::T) where T <: AbstractStabilizerCode = Z_minimum_distance_lower_bound(CSSTrait(T), S)
Z_minimum_distance_lower_bound(::IsCSS, S::AbstractStabilizerCode) = get(S.cache, :l_bound_dz, missing)
Z_minimum_distance_lower_bound(::IsNotCSS, S::AbstractStabilizerCode) = error("Only valid for CSS codes")

"""
    Z_minimum_distance_upper_bound(S::AbstractStabilizerCode)
Return the currently stored upper bound on the minimum `Z`-distance.
"""
Z_minimum_distance_upper_bound(S::T) where T <: AbstractStabilizerCode = Z_minimum_distance_upper_bound(CSSTrait(T), S)
Z_minimum_distance_upper_bound(::IsCSS, S::AbstractStabilizerCode) = get(S.cache, :u_bound_dz, missing)
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
    
    curr_u_bound = get(S.cache, :u_bound, S.n)
    curr_l_bound = get(S.cache, :l_bound, 1)
    
    curr_u_bound < d && (@warn "The distance set is greater than the current upper bound of $curr_u_bound")
    d < curr_l_bound && (@warn "The distance set is less than the current lower bound of $curr_l_bound")
    
    S.cache[:d] = d
    S.cache[:l_bound] = d
    S.cache[:u_bound] = d

    if CSSTrait(typeof(S)) == IsCSS()
        !haskey(S.cache, :dx) && (S.cache[:l_bound_dx] = d)
        !haskey(S.cache, :dz) && (S.cache[:l_bound_dz] = d)
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
    
    curr_u_bound_dx = get(S.cache, :u_bound_dx, S.n)
    curr_l_bound_dx = get(S.cache, :l_bound_dx, 1)
    
    curr_u_bound_dx < d && (@warn "The distance set is greater than the current upper bound of $curr_u_bound_dx")
    d < curr_l_bound_dx && (@warn "The distance set is less than the current lower bound of $curr_l_bound_dx")
    
    S.cache[:dx] = d
    S.cache[:l_bound_dx] = d
    S.cache[:u_bound_dx] = d

    if haskey(S.cache, :dz)
        S.cache[:d] = min(d, S.cache[:dz])
        S.cache[:u_bound] = S.cache[:d]
        S.cache[:l_bound] = S.cache[:d]
    else
        curr_u_bound = get(S.cache, :u_bound, S.n)
        d < curr_u_bound && (S.cache[:u_bound] = d)
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
    
    curr_u_bound_dz = get(S.cache, :u_bound_dz, S.n)
    curr_l_bound_dz = get(S.cache, :l_bound_dz, 1)
    
    curr_u_bound_dz < d && (@warn "The distance set is greater than the current upper bound of $curr_u_bound_dz")
    d < curr_l_bound_dz && (@warn "The distance set is less than the current lower bound of $curr_l_bound_dz")
    
    S.cache[:dz] = d
    S.cache[:l_bound_dz] = d
    S.cache[:u_bound_dz] = d

    if haskey(S.cache, :dx)
        S.cache[:d] = min(S.cache[:dx], d)
        S.cache[:u_bound] = S.cache[:d]
        S.cache[:l_bound] = S.cache[:d]
    else
        curr_u_bound = get(S.cache, :u_bound, S.n)
        d < curr_u_bound && (S.cache[:u_bound] = d)
    end
    return nothing
end
set_Z_minimum_distance!(::IsNotCSS, S::AbstractStabilizerCode, d::Int) =
    error("Only valid for CSS codes")


#############################
     # general functions
#############################

function _logicals(stabs::T, dual_gens::T, logs_alg::Symbol = :sys_eqs) where {T <: CTMatrixTypes}
    logs_alg ∈ (:sys_eqs,) || throw(ArgumentError("Unrecognized logicals algorithm. Use :sys_eqs."))

    L = _quotient_space(dual_gens, stabs, logs_alg)
    logs = _make_pairs(L)
    # verify
    n = div(ncols(L), 2)
    logs_mat = reduce(vcat, [reduce(vcat, logs[i]) for i in 1:length(logs)])
    are_symplectic_orthogonal(stabs, logs_mat) || error("Computed logicals do not commute with the codespace.")
    prod = hcat(logs_mat[:, n + 1:end], -logs_mat[:, 1:n]) * transpose(logs_mat)
    sum(_Flint_matrix_to_Julia_int_matrix(prod), dims = 1) == ones(Int, 1, size(prod, 1)) ||
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
    C_X = LinearCode(X_stabilizers(S)) # C_2
    C_Z = LinearCode(Z_stabilizers(S)) # dual(C_1)
    
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
        Threads.@threads for m in 1:num_thrds
            c = deepcopy(z)
            for u in GrayCode(k, r, Int[], mutate = true)
                if flag[]
                    LinearAlgebra.mul!(c, G, u)
                    wt_c = 0
                    supp_c = Int[]
                    @inbounds for j in 1:n
                        c[j] % p != 0 && (wt_c += 1; append!(supp_c, j);)
                    end

                    if !iseven(wt_c)
                        Threads.atomic_cas!(flag, true, false)
                        verbose && println("not even weight")
                        break
                    end
                    verbose && println("shortening on ", setdiff(1:n, supp_c))
                    
                    temp = shorten(C_Z, setdiff(1:n, supp_c))
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

"""
    is_triorthogonal(S::AbstractStabilizerCodeCSS ; verbose::Bool=false)

Return `true` if the CSS code `S` is triorthogonal.
Cached to avoid recomputing on subsequent checks.
"""
function is_triorthogonal(S::AbstractStabilizerCodeCSS; verbose::Bool=false)
    haskey(S.cache, :is_triorthogonal) && return S.cache[:is_triorthogonal]
    
    result = is_triorthogonal(X_stabilizers(S), verbose)
    S.cache[:is_triorthogonal] = result
    return result
end
