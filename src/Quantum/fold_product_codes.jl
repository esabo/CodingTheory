# Copyright (c) 2022 - 2026 Eric Sabo, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
$(TYPEDSIGNATURES)

Return the single-parity-check `D`-fold product code.

# Note
- This is defined in https://arxiv.org/abs/2209.13474

# Example

[512, 174, 8]] Symmetric 2-fold product CSS code from [ostrev2024classical](@cite)

```jldoctest
julia> using CodingTheory, Oscar;

julia> F = Oscar.Nemo.Native.GF(2);

julia> h = matrix(F, [1 1]);

julia> id = identity_matrix(F, 2);

julia> H_X = vcat(
             h ⊗ h ⊗ h ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id,
             id ⊗ id ⊗ id ⊗ h ⊗ h ⊗ h ⊗ id ⊗ id ⊗ id,
             id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ h ⊗ h ⊗ h);

julia> H_Z = vcat(
             h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id,
             id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id,
             id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h);

julia> code = SPCDFoldProductCode(3);

julia> length(code), dimension(code)
(512, 174)
```
"""
function SPCDFoldProductCode(D::Int, s::Int = 1)
    D >= 2 || throw(DomainError(D, "The fold count must be at least two."))
    s > 0 || throw(DomainError(s, "The scale must be positive."))
    vec_S = AbstractSubsystemCode[]
    for i in 1:D
        for l in 1:D
            if l == i
                push!(vec_S, CSSCode(RepetitionCode(2, 2 * s)))
            else
                push!(vec_S, CSSCode(RepetitionCode(2, 2)))
            end
        end
    end
    
    S = symmetric_product(vec_S)
    
    # Inject exact distance boundaries manually 
    dist = 2^D
    S.d = dist
    S.l_bound = dist
    S.u_bound = dist
    S.cache[:pure] = true
    
    return S
end
"""
$(TYPEDSIGNATURES)

Return the single-parity-check ``D``-fold product code with scale `s`.
"""
SingleParityCheckDFoldProductCode(D::Int, s::Int = 1) = SPCDFoldProductCode(D, s)

#############################
      # getter functions
#############################

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

# TODO: present the stabilizers in the docs and mention how to switch X and Z by
# using the switch on the inputs beforehand
"""
$(TYPEDSIGNATURES)

Return the asymmetric 2-fold product quantum CSS code of the CSS codes `S1` and `S2`.

# Note
- This is defined in https://arxiv.org/abs/2209.13474
"""
function asymmetric_product(::IsCSS, S1::AbstractSubsystemCode, S2::AbstractSubsystemCode;
    char_vec::Union{Vector{zzModRingElem}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm)
    F = S1.F
    F == S2.F || throw(ArgumentError("Base rings must match."))
    
    # Calculate physical qubits O(1)
    n_new = S1.n * S2.n
    cache = _family_cache(n_new; logs_alg=logs_alg)
    result = AsymmetricProductCode(F, S1, S2, n_new, 0,
        _family_char_vec(char_vec, F, n_new), cache)
    result.k = n_new - rank(X_stabilizers(result)) - rank(Z_stabilizers(result))
    return _seed_quantum_singleton_bound!(result)
end
asymmetric_product(S1::T, S2::T; kwargs...) where {T <: AbstractSubsystemCode} =
    asymmetric_product(CSSTrait(T), S1, S2; kwargs...)
asymmetric_product(::IsNotCSS, S1::AbstractSubsystemCode, S2::AbstractSubsystemCode; kwargs...) =
    error("Only valid for CSS codes.")

"""
$(TYPEDSIGNATURES)

Return the symmetric `D`-fold product quantum CSS code, where `D` is
the square-root of the length of the vector of CSS codes `vec_S`.

# Note
- This is defined in https://arxiv.org/abs/2209.13474
"""
function symmetric_product(::IsCSS, vec_S::Vector{T};
    char_vec::Union{Vector{zzModRingElem}, Missing}=missing,
    logs_alg::Symbol=:stnd_frm) where {T <: AbstractSubsystemCode}
    length(vec_S) >= 4 || throw(DomainError("The length of the input vector must be at least 4"))
    D_float = sqrt(length(vec_S))
    isinteger(D_float) ? (D = Int(D_float)) : throw(ArgumentError("The number of CSS codes must be D^2"))
    
    F = vec_S[1].F
    n_new = 1
    
    # Calculate physical qubits O(1) by multiplying all lengths
    for S in vec_S
        S.F == F || throw(ArgumentError("All codes must share the same base ring."))
        n_new *= S.n
    end
    
    cache = _family_cache(n_new; logs_alg=logs_alg)
    result = SymmetricProductCode(F, vec_S, D, n_new, 0,
        _family_char_vec(char_vec, F, n_new), cache)
    result.k = n_new - rank(X_stabilizers(result)) - rank(Z_stabilizers(result))
    return _seed_quantum_singleton_bound!(result)
end
function symmetric_product(vec_S::Vector{T}; kwargs...) where {T <: AbstractSubsystemCode}
    isempty(vec_S) && throw(ArgumentError("Input vector of CSS codes cannot be empty"))
    for S in vec_S
        if CSSTrait(typeof(S)) == IsNotCSS()
            return symmetric_product(IsNotCSS(), vec_S; kwargs...)
        end
    end
    return symmetric_product(IsCSS(), vec_S; kwargs...)
end
symmetric_product(::IsNotCSS, vec_S::Vector{T}; kwargs...) where {T <: AbstractSubsystemCode} =
    error("Only valid for CSS codes.")

function X_stabilizers(S::AsymmetricProductCode)
    haskey(S.cache, :H_X) && return S.cache[:H_X]
    F = S.F
    H_X = vcat(kronecker_product(X_stabilizers(S.S1), identity_matrix(F, S.S2.n)), 
               kronecker_product(identity_matrix(F, S.S1.n), X_stabilizers(S.S2)))
    S.cache[:H_X] = H_X
    return H_X
end

function Z_stabilizers(S::AsymmetricProductCode)
    haskey(S.cache, :H_Z) && return S.cache[:H_Z]
    H_Z = kronecker_product(Z_stabilizers(S.S1), Z_stabilizers(S.S2))
    S.cache[:H_Z] = H_Z
    return H_Z
end

function stabilizers(S::AsymmetricProductCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    H_X = X_stabilizers(S)
    H_Z = Z_stabilizers(S)
    F = S.F
    stabs = vcat(hcat(H_X, zero_matrix(F, nrows(H_X), S.n)),
                 hcat(zero_matrix(F, nrows(H_Z), S.n), H_Z))
    S.cache[:stabilizers] = stabs
    return stabs
end

# --- Symmetric Product Accessors ---
function X_stabilizers(S::SymmetricProductCode)
    haskey(S.cache, :H_X) && return S.cache[:H_X]
    D = S.D
    F = S.F
    vec_S = S.vec_S
    
    H_X = nothing
    for j in 0:D - 1
        temp_row = nothing
        for l in 1:D^2
            mat = (j * D + 1 <= l <= (j + 1) * D) ? X_stabilizers(vec_S[l]) : identity_matrix(F, vec_S[l].n)
            temp_row = l == 1 ? mat : kronecker_product(temp_row, mat)
        end
        H_X = j == 0 ? temp_row : vcat(H_X, temp_row)
    end
    S.cache[:H_X] = H_X
    return H_X
end

function Z_stabilizers(S::SymmetricProductCode)
    haskey(S.cache, :H_Z) && return S.cache[:H_Z]
    D = S.D
    F = S.F
    vec_S = S.vec_S
    
    H_Z = nothing
    for j in 0:D - 1
        temp_row = nothing
        for l in 1:D^2
            mat = ((l - 1) % D == j) ? Z_stabilizers(vec_S[l]) : identity_matrix(F, vec_S[l].n)
            temp_row = l == 1 ? mat : kronecker_product(temp_row, mat)
        end
        H_Z = j == 0 ? temp_row : vcat(H_Z, temp_row)
    end
    S.cache[:H_Z] = H_Z
    return H_Z
end

function stabilizers(S::SymmetricProductCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    H_X = X_stabilizers(S)
    H_Z = Z_stabilizers(S)
    F = S.F
    stabs = vcat(hcat(H_X, zero_matrix(F, nrows(H_X), S.n)),
                 hcat(zero_matrix(F, nrows(H_Z), S.n), H_Z))
    S.cache[:stabilizers] = stabs
    return stabs
end
