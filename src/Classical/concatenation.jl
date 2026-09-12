# Copyright (c) 2023 - 2026 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.
#############################
        # Classical
#############################

# ==============================================================================
# CONCATENATED CODE LAZY GETTERS
# ==============================================================================

# Internal helper to rigorously compute the nullspace basis 
# using the native no-column-swap RREF algorithm.
function _parity_check_from_G(G::CTMatrixTypes)
    R = _rref_no_col_swap(G)
    nc = ncols(R)
    non_pivots = _rref_non_pivot_cols(R, :nsp)
    pivots = sort!(setdiff(1:nc, non_pivots))
    k = length(pivots)
    
    H = zero_matrix(base_ring(G), nc - k, nc)
    for (idx, np) in enumerate(non_pivots)
        H[idx, np] = 1
        for i in 1:k
            H[idx, pivots[i]] = -R[i, np]
        end
    end
    return H
end

function generator_matrix(C::ConcatenatedCode, stand_form::Bool = false)
    cache = getfield(C, :cache)
    if !haskey(cache, :G)
        G_out = generator_matrix(C.C_out)
        
        # Apply standard form permutation if it exists to align the outer code
        P_out = standard_form_permutation(C.C_out)
        ismissing(P_out) || (G_out = G_out * P_out)
        
        if C.type == :expanded
            G_out = expand_matrix(G_out, C.F, C.basis)
        else
            G_out = change_base_ring(C.F, G_out)
        end
        
        G_in = generator_matrix(C.C_in)
        P_in = standard_form_permutation(C.C_in)
        ismissing(P_in) || (G_in = G_in * P_in)
        
        cache[:G] = _concatenated_generator_matrix(G_out, G_in)
    end
    if stand_form
        if !haskey(cache, :G_stand)
            G_stand, H_stand, P, _ = _standard_form(cache[:G])
            cache[:G_stand] = G_stand
            cache[:H_stand] = H_stand
            cache[:P_stand] = P
        end
        return cache[:G_stand]
    end
    return cache[:G]
end

function parity_check_matrix(C::ConcatenatedCode, stand_form::Bool = false)
    cache = getfield(C, :cache)
    if !haskey(cache, :H)
        G_mat = generator_matrix(C)
        cache[:H] = _parity_check_from_G(G_mat)
    end
    if stand_form
        generator_matrix(C, true)
        return cache[:H_stand]
    end
    return cache[:H]
end

# Internal function for concatenated generator assembly
function _concatenated_generator_matrix(A::T, B::T) where T <: CTMatrixTypes
    nr_A, nc_A = size(A)
    nr_B, nc_B = size(B)
    t = div(nc_A, nr_B)
    M = zero_matrix(base_ring(A), nr_A, t * nc_B)
    for i in 1:t
        M[:, nc_B * (i - 1) + 1: nc_B * i] = view(A, :, nr_B * (i - 1) + 1:nr_B * i) * B
    end
    return M
end

# ==============================================================================
# CONSTRUCTORS
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the single-level concatenation of `C_out` and `C_in`.
Evaluates lazily without eagerly building the generator matrices.
"""
function concatenate(C_out::AbstractLinearCode, C_in::AbstractLinearCode)
    F_out = C_out.F
    F_in = C_in.F
    β, λ = missing, missing
    type = :same

    if Int(order(F_out)) != Int(order(F_in))
        flag, deg = is_extension(F_out, F_in)
        flag || throw(ArgumentError("Galois concatenation requires the outer code to be over an extension field of the inner code"))
        deg % C_in.k == 0 || C_out.n % C_in.k == 0 || throw(ArgumentError("Inner dimension must divide outer length or extension degree"))
        
        β, λ = primitive_basis(F_out, F_in)
        type = :expanded
    else
        C_out.n % C_in.k == 0 || throw(ArgumentError("Inner dimension must divide outer length"))
        type = :same
    end
    
    n_new = C_in.n * div(C_out.n, C_in.k)
    k_new = C_out.k
    
    # O(1) Distance bounds calculation
    if ismissing(C_out.d) || ismissing(C_in.d)
        d_new = missing
        lb = C_out.l_bound * C_in.l_bound
    else
        d_new = C_out.d * C_in.d
        lb = d_new
    end
    
    cache = Dict{Symbol, Any}()
    # Field order in struct is (C_in, C_out, ...)
    return ConcatenatedCode(C_in, C_out, type, β, λ, F_in, n_new, k_new, d_new, lb, n_new, cache)
end
∘(C_out::AbstractLinearCode, C_in::AbstractLinearCode) = concatenate(C_out, C_in)

# ==============================================================================
# MULTILEVEL CONCATENATION LAZY GETTERS
# ==============================================================================

function generator_matrix(C::MultilevelConcatenatedCode, stand_form::Bool = false)
    cache = getfield(C, :cache)
    if !haskey(cache, :G)
        F = C.F
        n_in = C.C_ins[1].n
        n_out = div(C.C_outs[1].n, C.C_ins[1].k)
        
        # 1. Expand and collect all outer generator matrices
        G_outs_expanded = []
        for i in 1:length(C.C_outs)
            G_curr = generator_matrix(C.C_outs[i])
            P_curr = standard_form_permutation(C.C_outs[i])
            ismissing(P_curr) || (G_curr = G_curr * P_curr)
            
            if C.types[i] == :expanded
                G_curr = expand_matrix(G_curr, F, C.bases[i])
            else
                G_curr = change_base_ring(F, G_curr)
            end
            push!(G_outs_expanded, G_curr)
        end
        G1 = reduce(direct_sum, G_outs_expanded)
        
        # 2. Collect inner block quotient spaces B_i = G_i / G_{i-1}
        B = [generator_matrix(C.C_ins[1])]
        for i in 2:length(C.C_ins)
            Gi = generator_matrix(C.C_ins[i])
            Gim1 = generator_matrix(C.C_ins[i - 1])
            push!(B, _quotient_space(Gim1, Gi))
        end
        
        # 3. Assemble the block matrix G2
        G2 = zero_matrix(F, ncols(G1), n_in * n_out)
        z = 1
        for i in 1:length(C.C_ins)
            rows_B = nrows(B[i])
            for j in 0:(n_out - 1)
                r_start = z
                r_end = z + rows_B - 1
                c_start = j * n_in + 1
                c_end = (j + 1) * n_in
                
                # Direct assignment for speed
                G2[r_start:r_end, c_start:c_end] = B[i]
                z += rows_B
            end
        end
        
        # 4. The final concatenated generator matrix
        cache[:G] = G1 * G2
    end
    
    if stand_form
        if !haskey(cache, :G_stand)
            G_stand, H_stand, P, _ = _standard_form(cache[:G])
            cache[:G_stand] = G_stand
            cache[:H_stand] = H_stand
            cache[:P_stand] = P
        end
        return cache[:G_stand]
    end
    return cache[:G]
end

function parity_check_matrix(C::MultilevelConcatenatedCode, stand_form::Bool = false)
    cache = getfield(C, :cache)
    if !haskey(cache, :H)
        G_mat = generator_matrix(C)
        cache[:H] = _parity_check_from_G(G_mat)
    end
    if stand_form
        generator_matrix(C, true)
        return cache[:H_stand]
    end
    return cache[:H]
end

# ==============================================================================
# MULTILEVEL CONSTRUCTORS
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the generalized concatenation of a list of outer codes and a nested list of inner codes.
Evaluates lazily without building the enormous generator matrices.
"""
function concatenate(outers_unexpanded::Vector{T}, inners::Vector{T}) where T <: AbstractLinearCode
    isempty(outers_unexpanded) && throw(ArgumentError("List of codes cannot be empty"))
    length(outers_unexpanded) == length(inners) || throw(ArgumentError("Must have the same number of inner and outer codes"))
    
    for i in 2:length(inners)
        inners[i - 1] ⊆ inners[i] || throw(ArgumentError("The inner subcodes must be in a nested sequence (C_{i-1} ⊆ C_i)"))
    end
    
    F = first(inners).F
    n_in = first(inners).n
    ord_F = Int(order(F))

    # O(1) dimension and parameter validation
    n_out = divexact(outers_unexpanded[1].n, inners[1].k)
    for i in 2:length(outers_unexpanded)
        n_out == divexact(outers_unexpanded[i].n, inners[i].k - inners[i - 1].k) || throw(ArgumentError("The outer matrices are not of the correct size"))
    end

    bases = Union{Vector{<:CTFieldElem}, Missing}[missing for _ in eachindex(outers_unexpanded)]
    dual_bases = Union{Vector{<:CTFieldElem}, Missing}[missing for _ in eachindex(outers_unexpanded)]
    types = Symbol[:same for _ in eachindex(outers_unexpanded)]
    
    for (i, C_out) in enumerate(outers_unexpanded)
        # If the fields are the same, NO expansion is needed. Just map it natively.
        if Int(order(C_out.F)) == ord_F
            types[i] = :same
        else
            flag, _ = is_subfield(F, C_out.F)
            flag || throw(ArgumentError("Cannot connect outer code $i field to inner code field"))
            bases[i], dual_bases[i] = primitive_basis(C_out.F, F)
            types[i] = :expanded
        end
    end

    n_new = n_in * n_out
    k_new = sum(C.k for C in outers_unexpanded)
    
    # Distance bounds calculation
    d_new = missing
    lb = 1
    if all(isequal(:expanded), types) && all(!ismissing(inners[i].d) && !ismissing(outers_unexpanded[i].d) for i in eachindex(inners))
        d_new = minimum(inners[i].d * outers_unexpanded[i].d for i in eachindex(inners))
        lb = d_new
    else
        lb = minimum(inners[i].l_bound * outers_unexpanded[i].l_bound for i in eachindex(inners))
    end
    
    cache = Dict{Symbol, Any}()
    # FIX: Correctly maps (C_outs, C_ins) without swapping!
    return MultilevelConcatenatedCode(outers_unexpanded, inners, types, bases, dual_bases, F, n_new, k_new, d_new, lb, n_new, cache)
end
"""
$(TYPEDSIGNATURES)

Return the generalized concatenation of `outers` with the nested sequence
`inners`. This is an alias for `concatenate(outers, inners)`; the result has
length ``n_{\\mathrm{in}} n_{\\mathrm{out}}`` and dimension equal to the sum
of the outer-code dimensions.
"""
multilevel_concatenation(outers::Vector{T}, inners::Vector{T}) where T <: AbstractLinearCode = concatenate(outers, inners)

# ==============================================================================
# GETTERS & UTILITIES
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the inner code(s) of the concatenation.
"""
inner_code(C::ConcatenatedCode) = C.C_in
inner_code(C::MultilevelConcatenatedCode) = C.C_ins

"""
$(TYPEDSIGNATURES)

Return the outer code(s) of the concatenation.
"""
outer_code(C::ConcatenatedCode) = C.C_out
outer_code(C::MultilevelConcatenatedCode) = C.C_outs

"""
$(TYPEDSIGNATURES)

Return the basis (or list of bases) used to expand the outer code(s), if applicable.
"""
expansion_basis(C::ConcatenatedCode) = C.basis
expansion_basis(C::MultilevelConcatenatedCode) = C.bases

"""
$(TYPEDSIGNATURES)

Return the dual basis (or list of dual bases) used to expand the outer code(s), if applicable.
"""
expansion_dual_basis(C::ConcatenatedCode) = C.dual_basis
expansion_dual_basis(C::MultilevelConcatenatedCode) = C.dual_bases

"""
$(TYPEDSIGNATURES)

Return the type(s) of concatenation used (`:same` or `:expanded`).
"""
concatenation_type(C::ConcatenatedCode) = C.type
concatenation_type(C::MultilevelConcatenatedCode) = C.types

# ==============================================================================
# ENCODING
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the encoding of `v` into `C`, where `v` is either a valid input for the outer code or the full code.
"""
function encode(C::ConcatenatedCode, v::Union{CTMatrixTypes, Vector{Int}})
    w = isa(v, Vector{Int}) ? matrix(C.C_out.F, 1, length(v), v) : v
    
    nr_w, nc_w = size(w)
    (nr_w != 1 && nc_w == 1) && (w = transpose(w))
    nc_w = ncols(w)
    
    if nc_w == C.C_out.k
        base_ring(w) == C.C_out.F || throw(ArgumentError("Vector must have the same base ring as the outer code."))
        
        G_out = generator_matrix(C.C_out)
        c_out = w * G_out
        
        if C.type == :expanded
            D = _expansion_dict(C.C_out.F, C.C_in.F, C.dual_basis)
            c_out = _expand_matrix(c_out, D, div(degree(C.C_out.F), degree(C.C_in.F)))
        else
            c_out = change_base_ring(C.C_in.F, c_out)
        end
        
        Gin = generator_matrix(C.C_in)
        
        k_in = C.C_in.k
        n_in = C.C_in.n
        t = div(ncols(c_out), k_in)
        
        c_final = zero_matrix(C.C_in.F, 1, t * n_in)
        for i in 1:t
            # Using matrix() to realize the view into a concrete fpMatrix to prevent SubMat multiplication issues
            block = matrix(C.C_in.F, 1, k_in, [c_out[1, (i-1)*k_in + c] for c in 1:k_in])
            c_final[1:1, (i-1)*n_in + 1 : i*n_in] = block * Gin
        end
        
        return c_final
        
    elseif nc_w == C.k
        base_ring(w) == C.F || throw(ArgumentError("Vector must have the same base ring as the code."))
        G = generator_matrix(C)
        return w * G
    else
        throw(ArgumentError("Vector has incorrect dimension; expected $(C.C_out.k) or $(C.k)."))
    end
end
