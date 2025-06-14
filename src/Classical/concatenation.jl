# Copyright (c) 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
# Classical
#############################

#############################
# constructors
#############################

# TODO: give control over expansion basis
"""
    ∘(C_out::AbstractLinearCode, C_in::AbstractLinearCode)
    concatenate(C_out::AbstractLinearCode, C_in::AbstractLinearCode)

Return the concatenation of `C_out` and `C_in`.
"""
function concatenate(C_out::AbstractLinearCode, C_in::AbstractLinearCode)
    F_out = C_out.F
    F_in = C_in.F
    β, λ = missing, missing
    if Int(order(F_out)) != Int(order(F_in))
        flag, deg = is_extension(F_out, F_in)
        flag || throw(
            ArgumentError(
                "Galois concatenation requires the outer code to be over an extension field of the inner code",
            ),
        )
        deg % C_in.k == 0 ||
            C_out.n % C_in.k == 0 ||
            throw(
                ArgumentError(
                    "Inner dimension must divide outer length or extension degree",
                ),
            )
        G_out = generator_matrix(C_out, true)
        ismissing(C_out.P_stand) || (G_out = G_out * C_out.P_stand)

        β, λ = primitive_basis(F_out, F_in)
        D = _expansion_dict(F_out, F_in, λ)
        G_out = _expand_matrix(G_out, D, deg)
        type = :expanded
    else
        C_out.n % C_in.k == 0 ||
            throw(ArgumentError("Inner dimension must divide outer length"))

        F_out == F_in || (C_out = change_field(C_out, F_in);)
        G_out = generator_matrix(C_out, true)
        ismissing(C_out.P_stand) || (G_out = G_out * C_out.P_stand)
        type = :same
    end

    Gin = generator_matrix(C_in, true)
    ismissing(C_in.P_stand) || (Gin = Gin * C_in.P_stand)
    G = _concatenated_generator_matrix(G_out, Gin)
    G_stand, H_stand, P, k = _standard_form(G)
    H = ismissing(P) ? H_stand : H_stand * P
    ub1, _ = _min_wt_row(G)
    ub2, _ = _min_wt_row(G_stand)
    ub = min(ub1, ub2)

    C = ConcatenatedCode(
        C_out,
        C_in,
        type,
        β,
        λ,
        F_in,
        ncols(G),
        k,
        missing,
        1,
        ub,
        G,
        H,
        G_stand,
        H_stand,
        P,
        missing,
    )
    # TODO: distance check here
    # for a lower bound on the distance, count the number of pieces of G_out in _concatenated_generator_matrix
    # with full rank and multiply by the distance of the inner code

    return C
end
∘(C_out::AbstractLinearCode, C_in::AbstractLinearCode) = concatenate(C_out, C_in)

function concatenate(
    outers_unexpanded::Vector{T},
    inners::Vector{T},
) where {T<:AbstractLinearCode}
    isempty(outers_unexpanded) && throw(ArgumentError("List of codes cannot be empty"))
    length(outers_unexpanded) == length(inners) ||
        throw(ArgumentError("Must have the same number of inner and outer codes"))
    for i = 2:length(inners)
        inners[i-1] ⊆ inners[i] || throw(
            ArgumentError("The inner subcodes must be in a decreasing nested sequence"),
        )
    end
    F = first(inners).F
    n_in = first(inners).n

    outers = copy(outers_unexpanded)
    β = Union{Vector{<:CTFieldElem},Missing}[missing for _ in eachindex(outers)]
    λ = Union{Vector{<:CTFieldElem},Missing}[missing for _ in eachindex(outers)]
    type = [:same for i in eachindex(outers)]
    ord_F = Int(order(F))
    for (i, C_out) in enumerate(outers)
        if Int(order(C_out.F)) == ord_F
            # it was either pre-expanded or just doesn't need expansion
            C_out.F != F || (outers[i] = change_field(C_out, F))

            # if it could have been "expanded" without changing
            # anything, list as expanded so that the distance
            # calculation at the end knows
            if (i == 1 && inners[i].k == 1) || (i > 1 && inners[i].k - inners[i-1].k == 1)
                type[i] = :expanded
            end
        elseif is_subfield(F, C_out.F)[1]
            β[i], λ[i] = primitive_basis(outers[i].F, F)
            outers[i] = expanded_code(outers[i], F, β[i])
            type[i] = :expanded
        else
            throw(ArgumentError("Cannot connect outer code $i field to inner code field"))
        end
    end

    # Are the outer matrices the right size?
    n_out = divexact(outers[1].n, inners[1].k)
    for i = 2:length(outers)
        n_out == divexact(outers[i].n, inners[i].k - inners[i-1].k) ||
            throw(ArgumentError("The outer matrices are not of the correct size"))
    end

    B = [generator_matrix(inners[1])]
    for i = 2:length(inners)
        Gi = generator_matrix(inners[i])
        Gim1 = generator_matrix(inners[i-1])
        push!(B, _quotient_space(Gi, Gim1, :VS))
    end

    G1 = reduce(direct_sum, generator_matrix(C) for C in outers)
    G2 = zero_matrix(F, ncols(G1), n_in * n_out)
    z = 1
    for i in eachindex(inners)
        for j = 0:(n_out-1)
            rows = range(z, z + size(B[i], 1) - 1)
            cols = range(j * n_in + 1, (j + 1) * n_in)
            G2[rows, cols] = B[i]
            z += nrows(B[i])
        end
    end

    G = G1 * G2
    G_stand, H_stand, P, k = _standard_form(G)
    H = ismissing(P) ? H_stand : H_stand * P

    d = if all(isequal(:expanded), type)
        # if any of these distances are missing, it correctly results in missing
        reduce(min, inners[i].d * outers_unexpanded[i].d for i in eachindex(inners))
    else
        missing
    end
    lb = ismissing(d) ? 1 : d
    ub1, _ = _min_wt_row(G)
    ub2, _ = _min_wt_row(G_stand)
    ub = ismissing(d) ? min(ub1, ub2) : d

    return ConcatenatedCode(
        outers_unexpanded,
        inners,
        type,
        β,
        λ,
        F,
        ncols(G),
        k,
        d,
        lb,
        ub,
        G,
        H,
        G_stand,
        H_stand,
        P,
        missing,
    )
end
multilevel_concatenation(
    outers::Vector{T},
    inners::Vector{T},
) where {T<:AbstractLinearCode} = concatenate(outers, inners)
# cascade?

# Eric had written this so it's in the way explained by the book but technically it's equivalent to the above
# so it was never finished
# function generalized_concatenation(outers::Vector{T}, inners::Vector{T}) where T <: AbstractLinearCode
#     isempty(outers) || isempty(inners) && throw(ArgumentError("List of codes cannot be empty"))
#     for i in 1:length(inners) - 1
#         inners[i + 1] ⊆ inners[i] || throw(ArgumentError("The inner subcodes must be in a decreasing nested sequence"))
#     end
#     F = first(inners).F
#     n_in = first(inners).n
#     iszero(generator_matrix(inners[end])) || push!(inners, ZeroCode(F, n_in))

#     # ord_F = Int(order(F))
#     # for (i, C_out) in enumerate(outers)
#     #     if Int(order(C_out.F)) == ord_F
#     #         C_out.F != F || (outers[i] = change_field(C_out, F);)
#     #     elseif is_subfield(F, C_out.F)
#     #         # TODO: expansion step here
#     #     else
#     #         throw(ArgumentError("Cannot connect outer code $i field to inner code field"))
#     #     end
#     # end

#     G_in_part = matrix(F, 0, n_in, [])
#     H_in_part = matrix(F, 0, n_in, [])
#     G_part_locs = Vector{Vector{Int}}()
#     H_part_locs = Vector{Vector{Int}}()
#     for i in 1:length(inners) - 1
#         if i != length(inners) - 1
#             Gi = generator_matrix(inners[i])
#             Gip1 = generator_matrix(inners[i + 1])
#             new_rows = _quotient_space(Gi, Gip1, :VS)
#             G_in_part = vcat(G_in_part, new_rows)
#             isempty(G_part_locs) ? (push!(G_part_locs, [1, nrows(new_rows)]);) : (push!(G_part_locs, [G_part_locs[end][2] + 1, G_part_locs[end][2] + nrows(new_rows)]);)
#             # println("G")
#             # display(G_in_part)
#             # println(" ")
#         end
#         Hi = parity_check_matrix(inners[i])
#         Hip1 = parity_check_matrix(inners[i + 1])
#         new_rows = _quotient_space(Hip1, Hi, :VS)
#         H_in_part = vcat(H_in_part, new_rows)
#         isempty(H_part_locs) ? (push!(H_part_locs, [1, nrows(new_rows)]);) : (push!(H_part_locs, [H_part_locs[end][2] + 1, H_part_locs[end][2] + nrows(new_rows)]);)
#         # println("H")
#         # display(H_in_part)
#         # println(" ")
#     end

#     # TODO: finish
#     # now to check to make sure outer code dimensions are equal to H_part_locs length (label sizes)
#     for C_out in outers

#     end

#     return G_in_part, H_in_part, G_part_locs, H_part_locs
# end
# Blokh_Zyablov_concatenation(outers::Vector{T}, inners::Vector{T}) where T <: AbstractLinearCode = generalized_concatenation(outers, inners)

#############################
# getter functions
#############################

"""
    inner_code(C::AbstractConcatenatedCode)

Return the inner code of the concatenation.
"""
inner_code(C::AbstractConcatenatedCode) = C.C_in

"""
    outer_code(C::AbstractConcatenatedCode)

Return the outer code of the concatenation.
"""
outer_code(C::AbstractConcatenatedCode) = C.C_out

"""
    expansion_basis(C::AbstractConcatenatedCode)

Return the basis used to expanded the outer code, if it exists; otherwise return `missing`.
"""
expansion_basis(C::AbstractConcatenatedCode) = C.basis

"""
    expansion_dual_basis(C::AbstractConcatenatedCode)

Return the dual basis used to expanded the outer code, if it exists; otherwise return `missing`.
"""
expansion_dual_basis(C::AbstractConcatenatedCode) = C.dual_basis

"""
    concatenation_type(C::AbstractConcatenatedCode)

Return `:expanded`, `:same`, or `:generalized` depending on the type of concatenation.
"""
concatenation_type(C::AbstractConcatenatedCode) = C.type

#############################
# setter functions
#############################

#############################
# general functions
#############################

function _concatenated_generator_matrix(A::T, B::T) where {T<:CTMatrixTypes}
    nr_A, nc_A = size(A)
    nr_B, nc_B = size(B)
    t = div(nc_A, nr_B)
    M = zero_matrix(base_ring(A), nr_A, t * nc_B)
    for i = 1:t
        M[:, (nc_B*(i-1)+1):(nc_B*i)] = view(A, :, (nr_B*(i-1)+1):(nr_B*i)) * B
    end
    return M
end

# TODO: untested, little chance this works without error
"""
    encode(C::AbstractConcatenatedCode, v::Union{CTMatrixTypes, Vector{Int}})

Return the encoding of `v` into `C`, where `v` is either a valid input for the outer code or the full code.
"""
function encode(C::AbstractConcatenatedCode, v::Union{CTMatrixTypes,Vector{Int}})
    if typeof(v) <: CTMatrixTypes
        nr_v, nc_v = size(v)
        nr_v == 1 || nc_v == 1 || throw(ArgumentError("Vector has incorrect dimension"))
        nr_v != 1 && nc_v == 1 ? (w = transpose(v);) : (w = v;)
        nc_w = ncols(w)
        if nc_w == C.C_out.k
            # TODO: should check order and then convert if they are the same but different pointers
            base_ring(w) == C.C_out.F || throw(
                ArgumentError("Vector must have the same base ring as the outer code."),
            )
            G_out = generator_matrix(C.C_out, true)
            ismissing(C.C_out.P_stand) || (G_out = G_out * C.C_out.P_stand)
            temp = w * G_out
            if C.type == :expanded
                D = _expansion_dict(C.C_out.F, C.C_in.F, C.dual_basis)
                temp = _expand_matrix(temp, D, div(degree(C.C_out.F), degree(C.C_in.F)))
                # TODO: this is no longer a vector, only want 1 * temp row instead of full basis?
                temp = temp[1, :]
            end
            # should automatically now be in field of inner code
            Gin = generator_matrix(C.C_in, true)
            ismissing(C.C_in.P_stand) || (Gin = Gin * C.C_in.P_stand)
            return temp * Gin
        elseif nc_w == C.k
            base_ring(w) == C.F ||
                throw(ArgumentError("Vector must have the same base ring as the code."))
            return w * C.G
        else
            throw(ArgumentError("Vector has incorrect dimension"))
        end
    else
        len = length(v)
        if len == C.C_out.k
            w = matrix(C.C_out.F, 1, len, v)
            G_out = generator_matrix(C.C_out, true)
            ismissing(C.C_out.P_stand) || (G_out = G_out * C.C_out.P_stand)
            temp = w * G_out
            if C.type == :expanded
                D = _expansion_dict(C.C_out.F, C.C_in.F, C.dual_basis)
                temp = _expand_matrix(temp, D, div(degree(C.C_out.F), degree(C.C_in.F)))
                # TODO: this is no longer a vector, only want 1 * temp row instead of full basis?
                temp = temp[1, :]
            end
            # should automatically now be in field of inner code
            Gin = generator_matrix(C.C_in, true)
            ismissing(C.C_in.P_stand) || (Gin = Gin * C.C_in.P_stand)
            return temp * Gin
        elseif len == C.k
            return matrix(C.F, 1, len, v) * C.G
        else
            throw(ArgumentError("Vector has incorrect dimension"))
        end
    end
end

# permute?

#############################
# Quantum
#############################

# something like this
# function concatenatedcode(S1::CSSCode, S2::CSSCode)
#     # need them to be F-linear
#     Xstabs = vcat(Xstabilizers(S1) ⊗ I(S2), Xlogicals(S1) ⊗ Xstabilizers(S2))
#     Zstabs = vcat(Zstabilizers(S1) ⊗ I(S2), Zlogicals(S1) ⊗ Zstabilizers(S2))
#
#     this should give [[n1 n2, k1, k2, (d1x d2x, d1z d2z)]]
# end
