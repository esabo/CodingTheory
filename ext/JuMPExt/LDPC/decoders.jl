# Copyright (c) 2023 - 2024 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # LP Decoders
#############################

function _init_LP_decoder_LDPC(H::Union{CodingTheory.CTMatrixTypes, AbstractMatrix{<:Number}})
    check_adj_list, _ = CodingTheory._node_adjacencies(H)
    subsets = Vector{Vector{Vector{Int}}}()
    nr, nc = size(H)
    hasi = [[Vector{Int}() for _ in 1:nc] for _ in 1:nr]
    wmap = zeros(Int, nr)
    curr = 1
    for (j, cn) in enumerate(check_adj_list)
        wmap[j] = curr
        inner_subsets = Vector{Vector{Int}}()
        for S in powerset(cn)
            if iseven(length(S)) # && !isempty(S)
                push!(inner_subsets, S)
                for i in S
                    push!(hasi[j][i], curr)
                end
                curr += 1
            end
        end
        push!(subsets, inner_subsets)
    end

    model = Model(GLPK.Optimizer)
    @variable(model, 0 <= f[1:nc] <= 1)
    @variable(model, 0 <= w[1:curr - 1] <= 1)
    for i in 1:nr
        if i != nr
            @constraint(model, sum(w[wmap[i]:wmap[i + 1] - 1]) == 1)
        else
            @constraint(model, sum(w[wmap[i]:end]) == 1)
        end
    end
    for j in 1:nr
        for i in 1:nc
            if !isempty(hasi[j][i])
                @constraint(model, f[i] == sum(w[hasi[j][i]]))
            end
        end
    end
    @objective(model, Min, sum(0 * f[i] for i in 1:nc))
    return model
end
_init_LP_decoder_LDPC(C::AbstractLinearCode) = _init_LP_decoder_LDPC(parity_check_matrix(C))

function _LP_decoder_LDPC(model::JuMP.Model, v::Union{CodingTheory.CTMatrixTypes,
    Vector{<:Integer}}, Ch::BinarySymmetricChannel)

    γ = CodingTheory._channel_init_BSC(isa(v, Vector) ? v : Int.(data.(v))[:], Ch.param)
    @objective(model, Min, dot(γ, model[:f]))
    optimize!(model)
    termination_status(model) == MOI.INFEASIBLE && throw(DomainError("No solution exists"))
    @assert termination_status(model) == MOI.OPTIMAL "Didn't find an optimal point"
    w = value.(model[:f])
    # all(isinteger(x) for x in w) || @warn "Solution is not integral"
    return w
end

function CodingTheory.LP_decoder_LDPC(H::Union{CodingTheory.CTMatrixTypes, AbstractMatrix{<:Number}}, v::Union{CodingTheory.CTMatrixTypes, Vector{<:Integer}}, Ch::BinarySymmetricChannel)
    
    model = _init_LP_decoder_LDPC(H)
    return _LP_decoder_LDPC(model, v, Ch)
end
CodingTheory.LP_decoder_LDPC(C::AbstractLinearCode, v::Union{CodingTheory.CTMatrixTypes,
    Vector{<:Integer}}, Ch::BinarySymmetricChannel) = CodingTheory.LP_decoder_LDPC(
        parity_check_matrix(C), v, Ch)
