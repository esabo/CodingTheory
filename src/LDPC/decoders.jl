# Copyright (c) 2023 - 2024 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # MP Decoders
#############################

# Example of using Gallager A and B (out should end up [1 1 1 0 0 0 0]
# H = matrix(GF(2), [1 1 0 1 1 0 0; 1 0 1 1 0 1 0; 0 1 1 1 0 0 1]);
# v = matrix(GF(2), 7, 1, [1, 1, 0, 0, 0, 0, 0]);
# flag, out, iter, vtoc, ctov = Gallager_A(H, v, 100);
# flag, out, iter, vtoc, ctov = Gallager_B(H, v, 100);
# nm = MPNoiseModel(:BSC, 1/7);
# flag, out, iter, vtoc, ctov = sum_product(H, v, nm, 100);

struct MPNoiseModel
    type::Symbol
    cross_over_prob::Union{Float64, Missing}
    sigma::Union{Float64, Missing}
end

function MPNoiseModel(type::Symbol, x::Float64)
    if type == :BSC
        return MPNoiseModel(type, x, missing)
    elseif type == :BAWGNC
        return MPNoiseModel(type, missing, x)
    else
        throw(ArgumentError("Unsupported noise model type $type"))
    end
end

function Gallager_A(H::T, v::T; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} =
    missing, schedule::Symbol = :flooding) where T <: CTMatrixTypes

    schedule ∈ (:flooding, :serial, :custom) || throw(ArgumentError("Unknown schedule type"))

    H_Int, w, var_adj_list, check_adj_list, chn_inits = _message_passing_init(H, v, missing,
        max_iter, :A, 2, chn_inits)
    if schedule == :flooding
        return _message_passing_flooding(H_Int, w, chn_inits, _Gallager_A_check_node_message,
            var_adj_list, check_adj_list, max_iter, :A)
    end
end

function Gallager_B(H::T, v::T; max_iter::Int = 100, threshold::Int = 2, chn_inits::Union{Missing,
    Vector{Float64}} = missing, schedule::Symbol = :flooding) where T <: CTMatrixTypes
    
    H_Int, w, var_adj_list, check_adj_list, chn_inits = _message_passing_init(H, v, missing,
        max_iter, :B, threshold, chn_inits)
    if schedule == :flooding
        return _message_passing_flooding(H_Int, w, chn_inits, _Gallager_B_check_node_message,
            var_adj_list, check_adj_list, max_iter, :B, threshold)
    end
end

function sum_product(H::S, v::T, chn::MPNoiseModel; max_iter::Int = 100, chn_inits::Union{Missing,
    Vector{Float64}} = missing, schedule::Symbol = :flooding) where {S <: CTMatrixTypes, T <:
    Union{Vector{<:Real}, CTMatrixTypes}}

    H_Int, w, var_adj_list, check_adj_list, chn_inits = _message_passing_init(H, v, chn, max_iter,
        :SP, 2, chn_inits)
    if schedule == :flooding
        return _message_passing_flooding(H_Int, w, chn_inits, _SP_check_node_message,
            var_adj_list, check_adj_list, max_iter, :SP)
    elseif schedule == :serial
        return _message_passing_serial(H_Int, w, chn_inits, _SP_check_node_message,
            var_adj_list, check_adj_list, max_iter, :SP)
    end
end

# believe this can be merged into an optional argument of the above but keep for now so as not to break working code
function sum_product_decimation(H::S, v::T, chn::MPNoiseModel, decimated_bits::Vector{Tuple{Int, 
    Int}}; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} = missing, 
    schedule::Symbol = :flooding) where {S <: CTMatrixTypes, T <: Union{Vector{<:Real},
    CTMatrixTypes}}

    # what do I want out of this?
    # do I want to completely delete the columns in the decimated_bits?
    H_Int, w, var_adj_list, check_adj_list, chn_inits = _message_passing_init(H, v, chn,
        max_iter, :SP, 2, chn_inits)
    if schedule == :flooding
        return _message_passing_flooding_decimation(H_Int, w, decimated_bits, chn_inits,
            _SP_check_node_message, var_adj_list, check_adj_list, max_iter, :SP)
    elseif schedule == :serial
        return _message_passing_serial_decimation(H_Int, w, decimated_bits, chn_inits,
            _SP_check_node_message, var_adj_list, check_adj_list, max_iter, :SP)
    end
end

function sum_product_syndrome(H::S, syndrome::T, chn::MPNoiseModel; max_iter::Int = 100, 
    chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where {S <: 
    CTMatrixTypes, T <: Union{Vector{<:Real}, CTMatrixTypes}}

    H_Int, syn, var_adj_list, check_adj_list, chn_inits = _message_passing_init_syndrome(H,
        syndrome, chn, max_iter, chn_inits, :SP)
    if schedule == :flooding
        return _message_passing_flooding_syndrome(H_Int, syn, chn_inits, _SP_check_node_message,
            var_adj_list, check_adj_list, max_iter)
    elseif schedule == :serial
        return _message_passing_serial(H_Int, w, chn_inits, _SP_check_node_message,
            var_adj_list, check_adj_list, max_iter, :SP)
    end
end

function sum_product_box_plus(H::S, v::T, chn::MPNoiseModel; max_iter::Int = 100,
    chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where {S <:
    CTMatrixTypes, T <: Union{Vector{<:Real}, CTMatrixTypes}}

    H_Int, w, var_adj_list, check_adj_list, chn_inits = _message_passing_init(H, v, chn, max_iter,
        :SP, 2, chn_inits)
    if schedule == :flooding
        return _message_passing_flooding(H_Int, w, chn_inits, _SP_check_node_message_box_plus,
            var_adj_list, check_adj_list, max_iter, :SP)
    end
end

# default on attenuation?
function min_sum(H::S, v::T, chn::MPNoiseModel; max_iter::Int = 100, attenuation::Float64 =
    0.5, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where
    {S <: CTMatrixTypes, T <: Vector{<:AbstractFloat}}

    H_Int, w, var_adj_list, check_adj_list, chn_inits = _message_passing_init(H, v, chn, max_iter,
        :MS, 2, chn_inits)
    if schedule == :flooding
        return _message_passing_flooding(H_Int, w, chn_inits, _MS_check_node_message,
            var_adj_list, check_adj_list, max_iter, :MS, attenuation)
    end
end

function _message_passing_init(H::S, v::T, chn::Union{Missing, MPNoiseModel}, max_iter::Int,
    kind::Symbol, Bt::Int, chn_inits::Union{Missing, Vector{Float64}}) where {S <: CTMatrixTypes,
    T <: Union{Vector{<:Real}, CTMatrixTypes}}

    kind ∈ (:SP, :MS, :A, :B) || throw(ArgumentError("Unknown value for parameter kind"))
    kind ∈ (:SP, :MS) && ismissing(chn) && throw(ArgumentError(":SP and :MS require a noise model"))
    Int(order(base_ring(H))) == 2 ||
        throw(ArgumentError("Currently only implemented for binary codes"))
    num_check, num_var = size(H)
    num_check > 0 && num_var > 0 || throw(ArgumentError("Input matrix of improper dimension"))
    length(v) == num_var || throw(ArgumentError("Vector has incorrect dimension"))
    (kind == :B && !(1 <= Bt <= num_check)) &&
        throw(DomainError("Improper threshold for Gallager B"))
    2 <= max_iter || throw(DomainError("Number of maximum iterations must be at least two"))
    kind ∈ (:SP, :MS) && chn.type == :BAWGNC && !isa(v, Vector{<:AbstractFloat}) &&
        throw(DomainError("Received message should be a vector of floats for BAWGNC."))
    kind ∈ (:SP, :MS) && chn.type == :BSC && !isa(v, Vector{Int}) && !isa(v, CTMatrixTypes) &&
        throw(DomainError("Received message should be a vector of Ints for BSC."))
    
    H_Int = FpmattoJulia(H)
    w = if T <: CTMatrixTypes
        Int.(data.(v)[:])
    else
        copy(v)
    end
    check_adj_list = [[] for _ in 1:num_check]
    var_adj_list = [[] for _ in 1:num_var]

    for r in 1:num_check
        for c in 1:num_var
            if !iszero(H_Int[r, c])
                push!(check_adj_list[r], c)
                push!(var_adj_list[c], r)
            end
        end
    end

    if ismissing(chn_inits) && kind ∈ (:SP, :MS)
        chn_inits = if chn.type == :BSC
            _channel_init_BSC(w, chn.cross_over_prob)
        elseif chn.type == :BAWGNC && kind == :SP
            _channel_init_BAWGNC_SP(w, chn.sigma)
        elseif chn.type == :BAWGNC && kind == :MS
            _channel_init_BAWGNC_MS(w)
        end
    end

    # R = kind ∈ (:A, :B) ? Int : Float64
    # check_to_var_messages = zeros(R, num_check, num_var, max_iter)
    # var_to_check_messages = zeros(R, num_var, num_check, max_iter)

    return H_Int, w, var_adj_list, check_adj_list, chn_inits#, check_to_var_messages, var_to_check_messages
end

function _message_passing_init_syndrome(H::S, syndrome::T, chn::Union{Missing, MPNoiseModel}, 
    max_iter::Int, chn_inits::Union{Missing, Vector{Float64}}, kind::Symbol) where {S <:
    CTMatrixTypes, T <: Union{Vector{<:Real}, CTMatrixTypes}}

    Int(order(base_ring(H))) == 2 ||
        throw(ArgumentError("Currently only implemented for binary codes"))
    num_check, num_var = size(H)
    num_check > 0 && num_var > 0 || throw(ArgumentError("Input matrix of improper dimension"))
    length(syndrome) == num_check || throw(ArgumentError("Syndrome has incorrect dimension"))
    2 <= max_iter || throw(DomainError("Number of maximum iterations must be at least two"))
    
    H_Int = FpmattoJulia(H)
    syn = if T <: CTMatrixTypes
        Int.(data.(syndrome)[:])
    else
        copy(syndrome)
    end
    check_adj_list = [[] for _ in 1:num_check]
    var_adj_list = [[] for _ in 1:num_var]

    for r in 1:num_check
        for c in 1:num_var
            if !iszero(H_Int[r, c])
                push!(check_adj_list[r], c)
                push!(var_adj_list[c], r)
            end
        end
    end

    # TODO: this part needs to be redone
    # if ismissing(chn_inits)
    #     chn_inits = if chn.type == :BSC
    #         _channel_init_BSC(w, chn.cross_over_prob)
    #     elseif chn.type == :BAWGNC && kind == :SP
    #         _channel_init_BAWGNC_SP(w, chn.sigma)
    #     elseif chn.type == :BAWGNC && kind == :MS
    #         _channel_init_BAWGNC_MS(w)
    #     end
    # end

    chn_inits = [log((1 - chn.cross_over_prob) / chn.cross_over_prob) for _ in 1:num_var]

    # R = kind ∈ (:A, :B) ? Int : Float64
    # check_to_var_messages = zeros(R, num_check, num_var, max_iter)
    # var_to_check_messages = zeros(R, num_var, num_check, max_iter)

    return H_Int, syn, var_adj_list, check_adj_list, chn_inits#, check_to_var_messages, var_to_check_messages
end

function _message_passing_flooding(H::Matrix{UInt64}, w::Vector{T}, chn_inits::Union{Missing,
    Vector{Float64}}, c_to_v_mess::Function, var_adj_list::Vector{Vector{Any}},
    check_adj_list::Vector{Vector{Any}}, max_iter::Int, kind::Symbol, Bt::Int = 2,
    attenuation::Float64 = 0.5) where T <: Union{Int, AbstractFloat}

    # the inclusion of the kind statements add less than a microsecond
    num_check, num_var = size(H)
    # could reduce this stuff to UInt8 if needed
    S = kind ∈ (:A, :B) ? Int : Float64
    current_bits = zeros(Int, num_var)
    if kind in (:SP, :MS)
        totals = zeros(S, 1, num_var)
    end
    syn = zeros(Int, num_check)
    # max_iter += 1 # probably should copy this
    check_to_var_messages = zeros(S, num_check, num_var, 2)
    var_to_check_messages = zeros(S, num_var, num_check, 2)
    
    # remove and accept as input
    if kind ∈ (:SP, :MS)
        @inbounds for vn in 1:num_var
            var_to_check_messages[vn, var_adj_list[vn], 1] .= chn_inits[vn]
        end
    elseif kind ∈ (:A, :B)
        @inbounds for vn in 1:num_var
            var_to_check_messages[vn, var_adj_list[vn], :] .= w[vn]
        end
    end

    iter = 1
    curr_iter = 1
    prev_iter = 2
    # does this propagate correctly?
    @inbounds while iter ≤ max_iter
        # variable node is already done for first iteration, so start with check nodes
        @simd for cn in 1:num_check
            for v1 in check_adj_list[cn]
                check_to_var_messages[cn, v1, curr_iter] = c_to_v_mess(cn, v1, iter, check_adj_list,
                    var_to_check_messages, attenuation)
            end
        end

        # one full iteration done, check if converged
        if kind ∈ (:SP, :MS)
            @simd for vn in 1:num_var
                totals[vn] = chn_inits[vn]
                for c in var_adj_list[vn]
                    totals[vn] += check_to_var_messages[c, vn, curr_iter]
                end
                current_bits[vn] = totals[vn] >= 0 ? 0 : 1
            end
        elseif kind ∈ (:A, :B)
            @simd for i in 1:num_var
                len = length(var_adj_list[i])
                one_count = count(isone, view(check_to_var_messages, var_adj_list[i], i, curr_iter))
                d = fld(len, 2)
                current_bits[i] = one_count + (isone(w[i]) && iseven(len)) > d
            end
        end

        LinearAlgebra.mul!(syn, H, current_bits)
        iszero(syn .% 2) && return true, current_bits, iter, var_to_check_messages,
            check_to_var_messages
        iter += 1
        temp = curr_iter
        curr_iter = prev_iter
        prev_iter = temp

        if iter <= max_iter
            for vn in 1:num_var
                for c1 in var_adj_list[vn]
                    if kind ∈ (:SP, :MS)
                        var_to_check_messages[vn, c1, curr_iter] = totals[vn] -
                            check_to_var_messages[c1, vn, prev_iter]
                    elseif kind == :A && length(var_adj_list[vn]) > 1
                        if all(!Base.isequal(w[vn]), check_to_var_messages[c2, vn, prev_iter] for
                            c2 in var_adj_list[vn] if c1 != c2)

                            var_to_check_messages[vn, c1, curr_iter] ⊻= 1 
                        end
                    elseif kind == :B && length(var_adj_list[vn]) >= Bt
                        if count(!Base.isequal(w[vn]), check_to_var_messages[c2, vn, prev_iter] for 
                            c2 in var_adj_list[vn] if c1 != c2) >= Bt

                            var_to_check_messages[vn, c1, curr_iter] ⊻= 1 
                        end
                    end
                end
            end
        end
    end

    return false, curr, iter, var_to_check_messages, check_to_var_messages
end

function _message_passing_serial(H::Matrix{UInt64}, w::Vector{T}, chn_inits::Union{Missing,
    Vector{Float64}}, c_to_v_mess::Function, var_adj_list::Vector{Vector{Any}},
    check_adj_list::Vector{Vector{Any}}, max_iter::Int, kind::Symbol, Bt::Int = 2,
    attenuation::Float64 = 0.5) where T <: Union{Int, AbstractFloat}

    # the inclusion of the kind statements add less than a microsecond
    num_check, num_var = size(H)
    current_bits = zeros(Int, num_var)
    totals = zeros(Float64, 1, num_var)
    syn = zeros(Int, num_check)
    # max_iter += 1 # probably should copy this
    check_to_var_messages = zeros(Float64, num_check, num_var, 2)
    var_to_check_messages = zeros(Float64, num_var, num_check, 2)
    
    # remove and accept as input
    @inbounds for vn in 1:num_var
        var_to_check_messages[vn, var_adj_list[vn], 1] .= chn_inits[vn]
    end

    iter = 1
    curr_iter = 1
    prev_iter = 2
    # ϕ(x) = -log(tanh(0.5 * x))
    # does this propagate correctly?
    @inbounds while iter ≤ max_iter
        # for cn in 1:num_check
        #     S = 0.0
        #     sgn = 1
        #     for v in check_adj_list[cn]
        #         x = var_to_check_messages[v, cn, curr_iter]
        #         if x > 0
        #             S += ϕ(x)
        #         else
        #             S += ϕ(-x)
        #             sgn *= -1
        #         end
        #     end
        #     S = sgn * S

        #     for v in check_adj_list[cn]
        #         Q_temp = var_to_check_messages[v, cn, curr_iter] - check_to_var_messages[cn, v, curr_iter]
        #         if Q_temp > 0
        #             temp = S - ϕ(Q_temp)
        #         else
        #             temp = S + ϕ(-Q_temp)
        #         end
        #         if temp > 0
        #             check_to_var_messages[cn, v, curr_iter] = ϕ(temp)
        #         else
        #             check_to_var_messages[cn, v, curr_iter] = -ϕ(-temp)
        #         end
        #         var_to_check_messages[v, cn, curr_iter] = Q_temp + check_to_var_messages[cn, v, curr_iter]
        #     end

        for c in 1:num_check
            for v in check_adj_list[c]
                check_to_var_messages[c, v, curr_iter] = c_to_v_mess(c, v, curr_iter,
                    check_adj_list, var_to_check_messages, attenuation)
            end
        end
        for v in 1:num_var
            for c in var_adj_list[v]
                var_to_check_messages[v, c, curr_iter] = var_to_check_messages[v, c, prev_iter] + check_to_var_messages[c, v, curr_iter]
            end
        end
        # - check_to_var_messages[c, v, prev_iter] 

        # one full iteration done, check if converged
        @simd for vn in 1:num_var
            totals[vn] = chn_inits[vn]
            for c in var_adj_list[vn]
                totals[vn] += check_to_var_messages[c, vn, curr_iter]
            end
            curr[vn] = totals[vn] >= 0 ? 0 : 1
        end

        LinearAlgebra.mul!(syn, H, curr)
        iszero(syn .% 2) && return true, curr, iter, var_to_check_messages, check_to_var_messages
        iter += 1
        temp = curr_iter
        curr_iter = prev_iter
        prev_iter = temp
    end
    return false, curr, iter, var_to_check_messages, check_to_var_messages
end

function _message_passing_flooding_syndrome(H::Matrix{UInt64}, syndrome::Vector{T},
    chn_inits::Union{Missing, Vector{Float64}}, c_to_v_mess::Function, 
    var_adj_list::Vector{Vector{Any}}, check_adj_list::Vector{Vector{Any}}, max_iter::Int,
    attenuation::Float64 = 0.5) where T <: Union{Int, AbstractFloat}

    # the inclusion of the kind statements add less than a microsecond
    num_check, num_var = size(H)
    current_bits = zeros(Int, num_var)
    totals = zeros(Float64, 1, num_var)
    syn = zeros(Int, num_check)
    # max_iter += 1 # probably should copy this
    check_to_var_messages = zeros(Float64, num_check, num_var, 2)
    var_to_check_messages = zeros(Float64, num_var, num_check, 2)
    
    # remove and accept as input
    @inbounds for vn in 1:num_var
        var_to_check_messages[vn, var_adj_list[vn], 1] .= chn_inits[vn]
        if chn_inits[vn] > 0
            current_bits[vn] = 1
        end
    end
    # display(var_to_check_messages

    iter = 1
    curr_iter = 1
    prev_iter = 2
    # does this propagate correctly?
    # @inbounds 
    while iter ≤ max_iter
        # variable node is already done for first iteration, so start with check nodes
        for cn in 1:num_check
            count = sum(current_bits[check_adj_list[cn]])
            for v1 in check_adj_list[cn]
                check_to_var_messages[cn, v1, curr_iter] = 0.0
                # here I want to loop through all possible options of errors
                for e in 0:1
                    # check if this e with the other values gives the correct syndrome bit
                    if (count + e) % 2 == syndrome[cn]
                        check_to_var_messages[cn, v1, curr_iter] += c_to_v_mess(cn, v1, prev_iter,
                            check_adj_list, var_to_check_messages, attenuation)
                    end
                end
            end
        end

        # one full iteration done, check if converged
        @simd for vn in 1:num_var
            totals[vn] = chn_inits[vn]
            for c in var_adj_list[vn]
                totals[vn] += check_to_var_messages[c, vn, curr_iter]
            end
            current_bits[vn] = totals[vn] >= 0 ? 0 : 1
        end

        LinearAlgebra.mul!(syn, H, current_bits)
        syn .% 2 == syndrome && return true, current_bits, iter, var_to_check_messages, check_to_var_messages
        iter += 1
        temp = curr_iter
        curr_iter = prev_iter
        prev_iter = temp

        if iter <= max_iter
            for vn in 1:num_var
                for c1 in var_adj_list[vn]
                    var_to_check_messages[vn, c1, curr_iter] = totals[vn] -
                        check_to_var_messages[c1, vn, prev_iter]
                end
            end
        end
    end

    return false, current_bits, iter, var_to_check_messages, check_to_var_messages
end

# this is the working version from master to benchmark against
function _message_passing_old(H::Matrix{UInt64}, w::Vector{T}, chn::Union{Missing, MPNoiseModel},
    c_to_v_mess::Function, var_adj_list::Vector{Vector{Any}}, check_adj_list::Vector{Vector{Any}},
    max_iter::Int, kind::Symbol, Bt::Int = 2, attenuation::Float64 = 0.5) where T <: Union{Int,
    AbstractFloat}

    num_check, num_var = size(H)
    S = kind ∈ (:A, :B) ? Int : Float64
    curr = zeros(Int, num_var)
    if kind in (:SP, :MS)
        totals = zeros(S, 1, num_var)
    end
    syn = zeros(Int, num_check)
    max_iter += 1 # probably should copy this
    check_to_var_messages = zeros(S, num_check, num_var, max_iter)
    var_to_check_messages = zeros(S, num_var, num_check, max_iter)
    
    iter = 1
    if kind in (:SP, :MS)
        chn_inits = if chn.type == :BSC
            _channel_init_BSC(w, chn.cross_over_prob)
        elseif chn.type == :BAWGNC && kind == :SP
            _channel_init_BAWGNC_SP(w, chn.sigma)
        elseif chn.type == :BAWGNC && kind == :MS
            _channel_init_BAWGNC_MS(w)
        end
        for vn in 1:num_var
            var_to_check_messages[vn, var_adj_list[vn], 1] .= chn_inits[vn]
        end
    elseif kind in (:A, :B)
        for vn in 1:num_var
            var_to_check_messages[vn, var_adj_list[vn], :] .= w[vn]
        end
    end

    while iter < max_iter
        for cn in 1:num_check
            for v1 in check_adj_list[cn]
                check_to_var_messages[cn, v1, iter] = c_to_v_mess(cn, v1, iter, check_adj_list,
                    var_to_check_messages, attenuation)
            end
        end

        if kind in (:SP, :MS)
            for vn in 1:num_var
                totals[vn] = chn_inits[vn]
                for c in var_adj_list[vn]
                    totals[vn] += check_to_var_messages[c, vn, iter]
                end
            end
        end

        if kind in (:SP, :MS)
            @simd for i in 1:num_var
                curr[i] = totals[i] >= 0 ? 0 : 1
            end
        elseif kind in (:A, :B)
            @simd for i in 1:num_var
                len = length(var_adj_list[i])
                one_count = count(isone, view(check_to_var_messages, var_adj_list[i], i, iter))
                d = fld(len, 2)
                curr[i] = one_count + (isone(w[i]) && iseven(len)) > d
            end
        end

        LinearAlgebra.mul!(syn, H, curr)
        # @show curr
        # @show syn .% 2
        iszero(syn .% 2) && return true, curr, iter, var_to_check_messages, check_to_var_messages
        iter += 1

        if iter <= max_iter
            for vn in 1:num_var
                for c1 in var_adj_list[vn]
                    if kind in (:SP, :MS)
                        var_to_check_messages[vn, c1, iter] = totals[vn] -
                            check_to_var_messages[c1, vn, iter - 1]
                    elseif kind == :A && length(var_adj_list[vn]) > 1
                        if all(!Base.isequal(w[vn]), check_to_var_messages[c2, vn, iter - 1] for c2 
                            in var_adj_list[vn] if c1 != c2)

                            var_to_check_messages[vn, c1, iter] ⊻= 1 
                        end
                    elseif kind == :B && length(var_adj_list[vn]) >= Bt
                        if count(!Base.isequal(w[vn]), check_to_var_messages[c2, vn, iter - 1] for 
                            c2 in var_adj_list[vn] if c1 != c2) >= Bt

                            var_to_check_messages[vn, c1, iter] ⊻= 1 
                        end
                    end
                end
            end
        end
    end

    return false, curr, iter, var_to_check_messages, check_to_var_messages
end

function _channel_init_BSC(v::Vector{T}, p::Float64) where T <: Integer
    temp = log((1 - p) / p)
    chn_init = zeros(Float64, length(v))
    for i in eachindex(v)
        @inbounds chn_init[i] = (-1)^v[i] * temp
    end
    return chn_init
end

function _channel_init_BAWGNC_SP(v::Vector{T}, σ::Float64) where T <: AbstractFloat
    temp = 2 / σ^2
    chn_init = zeros(Float64, length(v))
    for i in eachindex(v)
        @inbounds chn_init[i] = temp * v[i]
    end
    return chn_init
end

_channel_init_BAWGNC_MS(v::Vector{T}) where T <: AbstractFloat = v

function _SP_check_node_message(cn::Int, v1::Int, iter, check_adj_list, var_to_check_messages,
    atten = missing)

    ϕ(x) = -log(tanh(0.5 * x))
    temp = 0.0
    s = 1
    @inbounds for v2 in check_adj_list[cn]
        if v2 != v1
            x = var_to_check_messages[v2, cn, iter]
            # Note that x should never be 0 unless there is an erasure.
            # For now, this is coded as if there will never be an erasure.
            # This will simply error if x == 0.
            if x > 0
                temp += ϕ(x)
            else
                temp += ϕ(-x)
                s *= -1
            end
        end
    end
    return s * ϕ(temp)
end

⊞(a, b) = log((1 + exp(a + b)) / (exp(a) + exp(b)))
⊞(a...) = reduce(⊞, a...)
function _SP_check_node_message_box_plus(cn::Int, v1::Int, iter, check_adj_list,
    var_to_check_messages, atten = missing)

    @inbounds ⊞(var_to_check_messages[v2, cn, iter] for v2 in check_adj_list[cn] if v2 != v1)
end

function _MS_check_node_message(cn::Int, v1::Int, iter, check_adj_list, var_to_check_messages,
    attenuation::Float64 = 0.5)

    temp = var_to_check_messages[check_adj_list[cn][1], cn, iter]
    s = 1
    @inbounds for v2 in check_adj_list[cn]
        if v2 != v1
            x = var_to_check_messages[v2, cn, iter]
            # Note that x should never be 0
            if x > 0
                temp > x && (temp = x;)
            else
                temp > -x && (temp = -x;)
                s *= -1
            end
        end
    end
    return s * attenuation * temp
end

function _Gallager_A_check_node_message(cn::Int, v1::Int, iter::Int, check_adj_list,
    var_to_check_messages, atten = missing)

    @inbounds reduce(⊻, var_to_check_messages[v, cn, iter] for v in check_adj_list[cn] if v != v1)
end
_Gallager_B_check_node_message(cn::Int, v1::Int, iter::Int, check_adj_list, var_to_check_messages,
    atten = missing) = _Gallager_A_check_node_message(cn, v1, iter, check_adj_list,
    var_to_check_messages, atten)

# Mansour, Shanbhag, "Turbo Decoder Architectures for Low-Density Parity-Check Codes" (2002)

function find_MP_schedule(H::CTMatrixTypes)
    num_check, num_var = size(H)
    num_check > 0 && num_var > 0 || throw(ArgumentError("Input matrix of improper dimension"))

    check_adj_list = [[] for _ in 1:num_check]
    for r in 1:num_check
        for c in 1:num_var
            iszero(H[r, c]) || push!(check_adj_list[r], c)
        end
    end

    sched_list = [[1]]
    for cn in 2:num_check
        found = false
        for sched in sched_list
            if !any(x ∈ check_adj_list[y] for y in sched for x ∈ check_adj_list[cn])
                push!(sched, cn)
                sort!(sched_list, lt=(x, y) -> length(x) < length(y))
                found = true
                break
            end
        end
        !found && push!(sched_list, [cn])
    end
    
    return sched_list
end

#############################
        # LP Decoders
#############################

# function _init_LP_decoder_LDPC end

# function _LP_decoder_LDPC end

# TODO: docstring and in extension
"""
    LP_decoder_LDPC(H::Union{CTMatrixTypes, AbstractMatrix{<:Number}}, v::Union{CTMatrixTypes, Vector{<:Integer}}, Ch::BinarySymmetricChannel)
    LP_decoder_LDPC(C::AbstractLinearCode, v::Union{CTMatrixTypes, Vector{<:Integer}}, Ch::BinarySymmetricChannel)

Return

# Note
- Run `using JuMP, GLPK` to activate this extension.
"""
function LP_decoder_LDPC end

#############################
          # Methods
#############################

function _channel_to_SNR(chn::MPNoiseModel)
    if cnh.type == :BAWGNC
        -10 * log10(chn.sigma^2)
    else
        throw(ArgumentError("Only supports BAWGNC currently."))
    end
end

function _channel_to_SNR(type::Symbol, sigma::Real)
    if type == :BAWGNC
        -10 * log10(sigma^2)
    else
        throw(ArgumentError("Only supports BAWGNC currently."))
    end
end

#############################
        # Simulations
#############################

# function decodersimulation(H::CTMatrixTypes, decoder::Symbol, noisetype::Symbol,
#                            noise::Union{Vector{<:Real}, AbstractRange{<:Real}},
#                            max_iter::Int = 100, numruns::Int = 100000, seed::Union{Int, Nothing} = nothing,
#                            attenuation::Float64 = 0.5)

#     decodersimulation(H, decoder, noisetype, noise, max_iter,
#                       [numruns for n in noise], seed, attenuation)
# end

# function decodersimulation(H::CTMatrixTypes, decoder::Symbol, noisetype::Symbol,
#                            noise::Union{Vector{<:Real}, AbstractRange{<:Real}},
#                            max_iter::Int = 100, numruns::Vector{Int} = [100000 for n in noise],
#                            seed::Union{Int, Nothing} = nothing, attenuation::Float64 = 0.5)

#     decoder in (:A, :B, :SP, :MS) || throw(ArgumentError("Unsupported decoder"))
#     noisetype in (:BSC, :BAWGNC) || throw(ArgumentError("Only supports BSC and BAWGNC"))
#     decoder in (:A, :B) && noisetype == :BAWGNC && throw(ArgumentError("BAWGNC not supported for Gallager decoders."))
#     0 <= minimum(noise) || throw(ArgumentError("Must have non-negative noise"))
#     maximum(noise) > 1 && noisetype == :BSC && throw(ArgumentError("Crossover probability must be in the range [0,1]"))

#     # We use an explicit pseudoRNG with the given seed. (if `nothing`, it's just random)
#     # Note that Xoshiro is default in Julia. Threading breaks reproducibility.
#     rng = Xoshiro(seed)

#     FER = zeros(length(noise))
#     BER = zeros(length(noise))
#     n = ncols(H)
#     cnmsg = decoder == :SP ? _SP_check_node_message : _MS_check_node_message

#     for k in eachindex(noise)
#         chn = MPNoiseModel(noisetype, noise[k])
#         w = noisetype == :BSC ? zeros(Int, n) : ones(n)
#         H_Int, _, var_adj_list, check_adj_list = _message_passing_init(H, w, chn, max_iter, decoder, 2)
#         FEtotal = 0 # number of frame errors
#         BEtotal = 0 # number of bit errors
#         @threads for i in 1:numruns[k]
#         # for i in 1:numruns
#             for j in 1:n
#                 if noisetype == :BSC
#                     w[j] = Int(rand(rng) < chn.cross_over_prob)
#                 else # BAWGNC
#                     w[j] = 1.0 + randn(rng, Float64) * chn.sigma
#                 end
#             end
#             iszero(w) && continue
#             flag, curr, _, _, _ = _message_passing(H_Int, w, chn, cnmsg, var_adj_list, check_adj_list, max_iter, decoder, 2, attenuation)
#             if !(flag && iszero(curr))
#                 FEtotal += 1
#                 BEtotal += count(!iszero, curr)
#             end
#         end
#         FER[k] = FEtotal / numruns[k]
#         BER[k] = BEtotal / (numruns[k] * n)
#     end
#     return FER, BER
# end

# function decodersimulation(H::CTMatrixTypes, decoder::Symbol, noisetype::Symbol,
#                            noise::Union{Vector{T}, AbstractRange{T}} where T<:Real,
#                            max_iter::Int = 100, numruns::Int = 1000,
#                            seed::Union{Int, Nothing} = nothing, attenuation::Float64 = 1.0)

#     decoder in (:A, :B, :SP, :MS, :LP) || throw(ArgumentError("Unsupported decoder"))
#     noisetype == :BSC || throw(ArgumentError("Only supports BSC"))
#     0 <= minimum(noise) || throw(ArgumentError("Must have non-negative noise"))
#     maximum(noise) > 1 && noisetype == :BSC && throw(ArgumentError("Crossover probability must be in the range [0,1]"))

#     # we'll use an explicit pseudoRNG with the given seed. Note that Xoshiro is default in Julia.
#     # `seed == nothing` just gives a random choice, i.e., the default is not reproducible.
#     rng = Xoshiro(seed)

#     FER = zeros(length(noise))
#     BER = zeros(length(noise))
#     ε = zeros(length(noise))
#     n = ncols(H)
#     cnmsg = decoder == :SP ? _SP_check_node_message : _MS_check_node_message

#     model = decoder == :LP ? _initLPdecoderLDPC(H) : nothing

#     for k in eachindex(noise)

#         # p[i] is the probability of having i - 1 bit errors
#         temp = BigFloat(noise[k])
#         p = BigFloat[temp^i * (1 - temp)^(n - i) * binomial(big(n), big(i)) /
#                      sum(temp^j * (1 - temp)^(n - j) * binomial(big(n), big(j))
#                          for j in 0:n) for i in 0:n]
#         p_partialsum = [sum(p[j] for j in 1:i) for i in 1:length(p)]
#         maxnerr = max(findfirst(p_partialsum .>= 1 - BigFloat("1e-9")) - 1, 6)
#         # @show maxnerr
#         ε[k] = 1 - p_partialsum[maxnerr + 1]

#         chn = MPNoiseModel(noisetype, noise[k])
#         w = noisetype == :BSC ? zeros(Int, n) : ones(n)
#         if decoder == :LP
#             noisemodel = BSC(noise[k])
#         else
#             H_Int, _, var_adj_list, check_adj_list = _message_passing_init(H, w, chn, max_iter, decoder, 2)
#         end

#         FEtotal = zeros(Int, maxnerr)
#         BEtotal = zeros(Int, maxnerr)
#         FER[k] = p[1]
#         BER[k] = p[1]

#         for e in 1:maxnerr

#             # importance sampling:
#             numruns_for_e = Int(cld(numruns * p[e+1], sum(p[i] for i in 2:maxnerr+1)))

#             # naive sampling, still quite good:
#             # numruns_for_e = Int(cld(numruns, maxnerr))

#             for i in 1:numruns_for_e
#                 w = zeros(Int, n)
#                 w[shuffle(rng, 1:n)[1:e]] .= 1
#                 if decoder == :LP
#                     curr = _LPdecoderLDPC(model, w, noisemodel)
#                     flag = all(isinteger(x) for x in curr)
#                 else
#                     flag, curr, _, _, _ = _message_passing(H_Int, w, chn, cnmsg, var_adj_list, check_adj_list, max_iter, decoder, 2, attenuation)
#                 end
#                 if !(flag && iszero(curr))
#                     FEtotal[e] += 1
#                     BEtotal[e] += count(!iszero, curr)
#                 end
#             end
#             FER[k] += p[e+1] * (numruns_for_e - FEtotal[e]) / numruns_for_e
#             BER[k] += p[e+1] * (numruns_for_e * n - BEtotal[e]) / (numruns_for_e * n)
#         end
#         FER[k] = 1 - FER[k]
#         BER[k] = 1 - BER[k]
#     end
#     return FER, BER, ε
# end

# using Plots: plot, savefig, xticks, yticks, xticks!, yticks!
# function testsimulation(figfilename = "test.png")
#     # H = paritycheckmatrix(regularLDPCCode(500, 6, 3));
#     H = matrix(GF(2), 10, 20, [1 0 1 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0;
#                                0 1 0 1 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0;
#                                0 0 1 0 1 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0;
#                                1 0 0 1 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0;
#                                0 1 0 0 1 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0;
#                                1 1 0 0 0 0 1 0 1 0 0 0 0 0 0 1 0 0 0 0;
#                                0 1 1 0 0 0 0 1 0 1 0 0 0 0 0 0 1 0 0 0;
#                                0 0 1 1 0 1 0 0 1 0 0 0 0 0 0 0 0 1 0 0;
#                                0 0 0 1 1 0 1 0 0 1 0 0 0 0 0 0 0 0 1 0;
#                                1 0 0 0 1 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1]);

#     p = 10 .^ collect(-4:.2:-0.8);
#     @time F1, B1, e1 = CodingTheory.decodersimulation(H, :SP, :BSC, p, 10, 5000, 1, 1.0);
#     @time F2, B2, e2 = CodingTheory.decodersimulation(H, :MS, :BSC, p, 10, 5000, 1, 1.0);
#     @time F3, B3, e3 = CodingTheory.decodersimulation(H, :MS, :BSC, p, 10, 5000, 1, 0.5);
#     @time F4, B4, e4 = CodingTheory.decodersimulation(H, :MS, :BSC, p, 10, 5000, 1, 0.1);
#     @time F5, B5, e5 = CodingTheory.decodersimulation(H, :LP, :BSC, p, 10, 5000, 1);

#     plt = plot(log10.(p), log10.([F1 F2 F3 F4 F5]),
#                label = ["FER, SP" "FER, MS atten=1.0" "FER, MS atten=0.5" "FER, MS atten=0.1" "FER, LP"],
#                xlabel = "Crossover probability",
#                ylabel = "Error rate",
#                title = "BSC with a [20,10,5] code",
#                # xlims = (0, maximum(p) * 1.02),
#                # ylims = (0, max(maximum(FER), maximum(BER)) * 1.02),
#                # xticks = (-4:-1, ["1e-4", "1e-3", "1e-2", "1e-1"]),
#                # yticks = (-6:0, ["1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1", "1e0"]),
#                # yscale = :log,
#                marker = :dot);
#     xticks!(plt, (xticks(plt)[1][1], "1e" .* xticks(plt)[1][2]));
#     yticks!(plt, (yticks(plt)[1][1], "1e" .* yticks(plt)[1][2]));


#     # σ = 0.1:0.1:1
#     # FER, BER = decodersimulation(H, :SP, :BAWGNC, p, 100, 100, 123);
#     # SNR = CodingTheory._channel_to_SNR.(:BAWGNC, σ)
#     # plt = plot(SNR, [FER BER],
#     #            label = ["FER" "BER"],
#     #            xlabel = "Noise (dB)",
#     #            ylabel = "Error rate",
#     #            marker = :dot);


#     savefig(plt, figfilename);
# end


# function profiletest()
#     C = BCHCode(2,103,3);
#     H = paritycheckmatrix(C);
#     chn = MPNoiseModel(:BSC, 0.01);
#     v = zero_matrix(C.F, 1, 103);
#     v[1,1] = 1;
#     sum_product(H, v, chn);
#     Profile.clear()
#     @profile sum_product(H, v, chn)
# end
