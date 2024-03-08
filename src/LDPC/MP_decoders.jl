# Copyright (c) 2023 - 2024 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
         # Gallager
#############################

function _Gallager_A_check_node_message(cn::Int, v1::Int, iter::Int, check_adj_list,
    var_to_check_messages, atten = missing)

    @inbounds reduce(⊻, var_to_check_messages[v, cn, iter] for v in check_adj_list[cn] if v != v1)
end
_Gallager_B_check_node_message(cn::Int, v1::Int, iter::Int, check_adj_list, var_to_check_messages,
    atten = missing) = _Gallager_A_check_node_message(cn, v1, iter, check_adj_list,
    var_to_check_messages, atten)

# TODO: does it make sense to do a Gallager syndrome setup?
# TODO does it make sense to take channel inits here?
"""
    Gallager_A(H::T, v::T; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where T <: CTMatrixTypes

Run the Gallager-A decoder with the parity-check matrix `H` and received vector `v`.
"""
function Gallager_A(H::T, v::T; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} =
    missing, schedule::Symbol = :flooding) where T <: CTMatrixTypes

    schedule ∈ (:flooding, :serial, :custom) || throw(ArgumentError("Unknown schedule type"))

    H_Int, w, var_adj_list, check_adj_list, chn_inits = _message_passing_init(H, v, missing,
        max_iter, :A, 2, chn_inits)
    return _message_passing(H_Int, w, chn_inits, _Gallager_A_check_node_message,
        var_adj_list, check_adj_list, max_iter, :A, 0, 0.0, schedule)
end

# TODO: threshold in docstring
"""
    Gallager_B(H::T, v::T; max_iter::Int = 100, threshold::Int = 2, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where T <: CTMatrixTypes

Run the Gallager-B decoder with the parity-check matrix `H` and received vector `v`.
"""
function Gallager_B(H::T, v::T; max_iter::Int = 100, threshold::Int = 2, chn_inits::Union{Missing,
    Vector{Float64}} = missing, schedule::Symbol = :flooding) where T <: CTMatrixTypes
    
    H_Int, w, var_adj_list, check_adj_list, chn_inits = _message_passing_init(H, v, missing,
        max_iter, :B, threshold, chn_inits)
    return _message_passing(H_Int, w, chn_inits, _Gallager_B_check_node_message,
        var_adj_list, check_adj_list, max_iter, :B, threshold, 0.0, schedule)
end

#############################
        # Sum product
#############################

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

"""
    sum_product(H::S, v::T, chn::MPNoiseModel; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where {S <: CTMatrixTypes, T <: Union{Vector{<:Real}, CTMatrixTypes}}

Run the sum-product algorithm with the parity-check matrix `H`, received vector `v`, and channel
`chn`.

# Notes
* Use `chn_inits` to pass in soft information.
* The options for `schedule` are `:flooding` (parallel) or `:serial`.
"""
function sum_product(H::S, v::T, chn::MPNoiseModel; max_iter::Int = 100, chn_inits::Union{Missing,
    Vector{Float64}} = missing, schedule::Symbol = :flooding) where {S <: CTMatrixTypes, T <:
    Union{Vector{<:Real}, CTMatrixTypes}}

    # TODO: need to type check H and v and for sizes
    # TODO: check max_iter
    # TODO: check schedule
    # TODO: check chn_inits

    H_Int, w, var_adj_list, check_adj_list, chn_inits = _message_passing_init(H, v, chn, max_iter,
        :SP, 2, chn_inits)
    return _message_passing(H_Int, w, chn_inits, _SP_check_node_message, var_adj_list,
        check_adj_list, max_iter, :SP, 0, 0.0, schedule)
end

"""
    sum_product_box_plus(H::S, v::T, chn::MPNoiseModel; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where {S <: CTMatrixTypes, T <: Union{Vector{<:Real}, CTMatrixTypes}}

Run the sum-product box-plus algorithm with the parity-check matrix `H`, received vector `v`, and
channel `chn`.

# Notes
* Use `chn_inits` to pass in soft information.
* The options for `schedule` are `:flooding` (parallel) or `:serial`.
"""
function sum_product_box_plus(H::S, v::T, chn::MPNoiseModel; max_iter::Int = 100,
    chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where {S <:
    CTMatrixTypes, T <: Union{Vector{<:Real}, CTMatrixTypes}}

    H_Int, w, var_adj_list, check_adj_list, chn_inits = _message_passing_init(H, v, chn, max_iter,
        :SP, 2, chn_inits)
    return _message_passing(H_Int, w, chn_inits, _SP_check_node_message_box_plus,
        var_adj_list, check_adj_list, max_iter, :SP, 0, 0.0, schedule)
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
    return _message_passing_decimation(H_Int, w, decimated_bits, chn_inits, _SP_check_node_message,
        var_adj_list, check_adj_list, max_iter, :SP, 0, 0.0, schedule)
end

"""
    sum_product_syndrome(H::S, syndrome::T, chn::MPNoiseModel; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where {S <: CTMatrixTypes, T <: Union{Vector{<:Real}, CTMatrixTypes}}

Run the syndrome-based sum-product algorithm with the parity-check matrix `H`, syndrome `syndrome`,
and channel `chn`.

# Notes
* Use `chn_inits` to pass in soft information.
* The options for `schedule` are `:flooding` (parallel) or `:serial`.
"""
function sum_product_syndrome(H::S, syndrome::T, chn::MPNoiseModel; max_iter::Int = 100,
    chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where {S <:
    CTMatrixTypes, T <: Union{Vector{<:Real}, CTMatrixTypes}}

    H_Int, syn, var_adj_list, check_adj_list, chn_inits = _message_passing_init_syndrome(H,
        syndrome, chn, max_iter, chn_inits, :SP)
    return _message_passing_syndrome(H_Int, syn, chn_inits, _SP_check_node_message,
        var_adj_list, check_adj_list, max_iter, :SP, 0, 0.0, schedule)
end

# syndrome-based decimation

#############################
        # Min Sum
#############################

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

"""
    min_sum(H::S, v::T, chn::MPNoiseModel; max_iter::Int = 100, attenuation::Float64 = 0.5, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where {S <: CTMatrixTypes, T <: Vector{<:AbstractFloat}}

Run the min-sum algorithm with the parity-check matrix `H`, received vector `v`, and channel
`chn`.

# Notes
* Use `chn_inits` to pass in soft information.
* The options for `schedule` are `:flooding` (parallel) or `:serial`.
* Set the normalization constant with `attenuation`.
"""
function min_sum(H::S, v::T, chn::MPNoiseModel; max_iter::Int = 100, attenuation::Float64 =
    0.5, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where
    {S <: CTMatrixTypes, T <: Vector{<:AbstractFloat}}

    # TODO: check this init call
    H_Int, w, var_adj_list, check_adj_list, chn_inits = _message_passing_init(H, v, chn, max_iter,
        :MS, 2, chn_inits)
    return _message_passing(H_Int, w, chn_inits, _MS_check_node_message,
        var_adj_list, check_adj_list, max_iter, :MS, 0, attenuation, schedule)
end

"""
    min_sum_syndrome(H::S, syndrome::T, chn::MPNoiseModel; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where {S <: CTMatrixTypes, T <: Union{Vector{<:Real}, CTMatrixTypes}}

Run the syndrome-based min-sum algorithm with the parity-check matrix `H`, syndrome `syndrome`,
and channel `chn`.

# Notes
* Use `chn_inits` to pass in soft information.
* The options for `schedule` are `:flooding` (parallel) or `:serial`.
* Set the normalization constant with `attenuation`.
"""
function min_sum_syndrome(H::S, syndrome::T, chn::MPNoiseModel; max_iter::Int = 100,
    chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where {S <:
    CTMatrixTypes, T <: Union{Vector{<:Real}, CTMatrixTypes}}

    H_Int, syn, var_adj_list, check_adj_list, chn_inits = _message_passing_init_syndrome(H,
        syndrome, chn, max_iter, chn_inits, :MS)
    return _message_passing_syndrome(H_Int, syn, chn_inits, _MS_check_node_message,
        var_adj_list, check_adj_list, max_iter, :MS, 0, 0.0, schedule)
end

# min-sum with decimation
# syndrome-based min-sum with decimation
# min-sum with correction (and variants)

#############################
       # Initialization
#############################

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
    
    # TODO: replace all FpmattoJulia calls in this file with lift from up-deps branch
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

#############################
      # Message Passing
#############################

# TODO: speed-up by doing all these allocations in init?
function _message_passing(H::Matrix{UInt64}, w::Vector{T}, chn_inits::Union{Missing,
    Vector{Float64}}, c_to_v_mess::Function, var_adj_list::Vector{Vector{Any}},
    check_adj_list::Vector{Vector{Any}}, max_iter::Int, kind::Symbol, Bt::Int = 2,
    attenuation::Float64 = 0.5, schedule::Symbol) where T <: Union{Int, AbstractFloat}

    # the inclusion of the kind statements add less than a microsecond
    num_check, num_var = size(H)
    # could reduce this stuff to UInt8 if needed
    S = kind ∈ (:A, :B) ? Int : Float64
    current_bits = zeros(Int, num_var)
    if kind in (:SP, :MS)
        totals = zeros(S, num_var)
    end
    syn = zeros(Int, num_check)

    if schedule == :flooding
        check_to_var_messages = zeros(S, num_check, num_var, 2)
        var_to_check_messages = zeros(S, num_var, num_check, 2)
    else
        check_to_var_messages = zeros(S, num_check, num_var, 1)
        var_to_check_messages = zeros(S, num_var, num_check, 1)
    end
    
    # first iteration for variable nodes - set to channel initialization
    if kind ∈ (:SP, :MS)
        @inbounds for v in 1:num_var
            var_to_check_messages[v, var_adj_list[v], 1] .= chn_inits[v]
        end
    elseif kind ∈ (:A, :B)
        # TODO: remove and set this to chn_inits
        @inbounds for v in 1:num_var
            var_to_check_messages[v, var_adj_list[v], :] .= w[v]
        end
    end

    iter = 1
    curr_iter = 1
    schedule == :flooding ? (prev_iter = 2;) : (prev_iter = 1;)
    # does this propagate correctly?
    @inbounds while iter ≤ max_iter
        # variable node is already done for first iteration, so start with check nodes
        @simd for c in 1:num_check
            for v in check_adj_list[c]
                check_to_var_messages[c, v, curr_iter] = c_to_v_mess(c, v, curr_iter, check_adj_list,
                    var_to_check_messages, attenuation)
            end
        end

        # one full iteration done, check if converged
        if kind ∈ (:SP, :MS)
            @simd for v in 1:num_var
                totals[v] = chn_inits[v]
                for c in var_adj_list[v]
                    totals[v] += check_to_var_messages[c, v, curr_iter]
                end
                current_bits[v] = totals[v] >= 0 ? 0 : 1
            end
        elseif kind ∈ (:A, :B)
            @simd for v in 1:num_var
                len = length(var_adj_list[v])
                one_count = count(isone, view(check_to_var_messages, var_adj_list[v], v, curr_iter))
                d = fld(len, 2)
                current_bits[v] = one_count + (isone(w[v]) && iseven(len)) > d
            end
        end

        LinearAlgebra.mul!(syn, H, current_bits)
        iszero(syn .% 2) && return true, current_bits, iter, var_to_check_messages,
            check_to_var_messages
        iter += 1
        if schedule == :flooding
            temp = curr_iter
            curr_iter = prev_iter
            prev_iter = temp
        end

        if iter <= max_iter
            for v in 1:num_var
                for c in var_adj_list[v]
                    if kind ∈ (:SP, :MS)
                        # this includes the channel inputs in total
                        var_to_check_messages[v, c, curr_iter] = totals[v] -
                            check_to_var_messages[c, v, prev_iter]
                    elseif kind == :A && length(var_adj_list[v]) > 1
                        if all(!Base.isequal(w[v]), check_to_var_messages[c2, v, prev_iter] for
                            c2 in var_adj_list[v] if c != c2)

                            var_to_check_messages[v, c, curr_iter] ⊻= 1
                        end
                    elseif kind == :B && length(var_adj_list[v]) >= Bt
                        if count(!Base.isequal(w[v]), check_to_var_messages[c2, v, prev_iter] for 
                            c2 in var_adj_list[v] if c != c2) >= Bt

                            var_to_check_messages[v, c, curr_iter] ⊻= 1
                        end
                    end
                end
            end
        end
    end

    return false, current_bits, iter, var_to_check_messages, check_to_var_messages
end

# merged into above
# function _message_passing_serial(H::Matrix{UInt64}, w::Vector{T}, chn_inits::Union{Missing,
#     Vector{Float64}}, c_to_v_mess::Function, var_adj_list::Vector{Vector{Any}},
#     check_adj_list::Vector{Vector{Any}}, max_iter::Int, kind::Symbol, Bt::Int = 2,
#     attenuation::Float64 = 0.5) where T <: Union{Int, AbstractFloat}

#      # the inclusion of the kind statements add less than a microsecond
#      num_check, num_var = size(H)
#      # could reduce this stuff to UInt8 if needed
#      S = kind ∈ (:A, :B) ? Int : Float64
#      current_bits = zeros(Int, num_var)
#      if kind in (:SP, :MS)
#          totals = zeros(S, num_var)
#      end
#      syn = zeros(Int, num_check)
#      # max_iter += 1 # probably should copy this
#      check_to_var_messages = zeros(S, num_check, num_var, 1)
#      var_to_check_messages = zeros(S, num_var, num_check, 1)
    
#     # first iteration for variable nodes - set to channel initialization
#     if kind ∈ (:SP, :MS)
#         @inbounds for v in 1:num_var
#             var_to_check_messages[v, var_adj_list[v], 1] .= chn_inits[v]
#         end
#     elseif kind ∈ (:A, :B)
#         @inbounds for v in 1:num_var
#             var_to_check_messages[v, var_adj_list[v], 1] .= w[v]
#         end
#     end

#     iter = 1
#     # does this propagate correctly?
#     @inbounds while iter ≤ max_iter
#         for c in 1:num_check
#             for v in check_adj_list[c]
#                 check_to_var_messages[c, v, 1] = c_to_v_mess(c, v, 1,
#                     check_adj_list, var_to_check_messages, attenuation)
#             end
#         end

#         # one full iteration done, check if converged
#         if kind ∈ (:SP, :MS)
#             @simd for v in 1:num_var
#                 totals[v] = chn_inits[v]
#                 for c in var_adj_list[v]
#                     totals[v] += check_to_var_messages[c, v, 1]
#                 end
#                 current_bits[v] = totals[v] >= 0 ? 0 : 1
#             end
#         elseif kind ∈ (:A, :B)
#             @simd for v in 1:num_var
#                 len = length(var_adj_list[v])
#                 one_count = count(isone, view(check_to_var_messages, var_adj_list[v], v, 1))
#                 d = fld(len, 2)
#                 current_bits[v] = one_count + (isone(w[v]) && iseven(len)) > d
#             end
#         end

#         LinearAlgebra.mul!(syn, H, current_bits)
#         iszero(syn .% 2) && return true, current_bits, iter, var_to_check_messages, check_to_var_messages
#         iter += 1

#         if iter <= max_iter
#             for v in 1:num_var
#                 for c in var_adj_list[v]
#                     if kind ∈ (:SP, :MS)
#                         # this includes the channel inputs in total
#                         var_to_check_messages[v, c, 1] = totals[v] -
#                             check_to_var_messages[c, v, 1]
#                     elseif kind == :A && length(var_adj_list[v]) > 1
#                         if all(!Base.isequal(w[v]), check_to_var_messages[c2, v, 1] for
#                             c2 in var_adj_list[v] if c != c2)

#                             var_to_check_messages[v, c, 1] ⊻= 1 
#                         end
#                     elseif kind == :B && length(var_adj_list[v]) >= Bt
#                         if count(!Base.isequal(w[v]), check_to_var_messages[c, v, 1] for 
#                             c2 in var_adj_list[v] if c != c2) >= Bt

#                             var_to_check_messages[v, c, 1] ⊻= 1 
#                         end
#                     end
#                 end
#             end
#         end
#     end
#     return false, current_bits, iter, var_to_check_messages, check_to_var_messages
# end
 
function _message_passing_syndrome(H::Matrix{UInt64}, syndrome::Vector{T},
    chn_inits::Union{Missing, Vector{Float64}}, c_to_v_mess::Function, 
    var_adj_list::Vector{Vector{Any}}, check_adj_list::Vector{Vector{Any}}, max_iter::Int,
    kind::Symbol, Bt::Int = 2, attenuation::Float64 = 0.5, schedule::Symbol) where T <: Union{Int,
    AbstractFloat}

    # the inclusion of the kind statements add less than a microsecond
    num_check, num_var = size(H)
    # could reduce this stuff to UInt8 if needed
    S = kind ∈ (:A, :B) ? Int : Float64
    current_bits = zeros(Int, num_var)
    if kind in (:SP, :MS)
        totals = zeros(S, num_var)
    end
    syn = zeros(Int, num_check)

    if schedule == :flooding
        check_to_var_messages = zeros(S, num_check, num_var, 2)
        var_to_check_messages = zeros(S, num_var, num_check, 2)
    else
        check_to_var_messages = zeros(S, num_check, num_var, 1)
        var_to_check_messages = zeros(S, num_var, num_check, 1)
    end
    
    # first iteration for variable nodes - set to channel initialization
    if kind ∈ (:SP, :MS)
        @inbounds for v in 1:num_var
            var_to_check_messages[v, var_adj_list[v], 1] .= chn_inits[v]
        end
    elseif kind ∈ (:A, :B)
        # TODO: remove and set this to chn_inits
        @inbounds for v in 1:num_var
            var_to_check_messages[v, var_adj_list[v], :] .= w[v]
        end
    end

    iter = 1
    curr_iter = 1
    schedule == :flooding ? (prev_iter = 2;) : (prev_iter = 1;)
    # does this propagate correctly?
    @inbounds while iter ≤ max_iter
        # variable node is already done for first iteration, so start with check nodes
        for c in 1:num_check
            count = sum(current_bits[check_adj_list[c]])
            for v in check_adj_list[c]
                check_to_var_messages[c, v, curr_iter] = 0.0
                count_v = count - current_bits[v]
                # here I want to loop through all possible options of errors
                for e in 0:1
                    # check if combination of v = e with the other bits on this check gives the
                    # correct syndrome bit
                    if (count_v + e) % 2 == syndrome[c]
                        # if it produces the syndrome bit, it contributes to Pr(e | s);
                        # otherwise, Pr(e | s) = 0
                        check_to_var_messages[c, v, curr_iter] += (-1)^e * c_to_v_mess(c, v, curr_iter,
                            check_adj_list, var_to_check_messages, attenuation)
                    end
                end
            end
        end

        # one full iteration done, check if converged
        if kind ∈ (:SP, :MS)
            @simd for v in 1:num_var
                totals[v] = chn_inits[v]
                for c in var_adj_list[v]
                    totals[v] += check_to_var_messages[c, v, curr_iter]
                end
                current_bits[v] = totals[v] >= 0 ? 0 : 1
            end
        elseif kind ∈ (:A, :B)
            @simd for v in 1:num_var
                len = length(var_adj_list[v])
                one_count = count(isone, view(check_to_var_messages, var_adj_list[v], v, curr_iter))
                d = fld(len, 2)
                current_bits[v] = one_count + (isone(w[v]) && iseven(len)) > d
            end
        end

        LinearAlgebra.mul!(syn, H, current_bits)
        syn .% 2 == syndrome && return true, current_bits, iter, var_to_check_messages, check_to_var_messages
        iter += 1
        if schedule == :flooding
            temp = curr_iter
            curr_iter = prev_iter
            prev_iter = temp
        end

        if iter <= max_iter
            for v in 1:num_var
                for c in var_adj_list[v]
                    if kind ∈ (:SP, :MS)
                        # this includes the channel inputs in total
                        var_to_check_messages[v, c, curr_iter] = totals[v] -
                            check_to_var_messages[c, v, prev_iter]
                    elseif kind == :A && length(var_adj_list[v]) > 1
                        if all(!Base.isequal(w[v]), check_to_var_messages[c2, v, prev_iter] for
                            c2 in var_adj_list[v] if c != c2)

                            var_to_check_messages[v, c, curr_iter] ⊻= 1
                        end
                    elseif kind == :B && length(var_adj_list[v]) >= Bt
                        if count(!Base.isequal(w[v]), check_to_var_messages[c2, v, prev_iter] for 
                            c2 in var_adj_list[v] if c != c2) >= Bt

                            var_to_check_messages[v, c, curr_iter] ⊻= 1
                        end
                    end
                end
            end
        end
    end

    return false, current_bits, iter, var_to_check_messages, check_to_var_messages
end

# # this is the working version from master to benchmark against
# function _message_passing_old(H::Matrix{UInt64}, w::Vector{T}, chn::Union{Missing, MPNoiseModel},
#     c_to_v_mess::Function, var_adj_list::Vector{Vector{Any}}, check_adj_list::Vector{Vector{Any}},
#     max_iter::Int, kind::Symbol, Bt::Int = 2, attenuation::Float64 = 0.5) where T <: Union{Int,
#     AbstractFloat}

#     num_check, num_var = size(H)
#     S = kind ∈ (:A, :B) ? Int : Float64
#     curr = zeros(Int, num_var)
#     if kind in (:SP, :MS)
#         totals = zeros(S, 1, num_var)
#     end
#     syn = zeros(Int, num_check)
#     max_iter += 1 # probably should copy this
#     check_to_var_messages = zeros(S, num_check, num_var, max_iter)
#     var_to_check_messages = zeros(S, num_var, num_check, max_iter)
    
#     iter = 1
#     if kind in (:SP, :MS)
#         chn_inits = if chn.type == :BSC
#             _channel_init_BSC(w, chn.cross_over_prob)
#         elseif chn.type == :BAWGNC && kind == :SP
#             _channel_init_BAWGNC_SP(w, chn.sigma)
#         elseif chn.type == :BAWGNC && kind == :MS
#             _channel_init_BAWGNC_MS(w)
#         end
#         for vn in 1:num_var
#             var_to_check_messages[vn, var_adj_list[vn], 1] .= chn_inits[vn]
#         end
#     elseif kind in (:A, :B)
#         for vn in 1:num_var
#             var_to_check_messages[vn, var_adj_list[vn], :] .= w[vn]
#         end
#     end

#     while iter < max_iter
#         for cn in 1:num_check
#             for v1 in check_adj_list[cn]
#                 check_to_var_messages[cn, v1, iter] = c_to_v_mess(cn, v1, iter, check_adj_list,
#                     var_to_check_messages, attenuation)
#             end
#         end

#         if kind in (:SP, :MS)
#             for vn in 1:num_var
#                 totals[vn] = chn_inits[vn]
#                 for c in var_adj_list[vn]
#                     totals[vn] += check_to_var_messages[c, vn, iter]
#                 end
#             end
#         end

#         if kind in (:SP, :MS)
#             @simd for i in 1:num_var
#                 curr[i] = totals[i] >= 0 ? 0 : 1
#             end
#         elseif kind in (:A, :B)
#             @simd for i in 1:num_var
#                 len = length(var_adj_list[i])
#                 one_count = count(isone, view(check_to_var_messages, var_adj_list[i], i, iter))
#                 d = fld(len, 2)
#                 curr[i] = one_count + (isone(w[i]) && iseven(len)) > d
#             end
#         end

#         LinearAlgebra.mul!(syn, H, curr)
#         # @show curr
#         # @show syn .% 2
#         iszero(syn .% 2) && return true, curr, iter, var_to_check_messages, check_to_var_messages
#         iter += 1

#         if iter <= max_iter
#             for vn in 1:num_var
#                 for c1 in var_adj_list[vn]
#                     if kind in (:SP, :MS)
#                         # this includes the channel inputs in total
#                         var_to_check_messages[vn, c1, iter] = totals[vn] -
#                             check_to_var_messages[c1, vn, iter - 1]
#                     elseif kind == :A && length(var_adj_list[vn]) > 1
#                         if all(!Base.isequal(w[vn]), check_to_var_messages[c2, vn, iter - 1] for c2 
#                             in var_adj_list[vn] if c1 != c2)

#                             var_to_check_messages[vn, c1, iter] ⊻= 1 
#                         end
#                     elseif kind == :B && length(var_adj_list[vn]) >= Bt
#                         if count(!Base.isequal(w[vn]), check_to_var_messages[c2, vn, iter - 1] for 
#                             c2 in var_adj_list[vn] if c1 != c2) >= Bt

#                             var_to_check_messages[vn, c1, iter] ⊻= 1 
#                         end
#                     end
#                 end
#             end
#         end
#     end

#     return false, curr, iter, var_to_check_messages, check_to_var_messages
# end


#############################
          # Methods
#############################

# Mansour, Shanbhag, "Turbo Decoder Architectures for Low-Density Parity-Check Codes" (2002)

# TODO: figure out and describe what this actually does - probably layered/semi-serial schedule
# rename and add docstring
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
    for c in 2:num_check
        found = false
        for sched in sched_list
            if !any(x ∈ check_adj_list[y] for y in sched for x ∈ check_adj_list[c])
                push!(sched, c)
                sort!(sched_list, lt = (x, y) -> length(x) < length(y))
                found = true
                break
            end
        end
        !found && push!(sched_list, [c])
    end
    
    return sched_list
end

# Example of using Gallager A and B (out should end up [1 1 1 0 0 0 0]
# H = matrix(GF(2), [1 1 0 1 1 0 0; 1 0 1 1 0 1 0; 0 1 1 1 0 0 1]);
# v = matrix(GF(2), 7, 1, [1, 1, 0, 0, 0, 0, 0]);
# flag, out, iter, vtoc, ctov = Gallager_A(H, v, 100);
# flag, out, iter, vtoc, ctov = Gallager_B(H, v, 100);
# nm = MPNoiseModel(:BSC, 1/7);
# flag, out, iter, vtoc, ctov = sum_product(H, v, nm, 100);
