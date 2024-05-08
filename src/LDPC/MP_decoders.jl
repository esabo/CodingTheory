# Copyright (c) 2023 - 2024 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
         # Gallager
#############################

function _Gallager_A_check_node_message(c::Int, v::Int, iter::Int, check_adj_list::Vector{Vector{Int}},
    var_to_check_messages::Array{Int, 3}, attenuation::Float64)

    @inbounds reduce(⊻, var_to_check_messages[v2, c, iter] for v2 in check_adj_list[c] if v2 != v)
end
_Gallager_B_check_node_message(c::Int, v::Int, iter::Int, check_adj_list::Vector{Vector{Int}}, var_to_check_messages,
    attenuation::Float64) = _Gallager_A_check_node_message(c, v, iter, check_adj_list,
    var_to_check_messages, attenuation)

"""
    Gallager_A(H::T, v::T; max_iter::Int = 100, schedule::Symbol = :parallel) where T <: CTMatrixTypes

Run the Gallager-A decoder with the parity-check matrix `H` and received vector `v`.
"""
function Gallager_A(H::T, v::T; max_iter::Int = 100, schedule::Symbol = :parallel) where T <:
    CTMatrixTypes

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(v) ≠ (nc, 1) && size(v) ≠ (1, nc)) && throw(ArgumentError("Vector has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial) || throw(ArgumentError("Unknown schedule algorithm"))
    schedule == :parallel && (schedule = :flooding;)

    # initialization - do these outside to reduce allocations when looped
    H_Int, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages, var_to_check_messages,
        current_bits, syn = _message_passing_init_Int(H, v, max_iter, :A, 2, schedule, Int[])

    return _message_passing_Int(H_Int, missing, chn_inits_2, _Gallager_A_check_node_message,
        var_adj_list, check_adj_list, max_iter, :A, schedule, current_bits, syn,
        check_to_var_messages, var_to_check_messages, 0)
end

# TODO: threshold in docstring
"""
    Gallager_B(H::T, v::T; max_iter::Int = 100, threshold::Int = 2, schedule::Symbol = :parallel) where T <: CTMatrixTypes

Run the Gallager-B decoder with the parity-check matrix `H` and received vector `v`.
"""
function Gallager_B(H::T, v::T; max_iter::Int = 100, threshold::Int = 2, schedule::Symbol =
    :flooding) where T <: CTMatrixTypes
    
    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(v) ≠ (nc, 1) && size(v) ≠ (1, nc)) && throw(ArgumentError("Vector has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial) || throw(ArgumentError("Unknown schedule algorithm"))
    schedule == :parallel && (schedule = :flooding;)

    # initialization - do these outside to reduce allocations when looped
    H_Int, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages, var_to_check_messages,
        current_bits, syn = _message_passing_init_Int(H, v, max_iter, :B, threshold, schedule,
        Int[])

    return _message_passing_Int(H_Int, missing, chn_inits_2, _Gallager_B_check_node_message,
        var_adj_list, check_adj_list, max_iter, :B, schedule, current_bits, syn,
        check_to_var_messages, var_to_check_messages, threshold)
end

#############################
        # Sum product
#############################

function _SP_check_node_message(c::Int, v::Int, iter::Int, check_adj_list::Vector{Vector{Int}}, var_to_check_messages::Array{Float64, 3},
    attenuation::Float64)

    ϕ(x) = -log(tanh(0.5 * x))
    temp = 0.0
    s = 1
    @inbounds for v2 in check_adj_list[c]
        if v2 != v
            x = var_to_check_messages[v2, c, iter]
            if x >= 0
                temp += ϕ(x)
            else
                temp += ϕ(-x)
                s *= -1
            end
        end
    end
    return s * ϕ(temp)
end

⊞(a::Float64, b::Float64) = log((1 + exp(a + b)) / (exp(a) + exp(b)))
⊞(a...) = reduce(⊞, a...)
function _SP_check_node_message_box_plus(c::Int, v::Int, iter::Int, check_adj_list::Vector{Vector{Int}},
    var_to_check_messages::Array{Float64, 3}, attenuation::Float64)

    @inbounds ⊞(var_to_check_messages[v2, c, iter] for v2 in check_adj_list[c] if v2 != v)
end

"""
    sum_product(H::T, v::T, chn::MPNoiseModel; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :parallel, rand_sched::Bool = false, erasures::Vector{Int} = Int[]) where T <: CTMatrixTypes

Run the sum-product algorithm with the parity-check matrix `H`, received vector `v`, and channel
`chn`.

# Notes
- Use `chn_inits` to pass in soft information.
- The options for `schedule` are `:parallel` (`:flooding`), `:serial`, or `:layered` (`:semiserial`).
"""
function sum_product(H::T, v::T, chn::MPNoiseModel; max_iter::Int = 100, chn_inits::Union{Missing,
    Vector{Float64}} = missing, schedule::Symbol = :parallel, rand_sched::Bool = false,
    erasures::Vector{Int} = Int[]) where T <: CTMatrixTypes

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(v) ≠ (nc, 1) && size(v) ≠ (1, nc)) && throw(ArgumentError("Vector has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial, :layered, :semiserial) || 
        throw(ArgumentError("Unknown schedule algorithm"))
    schedule == :flooding && (schedule = :parallel;)
    schedule == :semiserial && (schedule = :layered;)

    # initialization - do these outside to reduce allocations when looped
    layers = layered_schedule(H, schedule = schedule, random = rand_sched)
    H_Int, _, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages,
        var_to_check_messages, current_bits, totals, syn = _message_passing_init(H, v, chn,
        max_iter, :SP, chn_inits, schedule, erasures)

    # return _message_passing(H_Int, missing, chn_inits_2, _SP_check_node_message, var_adj_list,
    #     check_adj_list, max_iter, schedule, current_bits, totals, syn, check_to_var_messages,
    #     var_to_check_messages, 0.0)

    return _message_passing_layered(H_Int, missing, chn_inits_2, _SP_check_node_message, 
        var_adj_list, check_adj_list, max_iter, schedule, current_bits, totals, syn,
        check_to_var_messages, var_to_check_messages, 0.0, layers)
end

"""
    sum_product_box_plus(H::T, v::T, chn::MPNoiseModel; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :parallel, rand_sched::Bool = false, erasures::Vector{Int} = Int[]) where T <: CTMatrixTypes

Run the sum-product box-plus algorithm with the parity-check matrix `H`, received vector `v`, and
channel `chn`.

# Notes
- Use `chn_inits` to pass in soft information.
- The options for `schedule` are `:parallel` (`:flooding`), `:serial`, or `:layered` (`:semiserial`).
"""
function sum_product_box_plus(H::T, v::T, chn::MPNoiseModel; max_iter::Int = 100,
    chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :parallel,
    rand_sched::Bool = false, erasures::Vector{Int} = Int[]) where T <: CTMatrixTypes

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(v) ≠ (nc, 1) && size(v) ≠ (1, nc)) && throw(ArgumentError("Vector has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial, :layered, :semiserial) || 
        throw(ArgumentError("Unknown schedule algorithm"))
    schedule == :flooding && (schedule = :parallel;)
    schedule == :semiserial && (schedule = :layered;)

    # initialization - do these outside to reduce allocations when looped
    layers = layered_schedule(H, schedule = schedule, random = rand_sched)
    H_Int, _, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages,
        var_to_check_messages, current_bits, totals, syn = _message_passing_init(H, v, chn,
        max_iter, :SP, chn_inits, schedule, erasures)

    # if schedule != :layered
    #     return _message_passing(H_Int, missing, chn_inits_2, _SP_check_node_message_box_plus,
    #         var_adj_list, check_adj_list, max_iter, schedule, current_bits, totals, syn,
    #         check_to_var_messages, var_to_check_messages, 0.0)
    # else
        return _message_passing_layered(H_Int, missing, chn_inits_2, 
            _SP_check_node_message_box_plus, var_adj_list, check_adj_list, max_iter, schedule,
            current_bits, totals, syn, check_to_var_messages, var_to_check_messages, 0.0, layers)

    # end
end

"""
    sum_product_syndrome(H::T, syndrome::T, chn::MPNoiseModel; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :parallel, rand_sched::Bool = false, erasures::Vector{Int} = Int[]) where T <: CTMatrixTypes

Run the syndrome-based sum-product algorithm with the parity-check matrix `H`, syndrome `syndrome`,
and channel `chn`.

# Notes
- Use `chn_inits` to pass in soft information.
- The options for `schedule` are `:parallel` (`:flooding`), `:serial`, or `:layered` (`:semiserial`).
"""
function sum_product_syndrome(H::T, syndrome::T, chn::MPNoiseModel; max_iter::Int = 100,
    chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :parallel,
    rand_sched::Bool = false, erasures::Vector{Int} = Int[]) where T <: CTMatrixTypes

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(syndrome) ≠ (nr, 1) && size(syndrome) ≠ (1, nr)) && throw(ArgumentError("Syndrome has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial, :layered, :semiserial) || 
        throw(ArgumentError("Unknown schedule algorithm"))
    schedule == :flooding && (schedule = :parallel;)
    schedule == :semiserial && (schedule = :layered;)

    # initialization - do these outside to reduce allocations when looped
    layers = layered_schedule(H, schedule = schedule, random = rand_sched)
    H_Int, _, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages,
        var_to_check_messages, current_bits, totals, syn = _message_passing_init(H,
        syndrome, chn, max_iter, :SP, chn_inits, schedule, erasures)
    syn_Int = _Flint_matrix_to_Julia_int_vector(syndrome)

    # return _message_passing(H_Int, syn_Int, chn_inits_2, _SP_check_node_message,
    #     var_adj_list, check_adj_list, max_iter, schedule, current_bits, totals, syn, 
    #     check_to_var_messages, var_to_check_messages, 0.0)
    return _message_passing_layered(H_Int, syn_Int, chn_inits_2, _SP_check_node_message,
        var_adj_list, check_adj_list, max_iter, schedule, current_bits, totals, syn, 
        check_to_var_messages, var_to_check_messages, 0.0, layers)
end

# believe this can be merged into an optional argument of the above but keep for now so as not to break working code
function sum_product_decimation(H::T, v::T, chn::MPNoiseModel, algorithm::Symbol;
    decimated_bits_values::Vector{Tuple{Int, S}} = Tuple{Int, S}[], max_iter::Int = 100,
    chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :parallel,
    rand_sched::Bool = false, guided_rounds::Int = 10, erasures::Vector{Int} = Int[]) where {T <:
    CTMatrixTypes, S <: CTFieldElem}

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(v) ≠ (nc, 1) && size(v) ≠ (1, nc)) && throw(ArgumentError("Vector has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial, :layered, :semiserial) || 
        throw(ArgumentError("Unknown schedule algorithm"))
    schedule == :flooding && (schedule = :parallel;)
    schedule == :semiserial && (schedule = :layered;)
    algorithm ∈ (:auto, :manual, :guided) || throw(ArgumentError("Unknown decimation algorithm"))
    (algorithm == :manual && !isempty(decimated_bits_values)) ||
        throw(ArgumentError("Manual decimation but no decimated bits and values provided"))
    # unclear how to interpret passed in values if auto or guided is set, so ignore
    (algorithm == :auto || algorithm == :guided) && (decimated_bits_values = Tuple{Int, S}[];)
    if algorithm == :guided
        guided_rounds > 0 || throw(DomainError("The number of rounds before decimation must be positive"))
    end

    # initialization - do these outside to reduce allocations when looped
    layers = layered_schedule(H, schedule = schedule, random = rand_sched)
    H_Int, w, var_adj_list, check_adj_list, chn_inits_2, decimated_bits, decimated_values,
        check_to_var_messages, var_to_check_messages, current_bits, totals, syn =
        _message_passing_init_decimation(H, v, chn, decimated_bits_values, max_iter, :SP, 2,
        chn_inits, schedule, erasures)

    # TODO layers
    return _message_passing_decimation(H_Int, w, chn_inits_2, _SP_check_node_message,
        var_adj_list, check_adj_list, max_iter, :SP, schedule, decimated_bits, decimated_values,
        current_bits, totals, syn, check_to_var_messages, var_to_check_messages, 0, 0.0, algorithm,
        guided_rounds)
end

# decimation, guided decimation, automatic decimation
# ensemble decoding
# syndrome-based decimation

#############################
        # Min Sum
#############################

box_plus_min(a::Float64, b::Float64) = sign(a) * sign(b) * min(abs(a), abs(b))
box_plus_min(a...) = reduce(box_plus_min, a...)
function _MS_check_node_message(c::Int, v::Int, iter::Int, check_adj_list::Vector{Vector{Int}},
    var_to_check_messages::Array{Float64, 3}, attenuation::Float64)

    @inbounds box_plus_min(var_to_check_messages[v2, c, iter] for v2 in check_adj_list[c] if
        v2 != v)
end

"""
    min_sum(H::T, v::T, chn::MPNoiseModel; max_iter::Int = 100, attenuation::Float64 = 0.5, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :parallel, rand_sched::Bool = false, erasures::Vector{Int} = Int[]) where T <: CTMatrixTypes

Run the min-sum algorithm with the parity-check matrix `H`, received vector `v`, and channel
`chn`.

# Notes
- Use `chn_inits` to pass in soft information.
- The options for `schedule` are `:parallel` (`:flooding`), `:serial`, or `:layered` (`:semiserial`).
- Set the normalization constant with `attenuation`.
"""
function min_sum(H::T, v::T, chn::MPNoiseModel; max_iter::Int = 100, attenuation::Float64 =
    0.5, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :parallel,
    rand_sched::Bool = false, erasures::Vector{Int} = Int[]) where T <: CTMatrixTypes

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(v) ≠ (nc, 1) && size(v) ≠ (1, nc)) && throw(ArgumentError("Vector has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial, :layered, :semiserial) || 
        throw(ArgumentError("Unknown schedule algorithm"))
    schedule == :flooding && (schedule = :parallel;)
    schedule == :semiserial && (schedule = :layered;)

    # initialization - do these outside to reduce allocations when looped
    layers = layered_schedule(H, schedule = schedule, random = rand_sched)
    H_Int, _, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages,
        var_to_check_messages, current_bits, totals, syn = _message_passing_init(H, v, chn,
        max_iter, :MS, chn_inits, schedule, erasures)

    # return _message_passing(H_Int, missing, chn_inits_2, _MS_check_node_message,
    #     var_adj_list, check_adj_list, max_iter, schedule, current_bits, totals, syn,
    #     check_to_var_messages, var_to_check_messages, attenuation)
    return _message_passing_layered(H_Int, missing, chn_inits_2, _MS_check_node_message,
        var_adj_list, check_adj_list, max_iter, schedule, current_bits, totals, syn,
        check_to_var_messages, var_to_check_messages, attenuation, layers)
end

"""
    min_sum_syndrome(H::T, syndrome::T, chn::MPNoiseModel; max_iter::Int = 100, attenuation::Float64 = 0.5, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :parallel, rand_sched::Bool = false, erasures::Vector{Int} = Int[]) where T <: CTMatrixTypes

Run the syndrome-based min-sum algorithm with the parity-check matrix `H`, syndrome `syndrome`,
and channel `chn`.

# Notes
- Use `chn_inits` to pass in soft information.
- The options for `schedule` are `:parallel` (`:flooding`), `:serial`, or `:layered` (`:semiserial`).
- Set the normalization constant with `attenuation`.
"""
function min_sum_syndrome(H::T, syndrome::T, chn::MPNoiseModel; max_iter::Int = 100,
    attenuation::Float64 = 0.5, chn_inits::Union{Missing, Vector{Float64}} = missing,
    schedule::Symbol = :parallel, rand_sched::Bool = false, erasures::Vector{Int} = Int[]) where
    T <: CTMatrixTypes

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(syndrome) ≠ (nr, 1) && size(syndrome) ≠ (1, nr)) && throw(ArgumentError("Syndrome has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial, :layered, :semiserial) || 
        throw(ArgumentError("Unknown schedule algorithm"))
    schedule == :flooding && (schedule = :parallel;)
    schedule == :semiserial && (schedule = :layered;)

    # initialization - do these outside to reduce allocations when looped
    layers = layered_schedule(H, schedule = schedule, random = rand_sched)
    H_Int, _, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages,
        var_to_check_messages, current_bits, totals, syn = _message_passing_init(H,
        syndrome, chn, max_iter, :MS, chn_inits, schedule, erasures)
    syn_Int = _Flint_matrix_to_Julia_int_vector(syndrome)

    # return _message_passing(H_Int, syn_Int, chn_inits_2, _MS_check_node_message,
    #     var_adj_list, check_adj_list, max_iter, schedule, current_bits, totals, syn, 
    #     check_to_var_messages, var_to_check_messages, attenuation)
    return _message_passing_layered(H_Int, syn_Int, chn_inits_2, _MS_check_node_message,
        var_adj_list, check_adj_list, max_iter, schedule, current_bits, totals, syn, 
        check_to_var_messages, var_to_check_messages, attenuation, layers)
end

# believe this can be merged into an optional argument of the above but keep for now so as not to break working code
function min_sum_decimation(H::T, v::T, chn::MPNoiseModel, algorithm::Symbol;
    decimated_bits_values::Vector{Tuple{Int, S}} = Tuple{Int, S}[], max_iter::Int = 100,
    attenuation::Float64 = 0.5, chn_inits::Union{Missing, Vector{Float64}} = missing, 
    schedule::Symbol = :parallel, rand_sched::Bool = false, guided_rounds::Int = 10,
    erasures::Vector{Int} = Int[]) where {T <: CTMatrixTypes, S <: CTFieldElem}

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(v) ≠ (nc, 1) && size(v) ≠ (1, nc)) && throw(ArgumentError("Vector has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial, :layered, :semiserial) || 
        throw(ArgumentError("Unknown schedule algorithm"))
    schedule == :flooding && (schedule = :parallel;)
    schedule == :semiserial && (schedule = :layered;)
    algorithm ∈ (:auto, :manual, :guided) || throw(ArgumentError("Unknown decimation algorithm"))
    (algorithm == :manual && !isempty(decimated_bits_values)) ||
        throw(ArgumentError("Manual decimation but no decimated bits and values provided"))
    # unclear how to interpret passed in values if auto or guided is set, so ignore
    (algorithm == :auto || algorithm == :guided) && (decimated_bits_values = Tuple{Int, S}[];)
    if algorithm == :guided
        guided_rounds > 0 || throw(DomainError("The number of rounds before decimation must be positive"))
    end

    # initialization - do these outside to reduce allocations when looped
    layers = layered_schedule(H, schedule = schedule, random = rand_sched)
    H_Int, w, var_adj_list, check_adj_list, chn_inits_2, decimated_bits, decimated_values,
        check_to_var_messages, var_to_check_messages, current_bits, totals, syn =
        _message_passing_init_decimation(H, v, chn, decimated_bits_values, max_iter, :MS, 2,
        chn_inits, schedule, erasures)

    # TODO layers
    return _message_passing_decimation(H_Int, w, chn_inits_2, _MS_check_node_message,
        var_adj_list, check_adj_list, max_iter, :MS, schedule, decimated_bits, decimated_values,
        current_bits, totals, syn, check_to_var_messages, var_to_check_messages, 0, 0.0, algorithm,
        guided_rounds)
end

# syndrome-based min-sum with decimation

#############################
  # Min Sum With Correction
#############################

function _min_sum_corr_s(a::Float64, b::Float64)
    if abs(a + b) < 2 && abs(a - b) > 2 * abs(a + b)
        return 0.5
    elseif abs(a - b) < 2 && abs(a + b) > 2 * abs(a - b)
        return -0.5
    else
        return 0
    end
end
box_plus_min_c(a::Float64, b::Float64) = sign(a) * sign(b) * min(abs(a), abs(b)) + _min_sum_corr_s(a, b)
box_plus_min_c(a...) = reduce(box_plus_min_c, a...)
function _MS_correction_check_node_message(c::Int, v::Int, iter::Int, check_adj_list::Vector{Vector{Int}},
    var_to_check_messages::Array{Float64, 3}, attenuation::Float64)

    @inbounds box_plus_min_c(var_to_check_messages[v2, c, iter] for v2 in check_adj_list[c] if
        v2 != v)
end

"""
    min_sum_with_correction(H::T, v::T, chn::MPNoiseModel; max_iter::Int = 100, attenuation::Float64 = 0.5, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :parallel, rand_sched::Bool = false, erasures::Vector{Int} = Int[]) where T <: CTMatrixTypes

Run the min-sum algorithm with the parity-check matrix `H`, received vector `v`, and channel
`chn`.

# Notes
- Use `chn_inits` to pass in soft information.
- The options for `schedule` are `:parallel` (`:flooding`), `:serial`, or `:layered` (`:semiserial`).
- Set the normalization constant with `attenuation`.
- A low-complexity approximation to the correction term is used.
"""
function min_sum_with_correction(H::T, v::T, chn::MPNoiseModel; max_iter::Int = 100,
    attenuation::Float64 = 0.5, chn_inits::Union{Missing, Vector{Float64}} = missing,
    schedule::Symbol = :parallel, rand_sched::Bool = false, erasures::Vector{Int} = Int[]) where
    T <: CTMatrixTypes

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(v) ≠ (nc, 1) && size(v) ≠ (1, nc)) && throw(ArgumentError("Vector has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial, :layered, :semiserial) || 
        throw(ArgumentError("Unknown schedule algorithm"))
    schedule == :flooding && (schedule = :parallel;)
    schedule == :semiserial && (schedule = :layered;)

    # initialization - do these outside to reduce allocations when looped
    layers = layered_schedule(H, schedule = schedule, random = rand_sched)
    H_Int, _, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages,
        var_to_check_messages, current_bits, totals, syn = _message_passing_init(H, v, chn,
        max_iter, :MS, chn_inits, schedule, erasures)

    # return _message_passing(H_Int, missing, chn_inits_2, _MS_correction_check_node_message,
    #     var_adj_list, check_adj_list, max_iter, schedule, current_bits, totals, syn,
    #     check_to_var_messages, var_to_check_messages, attenuation)
    return _message_passing_layered(H_Int, missing, chn_inits_2, _MS_correction_check_node_message,
        var_adj_list, check_adj_list, max_iter, schedule, current_bits, totals, syn,
        check_to_var_messages, var_to_check_messages, attenuation, layers)
end

"""
    min_sum_correction_syndrome(H::T, syndrome::T, chn::MPNoiseModel; max_iter::Int = 100, attenuation::Float64 = 0.5, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :parallel, rand_sched::Bool = false, erasures::Vector{Int} = Int[]) where T <: CTMatrixTypes

Run the syndrome-based min-sum-with-correction algorithm with the parity-check matrix `H`, syndrome 
`syndrome`, and channel `chn`.

# Notes
- Use `chn_inits` to pass in soft information.
- The options for `schedule` are `:parallel` (`:flooding`), `:serial`, or `:layered` (`:semiserial`).
- Set the normalization constant with `attenuation`.
- A low-complexity approximation to the correction term is used.
"""
function min_sum_with_correction_syndrome(H::T, syndrome::T, chn::MPNoiseModel; max_iter::Int = 100,
    attenuation::Float64 = 0.5, chn_inits::Union{Missing, Vector{Float64}} = missing,
    schedule::Symbol = :parallel, rand_sched::Bool = false, erasures::Vector{Int} = Int[]) where
    T <: CTMatrixTypes

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(syndrome) ≠ (nr, 1) && size(syndrome) ≠ (1, nr)) && throw(ArgumentError("Syndrome has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial, :layered, :semiserial) || 
        throw(ArgumentError("Unknown schedule algorithm"))
    schedule == :flooding && (schedule = :parallel;)
    schedule == :semiserial && (schedule = :layered;)

    # initialization - do these outside to reduce allocations when looped
    layers = layered_schedule(H, schedule = schedule, random = rand_sched)
    H_Int, _, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages,
        var_to_check_messages, current_bits, totals, syn = _message_passing_init(H,
        syndrome, chn, max_iter, :MS, chn_inits, schedule, erasures)
    syn_Int = _Flint_matrix_to_Julia_int_vector(syndrome)

    # return _message_passing(H_Int, syn_Int, chn_inits_2, _MS_correction_check_node_message,
    #     var_adj_list, check_adj_list, max_iter, schedule, current_bits, totals, syn, 
    #     check_to_var_messages, var_to_check_messages, attenuation)
    return _message_passing_layered(H_Int, syn_Int, chn_inits_2, _MS_correction_check_node_message,
        var_adj_list, check_adj_list, max_iter, schedule, current_bits, totals, syn, 
        check_to_var_messages, var_to_check_messages, attenuation, layers)
end

function min_sum_correction_decimation(H::T, v::T, chn::MPNoiseModel, algorithm::Symbol;
    decimated_bits_values::Vector{Tuple{Int, S}} = Tuple{Int, S}[], max_iter::Int = 100,
    attenuation::Float64 = 0.5, chn_inits::Union{Missing, Vector{Float64}} = missing, 
    schedule::Symbol = :parallel, rand_sched::Bool = false, guided_rounds::Int = 10,
    erasures::Vector{Int} = Int[]) where {T <: CTMatrixTypes, S <: CTFieldElem}

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(v) ≠ (nc, 1) && size(v) ≠ (1, nc)) && throw(ArgumentError("Vector has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial, :layered, :semiserial) || 
        throw(ArgumentError("Unknown schedule algorithm"))
    schedule == :flooding && (schedule = :parallel;)
    schedule == :semiserial && (schedule = :layered;)
    algorithm ∈ (:auto, :manual, :guided) || throw(ArgumentError("Unknown decimation algorithm"))
    (algorithm == :manual && !isempty(decimated_bits_values)) ||
        throw(ArgumentError("Manual decimation but no decimated bits and values provided"))
    # unclear how to interpret passed in values if auto or guided is set, so ignore
    (algorithm == :auto || algorithm == :guided) && (decimated_bits_values = Tuple{Int, S}[];)
    if algorithm == :guided
        guided_rounds > 0 || throw(DomainError("The number of rounds before decimation must be positive"))
    end

    # initialization - do these outside to reduce allocations when looped
    layers = layered_schedule(H, schedule = schedule, random = rand_sched)
    H_Int, w, var_adj_list, check_adj_list, chn_inits_2, decimated_bits, decimated_values,
        check_to_var_messages, var_to_check_messages, current_bits, totals, syn =
        _message_passing_init_decimation(H, v, chn, decimated_bits_values, max_iter, :MS, 2,
        chn_inits, schedule, erasures)

    # TODO layers
    return _message_passing_decimation(H_Int, w, chn_inits_2, _MS_correction_check_node_message,
        var_adj_list, check_adj_list, max_iter, :MS, schedule, decimated_bits, decimated_values,
        current_bits, totals, syn, check_to_var_messages, var_to_check_messages, 0, attenuation,
        algorithm, guided_rounds)
end


#############################
       # Initialization
#############################

function _channel_init_BSC(v::Vector{<: Integer}, p::Float64)
    temp = log((1 - p) / p)
    chn_init = zeros(Float64, length(v))
    for i in eachindex(v)
        @inbounds chn_init[i] = (-1)^v[i] * temp
    end
    return chn_init
end

function _channel_init_BAWGNC_SP(v::Vector{<: AbstractFloat}, σ::Float64)
    temp = 2 / σ^2
    chn_init = zeros(Float64, length(v))
    for i in eachindex(v)
        @inbounds chn_init[i] = temp * v[i]
    end
    return chn_init
end

_channel_init_BAWGNC_MS(v::Vector{<: AbstractFloat}) = v

function _message_passing_init(H::Union{Matrix{S}, T}, v::Union{Vector{S}, Vector{AbstractFloat},
    T}, chn::MPNoiseModel, max_iter::Int, kind::Symbol, chn_inits::Union{Missing,
    Vector{Float64}}, schedule::Symbol, erasures::Vector{Int}) where {S <: Integer,
    T <: CTMatrixTypes}

    kind ∈ (:SP, :MS) || throw(ArgumentError("Unknown value for parameter kind"))
    # will this work?
    if isa(H, CTMatrixTypes)
        Int(order(base_ring(H))) == 2 ||
        throw(ArgumentError("Currently only implemented for binary codes"))
        H_Int = _Flint_matrix_to_Julia_int_matrix(H)
    else
        H_Int = H
    end
    num_check, num_var = size(H_Int)
    num_check > 0 && num_var > 0 || throw(ArgumentError("Input matrix of improper dimension"))
    2 <= max_iter || throw(DomainError("Number of maximum iterations must be at least two"))
    
    len_v = length(v)
    if len_v == num_var
        syndrome_based = false
        kind ∈ (:SP, :MS) && chn.type == :BAWGNC && !isa(v, Vector{<: AbstractFloat}) &&
        throw(DomainError("Received message should be a vector of floats for BAWGNC."))
        kind ∈ (:SP, :MS) && chn.type == :BSC && !isa(v, Vector{Int}) && !isa(v, CTMatrixTypes) &&
        throw(DomainError("Received message should be a vector of Ints for BSC."))
    elseif len_v == num_check
        syndrome_based = true
    else
        throw(ArgumentError("Vector has incorrect dimension"))
    end
    
    check_adj_list = [Int[] for _ in 1:num_check]
    var_adj_list = [Int[] for _ in 1:num_var]
    for r in 1:num_check
        for c in 1:num_var
            if !iszero(H_Int[r, c])
                push!(check_adj_list[r], c)
                push!(var_adj_list[c], r)
            end
        end
    end

    if !syndrome_based && ismissing(chn_inits)
        chn_inits_2 = if chn.type == :BSC
            typeof(v) <: CTMatrixTypes ? (w = _Flint_matrix_to_Julia_int_vector(v);) : (w = w;)
            _channel_init_BSC(w, chn.cross_over_prob)
        elseif chn.type == :BAWGNC && kind == :SP
            _channel_init_BAWGNC_SP(v, chn.sigma)
        elseif chn.type == :BAWGNC && kind == :MS
            _channel_init_BAWGNC_MS(v)
        end
    elseif syndrome_based && ismissing(chn_inits)
        chn_inits_2 = [log((1 - chn.cross_over_prob) / chn.cross_over_prob) for _ in 1:num_var]
    else
        # println(length(chn_inits), ", ", num_var)
        length(chn_inits) ≠ num_var && throw(ArgumentError("Channel inputs has wrong size"))
        chn_inits_2 = chn_inits
    end

    # could reduce this stuff to UInt8 if needed
    current_bits = zeros(Int, num_var)
    totals = zeros(Float64, num_var)
    syn = zeros(Int, num_check)

    # remove for combing everything into a single layered system
    # if schedule == :serial
    #     check_to_var_messages = zeros(Float64, num_check, num_var, 1)
    #     var_to_check_messages = zeros(Float64, num_var, num_check, 1)
    # else
        check_to_var_messages = zeros(Float64, num_check, num_var, 2)
        var_to_check_messages = zeros(Float64, num_var, num_check, 2)
    # end

    if !isempty(erasures)
        all(1 ≤ bit ≤ num_var for bit in erasures) ||
            throw(ArgumentError("Invalid bit index in erasures"))
        @inbounds for i in erasures
            chn_inits_2[i] = 0.0
        end
    end

    return H_Int, syndrome_based, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages,
        var_to_check_messages, current_bits, totals, syn
end

function _message_passing_init_Int(H::Union{Matrix{S}, T}, v::Union{Vector{S},
    Vector{AbstractFloat}, T}, max_iter::Int, kind::Symbol, Bt::Int, schedule::Symbol,
    erasures::Vector{Int}) where {S <: Integer, T <: CTMatrixTypes}

    kind ∈ (:A, :B) || throw(ArgumentError("Unknown value for parameter kind"))
    # will this work?
    if isa(H, CTMatrixTypes)
        Int(order(base_ring(H))) == 2 ||
        throw(ArgumentError("Currently only implemented for binary codes"))
        H_Int = _Flint_matrix_to_Julia_int_matrix(H)
    else
        H_Int = H
    end
    num_check, num_var = size(H_Int)
    num_check > 0 && num_var > 0 || throw(ArgumentError("Input matrix of improper dimension"))
    (kind == :B && !(1 <= Bt <= num_check)) &&
        throw(DomainError("Improper threshold for Gallager B"))
    2 <= max_iter || throw(DomainError("Number of maximum iterations must be at least two"))
    
    len_v = length(v)
    if len_v == num_var
        syndrome_based = false
    elseif len_v == num_checks
        syndrome_based = true
    else
        throw(ArgumentError("Vector has incorrect dimension"))
    end
    
    check_adj_list = [Int[] for _ in 1:num_check]
    var_adj_list = [Int[] for _ in 1:num_var]
    for r in 1:num_check
        for c in 1:num_var
            if !iszero(H_Int[r, c])
                push!(check_adj_list[r], c)
                push!(var_adj_list[c], r)
            end
        end
    end

    # # TODO: what are proper inits for syndrome and erasures here?
    # if !syndrome_based && ismissing(chn_inits)
        chn_inits = vec(_Flint_matrix_to_Julia_int_matrix(v))
    # elseif syndrome_based && ismissing(chn_inits)
    #     chn_inits = [log((1 - chn.cross_over_prob) / chn.cross_over_prob) for _ in 1:num_var]
    # end

    # could reduce this stuff to UInt8 if needed
    current_bits = zeros(Int, num_var)
    syn = zeros(Int, num_check)
    if schedule == :serial
        check_to_var_messages = zeros(Int, num_check, num_var, 1)
        var_to_check_messages = zeros(Int, num_var, num_check, 1)
    else
        check_to_var_messages = zeros(Int, num_check, num_var, 2)
        var_to_check_messages = zeros(Int, num_var, num_check, 2)
    end

    # # TODO
    # if !isempty(erasures)
    #     all(1 ≤ bit ≤ num_var for bit in erasures) ||
    #         throw(ArgumentError("Invalid bit index in erasures"))
    #     @inbounds for i in erasures
    #         chn_inits[i] = 0.0
    #     end
    # end

    return H_Int, var_adj_list, check_adj_list, chn_inits, check_to_var_messages,
        var_to_check_messages, current_bits, syn
end

function _message_passing_init_decimation(H::T, v::T, chn::Union{Missing, MPNoiseModel},
    decimated_bits_values::Vector{Tuple{Int, S}}, max_iter::Int, kind::Symbol, Bt::Int,
    chn_inits::Union{Missing, Vector{Float64}}, schedule::Symbol) where {S <: CTFieldElem,
    T <: CTMatrixTypes}

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
    
    H_Int = _Flint_matrix_to_Julia_int_matrix(H)
    w = vec(_Flint_matrix_to_Julia_int_matrix(v))
    check_adj_list = [Int[] for _ in 1:num_check]
    var_adj_list = [Int[] for _ in 1:num_var]

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

    if !isempty(decimated_bits_values)
        decimated_bits = getindex.(decimated_bits_values, 1)
        all(1 ≤ bit ≤ num_var for bit in decimated_bits) ||
            throw(ArgumentError("Invalid bit index in decimated_bits_values"))
        # TODO: not safe but use for now
        decimated_values = [iszero(pair[2]) ? 0 : 1 for pair in decimated_bits_values]
        all(iszero(val) || isone(val) for val in decimated_values) ||
            throw(ArgumentError("Invalid bit value in decimated_bits_values"))
        if kind ∈ (:SP, :MS)
            @inbounds for (i, b) in enumerate(decimated_bits)
                if iszero(decimated_values[i])
                    chn_inits[b] = 255
                else
                    chn_inits[b] = -255
                end
            end
        end
    else
        decimated_bits = Int[]
        decimated_values = Int[]
    end

    # could reduce this stuff to UInt8 if needed
    current_bits = zeros(Int, num_var)
    @inbounds for (i, v) in enumerate(decimated_bits)
        current_bits[v] = decimated_values[i]
    end

    totals = zeros(Float64, num_var)
    syn = zeros(Int, num_check)
    R = kind ∈ (:A, :B) ? Int : Float64
    if schedule == :serial
        check_to_var_messages = zeros(R, num_check, num_var, 1)
        var_to_check_messages = zeros(R, num_var, num_check, 1)
    else
        check_to_var_messages = zeros(R, num_check, num_var, 2)
        var_to_check_messages = zeros(R, num_var, num_check, 2)
    end

    if !isempty(erasures)
        all(1 ≤ bit ≤ num_var for bit in erasures) ||
            throw(ArgumentError("Invalid bit index in erasures"))
        @inbounds for i in erasures
            chn_inits[i] = 0.0
        end
    end

    return H_Int, w, var_adj_list, check_adj_list, chn_inits, decimated_bits,
        decimated_values, check_to_var_messages, var_to_check_messages, current_bits,
        totals, syn
end

#############################
      # Message Passing
#############################

function _message_passing(H::Matrix{T}, syndrome::Union{Missing, Vector{T}},
    chn_inits::Vector{Float64}, c_to_v_mess::Function, var_adj_list::Vector{Vector{Int}},
    check_adj_list::Vector{Vector{Int}}, max_iter::Int, schedule::Symbol,
    current_bits::Vector{Int}, totals::Vector{Float64}, syn::Vector{Int},
    check_to_var_messages::Array{Float64, 3}, var_to_check_messages::Array{Float64, 3},
    attenuation::Float64) where T <: Integer

    num_check, num_var = size(H)

    # first iteration for variable nodes - set to channel initialization
    @inbounds for v in 1:num_var
        @simd for c in var_adj_list[v]
            var_to_check_messages[v, c, 1] = chn_inits[v]
        end
    end

    @inbounds @simd for c in 1:num_check
        for v in check_adj_list[c]
            if !ismissing(syndrome)
                check_to_var_messages[c, v, 1] = (-1)^syndrome[c] * c_to_v_mess(c, v, 1,
                    check_adj_list, var_to_check_messages, attenuation)
            else
                check_to_var_messages[c, v, 1] = c_to_v_mess(c, v, 1, check_adj_list,
                    var_to_check_messages, attenuation)
            end
        end
    end

    # one full iteration done, check if converged
    @simd for v in 1:num_var
        totals[v] = chn_inits[v]
        for c in var_adj_list[v]
            totals[v] += check_to_var_messages[c, v, 1]
        end
        current_bits[v] = totals[v] >= 0 ? 0 : 1
    end
    LinearAlgebra.mul!(syn, H, current_bits)
    if !ismissing(syndrome)
        all(syn[i] .% 2 == syndrome[i] for i in 1:num_check) && return true, current_bits, 1 
    else
        all(iszero(syn[i] .% 2) for i in 1:num_check) && return true, current_bits, 1 
    end

    iter = 2
    # unclear if this actually handles serial properly
    schedule == :parallel ? (curr_iter = 2;) : (curr_iter = 1;)
    prev_iter = 1
    # does this propagate correctly?
    @inbounds while iter ≤ max_iter
        @simd for v in 1:num_var
            for c in var_adj_list[v]
                # this includes the channel inputs in total
                var_to_check_messages[v, c, curr_iter] = totals[v] -
                    check_to_var_messages[c, v, prev_iter]
            end
        end

        @simd for c in 1:num_check
            for v in check_adj_list[c]
                if !ismissing(syndrome)
                    check_to_var_messages[c, v, curr_iter] = (-1)^syndrome[c] * c_to_v_mess(c, v,
                    curr_iter, check_adj_list, var_to_check_messages, attenuation)
                else
                    check_to_var_messages[c, v, curr_iter] = c_to_v_mess(c, v,
                    curr_iter, check_adj_list, var_to_check_messages, attenuation)
                end
            end
        end
    
        # iteration done, check if converged
        @simd for v in 1:num_var
            totals[v] = chn_inits[v]
            for c in var_adj_list[v]
                totals[v] += check_to_var_messages[c, v, curr_iter]
            end
            current_bits[v] = totals[v] >= 0 ? 0 : 1
        end
    
        LinearAlgebra.mul!(syn, H, current_bits)
        if !ismissing(syndrome)
            all(syn[i] .% 2 == syndrome[i] for i in 1:num_check) && return true, current_bits, iter
        else
            all(iszero(syn[i] .% 2) for i in 1:num_check) && return true, current_bits, iter
        end

        if schedule == :flooding
            temp = curr_iter
            curr_iter = prev_iter
            prev_iter = temp
        end
        iter += 1
    end

    return false, current_bits, iter
end

function _message_passing_layered(H::Matrix{T}, syndrome::Union{Missing, Vector{T}},
    chn_inits::Vector{Float64}, c_to_v_mess::Function, var_adj_list::Vector{Vector{Int}},
    check_adj_list::Vector{Vector{Int}}, max_iter::Int, schedule::Symbol,
    current_bits::Vector{Int}, totals::Vector{Float64}, syn::Vector{Int},
    check_to_var_messages::Array{Float64, 3}, var_to_check_messages::Array{Float64, 3},
    attenuation::Float64, layers::Vector{Vector{Int}}) where T <: Integer

    # first iteration for variable nodes - set to channel initialization
    num_check, num_var = size(H)
    @inbounds for v in 1:num_var
        @simd for c in var_adj_list[v]
            var_to_check_messages[v, c, 2] = chn_inits[v]
        end
    end

    iter = 1
    curr_iter = 2
    prev_iter = 1
    # does this propagate correctly?
    @inbounds while iter < max_iter
        for layer in layers
            @simd for c in layer
                for v in check_adj_list[c]
                    if !ismissing(syndrome)
                        check_to_var_messages[c, v, curr_iter] = (-1)^syndrome[c] * c_to_v_mess(c, v,
                        curr_iter, check_adj_list, var_to_check_messages, attenuation)
                    else
                        check_to_var_messages[c, v, curr_iter] = c_to_v_mess(c, v,
                        curr_iter, check_adj_list, var_to_check_messages, attenuation)
                    end
                end
            end

            # TODO the only values that should be changing here are the ones connected to the check nodes in the layer
            @simd for v in 1:num_var
                totals[v] = chn_inits[v]
                for c in var_adj_list[v]
                    totals[v] += check_to_var_messages[c, v, curr_iter]
                end
                current_bits[v] = totals[v] >= 0 ? 0 : 1
            end

            @simd for v in 1:num_var
                for c in var_adj_list[v]
                    # this includes the channel inputs in total
                    # TODO: here prev was changed to curr
                    var_to_check_messages[v, c, curr_iter] = totals[v] -
                        check_to_var_messages[c, v, curr_iter]
                end
            end

            # switch these current values to previous values so the next layer can use them
            @inbounds @simd for c in layer
                for v in check_adj_list[c]
                    check_to_var_messages[c, v, prev_iter] = check_to_var_messages[c, v, curr_iter]
                end
            end

            temp = curr_iter
            curr_iter = prev_iter
            prev_iter = temp
        end

        # iteration done, check if converged
        LinearAlgebra.mul!(syn, H, current_bits)
        if !ismissing(syndrome)
            all(syn[i] % 2 == syndrome[i] for i in 1:num_check) && return true, current_bits, iter
        else
            all(iszero(syn[i] % 2) for i in 1:num_check) && return true, current_bits, iter
        end

        iter += 1
    end
    return false, current_bits, iter
end

function _MP_test()
    # TODO
    # 1. setup init function
    # 2. pass in char_adj_list, size of H, max_iter, phi
    # 3. estabilish how phi and phi^{-1} behave wrt to above functions

    check_to_var_messages::Vector{Vector{Float64}} = [zeros(Float64, length(check_adj_list[c])) for c in eachindex(check_adj_list)]
    var_to_check_messages = zeros(Float64, num_var)
    var_to_check_messages .= chn_inits

    iter = 1
    @inbounds while iter < max_iter
        # here, S is totals, current bits threshold this
        S_flag = true
        S::Float64 = 0.0
        for c in 1:num_check
            for v in check_adj_list[c]
                # ⊞
                S += phi(var_to_check_messages[v] - check_to_var_messages[c][v])
            end

            if S_flag && S < 0
                S_flag = false
            end

            for v in check_adj_list[c]
                Q_temp = var_to_check_messages[v] - check_to_var_messages[c][v]
                check_to_var_messages[c][v] = phi_inv(S - phi(Q_temp))
                var_to_check_messages[v] = Q_temp + check_to_var_messages[c][v]
            end
        end

        # how to extract the error here?
        S_flag && (return true, Int[var_to_check_messages[:] >= 0 ? 0 : 1], iter;)
        iter += 1
    end
    return false, current_bits, iter
end

# significant speedups seperating the float and int code
function _message_passing_Int(H::Matrix{T}, syndrome::Union{Missing, Vector{T}},
    chn_inits::Vector{Int}, c_to_v_mess::Function, var_adj_list::Vector{Vector{Int}},
    check_adj_list::Vector{Vector{Int}}, max_iter::Int, kind::Symbol, schedule::Symbol,
    current_bits::Vector{Int}, syn::Vector{Int}, check_to_var_messages::Array{Int, 3},
    var_to_check_messages::Array{Int, 3}, Bt::Int) where T <: Integer

    num_check, num_var = size(H)

    # first iteration for variable nodes - set to channel initialization
    @inbounds for v in 1:num_var
        @simd for c in var_adj_list[v]
            var_to_check_messages[v, c, 1] = v
        end
    end

    @simd for c in 1:num_check
        for v in check_adj_list[c]
            if !ismissing(syndrome)
                check_to_var_messages[c, v, 1] = (-1)^syndrome[c] * c_to_v_mess(c, v, 1,
                    check_adj_list, var_to_check_messages, 0.0)
            else
                check_to_var_messages[c, v, 1] = c_to_v_mess(c, v, 1, check_adj_list,
                    var_to_check_messages, 0.0)
            end
        end
    end

    # one full iteration done, check if converged
    curr_iter = 1
    @simd for v in 1:num_var
        len = length(var_adj_list[v])
        one_count = count(isone, view(check_to_var_messages, var_adj_list[v], v, curr_iter))
        d = fld(len, 2)
        current_bits[v] = one_count + (isone(chn_inits[v]) && iseven(len)) > d
    end
    LinearAlgebra.mul!(syn, H, current_bits)
    if !ismissing(syndrome)
        all(syn[i] .% 2 == syndrome[i] for i in 1:num_check) && return true, current_bits, 1 
    else
        all(iszero(syn[i] .% 2) for i in 1:num_check) && return true, current_bits, 1 
    end

    iter = 2
    schedule == :parallel ? (curr_iter = 2;) : (curr_iter = 1;)
    prev_iter = 1
    # does this propagate correctly?
    @inbounds while iter ≤ max_iter
        @simd for v in 1:num_var
            for c in var_adj_list[v]
                if kind == :A && length(var_adj_list[v]) > 1
                    if all(!Base.isequal(chn_inits[v]), check_to_var_messages[c2, v, prev_iter] for
                        c2 in var_adj_list[v] if c != c2)

                        var_to_check_messages[v, c, curr_iter] ⊻= 1
                    end
                elseif kind == :B && length(var_adj_list[v]) >= Bt
                    if count(!Base.isequal(chn_inits[v]), check_to_var_messages[c2, v, prev_iter] for 
                        c2 in var_adj_list[v] if c != c2) >= Bt

                        var_to_check_messages[v, c, curr_iter] ⊻= 1
                    end
                end
            end
        end

        @simd for c in 1:num_check
            for v in check_adj_list[c]
                if !ismissing(syndrome)
                    check_to_var_messages[c, v, 1] = (-1)^syndrome[c] * c_to_v_mess(c, v, curr_iter,
                        check_adj_list, var_to_check_messages, 0.0)
                else
                    check_to_var_messages[c, v, 1] = c_to_v_mess(c, v, curr_iter, check_adj_list,
                        var_to_check_messages, 0.0)
                end
            end
        end
    
        # iteration done, check if converged
        @simd for v in 1:num_var
            len = length(var_adj_list[v])
            one_count = count(isone, view(check_to_var_messages, var_adj_list[v], v, curr_iter))
            d = fld(len, 2)
            current_bits[v] = one_count + (isone(chn_inits[v]) && iseven(len)) > d
        end
    
        LinearAlgebra.mul!(syn, H, current_bits)
        if !ismissing(syndrome)
            all(syn[i] .% 2 == syndrome[i] for i in 1:num_check) && return true, current_bits,
                curr_iter
        else
            all(iszero(syn[i] .% 2) for i in 1:num_check) && return true, current_bits, curr_iter
        end
        
        if schedule == :flooding
            temp = curr_iter
            curr_iter = prev_iter
            prev_iter = temp
        end
        iter += 1
    end

    return false, current_bits, iter
end

function _message_passing_decimation(H::Matrix{T}, w::Vector{T}, chn_inits::Union{Missing,
    Vector{Float64}}, c_to_v_mess::Function, var_adj_list::Vector{Vector{Int}},
    check_adj_list::Vector{Vector{Int}}, max_iter::Int, kind::Symbol, schedule::Symbol,
    decimated_bits::Vector{Int}, decimated_values::Vector{Int}, current_bits::Vector{Int},
    totals::Union{Vector{Int}, Vector{Float64}}, syn::Vector{Int},
    check_to_var_messages::Union{Array{Float64, 3}, Array{Int, 3}},
    var_to_check_messages::Union{Array{Float64, 3}, Array{Int, 3}}, Bt::Int,
    attenuation::Float64, algorithm::Symbol, guided_rounds::Int) where T <: Integer

    # the inclusion of the kind statements add less than a microsecond
    num_check, num_var = size(H)
    if !isempty(decimated_bits)
        @inbounds for (i, v) in enumerate(decimated_bits)
            current_bits[v] = decimated_values[i]
        end
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
    schedule == :parallel ? (prev_iter = 2;) : (prev_iter = 1;)
    # does this propagate correctly?
    @inbounds while iter ≤ max_iter
        # variable node is already done for first iteration, so start with check nodes
        @simd for c in 1:num_check
            for v in check_adj_list[c]
                if v ∉ decimated_bits
                    check_to_var_messages[c, v, curr_iter] = c_to_v_mess(c, v, curr_iter,
                        check_adj_list, var_to_check_messages, attenuation)
                end
            end
        end

        # one full iteration done, check if converged
        if kind ∈ (:SP, :MS)
            @simd for v in 1:num_var
                if v ∉ decimated_bits
                    totals[v] = chn_inits[v]
                    for c in var_adj_list[v]
                        totals[v] += check_to_var_messages[c, v, curr_iter]
                    end
                    current_bits[v] = totals[v] >= 0 ? 0 : 1
                end
            end
        elseif kind ∈ (:A, :B)
            @simd for v in 1:num_var
                if v ∉ decimated_bits
                    len = length(var_adj_list[v])
                    one_count = count(isone, view(check_to_var_messages, var_adj_list[v], v, curr_iter))
                    d = fld(len, 2)
                    current_bits[v] = one_count + (isone(w[v]) && iseven(len)) > d
                end
            end
        end

        LinearAlgebra.mul!(syn, H, current_bits)
        iszero(syn .% 2) && return true, current_bits, iter, var_to_check_messages,
            check_to_var_messages

        if algorithm == :guided && iszero(iter % guided_rounds)
            val, index = findmax(totals)
            if val ≥ 0
                push!(decimated_bits, index)
                # push!(decimated_values, 0)
                chn_inits[v] = 255
            else
                push!(decimated_bits, index)
                # push!(decimated_values, 1)
                chn_inits[v] = -255
            end
        end

        iter += 1
        if schedule == :flooding
            temp = curr_iter
            curr_iter = prev_iter
            prev_iter = temp
        end

        if iter <= max_iter
            for v in 1:num_var
                if v ∉ decimated_bits
                    for c in var_adj_list[v]
                        if kind ∈ (:SP, :MS)
                            # this includes the channel inputs in total
                            var_to_check_messages[v, c, curr_iter] = totals[v] -
                                check_to_var_messages[c, v, prev_iter]
                            if algorithm == :auto
                                if var_to_check_messages[v, c, curr_iter] > 8
                                    push!(decimated_bits, v)
                                    # push!(decimated_values, 0)
                                    chn_inits[v] = 255
                                elseif var_to_check_messages[v, c, curr_iter] < -8
                                    push!(decimated_bits, v)
                                    # push!(decimated_values, 1)
                                    chn_inits[v] = -255
                                end
                            end
                        elseif kind == :A && length(var_adj_list[v]) > 1
                            if all(!Base.isequal(w[v]), check_to_var_messages[c2, v, prev_iter] for
                                c2 in var_adj_list[v] if c != c2)

                                var_to_check_messages[v, c, curr_iter] ⊻= 1
                                # TODO: how to auto this?
                            end
                        elseif kind == :B && length(var_adj_list[v]) >= Bt
                            if count(!Base.isequal(w[v]), check_to_var_messages[c2, v, prev_iter] for 
                                c2 in var_adj_list[v] if c != c2) >= Bt

                                var_to_check_messages[v, c, curr_iter] ⊻= 1
                                # TODO: how to auto this?
                            end
                        end
                    end
                end
            end
        end
    end

    return false, current_bits, iter, var_to_check_messages, check_to_var_messages
end

#############################
          # Methods
#############################

# Mansour, Shanbhag, "Turbo Decoder Architectures for Low-Density Parity-Check Codes" (2002)

# A layer is a collection of check-nodes such that any two check-nodes have no neighbouring variable-node in common.
# TODO latex the H
"""
    layered_schedule(H::CTMatrixTypes; schedule::Symbol = :layered, random::Bool = false)

Return a layered schedule for the parity-check matrix `H`. If `schedule` is `:parallel` or
`:serial`, layers representing these two extreme cases are returned. If `random` is `true`, the
schedule is shuffled.
"""
function layered_schedule(H::CTMatrixTypes; schedule::Symbol = :layered, random::Bool = false)
    num_check, num_var = size(H)
    num_check > 0 && num_var > 0 || throw(ArgumentError("Input matrix of improper dimension"))
    schedule ∈ (:flooding, :parallel, :serial, :layered, :semiserial) || 
        throw(ArgumentError("Unknown schedule algorithm"))
    schedule == :flooding && (schedule = :parallel;)
    schedule == :semiserial && (schedule = :layered;)

    if schedule == :layered
        check_adj_list = [Int[] for _ in 1:num_check]
        for r in 1:num_check
            for c in 1:num_var
                iszero(H[r, c]) || push!(check_adj_list[r], c)
            end
        end

        sched_list = [[1]]
        list = collect(2:num_check)
        random && shuffle!(list)
        for c in list
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
        random && shuffle!(sched_list)
    elseif schedule == :parallel
        sched_list = [collect(1:num_check)]
        random && shuffle!(sched_list[1])
    else
        # serial
        sched_list = [[i] for i in 1:num_check]
        random && shuffle!(sched_list)
    end
    return sched_list
end

# ref: Layered Decoding of Quantum LDPC Codes
function balance_of_layered_schedule(sch::Vector{Vector{Int}})
    is_empty(sch) && throw(ArgumentError("Schedule cannot be empty"))
    any(x -> is_empty(x), sch) && throw(ArgumentError("Schedule cannot contain an empty layer"))

    len = sch[1]
    all(x -> length(x) == len, sch) && return 1
    γ = 0.0
    for L_i in sch
        for L_j in sch
            temp = length(L_i) / length(L_j)
            temp > γ && (γ = temp;)
        end
    end
    return γ
end

# output should end up [1 1 1 0 0 0 0]

# H = matrix(GF(2), [1 1 0 1 1 0 0; 1 0 1 1 0 1 0; 0 1 1 1 0 0 1]);
# v = matrix(GF(2), 7, 1, [1, 1, 0, 0, 0, 0, 0]);
# syn = H * v;
# nm = MPNoiseModel(:BSC, 1/7);
# decimated_bits_values = [(1, base_ring(v)(1))];
# flag, out, iter = sum_product(H, v, nm); flag
# flag, out, iter = sum_product_syndrome(H, syn, nm); flag
# flag, out, iter = CodingTheory.sum_product_decimation(H, v, nm, decimated_bits_values); flag

# TODO: all of the min-sum versions fail
# flag, out, iter = min_sum(H, v, nm)
# flag, out, iter = min_sum_syndrome(H, syn, nm); flag
# flag, out, iter = CodingTheory.min_sum_decimation(H, v, nm, decimated_bits_values); flag

# flag, out, iter = Gallager_A(H, v); flag
# flag, out, iter = Gallager_B(H, v); flag
# TODO: B fails
