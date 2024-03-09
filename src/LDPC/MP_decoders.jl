# Copyright (c) 2023 - 2024 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
         # Gallager
#############################

function _Gallager_A_check_node_message(c::Int, v::Int, iter::Int, check_adj_list,
    var_to_check_messages, attenuation::Float64)

    @inbounds reduce(⊻, var_to_check_messages[v2, c, iter] for v2 in check_adj_list[c] if v2 != v)
end
_Gallager_B_check_node_message(c::Int, v::Int, iter::Int, check_adj_list, var_to_check_messages,
    attenuation::Float64) = _Gallager_A_check_node_message(c, v, iter, check_adj_list,
    var_to_check_messages, attenuation)

# TODO: does it make sense to do a Gallager syndrome setup?
# TODO does it make sense to take channel inits here?
"""
    Gallager_A(H::T, v::T; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where T <: CTMatrixTypes

Run the Gallager-A decoder with the parity-check matrix `H` and received vector `v`.
"""
function Gallager_A(H::T, v::T; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} =
    missing, schedule::Symbol = :flooding) where T <: CTMatrixTypes

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(v) ≠ (nc, 1) && size(v) ≠ (1, nc)) && throw(ArgumentError("Vector has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial) || throw(ArgumentError("Unknown schedule type"))
    schedule == :parallel && (schedule = :flooding;)

    # initialization - do these outside to reduce allocations when looped
    H_Int, w, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages,
        var_to_check_messages, current_bits, totals, syn = _message_passing_init(H, v, missing,
        max_iter, :A, 2, chn_inits, schedule)

    # check inits here
    return _message_passing(H_Int, w, chn_inits_2, _Gallager_A_check_node_message,
        var_adj_list, check_adj_list, max_iter, :A, schedule, current_bits, totals, syn,
        check_to_var_messages, var_to_check_messages, 0, 0.0)
end

# TODO: threshold in docstring
"""
    Gallager_B(H::T, v::T; max_iter::Int = 100, threshold::Int = 2, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where T <: CTMatrixTypes

Run the Gallager-B decoder with the parity-check matrix `H` and received vector `v`.
"""
function Gallager_B(H::T, v::T; max_iter::Int = 100, threshold::Int = 2, chn_inits::Union{Missing,
    Vector{Float64}} = missing, schedule::Symbol = :flooding) where T <: CTMatrixTypes
    
    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(v) ≠ (nc, 1) && size(v) ≠ (1, nc)) && throw(ArgumentError("Vector has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial) || throw(ArgumentError("Unknown schedule type"))
    schedule == :parallel && (schedule = :flooding;)

    # initialization - do these outside to reduce allocations when looped
    H_Int, w, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages,
        var_to_check_messages, current_bits, totals, syn = _message_passing_init(H, v, missing,
        max_iter, :B, threshold, chn_inits, schedule)

    # check inits here
    return _message_passing(H_Int, w, chn_inits_2, _Gallager_B_check_node_message,
        var_adj_list, check_adj_list, max_iter, :B, schedule, current_bits, totals, syn,
        check_to_var_messages, var_to_check_messages, threshold, 0.0)
end

#############################
        # Sum product
#############################

function _SP_check_node_message(c::Int, v::Int, iter::Int, check_adj_list, var_to_check_messages,
    attenuation::Float64)

    ϕ(x) = -log(tanh(0.5 * x))
    temp = 0.0
    s = 1
    @inbounds for v2 in check_adj_list[c]
        if v2 != v
            x = var_to_check_messages[v2, c, iter]
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
function _SP_check_node_message_box_plus(c::Int, v::Int, iter::Int, check_adj_list,
    var_to_check_messages, attenuation::Float64)

    @inbounds ⊞(var_to_check_messages[v2, c, iter] for v2 in check_adj_list[c] if v2 != v)
end

"""
    sum_product(H::T, v::T, chn::MPNoiseModel; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where T <: CTMatrixTypes

Run the sum-product algorithm with the parity-check matrix `H`, received vector `v`, and channel
`chn`.

# Notes
* Use `chn_inits` to pass in soft information.
* The options for `schedule` are `:flooding` (`:parallel`) or `:serial`.
"""
function sum_product(H::T, v::T, chn::MPNoiseModel; max_iter::Int = 100, chn_inits::Union{Missing,
    Vector{Float64}} = missing, schedule::Symbol = :flooding) where T <: CTMatrixTypes

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(v) ≠ (nc, 1) && size(v) ≠ (1, nc)) && throw(ArgumentError("Vector has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial) || throw(ArgumentError("Unknown schedule type"))
    schedule == :parallel && (schedule = :flooding;)

    # initialization - do these outside to reduce allocations when looped
    H_Int, w, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages,
        var_to_check_messages, current_bits, totals, syn = _message_passing_init(H, v, chn,
        max_iter, :SP, 2, chn_inits, schedule)

    return _message_passing(H_Int, w, chn_inits_2, _SP_check_node_message, var_adj_list,
        check_adj_list, max_iter, :SP, schedule, current_bits, totals, syn, check_to_var_messages,
        var_to_check_messages, 0, 0.0)
end

"""
    sum_product_box_plus(H::T, v::T, chn::MPNoiseModel; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where T <: CTMatrixTypes

Run the sum-product box-plus algorithm with the parity-check matrix `H`, received vector `v`, and
channel `chn`.

# Notes
* Use `chn_inits` to pass in soft information.
* The options for `schedule` are `:flooding` (`:parallel`) or `:serial`.
"""
function sum_product_box_plus(H::T, v::T, chn::MPNoiseModel; max_iter::Int = 100,
    chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where T <:
    CTMatrixTypes

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(v) ≠ (nc, 1) && size(v) ≠ (1, nc)) && throw(ArgumentError("Vector has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial) || throw(ArgumentError("Unknown schedule type"))
    schedule == :parallel && (schedule = :flooding;)

    # initialization - do these outside to reduce allocations when looped
    H_Int, w, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages,
        var_to_check_messages, current_bits, totals, syn = _message_passing_init(H, v, chn,
        max_iter, :SP, 2, chn_inits, schedule)

    return _message_passing(H_Int, w, chn_inits_2, _SP_check_node_message_box_plus, var_adj_list,
        check_adj_list, max_iter, :SP, schedule, current_bits, totals, syn, check_to_var_messages,
        var_to_check_messages, 0, 0.0)
end

"""
    sum_product_syndrome(H::T, syndrome::T, chn::MPNoiseModel; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where T <: CTMatrixTypes

Run the syndrome-based sum-product algorithm with the parity-check matrix `H`, syndrome `syndrome`,
and channel `chn`.

# Notes
* Use `chn_inits` to pass in soft information.
* The options for `schedule` are `:flooding` (`:parallel`) or `:serial`.
"""
function sum_product_syndrome(H::T, syndrome::T, chn::MPNoiseModel; max_iter::Int = 100,
    chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where T <:
    CTMatrixTypes

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(syndrome) ≠ (nr, 1) && size(syndrome) ≠ (1, nr)) && throw(ArgumentError("Syndrome has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial) || throw(ArgumentError("Unknown schedule type"))
    schedule == :parallel && (schedule = :flooding;)

    # initialization - do these outside to reduce allocations when looped
    H_Int, syn_Int, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages,
        var_to_check_messages, current_bits, totals, syn = _message_passing_init_syndrome(H,
        syndrome, chn, max_iter, chn_inits, :SP, schedule)

    return _message_passing_syndrome(H_Int, syn_Int, chn_inits_2, _SP_check_node_message,
        var_adj_list, check_adj_list, max_iter, :SP, schedule, current_bits, totals, syn, 
        check_to_var_messages, var_to_check_messages, 0, 0.0)
end

# believe this can be merged into an optional argument of the above but keep for now so as not to break working code
function sum_product_decimation(H::T, v::T, chn::MPNoiseModel, decimated_bits_values::Vector{Tuple{
    Int, S}}; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} = missing, 
    schedule::Symbol = :flooding) where {T <: CTMatrixTypes, S <: CTFieldElem}

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(v) ≠ (nc, 1) && size(v) ≠ (1, nc)) && throw(ArgumentError("Vector has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial) || throw(ArgumentError("Unknown schedule type"))
    schedule == :parallel && (schedule = :flooding;)

    # initialization - do these outside to reduce allocations when looped
    H_Int, w, var_adj_list, check_adj_list, chn_inits_2, decimated_bits, decimated_values,
        check_to_var_messages, var_to_check_messages, current_bits, totals, syn =
        _message_passing_init_decimation(H, v, chn, decimated_bits_values, max_iter, :SP, 2,
        chn_inits, schedule)

    return _message_passing_decimation(H_Int, w, chn_inits_2, _SP_check_node_message,
        var_adj_list, check_adj_list, max_iter, :SP, schedule, decimated_bits, decimated_values,
        current_bits, totals, syn, check_to_var_messages, var_to_check_messages, 0, 0.0)
end

# decimation, guided decimation, automatic decimation
# ensemble decoding
# syndrome-based decimation

#############################
        # Min Sum
#############################

function _MS_check_node_message(c::Int, v::Int, iter::Int, check_adj_list, var_to_check_messages,
    attenuation::Float64)

    @inbounds temp = var_to_check_messages[check_adj_list[c][1], c, iter]
    s = 1
    @inbounds for v2 in check_adj_list[c]
        if v2 != v
            x = var_to_check_messages[v2, c, iter]
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
    min_sum(H::T, v::T, chn::MPNoiseModel; max_iter::Int = 100, attenuation::Float64 = 0.5, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where T <: CTMatrixTypes

Run the min-sum algorithm with the parity-check matrix `H`, received vector `v`, and channel
`chn`.

# Notes
* Use `chn_inits` to pass in soft information.
* The options for `schedule` are `:flooding` (`:parallel`) or `:serial`.
* Set the normalization constant with `attenuation`.
"""
function min_sum(H::T, v::T, chn::MPNoiseModel; max_iter::Int = 100, attenuation::Float64 =
    0.5, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where
    T <: CTMatrixTypes

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(v) ≠ (nc, 1) && size(v) ≠ (1, nc)) && throw(ArgumentError("Vector has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial) || throw(ArgumentError("Unknown schedule type"))
    schedule == :parallel && (schedule = :flooding;)

    # initialization - do these outside to reduce allocations when looped
    # TODO: check this init call
    H_Int, w, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages,
        var_to_check_messages, current_bits, totals, syn = _message_passing_init(H, v, chn,
        max_iter, :MS, 2, chn_inits, schedule)

    return _message_passing(H_Int, w, chn_inits_2, _MS_check_node_message, var_adj_list,
        check_adj_list, max_iter, :MS, schedule, current_bits, totals, syn, check_to_var_messages,
        var_to_check_messages, 0, attenuation)
end

"""
    min_sum_syndrome(H::T, syndrome::T, chn::MPNoiseModel; max_iter::Int = 100, attenuation::Float64 = 0.5, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol = :flooding) where T <: CTMatrixTypes

Run the syndrome-based min-sum algorithm with the parity-check matrix `H`, syndrome `syndrome`,
and channel `chn`.

# Notes
* Use `chn_inits` to pass in soft information.
* The options for `schedule` are `:flooding` (`:parallel`) or `:serial`.
* Set the normalization constant with `attenuation`.
"""
function min_sum_syndrome(H::T, syndrome::T, chn::MPNoiseModel; max_iter::Int = 100,
    attenuation::Float64 = 0.5, chn_inits::Union{Missing, Vector{Float64}} = missing,
    schedule::Symbol = :flooding) where T <: CTMatrixTypes

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(syndrome) ≠ (nr, 1) && size(syndrome) ≠ (1, nr)) && throw(ArgumentError("Syndrome has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial) || throw(ArgumentError("Unknown schedule type"))
    schedule == :parallel && (schedule = :flooding;)

    # initialization - do these outside to reduce allocations when looped
    H_Int, syn_Int, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages,
        var_to_check_messages, current_bits, totals, syn = _message_passing_init_syndrome(H,
        syndrome, chn, max_iter, chn_inits, :MS, schedule)

    return _message_passing_syndrome(H_Int, syn_Int, chn_inits_2, _MS_check_node_message,
        var_adj_list, check_adj_list, max_iter, :MS, schedule, current_bits, totals, syn, 
        check_to_var_messages, var_to_check_messages, 0, attenuation)
end

# believe this can be merged into an optional argument of the above but keep for now so as not to break working code
function min_sum_decimation(H::T, v::T, chn::MPNoiseModel, decimated_bits_values::Vector{Tuple{
    Int, S}}; max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} = missing, 
    schedule::Symbol = :flooding) where {T <: CTMatrixTypes, S <: CTFieldElem}

    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    (size(v) ≠ (nc, 1) && size(v) ≠ (1, nc)) && throw(ArgumentError("Vector has incorrect dimension"))
    # do we want to flip it if necessary?
    max_iter > 0 || throw(DomainError("max_iter must be a positive integer"))
    schedule ∈ (:flooding, :parallel, :serial) || throw(ArgumentError("Unknown schedule type"))
    schedule == :parallel && (schedule = :flooding;)

    # initialization - do these outside to reduce allocations when looped
    H_Int, w, var_adj_list, check_adj_list, chn_inits_2, decimated_bits, decimated_values,
        check_to_var_messages, var_to_check_messages, current_bits, totals, syn =
        _message_passing_init_decimation(H, v, chn, decimated_bits_values, max_iter, :MS, 2,
        chn_inits, schedule)

    return _message_passing_decimation(H_Int, w, chn_inits_2, _MS_check_node_message,
        var_adj_list, check_adj_list, max_iter, :MS, schedule, decimated_bits, decimated_values,
        current_bits, totals, syn, check_to_var_messages, var_to_check_messages, 0, 0.0)
end

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

function _message_passing_init(H::T, v::T, chn::Union{Missing, MPNoiseModel}, max_iter::Int,
    kind::Symbol, Bt::Int, chn_inits::Union{Missing, Vector{Float64}}, schedule::Symbol) where T <:
    CTMatrixTypes

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

    # could reduce this stuff to UInt8 if needed
    current_bits = zeros(Int, num_var)
    totals = zeros(Float64, num_var)
    syn = zeros(Int, num_check)
    R = kind ∈ (:A, :B) ? Int : Float64
    if schedule == :flooding
        check_to_var_messages = zeros(R, num_check, num_var, 2)
        var_to_check_messages = zeros(R, num_var, num_check, 2)
    else
        check_to_var_messages = zeros(R, num_check, num_var, 1)
        var_to_check_messages = zeros(R, num_var, num_check, 1)
    end

    return H_Int, w, var_adj_list, check_adj_list, chn_inits, check_to_var_messages,
        var_to_check_messages, current_bits, totals, syn
end

function _message_passing_init_syndrome(H::T, syndrome::T, chn::Union{Missing, MPNoiseModel}, 
    max_iter::Int, chn_inits::Union{Missing, Vector{Float64}}, kind::Symbol, schedule::Symbol
    ) where T <: CTMatrixTypes

    Int(order(base_ring(H))) == 2 ||
        throw(ArgumentError("Currently only implemented for binary codes"))
    num_check, num_var = size(H)
    num_check > 0 && num_var > 0 || throw(ArgumentError("Input matrix of improper dimension"))
    length(syndrome) == num_check || throw(ArgumentError("Syndrome has incorrect dimension"))
    2 <= max_iter || throw(DomainError("Number of maximum iterations must be at least two"))
    
    H_Int = _Flint_matrix_to_Julia_int_matrix(H)
    syndrome_Int = vec(_Flint_matrix_to_Julia_int_matrix(syndrome))
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

   # could reduce this stuff to UInt8 if needed
   current_bits = zeros(Int, num_var)
   totals = zeros(Float64, num_var)
   syn = zeros(Int, num_check)
   R = kind ∈ (:A, :B) ? Int : Float64
   if schedule == :flooding
       check_to_var_messages = zeros(R, num_check, num_var, 2)
       var_to_check_messages = zeros(R, num_var, num_check, 2)
   else
       check_to_var_messages = zeros(R, num_check, num_var, 1)
       var_to_check_messages = zeros(R, num_var, num_check, 1)
   end

   return H_Int, syndrome_Int, var_adj_list, check_adj_list, chn_inits, check_to_var_messages,
       var_to_check_messages, current_bits, totals, syn
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

    # could reduce this stuff to UInt8 if needed
    current_bits = zeros(Int, num_var)
    @inbounds for (i, v) in enumerate(decimated_bits)
        current_bits[v] = decimated_values[i]
    end

    totals = zeros(Float64, num_var)
    syn = zeros(Int, num_check)
    R = kind ∈ (:A, :B) ? Int : Float64
    if schedule == :flooding
        check_to_var_messages = zeros(R, num_check, num_var, 2)
        var_to_check_messages = zeros(R, num_var, num_check, 2)
    else
        check_to_var_messages = zeros(R, num_check, num_var, 1)
        var_to_check_messages = zeros(R, num_var, num_check, 1)
    end

    return H_Int, w, var_adj_list, check_adj_list, chn_inits, decimated_bits,
        decimated_values, check_to_var_messages, var_to_check_messages, current_bits,
        totals, syn
end

#############################
      # Message Passing
#############################

function _message_passing(H::Matrix{T}, w::Vector{T}, chn_inits::Union{Missing,
    Vector{Float64}}, c_to_v_mess::Function, var_adj_list::Vector{Vector{Any}},
    check_adj_list::Vector{Vector{Any}}, max_iter::Int, kind::Symbol, schedule::Symbol,
    current_bits::Vector{Int}, totals::Union{Vector{Int}, Vector{Float64}}, syn::Vector{Int},
    check_to_var_messages::Union{Array{Float64, 3}, Array{Int, 3}},
    var_to_check_messages::Union{Array{Float64, 3}, Array{Int, 3}}, Bt::Int,
    attenuation::Float64) where T <: Integer

    # the inclusion of the kind statements add less than a microsecond
    num_check, num_var = size(H)
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
    # @inbounds
    while iter ≤ max_iter
        # variable node is already done for first iteration, so start with check nodes
        @simd for c in 1:num_check
            for v in check_adj_list[c]
                check_to_var_messages[c, v, curr_iter] = c_to_v_mess(c, v, curr_iter, 
                    check_adj_list, var_to_check_messages, attenuation)
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
 
function _message_passing_syndrome(H::Matrix{T}, syndrome::Vector{T}, chn_inits::Union{
    Missing, Vector{Float64}}, c_to_v_mess::Function, var_adj_list::Vector{Vector{Any}},
    check_adj_list::Vector{Vector{Any}}, max_iter::Int, kind::Symbol, schedule::Symbol,
    current_bits::Vector{Int}, totals::Union{Vector{Int}, Vector{Float64}}, syn::Vector{Int},
    check_to_var_messages::Union{Array{Float64, 3}, Array{Int, 3}},
    var_to_check_messages::Union{Array{Float64, 3}, Array{Int, 3}}, Bt::Int,
    attenuation::Float64) where T <: Integer

    # the inclusion of the kind statements add less than a microsecond
    num_check, num_var = size(H)
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

function _message_passing_decimation(H::Matrix{T}, w::Vector{T}, chn_inits::Union{Missing,
    Vector{Float64}}, c_to_v_mess::Function, var_adj_list::Vector{Vector{Any}},
    check_adj_list::Vector{Vector{Any}}, max_iter::Int, kind::Symbol, schedule::Symbol,
    decimated_bits::Vector{Int}, decimated_values::Vector{Int}, current_bits::Vector{Int},
    totals::Union{Vector{Int}, Vector{Float64}}, syn::Vector{Int},
    check_to_var_messages::Union{Array{Float64, 3}, Array{Int, 3}},
    var_to_check_messages::Union{Array{Float64, 3}, Array{Int, 3}}, Bt::Int,
    attenuation::Float64) where T <: Integer

    # the inclusion of the kind statements add less than a microsecond
    num_check, num_var = size(H)
    @inbounds for (i, v) in enumerate(decimated_bits)
        current_bits[v] = decimated_values[i]
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
    end

    return false, current_bits, iter, var_to_check_messages, check_to_var_messages
end

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

# output should end up [1 1 1 0 0 0 0]

# H = matrix(GF(2), [1 1 0 1 1 0 0; 1 0 1 1 0 1 0; 0 1 1 1 0 0 1]);
# v = matrix(GF(2), 7, 1, [1, 1, 0, 0, 0, 0, 0]);
# syn = H * v;
# nm = MPNoiseModel(:BSC, 1/7);
# decimated_bits_values = [(1, base_ring(v)(1))];
# flag, out, iter, vtoc, ctov = sum_product(H, v, nm); flag
# flag, out, iter, vtoc, ctov = sum_product_syndrome(H, syn, nm); flag
# flag, out, iter, vtoc, ctov = CodingTheory.sum_product_decimation(H, v, nm, decimated_bits_values); flag

# TODO: all of the min-sum versions fail
# flag, out, iter, vtoc, ctov = min_sum(H, v, nm); flag
# flag, out, iter, vtoc, ctov = min_sum_syndrome(H, syn, nm); flag
# flag, out, iter, vtoc, ctov = CodingTheory.min_sum_decimation(H, v, nm, decimated_bits_values); flag

# flag, out, iter, vtoc, ctov = Gallager_A(H, v); flag
# flag, out, iter, vtoc, ctov = Gallager_B(H, v); flag
# TODO: B fails
