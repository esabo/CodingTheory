# Copyright (c) 2023 - 2024 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
# Noise Model
#############################

struct MPNoiseModel
    type::Symbol
    cross_over_prob::Union{Float64,Missing}
    sigma::Union{Float64,Missing}
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

#############################
# Importance Sampling
#############################

subset_probability(n::Int, wt::Int, p::Float64) = Float64(
    binomial(BigInt(n), BigInt(wt)) * (BigFloat(p)^wt) * ((1 - BigFloat(p))^(n - wt)),
)

function find_subsets(n::Int, p::Float64, tol::Float64)
    subsets = Vector{Tuple{Int,Float64}}()
    for i = 0:n
        prob = subset_probability(n, i, p)
        if prob > tol
            push!(subsets, (i, prob))
        else
            return subsets
        end
    end
    return subsets
end

function min_max_subsets(n::Int, p_arr::Vector{Float64}, tol::Float64)
    min_subset = 100.0
    max_subset = 0.0
    for p in p_arr
        subsets = find_subsets(n, p, tol)
        if isempty(subsets)
            return -1, -1
        else
            if subsets[1][1] < min_subset
                min_subset = subsets[1][1]
            end
            if subsets[end][1] > max_subset
                max_subset = subsets[end][1]
            end
        end
        # println("p = ", p)
        # println("subsets: ", subsets)
    end
    return min_subset, max_subset
end

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

# function decoder_simulation(H::CTMatrixTypes, decoder::Symbol, noise_type::Symbol,
#     noise::Union{Vector{<:Real}, AbstractRange{<:Real}}, max_iter::Int = 100, num_runs::Union{Int,
#     Vector{Int}} = [100000 for n in noise], seed::Union{Int, Nothing} = nothing)

# decoder ∈ (:A, :B, :SP, :MS) || throw(ArgumentError("Unsupported decoder"))
# noise_type ∈ (:BSC, :BAWGNC) || throw(ArgumentError("Only supports BSC and BAWGNC"))
# decoder ∈ (:A, :B) && noise_type == :BAWGNC && throw(ArgumentError("BAWGNC not supported for Gallager decoders."))
# 0 <= minimum(noise) || throw(ArgumentError("Must have non-negative noise"))
# maximum(noise) > 1 && noise_type == :BSC && throw(ArgumentError("Crossover probability must be in the range [0, 1]"))

function decoders_test(H::CTMatrixTypes; verbose::Bool = true)
    # initial parameters
    # noise is assumed to be sorted
    # noise = [5e-4, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05]
    left = 1e-4
    right = 0.04
    intervals = 15
    Δ = (right - left) / (intervals - 1)
    noise = collect(left:Δ:right)
    len_noise = length(noise)
    noise_model = :BSC
    num_runs = 15000
    nr, n = size(H)
    # it suffices to decode the zero codeword
    v = zero_matrix(base_ring(H), 1, n)
    max_iter = 32
    dist = Uniform(0.0, 1.0)
    algorithm = :manual
    guided_rounds = 10
    attenuation = 0.625
    tolerance = 1e-9

    num_threads = Threads.nthreads()
    runs_per_thread = cld(num_runs, num_threads)
    new_num_runs = runs_per_thread * num_threads
    verbose && println(
        "Number of threads: $num_threads, runs per thread: $runs_per_thread, new number of runs: $new_num_runs",
    )

    # flooding
    # SP
    FER_SP = zeros(Float64, len_noise)
    # BER_SP = zeros(Float64, len_noise)
    # syndrome-based SP
    FER_SP_syn = zeros(Float64, len_noise)
    # BER_SP_syn = zeros(Float64, len_noise)
    # SP with decimation
    FER_SP_dec = zeros(Float64, len_noise)
    # BER_SP_dec = zeros(Float64, len_noise)
    # MS
    FER_MS = zeros(Float64, len_noise)
    # BER_MS = zeros(Float64, len_noise)
    # syndrome-based MS
    FER_MS_syn = zeros(Float64, len_noise)
    # BER_MS_syn = zeros(Float64, len_noise)
    # MS with decimation
    FER_MS_dec = zeros(Float64, len_noise)
    # BER_MS_dec = zeros(Float64, len_noise)
    # MS with correction
    FER_MS_C = zeros(Float64, len_noise)
    # BER_MS_C = zeros(Float64, len_noise)
    # syndrome-based MS with correction
    FER_MS_C_syn = zeros(Float64, len_noise)
    # BER_MS_C_syn = zeros(Float64, len_noise)
    # MS with correction and decimation
    FER_MS_C_dec = zeros(Float64, len_noise)
    # BER_MS_C_dec = zeros(Float64, len_noise)

    # serial
    # SP
    FER_SP_s = zeros(Float64, len_noise)
    # BER_SP = zeros(Float64, len_noise)
    # syndrome-based SP
    FER_SP_syn_s = zeros(Float64, len_noise)
    # BER_SP_syn = zeros(Float64, len_noise)
    # SP with decimation
    FER_SP_dec_s = zeros(Float64, len_noise)
    # BER_SP_dec = zeros(Float64, len_noise)
    # MS
    FER_MS_s = zeros(Float64, len_noise)
    # BER_MS = zeros(Float64, len_noise)
    # syndrome-based MS
    FER_MS_syn_s = zeros(Float64, len_noise)
    # BER_MS_syn = zeros(Float64, len_noise)
    # MS with decimation
    FER_MS_dec_s = zeros(Float64, len_noise)
    # BER_MS_dec = zeros(Float64, len_noise)
    # MS with correction
    FER_MS_C_s = zeros(Float64, len_noise)
    # BER_MS_C = zeros(Float64, len_noise)
    # syndrome-based MS with correction
    FER_MS_C_syn_s = zeros(Float64, len_noise)
    # BER_MS_C_syn = zeros(Float64, len_noise)
    # MS with correction and decimation
    FER_MS_C_dec_s = zeros(Float64, len_noise)
    # BER_MS_C_dec = zeros(Float64, len_noise)

    # flooding
    local_counts_SP = zeros(Int, num_threads)
    local_counts_SP_syn = zeros(Int, num_threads)
    local_counts_SP_dec = zeros(Int, num_threads)
    local_counts_MS = zeros(Int, num_threads)
    local_counts_MS_syn = zeros(Int, num_threads)
    local_counts_MS_dec = zeros(Int, num_threads)
    local_counts_MS_C = zeros(Int, num_threads)
    local_counts_MS_C_syn = zeros(Int, num_threads)
    local_counts_MS_C_dec = zeros(Int, num_threads)

    # serial
    local_counts_SP_s = zeros(Int, num_threads)
    local_counts_SP_syn_s = zeros(Int, num_threads)
    local_counts_SP_dec_s = zeros(Int, num_threads)
    local_counts_MS_s = zeros(Int, num_threads)
    local_counts_MS_syn_s = zeros(Int, num_threads)
    local_counts_MS_dec_s = zeros(Int, num_threads)
    local_counts_MS_C_s = zeros(Int, num_threads)
    local_counts_MS_C_syn_s = zeros(Int, num_threads)
    local_counts_MS_C_dec_s = zeros(Int, num_threads)

    # check which should be importance sampled
    imp_sam_ind = len_noise
    sort!(noise)
    for i = 1:len_noise
        # first probability where the expected value of the error weight under
        # direct sampling is nontrivial in a way that isn't going to severely
        # under estimate the FER
        if n * noise[i] > 3
            imp_sam_ind = i - 1
            break
        end
    end

    # imp_sam_ind = 0
    if imp_sam_ind ≠ 0
        noise_impt = noise[1:imp_sam_ind]
        min_subset, max_subset = min_max_subsets(n, noise, tolerance)
        subsets = [i for i = min_subset:max_subset]

        # flooding
        subset_dic_SP = Dict{Int,Float64}()
        subset_dic_SP_syn = Dict{Int,Float64}()
        subset_dic_SP_dec = Dict{Int,Float64}()
        subset_dic_MS = Dict{Int,Float64}()
        subset_dic_MS_syn = Dict{Int,Float64}()
        subset_dic_MS_dec = Dict{Int,Float64}()
        subset_dic_MS_C = Dict{Int,Float64}()
        subset_dic_MS_C_syn = Dict{Int,Float64}()
        subset_dic_MS_C_dec = Dict{Int,Float64}()

        # serial
        subset_dic_SP_s = Dict{Int,Float64}()
        subset_dic_SP_syn_s = Dict{Int,Float64}()
        subset_dic_SP_dec_s = Dict{Int,Float64}()
        subset_dic_MS_s = Dict{Int,Float64}()
        subset_dic_MS_syn_s = Dict{Int,Float64}()
        subset_dic_MS_dec_s = Dict{Int,Float64}()
        subset_dic_MS_C_s = Dict{Int,Float64}()
        subset_dic_MS_C_syn_s = Dict{Int,Float64}()
        subset_dic_MS_C_dec_s = Dict{Int,Float64}()

        first_fail_flag = false
        second_fail_flag = false

        if verbose
            println("Starting importance sampling on $(length(noise_impt)) noise values.")
            println(
                "Minimum subset: $min_subset, maximum subset: $max_subset, tolerance: $tolerance",
            )
        end

        p = noise_impt[end]
        verbose && println("Starting importance sampling in noise range $noise_impt")
        # # initialize everything for p
        chn = MPNoiseModel(noise_model, p)
        init_0 = log((1 - p) / p)
        init_1 = log(p / (1 - p))

        sample_range = collect(1:n)
        for s in subsets
            if s == 0
                # flooding
                subset_dic_SP[0] = 0.0
                subset_dic_SP_syn[0] = 0.0
                subset_dic_SP_dec[0] = 0.0
                subset_dic_MS[0] = 0.0
                subset_dic_MS_syn[0] = 0.0
                subset_dic_MS_dec[0] = 0.0
                subset_dic_MS_C[0] = 0.0
                subset_dic_MS_C_syn[0] = 0.0
                subset_dic_MS_C_dec[0] = 0.0

                # serial
                subset_dic_SP_s[0] = 0.0
                subset_dic_SP_syn_s[0] = 0.0
                subset_dic_SP_dec_s[0] = 0.0
                subset_dic_MS_s[0] = 0.0
                subset_dic_MS_syn_s[0] = 0.0
                subset_dic_MS_dec_s[0] = 0.0
                subset_dic_MS_C_s[0] = 0.0
                subset_dic_MS_C_syn_s[0] = 0.0
                subset_dic_MS_C_dec_s[0] = 0.0
            elseif second_fail_flag
                # flooding
                subset_dic_SP[s] = 1.0
                subset_dic_SP_syn[s] = 1.0
                subset_dic_SP_dec[s] = 1.0
                subset_dic_MS[s] = 1.0
                subset_dic_MS_syn[s] = 1.0
                subset_dic_MS_dec[s] = 1.0
                subset_dic_MS_C[s] = 1.0
                subset_dic_MS_C_syn[s] = 1.0
                subset_dic_MS_C_dec[s] = 1.0

                # serial
                subset_dic_SP_s[s] = 1.0
                subset_dic_SP_syn_s[s] = 1.0
                subset_dic_SP_dec_s[s] = 1.0
                subset_dic_MS_s[s] = 1.0
                subset_dic_MS_syn_s[s] = 1.0
                subset_dic_MS_dec_s[s] = 1.0
                subset_dic_MS_C_s[s] = 1.0
                subset_dic_MS_C_syn_s[s] = 1.0
                subset_dic_MS_C_dec_s[s] = 1.0
            else
                verbose && println("Subset $s")
                Threads.@threads for th = 1:num_threads
                    # for th in 1:1
                    err = zeros(Int, n)
                    err_locs = zeros(Int, s)
                    syn_Int = zeros(Int, nr)
                    erasures = Int[]
                    chn_inits = zeros(Float64, n)
                    # SP and MS for BSC
                    H_Int,
                    _,
                    var_adj_list,
                    check_adj_list,
                    chn_inits,
                    check_to_var_messages_f,
                    var_to_check_messages_f,
                    current_bits,
                    totals,
                    syn = _message_passing_init(
                        H,
                        v,
                        chn,
                        max_iter,
                        :SP,
                        chn_inits,
                        :flooding,
                        erasures,
                    )
                    _,
                    _,
                    _,
                    _,
                    _,
                    check_to_var_messages_s,
                    var_to_check_messages_s,
                    _,
                    _,
                    _ = _message_passing_init(
                        H,
                        v,
                        chn,
                        max_iter,
                        :SP,
                        chn_inits,
                        :serial,
                        erasures,
                    )
                    # for syndrome-based we need to change the inits
                    chn_inits_syn = [log((1 - p) / p) for _ = 1:n]

                    for _ = 1:runs_per_thread
                        # sample
                        err[:] .= 0
                        chn_inits[:] .= init_0
                        sample!(sample_range, err_locs, replace = false)
                        @inbounds for e in err_locs
                            err[e] = 1
                            chn_inits[e] = init_1
                        end
                        LinearAlgebra.mul!(syn_Int, H_Int, err)
                        @inbounds @simd for i = 1:nr
                            syn_Int[i] %= 2
                        end
                        # pick a random bit to decimate
                        # bit = rand(1:n)
                        # decimated_bits = [bit]
                        # # for now, setup as a genie-assisted decoder
                        # decimated_values = [err[bit]]

                        # flooding
                        # run SP
                        flag, out, _ = _message_passing(
                            H_Int,
                            missing,
                            chn_inits,
                            _SP_check_node_message_box_plus,
                            var_adj_list,
                            check_adj_list,
                            max_iter,
                            :flooding,
                            current_bits,
                            totals,
                            syn,
                            check_to_var_messages_f,
                            var_to_check_messages_f,
                            0.0,
                        )
                        (!flag || !iszero(out)) && (local_counts_SP[th] += 1;)

                        # reset inputs for next run, but don't re-allocate new memory
                        check_to_var_messages_f[:, :, :] .= 0.0
                        var_to_check_messages_f[:, :, :] .= 0.0

                        # run syndrome-based SP
                        flag, out, _, = _message_passing(
                            H_Int,
                            syn_Int,
                            chn_inits_syn,
                            _SP_check_node_message_box_plus,
                            var_adj_list,
                            check_adj_list,
                            max_iter,
                            :flooding,
                            current_bits,
                            totals,
                            syn,
                            check_to_var_messages_f,
                            var_to_check_messages_f,
                            0.0,
                        )
                        (!flag || out ≠ err) && (local_counts_SP_syn[th] += 1;)

                        # reset inputs for next run, but don't re-allocate new memory
                        check_to_var_messages_f[:, :, :] .= 0.0
                        var_to_check_messages_f[:, :, :] .= 0.0

                        # run MS
                        flag, out, _ = _message_passing(
                            H_Int,
                            missing,
                            chn_inits,
                            _MS_check_node_message,
                            var_adj_list,
                            check_adj_list,
                            max_iter,
                            :flooding,
                            current_bits,
                            totals,
                            syn,
                            check_to_var_messages_f,
                            var_to_check_messages_f,
                            attenuation,
                        )
                        (!flag && !iszero(out)) && (local_counts_MS[th] += 1;)

                        # reset inputs for next run, but don't re-allocate new memory
                        check_to_var_messages_f[:, :, :] .= 0.0
                        var_to_check_messages_f[:, :, :] .= 0.0

                        # run syndrome-based MS
                        flag, out, _, = _message_passing(
                            H_Int,
                            syn_Int,
                            chn_inits_syn,
                            _MS_check_node_message,
                            var_adj_list,
                            check_adj_list,
                            max_iter,
                            :flooding,
                            current_bits,
                            totals,
                            syn,
                            check_to_var_messages_f,
                            var_to_check_messages_f,
                            0.0,
                        )
                        (!flag && out ≠ err) && (local_counts_MS_syn[th] += 1;)

                        # reset inputs for next run, but don't re-allocate new memory
                        check_to_var_messages_f[:, :, :] .= 0.0
                        var_to_check_messages_f[:, :, :] .= 0.0

                        # run MS with correction
                        flag, out, _ = _message_passing(
                            H_Int,
                            missing,
                            chn_inits,
                            _MS_correction_check_node_message,
                            var_adj_list,
                            check_adj_list,
                            max_iter,
                            :flooding,
                            current_bits,
                            totals,
                            syn,
                            check_to_var_messages_f,
                            var_to_check_messages_f,
                            attenuation,
                        )
                        (!flag || !iszero(out)) && (local_counts_MS_C[th] += 1;)

                        # reset inputs for next run, but don't re-allocate new memory
                        check_to_var_messages_f[:, :, :] .= 0.0
                        var_to_check_messages_f[:, :, :] .= 0.0

                        # run syndrome-based MS with correction
                        flag, out, _, = _message_passing(
                            H_Int,
                            syn_Int,
                            chn_inits_syn,
                            _MS_correction_check_node_message,
                            var_adj_list,
                            check_adj_list,
                            max_iter,
                            :flooding,
                            current_bits,
                            totals,
                            syn,
                            check_to_var_messages_f,
                            var_to_check_messages_f,
                            attenuation,
                        )
                        (!flag || out ≠ err) && (local_counts_MS_C_syn[th] += 1;)

                        # reset inputs for next run, but don't re-allocate new memory
                        check_to_var_messages_f[:, :, :] .= 0.0
                        var_to_check_messages_f[:, :, :] .= 0.0

                        # serial
                        # run SP
                        flag, out, _ = _message_passing(
                            H_Int,
                            missing,
                            chn_inits,
                            _SP_check_node_message_box_plus,
                            var_adj_list,
                            check_adj_list,
                            max_iter,
                            :serial,
                            current_bits,
                            totals,
                            syn,
                            check_to_var_messages_s,
                            var_to_check_messages_s,
                            0.0,
                        )
                        (!flag || !iszero(out)) && (local_counts_SP_s[th] += 1;)

                        # reset inputs for next run, but don't re-allocate new memory
                        check_to_var_messages_s[:, :, :] .= 0.0
                        var_to_check_messages_s[:, :, :] .= 0.0

                        # run syndrome-based SP
                        flag, out, _, = _message_passing(
                            H_Int,
                            syn_Int,
                            chn_inits_syn,
                            _SP_check_node_message_box_plus,
                            var_adj_list,
                            check_adj_list,
                            max_iter,
                            :serial,
                            current_bits,
                            totals,
                            syn,
                            check_to_var_messages_s,
                            var_to_check_messages_s,
                            0.0,
                        )
                        (!flag || out ≠ err) && (local_counts_SP_syn_s[th] += 1;)

                        # reset inputs for next run, but don't re-allocate new memory
                        check_to_var_messages_s[:, :, :] .= 0.0
                        var_to_check_messages_s[:, :, :] .= 0.0

                        # run MS
                        flag, out, _ = _message_passing(
                            H_Int,
                            missing,
                            chn_inits,
                            _MS_check_node_message,
                            var_adj_list,
                            check_adj_list,
                            max_iter,
                            :serial,
                            current_bits,
                            totals,
                            syn,
                            check_to_var_messages_s,
                            var_to_check_messages_s,
                            attenuation,
                        )
                        (!flag && !iszero(out)) && (local_counts_MS_s[th] += 1;)

                        # reset inputs for next run, but don't re-allocate new memory
                        check_to_var_messages_s[:, :, :] .= 0.0
                        var_to_check_messages_s[:, :, :] .= 0.0

                        # run syndrome-based MS
                        flag, out, _, = _message_passing(
                            H_Int,
                            syn_Int,
                            chn_inits_syn,
                            _MS_check_node_message,
                            var_adj_list,
                            check_adj_list,
                            max_iter,
                            :serial,
                            current_bits,
                            totals,
                            syn,
                            check_to_var_messages_s,
                            var_to_check_messages_s,
                            0.0,
                        )
                        (!flag && out ≠ err) && (local_counts_MS_syn_s[th] += 1;)

                        # reset inputs for next run, but don't re-allocate new memory
                        check_to_var_messages_s[:, :, :] .= 0.0
                        var_to_check_messages_s[:, :, :] .= 0.0

                        # run MS with correction
                        flag, out, _ = _message_passing(
                            H_Int,
                            missing,
                            chn_inits,
                            _MS_correction_check_node_message,
                            var_adj_list,
                            check_adj_list,
                            max_iter,
                            :serial,
                            current_bits,
                            totals,
                            syn,
                            check_to_var_messages_s,
                            var_to_check_messages_s,
                            attenuation,
                        )
                        (!flag || !iszero(out)) && (local_counts_MS_C_s[th] += 1;)

                        # reset inputs for next run, but don't re-allocate new memory
                        check_to_var_messages_s[:, :, :] .= 0.0
                        var_to_check_messages_s[:, :, :] .= 0.0

                        # run syndrome-based MS with correction
                        flag, out, _, = _message_passing(
                            H_Int,
                            syn_Int,
                            chn_inits_syn,
                            _MS_correction_check_node_message,
                            var_adj_list,
                            check_adj_list,
                            max_iter,
                            :serial,
                            current_bits,
                            totals,
                            syn,
                            check_to_var_messages_s,
                            var_to_check_messages_s,
                            attenuation,
                        )
                        (!flag || out ≠ err) && (local_counts_MS_C_syn_s[th] += 1;)

                        # reset inputs for next run, but don't re-allocate new memory
                        check_to_var_messages_s[:, :, :] .= 0.0
                        var_to_check_messages_s[:, :, :] .= 0.0
                    end
                end
                verbose && println(local_counts_SP)

                # flooding
                subset_dic_SP[s] = sum(local_counts_SP) / new_num_runs
                subset_dic_SP_syn[s] = sum(local_counts_SP_syn) / new_num_runs
                subset_dic_SP_dec[s] = sum(local_counts_SP_dec) / new_num_runs
                subset_dic_MS[s] = sum(local_counts_MS) / new_num_runs
                subset_dic_MS_syn[s] = sum(local_counts_MS_syn) / new_num_runs
                subset_dic_MS_dec[s] = sum(local_counts_MS_dec) / new_num_runs
                subset_dic_MS_C[s] = sum(local_counts_MS_C) / new_num_runs
                subset_dic_MS_C_syn[s] = sum(local_counts_MS_C_syn) / new_num_runs
                subset_dic_MS_C_dec[s] = sum(local_counts_MS_C_dec) / new_num_runs

                # serial
                subset_dic_SP_s[s] = sum(local_counts_SP_s) / new_num_runs
                subset_dic_SP_syn_s[s] = sum(local_counts_SP_syn_s) / new_num_runs
                subset_dic_SP_dec_s[s] = sum(local_counts_SP_dec_s) / new_num_runs
                subset_dic_MS_s[s] = sum(local_counts_MS_s) / new_num_runs
                subset_dic_MS_syn_s[s] = sum(local_counts_MS_syn_s) / new_num_runs
                subset_dic_MS_dec_s[s] = sum(local_counts_MS_dec_s) / new_num_runs
                subset_dic_MS_C_s[s] = sum(local_counts_MS_C_s) / new_num_runs
                subset_dic_MS_C_syn_s[s] = sum(local_counts_MS_C_syn_s) / new_num_runs
                subset_dic_MS_C_dec_s[s] = sum(local_counts_MS_C_dec_s) / new_num_runs

                # short circuit the rest of the subsets if they are all going to fail
                if subset_dic_SP_s[s] ≥ 0.97
                    if first_fail_flag
                        if !second_fail_flag
                            second_fail_flag = true
                            verbose && println("Second short circuit flag set at subset $s")
                            verbose && println(
                                "Short circuiting importance sampling after subset $s",
                            )
                        end
                    else
                        first_fail_flag = true
                        verbose && println("First short circuit flag set at subset $s")
                    end
                end

                # flooding
                local_counts_SP[:] .= 0
                local_counts_SP_syn[:] .= 0
                local_counts_SP_dec[:] .= 0
                local_counts_MS[:] .= 0
                local_counts_MS_syn[:] .= 0
                local_counts_MS_dec[:] .= 0
                local_counts_MS_C[:] .= 0
                local_counts_MS_C_syn[:] .= 0
                local_counts_MS_C_dec[:] .= 0

                # serial
                local_counts_SP_s[:] .= 0
                local_counts_SP_syn_s[:] .= 0
                local_counts_SP_dec_s[:] .= 0
                local_counts_MS_s[:] .= 0
                local_counts_MS_syn_s[:] .= 0
                local_counts_MS_dec_s[:] .= 0
                local_counts_MS_C_s[:] .= 0
                local_counts_MS_C_syn_s[:] .= 0
                local_counts_MS_C_dec_s[:] .= 0
            end
        end

        for (i, p) in enumerate(noise_impt)
            subsets = find_subsets(n, p, tolerance)
            for s in subsets
                # flooding
                FER_SP[i] += s[2] * subset_dic_SP[s[1]]
                FER_SP_syn[i] += s[2] * subset_dic_SP_syn[s[1]]
                FER_SP_dec[i] += s[2] * subset_dic_SP_dec[s[1]]
                FER_MS[i] += s[2] * subset_dic_MS[s[1]]
                FER_MS_syn[i] += s[2] * subset_dic_MS_syn[s[1]]
                FER_MS_dec[i] += s[2] * subset_dic_MS_dec[s[1]]
                FER_MS_C[i] += s[2] * subset_dic_MS_C[s[1]]
                FER_MS_C_syn[i] += s[2] * subset_dic_MS_C_syn[s[1]]
                FER_MS_C_dec[i] += s[2] * subset_dic_MS_C_dec[s[1]]

                # serial
                FER_SP_s[i] += s[2] * subset_dic_SP_s[s[1]]
                FER_SP_syn_s[i] += s[2] * subset_dic_SP_syn_s[s[1]]
                FER_SP_dec_s[i] += s[2] * subset_dic_SP_dec_s[s[1]]
                FER_MS_s[i] += s[2] * subset_dic_MS_s[s[1]]
                FER_MS_syn_s[i] += s[2] * subset_dic_MS_syn_s[s[1]]
                FER_MS_dec_s[i] += s[2] * subset_dic_MS_dec_s[s[1]]
                FER_MS_C_s[i] += s[2] * subset_dic_MS_C_s[s[1]]
                FER_MS_C_syn_s[i] += s[2] * subset_dic_MS_C_syn_s[s[1]]
                FER_MS_C_dec_s[i] += s[2] * subset_dic_MS_C_dec_s[s[1]]
            end
        end
    end

    # imp_sam_ind = 3
    if imp_sam_ind ≠ len_noise
        verbose && println("Starting direct sampling")
        # direct sampling the rest of the points
        if imp_sam_ind ≠ 0
            # flooding
            local_counts_SP[:] .= 0
            local_counts_SP_syn[:] .= 0
            local_counts_SP_dec[:] .= 0
            local_counts_MS[:] .= 0
            local_counts_MS_syn[:] .= 0
            local_counts_MS_dec[:] .= 0
            local_counts_MS_C[:] .= 0
            local_counts_MS_C_syn[:] .= 0
            local_counts_MS_C_dec[:] .= 0

            # serial
            local_counts_SP_s[:] .= 0
            local_counts_SP_syn_s[:] .= 0
            local_counts_SP_dec_s[:] .= 0
            local_counts_MS_s[:] .= 0
            local_counts_MS_syn_s[:] .= 0
            local_counts_MS_dec_s[:] .= 0
            local_counts_MS_C_s[:] .= 0
            local_counts_MS_C_syn_s[:] .= 0
            local_counts_MS_C_dec_s[:] .= 0
        end

        for i = (imp_sam_ind+1):len_noise
            p = noise[i]
            println("Starting p = $p")
            # initialize everything for p
            chn = MPNoiseModel(noise_model, p)
            init_0 = log((1 - p) / p)
            init_1 = log(p / (1 - p))

            Threads.@threads for th = 1:num_threads
                err_dir = zeros(Int, n)
                syn_Int_dir = zeros(Int, nr)
                erasures_dir = Int[]
                chn_inits_dir = zeros(Float64, n)
                # SP and MS
                H_Int_dir,
                _,
                var_adj_list_dir,
                check_adj_list_dir,
                chn_inits_dir,
                check_to_var_messages_f_dir,
                var_to_check_messages_f_dir,
                current_bits_dir,
                totals_dir,
                syn_dir = _message_passing_init(
                    H,
                    v,
                    chn,
                    max_iter,
                    :SP,
                    chn_inits_dir,
                    :flooding,
                    erasures_dir,
                )
                _,
                _,
                _,
                _,
                _,
                check_to_var_messages_s_dir,
                var_to_check_messages_s_dir,
                _,
                _,
                _ = _message_passing_init(
                    H,
                    v,
                    chn,
                    max_iter,
                    :SP,
                    chn_inits_dir,
                    :serial,
                    erasures_dir,
                )
                # for syndrome-based we need to change the inits
                chn_inits_syn_dir = [init_0 for _ = 1:n]

                for _ = 1:runs_per_thread
                    # sample
                    @inbounds for j = 1:n
                        rand(dist) ≤ p ? (err_dir[j] = 1; chn_inits_dir[j] = init_1;) :
                        (err_dir[j] = 0; chn_inits_dir[j] = init_0;)
                    end
                    # err = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                    # chn_inits = [init_0 for _ in 1:n]
                    # chn_inits[81] = init_1
                    # th == 1 && println(err)
                    # println("thread: $th, weight: $wt")
                    LinearAlgebra.mul!(syn_Int_dir, H_Int_dir, err_dir)
                    @inbounds @simd for i = 1:nr
                        syn_Int_dir[i] %= 2
                    end
                    # pick a random bit to decimate
                    # bit = rand(1:n)
                    # decimated_bits = [bit]
                    # # for now, setup as a genie-assisted decoder
                    # decimated_values = [err[bit]]

                    # flooding
                    # run SP
                    flag_dir, out_dir, _ = _message_passing(
                        H_Int_dir,
                        missing,
                        chn_inits_dir,
                        _SP_check_node_message_box_plus,
                        var_adj_list_dir,
                        check_adj_list_dir,
                        max_iter,
                        :flooding,
                        current_bits_dir,
                        totals_dir,
                        syn_dir,
                        check_to_var_messages_f_dir,
                        var_to_check_messages_f_dir,
                        0.0,
                    )
                    (!flag_dir || !iszero(out_dir)) && (local_counts_SP[th] += 1;)
                    # if th == 1 && !iszero(out)
                    #     println("out = $out;")
                    #     println("err = $err;")
                    #     return
                    # end
                    # println(flag)
                    # println(out)
                    # return

                    # reset inputs for next run, but don't re-allocate new memory
                    check_to_var_messages_f_dir[:, :, :] .= 0.0
                    var_to_check_messages_f_dir[:, :, :] .= 0.0

                    # run syndrome-based SP
                    flag_dir, out_dir, _, = _message_passing(
                        H_Int_dir,
                        syn_Int_dir,
                        chn_inits_syn_dir,
                        _SP_check_node_message_box_plus,
                        var_adj_list_dir,
                        check_adj_list_dir,
                        max_iter,
                        :flooding,
                        current_bits_dir,
                        totals_dir,
                        syn_dir,
                        check_to_var_messages_f_dir,
                        var_to_check_messages_f_dir,
                        0.0,
                    )
                    (!flag_dir || out_dir ≠ err_dir) && (local_counts_SP_syn[th] += 1;)

                    # reset inputs for next run, but don't re-allocate new memory
                    check_to_var_messages_f_dir[:, :, :] .= 0.0
                    var_to_check_messages_f_dir[:, :, :] .= 0.0

                    # run MS
                    flag_dir, out_dir, _ = _message_passing(
                        H_Int_dir,
                        missing,
                        chn_inits_dir,
                        _MS_check_node_message,
                        var_adj_list_dir,
                        check_adj_list_dir,
                        max_iter,
                        :flooding,
                        current_bits_dir,
                        totals_dir,
                        syn_dir,
                        check_to_var_messages_f_dir,
                        var_to_check_messages_f_dir,
                        attenuation,
                    )
                    (!flag_dir && !iszero(out_dir)) && (local_counts_MS[th] += 1;)

                    # reset inputs for next run, but don't re-allocate new memory
                    check_to_var_messages_f_dir[:, :, :] .= 0.0
                    var_to_check_messages_f_dir[:, :, :] .= 0.0

                    # run syndrome-based MS
                    flag_dir, out_dir, _, = _message_passing(
                        H_Int_dir,
                        syn_Int_dir,
                        chn_inits_syn_dir,
                        _MS_check_node_message,
                        var_adj_list_dir,
                        check_adj_list_dir,
                        max_iter,
                        :flooding,
                        current_bits_dir,
                        totals_dir,
                        syn_dir,
                        check_to_var_messages_f_dir,
                        var_to_check_messages_f_dir,
                        0.0,
                    )
                    (!flag_dir && out_dir ≠ err_dir) && (local_counts_MS_syn[th] += 1;)

                    # reset inputs for next run, but don't re-allocate new memory
                    check_to_var_messages_f_dir[:, :, :] .= 0.0
                    var_to_check_messages_f_dir[:, :, :] .= 0.0

                    # run MS with correction
                    flag_dir, out_dir, _ = _message_passing(
                        H_Int_dir,
                        missing,
                        chn_inits_dir,
                        _MS_correction_check_node_message,
                        var_adj_list_dir,
                        check_adj_list_dir,
                        max_iter,
                        :flooding,
                        current_bits_dir,
                        totals_dir,
                        syn_dir,
                        check_to_var_messages_f_dir,
                        var_to_check_messages_f_dir,
                        attenuation,
                    )
                    (!flag_dir || !iszero(out_dir)) && (local_counts_MS_C[th] += 1;)

                    # reset inputs for next run, but don't re-allocate new memory
                    check_to_var_messages_f_dir[:, :, :] .= 0.0
                    var_to_check_messages_f_dir[:, :, :] .= 0.0

                    # run syndrome-based MS with correction
                    flag_dir, out_dir, _, = _message_passing(
                        H_Int_dir,
                        syn_Int_dir,
                        chn_inits_syn_dir,
                        _MS_correction_check_node_message,
                        var_adj_list_dir,
                        check_adj_list_dir,
                        max_iter,
                        :flooding,
                        current_bits_dir,
                        totals_dir,
                        syn_dir,
                        check_to_var_messages_f_dir,
                        var_to_check_messages_f_dir,
                        attenuation,
                    )
                    (!flag_dir || out_dir ≠ err_dir) && (local_counts_MS_C_syn[th] += 1;)

                    # reset inputs for next run, but don't re-allocate new memory
                    check_to_var_messages_f_dir[:, :, :] .= 0.0
                    var_to_check_messages_f_dir[:, :, :] .= 0.0

                    # serial
                    # run SP
                    flag_dir, out_dir, _ = _message_passing(
                        H_Int_dir,
                        missing,
                        chn_inits_dir,
                        _SP_check_node_message_box_plus,
                        var_adj_list_dir,
                        check_adj_list_dir,
                        max_iter,
                        :serial,
                        current_bits_dir,
                        totals_dir,
                        syn_dir,
                        check_to_var_messages_s_dir,
                        var_to_check_messages_s_dir,
                        0.0,
                    )
                    (!flag_dir || !iszero(out_dir)) && (local_counts_SP_s[th] += 1;)

                    # reset inputs for next run, but don't re-allocate new memory
                    check_to_var_messages_s_dir[:, :, :] .= 0.0
                    var_to_check_messages_s_dir[:, :, :] .= 0.0

                    # run syndrome-based SP
                    flag_dir, out_dir, _, = _message_passing(
                        H_Int_dir,
                        syn_Int_dir,
                        chn_inits_syn_dir,
                        _SP_check_node_message_box_plus,
                        var_adj_list_dir,
                        check_adj_list_dir,
                        max_iter,
                        :serial,
                        current_bits_dir,
                        totals_dir,
                        syn_dir,
                        check_to_var_messages_s_dir,
                        var_to_check_messages_s_dir,
                        0.0,
                    )
                    (!flag_dir || out_dir ≠ err_dir) && (local_counts_SP_syn_s[th] += 1;)

                    # reset inputs for next run, but don't re-allocate new memory
                    check_to_var_messages_s_dir[:, :, :] .= 0.0
                    var_to_check_messages_s_dir[:, :, :] .= 0.0

                    # run MS
                    flag_dir, out_dir, _ = _message_passing(
                        H_Int_dir,
                        missing,
                        chn_inits_dir,
                        _MS_check_node_message,
                        var_adj_list_dir,
                        check_adj_list_dir,
                        max_iter,
                        :serial,
                        current_bits_dir,
                        totals_dir,
                        syn_dir,
                        check_to_var_messages_s_dir,
                        var_to_check_messages_s_dir,
                        attenuation,
                    )
                    (!flag_dir && !iszero(out_dir)) && (local_counts_MS_s[th] += 1;)

                    # reset inputs for next run, but don't re-allocate new memory
                    check_to_var_messages_s_dir[:, :, :] .= 0.0
                    var_to_check_messages_s_dir[:, :, :] .= 0.0

                    # run syndrome-based MS
                    flag_dir, out_dir, _, = _message_passing(
                        H_Int_dir,
                        syn_Int_dir,
                        chn_inits_syn_dir,
                        _MS_check_node_message,
                        var_adj_list_dir,
                        check_adj_list_dir,
                        max_iter,
                        :serial,
                        current_bits_dir,
                        totals_dir,
                        syn_dir,
                        check_to_var_messages_s_dir,
                        var_to_check_messages_s_dir,
                        0.0,
                    )
                    (!flag_dir && out_dir ≠ err_dir) && (local_counts_MS_syn_s[th] += 1;)

                    # reset inputs for next run, but don't re-allocate new memory
                    check_to_var_messages_s_dir[:, :, :] .= 0.0
                    var_to_check_messages_s_dir[:, :, :] .= 0.0

                    # run MS with correction
                    flag_dir, out_dir, _ = _message_passing(
                        H_Int_dir,
                        missing,
                        chn_inits_dir,
                        _MS_correction_check_node_message,
                        var_adj_list_dir,
                        check_adj_list_dir,
                        max_iter,
                        :serial,
                        current_bits_dir,
                        totals_dir,
                        syn_dir,
                        check_to_var_messages_s_dir,
                        var_to_check_messages_s_dir,
                        attenuation,
                    )
                    (!flag_dir || !iszero(out_dir)) && (local_counts_MS_C_s[th] += 1;)

                    # reset inputs for next run, but don't re-allocate new memory
                    check_to_var_messages_s_dir[:, :, :] .= 0.0
                    var_to_check_messages_s_dir[:, :, :] .= 0.0

                    # run syndrome-based MS with correction
                    flag_dir, out_dir, _, = _message_passing(
                        H_Int_dir,
                        syn_Int_dir,
                        chn_inits_syn_dir,
                        _MS_correction_check_node_message,
                        var_adj_list_dir,
                        check_adj_list_dir,
                        max_iter,
                        :serial,
                        current_bits_dir,
                        totals_dir,
                        syn_dir,
                        check_to_var_messages_s_dir,
                        var_to_check_messages_s_dir,
                        attenuation,
                    )
                    (!flag_dir || out_dir ≠ err_dir) && (local_counts_MS_C_syn_s[th] += 1;)

                    # reset inputs for next run, but don't re-allocate new memory
                    check_to_var_messages_s_dir[:, :, :] .= 0.0
                    var_to_check_messages_s_dir[:, :, :] .= 0.0
                end
            end
            # println(local_counts_SP)
            # flooding
            @inbounds FER_SP[i] = sum(local_counts_SP) / new_num_runs
            # println("FER_SP[$i] = $(FER_SP[i])")
            # return
            @inbounds FER_SP_syn[i] = sum(local_counts_SP_syn) / new_num_runs
            @inbounds FER_SP_dec[i] = sum(local_counts_SP_dec) / new_num_runs
            @inbounds FER_MS[i] = sum(local_counts_MS) / new_num_runs
            @inbounds FER_MS_syn[i] = sum(local_counts_MS_syn) / new_num_runs
            @inbounds FER_MS_dec[i] = sum(local_counts_MS_dec) / new_num_runs
            @inbounds FER_MS_C[i] = sum(local_counts_MS_C) / new_num_runs
            @inbounds FER_MS_C_syn[i] = sum(local_counts_MS_C_syn) / new_num_runs
            @inbounds FER_MS_C_dec[i] = sum(local_counts_MS_C_dec) / new_num_runs

            # serial
            @inbounds FER_SP_s[i] = sum(local_counts_SP_s) / new_num_runs
            @inbounds FER_SP_syn_s[i] = sum(local_counts_SP_syn_s) / new_num_runs
            @inbounds FER_SP_dec_s[i] = sum(local_counts_SP_dec_s) / new_num_runs
            @inbounds FER_MS_s[i] = sum(local_counts_MS_s) / new_num_runs
            @inbounds FER_MS_syn_s[i] = sum(local_counts_MS_syn_s) / new_num_runs
            @inbounds FER_MS_dec_s[i] = sum(local_counts_MS_dec_s) / new_num_runs
            @inbounds FER_MS_C_s[i] = sum(local_counts_MS_C_s) / new_num_runs
            @inbounds FER_MS_C_syn_s[i] = sum(local_counts_MS_C_syn_s) / new_num_runs
            @inbounds FER_MS_C_dec_s[i] = sum(local_counts_MS_C_dec_s) / new_num_runs

            # flooding
            local_counts_SP[:] .= 0
            local_counts_SP_syn[:] .= 0
            local_counts_SP_dec[:] .= 0
            local_counts_MS[:] .= 0
            local_counts_MS_syn[:] .= 0
            local_counts_MS_dec[:] .= 0
            local_counts_MS_C[:] .= 0
            local_counts_MS_C_syn[:] .= 0
            local_counts_MS_C_dec[:] .= 0

            # serial
            local_counts_SP_s[:] .= 0
            local_counts_SP_syn_s[:] .= 0
            local_counts_SP_dec_s[:] .= 0
            local_counts_MS_s[:] .= 0
            local_counts_MS_syn_s[:] .= 0
            local_counts_MS_dec_s[:] .= 0
            local_counts_MS_C_s[:] .= 0
            local_counts_MS_C_syn_s[:] .= 0
            local_counts_MS_C_dec_s[:] .= 0

            println("Finished p = $p")
        end
    end
    return FER_SP,
    FER_SP_syn,
    FER_SP_dec,
    FER_MS,
    FER_MS_syn,
    FER_MS_dec,
    FER_MS_C,
    FER_MS_C_syn,
    FER_MS_C_dec,
    FER_SP_s,
    FER_SP_syn_s,
    FER_SP_dec_s,
    FER_MS_s,
    FER_MS_syn_s,
    FER_MS_dec_s,
    FER_MS_C_s,
    FER_MS_C_syn_s,
    FER_MS_C_dec_s,
    subset_dic_SP,
    subset_dic_SP_syn,
    subset_dic_SP_dec,
    subset_dic_MS,
    subset_dic_MS_syn,
    subset_dic_MS_dec,
    subset_dic_MS_C,
    subset_dic_MS_C_syn,
    subset_dic_MS_C_dec,
    subset_dic_SP_s,
    subset_dic_SP_syn_s,
    subset_dic_SP_dec_s,
    subset_dic_MS_s,
    subset_dic_MS_syn_s,
    subset_dic_MS_dec_s,
    subset_dic_MS_C_s,
    subset_dic_MS_C_syn_s,
    subset_dic_MS_C_dec_s
end

function single_decoder_test(H::CTMatrixTypes)
    # initial parameters
    noise = [0.01]
    # noise = [0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05]
    # left = 0.001
    # right = 0.04
    # intervals = 10
    # Δ = (right - left) / (intervals - 1)
    # noise = collect(left:Δ:right)
    len_noise = length(noise)
    noise_model = :BSC
    num_runs = 15000
    nr, n = size(H)
    # it suffices to decode the zero codeword
    v = zero_matrix(base_ring(H), 1, n)
    max_iter = 50
    schedule = :flooding
    dist = Uniform(0.0, 1.0)
    algorithm = :manual
    guided_rounds = 10
    attenuation = 0.5
    num_threads = Threads.nthreads()
    runs_per_thread = cld(num_runs, num_threads)
    new_num_runs = runs_per_thread * num_threads
    # println(num_threads, ", ", runs_per_thread, ", ", new_num_runs)

    FER = zeros(Float64, len_noise)
    for (i, p) in enumerate(noise)
        println("Starting p = $p")
        # initialize everything for p
        chn = MPNoiseModel(noise_model, p)
        init_0 = log((1 - p) / p)
        init_1 = log(p / (1 - p))
        local_counts = zeros(Int, num_threads)
        Threads.@threads for th = 1:num_threads
            err = zeros(Int, n)
            syn_Int = zeros(Int, nr)
            erasures = Int[]
            chn_inits = zeros(Float64, n)
            # SP/MS for BSC
            H_Int,
            _,
            var_adj_list,
            check_adj_list,
            chn_inits,
            check_to_var_messages,
            var_to_check_messages,
            current_bits,
            totals,
            syn = _message_passing_init(
                H,
                v,
                chn,
                max_iter,
                :SP,
                chn_inits,
                schedule,
                erasures,
            )
            # chn_inits_syn = [init_0 for _ in 1:n]

            for _ = 1:runs_per_thread
                # sample
                # @inbounds for j in 1:n
                #     rand(dist) ≤ p ? (err[j] = 1; chn_inits[j] = init_1;) : (err[j] = 0; chn_inits[j] = init_0;)
                # end
                err = zeros(Int, 254)
                err[81] = 1
                chn_inits = [init_0 for _ = 1:n]
                chn_inits[81] = init_1

                # @inbounds for j in 1:n
                #     rand(dist) ≤ p ? (err[j] = 1;) : (err[j] = 0;)
                # end
                # LinearAlgebra.mul!(syn_Int, H_Int, err)
                # @inbounds @simd for i in 1:nr
                #     syn_Int[i] %= 2
                # end
                # pick a random bit to decimate
                # bit = rand(1:n)
                # decimated_bits = [bit]
                # # for now, setup as a genie-assisted decoder
                # decimated_values = [err[bit]]

                # SP
                flag, out, _ = _message_passing(
                    H_Int,
                    missing,
                    chn_inits,
                    _SP_check_node_message_box_plus,
                    var_adj_list,
                    check_adj_list,
                    max_iter,
                    schedule,
                    current_bits,
                    totals,
                    syn,
                    check_to_var_messages,
                    var_to_check_messages,
                    0.0,
                )
                (!flag || !iszero(out)) && (local_counts[th] += 1;)
                println(flag)
                println(out)
                return
                # MS
                # flag, out, _ = _message_passing(H_Int, missing, chn_inits, _MS_check_node_message,
                #     var_adj_list, check_adj_list, max_iter, schedule, current_bits, totals, syn,
                #     check_to_var_messages, var_to_check_messages, attenuation)
                # (!flag || !iszero(out)) && (count += 1;)
                # SP syn
                # flag, out, _, = _message_passing(H_Int, syn_Int, chn_inits_syn,
                #     _SP_check_node_message_box_plus, var_adj_list, check_adj_list, max_iter, schedule,
                #     current_bits, totals, syn, check_to_var_messages, var_to_check_messages, 0.0)
                # (!flag || out ≠ err) && (count += 1;)
                # !flag && println("err = $err, out = $out, - $flag")
                # print("$iter, ")
                # run SP with decimation
                # flag, out, _ = _message_passing_decimation(H_Int, err, chn_inits,
                #     _SP_check_node_message_box_plus, var_adj_list, check_adj_list, max_iter, :SP, schedule,
                #     decimated_bits, decimated_values, current_bits, totals, syn, check_to_var_messages,
                #     var_to_check_messages, 0, 0.0, algorithm, guided_rounds)
                # (!flag || !iszero(out)) && (count += 1;)
                # run MS with decimation
                # flag, out, _= _message_passing_decimation(H_Int, err, chn_inits, 
                #     _MS_check_node_message, var_adj_list, check_adj_list, max_iter, :MS, schedule,
                #     decimated_bits, decimated_values, current_bits, totals, syn, check_to_var_messages, 
                #     var_to_check_messages, 0, attenuation, algorithm, guided_rounds)
                # (!flag || !iszero(out)) && (count += 1;)
                # run syndrome-based MS
                # flag, out, _= _message_passing(H_Int, syn_Int, chn_inits_syn, 
                #     _MS_check_node_message, var_adj_list, check_adj_list, max_iter, schedule, 
                #     current_bits, totals, syn, check_to_var_messages, var_to_check_messages, 
                #     attenuation)
                # (!flag && out ≠ err) && (count += 1;)

                # reset inputs for next run, but don't re-allocate new memory
                check_to_var_messages[:, :, :] .= 0.0
                var_to_check_messages[:, :, :] .= 0.0
            end
        end
        println(local_counts)
        @inbounds FER[i] = sum(local_counts) / new_num_runs
        println("FER = $(FER[i])")
        println("Finished p = $p")
    end
    return FER
end

#     for k in eachindex(noise)

#         """
#         first thing I want to do here is defind the important and direct sampling split
#         to do this, I need to find the number of subsets for each p
#         then I need to estimate the standard deviation at each point
#         then I need to determine the number of samples at each point for the confidence interval
#         then compute the speedup factor at the overlap to check it's a good transition?
#         """

#         # p[i] is the probability of having i - 1 bit errors
#         temp = BigFloat(noise[k])
#         p = BigFloat[temp^i * (1 - temp)^(n - i) * binomial(big(n), big(i)) / sum(temp^j * (1 -
#             temp)^(n - j) * binomial(big(n), big(j)) for j in 0:n) for i in 0:n]
#         p_partial_sum = [sum(p[j] for j in 1:i) for i in 1:length(p)]
#         max_num_err = max(findfirst(p_partial_sum .>= 1 - BigFloat("1e-9")) - 1, 6)
#         # @show max_num_err
#         ε[k] = 1 - p_partial_sum[max_num_err + 1]

#         """
#         having settled all of that, let's make a loop for the importance sampled points
#         do the MP init here
#         inside the loop, sample an error and run the same one for every decoder
#         save thread local totals
#         at the end, run a thread safe combination of local totals
#         save at the appropriate parts of FER/BER vectors
#         """


#         """
#         now do the same with the points from direct sampling
#         use the same inits and memory from above
#         repeat entire process
#         return an array for every decoder
#         """
#         chn = MPNoiseModel(noise_type, noise[k])
#         w = noise_type == :BSC ? zeros(Int, n) : ones(n)
#         if decoder == :LP
#             noisemodel = BSC(noise[k])
#         else
#             H_Int, _, var_adj_list, check_adj_list = _message_passing_init(H, w, chn, max_iter, decoder, 2)
#         end

#         FEtotal = zeros(Int, max_num_err)
#         BEtotal = zeros(Int, max_num_err)
#         FER[k] = p[1]
#         BER[k] = p[1]

#         for e in 1:max_num_err

#             # importance sampling:
#             # cld(x, y) = div(x, y, RoundUp)
#             num_runs_for_e = Int(cld(num_runs * p[e + 1], sum(p[i] for i in 2:max_num_err + 1)))

#             # direct sampling:
#             # num_runs_for_e = Int(cld(num_runs, max_num_err))

#             for i in 1:num_runs_for_e
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
#             FER[k] += p[e+1] * (num_runs_for_e - FEtotal[e]) / num_runs_for_e
#             BER[k] += p[e+1] * (num_runs_for_e * n - BEtotal[e]) / (num_runs_for_e * n)
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
#     @time F1, B1, e1 = CodingTheory.decoder_simulation(H, :SP, :BSC, p, 10, 5000, 1, 1.0);
#     @time F2, B2, e2 = CodingTheory.decoder_simulation(H, :MS, :BSC, p, 10, 5000, 1, 1.0);
#     @time F3, B3, e3 = CodingTheory.decoder_simulation(H, :MS, :BSC, p, 10, 5000, 1, 0.5);
#     @time F4, B4, e4 = CodingTheory.decoder_simulation(H, :MS, :BSC, p, 10, 5000, 1, 0.1);
#     @time F5, B5, e5 = CodingTheory.decoder_simulation(H, :LP, :BSC, p, 10, 5000, 1);

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
#     # FER, BER = decoder_simulation(H, :SP, :BAWGNC, p, 100, 100, 123);
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
