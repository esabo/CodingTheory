# Copyright (c) 2024 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
# Methods
#############################

function _message_passing_init_CSS(
    S::AbstractStabilizerCode,
    chn::Dict{Int,Float64},
    max_iter::Int,
    X_chn_inits::Union{Missing,Vector{Float64}},
    Z_chn_inits::Union{Missing,Vector{Float64}},
    schedule::Symbol,
    erasures::Vector{Int},
)

    # initial checks
    X_stabs = _Flint_matrix_to_Julia_int_matrix(X_stabilizers(S))
    Z_stabs = _Flint_matrix_to_Julia_int_matrix(Z_stabilizers(S))
    logs = _Flint_matrix_to_Julia_int_matrix(logicals_matrix(S))
    X_num_check, n = size(X_stabs)
    X_num_check > 0 && n > 0 || throw(ArgumentError("Input matrix of improper dimension"))
    Z_num_check = size(Z_stabs, 1)
    Z_num_check > 0 || throw(ArgumentError("Input matrix of improper dimension"))
    2 <= max_iter || throw(DomainError("Number of maximum iterations must be at least two"))

    X_check_adj_list = [Int[] for _ = 1:X_num_check]
    X_var_adj_list = [Int[] for _ = 1:n]
    for r = 1:X_num_check
        for c = 1:n
            if !iszero(X_stabs[r, c])
                push!(X_check_adj_list[r], c)
                push!(X_var_adj_list[c], r)
            end
        end
    end

    Z_check_adj_list = [Int[] for _ = 1:Z_num_check]
    Z_var_adj_list = [Int[] for _ = 1:n]
    for r = 1:Z_num_check
        for c = 1:n
            if !iszero(Z_stabs[r, c])
                push!(Z_check_adj_list[r], c)
                push!(Z_var_adj_list[c], r)
            end
        end
    end

    Pauli_weights = collect(values(chn))
    if ismissing(X_chn_inits)
        # CSS - so initial condition is log(p_I / p_X)
        X_chn_inits_2 = [log(Pauli_weights[1] / Pauli_weights[2]) for _ = 1:n]
    else
        length(X_chn_inits) ≠ n && throw(ArgumentError("Channel inputs has wrong size"))
        X_chn_inits_2 = X_chn_inits
    end
    if ismissing(Z_chn_inits)
        # CSS - so initial condition is log(p_I / p_Z)
        Z_chn_inits_2 = [log(Pauli_weights[1] / Pauli_weights[4]) for _ = 1:n]
    else
        length(Z_chn_inits) ≠ n && throw(ArgumentError("Channel inputs has wrong size"))
        Z_chn_inits_2 = Z_chn_inits
    end

    # could reduce this stuff to UInt8 if needed
    current_bits = zeros(Int, n)
    totals = zeros(Float64, n)
    # X error vector
    X_err = zeros(Int, n)
    # Z error vector
    Z_err = zeros(Int, n)
    # X syndrome vector
    X_syn = zeros(Int, size(X_stabs, 1))
    # Z syndrome vector
    Z_syn = zeros(Int, size(Z_stabs, 1))

    if schedule == :serial
        X_check_to_var_messages = zeros(Float64, X_num_check, n, 1)
        X_var_to_check_messages = zeros(Float64, n, X_num_check, 1)
        Z_check_to_var_messages = zeros(Float64, Z_num_check, n, 1)
        Z_var_to_check_messages = zeros(Float64, n, Z_num_check, 1)
    else
        X_check_to_var_messages = zeros(Float64, X_num_check, n, 2)
        X_var_to_check_messages = zeros(Float64, n, X_num_check, 2)
        Z_check_to_var_messages = zeros(Float64, Z_num_check, n, 2)
        Z_var_to_check_messages = zeros(Float64, n, Z_num_check, 2)
    end

    if !isempty(erasures)
        all(1 ≤ bit ≤ num_var for bit in erasures) ||
            throw(ArgumentError("Invalid bit index in erasures"))
        @inbounds for i in erasures
            chn_inits_2[i] = 0.0
        end
    end

    return X_stabs,
    Z_stabs,
    logs,
    X_err,
    Z_err,
    X_var_adj_list,
    X_check_adj_list,
    Z_var_adj_list,
    Z_check_adj_list,
    X_chn_inits_2,
    Z_chn_inits_2,
    X_check_to_var_messages,
    X_var_to_check_messages,
    Z_check_to_var_messages,
    Z_var_to_check_messages,
    current_bits,
    totals,
    X_syn,
    Z_syn
end

CSS_decoder_test(S::T) where {T<:AbstractStabilizerCode} = CSS_decoder_test(CSSTrait(T), S)
function CSS_decoder_test(::IsCSS, S::AbstractStabilizerCode)
    # noise details
    # noise = [0.01]
    noise = [0.02, 0.03, 0.04, 0.05]
    # noise = [0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05]
    # left = 0.001
    # right = 0.04
    # intervals = 10
    # Δ = (right - left) / (intervals - 1)
    # noise = collect(left:Δ:right)
    len_noise = length(noise)

    # BP details
    n = length(S)
    num_runs = 10000
    max_iter = 50
    # schedule = :parallel
    schedule = :serial
    erasures = Int[]

    # multi-threading details
    num_threads = Threads.nthreads()
    runs_per_thread = cld(num_runs, num_threads)
    new_num_runs = runs_per_thread * num_threads
    println(num_threads, ", ", runs_per_thread, ", ", new_num_runs)

    # a place to store results
    FER = zeros(Float64, len_noise)
    X_local_iters = [Int[] for _ = 1:num_threads]
    Z_local_iters = [Int[] for _ = 1:num_threads]
    for (i, p) in enumerate(noise)
        println("Starting p = $p")

        # setup noise model
        chn = Dict(0 => 1 - p, 1 => p / 3, 2 => p / 3, 3 => p / 3)
        Pauli_types = collect(keys(chn))
        Pauli_weights = Weights(collect(values(chn)))
        Pauli_op = [0]

        # a place to store results
        local_counts = zeros(Int, num_threads)
        Threads.@threads for th = 1:num_threads
            X_stabs_dir,
            Z_stabs_dir,
            logs_dir,
            X_err_dir,
            Z_err_dir,
            X_var_adj_list_dir,
            X_check_adj_list_dir,
            Z_var_adj_list_dir,
            Z_check_adj_list_dir,
            X_chn_inits_dir,
            Z_chn_inits_dir,
            X_check_to_var_messages_dir,
            X_var_to_check_messages_dir,
            Z_check_to_var_messages_dir,
            Z_var_to_check_messages_dir,
            current_bits_dir,
            totals_dir,
            X_syn_dir,
            Z_syn_dir = _message_passing_init_CSS(
                S,
                chn,
                max_iter,
                missing,
                missing,
                schedule,
                erasures,
            )

            nr_X = size(X_stabs_dir, 1)
            nr_Z = size(Z_stabs_dir, 1)
            nr_logs = size(logs_dir, 1)
            X_meas_syn_dir = zeros(Int, nr_X)
            Z_meas_syn_dir = zeros(Int, nr_Z)

            for _ = 1:runs_per_thread
                # sample
                @inbounds @simd for j = 1:n
                    sample!(Pauli_types, Pauli_weights, Pauli_op)
                    # println(Pauli_op)
                    if Pauli_op[1] == 0
                        # I
                        X_err_dir[j] = 0
                        Z_err_dir[j] = 0
                    elseif Pauli_op[1] == 1
                        # X
                        X_err_dir[j] = 1
                        Z_err_dir[j] = 0
                    elseif Pauli_op[1] == 2
                        # Y
                        X_err_dir[j] = 1
                        Z_err_dir[j] = 1
                    else
                        # Z
                        X_err_dir[j] = 0
                        Z_err_dir[j] = 1
                    end
                end
                # println(sum(X_err_dir))
                # println(sum(Z_err_dir))

                # TODO do some timing on this reset vs map or broadcast
                # X syndrome
                LinearAlgebra.mul!(X_meas_syn_dir, X_stabs_dir, Z_err_dir)
                @inbounds @simd for i = 1:nr_X
                    X_meas_syn_dir[i] %= 2
                end
                # println(X_syn_dir)

                # decode X
                X_flag_dir, X_out_dir, X_iter_dir, = _message_passing(
                    X_stabs_dir,
                    X_meas_syn_dir,
                    X_chn_inits_dir,
                    _SP_check_node_message_box_plus,
                    X_var_adj_list_dir,
                    X_check_adj_list_dir,
                    max_iter,
                    schedule,
                    current_bits_dir,
                    totals_dir,
                    X_syn_dir,
                    X_check_to_var_messages_dir,
                    X_var_to_check_messages_dir,
                    0.0,
                )
                # shouldn't matter if I mod 2 this now or later?
                Z_err_dir += X_out_dir
                X_flag_dir && push!(X_local_iters[th], X_iter_dir)
                # println(X_flag_dir, " ", X_out_dir, " ", X_iter_dir)

                # Z syndrome
                LinearAlgebra.mul!(Z_meas_syn_dir, Z_stabs_dir, X_err_dir)
                @inbounds @simd for i = 1:nr_Z
                    Z_meas_syn_dir[i] %= 2
                end
                # println(Z_syn_dir)

                # decode Z
                Z_flag_dir, Z_out_dir, Z_iter_dir, = _message_passing(
                    Z_stabs_dir,
                    Z_meas_syn_dir,
                    Z_chn_inits_dir,
                    _SP_check_node_message_box_plus,
                    Z_var_adj_list_dir,
                    Z_check_adj_list_dir,
                    max_iter,
                    schedule,
                    current_bits_dir,
                    totals_dir,
                    Z_syn_dir,
                    Z_check_to_var_messages_dir,
                    Z_var_to_check_messages_dir,
                    0.0,
                )
                X_err_dir += Z_out_dir
                Z_flag_dir && push!(Z_local_iters[th], Z_iter_dir)
                # println(Z_flag_dir, " ", Z_out_dir, " ", Z_iter_dir)
                # return

                if X_flag_dir && Z_flag_dir
                    # converged
                    # this implements !iszero(hcat(Z_err, -X_err) * transpose(logs)) but avoids the
                    # allocations and stops without having the compute the entire product when the
                    # answer is known
                    # temp = vcat(Z_err_dir, -X_err_dir)' * transpose(logs_dir)
                    # if !iszero(temp .% 2)
                    #     local_counts[th] += 1
                    # end
                    @inbounds for j = 1:nr_logs
                        # introduced a logical error
                        iseven(
                            dot(view(logs_dir, j, 1:n), Z_err_dir) -
                            dot(view(logs_dir, j, (n+1):(2*n)), X_err_dir),
                        ) || (local_counts[th] += 1; break;)
                    end
                else
                    # did not converge
                    local_counts[th] += 1
                end

                # reset inputs for next run, but don't re-allocate new memory
                X_check_to_var_messages_dir[:, :, :] .= 0.0
                X_var_to_check_messages_dir[:, :, :] .= 0.0
                Z_check_to_var_messages_dir[:, :, :] .= 0.0
                Z_var_to_check_messages_dir[:, :, :] .= 0.0
            end
        end
        println(local_counts)
        @inbounds FER[i] = sum(local_counts) / new_num_runs
        println("FER = $(FER[i])")
        println("Finished p = $p")
    end
    return FER,
    StatsBase.countmap(reduce(vcat, X_local_iters)),
    StatsBase.countmap(reduce(vcat, Z_local_iters))
end
CSS_decoder_test(::IsNotCSS, S::AbstractStabilizerCode) =
    throw(ArgumentError("CSS decoders are only valid for CSS codes"))

# using Bayes from X to update the Z priors
CSS_decoder_with_Bayes(S::T; verbose::Bool = true) where {T<:AbstractStabilizerCode} =
    CSS_decoder_with_Bayes(CSSTrait(T), S, verbose)
function CSS_decoder_with_Bayes(::IsCSS, S::AbstractStabilizerCode, verbose::Bool)
    # noise details
    # noise = [0.01]
    # these are the markers from the paper
    # noise = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14]
    # noise = [0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14]
    noise = 10 .^ (-3:0.25:-1);
    # left = 1e-2
    # right = 1e-1
    # intervals = 11
    # Δ = (right - left) / (intervals - 1)
    # noise = collect(left:Δ:right)
    len_noise = length(noise)

    # BP details
    n = length(S)
    num_runs = 20000
    max_iter = 32
    # schedule = :parallel
    # schedule = :serial
    schedule = :layered
    X_layers = layered_schedule(S.X_stabs)
    Z_layers = layered_schedule(S.Z_stabs)
    erasures = Int[]
    tolerance = 1e-12

    # multi-threading details
    num_threads = Threads.nthreads()
    runs_per_thread = cld(num_runs, num_threads)
    new_num_runs = runs_per_thread * num_threads
    verbose && println(
        "Number of threads: $num_threads, runs per thread: $runs_per_thread, new number of runs: $new_num_runs",
    )

    # preallocate places to store results
    FER = zeros(Float64, len_noise)
    local_log_failure_counts = zeros(Int, num_threads)
    local_conv_failure_counts = zeros(Int, num_threads)
    X_local_iters = [Dict{Int,Int}() for _ = 1:num_threads]
    Z_local_iters = [Dict{Int,Int}() for _ = 1:num_threads]
    for i = 1:num_threads
        sizehint!(X_local_iters[i], max_iter)
        sizehint!(Z_local_iters[i], max_iter)
        for j = 1:max_iter
            X_local_iters[i][j] = 0
            Z_local_iters[i][j] = 0
        end
    end

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

    # imp_sam_ind = 2
    if imp_sam_ind ≠ 0
        noise_impt = noise[1:imp_sam_ind]
        min_subset, max_subset = min_max_subsets(n, noise, tolerance)
        # TODO: fix
        min_subset = 0
        max_subset = 36
        subsets = [i for i = min_subset:max_subset]
        subset_dic = Dict{Int,Float64}()
        fail_flag = false

        if verbose
            println("Starting importance sampling on $(length(noise_impt)) noise values.")
            println(
                "Minimum subset: $min_subset, maximum subset: $max_subset, tolerance: $tolerance",
            )
        end

        p = noise_impt[end]
        verbose && println("Starting importance sampling in noise range $noise_impt")

        # setup noise model
        chn = Dict(0 => 1 - p, 1 => p / 3, 2 => p / 3, 3 => p / 3)
        Pauli_types = collect(keys(chn))
        Pauli_weights = collect(values(chn))
        PrII = Pauli_weights[1] / (Pauli_weights[1] + Pauli_weights[4])
        PrIX = Pauli_weights[2] / (Pauli_weights[1] + Pauli_weights[2])
        PrZI = Pauli_weights[4] / (Pauli_weights[1] + Pauli_weights[4])
        PrZX = Pauli_weights[3] / (Pauli_weights[2] + Pauli_weights[3])
        without_X_err = log(PrII / PrZI)
        with_X_err = log(PrIX / PrZX)

        sample_range = collect(1:n)
        for s in subsets
            if s == 0
                subset_dic[0] = 0.0
            elseif fail_flag
                subset_dic[s] = 1.0
            else
                verbose && println("Subset $s")
                Threads.@threads for th = 1:num_threads
                    X_stabs,
                    Z_stabs,
                    logs,
                    X_err,
                    Z_err,
                    X_var_adj_list,
                    X_check_adj_list,
                    Z_var_adj_list,
                    Z_check_adj_list,
                    X_chn_inits,
                    Z_chn_inits,
                    X_check_to_var_messages,
                    X_var_to_check_messages,
                    Z_check_to_var_messages,
                    Z_var_to_check_messages,
                    current_bits,
                    totals,
                    X_syn,
                    Z_syn = _message_passing_init_CSS(
                        S,
                        chn,
                        max_iter,
                        missing,
                        missing,
                        schedule,
                        erasures,
                    )

                    nr_X = size(X_stabs, 1)
                    nr_Z = size(Z_stabs, 1)
                    nr_logs = size(logs, 1)
                    X_meas_syn = zeros(Int, nr_X)
                    Z_meas_syn = zeros(Int, nr_Z)
                    err_locs = zeros(Int, s)

                    for _ = 1:runs_per_thread
                        # sample
                        X_err[:] .= 0
                        Z_err[:] .= 0
                        sample!(sample_range, err_locs, replace = false)
                        @inbounds for e in err_locs
                            # exclude I here
                            Pauli_op =
                                sample(Pauli_types[2:end], Weights(Pauli_weights[2:end]))
                            # println(Pauli_op)
                            if Pauli_op == 1
                                # X
                                X_err[e] = 1
                                Z_err[e] = 0
                            elseif Pauli_op == 2
                                # Y
                                X_err[e] = 1
                                Z_err[e] = 1
                            else
                                # Z
                                X_err[e] = 0
                                Z_err[e] = 1
                            end
                        end

                        # X syndrome
                        LinearAlgebra.mul!(X_meas_syn, X_stabs, Z_err)
                        @inbounds @simd for i = 1:nr_X
                            X_meas_syn[i] %= 2
                        end

                        # decode X
                        # X_flag, X_out, X_iter, = _message_passing(X_stabs, X_meas_syn,
                        #     X_chn_inits, _SP_check_node_message_box_plus, X_var_adj_list,
                        #     X_check_adj_list, max_iter, schedule, current_bits, totals,
                        #     X_syn, X_check_to_var_messages, X_var_to_check_messages, 0.0)
                        X_flag, X_out, X_iter, = _message_passing_layered(
                            X_stabs,
                            X_meas_syn,
                            X_chn_inits,
                            _SP_check_node_message_box_plus,
                            X_var_adj_list,
                            X_check_adj_list,
                            max_iter,
                            schedule,
                            current_bits,
                            totals,
                            X_syn,
                            X_check_to_var_messages,
                            X_var_to_check_messages,
                            0.0,
                            X_layers,
                        )
                        # shouldn't matter if I mod 2 this now or later
                        Z_err += X_out
                        X_flag && (X_local_iters[th][X_iter] += 1;)
                        # println(X_flag, " ", X_out, " ", X_iter)

                        # Z syndrome
                        LinearAlgebra.mul!(Z_meas_syn, Z_stabs, X_err)
                        @inbounds @simd for i = 1:nr_Z
                            Z_meas_syn[i] %= 2
                        end
                        # println(Z_syn)

                        # update priors on Z given the results from X
                        # always log(Pr(no error) / Pr(error)) so now
                        # log(Pr(I | I) / Pr(Z | I)) and log(Pr(I | X) / Pr(Z | X))
                        @inbounds @simd for i = 1:n
                            if X_out[i] == 0
                                Z_chn_inits[i] = without_X_err
                            else
                                Z_chn_inits[i] = with_X_err
                            end
                        end

                        # decode Z
                        # Z_flag, Z_out, Z_iter, = _message_passing(Z_stabs, Z_meas_syn,
                        #     Z_chn_inits, _SP_check_node_message_box_plus, Z_var_adj_list,
                        #     Z_check_adj_list, max_iter, schedule, current_bits, totals,
                        #     Z_syn, Z_check_to_var_messages, Z_var_to_check_messages, 0.0)
                        Z_flag, Z_out, Z_iter, = _message_passing_layered(
                            Z_stabs,
                            Z_meas_syn,
                            Z_chn_inits,
                            _SP_check_node_message_box_plus,
                            Z_var_adj_list,
                            Z_check_adj_list,
                            max_iter,
                            schedule,
                            current_bits,
                            totals,
                            Z_syn,
                            Z_check_to_var_messages,
                            Z_var_to_check_messages,
                            0.0,
                            Z_layers,
                        )
                        X_err += Z_out
                        Z_flag && (Z_local_iters[th][Z_iter] += 1;)
                        # println(Z_flag, " ", Z_out, " ", Z_iter)

                        if X_flag && Z_flag
                            # converged
                            # this implements !iszero(hcat(Z_err, -X_err) * transpose(logs)) but avoids the
                            # allocations and stops without having the compute the entire product when the
                            # answer is known
                            # temp = vcat(Z_err, -X_err)' * transpose(logs)
                            # if !iszero(temp .% 2)
                            #     local_counts[th] += 1
                            # end
                            # TODO keep track of the logical error rate on individual logical qubits
                            @inbounds for j = 1:nr_logs
                                # introduced a logical error
                                iseven(
                                    dot(view(logs, j, 1:n), Z_err) -
                                    dot(view(logs, j, (n+1):(2*n)), X_err),
                                ) || (local_log_failure_counts[th] += 1; break;)
                            end
                            # continue
                        else
                            # did not converge
                            local_conv_failure_counts[th] += 1
                        end

                        # reset inputs for next run, but don't re-allocate new memory
                        X_check_to_var_messages[:, :, :] .= 0.0
                        X_var_to_check_messages[:, :, :] .= 0.0
                        Z_check_to_var_messages[:, :, :] .= 0.0
                        Z_var_to_check_messages[:, :, :] .= 0.0
                    end
                end

                verbose && println("logical failures: $local_log_failure_counts")
                verbose && println("convergence failures: $local_conv_failure_counts")
                subset_dic[s] =
                    (sum(local_conv_failure_counts) + sum(local_log_failure_counts)) /
                    new_num_runs

                # short circuit the rest of the subsets if they are all going to fail
                if subset_dic[s] ≥ 0.97
                    fail_flag = true
                    verbose && println("Short circuit flag set at subset $s")
                end

                local_conv_failure_counts[:] .= 0
                local_log_failure_counts[:] .= 0
            end
        end

        for (i, p) in enumerate(noise_impt)
            subsets = find_subsets(n, p, tolerance)
            for s in subsets
                FER[i] += s[2] * subset_dic[s[1]]
            end
        end
        verbose && println(FER[1:imp_sam_ind])
    end

    if imp_sam_ind ≠ len_noise
        verbose && println("Starting direct sampling")
        # direct sampling the rest of the points
        if imp_sam_ind ≠ 0
            local_conv_failure_counts[:] .= 0
            local_log_failure_counts[:] .= 0
        end

        for i = (imp_sam_ind+1):len_noise
            p = noise[i]
            println("Starting p = $p")

            # setup noise model
            chn = Dict(0 => 1 - p, 1 => p / 3, 2 => p / 3, 3 => p / 3)
            Pauli_types = collect(keys(chn))
            Pauli_weights = collect(values(chn))
            PrII = Pauli_weights[1] / (Pauli_weights[1] + Pauli_weights[4])
            PrIX = Pauli_weights[2] / (Pauli_weights[1] + Pauli_weights[2])
            PrZI = Pauli_weights[4] / (Pauli_weights[1] + Pauli_weights[4])
            PrZX = Pauli_weights[3] / (Pauli_weights[2] + Pauli_weights[3])
            without_X_err = log(PrII / PrZI)
            with_X_err = log(PrIX / PrZX)

            Threads.@threads for th = 1:num_threads
                X_stabs_dir,
                Z_stabs_dir,
                logs_dir,
                X_err_dir,
                Z_err_dir,
                X_var_adj_list_dir,
                X_check_adj_list_dir,
                Z_var_adj_list_dir,
                Z_check_adj_list_dir,
                X_chn_inits_dir,
                Z_chn_inits_dir,
                X_check_to_var_messages_dir,
                X_var_to_check_messages_dir,
                Z_check_to_var_messages_dir,
                Z_var_to_check_messages_dir,
                current_bits_dir,
                totals_dir,
                X_syn_dir,
                Z_syn_dir = _message_passing_init_CSS(
                    S,
                    chn,
                    max_iter,
                    missing,
                    missing,
                    schedule,
                    erasures,
                )

                nr_X = size(X_stabs_dir, 1)
                nr_Z = size(Z_stabs_dir, 1)
                nr_logs = size(logs_dir, 1)
                X_meas_syn_dir = zeros(Int, nr_X)
                Z_meas_syn_dir = zeros(Int, nr_Z)

                for _ = 1:runs_per_thread
                    # sample
                    @inbounds @simd for j = 1:n
                        Pauli_op = sample(Pauli_types, Weights(Pauli_weights))
                        # println(Pauli_op)
                        if Pauli_op == 0
                            # I
                            X_err_dir[j] = 0
                            Z_err_dir[j] = 0
                        elseif Pauli_op == 1
                            # X
                            X_err_dir[j] = 1
                            Z_err_dir[j] = 0
                        elseif Pauli_op == 2
                            # Y
                            X_err_dir[j] = 1
                            Z_err_dir[j] = 1
                        else
                            # Z
                            X_err_dir[j] = 0
                            Z_err_dir[j] = 1
                        end
                    end
                    # println(sum(X_err_dir))
                    # println(sum(Z_err_dir))

                    # X syndrome
                    LinearAlgebra.mul!(X_meas_syn_dir, X_stabs_dir, Z_err_dir)
                    @inbounds @simd for i = 1:nr_X
                        X_meas_syn_dir[i] %= 2
                    end
                    # println(X_syn_dir)

                    # decode X
                    # X_flag_dir, X_out_dir, X_iter_dir, = _message_passing(X_stabs_dir, X_meas_syn_dir,
                    #     X_chn_inits_dir, _SP_check_node_message_box_plus, X_var_adj_list_dir,
                    #     X_check_adj_list_dir, max_iter, schedule, current_bits_dir, totals_dir,
                    #     X_syn_dir, X_check_to_var_messages_dir, X_var_to_check_messages_dir, 0.0)
                    X_flag_dir, X_out_dir, X_iter_dir, = _message_passing_layered(
                        X_stabs_dir,
                        X_meas_syn_dir,
                        X_chn_inits_dir,
                        _SP_check_node_message_box_plus,
                        X_var_adj_list_dir,
                        X_check_adj_list_dir,
                        max_iter,
                        schedule,
                        current_bits_dir,
                        totals_dir,
                        X_syn_dir,
                        X_check_to_var_messages_dir,
                        X_var_to_check_messages_dir,
                        0.0,
                        X_layers,
                    )
                    # shouldn't matter if I mod 2 this now or later?
                    Z_err_dir += X_out_dir
                    # X_flag_dir && push!(X_local_iters[th], X_iter_dir)
                    X_flag_dir && (X_local_iters[th][X_iter_dir] += 1;)
                    # println(X_flag_dir, " ", X_out_dir, " ", X_iter_dir)

                    # Z syndrome
                    LinearAlgebra.mul!(Z_meas_syn_dir, Z_stabs_dir, X_err_dir)
                    @inbounds @simd for i = 1:nr_Z
                        Z_meas_syn_dir[i] %= 2
                    end
                    # println(Z_syn_dir)

                    # update priors on Z given the results from X
                    # always log(Pr(no error) / Pr(error)) so now
                    # log((Pr(I | I) + Pr(I | X)) / (Pr(Z | I) + Pr(Z | X)))
                    @inbounds @simd for i = 1:n
                        if X_out_dir[i] == 0
                            Z_chn_inits_dir[i] = without_X_err
                        else
                            Z_chn_inits_dir[i] = with_X_err
                        end
                    end
                    # decode Z
                    # Z_flag_dir, Z_out_dir, Z_iter_dir, = _message_passing(Z_stabs_dir, Z_meas_syn_dir,
                    #     Z_chn_inits_dir, _SP_check_node_message_box_plus, Z_var_adj_list_dir,
                    #     Z_check_adj_list_dir, max_iter, schedule, current_bits_dir, totals_dir,
                    #     Z_syn_dir, Z_check_to_var_messages_dir, Z_var_to_check_messages_dir, 0.0)
                    Z_flag_dir, Z_out_dir, Z_iter_dir, = _message_passing_layered(
                        Z_stabs_dir,
                        Z_meas_syn_dir,
                        Z_chn_inits_dir,
                        _SP_check_node_message_box_plus,
                        Z_var_adj_list_dir,
                        Z_check_adj_list_dir,
                        max_iter,
                        schedule,
                        current_bits_dir,
                        totals_dir,
                        Z_syn_dir,
                        Z_check_to_var_messages_dir,
                        Z_var_to_check_messages_dir,
                        0.0,
                        Z_layers,
                    )
                    X_err_dir += Z_out_dir
                    # Z_flag_dir && push!(Z_local_iters[th], Z_iter_dir)
                    Z_flag_dir && (Z_local_iters[th][Z_iter_dir] += 1;)
                    # println(Z_flag_dir, " ", Z_out_dir, " ", Z_iter_dir)
                    # return

                    if X_flag_dir && Z_flag_dir
                        # converged
                        # this implements !iszero(hcat(Z_err, -X_err) * transpose(logs)) but avoids the
                        # allocations and stops without having the compute the entire product when the
                        # answer is known
                        # temp = vcat(Z_err_dir, -X_err_dir)' * transpose(logs_dir)
                        # if !iszero(temp .% 2)
                        #     local_counts[th] += 1
                        # end
                        # TODO keep track of the logical error rate on individual logical qubits
                        @inbounds for j = 1:nr_logs
                            # introduced a logical error
                            iseven(
                                dot(view(logs_dir, j, 1:n), Z_err_dir) -
                                dot(view(logs_dir, j, (n+1):(2*n)), X_err_dir),
                            ) || (local_log_failure_counts[th] += 1; break;)
                        end
                        # continue
                    else
                        # did not converge
                        local_conv_failure_counts[th] += 1
                    end

                    # reset inputs for next run, but don't re-allocate new memory
                    X_check_to_var_messages_dir[:, :, :] .= 0.0
                    X_var_to_check_messages_dir[:, :, :] .= 0.0
                    Z_check_to_var_messages_dir[:, :, :] .= 0.0
                    Z_var_to_check_messages_dir[:, :, :] .= 0.0
                end
            end
            verbose && println("logical failures: $local_log_failure_counts")
            verbose && println("convergence failures: $local_conv_failure_counts")
            @inbounds FER[i] =
                (sum(local_conv_failure_counts) + sum(local_log_failure_counts)) /
                new_num_runs
            verbose && println("FER = $(FER[i])")
            verbose && println("Finished p = $p")

            local_conv_failure_counts[:] .= 0
            local_log_failure_counts[:] .= 0
        end
    end
    # return FER, StatsBase.countmap(reduce(vcat, X_local_iters)),
    # StatsBase.countmap(reduce(vcat, Z_local_iters))
    return FER,
    foldl(mergewith!(+), X_local_iters; init = Dict{Int,Int}()),
    foldl(mergewith!(+), Z_local_iters; init = Dict{Int,Int}())
end
CSS_decoder_with_Bayes(::IsNotCSS, S::AbstractStabilizerCode, verbose::Bool) =
    throw(ArgumentError("CSS decoders are only valid for CSS codes"))

####
# M

CSS_metachecks_Mike(S::T; verbose::Bool = true) where {T<:AbstractStabilizerCode} =
    CSS_metachecks_Mike(CSSTrait(T), S, verbose)
function CSS_metachecks_Mike(::IsCSS, S::AbstractStabilizerCode, verbose::Bool)
    # noise details
    # these are the markers from the paper
    # noise = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14]
    noise = 10 .^ (-3:0.25:-1);
    len_noise = length(noise)

    # BP details
    n = length(S)
    num_runs = 20000
    max_iter = 32
    # schedule = :parallel
    # schedule = :serial
    schedule = :layered
    X_layers = layered_schedule(S.X_stabs)
    Z_layers = layered_schedule(S.Z_stabs)
    X_layers_meta = layered_schedule(X_metacheck(S))
    Z_layers_meta = layered_schedule(Z_metacheck(S))
    erasures = Int[]
    X_syn_erasures = Int[]
    Z_syn_erasures = Int[]
    tolerance = 1e-12

    # multi-threading details
    num_threads = Threads.nthreads()
    runs_per_thread = cld(num_runs, num_threads)
    new_num_runs = runs_per_thread * num_threads
    verbose && println(
        "Number of threads: $num_threads, runs per thread: $runs_per_thread, new number of runs: $new_num_runs",
    )

    # preallocate places to store results
    FER = zeros(Float64, len_noise)
    local_log_failure_counts = zeros(Int, num_threads)
    local_conv_failure_counts = zeros(Int, num_threads)
    X_local_iters = [Dict{Int,Int}() for _ = 1:num_threads]
    Z_local_iters = [Dict{Int,Int}() for _ = 1:num_threads]
    for i = 1:num_threads
        sizehint!(X_local_iters[i], max_iter)
        sizehint!(Z_local_iters[i], max_iter)
        for j = 1:max_iter
            X_local_iters[i][j] = 0
            Z_local_iters[i][j] = 0
        end
    end

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

    # imp_sam_ind = 2
    if imp_sam_ind ≠ 0
        noise_impt = noise[1:imp_sam_ind]
        min_subset, max_subset = min_max_subsets(n, noise, tolerance)
        # TODO: fix
        min_subset = 0
        max_subset = 36
        subsets = [i for i = min_subset:max_subset]
        subset_dic = Dict{Int,Float64}()
        fail_flag = false

        if verbose
            println("Starting importance sampling on $(length(noise_impt)) noise values.")
            println(
                "Minimum subset: $min_subset, maximum subset: $max_subset, tolerance: $tolerance",
            )
        end

        p = noise_impt[end]
        verbose && println("Starting importance sampling in noise range $noise_impt")

        # setup noise model
        chn = Dict(0 => 1 - p, 1 => p / 3, 2 => p / 3, 3 => p / 3)
        Pauli_types = collect(keys(chn))
        Pauli_weights = collect(values(chn))
        PrII = Pauli_weights[1] / (Pauli_weights[1] + Pauli_weights[4])
        PrIX = Pauli_weights[2] / (Pauli_weights[1] + Pauli_weights[2])
        PrZI = Pauli_weights[4] / (Pauli_weights[1] + Pauli_weights[4])
        PrZX = Pauli_weights[3] / (Pauli_weights[2] + Pauli_weights[3])
        without_X_err = log(PrII / PrZI)
        with_X_err = log(PrIX / PrZX)
        p_meas = 1e-3
        chn_X_meas = BSC(p_meas)
        chn_Z_meas = BSC(p_meas)
        syn_err_choices = [0, 1]
        syn_err_weights = Weights([1 - p_meas, p_meas])

        sample_range = collect(1:n)
        for s in subsets
            if s == 0
                subset_dic[0] = 0.0
            elseif fail_flag
                subset_dic[s] = 1.0
            else
                verbose && println("Subset $s")
                Threads.@threads for th = 1:num_threads
                    # code
                    X_stabs,
                    Z_stabs,
                    logs,
                    X_err,
                    Z_err,
                    X_var_adj_list,
                    X_check_adj_list,
                    Z_var_adj_list,
                    Z_check_adj_list,
                    X_chn_inits,
                    Z_chn_inits,
                    X_check_to_var_messages,
                    X_var_to_check_messages,
                    Z_check_to_var_messages,
                    Z_var_to_check_messages,
                    current_bits,
                    totals,
                    X_syn,
                    Z_syn = _message_passing_init_CSS(
                        S,
                        chn,
                        max_iter,
                        missing,
                        missing,
                        schedule,
                        erasures,
                    )

                    nr_X = size(X_stabs, 1)
                    nr_Z = size(Z_stabs, 1)
                    nr_logs = size(logs, 1)
                    X_meas_syn = zeros(Int, nr_X)
                    X_syn_err = zeros(Int, nr_X)
                    Z_meas_syn = zeros(Int, nr_Z)
                    Z_syn_err = zeros(Int, nr_Z)
                    err_locs = zeros(Int, s)

                    # metachecks
                    X_meta_Int,
                    _,
                    var_adj_list_X_meta,
                    check_adj_list_X_meta,
                    chn_inits_X_meta,
                    check_to_var_messages_X_meta,
                    var_to_check_messages_X_meta,
                    current_bits_X_meta,
                    totals_X_meta,
                    syn_X_meta = _message_passing_init(
                        X_metacheck(S),
                        X_meas_syn,
                        chn_X_meas,
                        max_iter,
                        :SP,
                        missing,
                        schedule,
                        X_syn_erasures,
                    )
                    Z_meta_Int,
                    _,
                    var_adj_list_Z_meta,
                    check_adj_list_Z_meta,
                    chn_inits_Z_meta,
                    check_to_var_messages_Z_meta,
                    var_to_check_messages_Z_meta,
                    current_bits_Z_meta,
                    totals_Z_meta,
                    syn_Z_meta = _message_passing_init(
                        Z_metacheck(S),
                        Z_meas_syn,
                        chn_Z_meas,
                        max_iter,
                        :SP,
                        missing,
                        schedule,
                        Z_syn_erasures,
                    )

                    for _ = 1:runs_per_thread
                        # sample
                        X_err[:] .= 0
                        Z_err[:] .= 0
                        sample!(sample_range, err_locs, replace = false)
                        @inbounds for e in err_locs
                            # exclude I here
                            Pauli_op =
                                sample(Pauli_types[2:end], Weights(Pauli_weights[2:end]))
                            # println(Pauli_op)
                            if Pauli_op == 1
                                # X
                                X_err[e] = 1
                                Z_err[e] = 0
                            elseif Pauli_op == 2
                                # Y
                                X_err[e] = 1
                                Z_err[e] = 1
                            else
                                # Z
                                X_err[e] = 0
                                Z_err[e] = 1
                            end
                        end

                        # X syndrome
                        LinearAlgebra.mul!(X_meas_syn, X_stabs, Z_err)
                        @inbounds @simd for i = 1:nr_X
                            X_meas_syn[i] %= 2
                        end

                        # add measurement errors
                        sample!(syn_err_choices, syn_err_weights, X_syn_err)
                        if !iszero(X_syn_err)
                            X_meas_syn += X_syn_err
                            X_meas_syn .%= 2

                            # we are justified putting this inside because it will always come out 
                            # zero otherwise
                            # TODO: this multiplication
                            if !iszero(X_metacheck * X_meas_syn)
                                # decode X syndrome
                                X_flag_meta, X_out_meta, X_iter_meta, =
                                    _message_passing_layered(
                                        X_meta_Int,
                                        X_meas_syn,
                                        chn_inits_X_meta,
                                        _SP_check_node_message_box_plus,
                                        var_adj_list_X_meta,
                                        check_adj_list_X_meta,
                                        max_iter,
                                        schedule,
                                        current_bits_X_meta,
                                        totals_X_meta,
                                        syn_X_meta,
                                        check_to_var_messages_X_meta,
                                        var_to_check_messages_X_meta,
                                        0.0,
                                        X_layers_meta,
                                    )

                                if X_flag_meta
                                    # converged
                                    # check if valid syndrome using element of homology
                                    # TODO this multiplication
                                    if !iszero(L_X * X_meas_syn)
                                        # detected an invalid syndrome, re-decode with L_X
                                        # TODO this matrix and these values
                                        # TODO unclear some of these need an L version
                                        X_flag_meta, X_out_meta, X_iter_meta, =
                                            _message_passing_layered(
                                                X_meta_L_Int,
                                                X_meas_syn,
                                                chn_inits_X_meta_L,
                                                _SP_check_node_message_box_plus,
                                                var_adj_list_X_meta_L,
                                                check_adj_list_X_meta_L,
                                                max_iter,
                                                schedule,
                                                current_bits_X_meta_L,
                                                totals_X_meta_L,
                                                syn_X_meta_L,
                                                check_to_var_messages_X_meta_L,
                                                var_to_check_messages_X_meta_L,
                                                0.0,
                                                X_layers_meta_L,
                                            )

                                        if !X_flag_meta
                                            # failed to converge
                                            # TODO set
                                            skip = true
                                        end
                                    end

                                    if !skip
                                        X_meas_syn += X_out_meta
                                        X_meas_syn .%= 2
                                    end
                                else
                                    # TODO should automatically retry with L?
                                    # failed to converge
                                    skip = true
                                end
                            end
                        end

                        # if syndrome is known to be bad, declare it a failure since we will be 
                        # forcing it to find an error which doesn't return to the codespace
                        if !skip
                            # decode X
                            X_flag, X_out, X_iter, = _message_passing_layered(
                                X_stabs,
                                X_meas_syn,
                                X_chn_inits,
                                _SP_check_node_message_box_plus,
                                X_var_adj_list,
                                X_check_adj_list,
                                max_iter,
                                schedule,
                                current_bits,
                                totals,
                                X_syn,
                                X_check_to_var_messages,
                                X_var_to_check_messages,
                                0.0,
                                X_layers,
                            )
                            # shouldn't matter if I mod 2 this now or later
                            Z_err += X_out
                            X_flag && (X_local_iters[th][X_iter] += 1; skip = true;)
                        end

                        if !skip
                            # Z syndrome
                            LinearAlgebra.mul!(Z_meas_syn, Z_stabs, X_err)
                            @inbounds @simd for i = 1:nr_Z
                                Z_meas_syn[i] %= 2
                            end

                            # add measurement errors
                            sample!(syn_err_choices, syn_err_weights, Z_syn_err)
                            if !iszero(Z_syn_err)
                                Z_meas_syn += Z_syn_err
                                Z_meas_syn .%= 2

                                # we are justified putting this inside because it will always come 
                                # out zero otherwise
                                # TODO: this multiplication
                                if !iszero(Z_metacheck * Z_meas_syn)
                                    # decode Z syndrome
                                    Z_flag_meta, Z_out_meta, Z_iter_meta, =
                                        _message_passing_layered(
                                            Z_meta_Int,
                                            Z_meas_syn,
                                            chn_inits_Z_meta,
                                            _SP_check_node_message_box_plus,
                                            var_adj_list_Z_meta,
                                            check_adj_list_Z_meta,
                                            max_iter,
                                            schedule,
                                            current_bits_Z_meta,
                                            totals_Z_meta,
                                            syn_Z_meta,
                                            check_to_var_messages_Z_meta,
                                            var_to_check_messages_Z_meta,
                                            0.0,
                                            Z_layers_meta,
                                        )

                                    if Z_flag_meta
                                        # converged
                                        # check if valid syndrome using element of homology
                                        # TODO this multiplication
                                        if !iszero(L_Z * Z_meas_syn)
                                            # detected an invalid syndrome, re-decode with L_X
                                            # TODO this matrix and these values
                                            # TODO unclear some of these need an L version
                                            Z_flag_meta, Z_out_meta, Z_iter_meta, =
                                                _message_passing_layered(
                                                    Z_meta_L_Int,
                                                    Z_meas_syn,
                                                    chn_inits_Z_meta_L,
                                                    _SP_check_node_message_box_plus,
                                                    var_adj_list_Z_meta_L,
                                                    check_adj_list_Z_meta_L,
                                                    max_iter,
                                                    schedule,
                                                    current_bits_Z_meta_L,
                                                    totals_Z_meta_L,
                                                    syn_Z_meta_L,
                                                    check_to_var_messages_Z_meta_L,
                                                    var_to_check_messages_Z_meta_L,
                                                    0.0,
                                                    Z_layers_meta_L,
                                                )

                                            if !X_flag_meta
                                                # failed to converge
                                                # TODO set
                                                skip = true
                                            end
                                        end

                                        if !skip
                                            X_meas_syn += X_out_meta
                                            X_meas_syn .%= 2
                                        end
                                    else
                                        # TODO should automatically retry with L?
                                        # failed to converge
                                        skip = true
                                    end
                                end
                            end

                            # if syndrome is known to be bad, declare it a failure since we will be 
                            # forcing it to find an error which doesn't return to the codespace
                            if !skip
                                # update priors on Z given the results from X
                                # always log(Pr(no error) / Pr(error)) so now
                                # log(Pr(I | I) / Pr(Z | I)) and log(Pr(I | X) / Pr(Z | X))
                                @inbounds @simd for i = 1:n
                                    if X_out[i] == 0
                                        Z_chn_inits[i] = without_X_err
                                    else
                                        Z_chn_inits[i] = with_X_err
                                    end
                                end

                                # decode Z
                                Z_flag, Z_out, Z_iter, = _message_passing_layered(
                                    Z_stabs,
                                    Z_meas_syn,
                                    Z_chn_inits,
                                    _SP_check_node_message_box_plus,
                                    Z_var_adj_list,
                                    Z_check_adj_list,
                                    max_iter,
                                    schedule,
                                    current_bits,
                                    totals,
                                    Z_syn,
                                    Z_check_to_var_messages,
                                    Z_var_to_check_messages,
                                    0.0,
                                    Z_layers,
                                )
                                X_err += Z_out
                                Z_flag && (Z_local_iters[th][Z_iter] += 1; skip = true)
                            end
                        end

                        if !skip
                            # converged
                            # this implements !iszero(hcat(Z_err, -X_err) * transpose(logs)) but avoids the
                            # allocations and stops without having the compute the entire product when the
                            # answer is known
                            # temp = vcat(Z_err, -X_err)' * transpose(logs)
                            # if !iszero(temp .% 2)
                            #     local_counts[th] += 1
                            # end
                            # TODO keep track of the logical error rate on individual logical qubits
                            @inbounds for j = 1:nr_logs
                                # introduced a logical error
                                iseven(
                                    dot(view(logs, j, 1:n), Z_err) -
                                    dot(view(logs, j, (n+1):(2*n)), X_err),
                                ) || (local_log_failure_counts[th] += 1; break;)
                            end
                            # continue
                        else
                            # did not converge
                            local_conv_failure_counts[th] += 1
                        end

                        # reset inputs for next run, but don't re-allocate new memory
                        X_check_to_var_messages[:, :, :] .= 0.0
                        X_var_to_check_messages[:, :, :] .= 0.0
                        Z_check_to_var_messages[:, :, :] .= 0.0
                        Z_var_to_check_messages[:, :, :] .= 0.0
                    end
                end

                verbose && println("logical failures: $local_log_failure_counts")
                verbose && println("convergence failures: $local_conv_failure_counts")
                subset_dic[s] =
                    (sum(local_conv_failure_counts) + sum(local_log_failure_counts)) /
                    new_num_runs

                # short circuit the rest of the subsets if they are all going to fail
                if subset_dic[s] ≥ 0.97
                    fail_flag = true
                    verbose && println("Short circuit flag set at subset $s")
                end

                local_conv_failure_counts[:] .= 0
                local_log_failure_counts[:] .= 0
            end
        end

        for (i, p) in enumerate(noise_impt)
            subsets = find_subsets(n, p, tolerance)
            for s in subsets
                FER[i] += s[2] * subset_dic[s[1]]
            end
        end
        verbose && println(FER[1:imp_sam_ind])
    end

    if imp_sam_ind ≠ len_noise
        verbose && println("Starting direct sampling")
        # direct sampling the rest of the points
        if imp_sam_ind ≠ 0
            local_conv_failure_counts[:] .= 0
            local_log_failure_counts[:] .= 0
        end

        for i = (imp_sam_ind+1):len_noise
            p = noise[i]
            println("Starting p = $p")

            # setup noise model
            chn = Dict(0 => 1 - p, 1 => p / 3, 2 => p / 3, 3 => p / 3)
            Pauli_types = collect(keys(chn))
            Pauli_weights = collect(values(chn))
            PrII = Pauli_weights[1] / (Pauli_weights[1] + Pauli_weights[4])
            PrIX = Pauli_weights[2] / (Pauli_weights[1] + Pauli_weights[2])
            PrZI = Pauli_weights[4] / (Pauli_weights[1] + Pauli_weights[4])
            PrZX = Pauli_weights[3] / (Pauli_weights[2] + Pauli_weights[3])
            without_X_err = log(PrII / PrZI)
            with_X_err = log(PrIX / PrZX)

            Threads.@threads for th = 1:num_threads
                X_stabs_dir,
                Z_stabs_dir,
                logs_dir,
                X_err_dir,
                Z_err_dir,
                X_var_adj_list_dir,
                X_check_adj_list_dir,
                Z_var_adj_list_dir,
                Z_check_adj_list_dir,
                X_chn_inits_dir,
                Z_chn_inits_dir,
                X_check_to_var_messages_dir,
                X_var_to_check_messages_dir,
                Z_check_to_var_messages_dir,
                Z_var_to_check_messages_dir,
                current_bits_dir,
                totals_dir,
                X_syn_dir,
                Z_syn_dir = _message_passing_init_CSS(
                    S,
                    chn,
                    max_iter,
                    missing,
                    missing,
                    schedule,
                    erasures,
                )

                nr_X = size(X_stabs_dir, 1)
                nr_Z = size(Z_stabs_dir, 1)
                nr_logs = size(logs_dir, 1)
                X_meas_syn_dir = zeros(Int, nr_X)
                Z_meas_syn_dir = zeros(Int, nr_Z)

                for _ = 1:runs_per_thread
                    # sample
                    @inbounds @simd for j = 1:n
                        Pauli_op = sample(Pauli_types, Weights(Pauli_weights))
                        # println(Pauli_op)
                        if Pauli_op == 0
                            # I
                            X_err_dir[j] = 0
                            Z_err_dir[j] = 0
                        elseif Pauli_op == 1
                            # X
                            X_err_dir[j] = 1
                            Z_err_dir[j] = 0
                        elseif Pauli_op == 2
                            # Y
                            X_err_dir[j] = 1
                            Z_err_dir[j] = 1
                        else
                            # Z
                            X_err_dir[j] = 0
                            Z_err_dir[j] = 1
                        end
                    end

                    # X syndrome
                    LinearAlgebra.mul!(X_meas_syn_dir, X_stabs_dir, Z_err_dir)
                    @inbounds @simd for i = 1:nr_X
                        X_meas_syn_dir[i] %= 2
                    end

                    # decode X
                    X_flag_dir, X_out_dir, X_iter_dir, = _message_passing_layered(
                        X_stabs_dir,
                        X_meas_syn_dir,
                        X_chn_inits_dir,
                        _SP_check_node_message_box_plus,
                        X_var_adj_list_dir,
                        X_check_adj_list_dir,
                        max_iter,
                        schedule,
                        current_bits_dir,
                        totals_dir,
                        X_syn_dir,
                        X_check_to_var_messages_dir,
                        X_var_to_check_messages_dir,
                        0.0,
                        X_layers,
                    )
                    # shouldn't matter if I mod 2 this now or later?
                    Z_err_dir += X_out_dir
                    X_flag_dir && (X_local_iters[th][X_iter_dir] += 1;)

                    # Z syndrome
                    LinearAlgebra.mul!(Z_meas_syn_dir, Z_stabs_dir, X_err_dir)
                    @inbounds @simd for i = 1:nr_Z
                        Z_meas_syn_dir[i] %= 2
                    end

                    # update priors on Z given the results from X
                    # always log(Pr(no error) / Pr(error)) so now
                    # log((Pr(I | I) + Pr(I | X)) / (Pr(Z | I) + Pr(Z | X)))   
                    @inbounds @simd for i = 1:n
                        if X_out_dir[i] == 0
                            Z_chn_inits_dir[i] = without_X_err
                        else
                            Z_chn_inits_dir[i] = with_X_err
                        end
                    end
                    # decode Z
                    Z_flag_dir, Z_out_dir, Z_iter_dir, = _message_passing_layered(
                        Z_stabs_dir,
                        Z_meas_syn_dir,
                        Z_chn_inits_dir,
                        _SP_check_node_message_box_plus,
                        Z_var_adj_list_dir,
                        Z_check_adj_list_dir,
                        max_iter,
                        schedule,
                        current_bits_dir,
                        totals_dir,
                        Z_syn_dir,
                        Z_check_to_var_messages_dir,
                        Z_var_to_check_messages_dir,
                        0.0,
                        Z_layers,
                    )
                    X_err_dir += Z_out_dir
                    Z_flag_dir && (Z_local_iters[th][Z_iter_dir] += 1;)

                    if X_flag_dir && Z_flag_dir
                        # converged
                        # this implements !iszero(hcat(Z_err, -X_err) * transpose(logs)) but avoids the
                        # allocations and stops without having the compute the entire product when the
                        # answer is known
                        # temp = vcat(Z_err_dir, -X_err_dir)' * transpose(logs_dir)
                        # if !iszero(temp .% 2)
                        #     local_counts[th] += 1
                        # end
                        # TODO keep track of the logical error rate on individual logical qubits
                        @inbounds for j = 1:nr_logs
                            # introduced a logical error
                            iseven(
                                dot(view(logs_dir, j, 1:n), Z_err_dir) -
                                dot(view(logs_dir, j, (n+1):(2*n)), X_err_dir),
                            ) || (local_log_failure_counts[th] += 1; break;)
                        end
                        # continue
                    else
                        # did not converge
                        local_conv_failure_counts[th] += 1
                    end

                    # reset inputs for next run, but don't re-allocate new memory
                    X_check_to_var_messages_dir[:, :, :] .= 0.0
                    X_var_to_check_messages_dir[:, :, :] .= 0.0
                    Z_check_to_var_messages_dir[:, :, :] .= 0.0
                    Z_var_to_check_messages_dir[:, :, :] .= 0.0
                end
            end
            verbose && println("logical failures: $local_log_failure_counts")
            verbose && println("convergence failures: $local_conv_failure_counts")
            @inbounds FER[i] =
                (sum(local_conv_failure_counts) + sum(local_log_failure_counts)) /
                new_num_runs
            verbose && println("FER = $(FER[i])")
            verbose && println("Finished p = $p")

            local_conv_failure_counts[:] .= 0
            local_log_failure_counts[:] .= 0
        end
    end
    return FER,
    foldl(mergewith!(+), X_local_iters; init = Dict{Int,Int}()),
    foldl(mergewith!(+), Z_local_iters; init = Dict{Int,Int}())
end

CSS_metachecks_Mike(::IsNotCSS, S::AbstractStabilizerCode, verbose::Bool) =
    throw(ArgumentError("CSS decoders are only valid for CSS codes"))

function _make_single_shot_Tanner_graph(H::T, M::T) where {T<:CTMatrixTypes}
    # internal function, so skip checks on correctness
    nc_H = ncols(H)
    nc_M = ncols(M)
    single_shot_matrix = _Flint_matrix_to_Julia_int_matrix(H ⊕ M)
    @inbounds @simd for r = (nc_H+1):(nc_H+nc_M)
        single_shot_matrix[r, r] = 1
    end
    # of form (H I; 0 M)
    return single_shot_matrix
end
