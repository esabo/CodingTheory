function CSS_decoder_simple(S::AbstractStabilizerCode)
    # BP details
    n = length(S)
    num_runs = 1000
    max_iter = 20
    schedule = :layered
    X_layers::Vector{Vector{Int}} = layered_schedule(X_stabilizers(S))
    Z_layers::Vector{Vector{Int}} = layered_schedule(Z_stabilizers(S))
    erasures = Int[]

    # noise details
    noise = 0.1
    # noise = 10 .^(-2:0.1:-1);
    len_noise = length(noise)

    # multi-threading details
    num_threads = Threads.nthreads()
    runs_per_thread = cld(num_runs, num_threads)
    new_num_runs = runs_per_thread * num_threads
    println("Number of threads: $num_threads, runs per thread: $runs_per_thread, new number of runs: $new_num_runs")

    # preallocate places to store results
    FER = zeros(Float64, len_noise)
    local_log_failure_counts = zeros(Int, num_threads)
    local_conv_failure_counts = zeros(Int, num_threads)
    X_local_iters = [zeros(Int, max_iter) for _ in 1:num_threads]
    Z_local_iters = [zeros(Int, max_iter) for _ in 1:num_threads]

    for i in 1:len_noise
        p = noise[i]
        println("Starting p = $p")

        # setup noise model
        chn = Dict(0 => 1 - p, 1 => p / 3, 2 => p / 3, 3 => p / 3)
        Pauli_types = Int[0, 1, 2, 3]
        Pauli_weights = Float64[1 - p, p / 3, p / 3, p / 3]
        PrII::Float64 = Pauli_weights[1] / (Pauli_weights[1] + Pauli_weights[4])
        PrIX::Float64 = Pauli_weights[2] / (Pauli_weights[1] + Pauli_weights[2])
        PrZI::Float64 = Pauli_weights[4] / (Pauli_weights[1] + Pauli_weights[4])
        PrZX::Float64 = Pauli_weights[3] / (Pauli_weights[2] + Pauli_weights[3])
        without_X_err::Float64 = log(PrII / PrZI)
        with_X_err::Float64 = log(PrIX / PrZX)
        Pauli_weights_W = Weights(Pauli_weights)
        Pauli_op::Vector{Int} = [0]

        Threads.@threads for th in 1:num_threads
            X_stabs, Z_stabs, logs, X_err, Z_err, X_var_adj_list,
                X_check_adj_list, Z_var_adj_list, Z_check_adj_list, X_chn_inits,
                Z_chn_inits, X_check_to_var_messages, X_var_to_check_messages,
                Z_check_to_var_messages, Z_var_to_check_messages, current_bits,
                totals, X_syn, Z_syn = CodingTheory._message_passing_init_CSS(S, chn, max_iter,
                missing, missing, schedule, erasures)

            nr_X = size(X_stabs, 1)
            nr_Z = size(Z_stabs, 1)
            nr_logs = size(logs, 1)
            X_meas_syn = zeros(Int, nr_X)
            Z_meas_syn = zeros(Int, nr_Z)
            skip = false

            for _ in 1:runs_per_thread
                # sample
                @inbounds @simd for j in 1:n
                    sample!(Pauli_types, Pauli_weights_W, Pauli_op)
                    # println(Pauli_op)
                    if Pauli_op[1] == 0
                        # I
                        X_err[j] = 0
                        Z_err[j] = 0
                    elseif Pauli_op[1] == 1
                        # X
                        X_err[j] = 1
                        Z_err[j] = 0
                    elseif Pauli_op[1] == 2
                        # Y
                        X_err[j] = 1
                        Z_err[j] = 1
                    else
                        # Z
                        X_err[j] = 0
                        Z_err[j] = 1
                    end
                end

                # bypass = false
                # X syndrome
                LinearAlgebra.mul!(X_meas_syn, X_stabs, Z_err)
                @inbounds @simd for i in 1:nr_X
                    X_meas_syn[i] %= 2
                end

                # decode X
                X_flag, X_out, X_iter = CodingTheory._message_passing_layered(X_stabs, X_meas_syn,
                    X_chn_inits, CodingTheory._SP_check_node_message_box_plus, X_var_adj_list,
                    X_check_adj_list, max_iter, schedule, current_bits, totals,
                    X_syn, X_check_to_var_messages, X_var_to_check_messages, 0.0, X_layers)

                
                if !X_flag
                    # did not converge
                    skip = true
                else
                    # did converge
                    # shouldn't matter if I mod 2 this now or later
                    # a non-allocating version of Z_err += X_out
                    # Z_err += X_out
                    @inbounds @simd for i in 1:n
                        Z_err[i] = Z_err[i] + X_out[i]
                    end
                    X_local_iters[th][X_iter] += 1
                end

                # skip Z if X failed because the result will not end up in the codespace
                # if !skip
                    # Z syndrome
                    # this allocates ~400 MB for 15k samples
                    LinearAlgebra.mul!(Z_meas_syn, Z_stabs, X_err)
                    @inbounds @simd for i in 1:nr_Z
                        Z_meas_syn[i] %= 2
                    end
                
                    # Bayes' Theorem: update priors on Z given the results from X
                    # always log(Pr(no error) / Pr(error)) so now
                    # log((Pr(I | I) + Pr(I | X)) / (Pr(Z | I) + Pr(Z | X))) etc
                
                    # this single loop allocates 0.6 GB for 15k samples
                    @inbounds @simd for i in 1:n
                        if X_out[i] == 0
                            Z_chn_inits[i] = without_X_err
                        else
                            Z_chn_inits[i] = with_X_err
                        end
                    end
                    
                    # if bypass
                    Z_flag, Z_out, Z_iter = CodingTheory._message_passing_layered(Z_stabs, Z_meas_syn,
                        Z_chn_inits, CodingTheory._SP_check_node_message_box_plus, Z_var_adj_list,
                        Z_check_adj_list, max_iter, schedule, current_bits, totals,
                        Z_syn, Z_check_to_var_messages, Z_var_to_check_messages, 0.0, Z_layers)

                    if !Z_flag
                        # did not converge
                        skip = true
                    else
                        # did converge
                        # shouldn't matter if I mod 2 this now or later
                        # a non-allocating version of X_err += Z_out
                        @inbounds @simd for i in 1:n
                            X_err[i] = X_err[i] + Z_out[i]
                        end
                        Z_local_iters[th][Z_iter] += 1
                    end
                # end

                if !skip
                    # converged, check for logical errors
                    @inbounds for j in 1:nr_logs
                        iseven(dot(view(logs, j, 1:n), Z_err) - dot(view(logs, j, n +
                            1:2 * n), X_err)) || (local_log_failure_counts[th] += 1; break;)
                    end
                else
                    # did not converge
                    local_conv_failure_counts[th] += 1
                end

                # # reset inputs for next run, but don't re-allocate new memory
                X_check_to_var_messages .= 0.0
                X_var_to_check_messages .= 0.0
                Z_check_to_var_messages .= 0.0
                Z_var_to_check_messages .= 0.0
                skip = false
            end
        end
        println("logical failures: $local_log_failure_counts")
        println("convergence failures: $local_conv_failure_counts")
        @inbounds FER[i] = (sum(local_conv_failure_counts) + sum(local_log_failure_counts)) / new_num_runs
        println("FER = $(FER[i])")
        println("Finished p = $p")

        local_conv_failure_counts[:] .= 0
        local_log_failure_counts[:] .= 0
    end

    # reduce iteration counts
    @inbounds for i in 2:num_threads
        @simd for j in 1:max_iter
            X_local_iters[1][j] += X_local_iters[i][j]
            Z_local_iters[1][j] += Z_local_iters[i][j]
        end
    end

    return FER, X_local_iters[1], Z_local_iters[1]
end

# 43.825219 seconds (39.97 M allocations: 1.635 GiB, 0.34% gc time)

# Data processing
# noise = 10 .^(-2:0.1:-1);
# fig = Figure()
# ax = Axis(fig[1, 1],
#     limits = ((noise[1], noise[end]), (1e-5, 1)),
#     xscale = log10,
#     yscale = log10,
#     yminorgridvisible = true,
#     xminorgridvisible = true,
#     title = "CSS Decoder - Simple",
#     xlabel = L"Depolarizing Probability $p$",
#     ylabel = "FER"
# )
# lines!(noise, FER, color = :red)
# scatter!(noise, FER, color = :red)
# current_figure()

# X_iters_15 = filter(((k, v), ) -> k <= 15, X_iters)
# Z_iters_15 = filter(((k, v), ) -> k <= 15, Z_iters)
# X_values = collect(values(X_iters_15))
# Z_values = collect(values(Z_iters_15))
# m = max(maximum(X_values), maximum(Z_values))
# fig2 = Figure()
# ax1 = Axis(fig2[1, 1], xlabel = "Iteration Number", ylabel = "Number Of Times Converged",
#         title = L"$X$ Convergence Distribution", limits = ((0, 16), (0, m + 500)))
# ax2 = Axis(fig2[1, 2], xlabel = "Iteration Number", ylabel = "Number Of Times Converged",
#     title = L"$Z$ Convergence Distribution", limits = ((0, 16), (0, m + 500)))
# barplot!(ax1, collect(keys(X_iters_15)), X_values, bar_width = 1)
# barplot!(ax2, collect(keys(Z_iters_15)), Z_values, bar_width = 1)


# X_iters_no_Bayes_15 = filter(((k, v), ) -> k <= 15, X_iters_no_Bayes)
# Z_iters_no_Bayes_15 = filter(((k, v), ) -> k <= 15, Z_iters_no_Bayes)
# X_values_no_Bayes = collect(values(X_iters_no_Bayes_15))
# Z_values_no_Bayes = collect(values(Z_iters_no_Bayes_15))
# m = max(maximum(X_values_no_Bayes), maximum(Z_values_no_Bayes))
# fig3 = Figure()
# ax1_no_Bayes = Axis(fig3[1, 1], xlabel = "Iteration Number", ylabel = "Number Of Times Converged",
#         title = L"$X$ Convergence Distribution", limits = ((0, 16), (0, m + 500)))
# ax2_no_Bayes = Axis(fig3[1, 2], xlabel = "Iteration Number", ylabel = "Number Of Times Converged",
#     title = L"$Z$ Convergence Distribution", limits = ((0, 16), (0, m + 500)))
# barplot!(ax1_no_Bayes, collect(keys(X_iters_no_Bayes_15)), X_values_no_Bayes, bar_width = 1)
# barplot!(ax2_no_Bayes, collect(keys(Z_iters_no_Bayes_15)), Z_values_no_Bayes, bar_width = 1)

# fig4 = Figure()
# ax4 = Axis(fig4[1, 1],
#     limits = ((noise[1], noise[end]), (1e-5, 1)),
#     xscale = log10,
#     yscale = log10,
#     yminorgridvisible = true,
#     xminorgridvisible = true,
#     title = "CSS Decoder - Simple",
#     xlabel = L"Depolarizing Probability $p$",
#     ylabel = "FER"
# )
# lines!(noise, FER, label = "Bayes", color = :red)
# scatter!(noise, FER, color = :red)
# lines!(noise, FER_no_Bayes, label = "No Bayes", color = :blue)
# scatter!(noise, FER_no_Bayes, color = :blue)
# axislegend(position = :rb)
# current_figure()

# Metachecks


function _make_single_shot_Tanner_graph(H::T, M::T) where T <: CodingTheory.CTMatrixTypes
    # internal function, so skip checks on correctness
    F = base_ring(H)
    Fone = F(1)
    nr_H, nc_H = size(H)
    nc_M = ncols(M)
    single_shot_matrix = CodingTheory.direct_sum(H, M)
    row = 1
    @inbounds @simd for c in nc_H + 1:nc_H + nc_M
        single_shot_matrix[row, c] = Fone
        row += 1
    end
    # of form (H I; 0 M)
    return single_shot_matrix
end

function CSS_metachecks(X_stabs::T, X_logs::T, X_meta::T, X_meta_L::T) where T <: CodingTheory.CTMatrixTypes
    nr, n = size(X_stabs)
    n == ncols(X_logs) || throw(ArgumentError("Logs are the wrong size"))
    ncols(X_meta) == nr || throw(ArgumentError("Metacheck is the wrong size"))
    iszero(X_meta * X_stabs) || throw(ArgumentError("Either metachecks or stabilizers wrong"))
    X_ss = _make_single_shot_Tanner_graph(X_stabs, X_meta)

    # noise details
    p_meas = 1e-3
    chn_meas =  MPNoiseModel(:BSC, p_meas)
    chn_inits_meta_in = Float64[log((1 - p_meas) / p_meas) for _ in 1:nr]
    syn_err_weights = Weights(Float64[1 - p_meas, p_meas])
    # noise = 10 .^(-3:.25:-1.3)
    noise = 10 .^(-1.5:.1:-1)
    len_noise = length(noise)

    # BP details
    num_runs = 10000
    max_iter = 20
    # 8 w/ measurement error, 1 perfect
    single_shot_rounds = 9
    schedule = :layered
    X_layers_stabs = layered_schedule(X_stabs)
    X_layers_meta = layered_schedule(X_meta)
    X_layers_meta_L = layered_schedule(X_meta_L)
    X_layers_ss = layered_schedule(X_ss)
    erasures = Int[]
    nr_logs = size(X_logs, 1)
    Pauli_types = Int[0, 1]
    v = zeros(Int, n)
    v2 = zeros(Int, nr)
    v3 = zeros(Int, n + nr)
    X_logs_Int = CodingTheory._Flint_matrix_to_Julia_int_matrix(X_logs)

    # multi-threading details
    num_threads = Threads.nthreads()
    runs_per_thread = cld(num_runs, num_threads)
    new_num_runs = runs_per_thread * num_threads
    println("Number of threads: $num_threads, runs per thread: $runs_per_thread, new number of runs: $new_num_runs")

    # preallocate places to store results
    FER_scheme_1 = zeros(Float64, len_noise)
    FER_scheme_2 = zeros(Float64, len_noise)
    local_log_failure_counts_scheme_1 = [zeros(Int, nr_logs) for _ in 1:num_threads]
    local_log_failure_counts_scheme_2 = [zeros(Int, nr_logs) for _ in 1:num_threads]
    local_failure_counts_scheme_1 = zeros(Int, num_threads)
    local_failure_counts_scheme_2 = zeros(Int, num_threads)
    for i in 1:len_noise
        p = noise[i]
        println("Starting p = $p")

        # setup noise model
        chn = MPNoiseModel(:BSC, p)
        Pauli_weights = Weights(Float64[1 - p, p])
        chn_inits_stabs_in = Float64[log((1 - p) / p) for _ in 1:n]
        chn_inits_ss_in = Float64[chn_inits_stabs_in; chn_inits_meta_in]

        Threads.@threads for th in 1:num_threads
            # Tanner graph for the stabilizers
            X_stabs_Int, _, var_adj_list_stabs, check_adj_list_stabs, chn_inits_stabs,
                check_to_var_messages_stabs, var_to_check_messages_stabs, current_bits_stabs,
                totals_stabs, syn_stabs = CodingTheory._message_passing_init(X_stabs, v, chn, max_iter, :SP,
                chn_inits_stabs_in, schedule, erasures)
            # Tanner graph for the metachecks
            X_meta_Int, _, var_adj_list_meta, check_adj_list_meta, chn_inits_meta,
                check_to_var_messages_meta, var_to_check_messages_meta, current_bits_meta,
                totals_meta, syn_meta = CodingTheory._message_passing_init(X_meta, v2, chn_meas, max_iter, :SP,
                chn_inits_meta_in, schedule, erasures)
            # Tanner graph for the second metacheck matrix
            X_meta_L_Int, _, var_adj_list_meta_L, check_adj_list_meta_L, chn_inits_meta_L,
                check_to_var_messages_meta_L, var_to_check_messages_meta_L, current_bits_meta_L,
                totals_meta_L, syn_meta_L = CodingTheory._message_passing_init(X_meta_L, v2, chn_meas, max_iter, :SP,
                chn_inits_meta_in, schedule, erasures)
            # Tanner graph for scheme 2
            X_ss_Int, _, var_adj_list_ss, check_adj_list_ss, chn_inits_ss, check_to_var_messages_ss,
                var_to_check_messages_ss, current_bits_ss, totals_ss, syn_ss =
                CodingTheory._message_passing_init(X_ss, v3, chn, max_iter, :SP, chn_inits_ss_in, schedule,
                erasures)
            X_meas_syn_scheme_1 = zeros(Int, nr)
            X_meas_syn_scheme_2 = zeros(Int, nr)
            X_syn_err = zeros(Int, nr)
            X_syn_temp = zeros(Int, nr)
            Z_err = zeros(Int, n)
            state_scheme_1 = zeros(Int, n)
            state_scheme_2 = zeros(Int, n)
            m = zeros(Int, size(X_meta_Int, 1))
            m_2 = zeros(Int, size(X_meta_Int, 1))
            m2 = zeros(Int, size(X_meta_L_Int, 1))
            logs_check = zeros(Int, nr_logs, 1)

            for _ in 1:runs_per_thread
                X_flag_stabs = false
                X_flag_ss = false
                X_out_stabs = Int[0]
                X_out_ss = Int[0]
                for ss_iter in 1:single_shot_rounds
                    # sample
                    sample!(Pauli_types, Pauli_weights, Z_err)
                    @inbounds @simd for j in i:n
                        state_scheme_1[j] += Z_err[j]
                        state_scheme_2[j] += Z_err[j]
                    end

                    # syndrome
                    LinearAlgebra.mul!(X_meas_syn_scheme_1, X_stabs_Int, state_scheme_1)
                    LinearAlgebra.mul!(X_meas_syn_scheme_2, X_stabs_Int, state_scheme_2)
                    @inbounds @simd for j in 1:nr
                        X_meas_syn_scheme_1[j] %= 2
                        X_meas_syn_scheme_2[j] %= 2
                    end

                    # last measurement is considered perfect
                    if ss_iter != single_shot_rounds
                        # add measurement errors
                        sample!(Pauli_types, syn_err_weights, X_syn_err)
                        @inbounds @simd for j in 1:nr
                            X_meas_syn_scheme_1[j] = (X_meas_syn_scheme_1[j] + X_syn_err[j]) % 2
                            X_meas_syn_scheme_2[j] = (X_meas_syn_scheme_2[j] + X_syn_err[j]) % 2
                        end
                    end

                    ###########
                    # Scheme 1
                    ###########

                    LinearAlgebra.mul!(m, X_meta_Int, X_meas_syn_scheme_1)
                    LinearAlgebra.mul!(m_2, X_meta_Int, X_meas_syn_scheme_2)
                    @inbounds @simd for j in 1:length(m)
                        m[j] %= 2
                        m_2[j] %= 2
                    end

                    if !iszero(m)
                        # decode with metachecks
                        X_flag_meta, X_out_meta, _, = CodingTheory._message_passing_layered(X_meta_Int,
                            X_meas_syn_scheme_1, chn_inits_meta,
                            CodingTheory._SP_check_node_message_box_plus, var_adj_list_meta,
                            check_adj_list_meta, max_iter, schedule, current_bits_meta,
                            totals_meta, syn_meta, check_to_var_messages_meta,
                            var_to_check_messages_meta, 0.0, X_layers_meta)
                        # println("meta flag: $X_flag_meta")

                        # reset inputs for next run, but don't re-allocate new memory
                        check_to_var_messages_meta .= 0.0
                        var_to_check_messages_meta .= 0.0

                        if X_flag_meta
                            @inbounds @simd for j in 1:nr
                                X_syn_temp[j] = (X_meas_syn_scheme_1[j] + X_out_meta[j]) % 2
                            end

                            # check if has element of homology
                            LinearAlgebra.mul!(m2, X_meta_L_Int, X_syn_temp)
                            @inbounds @simd for j in 1:length(m2)
                                m2[j] %= 2
                            end

                            if !iszero(m2)
                                # println("decoding other metachecks")
                                X_flag_meta_L, X_out_meta_L, _, = 
                                CodingTheory._message_passing_layered(X_meta_L_Int, X_meas_syn_scheme_1,
                                    chn_inits_meta_L, CodingTheory._SP_check_node_message_box_plus,
                                    var_adj_list_meta_L, check_adj_list_meta_L, max_iter, schedule,
                                    current_bits_meta_L, totals_meta_L, syn_meta_L,
                                    check_to_var_messages_meta_L, var_to_check_messages_meta_L,
                                    0.0, X_layers_meta_L)
                                # println("meta_L flag: $X_flag_meta_L")
    
                                # reset inputs for next run, but don't re-allocate new memory
                                check_to_var_messages_meta_L .= 0.0
                                var_to_check_messages_meta_L .= 0.0
    
                                # this is new and still doesn't seem to help
                                if X_flag_meta_L
                                    @inbounds @simd for j in i:nr
                                        X_meas_syn_scheme_1[j] += X_out_meta_L[j]
                                        X_meas_syn_scheme_1[j] %= 2
                                    end
                                else
                                    @inbounds @simd for j in i:nr
                                        X_meas_syn_scheme_1[j] = X_syn_temp[j]
                                    end
                                end    
                            else
                                @inbounds @simd for j in i:nr
                                    X_meas_syn_scheme_1[j] = X_syn_temp[j]
                                end
                            end
                        else
                            X_flag_meta_L, X_out_meta_L, _, = 
                            CodingTheory._message_passing_layered(X_meta_L_Int, X_meas_syn_scheme_1,
                                chn_inits_meta_L, CodingTheory._SP_check_node_message_box_plus,
                                var_adj_list_meta_L, check_adj_list_meta_L, max_iter, schedule,
                                current_bits_meta_L, totals_meta_L, syn_meta_L,
                                check_to_var_messages_meta_L, var_to_check_messages_meta_L,
                                0.0, X_layers_meta_L)

                            # reset inputs for next run, but don't re-allocate new memory
                            check_to_var_messages_meta_L .= 0.0
                            var_to_check_messages_meta_L .= 0.0

                            # this is new and still doesn't seem to help
                            if X_flag_meta_L
                                @inbounds @simd for j in i:nr
                                    X_meas_syn_scheme_1[j] += X_out_meta_L[j]
                                    X_meas_syn_scheme_1[j] %= 2
                                end
                            end
                        end
                    end

                    # decode stabilizers
                    X_flag_stabs, X_out_stabs, _, = CodingTheory._message_passing_layered(X_stabs_Int, 
                        X_meas_syn_scheme_1, chn_inits_stabs, CodingTheory._SP_check_node_message_box_plus,
                        var_adj_list_stabs, check_adj_list_stabs, max_iter, schedule,
                        current_bits_stabs, totals_stabs, syn_stabs, check_to_var_messages_stabs,
                        var_to_check_messages_stabs, 0.0, X_layers_stabs)

                    # reset inputs for next run, but don't re-allocate new memory
                    check_to_var_messages_stabs .= 0.0
                    var_to_check_messages_stabs .= 0.0

                    if X_flag_stabs
                        @inbounds @simd for j in 1:n
                            state_scheme_1[j] += X_out_stabs[j]
                        end
                    end

                    ###########
                    # Scheme 2
                    ###########

                    # decode single-shot matrix
                    X_flag_ss, X_out_ss, _, = CodingTheory._message_passing_layered(X_ss_Int, Int[X_meas_syn_scheme_2; m_2],
                        chn_inits_ss, CodingTheory._SP_check_node_message_box_plus, var_adj_list_ss,
                        check_adj_list_ss, max_iter, schedule, current_bits_ss, totals_ss,
                        syn_ss, check_to_var_messages_ss, var_to_check_messages_ss, 0.0, X_layers_ss)

                    # reset inputs for next run, but don't re-allocate new memory
                    check_to_var_messages_ss .= 0.0
                    var_to_check_messages_ss .= 0.0

                    # only take the bits corresponding to the data error
                    @inbounds @simd for j in 1:n
                        state_scheme_2[j] += X_out_ss[j]
                    end
                end

                ###########
                # Scheme 1
                ###########

                # check if scheme 1 converged - returned to the code subspace
                if X_flag_stabs
                    # converged
                    # check for logical error
                    logs_flag = false
                    LinearAlgebra.mul!(logs_check, X_logs_Int, state_scheme_1)
                    @inbounds for j in 1:nr_logs
                        logs_check[j, 1] %= 2
                        if logs_check[j, 1] == 1
                            logs_flag = true
                            local_log_failure_counts_scheme_1[th][j] += 1
                        end
                    end
                    if logs_flag
                        # introduced a logical error
                        local_failure_counts_scheme_1[th] += 1
                    end
                else
                    # did not converge
                    local_failure_counts_scheme_1[th] += 1
                end

                ###########
                # Scheme 2
                ###########

                # check if scheme 2 converged - returned to the code subspace
                if X_flag_ss
                    # converged
                    # check for logical error
                    logs_flag = false
                    LinearAlgebra.mul!(logs_check, X_logs_Int, state_scheme_2)
                    @inbounds for j in 1:nr_logs
                        logs_check[j, 1] %= 2
                        if logs_check[j, 1] == 1
                            logs_flag = true
                            local_log_failure_counts_scheme_2[th][j] += 1
                        end
                    end
                    if logs_flag
                        # introduced a logical error
                        local_failure_counts_scheme_2[th] += 1
                    end
                else
                    # did not converge
                    local_failure_counts_scheme_2[th] += 1
                end

                # reset inputs for next run, but don't re-allocate new memory
                state_scheme_1 .= 0
                state_scheme_2 .= 0
            end
        end
        
        # reduce iteration counts
        @inbounds for i in 2:num_threads
            @simd for j in 1:nr_logs
                local_log_failure_counts_scheme_1[1][j] += local_log_failure_counts_scheme_1[i][j]
                local_log_failure_counts_scheme_2[1][j] += local_log_failure_counts_scheme_2[i][j]
            end
        end

        println("logical failures 1: $(local_log_failure_counts_scheme_1[1])")
        println("logical failures 2: $(local_log_failure_counts_scheme_2[1])")
        println("failures 1: $local_failure_counts_scheme_1")
        println("failures 2: $local_failure_counts_scheme_2")
        @inbounds FER_scheme_1[i] = sum(local_failure_counts_scheme_1) / new_num_runs
        @inbounds FER_scheme_2[i] = sum(local_failure_counts_scheme_2) / new_num_runs
        println("FER 1 = $(FER_scheme_1[i])")
        println("FER 2 = $(FER_scheme_2[i])")
        println("Finished p = $p")

        # reset inputs for next run, but don't re-allocate new memory
        @inbounds for i in 1:num_threads
            local_log_failure_counts_scheme_1[i] .= 0
            local_log_failure_counts_scheme_2[i] .= 0
        end
        local_failure_counts_scheme_1 .= 0
        local_failure_counts_scheme_2 .= 0
    end
    return FER_scheme_1, FER_scheme_2
end
