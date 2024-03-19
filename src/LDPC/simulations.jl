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

function decoders_test(H::CTMatrixTypes)
    # initial parameters
    noise = [0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05]
    len_noise = length(noise)
    noise_model = :BSC
    num_runs = 10000
    seed = nothing
    # We use an explicit pseudoRNG with the given seed. (if `nothing`, it's just random)
    # Note that Xoshiro is default in Julia. Threading breaks reproducibility.
    rng = Xoshiro(seed)
    n = ncols(H)
    # it suffices to decode the zero codeword
    v = zero_matrix(base_ring(H), 1, n)
    max_iter = 10
    schedule = :flooding
    erasures = Int[]
    chn_inits = missing
    dist = Uniform(0.0, 1.0)
    err = zeros(Int, n)
    algorithm = :manual
    guided_rounds = 10
    attenuation = 1.0

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
    for (i, p) in enumerate(noise)
        println("Starting p = $p")
        # initialize everything for p
        chn = MPNoiseModel(noise_model, p)
        init_0 = log((1 - p) / p)
        init_1 = log(p / (1 - p))
        # SP and MS
        H_Int, v_Int, var_adj_list, check_adj_list, chn_inits, check_to_var_messages,
            var_to_check_messages, current_bits, totals, syn = _message_passing_init(H, v, chn,
            max_iter, :SP, 0, chn_inits, schedule, erasures)
        # for syndrome-based we need to change the inits
        chn_inits_syn = [log((1 - p) / p) for _ in 1:n]

        count_SP = 0
        count_SP_syn = 0
        count_SP_dec = 0
        count_MS = 0
        count_MS_syn = 0
        count_MS_dec = 0
        count_MS_C = 0
        count_MS_C_syn = 0
        count_MS_C_dec = 0

        for _ in 1:num_runs
            # sample
            @inbounds for j in 1:n
                rand(dist) ≤ p ? (err[j] = 1; chn_inits[j] = init_1;) : (err[j] = 0; chn_inits[j] = init_0;)
            end
            w = v_Int .+ err
            syn_Int = (H_Int * v_Int) .% 2
            # pick a random bit to decimate
            bit = rand(1:n)
            decimated_bits = [bit]
            # for now, setup as a genie-assisted decoder
            decimated_values = [err[bit]]

            # run SP
            flag, out, _, _, _ = _message_passing(H_Int, w, chn_inits, _SP_check_node_message_box_plus,
                var_adj_list, check_adj_list, max_iter, :SP, schedule, current_bits, totals, syn,
                check_to_var_messages, var_to_check_messages, 0, 0.0)
            (!flag || !iszero(out)) && (count_SP += 1;)

            # reset inputs for next run, but don't re-allocate new memory
            check_to_var_messages[:, :, :] .= 0.0
            var_to_check_messages[:, :, :] .= 0.0
            # current_bits[:] .= 0
            # totals[:] .= 0.0
            # syn[:] .= 0

            # run syndrome-based SP
            flag, out, _, _, _ = _message_passing_syndrome(H_Int, syn_Int, chn_inits_syn, _SP_check_node_message_box_plus,
            var_adj_list, check_adj_list, max_iter, :SP, schedule, current_bits, totals, syn, 
            check_to_var_messages, var_to_check_messages, 0, 0.0)
            (!flag || out ≠ err) && (count_SP_syn += 1;)

            # reset inputs for next run, but don't re-allocate new memory
            check_to_var_messages[:, :, :] .= 0.0
            var_to_check_messages[:, :, :] .= 0.0
            # current_bits[:] .= 0
            # totals[:] .= 0.0
            # syn[:] .= 0

            # run SP with decimation
            flag, out, _, _, _ = _message_passing_decimation(H_Int, w, chn_inits, _SP_check_node_message_box_plus,
            var_adj_list, check_adj_list, max_iter, :SP, schedule, decimated_bits, decimated_values,
            current_bits, totals, syn, check_to_var_messages, var_to_check_messages, 0, 0.0, algorithm,
            guided_rounds)
            (!flag || !iszero(out)) && (count_SP_dec += 1;)

            # reset inputs for next run, but don't re-allocate new memory
            check_to_var_messages[:, :, :] .= 0.0
            var_to_check_messages[:, :, :] .= 0.0
            # current_bits[:] .= 0
            # totals[:] .= 0.0
            # syn[:] .= 0

            # run MS
            flag, out, _, _, _ = _message_passing(H_Int, w, chn_inits, _MS_check_node_message, var_adj_list,
            check_adj_list, max_iter, :MS, schedule, current_bits, totals, syn, check_to_var_messages,
            var_to_check_messages, 0, attenuation)
            (!flag && !iszero(out)) && (count_MS += 1;)

            # reset inputs for next run, but don't re-allocate new memory
            check_to_var_messages[:, :, :] .= 0.0
            var_to_check_messages[:, :, :] .= 0.0
            # current_bits[:] .= 0
            # totals[:] .= 0.0
            # syn[:] .= 0

            # run syndrome-based MS
            flag, out, _, _, _ = _message_passing_syndrome(H_Int, w, chn_inits_syn, _MS_check_node_message,
            var_adj_list, check_adj_list, max_iter, :MS, schedule, current_bits, totals, syn, 
            check_to_var_messages, var_to_check_messages, 0, attenuation)
            (!flag && out ≠ err) && (count_MS_syn += 1;)

            # reset inputs for next run, but don't re-allocate new memory
            check_to_var_messages[:, :, :] .= 0.0
            var_to_check_messages[:, :, :] .= 0.0
            # current_bits[:] .= 0
            # totals[:] .= 0.0
            # syn[:] .= 0

            # run MS with decimation
            flag, out, _, _, _ = _message_passing_decimation(H_Int, w, chn_inits, _MS_check_node_message,
            var_adj_list, check_adj_list, max_iter, :MS, schedule, decimated_bits, decimated_values,
            current_bits, totals, syn, check_to_var_messages, var_to_check_messages, 0, 0.0, algorithm,
            guided_rounds)
            (!flag || !iszero(out)) && (count_MS_dec += 1;)

            # reset inputs for next run, but don't re-allocate new memory
            check_to_var_messages[:, :, :] .= 0.0
            var_to_check_messages[:, :, :] .= 0.0
            # current_bits[:] .= 0
            # totals[:] .= 0.0
            # syn[:] .= 0

            # run MS with correction
            flag, out, _, _, _ = _message_passing(H_Int, w, chn_inits, _MS_correction_check_node_message, var_adj_list,
            check_adj_list, max_iter, :MS, schedule, current_bits, totals, syn, check_to_var_messages,
            var_to_check_messages, 0, attenuation)
            (!flag || !iszero(out)) && (count_MS_C += 1;)

            # reset inputs for next run, but don't re-allocate new memory
            check_to_var_messages[:, :, :] .= 0.0
            var_to_check_messages[:, :, :] .= 0.0
            # current_bits[:] .= 0
            # totals[:] .= 0.0
            # syn[:] .= 0

            # run syndrome-based MS with correction
            flag, out, _, _, _ = _message_passing_syndrome(H_Int, syn_Int, chn_inits_syn, _MS_correction_check_node_message,
            var_adj_list, check_adj_list, max_iter, :MS, schedule, current_bits, totals, syn, 
            check_to_var_messages, var_to_check_messages, 0, attenuation)
            (!flag || out ≠ err) && (count_MS_C_syn += 1;)

            # reset inputs for next run, but don't re-allocate new memory
            check_to_var_messages[:, :, :] .= 0.0
            var_to_check_messages[:, :, :] .= 0.0
            # current_bits[:] .= 0
            # totals[:] .= 0.0
            # syn[:] .= 0

            # run MS with correction and decimation
            flag, out, _, _, _ = _message_passing_decimation(H_Int, w, chn_inits, _MS_check_node_message,
            var_adj_list, check_adj_list, max_iter, :MS, schedule, decimated_bits, decimated_values,
            current_bits, totals, syn, check_to_var_messages, var_to_check_messages, 0, 0.0, algorithm,
            guided_rounds)
            (!flag || !iszero(out)) && (count_MS_C_dec += 1;)

            # reset inputs for next run, but don't re-allocate new memory
            check_to_var_messages[:, :, :] .= 0.0
            var_to_check_messages[:, :, :] .= 0.0
            # current_bits[:] .= 0
            # totals[:] .= 0.0
            # syn[:] .= 0
        end

        FER_SP[i] = count_SP / num_runs
        println("FER_SP = $(FER_SP[i])")
        FER_SP_syn[i] = count_SP_syn / num_runs
        println("FER_SP_syn = $(FER_SP_syn[i])")
        FER_SP_dec[i] = count_SP_dec / num_runs
        println("FER_SP_dec = $(FER_SP_dec[i])")
        FER_MS[i] = count_MS / num_runs
        println("FER_MS = $(FER_MS[i])")
        FER_MS_syn[i] = count_MS_syn / num_runs
        println("FER_MS_syn = $(FER_MS_syn[i])")
        FER_MS_dec[i] = count_MS_dec / num_runs
        println("FER_MS_dec = $(FER_MS_dec[i])")
        FER_MS_C[i] = count_MS_C / num_runs
        println("FER_MS_C = $(FER_MS_C[i])")
        FER_MS_C_syn[i] = count_MS_C_syn / num_runs
        println("FER_MS_C_syn = $(FER_MS_C_syn[i])")
        FER_MS_C_dec[i] = count_MS_C_dec / num_runs
        println("FER_MS_C_dec = $(FER_MS_C_dec[i])")
        println("Finished p = $p")
    end
    return FER_SP, FER_SP_syn, FER_SP_dec, FER_MS, FER_MS_syn, FER_MS_dec, FER_MS_C, FER_MS_C_syn,
        FER_MS_C_dec
end

function single_decoder_test(H::CTMatrixTypes)
    # initial parameters
    noise = [0.02]
    # noise = [0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05]
    len_noise = length(noise)
    noise_model = :BSC
    num_runs = 5000
    nr, n = size(H)
    # it suffices to decode the zero codeword
    v = zero_matrix(base_ring(H), 1, n)
    syn_Int = zeros(Int, nr)
    max_iter = 50
    schedule = :flooding
    erasures = Int[]
    chn_inits = zeros(Float64, n)
    dist = Uniform(0.0, 1.0)
    err = zeros(Int, n)
    algorithm = :manual
    guided_rounds = 10
    attenuation = 0.5

    FER = zeros(Float64, len_noise)
    for (i, p) in enumerate(noise)
        # println("Starting p = $p")
        # initialize everything for p
        chn = MPNoiseModel(noise_model, p)
        init_0 = log((1 - p) / p)
        # init_1 = log(p / (1 - p))
        # SP
        # H_Int, v_Int, var_adj_list, check_adj_list, chn_inits, check_to_var_messages,
        #     var_to_check_messages, current_bits, totals, syn = _message_passing_init(H, v, chn,
        #     max_iter, :SP, 0, chn_inits, schedule, erasures)
        # MS
        H_Int, v_Int, var_adj_list, check_adj_list, chn_inits, check_to_var_messages,
            var_to_check_messages, current_bits, totals, syn = _message_passing_init(H, v, chn,
            max_iter, :MS, 0, chn_inits, schedule, erasures)
        chn_inits_syn = [init_0 for _ in 1:n]

        count = 0
        for _ in 1:num_runs
            # sample
            # err_wt = 0
            # @inbounds for j in 1:n
            #     rand(dist) ≤ p ? (err[j] = 1; chn_inits[j] = init_1;) : (err[j] = 0; chn_inits[j] = init_0;)
            # end
            @inbounds for j in 1:n
                rand(dist) ≤ p ? (err[j] = 1;) : (err[j] = 0;)
            end
            # w = v_Int .+ err
            # print("$err_wt, ")
            # syn_Int = (H_Int * err) .% 2
            LinearAlgebra.mul!(syn_Int, H_Int, err)
            @inbounds @simd for i in 1:nr
                syn_Int[i] %= 2
            end
            # pick a random bit to decimate
            # bit = rand(1:n)
            # decimated_bits = [bit]
            # # for now, setup as a genie-assisted decoder
            # decimated_values = [err[bit]]

            # SP
            # flag, out, iter, _, _ = _message_passing(H_Int, err, chn_inits, _SP_check_node_message_box_plus,
            #     var_adj_list, check_adj_list, max_iter, :SP, schedule, current_bits, totals, syn,
            #     check_to_var_messages, var_to_check_messages, 0, 0.0)
            # (!flag || !iszero(out)) && (count += 1;)
            # MS
            # flag, out, iter, _, _ = _message_passing(H_Int, err, chn_inits, _MS_check_node_message,
            #     var_adj_list, check_adj_list, max_iter, :MS, schedule, current_bits, totals, syn,
            #     check_to_var_messages, var_to_check_messages, 0, 0.0)
            # (!flag || !iszero(out)) && (count += 1;)
            # SP syn
            flag, out, _, = _message_passing_syndrome(H_Int, syn_Int, chn_inits_syn,
                _SP_check_node_message_box_plus, var_adj_list, check_adj_list, max_iter, :SP, schedule,
                current_bits, totals, syn, check_to_var_messages, var_to_check_messages, 0, 0.0)
            (!flag || out ≠ err) && (count += 1;)
            # !flag && println("err = $err, out = $out, - $flag")
            # print("$iter, ")
            # run SP with decimation
            # flag, out, _, _, _ = _message_passing_decimation(H_Int, err, chn_inits,
            #     _SP_check_node_message_box_plus, var_adj_list, check_adj_list, max_iter, :SP, schedule,
            #     decimated_bits, decimated_values, current_bits, totals, syn, check_to_var_messages,
            #     var_to_check_messages, 0, 0.0, algorithm, guided_rounds)
            # (!flag || !iszero(out)) && (count += 1;)
            # run MS with decimation
            # flag, out, _, _, _ = _message_passing_decimation(H_Int, err, chn_inits, 
            #     _MS_check_node_message, var_adj_list, check_adj_list, max_iter, :MS, schedule,
            #     decimated_bits, decimated_values, current_bits, totals, syn, check_to_var_messages, 
            #     var_to_check_messages, 0, attenuation, algorithm, guided_rounds)
            # (!flag || !iszero(out)) && (count += 1;)
            # run syndrome-based MS
            # flag, out, _, _, _ = _message_passing_syndrome(H_Int, syn_Int, chn_inits_syn, 
            #     _MS_check_node_message, var_adj_list, check_adj_list, max_iter, :MS, schedule, 
            #     current_bits, totals, syn, check_to_var_messages, var_to_check_messages, 0, 
            #     attenuation)
            # (!flag && out ≠ err) && (count += 1;)

            # reset inputs for next run, but don't re-allocate new memory
            check_to_var_messages[:, :, :] .= 0.0
            var_to_check_messages[:, :, :] .= 0.0
            # current_bits[:] .= 0
            # totals[:] .= 0.0
            # syn[:] .= 0

        end
        @inbounds FER[i] = count / num_runs
        # println("FER = $(FER[i])")
        # println("Finished p = $p")
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
