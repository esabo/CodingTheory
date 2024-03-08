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
