module JuMPExt

using CodingTheory
import CodingTheory: optimal_lambda, optimal_rho, optimal_lambda_and_rho, LP_decoder_LDPC, AbstractChannel, AbstractLinearCode, BinarySymmetricChannel, _minimum_distance_ILP, _minimum_distance_css_HiGHS, parity_check_matrix
using JuMP
# import JuMP: @variable, @constraint, @objective
using HiGHS
using Oscar
# import Oscar: order

include("Classical/min_dist_exact.jl")
include("Quantum/min_dist_exact.jl")
include("LDPC/decoders.jl")
include("LDPC/analysis.jl")

end
