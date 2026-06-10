module JuMPExt

using CodingTheory
import CodingTheory: optimal_lambda, optimal_rho, optimal_lambda_and_rho, LP_decoder_LDPC, AbstractLinearCode, BinarySymmetricChannel, _minimum_distance_ILP, parity_check_matrix
using JuMP
# import JuMP: @variable, @constraint, @objective
using GLPK
using Oscar
# import Oscar: order

include("Classical/min_dist_exact.jl")
include("LDPC/decoders.jl")
include("LDPC/analysis.jl")

end
