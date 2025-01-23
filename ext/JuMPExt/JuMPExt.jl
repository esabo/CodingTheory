module JuMPExt

import CodingTheory
import CodingTheory: optimal_lambda, optimal_rho, optimal_lambda_and_rho, LP_decoder_LDPC, AbstractLinearCode, BinarySymmetricChannel
import JuMP
import JuMP: @variable, @constraint, @objective
import GLPK
import Oscar

include("LDPC/decoders.jl")
include("LDPC/analysis.jl")

end
