using Test
using Oscar
using CodingTheory

# TODO: should setup test for traits
# TODO: add tests for _standardformstabilizer
# TODO: add tests for _logicalsstandardform

include("utils_test.jl")

include("Classical/linear_code_test.jl")
include("Classical/ReedMuller_test.jl")
include("Classical/misc_known_codes_test.jl")
include("Classical/cyclotomic_test.jl")
include("Classical/cyclic_code_test.jl")
# include("Classical/GeneralizedReedSolomon_test.jl")
include("Classical/quasi-cyclic_code_test.jl")
include("Classical/tilings_Tanner_test.jl")
include("Classical/concatenation_test.jl")

include("LDPC/codes_test.jl")

include("Quantum/stabilizer_code_test.jl")
# include("Quantum/quantum_MC_ids_test.jl")
include("Quantum/misc_known_codes_test.jl")
include("Quantum/product_codes_test.jl")
include("Quantum/subsystem_code_test.jl")
# # include("chain_complex_test.jl")
