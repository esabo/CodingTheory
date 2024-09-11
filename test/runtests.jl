using Oscar
using CodingTheory
using TestItemRunner

if get(ENV, "GPU_TESTS", "") != "true"
    println("skipping gpu tests (set GPU_TESTS=true to test gpu)")
end

# filter for the test
testfilter = ti -> begin
  exclude = Symbol[]
  if get(ENV,"JET_TEST","")!="true"
    push!(exclude, :jet)
  end
  if !(VERSION >= v"1.10")
    push!(exclude, :doctests)
    push!(exclude, :aqua)
  end

  if get(ENV, "GPU_TESTS", "")!="true"
    push!(exclude, :gpu)
  end

  if !(Base.Sys.islinux() & (Int===Int64))
    push!(exclude, :bitpack)
  end

  return all(!in(exclude), ti.tags)
end

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

@run_package_tests filter=testfilter

# TODO: should setup test for traits
# TODO: add tests for _standardformstabilizer
# TODO: add tests for _logicalsstandardform

# passes
# include("utils_test.jl")

# BUG go back and fix expanded codes here
# include("Classical/linear_code_test.jl")
# passes upto trellis
# include("Classical/ReedMuller_test.jl")
# works up to trellis code
# include("Classical/misc_known_codes_test.jl")
# passes
# include("Classical/cyclotomic_test.jl")
# passes up to min dist/wt enum
# include("Classical/cyclic_code_test.jl")
# passes
# include("Classical/quasi-cyclic_code_test.jl")
# passes
# include("Classical/tilings_Tanner_test.jl")
# BUG expanded codes above and min dist
# include("Classical/concatenation_test.jl")

# passes
# include("LDPC/codes_test.jl")

# passes
# include("Quantum/stabilizer_code_test.jl")
# include("Quantum/misc_known_codes_test.jl")
# include("Quantum/product_codes_test.jl")
# passed
# include("Quantum/subsystem_code_test.jl")

# currently skipping
# include("Classical/GeneralizedReedSolomon_test.jl")
# include("Quantum/quantum_MC_ids_test.jl")
# include("chain_complex_test.jl")
