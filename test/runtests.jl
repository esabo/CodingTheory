using Oscar
using CodingTheory
using TestItemRunner

if get(ENV, "GPU_TESTS", "") != "true"
    println("skipping gpu tests (set GPU_TESTS = true to test gpu)")
end


# filter for the test
test_filter = ti -> begin
    exclude = Symbol[]
    if get(ENV, "JET_TEST", "") != "true"
        push!(exclude, :jet)
    end
    if !(VERSION >= v"1.10")
        push!(exclude, :doctests)
        push!(exclude, :aqua)
    end

    if get(ENV, "GPU_TESTS", "") != "true"
        push!(exclude, :gpu)
    end

    if !(Base.Sys.islinux() & (Int === Int64))
        push!(exclude, :bitpack)
    end

    return all(!in(exclude), ti.tags)
end

# since unit tests need to be deterministic this saves us having to set it every time we test a random function 
println("Random.seed initialized to 0")
CodingTheory.Random.seed!(0)

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

@run_package_tests filter = test_filter

# TODO: should setup test for traits
# TODO: add tests for _standardformstabilizer
# TODO: add tests for _logicalsstandardform
