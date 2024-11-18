using Oscar, CodingTheory
using Combinatorics
using LinearAlgebra
using BenchmarkTools
using Random
using Dates
# using AllocCheck
# using Profile
CodingTheory.Random.seed!(0)

function white_n70_k35(use_mine, terse_description)
    F = Oscar.Nemo.Native.GF(2)
    mat = matrix(F, [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 0 0 0 0 1 1 0 1 1 1 0 1 0 1 1 1 1 0 1 1 0 1 1 0 0 1 0 0 1 1 0 1;
    0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 1 1 1 1 1 0 1 1 0 1 0 0 1 1 0 1 0 1 1 0 1 1 0 1 1 1 0 0 0 1 0;
    0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 1 1 0 0 1 1 0 0 0 0 0 1 0 0 1 1 1 0 1 1 0 1 0 0 1 1 1 1 0 0;
    0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 1 0 0 0 0 0 0 1 0 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 0 0 1;
    0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 1 0 1 0 0 0 0 1 0 0 0 1 0 1 1 0 1 0 0 0 1 1 0 1 0 0 0 0 0 0 1 0;
    0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 0 1 1 1 0 0 0 0 1 1 0 1 0 1 0 1 0 1 0 0 0 1 1 0 0 1 0 0 1 0 1 1;
    0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 1 0 0 0 0 1 1 0 0 1 0 0 1 1 1 1 1 0 0 1 1 1 0 1 1 1 0 0 0 0 1;
    0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 1 1 1 0 1 0 1 0 0 0 0 0 0 0 1 0 0 0 1 0 1 0 0 1 1 1 1 1 1 0;
    0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 1 1 0 1 1 1 0 1 0 0 0 1 1 1 1 0 0 1 0 0 0 1 0 1 1 1 1 1 0 1 0 1;
    0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 0 1 0 1 1 1 0 0 0 0 1 0 0 1 0 0 1 1 1 1 0 0 1 1 1 1 1 0;
    0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 0 1 0 1 1 1 0 0 0 0 1 0 0 1 0 0 1 1 1 1 0 0 1 1 1 1 1;
    0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 1 0 0 1 1 1 0 1 1 0 1 0 0 0 1 0 0 1 0 0 0 1 0 0 0 1 1 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 1 0 0 0 1 0 0 0 0 0 0 1 1 0 1 1 0 0 1 0 0 1 0 0 0 0 0 1 0 1 1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 1 1 0 0 1 0 1 0 0 1 0 1 1 1 0 1 0 0 1 0 0 1 1 0 0 0 0 0 1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 1 1 1 0 1 0 0 1 1 0 0 0 0 1 0 0 0 0 1 0 0 1 0 0 1 0 0 0 1 1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 1 1 0 0 0 1 0 0 0 1 1 1 1 1 0 1 0 0 1 1 0 0 1 1 1 1 1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 1 0 1 1 0 1 0 1 1 0 0 0 1 0 1 1 0 0 1 0 1 1 1 1 0 1 0 0 0 1 1 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 1 1 1 0 1 0 0 1 1 0 0 0 0 1 0 0 0 0 1 0 0 1 0 0 1 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 0 0 1 1 0 1 1 0 0 1 1 1 0 0 0 1 1 0 1 0 0 1 1 0 1 1 0 1 1 1 1 1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 1 0 1 1 0 1 1 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 1 1 1 0 0 1 1 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 1 0 1 1 0 1 0 1 1 0 0 0 1 0 1 1 0 0 1 0 1 1 1 1 0 1 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 0 1 0 0 0 0 0 0 0 1 1 1 1 0 1 1 1 1 1 0 1 1 1 1 0 0 1 1 0 1 1 1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 1 0 1 1 0 1 0 1 0 0 1 1 1 0 0 1 0 1 1 0 0 0 0 0 1 0 0 1 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 1 0 1 1 0 1 1 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 1 1 1 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 1 1 0 1 0 0 1 1 1 0 1 0 1 1 1 0 1 0 1 0 1 0 0 1 1 0 0 0 0 0 1 0 1 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 1 0 1 1 0 1 0 1 0 0 1 1 1 0 0 1 0 1 1 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 1 0 1 1 0 1 0 1 0 0 1 1 1 0 0 1 0 1 1 0 0 0 0 0 1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 1 0 0 0 1 1 1 0 0 0 1 1 1 1 0 0 0 1 0 0 0 1 0 0 0 1 1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 1 0 1 0 1 1 1 0 0 1 1 1 1 1 1 0 0 0 0 1 0 1 0 1 0 0 1 0 1 0 0 1 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 1 0 1 0 1 0 0 0 1 1 1 0 1 0 1 0 0 1 1 0 1 0 0 0 1 1 0 1 1 0 0 1 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 1 0 1 1 1 1 0 0 0 1 1 0 0 0 1 0 1 0 1 1 0 1 0 0 0 1 0 1 1 1 1 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 1 1 0 1 1 0 0 0 1 0 0 1 1 0 1 1 0 1 1 0 1 0 0 0 0 1 1 1 0 1 1 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 1 1 0 0 0 0 0 1 1 1 0 0 0 1 0 0 0 0 0 0 0 0 1 0 1 0 1 1 1 1 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 1 1 1 1 1 0 1 1 0 1 0 0 1 1 0 1 0 1 1 0 1 1 0 1 1 1 0 0 0 1 0 0 1 1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 1 1 0 1 1 0 0 0 1 0 0 1 1 0 1 1 0 1 1 0 1 0 0 0 0 1 1 1])

    num_thrds = Threads.nthreads()
    println("Benchmarking with $num_thrds threads.")
    if use_mine
        function_description = "_minimum_distance_BZ(codeN70K35,:Brouwer;false)"
        terse_description ? println(function_description) : println(benchmark_description(function_description))
        return @benchmark CodingTheory._minimum_distance_BZ(codeN70K35, info_set_alg = :Brouwer; verbose = false) setup = (codeN70K35=LinearCode($mat, false, false))
    else
        function_description = "_minimum_distance_enumeration_with_matrix_multiply(codeN70K35,:Brouwer;false)"
        terse_description ? println(function_description) : println(benchmark_description(function_description))
        return @benchmark CodingTheory._minimum_distance_enumeration_with_matrix_multiply(codeN70K35, info_set_alg = :Brouwer; verbose = false) setup = (codeN70K35=LinearCode($mat, false, false))
    end
end

function benchmark_description(function_description)
    return "@" * Dates.format(now(), dateformat"MM:SS") * " benchmark " * function_description * 
    "\nMost recent git commit\n" * current_git_hash() * "\ndate: " * Dates.format(now(), dateformat"y-m-d : H:M:S")
end

function current_git_hash()
    try
        # Run the Git command to get the latest commit hash
        hash = readchomp(`git rev-parse --short HEAD`)
        return hash
    catch
        # Return a default value if not in a Git repository or Git fails
        return "no-git-repo"
    end
end

use_new=false
t = white_n70_k35(use_new, true)
# t = test_white_105(use_mine)
t