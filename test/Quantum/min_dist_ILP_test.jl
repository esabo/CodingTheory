@testitem "Quantum/min_dist_ILP.jl" begin
    using CodingTheory, Oscar, SparseArrays
    try
        using HiGHS, JuMP
    catch
        @info "Skipping CSS ILP tests; JuMP/HiGHS are not available."
        return
    end

    @testset "Quotient-aware parity formulation" begin
        H = sparse([1, 1], [2, 3], [1, 1], 1, 3)
        logical_checks = sparse([1], [1], [1], 1, 3)
        d, witness = CodingTheory._minimum_distance(
            H, logical_checks; alg=:ILP, max_d=3,
            time_limit_sec=nothing, ilp_parity_cut_max_degree=10,
            ilp_optimizer=:HiGHS, ilp_threads=1, ilp_cyclic_period=3)
        @test d == 1
        @test witness == [1, 0, 0]
        @test_throws ArgumentError CodingTheory._minimum_distance(
            H, logical_checks; alg=:ILP, max_d=3, ilp_optimizer=:HiGHS,
            ilp_cyclic_period=2)
    end

    @testset "Auto selects the available ILP for large sparse kernels" begin
        H = spzeros(Int, 1, 96)
        H[1, 96] = 1
        logical_checks = spzeros(Int, 1, 96)
        logical_checks[1, 1] = 1
        d, witness = CodingTheory._minimum_distance(
            H, logical_checks; alg=:auto, max_d=1,
            time_limit_sec=nothing)
        @test d == 1
        @test witness[1] == 1
        @test sum(witness) == 1
    end
end
