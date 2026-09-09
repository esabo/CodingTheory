@testitem "Quantum/min_dist_ILP.jl" begin
    using CodingTheory, GLPK, JuMP, Oscar, SparseArrays

    @testset "Quotient-aware parity formulation" begin
        H = sparse([1, 1], [2, 3], [1, 1], 1, 3)
        logical_checks = sparse([1], [1], [1], 1, 3)
        d, witness = CodingTheory._minimum_distance(
            H, logical_checks; alg=:ILP, max_d=3,
            time_limit_sec=nothing, ilp_parity_cut_max_degree=10)
        @test d == 1
        @test witness == [1, 0, 0]
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
