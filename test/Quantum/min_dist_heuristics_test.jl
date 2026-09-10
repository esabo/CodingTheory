@testitem "Quantum/min_dist_heuristics.jl" begin
    using CodingTheory, Oscar

    F = Oscar.Nemo.Native.GF(2)
    H = matrix(F, [
        1 0 1 0 1 0 1;
        0 1 1 0 0 1 1;
        0 0 0 1 1 1 1
    ])

    function check_heuristic(alg::Symbol; kwargs...)
        S = CSSCode(H, H)
        set_X_minimum_distance_lower_bound!(S, 3)
        d, witness = heuristic_minimum_distance(
            S; which=:X, alg=alg, max_iters=500,
            pop_size=50, seed=20260909, kwargs...)
        @test d == 3
        @test wt(witness) == 3
        @test is_logical(S, witness)
        @test S.cache[:dx] == 3
    end

    @testset "Quotient-aware heuristic upper bounds" begin
        check_heuristic(:GGAOrder)
        check_heuristic(:NNCS)
        check_heuristic(:GA)
        check_heuristic(:ACO)
    end
end
