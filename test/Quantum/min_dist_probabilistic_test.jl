@testitem "Quantum/min_dist_probabilistic.jl" begin
    using CodingTheory, Oscar

    F = Oscar.Nemo.Native.GF(2)
    H = matrix(F, [
        1 0 1 0 1 0 1;
        0 1 1 0 0 1 1;
        0 0 0 1 1 1 1
    ])

    function check_isd_result(alg::Symbol; kwargs...)
        S = CSSCode(H, H)
        d, witness = probabilistic_minimum_distance(
            S; which=:X, alg=alg, max_weight=3,
            max_iters=500, seed=20260909, kwargs...)
        @test d == 3
        @test wt(witness) == d
        @test is_logical(S, witness)
        @test S.cache[:u_bound_dx] == 3
    end

    @testset "Quotient-aware binary ISD" begin
        check_isd_result(:Prange)
        check_isd_result(:LeeBrickell; p=2)
        check_isd_result(:Stern; p=1, l=2)
    end

    @testset "Information-set and automorphism preprocessors" begin
        S = CSSCode(H, H)
        identity_automorphism = [collect(1:S.n)]
        d, witness = probabilistic_minimum_distance(
            S; which=:Z, alg=:Prange, max_weight=3,
            max_iters=500, info_set_alg=:Brouwer,
            automorphisms=identity_automorphism, seed=7)
        @test d == 3
        @test wt(witness) == d
        @test is_logical(S, witness)
    end

    @testset "Failure is only an upper-bound miss" begin
        S = CSSCode(H, H)
        d, witness = probabilistic_minimum_distance(
            S; which=:X, alg=:LeeBrickell, p=2,
            max_weight=2, max_iters=5, seed=1)
        @test d == -1
        @test iszero(witness)
        @test !haskey(S.cache, :dx)
    end
end
