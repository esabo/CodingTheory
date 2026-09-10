@testitem "Quantum/min_dist_bounds.jl" begin
    using CodingTheory, Oscar

    F = Oscar.Nemo.Native.GF(2)
    H = matrix(F, [
        1 0 1 0 1 0 1;
        0 1 1 0 0 1 1;
        0 0 0 1 1 1 1
    ])
    S = CSSCode(H, H)
    X_logical = matrix(F, 1, 14, [
        1, 1, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0,
    ])

    @testset "Sector bounds close with a witness" begin
        @test is_logical(S, X_logical)
        set_X_minimum_distance_upper_bound!(S, 3, X_logical)
        set_X_minimum_distance_lower_bound!(S, 2)
        @test X_minimum_distance_lower_bound(S) == 2
        @test X_minimum_distance_upper_bound(S) == 3
        @test !haskey(S.cache, :dx)

        set_X_minimum_distance_lower_bound!(S, 3)
        @test S.cache[:dx] == 3
        @test S.cache[:X_minimum_distance_witness] == X_logical
    end

    @testset "Full lower bounds propagate to both sectors" begin
        T = CSSCode(H, H)
        set_minimum_distance_lower_bound!(T, 2)
        @test X_minimum_distance_lower_bound(T) == 2
        @test Z_minimum_distance_lower_bound(T) == 2
        @test minimum_distance_lower_bound(T) == 2
    end

    @testset "Exact search proves below a witnessed incumbent" begin
        T = CSSCode(H, H)
        set_X_minimum_distance_upper_bound!(T, 3, X_logical)
        d, witness = minimum_distance(T; which=:X, alg=:Gray, max_d=3)
        @test d == 3
        @test witness == X_logical
        @test X_minimum_distance_lower_bound(T) == 3
    end

    @testset "Full witnesses use Pauli weight" begin
        T = CSSCode(H, H)
        Y_logical = X_logical + matrix(F, 1, 14, [
            0, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 0, 0, 0, 0,
        ])
        set_minimum_distance_upper_bound!(T, 3, Y_logical)
        @test minimum_distance_upper_bound(T) == 3
    end

    @testset "Automorphism registration is validated" begin
        identity_permutation = collect(1:S.n)
        set_distance_automorphisms!(S, [identity_permutation])
        @test distance_automorphisms(S) == [identity_permutation]
        @test_throws ArgumentError set_distance_automorphisms!(
            S, [[1, 1, 3, 4, 5, 6, 7]])
    end
end
