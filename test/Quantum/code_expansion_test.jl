@testitem "Quantum/code_expansion.jl" begin
    using CodingTheory
    import Graphs
    using LinearAlgebra
    using SparseArrays

    @testset "Graph spectral expansion" begin
        G = Graphs.path_graph(4)
        @test Matrix(laplacian_matrix(G)) == [
            1 -1 0 0
            -1 2 -1 0
            0 -1 2 -1
            0 0 -1 1
        ]
        @test algebraic_connectivity(G) ≈ 2 - sqrt(2)
        @test normalized_spectral_gap(G) > 0
        @test length(fiedler_vector(G)) == 4
        @test estimated_edge_expansion(G) ≈ 0.5
        @test estimated_vertex_expansion(G) ≈ 0.5
        lower, upper = edge_expansion_bounds(G)
        @test lower <= 0.5 <= upper

        disconnected = Graphs.SimpleGraph(3)
        Graphs.add_edge!(disconnected, 1, 2)
        @test algebraic_connectivity(disconnected) ≈ 0 atol=1e-12
        @test normalized_spectral_gap(disconnected) ≈ 0 atol=1e-12

        C4 = Graphs.cycle_graph(4)
        @test nontrivial_adjacency_spectral_radius(C4) ≈ 0 atol=1e-12
        @test edge_expansion_bounds(Graphs.path_graph(2)) == (1.0, 2.0)
    end

    @testset "Exact bipartite expansion" begin
        H = Bool[
            1 1 0
            0 1 1
        ]
        @test bipartite_expansion_profile(H; max_subset_size=2) == [1.0, 1.0]
        @test bipartite_expansion_profile(
            sparse(H); max_subset_size=2) == [1.0, 1.0]
        @test is_expander(H, 2 // 3, 1)
        witness = expansion_witness(H, 2 // 3, 11 // 10)
        @test !isnothing(witness)
        @test length(witness) <= 2

        G = Graphs.SimpleGraph(5)
        for (left, right) in ((1, 4), (2, 4), (2, 5), (3, 5))
            Graphs.add_edge!(G, left, right)
        end
        @test is_bipartite_expander(G, 1:3, 4:5, 2 // 3, 1)
        @test is_left_right_expander(H, 1 // 3, 1, 1 // 2, 2)
        @test estimated_bipartite_vertex_expansion(
            H; max_subset_size=2) >= 1
        @test_throws DomainError is_expander(H, 1.1, 1)
        @test_throws ArgumentError is_expander([1 2], 1, 1)
    end

    @testset "Reduced syndrome profiles" begin
        H = Bool[
            1 1
            1 1
            0 1
        ]
        profile = confinement_profile(H; max_error_weight=2)
        @test profile == Dict(0 => 0, 1 => 2, 2 => 1, 3 => 1)
        @test deterministic_QLTC_soundness(H; max_error_weight=2) == 0.5
        @test verify_QLTC_soundness(H, 0.5; max_error_weight=2)
        @test !verify_QLTC_soundness(H, 0.51; max_error_weight=2)
        @test evaluate_single_shot_soundness(
            H, 2; max_error_weight=2) == Dict(0 => 0, 1 => 2)
        @test evaluate_confinement(H, 2) == profile
        @test !verify_confinement(H, 1; max_error_weight=2)
    end

    @testset "Cosystolic expansion" begin
        boundary = Bool[1 1]
        no_incoming = falses(2, 0)
        @test cosystolic_expansion(
            boundary, no_incoming; max_weight=1) == 1
        @test cosystolic_expansion(
            boundary, no_incoming; max_weight=2) == 0

        incoming = reshape(Bool[1, 1], 2, 1)
        @test cosystolic_expansion(
            boundary, incoming; max_weight=2) == 1
        @test_throws ArgumentError cosystolic_expansion(
            Bool[1 0], reshape(Bool[1, 1], 2, 1))
    end

    @testset "Quantum wrappers" begin
        S = Q422()
        @test length(fiedler_vector(S, :X)) > 0
        @test is_topologically_connected(S, :X)
        @test is_expander(S, 1 // 4, 1, :X)
        @test deterministic_QLTC_soundness(
            S, :X; max_error_weight=1) >= 1
        @test_throws ArgumentError is_expander(FiveQubitCode(), 1 // 5, 1, :X)
    end

    @testset "Sipser--Spielman consequences" begin
        H = Bool[
            1 0 1
            1 1 0
            0 1 1
        ]
        result = sipser_spielman_guarantees(H, 2 // 3, 1)
        @test result.degree == 2
        @test result.epsilon == 1 // 2
        @test !result.has_linear_distance
        @test_throws ArgumentError sipser_spielman_guarantees(
            H, 2 // 3, 2)
        @test_throws ArgumentError sipser_spielman_guarantees(
            Bool[1 1; 0 1], 1 // 2, 1)
    end
end
