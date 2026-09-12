@testitem "LDPC/cycles.jl" begin
    using Oscar, CodingTheory

    @testset "Adjacencies and Girth Basics" begin
        F = Oscar.Nemo.Native.GF(2)
        
        # Matrix with exactly one 4-cycle on (v1, v2) and (c1, c2)
        H_4 = matrix(F, [
            1 1 0 0;
            1 1 0 0;
            0 0 1 1;
        ])
        C_4 = LDPCCode(H_4)

        check_adj, var_adj = CodingTheory.node_adjacencies(C_4)
        @test length(check_adj) == 3
        @test length(var_adj) == 4
        @test check_adj[1] == [1, 2]
        @test var_adj[1] == [1, 2]

        @test girth(C_4) == 4
        
        # Test local girth (computation tree depth)
        @test local_girth(C_4, 1) == 4
        @test local_girth(C_4, 3) == -1  # Node 3 is part of a strict tree (no cycles reachable)
        @test local_girth(C_4, [1, 3]) == [4, -1]
        
        # Pure Tree matrix (no cycles)
        H_tree = matrix(F, [
            1 1 0;
            0 1 1
        ])
        C_tree = LDPCCode(H_tree)
        @test girth(C_tree) == -1
    end

    @testset "Cycle Enumeration and Distributions" begin
        F = Oscar.Nemo.Native.GF(2)
        
        # Graph with two disconnected 4-cycles
        H = matrix(F, [
            1 1 0 0;
            1 1 0 0;
            0 0 1 1;
            0 0 1 1
        ])
        C = LDPCCode(H)
        
        @test girth(C) == 4
        
        # Simple Cycles
        cycles = enumerate_simple_cycles(C; len=4)
        @test length(cycles) == 2
        
        dist = simple_cycle_length_distribution(C; len=4)
        @test dist[4] == 2
        @test count_simple_cycles(C; len=4) == 2
        
        @test average_simple_cycle_length(C; len=4) == 4.0
        @test median_simple_cycle_length(C; len=4) == 4.0
        @test mode_simple_cycle_length(C; len=4) == 4
        
        var_dist = simple_cycle_distribution_by_variable_node(C; len=4)
        @test var_dist[1] == 1 # v1 is involved in exactly one 4-cycle
        @test var_dist[3] == 1 # v3 is involved in exactly one 4-cycle
        
        # Short Cycles (For g=4, short cycles are defined as length g to 2g - 2 = 4 to 6)
        short_dist = short_cycle_length_distribution(C)
        @test short_dist[4] == 2
        @test count_short_cycles(C) == 2
        @test average_short_cycle_length(C) == 4.0
        @test median_short_cycle_length(C) == 4.0
        @test mode_short_cycle_length(C) == 4
        
        var_short_dist = short_cycle_distribution_by_variable_node(C)
        @test var_short_dist[2] == 1
    end

    @testset "Approximate Cycle EMD (ACE) Metrics" begin
        F = Oscar.Nemo.Native.GF(2)
        # Create a 4-cycle where the variable nodes have degree > 2 
        # so their ACE contribution (deg - 2) is positive.
        H_ace = matrix(F, [
            1 1 1 0;
            1 1 0 1;
            1 0 0 0;
            0 1 0 0
        ])
        C_ace = LDPCCode(H_ace)
        
        # v1 has degree 3 -> ACE contribution = 3 - 2 = 1
        # v2 has degree 3 -> ACE contribution = 3 - 2 = 1
        # The 4-cycle is on cols 1,2 and rows 1,2. Total ACE = 1 + 1 = 2.
        
        lens, ace_dists = CodingTheory._compute_ACE_distributions(C_ace)
        @test lens[1] == 4
        @test lens[2] == 4
        
        spectrum = ACE_spectrum(C_ace)
        @test haskey(spectrum, 4)
        
        # Test Statistical Arrays
        v1_aces = ACE_distribution(C_ace, 1)
        @test length(v1_aces) >= 1
        
        @test average_ACE_distribution(C_ace, 1) >= 0.0
        @test median_ACE_distribution(C_ace, 1) >= 0.0
        @test mode_ACE_distribution(C_ace, 1) >= 0
        
        # Check vectorized getters
        @test length(average_ACE_distribution(C_ace)) == C_ace.n
    end

    @testset "Cycle Removal (Socket-Swapping)" begin
        F = Oscar.Nemo.Native.GF(2)
        
        # Small bipartite graph with a known 4-cycle
        H = matrix(F, [
            1 1 1 0 0 0;
            1 1 0 1 0 0;
            0 0 1 0 1 1;
            0 0 0 1 1 1
        ])
        C = LDPCCode(H)
        @test girth(C) == 4
        
        # Execute BFS cycle removal targeting girth 6
        C_new = remove_cycles(C, 6; max_iters=200)
        
        # The girth should now be at least 6
        @test girth(C_new) >= 6
        
        # The algorithm strictly preserves the number of edges and dimensions
        @test density(C_new) == density(C)
        @test size(parity_check_matrix(C_new)) == size(parity_check_matrix(C))
        
        # Removing cycles from a graph that already meets the target should immediately return
        C_skip = remove_cycles(C_new, 4)
        @test girth(C_skip) >= 6
    end
end