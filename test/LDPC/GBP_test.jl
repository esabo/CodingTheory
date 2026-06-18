@testitem "LDPC/GBP.jl" begin
    using Oscar, CodingTheory

    @testset "Region Graph Construction and Validity" begin
        F = Oscar.Nemo.Native.GF(2)
        # Using a small parity check matrix (e.g., Hamming(7,4))
        H = matrix(F, [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ])
        
        # Build canonical region graph
        R = CodingTheory.canonical_region_graph(H)
        
        @test R isa CodingTheory.RegionGraph
        @test length(CodingTheory.regions(R)) > 0
        
        # A valid region graph MUST have exactly a sum of 1 for the overcounting 
        # numbers of all regions containing any specific variable
        @test CodingTheory.is_valid_region_graph(R)
        
        # Test basic getters
        reg = CodingTheory.regions(R)[1]
        @test !isempty(CodingTheory.id(reg))
        @test typeof(CodingTheory.overcounting_number(reg)) == Int
        
        # Test topological views
        @test length(collect(CodingTheory.base_regions(R))) > 0
        @test length(collect(CodingTheory.leaves(R))) > 0
    end

    @testset "Region Graph Topology Reductions" begin
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(F, [
            1 1 1 0 0 0;
            0 1 1 1 0 0;
            0 0 1 1 1 0;
            0 0 0 1 1 1
        ])
        
        R = CodingTheory.canonical_region_graph(H)
        orig_len = length(CodingTheory.regions(R))
        
        # 1. Remove Zero-Overcounting Regions
        R_no_zeros = CodingTheory.remove_zero_overcounting_numbers(R)
        @test CodingTheory.is_valid_region_graph(R_no_zeros)
        @test length(CodingTheory.regions(R_no_zeros)) <= orig_len
        
        # 2. Remove Generational Skips
        R_no_skips = CodingTheory.remove_generational_skips(R_no_zeros)
        @test CodingTheory.is_valid_region_graph(R_no_skips)
        
        # The topological schedule should still safely cover all remaining regions
        order = CodingTheory.message_passing_order(R_no_skips)
        @test length(order) == length(CodingTheory.regions(R_no_skips))
        
        # 3. Triangulation (Chordalization) for generalized base regions
        cliques = CodingTheory.triangulate_base_regions(H)
        @test length(cliques) > 0
        @test all(c -> typeof(c) == BitSet, cliques)
    end

    @testset "GBP Workspace and Message Passing" begin
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(F, [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ])
        
        R = CodingTheory.canonical_region_graph(H)
        R_opt = CodingTheory.remove_generational_skips(CodingTheory.remove_zero_overcounting_numbers(R))
        
        # Initialize the zero-allocation GBP workspace
        W = CodingTheory.init_gbp_workspace(R_opt, H)
        
        @test W.num_regions == length(CodingTheory.regions(R_opt))
        @test W.num_edges > 0
        @test length(W.edge_update_order) == 2 * W.num_edges # Full upward and downward sweeps
        
        # Test marginalized log-belief stride mapping sizes
        @test length(W.log_beliefs) == W.log_belief_offsets[end] - 1
    end

    @testset "GBP Decoding (Error Correction)" begin
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(F, [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ])
        
        R = CodingTheory.canonical_region_graph(H)
        W = CodingTheory.init_gbp_workspace(R, H)
        
        # Valid codeword: c = [1, 1, 1, 0, 0, 0, 0]
        # In BPSK: 0 -> +5.0, 1 -> -5.0
        # We inject a weak error at index 1 -> make it +1.0 (looks like a 0, but weak)
        llrs = Float64[1.0, -5.0, -5.0, 5.0, 5.0, 5.0, 5.0]
        expected_cw = UInt8[1, 1, 1, 0, 0, 0, 0]
        
        success, out_bits, iters = CodingTheory.gbp_decode!(W, R, H, llrs, max_iter=20, damping=0.5)
        
        @test success
        @test out_bits == expected_cw
        @test iters > 0
    end

    @testset "String Representations" begin
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(F, [1 1 0; 0 1 1])
        R = CodingTheory.canonical_region_graph(H)
        
        # Test Region
        reg = CodingTheory.regions(R)[1]
        out_reg = sprint(show, reg)
        @test contains(out_reg, "Region(id={")
        
        # Test RegionGraph (compact)
        out_rg = sprint(show, R)
        @test contains(out_rg, "RegionGraph(")
        
        # Test RegionGraph (MIME text/plain)
        out_rg_mime = sprint(show, MIME"text/plain"(), R)
        @test contains(out_rg_mime, "Base regions:")
        @test contains(out_rg_mime, "Leaf regions:")
    end
end