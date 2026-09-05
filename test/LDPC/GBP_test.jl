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

    # These pin the three properties that a correct region-based decoder must
    # have. An earlier engine passed every test above while failing all three: it
    # had no path for information to reach an outer region, so outer-region
    # beliefs were frozen at initialization, the iteration count did nothing, and
    # the result was one round of independent per-check decoding.
    @testset "Bethe region graph reduces to ordinary BP" begin
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(F, [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ])
        n = 7

        R = CodingTheory.bethe_region_graph(H)
        @test CodingTheory.is_valid_region_graph(R)
        # One region per check plus one per variable.
        @test length(CodingTheory.regions(R)) == 3 + n
        # Variable regions are singletons and carry counting number 1 - deg(v).
        singles = [r for r in CodingTheory.regions(R) if length(CodingTheory.id(r)) == 1]
        @test length(singles) == n
        @test all(r -> CodingTheory.overcounting_number(r) == 1 - length(CodingTheory.parents(r)),
                  singles)

        # With damping off, GBP on this graph is ordinary sum-product BP, so the
        # posteriors must agree to machine precision. Only runs where neither
        # decoder stopped early are compared, since otherwise the two would have
        # done different amounts of work.
        llrs = Float64[2.5, 1.0, 3.0, 0.7, 2.0, 1.5, 2.2]
        compared = 0
        worst = 0.0
        for (a, b) in ((1, 2), (2, 5), (3, 6), (1, 7))
            s = zeros(UInt8, 3)
            for c in 1:3
                s[c] = UInt8((Int(lift(ZZ, H[c, a])) + Int(lift(ZZ, H[c, b]))) % 2)
            end
            for it in 1:3
                W_bp = CodingTheory.init_soft_workspace(H)
                conv, _ = CodingTheory.decode!(W_bp, llrs; algorithm = :sum_product,
                                               syndrome = Int.(s), max_iter = it)
                conv && continue
                post = copy(W_bp.total_llrs)

                W_g = CodingTheory.init_gbp_workspace(R, H)
                ok, _, _ = CodingTheory.gbp_decode!(W_g, R, H, llrs; target_syndrome = s,
                                                    max_iter = it, damping = 1.0)
                ok && continue
                got = CodingTheory.gbp_marginal_llrs(W_g)

                compared += 1
                worst = max(worst, maximum(abs.(post .- got)))
            end
        end
        @test compared > 0          # guard against a vacuous pass
        @test worst < 1e-9
    end

    @testset "A single all-variable region gives exact marginals" begin
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(F, [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ])
        n = 7
        llrs = Float64[2.9, 2.2, 1.6, 2.5, 2.0, 1.2, 2.6]

        R = CodingTheory.region_graph_from_base_nodes([collect(1:n)])
        for (a, b) in ((2, 3), (1, 5), (4, 7))
            s = zeros(UInt8, 3)
            for c in 1:3
                s[c] = UInt8((Int(lift(ZZ, H[c, a])) + Int(lift(ZZ, H[c, b]))) % 2)
            end

            W = CodingTheory.init_gbp_workspace(R, H)
            CodingTheory.init_region_beliefs!(W, R, H, llrs, s)
            got = CodingTheory.extract_hard_decisions(W, R, n)

            # Brute-force bitwise MAP over the coset {e : H e = s}.
            lse(x, y) = x == -Inf ? y : (y == -Inf ? x : max(x, y) + log1p(exp(-abs(x - y))))
            l0 = fill(-Inf, n); l1 = fill(-Inf, n)
            for mask in 0:(2^n - 1)
                e = [(mask >> (j - 1)) & 1 for j in 1:n]
                ok = true
                for c in 1:3
                    par = 0
                    for v in 1:n
                        iszero(H[c, v]) || (par ⊻= e[v])
                    end
                    par == s[c] || (ok = false; break)
                end
                ok || continue
                w = -sum(llrs[j] for j in 1:n if e[j] == 1; init = 0.0)
                for j in 1:n
                    e[j] == 1 ? (l1[j] = lse(l1[j], w)) : (l0[j] = lse(l0[j], w))
                end
            end
            want = UInt8[l1[j] > l0[j] ? 0x01 : 0x00 for j in 1:n]
            @test got == want
        end
    end

    @testset "Outer-region beliefs move and iterations matter" begin
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(F, [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ])
        R = CodingTheory.canonical_region_graph(H)
        bases = [i for (i, r) in enumerate(CodingTheory.regions(R))
                 if isempty(CodingTheory.parents(r))]
        @test !isempty(bases)

        llrs = Float64[0.4, 0.5, 0.6, 0.45, 0.55, 0.5, 0.6]
        s = UInt8[1, 1, 0]

        W0 = CodingTheory.init_gbp_workspace(R, H)
        CodingTheory.init_region_beliefs!(W0, R, H, llrs, s)
        before = copy(W0.log_beliefs)

        W = CodingTheory.init_gbp_workspace(R, H)
        CodingTheory.gbp_decode!(W, R, H, llrs; target_syndrome = s, max_iter = 25,
                                 damping = 0.5)

        drift = 0.0
        for i in bases
            lo = W.log_belief_offsets[i]
            hi = W.log_belief_offsets[i + 1] - 1
            for k in lo:hi
                isfinite(before[k]) && isfinite(W.log_beliefs[k]) &&
                    (drift = max(drift, abs(before[k] - W.log_beliefs[k])))
            end
        end
        # An outer region must hear from the others through its descendants.
        @test drift > 1e-6
    end

    @testset "Marginal modes agree at a consistent fixed point" begin
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(F, [1 1 0 0; 0 1 1 0; 0 0 1 1])
        R = CodingTheory.bethe_region_graph(H)
        llrs = Float64[1.4, 1.1, 1.7, 1.2]
        s = UInt8[1, 0, 0]

        W = CodingTheory.init_gbp_workspace(R, H)
        ok, _, _ = CodingTheory.gbp_decode!(W, R, H, llrs; target_syndrome = s,
                                            max_iter = 400, damping = 0.5)
        a = CodingTheory.gbp_marginal_llrs(W; mode = :smallest)
        b = CodingTheory.gbp_marginal_llrs(W; mode = :counting)
        @test length(a) == 4 && length(b) == 4
        # Signs must agree even where the magnitudes differ.
        @test all(sign.(a) .== sign.(b))
        @test_throws ArgumentError CodingTheory.gbp_marginal_llrs(W; mode = :nonsense)
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