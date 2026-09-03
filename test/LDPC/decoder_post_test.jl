@testitem "LDPC/decoder_post.jl" begin
    using Oscar, CodingTheory

    @testset "Workspace Initializations" begin
        F = Oscar.Nemo.Native.GF(2)
        # Hamming(7,4) parity check matrix
        H = matrix(F, [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ])
        
        W_osd = CodingTheory.init_osd_workspace(H)
        @test W_osd.num_var == 7
        @test W_osd.num_check == 3
        @test size(W_osd.H_dense) == (3, 7)
        
        W_grand = CodingTheory.init_grand_workspace(H)
        @test W_grand.num_var == 7
        @test length(W_grand.var_to_checks) == 7
        
        # Every field has to be initialised on both constructor paths, so compare
        # the Flint method against the array method field by field rather than
        # spot-checking the ones we happen to remember.
        W_grand_plain = CodingTheory.init_grand_workspace(UInt8[
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ])
        for f in fieldnames(CodingTheory.GRANDWorkspace)
            @test getfield(W_grand, f) == getfield(W_grand_plain, f)
        end
        
        W_wbf = CodingTheory.init_wbf_workspace(H)
        @test W_wbf.num_check == 3
        @test length(W_wbf.var_to_checks) == 7
        @test length(W_wbf.chk_to_vars) == 3
    end

    @testset "Order Statistic Decoding (OSD)" begin
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(F, [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ])
        W = CodingTheory.init_osd_workspace(H)
        
        # Valid codeword: c = [1, 1, 1, 0, 0, 0, 0]
        # LLR mapping: < 0 is 1, > 0 is 0
        expected_cw = UInt8[1, 1, 1, 0, 0, 0, 0]
        
        # Scenario: Index 1 is actually a 1 (should be < 0), but the channel was noisy
        # so it's slightly > 0 (looks like a 0, but very weak/low reliability).
        llrs = Float64[1.0, -5.0, -5.0, 5.0, 5.0, 5.0, 5.0]
        
        # 1. Standard OSD (Order 0, 1, 2)
        cw_0 = CodingTheory.osd_decode!(W, llrs, method=:standard, order=0)
        @test cw_0 == expected_cw
        
        cw_1 = CodingTheory.osd_decode!(W, llrs, method=:standard, order=1)
        @test cw_1 == expected_cw
        
        cw_2 = CodingTheory.osd_decode!(W, llrs, method=:standard, order=2)
        @test cw_2 == expected_cw
        
        # 2. Combinatorial Sweep (OSD-CS)
        cw_cs = CodingTheory.osd_decode!(W, llrs, method=:cs, order=2, cs_lambda=4)
        @test cw_cs == expected_cw
    end

    @testset "Guessing Random Additive Noise Decoding (GRAND)" begin
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(F, [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ])
        W = CodingTheory.init_grand_workspace(H)
        
        # Scenario: Weak error on index 1 and 2
        # (True bits are 1, but LLRs are weakly positive)
        llrs = Float64[0.5, 0.8, -5.0, 5.0, 5.0, 5.0, 5.0]
        expected_cw = UInt8[1, 1, 1, 0, 0, 0, 0]
        
        # GRAND sweeps the least reliable bits and flips them to find a zero syndrome
        success, cw = CodingTheory.grand_decode!(W, llrs, max_weight=2, max_lrb=5)
        @test success
        @test cw == expected_cw
        
        # Test Syndrome Decoding: No LLR bias, just find the error pattern!
        neutral_llrs = zeros(Float64, 7) # All bits perfectly 0.0 reliability
        # Error on bit 1 -> target syndrome should be H[:, 1] = [1, 1, 0]
        target_syn = UInt8[1, 1, 0]
        expected_error = UInt8[1, 0, 0, 0, 0, 0, 0]
        
        success_syn, cw_syn = CodingTheory.grand_decode!(W, neutral_llrs, syndrome=target_syn, max_weight=1)
        @test success_syn
        @test cw_syn == expected_error
        
        # ------------------------------------------------------------------
        # Soft-cost ordering
        # ------------------------------------------------------------------
        num_check, num_var = size(H)
        H_bits = UInt8[iszero(H[c, v]) ? 0x00 : 0x01 for c in 1:num_check, v in 1:num_var]
        
        hard_decisions = l -> UInt8[x < 0.0 ? 0x01 : 0x00 for x in l]
        syndrome_of = cwd -> UInt8[reduce(⊻, H_bits[c, v] & cwd[v] for v in 1:num_var)
                                   for c in 1:num_check]
        soft_cost = (l, cwd) -> begin
            hard = hard_decisions(l)
            sum(abs(l[v]) for v in 1:num_var if cwd[v] != hard[v]; init = 0.0)
        end
        
        # Exhaustive minimum-soft-cost search over the same space, so the cases
        # below pin the maximum-likelihood property itself instead of one
        # hand-computed answer. Enumerating subsets by bitmask is fine at n = 7.
        brute_force_best = (l, target, max_weight) -> begin
            hard = hard_decisions(l)
            best_cost = Inf
            best_cw = nothing
            for mask in 0:(2^num_var - 1)
                count_ones(mask) > max_weight && continue
                cwd = copy(hard)
                cost = 0.0
                for v in 1:num_var
                    if isodd(mask >> (v - 1))
                        cwd[v] ⊻= 0x01
                        cost += abs(l[v])
                    end
                end
                if syndrome_of(cwd) == target && cost < best_cost
                    best_cost = cost
                    best_cw = cwd
                end
            end
            return best_cost, best_cw
        end
        
        # Soft cost, not Hamming weight, has to decide. The single flip of bit 4
        # matches this syndrome at cost 1.0, but flipping bits 1 and 7 matches it
        # at cost 0.5, so first-match-in-weight-order returns the wrong pattern.
        ordering_llrs = Float64[0.2, 5.0, 5.0, 1.0, 5.0, 5.0, 0.3]
        ordering_syn = UInt8[1, 1, 1]   # H[:, 4] == H[:, 1] ⊻ H[:, 7]
        success_ord, cw_ord = CodingTheory.grand_decode!(W, ordering_llrs,
                                                         syndrome=ordering_syn,
                                                         max_weight=2, max_lrb=num_var)
        @test success_ord
        @test cw_ord == UInt8[1, 0, 0, 0, 0, 0, 1]
        @test soft_cost(ordering_llrs, cw_ord) ≈ 0.5
        
        # A weight-3 pattern can win too: bits 5, 6 and 7 cost 0.6 against 1.0 for
        # the lone weight-1 match on bit 4. Capping max_weight still truncates the
        # search space, which is where the cheap weight-1 answer comes back.
        deep_llrs = Float64[9.0, 9.0, 9.0, 1.0, 0.1, 0.2, 0.3]
        success_3, cw_3 = CodingTheory.grand_decode!(W, deep_llrs, syndrome=UInt8[1, 1, 1],
                                                     max_weight=3, max_lrb=num_var)
        @test success_3
        @test cw_3 == UInt8[0, 0, 0, 0, 1, 1, 1]
        
        success_2, cw_2 = CodingTheory.grand_decode!(W, deep_llrs, syndrome=UInt8[1, 1, 1],
                                                     max_weight=2, max_lrb=num_var)
        @test success_2
        @test cw_2 == UInt8[0, 0, 0, 1, 0, 0, 0]
        
        # Every case must come back with the exhaustive minimum soft cost
        for (case_llrs, case_syn, case_weight) in (
                (Float64[0.5, 0.8, -5.0, 5.0, 5.0, 5.0, 5.0], zeros(UInt8, num_check), 2),
                (Float64[0.2, 5.0, 5.0, 1.0, 5.0, 5.0, 0.3],  UInt8[1, 1, 1],          2),
                (Float64[9.0, 9.0, 9.0, 1.0, 0.1, 0.2, 0.3],  UInt8[1, 1, 1],          3),
                (Float64[-0.4, 5.0, 5.0, -0.9, 5.0, 5.0, 5.0], zeros(UInt8, num_check), 2),
                (Float64[3.0, 0.7, 0.6, 0.5, 2.0, 1.0, 4.0],  UInt8[0, 1, 0],          3))
            ref_cost, ref_cw = brute_force_best(case_llrs, case_syn, case_weight)
            success_c, cw_c = CodingTheory.grand_decode!(W, case_llrs, syndrome=case_syn,
                                                         max_weight=case_weight,
                                                         max_lrb=num_var)
            @test success_c == (ref_cw !== nothing)
            @test syndrome_of(cw_c) == case_syn
            @test soft_cost(case_llrs, cw_c) ≈ ref_cost
        end
        
        # ------------------------------------------------------------------
        # Plain syndrome matching still works
        # ------------------------------------------------------------------
        
        # Base case: the hard decisions already satisfy the target syndrome
        success_base, cw_base = CodingTheory.grand_decode!(W, fill(5.0, num_var))
        @test success_base
        @test cw_base == zeros(UInt8, num_var)
        
        # No pattern in the search space matches. Bit 7 is the only single flip
        # with this syndrome and it sits outside the 3 least reliable bits, and no
        # pair drawn from those 3 reaches it either.
        narrow_llrs = Float64[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
        success_none, cw_none = CodingTheory.grand_decode!(W, narrow_llrs,
                                                           syndrome=UInt8[0, 0, 1],
                                                           max_weight=2, max_lrb=3)
        @test !success_none
        @test cw_none == zeros(UInt8, num_var)  # left at the hard decisions
        
        # Widening the search space finds it again
        success_wide, cw_wide = CodingTheory.grand_decode!(W, narrow_llrs,
                                                           syndrome=UInt8[0, 0, 1],
                                                           max_weight=1, max_lrb=num_var)
        @test success_wide
        @test cw_wide == UInt8[0, 0, 0, 0, 0, 0, 1]
        
        # ------------------------------------------------------------------
        # GRAND must agree with OSD wherever both reach the ML answer
        # ------------------------------------------------------------------
        W_osd = CodingTheory.init_osd_workspace(H)
        for agree_llrs in (Float64[1.0, -5.0, -5.0, 5.0, 5.0, 5.0, 5.0],
                           Float64[0.5, 0.8, -5.0, 5.0, 5.0, 5.0, 5.0],
                           Float64[-0.4, 5.0, 5.0, -0.9, 5.0, 5.0, 5.0])
            cw_osd = copy(CodingTheory.osd_decode!(W_osd, agree_llrs, method=:cs, order=2))
            success_g, cw_g = CodingTheory.grand_decode!(W, agree_llrs,
                                                         max_weight=3, max_lrb=num_var)
            @test success_g
            @test cw_g == cw_osd
        end
    end

    @testset "Residual Weighted Bit-Flipping (WBF)" begin
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(F, [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ])
        W = CodingTheory.init_wbf_workspace(H)
        
        # Scenario: One bit is in error, but we have LLR energies to guide the WBF
        llrs = Float64[1.5, -5.0, -5.0, 5.0, 5.0, 5.0, 5.0]
        expected_cw = UInt8[1, 1, 1, 0, 0, 0, 0]
        
        success, cw, iters = CodingTheory.wbf_decode!(W, llrs, alpha=0.5, max_iters=20)
        @test success
        @test cw == expected_cw
        @test iters > 0
        
        # Test Syndrome Decoding with WBF
        neutral_llrs = zeros(Float64, 7) 
        target_syn = UInt8[1, 1, 0] # Matches column 1
        expected_error = UInt8[1, 0, 0, 0, 0, 0, 0]
        
        success_syn, cw_syn, _ = CodingTheory.wbf_decode!(W, neutral_llrs, syndrome=target_syn, alpha=0.0)
        @test success_syn
        @test cw_syn == expected_error
    end
end