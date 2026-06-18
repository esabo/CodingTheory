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