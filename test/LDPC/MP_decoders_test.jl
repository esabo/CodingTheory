@testitem "LDPC/MP_decoders.jl" begin
    using Oscar, CodingTheory

    @testset "Workspace Initialization & Layered Schedules" begin
        F = Oscar.Nemo.Native.GF(2)
        # Standard Hamming(7,4) parity check matrix
        H = matrix(F, [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ])
        
        # Test Layered Schedule Generation
        sched_layered = CodingTheory.layered_schedule(H, schedule=:layered)
        @test length(sched_layered) > 0
        @test sum(length(layer) for layer in sched_layered) == 3 # All 3 rows accounted for
        
        sched_parallel = CodingTheory.layered_schedule(H, schedule=:parallel)
        @test length(sched_parallel) == 1
        @test length(sched_parallel[1]) == 3
        
        # Test Workspace allocations
        W_hard = CodingTheory.init_hard_workspace(H)
        @test W_hard.num_check == 3
        @test W_hard.num_var == 7
        @test W_hard.num_edges == 12 # Total ones in H
        
        W_soft = CodingTheory.init_soft_workspace(H, schedule=:layered)
        @test W_soft.num_edges == 12
        @test length(W_soft.layers) > 0
    end

    @testset "Hard Decision Decoding (Gallager)" begin
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(F, [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ])
        W_hard = CodingTheory.init_hard_workspace(H)
        
        # Valid codeword: c = [1, 1, 1, 0, 0, 0, 0]
        # Introduce an error at index 1 -> [0, 1, 1, 0, 0, 0, 0]
        received_bits = UInt8[0, 1, 1, 0, 0, 0, 0]
        expected_bits = UInt8[1, 1, 1, 0, 0, 0, 0]
        
        # Decode using Gallager logic (Bt = 2)
        success, out_bits, iters = CodingTheory.decode!(W_hard, received_bits, Bt=2)
        
        @test success
        @test out_bits == expected_bits
        @test iters > 0
    end

    @testset "Soft Decision Decoding (Belief Propagation)" begin
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(F, [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ])
        
        # Valid codeword: c = [1, 1, 1, 0, 0, 0, 0]
        # In BPSK: 0 -> +5.0, 1 -> -5.0 (High confidence)
        # We flip index 1 to +5.0 to simulate an error
        llrs = Float64[5.0, -5.0, -5.0, 5.0, 5.0, 5.0, 5.0]
        expected_bits = UInt8[1, 1, 1, 0, 0, 0, 0]
        
        # 1. Test Sum-Product (Exact Box-Plus)
        W_sp = CodingTheory.init_soft_workspace(H)
        success_sp, out_sp, _ = CodingTheory.decode!(W_sp, llrs, algorithm=:sum_product)
        @test success_sp
        @test out_sp == expected_bits
        
        # 2. Test Min-Sum
        W_ms = CodingTheory.init_soft_workspace(H)
        success_ms, out_ms, _ = CodingTheory.decode!(W_ms, llrs, algorithm=:min_sum)
        @test success_ms
        @test out_ms == expected_bits
        
        # 3. Test Offset Min-Sum with Layered Schedule
        W_oms = CodingTheory.init_soft_workspace(H, schedule=:layered)
        success_oms, out_oms, _ = CodingTheory.decode!(W_oms, llrs, algorithm=:offset_min_sum, schedule=:layered, offset=0.25)
        @test success_oms
        @test out_oms == expected_bits
    end

    @testset "Erasures and Decimation Hooks" begin
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(F, [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ])
        W = CodingTheory.init_soft_workspace(H)
        
        # Valid codeword: c = [1, 1, 1, 0, 0, 0, 0]
        # Send perfect LLRs, but erase index 1 and 2
        llrs = Float64[-5.0, -5.0, -5.0, 5.0, 5.0, 5.0, 5.0]
        expected_bits = UInt8[1, 1, 1, 0, 0, 0, 0]
        
        # Decode with erasures
        success, out, _ = CodingTheory.decode!(W, llrs, algorithm=:sum_product, erasures=[1, 2])
        @test success
        @test out == expected_bits
        
        # Decode with manual decimation (pinning bit 4 to 0)
        success_dec, out_dec, _ = CodingTheory.decode!(W, llrs, algorithm=:min_sum, decimated_bits_values=[(4, 0)])
        @test success_dec
        @test out_dec == expected_bits
        
        # Verify the decimation actually pinned the LLR
        @test W.is_decimated[4] == true
        @test W.channel_llrs[4] == 1000.0 # High confidence 0
    end

    @testset "Syndrome Decoding" begin
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(F, [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ])
        W = CodingTheory.init_soft_workspace(H)
        
        # Error vector: e = [1, 0, 0, 0, 0, 0, 0]
        # Syndrome: s = H * e = [1, 1, 0]
        expected_error = UInt8[1, 0, 0, 0, 0, 0, 0]
        target_syn = UInt8[1, 1, 0]
        
        # For syndrome decoding, the channel input is totally neutral (0.0)
        neutral_llrs = zeros(Float64, 7)
        
        success, out_error, _ = CodingTheory.decode!(W, neutral_llrs, algorithm=:min_sum_correction, syndrome=target_syn)
        
        @test success
        @test out_error == expected_error
    end
end