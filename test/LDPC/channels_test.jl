@testitem "LDPC/channels.jl" begin
    using Oscar, CodingTheory, Random

    @testset "Channel Getters" begin
        bec = BinaryErasureChannel(0.3)
        @test erasure_probability(bec) == 0.3
        
        bsc = BinarySymmetricChannel(0.1)
        @test crossover_probability(bsc) == 0.1
        
        bawgn = BAWGNChannel(2.0)
        @test standard_deviation(bawgn) == 2.0
        @test variance(bawgn) == 4.0
    end

    @testset "Channel Capacities" begin
        # Binary Erasure Channel (C = 1 - ε)
        @test capacity(BinaryErasureChannel(0.0)) == 1.0
        @test capacity(BinaryErasureChannel(1.0)) == 0.0
        @test capacity(BinaryErasureChannel(0.25)) == 0.75
        
        # Binary Symmetric Channel (C = 1 - H(p))
        @test capacity(BinarySymmetricChannel(0.0)) == 1.0
        @test capacity(BinarySymmetricChannel(1.0)) == 1.0 # Entropy is 0, Capacity is 1
        @test capacity(BinarySymmetricChannel(0.5)) == 0.0
        
        # Z-Channel
        @test capacity(ZChannel(0.0)) == 1.0
        @test capacity(ZChannel(1.0)) == 0.0
        @test 0.0 < capacity(ZChannel(0.5)) < 1.0
        
        # Binary AWGN Channel
        # Very low noise should approach C = 1.0
        @test capacity(BAWGNChannel(0.1)) ≈ 1.0 atol=1e-3
        # Very high noise should crash capacity towards 0.0
        @test capacity(BAWGNChannel(10.0)) < 0.1
        
        # Rayleigh Fading Channel
        @test capacity(RayleighFadingChannel(0.1)) > 0.9
        @test capacity(RayleighFadingChannel(10.0)) < 0.1
    end

    @testset "Transmit & LLR Calculations" begin
        # Common message vector
        x = [0, 1, 0, 1, 0, 0, 1, 1]
        
        # 1. Binary Erasure Channel
        bec = BinaryErasureChannel(0.5)
        y_bec = transmit(bec, x)
        @test all(v ∈ [-1, 0, 1] for v in y_bec)
        
        llr_bec = CodingTheory.llr(bec, y_bec)
        for (val, l) in zip(y_bec, llr_bec)
            if val == -1
                @test l == 0.0
            elseif val == 0
                @test l == 1e6
            elseif val == 1
                @test l == -1e6
            end
        end
        
        # 2. Binary Symmetric Channel
        p_crossover = 0.1
        bsc = BinarySymmetricChannel(p_crossover)
        y_bsc = transmit(bsc, x)
        @test all(v ∈ [0, 1] for v in y_bsc)
        
        llr_bsc = CodingTheory.llr(bsc, y_bsc)
        expected_mag = log((1.0 - p_crossover) / p_crossover)
        @test all(abs(l) ≈ expected_mag for l in llr_bsc)
        
        # 3. Z-Channel
        zc = ZChannel(0.5)
        y_zc = transmit(zc, x)
        @test all(v ∈ [0, 1] for v in y_zc)
        # Z-Channel strict property: 0s are NEVER flipped to 1s
        @test all(y_zc[i] == 0 for i in 1:length(x) if x[i] == 0)
        
        llr_zc = CodingTheory.llr(zc, y_zc)
        expected_ll_0 = -log(0.5)
        @test all(l == expected_ll_0 || l == -Inf for l in llr_zc)
        
        # 4. Binary AWGN Channel
        sigma_awgn = 1.5
        bawgn = BAWGNChannel(sigma_awgn)
        y_bawgn = transmit(bawgn, x)
        @test length(y_bawgn) == length(x)
        @test typeof(y_bawgn) == Vector{Float64}
        
        llr_bawgn = CodingTheory.llr(bawgn, y_bawgn)
        # Verify analytical identity LLR = 2y / σ^2
        @test all(llr_bawgn .≈ 2.0 .* y_bawgn ./ (sigma_awgn^2))
        
        # 5. Rayleigh Fading Channel
        sigma_rfc = 1.0
        rfc = RayleighFadingChannel(sigma_rfc)
        y_rfc, a_rfc = transmit(rfc, x)
        
        @test length(y_rfc) == length(x)
        @test length(a_rfc) == length(x)
        
        llr_rfc = CodingTheory.llr(rfc, (y_rfc, a_rfc))
        # Verify analytical fading identity LLR = 2ay / σ^2
        @test all(llr_rfc .≈ 2.0 .* a_rfc .* y_rfc ./ (sigma_rfc^2))
    end
    
    @testset "Show Methods" begin
        # Minimal checks to ensure string interpolation in show() doesn't crash
        bec = BinaryErasureChannel(0.15)
        bawgn = BAWGNChannel(2.0)
        
        @test contains(sprint(show, bec), "Binary erasure channel")
        @test contains(sprint(show, bawgn), "additive white Gaussian noise")
    end
end