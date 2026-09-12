@testitem "LDPC/analysis.jl" begin
    using Oscar, CodingTheory

    @testset "LDPC Ensemble Polynomials and Bounds" begin
        # Generate generic degree distribution polynomials from vectors
        # lambda(x) = 0.5*x + 0.5*x^2
        # rho(x) = x^3
        _, x = polynomial_ring(RealField(), :x)
        λ = 0.5 * x + 0.5 * x^2
        ρ = x^3
        
        E = LDPCEnsemble(λ, ρ)
        
        # Test average degrees and basic design rate calculations
        @test E.design_rate > 0.0
        
        # Multiplicative gap lower bound should strictly be between 0 and 1
        gap_lb = multiplicative_gap_lower_bound(E)
        @test 0.0 < gap_lb < 1.0
        
        # Check density lower bound on a BEC
        Ch_BEC = BinaryErasureChannel(0.5)
        d_lb = density_lower_bound(Ch_BEC, 0.8)
        @test d_lb > 0.0
        
        # Ensure domain constraint is caught
        @test_throws DomainError density_lower_bound(Ch_BEC, 1.5)
    end

    @testset "Density Evolution and Thresholds" begin
        _, x = polynomial_ring(RealField(), :x)
        # Use a classic (3,6) regular LDPC degree distribution
        λ_reg = x^2
        ρ_reg = x^5
        E_reg = LDPCEnsemble(λ_reg, ρ_reg)
        
        # Density evolution over BEC
        Ch_BEC = BinaryErasureChannel(0.4)
        evo_x, evo_y = density_evolution(E_reg, Ch_BEC)
        @test length(evo_x) > 1
        @test length(evo_y) == length(evo_x)
        
        # Test BEC analytical threshold
        thresh_bec = optimal_threshold(E_reg, BinaryErasureChannel)
        @test 0.0 < thresh_bec < 1.0
        
        # Multiplicative Gap mapping
        gap = multiplicative_gap(E_reg, Ch_BEC)
        @test gap >= 0.0 # Must be a non-negative gap
        
        # Test BAWGN Gaussian Approximation threshold
        thresh_awgn = optimal_threshold(E_reg, BAWGNChannel)
        @test thresh_awgn > 0.0
    end

    @testset "PEXIT Analysis on Protographs" begin
        # Example protograph base matrix
        B = [1 2; 2 1]
        
        # The AWGN PEXIT analysis takes the standard deviation vector of the *LLRs*.
        # A good channel (low noise) means HIGH LLR standard deviation (e.g., 4.0)
        sigma_ch_good = [4.0, 4.0]
        decoded_good, I_APP_good, iters_good = CodingTheory._PEXIT_AWGN(B, sigma_ch_good)
        @test decoded_good
        @test all(I_APP_good .> 0.999)
        @test iters_good < 100
        
        # A bad channel (high noise) means LOW LLR standard deviation (e.g., 0.4)
        sigma_ch_bad = [0.4, 0.4]
        decoded_bad, I_APP_bad, iters_bad = CodingTheory._PEXIT_AWGN(B, sigma_ch_bad)
        @test !decoded_bad
        @test any(I_APP_bad .< 0.999)
        
        # Test the threshold wrapper
        thresh_p = protograph_threshold(B)
        @test thresh_p > 0.0
        
        # Test specific column puncturing threshold (which weakens the code)
        thresh_punc = protograph_threshold(B, punctured=[true, false])
        @test thresh_punc < thresh_p
        
        # Catch array length mismatch
        @test_throws ArgumentError protograph_threshold(B, punctured=[true])
    end

    @testset "EXIT Chart Data Generation" begin
        _, x = polynomial_ring(RealField(), :x)
        λ = x^2
        ρ = x^5
        E = LDPCEnsemble(λ, ρ)
        
        Ch_BEC = BinaryErasureChannel(0.4)
        Ch_AWGN = BAWGNChannel(0.8)
        
        # Test BEC analytical curves
        vnd_x_b, vnd_y_b, cnd_x_b, cnd_y_b = EXIT_chart_data(E, Ch_BEC; pts=10)
        @test length(vnd_x_b) == 10
        @test vnd_x_b[1] == 0.0 && vnd_x_b[end] == 1.0 # MI limits
        
        # Test AWGN Gaussian Approximation curves
        vnd_x_a, vnd_y_a, cnd_x_a, cnd_y_a = EXIT_chart_data(E, Ch_AWGN; pts=10)
        @test length(vnd_x_a) == 10
        @test cnd_y_a[1] == 0.0 && cnd_y_a[end] == 1.0
        
        # Test Protograph averaged EXIT curves
        B = [2 1 0; 0 1 2]
        sigma_ch = [0.8, 0.8, 0.8]
        vnd_x_p, vnd_y_p, cnd_x_p, cnd_y_p = PEXIT_chart_data(B, sigma_ch; pts=10)
        @test length(vnd_x_p) == 10
        # A fully certain prior (I_A = 1.0) generally yields a perfect check transfer (I_E = 1.0)
        @test cnd_x_p[end] ≈ 1.0 atol=1e-3
    end
end