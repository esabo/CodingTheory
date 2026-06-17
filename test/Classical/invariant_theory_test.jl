@testitem "Classical/invariant_theory.jl" begin
    using Oscar, CodingTheory

    @testset "Distance Bounds" begin
        # Mallows-Sloane Bounds
        # Type II doubly-even codes only exist for lengths divisible by 8[cite: 88].
        @test Mallows_Sloane_bound(8, type=:TypeII) == 4
        @test Mallows_Sloane_bound(24, type=:TypeII) == 8
        @test_throws ArgumentError Mallows_Sloane_bound(12, type=:TypeII)

        # Type I singly-even codes exist for all even lengths[cite: 89].
        @test Mallows_Sloane_bound(8, type=:TypeI) == 4
        @test Mallows_Sloane_bound(22, type=:TypeI) == 6
        @test Mallows_Sloane_bound(24, type=:TypeI) == 8
        @test_throws ArgumentError Mallows_Sloane_bound(7, type=:TypeI)

        # Bachoc-Gaborit Bounds
        # If n = 22 mod 24, the bound is n/2 + 6[cite: 90].
        @test Bachoc_Gaborit_bound(22) == 17
        @test Bachoc_Gaborit_bound(46) == 29
        
        # Otherwise, the bound is n/2 + 4[cite: 90].
        @test Bachoc_Gaborit_bound(20) == 14
        @test Bachoc_Gaborit_bound(24) == 16
        @test_throws ArgumentError Bachoc_Gaborit_bound(9)
    end

    @testset "Gleason Generators" begin
        ϕ_2, ϕ_8, ϕ_24 = gleason_generators()
        
        # Verify the degrees and basic structures of the generators[cite: 91].
        @test total_degree(ϕ_2) == 2
        @test total_degree(ϕ_8) == 8
        @test total_degree(ϕ_24) == 24
        
        # ϕ_8 is the extended Hamming code generator[cite: 91].
        R_QQ, (x, y) = polynomial_ring(QQ, ["x", "y"])
        @test ϕ_8 == x^8 + 14 * x^4 * y^2 + y^8
    end

    @testset "Extremal Weight Enumerators & Sanity Checker" begin
        R_QQ, (x, y) = polynomial_ring(QQ, ["x", "y"])
        
        # Generate the extremal weight enumerator for the [24, 12, 8] Golay Code
        # Mathematically forces low-weight coefficients to 0[cite: 110, 114].
        W_24 = extremal_weight_enumerator(24, type=:TypeII)
        
        # The true Golay weight enumerator
        W_golay = x^24 + 759 * x^16 * y^8 + 2576 * x^12 * y^12 + 759 * x^8 * y^16 + y^24
        @test W_24 == W_golay
        
        # Validate the generated enumerator using the Sloane Sanity Checker
        # Checks if W is strictly in the Type II Gleason invariant ring[cite: 99, 100].
        is_valid, coeffs = is_valid_self_dual_enumerator(W_24, 24, type=:TypeII)
        @test is_valid == true
        # The basis coefficients must be integers[cite: 109].
        @test all(is_integer, coeffs)
        
        # Type I Extremal check
        W_8 = extremal_weight_enumerator(8, type=:TypeI)
        is_valid_8, _ = is_valid_self_dual_enumerator(W_8, 8, type=:TypeI)
        @test is_valid_8 == true
        
        # Ghost Code Rejection
        # Mutate the Golay enumerator so it is no longer algebraically valid
        W_ghost = W_golay + 5 * x^22 * y^2
        is_valid_ghost, ghost_coeffs = is_valid_self_dual_enumerator(W_ghost, 24, type=:TypeII)
        @test is_valid_ghost == false
        @test isempty(ghost_coeffs)
    end

    @testset "Shadow Transform" begin
        R_QQ, (x, y) = polynomial_ring(QQ, ["x", "y"])
        
        # Use the [8, 4, 4] extended Hamming code as our Type I base
        W_8 = x^8 + 14 * x^4 * y^4 + y^8
        
        # Transform the primal weight enumerator into its shadow W_S(x, y)[cite: 92, 95].
        W_shadow = shadow_transform(W_8, 8)
        
        # The shadow transform projects back to QQ[x, y] by extracting the real part[cite: 96, 98].
        @test parent(W_shadow) == R_QQ
        
        # Physical shadow weights must be rational/integer[cite: 96].
        # Verify that evaluating the shadow enumerator at x=1, y=1 gives 2^(n/2) = 16
        # (Total number of words in the shadow of an [8,4] code)
        shadow_eval = evaluate(W_shadow, [QQ(1), QQ(1)])
        @test shadow_eval == QQ(16)
        
        @test_throws ArgumentError shadow_transform(W_8, 7)
    end
end
