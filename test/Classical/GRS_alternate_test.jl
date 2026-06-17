@testitem "Classical/GRS_alternate.jl" begin
    using Oscar, CodingTheory

    @testset "GRS And Alternate Codes" begin
        # MacWilliams & Sloane, p. 335
        E = GF(8, :α)
        α = gen(E)
        γ = [α^i for i in 0:6]
        v = [E(1) for _ in 1:7]
        A = AlternateCode(GF(2), 2, v, γ)
        @test length(A) == 7
        @test dimension(A) == 3
        @test minimum_distance(A)[1] == 4

        v = γ
        A2 = AlternateCode(GF(2), 2, v, γ)
        @test length(A2) == 7
        @test dimension(A2) == 4
        @test minimum_distance(A2)[1] == 3

        # Ling & Xing, Example 9.2.4 (i), p. 193
        # To yield exactly the Hamming code, we use r = 1 and v = γ. 
        # This expands the single row H = [γ_1, ..., γ_n] into the full m x n parity check matrix.
        m_val = 3
        E_ham = GF(2^m_val, :α)
        α_ham = gen(E_ham)
        F_base = GF(2)
        n_ham = 2^m_val - 1
        γ_ham = [α_ham^i for i in 0:(n_ham - 1)]
        
        # r = 1, v = γ
        A_ham = AlternateCode(F_base, 1, γ_ham, γ_ham)
        
        @test length(A_ham) == 7
        @test dimension(A_ham) == 4
        @test minimum_distance(A_ham)[1] == 3

        # Ling & Xing, Example 9.2.4 (iii), p. 194
        E3 = GF(2^3, :α)
        n3 = 6
        α3 = gen(E3)
        v3 = [E3(1) for _ in 1:n3]
        γ3 = [α3, α3^2, α3^3, α3^4, α3^5, α3^6]
        A3 = AlternateCode(GF(2), 3, v3, γ3)
        @test length(A3) == 6
        @test dimension(A3) == 2
        @test minimum_distance(A3)[1] == 4
        
        @testset "GRS Casting (Goppa and Alternate)" begin
            E_cast = GF(8, :α)
            S_cast, z_cast = polynomial_ring(E_cast, :z)
            α_cast = gen(E_cast)
            F_cast = GF(2)
            
            # Setup a Goppa Code
            g_cast = z_cast^2 + z_cast + α_cast^3
            L_cast = [E_cast(0); [α_cast^i for i in 0:6 if !iszero(g_cast(α_cast^i))]]
            Γ = GoppaCode(F_cast, L_cast, g_cast)
            
            # Cast Goppa to GRS
            GRS_Γ = GeneralizedReedSolomonCode(Γ)
            @test length(GRS_Γ) == length(L_cast)
            @test dimension(GRS_Γ) == length(L_cast) - degree(g_cast)
            @test evaluation_points(GRS_Γ) == L_cast
            
            # Setup an Alternate Code
            γ_alt = [α_cast^i for i in 0:6]
            v_alt = [E_cast(1) for _ in 1:7]
            A_cast = AlternateCode(F_cast, 2, v_alt, γ_alt)
            
            # Cast Alternate to GRS
            GRS_A = GeneralizedReedSolomonCode(A_cast)
            @test length(GRS_A) == 7
        end
    end

    @testset "Srivastava codes" begin
        # MacWilliams & Sloane, Example, p. 358
        E = GF(2^6, :α)
        α = gen(E)
        a = [E(0), E(1), α^9, α^18, α^27, α^36, α^45, α^54]
        w = [α]
        z = [E(1) for _ in 1:8]
        F = GF(2)
        C = GeneralizedSrivastavaCode(F, a, w, z, 2)
        @test Int(order(extension_field(C))) == 64
        @test length(C) == 8
        @test dimension(C) == 2
        @test minimum_distance(C)[1] == 5
        
        # from Goppa_test.jl
        E2 = GF(8, :α)
        S, z2 = polynomial_ring(E2, :z)
        β = gen(E2)
        g = β^3 + z2 + z2^2
        L = [E2(0); [β^i for i in 0:6 if !iszero(g(β^i))]]
        C2 = GoppaCode(F, L, g)
        
        # broken because the perm function has been removed
        # flag, _ = are_permutation_equivalent(C, C2)
        # @test_broken flag

        # MacWilliams & Sloane, Problem (15), p. 359
        E4 = GF(2^4, :α)
        α4 = gen(E4)
        w4 = [E4(0), E4(1)]
        a4 = setdiff(collect(E4), w4)
        z4 = [E4(1) for _ in 1:length(a4)]
        C4 = GeneralizedSrivastavaCode(F, a4, w4, z4, 2)
        @test length(C4) == 14
        @test dimension(C4) == 6
        @test minimum_distance(C4)[1] == 5
        D4 = dual(C4)
        @test minimum_distance(D4)[1] == 4

        # MacWilliams & Sloane, Problem (16), p. 359
        E5 = GF(2^4, :α)
        α5 = gen(E5)
        w5 = [α5^-1, α5^-3]
        a5 = setdiff(collect(E5), [E5(0); w5])
        z5 = a5
        C5 = GeneralizedSrivastavaCode(F, a5, w5, z5, 2)
        @test length(C5) == 13
        @test dimension(C5) == 5
        @test minimum_distance(C5)[1] == 5
        D5 = dual(C5)
        # Problem 16 mathematically states that C^perp has d = 3
        @test minimum_distance(D5)[1] == 3
    end
end
