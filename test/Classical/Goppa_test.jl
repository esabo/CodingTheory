@testitem "Classical/Goppa.jl" begin
    using Oscar, CodingTheory

    @testset "Goppa Code" begin
        # Ling & Xing, Example 9.3.10 (iii), p. 200
        E = GF(8, :α)
        S, z = polynomial_ring(E, :z)
        α = gen(E)
        
        # Polynomial chosen to ensure no roots fall into our manual support L
        g = z^2 + z + α^3
        # Support avoids elements that evaluate to 0 under g
        L = [E(0); [α^i for i in 0:6 if !iszero(g(α^i))]]
        F = GF(2)
        
        C = GoppaCode(F, L, g)
        @test length(C) == length(L)
        @test dimension(C) > 0
        @test minimum_distance(C)[1] >= degree(g) + 1

        # MacWilliams & Sloane, p. 342
        # Setup using an irreducible polynomial over the field
        g_irr = z^2 + z + 1
        L_all = [E(0); [α^i for i in 0:6]] # Safe because g_irr has no roots in GF(8)
        C_mac = GoppaCode(F, L_all, g_irr)
        
        @test is_irreducible(C_mac)
        @test is_separable(C_mac)
        @test !is_cumulative(C_mac)
        @test nonzeros(C_mac) == L_all

        # Extension and Permutation tests
        C_ext = extend(C_mac)
        @test length(C_ext) == length(L_all) + 1
        
        C_ext_perm = permute_code(C_ext, collect(1:length(C_ext)))
        flag, C_cyc = is_cyclic(C_ext_perm)
        @test flag isa Bool

        # MacWilliams & Sloane, p. 343
        E5 = GF(2^5, :α)
        S5, z5 = polynomial_ring(E5, :z)
        α5 = gen(E5)
        L5 = [E5(0); [α5^i for i in 0:Int(order(E5)) - 2]]
        g5 = z5^3 + z5 + 1
        C5 = GoppaCode(GF(2), L5, g5)
        @test is_irreducible(C5)
        @test length(C5) == 32
        @test dimension(C5) == 17
        @test minimum_distance(C5)[1] >= 7

        # MacWilliams & Sloane, p. 344
        E4 = GF(2^4, :α)
        S4, z4 = polynomial_ring(E4, :z)
        α4 = gen(E4)
        g4 = z4^2 + z4 + α4^3
        L4 = [E4(0); [α4^i for i in 0:Int(order(E4)) - 2]]
        C4 = GoppaCode(GF(2), L4, g4)
        @test is_irreducible(C4)
        @test length(C4) == 16
        @test dimension(C4) == 8
        @test minimum_distance(C4)[1] >= 5

        # ==============================================================================
        # CRYPTOGRAPHIC / RANDOM GENERATOR TESTS
        # ==============================================================================
        @testset "Random Goppa Generation" begin
            F_base = GF(2)
            E_ext = GF(2^4, :α)
            n_len = 12
            t_deg = 2
            
            # Verify code constructs without an infinite loop or domain errors
            C_rand = RandomGoppaCode(F_base, E_ext, n_len, t_deg)
            @test length(C_rand) == n_len
            @test extension_field(C_rand) == E_ext
            @test degree(Goppa_polynomial(C_rand)) == t_deg
            
            # Error checking boundaries
            @test_throws DomainError RandomGoppaCode(F_base, E_ext, 100, t_deg) # n > |E|
            @test_throws DomainError RandomGoppaCode(F_base, E_ext, n_len, 0)    # t < 1
            @test_throws ArgumentError RandomGoppaCode(E_ext, F_base, n_len, t_deg) # Wrong extension order
        end
    end
end
