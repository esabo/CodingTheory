@testset "Classical/GeneralizedReedSolomon.jl" begin
    using Oscar, CodingTheory

    @testset "GRS And Alternate Codes" begin
        # the [q, k, q - k + 1] extended narrow-sense Reed-Solomon code over 𝔽_q is GRS and MDS

        # narrrow-sense RS codes are GRS codes with n = q - 1, γ_i = α^i, and v_i = 1 for 0 <= i <= n - 1

        # MacWilliams & Sloane, p. 335
        E = GF(8)
        α = gen(E)
        γ = [α^i for i in 0:6]
        v = [E(1) for _ in 1:7]
        A = AlternateCode(GF(2), 2, v, γ)
        @test length(A) == 7
        @test dimension(A) == 3
        @test minimum_distance(A) == 4

        v = γ
        A = AlternateCode(GF(2), 2, v, γ)
        @test length(A) == 7
        @test dimension(A) == 4
        @test minimum_distance(A) == 3

        # maybe problem on p. 338
        # would require a weird setup of scalars given their def of H based on g

        # Ling & Zing, Example 9.2.4 (i), p. 193
        # E = GF(2^m)
        # F = GF(2)
        # v = (1, α, α^2, ..., α^(2^m - 2)]
        # all non-zero
        # γ = collect(E)[2:end]
        # A = AlternateCode(F, 2^m - 2, v, γ)
        # change_base_ring to Oscar.Native.Nemo.GF(2)
        # should be equal to HammingCode(2, m)
        # m = 4
        # E = GF(2^m)
        # α = gen(E)
        # F = GF(2)
        # v = [α^i for i in 0:2^m - 2]
        # γ = [α^i for i in 1:2^m - 1]
        # A = AlternateCode(F, 2^m - 2, v, γ)
        # TODO not really understanding if I'm using notation right

        # Ling & Zing, Example 9.2.4 (ii), p. 194
        # BCH codes are alternate codes

        # Ling & Zing, Example 9.2.4 (iii), p. 194
        E = GF(2^3)
        n = 6
        α = gen(E)
        # α is root of α^3 + α + 1 = 0
        v = [E(1) for _ in 1:n]
        γ = [α, α^2, α^3, α^4, α^5, α^6]
        A = AlternateCode(GF(2), 3, v, γ)
        @test length(A) == 6
        @test dimension(A) == 2
        @test minimum_distance(A) == 4
        
        # TODO write tests for GRS(Γ), GRS(A), etc
    end
end
