@testitem "Classical/TwistedReedSolomon.jl" begin
    using Oscar, CodingTheory

    @testset "Twisted Reed-Solomon Codes" begin
        # Example 4: https://arxiv.org/abs/2107.06945
        F = GF(3, 2, :ω)
        ω = gen(F)
        k = 5
        α = collect(F)
        t = [2]
        h = [2]
        sqs = [i^2 for i in F]
        η = [setdiff(α, sqs)[1]]
        C = TwistedReedSolomonCode(k, α, t, h, η)
        G = zero_matrix(F, k, length(α))
        for c in 1:length(α)
            G[1, c] = α[c]^0
            G[2, c] = α[c]^1
            G[3, c] = α[c]^2 + η[1] * α[c]^6
            G[4, c] = α[c]^3
            G[5, c] = α[c]^4
        end
        @test G == generator_matrix(C)

        # Example 4 with more twists
        F_2 = GF(2, 3, :ω)
        ω_2 = gen(F_2)
        k_2 = 5
        α_2 = collect(F_2)
        t_2 = [1, 3, 3]
        h_2 = [4, 4, 2]
        η_2 = [F_2(0), ω_2, ω_2^2]
        C_2 = TwistedReedSolomonCode(k_2, α_2, t_2, h_2, η_2)
        G_2 = zero_matrix(F_2, k_2, length(α_2))
        for c in 1:length(α_2)
            G_2[1, c] = α_2[c]^0
            G_2[2, c] = α_2[c]^1
            G_2[3, c] = α_2[c]^2 + η_2[3] * α_2[c]^7
            G_2[4, c] = α_2[c]^3
            G_2[5, c] = α_2[c]^4 + η_2[1] * α_2[c]^5 + η_2[2] * α_2[c]^7
        end
        @test G_2 == generator_matrix(C_2)
    end

    @testset "Getters, Duals, and Parity Matrices" begin
        F = GF(2, 3, :a)
        a = gen(F)
        k = 3
        α = collect(F)
        t = [2]
        h = [1]
        η = [a]
        
        C = TwistedReedSolomonCode(k, α, t, h, η)
        
        # Test Parity Check Matrix Orthogonality
        G = generator_matrix(C)
        H = parity_check_matrix(C)
        @test iszero(G * transpose(H))
        
        # Getters
        @test twist_vector(C) == t
        @test hook_vector(C) == h
        @test coefficient_vector(C) == η
        @test number_of_twists(C) == 1
        
        # Dual Code tracking
        C_dual = dual(C)
        @test dimension(C_dual) == length(C) - dimension(C)
        @test twist_vector(C_dual) == k .- h
        @test hook_vector(C_dual) == (length(C) - k) .- t
        @test coefficient_vector(C_dual) == -η
    end

    @testset "Random Generation and Error Handling" begin
        F = GF(11)
        
        # Generate a valid random TRS code
        C_rand = RandomTwistedReedSolomonCode(F, 8, 4, 2)
        @test length(C_rand) == 8
        @test dimension(C_rand) == 4
        @test number_of_twists(C_rand) == 2
        
        # Domain errors
        @test_throws DomainError RandomTwistedReedSolomonCode(F, 12, 4, 1) # n > |F|
        @test_throws DomainError RandomTwistedReedSolomonCode(F, 8, 9, 1)  # k > n
        @test_throws DomainError RandomTwistedReedSolomonCode(F, 8, 4, 100) # Too many twists
    end
end