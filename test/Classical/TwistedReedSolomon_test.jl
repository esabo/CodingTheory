@testitem "Classical/TwistedReedSolomon.jl" begin
    using Oscar, CodingTheory

    @testset "Twisted Reed-Solomon Codes" begin
        # https://arxiv.org/abs/2107.06945
        # Example 4
        # η = 0 gives original RS code
        # TODO determine which RS codes with the below

        # Example 4
        F = GF(3, 2, :ω);
        ω = gen(F);
        k = 5;
        α = collect(F);
        t = [2];
        h = [2];
        sqs = [i^2 for i in F];
        η = [setdiff(α, sqs)[1]];
        C = TwistedReedSolomonCode(k, α, t, h, η);
        G = zero_matrix(F, k, length(α));
        for c in 1:length(α)
            G[1, c] = α[c]^0
            G[2, c] = α[c]^1
            G[3, c] = α[c]^2 + η[1] * α[c]^6
            G[4, c] = α[c]^3
            G[5, c] = α[c]^4
        end
        @test G == generator_matrix(C)

        # Example 4
        F = GF(2, 3, :ω);
        ω = gen(F);
        k = 5;
        α = collect(F);
        t = [1, 3, 3];
        h = [4, 4, 2];
        # paper does not specify, random suffices because here we are just matching the form of G
        η = [F(0), ω, ω^2];
        C = TwistedReedSolomonCode(k, α, t, h, η);
        G = zero_matrix(F, k, length(α));
        for c in 1:length(α)
            G[1, c] = α[c]^0
            G[2, c] = α[c]^1
            G[3, c] = α[c]^2 + η[3] * α[c]^7
            G[4, c] = α[c]^3
            G[5, c] = α[c]^4 + η[1] * α[c]^5 + η[2] * α[c]^7
        end
        @test G == generator_matrix(C)

        # BUG can't quite get these parameters to match mine
        # # https://arxiv.org/pdf/2211.06066
        # # this paper has shifted indices wrt the original definition
        # # Example 3.6
        # F = Oscar.Nemo.Native.GF(11)
        # l = 2
        # α = [F(1), F(2), F(3), F(5), F(6), F(8), F(9), F(10)]
        # k = 3
        # h = [k - l + i - 1 for i in 1:l]
        # t = [i for i in 1:l]
        # η = [F(0), F(0)]
        # C = TwistedReedSolomonCode(k, α, t, h, η);
        # @test length(C) == 8
        # @test dimension(C) == 3
        # # @test minimum_distance(C) == 6
        # # @test is_MDS(C)

        # η = [F(2), F(9)]
        # C = TwistedReedSolomonCode(k, α, t, h, η);
        # @test length(C) == 8
        # @test dimension(C) == 3
        # # @test minimum_distance(C) == 6
        # # @test is_MDS(C)

        # k = 4
        # h = [k - l + i - 1 for i in 1:l]
        # t = [i for i in 1:l]
        # η = [F(0), F(0)]
        # C = TwistedReedSolomonCode(k, α, t, h, η);
        # @test length(C) == 8
        # @test dimension(C) == 4
        # # @test minimum_distance(C) == 5
        # # @test is_MDS(C)

        # η = [F(4), F(4)]
        # C = TwistedReedSolomonCode(k, α, t, h, η);
        # @test length(C) == 8
        # @test dimension(C) == 4
        # # @test minimum_distance(C) == 5
        # # @test is_MDS(C)

        # η = [F(6), F(6)]
        # C = TwistedReedSolomonCode(k, α, t, h, η);
        # @test length(C) == 8
        # @test dimension(C) == 4
        # # @test minimum_distance(C) == 5
        # # @test is_MDS(C)

        # k = 5
        # h = [k - l + i - 1 for i in 1:l]
        # t = [i for i in 1:l]
        # η = [F(0), F(0)]
        # C = TwistedReedSolomonCode(k, α, t, h, η);
        # @test length(C) == 8
        # @test dimension(C) == 5
        # # @test minimum_distance(C) == 4
        # # @test is_MDS(C)

        # η = [F(9), F(10)]
        # C = TwistedReedSolomonCode(k, α, t, h, η);
        # @test length(C) == 8
        # @test dimension(C) == 5
        # # @test minimum_distance(C) == 4
        # # @test is_MDS(C)

        # # Example 3.7
        # F = Oscar.Nemo.Native.GF(13)
        # l = 3
        # α = [F(0), F(1), F(2), F(3), F(4), F(5), F(6), F(9), F(10), F(12)]
        # k = 5
        # h = [k - l + i - 1 for i in 1:l]
        # t = [i for i in 1:l]
        # η = [F(2), F(3), F(6)]
        # C = TwistedReedSolomonCode(k, α, t, h, η);
        # @test length(C) == 10
        # @test dimension(C) == 5
        # # @test minimum_distance(C) == 6
        # # @test is_MDS(C)

        # # this paper also does twisted-GRS codes
        # # the above examples are with v = 1 there
    end
end
