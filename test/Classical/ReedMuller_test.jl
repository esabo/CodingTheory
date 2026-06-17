@testitem "Classical/ReedMuller.jl" begin
    using Oscar, CodingTheory, Random

    @testset "Misc Known ReedMullerCodes" begin
        F = Oscar.Nemo.Native.GF(2)
        # Huffman, Pless, p. 34
        # identity used for RM(1, 1)
        @test CodingTheory._Reed_Muller_generator_matrix(1, 1, true) == matrix(F,
                [1 0;
                 0 1])
        @test generator_matrix(ReedMullerCode(1, 2, true)) == matrix(F,
                [1 0 1 0;
                 0 1 0 1;
                 0 0 1 1])
        @test generator_matrix(ReedMullerCode(1, 3, true)) == matrix(F,
                [1 0 1 0 1 0 1 0;
                 0 1 0 1 0 1 0 1;
                 0 0 1 1 0 0 1 1;
                 0 0 0 0 1 1 1 1])
        @test generator_matrix(ReedMullerCode(2, 3, true)) == matrix(F,
                [1 0 0 0 1 0 0 0;
                 0 1 0 0 0 1 0 0;
                 0 0 1 0 0 0 1 0;
                 0 0 0 1 0 0 0 1;
                 0 0 0 0 1 0 1 0;
                 0 0 0 0 0 1 0 1;
                 0 0 0 0 0 0 1 1])

        # Ling & Xing, p. 119
        # other sources, using [1 1; 0 1] for RM(1, 1)
        @test CodingTheory._Reed_Muller_generator_matrix(1, 1) == matrix(F,
                [1 1;
                 0 1])
        @test generator_matrix(ReedMullerCode(1, 2)) == matrix(F,
                [1 1 1 1;
                 0 1 0 1;
                 0 0 1 1])
        @test generator_matrix(ReedMullerCode(1, 3)) == matrix(F,
                [1 1 1 1 1 1 1 1;
                 0 1 0 1 0 1 0 1;
                 0 0 1 1 0 0 1 1;
                 0 0 0 0 1 1 1 1])
    end

    @testset "Parity Check and Orthogonality" begin
        # RM(r, m) has RM(m - r - 1, m) as its dual
        C = ReedMullerCode(1, 4)
        G = generator_matrix(C)
        H = parity_check_matrix(C)
        
        @test nrows(H) == length(C) - dimension(C)
        # Verify G * H^T = 0
        @test iszero(G * transpose(H))
        
        # Test standard form routing
        H_stand = parity_check_matrix(C, true)
        @test nrows(H_stand) == nrows(H)
    end

    @testset "self-dual property" begin
        # if m is odd and r = (m - 1)/2 then RM(r, m) = RM((m - 1)/2, m) is self-dual
        C = ReedMullerCode(2, 5)
        @test length(C) == 2^5
        @test is_self_dual(C)
        
        # RM(0, m) is the length 2^m repetition code
        @test are_equivalent(ReedMullerCode(0, 3), RepetitionCode(2, 8))
    end

    @testset "SimplexCode from ReedMullerCode" begin
        C = ReedMullerCode(1, 4)
        pC = puncture(C, [1])
        epC = even_subcode(pC)
        S = SimplexCode(2, 4)

        @test are_equivalent(epC, S)

        wt_dist = weight_distribution(C)
        expected = Dict(0 => 1, 8 => 30, 16 => 1)
        @test length(wt_dist) == length(expected)
        @test all(wt_dist[k] == v for (k, v) in expected)
    end

    @testset "Nested and Even-weight Property" begin
        m = rand(3:6)
        r = rand(1:m - 2)
        C = ReedMullerCode(r, m)
        C2 = ReedMullerCode(r + 1, m)
        @test C ⊆ C2

        # RM(m - 1, m) contains all vectors of even weight.
        C_even = ReedMullerCode(m - 1, m)
        wt_dist_even = weight_distribution(C_even)
        @test all(iseven(w) for w in keys(wt_dist_even))

        C_dist = ReedMullerCode(2, 5)
        wt_dist = weight_distribution(C_dist)
        expected = Dict(0 => 1, 8 => 620, 12 => 13888, 16 => 36518, 20 => 13888, 24 => 620, 32 => 1)
        @test length(wt_dist) == length(expected)
        @test all(wt_dist[k] == v for (k, v) in expected)
    end

    @testset "Getters and Random Cryptographic Generation" begin
        r, m = 2, 4
        C = ReedMullerCode(r, m)
        
        # Getters
        @test order(C) == r
        @test RM_r(C) == r
        @test RM_m(C) == m
        
        # Random Permuted Code
        C_perm = RandomPermutedReedMullerCode(r, m)
        @test length(C_perm) == length(C)
        @test dimension(C_perm) == dimension(C)
        
        # Random Boolean Functions
        bf = RandomBooleanFunction(r, m)
        @test size(bf) == (1, 2^m)
        @test base_ring(bf) == field(C)
        
        # Random Coset
        coset_rep = RandomRMCoset(r, m)
        @test size(coset_rep) == (1, 2^m)
        @test base_ring(coset_rep) == field(C)
    end
end