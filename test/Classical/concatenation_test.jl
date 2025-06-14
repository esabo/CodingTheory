@testitem "Classical/concatenation.jl" begin
    using Oscar, CodingTheory

    @testset "concatenate codes" begin
        # same field examples come from https://benjamin-hackl.at/downloads/BA-Thesis_Hackl.pdf
        # "Concatenated Error Correcting Codes: Galois and Binary Concatenation"
        F = Oscar.Nemo.Native.GF(2)
        Ham = matrix(
            F,
            4,
            8,
            [
                1 0 0 0 0 1 1 1;
                0 1 0 0 1 0 1 1;
                0 0 1 0 1 1 0 1;
                0 0 0 1 1 1 1 0
            ],
        )
        G_in = matrix(
            F,
            4,
            8,
            [
                1 0 0 0 0 1 0 0;
                0 1 0 0 0 0 1 0;
                0 0 1 0 1 0 0 1;
                0 0 0 1 1 1 1 0
            ],
        )
        # G_out should be the [8, 4, 4] extended binary Hamming code
        C_Ham = LinearCode(Ham)
        @test dimension(C_Ham) == 4
        @test minimum_distance(C_Ham) == 4
        # G_in should be [8, 4, 2]
        C_in = LinearCode(G_in)
        @test dimension(C_in) == 4
        @test minimum_distance(C_in) == 2
        C = concatenate(C_Ham, C_in)
        G_final = matrix(
            F,
            4,
            16,
            [
                1 0 0 0 0 1 0 0 0 1 1 1 0 1 0 1;
                0 1 0 0 0 0 1 0 1 0 1 1 0 0 1 1;
                0 0 1 0 1 0 0 1 1 1 0 1 1 0 0 0;
                0 0 0 1 1 1 1 0 1 1 1 0 1 1 1 1
            ],
        )
        @test generator_matrix(C) == G_final
        @test_broken minimum_distance(C) == 7

        C = concatenate(C_Ham, C_Ham)
        G_final = matrix(
            F,
            4,
            16,
            [
                1 0 0 0 0 1 1 1 0 1 1 1 1 0 0 0;
                0 1 0 0 1 0 1 1 1 0 1 1 0 1 0 0;
                0 0 1 0 1 1 0 1 1 1 0 1 0 0 1 0;
                0 0 0 1 1 1 1 0 1 1 1 0 0 0 0 1
            ],
        )
        @test generator_matrix(C) == G_final
        @test_broken minimum_distance(C) == 8

        G_out = matrix(
            F,
            4,
            12,
            [
                1 0 0 0 0 0 1 1 0 1 1 0
                0 1 0 0 1 0 0 0 1 0 0 1
                0 0 1 0 1 1 0 0 1 0 0 0
                0 0 0 1 0 0 1 1 0 0 1 1
            ],
        )
        C_out = LinearCode(G_out)
        @test dimension(C_out) == 4
        @test minimum_distance(C_out) == 4
        C = concatenate(C_out, C_Ham)
        G_final = matrix(
            F,
            4,
            24,
            [
                1 0 0 0 0 1 1 1 0 0 1 1 0 0 1 1 0 1 1 0 0 1 1 0;
                0 1 0 0 1 0 1 1 1 0 0 0 0 1 1 1 1 0 0 1 1 0 0 1;
                0 0 1 0 1 1 0 1 1 1 0 0 1 1 0 0 1 0 0 0 0 1 1 1;
                0 0 0 1 1 1 1 0 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1
            ],
        )
        @test generator_matrix(C) == G_final
        @test_broken minimum_distance(C) == 8

        G_out = matrix(
            F,
            4,
            12,
            [
                1 0 0 0 1 0 1 1 0 0 1 0;
                0 1 0 0 0 1 1 1 1 0 0 0;
                0 0 1 0 1 1 1 1 0 1 1 1;
                0 0 0 1 1 1 1 0 1 1 1 0
            ],
        )
        C_out = LinearCode(G_out)
        @test dimension(C_out) == 4
        @test minimum_distance(C_out) == 5
        C = concatenate(C_out, C_Ham)
        G_final = matrix(
            F,
            4,
            24,
            [
                1 0 0 0 0 1 1 1 1 0 1 1 0 1 0 0 0 0 1 0 1 1 0 1;
                0 1 0 0 1 0 1 1 0 1 1 1 1 0 0 0 1 0 0 0 0 1 1 1;
                0 0 1 0 1 1 0 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 0 0;
                0 0 0 1 1 1 1 0 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1
            ],
        )
        @test generator_matrix(C) == G_final
        @test_broken minimum_distance(C) == 12

        # claims this is ExtendedGolayCode(2) but I can't find the equivalence
        G_in = matrix(
            F,
            12,
            24,
            [
                1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 1 1 0 0 0 1 1;
                0 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 1 0 0 1 0;
                0 0 1 0 0 0 0 0 0 0 0 0 1 1 0 1 0 0 1 0 1 0 1 1;
                0 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1 1 0 1 1 0;
                0 0 0 0 1 0 0 0 0 0 0 0 1 1 0 0 1 1 0 1 1 0 0 1;
                0 0 0 0 0 1 0 0 0 0 0 0 0 1 1 0 0 1 1 0 1 1 0 1;
                0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 1 0 0 1 1 0 1 1 1;
                0 0 0 0 0 0 0 1 0 0 0 0 1 0 1 1 0 1 1 1 1 0 0 0;
                0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 1 1 0 1 1 1 1 0 0;
                0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 1 1 0 1 1 1 1 0;
                0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 1 1 0 0 0 1 1 0 1;
                0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 1 1 0 0 0 1 1 1
            ],
        )
        C_in = LinearCode(G_in)
        C = concatenate(C_out, C_in)
        G_final = matrix(
            F,
            4,
            24,
            [
                1 0 0 0 1 0 1 1 0 0 1 0 0 1 0 1 1 1 1 1 1 0 0 0
                0 1 0 0 0 1 1 1 1 0 0 0 0 1 0 0 0 0 0 0 1 1 0 0
                0 0 1 0 1 1 1 1 0 1 1 1 0 0 1 1 0 1 0 0 0 1 0 0
                0 0 0 1 1 1 1 0 1 1 1 0 1 0 0 1 0 0 0 1 1 0 1 0
            ],
        )
        @test generator_matrix(C) == G_final
        @test_broken minimum_distance(C) == 8

        # Bossert example 9.2
        A1 = extend(HammingCode(2, 3))
        A2 = RepetitionCode(4, 8)
        # A2_e = expanded_code(A2, Oscar.Nemo.Native.GF(2), primitive.basis(A2.F, Oscar.Nemo.Native.GF(2)))
        # outers = [A1, A2_e]
        inners = [RepetitionCode(2, 4), SPCCode(2, 4)]
        # C = concatenate(outers, inners)
        # @test C.n == 32
        # @test C.k == 6
        # should also work without expanding, and in that case it should know the distance:
        # outers = [A1, A2]
        # C = concatenate(outers, inners)
        # @test C.n == 32
        # @test C.k == 6
        # @test C.d == 16

        # Huffman and Pless, example 5.5.2
        # Note that when the inner code is an identity code, the expanded code should already be the final answer. Also test if the generalized concatenate gives the same answer as the usual concatenate
        # G = matrix(Oscar.Nemo.Native.GF(2), [1  0  0  0  0  0  1  0  0  1  0  1
        #                    0  1  0  0  0  0  0  1  1  1  1  1
        #                    0  0  1  0  0  0  0  1  1  0  0  1
        #                    0  0  0  1  0  0  1  1  0  1  1  1
        #                    0  0  0  0  1  0  0  1  0  1  1  0
        #                    0  0  0  0  0  1  1  1  1  1  0  1])
        # C = Hexacode()
        # Ce = expanded_code(C, Oscar.Nemo.Native.GF(2), basis(C.F, Oscar.Nemo.Native.GF(2)))
        # @test Ce.G == G
        # C1 = concatenate([Ce], [IdentityCode(2, 2)]) # Generalized concat, pre-expanded
        # @test C1.G == G
        # C2 = concatenate([C], [IdentityCode(2, 2)]) # Generalized concat, not pre-expanded
        # @test C2.G == G
        # C3 = concatenate(Ce, IdentityCode(2, 2)) # usual concat, pre-expanded
        # @ test C3.G == G
        # TODO: This one doesn't work (usual concat, not pre-expanded):
        # C4 = concatenate(C, IdentityCode(2, 2)) # usual concat, not pre-expanded
        # @test C4.G == G
    end
end
