@testitem "Code containment" begin
    using Oscar, CodingTheory

    # ⊆, ⊂, and is_subcode are all non-strict; ⊊ is the proper version
    @testset "Linear Codes" begin
        C = HammingCode(2, 3)
        D = RepetitionCode(2, 7)

        @test D ⊆ C
        @test D ⊂ C
        @test is_subcode(D, C)
        @test D ⊊ C
        @test !(C ⊆ D)
        @test !(C ⊊ D)

        @test C ⊆ C
        @test C ⊂ C
        @test is_subcode(C, C)
        @test !(C ⊊ C)

        # containment requires equal lengths
        @test !(RepetitionCode(2, 3) ⊆ C)
    end

    @testset "Cyclic Codes" begin
        # a larger designed distance gives a smaller code
        C = BCHCode(2, 15, 3)
        D = BCHCode(2, 15, 5)

        @test D ⊆ C
        @test D ⊂ C
        @test is_subcode(D, C)
        @test D ⊊ C
        @test !(C ⊆ D)

        # the cyclic methods agree with the general linear ones on equality
        @test C ⊆ C
        @test C ⊂ C
        @test is_subcode(C, C)
        @test !(C ⊊ C)
    end

    @testset "Dual Containment" begin
        # a self-orthogonal code is contained in its dual
        C = RepetitionCode(2, 4)
        @test C ⊆ dual(C)
        @test C ⊊ dual(C)

        # a self-dual code equals its dual
        E = ExtendedHammingCode(3)
        @test E ⊆ dual(E)
        @test dual(E) ⊆ E
        @test !(E ⊊ dual(E))
    end
end
