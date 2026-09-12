@testitem "tilings.jl" begin
    using CodingTheory

    C = CoxeterMatrix(3, [1, 3, 2, 1, 7, 1])
    @test size(C) == (3, 3)
    @test C[1, 2] == C[2, 1] == 3
    @test C[2, 3] == 7
    @test_throws BoundsError C[0, 1]
    @test_throws ArgumentError CoxeterMatrix(3, [1, 2])

    G = r_s_group(3, 7)
    @test G.cox_mat == C
    @test length(G.generators) == 3

    @test_throws ArgumentError tetrahedron_group([2, 3])
    @test_throws ArgumentError tetrahedron_group([1, 2, 2, 2, 2, 2])
    T = tetrahedron_group(fill(2, 6))
    @test size(T.cox_mat) == (4, 4)
end
