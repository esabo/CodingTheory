@testitem "Classical/quasi-cyclic_code.jl" begin
    using Oscar, CodingTheory

    @testset "QuasiCyclicCode" begin
        F = GF(2)
        v = matrix(F, 1, 8, [1, 0, 1, 1, 1, 0, 0, 0])
        v2 = matrix(F, 1, 8, [1, 1, 1, 0, 0, 0, 1, 0])
        C = QuasiCyclicCode([v, v2], 2, false)
        # v and v2 are shifts of each other
        @test C.k == 4
        v2 = matrix(F, 1, 8, [1, 0, 1, 0, 1, 0, 1, 0])
        C = QuasiCyclicCode([v, v2], 2, false)
        v = [matrix(F, 1, 4, [1, 0, 1, 1]), matrix(F, 1, 4, [0, 0, 0, 1]), matrix(F, 1, 4, [1, 1, 1, 1]), matrix(F, 1, 4, [0, 0, 0, 0])]
        C2 = QuasiCyclicCode(v, 2, true)
        @test are_equivalent(C, C2)
    end
end
