@testitem "Quantum/misc_known_codes.jl" begin
    using JLD2
    using CodingTheory

    @testset "Misc known Quantum codes" begin
        # TODO maybe we want to test to make sure we can compute these distances
        S = Q9143()
        @test S.n == 9
        @test S.k == 1
        @test S.r == 4
        @test S.d_dressed == 3
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasGauges()

        S = Q513()
        @test S.n == 5
        @test S.k == 1
        @test S.d == 3
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        S = Q713()
        @test S.n == 7
        @test S.k == 1
        @test S.d == 3
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        S = Q913()
        @test S.n == 9
        @test S.k == 1
        @test S.d == 3
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        S = Q412()
        @test S.n == 4
        @test S.k == 1
        @test S.d == 2
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        S = Q422()
        @test S.n == 4
        @test S.k == 2
        @test S.d == 2
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        S = Q511()
        @test S.n == 5
        @test S.k == 1
        @test S.d == 1
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        S = Q823()
        @test S.n == 8
        @test S.k == 2
        @test S.d == 3
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        S = Q1513()
        @test S.n == 15
        @test S.k == 1
        @test S.d == 3
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        S = Q1573()
        @test S.n == 15
        @test S.k == 7
        @test S.d == 3
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        S = RotatedSurfaceCode(5)
        @test S.n == 25
        @test S.k == 1
        @test S.d == 5
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        # TODO: fix bug
        # S = RotatedSurfaceCode(4)
        # @test S.n == 16
        # @test S.k == 1
        # @test S.d == 3
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        S = XZZXSurfaceCode(5)
        @test S.n == 25
        @test S.k == 1
        @test_broken S.d == 3
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 3
        S = TriangularColorCode488(d)
        # @test S.n == Int(d^2 // 2 + d - 1 // 2)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 5
        S = TriangularColorCode488(d)
        # @test S.n == Int(d^2 // 2 + d - 1 // 2)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 7
        S = TriangularColorCode488(d)
        # @test S.n == Int(d^2 // 2 + d - 1 // 2)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 9
        S = TriangularColorCode488(d)
        # @test S.n == Int(d^2 // 2 + d - 1 // 2)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 11
        S = TriangularColorCode488(d)
        # @test S.n == Int(d^2 // 2 + d - 1 // 2)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 13
        S = TriangularColorCode488(d)
        # @test S.n == Int(d^2 // 2 + d - 1 // 2)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 15
        S = TriangularColorCode488(d)
        # @test S.n == Int(d^2 // 2 + d - 1 // 2)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 17
        S = TriangularColorCode488(d)
        # @test S.n == Int(d^2 // 2 + d - 1 // 2)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 19
        S = TriangularColorCode488(d)
        # @test S.n == Int(d^2 // 2 + d - 1 // 2)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        # d = 21
        # S = TriangularColorCode488(d)
        # # @test S.n == Int(d^2 // 2 + d - 1 // 2)
        # @test S.k == 1
        # @test S.d == d
        # @test LogicalTrait(typeof(S)) == HasLogicals()
        # @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 3
        S = TriangularColorCode666(d)
        @test S.n == Int(3 * d^2 // 4 + 1 // 4)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 5
        S = TriangularColorCode666(d)
        @test S.n == Int(3 * d^2 // 4 + 1 // 4)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 7
        S = TriangularColorCode666(d)
        @test S.n == Int(3 * d^2 // 4 + 1 // 4)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 9
        S = TriangularColorCode666(d)
        @test S.n == Int(3 * d^2 // 4 + 1 // 4)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 11
        S = TriangularColorCode666(d)
        @test S.n == Int(3 * d^2 // 4 + 1 // 4)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 13
        S = TriangularColorCode666(d)
        @test S.n == Int(3 * d^2 // 4 + 1 // 4)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 15
        S = TriangularColorCode666(d)
        @test S.n == Int(3 * d^2 // 4 + 1 // 4)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 17
        S = TriangularColorCode666(d)
        @test S.n == Int(3 * d^2 // 4 + 1 // 4)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 19
        S = TriangularColorCode666(d)
        @test S.n == Int(3 * d^2 // 4 + 1 // 4)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        d = 21
        S = TriangularColorCode666(d)
        @test S.n == Int(3 * d^2 // 4 + 1 // 4)
        @test S.k == 1
        @test S.d == d
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        S = ToricCode(2)
        @test S.n == 2 * 2^2
        @test S.k == 2
        @test S.d == 2
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        S = ToricCode(3)
        @test S.n == 2 * 3^2
        @test S.k == 2
        @test S.d == 3
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        # BUG
        # S = PlanarSurfaceCode(3, 4)
        # @test S.n == (3 * 4 + 2 * 3)
        # @test S.k == 1
        # @test S.d == 3
        # @test LogicalTrait(typeof(S)) == HasLogicals()
        # @test GaugeTrait(typeof(S)) == HasNoGauges()

        S = HCode(8)
        @test S.n == 12
        @test S.k == 8
        @test S.d == 2
        @test LogicalTrait(typeof(S)) == HasLogicals()
        @test GaugeTrait(typeof(S)) == HasNoGauges()

        # TODO: need to remove quadratic
        # S = XYSurfaceCode(3, 4)
        # @test S.n == (3 * 4 + 2 * 3)
        # @test S.k == 1
        # @test_broken S.d == 3
        # @test LogicalTrait(typeof(S)) == HasLogicals()
        # @test GaugeTrait(typeof(S)) == HasNoGauges()
    end
end
