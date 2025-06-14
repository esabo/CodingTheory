@testitem "Quantum/GeneralizedToricCode.jl" begin
    using Oscar
    using CodingTheory

    @testset "Generalized Toric Codes" begin
        # Table 1 of https://arxiv.org/abs/2503.03827
        F = Oscar.Nemo.Native.GF(2);
        R, (x, y) = laurent_polynomial_ring(F, [:x, :y]);
        f = 1 + x + x * y;
        g = 1 + y + x * y;
        a1 = (0, 3);
        a2 = (2, 1);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 12
        @test dimension(C) == 4

        f = 1 + x + y;
        g = 1 + y + x;
        a1 = (0, 7);
        a2 = (1, 2);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 14
        @test dimension(C) == 6

        f = 1 + x + x * y;
        g = 1 + y + x * y;
        a1 = (0, 3);
        a2 = (3, 0);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 18
        @test dimension(C) == 4

        f = 1 + x + x * y;
        g = 1 + y + x * y;
        a1 = (0, 3);
        a2 = (4, 2);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 24
        @test dimension(C) == 4

        f = 1 + x + x^-1 * y;
        g = 1 + y + x * y;
        a1 = (0, 7);
        a2 = (2, 3);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 28
        @test dimension(C) == 6

        f = 1 + x + x^2;
        g = 1 + y + x^2;
        a1 = (0, 3);
        a2 = (5, 1);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 30
        @test dimension(C) == 4

        f = 1 + x + x^-1;
        g = 1 + y + y^-1;
        a1 = (0, 9);
        a2 = (2, 4);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 36
        @test dimension(C) == 4

        f = 1 + x + x * y;
        g = 1 + y + x * y^-1;
        a1 = (0, 7);
        a2 = (3, 2);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 42
        @test dimension(C) == 6

        f = 1 + x + x^2;
        g = 1 + y + x^2;
        a1 = (0, 3);
        a2 = (8, 1);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 48
        @test dimension(C) == 4

        f = 1 + x + x^-1;
        g = 1 + y + x^3 * y^2;
        a1 = (0, 3);
        a2 = (9, 0);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 54
        @test dimension(C) == 8

        f = 1 + x + y^-2;
        g = 1 + y + x^-2;
        a1 = (0, 7);
        a2 = (4, 3);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 56
        @test dimension(C) == 6

        f = 1 + x + y^-2;
        g = 1 + y + x^2;
        a1 = (0, 10);
        a2 = (3, 3);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 60
        @test dimension(C) == 8

        f = 1 + x + x^-1 * y;
        g = 1 + y + x^-1 * y^-1;
        a1 = (0, 31);
        a2 = (1, 13);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 62
        @test dimension(C) == 10

        f = 1 + x + x^-2 * y^-1;
        g = 1 + y + x^2 * y;
        a1 = (0, 3);
        a2 = (11, 2);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 66
        @test dimension(C) == 4

        f = 1 + x + x * y;
        g = 1 + y + x * y^-1;
        a1 = (0, 7);
        a2 = (5, 1);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 70
        @test dimension(C) == 6

        f = 1 + x + x^-1 * y^3;
        g = 1 + y + x^3 * y^-1;
        a1 = (0, 12);
        a2 = (3, 3);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 72
        @test dimension(C) == 8

        f = 1 + x + x^-2 * y^-1;
        g = 1 + y + x^2 * y;
        a1 = (0, 3);
        a2 = (13, 1);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 78
        @test dimension(C) == 4

        f = 1 + x + x^-2;
        g = 1 + y + x^-2 * y^2;
        a1 = (0, 14);
        a2 = (3, -6);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 84
        @test dimension(C) == 6

        f = 1 + x + x^-1 * y^-3;
        g = 1 + y + x^3 * y^-1;
        a1 = (0, 15);
        a2 = (3, -6);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 90
        @test dimension(C) == 8

        f = 1 + x + x^-2 * y;
        g = 1 + y + x * y^-2;
        a1 = (0, 12);
        a2 = (4, 2);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 96
        @test dimension(C) == 4

        f = 1 + x + x^-1 * y^2;
        g = 1 + y + x^-2 * y^-1;
        a1 = (0, 7);
        a2 = (7, 0);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 98
        @test dimension(C) == 6

        f = 1 + x + x^-3 * y;
        g = 1 + y + x^3 * y^2;
        a1 = (0, 3);
        a2 = (17, 2);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 102
        @test dimension(C) == 4

        f = 1 + x + x^-1 * y^-3;
        g = 1 + y + x^3 * y^-1;
        a1 = (0, 9);
        a2 = (6, 0);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 108
        @test dimension(C) == 8

        f = 1 + x + x^-1 * y^3;
        g = 1 + y + x^3 * y^-1;
        a1 = (0, 9);
        a2 = (6, 0);
        C = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(C) == 108
        @test dimension(C) == 8

    end
end
