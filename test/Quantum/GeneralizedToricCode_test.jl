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
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 12
        @test dimension(S) == 4
        
        f = 1 + x + y;
        g = 1 + y + x;
        a1 = (0, 7);
        a2 = (1, 2);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 14
        @test dimension(S) == 6

        f = 1 + x + x * y;
        g = 1 + y + x * y;
        a1 = (0, 3);
        a2 = (3, 0);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 18
        @test dimension(S) == 4
                   
        f = 1 + x + x * y;
        g = 1 + y + x * y;
        a1 = (0, 3);
        a2 = (4, 2);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 24
        @test dimension(S) == 4
                   
        f = 1 + x + x^-1 * y;
        g = 1 + y + x * y;
        a1 = (0, 7);
        a2 = (2, 3);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 28
        @test dimension(S) == 6
       
        f = 1 + x + x^2;
        g = 1 + y + x^2;
        a1 = (0, 3);
        a2 = (5, 1);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 30
        @test dimension(S) == 4     
        
        f = 1 + x + x^-1; 
        g = 1 + y + y^-1;
        a1 = (0, 9);
        a2 = (2, 4);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 36
        @test dimension(S) == 4 
                
        f = 1 + x + x * y;
        g = 1 + y + x * y^-1;
        a1 = (0, 7);
        a2 = (3, 2);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 42
        @test dimension(S) == 6
        
        f = 1 + x + x^2;
        g = 1 + y + x^2;
        a1 = (0, 3);
        a2 = (8, 1);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 48
        @test dimension(S) == 4
                   
        f = 1 + x + x^-1;
        g = 1 + y + x^3 * y^2;
        a1 = (0, 3);
        a2 = (9, 0);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 54
        @test dimension(S) == 8
                   
        f = 1 + x + y^-2;
        g = 1 + y + x^-2;
        a1 = (0, 7);
        a2 = (4, 3);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 56
        @test dimension(S) == 6
                   
        f = 1 + x + y^-2;
        g = 1 + y + x^2;
        a1 = (0, 10);
        a2 = (3, 3);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 60
        @test dimension(S) == 8
                   
        f = 1 + x + x^-1 * y;
        g = 1 + y + x^-1 * y^-1;
        a1 = (0, 31);
        a2 = (1, 13);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 62
        @test dimension(S) == 10
                   
        f = 1 + x + x^-2 * y^-1;
        g = 1 + y + x^2 * y;
        a1 = (0, 3);
        a2 = (11, 2);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 66
        @test dimension(S) == 4
                   
        f = 1 + x + x * y;
        g = 1 + y + x * y^-1;
        a1 = (0, 7);
        a2 = (5, 1);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 70
        @test dimension(S) == 6
                   
        f = 1 + x + x^-1 * y^3;
        g = 1 + y + x^3 * y^-1;
        a1 = (0, 12);
        a2 = (3, 3);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 72
        @test dimension(S) == 8
        
        f = 1 + x + x^-2 * y^-1;
        g = 1 + y + x^2 * y;
        a1 = (0, 3);
        a2 = (13, 1);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 78
        @test dimension(S) == 4
                   
        f = 1 + x + x^-2;
        g = 1 + y + x^-2 * y^2;
        a1 = (0, 14);
        a2 = (3, -6);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 84
        @test dimension(S) == 6
                   
        f = 1 + x + x^-1 * y^-3;
        g = 1 + y + x^3 * y^-1;
        a1 = (0, 15);
        a2 = (3, -6);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 90
        @test dimension(S) == 8
                   
        f = 1 + x + x^-2 * y;
        g = 1 + y + x * y^-2;
        a1 = (0, 12);
        a2 = (4, 2);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 96
        @test dimension(S) == 4

        f = 1 + x + x^-1 * y^2;
        g = 1 + y + x^-2 * y^-1;
        a1 = (0, 7);
        a2 = (7, 0);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 98
        @test dimension(S) == 6
                   
        f = 1 + x + x^-3 * y;
        g = 1 + y + x^3 * y^2;
        a1 = (0, 3);
        a2 = (17, 2);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 102
        @test dimension(S) == 4
                   
        f = 1 + x + x^-1 * y^-3;
        g = 1 + y + x^3 * y^-1;
        a1 = (0, 9);
        a2 = (6, 0);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 108
        @test dimension(S) == 8
                   
        f = 1 + x + x^-1 * y^3;
        g = 1 + y + x^3 * y^-1;
        a1 = (0, 9);
        a2 = (6, 0);
        S = CSSCode(FiniteGeneralizedToricCode(f, g, a1, a2))
        @test length(S) == 108
        @test dimension(S) == 8
           
    end
end
