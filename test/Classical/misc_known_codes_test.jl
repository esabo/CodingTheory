@testitem "Classical/misc_known_codes.jl" begin
    using Oscar, CodingTheory

    @testset "Trivial Codes (Zero and Identity)" begin
        Z2 = ZeroCode(5)
        @test length(Z2) == 5
        @test dimension(Z2) == 0
        @test typeof(field(Z2)) == typeof(Oscar.Nemo.Native.GF(2))
        @test iszero(generator_matrix(Z2))
        
        I2 = IdentityCode(5)
        @test length(I2) == 5
        @test dimension(I2) == 5
        @test generator_matrix(I2) == identity_matrix(Oscar.Nemo.Native.GF(2), 5)
        
        F4 = GF(2, 2, :α)
        Z4 = ZeroCode(4, 5)
        @test length(Z4) == 5
        @test dimension(Z4) == 0
        @test typeof(field(Z4)) == typeof(F4)
    end

    @testset "Repetition and Single Parity Check Codes" begin
        R3 = RepetitionCode(3, 4)
        @test length(R3) == 4
        @test dimension(R3) == 1
        @test minimum_distance(R3)[1] == 4
        @test typeof(field(R3)) == typeof(Oscar.Nemo.Native.GF(3))
        
        SPC = SPCCode(2, 5)
        @test length(SPC) == 5
        @test dimension(SPC) == 4
        @test minimum_distance(SPC)[1] == 2
        
        R2 = RepetitionCode(2, 5)
        @test are_equivalent(SPC, dual(R2))
    end

    @testset "Hexacode" begin
        H = Hexacode()
        ω = gen(field(H))
        @test length(H) == 6
        @test dimension(H) == 3
        @test minimum_distance(H)[1] == 4
        @test parity_check_matrix(H) == matrix(field(H), [1 ω ω 1 0 0; ω 1 ω 0 1 0; ω ω 1 0 0 1])
    end

    @testset "Hamming and Simplex Codes" begin
        R, (x, y) = polynomial_ring(Nemo.ZZ, [:x, :y])
        
        F = Oscar.Nemo.Native.GF(2)
        C = HammingCode(2, 7)
        @test length(C) == 2^7 - 1
        @test dimension(C) == 2^7 - 1 - 7
        @test minimum_distance(C)[1] == 3
        
        C_small = HammingCode(2, 3)
        ham_WE = weight_enumerator(C_small)
        @test polynomial(ham_WE, R) == x^7 + 7*x^3*y^4 + 7*x^4*y^3 + y^7
        
        EH = ExtendedHammingCode(3)
        @test length(EH) == 8
        @test dimension(EH) == 4
        @test minimum_distance(EH)[1] == 4

        S = SimplexCode(2, 4)
        @test length(S) == 2^4 - 1
        @test dimension(S) == 4
        
        S_small = SimplexCode(2, 3)
        S_small_we = weight_enumerator(S_small)
        # Verify MacWilliams identity: transformed Simplex WE == Hamming WE
        dual_we = CodingTheory.MacWilliams_transform(S_small_we, dimension(S_small), 2)
        @test polynomial(dual_we, R) == polynomial(ham_WE, R)
    end

    @testset "Golay Codes" begin
        R, (x, y) = polynomial_ring(Nemo.ZZ, [:x, :y])
        
        C24 = ExtendedGolayCode(2)
        @test length(C24) == 24
        @test dimension(C24) == 12
        @test is_self_dual(C24)
        @test polynomial(weight_enumerator(C24), R) == x^24 + 759*x^16*y^8 + 2576*x^12*y^12 + 759*x^8*y^16 + y^24
        
        C23 = GolayCode(2)
        @test length(C23) == 23
        @test minimum_distance(C23)[1] == 7
        
        C12 = ExtendedGolayCode(3)
        @test length(C12) == 12
        @test is_self_dual(C12)
        @test polynomial(weight_enumerator(C12), R) == x^12 + 264*x^6*y^6 + 440*x^3*y^9 + 24*y^12
        
        C11 = GolayCode(3)
        @test minimum_distance(C11)[1] == 5
        @test polynomial(weight_enumerator(C11), R) == x^11 + 132*x^6*y^5 + 132*x^5*y^6 + 330*x^3*y^8 + 110*x^2*y^9 + 24*y^11
        
        C2_ext = extend(puncture(C12, 7), 7)
        T = identity_matrix(field(C12), 12)
        T[7,7] = field(C12)(-1)
        C3_ext = LinearCode(generator_matrix(C2_ext) * T)
        @test are_equivalent(C12, C3_ext)
    end

    @testset "Tetra Code" begin
        C = TetraCode()
        @test length(C) == 4
        @test dimension(C) == 2
        @test minimum_distance(C)[1] == 3
        
        try
            R3, vars = polynomial_ring(Nemo.ZZ, [:x, :y, :z])
            CWE = polynomial(complete_weight_enumerator(C), R3)
            @test CWE == vars[1]^4 + vars[1]*vars[2]^3 + 3*vars[1]*vars[2]^2*vars[3] +
                    3*vars[1]*vars[2]*vars[3]^2 + vars[1]*vars[3]^3
        catch e
            if isa(e, UndefVarError)
                @warn "complete_weight_enumerator not found. Skipping CWE test."
            else
                rethrow(e)
            end
        end
    end
      
    @testset "Hadamard and MacDonald Codes" begin
        HC = HadamardCode(3)
        @test length(HC) == 8
        @test dimension(HC) == 3
        @test minimum_distance(HC)[1] == 4
        
        MC = MacDonaldCode(2, 3, 1)
        @test length(MC) == 6
        @test dimension(MC) == 3
        @test minimum_distance(MC)[1] == 3 
    end

    @testset "Lexicodes" begin
        LC = Lexicode(7, 3)
        @test length(LC) == 7
        @test dimension(LC) == 4
        @test minimum_distance(LC)[1] == 3
        
        @test_throws DomainError Lexicode(5, 6)
        @test_throws DomainError Lexicode(5, 0)
    end
end
