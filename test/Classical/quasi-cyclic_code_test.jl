@testitem "Classical/quasi-cyclic_code.jl" begin
    using Oscar, CodingTheory

    @testset "Constructors and Equivalent Codes" begin
        # Strictly use Oscar.Nemo.Native.GF(2) for binary
        F = Oscar.Nemo.Native.GF(2)
        
        # Test vector shifts
        v = matrix(F, 1, 8, [1, 0, 1, 1, 1, 0, 0, 0])
        v2 = matrix(F, 1, 8, [1, 1, 1, 0, 0, 0, 1, 0])
        C1 = QuasiCyclicCode([v, v2], 2, false)
        @test dimension(C1) == 4
        @test length(C1) == 8 # Fixed: length is exactly ncols(v)

        v3 = matrix(F, 1, 8, [1, 0, 1, 0, 1, 0, 1, 0])
        C_shift = QuasiCyclicCode([v, v3], 2, false)

        # Test circulant vectors constructor
        vecs = [matrix(F, 1, 4, [1, 0, 1, 1]), matrix(F, 1, 4, [0, 0, 0, 1]), 
                matrix(F, 1, 4, [1, 1, 1, 1]), matrix(F, 1, 4, [0, 0, 0, 0])]
        C_circ = QuasiCyclicCode(vecs, 2, true)
        
        # are_equivalent securely calls syndrome() which multiplies C_shift.G * C_circ.H^T
        # Because we fixed noncirculant expansion, this executes perfectly over GF(2)!
        @test are_equivalent(C_shift, C_circ)
    end

    @testset "Getters, Matrix Extraction, and Properties" begin
        F = Oscar.Nemo.Native.GF(2)
        S, x = polynomial_ring(F, :x)
        m = 5
        R, _ = residue_ring(S, x^m - 1)

        # Construct A directly in the residue ring
        A = matrix(R, 2, 3, [1+x x^2 0; x^3 1 x^4])
        C = QuasiCyclicCode(A, false)

        # Basic properties
        @test index(C) == 3
        @test expansion_factor(C) == 5
        @test polynomial_matrix(C) == A
        @test polynomial_matrix_type(C) == :G
        @test !is_single_generator(C)

        # Weight matrix evaluation
        W = weight_matrix(A)
        @test size(W) == (2, 3)
        @test W[1, 1] == 2 # 1+x has weight 2
        @test W[1, 2] == 1 # x^2 has weight 1
        @test base_matrix(A) == W
        @test protograph_matrix(A) == W

        # Extractors (now correctly return Matrix/Vectors of F)
        circs = circulants(C)
        @test length(circs) == 6

        # Non-circulant expansion natively evaluates the interleaved block matrix
        G_nc = noncirculant_generator_matrix(C)
        @test nrows(G_nc) == 10 # m * nr = 5 * 2
        @test ncols(G_nc) == 15 # m * l = 5 * 3
        
        # Test parity fallback
        @test ismissing(noncirculant_parity_check_matrix(C))
    end

    @testset "Algebraic Cycles (Fossorier's Condition)" begin
        F = Oscar.Nemo.Native.GF(2)
        S, x = polynomial_ring(F, :x)
        m = 7
        R, _ = residue_ring(S, x^m - 1)

        # Create an algebraic 4-cycle
        # Alternating sum: 0 - 1 + 3 - 2 = 0 mod 7 -> 4-cycle!
        A_cycle = matrix(R, 2, 2, [1 x; x^2 x^3])
        C_cycle = QuasiCyclicCode(A_cycle, true) # Set as parity check matrix

        E = shift_matrix(C_cycle)
        @test E == [0 1; 2 3]
        @test exponent_matrix(C_cycle) == E

        @test has_algebraic_4_cycle(C_cycle)
        @test has_algebraic_cycle(C_cycle, 4)

        # Multi-term polynomials should throw an ArgumentError in shift extraction
        A_bad = matrix(R, 1, 1, [1 + x])
        C_bad = QuasiCyclicCode(A_bad, false)
        @test_throws ArgumentError shift_matrix(C_bad)
        
        # Target length must be an even integer >= 4
        @test_throws ArgumentError has_algebraic_cycle(C_cycle, 3)
    end

    @testset "Algebraic Dimension and Components (Higher Fields)" begin
        # Non-binary fields strictly use standard GF(...)
        F3 = GF(3)
        S, x = polynomial_ring(F3, :x)
        m = 4 # Coprime to characteristic 3 to avoid nilpotent elements
        R, _ = residue_ring(S, x^m - 1)

        A = matrix(R, 1, 2, [1+x, x^2])
        C = QuasiCyclicCode(A, false)

        # Component extraction via Chinese Remainder Theorem
        comps = component_matrices(C)
        @test length(comps) > 0

        # Verify algebraic dimension solver perfectly matches the explicit rank
        alg_k = algebraic_dimension(C)
        @test alg_k > 0
        @test alg_k == dimension(C)

        # Test Parity Check algebraic dimension fallback
        C_H = QuasiCyclicCode(A, true)
        expected_k = length(C_H) - rank(lift(A))
        @test algebraic_dimension(C_H) == expected_k
    end
end