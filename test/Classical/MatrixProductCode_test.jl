@testitem "Classical/MatrixProductCode.jl" begin
    using Oscar, CodingTheory, Random

    @testset "Constructors, Dimensions, and Getters" begin
        F = Oscar.Nemo.Native.GF(2)
        # C1: [7, 4, 3] Hamming Code
        C1 = HammingCode(2, 3) 
        
        # C2: [7, 1, 7] Repetition Code
        G2 = matrix(F, [1 1 1 1 1 1 1])
        C2 = LinearCode(G2)
        
        # Defining matrix A: 2 rows (s=2), 3 columns (l=3)
        A = matrix(F, [1 0 1; 0 1 1])
        
        MPC = MatrixProductCode([C1, C2], A)
        
        # Mathematical checks
        # length = l * n_sub = 3 * 7 = 21
        @test length(MPC) == 21
        
        # dimension = sum(k_i) = 4 + 1 = 5
        @test dimension(MPC) == 5
        @test field(MPC) == F
        
        # Getters
        @test constituent_codes(MPC) == [C1, C2]
        @test defining_matrix(MPC) == A
        
        # Generator Matrix Cache Expansion
        G_mpc = generator_matrix(MPC)
        @test nrows(G_mpc) == 5
        @test ncols(G_mpc) == 21
        
        # Standard form generation (forces cache update and rref extraction)
        G_stand = generator_matrix(MPC, true)
        @test G_stand[:, 1:5] == identity_matrix(F, 5)
    end
    
    @testset "Exceptions and Error Handling" begin
        F = Oscar.Nemo.Native.GF(2)
        C1 = HammingCode(2, 3)
        C2 = LinearCode(matrix(F, [1 1 1 1 1 1 1]))
        
        # Empty code vector rejection
        A_valid = matrix(F, [1 0; 0 1])
        @test_throws ArgumentError MatrixProductCode(AbstractLinearCode[], A_valid)
        
        # Zero defining matrix rejection
        @test_throws ArgumentError MatrixProductCode([C1, C2], zero_matrix(F, 2, 3))
        
        # Mismatched rows (A has 1 row, but we pass 2 codes)
        A_bad_rows = matrix(F, [1 1 1])
        @test_throws ArgumentError MatrixProductCode([C1, C2], A_bad_rows)
        
        # Rank-deficient matrix A (Rank 1, but s=2)
        A_bad_rank = matrix(F, [1 1 1; 1 1 1])
        @test_throws ArgumentError MatrixProductCode([C1, C2], A_bad_rank)
        
        # Mismatched base rings
        C_F4 = LinearCode(matrix(GF(4, :ω), [1 0 0; 0 1 0]))
        @test_throws ArgumentError MatrixProductCode([C1, C_F4], matrix(F, [1 0; 0 1]))
        
        # Mismatched sub-code lengths (C1 has length 7, C3 has length 3)
        C3 = LinearCode(matrix(F, [1 1 1]))
        @test_throws ArgumentError MatrixProductCode([C1, C3], matrix(F, [1 0; 0 1]))
    end
    
    @testset "Random Matrix Product Code Generator" begin
        F = Oscar.Nemo.Native.GF(2)
        C1 = HammingCode(2, 3)
        C2 = LinearCode(matrix(F, [1 1 1 1 1 1 1]))
        
        # Target l columns
        l = 4
        MPC_rand = RandomMatrixProductCode([C1, C2], l)
        
        # Properties
        @test length(MPC_rand) == 28 # 4 * 7
        @test dimension(MPC_rand) == 5 # 4 + 1
        
        # Validate that the random A matrix is actually s x l and full row rank
        A_rand = defining_matrix(MPC_rand)
        @test nrows(A_rand) == 2
        @test ncols(A_rand) == 4
        @test rank(A_rand) == 2
        
        # Catch domain error if we request fewer columns than rows (which makes full-rank impossible)
        @test_throws DomainError RandomMatrixProductCode([C1, C2], 1)
    end
end
