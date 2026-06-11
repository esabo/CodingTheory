@testitem "Classical/trellis.jl" begin
    using Oscar, CodingTheory

    @testset "Trellis Construction and Invariants" begin
        F2 = Oscar.Nemo.Native.GF(2)
        
        # We use Hamming(7, 4) as our baseline. 
        # |C| = 2^4 = 16 codewords.
        G_hamming = matrix(F2, [
            1 0 0 0 0 1 1;
            0 1 0 0 1 0 1;
            0 0 1 0 1 1 0;
            0 0 0 1 1 1 1
        ])
        C = LinearCode(G_hamming)
        G_mat = Array(G_hamming)
        
        @testset "Trellis-Oriented Form (TOF) Uniqueness" begin
            G_tof = copy(G_mat)
            CodingTheory._make_trellis_oriented!(G_tof)
            L, R = CodingTheory._get_LR_indices(G_tof)
            
            # In Minimal Span Form, all non-trivial L and R indices must be strictly unique
            L_active = filter(x -> x <= C.n, L)
            R_active = filter(x -> x > 0, R)
            
            @test length(L_active) == length(unique(L_active))
            @test length(R_active) == length(unique(R_active))
        end
        
        @testset "Trellis Complexity & Wolf Bound" begin
            G_tof = copy(G_mat)
            CodingTheory._make_trellis_oriented!(G_tof)
            L, R = CodingTheory._get_LR_indices(G_tof)
            past, future = CodingTheory.past_future_profiles(L, R, C.n)
            
            V_profile = CodingTheory.vertex_counts(C.k, C.n, past, future)
            
            # Wolf Bound: dimension of V_i <= min(k, n-k)
            max_state_dim = min(C.k, C.n - C.k) 
            @test maximum(V_profile) <= max_state_dim
            
            # The Trellis must always start and merge at a single state (dimension 0)
            @test V_profile[1] == 0
            @test V_profile[end] == 0
        end
        
        @testset "CWE Path Conservation (|C| = q^k)" begin
            G_tof = copy(G_mat)
            CodingTheory._make_trellis_oriented!(G_tof)
            
            # 1. Test the Primal Generator Trellis Product
            cwe_dict_gen = CodingTheory._CWE_classical_TP_sectionalized(G_tof) 
            total_paths_gen = sum(values(cwe_dict_gen))
            
            # The sum of all compositions must exactly equal the total number of codewords
            @test total_paths_gen == 2^4
            
            # 2. Test the Syndrome Trellis (Parity-Check)
            H_mat = Array(parity_check_matrix(C))
            CodingTheory._make_trellis_oriented!(H_mat)
            cwe_dict_syn = CodingTheory._CWE_classical_syndrome_sectionalized(H_mat)
            
            total_paths_syn = sum(values(cwe_dict_syn))
            @test total_paths_syn == 2^4
        end

        @testset "Optimal Sectionalization" begin
            F2 = Oscar.Nemo.Native.GF(2)
            
            @testset "Linear Sectionalization" begin
                # Using Hamming(7,4) again
                G_lin = matrix(F2, [
                    1 0 0 0 0 1 1;
                    0 1 0 0 1 0 1;
                    0 0 1 0 1 1 0;
                    0 0 0 1 1 1 1
                ])
                CodingTheory._make_trellis_oriented!(G_lin)
                
                bounds = optimal_sectionalization(G_lin, 2, type=:linear)
                
                # Boundaries must start at 0, end at n, and be strictly increasing
                @test bounds[1] == 0
                @test bounds[end] == 7
                @test issorted(bounds, lt=<=) # strictly increasing
            end
            
            @testset "Quasi-Cyclic (QC) Sectionalization" begin
                # Mock a k=2, n=6 matrix (2 blocks of size p=3)
                G_qc = matrix(F2, [
                    1 1 0 1 1 0;
                    0 1 1 0 1 1
                ])
                CodingTheory._make_trellis_oriented!(G_qc)
                
                bounds = optimal_sectionalization(G_qc, 2, type=:QC, p=3)
                
                @test bounds[1] == 0
                @test bounds[end] == 6
                @test issorted(bounds, lt=<=)
                
                # Missing 'p' kwarg should throw an AssertionError
                @test_throws AssertionError optimal_sectionalization(G_qc, 2, type=:QC)
            end
            
            @testset "2D Cyclic Sectionalization" begin
                # Mock a k=4, n=12 matrix 
                # Let grid_x = 4, grid_y = 3 (4x3 = 12)
                # Let unit cell p_x = 2, p_y = 1
                G_2d = matrix(F2, [
                    1 1 1 0 0 0 0 0 0 0 0 0;
                    0 0 0 1 1 1 0 0 0 0 0 0;
                    0 0 0 0 0 0 1 1 1 0 0 0;
                    0 0 0 0 0 0 0 0 0 1 1 1
                ])
                CodingTheory._make_trellis_oriented!(G_2d)
                
                bounds = optimal_sectionalization(G_2d, 2, type=:twoD, p_x=2, p_y=1, grid_x=4, grid_y=3)
                
                @test bounds[1] == 0
                @test bounds[end] == 12
                @test issorted(bounds, lt=<=)
                
                # Missing kwargs should throw an AssertionError
                @test_throws AssertionError optimal_sectionalization(G_2d, 2, type=:twoD, p_x=2)
            end
            
            @testset "Router Error Handling" begin
                G_err = matrix(F2, [1 1; 0 1])
                # Unknown type should throw an ArgumentError
                @test_throws ArgumentError optimal_sectionalization(G_err, 2, type=:UnknownType)
            end
        end
    end
end

@testitem "Classical/trellis_distance.jl" begin
    using Oscar, CodingTheory

    @testset "Trellis Distance & Weight Distributions" begin
        F2 = Oscar.Nemo.Native.GF(2)
        
        # Hamming(7, 4) has a known weight enumerator: 
        # 1 word of weight 0, 7 of weight 3, 7 of weight 4, 1 of weight 7.
        G_hamming = matrix(F2, [
            1 0 0 0 0 1 1;
            0 1 0 0 1 0 1;
            0 0 1 0 1 1 0;
            0 0 0 1 1 1 1
        ])
        C = LinearCode(G_hamming)
        
        @testset "Weight Distribution" begin
            hwe_dict = weight_distribution_trellis(C)
            
            @test hwe_dict[0] == 1
            @test hwe_dict[3] == 7
            @test hwe_dict[4] == 7
            @test hwe_dict[7] == 1
            @test sum(values(hwe_dict)) == 16
        end
        
        @testset "Minimum Distance Extraction" begin
            C.d = missing # clear cache
            
            # Test the explicit trellis router (Pure Viterbi)
            d_trellis = minimum_distance(C, alg=:trellis, verbose=false)
            @test d_trellis == 3
            
            C.d = missing
            
            # Test the explicit hybrid bridge
            # (Forcing max_span=1 to ensure the BZ DFS bridge actually engages)
            d_hybrid = CodingTheory._minimum_distance_hybrid(C, max_span=1, verbose=false)
            @test d_hybrid == 3
        end
        
        @testset "MacWilliams Identity" begin
            # The dual of Hamming(7,4) is the Simplex(7,3) code.
            # Its weights are all exactly 4 (except the zero vector).
            C_dual = dual(C)
            hwe_dual = weight_distribution_trellis(C_dual)
            
            @test hwe_dual[0] == 1
            @test hwe_dual[4] == 7
            @test sum(values(hwe_dual)) == 8
            
            # Test the MacWilliams HWE transform back to primal
            primal_reconstructed = Macwilliams_HWE_transform(hwe_dual, C.n, C.k, 2)
            
            # It should perfectly reconstruct the Hamming(7,4) weights
            @test primal_reconstructed[3] == 7
            @test primal_reconstructed[4] == 7
        end
    end
end

@testitem "Classical/trellis_robust.jl" begin
    using Oscar, CodingTheory, Random

    @testset "Robust Trellis Invariants & Distance" begin
        F2 = Oscar.Nemo.Native.GF(2)
        F3 = Oscar.Nemo.Native.GF(3)
        
        @testset "TOF Uniqueness on Random Matrices" begin
            # Test over GF(2) and GF(3) to ensure pivot normalization works
            for F in [F2, F3]
                for _ in 1:10
                    # Random 10x20 matrix
                    G_rand = matrix(F, rand(0:(Int(order(F))-1), 10, 20))
                    G_arr = Array(G_rand)
                    
                    CodingTheory._make_trellis_oriented!(G_arr)
                    L, R = CodingTheory._get_LR_indices(G_arr)
                    
                    # Filter out the empty rows (L = n+1, R = 0)
                    L_active = filter(x -> x <= 20, L)
                    R_active = filter(x -> x > 0, R)
                    
                    @test length(L_active) == length(unique(L_active))
                    @test length(R_active) == length(unique(R_active))
                end
            end
        end
        
        @testset "CWE Path Conservation on Mid-Size Code" begin
            # Random 8x16 code over GF(2). |C| = 256.
            G_rand = matrix(F2, rand(0:1, 8, 16))
            # Ensure full rank for exact path counts
            while rank(G_rand) < 8
                G_rand = matrix(F2, rand(0:1, 8, 16))
            end
            
            C_rand = LinearCode(G_rand)
            
            # Primal Trellis Product
            cwe_primal = CodingTheory._CWE_classical_TP_sectionalized(Array(generator_matrix(C_rand)))
            @test sum(values(cwe_primal)) == 2^8
            
            # Dual Syndrome Trellis
            H_arr = Array(parity_check_matrix(C_rand))
            CodingTheory._make_trellis_oriented!(H_arr)
            cwe_dual = CodingTheory._CWE_classical_syndrome_sectionalized(H_arr)
            @test sum(values(cwe_dual)) == 2^8
        end
        
        @testset "Hybrid Bridge Stress Test (Hamming 15,11)" begin
            # Generate Hamming(15, 11). d = 3.
            H_15 = matrix(F2, [
                0 0 0 0 0 0 0 1 1 1 1 1 1 1 1;
                0 0 0 1 1 1 1 0 0 0 0 1 1 1 1;
                0 1 1 0 0 1 1 0 0 1 1 0 0 1 1;
                1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
            ])
            C_15 = LinearCode(H_15, true)
            
            # Pure Trellis
            d_pure = minimum_distance(C_15, alg=:trellis, verbose=false)
            @test d_pure == 3
            
            # Force the DFS to handle a large death zone (max_span = 3)
            # The peak span for this code is usually 6-8 depending on permutation.
            C_15.d = missing
            d_hybrid = CodingTheory._minimum_distance_hybrid(C_15, max_span=3, num_trials=5, verbose=false)
            @test d_hybrid == 3
        end
    end
end

@testitem "Classical/weight_enumerator_test.jl" begin
    using Oscar, CodingTheory

    @testset "Known Code Weight Enumerators (HWE)" begin
        
        @testset "Binary Repetition Code" begin
            # RepetitionCode(2, 5) -> [5, 1, 5] [cite: 6]
            C = RepetitionCode(2, 5)
            hwe = weight_distribution_trellis(C)
            
            # The only codewords are the all-zero and all-one vectors
            @test hwe[0] == 1
            @test hwe[5] == 1
            @test length(hwe) == 2
        end
        
        @testset "Binary Simplex Code" begin
            # SimplexCode(2, 3) -> [7, 3, 4] [cite: 19, 20, 25]
            C = SimplexCode(2, 3)
            hwe = weight_distribution_trellis(C)
            
            # All 2^3 - 1 = 7 non-zero codewords must have weight 2^(3-1) = 4 [cite: 25]
            @test hwe[0] == 1
            @test hwe[4] == 7
            @test length(hwe) == 2
            @test sum(values(hwe)) == 2^3
        end

        @testset "Tetra Code over GF(3)" begin
            # TetraCode() -> [4, 2, 3] over GF(3) [cite: 17, 18]
            C = TetraCode()
            hwe = weight_distribution_trellis(C)
            
            # From the theoretical complete weight enumerator provided [cite: 18]
            # It has exactly 1 word of weight 0, and 8 words of weight 3.
            @test hwe[0] == 1
            @test hwe[3] == 8
            @test length(hwe) == 2
            @test sum(values(hwe)) == 3^2
        end

        @testset "Extended Binary Golay Code" begin
            # ExtendedGolayCode(2) -> [24, 12, 8] [cite: 26, 38]
            C = ExtendedGolayCode(2)
            hwe = weight_distribution_trellis(C)
            
            # From the known weight enumerator polynomial for the Golay code [cite: 31]
            @test hwe[0] == 1
            @test hwe[8] == 759
            @test hwe[12] == 2576
            @test hwe[16] == 759
            @test hwe[24] == 1
            
            # Ensure no other weights leaked into the dictionary
            @test length(hwe) == 5
            @test sum(values(hwe)) == 2^12
        end

        @testset "First-Order Reed-Muller / Walsh-Hadamard" begin
            # RM(1, 3) is a [8, 4, 4] code. 
            # In your notes, these are equivalent to Walsh-Hadamard codes[cite: 40].
            # Generating it directly via a known matrix if ReedMullerCode isn't explicitly defined.
            F2 = Oscar.Nemo.Native.GF(2)
            G_rm = matrix(F2, [
                1 1 1 1 1 1 1 1;
                0 0 0 0 1 1 1 1;
                0 0 1 1 0 0 1 1;
                0 1 0 1 0 1 0 1
            ])
            C = LinearCode(G_rm)
            hwe = weight_distribution_trellis(C)
            
            # The weights for RM(1, m) are 0, 2^(m-1), and 2^m.
            @test hwe[0] == 1
            @test hwe[4] == 14
            @test hwe[8] == 1
            @test sum(values(hwe)) == 16
        end

        @testset "MDS Reed-Solomon Code (Theoretical Formula Match)" begin
            # Standard generic RS generator over GF(8) or GF(4) 
            # Let's use a simple [3, 2, 2] MDS code over GF(4) to test the mathematical identity
            F4 = GF(2, 2, :a)
            a = gen(F4)
            # Generator for a [n=3, k=2] RS code
            G_rs = matrix(F4, [
                1 1 1;
                1 a a^2
            ])
            C = LinearCode(G_rs)
            hwe = weight_distribution_trellis(C)
            
            # For an MDS code, the number of codewords of weight w >= d is governed by:
            # A_w = binom(n, w) * (q - 1) * sum( (-1)^j * binom(w-1, j) * q^(w - d - j) )
            function mds_weight(n, k, q, w)
                d = n - k + 1
                if w < d return 0 end
                
                term_sum = 0
                for j in 0:(w - d)
                    term_sum += (-1)^j * binomial(w - 1, j) * q^(w - d - j)
                end
                return binomial(n, w) * (q - 1) * term_sum
            end
            
            @test hwe[0] == 1
            @test hwe[2] == mds_weight(3, 2, 4, 2)
            @test hwe[3] == mds_weight(3, 2, 4, 3)
            @test sum(values(hwe)) == 4^2
        end
    end
end

