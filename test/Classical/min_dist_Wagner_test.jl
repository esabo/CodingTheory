@testitem "Classical/min_dist_Wagner.jl" begin
    using Oscar, CodingTheory

    @testset "Wagner Meet-in-the-Middle Minimum Distance Solvers" begin
        F2 = Oscar.Nemo.Native.GF(2)
        F3 = Oscar.Nemo.Native.GF(3)
        F5 = Oscar.Nemo.Native.GF(5)

        # Helper to strip the LinearCode constructor's cached brute-force distance
        function reset_bounds!(C::AbstractLinearCode)
            C.l_bound = 1
            C.u_bound = C.n + 1
            C.d = missing
        end

        @testset "Binary Wagner (GF(2))" begin
            # Test 1: Hamming(7, 4) Code (Standard MitM split)
            G_hamming = matrix(F2, [
                1 0 0 0 0 1 1;
                0 1 0 0 1 0 1;
                0 0 1 0 1 1 0;
                0 0 0 1 1 1 1
            ])
            C_bin = LinearCode(G_hamming)
            
            reset_bounds!(C_bin)
            d_bin, witness_bin = CodingTheory._minimum_distance_wagner_mitm_binary(C_bin, max_d=5, verbose=false)
            @test d_bin == 3
            @test iszero(parity_check_matrix(C_bin) * transpose(witness_bin))
            
            # Test 2: Master Routing
            reset_bounds!(C_bin)
            d_master, witness_master = minimum_distance(C_bin, alg=:Wagner, verbose=false)
            @test d_master == 3
            @test iszero(parity_check_matrix(C_bin) * transpose(witness_master))

            # Test 3: The Asymmetric Right-Side Edge Case
            # H = [1 1 1 1 1 1]. The split is L={1,2,3}, R={4,5,6}.
            # Codewords have weight 2. A valid codeword is [0 0 0 1 1 0].
            # This forces the algorithm to handle w_L = 0 and w_R = 2.
            H_asym = matrix(F2, 1, 6, [1 1 1 1 1 1])
            C_asym = dual(LinearCode(H_asym)) # Dual creates the code where H_asym is the parity-check
            
            reset_bounds!(C_asym)
            d_asym, witness_asym = CodingTheory._minimum_distance_wagner_mitm_binary(C_asym, max_d=4, verbose=false)
            @test d_asym == 2
            @test iszero(parity_check_matrix(C_asym) * transpose(witness_asym))

            # Test 4: Early Exit / max_d limit
            reset_bounds!(C_bin)
            d_miss, witness_miss = CodingTheory._minimum_distance_wagner_mitm_binary(C_bin, max_d=2, verbose=false)
            @test d_miss == -1
            @test iszero(witness_miss)

            # Test 5: A Dense [31, k] BCH Code (Stress test on larger code and pruning)
            C = BCHCode(2, 31, 7) 
            # Clear out any cached brute-force distances
            reset_bounds!(C)
            # println("Running Wagner on a dense [31, $(C.k)] BCH Code...")
            d_bch, witness_bch = CodingTheory._minimum_distance_wagner_mitm_binary(C, max_d=10, verbose=true)
            @test d_bch == 8
            @test iszero(parity_check_matrix(C) * transpose(witness_bch))
        end

        @testset "Non-Binary Wagner (GF(q))" begin
            # Test 5: Ternary Repetition Code [5, 1, 5] (Basic non-binary)
            G_rep3 = matrix(F3, 1, 5, [1 1 1 1 1])
            C_ternary = LinearCode(G_rep3)
            
            reset_bounds!(C_ternary)
            d_ternary, witness_ternary = CodingTheory._minimum_distance_wagner_mitm_nonbinary(C_ternary, max_d=6, verbose=false)
            @test d_ternary == 5
            @test iszero(parity_check_matrix(C_ternary) * transpose(witness_ternary))
            
            reset_bounds!(C_ternary)
            d_master_ternary, witness_master_ternary = minimum_distance(C_ternary, alg=:Wagner, verbose=false)
            @test d_master_ternary == 5
            @test iszero(parity_check_matrix(C_ternary) * transpose(witness_master_ternary))
            
            # Test 6: GF(5) Code (Stress testing scalar loops and additive inverses)
            # This is a small [5, 3, 3] MDS-like code over GF(5)
            # Evaluates the solver's ability to map elements and negate targets properly.
            H_gf5 = matrix(F5, [
                1 1 1 1 1;
                1 2 3 4 0
            ])
            C_gf5 = dual(LinearCode(H_gf5))
            
            reset_bounds!(C_gf5)
            d_gf5, witness_gf5 = CodingTheory._minimum_distance_wagner_mitm_nonbinary(C_gf5, max_d=4, verbose=false)
            @test d_gf5 == 3
            @test iszero(parity_check_matrix(C_gf5) * transpose(witness_gf5))

            # Test 7: Early Exit / max_d limit on non-binary
            reset_bounds!(C_ternary)
            d_miss_t, witness_miss_t = CodingTheory._minimum_distance_wagner_mitm_nonbinary(C_ternary, max_d=3, verbose=false)
            @test d_miss_t == -1
            @test iszero(witness_miss_t)

            # Test 8
            println("\n=== Non-Binary Wagner Stress Test ===")
            # Constructs a Reed-Solomon Code: Length 10, Dimension 5. 
            # Because it's MDS, the distance is exactly n - k + 1 = 6.
            C_rs = ReedSolomonCode(11, 6)
            reset_bounds!(C_rs)
            d_rs, witness_rs = CodingTheory._minimum_distance_wagner_mitm_nonbinary(C_rs, max_d=8, verbose=true)
            @test d_rs == 6
            @test iszero(parity_check_matrix(C_rs) * transpose(witness_rs))
        end
    end
end
