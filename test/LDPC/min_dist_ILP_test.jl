@testitem "LDPC/min_dist_ILP.jl" begin
    using Oscar, CodingTheory, JuMP, GLPK

    @testset "ILP Minimum Distance Solver" begin
        F = Oscar.Nemo.Native.GF(2)

        # Test 1: Sanity Check on Hamming(7, 4)
        # Small enough that the density of the parity-check matrix doesn't cause a branch-and-bound explosion.
        G_hamming = matrix(F, [
            1 0 0 0 0 1 1;
            0 1 0 0 1 0 1;
            0 0 1 0 1 1 0;
            0 0 0 1 1 1 1
        ])
        C = LinearCode(G_hamming)
        C.l_bound = 1
        C.u_bound = C.n + 1
        C.d = missing
        println("Running ILP solver on Hamming(7, 4)...")
        @test minimum_distance(C, alg = :ILP, verbose = true) == 3
        println("\n")

        # Test 2: The Repetition Code
        # This is a critical test. It ensures the sum(x) >= 1 constraint correctly 
        # forces the solver to find the all-ones vector instead of the all-zeros vector.
        G_rep = matrix(F, 1, 15, ones(Int, 15))
        C = LinearCode(G_rep)
        C.l_bound = 1
        C.u_bound = C.n + 1
        C.d = missing
        println("Running ILP solver on length 15 repetition code...")
        @test minimum_distance(C, alg = :ILP, verbose = true) == 15
        println("\n")

        # Test 3: Random Regular LDPC Code
        # We generate a deterministic, small (3, 6)-regular LDPC code.
        # This tests the ILP solver on the highly sparse matrices it is optimized for,
        # keeping `n` small enough (24) to ensure sub-second execution in the test suite.
        println("Running ILP solver on a random (3, 6)-regular LDPC code (n = 24)...")
        C = regular_LDPC_code(2, 24, 3, 6, seed = 42)
        C.l_bound = 1
        C.u_bound = C.n + 1
        C.d = missing
        d = minimum_distance(C, alg = :ILP, verbose = true)
        println("\n")
        @test d > 0
        C = LinearCode(C.H, true)
        C.l_bound = 1
        C.u_bound = C.n + 1
        C.d = missing
        @test d == minimum_distance(C, alg = :BZ, verbose = true)[1]

        # Test 4: Time Limit Trigger on Dense Code
        # We feed the solver a large, dense matrix where the LP relaxation is useless.
        # We enforce a 1.0 second time limit to ensure the timeout safety-valve works
        # and returns -1 instead of hanging the CI pipeline.
        println("Testing ILP time limit on a dense BCH code (n = 127)...")
        C = BCHCode(2, 127, 21, 1)
        C.l_bound = 1
        C.u_bound = C.n + 1
        C.d = missing
        d_timeout = minimum_distance(C, alg = :ILP, verbose = true, time_limit_sec = 1.0)
        # println("\n")
        # We expect the solver to hit the MOI.TIME_LIMIT status and safely return -1
        @test d_timeout == -1

        # the answer here is 4 because the stabilizers have wt 4 and this is a classical distance
        # S = ToricCode(13);
        # H_X = X_stabilizers(S);
        # C_X = LinearCode(H_X, true);
        # println("Running ILP on distance $(S.L) Toric Code X stabilizers (n = $(C_X.n))...")
        # d_exact = minimum_distance(C_X, alg = :ILP, verbose = true)
    end
end