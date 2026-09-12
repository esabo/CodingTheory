@testitem "Classical/min_dist_exact.jl" begin
    using Oscar, CodingTheory

    @testset "Exact Minimum Distance Solvers" begin
        F1 = Oscar.Nemo.Native.GF(2)
        F2 = GF(2)
        for F in (F1, F2)
            # Generator matrix for Hamming(7, 4)
            G = matrix(F, [
                1 0 0 0 0 1 1;
                0 1 0 0 1 0 1;
                0 0 1 0 1 1 0;
                0 0 0 1 1 1 1
            ])
            C = LinearCode(G)
            for info_set_alg ∈ (:auto, :Brouwer, :Zimmermann, :Bouyuklieva)
                C.l_bound = 1
                C.u_bound = C.n + 1
                C.d = missing
                println("Running Hamming code with information set algorithm $info_set_alg")
                d, witness = minimum_distance(C, alg = :BZ, info_set_alg = info_set_alg)
                @test minimum_distance(C)[1] == 3
                @test base_ring(witness) == C.F
                if !iszero(witness)
                    @test wt(witness) == C.d
                end
            end

            # Generator matrix for Extended Binary Golay
            G = matrix(F, [
                1 0 0 0 0 0 0 0 0 0 0 0  1 1 0 1 1 1 0 0 0 1 0 1;
                0 1 0 0 0 0 0 0 0 0 0 0  1 0 1 1 1 0 0 0 1 0 1 1;
                0 0 1 0 0 0 0 0 0 0 0 0  0 1 1 1 0 0 0 1 0 1 1 1;
                0 0 0 1 0 0 0 0 0 0 0 0  1 1 1 0 0 0 1 0 1 1 0 1;
                0 0 0 0 1 0 0 0 0 0 0 0  1 1 0 0 0 1 0 1 1 0 1 1;
                0 0 0 0 0 1 0 0 0 0 0 0  1 0 0 0 1 0 1 1 0 1 1 1;
                0 0 0 0 0 0 1 0 0 0 0 0  0 0 0 1 0 1 1 0 1 1 1 1;
                0 0 0 0 0 0 0 1 0 0 0 0  0 0 1 0 1 1 0 1 1 1 0 1;
                0 0 0 0 0 0 0 0 1 0 0 0  0 1 0 1 1 0 1 1 1 0 0 1;
                0 0 0 0 0 0 0 0 0 1 0 0  1 0 1 1 0 1 1 1 0 0 0 1;
                0 0 0 0 0 0 0 0 0 0 1 0  0 1 1 0 1 1 1 0 0 0 1 1;
                0 0 0 0 0 0 0 0 0 0 0 1  1 1 1 1 1 1 1 1 1 1 1 0
            ])
            C = LinearCode(G)
            for info_set_alg ∈ (:auto, :Brouwer, :Zimmermann, :Bouyuklieva)
                C.l_bound = 1
                C.u_bound = C.n + 1
                C.d = missing
                println("Running extended binary Golay code with information set algorithm $info_set_alg")
                d, witness = minimum_distance(C, alg = :BZ, info_set_alg = info_set_alg)
                @test minimum_distance(C)[1] == 8
                @test base_ring(witness) == C.F
                if !iszero(witness)
                    @test wt(witness) == C.d
                end
            end

            # generator matrix for repetition code
            G = matrix(F, 1, 15, ones(Int, 15))
            C = LinearCode(G)
            for info_set_alg ∈ (:auto, :Brouwer, :Zimmermann, :Bouyuklieva)
                C.l_bound = 1
                C.u_bound = C.n + 1
                C.d = missing
                println("Running repetition code with information set algorithm $info_set_alg")
                d, witness = minimum_distance(C, alg = :BZ, info_set_alg = info_set_alg)
                @test minimum_distance(C)[1] == 15
                @test base_ring(witness) == C.F
                if !iszero(witness)
                    @test wt(witness) == C.d
                end
            end
        end

        C = BCHCode(2, 127, 21, 1)
        info_set_alg = :Bouyuklieva
        # for info_set_alg ∈ (:Bouyuklieva, :Brouwer, :Zimmermann, :Chen) #  is too slow for this code
        C.l_bound = 1
        C.u_bound = C.n + 1
        C.d = missing
        println("Running BCH code with information set algorithm $info_set_alg")
        d, witness = minimum_distance(C, alg = :BZ, info_set_alg = info_set_alg, verbose = true)
        @test minimum_distance(C)[1] == 21
        @test base_ring(witness) == C.F
        if !iszero(witness)
            @test wt(witness) == C.d
        end
        println("\n")
        # end

        C = QuadraticResidueCode(2, 71)
        for info_set_alg ∈ (:auto, :Brouwer, :Zimmermann, :Bouyuklieva, :Chen)
            C.l_bound = 1
            C.u_bound = C.n + 1
            C.d = missing
            println("Running QR code with information set algorithm $info_set_alg")
            d, witness = minimum_distance(C, alg = :BZ, info_set_alg = info_set_alg, verbose = true)
            @test minimum_distance(C)[1] == 11
            @test base_ring(witness) == C.F
            if !iszero(witness)
                @test wt(witness) == C.d
            end
            println("\n")
        end

        C = ReedMullerCode(1, 6)
        for info_set_alg ∈ (:auto, :Brouwer, :Zimmermann, :Bouyuklieva)
            C.l_bound = 1
            C.u_bound = C.n + 1
            C.d = missing
            println("Running Reed-Muller code with information set algorithm $info_set_alg")
            d, witness = minimum_distance(C, alg = :BZ, info_set_alg = info_set_alg, verbose = true)
            @test minimum_distance(C)[1] == 32
            @test base_ring(witness) == C.F
            if !iszero(witness)
                @test wt(witness) == C.d
            end
            println("\n")
        end
    end
end
