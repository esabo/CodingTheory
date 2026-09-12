@testitem "Probabilistic minimum-distance searches" begin
    using Oscar, CodingTheory, Random

    seed_search!(seed) =
        Random.seed!(Random.default_rng(), rand(MersenneTwister(seed), UInt64))

    function check_upper_bound(C, true_distance, result)
        distance, witness = result
        entries = witness isa AbstractVector ? witness : vec(Array(witness))
        F = field(C)
        word = matrix(F, C.n, 1, F.(entries))

        @test distance >= true_distance
        @test distance == count(!iszero, entries)
        @test distance > 0
        @test iszero(parity_check_matrix(C) * word)
    end

    G = generator_matrix(HammingCode(2, 3))
    searches = [
        C -> probabilistic_minimum_distance_prange(C; confidence=0.99, verbose=false),
        C -> probabilistic_minimum_distance_lee_brickell(C; confidence=0.99, p=1,
            verbose=false),
        C -> probabilistic_minimum_distance_leon(C; confidence=0.99, p=1, l=1,
            verbose=false),
        C -> probabilistic_minimum_distance_stern_DOOM(C; confidence=0.99, p=1, l=1,
            verbose=false),
    ]

    for (seed, search) in enumerate(searches)
        C = LinearCode(G)
        seed_search!(100 + seed)
        check_upper_bound(C, 3, search(C))
    end
end

@testitem "Heuristic minimum-distance searches" begin
    using Oscar, CodingTheory, Random

    seed_search!(seed) =
        Random.seed!(Random.default_rng(), rand(MersenneTwister(seed), UInt64))

    function check_upper_bound(C, true_distance, result)
        distance, witness = result
        entries = vec(Array(witness))
        word = matrix(field(C), C.n, 1, entries)

        @test distance >= true_distance
        @test distance == count(!iszero, entries)
        @test distance > 0
        @test iszero(parity_check_matrix(C) * word)
    end

    G = generator_matrix(HammingCode(2, 3))
    searches = [
        C -> heuristic_minimum_distance_ga(C; max_gens=3, pop_size=10),
        C -> heuristic_minimum_distance_aco(C; y_max=3, m_ants=8),
        C -> heuristic_minimum_distance_irons(C; num_iters=20, p_max=2),
        C -> heuristic_minimum_distance_nncs(C),
        C -> heuristic_minimum_distance_gga_order(C; max_gens=3, pop_size=4),
    ]

    for (seed, search) in enumerate(searches)
        C = LinearCode(G)
        seed_search!(200 + seed)
        check_upper_bound(C, 3, search(C))
    end
end

@testitem "Minimum-distance search utilities" begin
    using Oscar, CodingTheory

    for n in 1:20, k in 0:n
        @test logbinomial(n, k) ≈ log(Float64(binomial(n, k))) atol=1e-12
    end
    @test logbinomial(12, 3) < logbinomial(13, 3)
    @test logbinomial(12, -1) == -Inf
    @test logbinomial(12, 13) == -Inf

    C = HammingCode(2, 3)
    G = generator_matrix(C)
    matrices, permutations, ranks = information_sets(G, :Edmonds)
    @test !isempty(matrices)
    @test length(matrices) == length(permutations) == length(ranks)
    for (permutation, r) in zip(permutations, ranks)
        @test rank(G[:, permutation[1:r]]) == r
    end
    @test heuristic_info_set_selection(C) == :Edmonds

    RM = ReedMullerCode(1, 3)
    H_RM = parity_check_matrix(RM)
    automorphisms = generate_automorphisms(RM)
    @test !isempty(automorphisms)
    for permutation in automorphisms
        @test sort(permutation) == collect(1:RM.n)
        @test iszero(H_RM * transpose(generator_matrix(RM)[:, permutation]))
    end

    distance, witness = minimum_distance_zssmp(RM; verbose=false)
    @test distance >= 4
    @test distance == count(!iszero, witness)
    @test iszero(H_RM * transpose(witness))
end
