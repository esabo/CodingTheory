@testitem "ISD attacks" begin
    using Oscar, CodingTheory, Random

    seed_search!(seed) =
        Random.seed!(Random.default_rng(), rand(MersenneTwister(seed), UInt64))

    function check_decoding_result(C, found, received, target_w)
        F = field(C)
        H = parity_check_matrix(C)
        syndrome_target = H * matrix(F, C.n, 1, F.(received))

        @test !isempty(found)
        @test all(e -> count(!iszero, e) <= target_w, found)
        @test all(e -> H * matrix(F, C.n, 1, F.(e)) == syndrome_target, found)
    end

    C = HammingCode(2, 3)
    received = [1, 0, 0, 0, 0, 0, 0]

    attacks = [
        () -> Prange_attack(C, 1; w_recv=received, max_iters=50),
        () -> Lee_Brickell_attack(C, 1; w_recv=received, p=1, max_iters=50),
        () -> Leon_attack(C, 1; w_recv=received, p=0, l=1, max_iters=50),
        () -> Canteaut_Chabaud_attack(C, 1; w_recv=received, p=0, l=1,
            max_iters=50),
    ]

    for (seed, attack) in enumerate(attacks)
        seed_search!(seed)
        check_decoding_result(C, attack(), received, 1)
    end

    seed_search!(20)
    decoded = syndrome_decode(C, received, 1; alg=:Prange, confidence=0.999999,
        verbose=false)
    check_decoding_result(C, decoded, received, 1)
end

@testitem "ISD deterministic bounds" begin
    using Oscar, CodingTheory

    n, k, w = 20, 10, 3
    success_probability = Float64(binomial(big(n - k), big(w))) /
        Float64(binomial(big(n), big(w)))
    expected(rate) = ceil(Int, log1p(-rate) / log1p(-success_probability))

    @test required_ISD_iterations(:Prange, n, k, w, 0.5) == expected(0.5)
    @test required_ISD_iterations(:Prange, n, k, w, 0.9) == expected(0.9)
    @test required_ISD_iterations(:Prange, n, k, w, 0.99) == expected(0.99)
    @test required_ISD_iterations(:Prange, n, k, w, 0.5) <
        required_ISD_iterations(:Prange, n, k, w, 0.99)

    # the bound is the largest d with sum_{i=0}^{d-2} C(n-1,i) (q-1)^i < q^(n-k)
    function GV_reference(n, k, q)
        target = big(q)^(n - k)
        d = 1
        while sum(binomial(big(n - 1), big(i)) * big(q - 1)^i for i in 0:d - 1) <
              target
            d += 1
        end
        return d
    end

    for (n, k, q) in ((7, 4, 2), (7, 1, 2), (7, 7, 2), (15, 11, 2), (23, 12, 2),
                      (9, 5, 3))
        @test Gilbert_Varshamov_bound(n, k, q) == GV_reference(n, k, q)
    end

    # the bound is met by perfect codes and never exceeds the true distance
    @test Gilbert_Varshamov_bound(7, 4, 2) == 3
    @test Gilbert_Varshamov_bound(7, 1, 2) == 7
    @test Gilbert_Varshamov_bound(7, 7, 2) == 1
    @test Gilbert_Varshamov_bound(7, 4, 2) <= minimum_distance(HammingCode(2, 3))[1]
end
