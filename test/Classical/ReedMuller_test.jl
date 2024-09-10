@testset "Classical/ReedMuller.jl" begin
    # using Oscar, CodingTheory

    F = GF(2)
    # Huffman, Pless, p. 34
    # identity used for RM(1, 1)
    @test CodingTheory._Reed_Muller_generator_matrix(1, 1, true) == matrix(F,
        [1 0;
         0 1]);
    @test generator_matrix(ReedMullerCode(1, 2, true)) == matrix(F,
        [1 0 1 0;
         0 1 0 1;
         0 0 1 1]);
    @test generator_matrix(ReedMullerCode(1, 3, true)) == matrix(F,
        [1 0 1 0 1 0 1 0;
         0 1 0 1 0 1 0 1;
         0 0 1 1 0 0 1 1;
         0 0 0 0 1 1 1 1])
    @test generator_matrix(ReedMullerCode(2, 3, true)) == matrix(F,
        [1 0 0 0 1 0 0 0;
         0 1 0 0 0 1 0 0;
         0 0 1 0 0 0 1 0;
         0 0 0 1 0 0 0 1;
         0 0 0 0 1 0 1 0;
         0 0 0 0 0 1 0 1;
         0 0 0 0 0 0 1 1])

    # Ling & Xing, p. 119
    # other sources, using [1 1; 0 1] for RM(1, 1)
    @test CodingTheory._Reed_Muller_generator_matrix(1, 1) == matrix(F,
        [1 1;
         0 1]);
    @test generator_matrix(ReedMullerCode(1, 2)) == matrix(F,
        [1 1 1 1;
         0 1 0 1;
         0 0 1 1]);
    @test generator_matrix(ReedMullerCode(1, 3)) == matrix(F,
        [1 1 1 1 1 1 1 1;
         0 1 0 1 0 1 0 1;
         0 0 1 1 0 0 1 1;
         0 0 0 0 1 1 1 1])

    # if m is odd and r = (m - 1)/2 then RM(r, m) = RM((m - 1)/2, m) is self-dual
    # random m
    C = ReedMullerCode(2, 5)
    # length 2^m
    @test length(C) == 2^5
    @test is_self_dual(C)
    # RM(0, m) is the length 2^m repetition code
    @test are_equivalent(ReedMullerCode(0, 3), RepetitionCode(2, 8))
    # puncturing RM(1, m) and taking the even subcode is the simplex code S_m
    # random parameters
    C = ReedMullerCode(1, 4)
    pC = puncture(C, [1])
    epC = even_subcode(pC)
    S = SimplexCode(2, 4)
    
    @test are_equivalent(epC, S)
    C.d = missing
    # BUG syndrome trellis still weird
    # @test minimum_distance(C) == 8

    # the weight distribution of RM(1, m) is [[0, 1], [2^(m - 1), 2^(m + 1) - 2], [2^m, 1]]
    C.weight_enum = missing
    wt_dist = weight_distribution(C, alg = :auto, compact = true)
    @test wt_dist == [(2^4, 1), (2^3, 2^5 - 2), (0, 1)]

    # Reed-Muller codes are nested
    m = rand(3:6)
    r = rand(1:m - 2)
    C = ReedMullerCode(r, m)
    C2 = ReedMullerCode(r + 1, m)
    @test C ⊆ C2

    # # all weights of RM(r, m) are multiples of 2^(Int(ceil(m / r) - 1)
    # sup = support(C)
    # flag = true
    # for i in sup
    #     if !iszero(i % 2^(Int(ceil(m / r) - 1)))
    #         flag = false
    #         break
    #     end
    # end
    # @test flag == true

    # RM(m - 1, m) contains all vectors of even weight
    C = ReedMullerCode(m - 1, m)
    sup = support(C)
    flag = true
    for i in sup
        if isodd(i)
            flag = false
            break
        end
    end
    @test flag == true

    C = ReedMullerCode(2, 5)
    @test weight_distribution(C, alg = :auto, compact = true) == [(32, 1), (24, 620),
        (20, 13888), (16, 36518), (12, 13888), (8, 620), (0, 1)]
end
