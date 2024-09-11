@testitem "Classical/misc_known_codes.jl" begin
    using Oscar, CodingTheory

    @testset "Hexacode" begin
        # Right now these are hardcoded in the function Hexacode(), but
        # testing in case we move to having some of it done automatically
        H = Hexacode()
        ω = gen(H.F)
        @test H.n == 6
        @test H.k == 3
        @test H.d == H.l_bound == H.u_bound == 4
        @test H.H == matrix(H.F, [1 ω ω 1 0 0; ω 1 ω 0 1 0; ω ω 1 0 0 1])
    end

    @testset "HammingCode and SimplexCode" begin
        R, (x, y) = polynomial_ring(Nemo.ZZ, [:x, :y])
        # Hamming codes
        # Tetra code is Hammingcode(3, 2)
        # random Hamming code
        F = Oscar.Nemo.Native.GF(2)
        C = HammingCode(2, 7)
        col = rand(1:length(C))
        # columns are 1, 2, ... 2^r - 1 written as binary numerals
        @test parity_check_matrix(C)[:, col:col] == matrix(F, length(C) -
                dimension(C), 1, reverse(digits(col, base=2, pad=7)))
        # should be [2^r - 1, 2^r - 1 - r, 3]
        @test length(C) == 2^7 - 1
        # TODO: fix Oscar imports and remove all CodingTheory.'s here
        @test CodingTheory.dimension(C) == 2^7 - 1 - 7
        C.d = missing
        #@test minimum_distance(C) == 3
        C = HammingCode(2, 3)
        ham_WE = weight_enumerator(C, type = :Hamming)
        @test polynomial(ham_WE) == x^7 + 7*x^3*y^4 + 7*x^4*y^3 + y^7
        n = length(C)
        C.weight_enum = missing
        @test polynomial(ham_WE) == divexact((x + y)^n + n*(x + y)^div(n - 1, 2)*(y - x)^div(n + 1, 2), n + 1)

        # simplex codes
        # random simplex code
        C = SimplexCode(2, 4)
        known = C.weight_enum
        # C.weight_enum = missing
        # HWEbf = weight_enumerator(C, type = :Hamming)
        # C.weight_enum = missing
        # HWEtrellis = weight_enumerator(C, type = :Hamming, "trellis")
        # @test CWEtoHWE(known) == HWEbf
        # @test HWEbf == HWEtrellis
        # all nonzero codewords have weights q^{r - 1}
        # flag = true
        # for exps in [collect(exponent_vectors(polynomial(HWEtrellis)))[i][1]
        #                 for i in 1:length(polynomial(HWEtrellis))]
        #         if !iszero(exps % 2^(4 - 1))
        #                 flag = false
        #                 break
        #         end
        # end
        # @test flag == true
        @test length(C) == 2^4 - 1
        @test CodingTheory.dimension(C) == 4
        C = SimplexCode(2, 3)
        @test MacWilliams_identity(C, weight_enumerator(C, type = :Hamming, alg = :bruteforce)) == ham_WE
    end

    @testset "Golay code" begin
        R, (x, y) = polynomial_ring(Nemo.ZZ, [:x, :y])
        # Golay codes
        # TODO: test extend for the ternary Golay code
        C = ExtendedGolayCode(2)
        @test is_self_dual(C)
        C.weight_enum = missing
        @test polynomial(weight_enumerator(C, type = :Hamming)) == y^24 + 759*x^8*y^16 + 2576*x^12*y^12 + 759*x^16*y^8 + x^24
        C = GolayCode(2)
        C.weight_enum = missing
        @test polynomial(weight_enumerator(C, type = :Hamming)) == y^23 + 253*x^7*y^16 +
                506*x^8*y^15 + 1288*x^11*y^12 + 1288*x^12*y^11 + 506*x^15*y^8 + 253*x^16*y^7 + x^23
        C = ExtendedGolayCode(3)
        @test is_self_dual(C)
        # well-known weight enumerators
        C.weight_enum = missing
        @test polynomial(weight_enumerator(C, type = :Hamming)) == y^12 + 264*x^6*y^6 + 440*x^9*y^3 + 24*x^12
        C = GolayCode(3)
        @test polynomial(weight_enumerator(C, type = :Hamming)) == y^11 + 132*x^5*y^6 + 132*x^6*y^5 + 330*x^8*y^3 + 110*x^9*y^2 + 24*x^11
        # cyclic code with generator polynomial g(x) = -1 + x^2 - x^3 + x^4 + x^5
        # and idempotent e(x) = -(x^2 + x^6 + x^7 + x^8 + x^10)
        # should be eqivalent to the [11, 6, 5] Golay code (maybe permutation?)

        # Huffman & Pless, p33, exercise 61d
        C = ExtendedGolayCode(3)
        C2 = extend(puncture(C, 7), 7)
        T = identity_matrix(C.F, 12)
        T[7,7] = C.F(-1)
        C3 = LinearCode(C2.G * T)
        @test are_equivalent(C, C3)
    end

    @testset "Tetra code" begin
        # tetra code
        C = TetraCode()
        C.weight_enum = missing
        CWE = polynomial(weight_enumerator(C, type = :complete))
        vars = gens(parent(CWE))
        @test CWE == vars[1]^4 + vars[1]*vars[2]^3 + 3*vars[1]*vars[2]^2*vars[3] +
                3*vars[1]*vars[2]*vars[3]^2 + vars[1]*vars[3]^3
    end
        # Hadamard code
        # the dual code of the Hamming code is the shortened Hadamard code
        # equivalent to RM(1, m)

end
