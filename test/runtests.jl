using Test

# these are subject to change as they develop, let me know when it breaks
# TODO: add hypergraph, LDPC
@testset "Types" begin
    using CodingTheory

    @test AbstractLinearCode <: AbstractCode
    @test AbstractCyclicCode <: AbstractLinearCode
    @test AbstractBCHCode <: AbstractCyclicCode
    @test AbstractReedSolomonCode <: AbstractBCHCode
    @test AbstractGeneralizedReedSolomonCode <: AbstractLinearCode
    @test AbstractAdditiveCode <: AbstractCode
    @test AbstractStabilizerCode <: AbstractAdditiveCode
    @test AbstractCSSCode <: AbstractStabilizerCode

    # check the hierarchy
    @test AbstractReedSolomonCode <: AbstractCyclicCode
    @test AbstractReedSolomonCode <: AbstractLinearCode
end

@testset "utils.jl" begin
    using CodingTheory

    # example: Betten et al
    # Golay code G_23
    qres, _ = quadraticresidues(2, 23)
    @test qres == [1, 2, 3, 4, 6, 8, 9, 12, 13, 16, 18]

    # example: Betten et al
    # ternary Golary code G_11
    qres, _ = quadraticresidues(3, 11)
    @test qres == [1, 3, 4, 5, 9]

end

@testset "linearcode.jl" begin
    using Nemo, CodingTheory

    F, _ = FiniteField(2, 1, "α");
    G = matrix(F, [1 0 0 0 0 1 1;
           0 1 0 0 1 0 1;
           0 0 1 0 1 1 0;
           0 0 0 1 1 1 1]);
    C = LinearCode(G);
    @test field(C) == F
    @test length(C) == 7
    @test rank(G) == dimension(C)
    @test cardinality(C) == BigInt(2)^4
    @test dimension(C) == 4
    @test rate(C) == 4 / 7
    @test ismissing(C.d)
    setminimumdistance!(C, 3)
    @test minimumdistance(C) == 3
    @test numbercorrectableerrors(C) == 1
    @test G == generatormatrix(C)
    @test G == originalgeneratormatrix(C)
    H = paritycheckmatrix(C)
    @test iszero(G * transpose(H))
    @test iszero(H * transpose(G))
    @test C ⊆ C
    D = dual(C)
    @test !(C ⊆ D)
    @test !issubcode(C, D)
    @test !isequivalent(C, D)
    @test isequivalent(C, C)
    @test !isselfdual(C)
    @test !isselforthogonal(C)
    cw = matrix(F, [1 0 0 0 0 1 1]);
    @test encode(C.G[:, 1], C) == cw
    # these v's are C.G[:, 1], just testing different formats
    v = [1, 0, 0, 0];
    @test encode(v, C) == cw
    v2 = [1; 0; 0; 0];
    @test encode(v2, C) == cw
    # this vector is the first row of the generator matrix and should
    # therefore have zero syndrome
    v = [1, 0, 0, 0, 0, 1, 1];
    @test iszero(syndrome(v, C))
    v = [1; 0; 0; 0; 0; 1; 1];
    @test iszero(syndrome(v, C))

    # lower rank test
    GandG = vcat(G, G);
    CGandG = LinearCode(G);
    @test rank(GandG) == dimension(CGandG)
    @test G == generatormatrix(CGandG)
    # this fails, is of wrong size, probably a transpose mistake
    # @test GandG == originalgeneratormatrix(CGandG)

    # puncturing examples from Huffman/Pless
    G = matrix(F, [1 1 0 0 0; 0 0 1 1 1])
    C = LinearCode(G)
    @test generatormatrix(puncture(C, [1])) == matrix(F, [1 0 0 0; 0 1 1 1])
    @test generatormatrix(puncture(C, [5])) == matrix(F, [1 1 0 0; 0 0 1 1])
    G = matrix(F, [1 0 0 0; 0 1 1 1])
    C = LinearCode(G)
    @test generatormatrix(puncture(C, [1])) == matrix(F, [1 1 1])
    @test generatormatrix(puncture(C, [4])) == matrix(F, [1 0 0; 0 1 1])

    # extending examples from Huffman/Pless
    C = TetraCode()
    exC = extend(C)
    @test generatormatrix(exC) == matrix(field(C), [1 0 1 1 0; 0 1 1 -1 -1])
    @test paritycheckmatrix(exC) == matrix(field(C), [1 1 1 1 1; -1 -1 1 0 0; -1 1 0 1 0])
    G = matrix(F, [1 1 0 0 1; 0 0 1 1 0])
    C = LinearCode(G)
    @test generatormatrix(extend(puncture(C, [5]))) == matrix(F, [1 1 0 0 0; 0 0 1 1 0])
    G = matrix(F, [1 0 0 1 1 1; 0 1 0 1 1 1; 0 0 1 1 1 1])
    C = LinearCode(G)
    @test generatormatrix(puncture(C, [5, 6])) == matrix(F, [1 0 0 1; 0 1 0 1; 0 0 1 1])

    # shortening examples from Huffman/Pless
    shC = shorten(C, [5, 6])
    shCtest = LinearCode(matrix(F, [1 0 1 0; 0 1 1 0]))
    @test isequivalent(shC, shCtest)

    # need to fix Hermitian stuff, can't find Hermitianconjugatematrix
    # C = Hexacode()
    # D = Hermitiandual(C)
    # @test isequivalent(C, D)



    # missing so far in tests:
    # expurgate, augment, lengthen, uuplusv, subcode,
    # juxtaposition, constructionX, constructionX3, upluswvpluswuplusvplusw,
    # expandedcode, entrywiseproductcode, evensubcode

    # subfield subcode example from Huffman/Pless
    K, ω = FiniteField(2, 2, "ω");
    C = Hexacode();
    dualbasis = [ω^2, K(1)]; # dual?
    CF2 = subfieldsubcode(C, F, dualbasis)
    @test isequivalent(CF2, RepetitionCode(2, 6))

    # # test Delsarte's theorem
    # @test isequivalent(CF2, dual(tracecode(dual(C), F)))

    # to test perms
    # C = HammingCode(2, 3)
    # S7 = SymmetricGroup(7);
    # σ = S7([3, 2, 1, 4, 5, 6, 7])
    # permutecode(C, σ)
    # # or
    # permutecode(C, [1, 3, 2, 4, 5, 6, 7])


    # "On the Schur Product of Vector Spaces over Finite Fields"
    # Christiaan Koster
    # Lemma 14: If C is cyclic and dim(C) > (1/2)(n + 1), then C * C = F^n

    # simplex code itself has dimension k(k + 1)/2
    #
end

@testset "ReedMuller.jl" begin
    using Nemo, CodingTheory

    F, _ = FiniteField(2, 1, "α");
    @test generatormatrix(ReedMullerCode(2, 1, 2)) == matrix(F,
        [1 0 1 0;
         0 1 0 1;
         0 0 1 1]);
    @test generatormatrix(ReedMullerCode(2, 1, 3)) == matrix(F,
        [1 0 1 0 1 0 1 0;
         0 1 0 1 0 1 0 1;
         0 0 1 1 0 0 1 1;
         0 0 0 0 1 1 1 1])
    @test generatormatrix(ReedMullerCode(2, 2, 3)) == matrix(F,
        [1 0 0 0 1 0 0 0;
         0 1 0 0 0 1 0 0;
         0 0 1 0 0 0 1 0;
         0 0 0 1 0 0 0 1;
         0 0 0 0 1 0 1 0;
         0 0 0 0 0 1 0 1;
         0 0 0 0 0 0 1 1])

    # if m is odd and r = (m - 1)/2 then RM(r, m) = RM((m - 1)/2, m) is self-dual
    # random m
    C = ReedMullerCode(2, 2, 5)
    # length 2^m
    @test length(C) == 2^5
    @test isselfdual(C)
    # RM(0, m) is the length 2^m repetition code
    @test isequivalent(ReedMullerCode(2, 0, 3), RepetitionCode(2, 8))
    # puncturing RM(1, m) and taking the even subcode is the simplex code S_m
    # random parameters
    C = ReedMullerCode(2, 1, 4)
    pC = puncture(C, [1])
    epC = evensubcode(pC)
    S = SimplexCode(2, 4)
    @test isequivalent(epC, S)
    C.d = missing
    @test minimumdistance(C) == 8

    # the weight distribution of RM(1, m) is [[0, 1], [2^(m - 1), 2^(m + 1) - 2], [2^m, 1]]
    C.weightenum = missing
    wtdist = weightdistribution(C, "auto", "compact")
    @test wtdist == [(2^4, 1), (2^3, 2^5 - 2), (0, 1)]

    # Reed-Muller codes are nested
    m = rand(3:6)
    r = rand(1:m - 2)
    C = ReedMullerCode(2, r, m)
    C2 = ReedMullerCode(2, r + 1, m)
    @test C ⊆ C2

    # all weights of RM(r, m) are multiples of 2^(Int(ceil(m / r) - 1)
    sup = support(C)
    flag = true
    for i in sup
        if !iszero(i % 2^(Int(ceil(m / r) - 1)))
            flag = false
            break
        end
    end
    @test flag == true

    # RM(m - 1, m) contains all vectors of even weight
    C = ReedMullerCode(2, m - 1, m)
    sup = support(C)
    flag = true
    for i in sup
        if isodd(i)
            flag = false
            break
        end
    end
    @test flag == true

    C = ReedMullerCode(2, 2, 5)
    @test weightdistribution(C, "auto", "compact") == [(32, 1), (24, 620),
        (20, 13888), (16, 36518), (12, 13888), (8, 620), (0, 1)]

end

@testset "miscknowncodes.jl" begin
    using Nemo, CodingTheory

    R, (x, y) = PolynomialRing(Nemo.ZZ, ["x", "y"])

    # Hamming codes
    # Tetra code is Hammingcode(3, 2)
    # random Hamming code
    F, _ = FiniteField(2, 1, "α")
    C = HammingCode(2, 7)
    col = rand(1:length(C))
    # columns are 1, 2, ... 2^r - 1 written as binary numerals
    @test paritycheckmatrix(C)[:, col] == matrix(F, length(C) -
        dimension(C), 1, reverse(digits(col, base=2, pad=7)))
    # should be [2^r - 1, 2^r - 1 - r, 3]
    @test length(C) == 2^7 - 1
    @test dimension(C) == 2^7 - 1 - 7
    C.d = missing
    @test minimumdistance(C) == 3
    C = HammingCode(2, 3)
    hamWE = weightenumerator(C, "Hamming")
    @test polynomial(hamWE) == x^7 + 7*x^3*y^4 + 7*x^4*y^3 + y^7
    n = length(C)
    C.weightenum = missing
    @test polynomial(hamWE) == divexact((x + y)^n + n*(x + y)^div(n - 1, 2)*(y - x)^div(n + 1, 2), n + 1)

    # simplex codes
    # random simplex code
    C = SimplexCode(2, 4)
    known = C.weightenum
    C.weightenum = missing
    HWEbf = weightenumerator(C, "Hamming")
    C.weightenum = missing
    HWEtrellis = weightenumerator(C, "Hamming", "trellis")
    @test CWEtoHWE(known) == HWEbf
    @test HWEbf == HWEtrellis
    # all nonzero codewords have weights q^{r - 1}
    flag = true
    for exps in [collect(exponent_vectors(polynomial(HWEtrellis)))[i][1]
            for i in 1:length(polynomial(HWEtrellis))]
        if !iszero(exps % 2^(4 - 1))
            flag = false
            break
        end
    end
    @test flag == true
    @test length(C) == 2^4 - 1
    @test dimension(C) == 4
    C = SimplexCode(2, 3)
    @test MacWilliamsIdentity(C, weightenumerator(C, "Hamming", "bruteforce")) == hamWE

    # Golay codes
    C = ExtendedGolayCode(2)
    @test isselfdual(C)
    C.weightenum = missing
    @test polynomial(weightenumerator(C, "Hamming")) == y^24 + 759*x^8*y^16 + 2576*x^12*y^12 + 759*x^16*y^8 + x^24
    C = GolayCode(2)
    C.weightenum = missing
    @test polynomial(weightenumerator(C, "Hamming")) == y^23 + 253*x^7*y^16 +
        506*x^8*y^15 + 1288*x^11*y^12 + 1288*x^12*y^11 + 506*x^15*y^8 + 253*x^16*y^7 + x^23
    C = ExtendedGolayCode(3)
    @test isselfdual(C)
    # well-known weight enumerators
    C.weightenum = missing
    @test polynomial(weightenumerator(C, "Hamming")) == y^12 + 264*x^6*y^6 + 440*x^9*y^3 + 24*x^12
    C = GolayCode(3)
    @test polynomial(weightenumerator(C, "Hamming")) == y^11 + 132*x^5*y^6 + 132*x^6*y^5 + 330*x^8*y^3 + 110*x^9*y^2 + 24*x^11
    # cyclic code with generator polynomial g(x) = -1 + x^2 - x^3 + x^4 + x^5
    # and idempotent e(x) = -(x^2 + x^6 + x^7 + x^8 + x^10)
    # should be eqivalent to the [11, 6, 5] Golay code (maybe permutation?)

    # tetra code
    C = TetraCode()
    C.weightenum = missing
    CWE = polynomial(weightenumerator(C, "complete"))
    vars = gens(parent(CWE))
    @test CWE == vars[1]^4 + vars[1]*vars[2]^3 + 3*vars[1]*vars[2]^2*vars[3] +
        3*vars[1]*vars[2]*vars[3]^2 + vars[1]*vars[3]^3
    # MacWilliams & Sloane have a different generator matrix for this
    # this test is false, as it makes sense to be, shoud be permutation equivalent though
    # test this later when implemented
    # F, _ = FiniteField(3, 1, "α")
    # G = matrix(F, [1 0 2 2 ; 0 1 2 1])
    # H = matrix(F, [1 1 1 0; 1 2 0 1])
    # C2 = LinearCode(F, 4, 2, 3, G, G, H, H, G, H, missing)
    # @test isequivalent(C, C2)

    # Hadamard code
    # the dual code of the Hamming code is the shortened Hadamard code
    # equivalent to RM(1, m)

end

@testset "cyclotomic.jl" begin
    using CodingTheory

    @test allcyclotomiccosets(2, 15, false) == [[0], [1, 2, 4, 8], [3, 6, 12, 9],
        [5, 10], [7, 14, 13, 11]]
    @test allcyclotomiccosets(3, 13, true) == [[0], [1, 3, 9], [2, 5, 6], [4, 10, 12] , [7, 8, 11]]

end

@testset "cycliccode.jl" begin
    using Nemo, CodingTheory

    # @test_throws ErrorException("There is no finite field of order 6.") CyclicCode(6, 9, [[0]])

    # this fails due to a column swap
    # CyclicCode(9, 14, definingset([1, 2, 3], 9, 14, false))


    # examples: Huffman & Pless
    cosets = definingset([1, 2, 3, 4, 5, 6], 2, 7, false)
    C = CyclicCode(2, 7, cosets)
    R = polynomialring(C)
    x = gen(R)
    @test dimension(C) == 1
    @test generatorpolynomial(C) == 1 + x + x^2 + x^3 + x^4 + x^5 + x^6
    @test idempotent(C) == 1 + x + x^2 + x^3 + x^4 + x^5 + x^6
    cosets = definingset([0, 1, 2, 4], 2, 7, false)
    C = CyclicCode(2, 7, cosets)
    @test dimension(C) == 3
    @test generatorpolynomial(C) == 1 + x^2 + x^3 + x^4
    @test idempotent(C) == 1 + x^3 + x^5 + x^6
    cosets = definingset([0, 3, 5, 6], 2, 7, false)
    C = CyclicCode(2, 7, cosets)
    @test dimension(C) == 3
    @test generatorpolynomial(C) == 1 + x + x^2 + x^4
    @test idempotent(C) == 1 + x + x^2 + x^4
    cosets = definingset([1, 2, 4], 2, 7, false)
    C = CyclicCode(2, 7, cosets)
    @test dimension(C) == 4
    @test generatorpolynomial(C) == 1 + x + x^3
    @test idempotent(C) == x + x^2 + x^4
    cosets = definingset([3, 5, 6], 2, 7, false)
    C = CyclicCode(2, 7, cosets)
    @test dimension(C) == 4
    @test generatorpolynomial(C) == 1 + x^2 + x^3
    @test idempotent(C) == x^3 + x^5 + x^6

    # the dual of a cyclic code is the complement code with multiplier -1 (p.146)
    # do Theorem 4.4.11 on page 147

    # a self-orthogonal binary cyclic code is doubly-even

    # example: Huffman & Pless
    C = BCHCode(3, 13, 2, 1)
    @test definingset(C) == [1, 3, 9]
    R = polynomialring(C)
    x = gen(R)
    @test generatorpolynomial(C) == 2 + x + x^2 + x^3
    @test dimension(C) == 10
    @test minimumdistance(C) == 3
    C = BCHCode(3, 13, 3, 1)
    @test definingset(C) == [1, 2, 3, 5, 6, 9]
    @test generatorpolynomial(C) == 1 + 2*x + x^2 + 2*x^3 + 2*x^4 + 2*x^5 + x^6
    @test dimension(C) == 7
    @test minimumdistance(C) == 4
    C = BCHCode(3, 13, 5, 1)
    @test definingset(C) == [1, 2, 3, 4, 5, 6, 9, 10, 12]
    @test generatorpolynomial(C) == 2 + 2*x^2 + 2*x^3 + x^5 + 2*x^7 + x^8 + x^9
    @test dimension(C) == 4
    @test minimumdistance(C) == 7
    @test dimension(C) >= length(C) - ord(length(C), 3)*(5 - 1)

    R, (x, y) = PolynomialRing(Nemo.ZZ, ["x", "y"])

    # example: MacWilliams & Sloane
    # any cyclic code over GF(2^m) of length 2^m + 1 is reversible

    # example: MacWilliams & Sloane
    C = BCHCode(2, 31, 5, 1)
    @test dimension(C) == 21
    @test minimumdistance(C) == 5
    @test polynomial(MacWilliamsIdentity(C, weightenumerator(C, "Hamming"))) == y^31 + 310*x^12*y^19 + 527*x^16*y^15 + 186*x^20*y^11

    # example: Huffman & Pless
    C = ReedSolomonCode(13, 5, 1)
    @test length(C) == 12
    @test dimension(C) == 8
    @test minimumdistance(C) == 5
    @test isMDS(C) == true
    @test definingset(C) == [1, 2, 3, 4]
    R = polynomialring(C)
    x = gen(R)
    @test generatorpolynomial(C) == 10 + 2*x + 7*x^2 + 9*x^3 + x^4
    D = dual(C)
    @test dimension(D) == 4
    @test minimumdistance(D) == 9
    @test isMDS(D) == true
    @test definingset(D) == [0, 1, 2, 3, 4, 5, 6, 7]
    @test generatorpolynomial(D) == 3 + 12*x + x^2 + 5*x^3 + 11*x^4 + 4*x^5 + 10*x^6 + 5*x^7 + x^8
    Cc = complement(C)
    @test length(Cc) == 12
    @test dimension(Cc) == 4
    @test minimumdistance(Cc) == 9
    @test definingset(Cc) == [0, 5, 6, 7, 8, 9, 10, 11]
    @test generatorpolynomial(Cc) == 9 + 6*x + 12*x^2 + 10*x^3 + 8*x^4 + 6*x^5 + 9*x^6 + 4*x^7 + x^8

    # example: Huffman & Pless
    C = ReedSolomonCode(16, 7, 1)
    @test length(C) == 15
    @test dimension(C) == 9
    @test minimumdistance(C) == 7
    @test definingset(C) == [1, 2, 3, 4, 5, 6]
    R = polynomialring(C)
    x = gen(R)
    α = primitiveroot(C)
    @test generatorpolynomial(C) == α^6 + α^9*x + α^6*x^2 + α^4*x^3 + α^14*x^4 + α^10*x^5 + x^6

    # example: MacWilliams & Sloane
    C = ReedSolomonCode(5, 3, 1)
    z = gen(polynomialring(C))
    @test generatorpolynomial(C) == z^2 + 4*z + 3

    # example: MacWilliams & Sloane
    C = ReedSolomonCode(8, 6)
    @test dimension(C) == 2
    z = gen(polynomialring(C))
    α = primitiveroot(C)
    @test idempotent(C) == α^4*z + α*z^2 + α^4*z^3 + α^2*z^4 + α^2*z^5 + α*z^6

    # example: MacWilliams & Sloane
    C = ReedSolomonCode(8, 3, 5)
    @test dimension(C) == 5
    z = gen(polynomialring(C))
    α = primitiveroot(C)
    @test generatorpolynomial(C) == α^4 + α*z + z^2
    # expand this code over F_2, is equivalent to the following BCH code
    # C2 = BCHCode(2, 21, 3, 1) # maybe not b = 1?
    # z2 = gen(polynomialring(C2))
    # @test generatorpolynomial(C2) == 1 + z2 + z2^2 + z2^4 + z2^6
    # @test isequivalent(expC, C2)

    # # example: MacWilliams & Sloane
    # # extended Reed-Solomon codes have distance d + 1
    # # TODO: fix extend here
    # # extC = extend(C)
    # # @test minimumdistance(extC) == minimumdistance(C) + 1

    # example: MacWilliams & Sloane
    # some [15, 6, 6] binary BCH code with min polys (-1, 0, 1) is reversible
    # C = BCHCode(2, 15, 6, 13)
    # println(C)
    # @test isreversible(C) == true

    # RS codes contain BCH codes
    C = ReedSolomonCode(16, 5)
    C2 = BCHCode(2, 15, 5)
    @test C2 ⊆ C
    @test C2 ⊂ C
    @test issubcode(C2, C)
    @test C == CyclicCode(16, 15, definingset([i for i = 0:(0 + 5 - 2)], 16, 15, false))
    @test C == BCHCode(16, 15, 5)
    @test designdistance(C) == 5
    @test isnarrowsense(C)
    @test isprimitive(C)

end
#
# @testset "GeneralizedReedSolomon.jl" begin
#     using CodingTheory
#
#     # the [q, k, q - k + 1] extended narrow-sense Reed-Solomon code over 𝔽_q is GRS and MDS
#
#     # narrrow-sense RS codes are GRS codes with n = q - 1, γ_i = α^i, and v_i = 1 for 0 <= i <= n - 1
# end

# @testset "quantumcode.jl" begin
#     using Nemo, CodingTheory
#
#     S = SteaneCode()
# end

@testset "LDPC.jl" begin
    using Nemo, CodingTheory

    # example from Ryan & Lin
    F, _ = FiniteField(2, 1, "α")
    H = matrix(F, [
        1 1 1 1 0 0 0 0 0 0;
        1 0 0 0 1 1 1 0 0 0;
        0 1 0 0 1 0 0 1 1 0;
        0 0 1 0 0 1 0 1 0 1;
        0 0 0 1 0 0 1 0 1 1])
    C = LDPCCode(H)
    @test originalparitycheckmatrix(C) == H
    @test columnrowbounds(C) == (2, 4)
    @test rate(C) == 3 / 5
    @test isregular(C)
    @test unique!(variabledegreedistribution(C)) == [2]
    @test unique(checkdegreedistribution(C)) == [4]
    R = parent(variabledegreepolynomial(C))
    x = gen(R)
    @test variabledegreepolynomial(C) == x
    @test checkdegreepolynomial(C) == x^3

end
