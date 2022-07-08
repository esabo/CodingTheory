using CodingTheory
using Test

# these are subject to change as they develop, let me know when it breaks
@testset "Types" begin
    @test AbstractLinearCode <: AbstractCode
    @test AbstractCyclicCode <: AbstractLinearCode
    @test AbstractBCHCode <: AbstractCyclicCode
    @test AbstractReedSolomonCode <: AbstractBCHCode
    @test AbstractAdditiveCode <: AbstractCode
    @test AbstractStabilizerCode <: AbstractAdditiveCode
    @test AbstractCSSCode <: AbstractStabilizerCode
end

# @testset "utils.jl" begin




# end

@testset "linearcode.jl" begin
    using Nemo
    using CodingTheory

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
    # cannot find below, don't want to provoke it just yet
    @test minimumdistance(C) == 3
    # @test !isMDS(C)
    # @test numbercorrectableerrors(C) == 1
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
end

@testset "ReedMuller.jl" begin
    using Nemo
    using CodingTheory

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
end

@testset "miscknowncodes.jl" begin
    using Nemo
    using CodingTheory

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
    @test minimumdistance(C) == 3

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
    for term in HWEtrellis.polynomial
        if !iszero(term[2] % 2^(4 - 1))
            flag = false
            break
        end
    end
    @test flag == true
    @test length(C) == 2^4 - 1
    @test dimension(C) == 4

    # Golay codes
    C = ExtendedGolayCode(3)
    # well-known weight enumerators
    @test weightenumerator(C, "Hamming").polynomial == [[1, 0, 12], [264, 6, 6], [440, 9, 3], [24, 12, 0]]
    C = GolayCode(3)
    @test weightenumerator(C, "Hamming").polynomial == [[1, 0, 11], [132, 5, 6], [132, 6, 5], [330, 8, 3], [110, 9, 2], [24, 11, 0]]
    # cyclic code with generator polynomial g(x) = -1 + x^2 - x^3 + x^4 + x^5
    # and idempotent e(x) = -(x^2 + x^6 + x^7 + x^8 + x^10)
    # should be eqivalent to the [11, 6, 5] Golay code (maybe permutation?)

end

@testset "cyclotomic.jl" begin
    using CodingTheory

    @test allcyclotomiccosets(2, 15, false) == [[0], [1, 2, 4, 8], [3, 6, 12, 9],
        [5, 10], [7, 14, 13, 11]]

end
