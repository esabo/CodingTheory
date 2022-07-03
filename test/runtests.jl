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
    @test length(C) == 7
    @test rank(G) == LinearCodeMod.dimension(C)
    @test LinearCodeMod.dimension(C) == 4
    @test ismissing(C.d)
    LinearCodeMod.setminimumdistance!(C, 3)
    # cannot find below, don't want to provoke it just yet
    @test minimumdistance(C) == 3
    # @test !LinearCodeMod.isMDS(C)
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
    @test iszero(LinearCodeMod.syndrome(v, C))
    v = [1; 0; 0; 0; 0; 1; 1];
    @test iszero(LinearCodeMod.syndrome(v, C))

    # lower rank test
    GandG = vcat(G, G);
    CGandG = LinearCode(G);
    @test rank(GandG) == LinearCodeMod.dimension(CGandG)
    @test G == generatormatrix(CGandG)
    # this fails, is of wrong size, probably a transpose mistake
    # @test GandG == originalgeneratormatrix(CGandG)
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
    # RM(0, m) is the length 2^m repetition code
    # puncturing RM(1, m) and taking the even subcode is the simplex code S_m of
    # length 2^m - 1


end

# tests for simplex codes
# all nonzero codewords have weights q^{r - 1}
# length (q^r - 1)/(q - 1)
# dim r

# tests for Golay codes
# all weights are divisible by four for nonpunctured
# Hamming weight enumerator of [12, 6, 6] is
# [[1, 0, 12], [264, 6, 6], [440, 9, 3], [24, 12, 0]]
# Hamming weight enumerator of [11, 6, 5] is
# [[1, 0, 11], [132, 5, 6], [132, 6, 5], [330, 8, 3], [110, 9, 2], [24, 11, 0]]
# cyclic code with generator polynomial g(x) = -1 + x^2 - x^3 + x^4 + x^5
# and idempotent e(x) = -(x^2 + x^6 + x^7 + x^8 + x^10)
# should be eqivalent to the [11, 6, 5] Golay code (maybe permutation?)
