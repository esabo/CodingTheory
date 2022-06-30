
using my package here
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

@testset "utils.jl" begin




end


@testset "linearcode.jl" begin
    # Hammming code
    F, _ = FiniteField(2, 1, "α");
    G = matrix(F, [1 0 0 0 0 1 1;
           0 1 0 0 1 0 1;
           0 0 1 0 1 1 0;
           0 0 0 1 1 1 1]);
    C = LinearCode(G);
    @test length(C) == 7
    @test rank(G) == dimension(C)
    @test dimension(C) == 4
    @test ismissing(minimumdistance(C))
    setminimumdistance!(C, 3)
    @test minimumdistance(C) == 3
    @test !isMDS(C)
    @test G == generatormatrix(C)
    @test G == originalgeneratormatrix(C)
    H = paritycheckmatrix(C)
    @test iszero(G * H')
    @test iszero(H * G')
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
    @test iszero(syndrome(v, C))
    v = [1; 0; 0; 0; 0; 1; 1];
    @test iszero(syndrome(v, C))

    # lower rank test
    GandG = vcat(G, G);
    CGandG = LinearCode(G);
    @test rank(GandG) == dimension(CGandG)
    @test G == generatormatrix(CGandG)
    @test GandG == originalgeneratormatrix(CGandG)




end
