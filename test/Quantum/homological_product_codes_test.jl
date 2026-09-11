@testitem "Quantum/homological_product_codes.jl" begin
    using CodingTheory, Oscar

    Q = homological_product(Q422(), Q422())
    @test Q isa AbstractStabilizerCodeCSS
    @test length(Q) == 16
    @test ncols(X_stabilizers(Q)) == 16
    @test iszero(X_stabilizers(Q) * transpose(Z_stabilizers(Q)))

    random = random_homological_product_code(3, 1, 3, 1)
    @test random isa AbstractStabilizerCodeCSS
    @test length(random) == 9
end
