@testitem "Quantum/fold_product_codes.jl" begin
    using CodingTheory, Oscar

    Q = SPCDFoldProductCode(2)
    @test Q isa AbstractStabilizerCodeCSS
    @test length(Q) == 16
    @test dimension(Q) == 2
    @test Q.d == 4
    @test iszero(X_stabilizers(Q) * transpose(Z_stabilizers(Q)))

    A = Q422()
    product = asymmetric_product(A, A)
    @test length(product) == 16
    @test ncols(X_stabilizers(product)) == 16
end
