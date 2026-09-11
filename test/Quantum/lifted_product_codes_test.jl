@testitem "Quantum/lifted_product_codes.jl" begin
    using CodingTheory, Oscar

    F = Oscar.Nemo.Native.GF(2)
    P, x = polynomial_ring(F, :x)
    R, _ = residue_ring(P, x^2 - 1)
    A = matrix(R, 1, 1, [R(1 + x)])
    Q = LiftedProductCode(A, R(1 + x))
    @test Q isa AbstractStabilizerCodeCSS
    @test ncols(X_stabilizers(Q)) == length(Q)
    @test iszero(X_stabilizers(Q) * transpose(Z_stabilizers(Q)))
    @test GeneralizedHypergraphProductCode(A, R(1 + x)) isa LiftedProductCode

    biased = BiasTailoredLiftedProductCode(A, A)
    @test biased isa AbstractStabilizerCode
    @test ncols(stabilizers(biased)) == 2 * length(biased)
end
