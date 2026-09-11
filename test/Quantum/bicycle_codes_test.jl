@testitem "Quantum/bicycle_codes.jl" begin
    using CodingTheory, Oscar

    F = Oscar.Nemo.Native.GF(2)
    P, x = polynomial_ring(F, :x)
    R, _ = residue_ring(P, x^7 - 1)
    Q = GeneralizedBicycleCode(R(1 + x + x^3), R(1 + x^2))
    @test Q isa AbstractStabilizerCodeCSS
    @test length(Q) == 14
    @test ncols(X_stabilizers(Q)) == 14
    @test iszero(X_stabilizers(Q) * transpose(Z_stabilizers(Q)))

    B = BicycleCode(residue_polynomial_to_circulant_matrix(R(1 + x)))
    @test B isa GeneralizedBicycleCode
    @test isempty(character_vector(B))
end
