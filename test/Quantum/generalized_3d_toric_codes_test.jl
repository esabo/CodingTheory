@testitem "Quantum/generalized_3d_toric_codes.jl" begin
    using CodingTheory, Oscar

    F = Oscar.Nemo.Native.GF(2)
    LR, (x, y, z) = laurent_polynomial_ring(F, [:x, :y, :z])
    infinite = Generalized3DToricCode(x + y, 1 + z)
    @test !(infinite isa AbstractSubsystemCode)

    finite = FiniteGeneralized3DToricCode(x + y, 1 + z, 2, 2, 2)
    @test finite isa AbstractStabilizerCodeCSS
    @test !finite.twisted
    @test length(finite) == 16
    @test ncols(X_stabilizers(finite)) == 16

    twisted = BBCode3D(x + y, 1 + z, (2, 1), (0, 2), 2)
    @test twisted.twisted
    @test length(twisted) == 16
end
