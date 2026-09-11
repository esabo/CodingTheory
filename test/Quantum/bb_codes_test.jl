@testitem "Quantum/BB_codes.jl" begin
    using CodingTheory, Oscar

    F = Oscar.Nemo.Native.GF(2)
    P, (x, y) = polynomial_ring(F, [:x, :y])
    R, _ = quo(P, ideal(P, [x^2 - 1, y^2 - 1]))
    standard = BBCode(R(x + y), R(1 + x * y))
    @test standard isa AbstractStabilizerCodeCSS
    @test !standard.twisted
    @test length(standard) == 8
    @test ncols(X_stabilizers(standard)) == 8

    LR, (u, v) = laurent_polynomial_ring(F, [:u, :v])
    infinite = InfiniteBBCode(u + v, 1 + u)
    @test !(infinite isa AbstractSubsystemCode)
    twisted = BBCode(u + v, 1 + u, (2, 1), (0, 2))
    @test twisted.twisted
    @test (length(twisted), dimension(twisted)) == (8, 2)

    U, z = polynomial_ring(F, :z)
    coprime = BBCode(1 + z, 1 + z^2, 5)
    @test (length(coprime), dimension(coprime)) == (10, 2)

    set_minimum_distance!(coprime, 2)
    @test coprime.d == 2
    @test coprime.l_bound == 2
    @test coprime.u_bound == 2
end
