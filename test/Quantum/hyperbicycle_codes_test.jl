@testitem "Quantum/hyperbicycle_codes.jl" begin
    using CodingTheory, Oscar

    F = Oscar.Nemo.Native.GF(2)
    h = matrix(F, 1, 2, [1, 1])
    css = HyperBicycleCodeCSS([h], [h], 1)
    @test css isa AbstractStabilizerCodeCSS
    @test length(css) == 4
    @test iszero(X_stabilizers(css) * transpose(Z_stabilizers(css)))

    one = matrix(F, 1, 1, [1])
    noncss = HyperBicycleCode([one], [one], 1)
    @test noncss isa AbstractStabilizerCode
    @test length(noncss) == 1
    @test ncols(stabilizers(noncss)) == 2

    @test_throws ArgumentError HyperBicycleCodeCSS(typeof(h)[], typeof(h)[], 1)
end
