@testitem "Quantum/generalized_shor_codes.jl" begin
    using CodingTheory, Oscar

    C = RepetitionCode(2, 2)
    Q = GeneralizedShorCode(C, C)
    @test Q isa AbstractSubsystemCode
    @test GaugeTrait(typeof(Q)) == HasGauges()
    @test (length(Q), dimension(Q)) == (4, 1)
    @test Q.r >= 0
    @test ncols(stabilizers(Q)) == 8
    @test ncols(gauges_matrix(Q)) == 8
    @test BaconCasaccinoConstruction(C, C) isa GeneralizedShorCode
end
