@testitem "Quantum/concatenated_codes.jl" begin
    using CodingTheory, Oscar

    Q = concatenate(Q422(), FiveQubitCode())
    @test Q isa AbstractStabilizerCode
    @test (length(Q), dimension(Q)) == (20, 2)
    @test ncols(stabilizers(Q)) == 40
    @test field(Q) == field(Q422())
end
