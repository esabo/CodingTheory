@testset "chaincomplex.jl" begin
    # TODO: better chain complex tests
    # these tests are simply to check that it runs at all
    @test !isnothing(tensor_product(ChainComplex(Q15RM()), ChainComplex(Q1573())))
    @test !isnothing(distance_balancing(Q15RM(), RepetitionCode(2,4)))
end
