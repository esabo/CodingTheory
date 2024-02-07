@testset "LDPC/codes.jl" begin
    using Oscar, CodingTheory

    # example from Ryan & Lin
    F = GF(2)
    H = matrix(F, [
        1 1 1 1 0 0 0 0 0 0;
        1 0 0 0 1 1 1 0 0 0;
        0 1 0 0 1 0 0 1 1 0;
        0 0 1 0 0 1 0 1 0 1;
        0 0 0 1 0 0 1 0 1 1])
    C = LDPCCode(H)
    @test column_row_bounds(C) == (2, 4)
    # @test rate(C) == 3 / 5
    @test is_regular(C)
    @test unique!(variable_degree_distribution(C)) == [2]
    @test unique(check_degree_distribution(C)) == [4]
    R = parent(variable_degree_polynomial(C))
    x = gen(R)
    @test variable_degree_polynomial(C) == x
    @test check_degree_polynomial(C) == x^3
end
