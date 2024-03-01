@testset "Quantum/stabilizer_code.jl" begin
    using Oscar, CodingTheory

    Q = Q1573()
    logs = logicals(Q)
    new_stab = logs[1][2] + logs[2][2]
    Q2 = augment(Q, new_stab, verbose = false)
    Q3 = expurgate(Q2, [9], verbose = false)
    @test are_equivalent(Q, Q3)

    # basic tests
    @test is_CSS_T_code(SteaneCode()) == false
    @test is_CSS_T_code(Q1513())
    @test is_CSS_T_code(Q1573()) == false
    @test is_CSS_T_code(SmallestInterestingColorCode())

    # examples in https://arxiv.org/abs/1910.09333
    # Example 2 - [[6, 2, 2]], transversal T implements logical identity
    G = matrix(GF(2), 4, 12, [
        1 1 1 1 1 1 0 0 0 0 0 0
        0 0 0 0 0 0 1 1 0 0 0 0
        0 0 0 0 0 0 0 0 1 1 0 0
        0 0 0 0 0 0 0 0 0 0 1 1])
    S = StabilizerCode(G)
    @test is_CSS_T_code(S)

    # Example 4/5 - [[16, 3, 2]], transversal T implements logical CCZ (up to logical Paulis)
    X = matrix(GF(2), 3, 16,[
        1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0;
        0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0;
        0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1])
    Z = matrix(GF(2), 10, 16, [
        1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0;
        0 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0;
        0 0 1 1 0 0 1 1 0 0 0 0 0 0 0 0;
        0 0 0 0 1 1 0 0 1 1 0 0 0 0 0 0;
        0 0 0 0 0 1 1 0 0 1 1 0 0 0 0 0;
        0 0 0 0 0 0 1 1 0 0 1 1 0 0 0 0;
        0 0 0 0 0 0 0 0 1 1 0 0 1 1 0 0;
        0 0 0 0 0 0 0 0 0 1 1 0 0 1 1 0;
        0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1;
        0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1])
    S = CSSCode(X, Z)
    @test is_CSS_T_code(S)

    # Example 6 - [[64, 15, 4]], transversal T realizes a logical diagonal gate at the 3rd level that is a product of elementary gates
    C1 = ReedMullerCode(2, 6)
    C2 = ReedMullerCode(1, 6)
    S = CSSCode(C1, C2)
    @test is_CSS_T_code(S)

    # Example 7 - [[128, 21, 4]], transversal T realizes a logical diagonal gate at the 3rd level that is a product of elementary gates
    C1 = ReedMullerCode(2, 7)
    C2 = ReedMullerCode(1, 7)
    S = CSSCode(C1, C2)
    @test is_CSS_T_code(S)
    @test is_triorthogonal(C1.G)

end
