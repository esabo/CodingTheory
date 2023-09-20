@testset "stabilizercode.jl" begin
    using Oscar, CodingTheory

    Q = Q1573()
    logs = logicals(Q)
    new_stab = logs[1][2] + logs[2][2]
    Q2 = augment(Q, new_stab, false)
    Q3 = expurgate(Q2, [9], false)
    @test are_equivalent(Q, Q3)
end
