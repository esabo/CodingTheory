@testitem "Shor-Laflamme Hamming enumerators" begin
    using Oscar, CodingTheory

    F = GF(2)
    bell = StabilizerCode(matrix(F, [
        1 1  0 0
        0 0  1 1
    ]); logs_alg=:sys_eqs)
    bell_SL = Shor_Laflamme_weight_enumerator(bell)
    @test bell_SL isa ShorLaflammeWeightEnumerator
    @test bell_SL.A.counts == Dict(0 => BigInt(1), 2 => BigInt(3))
    @test bell_SL.B.counts == bell_SL.A.counts
    @test weight_enumerator(bell; set=:A) === bell_SL.A
    @test weight_distribution(bell; set=:quotient) == Dict{Int, BigInt}()
    @test sum(values(bell_SL.A.counts)) == cardinality(bell)
    @test !haskey(bell.cache, :stabilizer_operators)

    five_qubit = FiveQubitCode()
    enumerator = SL_weight_enumerator(five_qubit)
    @test enumerator.A.counts ==
        Dict(0 => BigInt(1), 4 => BigInt(15))
    @test sum(values(enumerator.B.counts)) == BigInt(64)
    @test minimum(keys(weight_distribution(
        five_qubit; set=:quotient))) == 3
    @test five_qubit.cache[:d] == 3
    @test weight_distribution_array(
        five_qubit; set=:stabilizers)[5] == 15

    F4 = GF(4)
    additive = StabilizerCode(
        matrix(F4, [1 0]); logs_alg=:sys_eqs)
    additive_SL = Shor_Laflamme_weight_enumerator(additive)
    @test additive.k == 1 // 2
    @test additive_SL.A.counts ==
        Dict(0 => BigInt(1), 1 => BigInt(1))
    @test additive_SL.B.counts ==
        Dict(0 => BigInt(1), 1 => BigInt(7))

    too_large = random_stabilizer_code(F, 8, 0)
    @test_throws ArgumentError Shor_Laflamme_weight_enumerator(
        too_large; max_terms=100)
end
