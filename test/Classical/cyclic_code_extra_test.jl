@testitem "Cyclic code extras" begin
    using Oscar, CodingTheory, Random

    @testset "Cyclotomic data and cyclic invariants" begin
        C = BCHCode(2, 15, 5, 1)
        Z = CodingTheory.defining_set(C)
        @test BCH_offset(C) == 1
        @test CodingTheory.qcosets_reps(C) == sort(first.(CodingTheory.qcosets(C)))
        @test CodingTheory.dual_defining_set(Z, 15) ==
              sort([mod(-i, 15) for i in setdiff(0:14, Z)])

        R = polynomial_ring(C)
        x = gen(R)
        @test CodingTheory.generator_polynomial(C) *
              CodingTheory.parity_check_polynomial(C) == x^length(C) - 1
        @test HT_bound(C) >= BCH_bound(C)
        @test Roos_bound(C) >= HT_bound(C)

        @test CodingTheory.splitting_field(C) == parent(CodingTheory.primitive_root(C))
        @test Int(order(CodingTheory.splitting_field(C))) == 16
        expected_trace_reps = Int[]
        unseen = Set(setdiff(0:14, Z))
        while !isempty(unseen)
            r = minimum(unseen)
            Q = CodingTheory.cyclotomic_coset(r, 2, 15)
            push!(expected_trace_reps, first(Q))
            setdiff!(unseen, Q)
        end
        @test CodingTheory.trace_representation(C) == sort(expected_trace_reps)

        @test is_reversible(CyclicCode(2, 7, [1, 2, 3, 4, 5, 6]))
        @test !is_reversible(CyclicCode(2, 7, [1, 2, 4]))
        @test is_antiprimitive(BCHCode(2, 3, 2, 1))

        table, result = mktemp() do _, io
            value = redirect_stdout(io) do
                CodingTheory.qcoset_table(7, 7, 2)
            end
            flush(io)
            seekstart(io)
            return read(io, String), value
        end
        @test isnothing(result)
        @test occursin("n = 7", table)
        @test occursin("C_1 = {1, 2, 4}", table)
    end

    @testset "Mattson-Solomon transform" begin
        C = BCHCode(2, 7, 3, 1)
        E = CodingTheory.splitting_field(C)
        β = CodingTheory.primitive_root(C)
        row = [E(generator_matrix(C)[1, j]) for j in 1:length(C)]
        MS = MattsonSolomon_transform(row, β)
        @test inverse_MattsonSolomon_transform(MS, length(C), β) == row

        vanish = sort([mod(j, length(C)) for j in 1:length(C)
                       if iszero(coeff(MS, length(C) - j))])
        @test vanish == CodingTheory.defining_set(C)

        rng = MersenneTwister(2026)
        v = [E(rand(rng, 0:1)) for _ in 1:length(C)]
        @test inverse_MattsonSolomon_transform(
                  MattsonSolomon_transform(v, β), length(C), β) == v
    end

    @testset "Multipliers and constituents" begin
        C = CyclicCode(2, 7, [1, 2, 4])
        M = multiplier_group(C)
        @test M == [1, 2, 4]
        for a in M
            @test apply_multiplier(C, a) == C
        end
        @test apply_multiplier(C, 3) != C

        v = collect(1:7)
        @test apply_multiplier(v, 3) == [1, 6, 4, 2, 7, 5, 3]
        @test_throws ArgumentError apply_multiplier(v, 7)

        D = apply_multiplier(C, 3)
        equiv, a = is_multiplier_equivalent(C, D)
        @test equiv
        @test apply_multiplier(C, a) == D

        Hs, _ = multiplier_subgroup_Sn(C)
        @test Int(order(Hs)) == length(M)

        ambient = CodingTheory.ambient_constituents(2, 7)
        parts = CodingTheory.constituents(C)
        @test length(ambient) == 3
        @test sum(dimension, ambient) == 7
        @test all(is_irreducible, ambient)
        @test sum(dimension, parts) == dimension(C)
        @test all(P -> is_subcode(P, C), parts)
    end

    @testset "Fire, polyadic, and entrywise products" begin
        F = Oscar.Nemo.Native.GF(2)
        R, x = polynomial_ring(F, :x)
        p = x^3 + x + 1
        Cfire = FireCode(p, 2)
        @test length(Cfire) == 21
        @test dimension(Cfire) == 15
        Rfire = polynomial_ring(Cfire)
        y = gen(Rfire)
        @test CodingTheory.generator_polynomial(Cfire) == (y^3 + 1) * (y^3 + y + 1)
        @test first(is_cyclic(Cfire))

        P = PolyadicCodes(2, 31, 3)
        T = TriadicCodes(2, 31)
        @test length(P.codes) == 3
        @test dimension.(P.codes) == [21, 21, 21]
        @test CodingTheory.defining_set.(P.codes) ==
              CodingTheory.defining_set.(T.codes)
        sets = CodingTheory.defining_set.(P.codes)
        @test sort(reduce(vcat, sets)) == collect(1:30)
        @test all(i -> sort(mod.(P.multiplier .* sets[i], 31)) ==
                       sets[mod1(i + 1, 3)], 1:3)
        @test all(first(is_cyclic(Ci)) for Ci in P.codes)

        Q = TetradicCodes(3, 13)
        qsets = CodingTheory.defining_set.(Q.codes)
        @test length(Q.codes) == 4
        @test dimension.(Q.codes) == [10, 10, 10, 10]
        @test sort(reduce(vcat, qsets)) == collect(1:12)
        @test all(i -> sort(mod.(Q.multiplier .* qsets[i], 13)) ==
                       qsets[mod1(i + 1, 4)], 1:4)

        C1 = CyclicCode(2, 7, [1, 2, 4])
        C2 = CyclicCode(2, 7, [3, 5, 6])
        S = Schur_product_code(C1, C2)
        @test S == Hadamard_product_code(C1, C2)
        @test S == componentwise_product_code(C1, C2)
        @test S == CodingTheory.entrywise_product_code(C1, C2)
        N1 = setdiff(0:6, CodingTheory.defining_set(C1))
        N2 = setdiff(0:6, CodingTheory.defining_set(C2))
        expected_nonzeros = sort(unique(mod(a + b, 7) for a in N1 for b in N2))
        @test CodingTheory.defining_set(S) == setdiff(0:6, expected_nonzeros)
        @test is_subcode(C1, Schur_product_code(C1))
    end

    @testset "Multiplier Subgroup" begin
        C = CyclicCode(2, 7, [[1, 2, 4]])
        H, inc = multiplier_subgroup_Zn(C)

        # the subgroup sits inside the unit group of Z/nZ
        R, _ = residue_ring(Oscar.ZZ, C.n)
        U, f = unit_group(R)
        @test Oscar.domain(inc) === H
        @test Oscar.order(Oscar.codomain(inc)) == Oscar.order(U)

        # every multiplier of the code appears in the subgroup
        multipliers = CodingTheory.multiplier_group(C)
        recovered = Set(Oscar.lift(f(inc(h))) for h in H)
        @test Set(Oscar.ZZ(a) for a in multipliers) ⊆ recovered
    end
end
