@testitem "Quantum/subsystem_code.jl" begin
    using Oscar, CodingTheory

    @testset "Stabilizer Subsystem code" begin
        # Poulin, "Stabilizer Formalism for Operator Quantum Error Correction", (2008)
        # [[9, 1, 4, 3]] gauged Shor code
        S = ["XXXXXXIII", "XXXIIIXXX", "ZZIZZIZZI", "IZZIZZIZZ"]
        # these are the {X, Z} pairings
        G_ops = [
            "IZZIIIIII",
            "IIXIIIIIX",
            "IIIIZZIII",
            "IIIIIXIIX",
            "ZZIIIIIII",
            "XIIIIIXII",
            "IIIZZIIII",
            "IIIXIIXII",
        ]
        G = S ∪ G_ops
        L = ["ZZZZZZZZZ", "XXXXXXXXX"]
        Q = SubsystemCode(G)
        @test length(Q) == 9
        @test dimension(Q) == 1
        @test Q.r == 4
        @test_broken minimum_distance(Q) == 3
        @test LogicalTrait(typeof(Q)) == HasLogicals()
        @test GaugeTrait(typeof(Q)) == HasGauges()

        Q2 = SubsystemCode(S, L, G_ops)
        @test are_equivalent(Q, Q2)

        # # TODO: BaconShorCode
        Q3 = BaconShorCode(3, 3)

        F = Oscar.Nemo.Native.GF(2)
        A = matrix(F, 3, 3, ones(Int, 3, 3))
        Q4 = BravyiBaconShorCode(A)
        @test dimension(Q4) == rank(A)
        # # TODO: add a test here on the min distances
        # min d should be min(row wts, col wts)
        @test are_equivalent(Q3, Q4)
    end

    # # Klappenecker and Sarvepalli (2007) give a CSS construction equivalent to Bacon-Shor

end
