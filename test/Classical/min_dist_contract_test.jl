@testitem "Classical/min_dist_contract.jl" begin
    using Oscar, CodingTheory, SparseArrays

    function reset_bounds!(C::AbstractLinearCode)
        C.l_bound = 1
        C.u_bound = C.n + 1
        C.d = missing
        delete!(getfield(C, :cache), :minimum_distance_witness)
    end

    F = Oscar.Nemo.Native.GF(2)
    G = matrix(F, [
        1 0 0 0 0 1 1;
        0 1 0 0 1 0 1;
        0 0 1 0 1 1 0;
        0 0 0 1 1 1 1
    ])
    C = LinearCode(G)
    reset_bounds!(C)

    @testset "Code-level return contract" begin
        d, witness = minimum_distance(C, alg=:BZ, info_set_alg=:Brouwer)
        @test d isa Int
        @test d == 3
        @test size(witness) == (1, C.n)
        @test base_ring(witness) == C.F
        @test !iszero(witness)
        @test wt(witness) == 3
        @test iszero(parity_check_matrix(C) * transpose(witness))

        d2, witness2 = minimum_distance(C)
        @test d2 == 3
        @test witness2 == witness
    end

    @testset "Wagner miss still returns a matrix" begin
        reset_bounds!(C)
        d, witness = CodingTheory._minimum_distance_wagner_mitm_binary(C; max_d=2)
        @test d == -1
        @test size(witness) == (1, C.n)
        @test iszero(witness)
    end

    @testset "Matrix-level kernel on sparse H" begin
        H = sparse([1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3],
                   [2, 3, 4, 5, 1, 3, 4, 6, 1, 2, 4, 7],
                   ones(Int, 12), 3, 7)
        d, witness = CodingTheory._minimum_distance(H; alg=:Wagner, max_d=5)
        @test d == 3
        @test witness isa Vector{Int}
        @test length(witness) == 7
        @test sum(witness) == 3
        @test iszero((Array(H) * witness) .% 2)
    end

    @testset "Queue scheduler preserves its winning witness" begin
        C_queue = LinearCode(G)
        reset_bounds!(C_queue)
        d, witness = CodingTheory._minimum_distance_BZ_binary(
            C_queue; info_set_alg=:Brouwer, scheduler=:queue)
        @test d == 3
        @test wt(witness) == d
        @test iszero(parity_check_matrix(C_queue) * transpose(witness))
    end
end
