@testitem "Quantum/min_dist_exact.jl" begin
    using CodingTheory, Oscar, Random, SparseArrays

    F = Oscar.Nemo.Native.GF(2)

    function steane_code()
        H = matrix(F, [
            1 0 1 0 1 0 1;
            0 1 1 0 0 1 1;
            0 0 0 1 1 1 1
        ])
        return CSSCode(H, H)
    end

    function asymmetric_code()
        H_X = matrix(F, 1, 3, [1, 1, 0])
        H_Z = matrix(F, 1, 3, [0, 0, 1])
        return CSSCode(H_X, H_Z)
    end

    function check_witness(S, d, witness)
        @test size(witness) == (1, 2 * S.n)
        @test wt(witness) == d
        @test is_logical(S, witness)
    end

    @testset "Lazy logical detector" begin
        S = steane_code()
        @test !haskey(S.cache, :logs_mat)
        L = logicals_matrix(S)
        @test size(L) == (2 * S.k, 2 * S.n)
        @test haskey(S.cache, :logs_mat)
    end

    @testset "Gray and Wagner agree" begin
        for alg in (:Gray, :Wagner)
            S = steane_code()
            dx, wx = minimum_distance(S; which=:X, alg=alg)
            dz, wz = minimum_distance(S; which=:Z, alg=alg)
            @test (dx, dz) == (3, 3)
            check_witness(S, dx, wx)
            check_witness(S, dz, wz)

            # Cached calls preserve both the return shape and the witness.
            @test minimum_distance(S; which=:X, alg=alg) == (dx, wx)
            @test minimum_distance(S; which=:Z, alg=alg) == (dz, wz)
        end
    end

    @testset "Asymmetric X and Z distances" begin
        for alg in (:Gray, :Wagner)
            S = asymmetric_code()
            dx, wx = minimum_distance(S; which=:X, alg=alg)
            dz, wz = minimum_distance(S; which=:Z, alg=alg)
            @test dx == 1
            @test dz == 2
            check_witness(S, dx, wx)
            check_witness(S, dz, wz)

            S_full = asymmetric_code()
            d, witness = minimum_distance(S_full; which=:full, alg=alg)
            @test d == 1
            check_witness(S_full, d, witness)
            @test S_full.cache[:d] == 1
        end
    end

    @testset "Sparse quotient-aware matrix kernel" begin
        H = sparse([1], [3], [1], 1, 3)
        logical_checks = sparse([1], [1], [1], 1, 3)
        d, witness = CodingTheory._minimum_distance(
            H, logical_checks; alg=:Wagner, max_d=3)
        @test d == 1
        @test witness == [1, 0, 0]

        d_miss, witness_miss = CodingTheory._minimum_distance(
            H, logical_checks; alg=:Wagner, max_d=0)
        @test d_miss == -1
        @test iszero(witness_miss)

        H_X = sparse([1, 1], [1, 2], [1, 1], 1, 3)
        H_Z = sparse([1], [3], [1], 1, 3)
        S_sparse = CSSCode(H_X, H_Z)
        @test X_stabilizers(S_sparse) isa SparseMatrixCSC
        @test Z_stabilizers(S_sparse) isa SparseMatrixCSC
        @test minimum_distance(S_sparse; which=:X, alg=:Wagner)[1] == 1
        @test minimum_distance(S_sparse; which=:Z, alg=:Wagner)[1] == 2
    end

    @testset "Random small CSS cross-checks" begin
        Random.seed!(20260909)
        for _ in 1:5
            H_X = matrix(F, rand(0:1, 2, 7))
            while rank(H_X) != 2
                H_X = matrix(F, rand(0:1, 2, 7))
            end
            kernel_X = generator_matrix(LinearCode(H_X, true))
            H_Z = kernel_X[1:2, :]

            for which in (:X, :Z)
                S_gray = CSSCode(H_X, H_Z)
                S_wagner = CSSCode(H_X, H_Z)
                gray = minimum_distance(S_gray; which=which, alg=:Gray)
                wagner = minimum_distance(S_wagner; which=which, alg=:Wagner)
                @test gray[1] == wagner[1]
                check_witness(S_gray, gray...)
                check_witness(S_wagner, wagner...)
            end
        end
    end

    @testset "Input validation" begin
        S = asymmetric_code()
        @test_throws ArgumentError minimum_distance(S; which=:Y)
        @test_throws ArgumentError minimum_distance(S; which=:X, alg=:trellis)
    end
end
