@testitem "Quantum/weight_reduction.jl" begin
    using CodingTheory, Oscar, Random, SparseArrays

    F = Oscar.Nemo.Native.GF(2)
    H_X = matrix(F, 2, 4, [1, 1, 0, 0,
                           0, 1, 1, 0])
    H_Z = matrix(F, 1, 4, [0, 0, 0, 1])

    for method in (:Hastings, :reduced, :target)
        X, Z = copying(H_X, H_Z; method=method)
        @test ncols(X) == ncols(Z)
        @test iszero(X * transpose(Z))
        X_cone, Z_cone = copying_as_coning(
            H_X, H_Z; method=method, rng=Xoshiro(1))
        @test (X_cone, Z_cone) == (X, Z)
    end

    S = CSSCode(H_X, H_Z)
    copied = copying(S; method=:reduced)
    @test copied isa AbstractStabilizerCodeCSS
    @test iszero(X_stabilizers(copied) * transpose(Z_stabilizers(copied)))
    X_sparse, Z_sparse = copying(
        sparse(Array(H_X)), sparse(Array(H_Z)); method=:reduced)
    @test iszero(X_sparse * transpose(Z_sparse))

    X, Z = gauging(H_X, H_Z)
    @test iszero(X * transpose(Z))

    high_weight_X = matrix(F, 1, 4, [1, 1, 1, 1])
    X, Z = gauging_as_coning(
        high_weight_X, zero_matrix(F, 0, 4); rng=Xoshiro(1))
    @test iszero(X * transpose(Z))

    X, Z = thickening_and_choose_heights(H_X, H_Z, 2, [1])
    @test iszero(X * transpose(Z))

    @test_throws ArgumentError copying(H_X, matrix(F, 1, 3, [1, 0, 0]))
end
