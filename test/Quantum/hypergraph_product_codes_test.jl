@testitem "Quantum/hypergraph_product_codes.jl" begin
    using CodingTheory, Oscar

    C = RepetitionCode(2, 2)
    Q = HypergraphProductCode(C)
    @test Q isa AbstractHypergraphProductCode
    @test (length(Q), dimension(Q)) == (5, 1)
    @test iszero(X_stabilizers(Q) * transpose(Z_stabilizers(Q)))
    @test length(Q.C1T) == nrows(parity_check_matrix(C))
    @test Q.u_bound == 2
    @test X_minimum_distance_upper_bound(Q) == 2

    H = parity_check_matrix(C)
    @test length(HypergraphProductCode(H)) == 5
    lx, lz = Quintavalle_basis(Q)
    @test size(lx) == size(lz) == (1, 5)

    exact_Q = HypergraphProductCode(C)
    exact_Q.cache[:l_bound_dx] = 1
    exact_Q.cache[:u_bound_dx] = exact_Q.n
    d, witness = minimum_distance(exact_Q; which=:X, alg=:Gray, max_d=3)
    @test d == 2
    @test wt(witness) == d
    isd_Q = HypergraphProductCode(C)
    isd_Q.cache[:l_bound_dx] = 1
    isd_Q.cache[:u_bound_dx] = isd_Q.n
    d_isd, witness_isd = probabilistic_minimum_distance(
        isd_Q; which=:X, alg=:Prange, max_weight=3, max_iters=100, seed=1)
    @test d_isd == 2
    @test wt(witness_isd) == d_isd
end
