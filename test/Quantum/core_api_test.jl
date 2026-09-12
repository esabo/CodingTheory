@testitem "Quantum core structural APIs and constructors" begin
    using Oscar, CodingTheory, Random

    F = GF(2)
    H = matrix(F, [
        0 0 0 1 1 1 1
        0 1 1 0 0 1 1
        1 0 1 0 1 0 1
    ])
    S = CSSCode(H, H)

    @test symplectic_weight(matrix(F, [1 1])) == 1
    @test nrows(normalizer_matrix(S)) == S.n + S.k
    stab = stabilizers(S)[1:1, :]
    @test is_stabilizer(S, stab)
    @test is_normalizer(S, stab)
    @test !is_logical(S, stab)

    @test minimum_stabilizer_weight(S; alg=:bruteforce) == 4
    @test is_pure(S; distance=3, alg=:bruteforce)

    X_degenerate = hcat(H, zero_matrix(F, 3, 1))
    Z_degenerate = vcat(
        X_degenerate, matrix(F, [0 0 0 0 0 0 0 1]))
    degenerate = CSSCode(X_degenerate, Z_degenerate)
    @test minimum_stabilizer_weight(
        degenerate; alg=:bruteforce) == 1
    @test is_degenerate(
        degenerate; distance=3, alg=:bruteforce)

    random_stabilizer =
        random_stabilizer_code(MersenneTwister(1), F, 6, 2)
    @test (random_stabilizer.n, random_stabilizer.k) == (6, 2)
    @test rank(stabilizers(random_stabilizer)) == 4
    @test are_symplectic_orthogonal(
        stabilizers(random_stabilizer), stabilizers(random_stabilizer))

    trivial = random_stabilizer_code(MersenneTwister(2), F, 4, 4)
    @test trivial.k == 4
    @test rank(stabilizers(trivial)) == 0

    random_subsystem =
        random_subsystem_code(MersenneTwister(3), F, 6, 2, 2)
    @test (random_subsystem.n, random_subsystem.k,
        random_subsystem.r) == (6, 2, 2)
    @test rank(stabilizers(random_subsystem)) == 2
    @test rank(gauge_group(random_subsystem)) == 6

    F4 = GF(4, :b)
    qary = random_subsystem_code(MersenneTwister(4), F4, 5, 2, 1)
    @test (qary.k, qary.r) == (2, 1)
    @test CodingTheory._additive_rank(stabilizers(qary), F4) == 5
    @test CodingTheory._additive_rank(gauge_group(qary), F4) == 7

    half_integral =
        random_stabilizer_code(MersenneTwister(5), F4, 2, 1 // 2)
    @test half_integral.k == 1 // 2
    @test CodingTheory._additive_rank(
        stabilizers(half_integral), F4) == 3
    @test are_symplectic_orthogonal(
        stabilizers(half_integral), stabilizers(half_integral))

    E = GF(4, :a)
    C = LinearCode(matrix(E, [1 1 0]))
    quadratic = StabilizerCode(C, F)
    @test (quadratic.n, quadratic.k) == (3, 1)
    @test rank(stabilizers(quadratic)) == 2
    @test are_symplectic_orthogonal(
        stabilizers(quadratic), stabilizers(quadratic))

    C_gauge = LinearCode(matrix(E, [
        1 1 0
        0 0 1
    ]))
    quadratic_subsystem = SubsystemCode(C_gauge, F)
    @test (quadratic_subsystem.n, quadratic_subsystem.k,
        quadratic_subsystem.r) == (3, 0, 1)
    @test rank(stabilizers(quadratic_subsystem)) == 2
    @test rank(gauge_group(quadratic_subsystem)) == 4
end

@testitem "Sparse quantum constructors preserve storage" begin
    using Oscar, CodingTheory, SparseArrays

    H = sparse([
        1 1 0 0
        0 0 1 1
    ])
    S = CSSCode(H, H)
    @test stabilizers(S) isa SparseMatrixCSC
    @test X_stabilizers(S) isa SparseMatrixCSC
    @test Z_stabilizers(S) isa SparseMatrixCSC
    @test S.k == 0
    @test are_symplectic_orthogonal(
        stabilizers(S), stabilizers(S))
    @test stabilizer_weights(S) == fill(2, 4)
    @test qubit_degrees(S) == fill(2, 4)

    F = GF(2)
    H_oscar = sparse_matrix(F, Matrix(H))
    @test CodingTheory._quantum_degree_distribution(H_oscar, 4) ==
        (fill(1, 4), fill(2, 2))
    S_oscar = CSSCode(H_oscar, H_oscar)
    @test stabilizers(S_oscar) isa SparseMatrixCSC
    @test stabilizer_weights(S_oscar) == fill(2, 4)

    F4 = GF(4)
    additive_generator =
        sparse_matrix(F4, [F4(1) F4(0)])
    qary_sparse = StabilizerCode(
        additive_generator; logs_alg=:sys_eqs)
    @test stabilizers(qary_sparse) isa SMat
    @test qary_sparse.k == 1 // 2
    @test length(logicals(qary_sparse)) == 1

    G = sparse([
        1 1 0 0  0 0 0 0
        0 0 1 1  0 0 0 0
        0 0 0 0  1 0 1 0
        0 0 0 0  0 1 0 1
    ])
    Q = SubsystemCode(G)
    @test stabilizers(Q) isa SparseMatrixCSC
    @test gauge_group(Q) isa SparseMatrixCSC
    @test normalizer_matrix(Q) isa SparseMatrixCSC
end

@testitem "Quantum CSS subsystem and gauge fixing" begin
    using Oscar, CodingTheory

    F = GF(2)
    X_gauges = matrix(F, [
        1 1 0 0
        0 0 1 1
    ])
    Z_gauges = matrix(F, [
        1 0 1 0
        0 1 0 1
    ])
    Q = CSSSubsystemCode(X_gauges, Z_gauges)

    @test is_CSS(Q)
    @test (Q.n, Q.k, Q.r) == (4, 1, 1)
    @test rank(stabilizers(Q)) == 2
    @test rank(gauge_group(Q)) == 4
    @test is_gauge(Q, gauges_matrix(Q)[1:1, :])

    fixed_X = fix_all_gauges(Q; choice=:X)
    fixed_Z = fix_all_gauges(Q; choice=:Z)
    for fixed in (fixed_X, fixed_Z)
        @test fixed isa AbstractStabilizerCode
        @test is_CSS(fixed)
        @test fixed.k == 1
        @test rank(stabilizers(fixed)) == 3
        @test are_symplectic_orthogonal(
            stabilizers(fixed), stabilizers(fixed))
        @test logicals_matrix(fixed) == logicals_matrix(Q)
        @test haskey(fixed.cache, :gauge_fixed_from)
    end
end

@testitem "Quantum LDPC presentation statistics" begin
    using Oscar, CodingTheory

    F = GF(2)
    H = matrix(F, [
        0 0 0 1 1 1 1
        0 1 1 0 0 1 1
        1 0 1 0 1 0 1
    ])
    S = CSSCode(H, H)

    @test stabilizer_weights(S) == fill(4, 6)
    @test qubit_degrees(S) == [2, 2, 4, 2, 4, 4, 6]
    @test X_qubit_degrees(S) == [1, 1, 2, 1, 2, 2, 3]
    @test check_weights(S) == (4, 3, 4, 3)
    @test quantum_LDPC_parameters(S) ==
        (max_generator_weight=4, max_qubit_degree=6)
    @test is_quantum_LDPC(
        S; max_generator_weight=4, max_qubit_degree=6)
    @test is_LDPC(S; check_bound=4, column_bound=3)
    @test !is_LDPC(S; check_bound=3, column_bound=3)

    @test variable_degree_distribution(H) ==
        variable_degree_distribution(LDPCCode(H))
    @test check_degree_distribution(H) ==
        check_degree_distribution(LDPCCode(H))
    @test column_row_bounds(H) == (3, 4)
    @test check_weights(H) == (4, 3)
end
