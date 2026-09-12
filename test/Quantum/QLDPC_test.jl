@testitem "Quantum LDPC matrix statistics" begin
    using Oscar, CodingTheory

    F = GF(2)
    H_X = matrix(F, [
        1 1 1 0 0 0
        0 0 1 1 1 1
    ])
    H_Z = matrix(F, [
        1 1 0 0 0 0
        1 0 1 1 0 0
        0 0 0 0 1 1
    ])
    S = CSSCode(H_X, H_Z)

    matrix_degrees(H) = (
        [count(!iszero, H[:, j]) for j in 1:ncols(H)],
        [count(!iszero, H[i, :]) for i in 1:nrows(H)],
    )
    x_columns, x_rows = matrix_degrees(H_X)
    z_columns, z_rows = matrix_degrees(H_Z)
    all_columns = x_columns + z_columns
    all_rows = vcat(x_rows, z_rows)
    x_edges = sum(x_rows)
    z_edges = sum(z_rows)

    @test x_columns == [1, 1, 2, 1, 1, 1]
    @test x_rows == [3, 4]
    @test z_columns == [2, 1, 1, 1, 1, 1]
    @test z_rows == [2, 3, 2]

    @test X_qubit_degrees(S) == x_columns
    @test Z_qubit_degrees(S) == z_columns
    @test qubit_degrees(S) == all_columns
    @test qubit_degree_distribution(S) == all_columns
    @test X_stabilizer_weights(S) == x_rows
    @test Z_stabilizer_weights(S) == z_rows
    @test stabilizer_weights(S) == all_rows
    @test stabilizer_weight_distribution(S) == all_rows
    @test generator_weights(S) == all_rows

    @test X_variable_degree_distribution(S) == x_columns
    @test Z_variable_degree_distribution(S) == z_columns
    @test X_check_degree_distribution(S) == x_rows
    @test Z_check_degree_distribution(S) == z_rows
    @test X_degree_distributions(S) == (x_columns, x_rows)
    @test Z_degree_distributions(S) == (z_columns, z_rows)

    @test X_column_bound(S) == maximum(x_columns)
    @test Z_column_bound(S) == maximum(z_columns)
    @test X_row_bound(S) == maximum(x_rows)
    @test Z_row_bound(S) == maximum(z_rows)
    @test X_column_row_bounds(S) ==
        (maximum(x_columns), maximum(x_rows))
    @test Z_column_row_bounds(S) ==
        (maximum(z_columns), maximum(z_rows))
    @test X_limited(S) == max(maximum(x_columns), maximum(x_rows))
    @test Z_limited(S) == max(maximum(z_columns), maximum(z_rows))

    @test X_density(S) == x_edges / length(H_X)
    @test Z_density(S) == z_edges / length(H_Z)
    @test num_edges(H_X) == x_edges
    @test num_edges(H_Z) == z_edges
    @test num_edges(S) == x_edges + z_edges
    @test check_weights(S) == (
        maximum(x_rows), maximum(x_columns),
        maximum(z_rows), maximum(z_columns))
    @test maximum_qubit_degree(S) == maximum(all_columns)
    @test minimum_qubit_degree(S) == minimum(all_columns)
    @test maximum_stabilizer_weight(S) == maximum(all_rows)

    parameters = quantum_LDPC_parameters(S)
    @test parameters == (
        max_generator_weight=maximum(all_rows),
        max_qubit_degree=maximum(all_columns))
    @test is_quantum_LDPC(S;
        max_generator_weight=maximum(all_rows),
        max_qubit_degree=maximum(all_columns))
    @test !is_quantum_LDPC(S;
        max_generator_weight=maximum(all_rows) - 1,
        max_qubit_degree=maximum(all_columns))
    @test !is_quantum_LDPC(S;
        max_generator_weight=maximum(all_rows),
        max_qubit_degree=maximum(all_columns) - 1)

    @test is_X_LDPC(S;
        check_bound=maximum(x_rows), column_bound=maximum(x_columns))
    @test !is_X_LDPC(S;
        check_bound=maximum(x_rows) - 1, column_bound=maximum(x_columns))
    @test !is_X_LDPC(S;
        check_bound=maximum(x_rows), column_bound=maximum(x_columns) - 1)
    @test is_Z_LDPC(S;
        check_bound=maximum(z_rows), column_bound=maximum(z_columns))
    @test !is_Z_LDPC(S;
        check_bound=maximum(z_rows) - 1, column_bound=maximum(z_columns))
    @test !is_Z_LDPC(S;
        check_bound=maximum(z_rows), column_bound=maximum(z_columns) - 1)
    @test is_LDPC(S;
        check_bound=max(maximum(x_rows), maximum(z_rows)),
        column_bound=max(maximum(x_columns), maximum(z_columns)))
    @test !is_LDPC(S;
        check_bound=max(maximum(x_rows), maximum(z_rows)) - 1,
        column_bound=max(maximum(x_columns), maximum(z_columns)))
    @test !is_LDPC(S;
        check_bound=max(maximum(x_rows), maximum(z_rows)),
        column_bound=max(maximum(x_columns), maximum(z_columns)) - 1)

    C_X, C_Z = LDPC_codes(S)
    @test parity_check_matrix(C_X) == H_X
    @test parity_check_matrix(C_Z) == H_Z
    @test parity_check_matrix(X_LDPC_code(S)) == H_X
    @test parity_check_matrix(Z_LDPC_code(S)) == H_Z

    H_regular = matrix(F, [
        1 1 0 0
        0 0 1 1
    ])
    regular = CSSCode(H_regular, H_regular)
    @test X_is_regular(regular)
    @test Z_is_regular(regular)
    @test !X_is_regular(S)
    @test !Z_is_regular(S)

    non_css = FiveQubitCode()
    @test_throws ErrorException X_qubit_degrees(non_css)
    @test_throws ErrorException Z_stabilizer_weights(non_css)
    @test_throws ErrorException X_degree_distributions(non_css)
    @test_throws ErrorException Z_density(non_css)
    @test_throws ErrorException X_LDPC_code(non_css)
    @test_throws ErrorException LDPC_codes(non_css)
end

@testitem "Quantum LDPC subsystem weights" begin
    using Oscar, CodingTheory

    pauli_weights(M, n) = [
        count(q -> !iszero(M[i, q]) || !iszero(M[i, n + q]), 1:n)
        for i in 1:nrows(M)
    ]
    pauli_degrees(M, n) = [
        count(i -> !iszero(M[i, q]) || !iszero(M[i, n + q]), 1:nrows(M))
        for q in 1:n
    ]

    Q = BaconShorCode(2, 3)
    gauges = gauges_matrix(Q)
    group = gauge_group(Q)
    stabilizer_matrix = stabilizers(Q)
    expected_gauge_weights = pauli_weights(gauges, Q.n)
    expected_group_weights = pauli_weights(group, Q.n)
    expected_stabilizer_weights = pauli_weights(stabilizer_matrix, Q.n)
    expected_gauge_degrees = pauli_degrees(gauges, Q.n)

    @test gauge_weights(Q) == expected_gauge_weights
    @test gauge_group_weights(Q) == expected_group_weights
    @test generator_weights(Q; generators=:gauges) ==
        expected_gauge_weights
    @test generator_weights(Q; generators=:gauge_group) ==
        expected_group_weights
    @test generator_weights(Q; generators=:stabilizers) ==
        expected_stabilizer_weights
    @test qubit_degrees(Q; generators=:gauges) == expected_gauge_degrees
    @test quantum_LDPC_parameters(Q; generators=:gauges) == (
        max_generator_weight=maximum(expected_gauge_weights),
        max_qubit_degree=maximum(expected_gauge_degrees))
    @test is_quantum_LDPC(Q;
        max_generator_weight=maximum(expected_gauge_weights),
        max_qubit_degree=maximum(expected_gauge_degrees),
        generators=:gauges)
    @test !is_quantum_LDPC(Q;
        max_generator_weight=maximum(expected_gauge_weights) - 1,
        max_qubit_degree=maximum(expected_gauge_degrees),
        generators=:gauges)
end
