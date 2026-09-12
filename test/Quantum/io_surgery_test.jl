@testitem "Quantum portable I/O round trips" begin
    using Oscar, CodingTheory, SparseArrays

    H = sparse([
        1 1 0 0
        0 0 1 1
    ])
    S = CSSCode(H, H; logs_alg=:sys_eqs)
    S.cache[:d] = 2
    SL_weight_enumerator(S)

    payload = quantum_code_data(S)
    restored = quantum_code_from_data(payload)
    @test (restored.n, restored.k) == (S.n, S.k)
    @test is_CSS(restored)
    @test CodingTheory._additive_row_spaces_equal(
        stabilizers(restored), stabilizers(S), S.F)
    @test restored.cache[:d] == 2
    @test restored.cache[:SL_weight_enum].A.counts ==
        S.cache[:SL_weight_enum].A.counts
    @test stabilizers(restored) isa SparseMatrixCSC
    uncached = quantum_code_from_data(payload; restore_cache=false)
    @test !haskey(uncached.cache, :d)
    @test !haskey(uncached.cache, :SL_weight_enum)

    tampered_cache = deepcopy(payload)
    tampered_cache["cache"]["d"] = 1
    @test_throws ArgumentError quantum_code_from_data(tampered_cache)
    tampered_generators = deepcopy(payload)
    tampered_generators["generators"]["data"][1] ⊻= 1
    @test_throws ArgumentError quantum_code_from_data(tampered_generators)

    mktempdir() do directory
        path = joinpath(directory, "code.toml")
        @test save_code(path, S; type=:toml) == path
        loaded = load_quantum_code(path)
        @test CodingTheory._additive_row_spaces_equal(
            stabilizers(loaded), stabilizers(S), S.F)

        pauli_path = joinpath(directory, "code.pauli")
        @test_throws ArgumentError save_code(pauli_path, S; type=:pauli)
        write_pauli_strings(pauli_path, S)
        pauli_loaded = read_pauli_strings(pauli_path)
        @test CodingTheory._additive_row_spaces_equal(
            stabilizers(pauli_loaded), stabilizers(S), S.F)

        csv_path = joinpath(directory, "generators.csv")
        @test save_code(csv_path, S) == csv_path
        csv_rows = [
            parse.(Int, split(line, ','))
            for line in readlines(csv_path)
        ]
        @test reduce(vcat, permutedims.(csv_rows)) ==
            quantum_generator_array(S)

        signed_path = joinpath(directory, "signed.pauli")
        write(signed_path, "-XXII\n")
        @test_throws ArgumentError read_pauli_strings(signed_path)
    end

    @test pauli_strings(S) == ["XXII", "IIXX", "ZZII", "IIZZ"]

    G = sparse([
        1 1 0 0  0 0 0 0
        0 0 1 1  0 0 0 0
        0 0 0 0  1 0 1 0
        0 0 0 0  0 1 0 1
    ])
    Q = SubsystemCode(G)
    Q_loaded = quantum_code_from_data(quantum_code_data(Q))
    @test (Q_loaded.n, Q_loaded.k, Q_loaded.r) == (Q.n, Q.k, Q.r)
    @test CodingTheory._additive_row_spaces_equal(
        gauge_group(Q_loaded), gauge_group(Q), Q.F)

    F4 = GF(4)
    qary = StabilizerCode(
        sparse_matrix(F4, [F4(1) F4(0)]); logs_alg=:sys_eqs)
    qary_loaded = quantum_code_from_data(quantum_code_data(qary))
    @test qary_loaded.k == 1 // 2
    @test CodingTheory._additive_row_spaces_equal(
        stabilizers(qary_loaded), stabilizers(qary), qary.F)
    primitive_qary = StabilizerCode(
        matrix(F4, [gen(F4) F4(0)]); logs_alg=:sys_eqs)
    primitive_loaded =
        quantum_code_from_data(quantum_code_data(primitive_qary))
    @test CodingTheory._additive_row_spaces_equal(
        stabilizers(primitive_loaded), stabilizers(primitive_qary), F4)
    primitive_array = quantum_generator_array(primitive_qary)
    @test size(primitive_array) == (1, 4)
    @test eltype(primitive_array) == Int

    R4, _ = residue_ring(ZZ, 4)
    phased = StabilizerCode(
        matrix(S.F, [1 1 0 0]);
        char_vec=[R4(1), R4(0), R4(0), R4(0)],
        logs_alg=:sys_eqs)
    mktemp() do path, _
        @test_throws ArgumentError write_pauli_strings(path, phased)
    end
end

@testitem "Quantum new-codes-from-old operations" begin
    using Oscar, CodingTheory, SparseArrays

    F = GF(2)
    bell = StabilizerCode(matrix(F, [
        1 1  0 0
        0 0  1 1
    ]); logs_alg=:sys_eqs)

    doubled = quantum_direct_sum(bell, bell)
    @test (doubled.n, doubled.k) == (4, 0)
    @test nrows(stabilizers(doubled)) == 4
    @test are_symplectic_orthogonal(
        stabilizers(doubled), stabilizers(doubled))
    @test doubled ⊕ bell isa AbstractStabilizerCode

    bell.cache[:d] = 2
    transformed = local_fourier(bell, [1])
    restored = local_fourier(transformed, [1])
    @test CodingTheory._additive_row_spaces_equal(
        stabilizers(restored), stabilizers(bell), F)
    @test transformed.cache[:d] == 2

    punctured = puncture(bell, [2])
    @test punctured.n == 1
    @test GaugeTrait(typeof(punctured)) == HasGauges()
    @test punctured.r == 1
    @test puncture(bell, 2).n == 1

    shortened = shorten(bell, [2])
    @test shortened isa AbstractStabilizerCode
    @test (shortened.n, shortened.k) == (1, 1)
    @test shorten(bell, 2).n == 1

    measured = augment(
        bell, matrix(F, [0 0  1 0]); verbose=false)
    @test measured isa AbstractStabilizerCode
    @test measured.k == 0
    @test are_symplectic_orthogonal(
        stabilizers(measured), stabilizers(measured))
    measured_column = augment(
        bell, matrix(F, 4, 1, [0, 0, 1, 0]); verbose=false)
    @test measured_column.k == 0
    measured_sparse = augment(
        bell, sparse([0 0 1 0]); verbose=false)
    @test measured_sparse.k == 0

    expanded = expurgate(bell, 1; verbose=false)
    @test expanded isa AbstractStabilizerCode
    @test expanded.k == 1

    seed = StabilizerCode(matrix(F, [0 0  1 1]); logs_alg=:sys_eqs)
    gauged = gauge_code(seed, matrix(F, [1 0  0 0]))
    @test GaugeTrait(typeof(gauged)) == HasGauges()
    @test (gauged.n, gauged.k, gauged.r) == (2, 1, 1)

    F3 = GF(3)
    qutrit = StabilizerCode(
        matrix(F3, [1 0]); logs_alg=:sys_eqs)
    transformed_qutrit = foldl(
        (code, _) -> local_fourier(code, 1), 1:4; init=qutrit)
    @test CodingTheory._additive_row_spaces_equal(
        stabilizers(transformed_qutrit), stabilizers(qutrit), F3)
end
