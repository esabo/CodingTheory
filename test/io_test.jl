@testitem "Classical and LDPC CSV exports" begin
    using CodingTheory, Oscar, SparseArrays

    function read_integer_csv(path)
        rows = [
            parse.(Int, split(line, ','))
            for line in readlines(path)
        ]
        return reduce(vcat, permutedims.(rows))
    end

    C = HammingCode(2, 3)
    generator_values = code_matrix_array(C)
    @test size(generator_values) == (C.k, C.n)
    G = generator_matrix(C)
    @test generator_values == [
        Int(lift(ZZ, G[r, c]))
        for r in 1:nrows(G), c in 1:ncols(G)
    ]

    F = GF(2)
    H = sparse([
        1 1 0 1
        0 1 1 1
    ])
    L = LDPCCode(matrix(F, Matrix(H)))
    @test code_matrix_array(L) == Matrix(H)
    L_generator = code_matrix_array(L; representation=:generator)
    @test size(L_generator, 2) == L.n
    @test iszero(matrix(F, L_generator) *
        transpose(parity_check_matrix(L)))

    F4 = GF(4)
    extension_code = LinearCode(matrix(F4, [one(F4) gen(F4)]))
    @test size(code_matrix_array(extension_code)) == (1, 4)

    mktempdir() do directory
        classical_path = joinpath(directory, "classical.csv")
        ldpc_path = joinpath(directory, "ldpc.csv")
        @test save_code(classical_path, C; type=:csv) == classical_path
        @test save_code(ldpc_path, L) == ldpc_path
        @test read_integer_csv(classical_path) == generator_values
        @test read_integer_csv(ldpc_path) == Matrix(H)
    end

    @test_throws ArgumentError save_code(
        "unsupported.toml", C; type=:toml)
end
