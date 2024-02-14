@testset "utils.jl" begin
    # using Oscar, CodingTheory

    # NOTE: circshift is currently commented out, might be deleted in the future
    # F = GF(2,1)
    # v = matrix(F, 1, 8, [1, 0, 1, 1, 1, 0, 0, 0])
    # @test circshift(v, 1) == matrix(F, 1, 8, [0, 1, 0, 1, 1, 1, 0, 0])
    # @test circshift(v, 2) == matrix(F, 1, 8, [0, 0, 1, 0, 1, 1, 1, 0])
    # @test circshift(v, -3) == matrix(F, 1, 8, [1, 1, 0, 0, 0, 1, 0, 1])
    # @test circshift(v, 10) == matrix(F, 1, 8, [0, 0, 1, 0, 1, 1, 1, 0])
    # @test circshift(v, -11) == matrix(F, 1, 8, [1, 1, 0, 0, 0, 1, 0, 1])

    GF4 = GF(2, 2, :ω)
    GF2 = GF(2)

    M_min_wt = [1 0 0 1 0
              1 1 0 0 0
              0 0 1 0 0
              1 1 1 1 1
              0 1 0 1 0]
    M_min_wt2 = matrix(GF2, M_min_wt)
    M_min_wt3 = matrix(GF4, M_min_wt)
    @test CodingTheory._min_wt_row(M_min_wt) == (1, 3)
    @test CodingTheory._min_wt_col(M_min_wt) == (1, 5)
    @test CodingTheory._min_wt_row(M_min_wt2) == (1, 3)
    @test CodingTheory._min_wt_col(M_min_wt2) == (1, 5)
    @test CodingTheory._min_wt_row(M_min_wt3) == (1, 3)
    @test CodingTheory._min_wt_col(M_min_wt3) == (1, 5)

    v = [0, 1, 1, 0, 1, 1]
    w = [0, 1, 0, 1, 1, 1]
    v2 = matrix(GF2, 1, 6, v)
    w2 = matrix(GF2, 1, 6, w)
    v3 = matrix(GF4, 1, 6, v)
    w3 = matrix(GF4, 1, 6, w)
    @test Hamming_distance(v, w) == 2
    @test Hamming_distance(v2, w2) == 2
    @test Hamming_distance(v3, w3) == 2

    @test symplectic_inner_product(v2, w2) == 1
    @test symplectic_inner_product(v3, w3) == 1
    @test symplectic_inner_product(v3, v3) == 0
    @test symplectic_inner_product(w3, w3) == 0
    @test are_symplectic_orthogonal(v3, v3)
    @test are_symplectic_orthogonal(w3, w3)

    ω = gen(GF4)
    hexacode = matrix(GF4, [1 0 0 1 ω ω; 0 1 0 ω 1 ω; 0 0 1 ω ω 1])
    @test Hermitian_inner_product(hexacode[1, :], matrix(GF4, [1 0 0 1 1 0])) == ω
    @test Hermitian_inner_product(hexacode[1, :], hexacode[2, :]) == 0
    @test iszero(Hermitian_conjugate_matrix(hexacode) * transpose(hexacode))

    # _remove_empty
    M = ones(Int, rand(20:30), rand(20:30))
    row_index = rand(1:size(M, 1))
    colindex = rand(1:size(M, 2))
    for i in axes(M, 1)
        M[i, colindex] = 0
    end
    for j in axes(M, 2)
        M[row_index, j] = 0
    end
    M2 = matrix(GF2, M)
    M3 = matrix(GF4, M)
    M2_rem_row = CodingTheory._remove_empty(M2, :rows)
    M2_rem_col = CodingTheory._remove_empty(M2, :cols)
    M3_rem_row = CodingTheory._remove_empty(M3, :rows)
    M3_rem_col = CodingTheory._remove_empty(M3, :cols)
    @test !any(iszero(M2_rem_row[i, :]) for i in axes(M2_rem_row, 1))
    @test !any(iszero(M2_rem_col[:, j]) for j in axes(M2_rem_col, 2))
    @test !any(iszero(M3_rem_row[i, :]) for i in axes(M3_rem_row, 1))
    @test !any(iszero(M3_rem_col[:, j]) for j in axes(M3_rem_col, 2))

    # TODO: _rref_no_col_swap and _rref_col_swap - come back to when going over weightdist.jl

    # digits_to_int
    @test all(d == digits(d, base = 2, pad = 15) |> reverse |> digits_to_int for d in rand(0:2^15, 100))

    # _concat
    locations = [0 1; 1 1]
    M1 = matrix(GF2, ones(Int, 3, 2))
    M2 = matrix(GF4, ones(Int, 3, 2))
    @test CodingTheory._concat(locations, M1) == matrix(GF2, [0 0 1 1; 0 0 1 1; 0 0 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1])
    @test CodingTheory._concat(locations, M2) == matrix(GF4, [0 0 1 1; 0 0 1 1; 0 0 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1])

    # TODO: pseudoinverse test

    # Tri-orthogonal matrix from Bravyi and Haah 2012, equation 3
    M_tri_orth = [1 1 1 1 1 1 1 0 0 0 0 0 0 0
                0 0 0 0 0 0 0 1 1 1 1 1 1 1
                1 0 1 0 1 0 1 1 0 1 0 1 0 1
                0 1 1 0 0 1 1 0 1 1 0 0 1 1
                0 0 0 1 1 1 1 0 0 0 1 1 1 1]
    @test is_triorthogonal(M_tri_orth) # test for Matrix{Int}
    @test is_triorthogonal(matrix(GF2, M_tri_orth)) # test for fpMatrix

    # example: Betten et al
    # Golay code G_23
    qres, _ = quadratic_residues(2, 23)
    @test qres == [1, 2, 3, 4, 6, 8, 9, 12, 13, 16, 18]

    # example: Betten et al
    # ternary Golary code G_11
    qres, _ = quadratic_residues(3, 11)
    @test qres == [1, 3, 4, 5, 9]

    F = GF(2)
    E = GF(2, 3, :α)
    α = gen(E)
    flag, _ = is_extension(E, F)
    @test flag
    basis = [α^3, α^5, α^6];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis, _ = primitive_basis(E, F)
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis = normal_basis(E, F)
    flag, _ = is_basis(E, F, basis[1])
    @test flag
    @test verify_dual_basis(E, F, basis[1], basis[2])

    E = GF(2, 4, :α)
    α = gen(E)
    flag, _ = is_extension(E, F)
    @test flag
    basis = [α^3, α^6, α^9, α^12];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis = [α^7, α^11, α^13, α^14]
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis, _ = primitive_basis(E, F)
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis, _ = normal_basis(E, F)
    flag, _ = is_basis(E, F, basis)
    @test flag

    F = GF(3)
    E = GF(3, 2, :α)
    α = gen(E)
    flag, _ = is_extension(E, F)
    @test flag
    basis = [α, α^3];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis = [α^5, α^7];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis, _ = primitive_basis(E, F)
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis, _ = normal_basis(E, F)
    flag, _ = is_basis(E, F, basis)
    @test flag

    E = GF(3, 3, :α);
    α = gen(E)
    flag, _ = is_extension(E, F)
    @test flag
    basis = [α^2, α^6, α^18];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis = [α^4, α^10, α^12];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis = [α^5, α^15, α^19];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis = [α^7, α^11, α^21];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis = [α^8, α^20, α^24];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis = [α^17, α^23, α^25];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis, _ = primitive_basis(E, F)
    flag, _ = is_basis(E, F, basis)
    @test flag
    @test is_primitive_basis(E, F, basis)
    basis, _ = normal_basis(E, F)
    flag, _ = is_basis(E, F, basis)
    @test flag
    @test is_normal_basis(E, F, basis)

    F = GF(5)
    E = GF(5, 2, :α);
    α = gen(E)
    flag, _ = is_extension(E, F)
    @test flag
    basis = [α, α^5];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis = [α^2, α^10];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis = [α^4, α^20];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis = [α^7, α^11];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis = [α^8, α^16];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis = [α^13, α^17];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis = [α^14, α^22];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis = [α^19, α^23];
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis, _ = primitive_basis(E, F)
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis, _ = normal_basis(E, F)
    flag, _ = is_basis(E, F, basis)
    @test flag
    basis2 = [α * basis[i] for i in 1:2]
    @test are_equivalent_basis(basis, basis2)

    # TODO: work these in
    # F8 = GF(2, 3, :α)
    # β = [F8(1), α, α^6]
    # λ = dual_basis(F8, F, β)
    # D = CodingTheory._expansion_dict(F8, F, λ)
    # D2 = typeof(D)()
    # D2[F8(0)] = zero_matrix(F8, 1, 3)
    # D2[F8(1)] = matrix(F8, 1, 3, [1, 0, 0])
    # D2[α] = matrix(F8, 1, 3, [0, 1, 0])
    # D2[α^2] = matrix(F8, 1, 3, [1, 0, 1])
    # D2[α^3] = matrix(F8, 1, 3, [1, 1, 0])
    # D2[α^4] = matrix(F8, 1, 3, [1, 1, 1])
    # D2[α^5] = matrix(F8, 1, 3, [0, 1, 1])
    # D2[α^6] = matrix(F8, 1, 3, [0, 0, 1])
    # @test D == D2

    # β = [α^3, α^6, α^5]
    # λ = dual_basis(F8, F, β)
    # D = CodingTheory._expansion_dict(F8, F, λ)
    # D2 = typeof(D)()
    # D2[F8(0)] = zero_matrix(F8, 1, 3)
    # D2[F8(1)] = matrix(F8, 1, 3, [1, 1, 1])
    # D2[α] = matrix(F8, 1, 3, [0, 1, 1])
    # D2[α^2] = matrix(F8, 1, 3, [1, 0, 1])
    # D2[α^3] = matrix(F8, 1, 3, [1, 0, 0])
    # D2[α^4] = matrix(F8, 1, 3, [1, 1, 0])
    # D2[α^5] = matrix(F8, 1, 3, [0, 0, 1])
    # D2[α^6] = matrix(F8, 1, 3, [0, 1, 0])
    # @test D == D2


    F = GF(2)
    flag, _ = is_extension(E, F)
    @test !flag

    F = GF(2)
    S, x = PolynomialRing(F, :x)
    l = 3
    R = residue_ring(S, x^l - 1)
    A = matrix(R, 2, 3, [1, 0, 1 + x^2, 1 + x, 1 + x + x^2, x^2])
    @test lift(A) == matrix(F, 6, 9,
        [1, 0, 0, 0, 0, 0, 1, 1, 0,
        0, 1, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 1, 0, 0, 0, 1, 0, 1,
        1, 0, 1, 1, 1, 1, 0, 1, 0,
        1, 1, 0, 1, 1, 1, 0, 0, 1,
        0, 1, 1, 1, 1, 1, 1, 0, 0])
    @test weight_matrix(A) == [1 0 2; 2 3 1]
end
