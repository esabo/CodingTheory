@testitem "Linear code duality" begin
    using Oscar, CodingTheory, Random

    euclidean_gram(G) = G * transpose(G)
    function galois_gram(G, ell)
        p = Int(characteristic(base_ring(G)))
        return G * transpose(G .^ (p^ell))
    end
    rows_are_orthogonal(G, ell=0) = iszero(galois_gram(G, ell))

    @testset "Classical Euclidean examples" begin
        ext_hamming = ExtendedHammingCode(3)
        golay = ExtendedGolayCode(2)
        hamming = HammingCode(2, 3)
        tetra = TetraCode()

        for C in (ext_hamming, golay, tetra)
            G = generator_matrix(C)
            @test rows_are_orthogonal(G)
            @test dimension(C) == length(C) ÷ 2
            @test is_weakly_self_dual(C)
            @test is_Euclidean_self_orthogonal(C)
            @test is_self_dual(C)
            @test are_equivalent(Euclidean_dual(C), dual(C))
        end

        G_hamming = generator_matrix(hamming)
        H_hamming = parity_check_matrix(hamming)
        @test iszero(G_hamming * transpose(H_hamming))
        @test all(iszero(sum(H_hamming[i, t] * H_hamming[j, t]
                             for t in axes(H_hamming, 2)))
                  for i in axes(H_hamming, 1), j in axes(H_hamming, 1))
        @test is_dual_containing(hamming)
        @test is_Euclidean_dual_containing(hamming)
        @test !is_weakly_self_dual(hamming)

        @test is_doubly_even(golay)
        @test !is_triply_even(golay)
        @test is_triply_even(RepetitionCode(2, 8))
        @test contains_self_dual_subcode(ext_hamming)
        @test contains_self_dual_subcode(tetra)
        @test !contains_self_dual_subcode(hamming)
    end

    @testset "Hermitian and Galois forms" begin
        hexacode = Hexacode()
        G_hex = generator_matrix(hexacode)
        @test iszero(galois_gram(G_hex, 1))
        @test !iszero(euclidean_gram(G_hex))
        @test is_Hermitian_self_orthogonal(hexacode)
        @test is_Hermitian_weakly_self_dual(hexacode)
        @test is_Hermitian_self_dual(hexacode)
        @test is_Hermitian_dual_containing(hexacode)
        @test are_equivalent(Hermitian_dual(hexacode), l_Galois_dual(hexacode, 1))
        @test !is_self_orthogonal(hexacode)

        F4 = field(hexacode)
        ω = gen(F4)
        C_euclidean = LinearCode(matrix(F4, [1 ω ω + 1]))
        C_hermitian = LinearCode(matrix(F4, [1 ω]))

        @test rows_are_orthogonal(generator_matrix(C_euclidean), 0)
        @test !rows_are_orthogonal(generator_matrix(C_euclidean), 1)
        @test !rows_are_orthogonal(generator_matrix(C_hermitian), 0)
        @test rows_are_orthogonal(generator_matrix(C_hermitian), 1)

        @test is_l_Galois_self_orthogonal(C_euclidean, 0)
        @test is_l_Galois_weakly_self_dual(C_euclidean, 0)
        @test !is_l_Galois_self_orthogonal(C_euclidean, 1)
        @test !is_l_Galois_self_orthogonal(C_hermitian, 0)
        @test is_l_Galois_self_orthogonal(C_hermitian, 1)
        @test is_Hermitian_self_orthogonal(C_hermitian)
        @test is_l_Galois_self_dual(C_hermitian, 1)
        @test is_l_Galois_dual_containing(hexacode, 1)
        @test is_l_Galois_self_dual(hexacode, 1)

        @test are_equivalent(l_Galois_dual(C_euclidean, 0),
                             Euclidean_dual(C_euclidean))
        @test are_equivalent(l_Galois_dual(C_hermitian, 1),
                             Hermitian_dual(C_hermitian))
    end

    @testset "Hulls from independent Gram ranks" begin
        F3 = GF(3)
        C = LinearCode(matrix(F3, [1 0 0; 0 1 0]))
        G = generator_matrix(C)
        expected_dim = dimension(C) - rank(euclidean_gram(G))
        C_hull, dim_hull = Euclidean_hull(C)
        @test expected_dim == 0
        @test dim_hull == expected_dim
        @test ismissing(C_hull)

        F4 = GF(2, 2, :ω)
        ω = gen(F4)
        C4 = LinearCode(matrix(F4, [1 ω ω + 1]))
        G4 = generator_matrix(C4)
        expected = dimension(C4) - rank(galois_gram(G4, 1))
        H4, actual = l_Galois_hull(C4, 1)
        @test expected == 0
        @test actual == expected
        @test ismissing(H4)
    end

    @testset "Code metadata and representations" begin
        tetra = TetraCode()
        @test genus(tetra) == 0
        @test ismissing(genus(LinearCode(matrix(GF(2), [1 0 1]))))

        V, embedding = vector_space(tetra)
        @test AbstractAlgebra.dim(V) == dimension(tetra)
        @test all(parent(embedding(g)) == codomain(embedding) for g in gens(V))

        F2 = GF(2)
        G = matrix(F2, [0 1 1 0; 0 0 1 1])
        C = LinearCode(G)
        G_standard = generator_matrix(C, true)
        P = standard_form_permutation(C)
        @test rank(P) == ncols(P)
        @test all(Hamming_weight(P[i:i, :]) == 1 for i in axes(P, 1))
        @test all(Hamming_weight(P[:, j:j]) == 1 for j in axes(P, 2))
        @test CodingTheory._has_equivalent_row_spaces(G, G_standard * P)

        C_changed = deepcopy(C)
        F4 = GF(2, 2, :a)
        @test isnothing(change_field!(C_changed, F4))
        @test field(C_changed) == F4
        @test base_ring(generator_matrix(C_changed)) == F4

        info1 = random_information_set(HammingCode(2, 3);
                                       rng=MersenneTwister(2026))
        info2 = random_information_set(HammingCode(2, 3);
                                       rng=MersenneTwister(2026))
        @test info1 == info2
        @test rank(generator_matrix(HammingCode(2, 3))[:, info1]) == 4
    end

    @testset "Nonzero Hulls" begin
        # a self-orthogonal code is its own hull
        C = RepetitionCode(2, 4)
        Hl, d = Euclidean_hull(C)
        @test d == 1
        @test (Hl.n, Hl.k) == (4, 1)
        @test Hl ⊆ C
        @test Hl ⊆ dual(C)

        # the hull of the Hamming code is the simplex code it contains
        C = HammingCode(2, 3)
        Hl, d = Euclidean_hull(C)
        @test d == 3
        @test (Hl.n, Hl.k) == (7, 3)
        @test Hl ⊆ C
        @test Hl ⊆ dual(C)
        @test Hl ⊊ C

        # the hull is trivial exactly for linear complementary dual codes
        @test !is_LCD(HammingCode(2, 3))
        Hl, d = Hermitian_hull(HammingCode(4, 2))
        @test d == 2
        @test Hl ⊆ Hermitian_dual(HammingCode(4, 2))

        # the all-ones word has coordinate sum 1 over GF(4), so it is not in the
        # 1-Galois dual and the hull is trivial
        Hl, d = l_Galois_hull(RepetitionCode(4, 3), 1)
        @test d == 0
        @test ismissing(Hl)
        @test is_l_Galois_LCD(RepetitionCode(4, 3), 1)
    end

    @testset "Characteristic Polynomial" begin
        # the product runs over the nonzero dual weights, which for the [7,4,3]
        # Hamming code is the single weight 4 of the simplex code
        p = CodingTheory.characteristic_polynomial(HammingCode(2, 3))
        _, x = polynomial_ring(Oscar.QQ, :x)
        @test p == 8 - 2 * x

        C = RepetitionCode(2, 3)
        q_pow = 2^(C.n - C.k)
        weights = sort!([w for (w, c) in weight_distribution(dual(C)) if
                         !iszero(c) && w > 0])
        @test CodingTheory.characteristic_polynomial(C) ==
            q_pow * prod(1 - x // j for j in weights)
    end
end
