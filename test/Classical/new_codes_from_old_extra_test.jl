@testitem "New codes from old extras" begin
    using Oscar, CodingTheory, Random

    function manual_kron(A, B)
        F = base_ring(A)
        K = zero_matrix(F, nrows(A) * nrows(B), ncols(A) * ncols(B))
        for i in 1:nrows(A), j in 1:ncols(A), r in 1:nrows(B), s in 1:ncols(B)
            K[(i - 1) * nrows(B) + r, (j - 1) * ncols(B) + s] = A[i, j] * B[r, s]
        end
        return K
    end

    @testset "Direct, product, tensor, and Plotkin structures" begin
        C1 = HammingCode(2, 3)       # [7, 4, 3]
        C2 = RepetitionCode(2, 3)    # [3, 1, 3]
        F = field(C1)
        G1, G2 = generator_matrix(C1), generator_matrix(C2)
        H1, H2 = parity_check_matrix(C1), parity_check_matrix(C2)

        S = CodingTheory.direct_sum(C1, C2)
        expected_sum = vcat(
            hcat(G1, zero_matrix(F, dimension(C1), length(C2))),
            hcat(zero_matrix(F, dimension(C2), length(C1)), G2))
        @test (length(S), dimension(S)) == (10, 5)
        @test generator_matrix(S) == expected_sum

        P = CodingTheory.product_code(C1, C2)
        DP = CodingTheory.direct_product(C1, C2)
        expected_product = manual_kron(G1, G2)
        @test (length(P), dimension(P)) == (21, 4)
        @test generator_matrix(P) == expected_product
        @test generator_matrix(DP) == expected_product
        @test generator_matrix(C1 × C2) == expected_product

        K = CodingTheory.kron(C1, C2)
        expected_check = manual_kron(H1, H2)
        @test (length(K), dimension(K)) == (21, 15)
        @test parity_check_matrix(K) == expected_check
        @test iszero(generator_matrix(K) * transpose(expected_check))

        R = RepetitionCode(2, 7)
        Plot = Plotkin_construction(C1, R)
        expected_plotkin = vcat(
            hcat(G1, G1),
            hcat(zero_matrix(F, dimension(R), length(C1)), generator_matrix(R)))
        expected_plotkin_check = vcat(
            hcat(H1, zero_matrix(F, nrows(H1), ncols(H1))),
            hcat(-parity_check_matrix(R), parity_check_matrix(R)))
        @test (length(Plot), dimension(Plot)) == (14, 5)
        @test generator_matrix(Plot) == expected_plotkin
        @test parity_check_matrix(Plot) == expected_plotkin_check
        @test iszero(expected_plotkin * transpose(expected_plotkin_check))
    end

    @testset "Entrywise products" begin
        C = HammingCode(2, 3)
        R = RepetitionCode(2, 7)
        E = CodingTheory.entrywise_product_code(C, R)
        @test dimension(E) == dimension(C)
        @test E ⊆ C && C ⊆ E
        @test Schur_product_code(C, R) ⊆ E && E ⊆ Schur_product_code(C, R)
        @test Hadamard_product_code(C, R) ⊆ E && E ⊆ Hadamard_product_code(C, R)
        @test componentwise_product_code(C, R) ⊆ E &&
              E ⊆ componentwise_product_code(C, R)

        Sq = CodingTheory.entrywise_product_code(C, C)
        G = generator_matrix(C)
        for i in 1:nrows(G)
            square = [G[i, j]^2 for j in 1:ncols(G)]
            @test square ∈ Sq
        end
    end

    @testset "Residue and Y constructions" begin
        C = HammingCode(2, 3) # [7, 4, 3]
        c = vec(Array(generator_matrix(C)[1, :]))
        w = wt(c)
        A = construction_A(C, c)
        Res = residue_code(C, c)
        @test (length(A), dimension(A)) == (7 - w, 3)
        @test generator_matrix(A) == generator_matrix(Res)
        @test all(iszero, generator_matrix(A) * transpose(parity_check_matrix(A)))

        h = vec(Array(parity_check_matrix(C)[1, :]))
        s = wt(h)
        B = construction_B(C, h)
        Y = construction_Y(C, h)
        Y1 = construction_Y1(C, h)
        @test length(B) == 7 - s
        @test dimension(B) >= 4 - s + 1
        @test generator_matrix(B) == generator_matrix(Y)
        @test generator_matrix(B) == generator_matrix(Y1)

        B2 = construction_B2(C, h, 1)
        @test length(B2) == 7 - s
        @test dimension(B2) <= dimension(C)
        @test minimum_distance_lower_bound(B2) >= minimum_distance_lower_bound(C) - 2
    end

    @testset "Quotients, complements, and intermediate subcodes" begin
        C = HammingCode(2, 3)
        D = subcode(C, 2)
        Q = CodingTheory.quotient(C, D)
        Comp = CodingTheory.code_complement(D, C)
        @test (length(Q), dimension(Q)) == (7, 2)
        @test generator_matrix(Q) == generator_matrix(Comp)
        @test all(vec(Array(generator_matrix(Q)[i, :])) ∈ C for i in 1:dimension(Q))
        @test rank(vcat(generator_matrix(D), generator_matrix(Q))) == dimension(C)

        Mid = subcode_of_dimension_between_codes(C, D, 3)
        @test (length(Mid), dimension(Mid)) == (7, 3)
        @test D ⊆ Mid
        @test Mid ⊆ C
        @test subcode_of_dimension_between_codes(C, D, 2) === D
        @test subcode_of_dimension_between_codes(C, D, 4) === C
    end

    @testset "Extensions, permutations, and triply-even subcodes" begin
        C = HammingCode(2, 3)
        E = even_extension(C)
        GE = generator_matrix(E)
        @test (length(E), dimension(E)) == (8, 4)
        @test GE[:, 1:7] == generator_matrix(C)
        @test all(iszero(sum(GE[i, j] for j in 1:8)) for i in 1:nrows(GE))

        σ = reverse(collect(1:7))
        P = permute_code(C, σ)
        @test (length(P), dimension(P)) == (7, 4)
        @test generator_matrix(P) == generator_matrix(C)[:, σ]

        T = triply_even_subcode(ExtendedHammingCode(3))
        @test !ismissing(T)
        GT = generator_matrix(T)
        for mask in 0:(2^dimension(T) - 1)
            word = zero_matrix(field(T), 1, length(T))
            for i in 1:dimension(T)
                isodd(mask >> (i - 1)) && (word += GT[i:i, :])
            end
            @test wt(word) % 8 == 0
        end
    end
end
