@testitem "Utils extra helpers" begin
    using Oscar, CodingTheory, Graphs

    @testset "Weights, distances, digits, and binomials" begin
        v = [0, -2, 0, 5, 7, 0]
        w = [1, -2, 3, 0, 7, 0]
        @test Hamming_weight(v) == 3
        @test CodingTheory.weight(v) == 3
        @test dist(v, w) == 3
        @test distance(v, w) == 3

        F3 = GF(3)
        vf = matrix(F3, [0 1 2 0])
        wf = matrix(F3, [1 1 0 0])
        @test Hamming_weight(vf) == 2
        @test CodingTheory.weight(vf) == 2
        @test dist(vf, wf) == 2
        @test distance(vf, wf) == 2

        @test digits_to_int([1, 0, 1, 1], 2) == 11
        @test digits_to_int([2, 1, 0, 2], 3) == 65
        for (n, base) in ((37, 2), (314, 5), (2026, 10))
            @test digits_to_int(reverse(digits(n; base=base)), base) == n
        end

        @test extended_binomial(12, 5) == UInt128(binomial(12, 5))
        @test extended_binomial(5, 8) == UInt128(0)
    end

    @testset "Finite-field bases, traces, and expansion" begin
        F2 = GF(2)
        E = GF(2, 3, :α)
        α = gen(E)
        basis = [one(E), α, α^2]
        complementary = complementary_basis(E, F2, basis)

        @test verify_complementary_basis(E, F2, basis, complementary)
        @test all(CodingTheory.tr(basis[i] * complementary[j], F2) ==
                  (i == j ? one(E) : zero(E))
                  for i in eachindex(basis), j in eachindex(complementary))

        normal, normal_dual = normal_basis(E, F2)
        @test normal == normal_dual
        @test is_self_dual_basis(E, F2, normal)

        @test CodingTheory.tr(one(E), F2) == one(E)
        @test CodingTheory.tr(α, F2) == α + α^2 + α^4

        M = matrix(E, [one(E) α α^2; α^2 + one(E) zero(E) α + one(E)])
        expanded = expand_matrix(M, F2, basis)
        @test size(expanded) == (2, 9)
        @test expanded == matrix(F2, [
            1 0 0  0 1 0  0 0 1;
            1 0 1  0 0 0  1 1 0
        ])
    end

    @testset "Supports and triangular reduction" begin
        F2 = GF(2)
        M = matrix(F2, [1 0 1 0; 0 1 1 0; 0 0 0 1])
        @test row_supports(M) == [[1, 3], [2, 3], [4]]

        S = matrix(F2, [1 0 0 1 0 1; 0 1 0 0 1 0])
        @test row_supports_symplectic(S) ==
              [([1], [1, 3]), ([2], [2])]
        @test_throws ArgumentError row_supports_symplectic(M[:, 1:3])

        A = matrix(F2, [1 1 0; 0 1 1])
        ker_basis, image_complement = strongly_lower_triangular_reduction(A)
        @test ker_basis == matrix(F2, 3, 1, [1, 1, 1])
        @test image_complement == matrix(F2, 3, 1, [0, 0, 1])
        @test iszero(A * ker_basis)
        @test rank(hcat(transpose(A), image_complement)) == 3
    end

    @testset "Graph and group-algebra helpers" begin
        path = Graphs.path_graph(3)
        expected = [
            0 0 1 1 0;
            0 0 0 1 1;
            1 0 0 0 0;
            1 1 0 0 0;
            0 1 0 0 0
        ]
        incidence, left, right = edge_vertex_incidence_matrix(path)
        @test incidence == expected
        @test left == [1, 2]
        @test right == [3, 4, 5]

        incidence_graph, graph_left, graph_right =
            edge_vertex_incidence_graph(path)
        @test Matrix(Graphs.adjacency_matrix(incidence_graph)) == expected
        @test (graph_left, graph_right) == (left, right)
        @test is_valid_bipartition(incidence_graph, left, right)
        @test !is_valid_bipartition(incidence_graph, [1, 3], [2, 4, 5])

        F2 = Oscar.Nemo.Native.GF(2)
        group = abelian_group([3])
        algebra = group_algebra(F2, group)
        x = algebra(1) + algebra(gens(group)[1])
        expected_col = matrix(F2, [1 0 1; 1 1 0; 0 1 1])
        expected_row = matrix(F2, [1 1 0; 0 1 1; 1 0 1])
        @test group_algebra_element_to_circulant_matrix(x) == expected_col
        @test group_algebra_element_to_circulant_matrix(x, :row) == expected_row
    end

    @testset "alist loading and printing" begin
        alist = """
        2 3
        2 2
        2 2
        1 1 2
        1 3
        2 3
        1
        2
        1 2
        """
        mktemp() do path, io
            write(io, alist)
            close(io)
            @test load_alist(path) == [1 0 1; 0 1 1]
        end

        mktemp() do path, io
            redirect_stdout(io) do
                @test isnothing(print_string_array(["IX", "YZ"], true))
            end
            close(io)
            @test read(path, String) == " X\nYZ\n"
        end

        mktemp() do path, io
            redirect_stdout(io) do
                print_string_array(["IX", "YZ"])
            end
            close(io)
            @test read(path, String) == "IX\nYZ\n"
        end
    end

    @testset "extended_binomial" begin
        @test CodingTheory.extended_binomial(5, 2) == 10
        @test CodingTheory.extended_binomial(2, 5) == 0
        @test CodingTheory.extended_binomial(10, 0) == 1

        # the intermediate overflows Int64 while the result still fits UInt128
        @test CodingTheory.extended_binomial(100, 50) ==
            UInt128(binomial(big(100), big(50)))
        @test CodingTheory.extended_binomial(120, 60) ==
            UInt128(binomial(big(120), big(60)))
    end

    @testset "extract_bipartition" begin
        Grphs = CodingTheory.Graphs

        # the returned vectors hold vertex indices, one side each
        left, right = CodingTheory.extract_bipartition(Grphs.path_graph(3))
        @test sort(left) == [1, 3]
        @test sort(right) == [2]

        left, right = CodingTheory.extract_bipartition(
            Grphs.complete_bipartite_graph(2, 3))
        @test sort(left) == [1, 2]
        @test sort(right) == [3, 4, 5]

        # every vertex lands on exactly one side, and no edge stays within a side
        G = Grphs.complete_bipartite_graph(3, 4)
        left, right = CodingTheory.extract_bipartition(G)
        @test sort(vcat(left, right)) == collect(1:Grphs.nv(G))
        @test isempty(intersect(left, right))
        for e in Grphs.edges(G)
            @test (Grphs.src(e) in left) != (Grphs.dst(e) in left)
        end

        @test_throws ArgumentError CodingTheory.extract_bipartition(
            Grphs.cycle_graph(3))
    end
end
