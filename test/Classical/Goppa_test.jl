@testitem "Classical/Goppa.jl" begin
    using Oscar, CodingTheory

    @testset "Goppa Code" begin
        # Ling & Xing, Example 9.3.10 (iii), p. 200
        E = GF(8, :α)
        S, z = polynomial_ring(E, :z)
        α = gen(E)
        # α satisfies α^3 + α + 1
        # julia> minimal_polynomial(α)
        # x^3 + x + 1
        g = α^3 + z + z^2
        L = [E(0); [α^i for i = 0:6]]
        F = GF(2)
        C = GoppaCode(F, L, g)
        @test length(C) == 8
        @test dimension(C) == 2
        @test minimum_distance(C) == 5

        # Ling & Xing, Corollary 9.3.4, p. 198
        n = length(L)
        # v = [g(a)^(-1) for a in L]
        v = [g(L[i]) * prod(L[i] - L[j] for j = 1:n if i ≠ j)^-1 for i = 1:n]
        t = degree(g)
        C2 = AlternateCode(F, n - t, v, L)
        # TODO I highly suspect this result comes from using different bases and scalars throughout different books and the expansion is highly basis dependent
        @test_broken are_equivalent(C, C2)

        # MacWilliams & Sloane, p. 342
        g = z^2 + z + 1
        C = GoppaCode(F, L, g)
        @test is_irreducible(C)
        C_ext = extend(C)
        @test length(C_ext) == 9
        @test dimension(C_ext) == 2
        @test minimum_distance(C_ext) == 6
        C_ext_perm = permute_code(C_ext, [2, 6, 8, 9, 4, 7, 1, 3, 5])
        flag, C_cyc = is_cyclic(C_ext_perm)
        @test flag

        # TODO need permutation in M&S for cyclic test
        # MacWilliams & Sloane, p. 355
        # g = quadratic with distinct roots
        # L = set diff GF(2^m) minus those
        # extend GoppaCode
        # is cyclic
        # E = GF(2^5, :α)
        # S, z = polynomial_ring(E, :z)
        # α = gen(E)
        # g = z^2 + z + 1
        # L = [E(0); [α^i for i in 0:Int(order(E)) - 2]]
        # C = GoppaCode(F, L, g)
        # @test is_irreducible(C)
        # C_ext = extend(C)
        # @test length(C_ext) == 33
        # @test dimension(C_ext) == 22
        # @test minimum_distance(C_ext) == 6
        # C_ext_perm = permute_code(C_ext, need)
        # flag, C_cyc = is_cyclic(C_ext_perm)
        # @test flag
        # # which ring is this even in?
        # @test generator_polynomial(C_cyc) == 1
        # 1 + x^2 + x^5 + x^6 + x^9 + x^11

        # MacWilliams & Sloane, p. 343
        E = GF(2^5, :α)
        S, z = polynomial_ring(E, :z)
        α = gen(E)
        L = [E(0); [α^i for i = 0:(Int(order(E))-2)]]
        g = z^3 + z + 1
        C = GoppaCode(F, L, g)
        @test is_irreducible(C)
        @test length(C) == 32
        @test dimension(C) == 17
        @test minimum_distance(C) == 7

        # MacWilliams & Sloane, p. 344
        E = GF(2^4, :α)
        S, z = polynomial_ring(E, :z)
        α = gen(E)
        g = z^2 + z + α^3
        L = [E(0); [α^i for i = 0:(Int(order(E))-2)]]
        C = GoppaCode(GF(2), L, g)
        @test is_irreducible(C)
        @test length(C) == 16
        @test dimension(C) == 8
        @test minimum_distance(C) == 5

        # MacWilliams & Sloane, p. 349
        E = GF(2^4, :α)
        S, z = polynomial_ring(E, :z)
        α = gen(E)
        g = z^3 + z + 1
        L = [α^i for i = 0:14]
        C = GoppaCode(GF(2), L, g)
        @test length(C) == 15
        @test dimension(C) == 3
        @test minimum_distance(C) == 7

        # is equivalent to narrow-sense BCH code
        # Ling & Xing, p. 200
        # MacWilliams & Sloane, p. 345
        # Huffman & Pless, p. 522
        # the narrow-sense BCH code of length n and designed distance δ is the Goppa code with L = [1, β^(-1), β^(-2), ..., β^(1 - n)] and g = x^(δ - 1)
        # here β is the primitive root used in the BCH code
        # TODO I may not be fully understanding the variables being used here as they change from source to source
        # C = BCHCode(2, 15, 5)
        # α = primitive_root(C)
        # L = [α^(1 - i) for i in 1:C.n]
        # E = extension_field(C)
        # S, z = polynomial_ring(E, :z)
        # g = z^(design_distance(C) - 1)
        # C2 = GoppaCode(GF(2), L, g)
        # @test_broken are_equivalent(C, C2)
    end
end
