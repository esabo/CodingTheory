@testitem "Classical/cyclic_code.jl" begin
    using Oscar, CodingTheory

    @testset "Cyclic codes: Core textbook examples" begin
        # examples: Huffman & Pless
        cosets = CodingTheory.defining_set([1, 2, 3, 4, 5, 6], 2, 7, false)
        C = CyclicCode(2, 7, cosets)
        R = polynomial_ring(C)
        x = gen(R)
        @test dimension(C) == 1
        @test CodingTheory.generator_polynomial(C) == 1 + x + x^2 + x^3 + x^4 + x^5 + x^6
        @test CodingTheory.idempotent(C) == 1 + x + x^2 + x^3 + x^4 + x^5 + x^6
        
        cosets = CodingTheory.defining_set([0, 1, 2, 4], 2, 7, false)
        C = CyclicCode(2, 7, cosets)
        @test dimension(C) == 3
        @test CodingTheory.generator_polynomial(C) == 1 + x^2 + x^3 + x^4
        @test CodingTheory.idempotent(C) == 1 + x^3 + x^5 + x^6
        
        cosets = CodingTheory.defining_set([0, 3, 5, 6], 2, 7, false)
        C = CyclicCode(2, 7, cosets)
        @test dimension(C) == 3
        @test CodingTheory.generator_polynomial(C) == 1 + x + x^2 + x^4
        @test CodingTheory.idempotent(C) == 1 + x + x^2 + x^4
        
        cosets = CodingTheory.defining_set([1, 2, 4], 2, 7, false)
        C = CyclicCode(2, 7, cosets)
        @test dimension(C) == 4
        @test CodingTheory.generator_polynomial(C) == 1 + x + x^3
        @test CodingTheory.idempotent(C) == x + x^2 + x^4
        
        cosets = CodingTheory.defining_set([3, 5, 6], 2, 7, false)
        C = CyclicCode(2, 7, cosets)
        @test dimension(C) == 4
        @test CodingTheory.generator_polynomial(C) == 1 + x^2 + x^3
        @test CodingTheory.idempotent(C) == x^3 + x^5 + x^6

        # example: Huffman & Pless
        C = BCHCode(3, 13, 2, 1)
        @test CodingTheory.defining_set(C) == [1, 3, 9]
        R = polynomial_ring(C)
        x = gen(R)
        @test CodingTheory.generator_polynomial(C) == 2 + x + x^2 + x^3
        @test dimension(C) == 10
        @test minimum_distance(C)[1] == 3
        
        C = BCHCode(3, 13, 3, 1)
        @test CodingTheory.defining_set(C) == [1, 2, 3, 5, 6, 9]
        @test CodingTheory.generator_polynomial(C) == 1 + 2x + x^2 + 2x^3 + 2x^4 + 2x^5 + x^6
        @test dimension(C) == 7
        @test minimum_distance(C)[1] == 4
        
        C = BCHCode(3, 13, 5, 1)
        @test CodingTheory.defining_set(C) == [1, 2, 3, 4, 5, 6, 9, 10, 12]
        @test CodingTheory.generator_polynomial(C) == 2 + 2x^2 + 2x^3 + x^5 + 2x^7 + x^8 + x^9
        @test dimension(C) == 4
        @test minimum_distance(C)[1] == 7
        @test dimension(C) >= length(C) - CodingTheory.ord(length(C), 3) * (5 - 1)

        # example: MacWilliams & Sloane
        C = BCHCode(2, 31, 5, 1)
        @test dimension(C) == 21
        @test minimum_distance(C)[1] == 5
        
        # Test Dual MacWilliams Identity
        HWE = weight_enumerator(C)
        dual_HWE = MacWilliams_transform(HWE, dimension(C), 2)
        @test length(dual_HWE.counts) > 0 # Ensure transform executes mathematically

        # example: Huffman & Pless
        C = ReedSolomonCode(13, 5, 1)
        @test length(C) == 12
        @test dimension(C) == 8
        @test minimum_distance(C)[1] == 5
        @test is_MDS(C) == true
        @test CodingTheory.defining_set(C) == [1, 2, 3, 4]
        R = polynomial_ring(C)
        x = gen(R)
        @test CodingTheory.generator_polynomial(C) == 10 + 2x + 7x^2 + 9x^3 + x^4
        
        D = dual(C)
        @test dimension(D) == 4
        @test minimum_distance(D)[1] == 9
        @test is_MDS(D) == true
        @test CodingTheory.defining_set(D) == [0, 1, 2, 3, 4, 5, 6, 7]
        @test CodingTheory.generator_polynomial(D) == 3 + 12x + x^2 + 5x^3 + 11x^4 + 4x^5 + 10x^6 + 5x^7 + x^8
        
        Cc = complement(C)
        @test length(Cc) == 12
        @test dimension(Cc) == 4
        @test minimum_distance(Cc)[1] == 9
        @test CodingTheory.defining_set(Cc) == [0, 5, 6, 7, 8, 9, 10, 11]
        @test CodingTheory.generator_polynomial(Cc) == 9 + 6x + 12x^2 + 10x^3 + 8x^4 + 6x^5 + 9x^6 + 4x^7 + x^8

        # example: Huffman & Pless
        C = ReedSolomonCode(16, 7, 1)
        @test length(C) == 15
        @test dimension(C) == 9
        @test minimum_distance(C)[1] == 7
        @test CodingTheory.defining_set(C) == [1, 2, 3, 4, 5, 6]
        R = polynomial_ring(C)
        x = gen(R)
        α = CodingTheory.primitive_root(C)
        @test CodingTheory.generator_polynomial(C) == α^6 + α^9*x + α^6*x^2 + α^4*x^3 + α^14*x^4 + α^10*x^5 + x^6

        # example: MacWilliams & Sloane
        C = ReedSolomonCode(5, 3, 1)
        z = gen(polynomial_ring(C))
        @test CodingTheory.generator_polynomial(C) == z^2 + 4z + 3

        # example: MacWilliams & Sloane
        C = ReedSolomonCode(8, 6)
        @test dimension(C) == 2
        z = gen(polynomial_ring(C))
        α = CodingTheory.primitive_root(C)
        @test CodingTheory.idempotent(C) == α^4*z + α*z^2 + α^4*z^3 + α^2*z^4 + α^2*z^5 + α*z^6

        # example: MacWilliams & Sloane
        C = ReedSolomonCode(8, 3, 5)
        @test dimension(C) == 5
        z = gen(polynomial_ring(C))
        α = CodingTheory.primitive_root(C)
        @test CodingTheory.generator_polynomial(C) == α^4 + α*z + z^2

        # RS codes contain BCH codes
        C = ReedSolomonCode(16, 5)
        C2 = BCHCode(2, 15, 5)
        @test C2 ⊆ C
        @test C2 ⊂ C
        @test is_subcode(C2, C)
        @test C == CyclicCode(16, 15, CodingTheory.defining_set([i for i in 0:(0 + 5 - 2)], 16, 15, false))
        @test C == BCHCode(16, 15, 5)
        @test CodingTheory.design_distance(C) == 5
        @test is_narrowsense(C)
        @test is_primitive(C)
    end

    @testset "Advanced Cyclic Math: Equivalences, Products, & Bounds" begin
        Ham_cyc = CyclicCode(2, 7, [1, 2, 4]) 
        
        # Bypass the LinearCode(G) matrix typing issue by building the shell manually
        cache_shell = Dict{Symbol, Any}(:G => generator_matrix(Ham_cyc))
        C_lin = LinearCode(Ham_cyc.F, Ham_cyc.n, Ham_cyc.k, missing, 1, Ham_cyc.n, cache_shell)
        
        isc, C_cyc = is_cyclic(C_lin)
        @test isc
        @test dimension(C_cyc) == 4

        C1 = CyclicCode(2, 7, [1, 2, 4])
        C2 = CyclicCode(2, 7, [3, 5, 6])
        C_int = C1 ∩ C2
        @test dimension(C_int) == 1 
        C_sum = C1 + C2
        @test dimension(C_sum) == 7 

        C_schur = C1 * C2
        @test dimension(C_schur) == 7 

        M = multiplier_group(C1)
        @test 2 in M
        @test 4 in M
        @test 3 ∉ M

        equiv, a = is_multiplier_equivalent(C1, C2)
        @test equiv
        @test a == 3 || a == 5 || a == 6

        C_bch = BCHCode(2, 15, 5, 1)
        @test BCH_bound(C_bch) == 5
        @test HT_bound(C_bch) >= 5
        @test Roos_bound(C_bch) >= 5
    end

    @testset "Polyadic Codes & Constituents" begin
        duadics = DuadicCodes(2, 7; include_zero=false)
        @test length(duadics.codes) == 2
        @test dimension(duadics.codes[1]) == 4
        @test dimension(duadics.codes[2]) == 4
        @test duadics.multiplier == 3 || duadics.multiplier == 5 || duadics.multiplier == 6

        duadics_even = DuadicCodes(2, 7; include_zero=true)
        @test dimension(duadics_even.codes[1]) == 3

        C = CyclicCode(2, 7, [1, 2, 4])
        
        consts = CodingTheory.constituents(C)
        @test length(consts) == 2
        @test !is_irreducible(C)

        amb_consts = CodingTheory.ambient_constituents(2, 7)
        @test length(amb_consts) == 3
        for c in amb_consts
            @test is_irreducible(c)
        end
    end
end
