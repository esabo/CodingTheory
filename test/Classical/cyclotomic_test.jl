@testitem "Classical/cyclotomic.jl" begin
    using CodingTheory, Oscar

    @testset "Multiplicative Order and Basic Cosets" begin
        # ord(n, q)
        @test CodingTheory.ord(7, 2) == 3
        @test CodingTheory.ord(15, 2) == 4
        @test_throws ArgumentError CodingTheory.ord(6, 2) # Not coprime

        # Single coset generation
        @test CodingTheory.cyclotomic_coset(1, 2, 7) == [1, 2, 4]
        @test CodingTheory.cyclotomic_coset(3, 2, 7) == [3, 5, 6]

        # ALL YOUR ORIGINAL TEXTBOOK EXAMPLES ARE PRESERVED HERE
        @test CodingTheory.all_cyclotomic_cosets(2, 15, to_sort = false) == 
            [[0], [1, 2, 4, 8], [3, 6, 12, 9], [5, 10], [7, 14, 13, 11]]
        @test CodingTheory.all_cyclotomic_cosets(3, 13, to_sort = true) == 
            [[0], [1, 3, 9], [2, 5, 6], [4, 10, 12] , [7, 8, 11]]
            
        @test_throws DomainError CodingTheory.all_cyclotomic_cosets(4, 6)
    end

    @testset "Complement, Dual, and Pairings" begin
        # For q=2, n=7: Complete set is [[0], [1, 2, 4], [3, 5, 6]]
        
        # Complement of just the zero coset
        @test CodingTheory.complement_qcosets(2, 7, [[0]]) == [[1, 2, 4], [3, 5, 6]]
        
        # Dual qcosets = Negation mod n of the complement
        # Comp of [[1, 2, 4]] is [[0], [3, 5, 6]].
        # Negation of [3, 5, 6] mod 7 is [4, 2, 1], sorted -> [1, 2, 4].
        # So the dual of [[1, 2, 4]] is [[0], [1, 2, 4]]
        @test CodingTheory.dual_qcosets(2, 7, [[1, 2, 4]]) == [[0], [1, 2, 4]]

        # Pairings match a coset with its negation
        pairs, reps = CodingTheory.qcoset_pairings(2, 7)
        @test ([0], [0]) in pairs
        
        # Since qcoset_pairings(2, 7) defaults to to_sort=false internally,
        # the coset for 3 is generated sequentially as [3, 6, 5]
        @test ([1, 2, 4], [3, 6, 5]) in pairs || ([3, 6, 5], [1, 2, 4]) in pairs
    end

    @testset "Field Conjugates and Minimal Polynomials" begin
        # Setup GF(2^3)
        E, a = finite_field(2, 3, "a")
        
        # x and y are conjugates if y = x^(q^i)
        @test CodingTheory.are_conjugates(a, a^2, 2)
        @test CodingTheory.are_conjugates(a, a^4, 2)
        @test !CodingTheory.are_conjugates(a, a^3, 2) # a^3 is in the 3-coset
        
        # Minimal Polynomial Generation
        C1 = CodingTheory.cyclotomic_coset(1, 2, 7)
        M1 = CodingTheory.minimal_polynomial(C1, a)
        
        # The polynomial degree should match the size of the cyclotomic coset
        @test degree(M1) == 3
        
        # Mathematically, the coefficients must fall entirely within the base field GF(2).
        # Which means they evaluate to exactly E(0) or E(1).
        for c in coefficients(M1)
            @test c == E(1) || c == E(0)
        end
    end
end
