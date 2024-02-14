@testset "Classical/cyclic_code.jl" begin
    # using Oscar, CodingTheory

    # examples: Huffman & Pless
    cosets = defining_set([1, 2, 3, 4, 5, 6], 2, 7, false)
    C = CyclicCode(2, 7, cosets)
    R = polynomial_ring(C)
    x = gen(R)
    @test CodingTheory.dimension(C) == 1
    @test generator_polynomial(C) == 1 + x + x^2 + x^3 + x^4 + x^5 + x^6
    @test idempotent(C) == 1 + x + x^2 + x^3 + x^4 + x^5 + x^6
    cosets = defining_set([0, 1, 2, 4], 2, 7, false)
    C = CyclicCode(2, 7, cosets)
    @test CodingTheory.dimension(C) == 3
    @test generator_polynomial(C) == 1 + x^2 + x^3 + x^4
    @test idempotent(C) == 1 + x^3 + x^5 + x^6
    cosets = defining_set([0, 3, 5, 6], 2, 7, false)
    C = CyclicCode(2, 7, cosets)
    @test CodingTheory.dimension(C) == 3
    @test generator_polynomial(C) == 1 + x + x^2 + x^4
    @test idempotent(C) == 1 + x + x^2 + x^4
    cosets = defining_set([1, 2, 4], 2, 7, false)
    C = CyclicCode(2, 7, cosets)
    @test CodingTheory.dimension(C) == 4
    @test generator_polynomial(C) == 1 + x + x^3
    @test idempotent(C) == x + x^2 + x^4
    cosets = defining_set([3, 5, 6], 2, 7, false)
    C = CyclicCode(2, 7, cosets)
    @test CodingTheory.dimension(C) == 4
    @test generator_polynomial(C) == 1 + x^2 + x^3
    @test idempotent(C) == x^3 + x^5 + x^6

    # the dual of a cyclic code is the complement code with multiplier -1 (p.146)
    # do Theorem 4.4.11 on page 147

    # a self-orthogonal binary cyclic code is doubly-even

    # example: Huffman & Pless
    C = BCHCode(3, 13, 2, 1)
    @test defining_set(C) == [1, 3, 9]
    R = polynomial_ring(C)
    x = gen(R)
    @test generator_polynomial(C) == 2 + x + x^2 + x^3
    @test CodingTheory.dimension(C) == 10
    @test minimum_distance(C) == 3
    C = BCHCode(3, 13, 3, 1)
    @test defining_set(C) == [1, 2, 3, 5, 6, 9]
    @test generator_polynomial(C) == 1 + 2*x + x^2 + 2*x^3 + 2*x^4 + 2*x^5 + x^6
    @test CodingTheory.dimension(C) == 7
    @test minimum_distance(C) == 4
    C = BCHCode(3, 13, 5, 1)
    @test defining_set(C) == [1, 2, 3, 4, 5, 6, 9, 10, 12]
    @test generator_polynomial(C) == 2 + 2*x^2 + 2*x^3 + x^5 + 2*x^7 + x^8 + x^9
    @test CodingTheory.dimension(C) == 4
    @test minimum_distance(C) == 7
    @test CodingTheory.dimension(C) >= length(C) - ord(length(C), 3)*(5 - 1)

    R, (x, y) = PolynomialRing(Nemo.ZZ, (:x, :y))

    # example: MacWilliams & Sloane
    # any cyclic code over GF(2^m) of length 2^m + 1 is reversible

    # example: MacWilliams & Sloane
    C = BCHCode(2, 31, 5, 1)
    @test CodingTheory.dimension(C) == 21
    @test minimum_distance(C) == 5
    @test polynomial(MacWilliams_identity(C, weight_enumerator(C, :Hamming))) == y^31 + 310*x^12*y^19 + 527*x^16*y^15 + 186*x^20*y^11

    # example: Huffman & Pless
    C = ReedSolomonCode(13, 5, 1)
    @test length(C) == 12
    @test CodingTheory.dimension(C) == 8
    @test minimum_distance(C) == 5
    # @test isMDS(C) == true
    @test defining_set(C) == [1, 2, 3, 4]
    R = polynomial_ring(C)
    x = gen(R)
    @test generator_polynomial(C) == 10 + 2*x + 7*x^2 + 9*x^3 + x^4
    D = dual(C)
    @test CodingTheory.dimension(D) == 4
    @test minimum_distance(D) == 9
    # @test isMDS(D) == true
    @test defining_set(D) == [0, 1, 2, 3, 4, 5, 6, 7]
    @test generator_polynomial(D) == 3 + 12*x + x^2 + 5*x^3 + 11*x^4 + 4*x^5 + 10*x^6 + 5*x^7 + x^8
    Cc = complement(C)
    @test length(Cc) == 12
    @test CodingTheory.dimension(Cc) == 4
    @test minimum_distance(Cc) == 9
    @test defining_set(Cc) == [0, 5, 6, 7, 8, 9, 10, 11]
    @test generator_polynomial(Cc) == 9 + 6*x + 12*x^2 + 10*x^3 + 8*x^4 + 6*x^5 + 9*x^6 + 4*x^7 + x^8

    # example: Huffman & Pless
    C = ReedSolomonCode(16, 7, 1)
    @test length(C) == 15
    @test CodingTheory.dimension(C) == 9
    @test minimum_distance(C) == 7
    @test defining_set(C) == [1, 2, 3, 4, 5, 6]
    R = polynomial_ring(C)
    x = gen(R)
    α = primitive_root(C)
    @test generator_polynomial(C) == α^6 + α^9*x + α^6*x^2 + α^4*x^3 + α^14*x^4 + α^10*x^5 + x^6

    # example: MacWilliams & Sloane
    C = ReedSolomonCode(5, 3, 1)
    z = gen(polynomial_ring(C))
    @test generator_polynomial(C) == z^2 + 4*z + 3

    # example: MacWilliams & Sloane
    C = ReedSolomonCode(8, 6)
    @test CodingTheory.dimension(C) == 2
    z = gen(polynomial_ring(C))
    α = primitive_root(C)
    @test idempotent(C) == α^4*z + α*z^2 + α^4*z^3 + α^2*z^4 + α^2*z^5 + α*z^6

    # example: MacWilliams & Sloane
    C = ReedSolomonCode(8, 3, 5)
    @test CodingTheory.dimension(C) == 5
    z = gen(polynomial_ring(C))
    α = primitive_root(C)
    @test generator_polynomial(C) == α^4 + α*z + z^2
    # expand this code over F_2, is equivalent to the following BCH code
    # C2 = BCHCode(2, 21, 3, 1) # maybe not b = 1?
    # z2 = gen(polynomial_ring(C2))
    # @test generator_polynomial(C2) == 1 + z2 + z2^2 + z2^4 + z2^6
    # @test are_equivalent(expC, C2)

    # # example: MacWilliams & Sloane
    # # extended Reed-Solomon codes have distance d + 1
    # # TODO: fix extend here
    # # extC = extend(C)
    # # @test minimum_distance(extC) == minimum_distance(C) + 1

    # example: MacWilliams & Sloane
    # some [15, 6, 6] binary BCH code with min polys (-1, 0, 1) is reversible
    # C = BCHCode(2, 15, 6, 13)
    # println(C)
    # @test isreversible(C) == true

    # RS codes contain BCH codes
    C = ReedSolomonCode(16, 5)
    C2 = BCHCode(2, 15, 5)
    @test C2 ⊆ C
    @test C2 ⊂ C
    @test is_subcode(C2, C)
    @test C == CyclicCode(16, 15, defining_set([i for i = 0:(0 + 5 - 2)], 16, 15, false))
    @test C == BCHCode(16, 15, 5)
    @test design_distance(C) == 5
    @test is_narrowsense(C)
    @test is_primitive(C)

    # is_cyclic - true parameter also tests cyclic code constructor given generator polynomial
    # C = ReedSolomonCode(7, 3)
    # H = HammingCode(2, 3)
    # @test is_cyclic(H, false) == false
    # _, C2 = is_cyclic(C, true) # this true is construct, can do an are_equivalent here

end
