@testset "quantumproductcodes.jl" begin
    using Oscar, CodingTheory

    # Degenerate Quantum LDPC Codes With Good Finite Length Performance
    # Example A1
    F = GF(2)
    S, x = PolynomialRing(F, :x)
    l = 127
    R = residue_ring(S, x^l - 1)
    a = 1 + x^15 + x^20 + x^28 + x^66
    b = 1 + x^58 + x^59 + x^100 + x^121
    g = gcd(a, b, x^l - 1)
    aR = R(a)
    bR = R(b)
    Q = GeneralizedBicycleCode(aR, bR)
    @test length(Q) == 254
    @test CodingTheory.dimension(Q) == 28
    @test LogicalTrait(typeof(Q)) == HasLogicals()
    @test GaugeTrait(typeof(Q)) == HasNoGauges()

    # Example A2
    l = 24
    R = residue_ring(S, x^l - 1)
    a = 1 + x^2 + x^8 + x^15
    b = 1 + x^2 + x^12 + x^17
    Q = GeneralizedBicycleCode(R(a), R(b))
    @test length(Q) == 48
    @test CodingTheory.dimension(Q) == 6
    @test LogicalTrait(typeof(Q)) == HasLogicals()
    @test GaugeTrait(typeof(Q)) == HasNoGauges()

    # Example B1
    l = 63
    R = residue_ring(S, x^l - 1)
    A = matrix(R, 7, 7,
	    [x^27, 0, 0, 1, x^18, x^27, 1,
	     1, x^27, 0, 0, 1, x^18, x^27,
	     x^27, 1, x^27, 0, 0, 1, x^18,
	     x^18, x^27, 1, x^27, 0, 0, 1,
	     1, x^18, x^27, 1, x^27, 0, 0,
	     0, 1, x^18, x^27, 1, x^27, 0,
	     0, 0, 1, x^18, x^27, 1, x^27])
    b = R(1 + x + x^6)
    Q = LiftedGeneralizedHypergraphProductCode(A, b)
    @test length(Q) == 882
    @test CodingTheory.dimension(Q) == 48
    @test LogicalTrait(typeof(Q)) == HasLogicals()
    @test GaugeTrait(typeof(Q)) == HasNoGauges()

    # HGP codes are (l, q)-QLDPC code
    # l, q = columnrowweights(Sq2)

    # single code HGP tests
    # length C.n^2 + (n^T)^2
    # dimension C.k^2 + (k^T)^2,
    # minimum distance min(d, d^T)
    # stabilizers have row weights of the form i + j, where i and j are the
    # row and column weights of the H, respecitvely

    # two code HGP tests
    # [[(C1.n)^2 + (C2^T.n)^2, (C1.k)^2 + (C2^T.k)^2, min(d, d^T)]]
    # dX = min(d^T_1, d_2), dZ = min(d1, d^T_2), d = min(dX, dZ)

    # GeneralizedShorCode
    # [[n1 * n2, k1 * k2, min(d1, d2)]]
end
