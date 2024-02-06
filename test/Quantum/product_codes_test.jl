@testset "Quantum/product_codes.jl" begin
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
    Q = LiftedProductCode(A, b)
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

    # Quintavalle_basis
    Fone = F(1)
    for r in 2:4
        C = HammingCode(2, r)
        F = C.F
        Fone = F(1)
        H = parity_check_matrix(C)
        nr_H = nrows(H)

        # Add some degerate rows so the dual code has k > 0
        H = vcat(H, H[rand(1:nr_H), :] + H[rand(1:nr_H), :],
            H[rand(1:nr_H), :] + H[rand(1:nr_H), :])
        HGP = HypergraphProductCode(LinearCode(H, true))

        lx, lz = Quintavalle_basis(HGP)
        # Check the logical operators commute with the stabilizers
        @test iszero(HGP.X_stabs * transpose(lz))
        @test iszero(HGP.Z_stabs * transpose(lx))

        # Overlap should be zero in all cases except i == ii
        zero_sum_flag = true
        one_sum_flag = true
        count_zero_flag = true
        count_one_flag = true
        for i in 1:nrows(lx)
            for ii in 1:nrows(lz)
                if i != ii
                    iszero(sum(lx[i, :] .* lz[ii, :])) || (zero_sum_flag = false)
                    iszero(length(findall(x -> x == Fone, lx[i, :] .* lz[ii, :]))) || (count_zero_flag = false;)
                else
                    isone(sum(lx[i, :] .* lz[ii, :])) || (one_sum_flag = false;)
                    isone(length(findall(x -> x == Fone, lx[i, :] .* lz[ii, :]))) || (count_one_flag = false;)
                end
            end
        end
        @test zero_sum_flag
        @test one_sum_flag
        @test count_zero_flag
        @test count_one_flag

        # Check the logical operators have weight >= code distance
        weight_flag = true
        for i in 1:HGP.k
            (wt(lx[i, :]) < HGP.d || wt(lz[i, :]) < HGP.d) && (weight_flag = false;)
        end
        @test weight_flag
    end

    H1 = parity_check_matrix(HammingCode(2, 2))
    H2 = parity_check_matrix(HammingCode(2, 3))
    H1 = vcat(H1, H1[1, :] + H1[2, :])
    H2 = vcat(H2, H2[1, :] + H2[2, :])
    HGP = HypergraphProductCode(LinearCode(H1, true), LinearCode(H2, true))

    lx, lz = Quintavalle_basis(HGP)
    # Check the logical operators commute with the stabilizers
    @test iszero(HGP.X_stabs * transpose(lz))
    @test iszero(HGP.Z_stabs * transpose(lx))

    # Overlap should be zero in all cases except i == ii
    zero_sum_flag = true
    one_sum_flag = true
    count_zero_flag = true
    count_one_flag = true
    for i in 1:nrows(lx)
        for ii in 1:nrows(lz)
            if i != ii
                iszero(sum(lx[i, :] .* lz[ii, :])) || (zero_sum_flag = false)
                iszero(length(findall(x -> x == Fone, lx[i, :] .* lz[ii, :]))) || (count_zero_flag = false)
            else
                isone(sum(lx[i, :] .* lz[ii, :])) || (one_sum_flag = false)
                isone(length(findall(x -> x == Fone, lx[i, :] .* lz[ii, :]))) || (count_one_flag = false)
            end
        end
    end
    @test zero_sum_flag
    @test one_sum_flag
    @test count_zero_flag
    @test count_one_flag

    # Check the logical operators have weight >= code distance
    weight_flag = true
    for i in 1:HGP.k
        (wt(lx[i, :]) < HGP.d || wt(lz[i, :]) < HGP.d) && (weight_flag = false)
    end
    @test weight_flag

    # product codes
    F = GF(2);
    h = matrix(F, [1 1]);
    id = identity_matrix(F, 2);
    H_X = vcat(
        h ⊗ h ⊗ h ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id,
        id ⊗ id ⊗ id ⊗ h ⊗ h ⊗ h ⊗ id ⊗ id ⊗ id,
        id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ h ⊗ h ⊗ h);
    H_Z = vcat(
        h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id,
        id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id,
        id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h);
    SPCtest = CodingTheory.SPCDFoldProductCode(3);
    @test length(SPCtest) == 512
    @test dimension(SPCtest) == 174
    @test H_X = SPCtest.X_stabs
    @test H_Z = SPCtest.Z_stabs

end
