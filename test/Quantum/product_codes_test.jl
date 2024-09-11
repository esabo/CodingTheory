@testitem "Quantum/product_codes.jl" begin
    using Oscar
    using CodingTheory

    @testset "GeneralizedBicycleCode" begin
        # Degenerate Quantum LDPC Codes With Good Finite Length Performance
        # Example A1
        F = Oscar.Nemo.Native.GF(2)
        S, x = polynomial_ring(F, :x)
        l = 127
        R, _ = residue_ring(S, x^l - 1)
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

        # # Example A2
        l = 24
        R, _ = residue_ring(S, x^l - 1)
        a = 1 + x^2 + x^8 + x^15
        b = 1 + x^2 + x^12 + x^17
        Q = GeneralizedBicycleCode(R(a), R(b))
        @test length(Q) == 48
        @test CodingTheory.dimension(Q) == 6
        @test LogicalTrait(typeof(Q)) == HasLogicals()
        @test GaugeTrait(typeof(Q)) == HasNoGauges()
    end

    @testset "LiftedProductCode" begin
        # Degenerate Quantum LDPC Codes With Good Finite Length Performance
        # Example B1
        F = Oscar.Nemo.Native.GF(2)
        S, x = polynomial_ring(F, :x)
        l = 63
        R, _ = residue_ring(S, x^l - 1)
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
    end

    # @testset "HypergraphProductCode" begin
        # # HGP codes are (l, q)-QLDPC code
        # l, q = columnrowweights(Sq2)

        # # single code HGP tests
        # # length C.n^2 + (n^T)^2
        # # dimension C.k^2 + (k^T)^2,
        # # minimum distance min(d, d^T)
        # # stabilizers have row weights of the form i + j, where i and j are the
        # # row and column weights of the H, respecitvely

        # # two code HGP tests
        # # [[(C1.n)^2 + (C2^T.n)^2, (C1.k)^2 + (C2^T.k)^2, min(d, d^T)]]
        # # dX = min(d^T_1, d_2), dZ = min(d1, d^T_2), d = min(dX, dZ)

        # # GeneralizedShorCode
        # # [[n1 * n2, k1 * k2, min(d1, d2)]]

        # # Quintavalle_basis
        # Fone = F(1)
        # for r in 2:4
        #     C = HammingCode(2, r)
        #     F = C.F
        #     Fone = F(1)
        #     H = parity_check_matrix(C)
        #     nr_H = nrows(H)

        #     # Add some degerate rows so the dual code has k > 0
        #     rand1 = rand(1:nr_H)
        #     rand2 = rand(1:nr_H)
        #     rand3 = rand(1:nr_H)
        #     rand4 = rand(1:nr_H)
        #     H = vcat(H, H[rand1:rand1, :] + H[rand2:rand2, :],
        #              H[rand3:rand3, :] + H[rand4:rand4, :])
        #     HGP = HypergraphProductCode(LinearCode(H, true))

        #     lx, lz = Quintavalle_basis(HGP)
        #     # Check the logical operators commute with the stabilizers
        #     @test iszero(HGP.X_stabs * transpose(lz))
        #     @test iszero(HGP.Z_stabs * transpose(lx))

        #     # Overlap should be zero in all cases except i == ii
        #     zero_sum_flag = true
        #     one_sum_flag = true
        #     count_zero_flag = true
        #     count_one_flag = true
        #     for i in 1:nrows(lx)
        #         for ii in 1:nrows(lz)
        #             if i != ii
        #                 iszero(sum(lx[i:i, :] .* lz[ii:ii, :])) || (zero_sum_flag = false)
        #                 iszero(length(findall(x -> x == Fone, lx[i:i, :] .* lz[ii:ii, :]))) || (count_zero_flag = false;)
        #             else
        #                 isone(sum(lx[i:i, :] .* lz[ii:ii, :])) || (one_sum_flag = false;)
        #                 isone(length(findall(x -> x == Fone, lx[i:i, :] .* lz[ii:ii, :]))) || (count_one_flag = false;)
        #             end
        #         end
        #     end
        #     @test zero_sum_flag
        #     @test one_sum_flag
        #     @test count_zero_flag
        #     @test count_one_flag

        #     # Check the logical operators have weight >= code distance
              # weight_flag = true
              # for i in 1:HGP.k
              #    (wt(lx[i:i, :]) < HGP.d || wt(lz[i:i, :]) < HGP.d) && (weight_flag = false;)
              # end
              # @test weight_flag
        # end

        # H1 = parity_check_matrix(HammingCode(2, 2))
        # H2 = parity_check_matrix(HammingCode(2, 3))
        # H1 = vcat(H1, H1[1:1, :] + H1[2:2, :])
        # H2 = vcat(H2, H2[1:1, :] + H2[2:2, :])
        # HGP = HypergraphProductCode(LinearCode(H1, true), LinearCode(H2, true))

        # lx, lz = Quintavalle_basis(HGP)
        # Check the logical operators commute with the stabilizers
        # @test iszero(HGP.X_stabs * transpose(lz))
        # @test iszero(HGP.Z_stabs * transpose(lx))

        # Overlap should be zero in all cases except i == ii
        # zero_sum_flag = true
        # one_sum_flag = true
        # count_zero_flag = true
        # count_one_flag = true
        # for i in 1:nrows(lx)
        #     for ii in 1:nrows(lz)
        #         if i != ii
        #             iszero(sum(lx[i, :] .* lz[ii, :])) || (zero_sum_flag = false)
        #             iszero(length(findall(x -> x == Fone, lx[i, :] .* lz[ii, :]))) || (count_zero_flag = false)
        #         else
        #             isone(sum(lx[i, :] .* lz[ii, :])) || (one_sum_flag = false)
        #             isone(length(findall(x -> x == Fone, lx[i, :] .* lz[ii, :]))) || (count_one_flag = false)
        #         end
        #     end
        # end
        # @test zero_sum_flag
        # @test one_sum_flag
        # @test count_zero_flag
        # @test count_one_flag

        # # Check the logical operators have weight >= code distance
        # weight_flag = true
        # for i in 1:HGP.k
        #     (wt(lx[i, :]) < HGP.d || wt(lz[i, :]) < HGP.d) && (weight_flag = false)
        # end
        # @test weight_flag

        # # [[400,16,6]] code from Table 1 of https://doi.org/10.1103/PhysRevResearch.2.043423
        # H = matrix(
        #         Oscar.Nemo.Native.GF(2),
        #         [
        #                 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0;
        #                 0 0 1 0 0 0 1 1 0 0 0 0 1 0 0 0;
        #                 0 0 0 1 1 0 1 0 0 0 0 0 0 1 0 0;
        #                 0 0 0 0 0 1 0 0 1 1 0 0 0 0 0 1;
        #                 0 1 0 0 0 0 0 1 1 0 0 1 0 0 0 0;
        #                 0 0 0 0 0 0 0 0 1 0 0 0 1 1 1 0;
        #                 1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 1;
        #                 0 0 0 1 0 1 0 0 0 0 1 0 1 0 0 0;
        #                 0 0 1 1 0 0 0 0 0 0 0 1 0 0 0 1;
        #                 0 0 0 0 1 0 0 0 0 1 1 1 0 0 0 0;
        #                 0 1 0 0 0 0 1 0 0 0 1 0 0 0 1 0;
        #                 1 0 1 0 0 0 0 0 0 1 0 0 0 0 1 0
        #         ]
        # )
        # HGP = HypergraphProductCode(LinearCode(H, true))

        # lx, lz = Quintavalle_basis(HGP)
        # # Check the logical operators commute with the stabilizers
        # @test iszero(HGP.X_stabs * transpose(lz))
        # @test iszero(HGP.Z_stabs * transpose(lx))

        # # Overlap should be zero in all cases except i == ii
        # zero_sum_flag = true
        # one_sum_flag = true
        # count_zero_flag = true
        # count_one_flag = true
        # for i in 1:nrows(lx)
        #     for ii in 1:nrows(lz)
        #         if i != ii
        #             iszero(sum(lx[i, :] .* lz[ii, :])) || (zero_sum_flag = false)
        #             iszero(length(findall(x -> x == Fone, lx[i, :] .* lz[ii, :]))) || (count_zero_flag = false)
        #         else
        #             isone(sum(lx[i, :] .* lz[ii, :])) || (one_sum_flag = false)
        #             isone(length(findall(x -> x == Fone, lx[i, :] .* lz[ii, :]))) || (count_one_flag = false)
        #         end
        #     end
        # end
        # @test zero_sum_flag
        # @test one_sum_flag
        # @test count_zero_flag
        # @test count_one_flag

        # # Check the logical operators have weight >= code distance
        # weight_flag = true
        # for i in 1:HGP.k
        #     (wt(lx[i, :]) < HGP.d || wt(lz[i, :]) < HGP.d) && (weight_flag = false)
        # end
        # @test weight_flag

        # # product codes
        # F = Oscar.Nemo.Native.GF(2)
        # h = matrix(F, [1 1])
        # id = identity_matrix(F, 2)
        # H_X = vcat(
        #     h ⊗ h ⊗ h ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id,
        #     id ⊗ id ⊗ id ⊗ h ⊗ h ⊗ h ⊗ id ⊗ id ⊗ id,
        #     id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ h ⊗ h ⊗ h)
        # H_Z = vcat(
        #     h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id,
        #     id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id,
        #     id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h)
        # SPCtest = CodingTheory.SPCDFoldProductCode(3)
        # @test length(SPCtest) == 512
        # @test dimension(SPCtest) == 174
        # @test H_X == SPCtest.X_stabs
        # @test H_Z == SPCtest.Z_stabs
    # end

    @testset "BivariateBicycle Code" begin
        # bivariate bicycle codes
        S, (x, y) = polynomial_ring(Oscar.Nemo.Native.GF(2), [:x, :y])

        # Table 3 of https://arxiv.org/pdf/2308.07915
        # [[72, 12, 6]]
        l = 6
        m = 6
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^3 + y + y^2)
        b = R(y^3 + x + x^2)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 72
        @test dimension(Q) == 12
        # @test minimum_distance(S) == 6

        # [[90, 8, 10]]
        l = 15
        m = 3
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^9 + y + y^2)
        b = R(1 + x^2 + x^7)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 90
        @test dimension(Q) == 8
        # @test minimum_distance(S) == 10

        # [[108, 8, 10]]
        l = 9
        m = 6
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^3 + y + y^2)
        b = R(y^3 + x + x^2)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 108
        @test dimension(Q) == 8
        # @test minimum_distance(S) == 10

        # [[144, 12, 12]]
        l = 12
        m = 6
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^3 + y + y^2)
        b = R(y^3 + x + x^2)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 144
        @test dimension(Q) == 12
        # @test minimum_distance(S) == 12

        # [[288, 12, 18]]
        l = 12
        m = 12
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^3 + y^2 + y^7)
        b = R(y^3 + x + x^2)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 288
        @test dimension(Q) == 12
        # @test minimum_distance(S) == 18

        # [[360, 12, ≤24]]
        l = 30
        m = 6
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^9 + y + y^2)
        b = R(y^3 + x^25 + x^26)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 360
        @test dimension(Q) == 12

        # [[756, 16, ≤34]]
        l = 21
        m = 18
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^3 + y^10 + y^17)
        b = R(y^5 + x^3 + x^19)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 756
        @test dimension(Q) == 16

        # Table 1 of https://arxiv.org/pdf/2408.10001
        # [[54, 8, 6]]
        l = 3
        m = 9
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(1 + y^2 + y^4)
        b = R(y^3 + x + x^2)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 54
        @test dimension(Q) == 8
        # @test minimum_distance(S) == 6

        # [[98, 6, 12]]
        l = 7
        m = 7
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^3 + y^5 + y^6)
        b = R(y^2 + x^3 + x^5)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 98
        @test dimension(Q) == 6
        # @test minimum_distance(S) == 12

        # [[126, 8, 10]]
        l = 3
        m = 21
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(1 + y^2 + y^10)
        b = R(y^3 + x + x^2)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 126
        @test dimension(Q) == 8
    # @test minimum_distance(S) == 10

        # [[150, 16, 8]]
        l = 5
        m = 15
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(1 + y^6 + y^8)
        b = R(y^5 + x + x^4)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 150
        @test dimension(Q) == 16
        # @test minimum_distance(S) == 8

        # [[162, 8, 14]]
        l = 3
        m = 27
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(1 + y^10 + y^14)
        b = R(y^12 + x + x^2)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 162
        @test dimension(Q) == 8
        # @test minimum_distance(S) == 14

        # [[180, 8, 16]]
        l = 6
        m = 15
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^3 + y + y^2)
        b = R(y^6 + x^4 + x^5)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 180
        @test dimension(Q) == 8
        # @test minimum_distance(S) == 16
    end
end
