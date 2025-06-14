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
        # g = gcd(a, b, x^l - 1)
        Q = GeneralizedBicycleCode(R(a), R(b))
        @test length(Q) == 254
        @test dimension(Q) == 28
        @test LogicalTrait(typeof(Q)) == HasLogicals()
        @test GaugeTrait(typeof(Q)) == HasNoGauges()

        # Example A2
        l = 63
        R, _ = residue_ring(S, x^l - 1)
        a = 1 + x + x^14 + x^16 + x^22
        b = 1 + x^3 + x^13 + x^20 + x^42
        Q = GeneralizedBicycleCode(R(a), R(b))
        @test length(Q) == 126
        @test dimension(Q) == 28

        # Example A3
        l = 24
        R, _ = residue_ring(S, x^l - 1)
        a = 1 + x^2 + x^8 + x^15
        b = 1 + x^2 + x^12 + x^17
        Q = GeneralizedBicycleCode(R(a), R(b))
        @test length(Q) == 48
        @test dimension(Q) == 6

        # Example A4
        l = 23
        R, _ = residue_ring(S, x^l - 1)
        a = 1 + x^5 + x^8 + x^12
        b = 1 + x + x^5 + x^7
        Q = GeneralizedBicycleCode(R(a), R(b))
        @test length(Q) == 46
        @test dimension(Q) == 2

        # Example A5
        l = 90
        R, _ = residue_ring(S, x^l - 1)
        a = 1 + x^28 + x^80 + x^89
        b = 1 + x^2 + x^21 + x^25
        Q = GeneralizedBicycleCode(R(a), R(b))
        @test length(Q) == 180
        @test dimension(Q) == 10

        # Example A6
        l = 450
        R, _ = residue_ring(S, x^l - 1)
        a = 1 + x^97 + x^372 + x^425
        b = 1 + x^50 + x^265 + x^390
        Q = GeneralizedBicycleCode(R(a), R(b))
        @test length(Q) == 900
        @test dimension(Q) == 50
        # d = 15

        # Figure 9
        # [[126, 12, d < 11]]
        l = 63
        R, _ = residue_ring(S, x^l - 1)
        a = 1 + x^43 + x^37
        b = 1 + x^59 + x^31
        Q = GeneralizedBicycleCode(R(a), R(b))
        @test length(Q) == 126
        @test dimension(Q) == 12

        # Figure 9
        # [[254, 14, d < 17]]
        l = 127
        R, _ = residue_ring(S, x^l - 1)
        a = 1 + x^18 + x^53
        b = 1 + x^12 + x^125
        Q = GeneralizedBicycleCode(R(a), R(b))
        @test length(Q) == 254
        @test dimension(Q) == 14

        # Figure 9
        # [[510, 16, d < 19]]
        l = 255
        R, _ = residue_ring(S, x^l - 1)
        a = 1 + x^250 + x^133
        b = 1 + x^41 + x^157
        Q = GeneralizedBicycleCode(R(a), R(b))
        @test length(Q) == 510
        @test dimension(Q) == 16

        # Figure 9
        # [[1022, 18]]
        l = 511
        R, _ = residue_ring(S, x^l - 1)
        a = 1 + x^244 + x^157
        b = 1 + x^51 + x^454
        Q = GeneralizedBicycleCode(R(a), R(b))
        @test length(Q) == 1022
        @test dimension(Q) == 18

        # probably takes too long to compute
        # # Figure 9
        # # [[8190, 24]]
        # l = 4095
        # R, _ = residue_ring(S, x^l - 1)
        # a = 1 + x^2083 + x^2212
        # b = 1 + x^1802 + x^3220
        # Q = GeneralizedBicycleCode(R(a), R(b))
        # @test length(Q) == 8190
        # @test dimension(Q) == 24
    end

    @testset "LiftedProductCode" begin
        # Degenerate Quantum LDPC Codes With Good Finite Length Performance
        # Example B1
        F = Oscar.Nemo.Native.GF(2)
        S, x = polynomial_ring(F, :x)
        l = 63
        R, _ = residue_ring(S, x^l - 1)
        A = matrix(
            R,
            7,
            7,
            [
                x^27 0 0 0 0 1 x^54;
                x^54 x^27 0 0 0 0 1;
                1 x^54 x^27 0 0 0 0;
                0 1 x^54 x^27 0 0 0;
                0 0 1 x^54 x^27 0 0;
                0 0 0 1 x^54 x^27 0;
                0 0 0 0 1 x^54 x^27
            ],
        )
        b = R(1 + x + x^6)
        Q = LiftedProductCode(A, b)
        @test length(Q) == 882
        @test dimension(Q) == 24
        @test LogicalTrait(typeof(Q)) == HasLogicals()
        @test GaugeTrait(typeof(Q)) == HasNoGauges()

        # Example B2
        A = matrix(
            R,
            7,
            7,
            [
                x^27 0 0 1 x^18 x^27 1;
                1 x^27 0 0 1 x^18 x^27;
                x^27 1 x^27 0 0 1 x^18;
                x^18 x^27 1 x^27 0 0 1;
                1 x^18 x^27 1 x^27 0 0;
                0 1 x^18 x^27 1 x^27 0;
                0 0 1 x^18 x^27 1 x^27
            ],
        )
        Q = LiftedProductCode(A, b)
        @test length(Q) == 882
        @test dimension(Q) == 48

        # Example B3
        # Example B3
        l = 127
        R, _ = residue_ring(S, x^l - 1)
        A = matrix(
            R,
            5,
            5,
            [
                1 0 x^51 x^52 0;
                0 1 0 x^111 x^20;
                1 0 x^98 0 x^122;
                1 x^80 0 x^119 0;
                0 1 x^5 0 x^106
            ],
        )
        b = R(1 + x + x^7)
        Q = LiftedProductCode(A, b)
        @test length(Q) == 1270
        @test dimension(Q) == 28
    end

    @testset "HypergraphProductCode" begin
        # HGP codes are (l, q)-QLDPC code
        # l, q = column_row_weights(Sq2)

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

        # Degenerate Quantum LDPC Codes With Good Finite Length Performance
        # # Example C1
        # l = 63
        # R, _ = residue_ring(S, x^l - 1)
        # h = R(1 + x^3 + x^34 + x^41 + x^57)
        # A = residue_polynomial_to_circulant_matrix(h)
        # Q = HypergraphProductCode(A, A)
        # @test length(Q) == 7938
        # @test dimension(Q) == 578
        # takes awhile to run

        # Example C2
        F, x = polynomial_ring(Oscar.Nemo.Native.GF(2), :x)
        l = 31
        R, = residue_ring(F, x^l - 1)
        h = R(1 + x^2 + x^5)
        A = residue_polynomial_to_circulant_matrix(h)
        Q = HypergraphProductCode(A, A)
        @test length(Q) == 1922
        @test dimension(Q) == 50
        # [[1922, 50]]_2 subsystem code

        # Quintavalle_basis
        Fone = F(1)
        for r = 2:4
            C = HammingCode(2, r)
            Fone = F(1)
            H = parity_check_matrix(C)
            nr_H = nrows(H)

            # Add some degerate rows so the dual code has k > 0
            rand1 = rand(1:nr_H)
            rand2 = rand(1:nr_H)
            rand3 = rand(1:nr_H)
            rand4 = rand(1:nr_H)
            H = vcat(
                H,
                H[rand1:rand1, :] + H[rand2:rand2, :],
                H[rand3:rand3, :] + H[rand4:rand4, :],
            )
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
            for i = 1:nrows(lx)
                for ii = 1:nrows(lz)
                    if i != ii
                        iszero(sum(lx[i, :] .* lz[ii, :])) || (zero_sum_flag = false)
                        iszero(length(findall(x -> x == Fone, lx[i, :] .* lz[ii, :]))) ||
                            (count_zero_flag = false;)
                    else
                        isone(sum(lx[i, :] .* lz[ii, :])) || (one_sum_flag = false;)
                        isone(length(findall(x -> x == Fone, lx[i, :] .* lz[ii, :]))) ||
                            (count_one_flag = false;)
                    end
                end
            end
            @test zero_sum_flag
            @test one_sum_flag
            @test count_zero_flag
            @test count_one_flag

            # Check the logical operators have weight >= code distance
            weight_flag = true
            for i = 1:HGP.k
                (wt(lx[i:i, :]) < HGP.d || wt(lz[i:i, :]) < HGP.d) && (weight_flag = false;)
            end
            @test weight_flag
        end

        H1 = parity_check_matrix(HammingCode(2, 2))
        H2 = parity_check_matrix(HammingCode(2, 3))
        H1 = vcat(H1, H1[1:1, :] + H1[2:2, :])
        H2 = vcat(H2, H2[1:1, :] + H2[2:2, :])
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
        for i = 1:nrows(lx)
            for ii = 1:nrows(lz)
                if i != ii
                    iszero(sum(lx[i, :] .* lz[ii, :])) || (zero_sum_flag = false)
                    iszero(length(findall(x -> x == Fone, lx[i, :] .* lz[ii, :]))) ||
                        (count_zero_flag = false)
                else
                    isone(sum(lx[i, :] .* lz[ii, :])) || (one_sum_flag = false)
                    isone(length(findall(x -> x == Fone, lx[i, :] .* lz[ii, :]))) ||
                        (count_one_flag = false)
                end
            end
        end
        @test zero_sum_flag
        @test one_sum_flag
        @test count_zero_flag
        @test count_one_flag

        # Check the logical operators have weight >= code distance
        weight_flag = true
        for i = 1:HGP.k
            (wt(lx[i, :]) < HGP.d || wt(lz[i, :]) < HGP.d) && (weight_flag = false)
        end
        @test weight_flag

        # [[400,16,6]] code from Table 1 of https://doi.org/10.1103/PhysRevResearch.2.043423
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(
            F,
            12,
            16,
            [
                1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0;
                0 0 1 0 0 0 1 1 0 0 0 0 1 0 0 0;
                0 0 0 1 1 0 1 0 0 0 0 0 0 1 0 0;
                0 0 0 0 0 1 0 0 1 1 0 0 0 0 0 1;
                0 1 0 0 0 0 0 1 1 0 0 1 0 0 0 0;
                0 0 0 0 0 0 0 0 1 0 0 0 1 1 1 0;
                1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 1;
                0 0 0 1 0 1 0 0 0 0 1 0 1 0 0 0;
                0 0 1 1 0 0 0 0 0 0 0 1 0 0 0 1;
                0 0 0 0 1 0 0 0 0 1 1 1 0 0 0 0;
                0 1 0 0 0 0 1 0 0 0 1 0 0 0 1 0;
                1 0 1 0 0 0 0 0 0 1 0 0 0 0 1 0
            ],
        )
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
        for i = 1:nrows(lx)
            for ii = 1:nrows(lz)
                if i != ii
                    iszero(sum(lx[i, :] .* lz[ii, :])) || (zero_sum_flag = false)
                    iszero(length(findall(x -> x == Fone, lx[i, :] .* lz[ii, :]))) ||
                        (count_zero_flag = false)
                else
                    isone(sum(lx[i, :] .* lz[ii, :])) || (one_sum_flag = false)
                    isone(length(findall(x -> x == Fone, lx[i, :] .* lz[ii, :]))) ||
                        (count_one_flag = false)
                end
            end
        end
        @test zero_sum_flag
        @test one_sum_flag
        @test count_zero_flag
        @test count_one_flag

        # Check the logical operators have weight >= code distance
        weight_flag = true
        for i = 1:HGP.k
            (wt(lx[i, :]) < HGP.d || wt(lz[i, :]) < HGP.d) && (weight_flag = false)
        end
        @test weight_flag

        # TODO codes from "LRESC" paper, make tests (long range _ surface code?)
        # # the dual of any [3, 2, 2] code is a [3, 1, 3] repetition code
        # # but we'll build this the GAP way anyway
        # C_GAP = GAP.Globals.BestKnownLinearCode(3, 2, GAP.Globals.GF(2));
        # G = GAP.Globals.GeneratorMat(C_GAP);
        # y = [GAP.Globals.Int(G[i, j]) for i in 1:2, j in 1:3];
        # C32 = LinearCode(y, 2);
        # R4 = RepetitionCode(2, 4);
        # seed1 = concatenate(C32, R4)
        # LRESC1 = HypergraphProductCode(parity_check_matrix(seed1), parity_check_matrix(seed1))

        # # unclear as to which [6, 2, 4] code they used here
        # C_GAP = GAP.Globals.BestKnownLinearCode(6, 2, GAP.Globals.GF(2));
        # G = GAP.Globals.GeneratorMat(C_GAP);
        # y = [GAP.Globals.Int(G[i, j]) for i in 1:2, j in 1:6];
        # C62 = LinearCode(y, 2)
        # R2 = RepetitionCode(2, 2);
        # seed2 = concatenate(C62, R2)
        # LRESC2 = HypergraphProductCode(parity_check_matrix(seed2), parity_check_matrix(seed2))

        # # they either could have used this code or the extended Hamming code below
        # C_GAP = GAP.Globals.BestKnownLinearCode(8, 4, GAP.Globals.GF(2));
        # G = GAP.Globals.GeneratorMat(C_GAP);
        # y = [GAP.Globals.Int(G[i, j]) for i in 1:4, j in 1:8];
        # C84 = LinearCode(y, 2)
        # R3 = RepetitionCode(2, 3);
        # seed3 = concatenate(C84, R3)
        # LRESC3 = HypergraphProductCode(parity_check_matrix(seed3), parity_check_matrix(seed3))

        # ext_Ham = extend(HammingCode(2, 3))
        # seed3_alt = concatenate(ext_Ham, R3)
        # LRESC3_alt = HypergraphProductCode(parity_check_matrix(seed3_alt), parity_check_matrix(seed3_alt))

        # julia> are_permutation_equivalent(C84, ext_Ham)
        # (false, missing)
    end

    @testset "Product Codes" begin
        # [512, 174, 8]] SPCDFoldProductCode from Section 4 of https://arxiv.org/pdf/2209.13474
        F = Oscar.Nemo.Native.GF(2)
        h = matrix(F, [1 1])
        id = identity_matrix(F, 2)
        H_X = vcat(
            h ⊗ h ⊗ h ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id,
            id ⊗ id ⊗ id ⊗ h ⊗ h ⊗ h ⊗ id ⊗ id ⊗ id,
            id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ h ⊗ h ⊗ h,
        )
        H_Z = vcat(
            h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id,
            id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id,
            id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h,
        )
        SPCtest = SPCDFoldProductCode(3)
        @test length(SPCtest) == 512
        @test dimension(SPCtest) == 174
        # TODO compute this directly
        @test_broken minimum_distance(SPCtest) == 8
        @test H_X == SPCtest.X_stabs
        @test H_Z == SPCtest.Z_stabs
    end

    @testset "BivariateBicycleCode" begin
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
        @test_broken minimum_distance(Q) == 6

        # [[90, 8, 10]]
        l = 15
        m = 3
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^9 + y + y^2)
        b = R(1 + x^2 + x^7)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 90
        @test dimension(Q) == 8
        @test_broken minimum_distance(S) == 10

        # [[108, 8, 10]]
        l = 9
        m = 6
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^3 + y + y^2)
        b = R(y^3 + x + x^2)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 108
        @test dimension(Q) == 8
        @test_broken minimum_distance(S) == 10

        # [[144, 12, 12]]
        l = 12
        m = 6
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^3 + y + y^2)
        b = R(y^3 + x + x^2)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 144
        @test dimension(Q) == 12
        @test_broken minimum_distance(S) == 12

        # [[288, 12, 18]]
        l = 12
        m = 12
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^3 + y^2 + y^7)
        b = R(y^3 + x + x^2)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 288
        @test dimension(Q) == 12
        @test_broken minimum_distance(S) == 18

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
        @test_broken minimum_distance(S) == 6

        # [[98, 6, 12]]
        l = 7
        m = 7
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^3 + y^5 + y^6)
        b = R(y^2 + x^3 + x^5)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 98
        @test dimension(Q) == 6
        @test_broken minimum_distance(S) == 12

        # [[126, 8, 10]]
        l = 3
        m = 21
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(1 + y^2 + y^10)
        b = R(y^3 + x + x^2)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 126
        @test dimension(Q) == 8
        @test_broken minimum_distance(S) == 10

        # [[150, 16, 8]]
        l = 5
        m = 15
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(1 + y^6 + y^8)
        b = R(y^5 + x + x^4)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 150
        @test dimension(Q) == 16
        @test_broken minimum_distance(S) == 8

        # [[162, 8, 14]]
        l = 3
        m = 27
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(1 + y^10 + y^14)
        b = R(y^12 + x + x^2)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 162
        @test dimension(Q) == 8
        @test_broken minimum_distance(S) == 14

        # [[180, 8, 16]]
        l = 6
        m = 15
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^3 + y + y^2)
        b = R(y^6 + x^4 + x^5)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 180
        @test dimension(Q) == 8

        # Table 1 of https://arxiv.org/pdf/2404.17676
        # [[72, 8, 6]]
        l = 12
        m = 3
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^9 + y^1 + y^2)
        b = R(1 + x^1 + x^11)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 72
        @test dimension(Q) == 8
        @test_broken minimum_distance(S) == 6

        # [[90, 8, 6]]
        l = 9
        m = 5
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^8 + y^4 + y^1)
        b = R(y^5 + x^8 + x^7)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 90
        @test dimension(Q) == 8
        @test_broken minimum_distance(S) == 6

        # [[120, 8, 8]]
        l = 12
        m = 5
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^10 + y^4 + y^1)
        b = R(1 + x^1 + x^2)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 120
        @test dimension(Q) == 8
        @test_broken minimum_distance(S) == 8

        # [[150, 8, 8]]
        l = 15
        m = 5
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^5 + y^2 + y^3)
        b = R(y^2 + x^7 + x^6)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 150
        @test dimension(Q) == 8
        @test_broken minimum_distance(S) == 8

        # [[196, 12, 8]]
        l = 14
        m = 7
        R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
        a = R(x^6 + y^5 + y^6)
        b = R(1 + x^4 + x^13)
        Q = BivariateBicycleCode(a, b)
        @test length(Q) == 196
        @test dimension(Q) == 12
        @test_broken minimum_distance(S) == 8
    end

    @testset "CoprimeBivariateBicycleCode" begin
        # coprime bivariate bicycle codes
        S, P = polynomial_ring(Oscar.Nemo.Native.GF(2), :P)

        # Table 2 of https://arxiv.org/pdf/2408.10001v1
        # [[30, 4, 6]]
        l = 3
        m = 5
        R, _ = residue_ring(S, P^(l * m) - 1)
        a = R(1 + P + P^2)
        b = R(P + P^3 + P^8)
        Q = CoprimeBivariateBicycleCode(a, b)
        @test length(Q) == 30
        @test dimension(Q) == 4
        @test_broken minimum_distance(Q) == 6

        # [[42, 6, 6]]
        l = 3
        m = 7
        R, _ = residue_ring(S, P^(l * m) - 1)
        a = R(1 + P^2 + P^3)
        b = R(P + P^3 + P^11)
        Q = CoprimeBivariateBicycleCode(a, b)
        @test length(Q) == 42
        @test dimension(Q) == 6
        @test_broken minimum_distance(Q) == 6

        # [[70, 6, 8]]
        l = 5
        m = 7
        R, _ = residue_ring(S, P^(l * m) - 1)
        a = R(1 + P + P^5)
        b = R(1 + P + P^12)
        Q = CoprimeBivariateBicycleCode(a, b)
        @test length(Q) == 70
        @test dimension(Q) == 6
        @test_broken minimum_distance(Q) == 8

        # [[108, 12, 6]]
        l = 2
        m = 27
        R, _ = residue_ring(S, P^(l * m) - 1)
        a = R(P^2 + P^5 + P^44)
        b = R(P^8 + P^14 + P^47)
        Q = CoprimeBivariateBicycleCode(a, b)
        @test length(Q) == 108
        @test dimension(Q) == 12
        @test_broken minimum_distance(Q) == 6

        # [[126, 12, 10]]
        l = 7
        m = 9
        R, _ = residue_ring(S, P^(l * m) - 1)
        a = R(1 + P + P^58)
        b = R(P^3 + P^16 + P^44)
        Q = CoprimeBivariateBicycleCode(a, b)
        @test length(Q) == 126
        @test dimension(Q) == 12
        @test_broken minimum_distance(Q) == 10
    end

    @testset "BiasTailoredLiftedProductCode" begin
        # [[882, 24, d ≤ 24]] from Appendix B of https://arxiv.org/pdf/2202.01702
        F = Oscar.Nemo.Native.GF(2)
        S, x = polynomial_ring(F, :x)
        l = 63
        R, _ = residue_ring(S, x^l - 1)
        A1 = matrix(R, 1, 1, [1 + x^1 + x^6])
        A2 = matrix(
            R,
            7,
            7,
            [
                x^36 0 0 0 0 1 x^9;
                x^9 x^36 0 0 0 0 1;
                1 x^9 x^36 0 0 0 0;
                0 1 x^9 x^36 0 0 0;
                0 0 1 x^9 x^36 0 0;
                0 0 0 1 x^9 x^36 0;
                0 0 0 0 1 x^9 x^36
            ],
        )
        Q = BiasTailoredLiftedProductCode(A1, A2)
        @test length(Q) == 882
        @test dimension(Q) == 24
        @test_broken minimum_distance(Q) == 24

        # [[416, 18, d ≤ 20]] from Example 4.1 of https://arxiv.org/pdf/2202.01702
        S, x = polynomial_ring(F, :x)
        l = 13
        R, _ = residue_ring(S, x^l - 1)
        A1 = matrix(
            R,
            4,
            4,
            [
                1 x^11 x^7 x^12;
                x^1 x^8 x^2 x^8;
                x^11 1 x^4 x^8;
                x^6 x^1 x^4 x^12
            ],
        )
        A2 = A1
        Q = BiasTailoredLiftedProductCode(A1, A2)
        @test length(Q) == 416
        @test dimension(Q) == 18
        @test_broken minimum_distance(Q) == 20
    end
end
