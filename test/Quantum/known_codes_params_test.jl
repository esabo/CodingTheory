@testitem "Known quantum code parameters" begin
    using Oscar, CodingTheory, Random

    function check_stabilizer(S, n, k; css::Bool)
        @test (S.n, S.k) == (n, k)
        @test rank(stabilizers(S)) == n - k
        @test are_symplectic_orthogonal(stabilizers(S), stabilizers(S))
        @test is_CSS(S) == css
    end

    function check_subsystem(S, n, k, r; css::Bool)
        @test (S.n, S.k, S.r) == (n, k, r)
        @test rank(stabilizers(S)) == n - k - r
        @test are_symplectic_orthogonal(stabilizers(S), stabilizers(S))
        @test is_CSS(S) == css
        @test rank(gauge_group(S)) == n - k + r
        @test rank(gauge_group(S)) > rank(stabilizers(S))
    end

    gauged_shor = GaugedShorCode()
    check_subsystem(gauged_shor, 9, 1, 4; css=true)
    @test gauged_shor.d_dressed == 3

    F = GF(2)
    for A in (
        matrix(F, [1 1; 1 1]),
        matrix(F, [1 1 0; 1 1 1; 0 1 1]),
    )
        generalized = GeneralizedBaconShorCode(A)
        n = count(!iszero, A)
        k = rank(A)
        r = n + k - nrows(A) - ncols(A)
        check_subsystem(generalized, n, k, r; css=true)
        @test generalized.d_dressed == 2
    end

    for (constructor, n, k, d) in (
        (QC6, 6, 2, 2),
        (ShorCode, 9, 1, 3),
        (QC4, 4, 2, 2),
        (Q832, 8, 3, 2),
    )
        S = constructor()
        check_stabilizer(S, n, k; css=true)
        @test minimum_distance(S)[1] == d
    end

    # GrossCode is comparatively expensive, so exercise its unique fixed size once.
    gross = GrossCode()
    check_stabilizer(gross, 144, 12; css=true)
    @test minimum_distance(gross)[1] == 12

    for L in (2, 3)
        xcube = XCubeModel(L)
        check_stabilizer(xcube, 3 * L^3, 6 * L - 3; css=true)
        @test minimum_distance(xcube)[1] == L
    end

    for L in (3, 6)
        toric_color = ToricColorCode666(L)
        check_stabilizer(toric_color, 2 * L^2, 4; css=true)
    end

    cleve_gottesman = CleveGottesmanCode()
    check_stabilizer(cleve_gottesman, 8, 3; css=false)
    @test cleve_gottesman.d == 3

    for d in (3, 4)
        twist = TwistDefectSurfaceCode(d)
        check_stabilizer(twist, d^2, 2; css=false)
    end

    for d in (2, 3)
        num_vertical_edges = cld(d^2, 2)
        heavy_hex = HeavyHexCode(d)
        check_subsystem(
            heavy_hex,
            2 * d^2 + num_vertical_edges,
            num_vertical_edges,
            2 * num_vertical_edges;
            css=true,
        )
    end

    for d in (2, 3)
        heavy_square = HeavySquareCode(d)
        check_subsystem(heavy_square, 3 * d^2, d^2, 2 * d^2; css=true)
    end

    # the triangular lattice places a qubit on each of the 3L^2 edges and, being
    # on a torus, encodes two logical qubits
    for L in (2, 3)
        triangular = TriangularSurfaceCode(L)
        check_stabilizer(triangular, 3 * L^2, 2; css=true)
    end

    # the 4D toric code on a 4-torus encodes binomial(4, 2) logical qubits and
    # carries metachecks, which must annihilate their stabilizer matrices
    toric_4D = ToricCode4D(2)
    check_stabilizer(toric_4D, 96, 6; css=true)
    @test iszero(X_metacheck(toric_4D) * X_stabilizers(toric_4D))
    @test iszero(Z_metacheck(toric_4D) * Z_stabilizers(toric_4D))
    @test nrows(X_metacheck(toric_4D)) > 0

    planar_3D = PlanarSurfaceCode3D(2)
    check_stabilizer(planar_3D, 36, 15; css=true)

    toric_3D = ToricCode3D(2)
    check_stabilizer(toric_3D, 24, 3; css=true)

    check_subsystem(LocalBravyiBaconShorCode(matrix(F, [1 1; 1 1])), 4, 1, 1;
        css=true)
    check_subsystem(AugmentedBravyiBaconShorCode(matrix(F, [1 1; 1 1])), 4, 1, 1;
        css=true)
    check_subsystem(NappPreskill3DCode(2, 2, 2), 8, 1, 3; css=true)
    check_subsystem(SubsystemToricCode(2, 2), 12, 2, 4; css=true)
end
