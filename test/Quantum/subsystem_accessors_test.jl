@testitem "Quantum subsystem accessors" begin
    using Oscar, CodingTheory, Random

    F = GF(2)

    # Constructors, aliases, and the algebraic meaning of bare and dressed operators.
    XG = matrix(F, [
        1 1 0 0
        0 0 1 1
    ])
    ZG = matrix(F, [
        1 0 1 0
        0 1 0 1
    ])
    Q = SubsystemCodeCSS(XG, ZG)
    Q_alias = CSSSubsystemCode(XG, ZG)
    @test Q isa SubsystemCodeCSS
    @test are_equivalent(Q, Q_alias)
    @test (num_qubits(Q), Q.k, Q.r) == (4, 1, 1)
    @test num_X_stabs(Q) + num_Z_stabs(Q) ==
        nrows(X_stabilizers(Q)) + nrows(Z_stabilizers(Q)) ==
        nrows(stabilizers(Q))
    @test X_signs(Q) == CodingTheory.signs(Q)[1:num_X_stabs(Q)]
    @test Z_signs(Q) ==
        CodingTheory.signs(Q)[num_X_stabs(Q) + 1:end]

    L = logicals_matrix(Q)
    Gops = gauges_matrix(Q)
    Ggroup = gauge_group(Q)
    @test bare(Q) == bare_logicals(Q) == logical_operators(Q) == logicals(Q)
    @test gauges(Q) == gauge_operators(Q)
    @test gauge_operators_matrix(Q) == Gops
    @test gauge_group_matrix(Q) == gauge_generators_matrix(Q) ==
        gauge_group_generators_matrix(Q) == Ggroup
    @test are_symplectic_orthogonal(Ggroup, L)
    @test all(is_bare_normalizer(Q, L[i:i, :]) for i in 1:nrows(L))
    @test all(is_bare_logical(Q, L[i:i, :]) for i in 1:nrows(L))
    @test all(are_symplectic_orthogonal(stabilizers(Q), op)
        for pair in dressed(Q) for op in pair)
    @test dressed(Q) == dressed_logicals(Q) == dressed_operators(Q) ==
        logicals(Q) ∪ gauges(Q)
    @test any(!is_bare_normalizer(Q, Gops[i:i, :])
        for i in 1:nrows(Gops))
    @test bare_normalizer_matrix(Q) == gauge_centralizer_matrix(Q)
    @test stabilizer_centralizer_matrix(Q) == normalizer_matrix(Q)
    @test are_symplectic_orthogonal(
        gauge_centralizer_matrix(Q), Ggroup)
    @test are_symplectic_orthogonal(
        stabilizer_centralizer_matrix(Q), stabilizers(Q))

    # CSS syndrome conventions: X checks see Z errors and vice versa.
    x = matrix(F, [1 0 1 0])
    z = matrix(F, [0 1 0 1])
    @test X_syndrome(Q, z) == X_stabilizers(Q) * transpose(z)
    @test Z_syndrome(Q, transpose(x)) == Z_stabilizers(Q) * transpose(x)

    # Standard-form accessors are exact named slices and reassemble the matrix.
    R = random_subsystem_code(MersenneTwister(1201), F, 6, 1, 2)
    Hsf = stabilizers(R; standform=true)
    sr, sk = R.cache[:stand_r], R.cache[:stand_k]
    A, A1, A2 = standard_form_A(R), standard_form_A1(R), standard_form_A2(R)
    B = standard_form_B(R)
    C1, C2 = standard_form_C1(R), standard_form_C2(R)
    D, E = standard_form_D(R), standard_form_E(R)
    @test size(A) == (sr, R.n - sr)
    @test size(A1) == (sr, R.n - sk - sr)
    @test size(A2) == (sr, sk)
    @test size(B) == (sr, sr)
    @test size(C1) == (sr, R.n - sk - sr)
    @test size(C2) == (sr, sk)
    @test size(D) == (R.n - sk - sr, sr)
    @test size(E) == (R.n - sk - sr, sk)
    @test A == hcat(A1, A2)
    top = hcat(identity_matrix(F, sr), A1, A2, B, C1, C2)
    bottom_rows = R.n - sk - sr
    bottom = hcat(
        zero_matrix(F, bottom_rows, sr),
        zero_matrix(F, bottom_rows, R.n - sk - sr),
        zero_matrix(F, bottom_rows, sk),
        D,
        identity_matrix(F, bottom_rows),
        E)
    @test Hsf == vcat(top, bottom)

    # Mutating setters accept equivalent overcomplete presentations.
    S = BaconShorCode(2, 2)
    H0 = stabilizers(S)
    Hdup = vcat(H0, H0)
    set_stabilizers!(S, Hdup)
    @test stabilizers(S) == Hdup

    S = BaconShorCode(2, 2)
    X0, Z0 = X_stabilizers(S), Z_stabilizers(S)
    Xdup, Zdup = vcat(X0, X0), vcat(Z0, Z0)
    set_X_stabilizers!(S, Xdup)
    set_Z_stabilizers!(S, Zdup)
    @test X_stabilizers(S) == Xdup
    @test Z_stabilizers(S) == Zdup

    L0 = logicals_matrix(S)
    Lswap = vcat(L0[2:2, :], L0[1:1, :])
    set_logicals!(S, Lswap)
    @test logicals_matrix(S) == Lswap

    R4, _ = residue_ring(ZZ, 4)
    char_vec = [R4(i == 1 ? 2 : 0) for i in 1:2S.n]
    @test isempty(character_vector(S))
    @test_throws ArgumentError set_signs!(S, [R4(0)])
    set_signs!(S, typeof(R4(0))[])
    @test isempty(character_vector(S))
    @test all(iszero, CodingTheory.signs(S))

    # Valid metachecks annihilate their check matrices; invalid ones are rejected.
    MX = matrix(S.F, [1 1])
    MZ = matrix(S.F, [1 1])
    @test iszero(MX * X_stabilizers(S))
    @test iszero(MZ * Z_stabilizers(S))
    set_X_metacheck!(S, MX)
    set_Z_metacheck!(S, MZ)
    @test X_metacheck(S) == MX
    @test Z_metacheck(S) == MZ
    @test_throws ErrorException set_X_metacheck!(S, matrix(S.F, [1 0]))
    @test_throws ErrorException set_Z_metacheck!(S, matrix(S.F, [0 1]))
    @test_throws ErrorException metacheck(S)
    @test_throws ErrorException set_metacheck!(S, zero_matrix(F, 1, nrows(stabilizers(S))))

    N = FiveQubitCode()
    MN = zero_matrix(N.F, 1, nrows(stabilizers(N)))
    @test iszero(MN * stabilizers(N))
    set_metacheck!(N, MN)
    @test metacheck(N) == MN
    @test_throws ErrorException set_metacheck!(
        N, matrix(N.F, [1 0 0 0]))
    @test_throws ErrorException X_metacheck(N)
    @test_throws ErrorException Z_metacheck(N)
    @test_throws ErrorException set_X_metacheck!(N, zero_matrix(F, 1, 4))
    @test_throws ErrorException set_Z_metacheck!(N, zero_matrix(F, 1, 4))

    # Gauge/logical promotion conserves k+r, with copy and mutating semantics.
    P = BaconShorCode(2, 2)
    total = P.k + P.r
    Pgl = promote_gauges_to_logical(P, [1])
    @test (P.k, P.r) == (1, 1)
    @test (Pgl.k, Pgl.r, Pgl.k + Pgl.r) == (2, 0, total)
    promote_gauges_to_logical!(P, [1])
    @test (P.k, P.r, P.k + P.r) == (2, 0, total)
    Plg = promote_logicals_to_gauge(P, [1])
    @test (P.k, P.r) == (2, 0)
    @test (Plg.k, Plg.r, Plg.k + Plg.r) == (1, 1, total)
    promote_logicals_to_gauge!(P, [1])
    @test (P.k, P.r, P.k + P.r) == (1, 1, total)

    P = BaconShorCode(2, 2)
    oldL = logicals(P)[1]
    swap_X_Z_logicals!(P, [1])
    @test logicals(P)[1] == reverse(oldL)
    oldG = gauges(P)[1]
    swap_X_Z_gauge_operators!(P, [1])
    @test gauges(P)[1] == reverse(oldG)

    fixed = fix_gauge(BaconShorCode(2, 2), 1, :X)
    @test fixed.k == 1
    @test nrows(stabilizers(fixed)) == 1 + nrows(stabilizers(BaconShorCode(2, 2)))
    @test is_stabilizer(fixed, gauges(BaconShorCode(2, 2))[1][1])
    @test_throws ArgumentError fix_gauge(BaconShorCode(2, 2), 1, :Y)

    P = BaconShorCode(2, 2)
    Hbefore, Lbefore, Gbefore = stabilizers(P),
        logicals_matrix(P), gauges_matrix(P)
    permute_code!(P, [2, 1, 4, 3])
    permute_code!(P, [2, 1, 4, 3])
    @test stabilizers(P) == Hbefore
    @test logicals_matrix(P) == Lbefore
    @test gauges_matrix(P) == Gbefore

    # The Bacon-Shor gauge group is generated by weight-two operators.
    @test minimum_gauge_weight(BaconShorCode(2, 2); alg=:bruteforce) == 2

    # Stored distance APIs round-trip values without invoking distance solvers.
    B = BaconShorCode(3, 3)
    set_bare_minimum_distance!(B, 4)
    set_dressed_minimum_distance!(B, 2)
    @test (bare_minimum_distance_lower_bound(B),
        bare_minimum_distance_upper_bound(B)) == (4, 4)
    @test (dressed_minimum_distance_lower_bound(B),
        dressed_minimum_distance_upper_bound(B)) == (2, 2)
    set_bare_X_minimum_distance!(B, 4)
    set_bare_Z_minimum_distance!(B, 4)
    set_dressed_X_minimum_distance!(B, 2)
    set_dressed_Z_minimum_distance!(B, 3)
    @test (bare_X_minimum_distance_lower_bound(B),
        bare_X_minimum_distance_upper_bound(B)) == (4, 4)
    @test (bare_Z_minimum_distance_lower_bound(B),
        bare_Z_minimum_distance_upper_bound(B)) == (4, 4)
    @test (dressed_X_minimum_distance_lower_bound(B),
        dressed_X_minimum_distance_upper_bound(B)) == (2, 2)
    @test (dressed_Z_minimum_distance_lower_bound(B),
        dressed_Z_minimum_distance_upper_bound(B)) == (3, 3)
    @test relative_distance(B) ≈ 2 / 9

    @testset "Nonmutating Setters" begin
        S = BaconShorCode(2, 2)

        # the copy must share the base field, or Nemo rejects the matrices
        for (setter, getter) in ((set_stabilizers, stabilizers),
                                 (set_X_stabilizers, X_stabilizers),
                                 (set_Z_stabilizers, Z_stabilizers))
            T = setter(S, getter(S))
            @test T !== S
            @test (T.n, T.k) == (S.n, S.k)
            @test base_ring(stabilizers(T)) === base_ring(stabilizers(S))
        end
    end

    @testset "Character Vectors" begin
        S = BaconShorCode(2, 2)
        R, _ = residue_ring(Oscar.ZZ, 4)

        set_signs!(S, fill(R(0), 2 * S.n))
        @test length(S.char_vec) == 2 * S.n

        # the phases must live in Z/4 for a binary code
        R_wrong, _ = residue_ring(Oscar.ZZ, 3)
        @test_throws ArgumentError set_signs!(S, fill(R_wrong(0), 2 * S.n))
        @test_throws ArgumentError set_signs!(S, fill(R(0), S.n))
    end

    @testset "Syndromes of Symplectic Vectors" begin
        S = BaconShorCode(2, 2)

        # the X checks read the Z half of a full [X | Z] vector, and vice versa
        v = zero_matrix(S.F, 1, 2 * S.n)
        v[1, S.n + 1] = S.F(1)
        @test X_syndrome(S, v) == X_syndrome(S, v[:, S.n + 1:2 * S.n])
        @test iszero(Z_syndrome(S, v))

        w = zero_matrix(S.F, 1, 2 * S.n)
        w[1, 1] = S.F(1)
        @test Z_syndrome(S, w) == Z_syndrome(S, w[:, 1:S.n])
        @test iszero(X_syndrome(S, w))

        # column vectors are accepted too
        @test iszero(X_syndrome(S, zero_matrix(S.F, 2 * S.n, 1)))
        @test iszero(Z_syndrome(S, zero_matrix(S.F, 2 * S.n, 1)))
    end
end
