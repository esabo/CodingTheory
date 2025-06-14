@testitem "Classical/quasi-cyclic_code.jl" begin
    using Oscar, CodingTheory

    @testset "QuasiCyclicCode" begin
        F = Oscar.Nemo.Native.GF(2)
        v = matrix(F, 1, 8, [1, 0, 1, 1, 1, 0, 0, 0])
        v2 = matrix(F, 1, 8, [1, 1, 1, 0, 0, 0, 1, 0])
        C = QuasiCyclicCode([v, v2], 2, false)
        # v and v2 are shifts of each other
        @test C.k == 4
        v2 = matrix(F, 1, 8, [1, 0, 1, 0, 1, 0, 1, 0])
        C = QuasiCyclicCode([v, v2], 2, false)
        v = [
            matrix(F, 1, 4, [1, 0, 1, 1]),
            matrix(F, 1, 4, [0, 0, 0, 1]),
            matrix(F, 1, 4, [1, 1, 1, 1]),
            matrix(F, 1, 4, [0, 0, 0, 0]),
        ]
        C2 = QuasiCyclicCode(v, 2, true)
        @test are_equivalent(C, C2)

        # Quantum LDPC Codes with Almost Linear Minimum Distance
        # Example 4
        S, x = polynomial_ring(F, :x)
        l = 31
        R, _ = residue_ring(S, x^l - 1)
        A = matrix(
            R,
            3,
            5,
            [
                x x^2 x^4 x^8 x^16;
                x^5 x^10 x^20 x^9 x^18;
                x^25 x^19 x^7 x^14 x^28
            ],
        )
        # A = matrix(R, 3, 5,
        #     [x^(l - 1), x^(l - 2), x^(l - 4), x^(l - 8), x^(l - 16),
        #      x^(l - 5), x^(l - 10), x^(l - 20), x^(l - 9), x^(l - 18),
        #      x^(l - 25), x^(l - 19), x^(l - 7), x^(l - 14), x^(l - 28)])
        # A_tr = CodingTheory._CT_adjoint(A)
        C = QuasiCyclicCode(A, true)
        @test length(C) == 155
        # BUG I can't get this to work. This paper cites a previous paper which gives n, k but no A
        # So maybe this is incorrect. However, the ref uses right shifts whereas I am using left
        # F = Oscar.Nemo.Native.GF(2)
        # S, x = polynomial_ring(F, :x)
        # l = 5
        # R, _ = residue_ring(S, x^l - 1)
        # residue_polynomial_to_circulant_matrix(R(x))
        # residue_polynomial_to_circulant_matrix(R(x^2))
        # which is probably the difference. I tried accounting for this by shifting the exponents
        # but it doesn't seem to help
        @test_broken dimension(C) == 64
        @test index(C) == 5

        # TODO I have the following code for something, check the papers and make tests out of them
        # A_tr = CodingTheory._CT_adjoint(A)
        # LiftedProductCode(A, A_tr)


        # Bias-Tailored Quantum LDPC Codes
        # Example 2.2
        # l = 3
        # R, _ = residue_ring(S, x^l - 1)
        # A = matrix(R, 2, 3, [x + x^2, 1, 0, 0, 1 + x^2, x^2])
        # lift(A)
        # weightmatrix(A)

        # # julia> lift(A)
        # # [0   1   1   1   0   0   0   0   0]
        # # [1   0   1   0   1   0   0   0   0]
        # # [1   1   0   0   0   1   0   0   0]
        # # [0   0   0   1   1   0   0   1   0]
        # # [0   0   0   0   1   1   0   0   1]
        # # [0   0   0   1   0   1   1   0   0]
        # #
        # # julia> weightmatrix(A)
        # # 2×3 Matrix{Int64}:
        # #  2  1  0
        # #  0  2  1

        # # Example 3.3
        # l = 13
        # R, _ = residue_ring(S, x^l - 1)
        # A = matrix(R, 4, 4,
        #     [1, x^2, x^6, x,
        #      x^12, x^5, x^12, x^5,
        #      x^2, 1, x^9, x^5,
        #      x^7, x^11, x^9, x])
        # A_tr = CodingTheory._CT_adjoint(A)
        # code = LiftedQuasiCyclicLiftedProductCode(A, A_tr)

        # # TODO should remove from classical but put where?
        # # TODO this example was never finished because of the original paper? or is wrong?
        # # Quasi-cyclic constructions of quantum codes
        # # Example 1
        # n = 151
        # q = 2
        # deg = ord(n, q)
        # E = GF(2, 15, :α)
        # S, x = polynomial_ring(E, :x)
        # β = α^(div(q^deg - 1, n))
        # def_set_f = cyclotomic_coset(2, q, n)
        # f = CodingTheory._generator_polynomial(S, β, def_set_f)
        # def_set_g = [def_set_f; cyclotomic_coset(10, q, n)]
        # g = CodingTheory._generator_polynomial(S, β, def_set_g)
        # h = x + 1
        # R, _ = residue_ring(S, x^n - 1)
        # G1 = CodingTheory._generator_matrix(E, n, n - length(def_set_f), f);
        # G3 = CodingTheory._generator_matrix(E, n, n - length(def_set_f) - 1, f * h);
        # G2 = CodingTheory._generator_matrix(E, n, n - length(def_set_g), g);
        # # G1 = polytocircmatrix(R(f));
        # # G3 = polytocircmatrix(R(f * h));
        # # G2 = polytocircmatrix(R(g));
        # # cannot concatenate G1 and G3, although in the paper G1 and G2 are defined
        # # but G3 is never actually specified
        # stabs = vcat(hcat(G1, G3), hcat(zero(G2), G2))
        # Q = QuantumCode(change_base_ring(F, stabs), true)
        # # @test length(Q) ==
        # # @test dimension(Q) == 
    end
end
