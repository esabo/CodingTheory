# misc codes from the literature
@testitem "Quantum/misc_lit_codes.jl" begin
    using Oscar
    using CodingTheory

    @testset "Misc Codes From Literature" begin
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





        # radial codes
        # using StatsBase, Graphs, Oscar, CodingTheory
        # F = Oscar.Nemo.Native.GF(2)
        # S, x = polynomial_ring(F, :x)
        # l = 5
        # R, _ = residue_ring(S, x^l - 1)
        # A1 = matrix(R, 3, 3,
        #     [x^3, x^2, x,
        #     x^4, x^1, x^4,
        #     x, x^2, x^3])
        # A2 = matrix(R, 3, 3,
        #     [x^3, x^3, 1,
        #     x, 1, x,
        #     x^4, x^2, 1])
        # code = LiftedProductCode(A1, A2)
        # # [[90, 8]]_2 CSS stabilizer code
        # # julia> check_weights(code)
        # # (w_X, q_X, w_Z, q_Z) = (6, 3, 6, 3)
        # L_X = LDPCCode(code.X_stabs);
        # G_X, _, _ = Tanner_graph(L_X);
        # D_X = DiGraph(G_X);
        # L_X_cycles = CodingTheory._modified_hawick_james(D_X, 12);
        # countmap(length.(L_X_cycles))
        # # Dict{Int64, Int64} with 3 entries:
        # #   6  => 48
        # #   10 => 4842
        # #   8  => 777
        # L_Z = LDPCCode(code.Z_stabs);
        # G_Z, _, _ = Tanner_graph(L_Z);
        # D_Z = DiGraph(G_Z);
        # L_Z_cycles = CodingTheory._modified_hawick_james(D_Z, 12);
        # countmap(length.(L_Z_cycles))
        # # Dict{Int64, Int64} with 3 entries:
        # #   6  => 20
        # #   10 => 4074
        # #   8  => 678

        # l = 11
        # R, _ = residue_ring(S, x^l - 1)
        # A1 = matrix(R, 4, 4,
        #     [x^10, x^10, x, x^6,
        #     x^4, x^7, x^5, x^2,
        #     x^8, x^10, x^6, x^9,
        #     x, x^6, 1, x^6])
        # A2 = matrix(R, 4, 4,
        #     [x^9, x^5, x^8, x^3,
        #     x^5, x^4, x, 1,
        #     1, x^4, x^6, x^10,
        #     x^2, x^8, x^4, x^2])
        # code = LiftedProductCode(A1, A2)
        # # [[352, 18]]_2 CSS stabilizer code
        # # julia> check_weights(code)
        # # (w_X, q_X, w_Z, q_Z) = (8, 4, 8, 4)
        # L_X = LDPCCode(code.X_stabs);
        # G_X, _, _ = Tanner_graph(L_X);
        # D_X = DiGraph(G_X);
        # L_X_cycles = CodingTheory._modified_hawick_james(D_X, 12);
        # countmap(length.(L_X_cycles))
        # # Dict{Int64, Int64} with 3 entries:
        # #   6  => 96
        # #   10 => 61668
        # #   8  => 3970
        # L_Z = LDPCCode(code.Z_stabs);
        # G_Z, _, _ = Tanner_graph(L_Z);
        # D_Z = DiGraph(G_Z);
        # L_Z_cycles = CodingTheory._modified_hawick_james(D_Z, 12);
        # countmap(length.(L_Z_cycles))
        # # Dict{Int64, Int64} with 3 entries:
        # #   6  => 72
        # #   10 => 60787
        # #   8  => 3954
    end
end
