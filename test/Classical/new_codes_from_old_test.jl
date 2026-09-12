@testitem "Classical/new_codes_from_old.jl" begin
    using Oscar, CodingTheory

    @testset "Plotkin and u+v Constructions" begin
        F = Oscar.Nemo.Native.GF(2)
        C1 = HammingCode(2, 3)       # [7, 4, 3]
        C2 = RepetitionCode(2, 7)    # [7, 1, 7]
        
        C_plot = u_u_plus_v(C1, C2)
        @test length(C_plot) == 14
        @test dimension(C_plot) == 5
        
        # Test lazy alias
        @test are_equivalent(C_plot, Plotkin_construction(C1, C2))
        
        C_uvw = u_plus_w_v_plus_w_u_plus_v_plus_w(C1, C2)
        @test length(C_uvw) == 21
        # dimension is C1.k + C1.k + C2.k = 4 + 4 + 1 = 9
        @test dimension(C_uvw) == 9 
        
        C_bad = RepetitionCode(2, 5)
        @test_throws ArgumentError u_u_plus_v(C1, C_bad) # length mismatch
        @test_throws ArgumentError u_plus_w_v_plus_w_u_plus_v_plus_w(C1, C_bad)
    end

    @testset "Brouwer's Recursive Constructions" begin
        F = Oscar.Nemo.Native.GF(2)
        C = HammingCode(2, 3) # [7, 4, 3]
        
        # Construction A: Need codeword c in C
        G = generator_matrix(C)
        c = vec(Array(G[1, :])) 
        w = wt(c)
        
        C_A = construction_A(C, c)
        @test length(C_A) == 7 - w
        @test dimension(C_A) == 3
        @test_throws ArgumentError construction_A(C, zeros(Int, 7)) # cannot be zero vector
        
        # Construction B: Need dual codeword h in dual(C)
        H = parity_check_matrix(C)
        h = vec(Array(H[1, :]))
        s = wt(h)
        
        C_B = construction_B(C, h)
        @test length(C_B) == 7 - s
        @test dimension(C_B) >= 4 - s # typically drops by s-1
        
        # Construction B2: Needs binary code, dual codeword h, and 2j + 1 < wt(h)
        # wt(h) of a Hamming parity row is 4. Thus 2(1) + 1 = 3 < 4 is valid.
        C_B2 = construction_B2(C, h, 1)
        @test length(C_B2) == 7 - s
        
        # Construction X
        C3 = LinearCode(matrix(F, [1 0 1; 0 1 1])) # [3, 2]
        C2 = C # [7, 4, 3]
        C1 = subcode(C2, 2) # [7, 2]
        # C1 is a subcode of C2. C2.k (4) == C1.k (2) + C3.k (2).
        C_X = construction_X(C1, C2, C3)
        @test length(C_X) == 10
        @test dimension(C_X) == 4
        
        # Construction X3
        C4 = RepetitionCode(2, 3) # [3, 1]
        C5 = LinearCode(matrix(F, [1 0 0; 0 1 0])) # [3, 2]
        # Make a chain C_sub1 ⊂ C_sub2 ⊂ C
        C_sub2 = subcode(C, 3)
        C_sub1 = subcode(C_sub2, 1)
        # Check requirements: C.k(4) == C_sub2.k(3) + C4.k(1). C_sub2.k(3) == C_sub1.k(1) + C5.k(2).
        C_X3 = construction_X3(C_sub1, C_sub2, C, C4, C5)
        @test length(C_X3) == 13 # C1.n(7) + C4.n(3) + C5.n(3)
    end

    @testset "Operations: Products, Sums, Juxtaposition, Quotients" begin
        F = Oscar.Nemo.Native.GF(2)
        C1 = HammingCode(2, 3)       # [7, 4, 3]
        C2 = RepetitionCode(2, 3)    # [3, 1, 3]
        
        C_sum = C1 ⊕ C2
        @test length(C_sum) == 10
        @test dimension(C_sum) == 5
        
        C_prod = C1 × C2
        @test length(C_prod) == 21
        @test dimension(C_prod) == 4
        
        C_tensor = C1 ⊗ C2
        @test length(C_tensor) == 21
        @test dimension(C_tensor) == 21 - ((7-4) * (3-1)) # 21 - 6 = 15
        
        # Schur Product requires equal length
        C3 = RepetitionCode(2, 7)
        C_schur = Schur_product_code(C1, C3)
        @test length(C_schur) == 7
        
        # Juxtaposition requires equal dimension
        C1_sub = subcode(C1, 1)
        C_jux = juxtaposition(C1_sub, C2)
        @test length(C_jux) == 10
        @test dimension(C_jux) == 1
        
        # Quotient
        C_quo = C1 / C1_sub
        @test dimension(C_quo) == 3
        @test length(C_quo) == 7
    end

    @testset "Structural Modifications" begin
        F = Oscar.Nemo.Native.GF(2)
        C = HammingCode(2, 3) # [7, 4, 3]
        
        # Transpose
        C_trans = transpose(C)
        @test length(C_trans) == 3
        @test dimension(C_trans) == 3 - 3 # k_new = n_new - nrows(H_trans)
        
        # Permute
        C_perm = permute_code(C, [2, 1, 3, 4, 5, 6, 7])
        @test length(C_perm) == 7
        @test dimension(C_perm) == 4
        
        # Modifications
        C_ext = extend(C)
        @test length(C_ext) == 8
        @test dimension(C_ext) == 4
        
        C_punc = puncture(C, 1)
        @test length(C_punc) == 6
        @test dimension(C_punc) >= 3 # May drop
        
        C_exp = expurgate(C, 1)
        @test length(C_exp) == 7
        @test dimension(C_exp) == 3
        
        C_short = shorten(C, 1)
        @test length(C_short) == 6
        @test dimension(C_short) >= 3
        
        # Augment & Lengthen
        C_aug = augment(C, matrix(F, [1 1 1 1 1 1 1]))
        @test length(C_aug) == 7
        C_len = lengthen(C)
        @test length(C_len) == 8
    end

    # @testset "Subfield and Trace Codes" begin
    #     # Needs an extension field to test properly
    #     F4 = GF(2, 2, :a)
    #     a = gen(F4)
        
    #     # Simple code over GF(4)
    #     C4 = LinearCode(matrix(F4, [1 a a^2]))
        
    #     # Expand over GF(2)
    #     F2 = Oscar.Nemo.Native.GF(2)
    #     basis = [F4(1), a]
        
    #     C_exp = expanded_code(C4, F2, basis)
    #     @test length(C_exp) == 3
    #     @test dimension(C_exp) == 2
    #     @test field(C_exp) == F2
        
    #     # Subfield subcode / Trace code
    #     C_subf = subfield_subcode(C4, F2, basis)
    #     C_trace = trace_code(C4, F2, basis)
    #     @test field(C_subf) == F2
    #     @test field(C_trace) == F2
    # end

    @testset "Structural Parity (Even, Doubly, Triply)" begin
        F = Oscar.Nemo.Native.GF(2)
        # Extended Hamming [8, 4, 4] is a famous doubly-even code
        C = ExtendedHammingCode(3) 
        
        C_even = even_subcode(C)
        @test dimension(C_even) == 4 # already even
        
        C_de = doubly_even_subcode(C)
        @test dimension(C_de) == 4 # already doubly even
        
        C_te = triply_even_subcode(C)
        @test !ismissing(C_te)
        
        # Not all codes have doubly/triply even subcodes
        C_odd = RepetitionCode(2, 3) # [3, 1, 3] -> weight is 3 (odd)
        @test ismissing(even_subcode(C_odd))
    end
end