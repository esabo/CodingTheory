@testitem "Classical/concatenation.jl" begin
    using Oscar, CodingTheory

    @testset "Standard Concatenation & Distance" begin
        F = Oscar.Nemo.Native.GF(2)
        Ham = matrix(F, 4, 8, [
                1 0 0 0 0 1 1 1;
                0 1 0 0 1 0 1 1;
                0 0 1 0 1 1 0 1;
                0 0 0 1 1 1 1 0])
        G_in = matrix(F, 4, 8, [
                1 0 0 0 0 1 0 0;
                0 1 0 0 0 0 1 0;
                0 0 1 0 1 0 0 1;
                0 0 0 1 1 1 1 0])
        
        # G_out should be the [8, 4, 4] extended binary Hamming code
        C_Ham = LinearCode(Ham)
        @test dimension(C_Ham) == 4
        @test minimum_distance(C_Ham)[1] == 4
        
        # G_in should be [8, 4, 2]
        C_in = LinearCode(G_in)
        @test dimension(C_in) == 4
        @test minimum_distance(C_in)[1] == 2
        
        C = concatenate(C_Ham, C_in)
        G_final = matrix(F, 4, 16, [
                1 0 0 0 0 1 0 0 0 1 1 1 0 1 0 1;
                0 1 0 0 0 0 1 0 1 0 1 1 0 0 1 1;
                0 0 1 0 1 0 0 1 1 1 0 1 1 0 0 0;
                0 0 0 1 1 1 1 0 1 1 1 0 1 1 1 1])
        
        @test generator_matrix(C) == G_final
        @test minimum_distance(C)[1] == 8

        C2 = concatenate(C_Ham, C_Ham)
        G_final2 = matrix(F, 4, 16, [
                1 0 0 0 0 1 1 1 0 1 1 1 1 0 0 0;
                0 1 0 0 1 0 1 1 1 0 1 1 0 1 0 0;
                0 0 1 0 1 1 0 1 1 1 0 1 0 0 1 0;
                0 0 0 1 1 1 1 0 1 1 1 0 0 0 0 1])
        @test generator_matrix(C2) == G_final2
        @test minimum_distance(C2)[1] == 16 
    end

    @testset "Advanced Concatenation API & Orthogonality" begin
        F = Oscar.Nemo.Native.GF(2)
        
        # Setup generic [4, 2, 2] code for identical concatenation
        G_42 = matrix(F, 2, 4, [1 0 1 0; 0 1 0 1])
        C_42 = LinearCode(G_42)
        
        # 1. Operator Aliases & Getters
        C_cat = C_42 ∘ C_42 
        
        @test concatenation_type(C_cat) == :same
        @test inner_code(C_cat) === C_42
        @test outer_code(C_cat) === C_42
        
        # 2. Parity Check Orthogonality
        G_cat = generator_matrix(C_cat)
        H_cat = parity_check_matrix(C_cat)
        @test iszero(G_cat * transpose(H_cat))
        
        # 3. Encoding to Codespace
        msg_out = matrix(F, 1, 2, [1, 1])
        encoded_word = encode(C_cat, msg_out)
        
        @test ncols(encoded_word) == 8 # n_in * (n_out / k_in) = 4 * (4 / 2) = 8
        @test iszero(H_cat * transpose(encoded_word))
    end

    @testset "Multilevel Concatenation" begin
        F = Oscar.Nemo.Native.GF(2)
        
        # 1. Inner codes: C_in_1 ⊆ C_in_2
        C_in_1 = RepetitionCode(2, 4) # [4, 1, 4]
        C_in_2 = SPCCode(2, 4)        # [4, 3, 2]
        @test C_in_1 ⊆ C_in_2
        
        # 2. Outer codes
        # Target n_out = 3
        # C_out_1 length must be n_out * k_in_1 = 3 * 1 = 3
        C_out_1 = RepetitionCode(2, 3) 
        
        # C_out_2 length must be n_out * (k_in_2 - k_in_1) = 3 * (3 - 1) = 6
        G_63 = matrix(F, 3, 6, [1 0 0 1 1 0; 0 1 0 0 1 1; 0 0 1 1 0 1])
        C_out_2 = LinearCode(G_63)
        
        outers = [C_out_1, C_out_2]
        inners = [C_in_1, C_in_2]
        
        # 3. Execute Multilevel Build
        C_multi = multilevel_concatenation(outers, inners)
        
        @test C_multi.n == 12 # n_in * n_out = 4 * 3
        @test C_multi.k == 4  # k_out_1 + k_out_2 = 1 + 3
        
        # 4. Validate Block Quotient Matrices & Orthogonality
        G_multi = generator_matrix(C_multi)
        H_multi = parity_check_matrix(C_multi)
        @test iszero(G_multi * transpose(H_multi))
        
        # Test getters
        @test inner_code(C_multi) === inners
        @test outer_code(C_multi) === outers
        @test concatenation_type(C_multi) == [:same, :same]
    end
end
