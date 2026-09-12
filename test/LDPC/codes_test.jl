@testitem "LDPC/codes.jl" begin
    using Oscar, CodingTheory

    @testset "Regular LDPC Code Construction" begin
        # example from Ryan & Lin
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(F, [
            1 1 1 1 0 0 0 0 0 0;
            1 0 0 0 1 1 1 0 0 0;
            0 1 0 0 1 0 0 1 1 0;
            0 0 1 0 0 1 0 1 0 1;
            0 0 0 1 0 0 1 0 1 1])
        C = LDPCCode(H)
        
        # In a generic setup, these bounds might be accessed via internal properties or methods
        @test CodingTheory.column_row_bounds(C) == (2, 4)
        @test is_regular(C)
        
        v_degs = unique!(sort!(CodingTheory.variable_degree_distribution(C)))
        c_degs = unique!(sort!(CodingTheory.check_degree_distribution(C)))
        
        @test v_degs == [2]
        @test c_degs == [4]
        
        R = parent(variable_degree_polynomial(C))
        x = gen(R)
        @test variable_degree_polynomial(C) == x
        @test check_degree_polynomial(C) == x^3
    end

    @testset "Irregular LDPC Code Construction" begin
        F = Oscar.Nemo.Native.GF(2)
        # Custom irregular parity-check matrix
        # Columns: four of degree 2, four of degree 1
        # Rows: one of degree 2, two of degree 3, one of degree 4
        # Total edges = 12
        H_irreg = matrix(F, [
            1 1 0 0 0 0 0 0;
            1 0 1 1 0 0 0 0;
            0 1 1 0 1 0 0 0;
            0 0 0 1 0 1 1 1
        ])
        C_irreg = LDPCCode(H_irreg)
        
        @test CodingTheory.column_row_bounds(C_irreg) == (2, 4)
        @test !is_regular(C_irreg)
        
        R = parent(variable_degree_polynomial(C_irreg))
        x = gen(R)
        
        # Edge-perspective variable degree polynomial λ(x):
        # 4 edges from deg-1 nodes (4/12 = 1/3) -> coefficient for x^0
        # 8 edges from deg-2 nodes (8/12 = 2/3) -> coefficient for x^1
        expected_v_poly = (1//3) + (2//3) * x
        @test variable_degree_polynomial(C_irreg) == expected_v_poly
        
        # Edge-perspective check degree polynomial ρ(x):
        # 2 edges from deg-2 nodes (2/12 = 1/6) -> coefficient for x^1
        # 6 edges from deg-3 nodes (6/12 = 1/2) -> coefficient for x^2
        # 4 edges from deg-4 nodes (4/12 = 1/3) -> coefficient for x^3
        expected_c_poly = (1//6)*x + (1//2)*x^2 + (1//3)*x^3
        @test check_degree_polynomial(C_irreg) == expected_c_poly
    end

    @testset "Degree-Zero Nodes" begin
        F = Oscar.Nemo.Native.GF(2)
        R, x = polynomial_ring(Oscar.Nemo.QQ, :x)

        # an unused column has degree zero and contributes nothing to λ(x)
        C = LDPCCode(matrix(F, [1 1 0]))
        @test variable_degree_polynomial(C) == one(R)
        @test check_degree_polynomial(C) == x
        @test CodingTheory.column_row_bounds(C) == (1, 2)
        @test !is_regular(C)

        # a code with no edges at all has zero degree polynomials
        C_empty = LDPCCode(zero_matrix(F, 2, 3))
        @test iszero(variable_degree_polynomial(C_empty))
        @test iszero(check_degree_polynomial(C_empty))
    end

    @testset "String Representations" begin
        F = Oscar.Nemo.Native.GF(2)
        H = matrix(F, [1 1 0; 0 1 1])
        C = LDPCCode(H)
        
        # Capture the output of the custom show method
        output = sprint(show, C)
        
        @test contains(output, "irregular") || contains(output, "regular")
        @test contains(output, "density")
        @test contains(output, "Variable degree polynomial:")
        @test contains(output, "Check degree polynomial:")
        # @test contains(output, "Parity-check matrix:")
    end
end