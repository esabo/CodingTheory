@testitem "Classical/designs.jl" begin
    using Oscar, CodingTheory
    
    # Test Assmus-Mattson on the Hamming [7, 4, 3] code
    # The Hamming code is known to form a 2-design
    F = GF(2)
    H = matrix(F, [[1,0,1,1,1,0,0], [0,1,0,1,1,1,0], [0,0,1,0,1,1,1]])
    C = LinearCode(H)
    @test is_design_holder(C, 2) == true
    @test design_strength(C) >= 2

    # Test Reed-Muller 3-design property
    # RM(1, 3) is a [8, 4, 4] code
    RM = ReedMullerCode(1, 3)
    @test is_design_holder(RM, 3) == true
    @test minimum_weight_blocks(RM) == 14 # Number of weight-4 codewords
end