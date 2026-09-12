@testitem "Classical/Gleason.jl" begin
    using Oscar, CodingTheory
    
    @testset "Binary Extremal: Extended Golay" begin
        # Extended Binary Golay Code: [24, 12, 8]
        # This is a Type II (Doubly-Even, Self-Dual) binary code
        C_ext_bin = ExtendedGolayCode(2)
        
        # The Gleason bound for n = 24 is 4 * floor(24/24) + 4 = 8
        @test Gleason_bound(C_ext_bin) == 8
        
        # Since d = 8, it achieves the bound
        @test is_extremal(C_ext_bin)
    end
    
    @testset "Ternary Extremal: Extended Golay" begin
        # Extended Ternary Golay Code: [12, 6, 6]
        # This is a Ternary Self-Dual code
        C_ext_tern = ExtendedGolayCode(3)
        
        # The Gleason bound for n = 12 is 3 * floor(12/12) + 3 = 6
        @test Gleason_bound(C_ext_tern) == 6
        
        # Since d = 6, it achieves the bound
        @test is_extremal(C_ext_tern)
    end
    
    @testset "Non-Self-Dual Rejection" begin
        # The standard Golay [23, 12, 7] is NOT self-dual
        C_bin = GolayCode(2)
        
        # The Gleason bound should return missing for non-self-dual codes
        @test ismissing(Gleason_bound(C_bin))
        @test is_extremal(C_bin) == false
    end
end
