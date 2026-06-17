@testitem "Classical/Gleason.jl" begin
    using Oscar, CodingTheory
    
    # Test Extremal Type II Code (Extended Hamming [8, 4, 4] is Type I, 
    # but let's check a known extremal Type II, e.g., [24, 12, 8] Golay code)
    # Note: Requires a constructor for GolayCode or manual setup
    # Here we test the bound logic:
    n = 24
    # Placeholder for a self-dual code object
    # C_golay = GolayCode() 
    # @test Gleason_bound(C_golay) == 8
    
    # Test Quaternary Hermitian Self-Dual
    # n=6, d=2+2=4?
    # Logic: 2 * fld(6, 6) + 2 = 4
    @test Gleason_bound(DummySelfDual(GF(4), 6)) == 4 # Assumes a helper or direct check
end
