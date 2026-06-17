@testitem "Classical/Gabidulin.jl" begin
    using Oscar, CodingTheory
    
    F = GF(2)
    E = GF(2^4, :α)
    n, k, s = 4, 2, 1
    
    C = RandomGabidulinCode(F, E, n, k, s)
    @test length(C) == n
    @test dimension(C) == k
    
    # Test dual property: (C^⊥)^⊥ = C
    C_dual = dual(C)
    @test dimension(C_dual) == n - k
    C_double_dual = dual(C_dual)
    @test dimension(C_double_dual) == k
end