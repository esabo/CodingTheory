@testitem "Quantum/homological_measurements.jl" begin
    using Oscar
    using CodingTheory

    @testset "Homological Measurements" begin
        n = rand(4:2:30)
        M = zeros(Int, n, n)
        for i in 1:n - 1
            M[i, i] = M[i, i + 1] = 1
        end
        M[n, 1] = M[n, n] = 1
        @test Cheeger_constant(M) == 4 / n
    end
end
