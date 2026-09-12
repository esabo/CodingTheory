@testitem "Quantum/homological_measurements.jl" begin
    using Oscar
    using CodingTheory

    @testset "Homological Measurements" begin
        n = 8
        M = zeros(Int, n, n)
        for i in 1:n - 1
            M[i, i] = M[i, i + 1] = 1
        end
        M[n, 1] = M[n, n] = 1
        @test Cheeger_constant(M) == 4 / n

        F = Oscar.Nemo.Native.GF(2)
        S = Q422()
        L = matrix(F, 1, 2 * length(S), [1, 1, 0, 0, 0, 0, 0, 0])
        @test is_logical(S, L)
        measured = homological_measurement(S, L; style=:Cohen, r=1)
        @test measured isa AbstractStabilizerCode
        @test length(measured) > length(S)
        @test dimension(measured) == dimension(S) - 1

        L_Z = matrix(F, 1, 2 * length(S), [0, 0, 0, 0, 1, 1, 0, 0])
        @test is_logical(S, L_Z)
        measured_Z = homological_measurement(S, L_Z; style=:Cohen, r=1)
        @test dimension(measured_Z) == dimension(S) - 1
    end
end
