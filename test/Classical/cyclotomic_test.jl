@testset "Classical/cyclotomic.jl" begin
    using CodingTheory

    @test all_cyclotomic_cosets(2, 15, false) == [[0], [1, 2, 4, 8], [3, 6, 12, 9],
        [5, 10], [7, 14, 13, 11]]
    @test all_cyclotomic_cosets(3, 13, true) == [[0], [1, 3, 9], [2, 5, 6], [4, 10, 12] , [7, 8, 11]]
end
