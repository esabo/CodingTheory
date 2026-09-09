@testitem "Classical/linear_code.jl" begin
    using Oscar, CodingTheory, Random

    @testset "Basic Linear Code Properties" begin
        F = Oscar.Nemo.Native.GF(2)
        G = matrix(F, [1 0 0 0 0 1 1;
                       0 1 0 0 1 0 1;
                       0 0 1 0 1 1 0;
                       0 0 0 1 1 1 1])
        C = LinearCode(G)
        
        @test field(C) == F
        @test length(C) == 7
        @test rank(G) == dimension(C)
        @test cardinality(C) == BigInt(2)^4
        @test CodingTheory.dimension(C) == 4
        @test rate(C) == 4 / 7
        
        @test minimum_distance(C)[1] == 3
        @test number_correctable_errors(C) == 1
        @test G == generator_matrix(C)
    
        H = parity_check_matrix(C)
        @test iszero(G * transpose(H))
        @test iszero(H * transpose(G))
        @test C ⊆ C
        
        D = dual(C)
        @test !(C ⊆ D)
        @test !is_subcode(C, D)
        @test !are_equivalent(C, D)
        @test are_equivalent(C, C)
        @test !is_self_dual(C)
        @test !is_self_orthogonal(C)
        
        cw = matrix(F, [1 0 0 0 0 1 1])
        @test encode(C, generator_matrix(C)[:, 1:1]) == cw
        
        v = [1, 0, 0, 0]
        @test encode(C, v) == cw
        v2 = [1; 0; 0; 0]
        @test encode(C, v2) == cw
        
        v_cw = [1, 0, 0, 0, 0, 1, 1]
        @test iszero(syndrome(C, v_cw))
        v_cw2 = [1; 0; 0; 0; 0; 1; 1]
        @test iszero(syndrome(C, v_cw2))
        
        @test !is_overcomplete(HammingCode(2, 3))
        @test !is_overcomplete(HammingCode(2, 3), :H)

        # Lower rank / redundant rows test
        # The constructor explicitly preserves the user's input topology!
        G_and_G = vcat(G, G)
        C_G_and_G = LinearCode(G_and_G)
        @test dimension(C_G_and_G) == rank(G)
        @test nrows(generator_matrix(C_G_and_G)) == 8 # Retained all rows!

        # Information Sets
        Gstd = matrix(F, [1 0; 0 1])
        Cstd = LinearCode(Gstd)
        @test information_set(Cstd) == [1, 2]

        G2 = matrix(F, [1 1 0 0 0; 0 0 1 1 1])
        C2 = LinearCode(G2)
        @test information_set(C2) == [1, 3]

        C_ham = HammingCode(2, 3)
        @test information_set(C_ham) == [1, 2, 3, 4]

        C3 = deepcopy(C_ham)
        G3 = generator_matrix(C3)
        temp_col = G3[:, 1]
        G3[:, 1] = G3[:, 5]
        G3[:, 5] = temp_col
        C3.cache[:G] = G3 
        @test information_set(C3) == [1, 2, 3, 5]
    end

    @testset "Constructors & Orthogonality Verification" begin
        F = GF(2)
        G = matrix(F, [1 0 1; 0 1 1])
        H = matrix(F, [1 1 1])
        
        C_GH = LinearCode(G, H, check_orthogonality=true)
        @test dimension(C_GH) == 2
        @test length(C_GH) == 3
        
        H_bad = matrix(F, [1 0 0])
        @test_throws ArgumentError LinearCode(G, H_bad, check_orthogonality=true)
        
        Gs = [[1, 0, 1], [0, 1, 1]]
        C_vecs = LinearCode(Gs, 2)
        @test dimension(C_vecs) == 2
        
        G1 = matrix(F, [1 0 1])
        G2 = matrix(F, [0 1 1])
        C_mats = LinearCode([G1, G2])
        @test dimension(C_mats) == 2
    end

    @testset "Bounds and Setters" begin
        F = GF(2)
        G = matrix(F, [1 0 1; 0 1 1])
        C = LinearCode(G)
        
        set_distance_upper_bound!(C, 2)
        @test minimum_distance_upper_bound(C) == 2
        
        set_distance_lower_bound!(C, 2)
        @test C.d == 2
        
        C_F4 = change_field(C, GF(4, :ω))
        @test Int(order(field(C_F4))) == 4
    end

    @testset "Duals, Hulls, and LCD Properties" begin
        F = GF(2)
        G = matrix(F, [1 0 1; 0 1 1])
        C_mds = LinearCode(G)
        @test is_MDS(C_mds)
        
        C_lcd = LinearCode(matrix(F, [1 0 0; 0 1 0]))
        hull_C, dim_hull = hull(C_lcd)
        @test dim_hull == 0
        @test is_LCD(C_lcd)
        
        F4 = GF(4, :ω)
        G4 = matrix(F4, [1 0 1; 0 1 1])
        C4 = LinearCode(G4)
        _, dim_H_hull = Hermitian_hull(C4)
        @test dim_H_hull >= 0
        @test is_Hermitian_LCD(C4) == (dim_H_hull == 0)
        
        C_gal = l_Galois_dual(C4, 1)
        @test dimension(C_gal) == 1
        _, dim_gal_hull = l_Galois_hull(C4, 1)
        @test is_l_Galois_LCD(C4, 1) == (dim_gal_hull == 0)
    end

    @testset "Even-ness Properties" begin
        F = GF(2)
        G_RM = matrix(F, [
            1 1 1 1 1 1 1 1;
            0 1 0 1 0 1 0 1;
            0 0 1 1 0 0 1 1;
            0 0 0 0 1 1 1 1
        ])
        C_RM = LinearCode(G_RM)
        # Explicit scope for custom even-ness methods
        @test CodingTheory.is_even(C_RM)
        @test CodingTheory.is_doubly_even(C_RM)
    end

    @testset "Puncturing, Extending, and Shortening" begin
        F = Oscar.Nemo.Native.GF(2)
        G = matrix(F, [1 1 0 0 0; 0 0 1 1 1])
        C = LinearCode(G)
        @test generator_matrix(puncture(C, [1])) == matrix(F, [1 0 0 0; 0 1 1 1])
        @test generator_matrix(puncture(C, [5])) == matrix(F, [1 1 0 0; 0 0 1 1])
        
        C_tetra = TetraCode()
        exC = extend(C_tetra)
        @test generator_matrix(exC) == matrix(field(C_tetra), [1 0 1 1 0; 0 1 1 -1 -1])
        
        shC = shorten(C, [4, 5])
        @test length(shC) == 3
    end

    # @testset "Permutations of LinearCode" begin
    #     C = HammingCode(2, 3)
    #     C1 = permute_code(C, [2,1,3,4,5,6,7])
    #     # P = (1, 2) permutation
    #     P = matrix(GF(2), [0 1 0 0 0 0 0; 
    #                        1 0 0 0 0 0 0; 
    #                        0 0 1 0 0 0 0;
    #                        0 0 0 1 0 0 0;
    #                        0 0 0 0 1 0 0;
    #                        0 0 0 0 0 1 0;
    #                        0 0 0 0 0 0 1])
    #     @test are_equivalent(permute_code(C, P), C1)
    # end

    @testset "Random Linear Code Functions" begin
        C_ham = HammingCode(2, 3) 
        pivs = random_information_set(C_ham)
        mat = generator_matrix(C_ham)[:, pivs]
        @test nrows(mat) == 4
        @test ncols(mat) == 4
        @test rank(mat) == 4

        C_ham4 = HammingCode(2, 4) 
        for i in 1:5
            rng = Random.seed!(i)
            infoset = random_information_set(C_ham4, rng = rng) 
            @test det(generator_matrix(C_ham4)[:, infoset]) != 0 
        end

        p = 2; n = 7; k = 4
        rng = Random.seed!(0)
        C_rand = random_linear_code(p, n, k, rng = rng)
        @test C_rand.n == n
        @test C_rand.k == k
    
        rng = Random.seed!(0)
        C_rand2 = random_linear_code(p^2, n, k, rng = rng)
        @test C_rand2.n == n
        
        rng_from_field = Random.seed!(0)
        C_rand3 = random_linear_code(GF(p, 2, :x), n, k, rng = rng_from_field)
        @test generator_matrix(C_rand2) == generator_matrix(C_rand3)
    end
    
    @testset "REPL Display Coverage (Silent)" begin
        io = IOBuffer()
        C = HammingCode(2, 3)
        show(io, C)
        show(io, MIME("text/plain"), C)
        
        # Fixed: raise gen(E) to powers, not the field E itself!
        E = GF(8, :α)
        α = gen(E)
        v = [E(1) for _ in 1:7]
        γ = [α^i for i in 0:6]
        
        C_GRS = GeneralizedReedSolomonCode(E, 7, 3, 5, 5, 5, v, v, γ, Dict{Symbol,Any}())
        show(io, C_GRS)
        @test true 
    end

    @testset "Sparse parity-check construction" begin
        using SparseArrays

        # Hamming H = [P^T | I] for the generator used above, plus a repeated row.
        H = sparse([1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4],
                   [2, 3, 4, 5, 1, 3, 4, 6, 1, 2, 4, 7, 2, 3, 4, 5],
                   ones(Int, 16), 4, 7)
        C = LinearCode(H, true)
        @test C.n == 7
        @test C.k == 4
        @test size(parity_check_matrix(C), 1) == 4
        G = generator_matrix(C)
        @test size(G, 1) == 4
        @test iszero(G * transpose(matrix(C.F, Array(H))))
    end
end
