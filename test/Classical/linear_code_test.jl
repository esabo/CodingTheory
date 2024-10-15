@testitem "Classical/linear_code.jl" begin
    using Oscar, CodingTheory

    @testset "Linear Code" begin
        F = Oscar.Nemo.Native.GF(2)
        G = matrix(F, [1 0 0 0 0 1 1;
                   0 1 0 0 1 0 1;
                   0 0 1 0 1 1 0;
                   0 0 0 1 1 1 1]);
        C = LinearCode(G);
        @test field(C) == F
        @test length(C) == 7
        @test rank(G) == dimension(C)
        @test cardinality(C) == BigInt(2)^4
        @test CodingTheory.dimension(C) == 4
        @test rate(C) == 4 / 7
        @test minimum_distance(C) == 3
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
        cw = matrix(F, [1 0 0 0 0 1 1]);
        @test encode(C, C.G[:, 1:1]) == cw
        # these v's are C.G[:, 1], just testing different formats
        v = [1, 0, 0, 0];
        @test encode(C, v) == cw
        v2 = [1; 0; 0; 0];
        @test encode(C, v2) == cw
        # this vector is the first row of the generator matrix and should
        # therefore have zero syndrome
        v = [1, 0, 0, 0, 0, 1, 1];
        @test iszero(syndrome(C, v))
        v = [1; 0; 0; 0; 0; 1; 1];
        @test iszero(syndrome(C, v))
        @test is_overcomplete(ZeroCode(5))
        @test is_overcomplete(IdentityCode(5), :H)
        @test !is_overcomplete(HammingCode(2,3))
        @test !is_overcomplete(HammingCode(2,3), :H)

        # lower rank test
        G_and_G = vcat(G, G);
        C_G_and_G = LinearCode(G);
        @test rank(G_and_G) == dimension(C_G_and_G)
        @test G == generator_matrix(C_G_and_G)

        # information set tests:
        F = Oscar.Nemo.Native.GF(2)

        # a code with G in standard form 
        Gstd = matrix(F, [1 0; 0 1])
        Cstd = LinearCode(Gstd);
        pivs = information_set(Cstd)
        @test pivs == [1,2]

        # codes with nontrivial pivots:
        G = matrix(F, [1 1 0 0 0; 0 0 1 1 1])
        C = LinearCode(G);
        pivs = information_set(C)
        @test pivs == [1,3]

        C_ham = HammingCode(2, 3)
        pivs = information_set(C_ham)
        expected_info_mat = matrix(F, [1 0 0 0; 1 1 0 0; 0 1 1 0; 1 1 0 1])
        @test pivs == [1, 2, 3, 4]

        C = deepcopy(C_ham)
        C.G[:,1] = C_ham.G[:, 5]
        C.G[:,5] = C_ham.G[:, 1]
        pivs = information_set(C)
        @test pivs == [1, 2, 3, 5]
    end

    @testset "Puncturing examples" begin
        F = Oscar.Nemo.Native.GF(2)
        # puncturing examples from Huffman/Pless
        G = matrix(F, [1 1 0 0 0; 0 0 1 1 1])
        C = LinearCode(G)
        @test generator_matrix(puncture(C, [1])) == matrix(F, [1 0 0 0; 0 1 1 1])
        @test generator_matrix(puncture(C, [5])) == matrix(F, [1 1 0 0; 0 0 1 1])
        G = matrix(F, [1 0 0 0; 0 1 1 1])
        C = LinearCode(G)
        @test generator_matrix(puncture(C, [1])) == matrix(F, [1 1 1])
        @test generator_matrix(puncture(C, [4])) == matrix(F, [1 0 0; 0 1 1])
    end

    @testset "Extending and Shortening Property" begin
        F = Oscar.Nemo.Native.GF(2)
        # extending examples from Huffman/Pless
        C = TetraCode()
        exC = extend(C)
        @test generator_matrix(exC) == matrix(field(C), [1 0 1 1 0; 0 1 1 -1 -1])
        @test parity_check_matrix(exC) == matrix(field(C), [1 1 1 1 1; -1 -1 1 0 0; -1 1 0 1 0])
        G = matrix(F, [1 1 0 0 1; 0 0 1 1 0])
        C = LinearCode(G)
        @test generator_matrix(extend(puncture(C, [5]))) == matrix(F, [1 1 0 0 0; 0 0 1 1 0])
        G = matrix(F, [1 0 0 1 1 1; 0 1 0 1 1 1; 0 0 1 1 1 1])
        C = LinearCode(G)
        @test generator_matrix(puncture(C, [5, 6])) == matrix(F, [1 0 0 1; 0 1 0 1; 0 0 1 1])
        # shortening examples from Huffman/Pless
        shC = shorten(C, [5, 6])
        shCtest = LinearCode(matrix(F, [1 0 1 0; 0 1 1 0]))
        @test are_equivalent(shC, shCtest)

        C = Hexacode()
        D = Hermitian_dual(C)
        @test are_equivalent(C, D)
        # verify our definition of Hermitian_dual is equivalent:
        @test are_equivalent(D, LinearCode(Hermitian_conjugate_matrix(generator_matrix(dual(C)))))

        C = HammingCode(2, 3)
        C2 = LinearCode(words(C))
        @test are_equivalent(C, C2)
    end

        # missing so far in tests:
        # expurgate, augment, lengthen, uuplusv, subcode,
        # juxtaposition, constructionX, constructionX3, upluswvpluswuplusvplusw,
        # expandedcode, entrywiseproductcode, evensubcode
    @testset "Subfield Subcode" begin
        # subfield subcode example from Huffman/Pless
        K = GF(2, 2, :ω);
        ω = gen(K)
        C = Hexacode();
        dual_basis = [ω^2, K(1)]; # dual?
        #CF2 = subfield_subcode(C, F, dual_basis)
        #@test are_equivalent(CF2, RepetitionCode(2, 6))
    end

    @testset "Permutations of LinearCode" begin
        # test permutations
        C = HammingCode(2, 3)
        #S7 = SymmetricGroup(7);
        σ = [3, 2, 1, 4, 5, 6, 7] # this is the permutation (1, 3)
        C1 = permute_code(C, σ)
        C2 = permute_code(C, perm(σ))
        C3 = permute_code(C, Perm(σ))
        #C4 = permute_code(C, S7(σ))
        #@test C1.G == C2.G == C3.G == C4.G == C.G[:, σ]
        C1 = permute_code(C, [2,1,3,4,5,6,7])
        C2 = permute_code(C, [1,4,3,2,5,6,7])
        flag, P = are_permutation_equivalent(C1, C2)
        @test flag
        @test are_equivalent(permute_code(C1, P), C2)

        # Huffman and Pless Ex 1.6.1
        F = Oscar.Nemo.Native.GF(2)
        G1 = matrix(F, [1 1 0 0 0 0; 0 0 1 1 0 0; 0 0 0 0 1 1])
        print(typeof(G1))
        nc=ncols(G1)
        C1 = LinearCode(G1)
        G2 = matrix(F, [1 0 0 0 0 1; 0 0 1 1 0 0; 0 1 0 0 1 0])
        C2 = LinearCode(G2)
        G3 = matrix(F, [1 1 0 0 0 0; 1 0 1 0 0 0; 1 1 1 1 1 1])
        C3 = LinearCode(G3)
        LinearCode(words(C1))
        bool1, permutation1 = CodingTheory.are_perm_equivalent_exhaustive_search(C1, C2)
        @test bool1
        @test !(permutation1 === missing)
        bool2, permutation2 = CodingTheory.are_perm_equivalent_exhaustive_search(C1, C3)
        @test !bool2
        @test permutation2 === missing
 
    end
        # "On the Schur Product of Vector Spaces over Finite Fields"
        # Christiaan Koster
        # Lemma 14: If C is cyclic and dim(C) > (1/2)(n + 1), then C * C = F^n

        # simplex code itself has dimension k(k + 1)/2
        #

       

end
