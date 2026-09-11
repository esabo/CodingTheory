@testitem "Quantum/bounds.jl" begin
    using CodingTheory, Oscar, JuMP, Tulip

    @testset "Quantum Singleton and MDS bounds" begin
        @test quantum_Singleton_bound(5, 1) == 3
        @test quantum_Singleton_bound(9, 1; r=4) == 3
        @test quantum_Singleton_bound(3, 1 // 2) == 2
        @test quantum_Singleton_bound(
            3, 1 // 2; r=1, field_degree=2) == 2
        @test_throws DomainError quantum_Singleton_bound(5, 0)
        @test_throws DomainError quantum_Singleton_bound(5, 3; r=3)

        five_qubit = FiveQubitCode()
        @test Singleton_bound(five_qubit) == 3
        @test is_quantum_MDS(five_qubit)
        @test is_MDS(five_qubit)
        @test_throws DomainError set_minimum_distance!(five_qubit, 4)

        bacon_shor = BaconShorCode(3)
        @test quantum_Singleton_bound(bacon_shor) == 3
        @test is_quantum_MDS(bacon_shor)
    end

    @testset "Pure quantum Hamming bound" begin
        @test quantum_Hamming_volume(5, 0, 2) == 1
        @test quantum_Hamming_volume(5, 1, 2) == 16
        @test satisfies_quantum_Hamming_bound(5, 1, 3, 2)
        @test !satisfies_quantum_Hamming_bound(5, 1, 5, 2)
        @test quantum_Hamming_bound(5, 1, 2) == 4
        @test quantum_Hamming_bound(3, 1 // 2, 4; r=1) == 2

        five_qubit = FiveQubitCode()
        @test_throws ArgumentError quantum_Hamming_bound(five_qubit)
        five_qubit.cache[:pure] = true
        @test quantum_Hamming_bound(five_qubit) == 4

        # The impure [[9,1,4,3]] Bacon--Shor code beats the pure bound.
        bacon_shor = BaconShorCode(3)
        @test_throws ArgumentError quantum_Hamming_bound(bacon_shor)
        @test quantum_Hamming_bound(bacon_shor; assume_pure=true) == 2
    end

    @testset "Additive stabilizer Gilbert--Varshamov benchmark" begin
        five_qubit = FiveQubitCode()
        bounds_before = (
            five_qubit.cache[:l_bound], five_qubit.cache[:u_bound],
            five_qubit.cache[:d])
        @test quantum_Gilbert_Varshamov_bound(5, 1, 2) == 2
        @test quantum_Gilbert_Varshamov_bound(five_qubit) == 2
        @test (
            five_qubit.cache[:l_bound], five_qubit.cache[:u_bound],
            five_qubit.cache[:d]) == bounds_before

        # q=4 and k=1/2 still describe an integral protected space.
        @test quantum_Gilbert_Varshamov_bound(3, 1 // 2, 4) >= 1
        @test quantum_Gilbert_Varshamov_bound(
            3, 1 // 2, 4; r=1) >= 1
        @test quantum_Gilbert_Varshamov_exists(
            7, 1, 2, 2; variant=:linear)
        @test quantum_Gilbert_Varshamov_exists(
            7, 3, 2, 2; variant=:pure_feng_ma)
        @test_throws DomainError quantum_Gilbert_Varshamov_bound(
            3, 1 // 3, 4)
    end

    @testset "Construction-specific distance metadata" begin
        simplex = CodingTheory.dual(HammingCode(2, 3))
        hamming = HammingCode(2, 3)
        css = CSSCode(hamming, simplex)
        @test minimum_distance_lower_bound(css) == 1
        @test X_minimum_distance_lower_bound(css) == 1
        @test Z_minimum_distance_lower_bound(css) == 3
        @test minimum_distance_upper_bound(css) == 4
    end

    @testset "Low-check-weight bounds" begin
        @test quantum_stabilizer_generator_weight_lower_bound(100, 75) == 8
        @test quantum_stabilizer_generator_weight_lower_bound(
            12, 7; min_distance=2) == 6
        @test quantum_check_weight_dimension_bound(100, 8) == 75
        @test quantum_low_weight_stabilizer_distance_bound(100, 20, 2) == 1
        @test quantum_low_weight_stabilizer_distance_bound(100, 20, 3) == 2
        @test quantum_low_weight_stabilizer_distance_bound(100, 26, 3) == 1
        @test ismissing(
            quantum_low_weight_stabilizer_distance_bound(100, 20, 4))
        @test quantum_CSS_subsystem_weight_two_distance_bound(100, 4) == 10
        @test satisfies_quantum_CSS_subsystem_weight_two_bounds(
            100, 4, 10, 10)
        @test !satisfies_quantum_CSS_subsystem_weight_two_bounds(
            100, 4, 11, 10)
        @test quantum_stabilizer_check_weight_existence_bound(40, 10, 2) == 3
        @test quantum_stabilizer_check_weight_existence_bound(90, 10, 3) == 4
        @test ismissing(
            quantum_stabilizer_check_weight_existence_bound(39, 10, 2))
        @test quantum_stabilizer_group_average_weight(5) == 15 // 4
        @test quantum_stabilizer_group_total_weight(5, 1) == 60
    end

    @testset "Arbitrary-precision combinatorics" begin
        volume = quantum_Hamming_volume(200, 20, 2)
        @test volume isa BigInt
        @test volume > BigInt(10)^23
        @test quantum_check_weight_dimension_bound(
            BigInt(10)^30, BigInt(10)^15) isa BigInt

        # This is K_40(0; 80, 2), exercising coefficients above 10^23.
        binary_K = quantum_Krawtchouk_matrix(80; alphabet_size=2)
        @test binary_K[41, 1] == binomial(BigInt(80), 40)
        @test binary_K[41, 1] > BigInt(10)^23
        binary_K[1, 1] = 2
        @test quantum_Krawtchouk_matrix(
            80; alphabet_size=2)[1, 1] == 1
    end

    @testset "Arbitrary-precision weight-enumerator LPs" begin
        stabilizer_result =
            quantum_weight_enumerator_LP(5, 1, 3; precisions=(128,))
        @test stabilizer_result.status == :feasible
        @test stabilizer_result.enumerator !== nothing

        hooked_result = quantum_weight_enumerator_LP(
            5, 1, 3;
            precisions=(128,),
            model_hook=(model, A) -> JuMP.@constraint(model, A[5] == 15),
        )
        @test hooked_result.status == :feasible

        css_result = quantum_CSS_weight_enumerator_LP(
            7, 4, 4, 3; check_weight=4, precisions=(128,))
        @test css_result.status == :feasible
        @test css_result.enumerator !== nothing

        # Wei et al.'s finite example W_opt(12, 7, 2) >= 6.
        generator_bound = quantum_stabilizer_generator_weight_LP_bound(
            12, 7, 2;
            max_weight=6, precisions=(128,))
        @test generator_bound.lower_bound == 6

        raw_one = Dict((10, 3, 4) => 2, (11, 3, 4) => 3)
        raw_two = Dict((9, 3, 4) => 1, (10, 3, 4) => 2)
        closed = quantum_check_weight_LP_postprocess(raw_one, raw_two)
        @test closed[(10, 3, 4)] == 2
    end
end
