@testitem "LDPC/decoders.jl" begin
    using Oscar, CodingTheory, JuMP, GLPK

    @testset "LP Decoder Initialization" begin
        # Standard Hamming(7,4) parity-check matrix
        # We use a standard Julia matrix since LP_decoder_LDPC natively takes AbstractMatrix{<:Number}
        H = [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ]
        
        # 1. Initialize from raw matrix
        model = CodingTheory._init_LP_decoder_LDPC(H)
        @test typeof(model) == JuMP.Model
        
        # Ensure the variable `f` for the bit values was successfully registered
        dict = object_dictionary(model)
        @test haskey(dict, :f)
        @test length(dict[:f]) == 7 # n = 7
        
        # 2. Initialize from an Oscar LinearCode
        F = Oscar.Nemo.Native.GF(2)
        H_oscar = matrix(F, H)
        C = LinearCode(H_oscar)
        
        model_code = CodingTheory._init_LP_decoder_LDPC(C)
        @test typeof(model_code) == JuMP.Model
        @test length(model_code[:f]) == 7
    end

    @testset "LP Decoding (Binary Symmetric Channel)" begin
        H = [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ]
        
        # Scenario: All-zeros codeword transmitted.
        # Received vector has an error at index 1.
        received = [1, 0, 0, 0, 0, 0, 0]
        
        Ch_bsc = BinarySymmetricChannel(0.1) # 10% crossover probability
        
        f_out = CodingTheory.LP_decoder_LDPC(H, received, Ch_bsc)
        
        # The LP solver should correct the error and pull the fractional bits back to exactly 0.0
        @test length(f_out) == 7
        @test all(x -> isapprox(x, 0.0, atol=1e-5), f_out)
        
        # Scenario: Error at index 4
        received_2 = [0, 0, 0, 1, 0, 0, 0]
        f_out_2 = CodingTheory.LP_decoder_LDPC(H, received_2, Ch_bsc)
        @test all(x -> isapprox(x, 0.0, atol=1e-5), f_out_2)
    end

    @testset "LP Decoding (Binary Erasure Channel)" begin
        H = [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ]
        
        # Scenario: All-zeros codeword transmitted.
        # Index 1 and 2 are erased (-1 represents an erasure in the classical channel mapper)
        received = [-1, -1, 0, 0, 0, 0, 0]
        
        Ch_bec = BinaryErasureChannel(0.3)
        
        f_out = CodingTheory.LP_decoder_LDPC(H, received, Ch_bec)
        
        # The LP solver should correctly infer that the erased bits must be 0
        @test length(f_out) == 7
        @test all(x -> isapprox(x, 0.0, atol=1e-5), f_out)
    end
    
    @testset "LP Decoding (All-Ones Codeword Recovery)" begin
        H = [
            1 1 0 1 1 0 0;
            1 0 1 1 0 1 0;
            0 1 1 1 0 0 1
        ]
        
        # The all-ones vector [1, 1, 1, 1, 1, 1, 1] is a valid codeword in Hamming(7,4).
        # We transmit it, but index 1 is flipped to 0.
        received = [0, 1, 1, 1, 1, 1, 1]
        
        Ch_bsc = BinarySymmetricChannel(0.1)
        f_out = CodingTheory.LP_decoder_LDPC(H, received, Ch_bsc)
        
        # The LP solver should correct index 1 back to 1.0
        @test all(x -> isapprox(x, 1.0, atol=1e-5), f_out)
    end
end