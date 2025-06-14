@testitem "LDPC/MP_decoders.jl" begin
    using Oscar, CodingTheory

    @testset "Message Passing" begin
        # not the most robust test, but I find that if something doesn't work everything fails
        F = Oscar.Nemo.Native.GF(2);
        H = matrix(F, [1 1 0 1 1 0 0; 1 0 1 1 0 1 0; 0 1 1 1 0 0 1]);
        v = matrix(F, 7, 1, [1, 1, 0, 0, 0, 0, 0]);
        correct_v = UInt8.([1, 1, 1, 0, 0, 0, 0]);
        correct_e = UInt8.([0, 0, 1, 0, 0, 0, 0]);
        syn = H * v;
        p = 1/7;
        nm = BSC(p);

        # basic cases
        flag, out, iter, _ = sum_product(H, v, nm);
        @test flag == true && out == correct_v
        flag, out, iter, _ = sum_product_box_plus(H, v, nm);
        @test flag == true && out == correct_v
        flag, out, iter, _ = sum_product_syndrome(H, syn, nm);
        @test flag == true && out == correct_e
        flag, out, iter, _ = min_sum(H, v, nm);
        @test flag == true && out == correct_v
        flag, out, iter, _ = min_sum_syndrome(H, syn, nm);
        @test flag == true && out == correct_e
        flag, out, iter, _ = min_sum_with_correction(H, v, nm);
        @test flag == true && out == correct_v
        flag, out, iter, _ = min_sum_with_correction_syndrome(H, syn, nm);
        @test flag == true && out == correct_e
        # flag, out, iter = Gallager_A(H, v);
        # @test flag == true && out == correct_v
        # flag, out, iter = Gallager_B(H, v);
        # @test flag == true && out == correct_v

        # all use the same init and loop functions so it suffices to test the options for a single function
        # some options
        flag, out, iter, _ = sum_product(H, v, nm, schedule = :parallel);
        @test flag == true && out == correct_v
        flag, out, iter, _ = sum_product(H, v, nm, schedule = :serial);
        @test flag == true && out == correct_v
        # flag, out, iter, _ = sum_product(H, v, nm, schedule = :layered);
        # @test flag == true && out == correct_v
        # flag, out, iter, _ = sum_product(H, v, nm, schedule = :layered, rand_sched = true);
        # @test flag == true && out == correct_v
        # flag, out, iter, _ = sum_product(H, v, nm, erasures = [rand(1:7)]);
        @test flag == true && out == correct_v
        flag, out, iter, _ = min_sum(H, v, nm, erasures = [rand(1:7)]);
        @test flag == true && out == correct_v
        # TODO this one fails for some reason
        flag, out, iter, _ = min_sum_with_correction(H, v, nm, erasures = [rand(1:7)]);
        @test_broken flag == true && out == correct_v
        # not particularly creative...
        temp = log((1 - p) / p);
        chn_inits = zeros(Float64, length(v));
        @inbounds for i = 1:nrows(v)
            iszero(v[i]) ? (chn_inits[i] = temp;) : (chn_inits[i] = -temp;)
        end
        flag, out, iter, _ = sum_product(H, v, nm, chn_inits = chn_inits);
        @test flag == true && out == correct_v
        attenuation = 0.6
        flag, out, iter, _ = min_sum(H, v, nm, attenuation = 0.6);
        @test flag == true && out == correct_v
        flag, out, iter, _ = min_sum_with_correction(H, v, nm, attenuation = 0.6);
        @test flag == true && out == correct_v

        # decimation
        # decimated_bits_values = [(1, base_ring(v)(1))];
        # flag, out, iter, _ = sum_product_decimation(H, v, nm, decimated_bits_values); flag
        # @test flag == true && out == correct_v
        # flag, out, iter, _ = min_sum_decimation(H, v, nm, decimated_bits_values);
        # @test flag == true && out == correct_v
        # flag, out, iter, _ = min_sum_correction_decimation(H, v, nm, decimated_bits_values);
        # @test flag == true && out == correct_v

        # other noise models
        # nm_BEC = BEC(p);
        # nm_G = BAWGNC(p);
    end
end
