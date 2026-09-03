@testitem "LDPC/MP_decoders.jl" begin
    using Random, Oscar, CodingTheory

    # Rebuild a dense matrix from a 1-based CSR pair, so `csr_of` can be checked by
    # round trip rather than against a hand-transcribed index list.
    function dense_of_csr(row_ptr, col_ind, num_check, num_var)
        H = zeros(UInt8, num_check, num_var)
        for c in 1:num_check
            for i in row_ptr[c]:(row_ptr[c + 1] - 1)
                H[c, col_ind[i]] = 0x01
            end
        end
        return H
    end

    syndrome_of(H, e) = UInt8[reduce(⊻, H[c, :] .& e) for c in 1:size(H, 1)]

    # Array dispersion of the exponent matrix `(i - 1) * (j - 1) mod L` into circulant
    # permutation matrices: a sparse, 4-cycle-free, column-weight-`rows` QC-LDPC code.
    function qc_ldpc(L, rows, cols)
        H = zeros(UInt8, rows * L, cols * L)
        for i in 1:rows, j in 1:cols
            shift = mod((i - 1) * (j - 1), L)
            for k in 0:(L - 1)
                H[(i - 1) * L + k + 1, (j - 1) * L + mod(k + shift, L) + 1] = 0x01
            end
        end
        return H
    end

    # All-but-one aggregates of a box-plus fold, computed three different ways.
    seq_agg(msgs, skip) = foldl(CodingTheory.boxplus_exact,
                                [m for (i, m) in enumerate(msgs) if i != skip])

    function fb_agg(msgs, skip)
        deg = length(msgs)
        pre = ones(Float64, deg)
        acc = 1.0e12
        for i in 1:deg
            pre[i] = acc
            acc = CodingTheory.boxplus_exact(acc, msgs[i])
        end
        acc = 1.0e12
        out = 0.0
        for i in deg:-1:1
            i == skip && (out = CodingTheory.boxplus_exact(pre[i], acc))
            acc = CodingTheory.boxplus_exact(acc, msgs[i])
        end
        return out
    end

    prod_agg(msgs, skip) =
        2 * atanh(prod(tanh(m / 2) for (i, m) in enumerate(msgs) if i != skip))

    hamming74 = UInt8[1 0 1 0 1 0 1; 0 1 1 0 0 1 1; 0 0 0 1 1 1 1]

    @testset "csr_of" begin
        H = UInt8[1 0 1 0 0; 0 0 0 0 0; 0 1 1 0 0]
        row_ptr, col_ind = CodingTheory.csr_of(H)

        # 1-based, one entry per row boundary, indices ascending inside a row.
        @test row_ptr == [1, 3, 3, 5]
        @test col_ind == [1, 3, 2, 3]
        @test length(row_ptr) == size(H, 1) + 1
        @test row_ptr[end] == length(col_ind) + 1
        @test all(issorted(col_ind[row_ptr[c]:(row_ptr[c + 1] - 1)]) for c in 1:size(H, 1))

        # An all-zero row is an empty CSR row and an all-zero column simply never
        # appears; both round trip.
        @test dense_of_csr(row_ptr, col_ind, size(H)...) == H

        # Element type is irrelevant: only `iszero` is consulted.
        @test CodingTheory.csr_of(Float64.(H)) == (row_ptr, col_ind)
        @test CodingTheory.csr_of(Int.(H)) == (row_ptr, col_ind)

        # Flint-native form agrees with the plain-array one.
        F = Oscar.Nemo.Native.GF(2)
        @test CodingTheory.csr_of(matrix(F, Int.(H))) == (row_ptr, col_ind)

        # A denser, less symmetric case.
        rng = MersenneTwister(2718)
        H2 = UInt8.(rand(rng, Bool, 9, 14))
        rp2, ci2 = CodingTheory.csr_of(H2)
        @test dense_of_csr(rp2, ci2, size(H2)...) == H2
        @test length(ci2) == count(!iszero, H2)
    end

    @testset "Schedule construction" begin
        H = qc_ldpc(7, 3, 5)
        num_check, num_var = size(H)
        layers = CodingTheory.layered_schedule(H)

        # A partition: every check exactly once.
        @test sort(reduce(vcat, layers)) == collect(1:num_check)
        @test all(!isempty, layers)

        # Independence: no two checks in a layer share a variable.
        for layer in layers, i in eachindex(layer), j in eachindex(layer)
            i < j || continue
            @test !any(!iszero(H[layer[i], v]) && !iszero(H[layer[j], v]) for v in 1:num_var)
        end

        # Flat form agrees with the nested form, and is what the workspace stores.
        row_ptr, col_ind = CodingTheory.csr_of(H)
        layer_ptr, layer_checks = CodingTheory.layered_schedule(row_ptr, col_ind,
                                                                num_check, num_var; base = 1)
        @test [layer_checks[layer_ptr[l]:(layer_ptr[l + 1] - 1)]
               for l in 1:(length(layer_ptr) - 1)] == layers

        # 0-based input, as scipy hands it over.
        @test CodingTheory.layered_schedule(row_ptr .- 1, col_ind .- 1, num_check, num_var;
                                            base = 0) == (layer_ptr, layer_checks)

        # A degree-0 check shares no variable with anything, so it joins layer one.
        H0 = UInt8[1 0 1 0; 0 0 0 0; 0 1 1 0]
        @test CodingTheory.layered_schedule(H0) == [[1, 2], [3]]

        F = Oscar.Nemo.Native.GF(2)
        @test CodingTheory.layered_schedule(matrix(F, Int.(H0))) == [[1, 2], [3]]

        # Serial: one check per layer.
        ser_ptr, ser_checks = CodingTheory.serial_schedule(num_check)
        @test ser_ptr == collect(1:(num_check + 1))
        @test ser_checks == collect(1:num_check)
        @test all(ser_ptr[l + 1] - ser_ptr[l] == 1 for l in 1:num_check)

        # Balance, on both forms it accepts.
        @test CodingTheory.balance_of_layered_schedule(layer_ptr) ==
              CodingTheory.balance_of_layered_schedule(layers)
        @test CodingTheory.balance_of_layered_schedule(layer_ptr) ==
              maximum(length, layers) / minimum(length, layers)
        @test CodingTheory.balance_of_layered_schedule(layer_ptr) >= 1.0
        @test CodingTheory.balance_of_layered_schedule(ser_ptr) == 1.0
        @test CodingTheory.balance_of_layered_schedule([[1, 2], [3]]) == 2.0

        @test_throws ArgumentError CodingTheory.balance_of_layered_schedule(Int[1])
        @test_throws ArgumentError CodingTheory.balance_of_layered_schedule(Vector{Int}[])
        @test_throws ArgumentError CodingTheory.balance_of_layered_schedule([[1, 2], Int[]])
        @test_throws ArgumentError CodingTheory.balance_of_layered_schedule([1, 3, 3, 5])
        @test_throws ArgumentError CodingTheory.layered_schedule([1], [1], 0, 4)
    end

    @testset "Box-plus kernels" begin
        # Hand-computed: base + log1p(exp(-|x + y|)) - log1p(exp(-|x - y|)).
        @test CodingTheory.boxplus_exact(1.0, 1.0) ≈ 1 + log1p(exp(-2.0)) - log(2.0)
        @test CodingTheory.boxplus_exact(1.0, 1.0) ≈ 0.4337808304830272
        @test CodingTheory.boxplus_exact(3.0, 4.0) ≈ 2.6876497789355516
        @test CodingTheory.boxplus_minsum(3.0, 4.0) == 3.0
        @test CodingTheory.boxplus_minsum_correction(3.0, 4.0) == 2.5

        # Sign handling: the sign is the product of the two signs.
        for (x, y) in ((2.0, 3.0), (-2.0, 3.0), (2.0, -3.0), (-2.0, -3.0))
            expected = sign(x) * sign(y)
            @test sign(CodingTheory.boxplus_exact(x, y)) == expected
            @test CodingTheory.boxplus_minsum(x, y) == expected * 2.0
            @test sign(CodingTheory.boxplus_minsum_correction(x, y)) == expected
        end
        @test CodingTheory.boxplus_exact(1.0, -1.0) == -CodingTheory.boxplus_exact(1.0, 1.0)
        @test CodingTheory.boxplus_minsum_correction(1.0, -1.0) == -0.5
        @test CodingTheory.boxplus_minsum_correction(1.0, 1.0) == 0.5

        # A zero argument annihilates: an uninformative input makes the pair
        # uninformative, for all three rules.
        for y in (-7.0, -0.25, 0.0, 0.25, 7.0)
            @test iszero(CodingTheory.boxplus_exact(0.0, y))
            @test iszero(CodingTheory.boxplus_exact(y, 0.0))
            @test iszero(CodingTheory.boxplus_minsum(0.0, y))
            @test iszero(CodingTheory.boxplus_minsum(y, 0.0))
            @test iszero(CodingTheory.boxplus_minsum_correction(0.0, y))
        end

        # Symmetry.
        rng = MersenneTwister(1234)
        for _ in 1:200
            x, y = 6 * randn(rng), 6 * randn(rng)
            @test CodingTheory.boxplus_exact(x, y) == CodingTheory.boxplus_exact(y, x)
            @test CodingTheory.boxplus_minsum(x, y) == CodingTheory.boxplus_minsum(y, x)
            @test CodingTheory.boxplus_minsum_correction(x, y) ==
                  CodingTheory.boxplus_minsum_correction(y, x)
        end

        # Agreement with the tanh form, over the range where that form is stable.
        rng = MersenneTwister(4321)
        for _ in 1:300
            x, y = 16 * (rand(rng) - 0.5), 16 * (rand(rng) - 0.5)
            @test CodingTheory.boxplus_exact(x, y) ≈ 2 * atanh(tanh(x / 2) * tanh(y / 2)) atol=1e-9
            # Min-sum is the magnitude the exact rule is bounded by.
            @test abs(CodingTheory.boxplus_exact(x, y)) <=
                  abs(CodingTheory.boxplus_minsum(x, y)) + 1e-12
            # The correction is exactly one of three offsets away from plain min-sum.
            delta = CodingTheory.boxplus_minsum_correction(x, y) -
                    CodingTheory.boxplus_minsum(x, y)
            @test minimum(abs(delta - c) for c in (0.0, 0.5, -0.5)) < 1e-12
        end

        # Large magnitudes: the tanh form has already saturated to Inf here, the
        # Jacobian-logarithm form has not.
        @test CodingTheory.boxplus_exact(50.0, 50.0) ≈ 50.0 - log(2.0)
        @test 2 * atanh(tanh(25.0) * tanh(25.0)) == Inf
        @test CodingTheory.boxplus_minsum(50.0, 50.0) == 50.0
        @test CodingTheory.boxplus_minsum_correction(50.0, 50.0) == 49.5
        @test isfinite(CodingTheory.boxplus_exact(1.0e300, 1.0e300))
        @test isfinite(CodingTheory.boxplus_exact(Inf, 1.0))
        @test CodingTheory.boxplus_exact(Inf, 1.0) == 1.0
        @test CodingTheory.boxplus_exact(-Inf, 1.0) == -1.0

        # `_LLR_ID` is the fold identity: it must return the other message untouched.
        for y in (-9.5, -1.0, 0.0, 0.5, 9.5)
            @test CodingTheory.boxplus_exact(CodingTheory._LLR_ID, y) == y
            @test CodingTheory.boxplus_minsum(CodingTheory._LLR_ID, y) == y
            @test CodingTheory.boxplus_minsum_correction(CodingTheory._LLR_ID, y) == y
        end
    end

    @testset "Check-node update equivalence" begin
        # The three routes to an all-but-one box-plus aggregate: the quadratic fold,
        # forward-backward accumulation, and the tanh product form.
        rng = MersenneTwister(97531)
        for _ in 1:100
            msgs = [3 * randn(rng) for _ in 1:7]
            for skip in eachindex(msgs)
                a = seq_agg(msgs, skip)
                @test a ≈ fb_agg(msgs, skip) atol=1e-9
                @test a ≈ prod_agg(msgs, skip) atol=1e-9
            end
        end

        # Same three routes as the decoder actually implements them, on a degree-6
        # check (at `_FB_MIN_DEGREE`, so forward-backward really runs) and on a
        # degree-4 one (below it, so the quadratic fold runs).
        for deg in (4, 6)
            H = zeros(UInt8, 1, deg)
            H[1, :] .= 0x01
            W = CodingTheory.init_soft_workspace(H)
            rng = MersenneTwister(24680)
            for _ in 1:100
                for e in eachindex(W.V2C)
                    W.V2C[e] = 4 * randn(rng)
                end
                # The min-sum family: the single pass is claimed bit-identical to
                # both folds, not merely close.
                for algo in (:min_sum, :normalized_min_sum, :offset_min_sum)
                    CodingTheory._check_update!(W, Val(algo), 1, 0.75, 0.5)
                    single_pass = copy(W.C2V)
                    CodingTheory._check_update!(Val(:forward_backward), W, Val(algo), 1,
                                                0.75, 0.5)
                    @test W.C2V == single_pass
                    CodingTheory._check_update!(Val(:sequential), W, Val(algo), 1, 0.75, 0.5)
                    @test W.C2V == single_pass
                end
                # Sum-product: the two folds and the product form agree numerically.
                CodingTheory._check_update!(Val(:forward_backward), W, Val(:sum_product), 1,
                                            0.75, 0.5)
                fb = copy(W.C2V)
                CodingTheory._check_update!(Val(:sequential), W, Val(:sum_product), 1,
                                            0.75, 0.5)
                @test W.C2V ≈ fb atol=1e-9
                CodingTheory._check_update!(Val(:product), W, Val(:sum_product_fast), 1,
                                            0.75, 0.5)
                @test W.C2V ≈ fb atol=1e-9
            end
        end

        # An exact zero among the inputs is the case the negative-count trick exists
        # for: a sign product would zero the aggregates the zero is excluded from.
        H = zeros(UInt8, 1, 6)
        H[1, :] .= 0x01
        W = CodingTheory.init_soft_workspace(H)
        W.V2C .= [0.0, -1.0, 2.0, 3.0, -4.0, 0.5]
        CodingTheory._check_update!(W, Val(:min_sum), 1, 0.75, 0.5)
        single_pass = copy(W.C2V)
        CodingTheory._check_update!(Val(:sequential), W, Val(:min_sum), 1, 0.75, 0.5)
        @test W.C2V == single_pass
        # Only the edge carrying the zero learns anything; every other aggregate is
        # annihilated by it.
        @test single_pass[1] == 0.5
        @test all(iszero, single_pass[2:end])

        # A degree-1 check determines its variable outright, in both syndromes.
        H1 = UInt8[1 0; 1 1]
        W1 = CodingTheory.init_soft_workspace(H1)
        CodingTheory.load_soft_channel!(W1, [1.0, 1.0]; syndrome = UInt8[0, 0])
        CodingTheory._check_update!(W1, Val(:sum_product), 1, 0.75, 0.5)
        @test W1.C2V[1] == CodingTheory._LLR_MAX
        CodingTheory.load_soft_channel!(W1, [1.0, 1.0]; syndrome = UInt8[1, 0])
        CodingTheory._check_update!(W1, Val(:sum_product), 1, 0.75, 0.5)
        @test W1.C2V[1] == -CodingTheory._LLR_MAX
    end

    @testset "Workspace construction and channel loading" begin
        H = qc_ldpc(7, 3, 5)
        num_check, num_var = size(H)

        W = CodingTheory.init_soft_workspace(H)
        @test W isa CodingTheory.SoftDecisionWorkspace{Float64}
        @test W.num_check == num_check
        @test W.num_var == num_var
        @test W.num_edges == count(!iszero, H)
        @test W.max_check_degree == maximum(sum(Int.(H), dims = 2))
        @test isempty(W.layer_ptr)          # flooding builds no partition
        @test W.chk_ptr == CodingTheory.csr_of(H)[1]
        @test W.edge_var == CodingTheory.csr_of(H)[2]
        # The variable-side gather list is the transpose of the same graph.
        for v in 1:num_var
            @test sort([W.edge_var[e] for e in W.var_edges[W.var_ptr[v]:(W.var_ptr[v + 1] - 1)]]) ==
                  fill(v, count(!iszero, H[:, v]))
        end

        Wl = CodingTheory.init_soft_workspace(H; schedule = :layered)
        @test !isempty(Wl.layer_ptr)
        @test sort(Wl.layer_checks) == collect(1:num_check)
        @test Wl.layer_ptr[end] == num_check + 1

        Ws = CodingTheory.init_soft_workspace(H; schedule = :serial)
        @test Ws.layer_ptr == collect(1:(num_check + 1))
        @test Ws.layer_checks == collect(1:num_check)

        # `:parallel` is a synonym for `:flooding`, `:semiserial` for `:layered`.
        @test isempty(CodingTheory.init_soft_workspace(H; schedule = :parallel).layer_ptr)
        @test CodingTheory.init_soft_workspace(H; schedule = :semiserial).layer_ptr ==
              Wl.layer_ptr

        # A caller-supplied partition overrides `schedule`.
        Wg = CodingTheory.init_soft_workspace(H; layer_ptr = [1, 4, num_check + 1],
                                              layer_checks = collect(num_check:-1:1))
        @test Wg.layer_ptr == [1, 4, num_check + 1]
        @test Wg.layer_checks == collect(num_check:-1:1)

        # 0-based compressed-sparse-row input, as scipy hands it over.
        row_ptr, col_ind = CodingTheory.csr_of(H)
        W0 = CodingTheory.init_soft_workspace(row_ptr .- 1, col_ind .- 1, num_check, num_var;
                                              base = 0, schedule = :layered)
        @test W0.chk_ptr == Wl.chk_ptr
        @test W0.edge_var == Wl.edge_var
        @test W0.layer_ptr == Wl.layer_ptr
        @test W0.layer_checks == Wl.layer_checks

        F = Oscar.Nemo.Native.GF(2)
        Wf = CodingTheory.init_soft_workspace(matrix(F, Int.(H)); schedule = :layered)
        @test Wf.chk_ptr == Wl.chk_ptr
        @test Wf.edge_var == Wl.edge_var
        @test Wf.layer_ptr == Wl.layer_ptr

        @test_throws ArgumentError CodingTheory.init_soft_workspace(H; schedule = :bogus)
        @test_throws ArgumentError CodingTheory.init_soft_workspace([1, 3], [1, 2], 3, 4;
                                                                    base = 1)
        @test_throws ArgumentError CodingTheory.init_soft_workspace([1, 3], [1, 9], 1, 4;
                                                                    base = 1)
        @test_throws ArgumentError CodingTheory.init_soft_workspace(Int[1], Int[], 0, 0)
        @test_throws ArgumentError CodingTheory.init_soft_workspace(
            H; layer_ptr = [1, num_check + 1], layer_checks = collect(1:(num_check - 1)))
        @test_throws ArgumentError CodingTheory.init_soft_workspace(
            H; layer_ptr = [1, 5, 5, num_check + 1], layer_checks = collect(1:num_check))
        @test_throws ArgumentError CodingTheory.init_soft_workspace(
            H; layer_ptr = [2, num_check + 1], layer_checks = collect(1:num_check))
        @test_throws ArgumentError CodingTheory.init_soft_workspace(
            H; layer_ptr = [1, num_check + 1], layer_checks = fill(1, num_check))

        # Channel loading.
        llr = collect(range(-2.0, 2.0; length = num_var))
        syn = zeros(UInt8, num_check)
        syn[2] = 0x01
        returned = CodingTheory.load_soft_channel!(W, llr; syndrome = syn, erasures = [1, 3],
                                                   decimated_bits = [5, 6],
                                                   decimated_values = [0, 1])
        @test returned === W
        @test W.target_syndrome == syn
        @test W.channel_llrs[1] == 0.0            # erased: neutral belief
        @test W.channel_llrs[3] == 0.0
        @test W.channel_llrs[5] == CodingTheory._LLR_MAX    # pinned to 0
        @test W.channel_llrs[6] == -CodingTheory._LLR_MAX   # pinned to 1
        @test W.is_decimated[5] && W.is_decimated[6]
        @test count(W.is_decimated) == 2
        @test W.channel_llrs[7] == llr[7]
        @test W.total_llrs == W.channel_llrs
        @test all(iszero, W.V2C) && all(iszero, W.C2V)
        @test all(==(0xFF), W.prev_bits_1) && all(==(0xFF), W.prev_bits_2)

        # A non-zero syndrome is normalized to 0x00/0x01, and reloading resets state.
        CodingTheory.load_soft_channel!(W, llr; syndrome = fill(3, num_check))
        @test all(==(0x01), W.target_syndrome)
        CodingTheory.load_soft_channel!(W, llr)
        @test all(iszero, W.target_syndrome)
        @test !any(W.is_decimated)
        @test W.channel_llrs == llr

        @test_throws ArgumentError CodingTheory.load_soft_channel!(W, llr[1:3])
        @test_throws ArgumentError CodingTheory.load_soft_channel!(W, llr;
                                                                   syndrome = zeros(UInt8, 2))
        @test_throws ArgumentError CodingTheory.load_soft_channel!(W, llr;
                                                                   decimated_bits = [1],
                                                                   decimated_values = Int[])
        # Flooding built no partition, so the layered engine has nothing to sweep.
        @test_throws ArgumentError CodingTheory.decode!(W, llr; schedule = :layered)
        @test_throws ArgumentError CodingTheory.decode!(W, llr; schedule = :serial)
    end

    @testset "Codeword decoding: zero syndrome" begin
        H = qc_ldpc(7, 3, 5)
        num_var = size(H, 2)
        p = 0.05
        llr = fill(log((1 - p) / p), num_var)

        for schedule in (:flooding, :layered, :serial)
            for algorithm in (:sum_product, :min_sum, :normalized_min_sum, :offset_min_sum,
                              :min_sum_correction)
                W = CodingTheory.init_soft_workspace(H; schedule = schedule)
                out = fill(0xFF, num_var)
                converged, iters = CodingTheory.decode!(W, llr; algorithm = algorithm,
                                                        schedule = schedule, out = out)
                # Nothing to correct: the very first hard decision already matches.
                @test converged
                @test iters == 1
                @test all(iszero, out)
                @test all(iszero, W.current_bits)
            end
        end

        # Same, with a genuine codeword's LLRs rather than a syndrome: c = 0 is a
        # codeword of every parity-check matrix, and so is any row-space complement.
        F = Oscar.Nemo.Native.GF(2)
        Hh = matrix(F, Int.(hamming74))
        W = CodingTheory.init_soft_workspace(Hh)
        # [1, 1, 1, 0, 0, 0, 0] is not in the kernel of this H; [1, 1, 0, 1, 0, 0, 1] is.
        codeword = UInt8[1, 1, 0, 1, 0, 0, 1]
        @test all(iszero, syndrome_of(hamming74, codeword))
        out = fill(0xFF, 7)
        converged, iters = CodingTheory.decode!(W, [c == 0x01 ? -5.0 : 5.0 for c in codeword];
                                                algorithm = :sum_product, out = out)
        @test converged
        @test iters == 1
        @test out == codeword
    end

    @testset "Syndrome decoding on a sparse LDPC code" begin
        H = qc_ldpc(13, 3, 5)
        num_check, num_var = size(H)
        @test all(==(3), sum(Int.(H), dims = 1))    # column weight 3
        p = 0.05
        llr = fill(log((1 - p) / p), num_var)

        for schedule in (:flooding, :layered)
            for algorithm in (:sum_product, :min_sum, :normalized_min_sum, :offset_min_sum,
                              :min_sum_correction)
                W = CodingTheory.init_soft_workspace(H; schedule = schedule)
                rng = MersenneTwister(20260902)
                out = zeros(UInt8, num_var)
                for trial in 1:40
                    e = zeros(UInt8, num_var)
                    e[rand(rng, 1:num_var)] ⊻= 0x01
                    trial > 20 && (e[rand(rng, 1:num_var)] ⊻= 0x01)
                    s = syndrome_of(H, e)
                    converged, iters = CodingTheory.decode!(W, llr; algorithm = algorithm,
                                                            schedule = schedule, syndrome = s,
                                                            max_iter = 30, out = out)
                    @test converged
                    @test 1 <= iters <= 30
                    # Syndrome satisfaction, not exact recovery: any coset
                    # representative is a legitimate answer.
                    @test syndrome_of(H, out) == s
                end
            end
        end

        # Erasures: zeroing the channel belief on the flipped positions still leaves
        # enough structure to satisfy the syndrome.
        W = CodingTheory.init_soft_workspace(H; schedule = :layered)
        e = zeros(UInt8, num_var)
        e[3] = 0x01
        e[41] = 0x01
        s = syndrome_of(H, e)
        out = zeros(UInt8, num_var)
        converged, _ = CodingTheory.decode!(W, llr; algorithm = :sum_product,
                                            schedule = :layered, syndrome = s,
                                            erasures = [3, 41], max_iter = 50, out = out)
        @test converged
        @test syndrome_of(H, out) == s

        # Decimation and oscillation hooks all run to a valid answer on this code.
        for decimation in (:none, :auto, :guided, :manual), oscillation in (:none, :active)
            Wd = CodingTheory.init_soft_workspace(H; schedule = :layered)
            outd = zeros(UInt8, num_var)
            converged, _ = CodingTheory.decode!(Wd, llr; algorithm = :sum_product,
                                                schedule = :layered, syndrome = s,
                                                decimation = decimation,
                                                oscillation = oscillation,
                                                dec_thresh = 8.0, dec_rounds = 3,
                                                max_iter = 50, out = outd)
            @test converged
            @test syndrome_of(H, outd) == s
        end

        # Manual decimation pins the bits it is told to, and the pinned belief
        # survives the decode.
        Wm = CodingTheory.init_soft_workspace(H)
        outm = zeros(UInt8, num_var)
        converged, _ = CodingTheory.decode!(Wm, llr; algorithm = :sum_product, syndrome = s,
                                            decimation = :manual, decimated_bits = [3, 41],
                                            decimated_values = [1, 1], max_iter = 50,
                                            out = outm)
        @test converged
        @test Wm.is_decimated[3] && Wm.is_decimated[41]
        @test Wm.channel_llrs[3] == -CodingTheory._LLR_MAX
        @test outm[3] == 0x01 && outm[41] == 0x01
        @test syndrome_of(H, outm) == s

        # Auto decimation actually pins something once beliefs pass the threshold.
        Wa = CodingTheory.init_soft_workspace(H)
        CodingTheory.decode!(Wa, llr; algorithm = :sum_product, syndrome = s,
                             decimation = :auto, dec_thresh = 1.0, max_iter = 20)
        @test any(Wa.is_decimated)
    end

    @testset "The out buffer and the return shape" begin
        H = qc_ldpc(7, 3, 5)
        num_var = size(H, 2)
        p = 0.05
        llr = fill(log((1 - p) / p), num_var)
        e = zeros(UInt8, num_var)
        e[9] = 0x01
        s = syndrome_of(H, e)

        W = CodingTheory.init_soft_workspace(H)
        result = CodingTheory.decode!(W, llr; algorithm = :sum_product, syndrome = s,
                                      max_iter = 30)
        # Exactly two values, so that nothing crosses a language boundary unasked.
        @test result isa Tuple{Bool, Int}
        @test length(result) == 2
        @test result[1]
        # The decisions are still reachable on the workspace.
        @test syndrome_of(H, W.current_bits) == s

        # `out` is filled, and matches the workspace's own copy.
        for buffer in (fill(0xFF, num_var), fill(-1, num_var), fill(Int32(-1), num_var))
            converged, _ = CodingTheory.decode!(W, llr; algorithm = :sum_product,
                                                syndrome = s, max_iter = 30, out = buffer)
            @test converged
            @test buffer == W.current_bits
            @test all(b -> b in (0, 1), buffer)
            @test syndrome_of(H, UInt8.(buffer)) == s
        end

        # A failing decode still fills `out` with the decisions it ended on, and
        # reports the iteration budget it exhausted.
        W_stall = CodingTheory.init_soft_workspace(hamming74; schedule = :layered)
        stall_out = fill(0xFF, 7)
        stall_result = CodingTheory.decode!(W_stall, fill(log(0.9 / 0.1), 7);
                                            algorithm = :min_sum, schedule = :layered,
                                            syndrome = UInt8[1, 0, 0], max_iter = 5,
                                            out = stall_out)
        @test stall_result == (false, 5)
        @test stall_out == W_stall.current_bits
        @test all(b -> b in (0x00, 0x01), stall_out)
    end

    @testset "One layer equals serial" begin
        # The layered inner loop is strictly sequential, so a partition into a single
        # layer sweeps the checks in exactly the order the serial schedule does. The
        # partition documents independence; it does not change the arithmetic.
        H = qc_ldpc(13, 3, 5)
        num_check, num_var = size(H)
        p = 0.05
        llr = fill(log((1 - p) / p), num_var)

        W_serial = CodingTheory.init_soft_workspace(H; schedule = :serial)
        W_one = CodingTheory.init_soft_workspace(H; layer_ptr = [1, num_check + 1],
                                                 layer_checks = collect(1:num_check))
        @test length(W_serial.layer_ptr) == num_check + 1
        @test length(W_one.layer_ptr) == 2

        rng = MersenneTwister(13579)
        for algorithm in (:sum_product, :min_sum, :normalized_min_sum, :offset_min_sum,
                          :min_sum_correction, :sum_product_fast, :min_sum_correction_fast)
            for _ in 1:8
                e = zeros(UInt8, num_var)
                e[rand(rng, 1:num_var)] ⊻= 0x01
                e[rand(rng, 1:num_var)] ⊻= 0x01
                s = syndrome_of(H, e)
                a = CodingTheory.decode!(W_serial, llr; algorithm = algorithm,
                                         schedule = :serial, syndrome = s, max_iter = 25)
                bits_a, llrs_a = copy(W_serial.current_bits), copy(W_serial.total_llrs)
                b = CodingTheory.decode!(W_one, llr; algorithm = algorithm,
                                         schedule = :layered, syndrome = s, max_iter = 25)
                @test a == b
                @test W_one.current_bits == bits_a
                @test W_one.total_llrs == llrs_a   # bit-identical, not merely close
            end
        end

        # `:semiserial` is routed to the same engine as `:serial`.
        e = zeros(UInt8, num_var)
        e[7] = 0x01
        s = syndrome_of(H, e)
        @test CodingTheory.decode!(W_serial, llr; schedule = :serial, syndrome = s,
                                   max_iter = 25) ==
              CodingTheory.decode!(W_serial, llr; schedule = :semiserial, syndrome = s,
                                   max_iter = 25)
    end

    @testset "Fast variants agree with their reference rules" begin
        # `:sum_product_fast` swaps the Jacobian logarithm for the tanh product form.
        # Same rule algebraically; equal here because no belief approaches the clamp.
        H = qc_ldpc(13, 3, 5)
        num_var = size(H, 2)
        p = 0.05
        llr = fill(log((1 - p) / p), num_var)

        W_ref = CodingTheory.init_soft_workspace(H)
        W_fast = CodingTheory.init_soft_workspace(H)
        rng = MersenneTwister(864209)
        for _ in 1:25
            e = zeros(UInt8, num_var)
            e[rand(rng, 1:num_var)] ⊻= 0x01
            e[rand(rng, 1:num_var)] ⊻= 0x01
            s = syndrome_of(H, e)
            ref = CodingTheory.decode!(W_ref, llr; algorithm = :sum_product, syndrome = s,
                                       max_iter = 30)
            bits_ref = copy(W_ref.current_bits)
            fast = CodingTheory.decode!(W_fast, llr; algorithm = :sum_product_fast,
                                        syndrome = s, max_iter = 30)
            @test ref == fast
            @test W_fast.current_bits == bits_ref
        end

        # `:min_sum_correction_fast` is deliberately a different decoder: it moves the
        # threshold test out of the fold. It must still be a working decoder.
        W_mscf = CodingTheory.init_soft_workspace(H; schedule = :layered)
        out = zeros(UInt8, num_var)
        rng = MersenneTwister(864209)
        for _ in 1:25
            e = zeros(UInt8, num_var)
            e[rand(rng, 1:num_var)] ⊻= 0x01
            e[rand(rng, 1:num_var)] ⊻= 0x01
            s = syndrome_of(H, e)
            converged, _ = CodingTheory.decode!(W_mscf, llr;
                                                algorithm = :min_sum_correction_fast,
                                                schedule = :layered, syndrome = s,
                                                max_iter = 30, out = out)
            @test converged
            @test syndrome_of(H, out) == s
        end
    end

    @testset "Parameters change the answer they are supposed to" begin
        H = qc_ldpc(7, 3, 5)
        num_var = size(H, 2)
        p = 0.05
        llr = fill(log((1 - p) / p), num_var)
        e = zeros(UInt8, num_var)
        e[2] = 0x01
        e[19] = 0x01
        s = syndrome_of(H, e)

        # `attenuation` scales the normalized min-sum aggregate; `offset` shifts the
        # offset min-sum one. Both must reach the posteriors.
        W = CodingTheory.init_soft_workspace(H)
        CodingTheory.decode!(W, llr; algorithm = :normalized_min_sum, syndrome = s,
                             attenuation = 0.75, max_iter = 1)
        a75 = copy(W.total_llrs)
        CodingTheory.decode!(W, llr; algorithm = :normalized_min_sum, syndrome = s,
                             attenuation = 0.5, max_iter = 1)
        @test W.total_llrs != a75
        # The default is 0.75.
        CodingTheory.decode!(W, llr; algorithm = :normalized_min_sum, syndrome = s,
                             max_iter = 1)
        @test W.total_llrs == a75

        CodingTheory.decode!(W, llr; algorithm = :offset_min_sum, syndrome = s, offset = 0.5,
                             max_iter = 1)
        b50 = copy(W.total_llrs)
        CodingTheory.decode!(W, llr; algorithm = :offset_min_sum, syndrome = s, offset = 0.25,
                             max_iter = 1)
        @test W.total_llrs != b50
        # The default is 0.5, which is also the default algorithm.
        CodingTheory.decode!(W, llr; algorithm = :offset_min_sum, syndrome = s, max_iter = 1)
        @test W.total_llrs == b50
        CodingTheory.decode!(W, llr; syndrome = s, max_iter = 1)
        @test W.total_llrs == b50

        # `max_iter` is honoured, and a failure reports it.
        W_hard = CodingTheory.init_soft_workspace(hamming74; schedule = :layered)
        hllr = fill(log((1 - 0.1) / 0.1), 7)
        for limit in (3, 11)
            @test CodingTheory.decode!(W_hard, hllr; algorithm = :min_sum,
                                       schedule = :layered, syndrome = UInt8[1, 0, 0],
                                       max_iter = limit) == (false, limit)
        end
    end

    @testset "Min-sum stalls under a layered schedule on a dense high-rate matrix" begin
        # Documented degeneracy, not a defect. With a CONSTANT channel LLR every input
        # to an unsatisfied check shares one magnitude, so min-sum's outgoing magnitude
        # equals it exactly, the posterior lands on exactly 0.0, and the zeros then
        # propagate to every other edge. The state repeats from the second iteration.
        p = 0.1
        llr = fill(log((1 - p) / p), 7)
        stalling = UInt8[1, 0, 0]

        W = CodingTheory.init_soft_workspace(hamming74; schedule = :layered)
        @test CodingTheory.decode!(W, llr; algorithm = :min_sum, schedule = :layered,
                                   syndrome = stalling, max_iter = 100) == (false, 100)
        # `oscillation = :active` catches the two-cycle on the second iteration and
        # reports it with a negative iteration count.
        @test CodingTheory.decode!(W, llr; algorithm = :min_sum, schedule = :layered,
                                   syndrome = stalling, oscillation = :active,
                                   max_iter = 100) == (false, -2)
        # `:offset_min_sum` reaches the same fixed point; `max(0, |agg| - offset)`
        # maps the surviving magnitudes to zero.
        @test CodingTheory.decode!(W, llr; algorithm = :offset_min_sum, schedule = :layered,
                                   syndrome = stalling, oscillation = :active,
                                   max_iter = 100) == (false, -2)

        # Flooding is unaffected: it converges on every syndrome of this matrix, for
        # every algorithm.
        for algorithm in (:sum_product, :min_sum, :normalized_min_sum, :offset_min_sum,
                          :min_sum_correction)
            Wf = CodingTheory.init_soft_workspace(hamming74)
            out = zeros(UInt8, 7)
            for bits in 0:7
                s = UInt8[(bits >> k) & 1 for k in 0:2]
                converged, _ = CodingTheory.decode!(Wf, llr; algorithm = algorithm,
                                                    syndrome = s, max_iter = 50, out = out)
                @test converged
                @test syndrome_of(hamming74, out) == s
            end
        end

        # Sum-product under the same layered schedule still handles most syndromes,
        # and every syndrome it does converge on it converges on correctly.
        Wsp = CodingTheory.init_soft_workspace(hamming74; schedule = :layered)
        out = zeros(UInt8, 7)
        solved = 0
        for bits in 0:7
            s = UInt8[(bits >> k) & 1 for k in 0:2]
            converged, _ = CodingTheory.decode!(Wsp, llr; algorithm = :sum_product,
                                                schedule = :layered, syndrome = s,
                                                max_iter = 100, out = out)
            converged || continue
            solved += 1
            @test syndrome_of(hamming74, out) == s
        end
        @test solved >= 6

        # Breaking the tie in the channel LLRs removes the degeneracy entirely: the
        # stall is a property of the constant input, not of the schedule.
        jittered = [log((1 - p) / p) + 0.05 * v for v in 1:7]
        Wj = CodingTheory.init_soft_workspace(hamming74; schedule = :layered)
        outj = zeros(UInt8, 7)
        for bits in 0:7
            s = UInt8[(bits >> k) & 1 for k in 0:2]
            converged, _ = CodingTheory.decode!(Wj, jittered; algorithm = :min_sum,
                                                schedule = :layered, syndrome = s,
                                                max_iter = 100, out = outj)
            @test converged
            @test syndrome_of(hamming74, outj) == s
        end
    end
end
