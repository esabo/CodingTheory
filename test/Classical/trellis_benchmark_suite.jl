using Oscar, CodingTheory, Random, BenchmarkTools, Printf

function run_trellis_benchmark()
    F2 = Oscar.Nemo.Native.GF(2)
    q = 2
    
    println("==================================================")
    println("       CLASSICAL TRELLIS BENCHMARK SUITE          ")
    println("==================================================\n")
    
    test_sizes = [(10, 20), (15, 30), (20, 40), (25, 50)]
    
    for (k, n) in test_sizes
        @printf("Testing Code Size: [%d, %d]\n", n, k)
        
        G = matrix(F2, rand(0:1, k, n))
        while rank(G) < k
            G = matrix(F2, rand(0:1, k, n))
        end
        
        # 1. Profile Trellis Complexity
        M = Array(G)
        t_opt = @elapsed best_M, best_perm, peak_E = CodingTheory.optimize_trellis_permutation(M, 50)
        
        @printf("  -> Optimized Peak Active Span (E): %d (Opt Time: %.4fs)\n", peak_E, t_opt)
        
        # Precompute profiles for routing
        L, R = CodingTheory._get_LR_indices(best_M)
        past, future = CodingTheory.past_future_profiles(L, R, n)
        E_profile = [k - past[b] - future[b+1] for b in 1:n]
        
        # 2. Pure Viterbi Decoding (Min-Weight Engine)
        if peak_E <= 20 
            print("  -> Running Pure Viterbi Sectionalized... ")
            t_pure = @elapsed begin
                bounds = CodingTheory.optimal_sectionalization(best_M, q)
                # Swap to the fast min-weight integer engine; verbose=true to keep bench output clean
                d_pure = CodingTheory._min_weight_TP_Viterbi_sectionalized(best_M, bounds, true)
            end
            @printf("d = %d (Took: %.4f seconds)\n", d_pure, t_pure)
        else
            println("  -> Pure Viterbi Skipped (Peak E > 20 would exhaust memory).")
        end
        
        # 3. Hybrid BZ-DFS Bridge
        pinch_span = max(1, peak_E - 5) 
        print("  -> Running Hybrid Bridge (max_span = $pinch_span)... \n")
        t_hybrid = @elapsed begin
            B_L = findfirst(e -> e > pinch_span, E_profile) - 1
            B_R = findlast(e -> e > pinch_span, E_profile) 
            
            # Pass verbose=true to the boundary builds
            left_dict = CodingTheory._forward_trellis(best_M, B_L, q, 0, true)
            right_dict = CodingTheory._backward_trellis(best_M, n, B_R, q, true)
            
            # Pass verbose=true to see the ProgressMeter for the DFS Bridge
            d_hybrid = CodingTheory._BZ_middle_search(best_M, L, R, B_L, B_R, left_dict, right_dict, q, true)
        end
        @printf("d = %d (Took: %.4f seconds)\n", d_hybrid, t_hybrid)
        
        println("-"^50)
    end

    println("==================================================")
    println("       CLASSICAL TRELLIS EXTREME BENCHMARK        ")
    println("==================================================\n")
    
    test_sizes = [
        (10, 50, "Random Low Rate"),
        (25, 50, "Random Mid Rate"),
        (10, 128, "Random Low Rate (n=128)"),
        (118, 128, "Random High Rate (n=128) -> Auto-Routes to Syndrome")
    ]
    
    for (k, n, desc) in test_sizes
        @printf("Testing: %s [%d, %d]\n", desc, n, k)
        
        G = matrix(F2, rand(0:1, k, n))
        while rank(G) < k
            G = matrix(F2, rand(0:1, k, n))
        end
        C = LinearCode(G)
        
        # 1. Profile & Optimize Trellis Complexity
        is_generator = k <= n / 2
        mat = is_generator ? Array(generator_matrix(C)) : Array(parity_check_matrix(C))
        route_str = is_generator ? "Generator" : "Syndrome"
        
        t_opt = @elapsed best_M, best_perm, peak_E = CodingTheory.optimize_trellis_permutation(mat, 50)
        
        @printf("  -> Routed to: %s Matrix\n", route_str)
        @printf("  -> Optimized Peak Active Span (E): %d (Opt Time: %.4fs)\n", peak_E, t_opt)
        
        # 2. Pure Viterbi Decoding
        if peak_E <= 20
            print("  -> Running Pure Viterbi Sectionalized... ")
            t_pure = @elapsed begin
                if is_generator
                    bounds = CodingTheory.optimal_sectionalization(best_M, q)
                    d_pure = CodingTheory._min_weight_TP_Viterbi_sectionalized(best_M, bounds, true)
                else
                    # THE SYNDROME BYPASS: Force 1-column sections to prevent branch explosions
                    bounds = collect(0:n)
                    d_pure = CodingTheory._min_weight_syndrome_sectionalized(best_M, bounds, true)
                end
            end
            @printf("d = %d (Took: %.4f seconds)\n", d_pure, t_pure)
        else
            println("  -> Pure Viterbi Skipped (Peak E > 20 would exhaust memory).")
            
            # 3. Hybrid BZ-DFS Bridge
            # Note: If it hit this branch, we MUST use a Generator Matrix for the DFS bridge.
            if !is_generator
                mat_gen = Array(generator_matrix(C))
                best_M, _, peak_E_gen = CodingTheory.optimize_trellis_permutation(mat_gen, 50)
            end
            
            # Recompute profiles for the (potentially swapped) best_M
            L, R = CodingTheory._get_LR_indices(best_M)
            past, future = CodingTheory.past_future_profiles(L, R, n)
            active_rows = size(best_M, 1)
            E_profile = [active_rows - past[b] - future[b+1] for b in 1:n]
            
            pinch_span = 18 
            print("  -> Running Hybrid Bridge (max_span = $pinch_span)... \n")
            t_hybrid = @elapsed begin
                B_L = findfirst(e -> e > pinch_span, E_profile) - 1
                B_R = findlast(e -> e > pinch_span, E_profile)
                
                # Threads with verbose=true
                t_left = Threads.@spawn CodingTheory._forward_trellis(best_M, B_L, q, 0, true)
                t_right = Threads.@spawn CodingTheory._backward_trellis(best_M, n, B_R, q, true)
                left_dict = fetch(t_left)
                right_dict = fetch(t_right)
                
                # Execute DFS Bridge with verbose=true
                d_hybrid = CodingTheory._BZ_middle_search(best_M, L, R, B_L, B_R, left_dict, right_dict, q, true)
            end
            @printf("d = %d (Took: %.4f seconds)\n", d_hybrid, t_hybrid)
        end
        
        println("-"^50)
    end
end

run_trellis_benchmark()

# julia> run_trellis_benchmark()
# ==================================================
#        CLASSICAL TRELLIS BENCHMARK SUITE          
# ==================================================

# Testing Code Size: [20, 10]
#   -> Optimized Peak Active Span (E): 7 (Opt Time: 0.0002s)
#   -> Running Pure Viterbi Sectionalized... d = 2 (Took: 0.0013 seconds)
#   -> Running Hybrid Bridge (max_span = 2)...   -> Packing 16 columns into hardware registers...
# d = 2 (Took: 0.0635 seconds)
# --------------------------------------------------
# Testing Code Size: [30, 15]
#   -> Optimized Peak Active Span (E): 13 (Opt Time: 0.0005s)
#   -> Running Pure Viterbi Sectionalized... d = 4 (Took: 0.1819 seconds)
#   -> Running Hybrid Bridge (max_span = 8)...   -> Packing 14 columns into hardware registers...
# d = 4 (Took: 0.0012 seconds)
# --------------------------------------------------
# Testing Code Size: [40, 20]
#   -> Optimized Peak Active Span (E): 17 (Opt Time: 0.0012s)
#   -> Running Pure Viterbi Sectionalized... d = 5 (Took: 94.4809 seconds)
#   -> Running Hybrid Bridge (max_span = 12)...   -> Packing 16 columns into hardware registers...
# d = 5 (Took: 0.1389 seconds)
# --------------------------------------------------
# Testing Code Size: [50, 25]
#   -> Optimized Peak Active Span (E): 22 (Opt Time: 0.0024s)
#   -> Pure Viterbi Skipped (Peak E > 20 would exhaust memory).
#   -> Running Hybrid Bridge (max_span = 17)...   -> Packing 16 columns into hardware registers...
# Bridging Middle Trellis: 100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| Time: 0:00:01
# d = 6 (Took: 215.6349 seconds)
# --------------------------------------------------
# ==================================================
#        CLASSICAL TRELLIS EXTREME BENCHMARK        
# ==================================================

# Testing: Random Low Rate [50, 10]
#   -> Routed to: Generator Matrix
#   -> Optimized Peak Active Span (E): 10 (Opt Time: 0.0004s)
#   -> Running Pure Viterbi Sectionalized... d = 14 (Took: 0.0157 seconds)
# --------------------------------------------------
# Testing: Random Mid Rate [50, 25]
#   -> Routed to: Generator Matrix
#   -> Optimized Peak Active Span (E): 22 (Opt Time: 0.0027s)
#   -> Pure Viterbi Skipped (Peak E > 20 would exhaust memory).
#   -> Running Hybrid Bridge (max_span = 18)...   -> Packing 13 columns into hardware registers...
# d = 7 (Took: 1138.7888 seconds)
# --------------------------------------------------
# Testing: Random Low Rate (n=128) [128, 10]
#   -> Routed to: Generator Matrix
#   -> Optimized Peak Active Span (E): 10 (Opt Time: 0.0008s)
#   -> Running Pure Viterbi Sectionalized... d = 47 (Took: 0.0389 seconds)
# --------------------------------------------------
# Testing: Random High Rate (n=128) -> Auto-Routes to Syndrome [128, 118]
#   -> Routed to: Syndrome Matrix
#   -> Optimized Peak Active Span (E): 10 (Opt Time: 0.0007s)
#   -> Running Pure Viterbi Sectionalized... d = 2 (Took: 242.1872 seconds)
# --------------------------------------------------