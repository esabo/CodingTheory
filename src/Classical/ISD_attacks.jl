# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
          # Binary
#############################

"""
    _Stern_attack_binary(G::Matrix{Int}, target_w::Int; kwargs...)

Executes Stern's Information Set Decoding (ISD) algorithm optimized for binary (GF(2)) 
block codes. This algorithm uses a Meet-in-the-Middle hash collision technique 
to efficiently search for low-weight vectors in the code space.

### Dual Modes of Operation
1. **Minimum Distance Search (Default):** If `w_recv` is omitted (or set to the all-zero vector), 
   the algorithm searches for a valid non-zero codeword of exactly weight `target_w`.
2. **Syndrome Decoding (Error Correction):** If a noisy received word is passed to `w_recv`, 
   the algorithm searches for an error vector of weight `target_w` such that 
   `H * error == H * w_recv`.

### Arguments
* `G::Matrix{Int}`: The generator matrix of the binary linear block code.
* `target_w::Int`: The target Hamming weight of the codeword or error vector to find.

### Keyword Arguments
* `w_recv::Vector{Int} = zeros(Int, C.n)`: The received word. Defaults to the all-zero 
  vector for Minimum Distance searches.
* `p::Int = 2`: The combinatorial search weight per half of the Information Set. 
  Assumes exactly `2*p` errors land in the `k` information columns.
* `l::Int = 12`: The size of the collision window in the redundancy section. For binary codes, 
  `l = 12` perfectly balances the hash map size against the false-positive collision rate.
* `num_find::Int = 1`: The number of valid vectors to find before terminating the search.
* `max_iters::Int = 10000`: The maximum number of random matrix permutations to attempt 
  before giving up.

### Returns
* `Set{Vector{Int}}`: A set containing the discovered vectors of weight `target_w`. 
  Returns an empty set if `max_iters` is reached without finding `num_find` matches.
"""
function _Stern_attack_binary(G::Matrix{Int}, target_w::Int; w_recv::Vector{Int} = 
    zeros(Int, size(G, 2)), p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)

    k, n = size(G)
    @assert l <= 64 "For max performance, window size l should be <= 64 to fit in one CPU register."
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds
    remaining = max_iters % num_thrds
    
    keep_going = Threads.Atomic{Bool}(true)
    results_lock = Threads.SpinLock()
    
    # NEW: Initialize the collection Set
    found_vectors = Set{Vector{Int}}()
    
    Threads.@threads for th in 1:num_thrds
        iters_for_this_thread = thread_load + (th <= remaining ? 1 : 0)
        
        # Pre-allocate thread-local buffers to avoid GC pauses
        half_k = k ÷ 2
        X_cols = 1:half_k
        Y_cols = (half_k + 1):k
        tail_len = n - k - l
        num_tail_chunks = cld(tail_len, 64)
        
        tail_buffer = zeros(UInt64, num_tail_chunks)
        tail_rows = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
        for _ in 1:iters_for_this_thread
            if !keep_going[] break end
            
            # 1. Systematize [I_k | P]
            σ_loc = collect(1:n)
            shuffle!(σ_loc)
            
            G_loc = G[:, σ_loc]
            
            try
                _make_systematic_gf!(G_loc, σ_loc, k)
            catch
                continue # Skip if singular
            end
            
            w_loc = w_recv[σ_loc]
            
            # 2. 64-BIT PACKING
            window_rows = zeros(UInt64, k)
            for i in 1:k
                val = UInt64(0)
                for j in 1:l
                    if G_loc[i, k + j] == 1
                        val |= (UInt64(1) << (j - 1))
                    end
                end
                window_rows[i] = val
            end
            
            for i in 1:k
                for j in 1:tail_len
                    if G_loc[i, k + l + j] == 1
                        chunk = (j - 1) ÷ 64 + 1
                        bit   = (j - 1) % 64
                        tail_rows[i][chunk] |= (UInt64(1) << bit)
                    end
                end
            end
            
            # 3. COMPUTE TARGETS (w_recv offset logic)
            # target_window represents the syndrome of w_recv within the l-window
            target_window = UInt64(0)
            for j in 1:l
                if w_loc[k + j] == 1
                    target_window |= (UInt64(1) << (j - 1))
                end
            end
            for i in 1:k
                if w_loc[i] == 1
                    target_window ⊻= window_rows[i]
                end
            end
            
            target_tail = zeros(UInt64, num_tail_chunks)
            for j in 1:tail_len
                if w_loc[k + l + j] == 1
                    chunk = (j - 1) ÷ 64 + 1
                    bit   = (j - 1) % 64
                    target_tail[chunk] |= (UInt64(1) << bit)
                end
            end
            for i in 1:k
                if w_loc[i] == 1
                    target_tail .⊻= tail_rows[i]
                end
            end

            # 4. RECURSIVE HASH MAP BUILDER (X-Half)
            hash_X = Dict{UInt64, Vector{Vector{Int}}}()
            
            function _build_X!(depth, picked, current_val::UInt64, msg)
                if picked == p
                    if haskey(hash_X, current_val)
                        push!(hash_X[current_val], copy(msg))
                    else
                        hash_X[current_val] = [copy(msg)]
                    end
                    return
                end
                if depth > length(X_cols) || (length(X_cols) - depth + 1) < (p - picked) return end
                
                msg[depth] = 0
                _build_X!(depth+1, picked, current_val, msg)
                
                msg[depth] = 1
                _build_X!(depth+1, picked+1, current_val ⊻ window_rows[X_cols[depth]], msg)
            end
            
            _build_X!(1, 0, UInt64(0), zeros(Int, length(X_cols)))
            
            # 5. RECURSIVE COLLISION PROBER (Y-Half)
            function _probe_Y!(depth, picked, current_val::UInt64, msg_Y)
                if !keep_going[] return end
                
                if picked == p
                    search_val = target_window ⊻ current_val
                    
                    if haskey(hash_X, search_val)
                        for msg_X in hash_X[search_val]
                            
                            tail_buffer .= target_tail
                            for i in 1:length(msg_X)
                                if msg_X[i] == 1
                                    tail_buffer .⊻= tail_rows[X_cols[i]]
                                end
                            end
                            for i in 1:length(msg_Y)
                                if msg_Y[i] == 1
                                    tail_buffer .⊻= tail_rows[Y_cols[i]]
                                end
                            end
                            
                            wt_tail = sum(count_ones.(tail_buffer))
                            
                            # Check if the generated error vector perfectly matches the target weight
                            if 0 < wt_tail + 2*p <= target_w
                                lock(results_lock) do
                                    # Double check we haven't hit the limit while waiting for the lock
                                    if length(found_vectors) < num_find
                                        
                                        # Construct the pure ERROR vector
                                        e_loc = zeros(Int, n)
                                        
                                        for i in 1:length(msg_X)
                                            if msg_X[i] == 1 e_loc[X_cols[i]] = 1 end
                                        end
                                        for i in 1:length(msg_Y)
                                            if msg_Y[i] == 1 e_loc[Y_cols[i]] = 1 end
                                        end
                                        
                                        # The l-window error is 0 by definition of the MitM collision
                                        
                                        for j in 1:tail_len
                                            chunk = (j - 1) ÷ 64 + 1
                                            bit   = (j - 1) % 64
                                            if (tail_buffer[chunk] & (UInt64(1) << bit)) != 0
                                                e_loc[k + l + j] = 1
                                            end
                                        end
                                        
                                        # Un-permute and push to our Set
                                        push!(found_vectors, e_loc[invperm(σ_loc)])
                                        
                                        # Terminate early if we hit the requested number of vectors
                                        if length(found_vectors) >= num_find
                                            keep_going[] = false
                                        end
                                    end
                                end
                            end
                        end
                    end
                    return
                end
                
                if depth > length(Y_cols) || (length(Y_cols) - depth + 1) < (p - picked) return end
                
                msg_Y[depth] = 0
                _probe_Y!(depth+1, picked, current_val, msg_Y)
                
                msg_Y[depth] = 1
                _probe_Y!(depth+1, picked+1, current_val ⊻ window_rows[Y_cols[depth]], msg_Y)
            end
            
            _probe_Y!(1, 0, UInt64(0), zeros(Int, length(Y_cols)))
        end
    end
    
    return found_vectors
end

"""
    _Prange_attack_binary(G::Matrix{Int}, target_w::Int; w_recv::Vector{Int} = zeros(Int, size(G, 2)), num_find::Int = 1, max_iters::Int = 10000)

Executes Prange's Information Set Decoding algorithm for binary block codes.
"""
function _Prange_attack_binary(G::Matrix{Int}, target_w::Int; w_recv::Vector{Int} = zeros(Int,
    size(G, 2)), num_find::Int = 1, max_iters::Int = 10000)

    k, n = size(G)
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds
    remaining = max_iters % num_thrds
    
    keep_going = Threads.Atomic{Bool}(true)
    results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    is_min_dist = all(iszero, w_recv)
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            σ_loc = collect(1:n)
            shuffle!(σ_loc)
            G_loc = G[:, σ_loc]
            
            try
                _make_systematic_gf!(G_loc, σ_loc, k) # GF(2) version
            catch
                continue
            end
            
            if is_min_dist
                # MINIMUM DISTANCE MODE (p=1 Check)
                # The rows of G_sys are valid codewords with exactly one '1' in the info set.
                for i in 1:k
                    row_wt = sum(G_loc[i, j] for j in 1:n)
                    if 0 < row_wt <= target_w
                        lock(results_lock) do
                            if length(found_vectors) < num_find
                                e_loc = G_loc[i, :]
                                push!(found_vectors, e_loc[invperm(σ_loc)])
                                if length(found_vectors) >= num_find
                                    keep_going[] = false
                                end
                            end
                        end
                    end
                end
            else
                # SYNDROME DECODING MODE (p=0 Check)
                # We assume no errors in the information set.
                w_loc = w_recv[σ_loc]
                S = copy(w_loc[k+1:n])
                for i in 1:k
                    if w_loc[i] == 1
                        for j in 1:(n-k)
                            S[j] ⊻= G_loc[i, k+j]
                        end
                    end
                end
                
                if 0 < sum(S) <= target_w
                    lock(results_lock) do
                        if length(found_vectors) < num_find
                            e_loc = zeros(Int, n)
                            e_loc[k+1:n] .= S
                            push!(found_vectors, e_loc[invperm(σ_loc)])
                            if length(found_vectors) >= num_find
                                keep_going[] = false
                            end
                        end
                    end
                end
            end
        end
    end

    return found_vectors
end

"""
    _Lee_Brickell_attack_binary(G::Matrix{Int}, target_w::Int; w_recv::Vector{Int} = zeros(Int, size(G, 2)), p::Int = 2, num_find::Int = 1, max_iters::Int = 10000)

Executes the Lee-Brickell ISD algorithm for binary block codes.
"""
function _Lee_Brickell_attack_binary(G::Matrix{Int}, target_w::Int; w_recv::Vector{Int} =
    zeros(Int, size(G, 2)), p::Int = 2, num_find::Int = 1, max_iters::Int = 10000)

    k, n = size(G)
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds
    remaining = max_iters % num_thrds
    
    keep_going = Threads.Atomic{Bool}(true)
    results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        
        tail_len = n - k
        num_tail_chunks = cld(tail_len, 64)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            σ_loc = collect(1:n)
            shuffle!(σ_loc)
            G_loc = G[:, σ_loc]
            
            try 
                _make_systematic_gf!(G_loc, σ_loc, k) 
            catch
                continue 
            end
            
            w_loc = w_recv[σ_loc]
            
            # 1. PACK THE PARITY TAIL
            tail_rows = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            for i in 1:k, j in 1:tail_len
                if G_loc[i, k + j] == 1
                    tail_rows[i][(j - 1) ÷ 64 + 1] |= (UInt64(1) << ((j - 1) % 64))
                end
            end
            
            # Syndrome of w_recv
            target_tail = zeros(UInt64, num_tail_chunks)
            for j in 1:tail_len
                if w_loc[k + j] == 1
                    target_tail[(j - 1) ÷ 64 + 1] |= (UInt64(1) << ((j - 1) % 64))
                end
            end
            for i in 1:k
                if w_loc[i] == 1 
                    target_tail .⊻= tail_rows[i] 
                end
            end

            # 2. RECURSIVE XOR TREE
            function _probe_LB!(depth, picked, current_tail::Vector{UInt64}, msg)
                if !keep_going[] return end
                
                if picked == p
                    wt = sum(count_ones.(current_tail))
                    if 0 < wt + p <= target_w
                        lock(results_lock) do
                            if length(found_vectors) < num_find
                                # Reconstruct the TRUE error vector
                                e_loc = zeros(Int, n)
                                e_loc[1:k] .= msg
                                
                                for j in 1:tail_len
                                    chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                    if (current_tail[chunk] & (UInt64(1) << bit)) != 0
                                        e_loc[k + j] = 1
                                    end
                                end
                                
                                push!(found_vectors, e_loc[invperm(σ_loc)])
                                
                                if length(found_vectors) >= num_find
                                    keep_going[] = false
                                end
                            end
                        end
                    end
                    return
                end
                
                if depth > k || (k - depth + 1) < (p - picked) return end
                
                # Branch 0
                msg[depth] = 0
                _probe_LB!(depth+1, picked, current_tail, msg)
                
                # Branch 1 (Native SIMD XOR)
                msg[depth] = 1
                next_tail = current_tail .⊻ tail_rows[depth]
                _probe_LB!(depth+1, picked+1, next_tail, msg)
                msg[depth] = 0 # Backtrack
            end
            
            _probe_LB!(1, 0, copy(target_tail), zeros(Int, k))
        end
    end

    return found_vectors
end

"""
    _Leon_attack_binary(G::Matrix{Int}, target_w::Int; kwargs...)

Executes Leon's Information Set Decoding algorithm for binary block codes.
Uses a `Matrix{Int}` to pack the parity window and tail into `UInt64` registers 
for hardware-accelerated filtering and verification.
"""
function _Leon_attack_binary(G::Matrix{Int}, target_w::Int; 
                             w_recv::Vector{Int} = zeros(Int, size(G, 2)), 
                             p::Int = 2, l::Int = 12, 
                             num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    @assert l <= 64 "Window size l must be <= 64 for native register packing."
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds
    remaining = max_iters % num_thrds
    
    keep_going = Threads.Atomic{Bool}(true)
    results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        tail_len = n - k - l
        num_tail_chunks = cld(tail_len, 64)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            σ_loc = collect(1:n)
            shuffle!(σ_loc)
            G_loc = G[:, σ_loc]
            
            try 
                _make_systematic_gf!(G_loc, σ_loc, k) 
            catch
                continue 
            end
            
            w_loc = w_recv[σ_loc]
            
            # 1. PACK THE L-WINDOW
            window_rows = zeros(UInt64, k)
            for i in 1:k, j in 1:l
                if G_loc[i, k + j] == 1 
                    window_rows[i] |= (UInt64(1) << (j - 1)) 
                end
            end
            
            # Target Window Syndrome = w_window ⊻ w_info * P_window
            target_window = UInt64(0)
            for j in 1:l
                if w_loc[k + j] == 1 target_window |= (UInt64(1) << (j - 1)) end
            end
            for i in 1:k
                if w_loc[i] == 1 target_window ⊻= window_rows[i] end
            end
            
            # 2. PACK THE REMAINING TAIL
            tail_rows = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            for i in 1:k, j in 1:tail_len
                if G_loc[i, k + l + j] == 1
                    tail_rows[i][(j - 1) ÷ 64 + 1] |= (UInt64(1) << ((j - 1) % 64))
                end
            end
            
            target_tail = zeros(UInt64, num_tail_chunks)
            for j in 1:tail_len
                if w_loc[k + l + j] == 1
                    target_tail[(j - 1) ÷ 64 + 1] |= (UInt64(1) << ((j - 1) % 64))
                end
            end
            for i in 1:k
                if w_loc[i] == 1 target_tail .⊻= tail_rows[i] end
            end

            # 3. RECURSIVE LEON FILTER
            function _probe_Leon!(depth, picked, current_win::UInt64, current_tail::Vector{UInt64}, msg)
                if !keep_going[] return end
                
                if picked == p
                    # LEON FILTER: Is the error in the window exactly 0?
                    if current_win == UInt64(0)
                        # Only check the tail if it passes the window filter!
                        wt = sum(count_ones.(current_tail))
                        
                        if 0 < wt + p <= target_w
                            lock(results_lock) do
                                if length(found_vectors) < num_find
                                    e_loc = zeros(Int, n)
                                    e_loc[1:k] .= msg
                                    
                                    # Window (k+1 : k+l) is guaranteed to be 0
                                    
                                    # Extract the tail
                                    for j in 1:tail_len
                                        chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                        if (current_tail[chunk] & (UInt64(1) << bit)) != 0
                                            e_loc[k + l + j] = 1
                                        end
                                    end
                                    
                                    push!(found_vectors, e_loc[invperm(σ_loc)])
                                    
                                    if length(found_vectors) >= num_find
                                        keep_going[] = false
                                    end
                                end
                            end
                        end
                    end
                    return
                end
                
                if depth > k || (k - depth + 1) < (p - picked) return end
                
                # Branch 0
                msg[depth] = 0
                _probe_Leon!(depth+1, picked, current_win, current_tail, msg)
                
                # Branch 1
                msg[depth] = 1
                next_win = current_win ⊻ window_rows[depth]
                next_tail = current_tail .⊻ tail_rows[depth]
                _probe_Leon!(depth+1, picked+1, next_win, next_tail, msg)
                msg[depth] = 0 # Backtrack
            end
            
            _probe_Leon!(1, 0, target_window, copy(target_tail), zeros(Int, k))
        end
    end

    return found_vectors
end

"""
    _Canteaut_Chabaud_attack_binary(G::Matrix{Int}, target_w::Int; w_recv::Vector{Int} = zeros(Int, size(G, 2)), p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)

Executes the Canteaut-Chabaud variant of ISD for binary block codes.
This algorithm replaces the heavy O(n^3) full matrix permutation with 
O(n^2) single-column swaps and rapid pivot updates.
"""
function _Canteaut_Chabaud_attack_binary(G::Matrix{Int}, target_w::Int; w_recv::Vector{Int} =
    zeros(Int, size(G, 2)), p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)

    k, n = size(G)
    @assert l <= 64 "For max performance, window size l should be <= 64."
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds
    remaining = max_iters % num_thrds
    
    keep_going = Threads.Atomic{Bool}(true)
    results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    Threads.@threads for th in 1:num_thrds
        iters_for_this_thread = thread_load + (th <= remaining ? 1 : 0)
        
        half_k = k ÷ 2
        X_cols = 1:half_k
        Y_cols = (half_k + 1):k
        
        # Thread-local state
        σ_loc = collect(1:n)
        G_loc = zeros(Int, k, n)
        w_loc = zeros(Int, n)
        msg_X = zeros(Int, length(X_cols))
        msg_Y = zeros(Int, length(Y_cols))
        cur_win = zeros(Int, l)
        hash_X = Dict{UInt64, Vector{UInt64}}()
        force_rebuild = true
        
        for iter in 1:iters_for_this_thread
            if !keep_going[] break end
            
            # --- MATRIX UPDATE PHASE ---
            if force_rebuild || iter % 20 == 0
                σ_loc .= 1:n
                shuffle!(σ_loc)
                G_loc .= G[:, σ_loc]
                
                try 
                    _make_systematic_gf!(G_loc, σ_loc, k)
                    w_loc .= w_recv[σ_loc]
                    force_rebuild = false
                catch
                    force_rebuild = true
                    continue
                end
            else
                valid_swap = false
                c_in, c_out = 0, 0
                
                for _ in 1:50 
                    c_in = rand(1:k)
                    c_out = rand(k+1:n)
                    if G_loc[c_in, c_out] == 1
                        valid_swap = true
                        break
                    end
                end
                
                if !valid_swap
                    force_rebuild = true
                    continue
                end
                
                for r in 1:k
                    G_loc[r, c_in], G_loc[r, c_out] = G_loc[r, c_out], G_loc[r, c_in]
                end
                σ_loc[c_in], σ_loc[c_out] = σ_loc[c_out], σ_loc[c_in]
                w_loc[c_in], w_loc[c_out] = w_loc[c_out], w_loc[c_in]
                
                @inbounds for r in 1:k
                    if r != c_in && G_loc[r, c_in] == 1
                        for c in 1:n
                            G_loc[r, c] ⊻= G_loc[c_in, c]
                        end
                    end
                end
            end
            
            # --- 64-BIT STERN COLLISION SEARCH ---
            # 1. Determine safe window and tail sizes
            eff_l = min(l, n - k, 64)
            tail_len = n - k - eff_l

            # 2. Pack the window into UInt64 scalars
            window_rows_uint = zeros(UInt64, k)
            for i in 1:k
                bits = UInt64(0)
                for j in 1:eff_l
                    if G_loc[i, k + j] == 1
                        bits |= (UInt64(1) << (j - 1))
                    end
                end
                window_rows_uint[i] = bits
            end

            # 3. Pack the target window into a single UInt64
            # This replaces your target_window[j] -= ... loop
            target_window_uint = UInt64(0)
            for j in 1:eff_l
                if w_loc[k + j] == 1
                    target_window_uint |= (UInt64(1) << (j - 1))
                end
            end

            # Apply the contribution of the information symbols (w_loc[1:k]) to the target window
            for i in 1:k
                if w_loc[i] == 1
                    target_window_uint ⊻= window_rows_uint[i]
                end
            end

            # 4. Handle the Tail (using bit-packing for the tail is also faster)
            # If tail_len is small, you can keep it as vectors, but bit-packing is better.
            tail_rows = [G_loc[i, (k + eff_l + 1):n] for i in 1:k]
            target_tail = [w_loc[k + eff_l + j] for j in 1:tail_len]

            for i in 1:k
                if w_loc[i] == 1
                    # XOR the vectors element-wise
                    for j in 1:tail_len
                        target_tail[j] ⊻= tail_rows[i][j]
                    end
                end
            end
            
            # Helper 1: Build the hash map for the binary case (64-bit bitsliced)
            function _build_X_binary!(depth::Int, picked::Int, p::Int, 
                X_cols::UnitRange{Int}, window_rows_uint::Vector{UInt64}, 
                current_val::UInt64, msg::UInt64, 
                hash_X::Dict{UInt64, Vector{UInt64}})

                if picked == p
                    push!(get!(hash_X, current_val, Vector{UInt64}()), msg)
                    return
                end
                
                if depth > length(X_cols) || (length(X_cols) - depth + 1) < (p - picked)
                    return
                end

                # Branch 1: Row NOT included
                _build_X_binary!(depth + 1, picked, p, X_cols, window_rows_uint, 
                                current_val, msg, hash_X)

                # Branch 2: Row included (XOR the window bits)
                # We use a bitmask for the message to track which rows in X_cols were used
                new_msg = msg | (UInt64(1) << (depth - 1))
                _build_X_binary!(depth + 1, picked + 1, p, X_cols, window_rows_uint, 
                                current_val ⊻ window_rows_uint[X_cols[depth]], new_msg, hash_X)
            end
            _build_X_binary!(1, 0, p, X_cols, window_rows_uint, UInt64(0), UInt64(0), hash_X)
            
            function _probe_Y_binary!(depth::Int, picked::Int, p::Int, 
                target_w::Int, target_window_uint::UInt64, 
                target_tail::Vector{Int}, tail_rows::Vector{Vector{Int}}, 
                Y_cols::UnitRange{Int}, X_cols::UnitRange{Int}, 
                window_rows_uint::Vector{UInt64}, current_val::UInt64, msg_Y::UInt64, 
                hash_X::Dict{UInt64, Vector{UInt64}}, found_vectors::Set{Vector{Int}}, 
                results_lock::Threads.SpinLock, keep_going::Threads.Atomic{Bool}, 
                num_find::Int, n::Int, σ_loc::Vector{Int})

                if !keep_going[] return end

                if picked == p
                    # search_val is the bit-packed window we need to find in hash_X
                    search_val = target_window_uint ⊻ current_val
                    
                    if haskey(hash_X, search_val)
                        for msg_X in hash_X[search_val]
                            # Check the tail weight
                            # Start with target_tail and XOR the rows corresponding to msg_X and msg_Y
                            tail_buffer = copy(target_tail)
                            
                            # XOR rows from X part (using bitmask msg_X)
                            for i in 1:length(X_cols)
                                if (msg_X >> (i - 1)) & 1 == 1
                                    tail_buffer .⊻= tail_rows[X_cols[i]]
                                end
                            end
                            
                            # XOR rows from Y part (using bitmask msg_Y)
                            for i in 1:length(Y_cols)
                                if (msg_Y >> (i - 1)) & 1 == 1
                                    tail_buffer .⊻= tail_rows[Y_cols[i]]
                                end
                            end
                            
                            # Binary weight = wt(msg_X) + wt(msg_Y) + wt(tail)
                            # wt(msg_X) and wt(msg_Y) are both p by definition here
                            if 0 < count(!iszero, tail_buffer) + 2*p <= target_w
                                lock(results_lock) do
                                    if length(found_vectors) < num_find
                                        # Construct the error vector
                                        e_loc = zeros(Int, n)
                                        for i in 1:length(X_cols) 
                                            if (msg_X >> (i - 1)) & 1 == 1 e_loc[X_cols[i]] = 1 end 
                                        end
                                        for i in 1:length(Y_cols) 
                                            if (msg_Y >> (i - 1)) & 1 == 1 e_loc[Y_cols[i]] = 1 end 
                                        end
                                        # Tail is at the end of the systematic form
                                        tail_start = n - length(tail_buffer)
                                        for j in 1:length(tail_buffer)
                                            if tail_buffer[j] == 1 e_loc[tail_start + j] = 1 end
                                        end
                                        
                                        push!(found_vectors, e_loc[invperm(σ_loc)])
                                        if length(found_vectors) >= num_find
                                            keep_going[] = false
                                        end
                                    end
                                end
                            end
                        end
                    end
                    return
                end

                if depth > length(Y_cols) || (length(Y_cols) - depth + 1) < (p - picked)
                    return
                end

                # Branch 1: Row not picked
                _probe_Y_binary!(depth + 1, picked, p, target_w, target_window_uint, 
                                target_tail, tail_rows, Y_cols, X_cols, window_rows_uint, 
                                current_val, msg_Y, hash_X, found_vectors, results_lock, 
                                keep_going, num_find, n, σ_loc)

                # Branch 2: Row picked
                new_msg_Y = msg_Y | (UInt64(1) << (depth - 1))
                _probe_Y_binary!(depth + 1, picked + 1, p, target_w, target_window_uint, 
                                target_tail, tail_rows, Y_cols, X_cols, window_rows_uint, 
                                current_val ⊻ window_rows_uint[Y_cols[depth]], new_msg_Y, 
                                hash_X, found_vectors, results_lock, keep_going, num_find, n, σ_loc)
            end
            _probe_Y_binary!(1, 0, p, target_w, target_window_uint, 
                 target_tail, tail_rows, Y_cols, X_cols, 
                 window_rows_uint, UInt64(0), UInt64(0), 
                 hash_X, found_vectors, results_lock, 
                 keep_going, num_find, n, σ_loc)
        end
    end
    
    return found_vectors
end

#############################
        # Non-Binary
#############################

"""
    _Stern_attack_nonbinary(G::CTMatrixTypes, target_w::Int; w_recv::Vector{Int} = zeros(Int, size(G, 2)), p::Int = 2, l::Int = 3, num_find::Int = 1, max_iters::Int = 10000)

Executes Stern's Information Set Decoding (ISD) algorithm for nonbinary (GF(q)) 
block codes. This version incorporates modular field arithmetic and scalar iteration 
into the Meet-in-the-Middle search.

### Dual Modes of Operation
1. **Minimum Distance Search (Default):** If `w_recv` is the all-zero vector, 
   searches the nullspace for a codeword of weight `target_w`.
2. **Syndrome Decoding (Error Correction):** If `w_recv` is provided, searches for an 
   error vector of weight `target_w` to correct the received word.

### Arguments
* `G::CTMatrixTypes`: The generator matrix of the nonbinary linear block code over GF(q).
* `target_w::Int`: The target Hamming weight of the codeword or error vector.

### Keyword Arguments
* `w_recv::Vector{Int} = zeros(Int, size(G, 2))`: The received word. Defaults to the all-zero vector.
* `p::Int = 2`: The search weight per half of the Information Set. Keep this small, 
  as the search space grows exponentially with the field size `q`.
* `l::Int = 3`: The size of the collision window. Defaults to a much smaller value than 
  the binary version (`l=3`) because the number of possible window states is `q^l`. 
  A large `l` will stall the hash map generation.
* `num_find::Int = 1`: The number of valid vectors to find before terminating.
* `max_iters::Int = 10000`: The maximum number of random permutations to attempt.

### Returns
* `Set{Vector{Int}}`: A set containing the discovered vectors of weight `target_w`. 
  Returns an empty set if no vectors are found within `max_iters`.
"""
function _Stern_attack_nonbinary(G::CTMatrixTypes, target_w::Int; w_recv::Vector{Int} = 
    zeros(Int, size(G, 2)), p::Int = 2, l::Int = 3, num_find::Int = 1, max_iters::Int = 10000)

    # ... [Keep your existing ArgumentErrors and Setup] ...
    k, n = size(G)
    F = base_ring(G)
    p_char = Int(characteristic(F))
    d_deg = degree(F)
    non_zeros = filter(!iszero, collect(F))
    
    # ... [Keep Threading Setup] ...
    found_vectors = Set{Vector{Int}}()
    results_lock = ReentrantLock()
    Threads.@threads for th in 1:num_thrds
        # ... [Keep Thread-local buffers] ...
        for _ in 1:iters_for_this_thread
            if !keep_going[] break end
            
            # 1. Systematize & Permute
            perm = collect(1:n)
            shuffle!(perm)
            Gp = G[:, perm]
            try 
                _make_systematic_gf!(Gp, perm, k)
            catch
                continue
            end
            
            P_window = view(Gp, :, k+1 : k+l)
            P_tail   = view(Gp, :, k+l+1 : n)
            
            # 2. Compute Target Syndrome
            w_loc_F = [F(w_recv[perm[i]]) for i in 1:n]
            # [Keep your existing s_window and s_tail logic here]
            
            # 3. Hash Map for X
            hash_X = Dict{Vector{typeof(zero(F))}, Vector{Tuple{Vector{Int}, Vector{typeof(zero(F))}}}}()
            for cols in combinations(X_cols, p), scalars in Iterators.product(fill(non_zeros, p)...)
                scalar_vec = collect(scalars)
                v_X = [zero(F) for _ in 1:l]
                for i in 1:p, j in 1:l; v_X[j] += scalar_vec[i] * P_window[cols[i], j] end
                
                # Using push! into a vector to handle collisions in the hash map
                push!(get!(hash_X, v_X, []), (cols, scalar_vec))
            end
            
            # 4. Collision Search for Y
            for cols in combinations(Y_cols, p), scalars in Iterators.product(fill(non_zeros, p)...)
                if !keep_going[] break end
                scalar_vec = collect(scalars)
                v_Y = [zero(F) for _ in 1:l]
                for i in 1:p, j in 1:l; v_Y[j] += scalar_vec[i] * P_window[cols[i], j] end
                
                target_v = [s_window[j] - v_Y[j] for j in 1:l]
                
                if haskey(hash_X, target_v)
                    for (cols_X, scalars_X) in hash_X[target_v]
                        tail_wt = 0
                        is_valid = true
                        e_tail = [zero(F) for _ in 1:(n-k-l)]
                        
                        for j in 1:(n-k-l)
                            tail_val = s_tail[j]
                            for i in 1:p
                                tail_val -= (scalars_X[i] * P_tail[cols_X[i], j] + 
                                            scalar_vec[i]  * P_tail[cols[i], j])
                            end
                            e_tail[j] = tail_val
                            if !iszero(tail_val) tail_wt += 1 end
                            
                            # EARLY ABORT (Remains target_w)
                            if (2*p + tail_wt) > target_w
                                is_valid = false; break
                            end
                        end
                        
                        # CHANGE 1: THE VACUUM INEQUALITY
                        # Use <= target_w and ensure it's not the zero vector
                        total_wt = 2*p + tail_wt
                        if is_valid && 0 < total_wt <= target_w
                            
                            lock(results_lock) do
                                if length(found_vectors) < num_find
                                    e_loc = [zero(F) for _ in 1:n]
                                    for i in 1:p
                                        e_loc[cols_X[i]] = scalars_X[i]
                                        e_loc[cols[i]]   = scalar_vec[i]
                                    end
                                    for j in 1:(n-k-l); e_loc[k + l + j] = e_tail[j] end
                                    
                                    inv_p = invperm(perm)
                                    
                                    # CHANGE 2: THE BRIDGE PACKING
                                    # Use your _pack_field_elem function for extension fields!
                                    e_orig_Int = [_pack_field_elem(e_loc[inv_p[i]], p_char, d_deg) for i in 1:n]
                                    
                                    push!(found_vectors, e_orig_Int)
                                    
                                    if length(found_vectors) >= num_find
                                        keep_going[] = false
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return found_vectors
end

"""
    _Prange_attack_nonbinary(G::CTMatrixTypes, target_w::Int; w_recv::Vector{Int} = zeros(Int, size(G, 2)), num_find::Int = 1, max_iters::Int = 10000)

Executes Prange's Information Set Decoding algorithm for nonbinary block codes.
"""
function _Prange_attack_nonbinary(G::CTMatrixTypes, target_w::Int; w_recv::Vector{Int} = zeros(Int,
    size(G, 2)), num_find::Int = 1, max_iters::Int = 10000)

    k, n = size(G)
    F = base_ring(G)
    p_char = Int(characteristic(F))
    d_deg = degree(F)
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds
    remaining = max_iters % num_thrds
    
    keep_going = Threads.Atomic{Bool}(true)
    results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    is_min_dist = all(iszero, w_recv)
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            perm = collect(1:n)
            shuffle!(perm)
            G_loc = G[:, perm]
            
            try 
                _make_systematic_gf!(G_loc, perm, k) 
            catch
                continue 
            end
            
            if is_min_dist
                # MINIMUM DISTANCE MODE (p=1 Check)
                for i in 1:k
                    row_wt = count(!iszero, [G_loc[i, j] for j in 1:n])
                    if 0 < row_wt <= target_w
                        lock(results_lock) do
                            if length(found_vectors) < num_find
                                e_loc = [G_loc[i, j] for j in 1:n]
                                inv_p = invperm(perm)
                                push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                                if length(found_vectors) >= num_find
                                    keep_going[] = false
                                end
                            end
                        end
                    end
                end
            else
                # SYNDROME DECODING MODE (p=0 Check)
                w_loc_F = [F(w_recv[perm[i]]) for i in 1:n]
                S = [w_loc_F[k+j] for j in 1:(n-k)]
                wt_S = 0
                
                for j in 1:(n-k)
                    for i in 1:k
                        S[j] -= w_loc_F[i] * G_loc[i, k+j]
                    end
                    if !iszero(S[j])
                        wt_S += 1
                    end
                end
                
                if 0 < wt_S <= target_w
                    lock(results_lock) do
                        if length(found_vectors) < num_find
                            e_loc = [zero(F) for _ in 1:n]
                            e_loc[k+1:n] .= S
                            inv_p = invperm(perm)
                            # CHANGE 2: THE BRIDGE PACKING
                            # Use your _pack_field_elem function for extension fields!
                            e_orig_Int = [_pack_field_elem(e_loc[inv_p[i]], p_char, d_deg) for i in 1:n]      
                            push!(found_vectors, e_orig_Int)
                            if length(found_vectors) >= num_find
                                keep_going[] = false
                            end
                        end
                    end
                end
            end
        end
    end

    return found_vectors
end

"""
    _Lee_Brickell_attack_nonbinary(G::CTMatrixTypes, target_w::Int; w_recv::Vector{Int} = zeros(Int, size(G, 2)), p::Int = 2, num_find::Int = 1, max_iters::Int = 10000)

Executes the Lee-Brickell ISD algorithm for nonbinary block codes.
Accepts an Oscar field matrix directly for seamless GF(q) arithmetic.
"""
function _Lee_Brickell_attack_nonbinary(G::CTMatrixTypes, target_w::Int; w_recv::Vector{Int} =
    zeros(Int, size(G, 2)), p::Int = 2, num_find::Int = 1, max_iters::Int = 10000)

    k, n = size(G)
    F = base_ring(G)
    p_char = Int(characteristic(F))
    d_deg = degree(F)
    non_zeros = filter(!iszero, collect(F))
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds
    remaining = max_iters % num_thrds
    
    keep_going = Threads.Atomic{Bool}(true)
    results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            perm = collect(1:n)
            shuffle!(perm)
            G_loc = G[:, perm]
            
            try 
                _make_systematic_gf!(G_loc, perm, k) 
            catch
                continue 
            end
            
            w_loc_F = [F(w_recv[perm[i]]) for i in 1:n]
            
            # Compute Syndrome S = w_parity - w_info * P
            S = [zero(F) for _ in 1:(n-k)]
            for j in 1:(n-k)
                S[j] = w_loc_F[k+j]
                for i in 1:k
                    S[j] -= w_loc_F[i] * G_loc[i, k+j]
                end
            end
            
            # Guess e_info of weight p
            for cols in combinations(1:k, p)
                if !keep_going[] break end
                for scalars in Iterators.product(fill(non_zeros, p)...)
                    scalar_vec = collect(scalars)
                    
                    tail_wt = 0
                    is_valid = true
                    
                    e_parity = [zero(F) for _ in 1:(n-k)]
                    
                    # Compute e_parity = S + e_info * P
                    for j in 1:(n-k)
                        val = S[j]
                        for i in 1:p
                            val += scalar_vec[i] * G_loc[cols[i], k+j]
                        end
                        e_parity[j] = val
                        
                        if !iszero(val)
                            tail_wt += 1
                        end
                        
                        # Early abort if weight exceeds target
                        if p + tail_wt > target_w
                            is_valid = false
                            break
                        end
                    end
                    
                    if is_valid && (0 < p + tail_wt <= target_w)
                        lock(results_lock) do
                            if length(found_vectors) < num_find
                                e_loc = [zero(F) for _ in 1:n]
                                
                                for i in 1:p
                                    e_loc[cols[i]] = scalar_vec[i]
                                end
                                for j in 1:(n-k)
                                    e_loc[k+j] = e_parity[j]
                                end
                                
                                inv_p = invperm(perm)
                                # CHANGE 2: THE BRIDGE PACKING
                                # Use your _pack_field_elem function for extension fields!
                                e_orig_Int = [_pack_field_elem(e_loc[inv_p[i]], p_char, d_deg) for i in 1:n]
                                
                                push!(found_vectors, e_orig_Int)
                                if length(found_vectors) >= num_find
                                    keep_going[] = false
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return found_vectors
end

"""
    _Leon_attack_nonbinary(G::CTMatrixTypes, target_w::Int; w_recv::Vector{Int} = zeros(Int, size(G, 2)), p::Int = 2, l::Int = 3, num_find::Int = 1, max_iters::Int = 10000)

Executes Leon's Information Set Decoding algorithm for nonbinary block codes.
"""
function _Leon_attack_nonbinary(G::CTMatrixTypes, target_w::Int; w_recv::Vector{Int} = zeros(Int,
    size(G, 2)), p::Int = 2, l::Int = 3, num_find::Int = 1, max_iters::Int = 10000)

    k, n = size(G)
    @assert l <= (n - k) "Window size l must be less than or equal to n - k"
    
    F = base_ring(G)
    p_char = Int(characteristic(F))
    d_deg = degree(F)
    non_zeros = filter(!iszero, collect(F))
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds
    remaining = max_iters % num_thrds
    
    keep_going = Threads.Atomic{Bool}(true)
    results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            perm = collect(1:n)
            shuffle!(perm)
            G_loc = G[:, perm]
            
            try 
                _make_systematic_gf!(G_loc, perm, k) 
            catch
                continue 
            end
            
            w_loc_F = [F(w_recv[perm[i]]) for i in 1:n]
            
            # Compute S = w_parity - w_info * P
            S = [zero(F) for _ in 1:(n-k)]
            for j in 1:(n-k)
                S[j] = w_loc_F[k+j]
                for i in 1:k
                    S[j] -= w_loc_F[i] * G_loc[i, k+j]
                end
            end
            
            for cols in combinations(1:k, p)
                if !keep_going[] break end
                for scalars in Iterators.product(fill(non_zeros, p)...)
                    scalar_vec = collect(scalars)
                    
                    # 1. LEON WINDOW FILTER
                    window_zero = true
                    for j in 1:l
                        val = S[j]
                        for i in 1:p
                            val += scalar_vec[i] * G_loc[cols[i], k+j]
                        end
                        if !iszero(val)
                            window_zero = false
                            break # Fails the window check, instantly abort!
                        end
                    end
                    
                    if window_zero
                        # 2. CHECK THE REST OF THE TAIL
                        tail_wt = 0
                        is_valid = true
                        e_tail = [zero(F) for _ in (l+1):(n-k)]
                        
                        for (idx, j) in enumerate((l+1):(n-k))
                            val = S[j]
                            for i in 1:p
                                val += scalar_vec[i] * G_loc[cols[i], k+j]
                            end
                            e_tail[idx] = val
                            
                            if !iszero(val)
                                tail_wt += 1
                            end
                            if p + tail_wt > target_w
                                is_valid = false
                                break
                            end
                        end
                        
                        if is_valid && (0 < p + tail_wt <= target_w)
                            lock(results_lock) do
                                if length(found_vectors) < num_find
                                    e_loc = [zero(F) for _ in 1:n]
                                    
                                    for i in 1:p
                                        e_loc[cols[i]] = scalar_vec[i]
                                    end
                                    
                                    # We know k+1 : k+l are zero. Populate the rest.
                                    for (idx, j) in enumerate((l+1):(n-k))
                                        e_loc[k+j] = e_tail[idx]
                                    end
                                    
                                    inv_p = invperm(perm)
                                    # CHANGE 2: THE BRIDGE PACKING
                                    # Use your _pack_field_elem function for extension fields!
                                    e_orig_Int = [_pack_field_elem(e_loc[inv_p[i]], p_char, d_deg) for i in 1:n]
                                    
                                    push!(found_vectors, e_orig_Int)
                                    if length(found_vectors) >= num_find
                                        keep_going[] = false
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return found_vectors
end

"""
    _Canteaut_Chabaud_attack_nonbinary(G::CTMatrixTypes, target_w::Int; w_recv::Vector{Int} = zeros(Int, size(G_stand, 2)), p::Int = 2, l::Int = 3, num_find::Int = 1, max_iters::Int = 10000)

Executes the Canteaut-Chabaud variant of ISD for nonbinary block codes.
Uses an Oscar field matrix to handle field inversions and row reductions 
natively during the single-column pivot updates.
"""
function _Canteaut_Chabaud_attack_nonbinary(G::CTMatrixTypes, target_w::Int; w_recv::Vector{Int} =
    zeros(Int, size(G, 2)), p::Int = 2, l::Int = 3, num_find::Int = 1, max_iters::Int = 10000)

    k, n = size(G)
    F = base_ring(G)
    p_char = Int(characteristic(F))
    d_deg = degree(F)
    non_zeros = filter(!iszero, collect(F))
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds
    remaining = max_iters % num_thrds
    
    keep_going = Threads.Atomic{Bool}(true)
    results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    Threads.@threads for th in 1:num_thrds
        iters_for_this_thread = thread_load + (th <= remaining ? 1 : 0)
        
        half_k = k ÷ 2
        X_cols = 1:half_k
        Y_cols = (half_k + 1):k
        tail_len = n - k - l
        
        # Thread-local state
        σ_loc = collect(1:n)
        G_loc = zeros(F, k, n)
        w_loc = zeros(F, n)
        msg_X = zeros(F, length(X_cols))
        msg_Y = zeros(F, length(Y_cols))
        cur_win = zeros(F, l)
        hash_X = Dict{Vector{T}, Vector{Vector{T}}}()
        force_rebuild = true
        
        for iter in 1:iters_for_this_thread
            if !keep_going[] break end
            
            # --- MATRIX UPDATE PHASE ---
            if force_rebuild || iter % 20 == 0
                σ_loc .= 1:n
                shuffle!(σ_loc)
                G_loc .= G[:, σ_loc]
                
                try 
                    _make_systematic_gf!(G_loc, σ_loc, k)
                    w_loc .= [F(w_recv[σ_loc[i]]) for i in 1:n]
                    force_rebuild = false
                catch
                    force_rebuild = true
                    continue
                end
            else
                valid_swap = false
                c_in, c_out = 0, 0
                
                for _ in 1:50 
                    c_in = rand(1:k)
                    c_out = rand(k+1:n)
                    if !iszero(G_loc[c_in, c_out])
                        valid_swap = true
                        break
                    end
                end
                
                if !valid_swap
                    force_rebuild = true
                    continue
                end
                
                # 1. Swap
                for r in 1:k
                    G_loc[r, c_in], G_loc[r, c_out] = G_loc[r, c_out], G_loc[r, c_in]
                end
                σ_loc[c_in], σ_loc[c_out] = σ_loc[c_out], σ_loc[c_in]
                w_loc[c_in], w_loc[c_out] = w_loc[c_out], w_loc[c_in]
                
                # 2. Normalize Pivot
                inv_pivot = inv(G_loc[c_in, c_in])
                @inbounds for c in 1:n
                    G_loc[c_in, c] *= inv_pivot
                end
                
                # 3. Row reduction
                @inbounds for r in 1:k
                    if r != c_in
                        factor = G_loc[r, c_in]
                        if !iszero(factor)
                            for c in 1:n
                                G_loc[r, c] -= factor * G_loc[c_in, c]
                            end
                        end
                    end
                end
            end
            
            # --- NONBINARY STERN COLLISION SEARCH ---
            eff_l = min(l, n - k)
            curr_tail_len = n - k - eff_l
            # Extract rows for window (l-bits) and tail (remaining bits)
            window_rows = [G_loc[i, (k+1):(k+eff_l)] for i in 1:k]
            tail_rows   = [G_loc[i, (k+eff_l+1):n] for i in 1:k]

            # Adjust target vectors
            target_window = [w_loc[k+j] for j in 1:eff_l]
            for i in 1:k, j in 1:eff_l
                target_window[j] -= w_loc[i] * window_rows[i][j]
            end

            target_tail = [w_loc[k+eff_l+j] for j in 1:curr_tail_len]
            for i in 1:k, j in 1:curr_tail_len
                target_tail[j] -= w_loc[i] * tail_rows[i][j]
            end

            # Reset thread-local collision state
            empty!(hash_X)
            fill!(msg_X, zero(F))
            fill!(cur_win, zero(F))
            
            # Execute Search
            _build_X_recursive!(1, 0, p, l, X_cols, non_zeros, window_rows, cur_win, msg_X, hash_X)
            
            fill!(cur_win, zero(F)) # Reset window for prober
            _probe_Y_recursive!(1, 0, p, l, target_w, target_window, target_tail, tail_len, 
                                tail_rows, Y_cols, X_cols, non_zeros, window_rows, cur_win, 
                                msg_Y, hash_X, found_vectors, results_lock, keep_going, 
                                num_find, n, σ_loc, p_char, d_deg, F)
        end
    end
    
    return found_vectors
end

# Helper 1: Build the hash map from the first half of the message
function _build_X_recursive!(depth::Int, picked::Int, p::Int, l::Int, 
    X_cols::UnitRange{Int}, non_zeros::Vector{T}, window_rows::Vector{Vector{T}}, 
    current_window::Vector{T}, msg::Vector{T}, 
    hash_X::Dict{Vector{T}, Vector{Vector{T}}}) where T

    if picked == p
        # Use get! to efficiently handle the bucket logic
        push!(get!(hash_X, copy(current_window), Vector{Vector{T}}()), copy(msg))
        return
    end
    
    if depth > length(X_cols) || (length(X_cols) - depth + 1) < (p - picked)
        return
    end

    # Case 1: Scalar is zero (index depth is NOT picked)
    msg[depth] = zero(parent(non_zeros[1]))
    _build_X_recursive!(depth + 1, picked, p, l, X_cols, non_zeros, window_rows, 
                        current_window, msg, hash_X)

    # Case 2: Scalar is non-zero (index depth IS picked)
    for scalar in non_zeros
        msg[depth] = scalar
        # In-place window update (Backtracking)
        @inbounds for j in 1:l
            current_window[j] += scalar * window_rows[X_cols[depth]][j]
        end
        
        _build_X_recursive!(depth + 1, picked + 1, p, l, X_cols, non_zeros, window_rows, 
                            current_window, msg, hash_X)
        
        # Backtrack: subtract to restore state
        @inbounds for j in 1:l
            current_window[j] -= scalar * window_rows[X_cols[depth]][j]
        end
    end
    msg[depth] = zero(parent(non_zeros[1]))
end

function _probe_Y_recursive!(depth::Int, picked::Int, p::Int, l::Int, target_w::Int, 
    target_window::Vector{T}, target_tail::Vector{T}, tail_len::Int, tail_rows::Vector{Vector{T}}, 
    Y_cols::UnitRange{Int}, X_cols::UnitRange{Int}, non_zeros::Vector{T}, window_rows::Vector{Vector{T}}, 
    current_window::Vector{T}, msg_Y::Vector{T}, hash_X::Dict{Vector{T}, Vector{Vector{T}}}, 
    found_vectors::Set{Vector{Int}}, results_lock::Threads.SpinLock, keep_going::Threads.Atomic{Bool}, 
    num_find::Int, n::Int, σ_loc::Vector{Int}, p_char::Int, d_deg::Int, F::Any) where T

    if !keep_going[] return end

    if picked == p
        # Compute the value we need from hash_X
        # search_val = target_window - current_window_Y
        search_val = [target_window[j] - current_window[j] for j in 1:l]
        
        if haskey(hash_X, search_val)
            for msg_X in hash_X[search_val]
                # Check tail weight
                # We start with the target_tail and subtract the contributions of msg_X and msg_Y
                tail_buffer = copy(target_tail)
                
                # Subtract contribution from X columns
                for i in 1:length(msg_X)
                    if !iszero(msg_X[i])
                        sc = msg_X[i]
                        rows = tail_rows[X_cols[i]]
                        @inbounds for j in 1:tail_len
                            tail_buffer[j] -= sc * rows[j]
                        end
                    end
                end
                
                # Subtract contribution from Y columns
                for i in 1:length(msg_Y)
                    if !iszero(msg_Y[i])
                        sc = msg_Y[i]
                        rows = tail_rows[Y_cols[i]]
                        @inbounds for j in 1:tail_len
                            tail_buffer[j] -= sc * rows[j]
                        end
                    end
                end
                
                # Total weight = weight(msg_X) + weight(msg_Y) + weight(tail)
                # Since msg_X and msg_Y both have weight p:
                if 0 < count(!iszero, tail_buffer) + 2*p <= target_w
                    lock(results_lock) do
                        if length(found_vectors) < num_find
                            e_loc = zeros(F, n)
                            for i in 1:length(msg_X) e_loc[X_cols[i]] = msg_X[i] end
                            for i in 1:length(msg_Y) e_loc[Y_cols[i]] = msg_Y[i] end
                            # Systematic form puts tail at the very end
                            for j in 1:tail_len e_loc[n - tail_len + j] = tail_buffer[j] end
                            
                            inv_p = invperm(σ_loc)
                            e_orig_Int = [_pack_field_elem(e_loc[inv_p[idx]], p_char, d_deg) for idx in 1:n]
                            push!(found_vectors, e_orig_Int)
                            if length(found_vectors) >= num_find
                                keep_going[] = false
                            end
                        end
                    end
                end
            end
        end
        return
    end

    if depth > length(Y_cols) || (length(Y_cols) - depth + 1) < (p - picked)
        return
    end

    # Branch 1: Zero scalar
    msg_Y[depth] = zero(F)
    _probe_Y_recursive!(depth + 1, picked, p, l, target_w, target_window, target_tail, tail_len, 
                        tail_rows, Y_cols, X_cols, non_zeros, window_rows, current_window, 
                        msg_Y, hash_X, found_vectors, results_lock, keep_going, num_find, 
                        n, σ_loc, p_char, d_deg, F)

    # Branch 2: Non-zero scalars
    for scalar in non_zeros
        msg_Y[depth] = scalar
        @inbounds for j in 1:l
            current_window[j] += scalar * window_rows[Y_cols[depth]][j]
        end
        _probe_Y_recursive!(depth + 1, picked + 1, p, l, target_w, target_window, target_tail, tail_len, 
                            tail_rows, Y_cols, X_cols, non_zeros, window_rows, current_window, 
                            msg_Y, hash_X, found_vectors, results_lock, keep_going, num_find, 
                            n, σ_loc, p_char, d_deg, F)
        # Backtrack
        @inbounds for j in 1:l
            current_window[j] -= scalar * window_rows[Y_cols[depth]][j]
        end
    end
    msg_Y[depth] = zero(F)
end

#############################
         # Dispatch
#############################

"""
    Stern_attack(C::AbstractLinearCode, target_w::Int; w_recv::Vector{Int} = zeros(Int, C.n), p::Int = 2, l::Int = (Int(order(C.F)) == 2 ? 12 : 3), num_find::Int = 1, max_iters::Int = 10000)

Executes Stern's Information Set Decoding (ISD) algorithm for linear 
block codes. This function incorporate the most valuable parts of Dumer’s refinement.

### Dual Modes of Operation
1. **Minimum Distance Search (Default):** If `w_recv` is the all-zero vector, 
   searches the nullspace for a codeword of weight `target_w`.
2. **Syndrome Decoding (Error Correction):** If `w_recv` is provided, searches for an 
   error vector of weight `target_w` to correct the received word.

### Arguments
* `C::AbstractLinearCode`: The linear block code over GF(q).
* `target_w::Int`: The target Hamming weight of the codeword or error vector.

### Keyword Arguments
* `w_recv::Vector{Int} = zeros(Int, C.n)`: The received word. Defaults to the all-zero vector.
* `p::Int = 2`: The search weight per half of the Information Set. Keep this small, 
  as the search space grows exponentially with the field size `q`.
* `l::Int = 12 or 3`: The size of the collision window. For binary codes, `l = 12` perfectly balances the 
  hash map size against the false-positive collision rate. The general case defaults to a much 
  smaller value than the binary version (`l = 3`) because the number of possible window states is 
  `q^l`. A large `l` will stall the hash map generation.
* `num_find::Int = 1`: The number of valid vectors to find before terminating.
* `max_iters::Int = 10000`: The maximum number of random permutations to attempt.

### Returns
* `Set{Vector{Int}}`: A set containing the discovered vectors of weight `target_w`. 
  Returns an empty set if no vectors are found within `max_iters`.
"""
function Stern_attack(C::AbstractLinearCode, target_w::Int; w_recv::Vector{Int} = zeros(Int, C.n),
    p::Int = 2, l::Int = (Int(order(C.F)) == 2 ? 12 : 3), num_find::Int = 1, max_iters::Int = 10000)

    if Int(order(C.F)) == 2
        return _Stern_attack_binary(_convert_binary_to_int_matrix(generator_matrix(C, true)),
            target_w; w_recv = w_recv, p = p, l = l, num_find = num_find, max_iters = max_iters)
    else
        return _Stern_attack_nonbinary(generator_matrix(C, true), target_w; w_recv = w_recv, p = p,
            l = l, num_find = num_find, max_iters = max_iters)
    end
end

"""
    Prange_attack(C::AbstractLinearCode, target_w::Int; w_recv::Vector{Int} = zeros(Int, C.n), num_find::Int = 1, max_iters::Int = 10000)

Executes Prange's foundational Information Set Decoding (ISD) algorithm.

### Dual Modes of Operation
1. **Minimum Distance Search (Default):** If `w_recv` is the all-zero vector, 
   the algorithm evaluates the rows of the systematized matrix (p=1) to find 
   codewords of weight `target_w`.
2. **Syndrome Decoding:** If `w_recv` is provided, the algorithm assumes no 
   errors exist in the information set (p=0) and checks if the parity syndrome 
   has weight `target_w`.

### Arguments
* `C::AbstractLinearCode`: The linear block code.
* `target_w::Int`: The target Hamming weight of the codeword or error vector.

### Keyword Arguments
* `w_recv::Vector{Int} = zeros(Int, C.n)`: The received word. 
* `num_find::Int = 1`: The number of valid vectors to find before terminating.
* `max_iters::Int = 10000`: Maximum random permutations to attempt.

### Returns
* `Set{Vector{Int}}`: A set containing the discovered vectors. Returns an empty 
  set if `max_iters` is reached without finding matches.
"""
function Prange_attack(C::AbstractLinearCode, target_w::Int; w_recv::Vector{Int} = zeros(Int, C.n),
    num_find::Int = 1, max_iters::Int = 10000)

    if Int(order(C.F)) == 2
        return _Prange_attack_binary(_convert_binary_to_int_matrix(generator_matrix(C, true)),
            target_w; w_recv=w_recv, num_find=num_find, max_iters=max_iters)
    else
        
        return _Prange_attack_nonbinary(generator_matrix(C, true), target_w; w_recv=w_recv,
            num_find=num_find, max_iters=max_iters)
    end
end

"""
    Lee_Brickell_attack(C::AbstractLinearCode, target_w::Int; w_recv::Vector{Int} = zeros(Int, C.n), p::Int = 2, num_find::Int = 1, max_iters::Int = 10000)

Executes the Lee-Brickell Information Set Decoding (ISD) algorithm.
Automatically strips the structural objects and routes to optimized binary 
or nonbinary engines based on the code's base field.

### Dual Modes of Operation
1. **Minimum Distance Search (Default):** If `w_recv` is the all-zero vector, 
   searches the nullspace for a codeword of weight `target_w` by enforcing exactly 
   `p` non-zero elements in the Information Set.
2. **Syndrome Decoding:** If `w_recv` is provided, searches for an error vector 
   of weight `target_w` to correct the received word.

### Arguments
* `C::AbstractLinearCode`: The linear block code.
* `target_w::Int`: The target Hamming weight of the codeword or error vector.

### Keyword Arguments
* `w_recv::Vector{Int} = zeros(Int, C.n)`: The received word. 
* `p::Int = 2`: The number of errors assumed to be in the Information Set.
* `num_find::Int = 1`: The number of valid vectors to find before terminating.
* `max_iters::Int = 10000`: Maximum random permutations to attempt.

### Returns
* `Set{Vector{Int}}`: A set containing the discovered vectors. Returns an empty 
  set if `max_iters` is reached without finding matches.
"""
function Lee_Brickell_attack(C::AbstractLinearCode, target_w::Int; w_recv::Vector{Int} = zeros(Int,
    C.n), p::Int = 2, num_find::Int = 1, max_iters::Int = 10000)

    if Int(order(C.F)) == 2
        return _Lee_Brickell_attack_binary(_convert_binary_to_int_matrix(generator_matrix(C,
            true)), target_w; w_recv=w_recv, p=p, num_find=num_find, max_iters=max_iters)
    else
        return _Lee_Brickell_attack_nonbinary(generator_matrix(C, true), target_w; w_recv=w_recv,
            p=p, num_find=num_find, max_iters=max_iters)
    end
end

"""
    Leon_attack(C::AbstractLinearCode, target_w::Int; w_recv::Vector{Int} = zeros(Int, C.n),p::Int = 2, l::Int = (Int(order(C.F)) == 2 ? 12 : 3), num_find::Int = 1, max_iters::Int = 10000)

Executes Leon's Information Set Decoding (ISD) algorithm. Leon's algorithm improves 
upon Lee-Brickell by enforcing a strict filter: it demands that exactly `0` errors 
occur within a designated window of length `l`. This aggressively prunes the search 
tree before evaluating the full parity tail.

### Dual Modes of Operation
1. **Minimum Distance Search (Default):** If `w_recv` is the all-zero vector, 
   searches the nullspace for a codeword of weight `target_w`.
2. **Syndrome Decoding:** If `w_recv` is provided, searches for an error vector 
   of weight `target_w` to correct the received word.

### Arguments
* `C::AbstractLinearCode`: The linear block code.
* `target_w::Int`: The target Hamming weight of the codeword or error vector.

### Keyword Arguments
* `w_recv::Vector{Int} = zeros(Int, C.n)`: The received word. 
* `p::Int = 2`: The number of errors assumed to be in the Information Set.
* `l::Int = (Int(order(C.F)) == 2 ? 12 : 3)`: The size of the strict zero-error window. 
* `num_find::Int = 1`: The number of valid vectors to find before terminating.
* `max_iters::Int = 10000`: Maximum random permutations to attempt.

### Returns
* `Set{Vector{Int}}`: A set containing the discovered vectors. Returns an empty 
  set if `max_iters` is reached without finding matches.
"""
function Leon_attack(C::AbstractLinearCode, target_w::Int; w_recv::Vector{Int} = zeros(Int, C.n),p::Int = 2, l::Int = (Int(order(C.F)) == 2 ? 12 : 3), num_find::Int = 1, max_iters::Int = 10000)

    if Int(order(C.F)) == 2        
        return _Leon_attack_binary(_convert_binary_to_int_matrix(generator_matrix(C, true)),
            target_w; w_recv=w_recv, p=p, l=l, num_find=num_find, max_iters=max_iters)
    else
        return _Leon_attack_nonbinary(generator_matrix(C, true), target_w; w_recv=w_recv, p=p, l=l,
        num_find=num_find, max_iters=max_iters)
    end
end

"""
    Canteaut_Chabaud_attack(C::AbstractLinearCode, target_w::Int; w_recv::Vector{Int} = zeros(Int, C.n), p::Int = 2, l::Int = (Int(order(C.F)) == 2 ? 12 : 3), num_find::Int = 1, max_iters::Int = 10000)

Executes the Canteaut-Chabaud variant of Information Set Decoding (ISD). 
This algorithm modifies Stern's approach by replacing the O(n^3) full Gaussian 
elimination at each iteration with an O(n^2) single-column swap and pivot update. 
It reconstructs the full matrix entirely only after a set number of fast iterations.

### Dual Modes of Operation
1. **Minimum Distance Search (Default):** If `w_recv` is the all-zero vector, 
   searches the nullspace for a codeword of weight `target_w`.
2. **Syndrome Decoding:** If `w_recv` is provided, searches for an error vector 
   of weight `target_w` to correct the received word.

### Arguments
* `C::AbstractLinearCode`: The linear block code.
* `target_w::Int`: The target Hamming weight of the codeword or error vector.

### Keyword Arguments
* `w_recv::Vector{Int} = zeros(Int, C.n)`: The received word. 
* `p::Int = 2`: The number of errors assumed to be in *each half* of the Information Set.
* `l::Int = (Int(order(C.F)) == 2 ? 12 : 3)`: The size of the collision window. 
* `num_find::Int = 1`: The number of valid vectors to find before terminating.
* `max_iters::Int = 10000`: Maximum random permutations/pivot updates to attempt.

### Returns
* `Set{Vector{Int}}`: A set containing the discovered vectors. Returns an empty 
  set if `max_iters` is reached without finding matches.
"""
function Canteaut_Chabaud_attack(C::AbstractLinearCode, target_w::Int; w_recv::Vector{Int} =
    zeros(Int, C.n), p::Int = 2, l::Int = (Int(order(C.F)) == 2 ? 12 : 3), num_find::Int = 1,
    max_iters::Int = 10000)

    if Int(order(C.F)) == 2
        return _Canteaut_Chabaud_attack_binary(_convert_binary_to_int_matrix(
            generator_matrix(C, true)), target_w; w_recv=w_recv, p=p, l=l,num_find=num_find,
            max_iters=max_iters)
    else
        return _Canteaut_Chabaud_attack_nonbinary(generator_matrix(C, true), target_w;
            w_recv=w_recv, p=p, l=l, num_find=num_find, max_iters=max_iters)
    end
end

"""
    ISD_failure_probability(alg::Symbol, n::Int, k::Int, w::Int, iters::Int; p::Int=2, l::Int=12)

Computes the probability that a specific ISD algorithm will fail to find a vector 
of weight `w` after `iters` random permutations.
Valid algorithms: `:Prange`, `:LeeBrickell`, `:Leon`, `:Stern`, `:CanteautChabaud`.
"""
function ISD_failure_probability(alg::Symbol, n::Int, k::Int, w::Int, iters::Int; p::Int=2, l::Int=12)
    
    # Helper for BigInt binomial to prevent overflow
    bn(n, k) = (k < 0 || k > n) ? big(0) : binomial(big(n), big(k))
    
    total_space = bn(n, w)
    
    if alg == :Prange
        success_space = bn(n - k, w)
        
    elseif alg == :LeeBrickell
        success_space = bn(k, p) * bn(n - k, w - p)
        
    elseif alg == :Leon
        success_space = bn(k, p) * bn(n - k - l, w - p)
        
    elseif alg == :Stern || alg == :CanteautChabaud
        k1 = k ÷ 2
        k2 = k - k1
        success_space = bn(k1, p) * bn(k2, p) * bn(n - k - l, w - 2*p)
        
    else
        throw(ArgumentError("Unknown algorithm: $alg"))
    end
    
    # Convert back to Float64 for the probability calculation
    P_success = Float64(success_space) / Float64(total_space)
    
    if P_success == 0.0
        return 1.0 # Impossible constraints (e.g., asking for 10 errors in a 5-bit tail)
    end
    
    P_failure = (1.0 - P_success)^iters
    return P_failure
end

"""
    required_ISD_iterations(alg::Symbol, n::Int, k::Int, w::Int, target_success_rate::Float64; kwargs...)

Returns the required `max_iters` to achieve a desired success rate (e.g., 0.99 for 99%).
"""
function required_ISD_iterations(alg::Symbol, n::Int, k::Int, w::Int, target_success_rate::Float64; p::Int=2, l::Int=12)
    # Get P_success by running 1 iteration and finding the failure rate
    P_fail_1 = ISD_failure_probability(alg, n, k, w, 1; p=p, l=l)
    P_success = 1.0 - P_fail_1
    
    if P_success == 0.0
        error("Constraints make success mathematically impossible.")
    end
    
    target_failure = 1.0 - target_success_rate
    required_N = log(target_failure) / log(1.0 - P_success)
    
    return ceil(Int, required_N)
end

function _ISD_single_iteration_prob(alg::Symbol, n::Int, k::Int, w::Int; p::Int=0, l::Int=0)
    bn(n_val, k_val) = (k_val < 0 || k_val > n_val) ? big(0) : binomial(big(n_val), big(k_val))
    
    total_space = bn(n, w)
    
    if alg == :Prange
        success_space = bn(n - k, w)
    elseif alg == :LeeBrickell
        success_space = bn(k, p) * bn(n - k, w - p)
    elseif alg == :Leon
        success_space = bn(k, p) * bn(n - k - l, w - p)
    elseif alg == :Stern || alg == :CanteautChabaud
        k1 = k ÷ 2
        k2 = k - k1
        success_space = bn(k1, p) * bn(k2, p) * bn(n - k - l, w - 2*p)
    else
        throw(ArgumentError("Unknown ISD algorithm: $alg"))
    end
    
    return Float64(success_space) / Float64(total_space)
end

"""
    Gilbert_Varshamov_bound(n::Int, k::Int, q::Int)

Calculates the Gilbert-Varshamov (GV) bound for a linear code. 
This represents the expected minimum distance of a random [n, k] code over GF(q).
"""
function Gilbert_Varshamov_bound(n::Int, k::Int, q::Int)
    target_vol = big(q)^(n - k)
    vol = big(1)
    d = 1
    
    while true
        # Compute the volume of the Hamming sphere incrementally
        term = binomial(big(n - 1), big(d - 1)) * big(q - 1)^(d - 1)
        if vol + term >= target_vol
            break
        end
        vol += term
        d += 1
    end
    
    return d
end

"""
    minimum_distance_estimate(C::AbstractLinearCode; kwargs...)

Probabilistically estimates the minimum distance of a linear block code by incrementally 
searching for codewords of weight `w = start_w, ..., max_w` using Information Set Decoding.

### Keyword Arguments
* `alg::Symbol = :Stern`: The algorithm to use (`:Prange`, `:LeeBrickell`, `:Leon`, `:Stern`, `:CanteautChabaud`).
* `target_prob::Float64 = 0.99`: The statistical confidence threshold for ruling out a weight.
* `use_gv_bound::Bool = false`: If true, intelligently skips lower weights and starts the 
  search near the theoretical Gilbert-Varshamov bound to save massive computation time.
* `start_w::Int = 1`: The starting weight (overridden if `use_gv_bound = true`).
* `max_w::Int = 0`: The upper bound for the search. Defaults to the block length `C.n`.
* `p::Int = 2`: The number of errors assumed in the Information Set.
* `l::Int = (Int(order(C.F)) == 2 ? 12 : 3)`: The collision window size.

### Returns
* `Tuple{Int, Set{Vector{Int}}}`: The estimated minimum distance and the set of found codewords.
* `nothing`: If no codeword is found.
"""
function minimum_distance_estimate(C::AbstractLinearCode; 
                                   alg::Symbol = :Stern, 
                                   target_prob::Float64 = 0.99, 
                                   use_gv_bound::Bool = false,
                                   start_w::Int = 1,
                                   max_w::Int = 0,
                                   p::Int = 2,
                                   l::Int = (Int(order(C.F)) == 2 ? 12 : 3))
    
    n, k = C.n, C.k
    q = Int(order(C.F))
    upper_bound = max_w > 0 ? max_w : n
    
    # 1. OPTIONAL GV BOUND OPTIMIZATION
    if use_gv_bound
        gv = Gilbert_Varshamov_bound(n, k, q)
        # Start 2 weights below the expected GV bound to catch slight statistical deviations
        start_w = max(1, gv - 2) 
        println("GV Bound estimated at d = $gv. Starting search safely at w = $start_w.")
    end
    
    for w in start_w:upper_bound
        # 2. ITERATION CALCULATION
        iters = 0
        try
            iters = required_ISD_iterations(alg, n, k, w, target_prob; p=p, l=l)
        catch
            println("Weight $w is mathematically impossible with p=$p, l=$l. Skipping.")
            continue
        end
        
        # Guard against 0 iterations and mathematically invalid negative numbers
        if iters <= 0
            iters = 1
        end
        
        # 3. SAFETY NET FOR CRYPTOGRAPHIC SIZES
        if iters > 10^8
            @warn "Weight $w requires $iters iterations to reach $(target_prob * 100)% confidence. This may stall the CPU."
            # Uncomment the next line to strictly prevent massive iteration hangs
            # continue 
        end
        
        println("Searching for weight $w... (Running $iters iterations)")
        
        # 4. DISPATCH
        found_vectors = Set{Vector{Int}}()
        
        if alg == :Prange
            found_vectors = Prange_attack(C, w; num_find=1, max_iters=iters)
        elseif alg == :LeeBrickell
            found_vectors = Lee_Brickell_attack(C, w; p=p, num_find=1, max_iters=iters)
        elseif alg == :Leon
            found_vectors = Leon_attack(C, w; p=p, l=l, num_find=1, max_iters=iters)
        elseif alg == :Stern
            found_vectors = Stern_attack(C, w; p=p, l=l, num_find=1, max_iters=iters)
        elseif alg == :CanteautChabaud
            found_vectors = Canteaut_Chabaud_attack(C, w; p=p, l=l, num_find=1, max_iters=iters)
        else
            throw(ArgumentError("Unknown algorithm: $alg."))
        end
        
        # 5. VERIFICATION
        if !isempty(found_vectors)
            println("✅ Minimum distance established at w = $w!")
            return w, found_vectors
        end
    end
    
    println("No codeword found up to weight $upper_bound.")
    return nothing
end
