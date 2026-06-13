using SpecialFunctions # Required for loggamma

"""
    logbinomial(n::Int, k::Int)

Safely computes log(n choose k) using the Gamma function to prevent BigInt overflow.
Returns -Inf if k < 0 or k > n.
"""
function logbinomial(n::Int, k::Int)
    if k < 0 || k > n return -Inf end
    if k == 0 || k == n return 0.0 end
    return loggamma(n + 1) - loggamma(k + 1) - loggamma(n - k + 1)
end

function _prange_succ_prob(n::Int, k::Int, w::Int)
    if w > n - k return 0.0 end
    return exp(logbinomial(n - k, w) - logbinomial(n, w))
end

function _lee_brickell_succ_prob(n::Int, k::Int, w::Int, p::Int)
    if p > w || p > k || (w - p) > (n - k) return 0.0 end
    return exp(logbinomial(k, p) + logbinomial(n - k, w - p) - logbinomial(n, w))
end

function _leon_succ_prob(n::Int, k::Int, w::Int, p::Int, l::Int)
    if p > w || p > k || (w - p) > (n - k - l) return 0.0 end
    return exp(logbinomial(k, p) + logbinomial(n - k - l, w - p) - logbinomial(n, w))
end

function _stern_succ_prob(n::Int, k::Int, w::Int, p::Int, l::Int)
    if 2p > w || 2p > k || (w - 2p) > (n - k - l) return 0.0 end
    k1 = k ÷ 2
    k2 = k - k1
    if p > k1 || p > k2 return 0.0 end
    return exp(logbinomial(k1, p) + logbinomial(k2, p) + logbinomial(n - k - l, w - 2p) - logbinomial(n, w))
end

"""
    _Stern_minimum_distance_binary(G::Matrix{Int}, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)

Dedicated Minimum Distance solver using Stern's algorithm. 
Strips out all syndrome offset overhead to achieve raw nullspace collision speed.
"""
function _Stern_minimum_distance_binary(G::Matrix{Int}, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)
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
        tail_len = n - k - l
        num_tail_chunks = cld(tail_len, 64)
        
        tail_buffer = zeros(UInt64, num_tail_chunks)
        tail_rows = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
        
        for _ in 1:iters_for_this_thread
            if !keep_going[] break end
            
            # 1. Systematize
            σ_loc = collect(1:n)
            shuffle!(σ_loc)
            G_loc = G[:, σ_loc]
            
            try
                _make_systematic_gf!(G_loc, σ_loc, k)
            catch
                continue
            end
            
            # 2. 64-BIT PACKING (No targets to offset!)
            window_rows = zeros(UInt64, k)
            for i in 1:k
                val = UInt64(0)
                for j in 1:l
                    if G_loc[i, k + j] == 1 val |= (UInt64(1) << (j - 1)) end
                end
                window_rows[i] = val
            end
            
            for i in 1:k
                for j in 1:tail_len
                    if G_loc[i, k + l + j] == 1
                        tail_rows[i][(j - 1) ÷ 64 + 1] |= (UInt64(1) << ((j - 1) % 64))
                    end
                end
            end
            
            # 3. RECURSIVE HASH MAP BUILDER (X-Half)
            hash_X = Dict{UInt64, Vector{Vector{Int}}}()
            
            function _build_X!(depth, picked, current_val::UInt64, msg)
                if picked == p
                    push!(get!(hash_X, current_val, Vector{Vector{Int}}()), copy(msg))
                    return
                end
                if depth > length(X_cols) || (length(X_cols) - depth + 1) < (p - picked) return end
                
                msg[depth] = 0
                _build_X!(depth+1, picked, current_val, msg)
                
                msg[depth] = 1
                _build_X!(depth+1, picked+1, current_val ⊻ window_rows[X_cols[depth]], msg)
            end
            
            _build_X!(1, 0, UInt64(0), zeros(Int, length(X_cols)))
            
            # 4. RECURSIVE COLLISION PROBER (Y-Half)
            function _probe_Y!(depth, picked, current_val::UInt64, msg_Y)
                if !keep_going[] return end
                
                if picked == p
                    # NULLSPACE CHECK: We just look for current_val natively!
                    if haskey(hash_X, current_val)
                        for msg_X in hash_X[current_val]
                            
                            # Build the tail purely from X and Y components
                            fill!(tail_buffer, UInt64(0))
                            for i in 1:length(msg_X)
                                if msg_X[i] == 1 tail_buffer .⊻= tail_rows[X_cols[i]] end
                            end
                            for i in 1:length(msg_Y)
                                if msg_Y[i] == 1 tail_buffer .⊻= tail_rows[Y_cols[i]] end
                            end
                            
                            wt_tail = sum(count_ones.(tail_buffer))
                            
                            if 0 < wt_tail + 2*p <= target_w
                                lock(results_lock) do
                                    if length(found_vectors) < num_find
                                        e_loc = zeros(Int, n)
                                        for i in 1:length(msg_X)
                                            if msg_X[i] == 1 e_loc[X_cols[i]] = 1 end
                                        end
                                        for i in 1:length(msg_Y)
                                            if msg_Y[i] == 1 e_loc[Y_cols[i]] = 1 end
                                        end
                                        
                                        for j in 1:tail_len
                                            chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                            if (tail_buffer[chunk] & (UInt64(1) << bit)) != 0
                                                e_loc[k + l + j] = 1
                                            end
                                        end
                                        
                                        push!(found_vectors, e_loc[invperm(σ_loc)])
                                        if length(found_vectors) >= num_find keep_going[] = false end
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
    probabilistic_minimum_distance_prange(C::AbstractLinearCode; confidence::Float64 = 0.99, verbose::Bool=true)
"""
function probabilistic_minimum_distance_prange(C::AbstractLinearCode; confidence::Float64 = 0.99, verbose::Bool=true)
    k, n = C.k, C.n
    q = Int(order(C.F))
    ϵ = 1.0 - confidence
    
    G = q == 2 ? _convert_binary_to_int_matrix(generator_matrix(C, true)) : generator_matrix(C, true)
    
    for w in 1:n
        p_succ = _prange_succ_prob(n, k, w)
        if p_succ <= 0.0 continue end
        
        req_iters = p_succ >= 1.0 ? 1 : Int(ceil(log(ϵ) / log(1.0 - p_succ)))
        verbose && println("Prange (w=$w): Requires $req_iters iterations...")
        
        # ADVANCED ROUTING LOGIC
        if q == 2
            found = _Prange_minimum_distance_binary(G, w; num_find=1, max_iters=req_iters)
        elseif q == 3
            found = _Prange_minimum_distance_GF3(G, w; num_find=1, max_iters=req_iters)
        elseif q == 4
            found = _Prange_minimum_distance_GF4(G, w; num_find=1, max_iters=req_iters)
        else
            found = _Prange_minimum_distance_nonbinary(G, w; num_find=1, max_iters=req_iters)
        end
        
        if !isempty(found)
            actual_w = sum(count.(!iszero, first(found))) 
            verbose && println("🎉 Success! Found codeword of weight $actual_w.")
            return actual_w, first(found)
        end
    end
    return n, zeros(Int, n)
end

function _Prange_minimum_distance_nonbinary(G::CTMatrixTypes, target_w::Int; num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    F = base_ring(G)
    p_char, d_deg = Int(characteristic(F)), degree(F)
    
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
            
            try _make_systematic_gf!(G_loc, perm, k) catch; continue end
            
            # Nullspace search: Every row is a valid codeword
            for i in 1:k
                row_wt = count(!iszero, [G_loc[i, j] for j in 1:n])
                if 0 < row_wt <= target_w
                    lock(results_lock) do
                        if length(found_vectors) < num_find
                            e_loc = [G_loc[i, j] for j in 1:n]
                            inv_p = invperm(perm)
                            push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                            if length(found_vectors) >= num_find keep_going[] = false end
                        end
                    end
                end
            end
        end
    end
    return found_vectors
end

function _Leon_minimum_distance_nonbinary(G::CTMatrixTypes, target_w::Int; p::Int = 2, l::Int = 3, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    F = base_ring(G)
    p_char, d_deg = Int(characteristic(F)), degree(F)
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
            
            try _make_systematic_gf!(G_loc, perm, k) catch; continue end
            
            for cols in combinations(1:k, p)
                if !keep_going[] break end
                for scalars in Iterators.product(fill(non_zeros, p)...)
                    scalar_vec = collect(scalars)
                    
                    # 1. NULLSPACE LEON WINDOW FILTER (Check if sum == 0 natively)
                    window_zero = true
                    for j in 1:l
                        val = zero(F)
                        for i in 1:p
                            val += scalar_vec[i] * G_loc[cols[i], k+j]
                        end
                        if !iszero(val)
                            window_zero = false
                            break
                        end
                    end
                    
                    if window_zero
                        # 2. CHECK THE REST OF THE TAIL
                        tail_wt = 0
                        is_valid = true
                        e_tail = [zero(F) for _ in (l+1):(n-k)]
                        
                        for (idx, j) in enumerate((l+1):(n-k))
                            val = zero(F)
                            for i in 1:p
                                val += scalar_vec[i] * G_loc[cols[i], k+j]
                            end
                            e_tail[idx] = val
                            
                            if !iszero(val) tail_wt += 1 end
                            if p + tail_wt > target_w
                                is_valid = false; break
                            end
                        end
                        
                        if is_valid && (0 < p + tail_wt <= target_w)
                            lock(results_lock) do
                                if length(found_vectors) < num_find
                                    e_loc = [zero(F) for _ in 1:n]
                                    for i in 1:p e_loc[cols[i]] = scalar_vec[i] end
                                    for (idx, j) in enumerate((l+1):(n-k)) e_loc[k+j] = e_tail[idx] end
                                    
                                    inv_p = invperm(perm)
                                    e_orig_Int = [_pack_field_elem(e_loc[inv_p[i]], p_char, d_deg) for i in 1:n]
                                    
                                    push!(found_vectors, e_orig_Int)
                                    if length(found_vectors) >= num_find keep_going[] = false end
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
    _Leon_minimum_distance_GF4(G::CTMatrixTypes, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)

Dedicated Bitsliced Minimum Distance solver for GF(4) codes using Leon's algorithm.
Packs the base elements into dual UInt64 registers for SIMD parallel execution.
"""
function _Leon_minimum_distance_GF4(G::CTMatrixTypes, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    @assert l <= 64 "Window size l must be <= 64 for native register packing."
    
    F = base_ring(G)
    ω = gen(F)
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds
    remaining = max_iters % num_thrds
    
    keep_going = Threads.Atomic{Bool}(true)
    results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}() # Still returning standard Int vectors for consistency
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        tail_len = n - k - l
        num_tail_chunks = cld(tail_len, 64)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            σ_loc = collect(1:n)
            shuffle!(σ_loc)
            G_loc = G[:, σ_loc]
            
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # --- BITSLICING COMPILER ---
            # Parse the Oscar GF(4) elements into their 'a' and 'b' bit layers
            win_a = zeros(UInt64, k); win_b = zeros(UInt64, k)
            tail_a = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            tail_b = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            
            for i in 1:k
                # 1. Pack the l-window
                for j in 1:l
                    val = G_loc[i, k + j]
                    if val == 1 
                        win_b[i] |= (UInt64(1) << (j - 1))
                    elseif val == ω 
                        win_a[i] |= (UInt64(1) << (j - 1))
                    elseif val == ω + 1
                        win_a[i] |= (UInt64(1) << (j - 1))
                        win_b[i] |= (UInt64(1) << (j - 1))
                    end
                end
                
                # 2. Pack the parity tail
                for j in 1:tail_len
                    val = G_loc[i, k + l + j]
                    chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                    
                    if val == 1 
                        tail_b[i][chunk] |= (UInt64(1) << bit)
                    elseif val == ω 
                        tail_a[i][chunk] |= (UInt64(1) << bit)
                    elseif val == ω + 1
                        tail_a[i][chunk] |= (UInt64(1) << bit)
                        tail_b[i][chunk] |= (UInt64(1) << bit)
                    end
                end
            end
            
            # --- RECURSIVE HARDWARE PROBER ---
            function _probe_Leon_GF4!(depth, picked, cur_win_a::UInt64, cur_win_b::UInt64, cur_tail_a::Vector{UInt64}, cur_tail_b::Vector{UInt64}, msg)
                if !keep_going[] return end
                
                if picked == p
                    # NULLSPACE FILTER: Both register layers must be exactly 0
                    if cur_win_a == UInt64(0) && cur_win_b == UInt64(0)
                        
                        # GF(4) Weight Calculation: count_ones(a OR b)
                        wt = sum(count_ones.(cur_tail_a .| cur_tail_b))
                        
                        if 0 < wt + p <= target_w
                            lock(results_lock) do
                                if length(found_vectors) < num_find
                                    e_loc = zeros(Int, n)
                                    e_loc[1:k] .= msg
                                    
                                    # Decode the bit-layers back into standard field integers
                                    for j in 1:tail_len
                                        chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                        bit_a = (cur_tail_a[chunk] >> bit) & 1
                                        bit_b = (cur_tail_b[chunk] >> bit) & 1
                                        
                                        if bit_a == 1 && bit_b == 0 e_loc[k + l + j] = 2      # ω
                                        elseif bit_a == 1 && bit_b == 1 e_loc[k + l + j] = 3  # ω^2
                                        elseif bit_a == 0 && bit_b == 1 e_loc[k + l + j] = 1  # 1
                                        end
                                    end
                                    
                                    push!(found_vectors, e_loc[invperm(σ_loc)])
                                    if length(found_vectors) >= num_find keep_going[] = false end
                                end
                            end
                        end
                    end
                    return
                end
                
                if depth > k || (k - depth + 1) < (p - picked) return end
                
                # Branch 0 (Scalar = 0)
                msg[depth] = 0
                _probe_Leon_GF4!(depth+1, picked, cur_win_a, cur_win_b, cur_tail_a, cur_tail_b, msg)
                
                # Branch 1 (Scalar = 1) -> Add directly
                msg[depth] = 1
                _probe_Leon_GF4!(depth+1, picked+1, 
                                cur_win_a ⊻ win_a[depth], cur_win_b ⊻ win_b[depth], 
                                cur_tail_a .⊻ tail_a[depth], cur_tail_b .⊻ tail_b[depth], msg)
                
                # Branch 2 (Scalar = ω) -> Swap logic (a+b, a)
                msg[depth] = 2
                _probe_Leon_GF4!(depth+1, picked+1, 
                                cur_win_a ⊻ win_a[depth] ⊻ win_b[depth], cur_win_b ⊻ win_a[depth], 
                                cur_tail_a .⊻ tail_a[depth] .⊻ tail_b[depth], cur_tail_b .⊻ tail_a[depth], msg)
                
                # Branch 3 (Scalar = ω^2) -> Swap logic (b, a+b)
                msg[depth] = 3
                _probe_Leon_GF4!(depth+1, picked+1, 
                                cur_win_a ⊻ win_b[depth], cur_win_b ⊻ win_a[depth] ⊻ win_b[depth], 
                                cur_tail_a .⊻ tail_b[depth], cur_tail_b .⊻ tail_a[depth] .⊻ tail_b[depth], msg)
                
                msg[depth] = 0
            end
            
            # Start execution with empty dual-registers
            _probe_Leon_GF4!(1, 0, UInt64(0), UInt64(0), zeros(UInt64, num_tail_chunks), zeros(UInt64, num_tail_chunks), zeros(Int, k))
        end
    end
    return found_vectors
end

"""
    probabilistic_minimum_distance_lee_brickell(G::Matrix{Int}; confidence::Float64 = 0.99, p::Int = 2, verbose::Bool=true)

Finds the minimum distance using the Lee-Brickell algorithm.
"""
function probabilistic_minimum_distance_lee_brickell(G::Matrix{Int}; confidence::Float64 = 0.99, p::Int = 2, verbose::Bool=true)
    k, n = size(G)
    ϵ = 1.0 - confidence
    
    for w in 1:n
        p_succ = _lee_brickell_succ_prob(n, k, w, p)
        if p_succ <= 0.0 continue end
        
        req_iters = p_succ >= 1.0 ? 1 : Int(ceil(log(ϵ) / log(1.0 - p_succ)))
        verbose && println("Lee-Brickell (w=$w, p=$p): Requires $req_iters iterations...")
        
        found = _Lee_Brickell_minimum_distance_binary(G, w; p=p, num_find=1, max_iters=req_iters)
        if !isempty(found)
            actual_w = sum(first(found))
            verbose && println("🎉 Success! Found codeword of weight $actual_w.")
            return actual_w, first(found)
        end
    end
    return n, zeros(Int, n)
end

"""
    probabilistic_minimum_distance_leon(C::AbstractLinearCode; confidence::Float64 = 0.99, p::Int = 2, l::Int = (Int(order(C.F)) == 2 ? 12 : 3), verbose::Bool=true)

Finds the minimum distance using Leon's algorithm, scaling iterations to hit the target confidence level.
Automatically routes to bitsliced hardware-accelerated engines for GF(2), GF(3), and GF(4) for maximum performance.
"""
function probabilistic_minimum_distance_leon(C::AbstractLinearCode; confidence::Float64 = 0.99, p::Int = 2, l::Int = (Int(order(C.F)) == 2 ? 12 : 3), verbose::Bool=true)
    k, n = C.k, C.n
    q = Int(order(C.F))
    ϵ = 1.0 - confidence
    
    # Extract the correct matrix type
    G = q == 2 ? _convert_binary_to_int_matrix(generator_matrix(C, true)) : generator_matrix(C, true)
    
    for w in 1:n
        p_succ = _leon_succ_prob(n, k, w, p, l)
        if p_succ <= 0.0 continue end
        
        req_iters = p_succ >= 1.0 ? 1 : Int(ceil(log(ϵ) / log(1.0 - p_succ)))
        verbose && println("Leon (w=$w, p=$p, l=$l): Requires $req_iters iterations for $(confidence*100)% confidence...")
        
        # ADVANCED ROUTING LOGIC
        if q == 2
            found = _Leon_minimum_distance_binary(G, w; p=p, l=l, num_find=1, max_iters=req_iters)
        elseif q == 3
            found = _Leon_minimum_distance_GF3(G, w; p=p, l=l, num_find=1, max_iters=req_iters)
        elseif q == 4
            found = _Leon_minimum_distance_GF4(G, w; p=p, l=l, num_find=1, max_iters=req_iters)
        else
            found = _Leon_minimum_distance_nonbinary(G, w; p=p, l=l, num_find=1, max_iters=req_iters)
        end
        
        if !isempty(found)
            # Count non-zeros to accurately reflect Hamming weight for non-binary elements
            actual_w = sum(count.(!iszero, first(found))) 
            verbose && println("🎉 Success! Found codeword of weight $actual_w.")
            return actual_w, first(found)
        end
    end
    
    return n, zeros(Int, n)
end

"""
    probabilistic_minimum_distance_stern(G::Matrix{Int}; confidence::Float64 = 0.99, p::Int = 2, l::Int = 12, verbose::Bool=true)

Finds the minimum distance using Stern's algorithm.
"""
function probabilistic_minimum_distance_stern(G::Matrix{Int}; confidence::Float64 = 0.99, p::Int = 2, l::Int = 12, verbose::Bool=true)
    k, n = size(G)
    ϵ = 1.0 - confidence
    
    for w in 1:n
        p_succ = _stern_succ_prob(n, k, w, p, l)
        if p_succ <= 0.0 continue end
        
        req_iters = p_succ >= 1.0 ? 1 : Int(ceil(log(ϵ) / log(1.0 - p_succ)))
        verbose && println("Stern (w=$w, p=$p, l=$l): Requires $req_iters iterations...")
        
        found = _Stern_minimum_distance_binary(G, w; p=p, l=l, num_find=1, max_iters=req_iters)
        if !isempty(found)
            actual_w = sum(first(found))
            verbose && println("🎉 Success! Found codeword of weight $actual_w.")
            return actual_w, first(found)
        end
    end
    return n, zeros(Int, n)
end

"""
    _Prange_minimum_distance_binary(G::Matrix{Int}, target_w::Int; num_find::Int = 1, max_iters::Int = 10000)

Dedicated Minimum Distance solver using Prange's algorithm.
Strips out syndrome offset overhead to natively check the rows of the systematized matrix.
"""
function _Prange_minimum_distance_binary(G::Matrix{Int}, target_w::Int; num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
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
            
            σ_loc = collect(1:n)
            shuffle!(σ_loc)
            G_loc = G[:, σ_loc]
            
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # Nullspace search: The rows of G_sys are natively valid codewords!
            for i in 1:k
                row_wt = sum(G_loc[i, j] for j in 1:n)
                if 0 < row_wt <= target_w
                    lock(results_lock) do
                        if length(found_vectors) < num_find
                            push!(found_vectors, G_loc[i, invperm(σ_loc)])
                            if length(found_vectors) >= num_find keep_going[] = false end
                        end
                    end
                end
            end
        end
    end
    return found_vectors
end

"""
    _Lee_Brickell_minimum_distance_binary(G::Matrix{Int}, target_w::Int; p::Int = 2, num_find::Int = 1, max_iters::Int = 10000)

Dedicated Minimum Distance solver using the Lee-Brickell algorithm.
"""
function _Lee_Brickell_minimum_distance_binary(G::Matrix{Int}, target_w::Int; p::Int = 2, num_find::Int = 1, max_iters::Int = 10000)
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
            
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # PACK THE PARITY TAIL (No syndrome offsets needed!)
            tail_rows = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            for i in 1:k, j in 1:tail_len
                if G_loc[i, k + j] == 1
                    tail_rows[i][(j - 1) ÷ 64 + 1] |= (UInt64(1) << ((j - 1) % 64))
                end
            end
            
            # RECURSIVE XOR TREE
            function _probe_LB!(depth, picked, current_tail::Vector{UInt64}, msg)
                if !keep_going[] return end
                
                if picked == p
                    wt = sum(count_ones.(current_tail))
                    if 0 < wt + p <= target_w
                        lock(results_lock) do
                            if length(found_vectors) < num_find
                                e_loc = zeros(Int, n)
                                e_loc[1:k] .= msg
                                
                                for j in 1:tail_len
                                    chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                    if (current_tail[chunk] & (UInt64(1) << bit)) != 0
                                        e_loc[k + j] = 1
                                    end
                                end
                                
                                push!(found_vectors, e_loc[invperm(σ_loc)])
                                if length(found_vectors) >= num_find keep_going[] = false end
                            end
                        end
                    end
                    return
                end
                
                if depth > k || (k - depth + 1) < (p - picked) return end
                
                # Branch 0
                msg[depth] = 0
                _probe_LB!(depth+1, picked, current_tail, msg)
                
                # Branch 1
                msg[depth] = 1
                next_tail = current_tail .⊻ tail_rows[depth]
                _probe_LB!(depth+1, picked+1, next_tail, msg)
                msg[depth] = 0
            end
            
            # Start recursion with an all-zero vector (the nullspace target)
            _probe_LB!(1, 0, zeros(UInt64, num_tail_chunks), zeros(Int, k))
        end
    end
    return found_vectors
end

"""
    _Leon_minimum_distance_binary(G::Matrix{Int}, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)

Dedicated Minimum Distance solver using Leon's algorithm.
"""
function _Leon_minimum_distance_binary(G::Matrix{Int}, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)
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
            
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # PACK THE L-WINDOW
            window_rows = zeros(UInt64, k)
            for i in 1:k, j in 1:l
                if G_loc[i, k + j] == 1 window_rows[i] |= (UInt64(1) << (j - 1)) end
            end
            
            # PACK THE REMAINING TAIL
            tail_rows = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            for i in 1:k, j in 1:tail_len
                if G_loc[i, k + l + j] == 1
                    tail_rows[i][(j - 1) ÷ 64 + 1] |= (UInt64(1) << ((j - 1) % 64))
                end
            end

            function _probe_Leon!(depth, picked, current_win::UInt64, current_tail::Vector{UInt64}, msg)
                if !keep_going[] return end
                
                if picked == p
                    # NULLSPACE FILTER: Is the window exactly 0?
                    if current_win == UInt64(0)
                        wt = sum(count_ones.(current_tail))
                        
                        if 0 < wt + p <= target_w
                            lock(results_lock) do
                                if length(found_vectors) < num_find
                                    e_loc = zeros(Int, n)
                                    e_loc[1:k] .= msg
                                    
                                    for j in 1:tail_len
                                        chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                        if (current_tail[chunk] & (UInt64(1) << bit)) != 0
                                            e_loc[k + l + j] = 1
                                        end
                                    end
                                    
                                    push!(found_vectors, e_loc[invperm(σ_loc)])
                                    if length(found_vectors) >= num_find keep_going[] = false end
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
                msg[depth] = 0
            end
            
            # Start recursion with all-zero buffers
            _probe_Leon!(1, 0, UInt64(0), zeros(UInt64, num_tail_chunks), zeros(Int, k))
        end
    end
    return found_vectors
end

"""
    _Stern_minimum_distance_GF4(G::CTMatrixTypes, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)

Dedicated Bitsliced Minimum Distance solver for GF(4) codes using Stern's algorithm.
"""
function _Stern_minimum_distance_GF4(G::CTMatrixTypes, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    @assert l <= 64 "Window size l must be <= 64 for native register packing."
    
    F = base_ring(G)
    ω = gen(F)
    p_char, d_deg = Int(characteristic(F)), degree(F)
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds
    remaining = max_iters % num_thrds
    
    keep_going = Threads.Atomic{Bool}(true)
    results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        
        half_k = k ÷ 2
        X_cols = 1:half_k
        Y_cols = (half_k + 1):k
        tail_len = n - k - l
        num_tail_chunks = cld(tail_len, 64)
        
        tail_buf_a = zeros(UInt64, num_tail_chunks)
        tail_buf_b = zeros(UInt64, num_tail_chunks)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            σ_loc = collect(1:n)
            shuffle!(σ_loc)
            G_loc = G[:, σ_loc]
            
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # --- BITSLICING COMPILER ---
            win_a = zeros(UInt64, k); win_b = zeros(UInt64, k)
            tail_a = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            tail_b = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            
            for i in 1:k
                for j in 1:l
                    val = G_loc[i, k + j]
                    if val == 1 
                        win_b[i] |= (UInt64(1) << (j - 1))
                    elseif val == ω 
                        win_a[i] |= (UInt64(1) << (j - 1))
                    elseif val == ω + 1
                        win_a[i] |= (UInt64(1) << (j - 1))
                        win_b[i] |= (UInt64(1) << (j - 1))
                    end
                end
                
                for j in 1:tail_len
                    val = G_loc[i, k + l + j]
                    chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                    if val == 1 
                        tail_b[i][chunk] |= (UInt64(1) << bit)
                    elseif val == ω 
                        tail_a[i][chunk] |= (UInt64(1) << bit)
                    elseif val == ω + 1
                        tail_a[i][chunk] |= (UInt64(1) << bit)
                        tail_b[i][chunk] |= (UInt64(1) << bit)
                    end
                end
            end
            
            # --- RECURSIVE BUILDER (X-Half) ---
            hash_X = Dict{Tuple{UInt64, UInt64}, Vector{Vector{Int}}}()
            
            function _build_X_GF4!(depth, picked, cur_win_a::UInt64, cur_win_b::UInt64, msg)
                if picked == p
                    push!(get!(hash_X, (cur_win_a, cur_win_b), Vector{Vector{Int}}()), copy(msg))
                    return
                end
                if depth > length(X_cols) || (length(X_cols) - depth + 1) < (p - picked) return end
                
                idx = X_cols[depth]
                msg[depth] = 0; _build_X_GF4!(depth+1, picked, cur_win_a, cur_win_b, msg)
                
                msg[depth] = 1; _build_X_GF4!(depth+1, picked+1, cur_win_a ⊻ win_a[idx], cur_win_b ⊻ win_b[idx], msg)
                
                msg[depth] = 2; _build_X_GF4!(depth+1, picked+1, cur_win_a ⊻ win_a[idx] ⊻ win_b[idx], cur_win_b ⊻ win_a[idx], msg)
                
                msg[depth] = 3; _build_X_GF4!(depth+1, picked+1, cur_win_a ⊻ win_b[idx], cur_win_b ⊻ win_a[idx] ⊻ win_b[idx], msg)
            end
            _build_X_GF4!(1, 0, UInt64(0), UInt64(0), zeros(Int, length(X_cols)))
            
            # --- RECURSIVE PROBER (Y-Half) ---
            function _probe_Y_GF4!(depth, picked, cur_win_a::UInt64, cur_win_b::UInt64, msg_Y)
                if !keep_going[] return end
                
                if picked == p
                    if haskey(hash_X, (cur_win_a, cur_win_b))
                        for msg_X in hash_X[(cur_win_a, cur_win_b)]
                            fill!(tail_buf_a, UInt64(0)); fill!(tail_buf_b, UInt64(0))
                            
                            # Accumulate Tail for X
                            for (i, v) in enumerate(msg_X)
                                if v != 0
                                    idx = X_cols[i]
                                    if v == 1
                                        tail_buf_a .⊻= tail_a[idx]; tail_buf_b .⊻= tail_b[idx]
                                    elseif v == 2
                                        tail_buf_a .⊻= tail_a[idx] .⊻ tail_b[idx]; tail_buf_b .⊻= tail_a[idx]
                                    elseif v == 3
                                        tail_buf_a .⊻= tail_b[idx]; tail_buf_b .⊻= tail_a[idx] .⊻ tail_b[idx]
                                    end
                                end
                            end
                            
                            # Accumulate Tail for Y
                            for (i, v) in enumerate(msg_Y)
                                if v != 0
                                    idx = Y_cols[i]
                                    if v == 1
                                        tail_buf_a .⊻= tail_a[idx]; tail_buf_b .⊻= tail_b[idx]
                                    elseif v == 2
                                        tail_buf_a .⊻= tail_a[idx] .⊻ tail_b[idx]; tail_buf_b .⊻= tail_a[idx]
                                    elseif v == 3
                                        tail_buf_a .⊻= tail_b[idx]; tail_buf_b .⊻= tail_a[idx] .⊻ tail_b[idx]
                                    end
                                end
                            end
                            
                            wt_tail = sum(count_ones.(tail_buf_a .| tail_buf_b))
                            if 0 < wt_tail + 2*p <= target_w
                                lock(results_lock) do
                                    if length(found_vectors) < num_find
                                        e_loc = zeros(F, n)
                                        for (i, v) in enumerate(msg_X) if v != 0 e_loc[X_cols[i]] = v == 1 ? F(1) : (v == 2 ? ω : ω + 1) end end
                                        for (i, v) in enumerate(msg_Y) if v != 0 e_loc[Y_cols[i]] = v == 1 ? F(1) : (v == 2 ? ω : ω + 1) end end
                                        
                                        for j in 1:tail_len
                                            chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                            bit_a = (tail_buf_a[chunk] >> bit) & 1
                                            bit_b = (tail_buf_b[chunk] >> bit) & 1
                                            if bit_a == 1 && bit_b == 0 e_loc[k + l + j] = ω
                                            elseif bit_a == 1 && bit_b == 1 e_loc[k + l + j] = ω + 1
                                            elseif bit_a == 0 && bit_b == 1 e_loc[k + l + j] = F(1)
                                            end
                                        end
                                        
                                        inv_p = invperm(σ_loc)
                                        push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                                        if length(found_vectors) >= num_find keep_going[] = false end
                                    end
                                end
                            end
                        end
                    end
                    return
                end
                
                if depth > length(Y_cols) || (length(Y_cols) - depth + 1) < (p - picked) return end
                
                idx = Y_cols[depth]
                msg_Y[depth] = 0; _probe_Y_GF4!(depth+1, picked, cur_win_a, cur_win_b, msg_Y)
                
                msg_Y[depth] = 1; _probe_Y_GF4!(depth+1, picked+1, cur_win_a ⊻ win_a[idx], cur_win_b ⊻ win_b[idx], msg_Y)
                
                msg_Y[depth] = 2; _probe_Y_GF4!(depth+1, picked+1, cur_win_a ⊻ win_a[idx] ⊻ win_b[idx], cur_win_b ⊻ win_a[idx], msg_Y)
                
                msg_Y[depth] = 3; _probe_Y_GF4!(depth+1, picked+1, cur_win_a ⊻ win_b[idx], cur_win_b ⊻ win_a[idx] ⊻ win_b[idx], msg_Y)
            end
            _probe_Y_GF4!(1, 0, UInt64(0), UInt64(0), zeros(Int, length(Y_cols)))
        end
    end
    return found_vectors
end

"""
    _Stern_minimum_distance_GF3(G::CTMatrixTypes, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)

Dedicated Bitsliced Minimum Distance solver for GF(3) codes. 
Uses a native boolean logic circuit to perform parallel modular addition across 64 elements simultaneously.
"""
function _Stern_minimum_distance_GF3(G::CTMatrixTypes, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    @assert l <= 64 "Window size l must be <= 64 for native register packing."
    
    F = base_ring(G)
    p_char, d_deg = Int(characteristic(F)), degree(F)
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds
    remaining = max_iters % num_thrds
    
    keep_going = Threads.Atomic{Bool}(true)
    results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    # Inline macro for native Modulo-3 bitwise addition
    @inline function add_mod3(AH::UInt64, AL::UInt64, BH::UInt64, BL::UInt64)
        SH = ((AH ⊻ BH) & ~(AL | BL)) | (AL & BL)
        SL = ((AL ⊻ BL) & ~(AH | BH)) | (AH & BH)
        return SH, SL
    end
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        
        half_k = k ÷ 2
        X_cols = 1:half_k
        Y_cols = (half_k + 1):k
        tail_len = n - k - l
        num_tail_chunks = cld(tail_len, 64)
        
        tail_buf_H = zeros(UInt64, num_tail_chunks)
        tail_buf_L = zeros(UInt64, num_tail_chunks)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            σ_loc = collect(1:n)
            shuffle!(σ_loc)
            G_loc = G[:, σ_loc]
            
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            win_H = zeros(UInt64, k); win_L = zeros(UInt64, k)
            tail_H = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            tail_L = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            
            for i in 1:k
                for j in 1:l
                    val = G_loc[i, k + j]
                    if val == 1 win_L[i] |= (UInt64(1) << (j - 1))
                    elseif val == 2 win_H[i] |= (UInt64(1) << (j - 1))
                    end
                end
                
                for j in 1:tail_len
                    val = G_loc[i, k + l + j]
                    chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                    if val == 1 tail_L[i][chunk] |= (UInt64(1) << bit)
                    elseif val == 2 tail_H[i][chunk] |= (UInt64(1) << bit)
                    end
                end
            end
            
            hash_X = Dict{Tuple{UInt64, UInt64}, Vector{Vector{Int}}}()
            
            function _build_X_GF3!(depth, picked, cur_H::UInt64, cur_L::UInt64, msg)
                if picked == p
                    push!(get!(hash_X, (cur_H, cur_L), Vector{Vector{Int}}()), copy(msg))
                    return
                end
                if depth > length(X_cols) || (length(X_cols) - depth + 1) < (p - picked) return end
                
                idx = X_cols[depth]
                msg[depth] = 0; _build_X_GF3!(depth+1, picked, cur_H, cur_L, msg)
                
                msg[depth] = 1
                nH, nL = add_mod3(cur_H, cur_L, win_H[idx], win_L[idx])
                _build_X_GF3!(depth+1, picked+1, nH, nL, msg)
                
                msg[depth] = 2
                nH, nL = add_mod3(cur_H, cur_L, win_L[idx], win_H[idx]) # Scalar=2 swaps H and L
                _build_X_GF3!(depth+1, picked+1, nH, nL, msg)
            end
            _build_X_GF3!(1, 0, UInt64(0), UInt64(0), zeros(Int, length(X_cols)))
            
            function _probe_Y_GF3!(depth, picked, cur_H::UInt64, cur_L::UInt64, msg_Y)
                if !keep_going[] return end
                
                if picked == p
                    # GF(3) Nullspace Trick: X + Y = 0 => X = 2Y. So lookup the negation (swap H and L)!
                    if haskey(hash_X, (cur_L, cur_H))
                        for msg_X in hash_X[(cur_L, cur_H)]
                            fill!(tail_buf_H, UInt64(0)); fill!(tail_buf_L, UInt64(0))
                            
                            for (i, v) in enumerate(msg_X)
                                if v != 0
                                    idx = X_cols[i]
                                    RH = v == 1 ? tail_H[idx] : tail_L[idx]
                                    RL = v == 1 ? tail_L[idx] : tail_H[idx]
                                    for c in 1:num_tail_chunks
                                        tail_buf_H[c], tail_buf_L[c] = add_mod3(tail_buf_H[c], tail_buf_L[c], RH[c], RL[c])
                                    end
                                end
                            end
                            
                            for (i, v) in enumerate(msg_Y)
                                if v != 0
                                    idx = Y_cols[i]
                                    RH = v == 1 ? tail_H[idx] : tail_L[idx]
                                    RL = v == 1 ? tail_L[idx] : tail_H[idx]
                                    for c in 1:num_tail_chunks
                                        tail_buf_H[c], tail_buf_L[c] = add_mod3(tail_buf_H[c], tail_buf_L[c], RH[c], RL[c])
                                    end
                                end
                            end
                            
                            wt_tail = sum(count_ones.(tail_buf_H .| tail_buf_L))
                            if 0 < wt_tail + 2*p <= target_w
                                lock(results_lock) do
                                    if length(found_vectors) < num_find
                                        e_loc = zeros(F, n)
                                        for (i, v) in enumerate(msg_X) if v != 0 e_loc[X_cols[i]] = F(v) end end
                                        for (i, v) in enumerate(msg_Y) if v != 0 e_loc[Y_cols[i]] = F(v) end end
                                        
                                        for j in 1:tail_len
                                            chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                            bit_H = (tail_buf_H[chunk] >> bit) & 1
                                            bit_L = (tail_buf_L[chunk] >> bit) & 1
                                            if bit_H == 1 e_loc[k + l + j] = F(2)
                                            elseif bit_L == 1 e_loc[k + l + j] = F(1)
                                            end
                                        end
                                        
                                        inv_p = invperm(σ_loc)
                                        push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                                        if length(found_vectors) >= num_find keep_going[] = false end
                                    end
                                end
                            end
                        end
                    end
                    return
                end
                
                if depth > length(Y_cols) || (length(Y_cols) - depth + 1) < (p - picked) return end
                
                idx = Y_cols[depth]
                msg_Y[depth] = 0; _probe_Y_GF3!(depth+1, picked, cur_H, cur_L, msg_Y)
                
                msg_Y[depth] = 1
                nH, nL = add_mod3(cur_H, cur_L, win_H[idx], win_L[idx])
                _probe_Y_GF3!(depth+1, picked+1, nH, nL, msg_Y)
                
                msg_Y[depth] = 2
                nH, nL = add_mod3(cur_H, cur_L, win_L[idx], win_H[idx])
                _probe_Y_GF3!(depth+1, picked+1, nH, nL, msg_Y)
            end
            _probe_Y_GF3!(1, 0, UInt64(0), UInt64(0), zeros(Int, length(Y_cols)))
        end
    end
    return found_vectors
end

"""
    _Prange_minimum_distance_GF3(G::CTMatrixTypes, target_w::Int; num_find::Int = 1, max_iters::Int = 10000)

Dedicated Minimum Distance solver for GF(3) codes using Prange's algorithm.
Utilizes zero-allocation row scanning, exploiting the fact that scalar multiples in GF(3) 
do not alter the Hamming weight of the basis vectors.
"""
function _Prange_minimum_distance_GF3(G::CTMatrixTypes, target_w::Int; num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    F = base_ring(G)
    p_char, d_deg = Int(characteristic(F)), degree(F)
    
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
            
            σ_loc = collect(1:n)
            shuffle!(σ_loc)
            G_loc = G[:, σ_loc]
            
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # --- ZERO-ALLOCATION NULLSPACE SEARCH ---
            # We simply check the Hamming weight of the k rows directly in memory.
            for i in 1:k
                row_wt = 0
                @inbounds for j in 1:n
                    if !iszero(G_loc[i, j])
                        row_wt += 1
                    end
                end
                
                if 0 < row_wt <= target_w
                    lock(results_lock) do
                        if length(found_vectors) < num_find
                            # Only allocate the vector if we found a valid hit!
                            e_loc = zeros(F, n)
                            for j in 1:n e_loc[j] = G_loc[i, j] end
                            
                            inv_p = invperm(σ_loc)
                            push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                            
                            if length(found_vectors) >= num_find keep_going[] = false end
                        end
                    end
                end
            end
        end
    end
    return found_vectors
end

"""
    _Prange_minimum_distance_GF4(G::CTMatrixTypes, target_w::Int; num_find::Int = 1, max_iters::Int = 10000)

Dedicated Minimum Distance solver for GF(4) codes using Prange's algorithm.
Utilizes zero-allocation row scanning, bypassing the need for complex bitsliced 
addition circuits since Prange only evaluates basis vectors directly.
"""
function _Prange_minimum_distance_GF4(G::CTMatrixTypes, target_w::Int; num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    F = base_ring(G)
    p_char, d_deg = Int(characteristic(F)), degree(F) # For GF(4), p_char = 2, d_deg = 2
    
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
            
            σ_loc = collect(1:n)
            shuffle!(σ_loc)
            G_loc = G[:, σ_loc]
            
            # The Gaussian Elimination dominates the runtime here
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # --- ZERO-ALLOCATION NULLSPACE SEARCH ---
            for i in 1:k
                row_wt = 0
                @inbounds for j in 1:n
                    if !iszero(G_loc[i, j])
                        row_wt += 1
                    end
                end
                
                if 0 < row_wt <= target_w
                    lock(results_lock) do
                        if length(found_vectors) < num_find
                            # Only allocate the vector if we found a valid hit!
                            e_loc = zeros(F, n)
                            @inbounds for j in 1:n e_loc[j] = G_loc[i, j] end
                            
                            inv_p = invperm(σ_loc)
                            
                            # Pack the GF(4) elements back into standard Julia integers
                            push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                            
                            if length(found_vectors) >= num_find keep_going[] = false end
                        end
                    end
                end
            end
        end
    end
    return found_vectors
end

"""
    _Leon_minimum_distance_GF3(G::CTMatrixTypes, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)

Dedicated Bitsliced Minimum Distance solver for GF(3) codes using Leon's algorithm.
Utilizes native boolean modulo-3 logic and deferred tail-evaluation to achieve maximum performance.
"""
function _Leon_minimum_distance_GF3(G::CTMatrixTypes, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    @assert l <= 64 "Window size l must be <= 64 for native register packing."
    
    F = base_ring(G)
    p_char, d_deg = Int(characteristic(F)), degree(F)
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds
    remaining = max_iters % num_thrds
    
    keep_going = Threads.Atomic{Bool}(true)
    results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    # Inline macro for native Modulo-3 bitwise addition
    @inline function add_mod3(AH::UInt64, AL::UInt64, BH::UInt64, BL::UInt64)
        SH = ((AH ⊻ BH) & ~(AL | BL)) | (AL & BL)
        SL = ((AL ⊻ BL) & ~(AH | BH)) | (AH & BH)
        return SH, SL
    end
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        
        tail_len = n - k - l
        num_tail_chunks = cld(tail_len, 64)
        
        # Thread-local buffers for deferred tail evaluation
        eval_tail_H = zeros(UInt64, num_tail_chunks)
        eval_tail_L = zeros(UInt64, num_tail_chunks)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            σ_loc = collect(1:n)
            shuffle!(σ_loc)
            G_loc = G[:, σ_loc]
            
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # --- BITSLICING COMPILER ---
            win_H = zeros(UInt64, k); win_L = zeros(UInt64, k)
            tail_H = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            tail_L = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            
            for i in 1:k
                # Pack the L-window
                for j in 1:l
                    val = G_loc[i, k + j]
                    if val == 1 win_L[i] |= (UInt64(1) << (j - 1))
                    elseif val == 2 win_H[i] |= (UInt64(1) << (j - 1))
                    end
                end
                
                # Pack the Parity Tail
                for j in 1:tail_len
                    val = G_loc[i, k + l + j]
                    chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                    if val == 1 tail_L[i][chunk] |= (UInt64(1) << bit)
                    elseif val == 2 tail_H[i][chunk] |= (UInt64(1) << bit)
                    end
                end
            end
            
            # --- RECURSIVE PROBER (DEFERRED TAIL EVALUATION) ---
            function _probe_Leon_GF3!(depth, picked, cur_win_H::UInt64, cur_win_L::UInt64, msg)
                if !keep_going[] return end
                
                if picked == p
                    # NULLSPACE FILTER: Only proceed if the l-window is completely 0
                    if cur_win_H == UInt64(0) && cur_win_L == UInt64(0)
                        
                        # Filter passed! Assemble the tail on the fly.
                        fill!(eval_tail_H, UInt64(0))
                        fill!(eval_tail_L, UInt64(0))
                        
                        for i in 1:k
                            v = msg[i]
                            if v != 0
                                # If v=1, add normally. If v=2, swap H and L before adding.
                                RH = v == 1 ? tail_H[i] : tail_L[i]
                                RL = v == 1 ? tail_L[i] : tail_H[i]
                                
                                for c in 1:num_tail_chunks
                                    eval_tail_H[c], eval_tail_L[c] = add_mod3(eval_tail_H[c], eval_tail_L[c], RH[c], RL[c])
                                end
                            end
                        end
                        
                        wt_tail = sum(count_ones.(eval_tail_H .| eval_tail_L))
                        
                        if 0 < wt_tail + p <= target_w
                            lock(results_lock) do
                                if length(found_vectors) < num_find
                                    e_loc = zeros(F, n)
                                    for i in 1:k 
                                        if msg[i] != 0 e_loc[i] = F(msg[i]) end 
                                    end
                                    
                                    # Decode the Tail
                                    for j in 1:tail_len
                                        chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                        bit_H = (eval_tail_H[chunk] >> bit) & 1
                                        bit_L = (eval_tail_L[chunk] >> bit) & 1
                                        if bit_H == 1 e_loc[k + l + j] = F(2)
                                        elseif bit_L == 1 e_loc[k + l + j] = F(1)
                                        end
                                    end
                                    
                                    inv_p = invperm(σ_loc)
                                    push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                                    if length(found_vectors) >= num_find keep_going[] = false end
                                end
                            end
                        end
                    end
                    return
                end
                
                if depth > k || (k - depth + 1) < (p - picked) return end
                
                # Branch 0 (Scalar = 0)
                msg[depth] = 0
                _probe_Leon_GF3!(depth+1, picked, cur_win_H, cur_win_L, msg)
                
                # Branch 1 (Scalar = 1)
                msg[depth] = 1
                nH, nL = add_mod3(cur_win_H, cur_win_L, win_H[depth], win_L[depth])
                _probe_Leon_GF3!(depth+1, picked+1, nH, nL, msg)
                
                # Branch 2 (Scalar = 2 -> swap H and L of the window row)
                msg[depth] = 2
                nH, nL = add_mod3(cur_win_H, cur_win_L, win_L[depth], win_H[depth])
                _probe_Leon_GF3!(depth+1, picked+1, nH, nL, msg)
                
                msg[depth] = 0
            end
            
            # Start search carrying ONLY the 64-bit window!
            _probe_Leon_GF3!(1, 0, UInt64(0), UInt64(0), zeros(Int, k))
        end
    end
    return found_vectors
end

"""
    _Lee_Brickell_minimum_distance_GF3(G::CTMatrixTypes, target_w::Int; p::Int = 2, num_find::Int = 1, max_iters::Int = 10000)

Dedicated Bitsliced Minimum Distance solver for GF(3) codes using Lee-Brickell's algorithm.
Uses native modulo-3 boolean logic to recursively compute the parity tail without allocations.
"""
function _Lee_Brickell_minimum_distance_GF3(G::CTMatrixTypes, target_w::Int; p::Int = 2, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    F = base_ring(G)
    p_char, d_deg = Int(characteristic(F)), degree(F)
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds
    remaining = max_iters % num_thrds
    
    keep_going = Threads.Atomic{Bool}(true)
    results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    @inline function add_mod3(AH::UInt64, AL::UInt64, BH::UInt64, BL::UInt64)
        SH = ((AH ⊻ BH) & ~(AL | BL)) | (AL & BL)
        SL = ((AL ⊻ BL) & ~(AH | BH)) | (AH & BH)
        return SH, SL
    end
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        tail_len = n - k
        num_tail_chunks = cld(tail_len, 64)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            σ_loc = collect(1:n)
            shuffle!(σ_loc)
            G_loc = G[:, σ_loc]
            
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # --- BITSLICING COMPILER ---
            tail_H = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            tail_L = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            
            for i in 1:k
                for j in 1:tail_len
                    val = G_loc[i, k + j]
                    chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                    if val == 1 tail_L[i][chunk] |= (UInt64(1) << bit)
                    elseif val == 2 tail_H[i][chunk] |= (UInt64(1) << bit)
                    end
                end
            end
            
            # --- RECURSIVE PROBER ---
            function _probe_LB_GF3!(depth, picked, cur_tail_H::Vector{UInt64}, cur_tail_L::Vector{UInt64}, msg)
                if !keep_going[] return end
                
                if picked == p
                    wt = sum(count_ones.(cur_tail_H .| cur_tail_L))
                    if 0 < wt + p <= target_w
                        lock(results_lock) do
                            if length(found_vectors) < num_find
                                e_loc = zeros(F, n)
                                for i in 1:k if msg[i] != 0 e_loc[i] = F(msg[i]) end end
                                
                                for j in 1:tail_len
                                    chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                    if ((cur_tail_H[chunk] >> bit) & 1) == 1 e_loc[k + j] = F(2)
                                    elseif ((cur_tail_L[chunk] >> bit) & 1) == 1 e_loc[k + j] = F(1)
                                    end
                                end
                                
                                inv_p = invperm(σ_loc)
                                push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                                if length(found_vectors) >= num_find keep_going[] = false end
                            end
                        end
                    end
                    return
                end
                
                if depth > k || (k - depth + 1) < (p - picked) return end
                
                # Branch 0
                msg[depth] = 0
                _probe_LB_GF3!(depth+1, picked, cur_tail_H, cur_tail_L, msg)
                
                # Branch 1
                msg[depth] = 1
                nH = similar(cur_tail_H); nL = similar(cur_tail_L)
                for c in 1:num_tail_chunks nH[c], nL[c] = add_mod3(cur_tail_H[c], cur_tail_L[c], tail_H[depth][c], tail_L[depth][c]) end
                _probe_LB_GF3!(depth+1, picked+1, nH, nL, msg)
                
                # Branch 2 (Swap H and L from the matrix row)
                msg[depth] = 2
                for c in 1:num_tail_chunks nH[c], nL[c] = add_mod3(cur_tail_H[c], cur_tail_L[c], tail_L[depth][c], tail_H[depth][c]) end
                _probe_LB_GF3!(depth+1, picked+1, nH, nL, msg)
                
                msg[depth] = 0
            end
            
            _probe_LB_GF3!(1, 0, zeros(UInt64, num_tail_chunks), zeros(UInt64, num_tail_chunks), zeros(Int, k))
        end
    end
    return found_vectors
end

"""
    _Lee_Brickell_minimum_distance_GF4(G::CTMatrixTypes, target_w::Int; p::Int = 2, num_find::Int = 1, max_iters::Int = 10000)

Dedicated Bitsliced Minimum Distance solver for GF(4) codes using Lee-Brickell's algorithm.
"""
function _Lee_Brickell_minimum_distance_GF4(G::CTMatrixTypes, target_w::Int; p::Int = 2, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    F = base_ring(G)
    ω = gen(F)
    p_char, d_deg = Int(characteristic(F)), degree(F)
    
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
            
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # --- BITSLICING COMPILER ---
            tail_a = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            tail_b = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            
            for i in 1:k
                for j in 1:tail_len
                    val = G_loc[i, k + j]
                    chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                    if val == 1 tail_b[i][chunk] |= (UInt64(1) << bit)
                    elseif val == ω tail_a[i][chunk] |= (UInt64(1) << bit)
                    elseif val == ω + 1
                        tail_a[i][chunk] |= (UInt64(1) << bit)
                        tail_b[i][chunk] |= (UInt64(1) << bit)
                    end
                end
            end
            
            # --- RECURSIVE PROBER ---
            function _probe_LB_GF4!(depth, picked, cur_tail_a::Vector{UInt64}, cur_tail_b::Vector{UInt64}, msg)
                if !keep_going[] return end
                
                if picked == p
                    wt = sum(count_ones.(cur_tail_a .| cur_tail_b))
                    if 0 < wt + p <= target_w
                        lock(results_lock) do
                            if length(found_vectors) < num_find
                                e_loc = zeros(F, n)
                                for i in 1:k if msg[i] != 0 e_loc[i] = msg[i] == 1 ? F(1) : (msg[i] == 2 ? ω : ω + 1) end end
                                
                                for j in 1:tail_len
                                    chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                    bit_a = (cur_tail_a[chunk] >> bit) & 1
                                    bit_b = (cur_tail_b[chunk] >> bit) & 1
                                    if bit_a == 1 && bit_b == 0 e_loc[k + j] = ω
                                    elseif bit_a == 1 && bit_b == 1 e_loc[k + j] = ω + 1
                                    elseif bit_a == 0 && bit_b == 1 e_loc[k + j] = F(1)
                                    end
                                end
                                
                                inv_p = invperm(σ_loc)
                                push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                                if length(found_vectors) >= num_find keep_going[] = false end
                            end
                        end
                    end
                    return
                end
                
                if depth > k || (k - depth + 1) < (p - picked) return end
                
                msg[depth] = 0; _probe_LB_GF4!(depth+1, picked, cur_tail_a, cur_tail_b, msg)
                
                msg[depth] = 1; _probe_LB_GF4!(depth+1, picked+1, cur_tail_a .⊻ tail_a[depth], cur_tail_b .⊻ tail_b[depth], msg)
                
                msg[depth] = 2; _probe_LB_GF4!(depth+1, picked+1, cur_tail_a .⊻ tail_a[depth] .⊻ tail_b[depth], cur_tail_b .⊻ tail_a[depth], msg)
                
                msg[depth] = 3; _probe_LB_GF4!(depth+1, picked+1, cur_tail_a .⊻ tail_b[depth], cur_tail_b .⊻ tail_a[depth] .⊻ tail_b[depth], msg)
                
                msg[depth] = 0
            end
            
            _probe_LB_GF4!(1, 0, zeros(UInt64, num_tail_chunks), zeros(UInt64, num_tail_chunks), zeros(Int, k))
        end
    end
    return found_vectors
end

"""
    probabilistic_minimum_distance_lee_brickell(C::AbstractLinearCode; confidence::Float64 = 0.99, p::Int = 2, verbose::Bool=true)

Finds the minimum distance using the Lee-Brickell algorithm, scaling iterations to hit the target confidence level.
Automatically routes to bitsliced hardware-accelerated engines for GF(2), GF(3), and GF(4).
"""
function probabilistic_minimum_distance_lee_brickell(C::AbstractLinearCode; confidence::Float64 = 0.99, p::Int = 2, verbose::Bool=true)
    k, n = C.k, C.n
    q = Int(order(C.F))
    ϵ = 1.0 - confidence
    
    G = q == 2 ? _convert_binary_to_int_matrix(generator_matrix(C, true)) : generator_matrix(C, true)
    
    for w in 1:n
        p_succ = _lee_brickell_succ_prob(n, k, w, p)
        if p_succ <= 0.0 continue end
        
        req_iters = p_succ >= 1.0 ? 1 : Int(ceil(log(ϵ) / log(1.0 - p_succ)))
        verbose && println("Lee-Brickell (w=$w, p=$p): Requires $req_iters iterations for $(confidence*100)% confidence...")
        
        # ADVANCED ROUTING LOGIC
        if q == 2
            found = _Lee_Brickell_minimum_distance_binary(G, w; p=p, num_find=1, max_iters=req_iters)
        elseif q == 3
            found = _Lee_Brickell_minimum_distance_GF3(G, w; p=p, num_find=1, max_iters=req_iters)
        elseif q == 4
            found = _Lee_Brickell_minimum_distance_GF4(G, w; p=p, num_find=1, max_iters=req_iters)
        else
            found = _Lee_Brickell_minimum_distance_nonbinary(G, w; p=p, num_find=1, max_iters=req_iters)
        end
        
        if !isempty(found)
            actual_w = sum(count.(!iszero, first(found)))
            verbose && println("🎉 Success! Found codeword of weight $actual_w.")
            return actual_w, first(found)
        end
    end
    
    return n, zeros(Int, n)
end

"""
    _MMT_minimum_distance_binary(G::Matrix{Int}, target_w::Int; p::Int = 4, l1::Int = 8, l2::Int = 8, num_find::Int = 1, max_iters::Int = 10000)

Dedicated Minimum Distance solver using the May-Meurer-Thomae (MMT) 4-way merge algorithm.
Splits the information set into 4 quarters and performs a 2-level collision search over two stacked windows (l1 and l2).
"""
function _MMT_minimum_distance_binary(G::Matrix{Int}, target_w::Int; p::Int = 4, l1::Int = 8, l2::Int = 8, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    l = l1 + l2
    @assert l <= 64 "Combined window size (l1 + l2) must be <= 64."
    @assert p >= 4 "Target information weight p must be at least 4 for a 4-way split."
    
    # Calculate base weights per quarter
    p1 = p ÷ 4; p2 = (p - p1) ÷ 3; p3 = (p - p1 - p2) ÷ 2; p4 = p - p1 - p2 - p3
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds
    remaining = max_iters % num_thrds
    
    keep_going = Threads.Atomic{Bool}(true)
    results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        
        # Quarter bounds
        q1 = k ÷ 4; q2 = k ÷ 2; q3 = 3 * k ÷ 4
        R1 = 1:q1; R2 = (q1+1):q2; R3 = (q2+1):q3; R4 = (q3+1):k
        
        tail_len = n - k - l
        num_tail_chunks = cld(tail_len, 64)
        tail_buf = zeros(UInt64, num_tail_chunks)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            σ_loc = collect(1:n)
            shuffle!(σ_loc)
            G_loc = G[:, σ_loc]
            
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # --- PACK THE WINDOWS AND TAIL ---
            win1_rows = zeros(UInt64, k)
            win2_rows = zeros(UInt64, k)
            tail_rows = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            
            for i in 1:k
                for j in 1:l1
                    if G_loc[i, k + j] == 1 win1_rows[i] |= (UInt64(1) << (j - 1)) end
                end
                for j in 1:l2
                    if G_loc[i, k + l1 + j] == 1 win2_rows[i] |= (UInt64(1) << (j - 1)) end
                end
                for j in 1:tail_len
                    if G_loc[i, k + l + j] == 1
                        tail_rows[i][(j - 1) ÷ 64 + 1] |= (UInt64(1) << ((j - 1) % 64))
                    end
                end
            end
            
            # --- THE 4-WAY STREAMING MERGE TREE ---
            
            # H1 maps: win1 -> [(win2, cols_used)]
            H1 = Dict{UInt64, Vector{Tuple{UInt64, Vector{Int}}}}()
            # H12 maps: win2 -> [cols_used] (Since win1 is guaranteed 0 here)
            H12 = Dict{UInt64, Vector{Vector{Int}}}()
            # H3 maps: win1 -> [(win2, cols_used)]
            H3 = Dict{UInt64, Vector{Tuple{UInt64, Vector{Int}}}}()
            
            # 1. Build H1 (Base Level 1)
            function _build_H1!(depth, picked, cur_w1, cur_w2, msg)
                if picked == p1
                    push!(get!(H1, cur_w1, Tuple{UInt64, Vector{Int}}[]), (cur_w2, copy(msg)))
                    return
                end
                if depth > length(R1) || (length(R1) - depth + 1) < (p1 - picked) return end
                
                _build_H1!(depth+1, picked, cur_w1, cur_w2, msg)
                idx = R1[depth]
                push!(msg, idx)
                _build_H1!(depth+1, picked+1, cur_w1 ⊻ win1_rows[idx], cur_w2 ⊻ win2_rows[idx], msg)
                pop!(msg)
            end
            _build_H1!(1, 0, UInt64(0), UInt64(0), Int[])
            
            # 2. Probe H1 with R2 to build H12 (Level 1 Merge)
            function _build_H12!(depth, picked, cur_w1, cur_w2, msg2)
                if picked == p2
                    if haskey(H1, cur_w1) # Collision on l1!
                        for (w2_1, msg1) in H1[cur_w1]
                            # XOR win2 of both halves, combine columns
                            push!(get!(H12, cur_w2 ⊻ w2_1, Vector{Int}[]), vcat(msg1, msg2))
                        end
                    end
                    return
                end
                if depth > length(R2) || (length(R2) - depth + 1) < (p2 - picked) return end
                
                _build_H12!(depth+1, picked, cur_w1, cur_w2, msg2)
                idx = R2[depth]
                push!(msg2, idx)
                _build_H12!(depth+1, picked+1, cur_w1 ⊻ win1_rows[idx], cur_w2 ⊻ win2_rows[idx], msg2)
                pop!(msg2)
            end
            _build_H12!(1, 0, UInt64(0), UInt64(0), Int[])
            
            # Free H1 from memory early! We don't need it anymore.
            empty!(H1) 
            
            # 3. Build H3 (Base Level 3)
            function _build_H3!(depth, picked, cur_w1, cur_w2, msg)
                if picked == p3
                    push!(get!(H3, cur_w1, Tuple{UInt64, Vector{Int}}[]), (cur_w2, copy(msg)))
                    return
                end
                if depth > length(R3) || (length(R3) - depth + 1) < (p3 - picked) return end
                
                _build_H3!(depth+1, picked, cur_w1, cur_w2, msg)
                idx = R3[depth]
                push!(msg, idx)
                _build_H3!(depth+1, picked+1, cur_w1 ⊻ win1_rows[idx], cur_w2 ⊻ win2_rows[idx], msg)
                pop!(msg)
            end
            _build_H3!(1, 0, UInt64(0), UInt64(0), Int[])
            
            # 4. Probe H3 with R4 -> Probe H12 -> Tail Evaluation (Level 2 Merge)
            function _probe_H34!(depth, picked, cur_w1, cur_w2, msg4)
                if !keep_going[] return end
                
                if picked == p4
                    if haskey(H3, cur_w1) # Collision on l1!
                        for (w2_3, msg3) in H3[cur_w1]
                            w2_34 = cur_w2 ⊻ w2_3
                            
                            if haskey(H12, w2_34) # Collision on l2!
                                for msg12 in H12[w2_34]
                                    
                                    # We survived both filters. Evaluate the tail.
                                    fill!(tail_buf, UInt64(0))
                                    for c in msg12 tail_buf .⊻= tail_rows[c] end
                                    for c in msg3  tail_buf .⊻= tail_rows[c] end
                                    for c in msg4  tail_buf .⊻= tail_rows[c] end
                                    
                                    wt_tail = sum(count_ones.(tail_buf))
                                    if 0 < wt_tail + p <= target_w
                                        lock(results_lock) do
                                            if length(found_vectors) < num_find
                                                e_loc = zeros(Int, n)
                                                for c in msg12 e_loc[c] = 1 end
                                                for c in msg3  e_loc[c] = 1 end
                                                for c in msg4  e_loc[c] = 1 end
                                                
                                                for j in 1:tail_len
                                                    chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                                    if (tail_buf[chunk] & (UInt64(1) << bit)) != 0
                                                        e_loc[k + l + j] = 1
                                                    end
                                                end
                                                
                                                push!(found_vectors, e_loc[invperm(σ_loc)])
                                                if length(found_vectors) >= num_find keep_going[] = false end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    return
                end
                if depth > length(R4) || (length(R4) - depth + 1) < (p4 - picked) return end
                
                _probe_H34!(depth+1, picked, cur_w1, cur_w2, msg4)
                idx = R4[depth]
                push!(msg4, idx)
                _probe_H34!(depth+1, picked+1, cur_w1 ⊻ win1_rows[idx], cur_w2 ⊻ win2_rows[idx], msg4)
                pop!(msg4)
            end
            _probe_H34!(1, 0, UInt64(0), UInt64(0), Int[])
        end
    end
    return found_vectors
end

function _mmt_succ_prob(n::Int, k::Int, w::Int, p::Int, l1::Int, l2::Int)
    l = l1 + l2
    if p > w || p > k || (w - p) > (n - k - l) return 0.0 end
    
    # Mirror the quartering logic of the engine
    p1 = p ÷ 4; p2 = (p - p1) ÷ 3; p3 = (p - p1 - p2) ÷ 2; p4 = p - p1 - p2 - p3
    q1 = k ÷ 4; q2 = k ÷ 2; q3 = 3 * k ÷ 4
    
    s1 = q1; s2 = q2 - q1; s3 = q3 - q2; s4 = k - q3
    
    if p1 > s1 || p2 > s2 || p3 > s3 || p4 > s4 return 0.0 end
    
    num_log = logbinomial(s1, p1) + logbinomial(s2, p2) + logbinomial(s3, p3) + logbinomial(s4, p4) + logbinomial(n - k - l, w - p)
    den_log = logbinomial(n, w)
    
    return exp(num_log - den_log)
end

"""
    probabilistic_minimum_distance_mmt(C::AbstractLinearCode; confidence::Float64 = 0.99, p::Int = 4, l1::Int = (Int(order(C.F)) == 2 ? 8 : 2), l2::Int = (Int(order(C.F)) == 2 ? 8 : 2), verbose::Bool=true)

Finds the minimum distance using the May-Meurer-Thomae (MMT) 4-way merge tree. 
Automatically routes to bitsliced hardware-accelerated engines for GF(2), GF(3), and GF(4).
"""
function probabilistic_minimum_distance_mmt(C::AbstractLinearCode; confidence::Float64 = 0.99, p::Int = 4, l1::Int = (Int(order(C.F)) == 2 ? 8 : 2), l2::Int = (Int(order(C.F)) == 2 ? 8 : 2), verbose::Bool=true)
    k, n = C.k, C.n
    q = Int(order(C.F))
    ϵ = 1.0 - confidence
    
    G = q == 2 ? _convert_binary_to_int_matrix(generator_matrix(C, true)) : generator_matrix(C, true)
    
    for w in 1:n
        p_succ = _mmt_succ_prob(n, k, w, p, l1, l2)
        if p_succ <= 0.0 continue end
        
        req_iters = p_succ >= 1.0 ? 1 : Int(ceil(log(ϵ) / log(1.0 - p_succ)))
        verbose && println("MMT (w=$w, p=$p, l1=$l1, l2=$l2): Requires $req_iters iterations for $(confidence*100)% confidence...")
        
        if q == 2
            found = _MMT_minimum_distance_binary(G, w; p=p, l1=l1, l2=l2, num_find=1, max_iters=req_iters)
        elseif q == 3
            found = _MMT_minimum_distance_GF3(G, w; p=p, l1=l1, l2=l2, num_find=1, max_iters=req_iters)
        elseif q == 4
            found = _MMT_minimum_distance_GF4(G, w; p=p, l1=l1, l2=l2, num_find=1, max_iters=req_iters)
        else
            found = _MMT_minimum_distance_nonbinary(G, w; p=p, l1=l1, l2=l2, num_find=1, max_iters=req_iters)
        end
        
        if !isempty(found)
            actual_w = sum(count.(!iszero, first(found)))
            verbose && println("🎉 Success! Found codeword of weight $actual_w.")
            return actual_w, first(found)
        end
    end
    
    return n, zeros(Int, n)
end

"""
    _MMT_minimum_distance_GF4(G::CTMatrixTypes, target_w::Int; p::Int = 4, l1::Int = 8, l2::Int = 8, num_find::Int = 1, max_iters::Int = 10000)
"""
function _MMT_minimum_distance_GF4(G::CTMatrixTypes, target_w::Int; p::Int = 4, l1::Int = 8, l2::Int = 8, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    l = l1 + l2
    @assert l <= 64 "Combined window size (l1 + l2) must be <= 64."
    
    F = base_ring(G); ω = gen(F)
    p_char, d_deg = Int(characteristic(F)), degree(F)
    
    p1 = p ÷ 4; p2 = (p - p1) ÷ 3; p3 = (p - p1 - p2) ÷ 2; p4 = p - p1 - p2 - p3
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds; remaining = max_iters % num_thrds
    keep_going = Threads.Atomic{Bool}(true); results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        
        q1 = k ÷ 4; q2 = k ÷ 2; q3 = 3 * k ÷ 4
        R1 = 1:q1; R2 = (q1+1):q2; R3 = (q2+1):q3; R4 = (q3+1):k
        
        tail_len = n - k - l
        num_tail_chunks = cld(tail_len, 64)
        tail_buf_a = zeros(UInt64, num_tail_chunks)
        tail_buf_b = zeros(UInt64, num_tail_chunks)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            σ_loc = collect(1:n); shuffle!(σ_loc); G_loc = G[:, σ_loc]
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # --- BITSLICING COMPILER ---
            w1_a = zeros(UInt64, k); w1_b = zeros(UInt64, k)
            w2_a = zeros(UInt64, k); w2_b = zeros(UInt64, k)
            tail_a = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            tail_b = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            
            for i in 1:k
                for j in 1:l1
                    val = G_loc[i, k + j]
                    if val == 1 w1_b[i] |= (UInt64(1) << (j - 1))
                    elseif val == ω w1_a[i] |= (UInt64(1) << (j - 1))
                    elseif val == ω + 1 w1_a[i] |= (UInt64(1) << (j - 1)); w1_b[i] |= (UInt64(1) << (j - 1)) end
                end
                for j in 1:l2
                    val = G_loc[i, k + l1 + j]
                    if val == 1 w2_b[i] |= (UInt64(1) << (j - 1))
                    elseif val == ω w2_a[i] |= (UInt64(1) << (j - 1))
                    elseif val == ω + 1 w2_a[i] |= (UInt64(1) << (j - 1)); w2_b[i] |= (UInt64(1) << (j - 1)) end
                end
                for j in 1:tail_len
                    val = G_loc[i, k + l + j]
                    chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                    if val == 1 tail_b[i][chunk] |= (UInt64(1) << bit)
                    elseif val == ω tail_a[i][chunk] |= (UInt64(1) << bit)
                    elseif val == ω + 1 tail_a[i][chunk] |= (UInt64(1) << bit); tail_b[i][chunk] |= (UInt64(1) << bit) end
                end
            end
            
            # --- 4-WAY MERGE TREE ---
            H1 = Dict{Tuple{UInt64, UInt64}, Vector{Tuple{UInt64, UInt64, Vector{Tuple{Int, Int}}}}}()
            H12 = Dict{Tuple{UInt64, UInt64}, Vector{Vector{Tuple{Int, Int}}}}()
            H3 = Dict{Tuple{UInt64, UInt64}, Vector{Tuple{UInt64, UInt64, Vector{Tuple{Int, Int}}}}}()
            
            function _build_base!(depth, picked, Range, p_tgt, cur_w1a, cur_w1b, cur_w2a, cur_w2b, msg, target_Dict)
                if picked == p_tgt
                    push!(get!(target_Dict, (cur_w1a, cur_w1b), []), (cur_w2a, cur_w2b, copy(msg)))
                    return
                end
                if depth > length(Range) || (length(Range) - depth + 1) < (p_tgt - picked) return end
                
                _build_base!(depth+1, picked, Range, p_tgt, cur_w1a, cur_w1b, cur_w2a, cur_w2b, msg, target_Dict)
                
                idx = Range[depth]
                push!(msg, (idx, 1)); _build_base!(depth+1, picked+1, Range, p_tgt, cur_w1a ⊻ w1_a[idx], cur_w1b ⊻ w1_b[idx], cur_w2a ⊻ w2_a[idx], cur_w2b ⊻ w2_b[idx], msg, target_Dict); pop!(msg)
                
                push!(msg, (idx, 2)); _build_base!(depth+1, picked+1, Range, p_tgt, cur_w1a ⊻ w1_a[idx] ⊻ w1_b[idx], cur_w1b ⊻ w1_a[idx], cur_w2a ⊻ w2_a[idx] ⊻ w2_b[idx], cur_w2b ⊻ w2_a[idx], msg, target_Dict); pop!(msg)
                
                push!(msg, (idx, 3)); _build_base!(depth+1, picked+1, Range, p_tgt, cur_w1a ⊻ w1_b[idx], cur_w1b ⊻ w1_a[idx] ⊻ w1_b[idx], cur_w2a ⊻ w2_b[idx], cur_w2b ⊻ w2_a[idx] ⊻ w2_b[idx], msg, target_Dict); pop!(msg)
            end
            
            _build_base!(1, 0, R1, p1, UInt64(0), UInt64(0), UInt64(0), UInt64(0), Tuple{Int, Int}[], H1)
            
            function _build_H12!(depth, picked, cur_w1a, cur_w1b, cur_w2a, cur_w2b, msg2)
                if picked == p2
                    if haskey(H1, (cur_w1a, cur_w1b))
                        for (w2a_1, w2b_1, msg1) in H1[(cur_w1a, cur_w1b)]
                            push!(get!(H12, (cur_w2a ⊻ w2a_1, cur_w2b ⊻ w2b_1), []), vcat(msg1, msg2))
                        end
                    end
                    return
                end
                if depth > length(R2) || (length(R2) - depth + 1) < (p2 - picked) return end
                
                _build_H12!(depth+1, picked, cur_w1a, cur_w1b, cur_w2a, cur_w2b, msg2)
                idx = R2[depth]
                push!(msg2, (idx, 1)); _build_H12!(depth+1, picked+1, cur_w1a ⊻ w1_a[idx], cur_w1b ⊻ w1_b[idx], cur_w2a ⊻ w2_a[idx], cur_w2b ⊻ w2_b[idx], msg2); pop!(msg2)
                push!(msg2, (idx, 2)); _build_H12!(depth+1, picked+1, cur_w1a ⊻ w1_a[idx] ⊻ w1_b[idx], cur_w1b ⊻ w1_a[idx], cur_w2a ⊻ w2_a[idx] ⊻ w2_b[idx], cur_w2b ⊻ w2_a[idx], msg2); pop!(msg2)
                push!(msg2, (idx, 3)); _build_H12!(depth+1, picked+1, cur_w1a ⊻ w1_b[idx], cur_w1b ⊻ w1_a[idx] ⊻ w1_b[idx], cur_w2a ⊻ w2_b[idx], cur_w2b ⊻ w2_a[idx] ⊻ w2_b[idx], msg2); pop!(msg2)
            end
            _build_H12!(1, 0, UInt64(0), UInt64(0), UInt64(0), UInt64(0), Tuple{Int, Int}[])
            empty!(H1) # Free Memory!
            
            _build_base!(1, 0, R3, p3, UInt64(0), UInt64(0), UInt64(0), UInt64(0), Tuple{Int, Int}[], H3)
            
            function _probe_H34!(depth, picked, cur_w1a, cur_w1b, cur_w2a, cur_w2b, msg4)
                if !keep_going[] return end
                if picked == p4
                    if haskey(H3, (cur_w1a, cur_w1b)) # Collision l1
                        for (w2a_3, w2b_3, msg3) in H3[(cur_w1a, cur_w1b)]
                            w2a_34, w2b_34 = cur_w2a ⊻ w2a_3, cur_w2b ⊻ w2b_3
                            
                            if haskey(H12, (w2a_34, w2b_34)) # Collision l2
                                for msg12 in H12[(w2a_34, w2b_34)]
                                    # Tail Evaluation
                                    fill!(tail_buf_a, UInt64(0)); fill!(tail_buf_b, UInt64(0))
                                    for sub_msg in (msg12, msg3, msg4)
                                        for (idx, v) in sub_msg
                                            if v == 1 tail_buf_a .⊻= tail_a[idx]; tail_buf_b .⊻= tail_b[idx]
                                            elseif v == 2 tail_buf_a .⊻= tail_a[idx] .⊻ tail_b[idx]; tail_buf_b .⊻= tail_a[idx]
                                            elseif v == 3 tail_buf_a .⊻= tail_b[idx]; tail_buf_b .⊻= tail_a[idx] .⊻ tail_b[idx]
                                            end
                                        end
                                    end
                                    
                                    wt_tail = sum(count_ones.(tail_buf_a .| tail_buf_b))
                                    if 0 < wt_tail + p <= target_w
                                        lock(results_lock) do
                                            if length(found_vectors) < num_find
                                                e_loc = zeros(F, n)
                                                for sub_msg in (msg12, msg3, msg4)
                                                    for (idx, v) in sub_msg e_loc[idx] = v == 1 ? F(1) : (v == 2 ? ω : ω + 1) end
                                                end
                                                for j in 1:tail_len
                                                    chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                                    bit_a = (tail_buf_a[chunk] >> bit) & 1
                                                    bit_b = (tail_buf_b[chunk] >> bit) & 1
                                                    if bit_a == 1 && bit_b == 0 e_loc[k + l + j] = ω
                                                    elseif bit_a == 1 && bit_b == 1 e_loc[k + l + j] = ω + 1
                                                    elseif bit_a == 0 && bit_b == 1 e_loc[k + l + j] = F(1)
                                                    end
                                                end
                                                inv_p = invperm(σ_loc)
                                                push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                                                if length(found_vectors) >= num_find keep_going[] = false end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    return
                end
                if depth > length(R4) || (length(R4) - depth + 1) < (p4 - picked) return end
                
                _probe_H34!(depth+1, picked, cur_w1a, cur_w1b, cur_w2a, cur_w2b, msg4)
                idx = R4[depth]
                push!(msg4, (idx, 1)); _probe_H34!(depth+1, picked+1, cur_w1a ⊻ w1_a[idx], cur_w1b ⊻ w1_b[idx], cur_w2a ⊻ w2_a[idx], cur_w2b ⊻ w2_b[idx], msg4); pop!(msg4)
                push!(msg4, (idx, 2)); _probe_H34!(depth+1, picked+1, cur_w1a ⊻ w1_a[idx] ⊻ w1_b[idx], cur_w1b ⊻ w1_a[idx], cur_w2a ⊻ w2_a[idx] ⊻ w2_b[idx], cur_w2b ⊻ w2_a[idx], msg4); pop!(msg4)
                push!(msg4, (idx, 3)); _probe_H34!(depth+1, picked+1, cur_w1a ⊻ w1_b[idx], cur_w1b ⊻ w1_a[idx] ⊻ w1_b[idx], cur_w2a ⊻ w2_b[idx], cur_w2b ⊻ w2_a[idx] ⊻ w2_b[idx], msg4); pop!(msg4)
            end
            _probe_H34!(1, 0, UInt64(0), UInt64(0), UInt64(0), UInt64(0), Tuple{Int, Int}[])
        end
    end
    return found_vectors
end

function _MMT_minimum_distance_GF3(G::CTMatrixTypes, target_w::Int; p::Int = 4, l1::Int = 8, l2::Int = 8, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    l = l1 + l2
    
    F = base_ring(G); p_char, d_deg = Int(characteristic(F)), degree(F)
    p1 = p ÷ 4; p2 = (p - p1) ÷ 3; p3 = (p - p1 - p2) ÷ 2; p4 = p - p1 - p2 - p3
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds; remaining = max_iters % num_thrds
    keep_going = Threads.Atomic{Bool}(true); results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    @inline function add_mod3(AH::UInt64, AL::UInt64, BH::UInt64, BL::UInt64)
        SH = ((AH ⊻ BH) & ~(AL | BL)) | (AL & BL)
        SL = ((AL ⊻ BL) & ~(AH | BH)) | (AH & BH)
        return SH, SL
    end
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        
        q1 = k ÷ 4; q2 = k ÷ 2; q3 = 3 * k ÷ 4
        R1 = 1:q1; R2 = (q1+1):q2; R3 = (q2+1):q3; R4 = (q3+1):k
        tail_len = n - k - l; num_tail_chunks = cld(tail_len, 64)
        tail_buf_H = zeros(UInt64, num_tail_chunks); tail_buf_L = zeros(UInt64, num_tail_chunks)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            σ_loc = collect(1:n); shuffle!(σ_loc); G_loc = G[:, σ_loc]
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # --- BITSLICING COMPILER ---
            w1_H = zeros(UInt64, k); w1_L = zeros(UInt64, k)
            w2_H = zeros(UInt64, k); w2_L = zeros(UInt64, k)
            tail_H = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            tail_L = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            
            for i in 1:k
                for j in 1:l1
                    v = G_loc[i, k + j]
                    if v == 1 w1_L[i] |= (UInt64(1) << (j - 1)) elseif v == 2 w1_H[i] |= (UInt64(1) << (j - 1)) end
                end
                for j in 1:l2
                    v = G_loc[i, k + l1 + j]
                    if v == 1 w2_L[i] |= (UInt64(1) << (j - 1)) elseif v == 2 w2_H[i] |= (UInt64(1) << (j - 1)) end
                end
                for j in 1:tail_len
                    v = G_loc[i, k + l + j]; chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                    if v == 1 tail_L[i][chunk] |= (UInt64(1) << bit) elseif v == 2 tail_H[i][chunk] |= (UInt64(1) << bit) end
                end
            end
            
            # --- 4-WAY MERGE TREE ---
            H1 = Dict{Tuple{UInt64, UInt64}, Vector{Tuple{UInt64, UInt64, Vector{Tuple{Int, Int}}}}}()
            H12 = Dict{Tuple{UInt64, UInt64}, Vector{Vector{Tuple{Int, Int}}}}()
            H3 = Dict{Tuple{UInt64, UInt64}, Vector{Tuple{UInt64, UInt64, Vector{Tuple{Int, Int}}}}}()
            
            function _build_base!(depth, picked, Range, p_tgt, cw1_H, cw1_L, cw2_H, cw2_L, msg, tgt_Dict)
                if picked == p_tgt
                    push!(get!(tgt_Dict, (cw1_H, cw1_L), []), (cw2_H, cw2_L, copy(msg)))
                    return
                end
                if depth > length(Range) || (length(Range) - depth + 1) < (p_tgt - picked) return end
                
                _build_base!(depth+1, picked, Range, p_tgt, cw1_H, cw1_L, cw2_H, cw2_L, msg, tgt_Dict)
                
                idx = Range[depth]
                n1_H, n1_L = add_mod3(cw1_H, cw1_L, w1_H[idx], w1_L[idx])
                n2_H, n2_L = add_mod3(cw2_H, cw2_L, w2_H[idx], w2_L[idx])
                push!(msg, (idx, 1)); _build_base!(depth+1, picked+1, Range, p_tgt, n1_H, n1_L, n2_H, n2_L, msg, tgt_Dict); pop!(msg)
                
                n1_H, n1_L = add_mod3(cw1_H, cw1_L, w1_L[idx], w1_H[idx])
                n2_H, n2_L = add_mod3(cw2_H, cw2_L, w2_L[idx], w2_H[idx])
                push!(msg, (idx, 2)); _build_base!(depth+1, picked+1, Range, p_tgt, n1_H, n1_L, n2_H, n2_L, msg, tgt_Dict); pop!(msg)
            end
            
            _build_base!(1, 0, R1, p1, UInt64(0), UInt64(0), UInt64(0), UInt64(0), Tuple{Int, Int}[], H1)
            
            function _build_H12!(depth, picked, cw1_H, cw1_L, cw2_H, cw2_L, msg2)
                if picked == p2
                    # GF(3) Negation Lookup: lookup (L, H)
                    if haskey(H1, (cw1_L, cw1_H))
                        for (w2H_1, w2L_1, msg1) in H1[(cw1_L, cw1_H)]
                            nH, nL = add_mod3(cw2_H, cw2_L, w2H_1, w2L_1)
                            push!(get!(H12, (nH, nL), []), vcat(msg1, msg2))
                        end
                    end
                    return
                end
                if depth > length(R2) || (length(R2) - depth + 1) < (p2 - picked) return end
                
                _build_H12!(depth+1, picked, cw1_H, cw1_L, cw2_H, cw2_L, msg2)
                idx = R2[depth]
                n1_H, n1_L = add_mod3(cw1_H, cw1_L, w1_H[idx], w1_L[idx])
                n2_H, n2_L = add_mod3(cw2_H, cw2_L, w2_H[idx], w2_L[idx])
                push!(msg2, (idx, 1)); _build_H12!(depth+1, picked+1, n1_H, n1_L, n2_H, n2_L, msg2); pop!(msg2)
                
                n1_H, n1_L = add_mod3(cw1_H, cw1_L, w1_L[idx], w1_H[idx])
                n2_H, n2_L = add_mod3(cw2_H, cw2_L, w2_L[idx], w2_H[idx])
                push!(msg2, (idx, 2)); _build_H12!(depth+1, picked+1, n1_H, n1_L, n2_H, n2_L, msg2); pop!(msg2)
            end
            _build_H12!(1, 0, UInt64(0), UInt64(0), UInt64(0), UInt64(0), Tuple{Int, Int}[])
            empty!(H1)
            
            _build_base!(1, 0, R3, p3, UInt64(0), UInt64(0), UInt64(0), UInt64(0), Tuple{Int, Int}[], H3)
            
            function _probe_H34!(depth, picked, cw1_H, cw1_L, cw2_H, cw2_L, msg4)
                if !keep_going[] return end
                if picked == p4
                    if haskey(H3, (cw1_L, cw1_H)) # Collision l1 lookup negation
                        for (w2H_3, w2L_3, msg3) in H3[(cw1_L, cw1_H)]
                            w2H_34, w2L_34 = add_mod3(cw2_H, cw2_L, w2H_3, w2L_3)
                            
                            if haskey(H12, (w2L_34, w2H_34)) # Collision l2 lookup negation
                                for msg12 in H12[(w2L_34, w2H_34)]
                                    fill!(tail_buf_H, UInt64(0)); fill!(tail_buf_L, UInt64(0))
                                    for sub_msg in (msg12, msg3, msg4)
                                        for (idx, v) in sub_msg
                                            RH = v == 1 ? tail_H[idx] : tail_L[idx]
                                            RL = v == 1 ? tail_L[idx] : tail_H[idx]
                                            for c in 1:num_tail_chunks
                                                tail_buf_H[c], tail_buf_L[c] = add_mod3(tail_buf_H[c], tail_buf_L[c], RH[c], RL[c])
                                            end
                                        end
                                    end
                                    
                                    wt_tail = sum(count_ones.(tail_buf_H .| tail_buf_L))
                                    if 0 < wt_tail + p <= target_w
                                        lock(results_lock) do
                                            if length(found_vectors) < num_find
                                                e_loc = zeros(F, n)
                                                for sub_msg in (msg12, msg3, msg4)
                                                    for (idx, v) in sub_msg e_loc[idx] = F(v) end
                                                end
                                                for j in 1:tail_len
                                                    chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                                    bit_H = (tail_buf_H[chunk] >> bit) & 1
                                                    bit_L = (tail_buf_L[chunk] >> bit) & 1
                                                    if bit_H == 1 e_loc[k + l + j] = F(2)
                                                    elseif bit_L == 1 e_loc[k + l + j] = F(1) end
                                                end
                                                inv_p = invperm(σ_loc)
                                                push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                                                if length(found_vectors) >= num_find keep_going[] = false end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    return
                end
                if depth > length(R4) || (length(R4) - depth + 1) < (p4 - picked) return end
                
                _probe_H34!(depth+1, picked, cw1_H, cw1_L, cw2_H, cw2_L, msg4)
                idx = R4[depth]
                n1_H, n1_L = add_mod3(cw1_H, cw1_L, w1_H[idx], w1_L[idx]); n2_H, n2_L = add_mod3(cw2_H, cw2_L, w2_H[idx], w2_L[idx])
                push!(msg4, (idx, 1)); _probe_H34!(depth+1, picked+1, n1_H, n1_L, n2_H, n2_L, msg4); pop!(msg4)
                
                n1_H, n1_L = add_mod3(cw1_H, cw1_L, w1_L[idx], w1_H[idx]); n2_H, n2_L = add_mod3(cw2_H, cw2_L, w2_L[idx], w2_H[idx])
                push!(msg4, (idx, 2)); _probe_H34!(depth+1, picked+1, n1_H, n1_L, n2_H, n2_L, msg4); pop!(msg4)
            end
            _probe_H34!(1, 0, UInt64(0), UInt64(0), UInt64(0), UInt64(0), Tuple{Int, Int}[])
        end
    end
    return found_vectors
end

function _MMT_minimum_distance_nonbinary(G::CTMatrixTypes, target_w::Int; p::Int = 4, l1::Int = 3, l2::Int = 3, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    l = l1 + l2
    
    F = base_ring(G); p_char, d_deg = Int(characteristic(F)), degree(F)
    non_zeros = filter(!iszero, collect(F))
    
    p1 = p ÷ 4; p2 = (p - p1) ÷ 3; p3 = (p - p1 - p2) ÷ 2; p4 = p - p1 - p2 - p3
    
    num_thrds = Threads.nthreads(); thread_load = max_iters ÷ num_thrds; remaining = max_iters % num_thrds
    keep_going = Threads.Atomic{Bool}(true); results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        q1 = k ÷ 4; q2 = k ÷ 2; q3 = 3 * k ÷ 4
        R1 = 1:q1; R2 = (q1+1):q2; R3 = (q2+1):q3; R4 = (q3+1):k
        tail_len = n - k - l
        
        for _ in 1:iters
            if !keep_going[] break end
            
            perm = collect(1:n); shuffle!(perm); G_loc = G[:, perm]
            try _make_systematic_gf!(G_loc, perm, k) catch; continue end
            
            w1_rows = [G_loc[i, (k+1):(k+l1)] for i in 1:k]
            w2_rows = [G_loc[i, (k+l1+1):(k+l)] for i in 1:k]
            tail_rows = [G_loc[i, (k+l+1):n] for i in 1:k]
            
            # --- 4-WAY MERGE TREE ---
            H1 = Dict{Vector{typeof(zero(F))}, Vector{Tuple{Vector{typeof(zero(F))}, Vector{Tuple{Int, typeof(zero(F))}}}}}()
            H12 = Dict{Vector{typeof(zero(F))}, Vector{Vector{Tuple{Int, typeof(zero(F))}}}}()
            H3 = Dict{Vector{typeof(zero(F))}, Vector{Tuple{Vector{typeof(zero(F))}, Vector{Tuple{Int, typeof(zero(F))}}}}}()
            
            function _build_base!(depth, picked, Range, p_tgt, cur_w1, cur_w2, msg, tgt_Dict)
                if picked == p_tgt
                    push!(get!(tgt_Dict, cur_w1, []), (cur_w2, copy(msg)))
                    return
                end
                if depth > length(Range) || (length(Range) - depth + 1) < (p_tgt - picked) return end
                
                _build_base!(depth+1, picked, Range, p_tgt, cur_w1, cur_w2, msg, tgt_Dict)
                
                idx = Range[depth]
                for scalar in non_zeros
                    n1 = cur_w1 .+ scalar .* w1_rows[idx]; n2 = cur_w2 .+ scalar .* w2_rows[idx]
                    push!(msg, (idx, scalar))
                    _build_base!(depth+1, picked+1, Range, p_tgt, n1, n2, msg, tgt_Dict)
                    pop!(msg)
                end
            end
            
            _build_base!(1, 0, R1, p1, zeros(F, l1), zeros(F, l2), Tuple{Int, typeof(zero(F))}[], H1)
            
            function _build_H12!(depth, picked, cur_w1, cur_w2, msg2)
                if picked == p2
                    # GF(q) Negation Lookup: lookup -cur_w1
                    neg_w1 = -cur_w1
                    if haskey(H1, neg_w1)
                        for (w2_1, msg1) in H1[neg_w1]
                            push!(get!(H12, cur_w2 .+ w2_1, []), vcat(msg1, msg2))
                        end
                    end
                    return
                end
                if depth > length(R2) || (length(R2) - depth + 1) < (p2 - picked) return end
                
                _build_H12!(depth+1, picked, cur_w1, cur_w2, msg2)
                idx = R2[depth]
                for scalar in non_zeros
                    push!(msg2, (idx, scalar))
                    _build_H12!(depth+1, picked+1, cur_w1 .+ scalar .* w1_rows[idx], cur_w2 .+ scalar .* w2_rows[idx], msg2)
                    pop!(msg2)
                end
            end
            _build_H12!(1, 0, zeros(F, l1), zeros(F, l2), Tuple{Int, typeof(zero(F))}[])
            empty!(H1)
            
            _build_base!(1, 0, R3, p3, zeros(F, l1), zeros(F, l2), Tuple{Int, typeof(zero(F))}[], H3)
            
            function _probe_H34!(depth, picked, cur_w1, cur_w2, msg4)
                if !keep_going[] return end
                if picked == p4
                    neg_w1 = -cur_w1
                    if haskey(H3, neg_w1)
                        for (w2_3, msg3) in H3[neg_w1]
                            w2_34 = cur_w2 .+ w2_3
                            neg_w2_34 = -w2_34
                            
                            if haskey(H12, neg_w2_34) 
                                for msg12 in H12[neg_w2_34]
                                    tail_buf = zeros(F, tail_len)
                                    for sub_msg in (msg12, msg3, msg4)
                                        for (idx, sc) in sub_msg
                                            tail_buf .+= sc .* tail_rows[idx]
                                        end
                                    end
                                    
                                    wt_tail = count(!iszero, tail_buf)
                                    if 0 < wt_tail + p <= target_w
                                        lock(results_lock) do
                                            if length(found_vectors) < num_find
                                                e_loc = zeros(F, n)
                                                for sub_msg in (msg12, msg3, msg4)
                                                    for (idx, sc) in sub_msg e_loc[idx] = sc end
                                                end
                                                e_loc[(k+l+1):n] .= tail_buf
                                                
                                                inv_p = invperm(perm)
                                                push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                                                if length(found_vectors) >= num_find keep_going[] = false end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    return
                end
                if depth > length(R4) || (length(R4) - depth + 1) < (p4 - picked) return end
                
                _probe_H34!(depth+1, picked, cur_w1, cur_w2, msg4)
                idx = R4[depth]
                for scalar in non_zeros
                    push!(msg4, (idx, scalar))
                    _probe_H34!(depth+1, picked+1, cur_w1 .+ scalar .* w1_rows[idx], cur_w2 .+ scalar .* w2_rows[idx], msg4)
                    pop!(msg4)
                end
            end
            _probe_H34!(1, 0, zeros(F, l1), zeros(F, l2), Tuple{Int, typeof(zero(F))}[])
        end
    end
    return found_vectors
end

"""
    _DOOM_Stern_minimum_distance_binary(G::Matrix{Int}, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)

Dedicated Minimum Distance solver using Sendrier's DOOM framework applied to Stern's algorithm.
Simultaneously targets all n columns of H to find an error vector of weight w-1.
"""
function _DOOM_Stern_minimum_distance_binary(G::Matrix{Int}, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)
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
        tail_len = n - k - l
        num_tail_chunks = cld(tail_len, 64)
        
        tail_buffer = zeros(UInt64, num_tail_chunks)
        
        for _ in 1:iters_for_this_thread
            if !keep_going[] break end
            
            σ_loc = collect(1:n)
            shuffle!(σ_loc)
            G_loc = G[:, σ_loc]
            
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # 1. PACK THE GENERATOR MATRIX
            win_rows = zeros(UInt64, k)
            tail_rows = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            
            for i in 1:k
                for j in 1:l
                    if G_loc[i, k + j] == 1 win_rows[i] |= (UInt64(1) << (j - 1)) end
                end
                for j in 1:tail_len
                    if G_loc[i, k + l + j] == 1
                        tail_rows[i][(j - 1) ÷ 64 + 1] |= (UInt64(1) << ((j - 1) % 64))
                    end
                end
            end
            
            # 2. DOOM: PACK THE N TARGET SYNDROMES (Columns of H)
            # For G = [I | P], H = [P^T | I]. 
            # We group the n targets by their l-window for instant hash lookups.
            grouped_targets = Dict{UInt64, Vector{Tuple{Int, Vector{UInt64}}}}()
            
            # Targets 1 to k (Rows of P)
            for i in 1:k
                push!(get!(grouped_targets, win_rows[i], []), (i, tail_rows[i]))
            end
            
            # Targets k+1 to n (Identity Matrix I)
            for j in 1:(n-k)
                t_win = UInt64(0)
                t_tail = zeros(UInt64, num_tail_chunks)
                
                if j <= l
                    t_win |= (UInt64(1) << (j - 1))
                else
                    idx = j - l
                    t_tail[(idx - 1) ÷ 64 + 1] |= (UInt64(1) << ((idx - 1) % 64))
                end
                push!(get!(grouped_targets, t_win, []), (k + j, t_tail))
            end
            
            # 3. RECURSIVE HASH MAP BUILDER (X-Half)
            hash_X = Dict{UInt64, Vector{Vector{Int}}}()
            function _build_X!(depth, picked, current_val::UInt64, msg)
                if picked == p
                    push!(get!(hash_X, current_val, Vector{Vector{Int}}()), copy(msg))
                    return
                end
                if depth > length(X_cols) || (length(X_cols) - depth + 1) < (p - picked) return end
                
                msg[depth] = 0; _build_X!(depth+1, picked, current_val, msg)
                msg[depth] = 1; _build_X!(depth+1, picked+1, current_val ⊻ win_rows[X_cols[depth]], msg)
            end
            _build_X!(1, 0, UInt64(0), zeros(Int, length(X_cols)))
            
            # 4. RECURSIVE DOOM PROBER (Y-Half)
            function _probe_Y_DOOM!(depth, picked, current_val::UInt64, msg_Y)
                if !keep_going[] return end
                
                if picked == p
                    # DOOM CHECK: Instead of cur_X ⊻ cur_Y = 0, we check all unique target windows!
                    for (t_win, target_list) in grouped_targets
                        search_val = current_val ⊻ t_win
                        
                        if haskey(hash_X, search_val)
                            for msg_X in hash_X[search_val]
                                for (t_idx, t_tail) in target_list
                                    
                                    fill!(tail_buffer, UInt64(0))
                                    for i in 1:length(msg_X) if msg_X[i] == 1 tail_buffer .⊻= tail_rows[X_cols[i]] end end
                                    for i in 1:length(msg_Y) if msg_Y[i] == 1 tail_buffer .⊻= tail_rows[Y_cols[i]] end end
                                    tail_buffer .⊻= t_tail # XOR the target tail!
                                    
                                    wt_tail = sum(count_ones.(tail_buffer))
                                    
                                    # Target weight is now (target_w - 1) because we are finding e'
                                    if wt_tail + 2*p == target_w - 1
                                        lock(results_lock) do
                                            if length(found_vectors) < num_find
                                                e_loc = zeros(Int, n)
                                                for i in 1:length(msg_X) if msg_X[i] == 1 e_loc[X_cols[i]] = 1 end end
                                                for i in 1:length(msg_Y) if msg_Y[i] == 1 e_loc[Y_cols[i]] = 1 end end
                                                
                                                for j in 1:tail_len
                                                    chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                                    if (tail_buffer[chunk] & (UInt64(1) << bit)) != 0 e_loc[k + l + j] = 1 end
                                                end
                                                
                                                # Add the target bit back in to complete the codeword!
                                                e_loc[t_idx] = 1 
                                                
                                                push!(found_vectors, e_loc[invperm(σ_loc)])
                                                if length(found_vectors) >= num_find keep_going[] = false end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    return
                end
                
                if depth > length(Y_cols) || (length(Y_cols) - depth + 1) < (p - picked) return end
                
                msg_Y[depth] = 0; _probe_Y_DOOM!(depth+1, picked, current_val, msg_Y)
                msg_Y[depth] = 1; _probe_Y_DOOM!(depth+1, picked+1, current_val ⊻ win_rows[Y_cols[depth]], msg_Y)
            end
            
            _probe_Y_DOOM!(1, 0, UInt64(0), zeros(Int, length(Y_cols)))
        end
    end
    
    return found_vectors
end

function _doom_stern_succ_prob(n::Int, k::Int, w::Int, p::Int, l::Int)
    k1 = k ÷ 2; k2 = k - k1
    if p > k1 || p > k2 return 0.0 end
    
    # Probability that the codeword has exactly 2p ones in the info set
    prob_2p = 0.0
    if (w - 2p) >= 1 && (w - 2p) <= (n - k - l)
        num_log1 = logbinomial(k1, p) + logbinomial(k2, p) + logbinomial(n - k - l, w - 2p)
        prob_2p = exp(num_log1 - logbinomial(n, w))
    end
    
    # Probability that the codeword has exactly 2p+1 ones in the info set
    prob_2p_plus_1 = 0.0
    if (w - 2p - 1) >= 0 && (w - 2p - 1) <= (n - k - l)
        # It could be p+1 in left, p in right OR p in left, p+1 in right
        num_log_L = logbinomial(k1, p+1) + logbinomial(k2, p) + logbinomial(n - k - l, w - 2p - 1)
        num_log_R = logbinomial(k1, p) + logbinomial(k2, p+1) + logbinomial(n - k - l, w - 2p - 1)
        prob_2p_plus_1 = exp(num_log_L - logbinomial(n, w)) + exp(num_log_R - logbinomial(n, w))
    end
    
    return prob_2p + prob_2p_plus_1
end

"""
    _DOOM_Stern_minimum_distance_GF4(G::CTMatrixTypes, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)
"""
function _DOOM_Stern_minimum_distance_GF4(G::CTMatrixTypes, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    F = base_ring(G); ω = gen(F)
    p_char, d_deg = Int(characteristic(F)), degree(F)
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds; remaining = max_iters % num_thrds
    keep_going = Threads.Atomic{Bool}(true); results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        
        half_k = k ÷ 2; X_cols = 1:half_k; Y_cols = (half_k + 1):k
        tail_len = n - k - l; num_tail_chunks = cld(tail_len, 64)
        tail_buf_a = zeros(UInt64, num_tail_chunks); tail_buf_b = zeros(UInt64, num_tail_chunks)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            σ_loc = collect(1:n); shuffle!(σ_loc); G_loc = G[:, σ_loc]
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # 1. PACK THE GENERATOR MATRIX
            win_a = zeros(UInt64, k); win_b = zeros(UInt64, k)
            tail_a = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            tail_b = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            
            for i in 1:k
                for j in 1:l
                    val = G_loc[i, k + j]
                    if val == 1 win_b[i] |= (UInt64(1) << (j - 1))
                    elseif val == ω win_a[i] |= (UInt64(1) << (j - 1))
                    elseif val == ω + 1 win_a[i] |= (UInt64(1) << (j - 1)); win_b[i] |= (UInt64(1) << (j - 1)) end
                end
                for j in 1:tail_len
                    val = G_loc[i, k + l + j]; chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                    if val == 1 tail_b[i][chunk] |= (UInt64(1) << bit)
                    elseif val == ω tail_a[i][chunk] |= (UInt64(1) << bit)
                    elseif val == ω + 1 tail_a[i][chunk] |= (UInt64(1) << bit); tail_b[i][chunk] |= (UInt64(1) << bit) end
                end
            end
            
            # 2. DOOM: PACK ALL SYNDROME TARGETS (α * H_i)
            grouped_targets = Dict{Tuple{UInt64, UInt64}, Vector{Tuple{Int, Int, Vector{UInt64}, Vector{UInt64}}}}()
            
            # Targets 1 to k (Rows of P) multiplied by scalars {1, 2, 3}
            for i in 1:k
                wa, wb, ta, tb = win_a[i], win_b[i], tail_a[i], tail_b[i]
                push!(get!(grouped_targets, (wa, wb), []), (i, 1, ta, tb)) # sc = 1
                push!(get!(grouped_targets, (wa ⊻ wb, wa), []), (i, 2, ta .⊻ tb, ta)) # sc = 2 (ω)
                push!(get!(grouped_targets, (wb, wa ⊻ wb), []), (i, 3, tb, ta .⊻ tb)) # sc = 3 (ω^2)
            end
            
            # Targets k+1 to n (Identity Matrix I) multiplied by scalars
            for j in 1:(n-k)
                for sc in 1:3
                    t_wa = UInt64(0); t_wb = UInt64(0)
                    t_ta = zeros(UInt64, num_tail_chunks); t_tb = zeros(UInt64, num_tail_chunks)
                    
                    if j <= l
                        if sc == 1 t_wb |= (UInt64(1) << (j - 1))
                        elseif sc == 2 t_wa |= (UInt64(1) << (j - 1))
                        elseif sc == 3 t_wa |= (UInt64(1) << (j - 1)); t_wb |= (UInt64(1) << (j - 1)) end
                    else
                        idx = j - l; chunk, bit = (idx - 1) ÷ 64 + 1, (idx - 1) % 64
                        if sc == 1 t_tb[chunk] |= (UInt64(1) << bit)
                        elseif sc == 2 t_ta[chunk] |= (UInt64(1) << bit)
                        elseif sc == 3 t_ta[chunk] |= (UInt64(1) << bit); t_tb[chunk] |= (UInt64(1) << bit) end
                    end
                    push!(get!(grouped_targets, (t_wa, t_wb), []), (k + j, sc, t_ta, t_tb))
                end
            end
            
            # 3. RECURSIVE HASH MAP BUILDER (X-Half)
            hash_X = Dict{Tuple{UInt64, UInt64}, Vector{Vector{Int}}}()
            function _build_X_GF4!(depth, picked, cur_a, cur_b, msg)
                if picked == p
                    push!(get!(hash_X, (cur_a, cur_b), []), copy(msg)); return
                end
                if depth > length(X_cols) || (length(X_cols) - depth + 1) < (p - picked) return end
                idx = X_cols[depth]
                msg[depth] = 0; _build_X_GF4!(depth+1, picked, cur_a, cur_b, msg)
                msg[depth] = 1; _build_X_GF4!(depth+1, picked+1, cur_a ⊻ win_a[idx], cur_b ⊻ win_b[idx], msg)
                msg[depth] = 2; _build_X_GF4!(depth+1, picked+1, cur_a ⊻ win_a[idx] ⊻ win_b[idx], cur_b ⊻ win_a[idx], msg)
                msg[depth] = 3; _build_X_GF4!(depth+1, picked+1, cur_a ⊻ win_b[idx], cur_b ⊻ win_a[idx] ⊻ win_b[idx], msg)
            end
            _build_X_GF4!(1, 0, UInt64(0), UInt64(0), zeros(Int, length(X_cols)))
            
            # 4. RECURSIVE DOOM PROBER (Y-Half)
            function _probe_Y_DOOM!(depth, picked, cur_a, cur_b, msg_Y)
                if !keep_going[] return end
                if picked == p
                    for ((t_wa, t_wb), target_list) in grouped_targets
                        search_a, search_b = cur_a ⊻ t_wa, cur_b ⊻ t_wb
                        if haskey(hash_X, (search_a, search_b))
                            for msg_X in hash_X[(search_a, search_b)]
                                for (t_idx, t_sc, t_ta, t_tb) in target_list
                                    fill!(tail_buf_a, UInt64(0)); fill!(tail_buf_b, UInt64(0))
                                    for (msg, cols) in ((msg_X, X_cols), (msg_Y, Y_cols))
                                        for (i, v) in enumerate(msg)
                                            if v != 0
                                                idx = cols[i]
                                                if v == 1 tail_buf_a .⊻= tail_a[idx]; tail_buf_b .⊻= tail_b[idx]
                                                elseif v == 2 tail_buf_a .⊻= tail_a[idx] .⊻ tail_b[idx]; tail_buf_b .⊻= tail_a[idx]
                                                elseif v == 3 tail_buf_a .⊻= tail_b[idx]; tail_buf_b .⊻= tail_a[idx] .⊻ tail_b[idx] end
                                            end
                                        end
                                    end
                                    tail_buf_a .⊻= t_ta; tail_buf_b .⊻= t_tb # Target tail XOR
                                    
                                    if sum(count_ones.(tail_buf_a .| tail_buf_b)) + 2*p == target_w - 1
                                        lock(results_lock) do
                                            if length(found_vectors) < num_find
                                                e_loc = zeros(F, n)
                                                for (msg, cols) in ((msg_X, X_cols), (msg_Y, Y_cols))
                                                    for (i, v) in enumerate(msg) if v != 0 e_loc[cols[i]] = v == 1 ? F(1) : (v == 2 ? ω : ω + 1) end end
                                                end
                                                for j in 1:tail_len
                                                    chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                                    bit_a = (tail_buf_a[chunk] >> bit) & 1; bit_b = (tail_buf_b[chunk] >> bit) & 1
                                                    if bit_a == 1 && bit_b == 0 e_loc[k + l + j] = ω
                                                    elseif bit_a == 1 && bit_b == 1 e_loc[k + l + j] = ω + 1
                                                    elseif bit_a == 0 && bit_b == 1 e_loc[k + l + j] = F(1) end
                                                end
                                                e_loc[t_idx] = t_sc == 1 ? F(1) : (t_sc == 2 ? ω : ω + 1)
                                                
                                                inv_p = invperm(σ_loc)
                                                push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                                                if length(found_vectors) >= num_find keep_going[] = false end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    return
                end
                if depth > length(Y_cols) || (length(Y_cols) - depth + 1) < (p - picked) return end
                idx = Y_cols[depth]
                msg_Y[depth] = 0; _probe_Y_DOOM!(depth+1, picked, cur_a, cur_b, msg_Y)
                msg_Y[depth] = 1; _probe_Y_DOOM!(depth+1, picked+1, cur_a ⊻ win_a[idx], cur_b ⊻ win_b[idx], msg_Y)
                msg_Y[depth] = 2; _probe_Y_DOOM!(depth+1, picked+1, cur_a ⊻ win_a[idx] ⊻ win_b[idx], cur_b ⊻ win_a[idx], msg_Y)
                msg_Y[depth] = 3; _probe_Y_DOOM!(depth+1, picked+1, cur_a ⊻ win_b[idx], cur_b ⊻ win_a[idx] ⊻ win_b[idx], msg_Y)
            end
            _probe_Y_DOOM!(1, 0, UInt64(0), UInt64(0), zeros(Int, length(Y_cols)))
        end
    end
    return found_vectors
end

"""
    _DOOM_Stern_minimum_distance_GF3(G::CTMatrixTypes, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)
"""
function _DOOM_Stern_minimum_distance_GF3(G::CTMatrixTypes, target_w::Int; p::Int = 2, l::Int = 12, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    F = base_ring(G); p_char, d_deg = Int(characteristic(F)), degree(F)
    
    num_thrds = Threads.nthreads(); thread_load = max_iters ÷ num_thrds; remaining = max_iters % num_thrds
    keep_going = Threads.Atomic{Bool}(true); results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    @inline function add_mod3(AH::UInt64, AL::UInt64, BH::UInt64, BL::UInt64)
        SH = ((AH ⊻ BH) & ~(AL | BL)) | (AL & BL)
        SL = ((AL ⊻ BL) & ~(AH | BH)) | (AH & BH)
        return SH, SL
    end
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        half_k = k ÷ 2; X_cols = 1:half_k; Y_cols = (half_k + 1):k
        tail_len = n - k - l; num_tail_chunks = cld(tail_len, 64)
        tail_buf_H = zeros(UInt64, num_tail_chunks); tail_buf_L = zeros(UInt64, num_tail_chunks)
        
        for _ in 1:iters
            if !keep_going[] break end
            σ_loc = collect(1:n); shuffle!(σ_loc); G_loc = G[:, σ_loc]
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            win_H = zeros(UInt64, k); win_L = zeros(UInt64, k)
            tail_H = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            tail_L = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            
            for i in 1:k
                for j in 1:l
                    if G_loc[i, k + j] == 1 win_L[i] |= (UInt64(1) << (j - 1))
                    elseif G_loc[i, k + j] == 2 win_H[i] |= (UInt64(1) << (j - 1)) end
                end
                for j in 1:tail_len
                    v = G_loc[i, k + l + j]; chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                    if v == 1 tail_L[i][chunk] |= (UInt64(1) << bit)
                    elseif v == 2 tail_H[i][chunk] |= (UInt64(1) << bit) end
                end
            end
            
            # DOOM TARGETS (sc=1 and sc=2 -> swapping H and L)
            grouped_targets = Dict{Tuple{UInt64, UInt64}, Vector{Tuple{Int, Int, Vector{UInt64}, Vector{UInt64}}}}()
            for i in 1:k
                push!(get!(grouped_targets, (win_H[i], win_L[i]), []), (i, 1, tail_H[i], tail_L[i]))
                push!(get!(grouped_targets, (win_L[i], win_H[i]), []), (i, 2, tail_L[i], tail_H[i]))
            end
            
            for j in 1:(n-k)
                for sc in 1:2
                    t_wH = UInt64(0); t_wL = UInt64(0)
                    t_tH = zeros(UInt64, num_tail_chunks); t_tL = zeros(UInt64, num_tail_chunks)
                    if j <= l
                        if sc == 1 t_wL |= (UInt64(1) << (j - 1)) else t_wH |= (UInt64(1) << (j - 1)) end
                    else
                        idx = j - l; chunk, bit = (idx - 1) ÷ 64 + 1, (idx - 1) % 64
                        if sc == 1 t_tL[chunk] |= (UInt64(1) << bit) else t_tH[chunk] |= (UInt64(1) << bit) end
                    end
                    push!(get!(grouped_targets, (t_wH, t_wL), []), (k + j, sc, t_tH, t_tL))
                end
            end
            
            hash_X = Dict{Tuple{UInt64, UInt64}, Vector{Vector{Int}}}()
            function _build_X_GF3!(depth, picked, cH, cL, msg)
                if picked == p push!(get!(hash_X, (cH, cL), []), copy(msg)); return end
                if depth > length(X_cols) || (length(X_cols) - depth + 1) < (p - picked) return end
                idx = X_cols[depth]
                msg[depth] = 0; _build_X_GF3!(depth+1, picked, cH, cL, msg)
                nH, nL = add_mod3(cH, cL, win_H[idx], win_L[idx]); msg[depth] = 1; _build_X_GF3!(depth+1, picked+1, nH, nL, msg)
                nH, nL = add_mod3(cH, cL, win_L[idx], win_H[idx]); msg[depth] = 2; _build_X_GF3!(depth+1, picked+1, nH, nL, msg)
            end
            _build_X_GF3!(1, 0, UInt64(0), UInt64(0), zeros(Int, length(X_cols)))
            
            function _probe_Y_DOOM!(depth, picked, cH, cL, msg_Y)
                if !keep_going[] return end
                if picked == p
                    for ((t_wH, t_wL), tgt_list) in grouped_targets
                        # X = T - Y = T + 2Y (which is T + (cL, cH))
                        search_H, search_L = add_mod3(t_wH, t_wL, cL, cH)
                        if haskey(hash_X, (search_H, search_L))
                            for msg_X in hash_X[(search_H, search_L)]
                                for (t_idx, t_sc, t_tH, t_tL) in tgt_list
                                    fill!(tail_buf_H, UInt64(0)); fill!(tail_buf_L, UInt64(0))
                                    for (msg, cols) in ((msg_X, X_cols), (msg_Y, Y_cols))
                                        for (i, v) in enumerate(msg)
                                            if v != 0
                                                idx = cols[i]
                                                RH = v == 1 ? tail_H[idx] : tail_L[idx]
                                                RL = v == 1 ? tail_L[idx] : tail_H[idx]
                                                for c in 1:num_tail_chunks tail_buf_H[c], tail_buf_L[c] = add_mod3(tail_buf_H[c], tail_buf_L[c], RH[c], RL[c]) end
                                            end
                                        end
                                    end
                                    # Subtract target tail: tail - T = tail + 2T = tail + (t_tL, t_tH)
                                    for c in 1:num_tail_chunks tail_buf_H[c], tail_buf_L[c] = add_mod3(tail_buf_H[c], tail_buf_L[c], t_tL[c], t_tH[c]) end
                                    
                                    if sum(count_ones.(tail_buf_H .| tail_buf_L)) + 2*p == target_w - 1
                                        lock(results_lock) do
                                            if length(found_vectors) < num_find
                                                e_loc = zeros(F, n)
                                                for (msg, cols) in ((msg_X, X_cols), (msg_Y, Y_cols))
                                                    for (i, v) in enumerate(msg) if v != 0 e_loc[cols[i]] = F(v) end end
                                                end
                                                for j in 1:tail_len
                                                    chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                                    if ((tail_buf_H[chunk] >> bit) & 1) == 1 e_loc[k + l + j] = F(2)
                                                    elseif ((tail_buf_L[chunk] >> bit) & 1) == 1 e_loc[k + l + j] = F(1) end
                                                end
                                                e_loc[t_idx] = F(t_sc)
                                                inv_p = invperm(σ_loc)
                                                push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                                                if length(found_vectors) >= num_find keep_going[] = false end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    return
                end
                if depth > length(Y_cols) || (length(Y_cols) - depth + 1) < (p - picked) return end
                idx = Y_cols[depth]
                msg_Y[depth] = 0; _probe_Y_DOOM!(depth+1, picked, cH, cL, msg_Y)
                nH, nL = add_mod3(cH, cL, win_H[idx], win_L[idx]); msg_Y[depth] = 1; _probe_Y_DOOM!(depth+1, picked+1, nH, nL, msg_Y)
                nH, nL = add_mod3(cH, cL, win_L[idx], win_H[idx]); msg_Y[depth] = 2; _probe_Y_DOOM!(depth+1, picked+1, nH, nL, msg_Y)
            end
            _probe_Y_DOOM!(1, 0, UInt64(0), UInt64(0), zeros(Int, length(Y_cols)))
        end
    end
    return found_vectors
end

function _DOOM_Stern_minimum_distance_nonbinary(G::CTMatrixTypes, target_w::Int; p::Int = 2, l::Int = 3, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    F = base_ring(G); p_char, d_deg = Int(characteristic(F)), degree(F)
    non_zeros = filter(!iszero, collect(F))
    
    num_thrds = Threads.nthreads(); thread_load = max_iters ÷ num_thrds; remaining = max_iters % num_thrds
    keep_going = Threads.Atomic{Bool}(true); results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        half_k = k ÷ 2; X_cols = 1:half_k; Y_cols = (half_k + 1):k
        tail_len = n - k - l
        
        for _ in 1:iters
            if !keep_going[] break end
            perm = collect(1:n); shuffle!(perm); G_loc = G[:, perm]
            try _make_systematic_gf!(G_loc, perm, k) catch; continue end
            
            win_rows = [G_loc[i, (k+1):(k+l)] for i in 1:k]
            tail_rows = [G_loc[i, (k+l+1):n] for i in 1:k]
            
            grouped_targets = Dict{Vector{typeof(zero(F))}, Vector{Tuple{Int, typeof(zero(F)), Vector{typeof(zero(F))}}}}()
            for i in 1:k
                for sc in non_zeros
                    push!(get!(grouped_targets, sc .* win_rows[i], []), (i, sc, sc .* tail_rows[i]))
                end
            end
            for j in 1:(n-k)
                for sc in non_zeros
                    t_win = zeros(F, l); t_tail = zeros(F, tail_len)
                    if j <= l t_win[j] = sc else t_tail[j - l] = sc end
                    push!(get!(grouped_targets, t_win, []), (k + j, sc, t_tail))
                end
            end
            
            hash_X = Dict{Vector{typeof(zero(F))}, Vector{Vector{Tuple{Int, typeof(zero(F))}}}}()
            function _build_X!(depth, picked, cur_w, msg)
                if picked == p push!(get!(hash_X, cur_w, []), copy(msg)); return end
                if depth > length(X_cols) || (length(X_cols) - depth + 1) < (p - picked) return end
                _build_X!(depth+1, picked, cur_w, msg)
                idx = X_cols[depth]
                for sc in non_zeros
                    push!(msg, (idx, sc)); _build_X!(depth+1, picked+1, cur_w .+ sc .* win_rows[idx], msg); pop!(msg)
                end
            end
            _build_X!(1, 0, zeros(F, l), Tuple{Int, typeof(zero(F))}[])
            
            function _probe_Y_DOOM!(depth, picked, cur_w, msg_Y)
                if !keep_going[] return end
                if picked == p
                    for (t_win, tgt_list) in grouped_targets
                        search_val = -(cur_w .+ t_win)
                        if haskey(hash_X, search_val)
                            for msg_X in hash_X[search_val]
                                for (t_idx, t_sc, t_tail) in tgt_list
                                    t_buf = zeros(F, tail_len)
                                    for msg in (msg_X, msg_Y)
                                        for (idx, sc) in msg t_buf .+= sc .* tail_rows[idx] end
                                    end
                                    t_buf .+= t_tail
                                    
                                    if count(!iszero, t_buf) + 2*p == target_w - 1
                                        lock(results_lock) do
                                            if length(found_vectors) < num_find
                                                e_loc = zeros(F, n)
                                                for msg in (msg_X, msg_Y)
                                                    for (idx, sc) in msg e_loc[idx] = sc end
                                                end
                                                e_loc[(k+l+1):n] .= t_buf; e_loc[t_idx] = t_sc
                                                inv_p = invperm(perm)
                                                push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                                                if length(found_vectors) >= num_find keep_going[] = false end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    return
                end
                if depth > length(Y_cols) || (length(Y_cols) - depth + 1) < (p - picked) return end
                _probe_Y_DOOM!(depth+1, picked, cur_w, msg_Y)
                idx = Y_cols[depth]
                for sc in non_zeros
                    push!(msg_Y, (idx, sc)); _probe_Y_DOOM!(depth+1, picked+1, cur_w .+ sc .* win_rows[idx], msg_Y); pop!(msg_Y)
                end
            end
            _probe_Y_DOOM!(1, 0, zeros(F, l), Tuple{Int, typeof(zero(F))}[])
        end
    end
    return found_vectors
end

"""
    probabilistic_minimum_distance_stern_DOOM(C::AbstractLinearCode; confidence::Float64 = 0.99, p::Int = 2, l::Int = (Int(order(C.F)) == 2 ? 12 : 3), verbose::Bool=true)

Finds the minimum distance using Sendrier's DOOM (Decoding One Out of Many) framework applied to Stern's algorithm. 
Automatically routes to bitsliced hardware-accelerated engines for GF(2), GF(3), and GF(4).

DOOM targets an error vector of weight w-1 by simultaneously checking all n columns of the parity-check matrix, 
drastically reducing the required number of iterations compared to standard Stern.
"""
function probabilistic_minimum_distance_stern_DOOM(C::AbstractLinearCode; confidence::Float64 = 0.99, p::Int = 2, l::Int = (Int(order(C.F)) == 2 ? 12 : 3), verbose::Bool=true)
    k, n = C.k, C.n
    q = Int(order(C.F))
    ϵ = 1.0 - confidence
    
    # Extract the correct matrix type
    G = q == 2 ? _convert_binary_to_int_matrix(generator_matrix(C, true)) : generator_matrix(C, true)
    
    for w in 1:n
        # Use the DOOM-specific probability calculator (which checks for 2p or 2p+1 errors in the info set)
        p_succ = _doom_stern_succ_prob(n, k, w, p, l)
        
        if p_succ <= 0.0 continue end
        
        req_iters = p_succ >= 1.0 ? 1 : Int(ceil(log(ϵ) / log(1.0 - p_succ)))
        verbose && println("DOOM-Stern (w=$w, p=$p, l=$l): Requires $req_iters iterations for $(confidence*100)% confidence...")
        
        # ADVANCED ROUTING LOGIC
        if q == 2
            found = _DOOM_Stern_minimum_distance_binary(G, w; p=p, l=l, num_find=1, max_iters=req_iters)
        elseif q == 3
            found = _DOOM_Stern_minimum_distance_GF3(G, w; p=p, l=l, num_find=1, max_iters=req_iters)
        elseif q == 4
            found = _DOOM_Stern_minimum_distance_GF4(G, w; p=p, l=l, num_find=1, max_iters=req_iters)
        else
            found = _DOOM_Stern_minimum_distance_nonbinary(G, w; p=p, l=l, num_find=1, max_iters=req_iters)
        end
        
        if !isempty(found)
            # Count non-zeros to accurately reflect Hamming weight
            actual_w = sum(count.(!iszero, first(found))) 
            verbose && println("🎉 Success! Found codeword of weight $actual_w.")
            return actual_w, first(found)
        end
    end
    
    return n, zeros(Int, n)
end

"""
    _BJMM_minimum_distance_binary(G::Matrix{Int}, target_w::Int; p::Int = 4, ϵ1::Int = 1, l1::Int = 10, l2::Int = 10, num_find::Int = 1, max_iters::Int = 10000)

Dedicated Minimum Distance solver using the Becker-Joux-May-Meurer (BJMM) algorithm.
Utilizes heavier base weights (p/4 + ϵ1) and exploits intentional 1-cancellation via bitwise XOR 
to explosively increase the representation count and filter size.
"""
function _BJMM_minimum_distance_binary(G::Matrix{Int}, target_w::Int; p::Int = 4, ϵ1::Int = 1, l1::Int = 10, l2::Int = 10, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    l = l1 + l2
    @assert l <= 64 "Combined window size (l1 + l2) must be <= 64."
    @assert p >= 4 "Target info weight p must be >= 4."
    
    # Base weight per quarter includes the ϵ1 overlap inflation
    base_wt = (p ÷ 4) + ϵ1
    # Target weight after Level 1 merge (p/2)
    lvl1_target_wt = p ÷ 2
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds
    remaining = max_iters % num_thrds
    
    keep_going = Threads.Atomic{Bool}(true)
    results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        
        # Quarter bounds
        q1 = k ÷ 4; q2 = k ÷ 2; q3 = 3 * k ÷ 4
        R1 = 1:q1; R2 = (q1+1):q2; R3 = (q2+1):q3; R4 = (q3+1):k
        
        tail_len = n - k - l
        num_tail_chunks = cld(tail_len, 64)
        tail_buf = zeros(UInt64, num_tail_chunks)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            σ_loc = collect(1:n)
            shuffle!(σ_loc)
            G_loc = G[:, σ_loc]
            
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # --- PACK WINDOWS AND TAIL ---
            win1_rows = zeros(UInt64, k)
            win2_rows = zeros(UInt64, k)
            info_rows = [zeros(UInt64, cld(k, 64)) for _ in 1:k] # For tracking column masks
            tail_rows = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            
            for i in 1:k
                info_rows[i][(i-1)÷64 + 1] |= (UInt64(1) << ((i-1)%64))
                for j in 1:l1
                    if G_loc[i, k + j] == 1 win1_rows[i] |= (UInt64(1) << (j - 1)) end
                end
                for j in 1:l2
                    if G_loc[i, k + l1 + j] == 1 win2_rows[i] |= (UInt64(1) << (j - 1)) end
                end
                for j in 1:tail_len
                    if G_loc[i, k + l + j] == 1
                        tail_rows[i][(j - 1) ÷ 64 + 1] |= (UInt64(1) << ((j - 1) % 64))
                    end
                end
            end
            
            # --- THE BJMM MERGE TREE DICTIONARIES ---
            # Storing (win2, info_mask, cols_used)
            H1 = Dict{UInt64, Vector{Tuple{UInt64, Vector{UInt64}, Vector{Int}}}}()
            H12 = Dict{UInt64, Vector{Tuple{Vector{UInt64}, Vector{Int}}}}()
            H3 = Dict{UInt64, Vector{Tuple{UInt64, Vector{UInt64}, Vector{Int}}}}()
            
            # 1. Build H1 (Base Level 1 with base_wt)
            function _build_H1!(depth, picked, cur_w1, cur_w2, cur_info, msg)
                if picked == base_wt
                    push!(get!(H1, cur_w1, Tuple{UInt64, Vector{UInt64}, Vector{Int}}[]), (cur_w2, copy(cur_info), copy(msg)))
                    return
                end
                if depth > length(R1) || (length(R1) - depth + 1) < (base_wt - picked) return end
                
                _build_H1!(depth+1, picked, cur_w1, cur_w2, cur_info, msg)
                idx = R1[depth]
                push!(msg, idx)
                next_info = copy(cur_info)
                next_info[(idx-1)÷64 + 1] |= (UInt64(1) << ((idx-1)%64))
                _build_H1!(depth+1, picked+1, cur_w1 ⊻ win1_rows[idx], cur_w2 ⊻ win2_rows[idx], next_info, msg)
                pop!(msg)
            end
            _build_H1!(1, 0, UInt64(0), UInt64(0), zeros(UInt64, cld(k, 64)), Int[])
            
            # 2. Probe H1 with R2 -> Build H12 (Level 1 Merge with OVERLAP CHECK)
            function _build_H12!(depth, picked, cur_w1, cur_w2, cur_info, msg2)
                if picked == base_wt
                    if haskey(H1, cur_w1) # Window 1 Collision!
                        for (w2_1, info_1, msg1) in H1[cur_w1]
                            
                            # CRITICAL BJMM OVERLAP CHECK: Calculate true weight of info-XOR sum
                            # Because of base_wt inflation, some 1s MUST cancel out here to hit lvl1_target_wt
                            info_wt = sum(count_ones.(info_1 .⊻ cur_info))
                            
                            if info_wt == lvl1_target_wt
                                push!(get!(H12, cur_w2 ⊻ w2_1, Tuple{Vector{UInt64}, Vector{Int}}[]), (info_1 .⊻ cur_info, vcat(msg1, msg2)))
                            end
                        end
                    end
                    return
                end
                if depth > length(R2) || (length(R2) - depth + 1) < (base_wt - picked) return end
                
                _build_H12!(depth+1, picked, cur_w1, cur_w2, cur_info, msg2)
                idx = R2[depth]
                push!(msg2, idx)
                next_info = copy(cur_info)
                next_info[(idx-1)÷64 + 1] |= (UInt64(1) << ((idx-1)%64))
                _build_H12!(depth+1, picked+1, cur_w1 ⊻ win1_rows[idx], cur_w2 ⊻ win2_rows[idx], next_info, msg2)
                pop!(msg2)
            end
            _build_H12!(1, 0, UInt64(0), UInt64(0), zeros(UInt64, cld(k, 64)), Int[])
            
            empty!(H1) # Clear L1 from memory immediately
            
            # 3. Build H3 (Base Level 3 with base_wt)
            function _build_H3!(depth, picked, cur_w1, cur_w2, cur_info, msg)
                if picked == base_wt
                    push!(get!(H3, cur_w1, Tuple{UInt64, Vector{UInt64}, Vector{Int}}[]), (cur_w2, copy(cur_info), copy(msg)))
                    return
                end
                if depth > length(R3) || (length(R3) - depth + 1) < (base_wt - picked) return end
                
                _build_H3!(depth+1, picked, cur_w1, cur_w2, cur_info, msg)
                idx = R3[depth]
                push!(msg, idx)
                next_info = copy(cur_info)
                next_info[(idx-1)÷64 + 1] |= (UInt64(1) << ((idx-1)%64))
                _build_H3!(depth+1, picked+1, cur_w1 ⊻ win1_rows[idx], cur_w2 ⊻ win2_rows[idx], next_info, msg)
                pop!(msg)
            end
            _build_H3!(1, 0, UInt64(0), UInt64(0), zeros(UInt64, cld(k, 64)), Int[])
            
            # 4. Probe H3 with R4 -> Probe H12 -> Level 2 Merge + Tail Check
            function _probe_H34!(depth, picked, cur_w1, cur_w2, cur_info, msg4)
                if !keep_going[] return end
                
                if picked == base_wt
                    if haskey(H3, cur_w1) # Window 1 Collision!
                        for (w2_3, info_3, msg3) in H3[cur_w1]
                            w2_34 = cur_w2 ⊻ w2_3
                            
                            if haskey(H12, w2_34) # Window 2 Collision!
                                info_34 = info_3 .⊻ cur_info
                                
                                for (info_12, msg12) in H12[w2_34]
                                    
                                    # CRITICAL BJMM FINAL WEIGHT CHECK: Do we hit exactly the target info weight p?
                                    if sum(count_ones.(info_12 .⊻ info_34)) == p
                                        
                                        # Map survived all filters! Evaluate tail.
                                        fill!(tail_buf, UInt64(0))
                                        for c in msg12 tail_buf .⊻= tail_rows[c] end
                                        for c in msg3  tail_buf .⊻= tail_rows[c] end
                                        for c in msg4  tail_buf .⊻= tail_rows[c] end
                                        
                                        wt_tail = sum(count_ones.(tail_buf))
                                        if 0 < wt_tail + p <= target_w
                                            lock(results_lock) do
                                                if length(found_vectors) < num_find
                                                    e_loc = zeros(Int, n)
                                                    
                                                    # Decode information columns via combined mask
                                                    final_info_mask = info_12 .⊻ info_34
                                                    for i in 1:k
                                                        if (final_info_mask[(i-1)÷64 + 1] & (UInt64(1) << ((i-1)%64))) != 0
                                                            e_loc[i] = 1
                                                        end
                                                    end
                                                    
                                                    # Decode parity columns
                                                    for j in 1:tail_len
                                                        chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                                        if (tail_buf[chunk] & (UInt64(1) << bit)) != 0
                                                            e_loc[k + l + j] = 1
                                                        end
                                                    end
                                                    
                                                    push!(found_vectors, e_loc[invperm(σ_loc)])
                                                    if length(found_vectors) >= num_find keep_going[] = false end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    return
                end
                if depth > length(R4) || (length(R4) - depth + 1) < (base_wt - picked) return end
                
                _probe_H34!(depth+1, picked, cur_w1, cur_w2, cur_info, msg4)
                idx = R4[depth]
                push!(msg4, idx)
                next_info = copy(cur_info)
                next_info[(idx-1)÷64 + 1] |= (UInt64(1) << ((idx-1)%64))
                _probe_H34!(depth+1, picked+1, cur_w1 ⊻ win1_rows[idx], cur_w2 ⊻ win2_rows[idx], next_info, msg4)
                pop!(msg4)
            end
            _probe_H34!(1, 0, UInt64(0), UInt64(0), zeros(UInt64, cld(k, 64)), Int[])
        end
    end
    return found_vectors
end

function _bjmm_succ_prob(n::Int, k::Int, w::Int, p::Int, ϵ1::Int, l1::Int, l2::Int)
    l = l1 + l2
    if p > w || p > k || (w - p) > (n - k - l) return 0.0 end
    
    base_wt = (p ÷ 4) + ϵ1
    q1 = k ÷ 4; q2 = k ÷ 2; q3 = 3 * k ÷ 4
    s1 = q1; s2 = q2 - q1; s3 = q3 - q2; s4 = k - q3
    
    if base_wt > s1 || base_wt > s2 || base_wt > s3 || base_wt > s4 return 0.0 end
    
    # Probability of capturing the initial support inflated layout
    num_log = logbinomial(s1, base_wt) + logbinomial(s2, base_wt) + logbinomial(s3, base_wt) + logbinomial(s4, base_wt) + logbinomial(n - k - l, w - p)
    den_log = logbinomial(n, w)
    
    return exp(num_log - den_log)
end

"""
    probabilistic_minimum_distance_bjmm(C::AbstractLinearCode; confidence::Float64 = 0.99, p::Int = 4, ϵ1::Int = 1, l1::Int = (Int(order(C.F)) == 2 ? 10 : 2), l2::Int = (Int(order(C.F)) == 2 ? 10 : 2), verbose::Bool=true)

Finds the minimum distance using the Becker-Joux-May-Meurer (BJMM) 3rd-Generation algorithm.
Automatically hooks up bitsliced SIMD hardware-acceleration engines for GF(2), GF(3), and GF(4).
"""
function probabilistic_minimum_distance_bjmm(C::AbstractLinearCode; confidence::Float64 = 0.99, p::Int = 4, ϵ1::Int = 1, l1::Int = (Int(order(C.F)) == 2 ? 10 : 2), l2::Int = (Int(order(C.F)) == 2 ? 10 : 2), verbose::Bool=true)
    k, n = C.k, C.n
    q = Int(order(C.F))
    ϵ = 1.0 - confidence
    
    G = q == 2 ? _convert_binary_to_int_matrix(generator_matrix(C, true)) : generator_matrix(C, true)
    
    for w in 1:n
        p_succ = _bjmm_succ_prob(n, k, w, p, ϵ1, l1, l2)
        if p_succ <= 0.0 continue end
        
        req_iters = p_succ >= 1.0 ? 1 : Int(ceil(log(ϵ) / log(1.0 - p_succ)))
        verbose && println("BJMM (w=$w, p=$p, ϵ1=$ϵ1, l1=$l1, l2=$l2): Requires $req_iters iterations for $(confidence*100)% confidence...")
        
        if q == 2
            found = _BJMM_minimum_distance_binary(G, w; p=p, ϵ1=ϵ1, l1=l1, l2=l2, num_find=1, max_iters=req_iters)
        elseif q == 3
            found = _BJMM_minimum_distance_GF3(G, w; p=p, ϵ1=ϵ1, l1=l1, l2=l2, num_find=1, max_iters=req_iters)
        elseif q == 4
            found = _BJMM_minimum_distance_GF4(G, w; p=p, ϵ1=ϵ1, l1=l1, l2=l2, num_find=1, max_iters=req_iters)
        else
            found = _BJMM_minimum_distance_nonbinary(G, w; p=p, ϵ1=ϵ1, l1=l1, l2=l2, num_find=1, max_iters=req_iters)
        end
        
        if !isempty(found)
            actual_w = sum(count.(!iszero, first(found)))
            verbose && println("🎉 Success! Found codeword of weight $actual_w.")
            return actual_w, first(found)
        end
    end
    return n, zeros(Int, n)
end

"""
    _BJMM_minimum_distance_GF4(G::CTMatrixTypes, target_w::Int; p::Int = 4, ϵ1::Int = 1, l1::Int = 10, l2::Int = 10, num_find::Int = 1, max_iters::Int = 10000)

Dedicated Bitsliced Minimum Distance solver for GF(4) codes using the BJMM algorithm.
Exploits native self-inverse cancellation (a = -a) using parallel UInt64 register layers.
"""
function _BJMM_minimum_distance_GF4(G::CTMatrixTypes, target_w::Int; p::Int = 4, ϵ1::Int = 1, l1::Int = 10, l2::Int = 10, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    l = l1 + l2
    @assert l <= 64 "Combined window size (l1 + l2) must be <= 64."
    
    F = base_ring(G); ω = gen(F)
    p_char, d_deg = Int(characteristic(F)), degree(F)
    
    base_wt = (p ÷ 4) + ϵ1
    lvl1_target_wt = p ÷ 2
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds; remaining = max_iters % num_thrds
    keep_going = Threads.Atomic{Bool}(true); results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        
        q1 = k ÷ 4; q2 = k ÷ 2; q3 = 3 * k ÷ 4
        R1 = 1:q1; R2 = (q1+1):q2; R3 = (q2+1):q3; R4 = (q3+1):k
        
        tail_len = n - k - l; num_tail_chunks = cld(tail_len, 64)
        tail_buf_a = zeros(UInt64, num_tail_chunks); tail_buf_b = zeros(UInt64, num_tail_chunks)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            σ_loc = collect(1:n); shuffle!(σ_loc); G_loc = G[:, σ_loc]
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # --- BITSLICING COMPILER ---
            w1_a = zeros(UInt64, k); w1_b = zeros(UInt64, k)
            w2_a = zeros(UInt64, k); w2_b = zeros(UInt64, k)
            tail_a = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            tail_b = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            
            for i in 1:k
                for j in 1:l1
                    val = G_loc[i, k + j]
                    if val == 1 w1_b[i] |= (UInt64(1) << (j - 1))
                    elseif val == ω w1_a[i] |= (UInt64(1) << (j - 1))
                    elseif val == ω + 1 w1_a[i] |= (UInt64(1) << (j - 1)); w1_b[i] |= (UInt64(1) << (j - 1)) end
                end
                for j in 1:l2
                    val = G_loc[i, k + l1 + j]
                    if val == 1 w2_b[i] |= (UInt64(1) << (j - 1))
                    elseif val == ω w2_a[i] |= (UInt64(1) << (j - 1))
                    elseif val == ω + 1 w2_a[i] |= (UInt64(1) << (j - 1)); w2_b[i] |= (UInt64(1) << (j - 1)) end
                end
                for j in 1:tail_len
                    val = G_loc[i, k + l + j]; chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                    if val == 1 tail_b[i][chunk] |= (UInt64(1) << bit)
                    elseif val == ω tail_a[i][chunk] |= (UInt64(1) << bit)
                    elseif val == ω + 1 tail_a[i][chunk] |= (UInt64(1) << bit); tail_b[i][chunk] |= (UInt64(1) << bit) end
                end
            end
            
            # BJMM Structures tracking (win2_a, win2_b, info_mask_a, info_mask_b, msg)
            # In GF4, we must track masks for both basis coordinates (a and b) to ensure exact weight evaluation
            H1 = Dict{Tuple{UInt64, UInt64}, Vector{Tuple{UInt64, UInt64, UInt64, UInt64, Vector{Tuple{Int, Int}}}}}()
            H12 = Dict{Tuple{UInt64, UInt64}, Vector{Tuple{UInt64, UInt64, Vector{Tuple{Int, Int}}}}}()
            H3 = Dict{Tuple{UInt64, UInt64}, Vector{Tuple{UInt64, UInt64, UInt64, UInt64, Vector{Tuple{Int, Int}}}}}()
            
            function _build_base!(depth, picked, Range, p_tgt, cw1a, cw1b, cw2a, cw2b, cia, cib, msg, tgt_Dict)
                if picked == p_tgt
                    push!(get!(tgt_Dict, (cw1a, cw1b), []), (cw2a, cw2b, cia, cib, copy(msg)))
                    return
                end
                if depth > length(Range) || (length(Range) - depth + 1) < (p_tgt - picked) return end
                
                _build_base!(depth+1, picked, Range, p_tgt, cw1a, cw1b, cw2a, cw2b, cia, cib, msg, tgt_Dict)
                
                idx = Range[depth]; bit_mask = (UInt64(1) << ((idx-1)%64)); chunk_offset = (idx-1)÷64
                # We assume small enough k for single UInt64 info masks here for speed, or mask arrays if k > 64
                
                push!(msg, (idx, 1))
                _build_base!(depth+1, picked+1, Range, p_tgt, cw1a ⊻ w1_a[idx], cw1b ⊻ w1_b[idx], cw2a ⊻ w2_a[idx], cw2b ⊻ w2_b[idx], cia, cib | bit_mask, msg, tgt_Dict)
                pop!(msg)
                
                push!(msg, (idx, 2))
                _build_base!(depth+1, picked+1, Range, p_tgt, cw1a ⊻ w1_a[idx] ⊻ w1_b[idx], cw1b ⊻ w1_a[idx], cw2a ⊻ w2_a[idx] ⊻ w2_b[idx], cw2b ⊻ w2_a[idx], cia | bit_mask, cib, msg, tgt_Dict)
                pop!(msg)
                
                push!(msg, (idx, 3))
                _build_base!(depth+1, picked+1, Range, p_tgt, cw1a ⊻ w1_b[idx], cw1b ⊻ w1_a[idx] ⊻ w1_b[idx], cw2a ⊻ w2_b[idx], cw2b ⊻ w2_a[idx] ⊻ w2_b[idx], cia | bit_mask, cib | bit_mask, msg, tgt_Dict)
                pop!(msg)
            end
            
            _build_base!(1, 0, R1, base_wt, UInt64(0), UInt64(0), UInt64(0), UInt64(0), UInt64(0), UInt64(0), Tuple{Int, Int}[], H1)
            
            function _build_H12!(depth, picked, cw1a, cw1b, cw2a, cw2b, cia, cib, msg2)
                if picked == base_wt
                    if haskey(H1, (cw1a, cw1b))
                        for (w2a_1, w2b_1, cia_1, cib_1, msg1) in H1[(cw1a, cw1b)]
                            # True GF(4) Hamming weight calculation of the XOR info masks sum
                            info_wt = count_ones((cia_1 ⊻ cia) | (cib_1 ⊻ cib))
                            if info_wt == lvl1_target_wt
                                push!(get!(H12, (cw2a ⊻ w2a_1, cur_w2b = cw2b ⊻ w2b_1), []), (cia_1 ⊻ cia, cib_1 ⊻ cib, vcat(msg1, msg2)))
                            end
                        end
                    end
                    return
                end
                if depth > length(R2) || (length(R2) - depth + 1) < (base_wt - picked) return end
                
                _build_H12!(depth+1, picked, cw1a, cw1b, cw2a, cw2b, cia, cib, msg2)
                idx = R2[depth]; bit_mask = (UInt64(1) << ((idx-1)%64))
                
                push!(msg2, (idx, 1)); _build_H12!(depth+1, picked+1, cw1a ⊻ w1_a[idx], cw1b ⊻ w1_b[idx], cw2a ⊻ w2_a[idx], cw2b ⊻ w2_b[idx], cia, cib | bit_mask, msg2); pop!(msg2)
                push!(msg2, (idx, 2)); _build_H12!(depth+1, picked+1, cw1a ⊻ w1_a[idx] ⊻ w1_b[idx], cw1b ⊻ w1_a[idx], cw2a ⊻ w2_a[idx] ⊻ w2_b[idx], cw2b ⊻ w2_a[idx], cia | bit_mask, cib, msg2); pop!(msg2)
                push!(msg2, (idx, 3)); _build_H12!(depth+1, picked+1, cw1a ⊻ w1_b[idx], cw1b ⊻ w1_a[idx] ⊻ w1_b[idx], cw2a ⊻ w2_b[idx], cw2b ⊻ w2_a[idx] ⊻ w2_b[idx], cia | bit_mask, cib | bit_mask, msg2); pop!(msg2)
            end
            _build_H12!(1, 0, UInt64(0), UInt64(0), UInt64(0), UInt64(0), UInt64(0), UInt64(0), Tuple{Int, Int}[])
            empty!(H1)
            
            _build_base!(1, 0, R3, base_wt, UInt64(0), UInt64(0), UInt64(0), UInt64(0), UInt64(0), UInt64(0), Tuple{Int, Int}[], H3)
            
            function _probe_H34!(depth, picked, cw1a, cw1b, cw2a, cw2b, cia, cib, msg4)
                if !keep_going[] return end
                if picked == base_wt
                    if haskey(H3, (cw1a, cw1b))
                        for (w2a_3, w2b_3, cia_3, cib_3, msg3) in H3[(cw1a, cw1b)]
                            w2a_34, w2b_34 = cw2a ⊻ w2a_3, cw2b ⊻ w2b_3
                            if haskey(H12, (w2a_34, w2b_34))
                                cia_34, cib_34 = cia_3 ⊻ cia, cib_3 ⊻ cib
                                for (cia_12, cib_12, msg12) in H12[(w2a_34, w2b_34)]
                                    if count_ones((cia_12 ⊻ cia_34) | (cib_12 ⊻ cib_34)) == p
                                        fill!(tail_buf_a, UInt64(0)); fill!(tail_buf_b, UInt64(0))
                                        for sub_msg in (msg12, msg3, msg4)
                                            for (idx, v) in sub_msg
                                                if v == 1 tail_buf_a .⊻= tail_a[idx]; tail_buf_b .⊻= tail_b[idx]
                                                elseif v == 2 tail_buf_a .⊻= tail_a[idx] .⊻ tail_b[idx]; tail_buf_b .⊻= tail_a[idx]
                                                elseif v == 3 tail_buf_a .⊻= tail_b[idx]; tail_buf_b .⊻= tail_a[idx] .⊻ tail_b[idx] end
                                            end
                                        end
                                        
                                        if sum(count_ones.(tail_buf_a .| tail_buf_b)) + p <= target_w
                                            lock(results_lock) do
                                                if length(found_vectors) < num_find
                                                    e_loc = zeros(F, n)
                                                    final_a, final_b = cia_12 ⊻ cia_34, cib_12 ⊻ cib_34
                                                    for i in 1:k
                                                        ba = (final_a >> ((i-1)%64)) & 1; bb = (final_b >> ((i-1)%64)) & 1
                                                        if ba == 1 && bb == 0 e_loc[i] = ω
                                                        elseif ba == 1 && bb == 1 e_loc[i] = ω + 1
                                                        elseif ba == 0 && bb == 1 e_loc[i] = F(1) end
                                                    end
                                                    for j in 1:tail_len
                                                        chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                                        bit_a = (tail_buf_a[chunk] >> bit) & 1; bit_b = (tail_buf_b[chunk] >> bit) & 1
                                                        if bit_a == 1 && bit_b == 0 e_loc[k + l + j] = ω
                                                        elseif bit_a == 1 && bit_b == 1 e_loc[k + l + j] = ω + 1
                                                        elseif bit_a == 0 && bit_b == 1 e_loc[k + l + j] = F(1) end
                                                    end
                                                    inv_p = invperm(σ_loc)
                                                    push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                                                    if length(found_vectors) >= num_find keep_going[] = false end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    return
                end
                if depth > length(R4) || (length(R4) - depth + 1) < (base_wt - picked) return end
                
                _probe_H34!(depth+1, picked, cw1a, cw1b, cw2a, cw2b, cia, cib, msg4)
                idx = R4[depth]; bit_mask = (UInt64(1) << ((idx-1)%64))
                
                push!(msg4, (idx, 1)); _probe_H34!(depth+1, picked+1, cw1a ⊻ w1_a[idx], cw1b ⊻ w1_b[idx], cw2a ⊻ w2_a[idx], cw2b ⊻ w2_b[idx], cia, cib | bit_mask, msg4); pop!(msg4)
                push!(msg4, (idx, 2)); _probe_H34!(depth+1, picked+1, cw1a ⊻ w1_a[idx] ⊻ w1_b[idx], cw1b ⊻ w1_a[idx], cw2a ⊻ w2_a[idx] ⊻ w2_b[idx], cw2b ⊻ w2_a[idx], cia | bit_mask, cib, msg4); pop!(msg4)
                push!(msg4, (idx, 3)); _probe_H34!(depth+1, picked+1, cw1a ⊻ w1_b[idx], cw1b ⊻ w1_a[idx] ⊻ w1_b[idx], cw2a ⊻ w2_b[idx], cw2b ⊻ w2_a[idx] ⊻ w2_b[idx], cia | bit_mask, cib | bit_mask, msg4); pop!(msg4)
            end
            _probe_H34!(1, 0, UInt64(0), UInt64(0), UInt64(0), UInt64(0), UInt64(0), UInt64(0), Tuple{Int, Int}[])
        end
    end
    return found_vectors
end

"""
    _BJMM_minimum_distance_GF3(G::CTMatrixTypes, target_w::Int; p::Int = 4, ϵ1::Int = 1, l1::Int = 10, l2::Int = 10, num_find::Int = 1, max_iters::Int = 10000)
"""
function _BJMM_minimum_distance_GF3(G::CTMatrixTypes, target_w::Int; p::Int = 4, ϵ1::Int = 1, l1::Int = 10, l2::Int = 10, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    l = l1 + l2
    
    F = base_ring(G); p_char, d_deg = Int(characteristic(F)), degree(F)
    base_wt = (p ÷ 4) + ϵ1; lvl1_target_wt = p ÷ 2
    
    num_thrds = Threads.nthreads()
    thread_load = max_iters ÷ num_thrds; remaining = max_iters % num_thrds
    keep_going = Threads.Atomic{Bool}(true); results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    @inline function add_mod3(AH::UInt64, AL::UInt64, BH::UInt64, BL::UInt64)
        SH = ((AH ⊻ BH) & ~(AL | BL)) | (AL & BL)
        SL = ((AL ⊻ BL) & ~(AH | BH)) | (AH & BH)
        return SH, SL
    end
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        
        q1 = k ÷ 4; q2 = k ÷ 2; q3 = 3 * k ÷ 4
        R1 = 1:q1; R2 = (q1+1):q2; R3 = (q2+1):q3; R4 = (q3+1):k
        tail_len = n - k - l; num_tail_chunks = cld(tail_len, 64)
        tail_buf_H = zeros(UInt64, num_tail_chunks); tail_buf_L = zeros(UInt64, num_tail_chunks)
        
        for _ in 1:iters
            if !keep_going[] break end
            
            σ_loc = collect(1:n); shuffle!(σ_loc); G_loc = G[:, σ_loc]
            try _make_systematic_gf!(G_loc, σ_loc, k) catch; continue end
            
            # --- BITSLICING COMPILER ---
            w1_H = zeros(UInt64, k); w1_L = zeros(UInt64, k)
            w2_H = zeros(UInt64, k); w2_L = zeros(UInt64, k)
            tail_H = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            tail_L = [zeros(UInt64, num_tail_chunks) for _ in 1:k]
            
            for i in 1:k
                for j in 1:l1
                    v = G_loc[i, k + j]
                    if v == 1 w1_L[i] |= (UInt64(1) << (j - 1)) elseif v == 2 w1_H[i] |= (UInt64(1) << (j - 1)) end
                end
                for j in 1:l2
                    v = G_loc[i, k + l1 + j]
                    if v == 1 w2_L[i] |= (UInt64(1) << (j - 1)) elseif v == 2 w2_H[i] |= (UInt64(1) << (j - 1)) end
                end
                for j in 1:tail_len
                    v = G_loc[i, k + l + j]; chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                    if v == 1 tail_L[i][chunk] |= (UInt64(1) << bit) elseif v == 2 tail_H[i][chunk] |= (UInt64(1) << bit) end
                end
            end
            
            # (win2_H, win2_L, info_H, info_L, msg)
            H1 = Dict{Tuple{UInt64, UInt64}, Vector{Tuple{UInt64, UInt64, UInt64, UInt64, Vector{Tuple{Int, Int}}}}}()
            H12 = Dict{Tuple{UInt64, UInt64}, Vector{Tuple{UInt64, UInt64, Vector{Tuple{Int, Int}}}}}()
            H3 = Dict{Tuple{UInt64, UInt64}, Vector{Tuple{UInt64, UInt64, UInt64, UInt64, Vector{Tuple{Int, Int}}}}}()
            
            function _build_base!(depth, picked, Range, p_tgt, cw1H, cw1L, cw2H, cw2L, ciH, ciL, msg, tgt_Dict)
                if picked == p_tgt
                    push!(get!(tgt_Dict, (cw1H, cw1L), []), (cw2H, cw2L, ciH, ciL, copy(msg)))
                    return
                end
                if depth > length(Range) || (length(Range) - depth + 1) < (p_tgt - picked) return end
                
                _build_base!(depth+1, picked, Range, p_tgt, cw1H, cw1L, cw2H, cw2L, ciH, ciL, msg, tgt_Dict)
                idx = Range[depth]; bit_mask = (UInt64(1) << ((idx-1)%64))
                
                n1H, n1L = add_mod3(cw1H, cw1L, w1_H[idx], w1_L[idx]); n2H, n2L = add_mod3(cw2H, cw2L, w2_H[idx], w2_L[idx])
                push!(msg, (idx, 1)); _build_base!(depth+1, picked+1, Range, p_tgt, n1H, n1L, n2H, n2L, ciH, ciL | bit_mask, msg, tgt_Dict); pop!(msg)
                
                n1H, n1L = add_mod3(cw1H, cw1L, w1_L[idx], w1_H[idx]); n2H, n2L = add_mod3(cw2H, cw2L, w2_L[idx], w2_H[idx])
                push!(msg, (idx, 2)); _build_base!(depth+1, picked+1, Range, p_tgt, n1H, n1L, n2H, n2L, ciH | bit_mask, ciL, msg, tgt_Dict); pop!(msg)
            end
            
            _build_base!(1, 0, R1, base_wt, UInt64(0), UInt64(0), UInt64(0), UInt64(0), UInt64(0), UInt64(0), Tuple{Int, Int}[], H1)
            
            function _build_H12!(depth, picked, cw1H, cw1L, cw2H, cw2L, ciH, ciL, msg2)
                if picked == base_wt
                    # Inverse Lookup for GF(3): X + Y = 0 => X = 2Y (swap H and L)
                    if haskey(H1, (cw1L, cw1H))
                        for (w2H_1, w2L_1, ciH_1, ciL_1, msg1) in H1[(cw1L, cw1H)]
                            # Add information masks modulo 3 natively to evaluate resulting weight
                            res_iH, res_iL = add_mod3(ciH_1, ciL_1, ciH, ciL)
                            if count_ones(res_iH | res_iL) == lvl1_target_wt
                                nH, nL = add_mod3(cw2H, cw2L, w2H_1, w2L_1)
                                push!(get!(H12, (nH, nL), []), (res_iH, res_iL, vcat(msg1, msg2)))
                            end
                        end
                    end
                    return
                end
                if depth > length(R2) || (length(R2) - depth + 1) < (base_wt - picked) return end
                
                _build_H12!(depth+1, picked, cw1H, cw1L, cw2H, cw2L, ciH, ciL, msg2)
                idx = R2[depth]; bit_mask = (UInt64(1) << ((idx-1)%64))
                
                n1H, n1L = add_mod3(cw1H, cw1L, w1_H[idx], w1_L[idx]); n2H, n2L = add_mod3(cw2H, cw2L, w2_H[idx], w2_L[idx])
                push!(msg2, (idx, 1)); _build_H12!(depth+1, picked+1, n1H, n1L, n2H, n2L, ciH, ciL | bit_mask, msg2); pop!(msg2)
                
                n1H, n1L = add_mod3(cw1H, cw1L, w1_L[idx], w1_H[idx]); n2H, n2L = add_mod3(cw2H, cw2L, w2_L[idx], w2_H[idx])
                push!(msg2, (idx, 2)); _build_H12!(depth+1, picked+1, n1H, n1L, n2H, n2L, ciH | bit_mask, ciL, msg2); pop!(msg2)
            end
            _build_H12!(1, 0, UInt64(0), UInt64(0), UInt64(0), UInt64(0), UInt64(0), UInt64(0), Tuple{Int, Int}[])
            empty!(H1)
            
            _build_base!(1, 0, R3, base_wt, UInt64(0), UInt64(0), UInt64(0), UInt64(0), UInt64(0), UInt64(0), Tuple{Int, Int}[], H3)
            
            function _probe_H34!(depth, picked, cw1H, cw1L, cw2H, cw2L, ciH, ciL, msg4)
                if !keep_going[] return end
                if picked == base_wt
                    if haskey(H3, (cw1L, cw1H))
                        for (w2H_3, w2L_3, ciH_3, ciL_3, msg3) in H3[(cw1L, cw1H)]
                            w2H_34, w2L_34 = add_mod3(cw2H, cw2L, w2H_3, w2L_3)
                            if haskey(H12, (w2L_34, w2H_34))
                                ciH_34, ciL_34 = add_mod3(ciH_3, ciL_3, ciH, ciL)
                                for (ciH_12, ciL_12, msg12) in H12[(w2L_34, w2H_34)]
                                    final_H, final_L = add_mod3(ciH_12, ciL_12, ciH_34, ciL_34)
                                    if count_ones(final_H | final_L) == p
                                        fill!(tail_buf_H, UInt64(0)); fill!(tail_buf_L, UInt64(0))
                                        for sub_msg in (msg12, msg3, msg4)
                                            for (idx, v) in sub_msg
                                                RH = v == 1 ? tail_H[idx] : tail_L[idx]
                                                RL = v == 1 ? tail_L[idx] : tail_H[idx]
                                                for c in 1:num_tail_chunks tail_buf_H[c], tail_buf_L[c] = add_mod3(tail_buf_H[c], tail_buf_L[c], RH[c], RL[c]) end
                                            end
                                        end
                                        
                                        if sum(count_ones.(tail_buf_H .| tail_buf_L)) + p <= target_w
                                            lock(results_lock) do
                                                if length(found_vectors) < num_find
                                                    e_loc = zeros(F, n)
                                                    for i in 1:k
                                                        if ((final_H >> ((i-1)%64)) & 1) == 1 e_loc[i] = F(2)
                                                        elseif ((final_L >> ((i-1)%64)) & 1) == 1 e_loc[i] = F(1) end
                                                    end
                                                    for j in 1:tail_len
                                                        chunk, bit = (j - 1) ÷ 64 + 1, (j - 1) % 64
                                                        if ((tail_buf_H[chunk] >> bit) & 1) == 1 e_loc[k + l + j] = F(2)
                                                        elseif ((tail_buf_L[chunk] >> bit) & 1) == 1 e_loc[k + l + j] = F(1) end
                                                    end
                                                    inv_p = invperm(σ_loc)
                                                    push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                                                    if length(found_vectors) >= num_find keep_going[] = false end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    return
                end
                if depth > length(R4) || (length(R4) - depth + 1) < (base_wt - picked) return end
                
                _probe_H34!(depth+1, picked, cw1H, cw1L, cw2H, cw2L, ciH, ciL, msg4)
                idx = R4[depth]; bit_mask = (UInt64(1) << ((idx-1)%64))
                
                n1H, n1L = add_mod3(cw1H, cw1L, w1_H[idx], w1_L[idx]); n2H, n2L = add_mod3(cw2H, cw2L, w2_H[idx], w2_L[idx])
                push!(msg4, (idx, 1)); _probe_H34!(depth+1, picked+1, n1H, n1L, n2H, n2L, ciH, ciL | bit_mask, msg4); pop!(msg4)
                
                n1H, n1L = add_mod3(cw1H, cw1L, w1_L[idx], w1_H[idx]); n2H, n2L = add_mod3(cw2H, cw2L, w2_L[idx], w2_H[idx])
                push!(msg4, (idx, 2)); _probe_H34!(depth+1, picked+1, n1H, n1L, n2H, n2L, ciH | bit_mask, ciL, msg4); pop!(msg4)
            end
            _probe_H34!(1, 0, UInt64(0), UInt64(0), UInt64(0), UInt64(0), UInt64(0), UInt64(0), Tuple{Int, Int}[])
        end
    end
    return found_vectors
end

function _BJMM_minimum_distance_nonbinary(G::CTMatrixTypes, target_w::Int; p::Int = 4, ϵ1::Int = 1, l1::Int = 3, l2::Int = 3, num_find::Int = 1, max_iters::Int = 10000)
    k, n = size(G)
    l = l1 + l2
    
    F = base_ring(G); p_char, d_deg = Int(characteristic(F)), degree(F)
    non_zeros = filter(!iszero, collect(F))
    
    base_wt = (p ÷ 4) + ϵ1; lvl1_target_wt = p ÷ 2
    
    num_thrds = Threads.nthreads(); thread_load = max_iters ÷ num_thrds; remaining = max_iters % num_thrds
    keep_going = Threads.Atomic{Bool}(true); results_lock = Threads.SpinLock()
    found_vectors = Set{Vector{Int}}()
    
    Threads.@threads for th in 1:num_thrds
        iters = thread_load + (th <= remaining ? 1 : 0)
        q1 = k ÷ 4; q2 = k ÷ 2; q3 = 3 * k ÷ 4
        R1 = 1:q1; R2 = (q1+1):q2; R3 = (q2+1):q3; R4 = (q3+1):k
        tail_len = n - k - l
        
        for _ in 1:iters
            if !keep_going[] break end
            
            perm = collect(1:n); shuffle!(perm); G_loc = G[:, perm]
            try _make_systematic_gf!(G_loc, perm, k) catch; continue end
            
            w1_rows = [G_loc[i, (k+1):(k+l1)] for i in 1:k]
            w2_rows = [G_loc[i, (k+l1+1):(k+l)] for i in 1:k]
            tail_rows = [G_loc[i, (k+l+1):n] for i in 1:k]
            
            # (win2, info_vector, msg)
            H1 = Dict{Vector{typeof(zero(F))}, Vector{Tuple{Vector{typeof(zero(F))}, Vector{typeof(zero(F))}, Vector{Tuple{Int, typeof(zero(F))}}}}}()
            H12 = Dict{Vector{typeof(zero(F))}, Vector{Tuple{Vector{typeof(zero(F))}, Vector{Tuple{Int, typeof(zero(F))}}}}}()
            H3 = Dict{Vector{typeof(zero(F))}, Vector{Tuple{Vector{typeof(zero(F))}, Vector{typeof(zero(F))}, Vector{Tuple{Int, typeof(zero(F))}}}}}()
            
            function _build_base!(depth, picked, Range, p_tgt, cur_w1, cur_w2, cur_info, msg, tgt_Dict)
                if picked == p_tgt
                    push!(get!(tgt_Dict, cur_w1, []), (cur_w2, copy(cur_info), copy(msg)))
                    return
                end
                if depth > length(Range) || (length(Range) - depth + 1) < (p_tgt - picked) return end
                
                _build_base!(depth+1, picked, Range, p_tgt, cur_w1, cur_w2, cur_info, msg, tgt_Dict)
                
                idx = Range[depth]
                for scalar in non_zeros
                    n1 = cur_w1 .+ scalar .* w1_rows[idx]; n2 = cur_w2 .+ scalar .* w2_rows[idx]
                    next_info = copy(cur_info); next_info[idx] = scalar
                    push!(msg, (idx, scalar))
                    _build_base!(depth+1, picked+1, Range, p_tgt, n1, n2, next_info, msg, tgt_Dict)
                    pop!(msg)
                end
            end
            
            _build_base!(1, 0, R1, base_wt, zeros(F, l1), zeros(F, l2), zeros(F, k), Tuple{Int, typeof(zero(F))}[], H1)
            
            function _build_H12!(depth, picked, cur_w1, cur_w2, cur_info, msg2)
                if picked == base_wt
                    neg_w1 = -cur_w1
                    if haskey(H1, neg_w1)
                        for (w2_1, info_1, msg1) in H1[neg_w1]
                            res_info = info_1 .+ cur_info
                            if count(!iszero, res_info) == lvl1_target_wt
                                push!(get!(H12, cur_w2 .+ w2_1, []), (res_info, vcat(msg1, msg2)))
                            end
                        end
                    end
                    return
                end
                if depth > length(R2) || (length(R2) - depth + 1) < (base_wt - picked) return end
                
                _build_H12!(depth+1, picked, cur_w1, cur_w2, cur_info, msg2)
                idx = R2[depth]
                for scalar in non_zeros
                    next_info = copy(cur_info); next_info[idx] = scalar
                    push!(msg2, (idx, scalar))
                    _build_H12!(depth+1, picked+1, cur_w1 .+ scalar .* w1_rows[idx], cur_w2 .+ scalar .* w2_rows[idx], next_info, msg2)
                    pop!(msg2)
                end
            end
            _build_H12!(1, 0, zeros(F, l1), zeros(F, l2), zeros(F, k), Tuple{Int, typeof(zero(F))}[])
            empty!(H1)
            
            _build_base!(1, 0, R3, base_wt, zeros(F, l1), zeros(F, l2), zeros(F, k), Tuple{Int, typeof(zero(F))}[], H3)
            
            function _probe_H34!(depth, picked, cur_w1, cur_w2, cur_info, msg4)
                if !keep_going[] return end
                if picked == base_wt
                    neg_w1 = -cur_w1
                    if haskey(H3, neg_w1)
                        for (w2_3, info_3, msg3) in H3[neg_w1]
                            w2_34 = cur_w2 .+ w2_3
                            neg_w2_34 = -w2_34
                            if haskey(H12, neg_w2_34)
                                info_34 = info_3 .+ cur_info
                                for (info_12, msg12) in H12[neg_w2_34]
                                    final_info = info_12 .+ info_34
                                    if count(!iszero, final_info) == p
                                        tail_buf = zeros(F, tail_len)
                                        for sub_msg in (msg12, msg3, msg4)
                                            for (idx, sc) in sub_msg tail_buf .+= sc .* tail_rows[idx] end
                                        end
                                        if count(!iszero, tail_buf) + p <= target_w
                                            lock(results_lock) do
                                                if length(found_vectors) < num_find
                                                    e_loc = zeros(F, n)
                                                    e_loc[1:k] .= final_info
                                                    e_loc[(k+l+1):n] .= tail_buf
                                                    inv_p = invperm(perm)
                                                    push!(found_vectors, [_pack_field_elem(e_loc[inv_p[j]], p_char, d_deg) for j in 1:n])
                                                    if length(found_vectors) >= num_find keep_going[] = false end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    return
                end
                if depth > length(R4) || (length(R4) - depth + 1) < (base_wt - picked) return end
                
                _probe_H34!(depth+1, picked, cur_w1, cur_w2, cur_info, msg4)
                idx = R4[depth]
                for scalar in non_zeros
                    next_info = copy(cur_info); next_info[idx] = scalar
                    push!(msg4, (idx, scalar))
                    _probe_H34!(depth+1, picked+1, cur_w1 .+ scalar .* w1_rows[idx], cur_w2 .+ scalar .* w2_rows[idx], next_info, msg4)
                    pop!(msg4)
                end
            end
            _probe_H34!(1, 0, zeros(F, l1), zeros(F, l2), zeros(F, k), Tuple{Int, typeof(zero(F))}[])
        end
    end
    return found_vectors
end
