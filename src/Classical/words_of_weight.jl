# Copyright (c) 2022 - 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

function _is_canonical(msg_packed::UInt128, auts::Vector{Vector{Int}}, depth::Int)
    if isempty(auts) return true end
    
    # msg_packed contains the bit values at indices 1 through depth.
    # Bits are stored such that the i-th bit is at (msg_packed >> (i - 1)) & 1.
    
    for σ in auts
        is_smaller_or_equal = true
        
        for i in 1:depth
            orig_bit = (msg_packed >> (i - 1)) & 1
            
            # Where does this index map to under the automorphism?
            mapped_idx = σ[i]
            
            # If the automorphism maps it outside our currently decided bits, 
            # we can't fully evaluate this permutation yet.
            if mapped_idx > depth
                continue
            end
            
            mapped_bit = (msg_packed >> (mapped_idx - 1)) & 1
            
            # Lexicographical comparison
            if mapped_bit < orig_bit
                is_smaller_or_equal = false
                break
            elseif mapped_bit > orig_bit
                break # This permutation is strictly greater, so current msg is safe
            end
        end
        
        if !is_smaller_or_equal
            return false # We found a permutation that yields a smaller representation
        end
    end
    
    return true
end

# GF(4)
function _is_canonical_q(msg_packed::UInt128, auts::Vector{Vector{Int}}, depth::Int)
    if isempty(auts) return true end
    
    for σ in auts
        is_smaller_or_equal = true
        
        for i in 1:depth
            orig_val = (msg_packed >> (2 * (i - 1))) & 3
            
            mapped_idx = σ[i]
            if mapped_idx > depth
                continue
            end
            
            mapped_val = (msg_packed >> (2 * (mapped_idx - 1))) & 3
            
            if mapped_val < orig_val
                is_smaller_or_equal = false
                break
            elseif mapped_val > orig_val
                break 
            end
        end
        
        if !is_smaller_or_equal
            return false 
        end
    end
    
    return true
end

# GF(3)
function _is_canonical_q(msg_L::UInt128, msg_H::UInt128, auts::Vector{Vector{Int}}, depth::Int)
    if isempty(auts) return true end
    
    for σ in auts
        is_smaller_or_equal = true
        
        for i in 1:depth
            # Recombine Low and High bits into an integer 0, 1, or 2
            orig_L = (msg_L >> (i - 1)) & 1
            orig_H = (msg_H >> (i - 1)) & 1
            orig_val = orig_L | (orig_H << 1) 
            
            mapped_idx = σ[i]
            if mapped_idx > depth
                continue
            end
            
            map_L = (msg_L >> (mapped_idx - 1)) & 1
            map_H = (msg_H >> (mapped_idx - 1)) & 1
            mapped_val = map_L | (map_H << 1)
            
            if mapped_val < orig_val
                is_smaller_or_equal = false
                break
            elseif mapped_val > orig_val
                break
            end
        end
        
        if !is_smaller_or_equal
            return false
        end
    end
    
    return true
end

"""
Calculates the Hamming weight of a GF(4) state packed into a UInt128.
A 2-bit chunk is non-zero if either its left bit OR right bit is 1.
"""
@inline function _weight_gf4(x::UInt128)
    # Mask out the odd bits and even bits
    odd_bits  = x & 0x55555555555555555555555555555555
    even_bits = (x >> 1) & 0x55555555555555555555555555555555
    # If a chunk had any 1s, the OR will have a 1 in that position.
    return count_ones(odd_bits | even_bits)
end

"""
    _expand_orbits(representatives::Vector{Vector{T}}, auts::Vector{Vector{Int}}, q::Int=2, mul_table::Matrix{Int}=Matrix{Int}(undef,0,0)) where T

Expands orbit representatives using coordinate permutations and scalar multiplications.
"""
function _expand_orbits(representatives::Vector{Vector{T}}, auts::Vector{Vector{Int}}, q::Int=2, mul_table::Matrix{Int}=Matrix{Int}(undef,0,0)) where T
    if isempty(auts) && q == 2
        return representatives
    end

    all_words = Set{Vector{T}}()
    
    for rep in representatives
        queue = [rep]
        push!(all_words, rep)
        
        head = 1
        while head <= length(queue)
            curr = queue[head]
            head += 1
            
            # 1. Apply Coordinate Permutations
            for σ in auts
                shifted = curr[σ] 
                if !(shifted in all_words)
                    push!(all_words, shifted)
                    push!(queue, shifted)
                end
            end
            
            # 2. Apply Scalar Multiples (for non-binary)
            if q > 2 && !isempty(mul_table)
                # (Scalar mapping logic goes here for GF(3)/GF(4))
            end
        end
    end
    
    return collect(all_words)
end

function _prepare_packed_brouwer_binary(G_sys::Matrix{T}) where T
    k, n = size(G_sys)
    n_tail = n - k
    A_raw = G_sys[:, (k + 1):n]
    
    # 1. Pack A_raw into UInt128 registers
    A_rows = zeros(UInt128, k)
    for row in 1:k
        val = UInt128(0)
        for col in 1:n_tail
            if !iszero(A_raw[row, col]) 
                val |= (UInt128(1) << (col - 1)) 
            end
        end
        A_rows[row] = val
    end
    
    # 2. Compute the Lower Bound Table (LBT)
    # lbt[i+1] bounds the minimum weight of a combination of the LAST i rows.
    lbt = zeros(Int, k + 1)
    
    for i in 1:k
        # We look at the submatrix of the last i rows
        # A simple greedy bound: the minimum weight of any non-zero row in that block.
        # (For production, you can use a mini-DFS here to get the exact minimum distance 
        # of the subcode, which prunes the main search even faster).
        min_wt = n_tail
        for row in (k - i + 1):k
            wt = count_ones(A_rows[row])
            if wt > 0 && wt < min_wt
                min_wt = wt
            end
        end
        lbt[i + 1] = min_wt == n_tail ? 0 : min_wt
    end
    
    # Ensure monotonicity
    for i in 2:(k + 1)
        if lbt[i] < lbt[i - 1]
            lbt[i] = lbt[i - 1]
        end
    end
    
    return A_rows, lbt
end

function _minimum_words_binary_recursive!(
    A_rows::Vector{UInt128}, 
    depth::Int, 
    picked::Int, 
    curr_tail::UInt128, 
    results::Vector{Tuple{UInt128, UInt128}}, 
    lbt::Vector{Int}, 
    msg_so_far::UInt128, 
    spawn_depth::Int, 
    auts::Vector{Vector{Int}}, 
    res_lock::ReentrantLock,
    global_min_d::Threads.Atomic{Int},
    k::Int
)
    # 1. DYNAMIC WEIGHT PRUNING
    tw = count_ones(curr_tail)
    current_min = global_min_d[]
    
    # Prune if the absolute minimum weight this branch can produce exceeds the known minimum
    if (picked + tw + lbt[k - depth + 1]) > current_min
        return 
    end

    # 2. SYMMETRY PRUNING
    if picked > 1 && !isempty(auts)
        if !_is_canonical(msg_so_far, auts, depth) return end
    end

    # 3. BASE CASE: Full codeword reached
    if depth == k
        total_w = picked + tw
        if total_w > 0 && total_w <= global_min_d[]
            lock(res_lock) do
                # Double-check inside the lock to avoid race conditions
                if total_w < global_min_d[]
                    empty!(results) # We found a new global minimum! Discard old words.
                    global_min_d[] = total_w
                    push!(results, (msg_so_far, curr_tail))
                elseif total_w == global_min_d[]
                    push!(results, (msg_so_far, curr_tail))
                end
            end
        end
        return
    end

    # 4. HARDWARE-ACCELERATED BRANCHING
    if depth < spawn_depth
        # Branch 1: Include bit
        t = Threads.@spawn _minimum_words_binary_recursive!(
            A_rows, depth + 1, picked + 1, curr_tail ⊻ A_rows[depth + 1], 
            results, lbt, msg_so_far | (UInt128(1) << depth), spawn_depth, auts, res_lock, global_min_d, k
        )
        
        # Branch 2: Skip bit
        _minimum_words_binary_recursive!(
            A_rows, depth + 1, picked, curr_tail, 
            results, lbt, msg_so_far, spawn_depth, auts, res_lock, global_min_d, k
        )
        wait(t)
    else
        # Serial processing
        _minimum_words_binary_recursive!(
            A_rows, depth + 1, picked + 1, curr_tail ⊻ A_rows[depth + 1], 
            results, lbt, msg_so_far | (UInt128(1) << depth), spawn_depth, auts, res_lock, global_min_d, k
        )
        _minimum_words_binary_recursive!(
            A_rows, depth + 1, picked, curr_tail, 
            results, lbt, msg_so_far, spawn_depth, auts, res_lock, global_min_d, k
        )
    end
end

function _words_of_minimum_weight_binary(C::AbstractLinearCode; expand::Bool=true, verbose::Bool=false)
    k, n = C.k, C.n
    @assert n - k <= 128 "Parity tail exceeds 128 bits; requires chunked bit-packing."
    
    # 1. Standardize Matrix & Prepare LBT
    G_sys = Array(generator_matrix(C))
    A_rows, lbt = _prepare_packed_brouwer_binary(G_sys)
    
    # 2. Symmetry handling (assume generate_automorphisms is defined)
    auts = Vector{Vector{Int}}() 
    
    # 3. Setup thread-safe dynamic threshold
    global_min_d = Threads.Atomic{Int}(n)
    raw_results = Vector{Tuple{UInt128, UInt128}}()
    res_lock = ReentrantLock()
    
    verbose && println("Executing dynamic minimum words search...")
    _minimum_words_binary_recursive!(
        A_rows, 0, 0, UInt128(0), 
        raw_results, lbt, UInt128(0), 5, auts, res_lock, global_min_d, k
    )

    final_d = global_min_d[]
    verbose && println("Search complete. Minimum weight: $final_d. Found $(length(raw_results)) representative(s).")

    # 4. Unpack to standard Vectors
    F = C.F
    T_zero, T_one = zero(F), one(F)
    final_words = Vector{Vector{typeof(T_zero)}}()
    
    for (msg, tail) in raw_results
        word = fill(T_zero, n)
        for i in 1:k
            if ((msg >> (i - 1)) & 1) == 1 word[i] = T_one end
        end
        for i in 1:(n - k)
            if ((tail >> (i - 1)) & 1) == 1 word[k + i] = T_one end
        end
        push!(final_words, word)
    end
    
    if expand
        return final_d, _expand_orbits(final_words, auts, 2)
    end
    
    return final_d, final_words
end

function _words_of_minimum_weight_quaternary_recursive!(
    A_scaled::Matrix{UInt128}, 
    depth::Int, 
    picked::Int, 
    curr_tail::UInt128, 
    results::Vector{Tuple{UInt128, UInt128}}, 
    lbt::Vector{Int}, 
    msg_so_far::UInt128, 
    spawn_depth::Int, 
    auts::Vector{Vector{Int}}, 
    res_lock::ReentrantLock,
    global_min_d::Threads.Atomic{Int},
    k::Int
)
    # 1. DYNAMIC WEIGHT PRUNING
    tw = _weight_gf4(curr_tail)
    
    # Prune if the absolute minimum weight this branch can produce exceeds the known minimum
    if (picked + tw + lbt[k - depth + 1]) > global_min_d[]
        return 
    end

    # 2. SYMMETRY PRUNING 
    if picked > 1 && !isempty(auts)
        if !_is_canonical_q(msg_so_far, auts, depth) return end
    end

    # 3. BASE CASE: Full codeword reached
    if depth == k
        total_w = picked + tw
        if total_w > 0 && total_w <= global_min_d[]
            lock(res_lock) do
                # Double-check inside the lock to prevent thread races
                if total_w < global_min_d[]
                    empty!(results) # Found a new global minimum! Discard heavier words.
                    global_min_d[] = total_w
                    push!(results, (msg_so_far, curr_tail))
                elseif total_w == global_min_d[]
                    push!(results, (msg_so_far, curr_tail))
                end
            end
        end
        return
    end

    # 4. HARDWARE-ACCELERATED BRANCHING
    if depth < spawn_depth
        tasks = Task[]
        # Branch 1, 2, 3: Include bit with non-zero scalar (1, 2, or 3)
        for alpha in 1:3
            t = Threads.@spawn _words_of_minimum_weight_quaternary_recursive!(
                A_scaled, depth + 1, picked + 1, curr_tail ⊻ A_scaled[alpha, depth + 1], 
                results, lbt, msg_so_far | (UInt128(alpha) << (2 * depth)), spawn_depth, auts, res_lock, global_min_d, k
            )
            push!(tasks, t)
        end
        # Branch 4: Skip bit (scalar = 0)
        _words_of_minimum_weight_quaternary_recursive!(
            A_scaled, depth + 1, picked, curr_tail, 
            results, lbt, msg_so_far, spawn_depth, auts, res_lock, global_min_d, k
        )
        for t in tasks wait(t) end
    else
        # Serial processing
        for alpha in 1:3
            _words_of_minimum_weight_quaternary_recursive!(
                A_scaled, depth + 1, picked + 1, curr_tail ⊻ A_scaled[alpha, depth + 1], 
                results, lbt, msg_so_far | (UInt128(alpha) << (2 * depth)), spawn_depth, auts, res_lock, global_min_d, k
            )
        end
        _words_of_minimum_weight_quaternary_recursive!(
            A_scaled, depth + 1, picked, curr_tail, 
            results, lbt, msg_so_far, spawn_depth, auts, res_lock, global_min_d, k
        )
    end
end

function _words_of_minimum_weight_quaternary(C::AbstractLinearCode; expand::Bool=true, verbose::Bool=false)
    k, n = C.k, C.n
    @assert n - k <= 64 "Parity tail exceeds 64 elements (128 bits); requires chunked bit-packing."
    
    F = C.F
    elements = collect(F)
    elem_to_u8 = Dict(elements[1] => 0, elements[2] => 1, elements[3] => 2, elements[4] => 3)
    GF4_MULT = UInt8[0 0 0 0; 0 1 2 3; 0 2 3 1; 0 3 1 2]
    
    # 1. Standardize Matrix
    G_sys = Array(generator_matrix(C)) 
    A_raw = G_sys[:, (k + 1):n]
    
    # 2. Pre-scale and bit-pack the A columns for alphas 1, 2, 3
    A_scaled = zeros(UInt128, 3, k)
    for row in 1:k
        for alpha in 1:3
            val = UInt128(0)
            for col in 1:(n - k)
                matrix_val = elem_to_u8[A_raw[row, col]]
                prod = GF4_MULT[alpha + 1, matrix_val + 1]
                if prod != 0
                    val |= (UInt128(prod) << (2 * (col - 1)))
                end
            end
            A_scaled[alpha, row] = val
        end
    end
    
    lbt = zeros(Int, k + 1) # Will be filled by your Brouwer LBT pre-computation
    auts = Vector{Vector{Int}}()
    
    # 3. Setup thread-safe dynamic threshold
    global_min_d = Threads.Atomic{Int}(n)
    raw_results = Vector{Tuple{UInt128, UInt128}}()
    res_lock = ReentrantLock()
    
    verbose && println("Executing dynamic GF(4) minimum words search...")
    _words_of_minimum_weight_quaternary_recursive!(
        A_scaled, 0, 0, UInt128(0), 
        raw_results, lbt, UInt128(0), 4, auts, res_lock, global_min_d, k
    )

    final_d = global_min_d[]
    verbose && println("Search complete. Minimum weight: $final_d. Found $(length(raw_results)) representative(s).")

    # 4. Unpack to standard Oscar Vectors
    final_words = Vector{Vector{typeof(elements[1])}}()
    for (msg, tail) in raw_results
        word = fill(elements[1], n)
        for i in 1:k
            val = (msg >> (2 * (i - 1))) & 3
            word[i] = elements[val + 1]
        end
        for i in 1:(n - k)
            val = (tail >> (2 * (i - 1))) & 3
            word[k + i] = elements[val + 1]
        end
        push!(final_words, word)
    end
    
    if expand && (!isempty(auts) || length(elements) > 2)
        # Note: GF(4) requires a scalar multiplication map here!
        return final_d, _expand_orbits(final_words, auts, 4, Matrix{Int}(undef,0,0)) 
    end
    
    return final_d, final_words
end

function _words_of_minimum_weight_nonbinary_recursive!(
    A_idx::Matrix{UInt8}, 
    depth::Int, 
    picked::Int, 
    curr_tail_idx::Vector{UInt8}, 
    results::Vector{Vector{UInt8}}, 
    lbt::Vector{Int}, 
    msg_so_far_idx::Vector{UInt8}, 
    spawn_depth::Int, 
    q::Int, 
    add_table::Matrix{UInt8}, 
    mul_table::Matrix{UInt8}, 
    auts::Vector{Vector{Int}}, 
    res_lock::ReentrantLock,
    global_min_d::Threads.Atomic{Int}
)
    k = size(A_idx, 2)
    n_tail = size(A_idx, 1)

    # 1. DYNAMIC WEIGHT PRUNING
    tw = 0
    @inbounds for i in 1:n_tail
        if curr_tail_idx[i] != 1 tw += 1 end
    end
    
    if (picked + tw + lbt[k - depth + 1]) > global_min_d[]
        return 
    end

    if picked > 1 && !isempty(auts)
        if !_is_canonical_q(msg_so_far_idx, auts, q, mul_table) return end
    end

    if depth == k
        total_w = picked + tw
        if total_w > 0 && total_w <= global_min_d[]
            full_word = vcat(msg_so_far_idx, curr_tail_idx)
            lock(res_lock) do
                if total_w < global_min_d[]
                    empty!(results)
                    global_min_d[] = total_w
                    push!(results, full_word) 
                elseif total_w == global_min_d[]
                    push!(results, full_word)
                end
            end
        end
        return
    end

    # 4. BRANCHING WITH STRICT BACKTRACKING
    for alpha_idx in 2:q
        msg_so_far_idx[depth + 1] = alpha_idx
        old_tail = copy(curr_tail_idx) 
        
        @inbounds for i in 1:n_tail
            prod_idx = mul_table[alpha_idx, A_idx[i, depth + 1]]
            curr_tail_idx[i] = add_table[curr_tail_idx[i], prod_idx]
        end

        if depth < spawn_depth
            t = Threads.@spawn _words_of_minimum_weight_nonbinary_recursive!(
                A_idx, depth + 1, picked + 1, copy(curr_tail_idx), 
                results, lbt, copy(msg_so_far_idx), spawn_depth, q, add_table, mul_table, auts, res_lock, global_min_d
            )
            wait(t) 
        else
            _words_of_minimum_weight_nonbinary_recursive!(
                A_idx, depth + 1, picked + 1, curr_tail_idx, 
                results, lbt, msg_so_far_idx, spawn_depth, q, add_table, mul_table, auts, res_lock, global_min_d
            )
        end
        
        @inbounds for i in 1:n_tail curr_tail_idx[i] = old_tail[i] end
    end

    msg_so_far_idx[depth + 1] = 1
    _words_of_minimum_weight_nonbinary_recursive!(
        A_idx, depth + 1, picked, curr_tail_idx, 
        results, lbt, msg_so_far_idx, spawn_depth, q, add_table, mul_table, auts, res_lock, global_min_d
    )
end

"""
    words_of_minimum_weight(C::AbstractLinearCode; expand::Bool=true, verbose::Bool=false)

Calculates the minimum weight and returns all words of that minimum weight in a single pass.
Automatically routes to zero-allocation hardware engines based on the field size.
"""
function words_of_minimum_weight(C::AbstractLinearCode; expand::Bool=true, verbose::Bool=false)
    # 1. Check the cache!
    # (Assuming you cache the words too. If not, at least check if C.d is known and maybe skip the search if you only want the words)

    k, n = C.k, C.n
    q = length(collect(C.F))
    
    # 2. Route to the hardware engines
    if q == 2
        min_d, raw_words = _words_of_minimum_weight_binary(C; expand=expand, verbose=verbose)
    elseif q == 3 && (n - k) <= 128
        min_d, raw_words = _words_of_minimum_weight_ternary(C; expand=expand, verbose=verbose)
    elseif q == 4 && (n - k) <= 64
        min_d, raw_words = _words_of_minimum_weight_quaternary(C; expand=expand, verbose=verbose)
    else
        verbose && println("Routing to Generic Non-Binary minimum words engine...")
        min_d, raw_words = _words_of_minimum_weight_nonbinary(C; expand=expand, verbose=verbose)
    end
    
    # 3. LOCK IN THE STRUCT PROPERTIES
    C.d = min_d
    C.l_bound = min_d
    C.u_bound = min_d
    
    # Optional: C.minimum_words = raw_words
    
    return min_d, raw_words
end

"""
Adds two GF(3) vectors packed into Low/High registers in parallel.
"""
@inline function _add_gf3_packed(La::UInt128, Ha::UInt128, Lb::UInt128, Hb::UInt128)
    L_out = (La & ~Lb & ~Hb) | (~La & ~Ha & Lb) | (Ha & Hb)
    H_out = (Ha & ~Hb & ~Lb) | (~Ha & ~La & Hb) | (La & Lb)
    return L_out, H_out
end

function _words_of_minimum_weight_ternary_recursive!(
    A_L::Vector{UInt128}, 
    A_H::Vector{UInt128}, 
    depth::Int, 
    picked::Int, 
    curr_tail_L::UInt128, 
    curr_tail_H::UInt128, 
    results::Vector{NTuple{4, UInt128}}, 
    lbt::Vector{Int}, 
    msg_L::UInt128, 
    msg_H::UInt128, 
    spawn_depth::Int, 
    auts::Vector{Vector{Int}}, 
    res_lock::ReentrantLock,
    global_min_d::Threads.Atomic{Int},
    k::Int
)
    # 1. DYNAMIC WEIGHT PRUNING
    # A position is non-zero if it has a 1 OR a 2
    tw = count_ones(curr_tail_L | curr_tail_H)
    
    if (picked + tw + lbt[k - depth + 1]) > global_min_d[]
        return 
    end

    # 2. SYMMETRY PRUNING
    # (Requires a split-register overloaded _is_canonical_q)
    if picked > 1 && !isempty(auts)
        if !_is_canonical_q(msg_L, msg_H, auts, depth) return end
    end

    # 3. BASE CASE
    if depth == k
        total_w = picked + tw
        if total_w > 0 && total_w <= global_min_d[]
            lock(res_lock) do
                if total_w < global_min_d[]
                    empty!(results)
                    global_min_d[] = total_w
                    push!(results, (msg_L, msg_H, curr_tail_L, curr_tail_H))
                elseif total_w == global_min_d[]
                    push!(results, (msg_L, msg_H, curr_tail_L, curr_tail_H))
                end
            end
        end
        return
    end

    # 4. HARDWARE-ACCELERATED BRANCHING
    row_L = A_L[depth + 1]
    row_H = A_H[depth + 1]

    if depth < spawn_depth
        tasks = Task[]
        
        # Branch 1: Alpha = 1
        nL1, nH1 = _add_gf3_packed(curr_tail_L, curr_tail_H, row_L, row_H)
        t1 = Threads.@spawn _words_of_minimum_weight_ternary_recursive!(
            A_L, A_H, depth + 1, picked + 1, nL1, nH1, 
            results, lbt, msg_L | (UInt128(1) << depth), msg_H, 
            spawn_depth, auts, res_lock, global_min_d, k
        )
        push!(tasks, t1)

        # Branch 2: Alpha = 2 (Scalar multiplication by 2 swaps the L and H registers!)
        nL2, nH2 = _add_gf3_packed(curr_tail_L, curr_tail_H, row_H, row_L)
        t2 = Threads.@spawn _words_of_minimum_weight_ternary_recursive!(
            A_L, A_H, depth + 1, picked + 1, nL2, nH2, 
            results, lbt, msg_L, msg_H | (UInt128(1) << depth), 
            spawn_depth, auts, res_lock, global_min_d, k
        )
        push!(tasks, t2)

        # Branch 3: Alpha = 0 (Skip bit)
        _words_of_minimum_weight_ternary_recursive!(
            A_L, A_H, depth + 1, picked, curr_tail_L, curr_tail_H, 
            results, lbt, msg_L, msg_H, 
            spawn_depth, auts, res_lock, global_min_d, k
        )
        
        for t in tasks wait(t) end
    else
        # Serial processing
        nL1, nH1 = _add_gf3_packed(curr_tail_L, curr_tail_H, row_L, row_H)
        _words_of_minimum_weight_ternary_recursive!(
            A_L, A_H, depth + 1, picked + 1, nL1, nH1, 
            results, lbt, msg_L | (UInt128(1) << depth), msg_H, 
            spawn_depth, auts, res_lock, global_min_d, k
        )

        nL2, nH2 = _add_gf3_packed(curr_tail_L, curr_tail_H, row_H, row_L)
        _words_of_minimum_weight_ternary_recursive!(
            A_L, A_H, depth + 1, picked + 1, nL2, nH2, 
            results, lbt, msg_L, msg_H | (UInt128(1) << depth), 
            spawn_depth, auts, res_lock, global_min_d, k
        )

        _words_of_minimum_weight_ternary_recursive!(
            A_L, A_H, depth + 1, picked, curr_tail_L, curr_tail_H, 
            results, lbt, msg_L, msg_H, 
            spawn_depth, auts, res_lock, global_min_d, k
        )
    end
end

function _words_of_minimum_weight_ternary(C::AbstractLinearCode; expand::Bool=true, verbose::Bool=false)
    k, n = C.k, C.n
    @assert n - k <= 128 "Parity tail exceeds 128 bits; requires chunked bit-packing."
    
    F = C.F
    elements = collect(F)
    elem_to_idx = Dict(elements[1] => 0, elements[2] => 1, elements[3] => 2)
    
    # 1. Standardize Matrix
    G_sys = Array(generator_matrix(C)) 
    A_raw = G_sys[:, (k + 1):n]
    
    # 2. Pack the A columns into split L/H registers
    A_L = zeros(UInt128, k)
    A_H = zeros(UInt128, k)
    
    for row in 1:k
        val_L = UInt128(0)
        val_H = UInt128(0)
        for col in 1:(n - k)
            matrix_val = elem_to_idx[A_raw[row, col]]
            if matrix_val == 1
                val_L |= (UInt128(1) << (col - 1))
            elseif matrix_val == 2
                val_H |= (UInt128(1) << (col - 1))
            end
        end
        A_L[row] = val_L
        A_H[row] = val_H
    end
    
    lbt = zeros(Int, k + 1) # Will be filled by your Brouwer LBT logic
    auts = Vector{Vector{Int}}()
    
    # 3. Setup thread-safe dynamic threshold
    global_min_d = Threads.Atomic{Int}(n)
    raw_results = Vector{NTuple{4, UInt128}}()
    res_lock = ReentrantLock()
    
    verbose && println("Executing dynamic GF(3) minimum words search...")
    _words_of_minimum_weight_ternary_recursive!(
        A_L, A_H, 0, 0, UInt128(0), UInt128(0), 
        raw_results, lbt, UInt128(0), UInt128(0), 
        4, auts, res_lock, global_min_d, k
    )

    final_d = global_min_d[]
    verbose && println("Search complete. Minimum weight: $final_d. Found $(length(raw_results)) representative(s).")

    # 4. Unpack to standard Oscar Vectors
    final_words = Vector{Vector{typeof(elements[1])}}()
    for (mL, mH, tL, tH) in raw_results
        word = fill(elements[1], n)
        
        # Unpack message
        for i in 1:k
            if ((mL >> (i - 1)) & 1) == 1
                word[i] = elements[2]
            elseif ((mH >> (i - 1)) & 1) == 1
                word[i] = elements[3]
            end
        end
        
        # Unpack parity tail
        for i in 1:(n - k)
            if ((tL >> (i - 1)) & 1) == 1
                word[k + i] = elements[2]
            elseif ((tH >> (i - 1)) & 1) == 1
                word[k + i] = elements[3]
            end
        end
        push!(final_words, word)
    end
    
    if expand && (!isempty(auts) || length(elements) > 2)
        # Note: GF(3) requires a scalar multiplication map here!
        return final_d, _expand_orbits(final_words, auts, 3, Matrix{Int}(undef,0,0)) 
    end
    
    return final_d, final_words
end

function _minimal_spanning_set_binary_recursive!(
    A_rows::Vector{UInt128}, 
    depth::Int, 
    picked::Int, 
    curr_tail::UInt128, 
    lbt::Vector{Int}, 
    msg_so_far::UInt128, 
    spawn_depth::Int, 
    auts::Vector{Vector{Int}}, 
    res_lock::ReentrantLock,
    global_B::Threads.Atomic{Int},
    candidate_pool::Vector{Tuple{Int, UInt128, UInt128}}, # (weight, msg, tail)
    k::Int
)
    # 1. DYNAMIC WEIGHT PRUNING
    tw = count_ones(curr_tail)
    
    # Prune if the absolute minimum weight exceeds our spanning threshold B
    if (picked + tw + lbt[k - depth + 1]) > global_B[]
        return 
    end

    # 2. SYMMETRY PRUNING
    if picked > 1 && !isempty(auts)
        if !_is_canonical(msg_so_far, auts, depth) return end
    end

    # 3. BASE CASE: Full codeword reached
    if depth == k
        total_w = picked + tw
        if total_w > 0 && total_w <= global_B[]
            lock(res_lock) do
                # Double check inside the lock
                if total_w <= global_B[]
                    push!(candidate_pool, (total_w, msg_so_far, curr_tail))
                    
                    # Greedily Recompute the Basis
                    sort!(candidate_pool, by = x -> x[1])
                    basis = UInt128[]
                    max_w_in_basis = 0
                    
                    for (w, m, _) in candidate_pool
                        v_work = m
                        # Fast O(k) linear independence check using only the message register
                        for b in basis
                            msb_idx = 127 - leading_zeros(b)
                            if ((v_work >> msb_idx) & 1) == 1
                                v_work ⊻= b
                            end
                        end
                        
                        if v_work != 0
                            push!(basis, v_work)
                            # Keep basis sorted so the highest MSB is first
                            sort!(basis, by = leading_zeros) 
                            max_w_in_basis = max(max_w_in_basis, w)
                        end
                    end
                    
                    # If we have a full span, shrink the global upper bound B
                    if length(basis) == k
                        if max_w_in_basis < global_B[]
                            global_B[] = max_w_in_basis
                            # Instantly prune the pool of strictly heavier words
                            filter!(x -> x[1] <= global_B[], candidate_pool)
                        end
                    end
                end
            end
        end
        return
    end

    # 4. HARDWARE-ACCELERATED BRANCHING
    if depth < spawn_depth
        t = Threads.@spawn _minimal_spanning_set_binary_recursive!(
            A_rows, depth + 1, picked + 1, curr_tail ⊻ A_rows[depth + 1], 
            lbt, msg_so_far | (UInt128(1) << depth), spawn_depth, auts, res_lock, global_B, candidate_pool, k
        )
        
        _minimal_spanning_set_binary_recursive!(
            A_rows, depth + 1, picked, curr_tail, 
            lbt, msg_so_far, spawn_depth, auts, res_lock, global_B, candidate_pool, k
        )
        wait(t)
    else
        # Serial processing
        _minimal_spanning_set_binary_recursive!(
            A_rows, depth + 1, picked + 1, curr_tail ⊻ A_rows[depth + 1], 
            lbt, msg_so_far | (UInt128(1) << depth), spawn_depth, auts, res_lock, global_B, candidate_pool, k
        )
        _minimal_spanning_set_binary_recursive!(
            A_rows, depth + 1, picked, curr_tail, 
            lbt, msg_so_far, spawn_depth, auts, res_lock, global_B, candidate_pool, k
        )
    end
end

function _minimal_spanning_set_binary(C::AbstractLinearCode; expand::Bool=true, verbose::Bool=false)
    k, n = C.k, C.n
    @assert n - k <= 128 "Parity tail exceeds 128 bits; requires chunked bit-packing."
    
    # 1. Standardize Matrix & Prepare LBT
    G_sys = Array(generator_matrix(C))
    A_rows, lbt = _prepare_packed_brouwer_binary(G_sys)
    auts = Vector{Vector{Int}}() 
    
    # 2. Setup thread-safe dynamic spanning threshold B
    global_B = Threads.Atomic{Int}(n)
    candidate_pool = Vector{Tuple{Int, UInt128, UInt128}}()
    res_lock = ReentrantLock()
    
    verbose && println("Executing dynamic minimal spanning set search...")
    _minimal_spanning_set_binary_recursive!(
        A_rows, 0, 0, UInt128(0), 
        lbt, UInt128(0), 5, auts, res_lock, global_B, candidate_pool, k
    )

    B_final = global_B[]
    min_d = isempty(candidate_pool) ? -1 : minimum(x -> x[1], candidate_pool)
    
    verbose && println("Search complete. Minimum weight: $min_d. Spanning threshold (B): $B_final.")
    verbose && println("Found $(length(candidate_pool)) representative(s) spanning the code.")

    # 3. Unpack to standard Vectors
    F = C.F
    T_zero, T_one = zero(F), one(F)
    final_words = Vector{Vector{typeof(T_zero)}}()
    
    for (_, msg, tail) in candidate_pool
        word = fill(T_zero, n)
        for i in 1:k
            if ((msg >> (i - 1)) & 1) == 1 word[i] = T_one end
        end
        for i in 1:(n - k)
            if ((tail >> (i - 1)) & 1) == 1 word[k + i] = T_one end
        end
        push!(final_words, word)
    end
    
    if expand
        return min_d, B_final, _expand_orbits(final_words, auts, 2)
    end
    
    return min_d, B_final, final_words
end

function _minimal_spanning_set_quaternary_recursive!(
    A_scaled::Matrix{UInt128}, 
    depth::Int, 
    picked::Int, 
    curr_tail::UInt128, 
    lbt::Vector{Int}, 
    msg_so_far::UInt128, 
    spawn_depth::Int, 
    auts::Vector{Vector{Int}}, 
    res_lock::ReentrantLock,
    global_B::Threads.Atomic{Int},
    candidate_pool::Vector{Tuple{Int, UInt128, UInt128}},
    k::Int,
    F, elements, elem_to_u8
)
    tw = _weight_gf4(curr_tail)
    
    if (picked + tw + lbt[k - depth + 1]) > global_B[]
        return 
    end

    if picked > 1 && !isempty(auts)
        if !_is_canonical_q(msg_so_far, auts, depth) return end
    end

    if depth == k
        total_w = picked + tw
        if total_w > 0 && total_w <= global_B[]
            lock(res_lock) do
                if total_w <= global_B[]
                    push!(candidate_pool, (total_w, msg_so_far, curr_tail))
                    
                    # Unpack messages to check true GF(4) rank
                    mat_rows = []
                    for (_, m, _) in candidate_pool
                        vec = fill(elements[1], k)
                        for i in 1:k
                            val = (m >> (2 * (i - 1))) & 3
                            vec[i] = elements[val + 1]
                        end
                        push!(mat_rows, vec)
                    end
                    
                    # If we hit full rank, shrink B and prune
                    if length(mat_rows) >= k && rank(matrix(F, mat_rows)) == k
                        max_w = maximum(x -> x[1], candidate_pool)
                        if max_w < global_B[]
                            global_B[] = max_w
                            filter!(x -> x[1] <= global_B[], candidate_pool)
                        end
                    end
                end
            end
        end
        return
    end

    if depth < spawn_depth
        tasks = Task[]
        for alpha in 1:3
            t = Threads.@spawn _minimal_spanning_set_quaternary_recursive!(
                A_scaled, depth + 1, picked + 1, curr_tail ⊻ A_scaled[alpha, depth + 1], 
                lbt, msg_so_far | (UInt128(alpha) << (2 * depth)), spawn_depth, auts, res_lock, global_B, candidate_pool, k, F, elements, elem_to_u8
            )
            push!(tasks, t)
        end
        _minimal_spanning_set_quaternary_recursive!(
            A_scaled, depth + 1, picked, curr_tail, 
            lbt, msg_so_far, spawn_depth, auts, res_lock, global_B, candidate_pool, k, F, elements, elem_to_u8
        )
        for t in tasks wait(t) end
    else
        for alpha in 1:3
            _minimal_spanning_set_quaternary_recursive!(
                A_scaled, depth + 1, picked + 1, curr_tail ⊻ A_scaled[alpha, depth + 1], 
                lbt, msg_so_far | (UInt128(alpha) << (2 * depth)), spawn_depth, auts, res_lock, global_B, candidate_pool, k, F, elements, elem_to_u8
            )
        end
        _minimal_spanning_set_quaternary_recursive!(
            A_scaled, depth + 1, picked, curr_tail, 
            lbt, msg_so_far, spawn_depth, auts, res_lock, global_B, candidate_pool, k, F, elements, elem_to_u8
        )
    end
end

function _minimal_spanning_set_quaternary(C::AbstractLinearCode; expand::Bool=true, verbose::Bool=false)
    k, n = C.k, C.n
    F = C.F
    elements = collect(F)
    elem_to_u8 = Dict(elements[1] => 0, elements[2] => 1, elements[3] => 2, elements[4] => 3)
    GF4_MULT = UInt8[0 0 0 0; 0 1 2 3; 0 2 3 1; 0 3 1 2]
    
    G_sys = Array(generator_matrix(C)) 
    A_raw = G_sys[:, (k + 1):n]
    
    A_scaled = zeros(UInt128, 3, k)
    for row in 1:k
        for alpha in 1:3
            val = UInt128(0)
            for col in 1:(n - k)
                matrix_val = elem_to_u8[A_raw[row, col]]
                prod = GF4_MULT[alpha + 1, matrix_val + 1]
                if prod != 0 val |= (UInt128(prod) << (2 * (col - 1))) end
            end
            A_scaled[alpha, row] = val
        end
    end
    
    lbt = zeros(Int, k + 1)
    auts = Vector{Vector{Int}}()
    
    global_B = Threads.Atomic{Int}(n)
    candidate_pool = Vector{Tuple{Int, UInt128, UInt128}}()
    res_lock = ReentrantLock()
    
    _minimal_spanning_set_quaternary_recursive!(
        A_scaled, 0, 0, UInt128(0), lbt, UInt128(0), 4, auts, res_lock, global_B, candidate_pool, k, F, elements, elem_to_u8
    )

    B_final = global_B[]
    min_d = isempty(candidate_pool) ? -1 : minimum(x -> x[1], candidate_pool)
    
    final_words = Vector{Vector{typeof(elements[1])}}()
    for (_, msg, tail) in candidate_pool
        word = fill(elements[1], n)
        for i in 1:k
            val = (msg >> (2 * (i - 1))) & 3
            word[i] = elements[val + 1]
        end
        for i in 1:(n - k)
            val = (tail >> (2 * (i - 1))) & 3
            word[k + i] = elements[val + 1]
        end
        push!(final_words, word)
    end
    
    if expand return min_d, B_final, _expand_orbits(final_words, auts, 4) end
    return min_d, B_final, final_words
end

function _minimal_spanning_set_ternary_recursive!(
    A_L::Vector{UInt128}, A_H::Vector{UInt128}, 
    depth::Int, picked::Int, 
    curr_tail_L::UInt128, curr_tail_H::UInt128, 
    lbt::Vector{Int}, 
    msg_L::UInt128, msg_H::UInt128, 
    spawn_depth::Int, auts::Vector{Vector{Int}}, 
    res_lock::ReentrantLock, global_B::Threads.Atomic{Int},
    candidate_pool::Vector{Tuple{Int, UInt128, UInt128, UInt128, UInt128}},
    k::Int, F, elements
)
    tw = count_ones(curr_tail_L | curr_tail_H)
    
    if (picked + tw + lbt[k - depth + 1]) > global_B[]
        return 
    end

    if picked > 1 && !isempty(auts)
        if !_is_canonical_q(msg_L, msg_H, auts, depth) return end
    end

    if depth == k
        total_w = picked + tw
        if total_w > 0 && total_w <= global_B[]
            lock(res_lock) do
                if total_w <= global_B[]
                    push!(candidate_pool, (total_w, msg_L, msg_H, curr_tail_L, curr_tail_H))
                    
                    mat_rows = []
                    for (_, mL, mH, _, _) in candidate_pool
                        vec = fill(elements[1], k)
                        for i in 1:k
                            if ((mL >> (i - 1)) & 1) == 1
                                vec[i] = elements[2]
                            elseif ((mH >> (i - 1)) & 1) == 1
                                vec[i] = elements[3]
                            end
                        end
                        push!(mat_rows, vec)
                    end
                    
                    if length(mat_rows) >= k && rank(matrix(F, mat_rows)) == k
                        max_w = maximum(x -> x[1], candidate_pool)
                        if max_w < global_B[]
                            global_B[] = max_w
                            filter!(x -> x[1] <= global_B[], candidate_pool)
                        end
                    end
                end
            end
        end
        return
    end

    row_L = A_L[depth + 1]
    row_H = A_H[depth + 1]

    if depth < spawn_depth
        tasks = Task[]
        
        nL1, nH1 = _add_gf3_packed(curr_tail_L, curr_tail_H, row_L, row_H)
        t1 = Threads.@spawn _minimal_spanning_set_ternary_recursive!(
            A_L, A_H, depth + 1, picked + 1, nL1, nH1, 
            lbt, msg_L | (UInt128(1) << depth), msg_H, spawn_depth, auts, res_lock, global_B, candidate_pool, k, F, elements
        )
        push!(tasks, t1)

        nL2, nH2 = _add_gf3_packed(curr_tail_L, curr_tail_H, row_H, row_L)
        t2 = Threads.@spawn _minimal_spanning_set_ternary_recursive!(
            A_L, A_H, depth + 1, picked + 1, nL2, nH2, 
            lbt, msg_L, msg_H | (UInt128(1) << depth), spawn_depth, auts, res_lock, global_B, candidate_pool, k, F, elements
        )
        push!(tasks, t2)

        _minimal_spanning_set_ternary_recursive!(
            A_L, A_H, depth + 1, picked, curr_tail_L, curr_tail_H, 
            lbt, msg_L, msg_H, spawn_depth, auts, res_lock, global_B, candidate_pool, k, F, elements
        )
        for t in tasks wait(t) end
    else
        nL1, nH1 = _add_gf3_packed(curr_tail_L, curr_tail_H, row_L, row_H)
        _minimal_spanning_set_ternary_recursive!(
            A_L, A_H, depth + 1, picked + 1, nL1, nH1, 
            lbt, msg_L | (UInt128(1) << depth), msg_H, spawn_depth, auts, res_lock, global_B, candidate_pool, k, F, elements
        )

        nL2, nH2 = _add_gf3_packed(curr_tail_L, curr_tail_H, row_H, row_L)
        _minimal_spanning_set_ternary_recursive!(
            A_L, A_H, depth + 1, picked + 1, nL2, nH2, 
            lbt, msg_L, msg_H | (UInt128(1) << depth), spawn_depth, auts, res_lock, global_B, candidate_pool, k, F, elements
        )

        _minimal_spanning_set_ternary_recursive!(
            A_L, A_H, depth + 1, picked, curr_tail_L, curr_tail_H, 
            lbt, msg_L, msg_H, spawn_depth, auts, res_lock, global_B, candidate_pool, k, F, elements
        )
    end
end

function _minimal_spanning_set_ternary(C::AbstractLinearCode; expand::Bool=true, verbose::Bool=false)
    k, n = C.k, C.n
    F = C.F
    elements = collect(F)
    elem_to_idx = Dict(elements[1] => 0, elements[2] => 1, elements[3] => 2)
    
    G_sys = Array(generator_matrix(C)) 
    A_raw = G_sys[:, (k + 1):n]
    
    A_L = zeros(UInt128, k)
    A_H = zeros(UInt128, k)
    for row in 1:k
        val_L = UInt128(0)
        val_H = UInt128(0)
        for col in 1:(n - k)
            matrix_val = elem_to_idx[A_raw[row, col]]
            if matrix_val == 1 val_L |= (UInt128(1) << (col - 1))
            elseif matrix_val == 2 val_H |= (UInt128(1) << (col - 1)) end
        end
        A_L[row] = val_L
        A_H[row] = val_H
    end
    
    lbt = zeros(Int, k + 1)
    auts = Vector{Vector{Int}}()
    
    global_B = Threads.Atomic{Int}(n)
    candidate_pool = Vector{Tuple{Int, UInt128, UInt128, UInt128, UInt128}}()
    res_lock = ReentrantLock()
    
    _minimal_spanning_set_ternary_recursive!(
        A_L, A_H, 0, 0, UInt128(0), UInt128(0), lbt, UInt128(0), UInt128(0), 4, auts, res_lock, global_B, candidate_pool, k, F, elements
    )

    B_final = global_B[]
    min_d = isempty(candidate_pool) ? -1 : minimum(x -> x[1], candidate_pool)
    
    final_words = Vector{Vector{typeof(elements[1])}}()
    for (_, mL, mH, tL, tH) in candidate_pool
        word = fill(elements[1], n)
        for i in 1:k
            if ((mL >> (i - 1)) & 1) == 1 word[i] = elements[2]
            elseif ((mH >> (i - 1)) & 1) == 1 word[i] = elements[3] end
        end
        for i in 1:(n - k)
            if ((tL >> (i - 1)) & 1) == 1 word[k + i] = elements[2]
            elseif ((tH >> (i - 1)) & 1) == 1 word[k + i] = elements[3] end
        end
        push!(final_words, word)
    end
    
    if expand return min_d, B_final, _expand_orbits(final_words, auts, 3) end
    return min_d, B_final, final_words
end

function _words_of_weight_range_binary_recursive!(
    A_rows::Vector{UInt128}, 
    depth::Int, 
    picked::Int, 
    curr_tail::UInt128, 
    results::Vector{Tuple{UInt128, UInt128}}, 
    lbt::Vector{Int}, 
    msg_so_far::UInt128, 
    spawn_depth::Int, 
    auts::Vector{Vector{Int}}, 
    w_min::Int, 
    w_max::Int, 
    k::Int, 
    res_lock::ReentrantLock
)
    # 1. STATIC WEIGHT PRUNING
    tw = count_ones(curr_tail)
    
    # Prune if the absolute minimum weight this branch can produce exceeds w_max
    if (picked + tw + lbt[k - depth + 1]) > w_max
        return 
    end

    # 2. SYMMETRY PRUNING
    # (Requires the UInt128 overloaded _is_canonical we built earlier)
    if picked > 1 && !isempty(auts)
        if !_is_canonical(msg_so_far, auts, depth) return end
    end

    # 3. BASE CASE: Full codeword reached
    if depth == k
        total_w = picked + tw
        if w_min <= total_w <= w_max
            lock(res_lock) do
                push!(results, (msg_so_far, curr_tail))
            end
        end
        return
    end

    # 4. HARDWARE-ACCELERATED BRANCHING
    if depth < spawn_depth
        # Branch 1: Include bit
        t = Threads.@spawn _words_of_weight_range_binary_recursive!(
            A_rows, depth + 1, picked + 1, curr_tail ⊻ A_rows[depth + 1], 
            results, lbt, msg_so_far | (UInt128(1) << depth), spawn_depth, auts, w_min, w_max, k, res_lock
        )
        
        # Branch 2: Skip bit
        _words_of_weight_range_binary_recursive!(
            A_rows, depth + 1, picked, curr_tail, 
            results, lbt, msg_so_far, spawn_depth, auts, w_min, w_max, k, res_lock
        )
        wait(t)
    else
        # Serial processing (pure register math)
        _words_of_weight_range_binary_recursive!(
            A_rows, depth + 1, picked + 1, curr_tail ⊻ A_rows[depth + 1], 
            results, lbt, msg_so_far | (UInt128(1) << depth), spawn_depth, auts, w_min, w_max, k, res_lock
        )
        
        _words_of_weight_range_binary_recursive!(
            A_rows, depth + 1, picked, curr_tail, 
            results, lbt, msg_so_far, spawn_depth, auts, w_min, w_max, k, res_lock
        )
    end
end

function _words_of_weight_range_trellis(G::Matrix{T}, boundaries::Vector{Int}, w_min::Int, w_max::Int, q::Int, verbose::Bool) where T
    k, n = size(G)
    F = parent(G[1, 1])
    elements = collect(F)
    T_zero = zero(F)
    
    L, R = _get_LR_indices(G)
    start_sets, working_sets, active_sets, _ = _build_trellis_sets(L, R, boundaries)

    num_sections = length(boundaries) - 1
    
    # Forward Pass: Store all reachable weights and the incoming branches
    # Map: State => Vector of (Incoming_State, Branch_Scalars, Accumulated_Weight)
    layers = Vector{Dict{Vector{T}, Vector{Tuple{Vector{T}, Vector{T}, Int}}}}()
    
    init_dict = Dict{Vector{T}, Vector{Tuple{Vector{T}, Vector{T}, Int}}}()
    init_dict[Vector{T}()] = [(Vector{T}(), Vector{T}(), 0)]
    push!(layers, init_dict)
    
    u_buf = fill(T_zero, k)
    
    p = verbose ? Progress(num_sections, 0.1, "Trellis Forward Sweep: ") : nothing
    
    for m in 1:num_sections
        next_layer = Dict{Vector{T}, Vector{Tuple{Vector{T}, Vector{T}, Int}}}()
        
        active_curr = active_sets[m]
        left_b, right_b = boundaries[m], boundaries[m+1]
        active_prev = m == 1 ? Int[] : active_sets[m-1]
        starting = start_sets[m]
        working = working_sets[m]
        
        for (prev_scalars, paths) in layers[end]
            # Get the unique accumulated weights reaching this prev state to avoid duplicate branch exploration
            reaching_wts = unique([p[3] for p in paths])
            
            for branch_scalars in Iterators.product(fill(elements, length(starting))...)
                fill!(u_buf, T_zero)
                @inbounds for (idx, row) in enumerate(active_prev) u_buf[row] = prev_scalars[idx] end
                @inbounds for (idx, row) in enumerate(starting)    u_buf[row] = branch_scalars[idx] end
                
                chunk_wt = 0
                for col in (left_b + 1):right_b
                    c_i = T_zero
                    @inbounds for row in working c_i += u_buf[row] * G[row, col] end
                    if !iszero(c_i) chunk_wt += 1 end
                end
                
                next_scalars = [u_buf[row] for row in active_curr]
                branch_vec = [b for b in branch_scalars]
                
                if !haskey(next_layer, next_scalars)
                    next_layer[next_scalars] = Vector{Tuple{Vector{T}, Vector{T}, Int}}()
                end
                
                # Add valid paths
                dest_list = next_layer[next_scalars]
                for wt in reaching_wts
                    new_wt = wt + chunk_wt
                    if new_wt <= w_max
                        push!(dest_list, (prev_scalars, branch_vec, new_wt))
                    end
                end
            end
        end
        push!(layers, next_layer)
        verbose && next!(p)
    end
    
    # Traceback Pass
    verbose && println("Executing Trellis Traceback...")
    final_state = Vector{T}()
    valid_words = Vector{Vector{T}}()
    
    if !haskey(layers[end], final_state)
        return valid_words
    end
    
    # Recursive traceback to rebuild the full vectors
    function _traceback(section::Int, curr_state::Vector{T}, current_word::Vector{T}, current_wt::Int)
        if section == 1
            if w_min <= current_wt <= w_max
                push!(valid_words, reverse(current_word))
            end
            return
        end
        
        for (prev_state, branch_scalars, prev_wt) in layers[section][curr_state]
            # To reconstruct the exact codeword, we multiply the branch scalars by the generator matrix
            # Since this is just a traceback reconstruction, doing it section by section is trivial
            _traceback(section - 1, prev_state, vcat(current_word, branch_scalars), prev_wt)
        end
    end
    
    # Initiate traceback from the valid final states
    for (prev_state, branch_scalars, prev_wt) in layers[end][final_state]
        _traceback(num_sections, prev_state, branch_scalars, prev_wt)
    end
    
    # Note: The traceback recovers the Information Vectors (message). 
    # We must multiply them by G to get the full codewords.
    return [msg * G for msg in valid_words]
end

"""
    words_of_weight(C::AbstractLinearCode, w_range::UnitRange{Int}; max_span::Int=20, expand::Bool=true, verbose::Bool=false)

Extracts all codewords whose weight falls within `w_range` (e.g., `1:4`).
Dynamically routes to a Trellis Traceback sweep for low-complexity codes, or a bit-packed DFS for dense codes.
"""
function words_of_weight(C::AbstractLinearCode, w_range::UnitRange{Int}; max_span::Int=20, expand::Bool=true, verbose::Bool=false)
    w_min, w_max = first(w_range), last(w_range)
    k, n = C.k, C.n
    q = Int(order(C.F))
    
    # 1. Profile Trellis Complexity
    G_mat = Array(generator_matrix(C))
    best_M, best_perm, peak_E = optimize_trellis_permutation(G_mat, 50)
    
    # 2. Heuristic Routing
    if peak_E <= max_span
        verbose && println("Trellis profile is thin (Peak E: $peak_E <= $max_span). Routing to Trellis Traceback Sweep...")
        bounds = optimal_sectionalization(best_M, q)
        
        # Underscored internal call
        raw_words = _words_of_weight_range_trellis(best_M, bounds, w_min, w_max, q, verbose)
        
        # Invert the column permutation applied during Trellis optimization
        inv_perm = invperm(best_perm)
        final_words = [w[inv_perm] for w in raw_words]
        
    else
        verbose && println("Trellis profile too dense (Peak E: $peak_E > $max_span). Routing to Hardware DFS...")
        if q == 2
            final_words = _words_of_weight_range_dfs_binary(C, w_min, w_max, verbose)
        elseif q == 3 && (n - k) <= 128
            final_words = _words_of_weight_range_dfs_ternary(C, w_min, w_max, verbose)
        elseif q == 4 && (n - k) <= 64
            final_words = _words_of_weight_range_dfs_quaternary(C, w_min, w_max, verbose)
        else
            final_words = _words_of_weight_range_dfs_nonbinary(C, w_min, w_max, verbose)
        end
    end
    
    # 3. Orbit Expansion
    if expand
        auts = generate_automorphisms(C) 
        return _expand_orbits(final_words, auts, q)
    end
    
    return final_words
end

"""
    words_of_weight(C::AbstractLinearCode, w::Int; kwargs...)

Fallback to extract codewords of a single specific weight `w`.
"""
words_of_weight(C::AbstractLinearCode, w::Int; kwargs...) = words_of_weight(C, w:w; kwargs...)

function _words_of_weight_range_dfs_binary(C::AbstractLinearCode, w_min::Int, w_max::Int, verbose::Bool)
    k, n = C.k, C.n
    @assert n - k <= 128 "Parity tail exceeds 128 bits; requires chunked bit-packing."
    
    # 1. Standardize Matrix & Prepare LBT
    G_sys = Array(generator_matrix(C))
    A_rows, lbt = _prepare_packed_brouwer_binary(G_sys)
    auts = Vector{Vector{Int}}() # Assume generate_automorphisms is handled
    
    # 2. Setup Threading
    raw_results = Vector{Tuple{UInt128, UInt128}}()
    res_lock = ReentrantLock()
    
    verbose && println("Executing Hardware DFS for weight range [$w_min, $w_max]...")
    _words_of_weight_range_binary_recursive!(
        A_rows, 0, 0, UInt128(0), raw_results, lbt, UInt128(0), 5, auts, w_min, w_max, k, res_lock
    )

    verbose && println("Search complete. Found $(length(raw_results)) representative(s).")

    # 3. Unpack to standard Vectors
    F = C.F
    T_zero, T_one = zero(F), one(F)
    final_words = Vector{Vector{typeof(T_zero)}}()
    
    for (msg, tail) in raw_results
        word = fill(T_zero, n)
        for i in 1:k
            if ((msg >> (i - 1)) & 1) == 1 word[i] = T_one end
        end
        for i in 1:(n - k)
            if ((tail >> (i - 1)) & 1) == 1 word[k + i] = T_one end
        end
        push!(final_words, word)
    end
    
    return final_words
end

function _words_of_weight_range_quaternary_recursive!(
    A_scaled::Matrix{UInt128}, depth::Int, picked::Int, 
    curr_tail::UInt128, results::Vector{Tuple{UInt128, UInt128}}, 
    lbt::Vector{Int}, msg_so_far::UInt128, spawn_depth::Int, 
    auts::Vector{Vector{Int}}, w_min::Int, w_max::Int, 
    k::Int, res_lock::ReentrantLock
)
    tw = _weight_gf4(curr_tail)
    
    if (picked + tw + lbt[k - depth + 1]) > w_max
        return 
    end

    if picked > 1 && !isempty(auts)
        if !_is_canonical_q(msg_so_far, auts, depth) return end
    end

    if depth == k
        total_w = picked + tw
        if w_min <= total_w <= w_max
            lock(res_lock) do
                push!(results, (msg_so_far, curr_tail))
            end
        end
        return
    end

    if depth < spawn_depth
        tasks = Task[]
        for alpha in 1:3
            t = Threads.@spawn _words_of_weight_range_quaternary_recursive!(
                A_scaled, depth + 1, picked + 1, curr_tail ⊻ A_scaled[alpha, depth + 1], 
                results, lbt, msg_so_far | (UInt128(alpha) << (2 * depth)), spawn_depth, auts, w_min, w_max, k, res_lock
            )
            push!(tasks, t)
        end
        _words_of_weight_range_quaternary_recursive!(
            A_scaled, depth + 1, picked, curr_tail, 
            results, lbt, msg_so_far, spawn_depth, auts, w_min, w_max, k, res_lock
        )
        for t in tasks wait(t) end
    else
        for alpha in 1:3
            _words_of_weight_range_quaternary_recursive!(
                A_scaled, depth + 1, picked + 1, curr_tail ⊻ A_scaled[alpha, depth + 1], 
                results, lbt, msg_so_far | (UInt128(alpha) << (2 * depth)), spawn_depth, auts, w_min, w_max, k, res_lock
            )
        end
        _words_of_weight_range_quaternary_recursive!(
            A_scaled, depth + 1, picked, curr_tail, 
            results, lbt, msg_so_far, spawn_depth, auts, w_min, w_max, k, res_lock
        )
    end
end

function _words_of_weight_range_dfs_quaternary(C::AbstractLinearCode, w_min::Int, w_max::Int, verbose::Bool)
    k, n = C.k, C.n
    F = C.F
    elements = collect(F)
    elem_to_u8 = Dict(elements[1] => 0, elements[2] => 1, elements[3] => 2, elements[4] => 3)
    GF4_MULT = UInt8[0 0 0 0; 0 1 2 3; 0 2 3 1; 0 3 1 2]
    
    G_sys = Array(generator_matrix(C)) 
    A_raw = G_sys[:, (k + 1):n]
    
    A_scaled = zeros(UInt128, 3, k)
    for row in 1:k
        for alpha in 1:3
            val = UInt128(0)
            for col in 1:(n - k)
                matrix_val = elem_to_u8[A_raw[row, col]]
                prod = GF4_MULT[alpha + 1, matrix_val + 1]
                if prod != 0 val |= (UInt128(prod) << (2 * (col - 1))) end
            end
            A_scaled[alpha, row] = val
        end
    end
    
    lbt = zeros(Int, k + 1) # Assumes Brouwer LBT precomp will be patched here
    auts = Vector{Vector{Int}}()
    
    raw_results = Vector{Tuple{UInt128, UInt128}}()
    res_lock = ReentrantLock()
    
    verbose && println("Executing Hardware DFS (GF(4)) for weight range [$w_min, $w_max]...")
    _words_of_weight_range_quaternary_recursive!(
        A_scaled, 0, 0, UInt128(0), raw_results, lbt, UInt128(0), 4, auts, w_min, w_max, k, res_lock
    )

    verbose && println("Search complete. Found $(length(raw_results)) representative(s).")

    final_words = Vector{Vector{typeof(elements[1])}}()
    for (msg, tail) in raw_results
        word = fill(elements[1], n)
        for i in 1:k
            val = (msg >> (2 * (i - 1))) & 3
            word[i] = elements[val + 1]
        end
        for i in 1:(n - k)
            val = (tail >> (2 * (i - 1))) & 3
            word[k + i] = elements[val + 1]
        end
        push!(final_words, word)
    end
    
    return final_words
end

function _words_of_weight_range_ternary_recursive!(
    A_L::Vector{UInt128}, A_H::Vector{UInt128}, depth::Int, picked::Int, 
    curr_tail_L::UInt128, curr_tail_H::UInt128, results::Vector{NTuple{4, UInt128}}, 
    lbt::Vector{Int}, msg_L::UInt128, msg_H::UInt128, spawn_depth::Int, 
    auts::Vector{Vector{Int}}, w_min::Int, w_max::Int, k::Int, res_lock::ReentrantLock
)
    tw = count_ones(curr_tail_L | curr_tail_H)
    
    if (picked + tw + lbt[k - depth + 1]) > w_max
        return 
    end

    if picked > 1 && !isempty(auts)
        if !_is_canonical_q(msg_L, msg_H, auts, depth) return end
    end

    if depth == k
        total_w = picked + tw
        if w_min <= total_w <= w_max
            lock(res_lock) do
                push!(results, (msg_L, msg_H, curr_tail_L, curr_tail_H))
            end
        end
        return
    end

    row_L = A_L[depth + 1]
    row_H = A_H[depth + 1]

    if depth < spawn_depth
        tasks = Task[]
        
        nL1, nH1 = _add_gf3_packed(curr_tail_L, curr_tail_H, row_L, row_H)
        t1 = Threads.@spawn _words_of_weight_range_ternary_recursive!(
            A_L, A_H, depth + 1, picked + 1, nL1, nH1, results, lbt, 
            msg_L | (UInt128(1) << depth), msg_H, spawn_depth, auts, w_min, w_max, k, res_lock
        )
        push!(tasks, t1)

        nL2, nH2 = _add_gf3_packed(curr_tail_L, curr_tail_H, row_H, row_L)
        t2 = Threads.@spawn _words_of_weight_range_ternary_recursive!(
            A_L, A_H, depth + 1, picked + 1, nL2, nH2, results, lbt, 
            msg_L, msg_H | (UInt128(1) << depth), spawn_depth, auts, w_min, w_max, k, res_lock
        )
        push!(tasks, t2)

        _words_of_weight_range_ternary_recursive!(
            A_L, A_H, depth + 1, picked, curr_tail_L, curr_tail_H, results, lbt, 
            msg_L, msg_H, spawn_depth, auts, w_min, w_max, k, res_lock
        )
        for t in tasks wait(t) end
    else
        nL1, nH1 = _add_gf3_packed(curr_tail_L, curr_tail_H, row_L, row_H)
        _words_of_weight_range_ternary_recursive!(
            A_L, A_H, depth + 1, picked + 1, nL1, nH1, results, lbt, 
            msg_L | (UInt128(1) << depth), msg_H, spawn_depth, auts, w_min, w_max, k, res_lock
        )

        nL2, nH2 = _add_gf3_packed(curr_tail_L, curr_tail_H, row_H, row_L)
        _words_of_weight_range_ternary_recursive!(
            A_L, A_H, depth + 1, picked + 1, nL2, nH2, results, lbt, 
            msg_L, msg_H | (UInt128(1) << depth), spawn_depth, auts, w_min, w_max, k, res_lock
        )

        _words_of_weight_range_ternary_recursive!(
            A_L, A_H, depth + 1, picked, curr_tail_L, curr_tail_H, results, lbt, 
            msg_L, msg_H, spawn_depth, auts, w_min, w_max, k, res_lock
        )
    end
end

function _words_of_weight_range_dfs_ternary(C::AbstractLinearCode, w_min::Int, w_max::Int, verbose::Bool)
    k, n = C.k, C.n
    F = C.F
    elements = collect(F)
    elem_to_idx = Dict(elements[1] => 0, elements[2] => 1, elements[3] => 2)
    
    G_sys = Array(generator_matrix(C)) 
    A_raw = G_sys[:, (k + 1):n]
    
    A_L = zeros(UInt128, k)
    A_H = zeros(UInt128, k)
    for row in 1:k
        val_L = UInt128(0)
        val_H = UInt128(0)
        for col in 1:(n - k)
            matrix_val = elem_to_idx[A_raw[row, col]]
            if matrix_val == 1 val_L |= (UInt128(1) << (col - 1))
            elseif matrix_val == 2 val_H |= (UInt128(1) << (col - 1)) end
        end
        A_L[row] = val_L
        A_H[row] = val_H
    end
    
    lbt = zeros(Int, k + 1)
    auts = Vector{Vector{Int}}()
    
    raw_results = Vector{NTuple{4, UInt128}}()
    res_lock = ReentrantLock()
    
    verbose && println("Executing Hardware DFS (GF(3)) for weight range [$w_min, $w_max]...")
    _words_of_weight_range_ternary_recursive!(
        A_L, A_H, 0, 0, UInt128(0), UInt128(0), raw_results, lbt, 
        UInt128(0), UInt128(0), 4, auts, w_min, w_max, k, res_lock
    )

    verbose && println("Search complete. Found $(length(raw_results)) representative(s).")

    final_words = Vector{Vector{typeof(elements[1])}}()
    for (mL, mH, tL, tH) in raw_results
        word = fill(elements[1], n)
        for i in 1:k
            if ((mL >> (i - 1)) & 1) == 1 word[i] = elements[2]
            elseif ((mH >> (i - 1)) & 1) == 1 word[i] = elements[3] end
        end
        for i in 1:(n - k)
            if ((tL >> (i - 1)) & 1) == 1 word[k + i] = elements[2]
            elseif ((tH >> (i - 1)) & 1) == 1 word[k + i] = elements[3] end
        end
        push!(final_words, word)
    end
    
    return final_words
end

function _words_of_weight_range_nonbinary_recursive!(
    A_idx::Matrix{UInt8}, depth::Int, picked::Int, 
    curr_tail_idx::Vector{UInt8}, results::Vector{Vector{UInt8}}, 
    lbt::Vector{Int}, msg_so_far_idx::Vector{UInt8}, spawn_depth::Int, 
    q::Int, add_table::Matrix{UInt8}, mul_table::Matrix{UInt8}, 
    auts::Vector{Vector{Int}}, w_min::Int, w_max::Int, res_lock::ReentrantLock
)
    k = size(A_idx, 2)
    n_tail = size(A_idx, 1)

    tw = 0
    @inbounds for i in 1:n_tail
        if curr_tail_idx[i] != 1 tw += 1 end
    end
    
    if (picked + tw + lbt[k - depth + 1]) > w_max
        return 
    end

    if picked > 1 && !isempty(auts)
        if !_is_canonical_q(msg_so_far_idx, auts, q, mul_table) return end
    end

    if depth == k
        total_w = picked + tw
        if w_min <= total_w <= w_max
            full_word = vcat(msg_so_far_idx, curr_tail_idx)
            lock(res_lock) do
                push!(results, full_word) 
            end
        end
        return
    end

    for alpha_idx in 2:q
        msg_so_far_idx[depth + 1] = alpha_idx
        old_tail = copy(curr_tail_idx) 
        
        @inbounds for i in 1:n_tail
            prod_idx = mul_table[alpha_idx, A_idx[i, depth + 1]]
            curr_tail_idx[i] = add_table[curr_tail_idx[i], prod_idx]
        end

        if depth < spawn_depth
            t = Threads.@spawn _words_of_weight_range_nonbinary_recursive!(
                A_idx, depth + 1, picked + 1, copy(curr_tail_idx), results, lbt, 
                copy(msg_so_far_idx), spawn_depth, q, add_table, mul_table, auts, w_min, w_max, res_lock
            )
            wait(t) 
        else
            _words_of_weight_range_nonbinary_recursive!(
                A_idx, depth + 1, picked + 1, curr_tail_idx, results, lbt, 
                msg_so_far_idx, spawn_depth, q, add_table, mul_table, auts, w_min, w_max, res_lock
            )
        end
        
        @inbounds for i in 1:n_tail curr_tail_idx[i] = old_tail[i] end
    end

    msg_so_far_idx[depth + 1] = 1
    _words_of_weight_range_nonbinary_recursive!(
        A_idx, depth + 1, picked, curr_tail_idx, results, lbt, 
        msg_so_far_idx, spawn_depth, q, add_table, mul_table, auts, w_min, w_max, res_lock
    )
end

function _words_of_weight_range_dfs_nonbinary(C::AbstractLinearCode, w_min::Int, w_max::Int, verbose::Bool)
    k, n = C.k, C.n
    F = C.F
    q_order = length(collect(F))
    
    add_t, mul_t, q = _prepare_gf_tables(F) # Assuming table builder exists
    elements = collect(F)
    elem_to_idx = Dict(el => i for (i, el) in enumerate(elements))
    
    G_sys = Array(generator_matrix(C))
    A_raw = G_sys[:, (k + 1):n]
    A_idx = [elem_to_idx[A_raw[i, j]] for i in 1:size(A_raw, 1), j in 1:size(A_raw, 2)]
    
    lbt = zeros(Int, k + 1)
    auts = Vector{Vector{Int}}()
    
    raw_results = Vector{Vector{UInt8}}()
    res_lock = ReentrantLock()
    
    verbose && println("Executing Generic Non-Binary DFS for weight range [$w_min, $w_max]...")
    _words_of_weight_range_nonbinary_recursive!(
        A_idx, 0, 0, ones(UInt8, n - k), raw_results, lbt, 
        ones(UInt8, k), 2, q, add_t, mul_t, auts, w_min, w_max, res_lock
    )

    verbose && println("Search complete. Found $(length(raw_results)) representative(s).")
    
    return [[elements[idx] for idx in word] for word in raw_results]
end

# TODO add C.d saving to all relevant functions

function Base.show(io::IO, hwe::HammingWeightEnumerator)
    terms = String[]
    # Sort by weight so the polynomial prints in standard order (e.g., x^n first)
    for w in sort(collect(keys(hwe.counts)))
        c = hwe.counts[w]
        if c == 0 continue end
        
        # Format coefficient
        coeff_str = c == 1 ? "" : "$c"
        
        # Format x term (x represents zeros)
        p_x = hwe.n - w
        x_str = p_x == 0 ? "" : (p_x == 1 ? "x" : "x^$p_x")
        
        # Format y term (y represents ones/non-zeros)
        y_str = w == 0 ? "" : (w == 1 ? "y" : "y^$w")
        
        # Join them cleanly
        term = join(filter(!isempty, [coeff_str, x_str, y_str]), "*")
        push!(terms, isempty(term) ? "1" : term)
    end
    
    print(io, isempty(terms) ? "0" : join(terms, " + "))
end

"""
    polynomial(hwe::HammingWeightEnumerator, R::Oscar.MPolyRing)

Converts the HammingWeightEnumerator into an actual Oscar polynomial.
Requires a bivariate polynomial ring `R`, e.g., `R, (x, y) = PolynomialRing(ZZ, ["x", "y"])`.
"""
function polynomial(hwe::HammingWeightEnumerator, R)
    gens = Oscar.gens(R)
    @assert length(gens) >= 2 "Ring must have at least two variables (x and y)."
    x, y = gens[1], gens[2]
    
    poly = R() # Initialize empty polynomial in the ring
    for (w, count) in hwe.counts
        poly += count * (x^(hwe.n - w)) * (y^w)
    end
    
    return poly
end

"""
    HammingWeightEnumerator(cwe::CompleteWeightEnumerator)

Reduces a Complete Weight Enumerator down to a Homogeneous Hamming Weight Enumerator.
"""
function HammingWeightEnumerator(cwe::CompleteWeightEnumerator)
    hwe_counts = Dict{Int, BigInt}()
    
    for (counts, coeff) in cwe.counts
        # The first element in the tuple corresponds to the '0' field element.
        # The Hamming weight is the total length minus the number of zeros.
        hamming_weight = sum(counts) - counts[1]
        
        # Aggregate the coefficients of all configurations sharing this weight
        hwe_counts[hamming_weight] = get(hwe_counts, hamming_weight, BigInt(0)) + coeff
    end
    
    return HammingWeightEnumerator(cwe.n, hwe_counts)
end

function Base.show(io::IO, cwe::CompleteWeightEnumerator)
    terms = String[]
    
    # Sort lexicographically descending for standard polynomial presentation
    for counts in sort(collect(keys(cwe.counts)), rev=true)
        c = cwe.counts[counts]
        if c == 0 continue end
        
        # Format the coefficient
        coeff_str = c == 1 ? "" : "$c"
        
        # Format the z_i variables (z0 is the zero element, z1 is the first non-zero, etc.)
        var_strs = String[]
        for (i, p) in enumerate(counts)
            if p == 0 continue end
            var_name = "z$(i - 1)"
            push!(var_strs, p == 1 ? var_name : "$(var_name)^$p")
        end
        
        # Join coefficient and variables cleanly
        term = join(filter(!isempty, vcat([coeff_str], var_strs)), "*")
        push!(terms, isempty(term) ? "1" : term)
    end
    
    print(io, isempty(terms) ? "0" : join(terms, " + "))
end

"""
    polynomial(cwe::CompleteWeightEnumerator, R)

Converts the CompleteWeightEnumerator into an Oscar multivariate polynomial. 
Requires an MPolyRing `R` with at least `q` variables (e.g., z0, z1, ..., z_{q-1}).
"""
function polynomial(cwe::CompleteWeightEnumerator, R)
    gens = Oscar.gens(R)
    @assert length(gens) >= cwe.q "Ring must have at least $(cwe.q) variables."
    
    poly = R() # Initialize empty polynomial
    for (counts, coeff) in cwe.counts
        term = R(coeff)
        for (i, p) in enumerate(counts)
            if p > 0
                term *= gens[i]^p
            end
        end
        poly += term
    end
    
    return poly
end

"""
    MacWilliams_transform(dual_hwe::HammingWeightEnumerator, k::Int, q::Int)

Applies the MacWilliams identity to a Dual HWE to obtain the Primal HWE.
Uses the highly optimized Krawtchouk polynomial evaluation over the internal dictionary.
"""
function MacWilliams_transform(dual_hwe::HammingWeightEnumerator, k::Int, q::Int)
    # Route directly to your existing combinatorial Krawtchouk engine
    primal_counts = MacWilliams_HWE_transform(dual_hwe.counts, dual_hwe.n, k, q)
    return HammingWeightEnumerator(dual_hwe.n, primal_counts)
end

"""
    weight_distribution(C::AbstractLinearCode; verbose::Bool=false)

Calculates or retrieves the cached Hamming weight distribution of the code.
Returns a `Dict{Int, BigInt}`.
"""
function weight_distribution(C::AbstractLinearCode; verbose::Bool=false)
    return get!(C.cache, :weight_dist) do
        k, n = C.k, C.n
        q = Int(order(C.F))
        
        if k > n / 2
            verbose && println("High-rate code detected. Computing Dual Trellis...")
            dual_C = dual(C)
            dual_counts = _weight_distribution_trellis(dual_C, num_trials=50, verbose=verbose)
            
            verbose && println("Applying MacWilliams Transform...")
            # Implicitly returned and cached:
            MacWilliams_HWE_transform(dual_counts, n, k, q) 
        else
            verbose && println("Routing to Primal Trellis Product...")
            # Implicitly returned and cached:
            _weight_distribution_trellis(C, num_trials=50, verbose=verbose)
        end
    end
end

"""
    complete_weight_distribution(C::AbstractLinearCode; verbose::Bool=false)

Calculates the Complete Weight Distribution (CWD) of the code.
Returns a `Dict{Tuple, BigInt}` mapping field element frequencies to codeword count.
"""
function complete_weight_distribution(C::AbstractLinearCode; verbose::Bool=false)
    k, n = C.k, C.n
    q = Int(order(C.F))
    
    # The MacWilliams Shortcut
    if k > n / 2
        verbose && println("High-rate code detected. Computing Dual CWE Trellis...")
        dual_C = dual(C)
        dual_counts = _complete_weight_enumerator_trellis(dual_C, num_trials=50, verbose=verbose)
        
        verbose && println("Applying CWE MacWilliams Transform...")
        return MacWilliams_CWE_transform(dual_counts, n, k, q) # Assuming this raw Dict transformer exists
    end
    
    # Primal Trellis Routing
    verbose && println("Routing to Primal CWE Trellis Product...")
    return _complete_weight_enumerator_trellis(C, num_trials=50, verbose=verbose)
end

"""
    weight_enumerator(C::AbstractLinearCode; verbose::Bool=false)

Retrieves the cached Hamming Weight Enumerator, or builds it from the distribution.
Returns a `HammingWeightEnumerator` struct.
"""
function weight_enumerator(C::AbstractLinearCode; verbose::Bool=false)
    return get!(C.cache, :weight_enum) do
        # This will INSTANTLY return if weight_distribution(C) was already called!
        counts = weight_distribution(C, verbose=verbose)
        
        # Build the wrapper struct, implicitly returning and caching it
        HammingWeightEnumerator(C.n, counts)
    end
end

function complete_weight_distribution(C::AbstractLinearCode; verbose::Bool=false)
    return get!(C.cache, :cwe_dist) do
        # We use your clever syndrome Trellis wrapper here!
        _complete_weight_enumerator_trellis(C, num_trials=50, verbose=verbose)
    end
end

function complete_weight_enumerator(C::AbstractLinearCode; verbose::Bool=false)
    return get!(C.cache, :cwe_enum) do
        counts = complete_weight_distribution(C, verbose=verbose)
        F = C.F
        q = Int(order(F))
        elements = collect(F)
        
        CompleteWeightEnumerator(C.n, q, elements, counts)
    end
end

"""
    weight_distribution_array(C::AbstractLinearCode; verbose::Bool=false)

Calculates the Hamming weight distribution of the code.
Returns a `Vector{BigInt}` of length `n + 1`, where the `i`-th element 
is the number of codewords of weight `i - 1`.
"""
function weight_distribution_array(C::AbstractLinearCode; verbose::Bool=false)
    # Fetch the raw dictionary
    counts_dict = weight_distribution(C; verbose=verbose)
    
    # Initialize a dense array of BigInts for weights 0 through n
    arr = zeros(BigInt, C.n + 1)
    
    for (w, count) in counts_dict
        arr[w + 1] = count
    end
    
    return arr
end

"""
    weight_plot(C::AbstractLinearCode; alg::Symbol = :auto)

Return a bar graph of the weight distribution of `C`.

# Note
- Run `using Makie` to activate this extension.
"""
function weight_plot end

# ==============================================================================
# WEIGHT ALGEBRA FOR COMPOSITE CODES
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Compute the exact Hamming Weight Enumerator for a Direct Sum code via discrete convolution.
"""
function weight_enumerator(C::DirectSumCode)
    cache = getfield(C, :cache)
    if !haskey(cache, :weight_enum)
        we1 = weight_enumerator(C.C1)
        we2 = weight_enumerator(C.C2)
        
        # If either constituent code doesn't have a known WE, we can't compute it
        (ismissing(we1) || ismissing(we2)) && return missing
        
        new_counts = Dict{Int, BigInt}()
        for (w1, count1) in we1.weights
            for (w2, count2) in we2.weights
                w_new = w1 + w2
                new_counts[w_new] = get(new_counts, w_new, BigInt(0)) + (count1 * count2)
            end
        end
        cache[:weight_enum] = HammingWeightEnumerator(C.n, new_counts)
    end
    return cache[:weight_enum]
end
