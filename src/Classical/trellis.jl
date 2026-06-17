# Copyright (c) 2022 - 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
     # Helper Functions
#############################

function _build_trellis_sets(L::Vector{Int}, R::Vector{Int}, boundaries::Vector{Int})
    num_sections = length(boundaries) - 1
    k = length(L)
    
    start_sets   = [Int[] for _ in 1:num_sections] 
    working_sets = [Int[] for _ in 1:num_sections] 
    active_sets  = [Int[] for _ in 1:num_sections] 
    ending_sets  = [Int[] for _ in 1:num_sections]
    
    for m in 1:num_sections
        left_b, right_b = boundaries[m], boundaries[m+1]
        for row in 1:k
            if L[row] <= right_b < R[row] push!(active_sets[m], row) end
            if left_b < L[row] <= right_b push!(start_sets[m], row) end
            if L[row] <= right_b && R[row] > left_b push!(working_sets[m], row) end
            # Track when a span strictly terminates inside this section
            if left_b < R[row] <= right_b push!(ending_sets[m], row) end
        end
    end
    
    return start_sets, working_sets, active_sets, ending_sets
end

"""
    _make_trellis_oriented!(G::Matrix{T}) where T

Transforms a GF(q) generator matrix into Trellis-Oriented Form (Minimal Span Form) in-place.
This guarantees that all rows have unique start indices and unique end indices, 
minimizing the peak state complexity of the Viterbi trellis.
"""
function _make_trellis_oriented!(G::Matrix{T}) where T
    k, n = size(G)
    
    # ---------------------------------------------------------
    # PHASE 1: Forward Elimination (Make Start Indices Unique)
    # ---------------------------------------------------------
    r = 1
    for c in 1:n
        # Find pivot
        pivot_row = 0
        for i in r:k
            if !iszero(G[i, c])
                pivot_row = i
                break
            end
        end
        
        if pivot_row == 0
            continue
        end
        
        # Swap pivot to current row
        if pivot_row != r
            for j in 1:n
                G[r, j], G[pivot_row, j] = G[pivot_row, j], G[r, j]
            end
        end
        
        # Normalize pivot to 1
        inv_p = inv(G[r, c])
        for j in c:n
            G[r, j] = G[r, j] * inv_p
        end
        
        # Eliminate all entries BELOW the pivot
        for i in r+1:k
            if !iszero(G[i, c])
                factor = G[i, c]
                for j in c:n
                    G[i, j] -= factor * G[r, j]
                end
            end
        end
        r += 1
        if r > k break end
    end
    
    # ---------------------------------------------------------
    # PHASE 2: Backward Span Reduction (Make End Indices Unique)
    # ---------------------------------------------------------
    while true
        # 1. Map current starts and ends
        starts = zeros(Int, k)
        ends   = zeros(Int, k)
        
        for i in 1:k
            s_idx = findfirst(!iszero, view(G, i, :))
            e_idx = findlast(!iszero, view(G, i, :))
            
            starts[i] = isnothing(s_idx) ? n + 1 : s_idx
            ends[i]   = isnothing(e_idx) ? 0     : e_idx
        end
        
        # 2. Search for a collision in the end indices
        collision_found = false
        for i in 1:k
            if ends[i] == 0 continue end # Ignore zero rows
            
            for j in i+1:k
                if ends[i] == ends[j]
                    
                    # COLLISION! We must use the row with the LARGER start index 
                    # to eliminate the end element of the row with the SMALLER start index.
                    # Because the source has zeros where the target starts, 
                    # the target's start index will remain completely untouched!
                    
                    if starts[i] < starts[j]
                        target, source = i, j
                    else
                        target, source = j, i
                    end
                    
                    # Eliminate the end element of the target row
                    factor = G[target, ends[target]] * inv(G[source, ends[source]])
                    
                    for c in 1:n
                        G[target, c] -= factor * G[source, c]
                    end
                    
                    collision_found = true
                    break
                end
            end
            if collision_found break end
        end
        
        # If all end indices are unique, the matrix is in Minimal Span Form!
        if !collision_found
            break
        end
    end
end

"""
    _make_trellis_oriented!(G::CTMatrixTypes)

Wrapper to natively handle Oscar matrix types. Extracts to a standard Julia Array 
for fast O(1) memory access during the heavy row-reduction, then mutates the original 
Oscar matrix in-place.
"""
function _make_trellis_oriented!(G::CTMatrixTypes)
    # Extract to fast standard array
    G_mat = Array(G)
    
    # Run the optimized engine
    _make_trellis_oriented!(G_mat)
    
    # Pack the results back into the Oscar matrix in-place
    for i in 1:nrows(G)
        for j in 1:ncols(G)
            G[i, j] = G_mat[i, j]
        end
    end
    
    return G
end

"""
    _weight_distribution_trellis(C::AbstractLinearCode; num_trials::Int=50, verbose::Bool=false)

Computes the Hamming weight distribution of a linear code.
Uses the primal generator Trellis Product for low-rate codes, and the dual generator 
Trellis Product (followed by the MacWilliams Identity) for high-rate codes.
"""
function _weight_distribution_trellis(C::AbstractLinearCode; num_trials::Int=50, verbose::Bool=false)
    k = C.k
    n = C.n
    q = Int(order(C.F))
    
    is_generator = k <= n / 2
    
    if is_generator
        verbose && println("Low-rate code. Optimizing primal Generator matrix...")
        mat = Array(generator_matrix(C))
    else
        verbose && println("High-rate code. Optimizing dual Parity-Check matrix...")
        mat = Array(parity_check_matrix(C))
    end
    
    # Optimize the matrix
    best_M, best_perm, best_peak_E = optimize_trellis_permutation(mat, num_trials)
    boundaries = optimal_sectionalization(best_M, q)
    
    verbose && println("Executing Generator Trellis Product engine...")
    # Pass verbose explicitly down here!
    CWE_dict = _CWE_classical_TP_sectionalized(best_M, boundaries, verbose)
    
    HWE_dict = _CWE_to_HWE_dict(CWE_dict)
    
    if is_generator
        verbose && println("Primal distribution computed directly.")
        return HWE_dict
    else
        verbose && println("Applying Krawtchouk-MacWilliams Transform to dual distribution...")
        primal_hwe = Macwilliams_HWE_transform(HWE_dict, n, k, q)
        return primal_hwe
    end
end

"""
    _distance_from_CWE_dict(CWE_dict::Dict{NTuple{Q, Int}, BigInt}) where Q

Extracts the minimum distance directly from a Complete Weight Enumerator dictionary.
Calculates the Hamming weight by summing the tuple and subtracting the count 
of the zero element (index 1). Returns the minimum strictly positive weight.
"""
function _distance_from_CWE_dict(CWE_dict::Dict{NTuple{Q, Int}, BigInt}) where Q
    min_d = typemax(Int)
    
    for counts in keys(CWE_dict)
        # Hamming weight = total length - number of zeros
        w = sum(counts) - counts[1]
        
        if w > 0 && w < min_d
            min_d = w
        end
    end
    
    return min_d == typemax(Int) ? -1 : min_d
end

"""
    _CWE_to_HWE_dict(CWE_dict::Dict{NTuple{Q, Int}, BigInt}) where Q

Reduces a Complete Weight Enumerator (CWE) dictionary down to a Homogeneous 
Weight Enumerator (HWE) dictionary. 

Returns a `Dict{Int, BigInt}` mapping the Hamming weight to the total number 
of codewords possessing that weight.
"""
function _CWE_to_HWE_dict(CWE_dict::Dict{NTuple{Q, Int}, BigInt}) where Q
    HWE_dict = Dict{Int, BigInt}()
    
    for (counts, coeff) in CWE_dict
        # The first element in the tuple counts the '0' field symbol.
        # Hamming weight is the total length minus the number of zeros.
        hamming_weight = sum(counts) - counts[1]
        
        # Aggregate the counts of all compositions that share this Hamming weight
        HWE_dict[hamming_weight] = get(HWE_dict, hamming_weight, BigInt(0)) + coeff
    end
    
    return HWE_dict
end

function _min_weight_TP_Viterbi_binary(G::Matrix{T}, boundaries::Vector{Int}, verbose::Bool=true) where T
    k, n = size(G)
    L, R = _get_LR_indices(G)
    start_sets, working_sets, active_sets, _ = _build_trellis_sets(L, R, boundaries)

    num_sections = length(boundaries) - 1
    prev_layer = Dict{Tuple{UInt128, Bool}, Int}((UInt128(0), false) => 0)
    
    p = verbose ? Progress(num_sections, 1, "Building Viterbi Trellis (Binary): ") : nothing
    
    # Precompute column masks for lightning fast GF(2) dot products
    M_cols = zeros(UInt128, n)
    for col in 1:n
        col_val = UInt128(0)
        for row in 1:k
            if !iszero(G[row, col]) col_val |= (UInt128(1) << (row - 1)) end
        end
        M_cols[col] = col_val
    end

    for m in 1:num_sections
        next_layer = Dict{Tuple{UInt128, Bool}, Int}()
        
        left_b, right_b = boundaries[m], boundaries[m+1]
        active_curr = active_sets[m]
        starting = start_sets[m]
        working = working_sets[m]
        
        # PERFECT MEMORY PRE-ALLOCATION
        sizehint!(next_layer, 1 << length(active_curr))
        
        # Keep mask to safely clear dead rows so paths merge correctly
        keep_mask = UInt128(0)
        for row in active_curr keep_mask |= (UInt128(1) << (row - 1)) end
        
        w_mask = UInt128(0)
        for row in working w_mask |= (UInt128(1) << (row - 1)) end
        
        num_branches = 1 << length(starting)
        deltas = zeros(UInt128, num_branches)
        for b in 0:(num_branches-1)
            d = UInt128(0)
            for (idx, row) in enumerate(starting)
                if ((b >> (idx - 1)) & 1) == 1 d |= (UInt128(1) << (row - 1)) end
            end
            deltas[b+1] = d
        end
        
        cols_in_chunk = (left_b + 1):right_b
        
        for ((u_prev, is_nonzero), prev_wt) in prev_layer
            @inbounds for b in 1:num_branches
                u_work = u_prev | deltas[b]
                
                chunk_wt = 0
                for col in cols_in_chunk
                    chunk_wt += count_ones(u_work & M_cols[col] & w_mask) & 1
                end
                
                new_wt = prev_wt + chunk_wt
                # b=1 means all zero branch. b>1 means non-zero branch
                new_is_nonzero = is_nonzero | (b > 1)
                u_next = u_work & keep_mask
                
                key = (u_next, new_is_nonzero)
                if !haskey(next_layer, key) || new_wt < next_layer[key]
                    next_layer[key] = new_wt
                end
            end
        end
        prev_layer = next_layer
        verbose && next!(p)
    end
    
    return get(prev_layer, (UInt128(0), true), -1)
end

function _min_weight_TP_Viterbi_sectionalized(G::Matrix{T}, boundaries::Vector{Int}=collect(0:size(G,2)), verbose::Bool=true) where T
    k, n = size(G)
    F = parent(G[1, 1])
    q = length(collect(F))
    
    # ADD THE FAST-PATH ROUTING HERE!
    if q == 2 && k <= 128
        return _min_weight_TP_Viterbi_binary(G, boundaries, verbose)
    end
    
    elements = collect(F)
    T_zero = zero(F)
    
    L, R = _get_LR_indices(G)
    start_sets, working_sets, active_sets, _ = _build_trellis_sets(L, R, boundaries)

    num_sections = length(boundaries) - 1
    prev_layer = Dict{Tuple{Vector{T}, Bool}, Int}((Vector{T}(), false) => 0)
    
    p = verbose ? Progress(num_sections, 0.1, "Building Viterbi Trellis: ") : nothing
    u_buf = fill(T_zero, k)
    
    for m in 1:num_sections
        next_layer = Dict{Tuple{Vector{T}, Bool}, Int}()
        
        active_curr = active_sets[m]
        sizehint!(next_layer, q^length(active_curr))
        
        left_b, right_b = boundaries[m], boundaries[m+1]
        active_prev = m == 1 ? Int[] : active_sets[m-1]
        starting    = start_sets[m]
        working     = working_sets[m]
        
        for ((prev_scalars, is_nonzero), wt) in prev_layer
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
                
                new_wt = wt + chunk_wt
                new_is_nonzero = is_nonzero || any(!iszero, branch_scalars)
                next_scalars = [u_buf[row] for row in active_curr]
                key = (next_scalars, new_is_nonzero)
                
                if !haskey(next_layer, key) || new_wt < next_layer[key]
                    next_layer[key] = new_wt
                end
            end
        end
        prev_layer = next_layer
        verbose && next!(p)
    end
    
    return get(prev_layer, (Vector{T}(), true), -1)
end

function _min_weight_syndrome_sectionalized(H::Matrix{T}, boundaries::Vector{Int}=collect(0:size(H,2)), verbose::Bool=true) where T
    
    r, _ = size(H)
    F = parent(H[1, 1])
    q = length(collect(F))
    
    if q == 2 && r <= 128
        return _min_weight_TP_Viterbi_binary(H, boundaries, verbose)
    end

    elements = collect(F)
    T_zero = zero(F)
    
    L, R = _get_LR_indices(H)
    _, _, active_sets, ending_sets = _build_trellis_sets(L, R, boundaries)

    num_sections = length(boundaries) - 1
    prev_layer = Dict{Tuple{Vector{T}, Bool}, Int}((Vector{T}(), false) => 0)
    
    current_syn_buf = fill(T_zero, r)
    p = verbose ? Progress(num_sections, 1, "Building Syndrome Min-Weight Trellis: ") : nothing

    for m in 1:num_sections
        next_layer = Dict{Tuple{Vector{T}, Bool}, Int}()
        
        active_curr = active_sets[m]
        sizehint!(next_layer, q^length(active_curr))
        
        left_b, right_b = boundaries[m], boundaries[m+1]
        chunk_len = right_b - left_b
        
        active_prev = m == 1 ? Int[] : active_sets[m-1]
        ending_curr = ending_sets[m]
        
        for chunk_symbols in Iterators.product(fill(elements, chunk_len)...)
            chunk_vec = [sym for sym in chunk_symbols]
            
            chunk_wt = count(!iszero, chunk_vec)
            H_chunk = view(H, :, (left_b + 1):right_b)
            partial_syn = H_chunk * chunk_vec
            
            for ((prev_syn, is_nonzero), prev_wt) in prev_layer
                
                fill!(current_syn_buf, T_zero)
                @inbounds for (idx, row) in enumerate(active_prev) current_syn_buf[row] = prev_syn[idx] end
                @inbounds for i in 1:r current_syn_buf[i] += partial_syn[i] end
                
                valid_branch = true
                for row in ending_curr
                    if !iszero(current_syn_buf[row])
                        valid_branch = false
                        break
                    end
                end
                if !valid_branch continue end
                
                next_syn = [current_syn_buf[row] for row in active_curr]
                new_wt = prev_wt + chunk_wt
                new_is_nonzero = is_nonzero || chunk_wt > 0
                key = (next_syn, new_is_nonzero)
                
                if !haskey(next_layer, key) || new_wt < next_layer[key]
                    next_layer[key] = new_wt
                end
            end
        end
        prev_layer = next_layer
        verbose && next!(p)
    end
   
    final_state = fill(T_zero, length(active_sets[end]))
    return get(prev_layer, (final_state, true), -1)
end

function _CWE_classical_TP_sectionalized(G::Matrix{T}, boundaries::Vector{Int}=collect(0:size(G,2)), verbose::Bool=true) where T

    k, _ = size(G)
    F = parent(G[1, 1])
    q = length(collect(F))
    
    if q == 2 && k <= 128
        return _CWE_classical_TP_binary(G, boundaries, verbose)
    elseif q == 3 && k <= 64
        return _CWE_classical_TP_ternary(G, boundaries, verbose)
    elseif q == 4 && k <= 64
        return _CWE_classical_TP_quaternary(G, boundaries, verbose)
    end

    elements = collect(F)
    elem_idx = Dict(elements[i] => i for i in 1:q)
    T_zero = zero(F)
    
    L, R = _get_LR_indices(G)
    start_sets, working_sets, active_sets, _ = _build_trellis_sets(L, R, boundaries)

    num_sections = length(boundaries) - 1
    init_comp = ntuple(i -> 0, q)
    prev_layer = Dict{Vector{T}, Dict{NTuple{q, Int}, BigInt}}(Vector{T}() => Dict(init_comp => BigInt(1)))
    
    p = verbose ? Progress(num_sections, 1, "Building TP CWE Trellis: ") : nothing
    u_buf = fill(T_zero, k)
    
    for m in 1:num_sections
        next_layer = Dict{Vector{T}, Dict{NTuple{q, Int}, BigInt}}()
        
        active_curr = active_sets[m]
        sizehint!(next_layer, q^length(active_curr))
        
        left_b, right_b = boundaries[m], boundaries[m+1]
        active_prev = m == 1 ? Int[] : active_sets[m-1]
        starting = start_sets[m]
        working = working_sets[m]
        
        for (prev_scalars, partial_CWE) in prev_layer
            for branch_scalars in Iterators.product(fill(elements, length(starting))...)
                
                fill!(u_buf, T_zero)
                @inbounds for (idx, row) in enumerate(active_prev) u_buf[row] = prev_scalars[idx] end
                @inbounds for (idx, row) in enumerate(starting)    u_buf[row] = branch_scalars[idx] end
                
                chunk_counts = zeros(Int, q)
                for col in (left_b + 1):right_b
                    c_i = T_zero
                    @inbounds for row in working c_i += u_buf[row] * G[row, col] end
                    chunk_counts[elem_idx[c_i]] += 1
                end
                
                chunk_tuple = ntuple(i -> chunk_counts[i], q)
                next_scalars = [u_buf[row] for row in active_curr]
                
                if !haskey(next_layer, next_scalars)
                    next_layer[next_scalars] = Dict{NTuple{q, Int}, BigInt}()
                end
                
                dest_dict = next_layer[next_scalars]
                for (prev_comp, count) in partial_CWE
                    new_comp = ntuple(i -> prev_comp[i] + chunk_tuple[i], q)
                    dest_dict[new_comp] = get(dest_dict, new_comp, BigInt(0)) + count
                end
            end
        end
        prev_layer = next_layer
        verbose && next!(p)
    end
    
    return prev_layer[Vector{T}()]
end

function _CWE_classical_TP_binary(G::Matrix{T}, boundaries::Vector{Int}, verbose::Bool=true) where T
    k, n = size(G)
    L, R = _get_LR_indices(G)
    start_sets, working_sets, active_sets, _ = _build_trellis_sets(L, R, boundaries)

    num_sections = length(boundaries) - 1
    init_comp = (0, 0) # (zeros, ones)
    
    # Dict Key is now UInt128!
    prev_layer = Dict{UInt128, Dict{NTuple{2, Int}, BigInt}}(UInt128(0) => Dict(init_comp => BigInt(1)))

    p = verbose ? Progress(num_sections, 0.1, "Building TP CWE Trellis (Binary): ") : nothing

    M_cols = zeros(UInt128, n)
    for col in 1:n
        col_val = UInt128(0)
        for row in 1:k
            if !iszero(G[row, col]) col_val |= (UInt128(1) << (row - 1)) end
        end
        M_cols[col] = col_val
    end

    for m in 1:num_sections
        next_layer = Dict{UInt128, Dict{NTuple{2, Int}, BigInt}}()

        active_curr = active_sets[m]
        sizehint!(next_layer, 1 << length(active_curr))

        left_b, right_b = boundaries[m], boundaries[m+1]
        starting = start_sets[m]
        working = working_sets[m]

        keep_mask = UInt128(0)
        for row in active_curr keep_mask |= (UInt128(1) << (row - 1)) end

        w_mask = UInt128(0)
        for row in working w_mask |= (UInt128(1) << (row - 1)) end

        num_branches = 1 << length(starting)
        deltas = zeros(UInt128, num_branches)
        for b in 0:(num_branches-1)
            d = UInt128(0)
            for (idx, row) in enumerate(starting)
                if ((b >> (idx - 1)) & 1) == 1 d |= (UInt128(1) << (row - 1)) end
            end
            deltas[b+1] = d
        end

        cols_in_chunk = (left_b + 1):right_b
        chunk_len = right_b - left_b

        for (u_prev, partial_CWE) in prev_layer
            @inbounds for b in 1:num_branches
                u_work = u_prev | deltas[b]

                ones_wt = 0
                for col in cols_in_chunk
                    ones_wt += count_ones(u_work & M_cols[col] & w_mask) & 1
                end
                zeros_wt = chunk_len - ones_wt
                chunk_tuple = (zeros_wt, ones_wt)

                u_next = u_work & keep_mask

                if !haskey(next_layer, u_next)
                    next_layer[u_next] = Dict{NTuple{2, Int}, BigInt}()
                end

                dest_dict = next_layer[u_next]
                for (prev_comp, count) in partial_CWE
                    new_comp = (prev_comp[1] + chunk_tuple[1], prev_comp[2] + chunk_tuple[2])
                    dest_dict[new_comp] = get(dest_dict, new_comp, BigInt(0)) + count
                end
            end
        end
        prev_layer = next_layer
        verbose && next!(p)
    end

    return prev_layer[UInt128(0)]
end

function _CWE_classical_TP_ternary(G::Matrix{T}, boundaries::Vector{Int}, verbose::Bool=true) where T
    k, n = size(G)
    F = parent(G[1,1])
    elements = collect(F)
    elem_to_u8 = Dict(elements[1] => UInt8(0), elements[2] => UInt8(1), elements[3] => UInt8(2))
    
    G_u8 = zeros(UInt8, k, n)
    for i in 1:k, j in 1:n G_u8[i,j] = elem_to_u8[G[i,j]] end
    
    L, R = _get_LR_indices(G)
    start_sets, working_sets, active_sets, _ = _build_trellis_sets(L, R, boundaries)

    num_sections = length(boundaries) - 1
    init_comp = (0, 0, 0)
    prev_layer = Dict{UInt128, Dict{NTuple{3, Int}, BigInt}}(UInt128(0) => Dict(init_comp => BigInt(1)))

    p = verbose ? Progress(num_sections, 0.1, "Building TP CWE Trellis (GF(3)): ") : nothing
    u_buf = zeros(UInt8, k)

    for m in 1:num_sections
        next_layer = Dict{UInt128, Dict{NTuple{3, Int}, BigInt}}()

        active_curr = active_sets[m]
        sizehint!(next_layer, 3^length(active_curr))

        left_b, right_b = boundaries[m], boundaries[m+1]
        active_prev = m == 1 ? Int[] : active_sets[m-1]
        starting = start_sets[m]
        working = working_sets[m]

        for (u_packed, partial_CWE) in prev_layer
            fill!(u_buf, UInt8(0))
            for (idx, row) in enumerate(active_prev)
                u_buf[row] = UInt8((u_packed >> (2 * (idx - 1))) & 3)
            end

            for branch_scalars in Iterators.product(fill((UInt8(0), UInt8(1), UInt8(2)), length(starting))...)
                @inbounds for (idx, row) in enumerate(starting) u_buf[row] = branch_scalars[idx] end
                
                c0, c1, c2 = 0, 0, 0
                for col in (left_b + 1):right_b
                    c_i = UInt32(0)
                    @inbounds for row in working c_i += u_buf[row] * G_u8[row, col] end
                    val = c_i % 3
                    if val == 0 c0 += 1 elseif val == 1 c1 += 1 else c2 += 1 end
                end
                chunk_tuple = (c0, c1, c2)

                u_next = UInt128(0)
                @inbounds for (idx, row) in enumerate(active_curr)
                    u_next |= (UInt128(u_buf[row]) << (2 * (idx - 1)))
                end

                if !haskey(next_layer, u_next)
                    next_layer[u_next] = Dict{NTuple{3, Int}, BigInt}()
                end
                
                dest_dict = next_layer[u_next]
                for (prev_comp, count) in partial_CWE
                    new_comp = (prev_comp[1] + chunk_tuple[1], prev_comp[2] + chunk_tuple[2], prev_comp[3] + chunk_tuple[3])
                    dest_dict[new_comp] = get(dest_dict, new_comp, BigInt(0)) + count
                end
            end
        end
        prev_layer = next_layer
        verbose && next!(p)
    end
    return prev_layer[UInt128(0)]
end

function _CWE_classical_TP_quaternary(G::Matrix{T}, boundaries::Vector{Int}, verbose::Bool=true) where T
    k, n = size(G)
    F = parent(G[1,1])
    elements = collect(F)
    elem_to_u8 = Dict(elements[1] => UInt8(0), elements[2] => UInt8(1), elements[3] => UInt8(2), elements[4] => UInt8(3))
    
    GF4_MULT = UInt8[0 0 0 0; 0 1 2 3; 0 2 3 1; 0 3 1 2]
    
    G_u8 = zeros(UInt8, k, n)
    for i in 1:k, j in 1:n G_u8[i,j] = elem_to_u8[G[i,j]] end
    
    L, R = _get_LR_indices(G)
    start_sets, working_sets, active_sets, _ = _build_trellis_sets(L, R, boundaries)

    num_sections = length(boundaries) - 1
    init_comp = (0, 0, 0, 0)
    prev_layer = Dict{UInt128, Dict{NTuple{4, Int}, BigInt}}(UInt128(0) => Dict(init_comp => BigInt(1)))

    p = verbose ? Progress(num_sections, 0.1, "Building TP CWE Trellis (GF(4)): ") : nothing
    u_buf = zeros(UInt8, k)

    for m in 1:num_sections
        next_layer = Dict{UInt128, Dict{NTuple{4, Int}, BigInt}}()

        active_curr = active_sets[m]
        sizehint!(next_layer, 4^length(active_curr))

        left_b, right_b = boundaries[m], boundaries[m+1]
        active_prev = m == 1 ? Int[] : active_sets[m-1]
        starting = start_sets[m]
        working = working_sets[m]

        for (u_packed, partial_CWE) in prev_layer
            fill!(u_buf, UInt8(0))
            for (idx, row) in enumerate(active_prev)
                u_buf[row] = UInt8((u_packed >> (2 * (idx - 1))) & 3)
            end

            for branch_scalars in Iterators.product(fill((UInt8(0), UInt8(1), UInt8(2), UInt8(3)), length(starting))...)
                @inbounds for (idx, row) in enumerate(starting) u_buf[row] = branch_scalars[idx] end
                
                c0, c1, c2, c3 = 0, 0, 0, 0
                for col in (left_b + 1):right_b
                    c_i = UInt8(0)
                    @inbounds for row in working c_i ⊻= GF4_MULT[u_buf[row] + 1, G_u8[row, col] + 1] end
                    if c_i == 0 c0 += 1 elseif c_i == 1 c1 += 1 elseif c_i == 2 c2 += 1 else c3 += 1 end
                end
                chunk_tuple = (c0, c1, c2, c3)

                u_next = UInt128(0)
                @inbounds for (idx, row) in enumerate(active_curr)
                    u_next |= (UInt128(u_buf[row]) << (2 * (idx - 1)))
                end

                if !haskey(next_layer, u_next)
                    next_layer[u_next] = Dict{NTuple{4, Int}, BigInt}()
                end
                dest_dict = next_layer[u_next]
                for (prev_comp, count) in partial_CWE
                    new_comp = (prev_comp[1] + chunk_tuple[1], prev_comp[2] + chunk_tuple[2], prev_comp[3] + chunk_tuple[3], prev_comp[4] + chunk_tuple[4])
                    dest_dict[new_comp] = get(dest_dict, new_comp, BigInt(0)) + count
                end
            end
        end
        prev_layer = next_layer
        verbose && next!(p)
    end
    return prev_layer[UInt128(0)]
end

function _optimal_sectionalization_linear(M::Matrix{T}, q::Int) where T
    k, n = size(M)
    
    # 1. Extract Left and Right Indices
    L = zeros(Int, k)
    R = zeros(Int, k)
    for row in 1:k
        l_idx = findfirst(!iszero, view(M, row, :))
        r_idx = findlast(!iszero, view(M, row, :))
        L[row] = isnothing(l_idx) ? n + 1 : l_idx
        R[row] = isnothing(r_idx) ? 0     : r_idx
    end
    
    # 2. THE OPTIMIZATION: Precompute `past` and `future`
    # past[u]: rows that end at or before boundary u
    # future[v]: rows that start strictly after boundary v
    past = zeros(Int, n + 1)
    future = zeros(Int, n + 1)
    for b in 0:n
        past[b + 1] = count(r -> r > 0 && r <= b, R)
        future[b + 1] = count(l -> l > b, L)
    end
    
    # 3. Flat DP Arrays (Avoids Vertex/Edge struct allocations)
    min_cost = fill(Inf, n + 1)
    min_cost[1] = 0.0 # Cost to reach boundary 0
    parent = zeros(Int, n + 1)
    
    # 4. O(n^2) DAG Traversal
    for v in 1:n
        v_idx = v + 1
        for u in 0:(v - 1)
            u_idx = u + 1
            
            # The brilliant O(1) intersection math
            # (Note: dim_ker is subtracted here for quantum codes)
            active_count = k - past[u_idx] - future[v_idx] 
            
            # Your safety cutoff heuristic
            if active_count > 50
                continue # Equivalent to setting edge weight to Inf
            end
            
            section_cost = Float64(q)^active_count
            
            # DP Relaxation
            if min_cost[u_idx] + section_cost < min_cost[v_idx]
                min_cost[v_idx] = min_cost[u_idx] + section_cost
                parent[v_idx] = u
            end
        end
    end
    
    # 5. Reconstruct Path
    boundaries = Int[]
    curr = n
    while curr != 0
        push!(boundaries, curr)
        curr = parent[curr + 1]
    end
    push!(boundaries, 0)
    
    return reverse(boundaries)
end

function past_future_profiles(L::Vector{Int}, R::Vector{Int}, n::Int)
    past = zeros(Int, n + 1)
    future = zeros(Int, n + 1)
    for b in 0:n
        # past[b+1]: number of rows that end at or before boundary b
        past[b + 1] = count(r -> r > 0 && r <= b, R)
        # future[b+1]: number of rows that start strictly after boundary b
        future[b + 1] = count(l -> l > b, L)
    end
    return past, future
end

"""
    vertex_counts(k::Int, n::Int, past::Vector{Int}, future::Vector{Int}, boundaries::Vector{Int}=collect(0:n))

Returns an array representing the number of active generators (vertex exponent) 
exactly at each chosen section boundary for the Trellis Product graph. 
Defaults to the unsectionalized step-by-step counts if `boundaries` is omitted.
"""
function vertex_counts(k::Int, n::Int, past::Vector{Int}, future::Vector{Int}, boundaries::Vector{Int}=collect(0:n))
    num_b = length(boundaries)
    V_sect = zeros(Int, num_b)
    
    for i in 1:num_b
        b_idx = boundaries[i] + 1
        V_sect[i] = k - past[b_idx] - future[b_idx]
    end
    
    return V_sect
end

"""
    edge_counts(k::Int, n::Int, past::Vector{Int}, future::Vector{Int}, boundaries::Vector{Int}=collect(0:n))

Returns an array representing the edge complexity (number of active generators) 
for each macro-section defined by the boundaries for the Trellis Product graph.
Defaults to the unsectionalized step-by-step counts if `boundaries` is omitted.
"""
function edge_counts(k::Int, n::Int, past::Vector{Int}, future::Vector{Int}, boundaries::Vector{Int}=collect(0:n))
    num_sections = length(boundaries) - 1
    E_sect = zeros(Int, num_sections)
    
    for i in 1:num_sections
        u_idx = boundaries[i] + 1
        v_idx = boundaries[i+1] + 1
        
        # The O(1) intersection math for the macro-edge span
        E_sect[i] = k - past[u_idx] - future[v_idx]
    end
    
    return E_sect
end

"""
    optimize_trellis_permutation(M::Matrix{T}, num_trials::Int=50) where T

Applies random column permutations to matrix M, transforms each into Trellis-Oriented 
Form (TOF), and profiles the maximum active edge span. Returns the permuted TOF 
matrix that yields the smallest peak trellis complexity.
"""
function optimize_trellis_permutation(M::Matrix{T}, num_trials::Int=50) where T
    k, n = size(M)
    
    # Baseline
    best_M = copy(M)
    _make_trellis_oriented!(best_M)
    L, R = _get_LR_indices(best_M)
    past, future = past_future_profiles(L, R, n)
    best_peak_E = maximum([k - past[b] - future[b+1] for b in 1:n])
    best_perm = collect(1:n)
    
    if num_trials == 0 return best_M, best_perm, best_peak_E end

    # Thread-safe storage
    results = Vector{Tuple{Matrix{T}, Vector{Int}, Int}}(undef, num_trials)
    
    Threads.@threads for i in 1:num_trials
        perm = randperm(n)
        M_shuffled = M[:, perm]
        _make_trellis_oriented!(M_shuffled)
        
        L_s, R_s = _get_LR_indices(M_shuffled)
        past_s, future_s = past_future_profiles(L_s, R_s, n)
        peak_E_s = maximum([k - past_s[b] - future_s[b+1] for b in 1:n])
        
        results[i] = (M_shuffled, perm, peak_E_s)
    end
    
    # Find the best result from all threads
    for (M_s, perm_s, peak_E_s) in results
        if peak_E_s < best_peak_E
            best_peak_E = peak_E_s
            best_M = M_s
            best_perm = perm_s
        end
    end
    
    return best_M, best_perm, best_peak_E
end

"""
    _get_LR_indices(M::Matrix{T}) where T

Returns the left (L) and right (R) bounding indices for each row of the matrix M.
L[i] is the column index of the first non-zero element in row i.
R[i] is the column index of the last non-zero element in row i.

If row i is all zeros, L[i] is set to n + 1 and R[i] is set to 0 
so that it does not artificially inflate the active span complexity.
"""
function _get_LR_indices(M::Matrix{T}) where T
    k, n = size(M)
    
    L = zeros(Int, k)
    R = zeros(Int, k)
    
    for row in 1:k
        l_idx = findfirst(!iszero, view(M, row, :))
        r_idx = findlast(!iszero, view(M, row, :))
        
        L[row] = isnothing(l_idx) ? n + 1 : l_idx
        R[row] = isnothing(r_idx) ? 0     : r_idx
    end
    
    return L, R
end
_get_LR_indices(M::CTMatrixTypes) = _get_LR_indices(Array(M))

"""
    _minimum_distance_trellis(C::AbstractLinearCode; num_trials::Int=50, verbose::Bool=false)

Computes the minimum distance by building the optimally permuted and sectionalized Trellis.
"""
function _minimum_distance_trellis(C::AbstractLinearCode; num_trials::Int=50, verbose::Bool=false)
    k, n, q = C.k, C.n, Int(order(C.F))
    
    is_generator = k <= n / 2
    mat = is_generator ? Array(generator_matrix(C)) : Array(parity_check_matrix(C))
    
    verbose && println("Optimizing trellis via column permutations...")
    best_M, best_perm, best_peak_E = optimize_trellis_permutation(mat, num_trials)
    
    verbose && println("Computing optimal Lafourcade-Vardy sectionalization boundaries...")
    boundaries = optimal_sectionalization(best_M, q)
    
    # DYNAMICALLY ROUTE TO MIN-WEIGHT (VITERBI) INSTEAD OF CWE
    if is_generator
        return _min_weight_TP_Viterbi_sectionalized(best_M, boundaries, verbose)
    else
        return _min_weight_syndrome_sectionalized(best_M, boundaries, verbose)
    end
end

"""
    _complete_weight_enumerator_trellis(C::AbstractLinearCode; num_trials::Int=50, verbose::Bool=false)

Computes the Complete Weight Enumerator using dynamic routing and optimal sectionalization.
"""
function _complete_weight_enumerator_trellis(C::AbstractLinearCode; num_trials::Int=50, verbose::Bool=false)
    k, n, q = C.k, C.n, Int(order(C.F))
    
    is_generator = k <= n / 2
    mat = is_generator ? Array(generator_matrix(C)) : Array(parity_check_matrix(C))
    
    verbose && println("Optimizing trellis via column permutations...")
    best_M, best_perm, best_peak_E = optimize_trellis_permutation(mat, num_trials)
    
    verbose && println("Computing optimal Lafourcade-Vardy sectionalization boundaries...")
    boundaries = optimal_sectionalization(best_M, q)
    
    return is_generator ? _CWE_classical_TP_sectionalized(best_M, boundaries, verbose) : _CWE_classical_syndrome_sectionalized(best_M, boundaries, verbose)
end

function _CWE_classical_syndrome_sectionalized(H::Matrix{T}, boundaries::Vector{Int}=collect(0:size(H,2)), verbose::Bool=true) where T

    r, n = size(H) 
    F = parent(H[1, 1])
    elements = collect(F)
    q = length(elements)
    
   if q == 2 && r <= 128
        return _CWE_classical_syndrome_binary(H, boundaries, verbose)
    elseif q == 3 && r <= 64
        return _CWE_classical_syndrome_ternary(H, boundaries, verbose)
    elseif q == 4 && r <= 64
        return _CWE_classical_syndrome_quaternary(H, boundaries, verbose)
    end

    elem_idx = Dict(elements[i] => i for i in 1:q)
    T_zero = zero(F)
    
    L, R = _get_LR_indices(H)
    _, _, active_sets, ending_sets = _build_trellis_sets(L, R, boundaries)

    num_sections = length(boundaries) - 1
    init_comp = ntuple(i -> 0, q)
    prev_layer = Dict{Vector{T}, Dict{NTuple{q, Int}, BigInt}}(Vector{T}() => Dict(init_comp => BigInt(1)))
    
    current_syn_buf = fill(T_zero, r)
    p = verbose ? Progress(num_sections, 1, "Building Syndrome CWE Trellis: ") : nothing

    for m in 1:num_sections
        next_layer = Dict{Vector{T}, Dict{NTuple{q, Int}, BigInt}}()
        
        active_curr = active_sets[m]
        sizehint!(next_layer, q^length(active_curr))
        
        left_b, right_b = boundaries[m], boundaries[m+1]
        chunk_len = right_b - left_b
        
        active_prev = m == 1 ? Int[] : active_sets[m-1]
        ending_curr = ending_sets[m]
        
        for chunk_symbols in Iterators.product(fill(elements, chunk_len)...)
            chunk_vec = [sym for sym in chunk_symbols]
            
            chunk_counts = zeros(Int, q)
            for sym in chunk_vec chunk_counts[elem_idx[sym]] += 1 end
            chunk_tuple = ntuple(i -> chunk_counts[i], q)
            
            H_chunk = view(H, :, (left_b + 1):right_b)
            partial_syn = H_chunk * chunk_vec
            
            for (prev_syn, partial_CWE) in prev_layer
                
                fill!(current_syn_buf, T_zero)
                @inbounds for (idx, row) in enumerate(active_prev) 
                    current_syn_buf[row] = prev_syn[idx] 
                end
                
                @inbounds for i in 1:r
                    current_syn_buf[i] += partial_syn[i]
                end
                
                valid_branch = true
                for row in ending_curr
                    if !iszero(current_syn_buf[row])
                        valid_branch = false
                        break
                    end
                end
                if !valid_branch continue end
                
                next_syn = [current_syn_buf[row] for row in active_curr]
                
                if !haskey(next_layer, next_syn)
                    next_layer[next_syn] = Dict{NTuple{q, Int}, BigInt}()
                end
                
                dest_dict = next_layer[next_syn]
                for (prev_comp, count) in partial_CWE
                    new_comp = ntuple(i -> prev_comp[i] + chunk_tuple[i], q)
                    dest_dict[new_comp] = get(dest_dict, new_comp, BigInt(0)) + count
                end
            end
        end
        prev_layer = next_layer
        verbose && next!(p)
    end
   
    final_state = fill(T_zero, length(active_sets[end]))
    return get(prev_layer, final_state, Dict{NTuple{q, Int}, BigInt}())
end

function _CWE_classical_syndrome_binary(H::Matrix{T}, boundaries::Vector{Int}, verbose::Bool=true) where T
    r, n = size(H) 
    L, R = _get_LR_indices(H)
    _, _, active_sets, ending_sets = _build_trellis_sets(L, R, boundaries)

    num_sections = length(boundaries) - 1
    init_comp = (0, 0)
    
    prev_layer = Dict{UInt128, Dict{NTuple{2, Int}, BigInt}}(UInt128(0) => Dict(init_comp => BigInt(1)))
    p = verbose ? Progress(num_sections, 0.1, "Building Syndrome CWE Trellis (Binary): ") : nothing

    M_cols = zeros(UInt128, n)
    for col in 1:n
        col_val = UInt128(0)
        for row in 1:r
            if !iszero(H[row, col]) col_val |= (UInt128(1) << (row - 1)) end
        end
        M_cols[col] = col_val
    end

    for m in 1:num_sections
        next_layer = Dict{UInt128, Dict{NTuple{2, Int}, BigInt}}()
        
        active_curr = active_sets[m]
        sizehint!(next_layer, 1 << length(active_curr))
        
        left_b, right_b = boundaries[m], boundaries[m+1]
        chunk_len = right_b - left_b
        
        ending_curr = ending_sets[m]
        end_mask = UInt128(0)
        for row in ending_curr end_mask |= (UInt128(1) << (row - 1)) end
        
        keep_mask = UInt128(0)
        for row in active_curr keep_mask |= (UInt128(1) << (row - 1)) end
        
        # 2^chunk_len possible error configurations in this block
        for chunk_bits in 0:((1 << chunk_len) - 1)
            
            ones_wt = count_ones(chunk_bits)
            zeros_wt = chunk_len - ones_wt
            chunk_tuple = (zeros_wt, ones_wt)
            
            # Compute partial syndrome for this chunk
            partial_syn = UInt128(0)
            for i in 1:chunk_len
                if ((chunk_bits >> (i - 1)) & 1) == 1
                    partial_syn ⊻= M_cols[left_b + i]
                end
            end
            
            for (u_prev, partial_CWE) in prev_layer
                u_work = u_prev ⊻ partial_syn
                
                # Prune if the syndrome doesn't clear out on terminating rows
                if (u_work & end_mask) != 0 continue end
                
                u_next = u_work & keep_mask
                
                if !haskey(next_layer, u_next)
                    next_layer[u_next] = Dict{NTuple{2, Int}, BigInt}()
                end
                
                dest_dict = next_layer[u_next]
                for (prev_comp, count) in partial_CWE
                    new_comp = (prev_comp[1] + chunk_tuple[1], prev_comp[2] + chunk_tuple[2])
                    dest_dict[new_comp] = get(dest_dict, new_comp, BigInt(0)) + count
                end
            end
        end
        prev_layer = next_layer
        verbose && next!(p)
    end
   
    return get(prev_layer, UInt128(0), Dict{NTuple{2, Int}, BigInt}())
end

function _CWE_classical_syndrome_ternary(H::Matrix{T}, boundaries::Vector{Int}, verbose::Bool=true) where T
    r, n = size(H) 
    F = parent(H[1,1])
    elements = collect(F)
    elem_to_u8 = Dict(elements[1] => UInt8(0), elements[2] => UInt8(1), elements[3] => UInt8(2))
    
    H_u8 = zeros(UInt8, r, n)
    for i in 1:r, j in 1:n H_u8[i,j] = elem_to_u8[H[i,j]] end
    
    L, R = _get_LR_indices(H)
    _, _, active_sets, ending_sets = _build_trellis_sets(L, R, boundaries)

    num_sections = length(boundaries) - 1
    init_comp = (0, 0, 0)
    prev_layer = Dict{UInt128, Dict{NTuple{3, Int}, BigInt}}(UInt128(0) => Dict(init_comp => BigInt(1)))

    p = verbose ? Progress(num_sections, 0.1, "Building Syndrome CWE Trellis (GF(3)): ") : nothing
    u_buf = zeros(UInt8, r)

    for m in 1:num_sections
        next_layer = Dict{UInt128, Dict{NTuple{3, Int}, BigInt}}()
        
        active_curr = active_sets[m]
        sizehint!(next_layer, 3^length(active_curr))
        
        left_b, right_b = boundaries[m], boundaries[m+1]
        chunk_len = right_b - left_b
        active_prev = m == 1 ? Int[] : active_sets[m-1]
        ending_curr = ending_sets[m]

        for chunk_symbols in Iterators.product(fill((UInt8(0), UInt8(1), UInt8(2)), chunk_len)...)
            c0, c1, c2 = 0, 0, 0
            for sym in chunk_symbols
                if sym == 0 c0 += 1 elseif sym == 1 c1 += 1 else c2 += 1 end
            end
            chunk_tuple = (c0, c1, c2)
            
            chunk_vec = [sym for sym in chunk_symbols]
            partial_syn = zeros(UInt8, r)
            for col in 1:chunk_len
                if chunk_vec[col] != 0
                    for row in 1:r
                        partial_syn[row] = (partial_syn[row] + chunk_vec[col] * H_u8[row, left_b + col]) % 3
                    end
                end
            end
            
            for (u_packed, partial_CWE) in prev_layer
                fill!(u_buf, UInt8(0))
                for (idx, row) in enumerate(active_prev)
                    u_buf[row] = UInt8((u_packed >> (2 * (idx - 1))) & 3)
                end
                
                valid_branch = true
                for row in 1:r
                    u_buf[row] = (u_buf[row] + partial_syn[row]) % 3
                end
                for row in ending_curr
                    if u_buf[row] != 0
                        valid_branch = false
                        break
                    end
                end
                if !valid_branch continue end
                
                u_next = UInt128(0)
                for (idx, row) in enumerate(active_curr)
                    u_next |= (UInt128(u_buf[row]) << (2 * (idx - 1)))
                end
                
                if !haskey(next_layer, u_next)
                    next_layer[u_next] = Dict{NTuple{3, Int}, BigInt}()
                end
                dest_dict = next_layer[u_next]
                for (prev_comp, count) in partial_CWE
                    new_comp = (prev_comp[1] + chunk_tuple[1], prev_comp[2] + chunk_tuple[2], prev_comp[3] + chunk_tuple[3])
                    dest_dict[new_comp] = get(dest_dict, new_comp, BigInt(0)) + count
                end
            end
        end
        prev_layer = next_layer
        verbose && next!(p)
    end
    return get(prev_layer, UInt128(0), Dict{NTuple{3, Int}, BigInt}())
end

function _CWE_classical_syndrome_quaternary(H::Matrix{T}, boundaries::Vector{Int}, verbose::Bool=true) where T
    r, n = size(H) 
    F = parent(H[1,1])
    elements = collect(F)
    elem_to_u8 = Dict(elements[1] => UInt8(0), elements[2] => UInt8(1), elements[3] => UInt8(2), elements[4] => UInt8(3))
    GF4_MULT = UInt8[0 0 0 0; 0 1 2 3; 0 2 3 1; 0 3 1 2]
    
    H_u8 = zeros(UInt8, r, n)
    for i in 1:r, j in 1:n H_u8[i,j] = elem_to_u8[H[i,j]] end
    
    L, R = _get_LR_indices(H)
    _, _, active_sets, ending_sets = _build_trellis_sets(L, R, boundaries)

    num_sections = length(boundaries) - 1
    init_comp = (0, 0, 0, 0)
    prev_layer = Dict{UInt128, Dict{NTuple{4, Int}, BigInt}}(UInt128(0) => Dict(init_comp => BigInt(1)))

    p = verbose ? Progress(num_sections, 0.1, "Building Syndrome CWE Trellis (GF(4)): ") : nothing

    for m in 1:num_sections
        next_layer = Dict{UInt128, Dict{NTuple{4, Int}, BigInt}}()
        
        active_curr = active_sets[m]
        sizehint!(next_layer, 4^length(active_curr))
        
        left_b, right_b = boundaries[m], boundaries[m+1]
        chunk_len = right_b - left_b
        active_prev = m == 1 ? Int[] : active_sets[m-1]
        ending_curr = ending_sets[m]
        
        end_mask = UInt128(0)
        for row in ending_curr end_mask |= (UInt128(3) << (2 * (row - 1))) end
        
        keep_mask = UInt128(0)
        for (idx, row) in enumerate(active_curr) keep_mask |= (UInt128(3) << (2 * (idx - 1))) end

        for chunk_symbols in Iterators.product(fill((UInt8(0), UInt8(1), UInt8(2), UInt8(3)), chunk_len)...)
            c0, c1, c2, c3 = 0, 0, 0, 0
            for sym in chunk_symbols
                if sym == 0 c0 += 1 elseif sym == 1 c1 += 1 elseif sym == 2 c2 += 1 else c3 += 1 end
            end
            chunk_tuple = (c0, c1, c2, c3)
            
            chunk_vec = [sym for sym in chunk_symbols]
            partial_syn_packed = UInt128(0)
            
            for row in 1:r
                c_i = UInt8(0)
                for col in 1:chunk_len
                    if chunk_vec[col] != 0
                        c_i ⊻= GF4_MULT[chunk_vec[col] + 1, H_u8[row, left_b + col] + 1]
                    end
                end
                if c_i != 0
                    # Map the row's syndrome value into the correct position in the active_prev layout
                    idx = findfirst(x -> x == row, active_prev)
                    if !isnothing(idx)
                        partial_syn_packed |= (UInt128(c_i) << (2 * (idx - 1)))
                    else
                        # If it affects a row that wasn't active previously, it's either an ending row or an active_curr row
                        partial_syn_packed |= (UInt128(c_i) << (2 * (row - 1))) # We'll re-mask this properly below
                    end
                end
            end
            
            for (u_packed, partial_CWE) in prev_layer
                # GF(4) Native Hardware Addition!
                u_work = u_packed ⊻ partial_syn_packed
                
                if (u_work & end_mask) != 0 continue end
                
                # Compress u_work back down into the active_curr layout
                u_next = UInt128(0)
                for (idx, row) in enumerate(active_curr)
                    # Extract the value from its true row position
                    val = (u_work >> (2 * (row - 1))) & 3
                    # Pack it into its dense active_curr index position
                    u_next |= (UInt128(val) << (2 * (idx - 1)))
                end
                
                if !haskey(next_layer, u_next)
                    next_layer[u_next] = Dict{NTuple{4, Int}, BigInt}()
                end
                dest_dict = next_layer[u_next]
                for (prev_comp, count) in partial_CWE
                    new_comp = (prev_comp[1] + chunk_tuple[1], prev_comp[2] + chunk_tuple[2], prev_comp[3] + chunk_tuple[3], prev_comp[4] + chunk_tuple[4])
                    dest_dict[new_comp] = get(dest_dict, new_comp, BigInt(0)) + count
                end
            end
        end
        prev_layer = next_layer
        verbose && next!(p)
    end
    return get(prev_layer, UInt128(0), Dict{NTuple{4, Int}, BigInt}())
end

"""
    Krawtchouk(i::Int, j::Int, n::Int, q::Int)

Evaluates the Krawtchouk polynomial P_i(j; n, q).
"""
function Krawtchouk(i::Int, j::Int, n::Int, q::Int)
    val = BigInt(0)
    for r in 0:i
        if r <= j && (i - r) <= (n - j)
            term = (BigInt(-1)^r) * (BigInt(q - 1)^(i - r)) * binomial(BigInt(j), BigInt(r)) * binomial(BigInt(n - j), BigInt(i - r))
            val += term
        end
    end
    return val
end

"""
    MacWilliams_HWE_transform(input_hwe::Dict{Int, BigInt}, n::Int, k_in::Int, q::Int)

Applies the MacWilliams identity to convert a Hamming weight distribution 
into its dual Hamming weight distribution using Krawtchouk polynomials.
"""
function MacWilliams_HWE_transform(input_hwe::Dict{Int, BigInt}, n::Int, k_in::Int, q::Int)
    output_hwe = Dict{Int, BigInt}()
    scaling_factor = BigInt(q)^k_in
    
    for i in 0:n
        A_i = BigInt(0)
        for (j, A_in_j) in input_hwe
            A_i += A_in_j * Krawtchouk(i, j, n, q)
        end
        
        @assert A_i % scaling_factor == 0 "MacWilliams transform yielded non-integer. Check inputs."
        
        actual_A_i = A_i ÷ scaling_factor
        if actual_A_i > 0
            output_hwe[i] = actual_A_i
        end
    end
    
    return output_hwe
end

function _BZ_middle_search(M::Matrix{T}, L::Vector{Int}, R::Vector{Int}, B_L::Int, B_R::Int, 
                           left_dict, right_dict, q::Int, verbose::Bool=true) where T
    k, _ = size(M)
    if q == 2 && k <= 128
        return _BZ_middle_search_binary(M, L, R, B_L, B_R, left_dict, right_dict, verbose)
    elseif q == 3
        return _BZ_middle_search_ternary(M, L, R, B_L, B_R, left_dict, right_dict, verbose)
    elseif q == 4
        return _BZ_middle_search_quaternary(M, L, R, B_L, B_R, left_dict, right_dict, verbose)
    else
        return _BZ_middle_search_nonbinary(M, L, R, B_L, B_R, left_dict, right_dict, q, verbose)
    end
end

function _BZ_middle_search_binary(M::Matrix{T}, L::Vector{Int}, R::Vector{Int}, B_L::Int, B_R::Int, 
                                  left_dict::Dict{Tuple{UInt128, Bool}, Int}, right_dict::Dict{Tuple{UInt128, Bool}, Int}, 
                                  verbose::Bool=true) where T
    k, n = size(M)
    
    A_L = [row for row in 1:k if L[row] <= B_L && R[row] > B_L]
    A_R = [row for row in 1:k if L[row] <= B_R && R[row] > B_R]
    
    # The magical bit-mask to instantly extract the right state in O(1)
    A_R_mask = UInt128(0)
    for row in A_R A_R_mask |= (UInt128(1) << (row - 1)) end
    
    min_R = isempty(right_dict) ? 0 : minimum(values(right_dict))
    
    middle_cols = collect((B_L + 1):B_R)
    sort!(middle_cols, by = c -> count(!iszero, view(M, :, c)), rev = true)
    
    num_steps = length(middle_cols)
    M_cols = zeros(UInt128, num_steps)
    working_masks = zeros(UInt128, num_steps)
    branch_deltas_arr = Vector{Vector{UInt128}}(undef, num_steps)
    
    branched_rows = Set{Int}()
    
    verbose && println("  -> Packing $(length(middle_cols)) columns into hardware registers...")
    for (step, col) in enumerate(middle_cols)
        active_here = [row for row in 1:k if L[row] <= col <= R[row]]
        
        mask = UInt128(0)
        for row in active_here mask |= (UInt128(1) << (row - 1)) end
        working_masks[step] = mask
        
        col_val = UInt128(0)
        for row in 1:k
            if !iszero(M[row, col]) col_val |= (UInt128(1) << (row - 1)) end
        end
        M_cols[step] = col_val
        
        s_rows = Int[]
        for row in active_here
            if !(row in A_L) && !(row in branched_rows)
                push!(s_rows, row)
                push!(branched_rows, row)
            end
        end
        
        num_b = 1 << length(s_rows)
        deltas = zeros(UInt128, num_b)
        for b in 0:(num_b - 1)
            delta = UInt128(0)
            for (idx, row) in enumerate(s_rows)
                if ((b >> (idx - 1)) & 1) == 1
                    delta |= (UInt128(1) << (row - 1))
                end
            end
            deltas[b + 1] = delta 
        end
        branch_deltas_arr[step] = deltas
    end

    global_min_d = Ref(typemax(Int))
    min_d_lock = ReentrantLock()
    
    sorted_left = sort(collect(left_dict), by = x -> x[2])
    num_paths = length(sorted_left)
    
    paths_done = Threads.Atomic{Int}(0)
    print_threshold = max(1, div(num_paths, 20))
    verbose && println("  -> Bridging Middle Trellis ($num_paths left paths)...")
    
    Threads.@threads for i in 1:num_paths
        # No loop required here! The state is already a raw UInt128!
        ((u_local, l_nonzero), l_wt) = sorted_left[i]
        
        if l_wt + min_R >= global_min_d[] 
            Threads.atomic_add!(paths_done, 1)
            continue 
        end
        
        function dfs(step::Int, current_wt::Int, is_nonzero::Bool, u::UInt128)
            if current_wt + min_R >= global_min_d[] return end
            
            if step > num_steps
                # INSTANT O(1) HARDWARE LOOKUP
                u_right = u & A_R_mask
                
                key_f = (u_right, false)
                if haskey(right_dict, key_f)
                    total_wt = current_wt + right_dict[key_f]
                    if is_nonzero && total_wt < global_min_d[]
                        lock(min_d_lock) do
                            if total_wt < global_min_d[] global_min_d[] = total_wt end
                        end
                    end
                end
                
                key_t = (u_right, true)
                if haskey(right_dict, key_t)
                    total_wt = current_wt + right_dict[key_t]
                    if total_wt < global_min_d[]
                        lock(min_d_lock) do
                            if total_wt < global_min_d[] global_min_d[] = total_wt end
                        end
                    end
                end
                return
            end
            
            @inbounds for (b, delta) in enumerate(branch_deltas_arr[step])
                new_u = u | delta
                col_wt = count_ones(new_u & M_cols[step] & working_masks[step]) & 1
                new_nonzero = is_nonzero | (b > 1) 
                
                dfs(step + 1, current_wt + col_wt, new_nonzero, new_u)
            end
        end
        
        dfs(1, l_wt, l_nonzero, u_local)
        
        # Native Atomic Progress Tracking
        c = Threads.atomic_add!(paths_done, 1) + 1
        if verbose && (c % print_threshold == 0 || c == num_paths)
            pct = round(Int, (c / num_paths) * 100)
            println("     -> DFS Progress: $pct% ($c / $num_paths paths)")
        end
    end
    
    return global_min_d[]
end
function _BZ_middle_search_ternary(M::Matrix{T}, L::Vector{Int}, R::Vector{Int}, B_L::Int, B_R::Int, 
                                   left_dict::Dict{Tuple{Vector{T}, Bool}, Int}, right_dict::Dict{Tuple{Vector{T}, Bool}, Int}, 
                                   verbose::Bool=true) where T
    k, n = size(M)
    F = parent(M[1,1])
    elements = collect(F)
    
    elem_to_u8 = Dict(elements[1] => UInt8(0), elements[2] => UInt8(1), elements[3] => UInt8(2))
    u8_to_elem = Dict(UInt8(0) => elements[1], UInt8(1) => elements[2], UInt8(2) => elements[3])
    
    M_u8 = zeros(UInt8, k, n)
    for i in 1:k, j in 1:n M_u8[i,j] = elem_to_u8[M[i,j]] end
    
    A_L = [row for row in 1:k if L[row] <= B_L && R[row] > B_L]
    A_R = [row for row in 1:k if L[row] <= B_R && R[row] > B_R]
    
    min_R = isempty(right_dict) ? 0 : minimum(values(right_dict))
    
    middle_cols = collect((B_L + 1):B_R)
    sort!(middle_cols, by = c -> count(!iszero, view(M_u8, :, c)), rev = true)
    
    starts_at = [Int[] for _ in 1:length(middle_cols)]
    active_at_step = [Int[] for _ in 1:length(middle_cols)]
    branched_rows = Set{Int}()
    
    verbose && println("  -> Unboxing $(length(middle_cols)) middle columns to native UInt8 for GF(3)...")
    for (step, col) in enumerate(middle_cols)
        active_here = [row for row in 1:k if L[row] <= col <= R[row]]
        active_at_step[step] = active_here
        for row in active_here
            if !(row in A_L) && !(row in branched_rows)
                push!(starts_at[step], row)
                push!(branched_rows, row)
            end
        end
    end

    global_min_d = Ref(typemax(Int))
    min_d_lock = ReentrantLock()
    sorted_left = sort(collect(left_dict), by = x -> x[2])
    num_paths = length(sorted_left)
    
    p = Progress(num_paths, 1, "Bridging Middle Trellis (GF(3)): ")
    
    Threads.@threads for i in 1:num_paths
        ((l_state, l_nonzero), l_wt) = sorted_left[i]
        
        if l_wt + min_R >= global_min_d[] 
            next!(p)
            continue 
        end
        
        u_local = zeros(UInt8, k)
        @inbounds for (idx, row) in enumerate(A_L) u_local[row] = elem_to_u8[l_state[idx]] end
        
        right_state_buffer = fill(elements[1], length(A_R))
        
        function dfs(step::Int, current_wt::Int, is_nonzero::Bool, u::Vector{UInt8})
            if current_wt + min_R >= global_min_d[] return end
            
            if step > length(middle_cols)
                @inbounds for (idx, row) in enumerate(A_R) right_state_buffer[idx] = u8_to_elem[u[row]] end
                
                key_f = (right_state_buffer, false)
                if haskey(right_dict, key_f)
                    total_wt = current_wt + right_dict[key_f]
                    if is_nonzero && total_wt < global_min_d[]
                        lock(min_d_lock) do
                            if total_wt < global_min_d[] global_min_d[] = total_wt end
                        end
                    end
                end
                
                key_t = (right_state_buffer, true)
                if haskey(right_dict, key_t)
                    total_wt = current_wt + right_dict[key_t]
                    if total_wt < global_min_d[]
                        lock(min_d_lock) do
                            if total_wt < global_min_d[] global_min_d[] = total_wt end
                        end
                    end
                end
                return
            end
            
            col = middle_cols[step]
            starting_rows = starts_at[step]
            working_rows = active_at_step[step]
            
            @inbounds for branch_scalars in Iterators.product(fill((UInt8(0), UInt8(1), UInt8(2)), length(starting_rows))...)
                for (idx, row) in enumerate(starting_rows) u[row] = branch_scalars[idx] end
                
                c_i = UInt32(0) 
                for row in working_rows
                    c_i += u[row] * M_u8[row, col]
                end
                col_wt = (c_i % 3) == 0 ? 0 : 1
                
                new_nonzero = is_nonzero
                for val in branch_scalars
                    if val != 0 new_nonzero = true end
                end
                
                dfs(step + 1, current_wt + col_wt, new_nonzero, u)
            end
        end
        
        dfs(1, l_wt, l_nonzero, u_local)
        next!(p)
    end
    
    return global_min_d[]
end

function _BZ_middle_search_quaternary(M::Matrix{T}, L::Vector{Int}, R::Vector{Int}, B_L::Int, B_R::Int, 
                                      left_dict::Dict{Tuple{Vector{T}, Bool}, Int}, right_dict::Dict{Tuple{Vector{T}, Bool}, Int}, 
                                      verbose::Bool=true) where T
    k, n = size(M)
    F = parent(M[1,1])
    elements = collect(F)
    
    elem_to_u8 = Dict(elements[1] => UInt8(0), elements[2] => UInt8(1), elements[3] => UInt8(2), elements[4] => UInt8(3))
    u8_to_elem = Dict(UInt8(0) => elements[1], UInt8(1) => elements[2], UInt8(2) => elements[3], UInt8(3) => elements[4])
    
    GF4_MULT = UInt8[
        0 0 0 0;
        0 1 2 3;
        0 2 3 1;
        0 3 1 2
    ]
    
    M_u8 = zeros(UInt8, k, n)
    for i in 1:k, j in 1:n M_u8[i,j] = elem_to_u8[M[i,j]] end
    
    A_L = [row for row in 1:k if L[row] <= B_L && R[row] > B_L]
    A_R = [row for row in 1:k if L[row] <= B_R && R[row] > B_R]
    
    min_R = isempty(right_dict) ? 0 : minimum(values(right_dict))
    
    middle_cols = collect((B_L + 1):B_R)
    sort!(middle_cols, by = c -> count(!iszero, view(M_u8, :, c)), rev = true)
    
    starts_at = [Int[] for _ in 1:length(middle_cols)]
    active_at_step = [Int[] for _ in 1:length(middle_cols)]
    branched_rows = Set{Int}()
    
    verbose && println("  -> Unboxing $(length(middle_cols)) middle columns to native UInt8 for GF(4)...")
    for (step, col) in enumerate(middle_cols)
        active_here = [row for row in 1:k if L[row] <= col <= R[row]]
        active_at_step[step] = active_here
        for row in active_here
            if !(row in A_L) && !(row in branched_rows)
                push!(starts_at[step], row)
                push!(branched_rows, row)
            end
        end
    end

    global_min_d = Ref(typemax(Int))
    min_d_lock = ReentrantLock()
    sorted_left = sort(collect(left_dict), by = x -> x[2])
    num_paths = length(sorted_left)
    
    p = Progress(num_paths, 1, "Bridging Middle Trellis (GF(4)): ")
    
    Threads.@threads for i in 1:num_paths
        ((l_state, l_nonzero), l_wt) = sorted_left[i]
        
        if l_wt + min_R >= global_min_d[] 
            next!(p)
            continue 
        end
        
        u_local = zeros(UInt8, k)
        @inbounds for (idx, row) in enumerate(A_L) u_local[row] = elem_to_u8[l_state[idx]] end
        
        right_state_buffer = fill(elements[1], length(A_R))
        
        function dfs(step::Int, current_wt::Int, is_nonzero::Bool, u::Vector{UInt8})
            if current_wt + min_R >= global_min_d[] return end
            
            if step > length(middle_cols)
                @inbounds for (idx, row) in enumerate(A_R) right_state_buffer[idx] = u8_to_elem[u[row]] end
                
                key_f = (right_state_buffer, false)
                if haskey(right_dict, key_f)
                    total_wt = current_wt + right_dict[key_f]
                    if is_nonzero && total_wt < global_min_d[]
                        lock(min_d_lock) do
                            if total_wt < global_min_d[] global_min_d[] = total_wt end
                        end
                    end
                end
                
                key_t = (right_state_buffer, true)
                if haskey(right_dict, key_t)
                    total_wt = current_wt + right_dict[key_t]
                    if total_wt < global_min_d[]
                        lock(min_d_lock) do
                            if total_wt < global_min_d[] global_min_d[] = total_wt end
                        end
                    end
                end
                return
            end
            
            col = middle_cols[step]
            starting_rows = starts_at[step]
            working_rows = active_at_step[step]
            
            @inbounds for branch_scalars in Iterators.product(fill((UInt8(0), UInt8(1), UInt8(2), UInt8(3)), length(starting_rows))...)
                for (idx, row) in enumerate(starting_rows) u[row] = branch_scalars[idx] end
                
                c_i = UInt8(0)
                for row in working_rows
                    c_i ⊻= GF4_MULT[u[row] + 1, M_u8[row, col] + 1] 
                end
                col_wt = c_i == 0 ? 0 : 1
                
                new_nonzero = is_nonzero
                for val in branch_scalars
                    if val != 0 new_nonzero = true end
                end
                
                dfs(step + 1, current_wt + col_wt, new_nonzero, u)
            end
        end
        
        dfs(1, l_wt, l_nonzero, u_local)
        next!(p)
    end
    
    return global_min_d[]
end

function _BZ_middle_search_nonbinary(M::Matrix{T}, L::Vector{Int}, R::Vector{Int}, B_L::Int, B_R::Int, left_dict, right_dict, q::Int) where T
    k, n = size(M)
    F = parent(M[1,1])
    elements = collect(F)
    
    A_L = [row for row in 1:k if L[row] <= B_L && R[row] > B_L]
    A_R = [row for row in 1:k if L[row] <= B_R && R[row] > B_R]
    
    min_R = isempty(right_dict) ? 0 : minimum(values(right_dict))
    
    middle_cols = collect((B_L + 1):B_R)
    sort!(middle_cols, by = c -> count(!iszero, view(M, :, c)), rev = true)
    
    starts_at = [Int[] for _ in 1:length(middle_cols)]
    active_at_step = [Int[] for _ in 1:length(middle_cols)]
    branched_rows = Set{Int}()
    
    for (step, col) in enumerate(middle_cols)
        active_here = [row for row in 1:k if L[row] <= col <= R[row]]
        active_at_step[step] = active_here
        for row in active_here
            if !(row in A_L) && !(row in branched_rows)
                push!(starts_at[step], row)
                push!(branched_rows, row)
            end
        end
    end

    global_min_d = Ref(typemax(Int))
    min_d_lock = ReentrantLock()
    sorted_left = sort(collect(left_dict), by = x -> x[2])
    num_paths = length(sorted_left)
    
    # Initialize the thread-safe progress meter
    p = Progress(num_paths, 1, "Bridging Middle Trellis: ")
    
    Threads.@threads for i in 1:num_paths
        ((l_state, l_nonzero), l_wt) = sorted_left[i]
        
        if l_wt + min_R >= global_min_d[] continue end
        
        u_local = fill(zero(F), k)
        @inbounds for (idx, row) in enumerate(A_L) u_local[row] = l_state[idx] end
        
        function dfs(step::Int, current_wt::Int, is_nonzero::Bool, u::Vector{T})
            if current_wt + min_R >= global_min_d[] return end
            
            if step > length(middle_cols)
                right_state = [u[row] for row in A_R]
                
                key_f = (right_state, false)
                if haskey(right_dict, key_f)
                    total_wt = current_wt + right_dict[key_f]
                    if is_nonzero && total_wt < global_min_d[]
                        lock(min_d_lock) do
                            if total_wt < global_min_d[] global_min_d[] = total_wt end
                        end
                    end
                end
                
                key_t = (right_state, true)
                if haskey(right_dict, key_t)
                    total_wt = current_wt + right_dict[key_t]
                    if total_wt < global_min_d[]
                        lock(min_d_lock) do
                            if total_wt < global_min_d[] global_min_d[] = total_wt end
                        end
                    end
                end
                return
            end
            
            col = middle_cols[step]
            starting_rows = starts_at[step]
            working_rows = active_at_step[step]
            
            @inbounds for branch_scalars in Iterators.product(fill(elements, length(starting_rows))...)
                for (idx, row) in enumerate(starting_rows) u[row] = branch_scalars[idx] end
                
                c_i = zero(F)
                for row in working_rows
                    c_i += u[row] * M[row, col]
                end
                col_wt = iszero(c_i) ? 0 : 1
                
                new_nonzero = is_nonzero
                for val in branch_scalars
                    if !iszero(val) new_nonzero = true end
                end
                
                dfs(step + 1, current_wt + col_wt, new_nonzero, u)
            end
        end
        
        dfs(1, l_wt, l_nonzero, u_local)
        # The DFS has completely finished for this left_path.
        # Atomically increment the progress bar.
        next!(p)
    end
    
    return global_min_d[]
end

"""
    _optimal_sectionalization_QC(M::Matrix{T}, q::Int, p::Int) where T

Finds the optimal periodic sectionalization for a Quasi-Cyclic code with block size `p`.
Evaluates the global cost of cuts within a single block of size `p` and tiles 
the optimal boundaries across the entire code of length `n`.

Drastically reduces DAG complexity from O(n^2) to O(p^2).
"""
function _optimal_sectionalization_QC(M::Matrix{T}, q::Int, p::Int) where T
    k, n = size(M)
    
    # Ensure the matrix length is a perfect multiple of the block size
    @assert n % p == 0 "Code length n ($n) must be divisible by block size p ($p)."
    m = div(n, p) # Number of repeating blocks
    
    # 1. Precompute global past/future profiles
    L, R = _get_LR_indices(M)
    past, future = past_future_profiles(L, R, n)
    
    # 2. Setup the mini-DAG for a single block (vertices 0 to p)
    min_cost = fill(Inf, p + 1)
    min_cost[1] = 0.0
    parent = zeros(Int, p + 1)
    
    # 3. O(p^2) DAG Traversal (with global block averaging)
    for v in 1:p
        for u in 0:(v - 1)
            
            total_periodic_cost = 0.0
            
            # Evaluate this relative cut across ALL m blocks in the code
            for j in 0:(m - 1)
                u_global = u + j * p
                v_global = v + j * p
                
                # Our beautiful O(1) intersection math
                active_count = k - past[u_global + 1] - future[v_global + 1]
                
                # Safety cutoff for massive overlaps
                if active_count > 50
                    total_periodic_cost = Inf
                    break
                end
                
                total_periodic_cost += Float64(q)^active_count
            end
            
            # Relax the edge on the mini-DAG
            if min_cost[u + 1] + total_periodic_cost < min_cost[v + 1]
                min_cost[v + 1] = min_cost[u + 1] + total_periodic_cost
                parent[v + 1] = u
            end
        end
    end
    
    # 4. Reconstruct the base path for the single block
    base_path = Int[]
    curr = p
    while curr != 0
        push!(base_path, curr)
        curr = parent[curr + 1]
    end
    push!(base_path, 0)
    reverse!(base_path)
    
    # 5. Tile the base path across the entire matrix
    boundaries = Int[]
    for j in 0:(m - 1)
        # We drop the final element of base_path to avoid duplicating block boundaries 
        # (e.g., end of block 1 is start of block 2)
        for b in base_path[1:end-1]
            push!(boundaries, b + j * p)
        end
    end
    push!(boundaries, n) # Cap off the final boundary
    
    return boundaries
end

"""
    _optimal_sectionalization_2d_cyclic(M::Matrix{T}, q::Int, p_x::Int, p_y::Int, grid_x::Int, grid_y::Int) where T

Finds the optimal periodic sectionalization for a 2D Cyclic Code mapped to a 1D matrix.
Evaluates cuts inside a single 1D unit cell of size `p_x` and tiles them across 
both the x and y dimensions of the grid.
"""
function _optimal_sectionalization_2d_cyclic(M::Matrix{T}, q::Int, p_x::Int, p_y::Int, grid_x::Int, grid_y::Int) where T
    k, n = size(M)
    @assert n == grid_x * grid_y "Matrix length n must equal grid_x * grid_y"
    
    L, R = _get_LR_indices(M)
    past, future = past_future_profiles(L, R, n)
    
    # Setup mini-DAG for just the inner unit cell (size p_x)
    min_cost = fill(Inf, p_x + 1)
    min_cost[1] = 0.0
    parent = zeros(Int, p_x + 1)
    
    # Number of tiles in X and Y
    tiles_x = div(grid_x, p_x)
    tiles_y = div(grid_y, p_y)
    
    for v in 1:p_x
        for u in 0:(v - 1)
            
            total_2d_cost = 0.0
            
            # Evaluate this cut across the entire 2D lattice!
            for y in 0:(tiles_y - 1)
                for x in 0:(tiles_x - 1)
                    # Map the 2D block coordinate back to the 1D matrix index
                    offset = y * (p_y * grid_x) + x * p_x
                    u_global = u + offset
                    v_global = v + offset
                    
                    active_count = k - past[u_global + 1] - future[v_global + 1]
                    
                    if active_count > 50
                        total_2d_cost = Inf
                        break
                    end
                    total_2d_cost += Float64(q)^active_count
                end
                if total_2d_cost == Inf break end
            end
            
            # Relax the edge on the unit-cell DAG
            if min_cost[u + 1] + total_2d_cost < min_cost[v + 1]
                min_cost[v + 1] = min_cost[u + 1] + total_2d_cost
                parent[v + 1] = u
            end
        end
    end
    
    # Reconstruct the base path for the unit cell
    base_path = Int[]
    curr = p_x
    while curr != 0
        push!(base_path, curr)
        curr = parent[curr + 1]
    end
    push!(base_path, 0)
    reverse!(base_path)
    
    # TILE ACROSS THE 1D FLATTENED MATRIX
    boundaries = Int[]
    # We tile the base pattern across every horizontal block
    for i in 0:(div(n, p_x) - 1)
        for b in base_path[1:end-1]
            push!(boundaries, b + i * p_x)
        end
    end
    push!(boundaries, n)
    
    return boundaries
end

"""
    optimal_sectionalization(M::Matrix{T}, q::Int; type::Symbol=:linear, max_width::Int=10, kwargs...) where T

Computes the optimal sectionalization bounds for the Viterbi trellis to minimize 
peak state and branch complexity. Enforces a strict `max_width` to prevent 
exponential branch evaluation traps in low-density or syndrome matrices.

# Types
* `:linear` - Standard optimal sectionalization for a generic code.
* `:QC` - Imposes periodicity for a Quasi-Cyclic code. Requires kwarg `p` (block size).
* `:twoD` - Evaluates a 2D grid mapped to 1D. Requires kwargs `p_x`, `p_y`, `grid_x`, `grid_y`.
"""
function optimal_sectionalization(M::Matrix{T}, q::Int; 
                                  type::Symbol=:linear, 
                                  max_width::Int=10, 
                                  p::Int=0, 
                                  p_x::Int=0, p_y::Int=0, 
                                  grid_x::Int=0, grid_y::Int=0) where T
    
    if type == :linear
        k, n = size(M)
        L, R = CodingTheory._get_LR_indices(M)
        
        # DP array: min max-complexity to reach column j
        cost = fill(typemax(Int), n + 1)
        parent = zeros(Int, n + 1)
        cost[1] = 0 # Base case: cost to reach column 0 is 0
        
        for j in 1:n
            # THE CAP: Look backwards to find the best starting column i, 
            # but NEVER look further back than max_width.
            min_i = max(0, j - max_width)
            
            for i in min_i:(j - 1)
                active_rows = 0
                for row in 1:k
                    if L[row] <= j && R[row] > i
                        active_rows += 1
                    end
                end
                
                # The section cost is dominated by q^(active_rows)
                # We just track the exponent to prevent integer overflow in the DP
                section_cost = active_rows
                path_cost = max(cost[i + 1], section_cost)
                
                if path_cost < cost[j + 1]
                    cost[j + 1] = path_cost
                    parent[j + 1] = i
                end
            end
        end
        
        # Backtrack to extract the optimal section boundaries
        bounds = Int[]
        curr = n
        while curr > 0
            push!(bounds, curr)
            curr = parent[curr + 1]
        end
        push!(bounds, 0)
        reverse!(bounds)
        
        return bounds

    elseif type == :QC
        p > 0 || throw(AssertionError("Block size `p` must be provided for Quasi-Cyclic sectionalization."))
        return _optimal_sectionalization_QC(M, q, p)

    elseif type == :twoD
        (p_x > 0 && p_y > 0 && grid_x > 0 && grid_y > 0) || throw(AssertionError("Parameters `p_x`, `p_y`, `grid_x`, `grid_y` must be provided for 2D sectionalization."))
        return _optimal_sectionalization_2d_cyclic(M, q, p_x, p_y, grid_x, grid_y)

    else
        throw(ArgumentError("Unknown sectionalization type: $type"))
    end
end

optimal_sectionalization(M::CTMatrixTypes, q::Int; kwargs...) = optimal_sectionalization(Array(M), q; kwargs...)

function _forward_trellis(M::Matrix{T}, B_L::Int, q::Int, block_size::Int=0, verbose::Bool=true) where T
    k, _ = size(M)
    if q == 2 && k <= 128
        return _forward_trellis_binary(M, B_L, block_size, verbose)
    elseif q == 3 && k <= 64
        return _forward_trellis_ternary(M, B_L, block_size, verbose)
    elseif q == 4 && k <= 64
        return _forward_trellis_quaternary(M, B_L, block_size, verbose)
    else
        return _forward_trellis_nonbinary(M, B_L, q, block_size, verbose)
    end
end

function _forward_trellis_binary(M::Matrix{T}, B_L::Int, block_size::Int=0, verbose::Bool=true) where T
    k, n = size(M)
    
    L, R = _get_LR_indices(M)
    prev_layer = Dict{Tuple{UInt128, Bool}, Int}((UInt128(0), false) => 0)
    
    active_sets = [Int[] for _ in 1:B_L]
    start_sets  = [Int[] for _ in 1:B_L]
    work_sets   = [Int[] for _ in 1:B_L]
    M_cols      = zeros(UInt128, B_L)
    work_masks  = zeros(UInt128, B_L)

    for col in 1:B_L
        mask = UInt128(0)
        col_val = UInt128(0)
        for row in 1:k
            if L[row] <= col < R[row] push!(active_sets[col], row) end
            if L[row] == col push!(start_sets[col], row) end
            if L[row] <= col <= R[row]
                push!(work_sets[col], row)
                mask |= (UInt128(1) << (row - 1))
            end
            if !iszero(M[row, col]) col_val |= (UInt128(1) << (row - 1)) end
        end
        work_masks[col] = mask
        M_cols[col] = col_val
    end

    p = verbose ? Progress(B_L, 0.1, "Building Forward Trellis: ") : nothing
    
    for col in 1:B_L
        next_layer = Dict{Tuple{UInt128, Bool}, Int}()
        
        active_curr = active_sets[col]
        sizehint!(next_layer, 1 << length(active_curr))
        
        starting = start_sets[col]
        num_branches = 1 << length(starting)
        deltas = zeros(UInt128, num_branches)
        for b in 0:(num_branches-1)
            d = UInt128(0)
            for (idx, row) in enumerate(starting)
                if ((b >> (idx - 1)) & 1) == 1 d |= (UInt128(1) << (row - 1)) end
            end
            deltas[b+1] = d
        end
        
        keep_mask = UInt128(0)
        for row in active_curr keep_mask |= (UInt128(1) << (row - 1)) end
        
        col_val = M_cols[col]
        w_mask = work_masks[col]
        
        for ((u_prev, is_nonzero), prev_wt) in prev_layer
            @inbounds for b in 1:num_branches
                delta = deltas[b]
                u_work = u_prev | delta
                
                col_wt = count_ones(u_work & col_val & w_mask) & 1
                new_wt = prev_wt + col_wt
                
                new_is_nonzero = is_nonzero | (b > 1)
                u_next = u_work & keep_mask 
                
                key = (u_next, new_is_nonzero)
                if !haskey(next_layer, key) || new_wt < next_layer[key]
                    next_layer[key] = new_wt
                end
            end
        end
        
        if col == block_size
            zero_key = (UInt128(0), false)
            if haskey(next_layer, zero_key) && next_layer[zero_key] == 0
                delete!(next_layer, zero_key)
            end
        end
        
        prev_layer = next_layer
        verbose && next!(p)
    end
    
    # DO NOT UNPACK! Return raw hardware registers.
    return prev_layer
end

function _forward_trellis_ternary(M::Matrix{T}, B_L::Int, block_size::Int=0, verbose::Bool=true) where T
    k, n = size(M)
    F = parent(M[1,1])
    elements = collect(F)
    T_zero = zero(F)
    
    elem_to_u8 = Dict(elements[1] => UInt8(0), elements[2] => UInt8(1), elements[3] => UInt8(2))
    u8_to_elem = Dict(UInt8(0) => elements[1], UInt8(1) => elements[2], UInt8(2) => elements[3])
    
    M_u8 = zeros(UInt8, k, n)
    for i in 1:k, j in 1:n M_u8[i,j] = elem_to_u8[M[i,j]] end
    
    L, R = _get_LR_indices(M)
    
    prev_layer = Dict{Tuple{UInt128, Bool}, Int}((UInt128(0), false) => 0)
    
    active_sets = [Int[] for _ in 1:B_L]
    start_sets  = [Int[] for _ in 1:B_L]
    work_sets   = [Int[] for _ in 1:B_L]

    for col in 1:B_L
        for row in 1:k
            if L[row] <= col < R[row] push!(active_sets[col], row) end
            if L[row] == col push!(start_sets[col], row) end
            if L[row] <= col <= R[row] push!(work_sets[col], row) end
        end
    end

    p = verbose ? Progress(B_L, 0.1, "Building Forward Trellis (GF(3)): ") : nothing
    u_buf = zeros(UInt8, k)
    
    for col in 1:B_L
        next_layer = Dict{Tuple{UInt128, Bool}, Int}()
        active_curr = active_sets[col]
        sizehint!(next_layer, 3^length(active_curr))
        
        active_prev = col == 1 ? Int[] : active_sets[col-1]
        starting = start_sets[col]
        working = work_sets[col]
        
        for ((u_packed, is_nonzero), prev_wt) in prev_layer
            
            fill!(u_buf, UInt8(0))
            # Shift hardware bits to unpack local states
            for (idx, row) in enumerate(active_prev)
                u_buf[row] = UInt8((u_packed >> (2 * (idx - 1))) & 3)
            end
            
            for branch_scalars in Iterators.product(fill((UInt8(0), UInt8(1), UInt8(2)), length(starting))...)
                @inbounds for (idx, row) in enumerate(starting) u_buf[row] = branch_scalars[idx] end
                
                c_i = UInt32(0)
                @inbounds for row in working
                    c_i += u_buf[row] * M_u8[row, col]
                end
                col_wt = (c_i % 3) == 0 ? 0 : 1
                new_wt = prev_wt + col_wt
                
                new_is_nonzero = is_nonzero
                for val in branch_scalars
                    if val != 0 new_is_nonzero = true end
                end
                
                # Shift bits to pack hardware dictionary key
                next_packed = UInt128(0)
                @inbounds for (idx, row) in enumerate(active_curr)
                    next_packed |= (UInt128(u_buf[row]) << (2 * (idx - 1)))
                end
                
                key = (next_packed, new_is_nonzero)
                if !haskey(next_layer, key) || new_wt < next_layer[key]
                    next_layer[key] = new_wt
                end
            end
        end
        
        if col == block_size
            zero_key = (UInt128(0), false)
            if haskey(next_layer, zero_key) && next_layer[zero_key] == 0
                delete!(next_layer, zero_key)
            end
        end
        
        prev_layer = next_layer
        verbose && next!(p)
    end
    
    # Translate hardware registers back into Oscar Vectors for the Middle Search
    final_dict = Dict{Tuple{Vector{T}, Bool}, Int}()
    active_final = active_sets[B_L]
    for ((u_packed, is_nz), wt) in prev_layer
        u_unpacked = fill(T_zero, length(active_final))
        for (idx, row) in enumerate(active_final)
            val_u8 = UInt8((u_packed >> (2 * (idx - 1))) & 3)
            u_unpacked[idx] = u8_to_elem[val_u8]
        end
        final_dict[(u_unpacked, is_nz)] = wt
    end
    
    return final_dict
end

function _forward_trellis_quaternary(M::Matrix{T}, B_L::Int, block_size::Int=0, verbose::Bool=true) where T
    k, n = size(M)
    F = parent(M[1,1])
    elements = collect(F)
    T_zero = zero(F)
    
    elem_to_u8 = Dict(elements[1] => UInt8(0), elements[2] => UInt8(1), elements[3] => UInt8(2), elements[4] => UInt8(3))
    u8_to_elem = Dict(UInt8(0) => elements[1], UInt8(1) => elements[2], UInt8(2) => elements[3], UInt8(3) => elements[4])
    
    GF4_MULT = UInt8[
        0 0 0 0;
        0 1 2 3;
        0 2 3 1;
        0 3 1 2
    ]
    
    M_u8 = zeros(UInt8, k, n)
    for i in 1:k, j in 1:n M_u8[i,j] = elem_to_u8[M[i,j]] end
    
    L, R = _get_LR_indices(M)
    
    prev_layer = Dict{Tuple{UInt128, Bool}, Int}((UInt128(0), false) => 0)
    
    active_sets = [Int[] for _ in 1:B_L]
    start_sets  = [Int[] for _ in 1:B_L]
    work_sets   = [Int[] for _ in 1:B_L]

    for col in 1:B_L
        for row in 1:k
            if L[row] <= col < R[row] push!(active_sets[col], row) end
            if L[row] == col push!(start_sets[col], row) end
            if L[row] <= col <= R[row] push!(work_sets[col], row) end
        end
    end

    p = verbose ? Progress(B_L, 0.1, "Building Forward Trellis (GF(4)): ") : nothing
    u_buf = zeros(UInt8, k)
    
    for col in 1:B_L
        next_layer = Dict{Tuple{UInt128, Bool}, Int}()
        
        active_curr = active_sets[col]
        sizehint!(next_layer, 4^length(active_curr))
        
        active_prev = col == 1 ? Int[] : active_sets[col-1]
        starting = start_sets[col]
        working = work_sets[col]
        
        for ((u_packed, is_nonzero), prev_wt) in prev_layer
            
            fill!(u_buf, UInt8(0))
            for (idx, row) in enumerate(active_prev)
                u_buf[row] = UInt8((u_packed >> (2 * (idx - 1))) & 3)
            end
            
            for branch_scalars in Iterators.product(fill((UInt8(0), UInt8(1), UInt8(2), UInt8(3)), length(starting))...)
                @inbounds for (idx, row) in enumerate(starting) u_buf[row] = branch_scalars[idx] end
                
                c_i = UInt8(0)
                @inbounds for row in working
                    c_i ⊻= GF4_MULT[u_buf[row] + 1, M_u8[row, col] + 1]
                end
                col_wt = c_i == 0 ? 0 : 1
                new_wt = prev_wt + col_wt
                
                new_is_nonzero = is_nonzero
                for val in branch_scalars
                    if val != 0 new_is_nonzero = true end
                end
                
                next_packed = UInt128(0)
                @inbounds for (idx, row) in enumerate(active_curr)
                    next_packed |= (UInt128(u_buf[row]) << (2 * (idx - 1)))
                end
                
                key = (next_packed, new_is_nonzero)
                if !haskey(next_layer, key) || new_wt < next_layer[key]
                    next_layer[key] = new_wt
                end
            end
        end
        
        if col == block_size
            zero_key = (UInt128(0), false)
            if haskey(next_layer, zero_key) && next_layer[zero_key] == 0
                delete!(next_layer, zero_key)
            end
        end
        
        prev_layer = next_layer
        verbose && next!(p)
    end
    
    final_dict = Dict{Tuple{Vector{T}, Bool}, Int}()
    active_final = active_sets[B_L]
    for ((u_packed, is_nz), wt) in prev_layer
        u_unpacked = fill(T_zero, length(active_final))
        for (idx, row) in enumerate(active_final)
            val_u8 = UInt8((u_packed >> (2 * (idx - 1))) & 3)
            u_unpacked[idx] = u8_to_elem[val_u8]
        end
        final_dict[(u_unpacked, is_nz)] = wt
    end
    
    return final_dict
end

function _forward_trellis_nonbinary(M::Matrix{T}, B_L::Int, q::Int, block_size::Int=0, verbose::Bool=true) where T
    k, n = size(M)
    F = parent(M[1, 1])
    elements = collect(F)
    T_zero = zero(F)
    
    L, R = _get_LR_indices(M)
    prev_layer = Dict{Tuple{Vector{T}, Bool}, Int}((Vector{T}(), false) => 0)
    
    active_sets = [Int[] for _ in 1:B_L]
    start_sets  = [Int[] for _ in 1:B_L]
    work_sets   = [Int[] for _ in 1:B_L]

    for col in 1:B_L
        for row in 1:k
            if L[row] <= col < R[row] push!(active_sets[col], row) end
            if L[row] == col push!(start_sets[col], row) end
            if L[row] <= col <= R[row] push!(work_sets[col], row) end
        end
    end

    p = verbose ? Progress(B_L, 1, "Building Forward Trellis (Non-Binary): ") : nothing
    u_buf = fill(T_zero, k)
    
    for col in 1:B_L
        next_layer = Dict{Tuple{Vector{T}, Bool}, Int}()
        
        active_curr = active_sets[col]
        sizehint!(next_layer, q^length(active_curr))
        
        active_prev = col == 1 ? Int[] : active_sets[col-1]
        starting = start_sets[col]
        working = work_sets[col]
        
        for ((prev_scalars, is_nonzero), prev_wt) in prev_layer
            for branch_scalars in Iterators.product(fill(elements, length(starting))...)
                
                fill!(u_buf, T_zero)
                @inbounds for (idx, row) in enumerate(active_prev) u_buf[row] = prev_scalars[idx] end
                @inbounds for (idx, row) in enumerate(starting)    u_buf[row] = branch_scalars[idx] end
                
                c_i = T_zero
                @inbounds for row in working
                    c_i += u_buf[row] * M[row, col]
                end
                col_wt = iszero(c_i) ? 0 : 1
                new_wt = prev_wt + col_wt
                
                new_is_nonzero = is_nonzero || any(!iszero, branch_scalars)
                next_scalars = [u_buf[row] for row in active_curr]
                
                key = (next_scalars, new_is_nonzero)
                if !haskey(next_layer, key) || new_wt < next_layer[key]
                    next_layer[key] = new_wt
                end
            end
        end
        
        if col == block_size
            zero_key = (fill(T_zero, length(active_curr)), false)
            if haskey(next_layer, zero_key) && next_layer[zero_key] == 0
                delete!(next_layer, zero_key)
            end
        end
        
        prev_layer = next_layer
        verbose && next!(p)
    end
    
    return prev_layer
end

function _backward_trellis(M::Matrix{T}, n::Int, B_R::Int, q::Int, verbose::Bool=true) where T
    k, _ = size(M)
    if q == 2 && k <= 128
        return _backward_trellis_binary(M, n, B_R, verbose)
    elseif q == 3 && k <= 64
        return _backward_trellis_ternary(M, n, B_R, verbose)
    elseif q == 4 && k <= 64
        return _backward_trellis_quaternary(M, n, B_R, verbose)
    else
        return _backward_trellis_nonbinary(M, n, B_R, q, verbose)
    end
end

function _backward_trellis_binary(M::Matrix{T}, n::Int, B_R::Int, verbose::Bool=true) where T
    k, _ = size(M)
    
    L, R = _get_LR_indices(M)
    prev_layer = Dict{Tuple{UInt128, Bool}, Int}((UInt128(0), false) => 0)
    
    active_sets   = [Int[] for _ in 1:n]
    start_bw_sets = [Int[] for _ in 1:n]
    work_sets     = [Int[] for _ in 1:n]
    M_cols        = zeros(UInt128, n)
    work_masks    = zeros(UInt128, n)

    for col in n:-1:(B_R + 1)
        mask = UInt128(0)
        col_val = UInt128(0)
        for row in 1:k
            if L[row] <= col < R[row] push!(active_sets[col], row) end
            if R[row] == col push!(start_bw_sets[col], row) end
            if L[row] <= col <= R[row] 
                push!(work_sets[col], row)
                mask |= (UInt128(1) << (row - 1))
            end
            if !iszero(M[row, col]) col_val |= (UInt128(1) << (row - 1)) end
        end
        work_masks[col] = mask
        M_cols[col] = col_val
    end
    
    active_at_BR = [row for row in 1:k if L[row] <= B_R < R[row]]

    p = verbose ? Progress(n - B_R, 0.1, "Building Backward Trellis: ") : nothing
    
    for col in n:-1:(B_R + 1)
        next_layer = Dict{Tuple{UInt128, Bool}, Int}()
        
        active_curr = col == B_R + 1 ? active_at_BR : active_sets[col-1]
        sizehint!(next_layer, 1 << length(active_curr))
        
        starting_bw = start_bw_sets[col]
        num_branches = 1 << length(starting_bw)
        deltas = zeros(UInt128, num_branches)
        for b in 0:(num_branches-1)
            d = UInt128(0)
            for (idx, row) in enumerate(starting_bw)
                if ((b >> (idx - 1)) & 1) == 1 d |= (UInt128(1) << (row - 1)) end
            end
            deltas[b+1] = d
        end
        
        keep_mask = UInt128(0)
        for row in active_curr keep_mask |= (UInt128(1) << (row - 1)) end
        
        col_val = M_cols[col]
        w_mask = work_masks[col]
        
        for ((u_prev, is_nonzero), prev_wt) in prev_layer
            @inbounds for b in 1:num_branches
                delta = deltas[b]
                u_work = u_prev | delta
                
                col_wt = count_ones(u_work & col_val & w_mask) & 1
                new_wt = prev_wt + col_wt
                
                u_next = u_work & keep_mask
                new_is_nonzero = is_nonzero | (b > 1) | (u_next > 0)
                
                key = (u_next, new_is_nonzero)
                if !haskey(next_layer, key) || new_wt < next_layer[key]
                    next_layer[key] = new_wt
                end
            end
        end
        prev_layer = next_layer
        verbose && next!(p)
    end
    
    # DO NOT UNPACK! Return raw hardware registers for the BZ Middle Search
    return prev_layer
end

function _backward_trellis_ternary(M::Matrix{T}, n::Int, B_R::Int, verbose::Bool=true) where T
    k, _ = size(M)
    F = parent(M[1, 1])
    elements = collect(F)
    T_zero = zero(F)
    
    elem_to_u8 = Dict(elements[1] => UInt8(0), elements[2] => UInt8(1), elements[3] => UInt8(2))
    u8_to_elem = Dict(UInt8(0) => elements[1], UInt8(1) => elements[2], UInt8(2) => elements[3])
    
    M_u8 = zeros(UInt8, k, n)
    for i in 1:k, j in 1:n M_u8[i,j] = elem_to_u8[M[i,j]] end
    
    L, R = _get_LR_indices(M)
    
    prev_layer = Dict{Tuple{UInt128, Bool}, Int}((UInt128(0), false) => 0)
    
    active_sets   = [Int[] for _ in 1:n]
    start_bw_sets = [Int[] for _ in 1:n]
    work_sets     = [Int[] for _ in 1:n]

    for col in n:-1:(B_R + 1)
        for row in 1:k
            if L[row] <= col < R[row] push!(active_sets[col], row) end
            if R[row] == col push!(start_bw_sets[col], row) end
            if L[row] <= col <= R[row] push!(work_sets[col], row) end
        end
    end
    
    active_at_BR = [row for row in 1:k if L[row] <= B_R < R[row]]

    p = verbose ? Progress(n - B_R, 0.1, "Building Backward Trellis (GF(3)): ") : nothing
    u_buf = zeros(UInt8, k)
    
    for col in n:-1:(B_R + 1)
        next_layer = Dict{Tuple{UInt128, Bool}, Int}()
        
        active_curr = col == B_R + 1 ? active_at_BR : active_sets[col-1]
        sizehint!(next_layer, 3^length(active_curr))
        
        active_prev = active_sets[col]
        starting_bw = start_bw_sets[col]
        working     = work_sets[col]
        
        for ((u_packed, is_nonzero), prev_wt) in prev_layer
            
            fill!(u_buf, UInt8(0))
            for (idx, row) in enumerate(active_prev)
                u_buf[row] = UInt8((u_packed >> (2 * (idx - 1))) & 3)
            end
            
            for branch_scalars in Iterators.product(fill((UInt8(0), UInt8(1), UInt8(2)), length(starting_bw))...)
                @inbounds for (idx, row) in enumerate(starting_bw) u_buf[row] = branch_scalars[idx] end
                
                c_i = UInt32(0)
                @inbounds for row in working
                    c_i += u_buf[row] * M_u8[row, col]
                end
                
                col_wt = (c_i % 3) == 0 ? 0 : 1
                new_wt = prev_wt + col_wt
                
                next_packed = UInt128(0)
                new_is_nonzero = is_nonzero
                
                for val in branch_scalars
                    if val != 0 new_is_nonzero = true end
                end
                
                @inbounds for (idx, row) in enumerate(active_curr)
                    val = u_buf[row]
                    if val != 0 new_is_nonzero = true end
                    next_packed |= (UInt128(val) << (2 * (idx - 1)))
                end
                
                key = (next_packed, new_is_nonzero)
                if !haskey(next_layer, key) || new_wt < next_layer[key]
                    next_layer[key] = new_wt
                end
            end
        end
        prev_layer = next_layer
        verbose && next!(p)
    end
    
    final_dict = Dict{Tuple{Vector{T}, Bool}, Int}()
    for ((u_packed, is_nz), wt) in prev_layer
        u_unpacked = fill(T_zero, length(active_at_BR))
        for (idx, row) in enumerate(active_at_BR)
            val_u8 = UInt8((u_packed >> (2 * (idx - 1))) & 3)
            u_unpacked[idx] = u8_to_elem[val_u8]
        end
        final_dict[(u_unpacked, is_nz)] = wt
    end
    
    return final_dict
end

function _backward_trellis_quaternary(M::Matrix{T}, n::Int, B_R::Int, verbose::Bool=true) where T
    k, _ = size(M)
    F = parent(M[1, 1])
    elements = collect(F)
    T_zero = zero(F)
    
    elem_to_u8 = Dict(elements[1] => UInt8(0), elements[2] => UInt8(1), elements[3] => UInt8(2), elements[4] => UInt8(3))
    u8_to_elem = Dict(UInt8(0) => elements[1], UInt8(1) => elements[2], UInt8(2) => elements[3], UInt8(3) => elements[4])
    
    GF4_MULT = UInt8[
        0 0 0 0;
        0 1 2 3;
        0 2 3 1;
        0 3 1 2
    ]
    
    M_u8 = zeros(UInt8, k, n)
    for i in 1:k, j in 1:n M_u8[i,j] = elem_to_u8[M[i,j]] end
    
    L, R = _get_LR_indices(M)
    
    prev_layer = Dict{Tuple{UInt128, Bool}, Int}((UInt128(0), false) => 0)
    
    active_sets   = [Int[] for _ in 1:n]
    start_bw_sets = [Int[] for _ in 1:n]
    work_sets     = [Int[] for _ in 1:n]

    for col in n:-1:(B_R + 1)
        for row in 1:k
            if L[row] <= col < R[row] push!(active_sets[col], row) end
            if R[row] == col push!(start_bw_sets[col], row) end
            if L[row] <= col <= R[row] push!(work_sets[col], row) end
        end
    end
    
    active_at_BR = [row for row in 1:k if L[row] <= B_R < R[row]]

    p = verbose ? Progress(n - B_R, 0.1, "Building Backward Trellis (GF(4)): ") : nothing
    u_buf = zeros(UInt8, k)
    
    for col in n:-1:(B_R + 1)
        next_layer = Dict{Tuple{UInt128, Bool}, Int}()
        
        active_curr = col == B_R + 1 ? active_at_BR : active_sets[col-1]
        sizehint!(next_layer, 4^length(active_curr))
        
        active_prev = active_sets[col]
        starting_bw = start_bw_sets[col]
        working     = work_sets[col]
        
        for ((u_packed, is_nonzero), prev_wt) in prev_layer
            
            fill!(u_buf, UInt8(0))
            for (idx, row) in enumerate(active_prev)
                u_buf[row] = UInt8((u_packed >> (2 * (idx - 1))) & 3)
            end
            
            for branch_scalars in Iterators.product(fill((UInt8(0), UInt8(1), UInt8(2), UInt8(3)), length(starting_bw))...)
                @inbounds for (idx, row) in enumerate(starting_bw) u_buf[row] = branch_scalars[idx] end
                
                c_i = UInt8(0)
                @inbounds for row in working
                    c_i ⊻= GF4_MULT[u_buf[row] + 1, M_u8[row, col] + 1]
                end
                
                col_wt = c_i == 0 ? 0 : 1
                new_wt = prev_wt + col_wt
                
                next_packed = UInt128(0)
                new_is_nonzero = is_nonzero
                
                for val in branch_scalars
                    if val != 0 new_is_nonzero = true end
                end
                
                @inbounds for (idx, row) in enumerate(active_curr)
                    val = u_buf[row]
                    if val != 0 new_is_nonzero = true end
                    next_packed |= (UInt128(val) << (2 * (idx - 1)))
                end
                
                key = (next_packed, new_is_nonzero)
                if !haskey(next_layer, key) || new_wt < next_layer[key]
                    next_layer[key] = new_wt
                end
            end
        end
        prev_layer = next_layer
        verbose && next!(p)
    end
    
    final_dict = Dict{Tuple{Vector{T}, Bool}, Int}()
    for ((u_packed, is_nz), wt) in prev_layer
        u_unpacked = fill(T_zero, length(active_at_BR))
        for (idx, row) in enumerate(active_at_BR)
            val_u8 = UInt8((u_packed >> (2 * (idx - 1))) & 3)
            u_unpacked[idx] = u8_to_elem[val_u8]
        end
        final_dict[(u_unpacked, is_nz)] = wt
    end
    
    return final_dict
end

function _backward_trellis_nonbinary(M::Matrix{T}, n::Int, B_R::Int, q::Int, verbose::Bool=true) where T
    k, _ = size(M)
    F = parent(M[1, 1])
    elements = collect(F)
    T_zero = zero(F)
    
    L, R = _get_LR_indices(M)
    
    prev_layer = Dict{Tuple{Vector{T}, Bool}, Int}((Vector{T}(), false) => 0)
    
    active_sets   = [Int[] for _ in 1:n]
    start_bw_sets = [Int[] for _ in 1:n]
    work_sets     = [Int[] for _ in 1:n]

    for col in n:-1:(B_R + 1)
        for row in 1:k
            if L[row] <= col < R[row] push!(active_sets[col], row) end
            if R[row] == col push!(start_bw_sets[col], row) end
            if L[row] <= col <= R[row] push!(work_sets[col], row) end
        end
    end
    
    active_at_BR = [row for row in 1:k if L[row] <= B_R < R[row]]

    p = verbose ? Progress(n - B_R, 0.1, "Building Backward Trellis (Non-Binary): ") : nothing
    u_buf = fill(T_zero, k)
    
    for col in n:-1:(B_R + 1)
        next_layer = Dict{Tuple{Vector{T}, Bool}, Int}()
        
        active_curr = col == B_R + 1 ? active_at_BR : active_sets[col-1]
        sizehint!(next_layer, q^length(active_curr))
        
        active_prev = active_sets[col]
        starting_bw = start_bw_sets[col]
        working     = work_sets[col]
        
        for ((prev_scalars, is_nonzero), prev_wt) in prev_layer
            for branch_scalars in Iterators.product(fill(elements, length(starting_bw))...)
                
                fill!(u_buf, T_zero)
                @inbounds for (idx, row) in enumerate(active_prev) u_buf[row] = prev_scalars[idx] end
                @inbounds for (idx, row) in enumerate(starting_bw) u_buf[row] = branch_scalars[idx] end
                
                c_i = T_zero
                @inbounds for row in working
                    c_i += u_buf[row] * M[row, col]
                end
                
                col_wt = iszero(c_i) ? 0 : 1
                new_wt = prev_wt + col_wt
                
                next_scalars = [u_buf[row] for row in active_curr]
                
                new_is_nonzero = is_nonzero || any(!iszero, branch_scalars) || any(!iszero, next_scalars)
                
                key = (next_scalars, new_is_nonzero)
                if !haskey(next_layer, key) || new_wt < next_layer[key]
                    next_layer[key] = new_wt
                end
            end
        end
        prev_layer = next_layer
        verbose && next!(p)
    end
    
    return prev_layer
end

"""
    _minimum_distance_hybrid(C::AbstractLinearCode; max_span::Int=15, num_trials::Int=50, block_size::Int=0, verbose::Bool=true)

Computes the minimum distance using the Meet-in-the-Middle BZ-Trellis bridge.
Intelligently routes between aggressive column permutations and natural Quasi-Cyclic 
shift-symmetry pruning based on whether the "Death Zone" can be eliminated.
"""
function _minimum_distance_hybrid(C::AbstractLinearCode; max_span::Int=15, num_trials::Int=50, block_size::Int=0, verbose::Bool=true)
    k, n = C.k, C.n
    q = Int(order(C.F))
    
    # The BZ-DFS bridge evaluates `u * M` directly.
    # We must strictly use the generator matrix to compute the primal distance.
    mat = Array(generator_matrix(C))
    best_M = mat

    if block_size > 0
        verbose && println("Quasi-Cyclic Block Size: $block_size provided. Evaluating routing strategy...")
        
        perm_M, _, perm_peak_E = optimize_trellis_permutation(mat, num_trials)
        
        if perm_peak_E <= max_span
            verbose && println("Strategy A: Permutations crushed the profile ($perm_peak_E <= $max_span). Bypassing symmetries for Pure Viterbi.")
            best_M = perm_M
            block_size = 0 
        else
            verbose && println("Strategy B: Permutations failed to clear the Death Zone ($perm_peak_E > $max_span).")
            verbose && println("Reverting to natural column order to exploit First Block symmetry pruning.")
            best_M = copy(mat)
            _make_trellis_oriented!(best_M)
        end
    else
        verbose && println("Optimizing Trellis profile via permutations...")
        best_M, _, _ = optimize_trellis_permutation(mat, num_trials)
    end
    
    L, R = _get_LR_indices(best_M)
    past, future = past_future_profiles(L, R, n)
    E_profile = [k - past[b] - future[b+1] for b in 1:n]
    peak_E = maximum(E_profile)
    
    if peak_E <= max_span
        verbose && println("Profile is thin. Routing to Pure Sectionalized Viterbi...")
        boundaries = optimal_sectionalization(best_M, q)
        
        # USE MIN WEIGHT ENGINE INSTEAD OF CWE
        return _min_weight_TP_Viterbi_sectionalized(best_M, boundaries, verbose)
    end
    
    B_L = findfirst(e -> e > max_span, E_profile) - 1
    B_R = findlast(e -> e > max_span, E_profile) 
    
    if block_size > 0 && B_L < block_size
        verbose && println("Extending Left Trellis to column $block_size to guarantee First Block Pruning.")
        B_L = block_size
    end
    
    verbose && println("Pinching TOF Matrix. Left: [1, $B_L], BZ: [$(B_L+1), $B_R], Right: [$(B_R+1), $n]")
    
    # Spawn the independent boundary trellises on separate threads
    t_left = Threads.@spawn _forward_trellis(best_M, B_L, q, block_size)
    t_right = Threads.@spawn _backward_trellis(best_M, n, B_R, q)
    
    # The main thread halts here until BOTH dictionaries are fully constructed
    left_dict = fetch(t_left)
    right_dict = fetch(t_right)
    
    verbose && println("Executing Bounded BZ-DFS Bridge...")
    d = _BZ_middle_search(best_M, L, R, B_L, B_R, left_dict, right_dict, q)
    
    return d
end
