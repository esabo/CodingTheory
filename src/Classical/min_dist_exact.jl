# Copyright (c) 2021 - 2026 Eric Sabo, Benjamin Ide, Michael Vasmer, David Marquis
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
  # Binary Helper Functions
#############################

function _convert_binary_to_int_matrix(A::CTMatrixTypes)
    nr, nc = size(A)
    B = zeros(Int, nr, nc)
    for r in 1:nr
        for c in 1:nc
            B[r, c] = is_zero(A[r, c]) ? 0 : 1
        end
    end
    return B
end

# auto generated docstring
"""
    _precompute_pruning_bounds_binary(A_raw::Matrix{Int}, n::Int, k::Int)

Precomputes the dual theoretical bounds required for the Unified Triangle Inequality 
pruning step in the Brouwer-Zimmermann recursive search over GF(2).

For any given number of remaining message bits `r`, this function calculates both the 
absolute minimum and absolute maximum parity weight those `r` bits could possibly generate. 

# Mathematical Bounds Computed:
1. **Griesmer Bound (`lbt`)**: 
   The theoretical lower bound. Calculated via a binary search over the Griesmer 
   theorem for binary codes. It dictates the minimum possible weight `W_rem` that 
   `r` rows of the parity matrix `A` *must* produce.
   
2. **Maximum Cancellation Bound (`max_canc`)**: 
   The theoretical upper bound. Calculated by multiplying `r` by the weight of the 
   heaviest single column in `A_raw`, capped at the total parity length `n - k`. 
   It dictates the absolute worst-case drop in weight via bitwise XOR cancellation.

These two arrays guarantee that the unknown future parity weight `W_rem` strictly 
satisfies: `lbt[r] <= W_rem <= max_canc[r]`.

# Arguments
- `A_raw::Matrix{Int}`: The `(n-k) × k` parity portion of the systematic generator matrix.
- `n::Int`: The total length of the code.
- `k::Int`: The dimension of the code.

# Returns
- `lbt::Vector{Int}`: Array of length `k+1`. `lbt[r+1]` is the Griesmer minimum weight for `r` bits.
- `max_canc::Vector{Int}`: Array of length `k+1`. `max_canc[r+1]` is the maximum possible weight for `r` bits.
"""
function _precompute_pruning_bounds_binary(A_raw::Matrix{Int}, k::Int)
    k_rows, n_tail = size(A_raw) 
    lbt = zeros(Int, k + 1)
    max_canc = zeros(Int, k + 1)
    
    # Weight of parity columns (how many 1s in each column)
    col_weights = [sum(view(A_raw, :, j)) for j in 1:size(A_raw, 2)]
    
    # We use the lightest column as the absolute minimum contribution 
    # of any non-zero combination in the parity part.
    min_col_wt = isempty(col_weights) ? 0 : minimum(col_weights)

    lbt[1] = 0
    for r in 1:k
        # Conservative: a message of weight r must produce 
        # AT LEAST 1 parity bit (unless the matrix is singular).
        # We use a floor of 1, or better, use the column density.
        lbt[r + 1] = (r > 0) ? 1 : 0 
    end
    
    # Max cancellation (upper bound) logic remains valid
    row_weights = [sum(view(A_raw, i, :)) for i in 1:k_rows]
    sort!(row_weights, rev=true)
    current_max_sum = 0
    for r in 1:k
        current_max_sum += row_weights[r]
        max_canc[r + 1] = min(n_tail, current_max_sum)
    end
    
    return lbt, max_canc
end

"""
    _pack_matrix_simd(A_raw::Matrix{Int})

Compresses a binary matrix into an array of UInt64 chunks for bitwise operations.
"""
function _pack_matrix_simd(A_raw::Matrix{Int})
    n_tail, k = size(A_raw)
    num_chunks = cld(n_tail, 64)
    A_packed = zeros(UInt64, num_chunks, k)
    
    for j in 1:k
        for i in 1:n_tail
            if A_raw[i, j] == 1
                chunk_idx = (i - 1) ÷ 64 + 1
                bit_idx = (i - 1) % 64
                A_packed[chunk_idx, j] |= (UInt64(1) << bit_idx)
            end
        end
    end
    return A_packed
end

"""
    _fast_simd_wt(v::Vector{UInt64})

Calculates the Hamming weight of a bit-packed parity vector using CPU POPCNT.
"""
@inline function _fast_simd_wt(v::Vector{UInt64})
    s = 0
    @inbounds @simd for i in eachindex(v) 
        s += count_ones(v[i]) 
    end
    return s
end

"""
    _precompute_weight2_table(A_packed::Matrix{UInt64})

The "Brouwer Trick": Precomputes the exact minimum parity weight for all pairs 
of remaining message bits. Replaces deep tree branching with a O(1) lookup.
"""
function _precompute_weight2_table(A_packed::Vector{Vector{UInt64}})
    k = length(A_packed)
    # table[i, j] stores count_ones(A_packed[i] ^ A_packed[j])
    table = zeros(Int, k, k)
    
    for i in 1:k
        for j in i+1:k
            w = 0
            @inbounds @simd for c in eachindex(A_packed[i])
                w += count_ones(A_packed[i][c] ⊻ A_packed[j][c])
            end
            table[i, j] = w
            table[j, i] = w
        end
    end
    return table
end

"""
    _precompute_refined_lbt(A_raw::Matrix{Int}, n::Int, k::Int)

Creates the Lower Bound Table (LBT). lbt[r] defines the absolute minimum 
parity weight that a message of weight 'r' MUST generate.
"""
function _precompute_refined_lbt(A_packed, n, k)
    msg_rows = size(A_packed, 2)
    lbt = zeros(Int, msg_rows + 1)
    
    # Calculate weights of the columns (parity bits) [cite: 24]
    # BUG FIX: We must use a conservative minimum, not a sorted sum 
    # that assumes we always pick the heaviest columns.
    col_weights = [count_ones(A_packed[j]) for j in 1:msg_rows]
    min_weight = minimum(col_weights) # [cite: 12]
    
    for r in 1:msg_rows
        # Griesmer floor: max d such that Σ ⌈d/2ⁱ⌉ <= (n-k) [cite: 25]
        low, high, griesmer = 0, (n - k), 0
        while low <= high
            mid = (low + high) ÷ 2
            g_sum, d_val = 0, Float64(mid)
            for i in 1:r
                term = ceil(Int, d_val)
                g_sum += term
                d_val /= 2
                term <= 1 && (g_sum += (r - i); break) # [cite: 26]
            end
            g_sum <= (n - k) ? (griesmer = mid; low = mid + 1) : (high = mid - 1) # [cite: 26]
        end
        
        # BUG FIX: The lower bound for 'r' bits cannot exceed r * min_weight.
        # If your A matrix is sparse, Griesmer is often too optimistic.
        lbt[r] = griesmer # [cite: 27]
    end
    return lbt
end

"""
    _map_automorphisms(auts::Vector{Vector{Int}}, col_perm::Vector{Int})

Filters and projects the full length-n permutations down to the k-length 
information space, keeping only those that stabilize the information set.
"""
function _map_automorphisms(auts::Vector{Vector{Int}}, col_perm::Vector{Int})
    internal_auts = Vector{Vector{Int}}()
    k = length(col_perm)
    
    # Fast lookup for whether a target index is inside our information set
    inv_perm = Dict(col_perm[i] => i for i in 1:k)
    
    for σ in auts
        stabilizes = true
        mapped_σ = zeros(Int, k)
        for i in 1:k
            orig_idx = col_perm[i]
            mapped_orig = σ[orig_idx]
            
            if haskey(inv_perm, mapped_orig)
                mapped_σ[i] = inv_perm[mapped_orig]
            else
                stabilizes = false # Permutation mapped a bit into the parity set
                break
            end
        end
        
        if stabilizes
            push!(internal_auts, mapped_σ)
        end
    end
    return internal_auts
end

"""
    _is_canonical(msg::Vector{Int}, auts::Vector{Vector{Int}})

Lexicographical Orbit Leader check. Ensures we only evaluate the "smallest" 
version of any symmetric codeword pattern.
"""
function _is_canonical(msg::Vector{Int}, auts::Vector{Vector{Int}})
    for σ in auts
        for i in eachindex(msg)
            mapped_idx = σ[i]
            # msg[i] is the current value, msg[mapped_idx] is the permuted value
            if msg[i] > msg[mapped_idx]
                return false # Not canonical, kill branch
            elseif msg[i] < msg[mapped_idx]
                break # Strictly less than this permutation, safe to proceed
            end
        end
    end
    return true
end

"""
    _generate_known_automorphisms(C::AbstractLinearCod)

Generates a set of permutation vectors (automorphisms) for known code families.
These vectors can be passed to the `minimum_distance_master` search engine.
"""
function _generate_known_automorphisms(C::AbstractLinearCode)
    if typeof(C) <: AbstractCyclicCode
        return _generate_cyclic_auts(C.n)
        # TODO make this type
    # elseif isa(C, CyclicCode2D)
    #     return _generate_2d_cyclic_auts(C.n1, C.n2)
        # TODO make this type
    elseif isa(C, ExtendedQRCode)
        # length is p + 1, so p = n-1
        return _generate_extended_qr_auts(C.n - 1)
        # TODO make this type
    elseif isa(C, ProductCode)
        # For a product code C1 x C2, where k = k1 * k2
        # Shift in C1 across all blocks of C2
        row_shift = zeros(Int, C.k)
        for i in 0:C.k1-1, j in 0:C.k2-1
            row_shift[i * C.k2 + j + 1] = ((i + 1) % C.k1) * C.k2 + j + 1
        end
        push!(auts, row_shift)

        # Shift in C2 across all blocks of C1
        col_shift = zeros(Int, C.k)
        for i in 0:C.k1-1, j in 0:C.k2-1
            col_shift[i * C.k2 + j + 1] = i * C.k2 + ((j + 1) % C.k2) + 1
        end
        push!(auts, col_shift)
    elseif isa(C, ReedMullerCode)
        return _generate_RM_auts(C.m)
    else
        return Vector{Vector{Int}}() # No known automorphisms for this code family
    end
end

"""
    _generate_cyclic_auts(n::Int)

Returns the permutation generator for the cyclic group of order n.
"""
_generate_cyclic_auts(n::Int) = [[(i % n) + 1 for i in 1:n]]

"""
    _generate_toric_auts(n1::Int, n2::Int)

Returns the permutation generators for a 2D cyclic code of size n1 × n2.
Assumes the n1 * n2 coordinates are flattened row-by-row.
"""
function _generate_2d_cyclic_auts(n1::Int, n2::Int)
    n = n1 * n2
    shift_x = zeros(Int, n)
    shift_y = zeros(Int, n)
    
    for y in 1:n2
        for x in 1:n1
            idx = (y - 1) * n1 + x
            
            # Shift along the x-axis (wrap around n1)
            idx_x = (y - 1) * n1 + (x % n1) + 1
            
            # Shift along the y-axis (wrap around n2)
            idx_y = (y % n2) * n1 + x
            
            shift_x[idx] = idx_x
            shift_y[idx] = idx_y
        end
    end
    
    return [shift_x, shift_y]
end

"""
    _generate_extended_qr_auts(p::Int)

Returns the two generators for PSL(2, p) which form the automorphism group 
of an Extended QR code of length p + 1.
Indices 1 to p represent field elements 0 to p-1. Index p+1 represents infinity.
"""
function _generate_extended_qr_auts(p::Int)
    # Generator 1: Translation x ↦ x + 1 (mod p), with ∞ ↦ ∞
    g1 = zeros(Int, p + 1)
    for x in 0:p-1
        g1[x + 1] = ((x + 1) % p) + 1
    end
    g1[p + 1] = p + 1 

    # Generator 2: Inversion x ↦ -1/x (mod p), with 0 ↔ ∞
    g2 = zeros(Int, p + 1)
    g2[1] = p + 1     # 0 maps to ∞
    g2[p + 1] = 1     # ∞ maps to 0
    for x in 1:p-1
        # Find modular inverse, then negate
        inv_x = invmod(x, p)
        minus_inv_x = (-inv_x + p) % p
        g2[x + 1] = minus_inv_x + 1
    end
    
    return [g1, g2]
end

"""
    _generate_RM_auts(m::Int)

Returns the generators for the General Affine Group GA(m, 2) acting on a 
Reed-Muller code of length n = 2^m.
"""
function _generate_RM_auts(m::Int)
    n = 2^m
    auts = Vector{Vector{Int}}()
    
    # 1. Affine Translations (m generators: flip the i-th bit of the coordinate index)
    for bit in 0:m-1
        t = zeros(Int, n)
        for x in 0:n-1
            t[x + 1] = (x ⊻ (1 << bit)) + 1
        end
        push!(auts, t)
    end
    
    # 2. Linear GL(m, 2) Generators: Adjacent bit swaps (Generates the Symmetric group on bits)
    for bit in 0:m-2
        s = zeros(Int, n)
        for x in 0:n-1
            b1 = (x >> bit) & 1
            b2 = (x >> (bit + 1)) & 1
            
            # Clear the two bits, then swap their positions
            new_x = x & ~(1 << bit) & ~(1 << (bit + 1))
            new_x |= (b2 << bit) | (b1 << (bit + 1))
            s[x + 1] = new_x + 1
        end
        push!(auts, s)
    end
    
    # 3. Linear GL(m, 2) Generators: Single Transvection (Add bit 0 to bit 1)
    if m > 1
        v = zeros(Int, n)
        for x in 0:n-1
            b0 = x & 1
            b1 = (x >> 1) & 1
            
            # XOR bit 0 into bit 1
            new_b1 = b1 ⊻ b0
            new_x = (x & ~(1 << 1)) | (new_b1 << 1)
            v[x + 1] = new_x + 1
        end
        push!(auts, v)
    end
    
    return auts
end

"""
    _score_info_set(A_raw::Matrix{Int})

Scores a parity matrix based on density and variance. Denser parity matrices 
force the BZ lookahead to prune branches higher in the tree.
"""
function _score_info_set(A_raw::Matrix{Int})
    k = size(A_raw, 2)
    col_weights = [sum(A_raw[:, j]) for j in 1:k]
    avg_weight = sum(col_weights) / k
    variance = sum((col_weights .- avg_weight).^2) / k
    
    # Reward high average weight, penalize high variance (we want uniformly dense columns)
    return avg_weight - (0.1 * variance) 
end

# TODO: this does not produce the optimal set of matrices
# see section 7.3 of White's thesis for comments on this
function information_sets(G::CTMatrixTypes, alg::Symbol = :Edmonds; permute::Bool = false, only_A::Bool = false)

    alg ∈ (:Brouwer, :Zimmermann, :White, :Chen, :Bouyuklieva, :Edmonds) || throw(ArgumentError("Unknown information set algorithm. Expected `:Brouwer`, `:Zimmermann`, `:White`, `:Chen`, `:Bouyuklieva`, or `:Edmonds`."))
    # TODO should rref to begin with and remove empty rows?
    nr, nc = size(G)
    gen_mats = Vector{}()
    perms = Vector{}()
    rnks = Vector{Int}()
    
    rnk = nr
    #TODO Brouw and Zimm should use the same code for rref then discard the last matrix
    if alg == :Brouwer
        for i in 0:Int(floor(nc / nr)) - 1
            start_ind = i * nr + 1
            rnk, Gi, Pi = _rref_col_swap(G, 1:nr, start_ind:nc)
            if rnk < nr # for Brouwer the Gi must all have full rank
                break
            end

            if only_A
                Ai = Gi[:, setdiff(1:nc, start_ind:(i + 1) * nr)]
                push!(gen_mats, Ai)
                push!(perms, Pi)
                push!(rnks, rnk)
            else
                if permute
                    # permute identities to the front
                    pivots = collect(start_ind:(i + 1) * nr)
                    σ = [pivots; setdiff(1:nc, pivots)]
                    Gi = Gi[:, σ]
                    Pi = Pi[:, σ]
                end
                push!(gen_mats, Gi)
                push!(perms, Pi)
                push!(rnks, rnk)
            end
        end
    elseif alg == :Zimmermann
        for i in 0:Int(floor(nc / nr))
            rnk, Gi, Pi = _rref_col_swap(G, 1:nr, i * rnk + 1:nc)
            # BUG this option seems really bad for type stability
            if ismissing(Pi)
                Pi = identity_matrix(base_ring(G), nc)
            end
            if only_A
                Ai = Gi[:, setdiff(1:nc, i * nr + 1:i * nr + rnk)]
                push!(gen_mats, Ai)
                push!(perms, Pi)
                push!(rnks, rnk)
            else
                if permute
                    # Ensure the end of the range doesn't exceed the number of columns
                    start_idx = i * nr + 1
                    end_idx = min((i + 1) * nr, nc)
    
                    pivots = collect(start_idx:end_idx)
    
                    # If the last block is smaller than nr, the logic still holds 
                    # because setdiff will handle the remaining indices correctly.
                    σ = [pivots; setdiff(1:nc, pivots)]
    
                    Gi = Gi[:, σ]
                    Pi = Pi[:, σ]
                end
                push!(gen_mats, Gi)
                push!(perms, Pi)
                push!(rnks, rnk)
            end
        end
    elseif alg == :White
        # TODO: this is not true when the parity-check matrix is true
        # the expansion factor of the code
        for i in 0:div(nc, nr) - 1
            # could use Gi here instead of G
            rnk, Gi, Pi = _rref_col_swap(G, 1:nr, i * nr + 1:(i + 1) * nr)
            # display(Gi)
            # println(rnk)
            push!(gen_mats, Gi)
            push!(perms, Pi)
            push!(rnks, rnk)
        end
    elseif alg == :Chen
        Gi, _, Pi, rnk = _standard_form(G)
        if only_A
            Ai = Gi[:, rnk + 1:nc]
            push!(gen_mats, Ai)
            push!(perms, Pi)
            push!(rnks, rnk)
        else
            push!(gen_mats, Gi)
            push!(perms, Pi)
            push!(rnks, rnk)
        end
    elseif alg == :Bouyuklieva
        remaining_cols = collect(1:nc)
        while !isempty(remaining_cols)
            # 1. Greedy Rank Check on remaining columns
            G_rem = G[:, remaining_cols]
            # Use your existing RREF helper
            rnk, G_rref, pivots = _rref_col_swap(G_rem, 1:nr, 1:size(G_rem, 2))
            
            # 2. Extract Absolute Indices
            set_indices = remaining_cols[pivots[1:rnk]]
            push!(rnks, rnk)
            
            # 3. Create the systematic generator matrix (G_T)
            # The paper requires columns in T to form an identity submatrix 
            other_cols = setdiff(1:nc, set_indices)
            σ = [set_indices; other_cols]
            
            # Construct the matrix Gp and the Permutation matrix Pp
            # Ensure the rnk x rnk identity is at the top-left
            Gp = G[:, σ]
            _make_systematic_gf!(Gp, collect(1:nc), rnk)
            
            # Create a proper permutation matrix (or vector) to match your BZ loop
            # If your BZ loop expects a matrix, use: Pp = _permutation_matrix(σ, nc)
            # If it expects a vector, just use σ
            push!(gen_mats, Gp)
            push!(perms, σ) 

            # 4. Enforce Disjointness [cite: 827]
            filter!(x -> x ∉ set_indices, remaining_cols)
            
            if rnk == 0 break end
        end
    end

    return gen_mats, perms, rnks
end

"""
    _generate_scored_zimmermann_mats(C::AbstractLinearCode, max_m::Int=3)

Generates 'm' overlapping systematic matrices to facilitate the 
Zimmermann lower-bound warp (l_bound = m * (w + 1)).
"""
function _generate_scored_zimmermann_mats(C, max_m::Int=3)
    k, n = size(C.G_stand)
    
    # 1. Generate multiple valid information sets and score them
    # (Assuming information_sets() is an existing library function)
    A_mats, perms, _ = information_sets(C.G_stand, :Zimmermann, permute=true, only_A=false)
    
    best_score = -1.0
    best_idx = 1
    for i in 1:min(length(A_mats), 10) # Check top 10 candidates
        score = _score_info_set(Int.(A_mats[i][:, k+1:n]'))
        if score > best_score
            best_score, best_idx = score, i
        end
    end
    
    matrices = [(G = A_mats[best_idx], perm = perms[best_idx])]
    curr_G = matrices[1].G
    
    # 2. Iteratively reduce the parity part to create overlapping sets
    for i in 2:max_m
        A = curr_G[:, k+1:n]
        if rank(A) == 0 break end
        
        try
            G_new, P_new, _ = information_sets(A, :Zimmermann, permute=true)
            push!(matrices, (G = G_new[1], perm = P_new[1]))
            curr_G = G_new[1]
        catch
            break
        end
    end
    
    return matrices
end

function _information_set_lower_bound(r::Int, n::Int, k::Int, l::Int, rank_defs::Vector{Int},
    info_set_alg::Symbol; even::Bool = false, doubly_even::Bool = false, triply_even::Bool = false)

    info_set_alg ∈ (:auto, :Brouwer, :Zimmermann, :White, :Chen, :Bouyuklieva, :Edmonds) || throw(ArgumentError("Unknown information set algorithm. Expected `:auto`, `:Brouwer`, `:Zimmermann`, `:White`, `:Chen`, `:Bouyuklieva`, or `:Edmonds`."))

    lower = 0
    if info_set_alg == :Brouwer
        lower = r * length(rank_defs)
    elseif info_set_alg == :Zimmermann
        h = length(rank_defs)
        lower = count(x -> x != 0, rank_defs) 
        for i in 1:h
            lower += maximum([0, r - rank_defs[i]]) 
        end
    elseif info_set_alg == :Chen
        lower = Int(ceil(n * r / k))
    elseif info_set_alg == :White
        lower = 0
        for i in 1:l
	        lower += Int(ceil(n * maximum([0, r - rank_defs[i]]) / (l * (k + rank_defs[i]))))
        end
    elseif info_set_alg == :Bouyuklieva
        # BB21 Logic based on Theorem 1: w <= sum(a_i) + t + r - 1 [cite: 863]
        # Assumptions for this dispatch:
        # rank_defs: contains the a_i values (rows combined per matrix) [cite: 829]
        # r: corresponds to the number of reduced sets currently active [cite: 863, 881]
        
        # t is the count of full systematic sets (where columns form identity Ik) [cite: 825, 828]
        # In our disjoint partitioning, this is the number of sets with size k [cite: 858, 879]
        # We determine t from the length of rank_defs minus our active reduced sets r
        # L = sum(a_i) - s + (k - rank(union_of_subcodes))
        # Since T1 is an info set, rank(union) = k.
        # L = sum(a_i) - length(a_values) + 1
        return sum(rank_defs) - length(rank_defs) + 1
    # elseif info_set_alg == :Edmonds
    #     continue
    end
    if lower > 0 && verbose
        println("Initial lower bound raised to: $lower")
    end

    (!triply_even && !doubly_even && even) && (lower += lower % 2;)
    (!triply_even && doubly_even) && (lower += 4 - lower % 4;)
    triply_even && (lower += 8 - lower % 8;)
    return lower
end

"""
    _partition_disjoint_systematic_sets(G::CTMatrixTypes)

Partitions coordinates {1...n} into disjoint systematic sets T_1...T_s.
Returns (gen_mats, perms, rnks) where perms are the disjoint column sets.
"""
function _partition_disjoint_systematic_sets(G::CTMatrixTypes)
    k, n = size(G)
    remaining_cols = collect(1:n)
    gen_mats = []
    perms = []
    rnks = Int[]

    while !isempty(remaining_cols)
        # 1. Work only on the subset of columns not yet assigned
        G_rem = G[:, remaining_cols]
        
        rnk, _, P = _rref_col_swap(G_rem, 1:k, 1:size(G_rem, 2))
        
        if rnk == 0 break end
        
        # Guard against missing permutation matrices if no swaps were needed
        if ismissing(P)
            P = identity_matrix(base_ring(G), size(G_rem, 2))
        end
        
        # 2. Extract pivot indices
        pivot_indices_in_rem = Int[]
        for i in 1:rnk
            for j in 1:size(P, 2)
                if !iszero(P[j, i]) 
                    push!(pivot_indices_in_rem, j)
                    break
                end
            end
        end
        
        set_indices = remaining_cols[pivot_indices_in_rem]
        push!(rnks, rnk)
        
        # 3. Create the Systematic Generator Matrix
        other_cols = setdiff(1:n, set_indices)
        σ = [set_indices; other_cols] 
        Gp = G[:, σ]
        
        # ROOT FIX:
        # Since σ already forced the independent pivot columns to the front, 
        # standard RREF will perfectly create [I_rnk | A] without needing to 
        # swap any columns, ensuring σ and Gp remain perfectly synced.
        _, Gp_rref = rref(Gp) 
        Gp .= Gp_rref # Mutate Gp in place to match the reduced form
        
        push!(gen_mats, Gp)
        push!(perms, σ) 

        # 4. Enforce DISJOINTNESS
        filter!(x -> x ∉ set_indices, remaining_cols)
    end
    
    return gen_mats, perms, rnks
end

# actually had to use AI to make this less good because it stayed too balanced and highly structured codes are most often unbalanced
function _greedy_increment_bb21(a_values, rnks, q)
    s = length(a_values)
    costs = fill(Inf, s + 1)
    
    # Calculate total weight currently searched
    total_w = sum(a_values)
    
    for i in 1:s
        base_cost = Float64(_calculate_U_size(rnks[i], a_values[i] + 1, q))
        
        # HEURISTIC: Penalize states that are "too balanced" 
        # to encourage exploring the lopsided distributions 
        # common in cyclic/BCH codes.
        avg_w = total_w / s
        imbalance = abs((a_values[i] + 1) - avg_w)
        
        # Reduce the cost of moves that increase imbalance
        # This makes [15, 6] look more attractive relative to [11, 10]
        costs[i] = base_cost / (1.0 + 0.5 * imbalance)
    end
    
    # DEBUG: See how the bias shifts the selection
    # println("Biased costs: $costs")
    
    best_idx = argmin(costs)
    new_a = copy(a_values)
    if best_idx <= s
        new_a[best_idx] += 1
    end

    return new_a
end

"""
    _calculate_U_size(rnk::Int, a::Int, q::Int)

Calculates the number of non-proportional codewords in U_i^(a).
Based on the formulas in Section 3 of Bouyuklieva & Bouyukliev (2021).
"""
function _calculate_U_size(rnk::Int, a::Int, q::Int)
    if a == 0
        # For reduced sets (i > t), |U_i^(0)| = (q^(k-rnk) - 1) / (q-1)
        # However, for the greedy selector, we only care about new combinations.
        return 0
    end

    # sum_{j=1}^{a} (q-1)^(j-1) * binomial(rnk, j)
    # Note: We divide the paper's formula by (q-1) to get non-proportional vectors
    total = big(0)
    for j in 1:a
        term = binomial(big(rnk), big(j)) * big(q - 1)^(j - 1)
        total += term
    end
    
    return total
end

"""
    _is_new_codeword(v, systematic_sets, a_values, current_idx)

Returns true if the codeword 'v' has not been generated in a previous 
step of the BB21 algorithm. [cite: 934]
"""
function _is_new_codeword(v, systematic_sets, a_values, current_idx)
    for j in 1:(current_idx - 1)
        # wt(v|Tj) < aj check ensures uniqueness [cite: 934]
        if count(!iszero, v[systematic_sets[j]]) < a_values[j]
            return false
        end
    end
    return true
end

"""
    heuristic_info_set_selection(C::AbstractLinearCode)

Selects the most mathematically efficient information set algorithm 
based on the code's length, dimension, and algebraic structure.
"""
function heuristic_info_set_selection(C::AbstractLinearCode)    
    # 1. BB21 Priority: Length not divisible by dimension
    # Significant reduction in codewords for n != tk.
    if C.n % C.k != 0
        return :Bouyuklieva
    end

    # 2. White Priority: Quasi-Cyclic Structure
    # If the code is QC, White's expansion-factor sets are usually superior.
    if isa(C, QuasiCyclicCode)
        return :White
    end

    # 3. Chen Priority: Cyclic Structure
    # Optimized for standard RREF forms in cyclic codes.
    if isa(C, CyclicCode)
        return :Chen
    end

    # 4. Brouwer Priority: High Symmetry
    # If there's a large automorphism group, Brouwer's disjoint blocks 
    # exploit the symmetry for faster coverage.
    if length(_generate_known_automorphisms(C)) > 100
        return :Brouwer
    end

    # 5. Default: Zimmermann
    # For random codes where n = tk, overlapping sets with scoring 
    # are the most robust fallback.
    return :Zimmermann
end

# Convert a Permutation Matrix to a Vector of indices
function _matrix_to_perm_vector(P::CTMatrixTypes)
    n = size(P, 1)
    p_vec = zeros(Int, n)
    for r in 1:n
        for c in 1:n
            if P[r, c] != 0  # Find the '1' in this row
                p_vec[r] = c
                break
            end
        end
    end
    return p_vec
end

#############################
# Non-Binary Helper Functions
#############################

"""
    _pack_field_elem(x, p::Int, d::Int)

Converts a finite field element (prime or extension) into a unique integer.
Treats the polynomial coefficients of `x` as a base-`p` integer.
"""
function _pack_field_elem(x, p::Int, d::Int)
    if d == 1
        # Prime field: simply lift the element to an integer
        return Int(lift(x))
    else
        # Extension field: pack coefficients into a base-p integer
        val = 0
        for i in 0:(d-1)
            # lift(coeff(x, i)) grabs the i-th coefficient and makes it a native Int
            val += Int(lift(coeff(x, i))) * (p^i)
        end
        return val
    end
end

"""
    _unpack_field_elem(val::Int, F, p::Int, d::Int)

Reconstructs an Oscar finite field element from its integer representation.
"""
function _unpack_field_elem(val::Int, F, p::Int, d::Int)
    if d == 1
        return F(val)
    else
        g = gen(F) # The primitive element/generator of the extension field
        res = zero(F)
        v = val
        for i in 0:(d-1)
            c = v % p
            res += F(c) * (g^i)
            v = div(v, p)
        end
        return res
    end
end

"""
    _precompute_pruning_bounds_nonbinary(A_idx::Matrix{Int}, k::Int, q::Int)

Calculates the absolute lower bound (`lbt`) and worst-case upper bound (`max_canc`) 
for the parity weight generated by 'r' remaining message symbols. 
Expects an integer-mapped matrix where index 1 represents the zero element.
"""
function _precompute_pruning_bounds_nonbinary(A_idx::Matrix{Int}, k::Int, q::Int)
    n_tail = size(A_idx, 1)
    
    lbt = zeros(Int, k + 1)
    max_canc = zeros(Int, k + 1)
    
    # CHANGED: Index 1 is the zero element, so we count entries that are NOT 1.
    col_weights = [count(!=(1), view(A_idx, :, j)) for j in 1:k]
    
    # --- LOWER BOUND TABLE (lbt) ---
    lbt[1] = 0
    lbt[2] = minimum(col_weights) 
    
    for r in 2:k
        # GF(q) Griesmer bound
        low, high, griesmer = 0, n_tail, 0
        while low <= high
            mid = (low + high) ÷ 2
            g_sum, d_val = 0, Float64(mid)
            for i in 1:r
                term = ceil(Int, d_val)
                g_sum += term
                d_val /= q # q-ary division
                term <= 1 && (g_sum += (r - i); break)
            end
            g_sum <= n_tail ? (griesmer = mid; low = mid + 1) : (high = mid - 1)
        end
        lbt[r + 1] = griesmer
    end
    
    # --- MAX CANCELLATION TABLE (max_canc) ---
    sort!(col_weights, rev=true) 
    
    max_canc[1] = 0
    current_max_sum = 0
    for r in 1:k
        current_max_sum += col_weights[r]
        max_canc[r + 1] = min(n_tail, current_max_sum)
    end
    
    return lbt, max_canc
end

"""
    _precompute_weight2_table_nonbinary(A_idx::Matrix{Int}, add_t::Matrix{Int}, mul_t::Matrix{Int}, q::Int)

Precomputes the absolute minimum parity weight generated by any pair of remaining 
message symbols over GF(q). This is the non-binary generalization of the "Brouwer Trick", 
allowing the recursive search engine to replace deep tree branching with an O(1) 
table lookup for the final two symbols.

# Mathematical Optimization
Instead of naively checking all (q-1)^2 possible non-zero combinations of two 
message symbols (m_i * col_i + m_j * col_j), this function exploits the fact that 
Hamming weight is invariant under scalar multiplication. It calculates the minimum 
weight of (col_i + γ * col_j) for all non-zero scalars γ ∈ GF(q), reducing the 
complexity from O(k^2 * q^2) down to O(k^2 * q).

# Arguments
- `A_idx::Matrix{Int}`: The (n-k) × k parity tail matrix. Entries must be 1-based integer indices corresponding to field elements, NOT the raw field elements themselves.
- `add_t::Matrix{Int}`: The q × q addition lookup table for the field.
- `mul_t::Matrix{Int}`: The q × q multiplication lookup table for the field.
- `q::Int`: The order of the finite field.

!!! warning "Zero-Index Assumption"
    This function strictly assumes that **Index 1** represents the additive 
    identity (the zero element) in both the `add_t` and `mul_t` lookup tables.

# Returns
- A `k × k` symmetric integer matrix where entry `[i, j]` contains the absolute 
  minimum weight generated by any non-zero linear combination of columns `i` and `j`.
"""
function _precompute_weight2_table_nonbinary(A_idx::Matrix{Int}, add_t::Matrix{Int}, mul_t::Matrix{Int}, q::Int)
    n_tail, k = size(A_idx)
    w2_min = fill(typemax(Int), k, k)
    
    for i in 1:k
        for j in i+1:k
            min_w = typemax(Int)
            
            # 2:q assumes index 1 is the additive identity (0)
            for γ in 2:q 
                w = 0
                @inbounds @simd for c in 1:n_tail
                    val_j = mul_t[γ, A_idx[c, j]]
                    val_sum = add_t[A_idx[c, i], val_j]
                    
                    # BRANCHLESS: Booleans evaluate to 1 (true) or 0 (false)
                    w += (val_sum != 1) 
                end
                
                if w < min_w
                    min_w = w
                end
            end
            w2_min[i, j] = min_w
            w2_min[j, i] = min_w
        end
    end

    return w2_min
end

function _make_systematic_gf!(M::CTMatrixTypes, perm::Vector{Int}, k::Int)
    n = size(M, 2)
    
    for i in 1:k
        # 1. Find a pivot in column i (from row i downwards)
        pivot_row = 0
        for r in i:k
            if !iszero(M[r, i])
                pivot_row = r
                break
            end
        end
        
        # 2. If no pivot is found, the information set is singular. 
        # We must swap a column from the parity section (j > k) into position i.
        if pivot_row == 0
            found_col = false
            for c in k+1:n
                for r in i:k
                    if !iszero(M[r, c])
                        # Swap columns in the matrix
                        for row in 1:k
                            M[row, i], M[row, c] = M[row, c], M[row, i]
                        end
                        # Track the permutation!
                        perm[i], perm[c] = perm[c], perm[i]
                        pivot_row = r
                        found_col = true
                        break
                    end
                end
                found_col && break
            end
            
            # If still no pivot, the entire matrix doesn't have rank k. 
            # (Should never happen for a valid generator matrix).
            @assert found_col "Matrix does not have full row rank!"
        end
        
        # 3. Swap the pivot row to the current row i
        if pivot_row != i
            for c in 1:n
                M[i, c], M[pivot_row, c] = M[pivot_row, c], M[i, c]
            end
        end
        
        # 4. Normalize the pivot row (make the diagonal element 1)
        # Assumes inv() is overloaded for your GF(q) elements
        pivot_inv = inv(M[i, i]) 
        for c in i:n
            M[i, c] = M[i, c] * pivot_inv
        end
        
        # 5. Eliminate all other entries in column i
        for r in 1:k
            if r != i && !iszero(M[r, i])
                factor = M[r, i]
                for c in i:n
                    M[r, c] = M[r, c] - factor * M[i, c]
                end
            end
        end
    end
end

"""
    _generate_scored_zimmermann_mats_nonbinary(G::Matrix, m::Int; pool_size::Int = 20)

Generates `m` highly-optimized Zimmermann matrices for the GF(q) BZ algorithm.
Explores a pool of random column permutations and selects the ones that yield 
the sparsest parity matrices (lowest Hamming weight) after reduction.
"""
function _generate_scored_zimmermann_mats_nonbinary(G::Matrix, m::Int; pool_size::Int = 20)
    k, n = size(G)
    
    # We always want the original, unpermuted standard form as our first matrix
    z_mats = []
    
    # Generate a pool of candidates
    candidates = []
    for _ in 1:pool_size
        # Create a random permutation
        perm = collect(1:n)
        shuffle!(perm)
        
        # Apply permutation to G
        Gp = G[:, perm]
        
        # Put into systematic form [I_k | P]
        _make_systematic_gf!(Gp, perm, k)
        
        # Score it! (Count non-zero elements in the parity tail P)
        # BZ runs exponentially faster if the parity matrix is sparse.
        tail_weight = count(!iszero, view(Gp, :, k+1:n))
        
        push!(candidates, (score = tail_weight, G = Gp, perm = perm))
    end
    
    # Sort candidates by score (lowest tail weight first)
    sort!(candidates, by = x -> x.score)
    
    # Select the best 'm' matrices
    # (Or take all of them if pool_size was smaller than m)
    num_to_take = min(m, length(candidates))
    for i in 1:num_to_take
        # Format it exactly as your master function expects
        push!(z_mats, (G = candidates[i].G, perm = candidates[i].perm))
    end
    
    return z_mats
end

"""
    _pack_matrix_gf3_bitsliced(A_raw::Matrix{Int})

Compresses a GF(3) matrix into orthogonal bitsliced UInt64 matrices.
Returns a tuple (H_packed, L_packed).
"""
function _pack_matrix_gf3_bitsliced(A_raw::Matrix{Int})
    n_tail, k = size(A_raw)
    num_chunks = cld(n_tail, 64)
    
    H_packed = zeros(UInt64, num_chunks, k)
    L_packed = zeros(UInt64, num_chunks, k)
    
    for j in 1:k
        for i in 1:n_tail
            val = A_raw[i, j]
            if val == 1
                chunk_idx = (i - 1) ÷ 64 + 1
                bit_idx = (i - 1) % 64
                L_packed[chunk_idx, j] |= (UInt64(1) << bit_idx)
            elseif val == 2
                chunk_idx = (i - 1) ÷ 64 + 1
                bit_idx = (i - 1) % 64
                H_packed[chunk_idx, j] |= (UInt64(1) << bit_idx)
            end
        end
    end
    
    return H_packed, L_packed
end

"""
    _fast_simd_wt_gf3(H::Vector{UInt64}, L::Vector{UInt64})

Calculates the Hamming weight of a bitsliced GF(3) parity vector.
"""
@inline function _fast_simd_wt_gf3(H::Vector{UInt64}, L::Vector{UInt64})
    w = 0
    @inbounds @simd for i in eachindex(H)
        w += count_ones(H[i] | L[i])
    end
    return w
end

"""
    _pack_matrix_gf4_bitsliced(A_raw::Matrix{Int})

Compresses a GF(4) matrix into orthogonal bitsliced UInt64 matrices.
Assumes field elements are mapped as 0=>0, 1=>1, 2=>ω, 3=>ω².
"""
function _pack_matrix_gf4_bitsliced(A_raw::Matrix{Int})
    n_tail, k = size(A_raw)
    num_chunks = cld(n_tail, 64)
    
    H_packed = zeros(UInt64, num_chunks, k)
    L_packed = zeros(UInt64, num_chunks, k)
    
    for j in 1:k
        for i in 1:n_tail
            val = A_raw[i, j]
            chunk_idx = (i - 1) ÷ 64 + 1
            bit_idx = (i - 1) % 64
            
            # Extract the lowest bit for L, and the second bit for H
            if (val & 1) != 0
                L_packed[chunk_idx, j] |= (UInt64(1) << bit_idx)
            end
            if (val & 2) != 0
                H_packed[chunk_idx, j] |= (UInt64(1) << bit_idx)
            end
        end
    end
    
    return H_packed, L_packed
end

"""
    _add_scaled_gf4_simd!(tail_H, tail_L, col_H, col_L, scalar::Int)

Adds a column scaled by a GF(4) element (1, 2=ω, 3=ω²) to the tail.
"""
@inline function _add_scaled_gf4_simd!(
    tail_H::Vector{UInt64}, tail_L::Vector{UInt64}, 
    col_H::AbstractVector{UInt64}, col_L::AbstractVector{UInt64},
    scalar::Int
)
    @inbounds @simd for i in eachindex(tail_H)
        CH, CL = col_H[i], col_L[i]
        
        if scalar == 1      # Multiply by 1
            tail_H[i] ⊻= CH
            tail_L[i] ⊻= CL
        elseif scalar == 2  # Multiply by ω
            tail_H[i] ⊻= (CH ⊻ CL)
            tail_L[i] ⊻= CH
        elseif scalar == 3  # Multiply by ω²
            tail_H[i] ⊻= CL
            tail_L[i] ⊻= (CH ⊻ CL)
        end
    end
end

"""
    _precompute_weight2_table_gf3(H::Matrix{UInt64}, L::Matrix{UInt64})

Precomputes the GF(3) Brouwer weight-2 table directly from the bitsliced matrices.
Uses the branchless boolean circuit to evaluate combinations in SIMD.
"""
function _precompute_weight2_table_gf3(H::Matrix{UInt64}, L::Matrix{UInt64})
    num_chunks, k = size(H)
    w2_min = fill(typemax(Int), k, k)
    
    for i in 1:k
        for j in i+1:k
            min_w = typemax(Int)
            
            # Gamma = 1 and Gamma = 2
            for γ in 1:2
                w = 0
                @inbounds @simd for c in 1:num_chunks
                    AH, AL = H[c, i], L[c, i]
                    
                    # If γ == 1, use col j. If γ == 2, swap H and L of col j.
                    BH = γ == 1 ? H[c, j] : L[c, j]
                    BL = γ == 1 ? L[c, j] : H[c, j]
                    
                    # GF(3) Boolean Adder Circuit
                    SL = AL ⊻ BL
                    SH = AH ⊻ BH
                    XL = SL ⊻ (AH & BH)
                    XH = SH ⊻ (AL & BL)
                    mask = ~(XL & XH)
                    
                    SumH = XH & mask
                    SumL = XL & mask
                    
                    w += count_ones(SumH | SumL)
                end
                if w < min_w
                    min_w = w
                end
            end
            w2_min[i, j] = min_w
            w2_min[j, i] = min_w
        end
    end
    return w2_min
end

"""
    _precompute_weight2_table_gf4(H::Matrix{UInt64}, L::Matrix{UInt64})

Precomputes the GF(4) Brouwer weight-2 table directly from the bitsliced matrices.
Evaluates col_i + γ * col_j for γ ∈ {1, ω, ω²}.
"""
function _precompute_weight2_table_gf4(H::Matrix{UInt64}, L::Matrix{UInt64})
    num_chunks, k = size(H)
    w2_min = fill(typemax(Int), k, k)
    
    for i in 1:k
        for j in i+1:k
            min_w = typemax(Int)
            
            # γ = 1 (1), γ = 2 (ω), γ = 3 (ω²)
            for γ in 1:3
                w = 0
                @inbounds @simd for c in 1:num_chunks
                    AH, AL = H[c, i], L[c, i]
                    BH, BL = H[c, j], L[c, j]
                    
                    # Apply GF(4) scalar multiplication
                    if γ == 1
                        CH, CL = BH, BL
                    elseif γ == 2
                        CH, CL = (BH ⊻ BL), BH
                    else # γ == 3
                        CH, CL = BL, (BH ⊻ BL)
                    end
                    
                    # Addition is just XOR
                    SumH = AH ⊻ CH
                    SumL = AL ⊻ CL
                    
                    w += count_ones(SumH | SumL)
                end
                if w < min_w
                    min_w = w
                end
            end
            w2_min[i, j] = min_w
            w2_min[j, i] = min_w
        end
    end
    return w2_min
end

"""
    _generate_field_tables_safe(F, q::Int)

Safely generates 1-based integer lookup tables for any generic field representation.
"""
function _generate_field_tables_safe(F, q::Int)
    # Most Julia algebra libraries support collect() for small finite fields.
    elements = collect(F) 
    
    # Safely find the additive identity (zero)
    zero_idx = findfirst(iszero, elements)
    if zero_idx != 1
        elements[1], elements[zero_idx] = elements[zero_idx], elements[1]
    end
    
    # Map the opaque field element to a standard 1-based Int
    elem_to_idx = Dict(elements[i] => i for i in 1:q)
    
    add_t = zeros(Int, q, q)
    mul_t = zeros(Int, q, q)
    
    for i in 1:q
        for j in 1:q
            # The dictionary safely handles whatever custom type the library uses
            add_t[i, j] = elem_to_idx[elements[i] + elements[j]]
            mul_t[i, j] = elem_to_idx[elements[i] * elements[j]]
        end
    end
    
    return add_t, mul_t, elem_to_idx
end

"""
    _add_gf3_simd!(tail_H, tail_L, col_H, col_L)

Performs parallel GF(3) addition of a column into the current tail using 
a branchless boolean logic circuit.
"""
@inline function _add_gf3_simd!(
    tail_H::Vector{UInt64}, tail_L::Vector{UInt64}, 
    col_H::AbstractVector{UInt64}, col_L::AbstractVector{UInt64}
)
    @inbounds @simd for i in eachindex(tail_H)
        AH, AL = tail_H[i], tail_L[i]
        BH, BL = col_H[i], col_L[i]

        # Standard binary addition (XOR)
        SL = AL ⊻ BL
        SH = AH ⊻ BH
        
        # Cross-carries
        XL = SL ⊻ (AH & BH)
        XH = SH ⊻ (AL & BL)
        
        # Annihilation mask (1+2=0)
        mask = ~(XL & XH)
        
        tail_H[i] = XH & mask
        tail_L[i] = XL & mask
    end
end

# #############################
#   # Enumeration Based Algs
# #############################

# """
#     minimum_distance_Gray(C::AbstractLinearCode; alg::Symbol = :Zimmermann, v::Bool = false)
# Return the minimum distance of `C` using a deterministic algorithm based on enumerating
# constant weight codewords of the binary reflected Gray code. If a word of minimum weight
# is found before the lower and upper bounds cross, it is returned; otherwise, the zero vector 
# is returned.

# show_progress will display a progress meter for each iteration of a weight that takes longer than 
#     a second
# """
# function minimum_distance_Gray(C::AbstractLinearCode; alg::Symbol = :auto, v::Bool = false, 
#     show_progress = true)

#     ord_F = Int(order(C.F))
#     ord_F == 2 || throw(ArgumentError("Currently only implemented for binary codes."))
#     #TODO add :Bouyuklieva, and :Edmonds

#     # We implement the following algorithms described in White's thesis:
#     # :Brouwer      Algo 2.2
#     # :Zimmermann   Algo 2.4 
#     # :Chen         Algo 2.6 
#     # :White        Algo 3.1 
#     alg ∈ (:auto, :Brouwer, :Zimmermann, :White, :Chen) || throw(ArgumentError("Unknown information set algorithm. Expected `:auto`, `:Brouwer`, `:Zimmermann`, `:White`, `:Chen`"))

#     if alg == :auto
#         if typeof(C) <: AbstractCyclicCode
#             v && println("Detected a cyclic code, using Chen's adaption.")
#             alg = :Chen
#             # TODO: fix this case
#         elseif typeof(C) <: AbstractQuasiCyclicCode
#             v && println("Detected a quasi-cyclic code, using White's adaption.")
#             alg = :White
#         else
#             v && println("Using Zimmermann's algorithm.")
#             alg = :Zimmermann
#         end
#     end
#     alg == :auto && throw(ErrorException("Could not determine minimum distance algo automatically"))

#     if alg in (:Brouwer, :Zimmermann) 
#         return _minimum_distance_BZ(C::AbstractLinearCode; info_set_alg = alg, verbose = v, show_progress = show_progress)
#     end
#     println("Warning: old enumeration algorithm selected. Performance will be slow") # TODO remove when all code updated
#     return _minimum_distance_enumeration_with_matrix_multiply(C::AbstractLinearCode; info_set_alg = alg)
# end

# function _minimum_distance_BZ(C::AbstractLinearCode; info_set_alg::Symbol = :Zimmermann,
#     verbose::Bool = false, dbg = Dict(), show_progress=false)
#     dbg_key_exit_r = "exit_r"

#     ord_F = Int(order(C.F))
#     ord_F == 2 || throw(ArgumentError("Currently only implemented for binary codes."))
#     C.k < 2^16 || throw(DomainError("The given linear code has length k >= 2^16 which is not supported"))
#     info_set_alg ∈ (:Brouwer, :Zimmermann) || throw(ArgumentError("Unknown information set algorithm. Expected `:Brouwer`, `:Zimmermann`"))

#     generator_matrix(C, true) # ensure G_stand exists
#     if _has_empty_vec(C.G, :cols) 
#         #TODO err string can instruct the user to construct a new code without 0 cols and tell them the function for that
#         throw(ArgumentError("Codes with standard form of generator matrix having 0 columns not supported")) 
#     end
#     # generate if not pre-stored
#     parity_check_matrix(C)

#     A_mats, perms_mats, rnks = information_sets(C.G, info_set_alg, permute = true, only_A = false)

#     A_mats = [deepcopy(_Flint_matrix_to_Julia_T_matrix(Ai, UInt16)') for Ai in A_mats]
#     # A_mats_trunc = () 
#     perms_mats = [deepcopy(_Flint_matrix_to_Julia_T_matrix(Pi, UInt16)') for Pi in perms_mats]
#     h = length(A_mats)
#     # println("Starting loop to refine upper bound. Initial upper bound ", C.u_bound, " num of mats is ", length(A_mats), " dimension ", size(A_mats[1]))
#     rank_defs = zeros(Int, h)

#     if length(keys(dbg)) > 0
#         println("Debug mode ON")
#     end

#     if haskey(dbg, dbg_key_exit_r)
#         verbose && println("dbg Dict: largest message weight searched stored @key=$dbg_key_exit_r")
#         dbg[dbg_key_exit_r] = -1
#     end

#     k, n = size(C.G)
#     A_mats_trunc = [Matrix{UInt16}(undef, k, n-k) for _ in 1:length(A_mats)]
#     for i in 1:size(A_mats, 1) 
#        A_mats_trunc[i] = deepcopy(A_mats[i][k+1 : n, :])
#     end

#     if info_set_alg == :Brouwer && rnks[h] != k
#         println("Rank of last matrix too small")
#         return
#     end
#     if verbose
#         print("Generated $h information sets with ranks: ")
#         for i in 1:h
#             i == h ? (println(rnks[i]);) : (print("$(rnks[i]), ");)
#             # will only be using the rank deficits here
#             # at the moment, the information sets are always disjoint so the relative
#             # rank is zero
#             # TODO huh? check this comment and setup properly
#             rank_defs[i] = C.k - rnks[i]
#         end
#     end
    
#     even_flag = false
#     doubly_even_flag = false
#     triply_even_flag = false
#     ord_F == 2 && (even_flag = is_even(C);)
#     even_flag && (doubly_even_flag = is_doubly_even(C);)
#     doubly_even_flag && (triply_even_flag = is_triply_even(C);)
#     if verbose
#         triply_even_flag && println("Detected a triply even code.")
#         (!triply_even_flag && doubly_even_flag) && println("Detected a doubly even code.")
#         (!triply_even_flag && !doubly_even_flag && even_flag) && println("Detected an even code.")
#     end

#     # initial_perm_ind will match the permutation we use for the 'found' vector if the found vector is nonzero. To simplify the code below we're going to choose an initial permutation arbitrarily.  
#     initial_perm_ind = 1 
#     # following loop is the r=1 case of the enumeration. We do this case here because we want to make a good guess at the terminating r before we start multiple threads
#     for (j, g) in enumerate(A_mats) # loop over the A_mats rather than the original G because it would add another case to deal with later 
#         # can make this faster with dots and views
#         w, i = _min_wt_col(g)
#         if w <= C.u_bound
#             found = g[:, i] 
#             C.u_bound = w
#             y = perms_mats[j] * found 
#         end
#     end

#     verbose && println("Current upper bound: $(C.u_bound)")
#     found = A_mats[1][:, 1]

#     l = 0
#     if verbose 
#         _, _, b_rnks = information_sets(C.G, :Brouwer, permute = true, only_A = false)
#         b_h = length(b_rnks)
#         b_lower_bounds = [_information_set_lower_bound(r+1, n, k, l, [k - 0 for i in 1:b_h], :Brouwer, even = even_flag, doubly_even = doubly_even_flag, triply_even = triply_even_flag) for r in 1:k-1]
#         b_r_term = findfirst(x -> x ≥ C.u_bound, b_lower_bounds)

#         # _, _, z_rnks = information_sets(G, :Zimmermann, permute = true, only_A = false)
#         # z_h = length(b_rnks)
#         # z_lower_bounds = [_information_set_lower_bound(r+1, n, k, l, [k - z_rnks[i] for i in 1:z_h], :Zimmermann, even = even_flag, doubly_even = doubly_even_flag, triply_even = triply_even_flag) for r in 1:k-1]
#         # z_r_term = findfirst(x -> x ≥ C.u_bound, z_lower_bounds)
#         # verbose && println("ranks: Brouwer $b_rnks Zimm $z_rnks")
#         # verbose && println("Predicted termination weight based on current upper bound: Brouwer $b_r_term Zimm $z_r_term")
#     end

#     #Note the r+1 here. 
#     lower_bounds_for_prediction = [_information_set_lower_bound(r+1, n, k, l, rank_defs, info_set_alg, even = even_flag, doubly_even = doubly_even_flag, triply_even = triply_even_flag) for r in 1:k-1]
#     r_term = findfirst(x -> x ≥ C.u_bound, lower_bounds_for_prediction)
#     if isnothing(r_term)
#         raise(DomainError("invalid termination r")) 
#     end
#     verbose && println("Predicted termination weight based on current upper bound: $r_term")

#     #In the main loop we check if lower bound > upper bound before we enumerate and so the lower bounds for the loop use r not r+1
#     lower_bounds = [_information_set_lower_bound(r, n, k, l, rank_defs, info_set_alg, even = even_flag, doubly_even = doubly_even_flag, triply_even = triply_even_flag) for r in 1:k-1]

#     predicted_work_factor = fld(n, k) * sum([binomial(k, i) for i in 1:r_term])
#     verbose && println("Predicted work factor: $predicted_work_factor")
#     if show_progress 
#         prog_bar = Progress(predicted_work_factor, dt=1.0, showspeed=true) # updates no faster than once every 1s
#     end
#     weight_sum_bound = min(2 * C.u_bound + 5, n-k)
#     verbose && println("Codeword weights initially checked on first $weight_sum_bound entries")

#     num_thrds = Threads.nthreads()
#     verbose && println("Number of threads ", num_thrds)
#     for r in 2:k
#         if r > 2^16
#             verbose && println("Warning: Reached an r larger than 2^16") 
#         end
#         C.l_bound < lower_bounds[r] && (C.l_bound = lower_bounds[r];)
#         # an even code can't have have an odd minimum weight
#         # (!triply_even_flag && !doubly_even_flag && even_flag) && (C.l_bound += C.l_bound % 2;)
#         # (!triply_even_flag && doubly_even_flag) && (C.l_bound += 4 - C.l_bound % 4;)
#         # triply_even_flag && (C.l_bound += 8 - C.l_bound % 8;)
#         if C.l_bound >= C.u_bound
#             dbg[dbg_key_exit_r] = r-1
#             break
#         end
#         verbose && println("r: $r")
#         verbose && println("Lower bound: $(C.l_bound)")
#         verbose && println("Upper bound: $(C.u_bound)")

#         if verbose
#             i_count = 0
#             for i in 1:h
#                 r - rank_defs[i] ≤ 0 && (i_count += 1;)
#             end
#             i_count > 0 && println("$i_count of the original $h information sets no longer contribute to the lower bound")
#         end
#         p = Int(characteristic(C.F))

#         uppers = [C.u_bound for _ in 1:num_thrds]
#         founds = [found for _ in 1:num_thrds]
#         exit_thread_indicator_vec = [initial_perm_ind for _ in 1:num_thrds]
#         keep_going = Threads.Atomic{Bool}(true)

#         bin = extended_binomial(C.k, r)

#         thrd_stop_msg = "Stopping current thread, main loop finished"

#         Threads.@threads for ind in 1:num_thrds 
#             len = (ind == num_thrds) ? bin - (num_thrds - 1) * fld(bin, num_thrds) : fld(bin, num_thrds)

#             # iteration begins with a single matrix multiplication of the generator matrix by first_vec
#             init_rank = 1 + (ind - 1) * fld(bin, num_thrds)
#             first_vec = zeros(Int, k)
#             if init_rank == 1
#                 for i in 1:r
#                     first_vec[i] = 1
#                 end
#             else
#                 CodingTheory._subset_unrank_to_vec!(init_rank, UInt64(r), first_vec)
#             end
 
#             # as in White Algo 7.1 we loop over matrices first 
#             for i in 1:h
#                 if keep_going[] == false
#                     verbose && println(thrd_stop_msg)
#                     break
#                 end

#                 c_itr = zeros(UInt16, C.n - C.k) 
#                 is_first = true
#                 curr_mat = A_mats_trunc[i]
#                 count = UInt128(0)

#                 for u in SubsetGrayCode(k, r, len, init_rank)
#                     if keep_going[] == false
#                         println(thrd_stop_msg)
#                         break
#                     end
#                     show_progress && ProgressMeter.next!(prog_bar) 
#                     if r - rank_defs[i] > 0
#                         if is_first 
#                             LinearAlgebra.mul!(c_itr, curr_mat, first_vec)
#                             @inbounds @simd for j in eachindex(c_itr) 
#                                 c_itr[j] %= p
#                             end
#                             is_first = false
#                         else
#                             for ci in u 
#                                 if ci != -1
#                                     @simd for i in eachindex(c_itr)
#                                         @inbounds c_itr[i] = xor(c_itr[i], curr_mat[i, ci])
#                                     end
#                                 end
#                             end
#                         end

#                         partial_weight = r + sum(view(c_itr, 1:weight_sum_bound))

#                         if uppers[ind] > partial_weight
#                             w = r + sum(c_itr) 
#                             verbose && @assert w != 0
#                             if uppers[ind] > w 
#                                 subset_vec_full = zeros(Int, k)
#                                 CodingTheory._subset_unrank_to_vec!(UInt128(init_rank + count), UInt64(r), subset_vec_full)

#                                 uppers[ind] = w 
#                                 founds[ind] = vcat(subset_vec_full, c_itr)
#                                 verbose && @assert size(founds[ind], 1) == C.n "found vector has length $(size(founds[ind], 1)) but should be n=$(C.n)"
#                                 exit_thread_indicator_vec[ind] = i 

#                                 println("Adjusting (local) upper bound: $w for c_itr=$(Int.(c_itr))")
#                                 if C.l_bound == uppers[ind]
#                                     println("early exit")
#                                     Threads.atomic_cas!(keep_going, true, false)
#                                 else
#                                     r_term = findfirst(x -> x ≥ C.u_bound, lower_bounds)
#                                     isnothing(r_term) && (r_term = k;)
#                                     verbose && println("Updated termination weight: $r_term")
#                                 end
#                             end
#                         end
#                     end
#                     count = add!(count, count, 1)
#                 end
#             end
#             loc = argmin(uppers) 
#             C.u_bound = uppers[loc]
#             found = founds[loc]
#             initial_perm_ind = exit_thread_indicator_vec[loc]
#         end
#     end

#     C.d = C.u_bound
#     y = matrix(C.F, 1, n, perms_mats[initial_perm_ind] * found) # weight(y) >= C.d, with equality not being the typical case
#     verbose && @assert iszero(C.H * transpose(y))
#     if dbg[dbg_key_exit_r] == -1
#         dbg[dbg_key_exit_r] = r
#     end 
#     show_progress && ProgressMeter.finish!(prog_bar)
#     verbose && println("Computation complete")
#     return C.u_bound, y
# end

#############################
     # Minimum Distance
#############################

# the recursion will never cause a stack overflow. The depth is strictly bounded by k, 
# and in practice, it usually terminates far shallower than k due to the Griesmer/Cancellation pruning, 
# the Automorphism pruning, and the Pigeonhole bound.
function _Brouwer_Zimmermann_binary_recursive!(
    A_packed::Vector{Vector{UInt64}},  # Vector of length k (each element is a row of chunks)
    r::Int, depth::Int, picked::Int, 
    curr_tail::Vector{UInt64},         # The current XOR sum of parity bits
    best_w::Threads.Atomic{Int}, best_msg::Vector{Int}, update_lock::Threads.SpinLock,
    lbt::Vector{Int}, max_canc::Vector{Int}, w2_table::Matrix{Int},
    keep_going::Threads.Atomic{Bool}, l_bound::Int, current_msg::Vector{Int},
    spawn_depth::Int, auts::Vector{Vector{Int}}
)
    !keep_going[] && return
    
    k = length(A_packed)
    tw = 0
    @inbounds @simd for c in eachindex(curr_tail)
        tw += count_ones(curr_tail[c])
    end
    
    # Base Case
    if picked == r
        if !_is_canonical(current_msg, auts)
            return
        end

        w = r + tw
        if w < best_w[] 
            lock(update_lock) do
                if w < best_w[]
                    Threads.atomic_xchg!(best_w, w)
                    copyto!(best_msg, current_msg)
                    for i in (depth + 1):length(best_msg)
                        best_msg[i] = 0
                    end
                    if w <= l_bound 
                        Threads.atomic_cas!(keep_going, true, false)
                    end
                end
            end
        end
        return
    end

    # Structural Pruning
    rem_to_pick = r - picked
    if depth >= k || (k - depth) < rem_to_pick
        return
    end

    if rem_to_pick == 2 && depth < k - 1
        min_w2 = typemax(Int)
        for i in (depth + 1):k
            for j in (i + 1):k
                combined_parity_wt = 0
                @inbounds @simd for c in eachindex(curr_tail)
                    combined_parity_wt += count_ones(curr_tail[c] ⊻ A_packed[i][c] ⊻ A_packed[j][c])
                end
                if (r + combined_parity_wt) < min_w2
                    min_w2 = r + combined_parity_wt
                end
            end
        end
        if min_w2 >= best_w[]
            return 
        end
    end

    # Weight-based Pruning
    min_rem_wt = lbt[rem_to_pick + 1]
    max_rem_wt = max_canc[rem_to_pick + 1]
    
    if tw < min_rem_wt
        min_possible_tw = min_rem_wt - tw
    elseif tw > max_rem_wt
        min_possible_tw = tw - max_rem_wt
    else
        min_possible_tw = 0
    end
    
    if (r + min_possible_tw) >= best_w[] 
        return 
    end

    # Branching
    if depth < spawn_depth
        # Parallel Branch: Include Row
        msg_inc = copy(current_msg)
        msg_inc[depth + 1] = 1
        tail_inc = copy(curr_tail)
        @inbounds @simd for c in eachindex(tail_inc)
            tail_inc[c] ⊻= A_packed[depth + 1][c]
        end
        
        t = Threads.@spawn _Brouwer_Zimmermann_binary_recursive!(
            A_packed, r, depth + 1, picked + 1, tail_inc, best_w, best_msg, update_lock, 
            lbt, max_canc, w2_table, keep_going, l_bound, msg_inc, spawn_depth, auts
        )
            
        # Local Branch: Exclude Row
        current_msg[depth + 1] = 0
        _Brouwer_Zimmermann_binary_recursive!(
            A_packed, r, depth + 1, picked, curr_tail, best_w, best_msg, update_lock, 
            lbt, max_canc, w2_table, keep_going, l_bound, current_msg, spawn_depth, auts
        )
        wait(t)
    else
        # Serial Branch (Backtracking)
        # Option 1: Exclude row
        current_msg[depth + 1] = 0
        _Brouwer_Zimmermann_binary_recursive!(
            A_packed, r, depth + 1, picked, curr_tail, 
            best_w, best_msg, update_lock, lbt, max_canc, w2_table, 
            keep_going, l_bound, current_msg, spawn_depth, auts
        )
        
        # Option 2: Include row
        current_msg[depth + 1] = 1
        @inbounds @simd for c in eachindex(curr_tail)
            curr_tail[c] ⊻= A_packed[depth + 1][c]
        end
        
        _Brouwer_Zimmermann_binary_recursive!(
            A_packed, r, depth + 1, picked + 1, curr_tail, 
            best_w, best_msg, update_lock, lbt, max_canc, w2_table, 
            keep_going, l_bound, current_msg, spawn_depth, auts
        )
        
        @inbounds @simd for c in eachindex(curr_tail)
            curr_tail[c] ⊻= A_packed[depth + 1][c] # Backtrack without allocating!
        end
    end
end

function _Brouwer_Zimmermann_gf3_recursive!(
    A_packed_H::Matrix{UInt64}, 
    A_packed_L::Matrix{UInt64}, 
    r::Int, depth::Int, picked::Int, 
    curr_tail_H::Vector{UInt64}, curr_tail_L::Vector{UInt64}, 
    best_w::Threads.Atomic{Int}, best_msg::Vector{Int}, 
    update_lock::Threads.SpinLock, lbt::Vector{Int}, max_canc::Vector{Int}, 
    w2_table::Matrix{Int}, keep_going::Threads.Atomic{Bool}, l_bound::Int, 
    current_msg::Vector{Int}, spawn_depth::Int, auts::Vector{Vector{Int}}
)
    if !keep_going[] return end
    k = size(A_packed_H, 2)
    tw = _fast_simd_wt_gf3(curr_tail_H, curr_tail_L)
    
    # 1. BASE CASE
    if picked == r
        # FIX: Automorphism pruning
        if !_is_canonical(current_msg, auts)
            return
        end

        w = r + tw
        if w < best_w[] 
            lock(update_lock) do
                if w < best_w[]
                    Threads.atomic_xchg!(best_w, w)
                    copyto!(best_msg, current_msg)
                    
                    # FIX: Clean up the un-picked suffix left by deep branches
                    for i in (depth + 1):length(best_msg)
                        best_msg[i] = 0
                    end
                    
                    w <= l_bound && Threads.atomic_cas!(keep_going, true, false)
                end
            end
        end
        return
    end

    # 2. PRUNING & LOOKAHEAD
    rem_to_pick = r - picked
    if depth >= k || (k - depth) < rem_to_pick return end
    min_possible_tw = tw < lbt[rem_to_pick+1] ? lbt[rem_to_pick+1] - tw : (tw > max_canc[rem_to_pick+1] ? tw - max_canc[rem_to_pick+1] : 0)
    if (r + min_possible_tw) >= best_w[] return end

    if rem_to_pick == 2 && depth < k - 1
        min_w2 = typemax(Int)
        for i in depth+1:k, j in i+1:k
            lb = abs(tw - w2_table[i,j])
            (picked + lb) < min_w2 && (min_w2 = picked + lb)
        end
        min_w2 >= best_w[] && return
    end

    col_H = view(A_packed_H, :, depth + 1)
    col_L = view(A_packed_L, :, depth + 1)

    # 3. BRANCHING
    if depth < spawn_depth
        # Spawn for each possible non-zero scalar + zero branch
        t0 = Threads.@spawn begin
            msg0 = copy(current_msg); msg0[depth+1] = 0
            _Brouwer_Zimmermann_gf3_recursive!(A_packed_H, A_packed_L, r, depth+1, picked, copy(curr_tail_H), copy(curr_tail_L), best_w, best_msg, update_lock, lbt, max_canc, w2_table, keep_going, l_bound, msg0, spawn_depth, auts)
        end
        t1 = Threads.@spawn begin
            msg1 = copy(current_msg); msg1[depth+1] = 1
            th1, tl1 = copy(curr_tail_H), copy(curr_tail_L)
            _add_gf3_simd!(th1, tl1, col_H, col_L)
            _Brouwer_Zimmermann_gf3_recursive!(A_packed_H, A_packed_L, r, depth+1, picked+1, th1, tl1, best_w, best_msg, update_lock, lbt, max_canc, w2_table, keep_going, l_bound, msg1, spawn_depth, auts)
        end
        t2 = Threads.@spawn begin
            msg2 = copy(current_msg); msg2[depth+1] = 2
            th2, tl2 = copy(curr_tail_H), copy(curr_tail_L)
            _add_gf3_simd!(th2, tl2, col_L, col_H) # Swap H/L for mul by 2
            _Brouwer_Zimmermann_gf3_recursive!(A_packed_H, A_packed_L, r, depth+1, picked+1, th2, tl2, best_w, best_msg, update_lock, lbt, max_canc, w2_table, keep_going, l_bound, msg2, spawn_depth, auts)
        end
        wait(t0); wait(t1); wait(t2)
    else
        # ZERO ALLOCATION SERIAL BLOCK
        current_msg[depth+1] = 0
        _Brouwer_Zimmermann_gf3_recursive!(A_packed_H, A_packed_L, r, depth+1, picked, curr_tail_H, curr_tail_L, best_w, best_msg, update_lock, lbt, max_canc, w2_table, keep_going, l_bound, current_msg, spawn_depth, auts)
        
        current_msg[depth+1] = 1
        _add_gf3_simd!(curr_tail_H, curr_tail_L, col_H, col_L)
        _Brouwer_Zimmermann_gf3_recursive!(A_packed_H, A_packed_L, r, depth+1, picked+1, curr_tail_H, curr_tail_L, best_w, best_msg, update_lock, lbt, max_canc, w2_table, keep_going, l_bound, current_msg, spawn_depth, auts)
        
        current_msg[depth+1] = 2
        _add_gf3_simd!(curr_tail_H, curr_tail_L, col_H, col_L) # Now state is +2
        _Brouwer_Zimmermann_gf3_recursive!(A_packed_H, A_packed_L, r, depth+1, picked+1, curr_tail_H, curr_tail_L, best_w, best_msg, update_lock, lbt, max_canc, w2_table, keep_going, l_bound, current_msg, spawn_depth, auts)
        
        _add_gf3_simd!(curr_tail_H, curr_tail_L, col_H, col_L) # Backtrack to 0
    end
end

function _Brouwer_Zimmermann_gf4_recursive!(
    A_packed_H::Matrix{UInt64}, 
    A_packed_L::Matrix{UInt64}, 
    r::Int, depth::Int, picked::Int, 
    curr_tail_H::Vector{UInt64}, curr_tail_L::Vector{UInt64}, 
    best_w::Threads.Atomic{Int}, best_msg::Vector{Int}, 
    update_lock::Threads.SpinLock, lbt::Vector{Int}, max_canc::Vector{Int}, 
    w2_table::Matrix{Int}, keep_going::Threads.Atomic{Bool}, l_bound::Int, 
    current_msg::Vector{Int}, spawn_depth::Int, auts::Vector{Vector{Int}}
)
    if !keep_going[] return end
    k = size(A_packed_H, 2)
    tw = _fast_simd_wt_gf3(curr_tail_H, curr_tail_L)
    
    if picked == r
        # FIX: Automorphism pruning
        if !_is_canonical(current_msg, auts)
            return
        end

        w = r + tw
        if w < best_w[] 
            lock(update_lock) do
                if w < best_w[]
                    Threads.atomic_xchg!(best_w, w)
                    copyto!(best_msg, current_msg)
                    
                    # FIX: Clean up the un-picked suffix
                    for i in (depth + 1):length(best_msg)
                        best_msg[i] = 0
                    end
                    
                    w <= l_bound && Threads.atomic_cas!(keep_going, true, false)
                end
            end
        end
        return
    end

    rem_to_pick = r - picked
    if depth >= k || (k - depth) < rem_to_pick return end
    min_possible_tw = tw < lbt[rem_to_pick+1] ? lbt[rem_to_pick+1] - tw : (tw > max_canc[rem_to_pick+1] ? tw - max_canc[rem_to_pick+1] : 0)
    if (r + min_possible_tw) >= best_w[] return end

    col_H = view(A_packed_H, :, depth + 1)
    col_L = view(A_packed_L, :, depth + 1)

    if depth < spawn_depth
        tasks = []
        # Val 0 Branch
        t0 = Threads.@spawn begin
            msg0 = copy(current_msg); msg0[depth+1] = 0
            _Brouwer_Zimmermann_gf4_recursive!(A_packed_H, A_packed_L, r, depth+1, picked, copy(curr_tail_H), copy(curr_tail_L), best_w, best_msg, update_lock, lbt, max_canc, w2_table, keep_going, l_bound, msg0, spawn_depth, auts)
        end
        push!(tasks, t0)
        # Non-zero branches
        for v in 1:3
            ti = Threads.@spawn begin
                msg_v = copy(current_msg); msg_v[depth+1] = v
                th_v, tl_v = copy(curr_tail_H), copy(curr_tail_L)
                _add_scaled_gf4_simd!(th_v, tl_v, col_H, col_L, v)
                _Brouwer_Zimmermann_gf4_recursive!(A_packed_H, A_packed_L, r, depth+1, picked+1, th_v, tl_v, best_w, best_msg, update_lock, lbt, max_canc, w2_table, keep_going, l_bound, msg_v, spawn_depth, auts)
            end
            push!(tasks, ti)
        end
        foreach(wait, tasks)
    else
        # SERIAL BLOCK (XOR Backtracking)
        current_msg[depth+1] = 0
        _Brouwer_Zimmermann_gf4_recursive!(A_packed_H, A_packed_L, r, depth+1, picked, curr_tail_H, curr_tail_L, best_w, best_msg, update_lock, lbt, max_canc, w2_table, keep_going, l_bound, current_msg, spawn_depth, auts)
        
        for v in 1:3
            current_msg[depth+1] = v
            _add_scaled_gf4_simd!(curr_tail_H, curr_tail_L, col_H, col_L, v)
            _Brouwer_Zimmermann_gf4_recursive!(A_packed_H, A_packed_L, r, depth+1, picked+1, curr_tail_H, curr_tail_L, best_w, best_msg, update_lock, lbt, max_canc, w2_table, keep_going, l_bound, current_msg, spawn_depth, auts)
            _add_scaled_gf4_simd!(curr_tail_H, curr_tail_L, col_H, col_L, v) # Backtrack
        end
    end
end

function _Brouwer_Zimmermann_nonbinary_recursive!(
    A_raw::Matrix{T}, 
    r::Int, depth::Int, picked::Int, 
    curr_tail::Vector{T}, 
    best_w::Threads.Atomic{Int}, best_msg::Vector{T}, 
    update_lock::Threads.SpinLock, lbt::Vector{Int}, max_canc::Vector{Int}, 
    keep_going::Threads.Atomic{Bool}, l_bound::Int, 
    current_msg::Vector{T}, spawn_depth::Int, auts::Vector{Vector{Int}}, 
    non_zero_elements::Vector{T}
) where T
    if !keep_going[] return end
    k = size(A_raw, 2)
    tw = count(!iszero, curr_tail)
    
    if picked == r
        # FIX: Automorphism pruning
        if !_is_canonical(current_msg, auts)
            return
        end

        w = r + tw
        if w < best_w[] 
            lock(update_lock) do
                if w < best_w[]
                    Threads.atomic_xchg!(best_w, w)
                    copyto!(best_msg, current_msg)
                    
                    # FIX: Clean up the un-picked suffix with generic field zero
                    for i in (depth + 1):length(best_msg)
                        best_msg[i] = zero(T)
                    end
                    
                    w <= l_bound && Threads.atomic_cas!(keep_going, true, false)
                end
            end
        end
        return
    end

    rem_to_pick = r - picked
    if depth >= k || (k - depth) < rem_to_pick return end
    min_possible_tw = tw < lbt[rem_to_pick+1] ? lbt[rem_to_pick+1] - tw : (tw > max_canc[rem_to_pick+1] ? tw - max_canc[rem_to_pick+1] : 0)
    if (r + min_possible_tw) >= best_w[] return end

    if depth < spawn_depth
        tasks = []
        # Zero branch
        t0 = Threads.@spawn begin
            msg0 = copy(current_msg); msg0[depth+1] = zero(parent(A_raw[1]))
            _Brouwer_Zimmermann_nonbinary_recursive!(A_raw, r, depth+1, picked, copy(curr_tail), best_w, best_msg, update_lock, lbt, max_canc, keep_going, l_bound, msg0, spawn_depth, auts, non_zero_elements)
        end
        push!(tasks, t0)
        # Non-zero branches
        for α in non_zero_elements
            ti = Threads.@spawn begin
                msg_α = copy(current_msg); msg_α[depth+1] = α
                tail_α = copy(curr_tail)
                @inbounds for i in eachindex(tail_α)
                    tail_α[i] += α * A_raw[i, depth+1]
                end
                _Brouwer_Zimmermann_nonbinary_recursive!(A_raw, r, depth+1, picked+1, tail_α, best_w, best_msg, update_lock, lbt, max_canc, keep_going, l_bound, msg_α, spawn_depth, auts, non_zero_elements)
            end
            push!(tasks, ti)
        end
        foreach(wait, tasks)
    else
        # SERIAL BLOCK
        current_msg[depth+1] = zero(parent(A_raw[1]))
        _Brouwer_Zimmermann_nonbinary_recursive!(A_raw, r, depth+1, picked, curr_tail, best_w, best_msg, update_lock, lbt, max_canc, keep_going, l_bound, current_msg, spawn_depth, auts, non_zero_elements)

        for α in non_zero_elements
            current_msg[depth+1] = α
            @inbounds for i in eachindex(curr_tail)
                curr_tail[i] += α * A_raw[i, depth+1]
            end
            _Brouwer_Zimmermann_nonbinary_recursive!(A_raw, r, depth+1, picked+1, curr_tail, best_w, best_msg, update_lock, lbt, max_canc, keep_going, l_bound, current_msg, spawn_depth, auts, non_zero_elements)
            @inbounds for i in eachindex(curr_tail)
                curr_tail[i] += (Int(characteristic(parent(A_raw[1]))) - 1) * α * A_raw[i, depth+1] # Backtrack
            end
        end
    end
end

function _minimum_distance_BZ_binary(C::AbstractLinearCode; info_set_alg::Symbol = :auto, verbose::Bool = false)

    !ismissing(C.d) && return C.d
    num_thrds = Threads.nthreads()
    G = generator_matrix(C, true) 
    k, n = size(G)

    if k > 0.75 * n && (2^(n - k) < 1e7)
        verbose && println("High-rate code: using dual weight enumerator.")
        HWE = weight_enumerator(C)
        return minimum(filter(x -> x != 0, [collect(exponent_vectors(HWE.polynomial))[i][1] for i in 1:length(HWE.polynomial)]))
    end

    # 2. Information Set Selection
    if info_set_alg == :auto
        info_set_alg = heuristic_info_set_selection(C)
    end
    
    verbose && println("Using information set algorithm: `$info_set_alg`")
    local z_mats, perms_mats, rnks
    if info_set_alg == :Bouyuklieva
        verbose && println("Using Bouyuklieva (BB21) disjoint partitioning.")
        z_mats_raw, perms_mats, rnks = _partition_disjoint_systematic_sets(G)

        t = count(x -> x == k, rnks)
        a_values = zeros(Int, length(rnks))
        a_values[1] = 1
    else
        verbose && println("Using $info_set_alg overlapping information sets.")
        z_mats_raw, perms_mats, rnks = information_sets(C.G, info_set_alg, permute = true)
        
        # FIX: Filter out zero-rank AND corrupted identity matrices
        valid_idx = Int[]
        for i in 1:length(z_mats_raw)
            curr_rnk = rnks[i]
            if curr_rnk > 0
                is_ident = true
                for r in 1:curr_rnk
                    for c in 1:curr_rnk
                        expected = r == c ? 1 : 0
                        if z_mats_raw[i][r, c] != expected
                            is_ident = false
                            break
                        end
                    end
                    !is_ident && break
                end
                if is_ident
                    push!(valid_idx, i)
                end
            end
        end
        
        z_mats_raw = z_mats_raw[valid_idx]
        perms_mats = perms_mats[valid_idx]
        rnks = rnks[valid_idx]
        
        t = length(rnks) 
        a_values = zeros(Int, length(rnks))
    end
    
    z_mats = [(G = _convert_binary_to_int_matrix(z_mats_raw[i]), perm = perms_mats[i]) for i in 1:length(z_mats_raw)]
    m = length(z_mats)

    current_upper_bound = ismissing(C.u_bound) ? (n + 1) : C.u_bound
    global_min_codeword = zeros(Int, n)
    
    verbose && println("Starting initial row-check (r = 1)...")
    for entry in z_mats
        if info_set_alg == :Bouyuklieva
            for i in 1:size(entry.G, 1)
                row_vec = vec(Array(entry.G[i, :]))
                row_wt = count(!iszero, row_vec)
                if row_wt > 0 && (row_wt < current_upper_bound || (row_wt == current_upper_bound && iszero(global_min_codeword)))
                    current_upper_bound = row_wt
                    reconstructed = zeros(Int, n)
                    for idx in 1:n
                        reconstructed[entry.perm[idx]] = row_vec[idx]
                    end
                    global_min_codeword = reconstructed
                end
            end
        else
            perm_vec = _matrix_to_perm_vector(entry.perm)
            for i in 1:size(entry.G, 1)
                row_vec = vec(Array(entry.G[i, :]))
                row_wt = count(!iszero, row_vec)
                if row_wt > 0 && (row_wt < current_upper_bound || (row_wt == current_upper_bound && iszero(global_min_codeword)))
                    current_upper_bound = row_wt
                    reconstructed = zeros(Int, n)
                    for idx in 1:n
                        reconstructed[perm_vec[idx]] = row_vec[idx]
                    end
                    global_min_codeword = reconstructed
                end
            end
        end
    end
    C.u_bound = current_upper_bound
    verbose && println("Initial row-check set bound to: $(C.u_bound)")

    l_win = 12
    eff_l = min(l_win, n - k)
    
    verbose && println("Launching ISD preprocessor (target < $(C.u_bound))...")
    stagnation_counter = 0
    while stagnation_counter < 3
        target = C.u_bound - 1
        found = Canteaut_Chabaud_attack(C, target; p=2, l=eff_l, max_iters=500)
        if !isempty(found)
            best_wt, local_best = minimum(x -> (count(!iszero, x), x), found)
            if best_wt > 0 && (best_wt < C.u_bound || (best_wt == C.u_bound && iszero(global_min_codeword)))
                C.u_bound = best_wt
                global_min_codeword = Int.(local_best)
                verbose && println("ISD preprocessor dropped bound to: ", C.u_bound)
            end
            stagnation_counter = 0
        else
            stagnation_counter += 1
        end
    end

    processed_configs = []
    auts = _generate_known_automorphisms(C)
    
    for i in 1:m
        entry = z_mats[i]
        current_rnk = rnks[i]
        A_raw = Matrix{Int64}(entry.G[1:current_rnk, current_rnk + 1:end])
        
        # NEW: Construct arbitrary length chunks
        num_chunks = cld(size(A_raw, 2), 64)
        A_packed = Vector{Vector{UInt64}}(undef, current_rnk)
        for r_idx in 1:current_rnk
            chunk_vec = zeros(UInt64, num_chunks)
            for j in 1:size(A_raw, 2)
                if A_raw[r_idx, j] == 1
                    chunk_idx = (j - 1) ÷ 64 + 1
                    bit_idx = (j - 1) % 64
                    chunk_vec[chunk_idx] |= (UInt64(1) << bit_idx)
                end
            end
            A_packed[r_idx] = chunk_vec
        end
        
        lbt, max_canc = _precompute_pruning_bounds_binary(A_raw, current_rnk)
        w2_table = _precompute_weight2_table(A_packed)
        
        if isempty(auts)
            internal_auts = Vector{Vector{Int}}()
        else
            if typeof(entry.perm) <: Vector
                col_perm = entry.perm[1:current_rnk]
            else
                perm_vec = _matrix_to_perm_vector(entry.perm)
                col_perm = perm_vec[1:current_rnk]
            end
            internal_auts = _map_automorphisms(auts, col_perm)
        end
        
        push!(processed_configs, (A = A_packed, lbt = lbt, max_canc = max_canc, 
                                  w2 = w2_table, rnk = current_rnk, auts = internal_auts))
    end

    keep_going = Threads.Atomic{Bool}(true)
    if info_set_alg == :Bouyuklieva
        while true
            active_reduced = count(idx -> idx > t && a_values[idx] > 0, 1:m)
            C.l_bound = sum(a_values) + m - 1
            if !keep_going[] || C.l_bound >= C.u_bound break end
            
            verbose && println("BB21 State $a_values | Bounds: [$(C.l_bound), $(C.u_bound)]")

            for j in 1:m
                a_j = a_values[j]
                if a_j == 0 continue end
                
                config = processed_configs[j]
                search_target = iszero(global_min_codeword) ? (C.u_bound + 1) : C.u_bound
                best_w = Threads.Atomic{Int}(search_target)
                best_msg = zeros(Int, config.rnk)
                update_lock = Threads.SpinLock()

                _Brouwer_Zimmermann_binary_recursive!(
                    config.A, a_j, 0, 0, zeros(UInt64, length(config.A[1])), 
                    best_w, best_msg, update_lock, config.lbt, config.max_canc, config.w2, 
                    keep_going, C.l_bound, zeros(Int, config.rnk), 5, config.auts
                )
                 
                if best_w[] < C.u_bound
                    C.u_bound = best_w[]
                    
                    diff = size(z_mats[j].G, 1) - length(best_msg)
                    full_msg = diff > 0 ? vcat(best_msg, zeros(Int, diff)) : best_msg

                    full_local_c = vec(Array((full_msg' * z_mats[j].G) .% 2))
                    reconstructed = zeros(Int, n)
                   
                    if typeof(perms_mats[j]) <: Vector{Int64}
                        for idx in 1:n
                            reconstructed[perms_mats[j][idx]] = full_local_c[idx]
                        end
                    else
                        perm_vec = _matrix_to_perm_vector(perms_mats[j])
                        for idx in 1:n
                            reconstructed[perm_vec[idx]] = full_local_c[idx]
                        end
                    end
                    
                    global_min_codeword = reconstructed
                    verbose && println("New minimum/witness found in BB21: $(C.u_bound)")
                end
            end
            a_values = _greedy_increment_bb21(a_values, rnks, 2)
        end
    else
        lower_bounds = [_information_set_lower_bound(r, n, k, 0, [0], :auto) for r in 1:k]
        for r in 2:k
            C.l_bound = max(C.l_bound, lower_bounds[r])
            if !keep_going[] || C.l_bound >= C.u_bound break end
            
            num_combinations = binomial(k, r)
            p = Progress(num_combinations * m; dt=0.5, desc="Weight $r search: ", color=:cyan)

            for j in 1:m
                config = processed_configs[j]
                search_target = iszero(global_min_codeword) ? (C.u_bound + 1) : C.u_bound
                best_w = Threads.Atomic{Int}(search_target)
                best_msg = zeros(Int, config.rnk)
                update_lock = Threads.SpinLock()

                _Brouwer_Zimmermann_binary_recursive!(
                    config.A, r, 0, 0, zeros(UInt64, length(config.A[1])), 
                    best_w, best_msg, update_lock, config.lbt, config.max_canc, config.w2, 
                    keep_going, C.l_bound, zeros(Int, config.rnk), 5, config.auts
                )
                
                next!(p, step=num_combinations)

                if best_w[] < C.u_bound
                    C.u_bound = best_w[]
                    
                    # FIX: Unpack using exact matrix multiplication (Identical to BB21)
                    diff = size(z_mats[j].G, 1) - length(best_msg)
                    full_msg = diff > 0 ? vcat(best_msg, zeros(Int, diff)) : best_msg

                    full_local_c = vec(Array((full_msg' * z_mats[j].G) .% 2))
                    reconstructed = zeros(Int, n)
                   
                    if typeof(perms_mats[j]) <: Vector{Int64}
                        for idx in 1:n
                            reconstructed[perms_mats[j][idx]] = full_local_c[idx]
                        end
                    else
                        perm_vec = _matrix_to_perm_vector(perms_mats[j])
                        for idx in 1:n
                            reconstructed[perm_vec[idx]] = full_local_c[idx]
                        end
                    end
                    
                    global_min_codeword = reconstructed
                    verbose && println("New minimum found at weight $r: $(C.u_bound)")
                end
            end
            finish!(p)
        end
    end

    C.d = C.u_bound
    s_glcw = sum(global_min_codeword)
    if s_glcw == C.d
        found_witness = true
    else
        found_witness = false
        s_glcw = 0
    end
    found_witness = s_glcw > 0 

    if found_witness
        if !ismissing(C.P_stand)
            y_std = matrix(C.F, 1, n, global_min_codeword)
            y_orig = y_std * C.P_stand
            global_min_codeword = vec(Array(y_orig))
            verbose && println("Applied P_stand to map witness back to original codespace.")
        end

        y_test = matrix(C.F, 1, n, global_min_codeword)
        if !iszero(parity_check_matrix(C) * transpose(y_test))
            verbose && println("Warning: Saved witness failed parity check! (Likely a corrupted permutation matrix from the BB21 partitioner).")
            verbose && println("Discarding corrupted witness...")
            found_witness = false
            s_glcw = 0
        end
    end

    if s_glcw == 0
        verbose && println("Search confirmed d = $(C.d) but no vector was saved. Launching targeted ISD attack to recover a witness...")
        
        found = Canteaut_Chabaud_attack(C, C.d; p=2, l=eff_l, max_iters=5000)
        if !isempty(found)
            global_min_codeword = [c == 1 ? one(C.F) : zero(C.F) for c in only(found)]
            println("Weight of y: ", sum(Int.(global_min_codeword)))
            found_witness = true
             verbose && println("Targeted ISD attack successfully recovered a minimum-weight codeword!")
        end
        
        if !found_witness
            y = zero_matrix(C.F, 1, n)
            if verbose 
                println("Warning: Targeted search failed to find the codeword in the allotted iterations. Returning a zero vector.")
            end
        else
            y = matrix(C.F, 1, n, global_min_codeword)
            @assert iszero(parity_check_matrix(C) * transpose(y)) "Verification failed: computed codeword is not in the codespace."
        end
    else
        y = matrix(C.F, 1, n, global_min_codeword)
        @assert iszero(parity_check_matrix(C) * transpose(y)) "Verification failed: computed codeword is not in the codespace."
    end
    
    return C.d, y
end
"""
    _reconstruct_codeword(msg_bits::Vector{Int}, packed_tail::Vector{UInt64}, 
                         perm::Vector{Int}, n::Int, k::Int)

Maps a systematic codeword back to original coordinates.
- msg_bits: The weight-r message (first k coordinates of G_sys)
- packed_tail: The XORed parity bits (the n-k coordinates)
- perm: The specific permutation [info_indices; parity_indices] for this set
"""
function _reconstruct_codeword(msg_bits::Vector{Int}, packed_tail::Vector{UInt64}, 
                              perm::Vector{Int}, n::Int, k::Int)
    # 1. Build the codeword in the LOCAL (systematic) order [m | p]
    local_c = zeros(Int, n)
    
    # Fill Information bits
    for i in 1:k
        local_c[i] = msg_bits[i]
    end
    
    # Fill Parity bits by unpacking the UInt64 tail
    # We only care about n-k bits
    tail_len = n - k
    for j in 1:tail_len
        # Bits are packed into UInt64 blocks of 64
        block_idx = div(j - 1, 64) + 1
        bit_pos = (j - 1) % 64
        
        # Extract the bit at bit_pos from the UInt64 block
        if (packed_tail[block_idx] >> bit_pos) & 1 == 1
            local_c[k + j] = 1
        end
    end
    
    # 2. Map back to ORIGINAL coordinates
    # local_c[i] belongs at original index perm[i]
    original_c = zeros(Int, n)
    for i in 1:n
        original_c[perm[i]] = local_c[i]
    end
    
    return original_c
end

function _minimum_distance_BZ_nonbinary(C::AbstractLinearCode; verbose::Bool = false)
    !ismissing(C.d) && return C.d
    
    G_stand = generator_matrix(C, true)
    k, n = size(G_stand)
    q = Int(order(C.F))
    
    # 1. INITIALIZE BOUNDS
    C.u_bound = ismissing(C.u_bound) ? (n + 1) : C.u_bound
    global_min_codeword = [zero(C.F) for _ in 1:n]

    p_char = Int(characteristic(C.F))
    d_deg = degree(C.F)

    current_upper_bound = ismissing(C.u_bound) ? (n + 1) : C.u_bound
    verbose && println("Initial upper bound: $current_upper_bound")
    
    stagnation_counter = 0
    best_cw = nothing 
    l_win = (p_char == 2 && d_deg == 1) ? 12 : 3
    
    while stagnation_counter < 3
        target = current_upper_bound - 1
        found = Stern_attack(C, target; p=2, l=l_win, max_iters=500)
        
        if !isempty(found)
            current_upper_bound, local_best = minimum(x -> (count(!iszero, x), x), found)
            best_cw = local_best
            verbose && println("Preprocessor dropped bound to: ", current_upper_bound)
            stagnation_counter = 0 
        else
            stagnation_counter += 1 
        end
    end
    
    C.u_bound = current_upper_bound
    
    if best_cw !== nothing
        global_min_codeword = [_unpack_field_elem(c, C.F, p_char, d_deg) for c in best_cw]
        verbose && println("ISD preprocessor lowered upper bound to: $(C.u_bound)")
    else
        verbose && println("Preprocessor could not lower the initial bound.")
    end

    # 3. SETUP MATRICES & SMART PACKING
    info_set_alg = heuristic_info_set_selection(C)
    local z_mats, perms_list, rnks
    if info_set_alg == :Bouyuklieva
        z_mats, perms_list, rnks = _partition_disjoint_systematic_sets(G_stand)
    else
        z_mats_raw = _generate_scored_zimmermann_mats_nonbinary(G_stand, 3; pool_size = 20)
        z_mats = [entry.G for entry in z_mats_raw]
        perms_list = [entry.perm for entry in z_mats_raw]
        rnks = [k for _ in 1:length(z_mats)] 
    end

    auts = _generate_known_automorphisms(C)
    m = length(z_mats)

    add_t, mul_t, elem_to_idx = _generate_field_tables_safe(C.F, q)
    
    processed_configs = []
    for i in 1:m
        entry_G = z_mats[i]
        
        # FIX: Ensure permutation is safely mapped to a vector to avoid Oscar type issues
        if typeof(perms_list[i]) <: Vector
            entry_perm = perms_list[i]
        else
            entry_perm = _matrix_to_perm_vector(perms_list[i])
        end
        
        A_raw = entry_G[:, k + 1:end]'
        A_idx = zeros(Int, size(A_raw))
        for idx in eachindex(A_raw)
            A_idx[idx] = elem_to_idx[A_raw[idx]]
        end
        
        lbt, max_canc = _precompute_pruning_bounds_nonbinary(A_idx, k, q)
        
        col_weights = [count(!iszero, view(A_raw, :, j)) for j in 1:k]
        min_parity_wt, j = findmin(col_weights)
        min_row_wt = 1 + min_parity_wt
        
        if min_row_wt <= C.u_bound
            C.u_bound = min_row_wt
            msg = [zero(C.F) for _ in 1:k]
            msg[j] = one(C.F)
            msg_mat = matrix(C.F, 1, k, msg)
            temp_word = vec(Array(msg_mat * entry_G))
            global_min_codeword = temp_word[invperm(entry_perm)]
        end
        
        # FIX: Empty check bypasses overhead for random codes
        if isempty(auts)
            internal_auts = Vector{Vector{Int}}()
        else
            internal_auts = _map_automorphisms(auts, entry_perm[1:k])
        end
        
        if q == 3
            H, L = _pack_matrix_gf3_bitsliced(A_idx)
            w2_table = _precompute_weight2_table_gf3(H, L)
            push!(processed_configs, (A_raw=A_raw, H=H, L=L, lbt=lbt, max_canc=max_canc, auts=internal_auts, w2=w2_table))
        elseif q == 4
            H, L = _pack_matrix_gf4_bitsliced(A_idx)
            w2_table = _precompute_weight2_table_gf4(H, L)
            push!(processed_configs, (A_raw=A_raw, H=H, L=L, lbt=lbt, max_canc=max_canc, auts=internal_auts, w2=w2_table))
        else
            w2_table = _precompute_weight2_table_nonbinary(A_idx, add_t, mul_t, q)
            push!(processed_configs, (A_raw=A_raw, lbt=lbt, max_canc=max_canc, auts=internal_auts, w2=w2_table))
        end
    end

    C.l_bound = 1
    keep_going = Threads.Atomic{Bool}(true)

    # 4. RECURSIVE SEARCH
    if info_set_alg == :Bouyuklieva
        t = count(x -> x == k, rnks) 
        a_values = zeros(Int, m); a_values[1] = 1 
        
        while true
            active_reduced = count(i -> i > t && a_values[i] > 0, 1:m)
            C.l_bound = _information_set_lower_bound(active_reduced, n, k, 0, a_values, :Bouyuklieva; even=is_even(C))
            
            if !keep_going[] || C.l_bound >= C.u_bound
                break
            end
            verbose && println("BB21 State $a_values | Bounds: [$(C.l_bound), $(C.u_bound)]")

            for j in 1:m
                a_j = a_values[j]
                if a_j == 0
                    continue
                end

                config = processed_configs[j]
                best_w = Threads.Atomic{Int}(C.u_bound)
                update_lock = Threads.SpinLock()
                
                best_msg_idx = zeros(Int, k)
                best_msg_field = [zero(C.F) for _ in 1:k]

                if q == 3
                    _Brouwer_Zimmermann_gf3_recursive!(
                        config.H, config.L, a_j, 0, 0, zeros(UInt64, size(config.H, 1)), zeros(UInt64, size(config.L, 1)), 
                        best_w, best_msg_idx, update_lock, config.lbt, config.max_canc, config.w2,
                        keep_going, C.l_bound, zeros(Int, k), 3, config.auts
                    )
                elseif q == 4
                    _Brouwer_Zimmermann_gf4_recursive!(
                        config.H, config.L, a_j, 0, 0, zeros(UInt64, size(config.H, 1)), zeros(UInt64, size(config.L, 1)), 
                        best_w, best_msg_idx, update_lock, config.lbt, config.max_canc, config.w2,
                        keep_going, C.l_bound, zeros(Int, k), 3, config.auts
                    )
                else
                    non_zero_elements = collect(C.F)[2:end]
                    _Brouwer_Zimmermann_nonbinary_recursive!(
                        config.A_raw, a_j, 0, 0, [zero(C.F) for _ in 1:size(config.A_raw, 1)], 
                        best_w, best_msg_field, update_lock, config.lbt, config.max_canc, 
                        keep_going, C.l_bound, [zero(C.F) for _ in 1:k], 3, config.auts, non_zero_elements
                    )
                end

                if best_w[] < C.u_bound
                    lock(update_lock) do
                        C.u_bound = best_w[]
                        
                        if q <= 4
                            actual_elements = collect(C.F)
                            z_idx = findfirst(iszero, actual_elements)
                            if z_idx != 1
                                actual_elements[1], actual_elements[z_idx] = actual_elements[z_idx], actual_elements[1]
                            end
                            
                            for i in 1:k
                                idx = best_msg_idx[i]
                                best_msg_field[i] = (idx == 0) ? actual_elements[1] : actual_elements[idx + 1]
                            end
                        end

                        msg_mat = matrix(C.F, 1, k, best_msg_field)
                        temp_word = vec(Array(msg_mat * z_mats[j]))
                        
                        # FIX: Correct scoping variable call
                        global_min_codeword = temp_word[invperm(perms_list[j])]
                    end
                end
            end
            a_values = _greedy_increment_bb21(a_values, rnks, q)
        end
    else
        lower_bounds = [_information_set_lower_bound(r, n, k, 0, [0], :auto) for r in 1:k]
        for r in 2:k
            C.l_bound = max(C.l_bound, lower_bounds[r])
            if !keep_going[] || C.l_bound >= C.u_bound break end
            
            verbose && println("Weight $r starting. Bounds: [$(C.l_bound), $(C.u_bound)]")

            for j in 1:m
                config = processed_configs[j]
                best_w = Threads.Atomic{Int}(C.u_bound)
                
                best_msg_idx = zeros(Int, k)
                best_msg_field = [zero(C.F) for _ in 1:k]
                
                update_lock = Threads.SpinLock()

                if q == 3
                    _Brouwer_Zimmermann_gf3_recursive!(
                        config.H, config.L, r, 0, 0, zeros(UInt64, size(config.H, 1)), zeros(UInt64, size(config.L, 1)), 
                        best_w, best_msg_idx, update_lock, config.lbt, config.max_canc, config.w2,
                        keep_going, C.l_bound, zeros(Int, k), 3, config.auts
                    )
                elseif q == 4
                    _Brouwer_Zimmermann_gf4_recursive!(
                        config.H, config.L, r, 0, 0, zeros(UInt64, size(config.H, 1)), zeros(UInt64, size(config.L, 1)), 
                        best_w, best_msg_idx, update_lock, config.lbt, config.max_canc, config.w2,
                        keep_going, C.l_bound, zeros(Int, k), 3, config.auts
                    )
                else
                    non_zero_elements = collect(C.F)[2:end] 
                    _Brouwer_Zimmermann_nonbinary_recursive!(
                        config.A_raw, r, 0, 0, [zero(C.F) for _ in 1:size(config.A_raw, 1)], 
                        best_w, best_msg_field, update_lock, config.lbt, config.max_canc, 
                        keep_going, C.l_bound, [zero(C.F) for _ in 1:k], 3, config.auts, non_zero_elements
                    )
                end
                
                if best_w[] < C.u_bound
                    C.u_bound = best_w[]
                    
                    if q <= 4
                        actual_elements = collect(C.F)
                        z_idx = findfirst(iszero, actual_elements)
                        if z_idx != 1
                            actual_elements[1], actual_elements[z_idx] = actual_elements[z_idx], actual_elements[1]
                        end
                        
                        for i in 1:k
                            idx = best_msg_idx[i]
                            best_msg_field[i] = (idx == 0) ? actual_elements[1] : actual_elements[idx + 1]
                        end
                    end

                    msg_mat = matrix(C.F, 1, k, best_msg_field)
                    temp_word = vec(Array(msg_mat * z_mats[j]))
                    
                    # FIX: Correct scoping variable call
                    global_min_codeword = temp_word[invperm(perms_list[j])]
                end
            end
            C.l_bound = max(C.l_bound, m * (r + 1))
        end
    end

    # 5. ENDGAME / TARGETED ISD
    C.d = C.u_bound
    s_glcw = count(!iszero, global_min_codeword)
    
    if s_glcw == 0
        verbose && println("Search confirmed d = $(C.d) but no vector was saved. Launching targeted Stern's attack...")
        
        found_witness = false 
        
        target_vecs = Stern_attack(C, C.d; p=p_char, l=l_win, num_find=1, max_iters=1000)
        if !isempty(target_vecs)
            raw_vec = only(target_vecs)
            global_min_codeword = [_unpack_field_elem(c, C.F, p_char, d_deg) for c in raw_vec]
            found_witness = true
            verbose && println("Targeted search successfully recovered a minimum-weight codeword!")
        end
        
        if !found_witness
            y = zero_matrix(C.F, 1, n)
            verbose && println("Warning: Targeted search failed. Returning a zero vector.")
        else
            y = matrix(C.F, 1, n, global_min_codeword) * C.P_stand
        end
    else
        y = matrix(C.F, 1, n, global_min_codeword) * C.P_stand
    end
    
    return C.d, y
end

"""
    _minimum_distance_wagner_mitm_binary(C::AbstractLinearCode; max_d::Int=20, verbose::Bool=false)

Optimized exact minimum distance solver using Wagner's Meet-in-the-Middle 
strictly for binary (GF(2)) codes. Eliminates scalar coefficient loops and 
uses modulo 2 arithmetic for direct hash collisions.
"""
function _minimum_distance_wagner_mitm_binary(C::AbstractLinearCode; max_d::Int=20, verbose::Bool=false)
    # FIX: Use the safe binary-to-int converter instead of Int.()
    H = _convert_binary_to_int_matrix(parity_check_matrix(C))
    r, n = size(H)
    
    mid = div(n, 2)
    H_L = view(H, :, 1:mid)
    H_R = view(H, :, (mid+1):n)
    
    function build_syndrome_table_bin(H_half, target_wt, offset)
        _, cols = size(H_half)
        table = Dict{Vector{Int}, Vector{Int}}()
        
        # FIX: Explicitly qualify Combinatorics
        for col_indices in Combinatorics.combinations(1:cols, target_wt)
            syn = zeros(Int, r)
            for c in col_indices
                for row in 1:r
                    syn[row] = (syn[row] + H_half[row, c]) % 2
                end
            end
            
            if !haskey(table, syn)
                table[syn] = [c + offset for c in col_indices]
            end
        end
        return table
    end

    verbose && println("Starting Optimized Binary Wagner MitM search...")

    for w in 1:max_d
        verbose && println("  Checking for codewords of total weight $w...")
        
        for w_L in 0:w
            w_R = w - w_L
            
            if w_L > mid || w_R > (n - mid)
                continue
            end
            
            left_table = build_syndrome_table_bin(H_L, w_L, 0)
            
            # FIX: Explicitly qualify Combinatorics
            for col_indices in Combinatorics.combinations(1:(n - mid), w_R)
                syn_R = zeros(Int, r)
                for c in col_indices
                    for row in 1:r
                        syn_R[row] = (syn_R[row] + H_R[row, c]) % 2
                    end
                end
                
                if haskey(left_table, syn_R)
                    verbose && println("Binary Collision found! Left wt: $w_L, Right wt: $w_R")
                    
                    witness = zero_matrix(C.F, 1, n)
                    
                    L_cols = left_table[syn_R]
                    for idx in L_cols
                        witness[1, idx] = one(C.F)
                    end
                    
                    for c in col_indices
                        witness[1, c + mid] = one(C.F)
                    end
                    
                    @assert iszero(parity_check_matrix(C) * transpose(witness)) "Wagner reconstructed a failed witness!"
                    
                    return w, witness
                end
            end
        end
    end
    
    verbose && println("No codewords found up to weight $max_d.")
    return -1, zero_matrix(C.F, 1, n)
end

"""
    _minimum_distance_wagner_mitm_nonbinary(C::AbstractLinearCode; max_d::Int=20, verbose::Bool=false)

Computes the exact minimum distance of a code using the Syndrome Meet-in-the-Middle 
(Wagner's) algorithm. Extremely fast for small d, but RAM usage explodes as d grows.
"""
function _minimum_distance_wagner_mitm_nonbinary(C::AbstractLinearCode; max_d::Int=20, verbose::Bool=false)
    H = Array(parity_check_matrix(C))
    r, n = size(H)
    F = parent(H[1,1])
    q = Int(order(F))
    
    nonzero_elements = filter(!iszero, collect(F))
    
    mid = div(n, 2)
    H_L = view(H, :, 1:mid)
    H_R = view(H, :, (mid+1):n)
    
    function build_syndrome_table(H_half, target_wt, offset)
        _, cols = size(H_half)
        table = Dict{Vector{typeof(zero(F))}, Tuple{Vector{Int}, Vector{typeof(zero(F))}}}()
        
        # FIX: Explicitly qualify Combinatorics
        for col_indices in Combinatorics.combinations(1:cols, target_wt)
            for scalars in Iterators.product(fill(nonzero_elements, target_wt)...)
                
                syn = fill(zero(F), r)
                for (i, c) in enumerate(col_indices)
                    for row in 1:r
                        syn[row] += scalars[i] * H_half[row, c]
                    end
                end
                
                if !haskey(table, syn)
                    shifted_cols = [c + offset for c in col_indices]
                    table[syn] = (shifted_cols, collect(scalars))
                end
            end
        end
        return table
    end

    verbose && println("Starting Wagner Syndrome MitM search...")

    for w in 1:max_d
        verbose && println("  Checking for codewords of total weight $w...")
        
        for w_L in 0:w
            w_R = w - w_L
            
            if w_L > mid || w_R > (n - mid)
                continue
            end
            
            left_table = build_syndrome_table(H_L, w_L, 0)
            
            # FIX: Explicitly qualify Combinatorics
            for col_indices in Combinatorics.combinations(1:(n - mid), w_R)
                for scalars in Iterators.product(fill(nonzero_elements, w_R)...)
                    
                    syn_R = fill(zero(F), r)
                    for (i, c) in enumerate(col_indices)
                        for row in 1:r
                            syn_R[row] += scalars[i] * H_R[row, c]
                        end
                    end
                    
                    target_syn = [-x for x in syn_R]
                    
                    if haskey(left_table, target_syn)
                        verbose && println("Collision found! Left wt: $w_L, Right wt: $w_R")
                        
                        witness = zero_matrix(F, 1, n)
                        
                        L_cols, L_scalars = left_table[target_syn]
                        for idx in 1:length(L_cols)
                            witness[1, L_cols[idx]] = L_scalars[idx]
                        end
                        
                        shifted_R_cols = [c + mid for c in col_indices]
                        for idx in 1:length(shifted_R_cols)
                            witness[1, shifted_R_cols[idx]] = scalars[idx]
                        end
                        
                        @assert iszero(parity_check_matrix(C) * transpose(witness)) "Wagner reconstructed a failed witness!"
                        
                        return w, witness
                    end
                end
            end
        end
    end
    
    verbose && println("No codewords found up to weight $max_d.")
    return -1, zero_matrix(F, 1, n)
end

"""
    _minimum_distance_ILP(C::AbstractLinearCode; verbose::Bool = false)

Return the minimum distance of the linear code using an integer linear programming approach.

# Note
- Run `using JuMP, GLPK` to activate this extension.
"""
function _minimum_distance_ILP end

"""
    minimum_distance(C::AbstractLinearCode; alg::Symbol = :trellis, sect::Bool = false, verbose::Bool = false)

Return the minimum distance of the linear code if known, otherwise computes it
using the algorithm of `alg`. If `alg = "trellis"`, the sectionalization flag
`sect` can be set to true to further compactify the reprsentation.
"""
function minimum_distance(C::AbstractLinearCode; alg::Symbol = :auto,
    info_set_alg::Symbol = :auto, auts::Vector{Vector{Int}} = [Int[]], verbose::Bool = false)

    !ismissing(C.d) && return C.d

    alg ∈ (:auto, :BZ, :trellis, :bruteforce, :wt_dist, :Leon, :Wagner, :ILP) ||
        throw(ArgumentError("Unexpected algorithm '$alg'."))
    info_set_alg ∈ (:auto, :Brouwer, :Zimmermann, :White, :Chen, :Bouyuklieva, :Edmonds) || throw(ArgumentError("Unknown information set algorithm. Expected `:auto`, `:Brouwer`, `:Zimmermann`, `:White`, `:Chen`, `:Bouyuklieva`, or `:Edmonds`."))
    
    k, n = C.k, C.n
    q = Int(order(C.F))

    if alg == :auto
        card_C = BigInt(q)^k
        card_D = BigInt(q)^(n - k)

        # 1. TRIVIAL FAST PATH: Primal Brute Force
        if card_C <= 1e6 # random cutoff
            verbose && println("Auto: Small cardinality ($card_C). Routing to Primal Brute Force.")
            C.weight_enum = _weight_enumerator_BF(C.G)
            HWE = CWE_to_HWE(C.weight_enum)
            C.d = minimum(filter(x -> x != 0, [collect(exponent_vectors(HWE.polynomial))[i][1]
                for i in 1:length(HWE.polynomial)]))
            return C.d
        end

        # 2. TRIVIAL FAST PATH: Dual Brute Force
        if card_D <= 1e6 # random cutoff
            verbose && println("Auto: Small dual cardinality ($card_D). Routing to Dual Brute Force.")
            D = dual(C)
            D.weight_enum = _weight_enumerator_BF(D.G)
            C.weight_enum = MacWilliams_identity(D, D.weight_enum)
            HWE = CWE_to_HWE(C.weight_enum)
            C.d = minimum(filter(x -> x != 0, [collect(exponent_vectors(HWE.polynomial))[i][1]
                for i in 1:length(HWE.polynomial)]))
            return C.d
        end

        # 3. SPARSITY TRAP (LDPC / Sparse Parity Matrices)
        H = Array(parity_check_matrix(C))
        density = count(!iszero, H) / length(H)
        if density < 0.10 && n <= 150
            verbose && println("Auto: High sparsity ($(round(density*100, digits=1))%). Routing to ILP.")
            d_ilp = _minimum_distance_ILP(C; verbose = verbose)
            if d_ilp > 0
                C.d = d_ilp
                return C.d
            end
            verbose && println("ILP stalled or failed. Falling back...")
        end

        # 4. BINARY PRE-FLIGHT CHECK (Wagner Meet-in-the-Middle)
        if q == 2 && n <= 128
            verbose && println("Auto: Binary code detected. Running Wagner MitM pre-flight (d <= 5)...")
            d_wagner, witness_wagner = _minimum_distance_wagner_mitm_binary(C; max_d = 5, verbose = false)
            if d_wagner != -1
                verbose && println("Auto: Wagner MitM caught early collision!")
                C.d = d_wagner
                return C.d, witness_wagner
            end
            verbose && println("Auto: No low-weight words found. Proceeding to deep search...")
        end

        # 5. DUAL TRELLIS vs HYBRID BZ
        if rate(C) > 0.5
            verbose && println("Auto: High rate (> 0.5). Routing to Dual Syndrome Trellis.")
            D = dual(C)
            # Assuming weight_enumerator_classical stores the result in D.weight_enum
            weight_enumerator_classical(syndrome_trellis(D, "primal", false), type = :CWE)
            C.weight_enum = MacWilliams_identity(D, D.weight_enum)
            HWE = CWE_to_HWE(C.weight_enum)
            C.d = minimum(filter(x -> x != 0, [collect(exponent_vectors(HWE.polynomial))[i][1]
                for i in 1:length(HWE.polynomial)]))
            return C.d
        else
            verbose && println("Auto: Moderate/Low rate. Routing to Hybrid BZ-Trellis Bridge.")
            # We replace your Gray code fallback with the state-of-the-art Hybrid BZ
            C.d = minimum_distance_hybrid(C; max_span = 15, verbose = verbose)
            return C.d
        end

    # --- EXPLICIT ROUTING BYPASSES AUTO ---
    elseif alg == :BZ
        if q == 2
            return _minimum_distance_BZ_binary(C; info_set_alg = info_set_alg, verbose = verbose)
        else
            return _minimum_distance_BZ_nonbinary(C; verbose = verbose)
        end
    elseif alg == :trellis
        weight_enumerator_classical(syndrome_trellis(C, "primal", false), type = :CWE)
        return C.d
    elseif alg == :bruteforce
        C.weight_enum = _weight_enumerator_BF(C.G)
        HWE = CWE_to_HWE(C.weight_enum)
        C.d = minimum(filter(x -> x != 0, [collect(exponent_vectors(HWE.polynomial))[i][1]
            for i in 1:length(HWE.polynomial)]))
        return C.d
    elseif alg == :wt_dist
        HWE = weight_enumerator(C, type = :Hamming, alg = alg)
        !ismissing(C.d) && return C.d
        C.d = minimum(filter(x -> x != 0, [collect(exponent_vectors(HWE.polynomial))[i][1]
            for i in 1:length(HWE.polynomial)]))
        return C.d
    elseif alg == :Wagner
        if q == 2
            return _minimum_distance_wagner_mitm_binary(C; verbose = verbose)
        else
            return _minimum_distance_wagner_mitm_nonbinary(C; verbose = verbose)
        end
    elseif alg == :ILP
        return _minimum_distance_ILP(C; verbose = verbose)
    # elseif alg == :Leon
    #     Leon(C)
    end
end
