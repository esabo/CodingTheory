# ==============================================================================
# ALGORITHMS: ORDER STATISTIC DECODING (OSD)
# ==============================================================================

"""
A zero-allocation workspace for Order Statistic Decoding (OSD).
Gaussian Elimination is performed in-place using a pre-allocated dense matrix.
"""
struct OSDWorkspace
    num_var::Int
    num_check::Int
    
    # Gaussian Elimination Matrix
    H_dense::Matrix{UInt8}
    H_work::Matrix{UInt8}
    
    # Sorting and Permutation Tracking
    reliabilities::Vector{Float64}
    perm::Vector{Int}
    inv_perm::Vector{Int}
    
    # Bit Tracking
    is_mrb::Vector{Bool}      # True if the bit belongs to the Most Reliable Basis
    pivot_map::Vector{Int}    # Maps a parity column to its pivot row
    
    # Codeword Generation
    hard_decisions::Vector{UInt8}
    candidate_cw::Vector{UInt8}
    best_cw::Vector{UInt8}
    
    # Error pattern generation for OSD-m
    test_pattern::Vector{UInt8}
end

"""
Initialize the OSD workspace. 
Requires the original parity-check matrix to allocate the dense GE workspace.
"""
function init_osd_workspace(H::AbstractMatrix)
    num_check, num_var = size(H)
    
    H_dense = zeros(UInt8, num_check, num_var)
    for c in 1:num_check
        for v in 1:num_var
            H_dense[c, v] = H[c, v] != 0 ? 0x01 : 0x00
        end
    end
    
    return OSDWorkspace(
        num_var, num_check,
        H_dense, zeros(UInt8, num_check, num_var),
        zeros(Float64, num_var), collect(1:num_var), zeros(Int, num_var),
        zeros(Bool, num_var), zeros(Int, num_var),
        zeros(UInt8, num_var), zeros(UInt8, num_var), zeros(UInt8, num_var),
        zeros(UInt8, num_var)
    )
end

"""
Evaluates a specific OSD bit-flip pattern. 
Re-encodes the LRB using the permuted syndrome and updates W.best_cw if it is the new minimum.
"""
@inline function _evaluate_pattern!(W::OSDWorkspace, flip_1::Int, flip_2::Int, current_min_dist::Float64)
    # Load base MRB decisions
    for i in 1:W.num_var
        W.candidate_cw[i] = W.hard_decisions[W.perm[i]]
    end
    
    # Apply requested flips
    if flip_1 > 0; W.candidate_cw[flip_1] ⊻= 0x01; end
    if flip_2 > 0; W.candidate_cw[flip_2] ⊻= 0x01; end
    
    # Re-encode LRB
    for col in 1:W.num_var
        if !W.is_mrb[col]
            row = W.pivot_map[col]
            parity_val = W.s_work[row] 
            
            @simd for m_col in 1:W.num_var
                if W.is_mrb[m_col]
                    parity_val ⊻= (W.H_work[row, m_col] & W.candidate_cw[m_col])
                end
            end
            W.candidate_cw[col] = parity_val
        end
    end
    
    # Calculate Euclidean distance
    dist = 0.0
    for i in 1:W.num_var
        orig_i = W.perm[i]
        if W.candidate_cw[i] != W.hard_decisions[orig_i]
            dist += W.reliabilities[orig_i]
        end
    end
    
    # Track the global best
    if dist < current_min_dist
        for i in 1:W.num_var
            W.best_cw[W.perm[i]] = W.candidate_cw[i]
        end
        return dist
    end
    
    return current_min_dist
end

# ==============================================================================
# STANDARD OSD SWEEPS
# ==============================================================================

@inline function _run_osd_sweeps!(::Val{:standard}, ::Val{0}, W, mrb_indices, cs_lambda)
    _evaluate_pattern!(W, 0, 0, Inf)
end

@inline function _run_osd_sweeps!(::Val{:standard}, ::Val{1}, W, mrb_indices, cs_lambda)
    min_dist = _evaluate_pattern!(W, 0, 0, Inf)
    for i in mrb_indices
        min_dist = _evaluate_pattern!(W, i, 0, min_dist)
    end
end

@inline function _run_osd_sweeps!(::Val{:standard}, ::Val{2}, W, mrb_indices, cs_lambda)
    min_dist = _evaluate_pattern!(W, 0, 0, Inf)
    for idx1 in 1:length(mrb_indices)
        i = mrb_indices[idx1]
        min_dist = _evaluate_pattern!(W, i, 0, min_dist)
        
        for idx2 in (idx1+1):length(mrb_indices)
            j = mrb_indices[idx2]
            min_dist = _evaluate_pattern!(W, i, j, min_dist)
        end
    end
end

# ==============================================================================
# COMBINATORIAL SWEEP (OSD-CS)
# ==============================================================================
# (Note: OSD-CS with order 0 or 1 is mathematically identical to Standard OSD)

@inline _run_osd_sweeps!(::Val{:cs}, ::Val{0}, W, mrb_indices, cs_lambda) = _run_osd_sweeps!(Val(:standard), Val(0), W, mrb_indices, cs_lambda)
@inline _run_osd_sweeps!(::Val{:cs}, ::Val{1}, W, mrb_indices, cs_lambda) = _run_osd_sweeps!(Val(:standard), Val(1), W, mrb_indices, cs_lambda)

@inline function _run_osd_sweeps!(::Val{:cs}, ::Val{2}, W, mrb_indices, cs_lambda)
    # Phase 1: Weight-1 sweep across the ENTIRE MRB
    min_dist = _evaluate_pattern!(W, 0, 0, Inf)
    for i in mrb_indices
        min_dist = _evaluate_pattern!(W, i, 0, min_dist)
    end
    
    # Phase 2: Weight-2 sweep restricted strictly to the least reliable `cs_lambda` bits
    cs_start = max(1, length(mrb_indices) - cs_lambda + 1)
    for idx1 in cs_start:length(mrb_indices)
        i = mrb_indices[idx1]
        for idx2 in (idx1+1):length(mrb_indices)
            j = mrb_indices[idx2]
            min_dist = _evaluate_pattern!(W, i, j, min_dist)
        end
    end
end

function _fast_osd!(W::OSDWorkspace, method::Val, order::Val, total_llrs::Vector{Float64}, cs_lambda::Int)
    
    # 1. Sort by Reliability & Copy Syndrome
    for v in 1:W.num_var
        W.reliabilities[v] = abs(total_llrs[v])
        W.hard_decisions[v] = total_llrs[v] < 0.0 ? 0x01 : 0x00
    end
    sortperm!(W.perm, W.reliabilities, rev=true)
    for i in 1:W.num_var; W.inv_perm[W.perm[i]] = i; end
    
    for c in 1:W.num_check
        for i in 1:W.num_var
            W.H_work[c, i] = W.H_dense[c, W.perm[i]]
        end
    end
    copyto!(W.s_work, syndrome)
    
    # 2. Quantum Gaussian Elimination (In-place)
    fill!(W.is_mrb, true)
    fill!(W.pivot_map, 0)
    
    pivot_row = 1
    for col in W.num_var:-1:1
        target_row = 0
        for r in pivot_row:W.num_check
            if W.H_work[r, col] == 0x01
                target_row = r; break
            end
        end
        
        if target_row > 0
            if target_row != pivot_row
                for i in 1:W.num_var
                    temp = W.H_work[pivot_row, i]
                    W.H_work[pivot_row, i] = W.H_work[target_row, i]
                    W.H_work[target_row, i] = temp
                end
                temp_s = W.s_work[pivot_row]
                W.s_work[pivot_row] = W.s_work[target_row]
                W.s_work[target_row] = temp_s
            end
            
            for r in 1:W.num_check
                if r != pivot_row && W.H_work[r, col] == 0x01
                    @simd for i in 1:W.num_var
                        W.H_work[r, i] ⊻= W.H_work[pivot_row, i]
                    end
                    W.s_work[r] ⊻= W.s_work[pivot_row]
                end
            end
            
            W.is_mrb[col] = false
            W.pivot_map[col] = pivot_row
            pivot_row += 1
            if pivot_row > W.num_check; break; end
        end
    end
    
    # 3. Extract sorted MRB indices
    mrb_indices = Int[]
    for i in 1:W.num_var
        if W.is_mrb[i]; push!(mrb_indices, i); end
    end
    
    # 4. Trigger the zero-cost Sweep Trait
    _run_osd_sweeps!(method, order, W, mrb_indices, cs_lambda)
    
    return W.best_cw
end

# ==============================================================================
# CLASSICAL OSD API (No syndrome passed)
# ==============================================================================
function osd_decode!(W::OSDWorkspace, total_llrs::Vector{Float64};
                     method::Symbol = :cs, 
                     order::Int = 2, 
                     cs_lambda::Int = 10)
    
    # In-place wipe. Zero allocations.
    fill!(W.s_work, 0x00)
    
    return _fast_osd!(W, Val(method), Val(order), total_llrs, cs_lambda)
end

# ==============================================================================
# QUANTUM / SYNDROME OSD API (Syndrome passed)
# ==============================================================================
function osd_decode!(W::OSDWorkspace, total_llrs::Vector{Float64}, syndrome::Vector{UInt8};
                     method::Symbol = :cs, 
                     order::Int = 2, 
                     cs_lambda::Int = 10)
    
    # In-place copy. Zero allocations.
    copyto!(W.s_work, syndrome)
    
    return _fast_osd!(W, Val(method), Val(order), total_llrs, cs_lambda)
end

# ==============================================================================
# ALGORITHMS: GUESSING RANDOM ADDITIVE NOISE DECODING (GRAND)
# ==============================================================================

struct GRANDWorkspace
    num_var::Int
    num_check::Int
    
    # Graph Mapping (Sparse is faster for incremental syndrome updates)
    var_to_checks::Vector{Vector{Int}}
    
    # Trackers
    reliabilities::Vector{Float64}
    perm::Vector{Int}
    candidate_cw::Vector{UInt8}
    
    # Syndrome Tracking
    target_syndrome::Vector{UInt8}
    current_syndrome::Vector{UInt8}
end

function init_grand_workspace(H::AbstractMatrix)
    num_check, num_var = size(H)
    var_to_checks = [Int[] for _ in 1:num_var]
    
    for c in 1:num_check
        for v in 1:num_var
            if !iszero(H[c, v])
                push!(var_to_checks[v], c)
            end
        end
    end
    
    return GRANDWorkspace(
        num_var, num_check, var_to_checks,
        zeros(Float64, num_var), collect(1:num_var), zeros(UInt8, num_var),
        zeros(UInt8, num_check), zeros(UInt8, num_check)
    )
end

"""
Executes Post-BP GRAND. 
Sweeps all weight-1, weight-2, and weight-3 error patterns across the `max_lrb` Least Reliable Bits.
"""
function grand_decode!(W::GRANDWorkspace, total_llrs::Vector{Float64};
                       syndrome::Vector{UInt8} = zeros(UInt8, W.num_check),
                       max_lrb::Int = 12, 
                       max_weight::Int = 3)
    
    # 1. Load Hard Decisions & Calculate Base Parity
    fill!(W.current_syndrome, 0x00)
    copyto!(W.target_syndrome, syndrome)
    
    for v in 1:W.num_var
        W.reliabilities[v] = abs(total_llrs[v])
        bit_val = total_llrs[v] < 0.0 ? 0x01 : 0x00
        W.candidate_cw[v] = bit_val
        
        if bit_val == 0x01
            for c in W.var_to_checks[v]
                W.current_syndrome[c] ⊻= 0x01
            end
        end
    end
    
    # Base case check: Did BP actually succeed and just fail to flag it?
    if W.current_syndrome == W.target_syndrome
        return true, W.candidate_cw
    end
    
    # 2. Sort to find the Least Reliable Bits (Ascending order)
    sortperm!(W.perm, W.reliabilities)
    search_space = min(max_lrb, W.num_var)
    lrb_indices = W.perm[1:search_space]
    
    # Helper to check syndrome incrementally
    @inline function _check_pattern(flips...)
        # Apply flips incrementally
        for v in flips
            W.candidate_cw[v] ⊻= 0x01
            for c in W.var_to_checks[v]
                W.current_syndrome[c] ⊻= 0x01
            end
        end
        
        success = W.current_syndrome == W.target_syndrome
        
        # Revert flips to leave workspace clean for next pattern
        if !success
            for v in flips
                W.candidate_cw[v] ⊻= 0x01
                for c in W.var_to_checks[v]
                    W.current_syndrome[c] ⊻= 0x01
                end
            end
        end
        return success
    end
    
    # 3. Weight-1 Sweep
    if max_weight >= 1
        for i in 1:search_space
            if _check_pattern(lrb_indices[i])
                return true, W.candidate_cw
            end
        end
    end
    
    # 4. Weight-2 Sweep
    if max_weight >= 2
        for i in 1:search_space
            for j in (i+1):search_space
                if _check_pattern(lrb_indices[i], lrb_indices[j])
                    return true, W.candidate_cw
                end
            end
        end
    end
    
    # 5. Weight-3 Sweep
    if max_weight >= 3
        for i in 1:search_space
            for j in (i+1):search_space
                for k in (j+1):search_space
                    if _check_pattern(lrb_indices[i], lrb_indices[j], lrb_indices[k])
                        return true, W.candidate_cw
                    end
                end
            end
        end
    end
    
    return false, W.candidate_cw
end

# ==============================================================================
# ALGORITHMS: RESIDUAL WEIGHTED BIT-FLIPPING (WBF)
# ==============================================================================

struct WBFWorkspace
    num_var::Int
    num_check::Int
    
    var_to_checks::Vector{Vector{Int}}
    chk_to_vars::Vector{Vector{Int}}
    
    candidate_cw::Vector{UInt8}
    mismatch_syndrome::Vector{UInt8} # 1 if check fails, 0 if it passes
end

function init_wbf_workspace(H::AbstractMatrix)
    num_check, num_var = size(H)
    var_to_checks = [Int[] for _ in 1:num_var]
    chk_to_vars = [Int[] for _ in 1:num_check]
    
    for c in 1:num_check
        for v in 1:num_var
            if !iszero(H[c, v])
                push!(var_to_checks[v], c)
                push!(chk_to_vars[c], v)
            end
        end
    end
    
    return WBFWorkspace(
        num_var, num_check, var_to_checks, chk_to_vars,
        zeros(UInt8, num_var), zeros(UInt8, num_check)
    )
end

"""
Executes Post-BP WBF. 
Iteratively flips the single bit with the highest energy score until the syndrome is zero.
"""
function wbf_decode!(W::WBFWorkspace, total_llrs::Vector{Float64};
                     syndrome::Vector{UInt8} = zeros(UInt8, W.num_check),
                     alpha::Float64 = 0.5,
                     max_iters::Int = 50)
    
    # 1. Initialize hard decisions and calculate the mismatch syndrome
    fill!(W.mismatch_syndrome, 0x00)
    
    for v in 1:W.num_var
        W.candidate_cw[v] = total_llrs[v] < 0.0 ? 0x01 : 0x00
        if W.candidate_cw[v] == 0x01
            for c in W.var_to_checks[v]
                W.mismatch_syndrome[c] ⊻= 0x01
            end
        end
    end
    
    # XOR with target to strictly track failures (1 = fail, 0 = pass)
    for c in 1:W.num_check
        W.mismatch_syndrome[c] ⊻= syndrome[c]
    end
    
    # 2. Iterative bit-flipping loop
    for iter in 1:max_iters
        # Check if all parity equations are satisfied
        is_valid = true
        for c in 1:W.num_check
            if W.mismatch_syndrome[c] == 0x01
                is_valid = false
                break
            end
        end
        
        if is_valid
            return true, W.candidate_cw, iter
        end
        
        # Find the bit with the maximum energy score
        best_v = -1
        max_score = -Inf
        
        for v in 1:W.num_var
            # Only evaluate bits connected to AT LEAST one failed parity check
            failed_checks_count = 0
            for c in W.var_to_checks[v]
                if W.mismatch_syndrome[c] == 0x01
                    failed_checks_count += 1
                end
            end
            
            if failed_checks_count > 0
                # Calculate energy: reward failed checks, penalize LLR confidence
                score = failed_checks_count - (alpha * abs(total_llrs[v]))
                
                if score > max_score
                    max_score = score
                    best_v = v
                end
            end
        end
        
        # If no bits are connected to failed checks (should be impossible if invalid), break
        if best_v == -1
            break
        end
        
        # 3. Flip the winning bit and incrementally update the mismatch syndrome
        W.candidate_cw[best_v] ⊻= 0x01
        for c in W.var_to_checks[best_v]
            W.mismatch_syndrome[c] ⊻= 0x01
        end
    end
    
    return false, W.candidate_cw, max_iters
end
