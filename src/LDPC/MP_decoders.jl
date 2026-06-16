# Copyright (c) 2023 - 2026 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

# In Julia, a mutable struct is allocated on the heap.
# This struct can be stack-allocated or entirely compiled away into CPU registers.
struct DecoderWorkspace{T}
    # 1. Node states
    channel_LLRs::Vector{T}  # Size: N (The raw channel input)
    total_LLRs::Vector{T}    # Size: N (The running posterior beliefs)
    current_bits::Vector{Int}# Size: N (The hard decisions)
    
    # 2. Edge tracking (Flattened 1D arrays)
    # Instead of an M x N matrix, we just store the exact edges.
    C2V::Vector{T}           # Size: E (Check to Variable messages)
    V2C::Vector{T}           # Size: E (Variable to Check messages)
    
    # Pre-calculated routing (where does edge 'i' go?)
    edge_to_var::Vector{Int} # Size: E
    edge_to_chk::Vector{Int} # Size: E
    
    # Optional: Layer boundaries for layered scheduling
    layer_ptrs::Vector{Int} 
end

struct HardDecisionWorkspace
    num_var::Int
    num_check::Int
    num_edges::Int
    
    V2C::Vector{UInt8}
    C2V::Vector{UInt8}
    
    var_to_edges::Vector{Vector{Int}}
    chk_to_edges::Vector{Vector{Int}}
    chk_to_vars::Vector{Vector{Int}}
    
    channel_bits::Vector{UInt8}
    current_bits::Vector{UInt8}
    
    # O(1) Masks
    is_decimated::Vector{Bool}
    is_erased::Vector{Bool}
end

# ==============================================================================
# 2. THE ONE-TIME RAM ALLOCATOR
# ==============================================================================

function init_hard_workspace(H::AbstractMatrix)
    num_check, num_var = size(H)
    
    var_to_edges = [Int[] for _ in 1:num_var]
    chk_to_edges = [Int[] for _ in 1:num_check]
    chk_to_vars  = [Int[] for _ in 1:num_check]
    
    edge_id = 1
    for c in 1:num_check
        for v in 1:num_var
            if !iszero(H[c, v])
                push!(chk_to_edges[c], edge_id)
                push!(var_to_edges[v], edge_id)
                push!(chk_to_vars[c], v)
                edge_id += 1
            end
        end
    end
    
    num_edges = edge_id - 1
    
    return HardDecisionWorkspace(
        num_var, num_check, num_edges,
        zeros(UInt8, num_edges), zeros(UInt8, num_edges),
        var_to_edges, chk_to_edges, chk_to_vars,
        zeros(UInt8, num_var), zeros(UInt8, num_var),
        zeros(Bool, num_var), zeros(Bool, num_var)
    )
end

@inline function load_hard_channel!(W::HardDecisionWorkspace, v_in::Vector{UInt8};
                                    erasures::Vector{Int} = Int[],
                                    decimated_bits_values::Vector{Tuple{Int, Int}} = Tuple{Int, Int}[])
    
    # 1. Reset masks
    fill!(W.is_decimated, false)
    fill!(W.is_erased, false)
    
    # 2. Load standard bits
    copyto!(W.channel_bits, v_in)
    
    # 3. Apply Erasures
    for v in erasures
        W.is_erased[v] = true
        W.channel_bits[v] = 0x00 # Fallback guess
    end
    
    # 4. Apply Manual Decimation
    for (v, bit_val) in decimated_bits_values
        W.is_decimated[v] = true
        W.channel_bits[v] = UInt8(bit_val)
    end
    
    # 5. Initialize V2C messages
    copyto!(W.current_bits, W.channel_bits)
    fill!(W.C2V, 0x00)
    
    for v in 1:W.num_var
        val = W.channel_bits[v]
        for e in W.var_to_edges[v]
            W.V2C[e] = val
        end
    end
end

struct SoftDecisionWorkspace{T <: AbstractFloat}
    num_var::Int
    num_check::Int
    num_edges::Int
    
    # 1D Edge Trackers (Size: E)
    V2C::Vector{T}
    C2V::Vector{T}
    
    # Edge Mapping: Which edges belong to which nodes?
    var_to_edges::Vector{Vector{Int}}
    chk_to_edges::Vector{Vector{Int}}
    chk_to_vars::Vector{Vector{Int}} # For syndrome checking and Layered updates
    
    # Node Trackers
    channel_llrs::Vector{T}
    total_llrs::Vector{T}
    current_bits::Vector{UInt8}
    is_decimated::Vector{Bool}

    target_syndrome::Vector{UInt8}

    # Optional Layered Architecture
    layers::Vector{Vector{Int}}

    # Oscillation Tracking
    prev_bits_1::Vector{UInt8}
    prev_bits_2::Vector{UInt8}
end

function init_soft_workspace(H::AbstractMatrix; schedule::Symbol=:flooding)
    num_check, num_var = size(H)
    
    var_to_edges = [Int[] for _ in 1:num_var]
    chk_to_edges = [Int[] for _ in 1:num_check]
    chk_to_vars  = [Int[] for _ in 1:num_check]
    
    edge_id = 1
    for c in 1:num_check
        for v in 1:num_var
            if !iszero(H[c, v])
                push!(chk_to_edges[c], edge_id)
                push!(var_to_edges[v], edge_id)
                push!(chk_to_vars[c], v) 
                
                edge_id += 1
            end
        end
    end
    
    num_edges = edge_id - 1
    
    # Pre-compute the layer schedule if requested, otherwise leave it empty
    layers = if schedule == :layered 
        layered_schedule(H, schedule=:layered) 
    else 
        Vector{Vector{Int}}() 
    end
    
    return SoftDecisionWorkspace{Float64}(
        num_var, num_check, num_edges,
        zeros(Float64, num_edges), zeros(Float64, num_edges),
        var_to_edges, chk_to_edges, chk_to_vars,
        zeros(Float64, num_var), zeros(Float64, num_var), zeros(UInt8, num_var),
        zeros(Bool, num_var), layers
    )
end

"""
Loads a new received vector into the workspace. 
Applies erasures (neutral LLRs) and manual decimation (pinned LLRs) allocation-free.
"""
@inline function load_soft_channel!(W::SoftDecisionWorkspace{Float64}, LLR_in::Vector{Float64};
                                    syndrome::Vector{UInt8} = zeros(UInt8, W.num_check),
                                    erasures::Vector{Int} = Int[],
                                    decimated_bits_values::Vector{Tuple{Int, Int}} = Tuple{Int, Int}[])
    
    # 1. Load Syndrome
    copyto!(W.target_syndrome, syndrome)
    
    # 1. Reset decimation mask
    fill!(W.is_decimated, false)
    
    # 2. Load standard LLRs
    copyto!(W.channel_llrs, LLR_in)
    
    # 3. Apply Erasures (Mathematically, an erasure has exactly 0.0 LLR)
    for v in erasures
        W.channel_llrs[v] = 0.0
    end
    
    # 4. Apply Manual Decimation (Overwrites erasures if they overlap)
    for (v, bit_val) in decimated_bits_values
        W.is_decimated[v] = true
        # Pin to a massive LLR: +1000.0 means 100% confident it's a 0, -1000.0 means 1.
        llr_val = iszero(bit_val) ? 1000.0 : -1000.0
        W.channel_llrs[v] = llr_val
    end
    
    # 5. Initialize Total LLRs and V2C messages for Iteration 1
    copyto!(W.total_llrs, W.channel_llrs)
    fill!(W.C2V, 0.0)
    
    for v in 1:W.num_var
        val = W.channel_llrs[v]
        for e in W.var_to_edges[v]
            W.V2C[e] = val
        end
    end
end

#############################
         # Gallager
#############################

"""
The universal entry point for Hard-Decision decoding (Gallager A and B).
For Gallager A, set `Bt` equal to `d_v - 1` (where `d_v` is the variable node degree).
"""
function decode!(W::HardDecisionWorkspace, bits_in::Vector{UInt8};
                 max_iter::Int = 100,
                 Bt::Int = 2,
                 syndrome::Vector{UInt8} = zeros(UInt8, W.num_check),
                 erasures::Vector{Int} = Int[],
                 decimated_bits_values::Vector{Tuple{Int, Int}} = Tuple{Int, Int}[])
    
    # 1. Snap the channel data, syndrome, and masks into place
    copyto!(W.target_syndrome, syndrome)
    load_hard_channel!(W, bits_in, erasures=erasures, decimated_bits_values=decimated_bits_values)
    
    # 2. Fire the branchless math engine
    return _fast_decode!(W, max_iter, Bt)
end

# ==============================================================================
# HARD DECISION ENGINE: GALLAGER B (AND A)
# ==============================================================================

function _fast_decode!(W::HardDecisionWorkspace, max_iter::Int, Bt::Int)
    
    for iter in 1:max_iter
        # ---------------------------------------------------------
        # 1. Check Node Update (The Parity Trick & Syndrome Injection)
        # ---------------------------------------------------------
        for c in 1:W.num_check
            # Initialize parity to the target syndrome (0x00 for standard decoding)
            total_parity = W.target_syndrome[c]
            
            for e in W.chk_to_edges[c]
                total_parity ⊻= W.V2C[e]
            end
            
            for e in W.chk_to_edges[c]
                W.C2V[e] = total_parity ⊻ W.V2C[e]
            end
        end
        
        # ---------------------------------------------------------
        # 2. Variable Node Update (Disagreement Trick & Decimation Guard)
        # ---------------------------------------------------------
        for v in 1:W.num_var
            # If decimated, node is permanently frozen to its pinned channel bit
            if W.is_decimated[v]
                for e in W.var_to_edges[v]
                    W.V2C[e] = W.channel_bits[v]
                end
                W.current_bits[v] = W.channel_bits[v]
                continue
            end
            
            y_v = W.channel_bits[v]
            disagreements = 0
            
            # Count incoming messages that disagree with the channel bit
            for e in W.var_to_edges[v]
                if W.C2V[e] != y_v
                    disagreements += 1
                end
            end
            
            # Compute outgoing V2C messages
            for e in W.var_to_edges[v]
                disagree_excluding_e = disagreements - (W.C2V[e] != y_v ? 1 : 0)
                if disagree_excluding_e >= Bt
                    W.V2C[e] = y_v ⊻ 0x01
                else
                    W.V2C[e] = y_v
                end
            end
            
            # Tentative hard decision
            if disagreements >= Bt
                W.current_bits[v] = y_v ⊻ 0x01
            else
                W.current_bits[v] = y_v
            end
        end
        
        # ---------------------------------------------------------
        # 3. Syndrome Convergence Check
        # ---------------------------------------------------------
        is_valid = true
        for c in 1:W.num_check
            syn = W.target_syndrome[c]
            for v in W.chk_to_vars[c]
                syn ⊻= W.current_bits[v]
            end
            
            if syn != 0x00
                is_valid = false
                break
            end
        end
        
        if is_valid
            return true, W.current_bits, iter
        end
    end
    
    return false, W.current_bits, max_iter
end

"""
The Exact Box-Plus Operator (Jacobian Logarithm) for Sum-Product.
Mathematically equivalent to the tanh rule, but numerically bulletproof against NaNs.
"""
@inline function boxplus_exact(x::Float64, y::Float64)
    # The Min-Sum base: sign(x)*sign(y)*min(|x|, |y|)
    base = sign(x) * sign(y) * min(abs(x), abs(y))
    
    # The Sum-Product correction factor: log(1 + e^-|x+y|) - log(1 + e^-|x-y|)
    # log1p(z) accurately computes log(1+z) even for extremely small z
    corr = log1p(exp(-abs(x + y))) - log1p(exp(-abs(x - y)))
    
    return base + corr
end

"""
The Min-Sum Box-Plus Operator.
"""
@inline function boxplus_minsum(x::Float64, y::Float64)
    return sign(x) * sign(y) * min(abs(x), abs(y))
end

"""
The Min-Sum Box-Plus Operator with a low-complexity correction term.
"""
@inline function boxplus_minsum_correction(x::Float64, y::Float64)
    # The standard Min-Sum base
    base = sign(x) * sign(y) * min(abs(x), abs(y))
    
    # Pre-compute to avoid redundant absolute value evaluations
    sum_abs = abs(x + y)
    diff_abs = abs(x - y)
    
    # Low-complexity approximation to the correction term
    corr = 0.0
    if sum_abs < 2.0 && diff_abs > 2.0 * sum_abs
        corr = 0.5
    elseif diff_abs < 2.0 && sum_abs > 2.0 * diff_abs
        corr = -0.5
    end
    
    return base + corr
end

function decode!(W::SoftDecisionWorkspace{Float64}, LLR_in::Vector{Float64};
                 algorithm::Symbol = :offset_min_sum,
                 schedule::Symbol = :layered,
                 decimation::Symbol = :none,
                 oscillation::Symbol = :none,
                 perturbation::Symbol = :none,  # NEW: Default to :none
                 max_iter::Int = 100,
                 attenuation::Float64 = 0.75,
                 offset::Float64 = 0.5,
                 dec_thresh::Float64 = 10.0,
                 dec_rounds::Int = 10,
                 perturbation_scale::Float64 = 0.0,
                 syndrome::Vector{UInt8} = zeros(UInt8, W.num_check),
                 erasures::Vector{Int} = Int[],
                 decimated_bits_values::Vector{Tuple{Int, Int}} = Tuple{Int, Int}[])
    
    copyto!(W.target_syndrome, syndrome)
    load_soft_channel!(W, LLR_in, erasures=erasures, decimated_bits_values=decimated_bits_values)
    
    target_schedule = schedule == :serial ? :layered : schedule
    
    return _fast_decode!(W, Val(algorithm), Val(target_schedule), Val(decimation), Val(oscillation), Val(perturbation),
                         max_iter, attenuation, offset, dec_thresh, dec_rounds, perturbation_scale)
end

# ==============================================================================
# SOFT DECISION ENGINE: FLOODING (PARALLEL SCHEDULE)
# ==============================================================================

function _fast_decode!(W::SoftDecisionWorkspace{Float64}, 
                       algo::Val, 
                       ::Val{:flooding}, 
                       decimation_type::Val, 
                       osc_type::Val, 
                       pert_type::Val, 
                       max_iter::Int, 
                       attenuation::Float64, 
                       offset::Float64, 
                       dec_thresh::Float64, 
                       dec_rounds::Int, 
                       perturbation_scale::Float64)
    
    @inbounds begin
        for iter in 1:max_iter
            # ---------------------------------------------------------
            # 1. Variable Node Update (V2C Calculation)
            # ---------------------------------------------------------
            for v in 1:W.num_var
                for e in W.var_to_edges[v]
                    W.V2C[e] = W.total_llrs[v] - W.C2V[e]
                end
            end
            
            # ---------------------------------------------------------
            # 2. Check Node Update (C2V Calculation)
            # ---------------------------------------------------------
            for c in 1:W.num_check
                for e_out in W.chk_to_edges[c]
                    agg = 0.0
                    first_val = true
                    for e_in in W.chk_to_edges[c]
                        if e_in != e_out
                            if first_val
                                agg = W.V2C[e_in]
                                first_val = false
                            else
                                agg = _apply_boxplus(algo, agg, W.V2C[e_in])
                            end
                        end
                    end
                    
                    # Apply attenuation / offset
                    agg = _apply_post_process(algo, agg, attenuation, offset)
                    
                    # Inject Syndrome Math: Flip sign if parity check requires an odd number of 1s
                    if W.target_syndrome[c] == 0x01
                        agg = -agg
                    end
                    
                    W.C2V[e_out] = agg
                end
            end
            
            # ---------------------------------------------------------
            # 3. Total LLR Update & Hard Decision
            # ---------------------------------------------------------
            for v in 1:W.num_var
                # Only update the total if it hasn't been pinned by decimation
                if !W.is_decimated[v]
                    tot = W.channel_llrs[v]
                    for e in W.var_to_edges[v]
                        tot += W.C2V[e]
                    end
                    W.total_llrs[v] = tot
                end
                
                W.current_bits[v] = W.total_llrs[v] < 0.0 ? 0x01 : 0x00
            end
            
            # ---------------------------------------------------------
            # 4. End-of-Iteration Hooks & Convergence Check
            # ---------------------------------------------------------
            _apply_decimation!(decimation_type, W, iter, dec_thresh, dec_rounds)
            
            is_valid = true
            for c in 1:W.num_check
                syn = W.target_syndrome[c]
                for v in W.chk_to_vars[c]
                    syn ⊻= W.current_bits[v]
                end
                if syn != 0x00
                    is_valid = false
                    break
                end
            end
            
            if is_valid
                return true, W.current_bits, iter
            end

            # ---------------------------------------------------------
            # OPTIONAL OSCILLATION DETECTION & RESCUE
            # ---------------------------------------------------------
            is_oscillating = _check_oscillation(osc_type, W)
            
            if is_oscillating
                if pert_type != Val{:none}
                    _apply_perturbation!(pert_type, W, iter, perturbation_scale, force=true)
                    _wipe_history!(osc_type, W)
                    
                elseif decimation_type != Val{:none} && decimation_type != Val{:manual}
                    _apply_decimation!(decimation_type, W, iter, dec_thresh, dec_rounds, force=true)
                    _wipe_history!(osc_type, W)
                    
                else
                    return false, W.current_bits, -iter 
                end
            else
                _update_history!(osc_type, W)
                
                # The compiler deletes this entirely if pert_type == Val{:none}
                _apply_perturbation!(pert_type, W, iter, perturbation_scale, force=false)
                _apply_decimation!(decimation_type, W, iter, dec_thresh, dec_rounds, force=false)
            end
        end
    end

    return false, W.current_bits, max_iter
end

# ==============================================================================
# SOFT DECISION ENGINE: LAYERED (& SERIAL) SCHEDULE
# ==============================================================================

function _fast_decode!(W::SoftDecisionWorkspace{Float64}, 
                       algo::Val, 
                       ::Val{:flooding}, 
                       decimation_type::Val, 
                       osc_type::Val, 
                       pert_type::Val, 
                       max_iter::Int, 
                       attenuation::Float64, 
                       offset::Float64, 
                       dec_thresh::Float64, 
                       dec_rounds::Int, 
                       perturbation_scale::Float64)
    
    @inbounds begin
        for iter in 1:max_iter
            # ---------------------------------------------------------
            # 1. Variable Node Update (V2C Calculation)
            # ---------------------------------------------------------
            for v in 1:W.num_var
                for e in W.var_to_edges[v]
                    W.V2C[e] = W.total_llrs[v] - W.C2V[e]
                end
            end
            
            # ---------------------------------------------------------
            # 2. Check Node Update (C2V Calculation)
            # ---------------------------------------------------------
            for c in 1:W.num_check
                for e_out in W.chk_to_edges[c]
                    agg = 0.0
                    first_val = true
                    for e_in in W.chk_to_edges[c]
                        if e_in != e_out
                            if first_val
                                agg = W.V2C[e_in]
                                first_val = false
                            else
                                agg = _apply_boxplus(algo, agg, W.V2C[e_in])
                            end
                        end
                    end
                    
                    # Apply attenuation / offset
                    agg = _apply_post_process(algo, agg, attenuation, offset)
                    
                    # Inject Syndrome Math: Flip sign if parity check requires an odd number of 1s
                    if W.target_syndrome[c] == 0x01
                        agg = -agg
                    end
                    
                    W.C2V[e_out] = agg
                end
            end
            
            # ---------------------------------------------------------
            # 3. Total LLR Update & Hard Decision (Split for SIMD)
            # ---------------------------------------------------------
            # Loop A: Indirect memory access (No SIMD)
            for v in 1:W.num_var
                # Only update the total if it hasn't been pinned by decimation
                if !W.is_decimated[v]
                    tot = W.channel_llrs[v]
                    for e in W.var_to_edges[v]
                        tot += W.C2V[e]
                    end
                    W.total_llrs[v] = tot
                end
            end
            
            # Loop B: Perfectly linear, contiguous sweep (SIMD applied)
            @simd for v in 1:W.num_var
                W.current_bits[v] = W.total_llrs[v] < 0.0 ? 0x01 : 0x00
            end
            
            # ---------------------------------------------------------
            # 4. End-of-Iteration Hooks & Convergence Check
            # ---------------------------------------------------------
            _apply_decimation!(decimation_type, W, iter, dec_thresh, dec_rounds)
            
            is_valid = true
            for c in 1:W.num_check
                syn = W.target_syndrome[c]
                for v in W.chk_to_vars[c]
                    syn ⊻= W.current_bits[v]
                end
                if syn != 0x00
                    is_valid = false
                    break
                end
            end
            
            if is_valid
                return true, W.current_bits, iter
            end

            # ---------------------------------------------------------
            # OPTIONAL OSCILLATION DETECTION & RESCUE
            # ---------------------------------------------------------
            is_oscillating = _check_oscillation(osc_type, W)
            
            if is_oscillating
                if pert_type != Val{:none}
                    _apply_perturbation!(pert_type, W, iter, perturbation_scale, force=true)
                    _wipe_history!(osc_type, W)
                    
                elseif decimation_type != Val{:none} && decimation_type != Val{:manual}
                    _apply_decimation!(decimation_type, W, iter, dec_thresh, dec_rounds, force=true)
                    _wipe_history!(osc_type, W)
                    
                else
                    return false, W.current_bits, -iter 
                end
            else
                _update_history!(osc_type, W)
                
                # The compiler deletes this entirely if pert_type == Val{:none}
                _apply_perturbation!(pert_type, W, iter, perturbation_scale, force=false)
                _apply_decimation!(decimation_type, W, iter, dec_thresh, dec_rounds, force=false)
            end
        end
    end
    
    return false, W.current_bits, max_iter
end

# ==============================================================================
# ZERO-OVERHEAD DECIMATION HOOKS
# ==============================================================================

"""
Auto Decimation: 
If a node's belief exceeds a threshold, permanently pin it.
"""
@inline function _apply_decimation!(::Val{:auto}, W::SoftDecisionWorkspace, iter::Int, threshold::Float64, rounds::Int)
    @inbounds for v in 1:W.num_var
        if !W.is_decimated[v]
            if W.total_llrs[v] > threshold
                W.is_decimated[v] = true
                W.total_llrs[v] = 1000.0   # Pin to strong 0
                W.channel_llrs[v] = 1000.0
            elseif W.total_llrs[v] < -threshold
                W.is_decimated[v] = true
                W.total_llrs[v] = -1000.0  # Pin to strong 1
                W.channel_llrs[v] = -1000.0
            end
        end
    end
end

"""
Guided Decimation: 
Every `rounds` iterations, find the highest confidence undeclared node and pin it.
"""
@inline function _apply_decimation!(::Val{:guided}, W::SoftDecisionWorkspace, iter::Int, threshold::Float64, rounds::Int)
    if iter % rounds == 0
        best_v = -1
        max_belief = -1.0
        
        # Find the node with the highest absolute belief that isn't already decimated
        @inbounds for v in 1:W.num_var
            if !W.is_decimated[v]
                abs_belief = abs(W.total_llrs[v])
                if abs_belief > max_belief
                    max_belief = abs_belief
                    best_v = v
                end
            end
        end
        
        # Decimate the winner
        if best_v != -1
            W.is_decimated[best_v] = true
            if W.total_llrs[best_v] >= 0
                W.total_llrs[best_v] = 1000.0
                W.channel_llrs[best_v] = 1000.0
            else
                W.total_llrs[best_v] = -1000.0
                W.channel_llrs[best_v] = -1000.0
            end
        end
    end
end

# ==============================================================================
# 1. THE BOX-PLUS FOLDING TRAITS (Inner Loop)
# ==============================================================================
@inline _apply_boxplus(::Val{:sum_product}, x::Float64, y::Float64) = boxplus_exact(x, y)
@inline _apply_boxplus(::Val{:min_sum}, x::Float64, y::Float64) = boxplus_minsum(x, y)
@inline _apply_boxplus(::Val{:min_sum_correction}, x::Float64, y::Float64) = boxplus_minsum_correction(x, y)

# Normalized Min-Sum uses the exact same folding logic as standard Min-Sum!
@inline _apply_boxplus(::Val{:normalized_min_sum}, x::Float64, y::Float64) = boxplus_minsum(x, y)
@inline _apply_boxplus(::Val{:offset_min_sum}, x::Float64, y::Float64) = boxplus_minsum(x, y)

# ==============================================================================
# 2. THE POST-PROCESSING TRAITS (Applied once at the end)
# ==============================================================================
# Base cases (Do nothing)
@inline _apply_post_process(::Val{:sum_product}, agg::Float64, alpha::Float64, beta::Float64) = agg
@inline _apply_post_process(::Val{:min_sum}, agg::Float64, alpha::Float64, beta::Float64) = agg
@inline _apply_post_process(::Val{:min_sum_correction}, agg::Float64, alpha::Float64, beta::Float64) = agg

# Normalized Min-Sum (Uses alpha)
@inline _apply_post_process(::Val{:normalized_min_sum}, agg::Float64, alpha::Float64, beta::Float64) = agg * alpha

# Offset Min-Sum (Uses beta)
@inline _apply_post_process(::Val{:offset_min_sum}, agg::Float64, alpha::Float64, beta::Float64) = sign(agg) * max(0.0, abs(agg) - beta)

# No decimation: The compiler completely erases this call.
@inline _apply_decimation!(::Val{:none}, W, iter, threshold, rounds) = nothing

"""
Manual Decimation: 
Handled during initialization. Nodes are pinned to +/- 1000.0 and skipped.
"""
@inline _apply_decimation!(::Val{:manual}, W, iter, threshold, rounds) = nothing

# ==============================================================================
# PERTURBATION TRAITS
# ==============================================================================

# When disabled, the compiler completely erases this call
@inline _apply_perturbation!(::Val{:none}, W, iter, scale; force=false) = nothing

# When active, perform the stochastic noise injection
@inline function _apply_perturbation!(::Val{:active}, W::SoftDecisionWorkspace, iter::Int, scale::Float64; force::Bool=false)
    if force || iter % 10 == 0
        # can't @simd with randn() unless we pre-allocate a noise buffer
        @inbounds for v in 1:W.num_var
            if !W.is_decimated[v]
                W.total_llrs[v] += scale * randn()
            end
        end
    end
end

# ==============================================================================
# OSCILLATION TRAITS
# ==============================================================================

# When disabled, compiler returns false and deletes the history updates
@inline _check_oscillation(::Val{:none}, W) = false
@inline _update_history!(::Val{:none}, W) = nothing
@inline _wipe_history!(::Val{:none}, W) = nothing

# When active, perform the memory checks and shifts
@inline _check_oscillation(::Val{:active}, W) = (W.current_bits == W.prev_bits_1 || W.current_bits == W.prev_bits_2)

@inline function _update_history!(::Val{:active}, W)
    copyto!(W.prev_bits_2, W.prev_bits_1)
    copyto!(W.prev_bits_1, W.current_bits)
end

@inline function _wipe_history!(::Val{:active}, W)
    fill!(W.prev_bits_1, 0xFF)
    fill!(W.prev_bits_2, 0xFF)
end

#############################
          # Methods
#############################

# Mansour, Shanbhag, "Turbo Decoder Architectures for Low-Density Parity-Check Codes" (2002)

# A layer is a collection of check-nodes such that any two check-nodes have no neighbouring variable-node in common.
# TODO latex the H
"""
    layered_schedule(H::CTMatrixTypes; schedule::Symbol = :layered, random::Bool = false)

Return a layered schedule for the parity-check matrix `H`. If `schedule` is `:parallel` or
`:serial`, layers representing these two extreme cases are returned. If `random` is `true`, the
schedule is shuffled.
"""
function layered_schedule(H::CTMatrixTypes; schedule::Symbol = :layered, random::Bool = false)
    num_check, num_var = size(H)
    num_check > 0 && num_var > 0 || throw(ArgumentError("Input matrix of improper dimension"))
    schedule ∈ (:flooding, :parallel, :serial, :layered, :semiserial) || 
        throw(ArgumentError("Unknown schedule algorithm"))
    schedule == :flooding && (schedule = :parallel;)
    schedule == :semiserial && (schedule = :layered;)

    if schedule == :layered
        check_adj_list = [Int[] for _ in 1:num_check]
        for r in 1:num_check
            for c in 1:num_var
                iszero(H[r, c]) || push!(check_adj_list[r], c)
            end
        end

        sched_list = [[1]]
        list = collect(2:num_check)
        random && shuffle!(list)
        for c in list
            found = false
            for sched in sched_list
                if !any(x ∈ check_adj_list[y] for y in sched for x ∈ check_adj_list[c])
                    push!(sched, c)
                    sort!(sched_list, lt = (x, y) -> length(x) < length(y))
                    found = true
                    break
                end
            end
            !found && push!(sched_list, [c])
        end
        random && shuffle!(sched_list)
    elseif schedule == :parallel
        sched_list = [collect(1:num_check)]
        random && shuffle!(sched_list[1])
    else
        # serial
        sched_list = [[i] for i in 1:num_check]
        random && shuffle!(sched_list)
    end
    return sched_list
end
# TODO LDPCCode version

# ref: Layered Decoding of Quantum LDPC Codes
function balance_of_layered_schedule(sch::Vector{Vector{Int}})
    is_empty(sch) && throw(ArgumentError("Schedule cannot be empty"))
    any(x -> is_empty(x), sch) && throw(ArgumentError("Schedule cannot contain an empty layer"))

    len = sch[1]
    all(x -> length(x) == len, sch) && return 1
    γ = 0.0
    for L_i in sch
        for L_j in sch
            temp = length(L_i) / length(L_j)
            temp > γ && (γ = temp;)
        end
    end
    return γ
end
