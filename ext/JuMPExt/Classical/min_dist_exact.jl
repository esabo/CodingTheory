# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
    _minimum_distance_ILP(C::AbstractLinearCode; verbose::Bool = false, time_limit_sec::Float64 = 300.0)

Return the exact minimum distance and witness codeword using an Integer Linear Programming (ILP) formulation.
Intercepts the NP-hard search space by mapping the parity constraints to integer multiples of `q`.
"""
function _minimum_distance_ILP(C::AbstractLinearCode; verbose::Bool = false, time_limit_sec::Float64 = 300.0)
    r, n = C.n - C.k, C.n
    q = Int(order(C.F))
    
    # Extract integer representation of the parity check matrix
    H_mat = parity_check_matrix(C)
    H = zeros(Int, r, n)
    for i in 1:r, j in 1:n
        H[i, j] = CodingTheory._is_binary(C) ? Int(H_mat[i, j]) : CodingTheory._field_elem_to_int(H_mat[i, j])
    end
    
    verbose && println("Formulating ILP model for [$n, $(C.k)] code over GF($q)...")
    model = Model(GLPK.Optimizer)
    set_optimizer_attribute(model, "tm_lim", round(Int, time_limit_sec * 1000)) 
    
    if !verbose
        set_silent(model)
    else
        unset_silent(model)
    end
    
    # --- VARIABLES ---
    @variable(model, x[1:n], Bin)
    @variable(model, 0 <= v[1:n] <= q - 1, Int)
    @variable(model, z[1:r], Int)
    
    # --- OBJECTIVE ---
    @objective(model, Min, sum(x[i] for i in 1:n))
    
    # --- CONSTRAINTS ---
    @constraint(model, sum(x[i] for i in 1:n) >= 1)
    
    for i in 1:n
        @constraint(model, v[i] <= (q - 1) * x[i])
        @constraint(model, v[i] >= x[i])
    end
    
    for j in 1:r
        @constraint(model, sum(H[j, i] * v[i] for i in 1:n) == q * z[j])
    end
    
    # --- SOLVE ---
    verbose && println("Handing off to GLPK solver (Time limit: $(time_limit_sec)s)...")
    optimize!(model)
    
    status = termination_status(model)
    
    if status == MOI.OPTIMAL
        d = Int(round(objective_value(model)))
        # Extract and cast the witness vector back to the finite field
        witness_raw = Int.(round.(value.(v)))
        witness = matrix(C.F, 1, n, [C.F(val) for val in witness_raw])
        
        verbose && println("GLPK found exact optimum: d = $d")
        return d, witness
    elseif status == MOI.TIME_LIMIT
        verbose && println("Solver hit the time limit of $(time_limit_sec)s before proving optimality.")
        return -1, zero_matrix(C.F, 1, n)
    else
        verbose && println("Solver terminated early with status: $status")
        return -1, zero_matrix(C.F, 1, n)
    end
end

"""
    _fractional_distance_LP(C::AbstractLinearCode; verbose::Bool = false)

Calculates a tight lower bound on the fractional distance of an LDPC code by minimizing 
the L1 norm over the facets of the fundamental polytope not containing the all-zero vertex (Burshtein & Goldenberg, 2011).
"""
function _fractional_distance_LP(C::AbstractLinearCode; verbose::Bool = false)
    r, n = C.n - C.k, C.n
    H = CodingTheory._convert_binary_to_int_matrix(parity_check_matrix(C))
    
    verbose && println("Formulating LP relaxation for fractional distance bounding...")
    
    # We compute d_frac^(1) by taking the minimum of the L1 norm over the facets {c_i = 1}
    min_fractional_dist = Float64(n + 1)
    
    for target_i in 1:n
        model = Model(GLPK.Optimizer)
        set_silent(model)
        
        # Continuous variables for the LP relaxation
        @variable(model, 0 <= c[1:n] <= 1)
        @objective(model, Min, sum(c[i] for i in 1:n))
        
        # Facet constraint: Force the target variable to 1
        @constraint(model, c[target_i] == 1)
        
        # Fundamental Polytope Constraints for regular/plain LDPC
        for j in 1:r
            neighbors = findall(x -> x != 0, H[j, :])
            # For each odd-sized subset of neighbors S
            for S_size in 1:2:length(neighbors)
                for S in Combinatorics.combinations(neighbors, S_size)
                    not_S = setdiff(neighbors, S)
                    @constraint(model, sum(c[i] for i in not_S) + sum(1 - c[i] for i in S) >= 1)
                end
            end
        end
        
        optimize!(model)
        
        if termination_status(model) == MOI.OPTIMAL
            val = objective_value(model)
            if val < min_fractional_dist
                min_fractional_dist = val
            end
        end
    end
    
    verbose && println("LP Relaxation found fractional distance lower bound: $min_fractional_dist")
    return min_fractional_dist
end

# ==============================================================================
# DISPATCHER EXTENSIONS
# ==============================================================================

"""
    minimum_distance_ilp(C::AbstractLinearCode; time_limit_sec::Float64 = 300.0, verbose::Bool = false)

Directly invoke the Integer Linear Programming (ILP) solver to find the exact minimum distance.
Bypasses the standard exact combinatorial solvers, which is highly recommended for exceptionally sparse LDPC codes.
"""
function minimum_distance_ilp(C::AbstractLinearCode; time_limit_sec::Float64 = 300.0, verbose::Bool = false)
    !ismissing(C.d) && return C.d, (isdefined(C, :witness) ? C.witness : zero_matrix(C.F, 1, C.n))
    
    d_ilp, witness_ilp = _minimum_distance_ILP(C; verbose=verbose, time_limit_sec=time_limit_sec)
    
    if d_ilp != -1
        C.d = d_ilp
        C.u_bound = d_ilp
        C.l_bound = d_ilp
        return d_ilp, witness_ilp
    else
        error("ILP Solver failed to find the exact minimum distance within the allotted constraints.")
    end
end

"""
    fractional_distance_bound(C::AbstractLinearCode; verbose::Bool = false)

Estimates the fractional distance by searching the fundamental polytope via Linear Programming (LP).
Returns a lower bound on the fractional distance, which mathematically guarantees the error-correction 
capability of an LP decoder on the given code.
"""
function fractional_distance_bound(C::AbstractLinearCode; verbose::Bool = false)
    if !CodingTheory._is_binary(C)
        error("Fractional distance LP bounding is currently only implemented for binary codes.")
    end
    return _fractional_distance_LP(C; verbose=verbose)
end
