# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
    _minimum_distance_ILP(C::AbstractLinearCode; verbose::Bool = false, time_limit_sec::Float64 = 300.0)

Return the minimum distance of the linear code using an integer linear programming approach.

# Note
- Run `using JuMP, GLPK` to activate this extension.
"""
function _minimum_distance_ILP(C::AbstractLinearCode; verbose::Bool = false, time_limit_sec::Float64 = 300.0)
    r, n = C.n - C.k, C.n
    q = Int(order(C.F))
    H = CodingTheory._convert_binary_to_int_matrix(parity_check_matrix(C))
    
    verbose && println("Formulating ILP model for [$n, $(C.k)] code over GF($q)...")
    
    # Initialize the GLPK model
    model = Model(GLPK.Optimizer)
    
    # Set the time limit so it doesn't hang forever on dense codes
    set_optimizer_attribute(model, "tm_lim", round(Int, time_limit_sec * 1000)) # GLPK uses milliseconds
    
    if !verbose
        set_silent(model)
    else
        unset_silent(model) # Explicitly tell GLPK to print branch-and-cut progress
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
    verbose && println("Model generated with $(num_variables(model)) variables and $(num_constraints(model, VariableRef, MOI.Integer)) integer constraints.")
    verbose && println("Handing off to GLPK solver (Time limit: $(time_limit_sec)s)...")
    
    optimize!(model)
    
    # --- PARSE RESULTS ---
    status = termination_status(model)
    
    if status == MOI.OPTIMAL
        d = Int(round(objective_value(model)))
        verbose && println("GLPK found exact optimum: d = $d")
        return d
    elseif status == MOI.TIME_LIMIT
        verbose && println("Solver hit the time limit of $(time_limit_sec)s before proving optimality.")
        return -1
    elseif status == MOI.INFEASIBLE
        verbose && println("Model is mathematically infeasible (should not happen for valid codes).")
        return -1
    else
        verbose && println("Solver terminated early with status: $status")
        return -1
    end
end
