# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
    _minimum_distance_ILP(C::AbstractLinearCode; verbose::Bool = false)

Return the minimum distance of the linear code using an integer linear programming approach.

# Note
- Run `using JuMP, GLPK` to activate this extension.
"""
function _minimum_distance_ILP(C::AbstractLinearCode; verbose::Bool = false)
    r, n = C.n - C.k, C.n
    q = Int(order(C.F))
    
    # We must lift the parity-check matrix out of the finite field into standard Integers
    # Assuming parity_check_matrix(C) returns a matrix we can cast to Int
    H = Int.(Array(parity_check_matrix(C))) 
    
    # Initialize the GLPK model
    model = Model(GLPK.Optimizer)
    if !verbose
        set_silent(model)
    end
    
    # --- VARIABLES ---
    # x[i] is the binary indicator for the Hamming weight (1 if non-zero, 0 if zero)
    @variable(model, x[1:n], Bin)
    
    # v[i] is the actual symbol value in the physical codeword
    @variable(model, 0 <= v[1:n] <= q - 1, Int)
    
    # z[j] is the unconstrained auxiliary multiplier to handle the modulo q constraint
    @variable(model, z[1:r], Int)
    
    # --- OBJECTIVE ---
    # We want the Minimum Hamming Weight
    @objective(model, Min, sum(x[i] for i in 1:n))
    
    # --- CONSTRAINTS ---
    # 1. The codeword cannot be the all-zero vector
    @constraint(model, sum(x[i] for i in 1:n) >= 1)
    
    # 2. Link the indicator variables x[i] to the values v[i]
    for i in 1:n
        # If x[i] == 0, then v[i] MUST be 0
        @constraint(model, v[i] <= (q - 1) * x[i])
        
        # If x[i] == 1, then v[i] MUST be at least 1
        @constraint(model, v[i] >= x[i])
    end
    
    # 3. The Parity-Check Equations (H * v = 0 mod q)
    for j in 1:r
        # Instead of modulo, we enforce that the dot product is exactly a multiple of q
        @constraint(model, sum(H[j, i] * v[i] for i in 1:n) == q * z[j])
    end
    
    # --- SOLVE ---
    verbose && println("Passing ILP model to GLPK solver...")
    optimize!(model)
    
    # --- PARSE RESULTS ---
    status = termination_status(model)
    
    if status == MOI.OPTIMAL
        d = Int(round(objective_value(model)))
        verbose && println("GLPK found exact optimum: d = $d")
        return d
    elseif status == MOI.INFEASIBLE
        verbose && println("Model is mathematically infeasible (should not happen for valid codes).")
        return -1
    else
        verbose && println("Solver terminated early with status: $status")
        return -1
    end
end
