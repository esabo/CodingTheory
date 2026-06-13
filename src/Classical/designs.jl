"""
    is_design_holder(C::AbstractLinearCode, t::Int; verbose::Bool=false)

Evaluates the Assmus-Mattson theorem to determine if the codewords of `C` 
form a `t`-design.
"""
function is_design_holder(C::AbstractLinearCode, t::Int; verbose::Bool=false)
    # Ensure minimum distance is known (assuming you have a minimum_distance getter)
    d = ismissing(C.d) ? minimum_distance(C) : C.d
    
    # Lazily fetch or compute the dual weight distribution
    dual_C = dual(C)
    dual_wd = weight_distribution(dual_C, verbose=false)
    
    # Count the number of non-zero weights in C^⟂ up to n - t
    w = 0
    for (wt, count) in dual_wd
        if wt > 0 && wt <= C.n - t && count > 0
            w += 1
        end
    end
    
    holds_design = w < (d - t)
    
    if verbose
        println("--- Assmus-Mattson Check (t = $t) ---")
        println("Primal min distance (d): $d")
        println("Non-zero dual weights ≤ $(C.n - t) (w): $w")
        println("Condition (w < d - t): $w < $(d - t) -> $(holds_design)")
    end
    
    return holds_design
end

"""
    design_strength(C::AbstractLinearCode; verbose::Bool=false)

Finds the maximum strength `t` for which the Assmus-Mattson theorem 
guarantees the code forms a `t`-design.
"""
function design_strength(C::AbstractLinearCode; verbose::Bool=false)
    d = ismissing(C.d) ? minimum_distance(C) : C.d
    max_t = 0
    
    # The theorem requires w < d - t, which inherently implies t < d
    for t in 1:(d - 1)
        if is_design_holder(C, t, verbose=false)
            max_t = t
        else
            break # If it fails for t, it will fail for t + 1
        end
    end
    
    verbose && println("Maximum guaranteed design strength: t = $max_t")
    return max_t
end
