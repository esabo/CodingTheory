"""
$(TYPEDSIGNATURES)

Return the theoretical maximum minimum distance for a self-dual code 
under Gleason's Theorems. Returns `missing` if the code does not fit a standard extremal classification.
"""
function Gleason_bound(C::AbstractLinearCode)
    if !is_self_dual(C)
        return missing
    end
    
    q = Int(order(C.F))
    n = C.n
    
    if q == 2
        if is_doubly_even(C)
            # Type II (Doubly-Even Binary Self-Dual)
            # Exists only if n is a multiple of 8
            return n % 8 == 0 ? 4 * fld(n, 24) + 4 : missing
        else
            # Type I (Singly-Even Binary Self-Dual)
            # Exists only if n is a multiple of 2
            return n % 2 == 0 ? 2 * fld(n, 8) + 2 : missing
        end
    elseif q == 3
        # Ternary Self-Dual
        # Exists only if n is a multiple of 4
        return n % 4 == 0 ? 3 * fld(n, 12) + 3 : missing
    elseif q == 4
        # Hermitian Self-Dual Quaternary (Quantum preparation)
        # Exists only if n is even
        return n % 2 == 0 ? 2 * fld(n, 6) + 2 : missing
    end
    
    return missing
end

"""
$(TYPEDSIGNATURES)

Return `true` if a self-dual code attains the maximum minimum distance allowed
by Gleason's theorems; otherwise, return `false`.
"""
function is_extremal(C::AbstractLinearCode; verbose::Bool=false)
    bound = Gleason_bound(C)
    
    if ismissing(bound)
        verbose && println("Code does not fit a standard extremal Gleason classification.")
        return false
    end
    
    d = ismissing(C.d) ? minimum_distance(C) : C.d
    is_ext = (d == bound)
    
    if verbose
        println("--- Extremality Check ---")
        println("Theoretical Gleason Bound: d_max = $bound")
        println("Actual Minimum Distance: d = $d")
        println("Is Extremal: $is_ext")
    end
    
    return is_ext
end