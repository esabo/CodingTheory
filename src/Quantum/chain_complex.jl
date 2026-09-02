import Oscar: rank, nrows, ncols, iszero, zero_matrix, base_ring
import Base: getindex

"""
    ChainComplex{T <: CTMatrixTypes}

A homological chain complex built on Oscar matrix types.
`degrees` specifies the homological grading (e.g., `[4, 3, 2]`).
`maps` contains the boundary matrices. 
"""
struct ChainComplex{T <: CTMatrixTypes}
    degrees::Vector{Int}
    maps::Vector{T}
    
    function ChainComplex(degrees::Vector{Int}, maps::Vector{T}) where {T <: CTMatrixTypes}
        if length(degrees) != length(maps) + 1
            throw(ArgumentError("A complex with $(length(degrees)) spaces must have exactly $(length(degrees)-1) boundary maps."))
        end
        
        for i in 1:(length(maps) - 1)
            # maps[i] is ∂_k : C_k -> C_{k-1}
            # maps[i+1] is ∂_{k-1} : C_{k-1} -> C_{k-2}
            if nrows(maps[i]) != ncols(maps[i+1])
                throw(DimensionMismatch("Target dimension of maps[$i] does not match domain of maps[$(i+1)]."))
            end
            
            # Homological condition: ∂_{k-1} ∘ ∂_k == 0
            if !iszero(maps[i+1] * maps[i])
                throw(ArgumentError("Homology condition failed: maps[$(i+1)] * maps[$i] != 0"))
            end
        end
        
        new{T}(degrees, maps)
    end
end

"""
    Base.getindex(C::ChainComplex, k::Int)

Returns the boundary map ∂_k originating from degree `k`.
For example, if degrees=[4, 3, 2] and maps=[M_4, M_3]:
- `C[4]` returns `M_4`
- `C[3]` returns `M_3`
- `C[2]` returns a dynamically sized zero-matrix.
"""
function Base.getindex(C::ChainComplex, k::Int)
    idx = findfirst(==(k), C.degrees)
    
    if idx === nothing
        throw(KeyError("Degree $k not found in the chain complex grading: $(C.degrees)"))
    elseif idx > length(C.maps)
        # Terminal degree: ∂_k is the zero map.
        # We fetch the base ring from the incoming map to build a consistent zero matrix.
        R = base_ring(C.maps[end])
        dim_C_k = nrows(C.maps[end])
        return zero_matrix(R, 0, dim_C_k)
    end
    
    return C.maps[idx]
end

"""
    betti_number(C::ChainComplex, k::Int)

Computes the dimension of the k-th homology group of the chain complex
using native Oscar `rank` calculations.
"""
function betti_number(C::ChainComplex, k::Int)
    idx = findfirst(==(k), C.degrees)
    if idx === nothing
        return 0 
    end
    
    # Outgoing map ∂_k : C_k -> C_{k-1}
    # If k is the terminal degree, this natively fetches the 0-matrix we defined above!
    ∂_out = C[k]
    nullity_out = ncols(∂_out) - rank(∂_out)
    
    # Incoming map ∂_{k+1} : C_{k+1} -> C_k
    if idx == 1
        # It's the highest degree; nothing maps into it.
        rank_in = 0
    else
        deg_in = C.degrees[idx - 1]
        ∂_in = C[deg_in]
        rank_in = rank(∂_in)
    end
    
    return nullity_out - rank_in
end

"""
    set_indices!(C::ChainComplex, new_degrees::Vector{Int})

Overwrites the homological grading of the chain complex in place.
The new grading must have the same length as the original.
"""
function set_indices!(C::ChainComplex, new_degrees::Vector{Int})
    if length(new_degrees) != length(C.degrees)
        throw(ArgumentError("New degrees vector must have length $(length(C.degrees)), got $(length(new_degrees))."))
    end
    
    # Mutate the array in place so we don't break immutability of the struct
    C.degrees .= new_degrees
    return C
end

"""
    shift_indices!(C::ChainComplex, shift::Int)

Shifts the homological grading of the complex by a constant offset.
For example, shifting `[2, 1, 0]` by `+2` results in `[4, 3, 2]`.
"""
function shift_indices!(C::ChainComplex, shift::Int)
    C.degrees .+= shift
    return C
end
