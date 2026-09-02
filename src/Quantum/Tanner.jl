"""
    Tanner_graph_X(S::AbstractSubsystemCode)

Return the bipartite `SimpleGraph` representing only the X-checks and qubits.
Results are cached in `S.cache[:Tanner_graph_X]`.
"""
function Tanner_graph_X(S::AbstractSubsystemCode)
    CSSTrait(typeof(S)) == IsCSS() || throw(ArgumentError("Tanner_graph_X is only defined for CSS codes."))
    haskey(S.cache, :Tanner_graph_X) && return S.cache[:Tanner_graph_X]
    
    H_X = X_stabilizers(S)
    nr, nc = size(H_X)
    total_V = nr + nc
    
    G = Grphs.SimpleGraph(total_V)
    for r in 1:nr
        for c in 1:nc
            if !iszero(H_X[r, c])
                Grphs.add_edge!(G, c, r + nc)
            end
        end
    end
    
    result = (G, collect(1:nc), collect(nc + 1 : total_V))
    S.cache[:Tanner_graph_X] = result
    return result
end

"""
    Tanner_graph_Z(S::AbstractSubsystemCode)

Return the bipartite `SimpleGraph` representing only the Z-checks and qubits.
Results are cached in `S.cache[:Tanner_graph_Z]`.
"""
function Tanner_graph_Z(S::AbstractSubsystemCode)
    CSSTrait(typeof(S)) == IsCSS() || throw(ArgumentError("Tanner_graph_Z is only defined for CSS codes."))
    haskey(S.cache, :Tanner_graph_Z) && return S.cache[:Tanner_graph_Z]
    
    H_Z = Z_stabilizers(S)
    nr, nc = size(H_Z)
    total_V = nr + nc
    
    G = Grphs.SimpleGraph(total_V)
    for r in 1:nr
        for c in 1:nc
            if !iszero(H_Z[r, c])
                Grphs.add_edge!(G, c, r + nc)
            end
        end
    end
    
    result = (G, collect(1:nc), collect(nc + 1 : total_V))
    S.cache[:Tanner_graph_Z] = result
    return result
end

"""
    Tanner_graph(S::AbstractSubsystemCode)

Return the `SimpleGraph` object representing the Tanner graph of the code `S`.
Automatically generates a tripartite graph (4-tuple) for CSS codes and a bipartite graph (3-tuple) for non-CSS codes.
Results are cached in `S.cache[:Tanner_graph]`.
"""
function Tanner_graph(S::AbstractSubsystemCode)
    haskey(S.cache, :Tanner_graph) && return S.cache[:Tanner_graph]
    
    if CSSTrait(typeof(S)) == IsCSS()
        H_X = X_stabilizers(S)
        H_Z = Z_stabilizers(S)
        
        m_X, n = size(H_X)
        m_Z, _ = size(H_Z)
        total_V = n + m_X + m_Z
        
        G = Grphs.SimpleGraph(total_V)
        
        for r in 1:m_X
            for c in 1:n
                !iszero(H_X[r, c]) && Grphs.add_edge!(G, c, n + r)
            end
        end
        
        for r in 1:m_Z
            for c in 1:n
                !iszero(H_Z[r, c]) && Grphs.add_edge!(G, c, n + m_X + r)
            end
        end
        
        result = (G, collect(1:n), collect(n + 1 : n + m_X), collect(n + m_X + 1 : total_V))
    else
        H = stabilizers(S)
        nr, nc = size(H)
        n = div(nc, 2)
        total_V = nr + n
        
        G = Grphs.SimpleGraph(total_V)
        
        for r in 1:nr
            for c in 1:n
                (!iszero(H[r, c]) || !iszero(H[r, c + n])) && Grphs.add_edge!(G, c, r + n)
            end
        end
        
        result = (G, collect(1:n), collect(n + 1 : total_V))
    end
    
    S.cache[:Tanner_graph] = result
    return result
end
