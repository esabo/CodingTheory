# Copyright (c) 2024 - 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
$(TYPEDSIGNATURES)

Return the bipartite adjacency lists `(check_adj, var_adj)` for the LDPC code.
Computes lazily and caches the result.
"""
function node_adjacencies(C::LDPCCode)
    if !haskey(C.cache, :var_adj) || !haskey(C.cache, :check_adj)
        nr, nc = size(C.H)
        var_adj = [Int[] for _ in 1:nc]
        check_adj = [Int[] for _ in 1:nr]
        
        # O(|E|) traversal based on whether the matrix is natively sparse
        if typeof(C.H) <: SMatElem
            for (r, row) in enumerate(C.H)
                for (c, val) in row
                    if !iszero(val)
                        push!(var_adj[c], r)
                        push!(check_adj[r], c)
                    end
                end
            end
        else
            for c in 1:nc
                for r in 1:nr
                    if !iszero(C.H[r, c])
                        push!(var_adj[c], r)
                        push!(check_adj[r], c)
                    end
                end
            end
        end
        
        C.cache[:var_adj] = var_adj
        C.cache[:check_adj] = check_adj
    end
    return C.cache[:check_adj], C.cache[:var_adj]
end

"""
$(TYPEDSIGNATURES)

Return the girth of the Tanner graph of `C`.
"""
function girth(C::LDPCCode)
    haskey(C.cache, :girth) && return C.cache[:girth]
    
    if haskey(C.cache, :short_cycle_dist)
        C.cache[:girth] = minimum(keys(C.cache[:short_cycle_dist]))
        return C.cache[:girth]
    end

    check_adj, var_adj = node_adjacencies(C)
    nr, nc = size(C.H)
    total_nodes = nr + nc
    
    # Thread-safe global minimum
    global_min_girth = Threads.Atomic{Int}(typemax(Int))
    
    # Pre-allocate EXACTLY one workspace per thread to avoid GC thrashing
    n_threads = Threads.nthreads()
    dists   = [fill(-1, total_nodes) for _ in 1:n_threads]
    parents = [fill(-1, total_nodes) for _ in 1:n_threads]
    queues  = [Vector{Int}(undef, total_nodes) for _ in 1:n_threads]
    
    Threads.@threads for root in 1:nc
        # Bipartite graphs cannot have cycles < 4. Stop everything if we found the absolute floor.
        global_min_girth[] == 4 && continue 
        
        tid = Threads.threadid()
        dist = dists[tid]
        parent = parents[tid]
        queue = queues[tid]
        
        # Reset only the workspace for this specific thread
        fill!(dist, -1)
        fill!(parent, -1)
        
        head = 1
        tail = 2
        queue[1] = root
        dist[root] = 0
        
        while head < tail
            u = queue[head]
            head += 1
            
            # Pruning: Stop if this tree's depth exceeds the globally found shortest cycle
            if dist[u] * 2 >= global_min_girth[]
                break
            end
            
            is_var = u <= nc
            neighbors = is_var ? check_adj[u] : var_adj[u - nc]
            
            for n_idx in neighbors
                v = is_var ? n_idx + nc : n_idx
                
                if v != parent[u]
                    if dist[v] == -1
                        dist[v] = dist[u] + 1
                        parent[v] = u
                        queue[tail] = v
                        tail += 1
                    else
                        # Cycle found!
                        cycle_len = dist[u] + dist[v] + 1
                        if cycle_len < global_min_girth[]
                            # Safely update the global minimum across all threads
                            Threads.atomic_min!(global_min_girth, cycle_len)
                        end
                    end
                end
            end
        end
    end
    
    final_girth = global_min_girth[]
    C.cache[:girth] = final_girth == typemax(Int) ? -1 : final_girth
    return C.cache[:girth]
end

"""
$(TYPEDSIGNATURES)

Return the local girth (computation tree depth) for the variable node `v` in the 
Tanner graph of `C`. 

# Notes
* This metric intrinsically captures the "lollipop graph" penalty. It returns the length 
  of the shortest message-passing loop that originates from and returns to `v`, calculating 
  the cycle length plus twice the length of the stem.
* Returns `-1` if the node is part of a tree structure (no cycles reachable).
"""
function local_girth(C::LDPCCode, v::Int)
    check_adj, var_adj = node_adjacencies(C)
    nr, nc = size(C.H)
    total_nodes = nr + nc
    
    1 <= v <= nc || throw(BoundsError("Variable node index must be between 1 and $nc"))
    
    # Pre-allocate zero-allocation BFS workspace
    dist = fill(-1, total_nodes)
    parent = fill(-1, total_nodes)
    queue = Vector{Int}(undef, total_nodes)
    
    head = 1
    tail = 2
    queue[1] = v
    dist[v] = 0
    
    while head < tail
        u = queue[head]
        head += 1
        
        is_var = u <= nc
        neighbors = is_var ? check_adj[u] : var_adj[u - nc]
        
        for n_idx in neighbors
            w = is_var ? n_idx + nc : n_idx
            
            if w != parent[u]
                if dist[w] == -1
                    dist[w] = dist[u] + 1
                    parent[w] = u
                    queue[tail] = w
                    tail += 1
                else
                    # The very first collision in a BFS guarantees the shortest 
                    # possible topological loop back to the root vertex v.
                    return dist[u] + dist[w] + 1
                end
            end
        end
    end
    
    return -1 # Node v is part of a pure tree
end

"""
$(TYPEDSIGNATURES)

Return the local girth for a specified list of variable nodes `vs`.
"""
local_girth(C::LDPCCode, vs::Vector{Int}) = [local_girth(C, v) for v in vs]

"""
$(TYPEDSIGNATURES)

Return the local girth for every variable node in the LDPC code `C`.
"""
local_girth(C::LDPCCode) = local_girth(C, collect(1:C.n))

"""
    computation_graph(C::AbstractLDPCCode, lvl::Int, v::Int, v_type::Symbol = :v)

Return a figure representing the expansion of the Tanner graph of `C` to level `lvl`
for node `v`. If `v_type` is `:v`, `v` is interpreted as a variable node; otherwise,
`v_type` is `:c` and `v` is interpreted as a check node.

# Note
- Run `using Makie` to activate this extension.
"""
function computation_graph end

"""
$(TYPEDSIGNATURES)

Attempt to structurally remove cycles of length strictly less than `target_girth` 
from the LDPC code `C` using the BFS Socket-Swapping algorithm.

# Notes
* `target_girth` must be an even integer >= 4.
* This algorithm exactly preserves the degree distributions (λ and ρ) of the original code.
* Returns a strictly new `LDPCCode` object. If the target girth cannot be reached 
  within `max_iters`, it returns the best-effort matrix achieved so far.
"""
function remove_cycles(C::LDPCCode, target_girth::Int; max_iters::Int=2000)
    iseven(target_girth) && target_girth >= 4 || throw(ArgumentError("Target girth must be an even integer >= 4"))
    
    # If the matrix is empty or a tree, just return a copy
    girth(C) == -1 && return LDPCCode(C.H)
    
    nr, nc = size(C.H)
    F = base_ring(C.H)
    
    # 1. Build mutable dictionary adjacency lists to track edge field values
    var_adj = [Dict{Int, typeof(F(1))}() for _ in 1:nc]
    check_adj = [Dict{Int, typeof(F(1))}() for _ in 1:nr]
    
    if typeof(C.H) <: SMatElem
        for (r, row) in enumerate(C.H)
            for (c, val) in row
                if !iszero(val)
                    var_adj[c][r] = val
                    check_adj[r][c] = val
                end
            end
        end
    else
        for c in 1:nc
            for r in 1:nr
                val = C.H[r, c]
                if !iszero(val)
                    var_adj[c][r] = val
                    check_adj[r][c] = val
                end
            end
        end
    end

    total_nodes = nr + nc
    dist = fill(-1, total_nodes)
    parent = fill(-1, total_nodes)
    queue = Vector{Int}(undef, total_nodes)
    
    # Internal fast BFS to find exactly one edge involved in a short cycle
    function _get_bad_edge()
        for root in 1:nc
            fill!(dist, -1)
            fill!(parent, -1)
            head = 1
            tail = 2
            queue[1] = root
            dist[root] = 0
            
            while head < tail
                u = queue[head]
                head += 1
                
                # Prune if we've searched deep enough
                if dist[u] * 2 >= target_girth
                    break
                end
                
                is_var = u <= nc
                neighbors = is_var ? keys(var_adj[u]) : keys(check_adj[u - nc])
                
                for n_idx in neighbors
                    v = is_var ? n_idx + nc : n_idx
                    
                    if v != parent[u]
                        if dist[v] == -1
                            dist[v] = dist[u] + 1
                            parent[v] = u
                            queue[tail] = v
                            tail += 1
                        else
                            # Collision! A cycle is found.
                            cycle_len = dist[u] + dist[v] + 1
                            if cycle_len < target_girth
                                # Return the edge connecting u to its parent
                                if is_var
                                    return (u, parent[u] - nc)
                                else
                                    return (parent[u], u - nc)
                                end
                            end
                        end
                    end
                end
            end
        end
        return nothing # Target girth successfully achieved!
    end

    # 2. Main Socket-Swapping Loop
    iters_used = 0
    for iter in 1:max_iters
        bad_edge = _get_bad_edge()
        isnothing(bad_edge) && break
        
        v1, c1 = bad_edge
        
        # Pick a random edge (v2, c2) to swap with
        swap_valid = false
        attempts = 0
        v2, c2 = 0, 0
        
        while !swap_valid && attempts < 50
            v2 = rand(1:nc)
            # Skip empty columns
            isempty(var_adj[v2]) && continue 
            
            c2 = rand(collect(keys(var_adj[v2])))
            
            # Ensure the swap won't create multi-edges or self-loops
            if v1 != v2 && c1 != c2 && !haskey(var_adj[v1], c2) && !haskey(var_adj[v2], c1)
                swap_valid = true
            end
            attempts += 1
        end
        
        if swap_valid
            val1 = var_adj[v1][c1]
            val2 = var_adj[v2][c2]
            
            # Disconnect old edges
            delete!(var_adj[v1], c1)
            delete!(check_adj[c1], v1)
            delete!(var_adj[v2], c2)
            delete!(check_adj[c2], v2)
            
            # Reconnect new swapped edges
            var_adj[v1][c2] = val1
            check_adj[c2][v1] = val1
            var_adj[v2][c1] = val2
            check_adj[c1][v2] = val2
        end
        iters_used += 1
    end
    
    if iters_used == max_iters
        @warn "remove_cycles: Reached max_iters ($max_iters) without achieving target girth $target_girth. Returning best-effort matrix."
    end
    
    # 3. Reassemble the modified sparse matrix
    I_idx = Int[]
    J_idx = Int[]
    V_val = typeof(F(1))[]
    
    for c in 1:nc
        for (r, val) in var_adj[c]
            push!(I_idx, r)
            push!(J_idx, c)
            push!(V_val, val)
        end
    end
    
    new_H = sparse_matrix(F, nr, nc, I_idx, J_idx, V_val)
    
    # Returning via LDPCCode instantly primes the new cache and distributions
    return LDPCCode(new_H)
end

"""
$(TYPEDSIGNATURES)

Attempt to structurally remove cycles of length up to `n_max` from the LDPC code `C`.
Returns a strictly **new** `LDPCCode` object.
"""
function remove_cycles(C::LDPCCode, n_max::Int)
    new_H = remove_cycles(parity_check_matrix(C), n_max)
    return LDPCCode(new_H)
end

#############################
       # simple cycles
#############################

function _circuit_recursive!(v1::Int, v2::Int, blocked::Vector{Bool}, B::Vector{Vector{Int}}, 
                             stack::Vector{Int}, cycles::Vector{Vector{Int}}, 
                             unique_cycles::Set{Vector{Int}}, max_len::Int, 
                             check_adj::Vector{Vector{Int}}, var_adj::Vector{Vector{Int}}, nc::Int)
    
    flag = false
    push!(stack, v2)
    blocked[v2] = true

    is_var = v2 <= nc
    neighbors = is_var ? check_adj[v2] : var_adj[v2 - nc]

    for n_idx in neighbors
        w = is_var ? n_idx + nc : n_idx
        (w < v1) && continue
        
        # Don't trivially backtrack
        length(stack) > 1 && w == stack[end-1] && continue

        if w == v1
            if length(stack) >= 4 # Bipartite graphs only have cycles >= 4
                cycle = copy(stack)
                sorted_cycle = sort(cycle)
                
                # O(1) hash lookup replaces the old O(N^2) loop
                if !(sorted_cycle in unique_cycles)
                    push!(unique_cycles, sorted_cycle)
                    push!(cycles, cycle)
                end
            end
            flag = true
        elseif !blocked[w] && length(stack) < max_len
            flag |= _circuit_recursive!(v1, w, blocked, B, stack, cycles, unique_cycles, max_len, check_adj, var_adj, nc)
        end
    end

    if flag
        _unblock!(v2, blocked, B)
    else
        for n_idx in neighbors
            w = is_var ? n_idx + nc : n_idx
            (w < v1) && continue
            if !(v2 in B[w])
                push!(B[w], v2)
            end
        end
    end
    
    pop!(stack)
    return flag
end

function _unblock!(v::Int, blocked::Vector{Bool}, B::Vector{Vector{Int}})
    blocked[v] = false
    Bv = B[v]
    while !isempty(Bv)
        w = pop!(Bv)
        if blocked[w]
            _unblock!(w, blocked, B)
        end
    end
end

#############################
       # simple cycles
#############################

"""
$(TYPEDSIGNATURES)

Return the unique simple cycles up to length `len` of the Tanner graph of `L`.
An empty `Vector{Vector{Int}}` is returned when there is no cycles.

# Note
- Simple cycles do not contain the same vertex twice.
- Cycles are returned as a vector of vertex indices, where the vertices are ordered left-to-right by
  columns of `parity_check_matrix(L)` then top-to-bottom by rows.
"""
function enumerate_simple_cycles(C::LDPCCode; len::Int = 16)
    len > 0 || throw(DomainError("Cycle length must be positive"))
    
    # Check the tree edge-case using the Batch 1 girth function
    g = girth(C)
    g == -1 && return Vector{Vector{Int}}()

    if !haskey(C.cache, :max_cyc_len)
        C.cache[:max_cyc_len] = 0
        C.cache[:simple_cycles] = Vector{Vector{Int}}()
    end

    # Return instantly if already computed
    if len <= C.cache[:max_cyc_len]
        return filter(x -> length(x) <= len, C.cache[:simple_cycles])
    end

    check_adj, var_adj = node_adjacencies(C)
    nr, nc = size(C.H)
    total_nodes = nr + nc

    cycles = Vector{Vector{Int}}()
    unique_cycles = Set{Vector{Int}}()
    
    # Start the search from previously found cycles to save time
    if !isempty(C.cache[:simple_cycles])
        cycles = copy(C.cache[:simple_cycles])
        for c in cycles
            push!(unique_cycles, sort(c))
        end
    end

    Threads.@threads for i in 1:nc
        blocked = fill(false, total_nodes)
        B = [Int[] for _ in 1:total_nodes]
        stack = Int[]
        _circuit_recursive!(i, i, blocked, B, stack, cycles, unique_cycles, len, check_adj, var_adj, nc)
    end

    C.cache[:max_cyc_len] = len
    C.cache[:simple_cycles] = cycles
    return cycles
end

"""
    simple_cycle_length_distribution(L::AbstractLDPCCode; len::Int = 16)

Return a dictionary of (length, count) pairs for the unique simple cycles up to length `len` of the
Tanner graph of `L`. If `len` is `-1`, then all simple cycles will be enumerated. An empty
dictionary is returned when there are no cycles.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
function simple_cycle_length_distribution(C::LDPCCode; len::Int = 16)
    cycles = enumerate_simple_cycles(C; len = len)
    isempty(cycles) && return Dict{Int, Int}()
    
    dist = Dict{Int, Int}()
    for c in cycles
        l = length(c)
        dist[l] = get(dist, l, 0) + 1
    end
    
    # Mutual resolution: update girth if we found something shorter
    min_len = minimum(keys(dist))
    if !haskey(C.cache, :girth) || min_len < C.cache[:girth]
        C.cache[:girth] = min_len
    end
    
    return dist
end

"""
    simple_cycle_length_distribution_plot(L::AbstractLDPCCode; len::Int = 16)

Return a bar graph and dictionary of (length, count) pairs for the unique simple cycles up to
length `len` of the Tanner graph of `L`. If `len` is `-1`, then all simple cycles will be
enumerated. An empty figure and dictionary are returned when there are no cycles.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
- Run `using Makie` to activate this extension.
"""
function simple_cycle_length_distribution_plot end

"""
    average_simple_cycle_length(L::AbstractLDPCCode; len::Int = 16)

Return the average cycle length of unique simple cycles up to length `len` of the Tanner graph of
`L`. If `len` is `-1`, then all simple cycles will be enumerated.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
function average_simple_cycle_length(C::LDPCCode; len::Int = 16)
    dist = simple_cycle_length_distribution(C; len = len)
    isempty(dist) && return NaN
    
    total_length = sum(k * v for (k, v) in dist)
    total_count = sum(values(dist))
    return total_length / total_count
end

"""
    median_simple_cycle_length(L::AbstractLDPCCode; len::Int = 16)

Return the median cycle length of unique simple cycles up to length `len` of the Tanner graph of
`L`. If `len` is `-1`, then all simple cycles will be enumerated.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
function median_simple_cycle_length(C::LDPCCode; len::Int = 16)
    dist = simple_cycle_length_distribution(C; len = len)
    isempty(dist) && return NaN
    
    counts = Int[]
    for k in sort(collect(keys(dist)))
        append!(counts, fill(k, dist[k]))
    end
    return median(counts)
end

"""
    mode_simple_cycle_length(L::AbstractLDPCCode; len::Int = 16)

Return the most common cycle length of unique simple cycles up to length `len` of the Tanner graph
of `L`. If `len` is `-1`, then all simple cycles will be enumerated.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
function mode_simple_cycle_length(C::LDPCCode; len::Int = 16)
    dist = simple_cycle_length_distribution(C; len = len)
    isempty(dist) && return NaN
    
    max_count = maximum(values(dist))
    for (k, v) in dist
        v == max_count && return k
    end
end

"""
    count_simple_cycles(L::AbstractLDPCCode; len::Int = 16)

Return the total number of unique simple cycles up to length `len` of the Tanner graph
of `L`. If `len` is `-1`, then all simple cycles will be enumerated.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
count_simple_cycles(C::LDPCCode; len::Int=16) = sum(values(simple_cycle_length_distribution(C; len=len)), init=0)

"""
    simple_cycle_distribution_by_variable_node(L::AbstractLDPCCode; len::Int = 16)

Return a dictionary of (node, count) pairs for the unique simple cycles up to length `len` of the
Tanner graph of `L`. If `len` is `-1`, then all simple cycles will be enumerated. An empty
dictionary is returned when there are no cycles.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
function simple_cycle_distribution_by_variable_node(C::LDPCCode; len::Int = 16)
    cycles = enumerate_simple_cycles(C; len = len)
    isempty(cycles) && return Dict{Int, Int}()
    
    dist = Dict{Int, Int}()
    nc = size(C.H, 2)
    for c in cycles
        for v in c
            if v <= nc # Only count variable nodes
                dist[v] = get(dist, v, 0) + 1
            end
        end
    end
    return dist
end

"""
    simple_cycle_distribution_by_variable_node_plot(L::AbstractLDPCCode; len::Int = 16)

Return bar graph and a dictionary of (node, count) pairs for the unique simple cycles up to length
`len` of the Tanner graph of `L`. If `len` is `-1`, then all simple cycles will be enumerated. An
empty figure and dictionary are returned when there are no cycles.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
    already cached.
- Run `using Makie` to activate this extension.
"""
function simple_cycle_distribution_by_variable_node_plot end

#############################
        # short cycles
#############################

"""
$(TYPEDSIGNATURES)

Return the unique short cycles (length between `g` and `2g - 2`) of the Tanner graph of `C`.
"""
function enumerate_short_cycles(C::LDPCCode)
    g = girth(C)
    g == -1 && return Vector{Vector{Int}}()
    return enumerate_simple_cycles(C; len = 2*g - 2)
end

"""
    short_cycle_length_distribution(L::AbstractLDPCCode; len::Int = 16)

Return a dictionary of (length, count) pairs for the unique short cycles up to length `len` of the
Tanner graph of `L`. If `len` is `-1`, then all short cycles will be enumerated. An empty
dictionary is returned when there are no cycles.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
function short_cycle_length_distribution(C::LDPCCode)
    g = girth(C)
    g == -1 && return Dict{Int, Int}()
    return simple_cycle_length_distribution(C; len = 2*g - 2)
end

"""
    short_cycle_length_distribution_plot(L::AbstractLDPCCode; len::Int = 16)

Return a bar graph and dictionary of (length, count) pairs for the unique short cycles up to
length `len` of the Tanner graph of `L`. If `len` is `-1`, then all short cycles will be
enumerated. An empty figure and dictionary are returned when there are no cycles.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
- Run `using Makie` to activate this extension.
"""
function short_cycle_length_distribution_plot end

"""
    average_short_cycle_length(L::AbstractLDPCCode; len::Int = 16)

Return the average cycle length of unique short cycles up to length `len` of the Tanner graph of
`L`. If `len` is `-1`, then all short cycles will be enumerated.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
function average_short_cycle_length(C::LDPCCode)
    dist = short_cycle_length_distribution(C)
    isempty(dist) && return NaN
    
    total_len = sum(k * v for (k, v) in dist)
    total_count = sum(values(dist))
    return total_len / total_count
end

"""
    median_short_cycle_length(L::AbstractLDPCCode; len::Int = 16)

Return the median cycle length of unique short cycles up to length `len` of the Tanner graph of
`L`. If `len` is `-1`, then all short cycles will be enumerated.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
function median_short_cycle_length(C::LDPCCode)
    dist = short_cycle_length_distribution(C)
    isempty(dist) && return NaN
    
    # Expand into a sorted array for median computation
    counts = Int[]
    for k in sort(collect(keys(dist)))
        append!(counts, fill(k, dist[k]))
    end
    return median(counts)
end

"""
    mode_short_cycle_length(L::AbstractLDPCCode; len::Int = 16)

Return the most common cycle length of unique short cycles up to length `len` of the Tanner graph
of `L`. If `len` is `-1`, then all short cycles will be enumerated.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
function mode_short_cycle_length(C::LDPCCode)
    dist = short_cycle_length_distribution(C)
    isempty(dist) && return NaN
    
    # Mode is just the key with the maximum value
    max_count = maximum(values(dist))
    for (k, v) in dist
        v == max_count && return k
    end
end

 """
    count_short_cycles(L::AbstractLDPCCode; len::Int = 16)

Return the total number of unique short cycles up to length `len` of the Tanner graph
of `L`. If `len` is `-1`, then all short cycles will be enumerated.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
count_short_cycles(C::LDPCCode) = sum(values(short_cycle_length_distribution(C)), init=0)

"""
    short_cycle_distribution_by_variable_node(L::AbstractLDPCCode; len::Int = 16)

Return a dictionary of (node, count) pairs for the unique short cycles up to length `len` of the
Tanner graph of `L`. If `len` is `-1`, then all short cycles will be enumerated. An empty
dictionary is returned when there are no cycles.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
function short_cycle_distribution_by_variable_node(C::LDPCCode)
    cycles = enumerate_short_cycles(C)
    isempty(cycles) && return Dict{Int, Int}()
    
    dist = Dict{Int, Int}()
    nc = size(C.H, 2)
    for c in cycles
        for v in c
            if v <= nc # Only count variable nodes
                dist[v] = get(dist, v, 0) + 1
            end
        end
    end
    return dist
end

"""
    short_cycle_distribution_by_variable_node_plot(L::AbstractLDPCCode; len::Int = 16)

Return bar graph and a dictionary of (node, count) pairs for the unique short cycles up to length
`len` of the Tanner graph of `L`. If `len` is `-1`, then all short cycles will be enumerated. An
empty figure and dictionary are returned when there are no cycles.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
    already cached.
- Run `using Makie` to activate this extension.
"""
function short_cycle_distribution_by_variable_node_plot end

#############################
           # ACE
#############################

"""
$(TYPEDSIGNATURES)

Internal parallel engine to compute the shortest cycle lengths and ACE distributions 
for every variable node in the graph simultaneously.
"""
function _compute_ACE_distributions(C::LDPCCode)
    if haskey(C.cache, :ACE_dists) && haskey(C.cache, :shortest_cycle_lens)
        return C.cache[:shortest_cycle_lens], C.cache[:ACE_dists]
    end

    check_adj, var_adj = node_adjacencies(C)
    nr, nc = size(C.H)
    total_nodes = nr + nc

    # Output Arrays
    shortest_lens = fill(typemax(Int), nc)
    ace_dists = [Int[] for _ in 1:nc]

    # Pre-allocate exactly one zero-allocation workspace per thread
    n_threads = Threads.nthreads()
    dists   = [fill(-1, total_nodes) for _ in 1:n_threads]
    parents = [fill(-1, total_nodes) for _ in 1:n_threads]
    ace_wts = [fill(0, total_nodes) for _ in 1:n_threads]
    queues  = [Vector{Int}(undef, total_nodes) for _ in 1:n_threads]

    Threads.@threads for root in 1:nc
        tid = Threads.threadid()
        dist = dists[tid]
        parent = parents[tid]
        ace_wt = ace_wts[tid]
        queue = queues[tid]

        fill!(dist, -1)
        fill!(parent, -1)

        head = 1
        tail = 2
        queue[1] = root
        dist[root] = 0
        
        # Local ACE of the root node
        root_deg = length(var_adj[root])
        ace_wt[root] = root_deg - 2

        min_len = typemax(Int)
        local_aces = Int[]

        while head < tail
            u = queue[head]
            head += 1

            # Prune search instantly once we exceed the shortest cycle length found for THIS node
            if dist[u] * 2 > min_len
                break
            end

            is_var = u <= nc
            neighbors = is_var ? var_adj[u] : check_adj[u - nc]

            for n_idx in neighbors
                v = is_var ? n_idx + nc : n_idx

                if v != parent[u]
                    if dist[v] == -1
                        dist[v] = dist[u] + 1
                        parent[v] = u
                        
                        # Accumulate ACE mathematically
                        if v <= nc
                            ace_wt[v] = ace_wt[u] + length(var_adj[v]) - 2
                        else
                            ace_wt[v] = ace_wt[u]
                        end
                        
                        queue[tail] = v
                        tail += 1
                    else
                        # Collision! A cycle is closed.
                        cycle_len = dist[u] + dist[v] + 1
                        if cycle_len <= min_len
                            min_len = cycle_len
                            # The exact ACE of the cycle avoids double-counting the root
                            c_ace = ace_wt[u] + ace_wt[v] - (root_deg - 2)
                            push!(local_aces, c_ace)
                        end
                    end
                end
            end
        end
        
        shortest_lens[root] = min_len == typemax(Int) ? -1 : min_len
        # The BFS explores symmetrically, so collisions are detected twice. `unique` instantly deduplicates.
        ace_dists[root] = unique(local_aces) 
    end

    C.cache[:shortest_cycle_lens] = shortest_lens
    C.cache[:ACE_dists] = ace_dists
    
    # Mutual Resolution: Update global girth if we found a new minimum
    true_girth = minimum(filter(x -> x != -1, shortest_lens))
    C.cache[:girth] = true_girth == typemax(Int) ? -1 : true_girth
    
    return shortest_lens, ace_dists
end

"""
$(TYPEDSIGNATURES)

Return the ACE spectrum of the Tanner graph of `C`.
Returns a `Dict{Int, Dict{Int, Int}}` mapping `Cycle Length -> (Minimum ACE -> Count)`.
"""
function ACE_spectrum(C::LDPCCode)
    lens, ace_dists = _compute_ACE_distributions(C)
    
    # Handle the tree edge-case gracefully
    girth(C) == -1 && return Dict{Int, Dict{Int, Int}}()
    
    spectrum = Dict{Int, Dict{Int, Int}}()
    
    for root in 1:C.n
        l = lens[root]
        l == -1 && continue
        
        if !haskey(spectrum, l)
            spectrum[l] = Dict{Int, Int}()
        end
        
        # The true spectrum groups nodes by the minimum ACE among their shortest cycles
        min_ace = minimum(ace_dists[root])
        spectrum[l][min_ace] = get(spectrum[l], min_ace, 0) + 1
    end
    
    return spectrum
end

"""
    ACE_spectrum_plot(C::AbstractLDPCCode)

Return an interactive figure and data for the ACE spectrum of the Tanner graph of `C`.

# Note
- Run `using Makie` to activate this extension.
"""
function ACE_spectrum_plot end


"""
$(TYPEDSIGNATURES)

Return the exact ACE arrays for the shortest cycles of the given variable node(s).
"""
function ACE_distribution(C::LDPCCode, vs::Vector{Int})
    _, ace_dists = _compute_ACE_distributions(C)
    return [ace_dists[v] for v in vs]
end
ACE_distribution(C::LDPCCode, v::Int) = ACE_distribution(C, [v])[1]
ACE_distribution(C::LDPCCode) = ACE_distribution(C, collect(1:C.n))

"""
$(TYPEDSIGNATURES)

Return the average ACE of the vertex `v` or vertices `vs` of the Tanner graph of `C`. If no
vertices are given, all vertices are computed (individually) by default.
"""
function average_ACE_distribution(C::LDPCCode, vs::Vector{Int})
    _, ace_dists = _compute_ACE_distributions(C)
    return [isempty(ace_dists[v]) ? NaN : mean(ace_dists[v]) for v in vs]
end
average_ACE_distribution(C::LDPCCode, v::Int) = average_ACE_distribution(C, [v])[1]
average_ACE_distribution(C::LDPCCode) = average_ACE_distribution(C, collect(1:C.n))

"""
$(TYPEDSIGNATURES)

Return the median ACE of the vertex `v` or vertices `vs` of the Tanner graph of `C`. If no
vertices are given, all vertices are computed (individually) by default.
"""
function median_ACE_distribution(C::LDPCCode, vs::Vector{Int})
    _, ace_dists = _compute_ACE_distributions(C)
    return [isempty(ace_dists[v]) ? NaN : median(ace_dists[v]) for v in vs]
end
median_ACE_distribution(C::LDPCCode, v::Int) = median_ACE_distribution(C, [v])[1]
median_ACE_distribution(C::LDPCCode) = median_ACE_distribution(C, collect(1:C.n))

"""
$(TYPEDSIGNATURES)

Return the mode ACE of the vertex `v` or vertices `vs` of the Tanner graph of `C`. If no
vertices are given, all vertices are computed (individually) by default.

# Note
- In case of ties, the smallest tied value is returned.
"""
function mode_ACE_distribution(C::LDPCCode, vs::Vector{Int})
    _, ace_dists = _compute_ACE_distributions(C)
    # Returns NaN if the node has no cycles, otherwise finds the most frequent ACE
    return [isempty(ace_dists[v]) ? NaN : mode(ace_dists[v]) for v in vs]
end
mode_ACE_distribution(C::LDPCCode, v::Int) = mode_ACE_distribution(C, [v])[1]
mode_ACE_distribution(C::LDPCCode) = mode_ACE_distribution(C, collect(1:C.n))
