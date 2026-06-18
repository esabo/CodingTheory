# # Copyright (c) 2026 Eric Sabo
# # All rights reserved.
# #
# # This source code is licensed under the BSD-style license found in the
# # LICENSE file in the root directory of this source tree.

"""
Expanded workspace to support ultra-fast BFS and ACE tracking without allocations.
"""
struct PEGWorkspace
    m::Int
    n::Int
    
    v_adj::Vector{Vector{Int}}
    c_adj::Vector{Vector{Int}}
    c_degrees::Vector{Int}
    v_degrees::Vector{Int}
    
    # BFS Trackers
    depths::Vector{Int}
    visited_c::Vector{Bool}
    visited_v::Vector{Bool}
    
    # ACE Trackers: The accumulated Extrinsic Message Degree along a path
    ace_path::Vector{Int}
    
    # Preallocated queue for BFS (prevents push!/pop! allocations)
    q::Vector{Int}
    
    function PEGWorkspace(n::Int, m::Int)
        new(
            m, n,
            [Int[] for _ in 1:n],
            [Int[] for _ in 1:m],
            zeros(Int, m),
            zeros(Int, n),
            fill(-1, m),
            fill(false, m),
            fill(false, n),
            zeros(Int, m),
            zeros(Int, max(n, m))
        )
    end
end

function _add_edge!(W::PEGWorkspace, v::Int, c::Int)
    push!(W.v_adj[v], c)
    push!(W.c_adj[c], v)
    W.c_degrees[c] += 1
    W.v_degrees[v] += 1
end

"""
$(TYPEDSIGNATURES)

Expand the subgraph from `root_v` using BFS and find the optimal check node to connect to,
utilizing both the PEG (max depth/min degree) and ACE (cycle EMD) criteria.
"""
function _expand_and_find_optimal_check!(W::PEGWorkspace, root_v::Int)
    # 1. Reset the BFS trackers
    fill!(W.depths, -1)
    fill!(W.visited_c, false)
    fill!(W.visited_v, false)
    fill!(W.ace_path, 0)
    
    W.visited_v[root_v] = true
    
    # 2. Initialize BFS from the root's current check nodes
    q_head = 1
    q_tail = 1
    
    for c in W.v_adj[root_v]
        W.q[q_tail] = c
        q_tail += 1
        W.visited_c[c] = true
        W.depths[c] = 1
        W.ace_path[c] = 0 # No ACE contribution from the root itself
    end
    
    max_depth = 1
    
    # 3. BFS Traversal
    while q_head < q_tail
        curr_c = W.q[q_head]
        q_head += 1
        
        curr_depth = W.depths[curr_c]
        max_depth = max(max_depth, curr_depth)
        
        # Traverse up to connected variable nodes
        for v in W.c_adj[curr_c]
            if !W.visited_v[v]
                W.visited_v[v] = true
                
                # Calculate the ACE contribution of this variable node
                # Extrinsic connections = (current edges) - 2 (the path in and out)
                v_contribution = max(0, W.v_degrees[v] - 2)
                current_ace = W.ace_path[curr_c] + v_contribution
                
                # Traverse back down to the next layer of check nodes
                for next_c in W.v_adj[v]
                    if !W.visited_c[next_c]
                        W.visited_c[next_c] = true
                        W.depths[next_c] = curr_depth + 1
                        W.ace_path[next_c] = current_ace
                        
                        W.q[q_tail] = next_c
                        q_tail += 1
                    end
                end
            end
        end
    end
    
    # 4. Filter and select the optimal check node
    # We want a check node that was either unreachable (depth == -1) 
    # or is at the maximum possible depth.
    candidate_check = -1
    best_degree = typemax(Int)
    best_ace = -1
    
    # First, check for unreachable check nodes (infinite girth)
    for c in 1:W.m
        if W.depths[c] == -1
            if W.c_degrees[c] < best_degree
                best_degree = W.c_degrees[c]
                candidate_check = c
            end
        end
    end
    
    # If we found an unreachable check node, return it immediately.
    if candidate_check != -1
        return candidate_check
    end
    
    # Otherwise, evaluate nodes at the maximum depth using PEG and ACE
    for c in 1:W.m
        if W.depths[c] == max_depth
            c_deg = W.c_degrees[c]
            c_ace = W.ace_path[c]
            
            # The PEG+ACE Selection Logic:
            # 1. Prefer higher ACE recovery (if creating a cycle, make sure it has high extrinsic connections)
            # 2. Break ties with lower check node degree (standard PEG)
            if c_ace > best_ace || (c_ace == best_ace && c_deg < best_degree)
                best_ace = c_ace
                best_degree = c_deg
                candidate_check = c
            end
        end
    end
    
    return candidate_check
end

"""
$(TYPEDSIGNATURES)

Generate a parity-check matrix using the Progressive Edge-Growth (PEG) algorithm
enhanced with Approximate Cycle EMD (ACE) maximization.

# Arguments
* `n::Int`: Number of variable nodes (columns).
* `m::Int`: Number of check nodes (rows).
* `deg_v::Vector{Int}`: Target degree for each variable node.
"""
function progressive_edge_growth(n::Int, m::Int, deg_v::Vector{Int})
    length(deg_v) == n || throw(ArgumentError("Length of deg_v must match n"))
    
    workspace = PEGWorkspace(n, m)
    
    # ---------------------------------------------------------
    # The "Page 10 Optimization": Fisher-Yates 0-degree pool
    # ---------------------------------------------------------
    # We shuffle the check nodes so that early random selections 
    # distribute the edges uniformly without needing tree traversal.
    available_c = collect(1:m)
    shuffle!(available_c)
    
    for j in 1:n
        for k in 1:deg_v[j]
            if k == 1
                # The first edge can never create a cycle. Connect to an isolated check node if possible.
                if !isempty(available_c)
                    c_opt = pop!(available_c)
                else
                    # Fallback: find the absolute lowest degree check node globally
                    c_opt = argmin(workspace.c_degrees)
                end
                _add_edge!(workspace, j, c_opt)
            else
                # Subsequent edges risk creating cycles. Run the fast PEG+ACE BFS engine.
                c_opt = _expand_and_find_optimal_check!(workspace, j)
                _add_edge!(workspace, j, c_opt)
                
                # Housekeeping: If the BFS engine happened to select a degree-0 check node
                # that was still in our shuffled pool, remove it so we don't double-dip.
                if workspace.c_degrees[c_opt] == 1
                    idx = findfirst(==(c_opt), available_c)
                    if !isnothing(idx)
                        deleteat!(available_c, idx)
                    end
                end
            end
        end
    end
    
    return _build_sparse_matrix(workspace)
end

"""
Converts the internal PEG workspace adjacency lists into a highly optimized 
Julia SparseMatrixCSC.
"""
function _build_sparse_matrix(W::PEGWorkspace)
    total_edges = sum(W.v_degrees)
    
    rows = Int[]
    cols = Int[]
    vals = Int[]
    
    sizehint!(rows, total_edges)
    sizehint!(cols, total_edges)
    sizehint!(vals, total_edges)
    
    # Iterate through the check nodes to populate the standard sparse coordinate format
    for c in 1:W.m
        for v in W.c_adj[c]
            push!(rows, c)
            push!(cols, v)
            push!(vals, 1) # Standard binary parity-check matrix
        end
    end
    
    # Create the m x n sparse matrix
    H_sparse = sparse(rows, cols, vals, W.m, W.n)
    
    # Note: If your core library expects a CTMatrix, you can easily wrap this output 
    # wherever generate_peg is called (e.g., CTMatrix(H_sparse)).
    return H_sparse
end

"""
$(TYPEDSIGNATURES)

Generate a Quasi-Cyclic shift matrix using the QC-PEG algorithm. 
Maximizes the lifted graph's girth while strictly maintaining a block-circulant structure.

# Arguments
* `B::Matrix{Int}`: The base matrix (protograph/macro-graph) where `1` indicates an edge exists, and `0` indicates no edge.
* `Z::Int`: The lifting factor (circulant size).

# Returns
* A matrix of the same size as `B` containing the optimal circulant shifts `p ∈ [0, Z-1]`.
  Null edges are represented as `-1`.
"""
function progressive_edge_growth_QC(B::Matrix{Int}, Z::Int)
    mb, nb = size(B)
    shifts = fill(-1, mb, nb)
    
    # Lifted graph adjacency lists
    # Variables: 1 to nb*Z  |  Checks: 1 to mb*Z
    v_adj = [Int[] for _ in 1:(nb * Z)]
    c_adj = [Int[] for _ in 1:(mb * Z)]
    
    # Preallocated BFS Trackers (sized for the full lifted graph)
    depths = fill(-1, mb * Z)
    visited_v = fill(false, nb * Z)
    q = zeros(Int, mb * Z)
    candidate_shifts = Int[]
    
    for j in 1:nb
        for i in 1:mb
            if B[i, j] == 1
                # 1. Run BFS from the 0-th lifted variable node of macro-column j
                v_root = (j - 1) * Z + 1
                
                fill!(depths, -1)
                fill!(visited_v, false)
                visited_v[v_root] = true
                
                q_head = 1
                q_tail = 1
                
                # Initialize BFS with the direct check node connections of v_root
                for c in v_adj[v_root]
                    q[q_tail] = c
                    q_tail += 1
                    depths[c] = 1
                end
                
                # Expand the BFS tree
                while q_head < q_tail
                    curr_c = q[q_head]
                    q_head += 1
                    curr_depth = depths[curr_c]
                    
                    for v in c_adj[curr_c]
                        if !visited_v[v]
                            visited_v[v] = true
                            for next_c in v_adj[v]
                                if depths[next_c] == -1
                                    depths[next_c] = curr_depth + 1
                                    q[q_tail] = next_c
                                    q_tail += 1
                                end
                            end
                        end
                    end
                end
                
                # 2. Evaluate all possible shifts p ∈ [0, Z-1]
                best_depth = -2
                empty!(candidate_shifts)
                
                for p in 0:(Z - 1)
                    # The physical check node in macro-row i corresponding to shift p
                    c_target = (i - 1) * Z + p + 1
                    d = depths[c_target]
                    
                    # -1 means unreachable: connecting this edge creates infinite local girth!
                    if d == -1
                        d = typemax(Int) 
                    end
                    
                    # Keep track of the shifts that yield the absolute maximum depth
                    if d > best_depth
                        best_depth = d
                        empty!(candidate_shifts)
                        push!(candidate_shifts, p)
                    elseif d == best_depth
                        push!(candidate_shifts, p)
                    end
                end
                
                # 3. Select a shift (breaking ties randomly to prevent structured stopping sets)
                chosen_p = rand(candidate_shifts)
                shifts[i, j] = chosen_p
                
                # 4. Apply the Z-circulant block to the physical lifted graph
                for k in 0:(Z - 1)
                    v_node = (j - 1) * Z + k + 1
                    c_node = (i - 1) * Z + ((k + chosen_p) % Z) + 1
                    push!(v_adj[v_node], c_node)
                    push!(c_adj[c_node], v_node)
                end
            end
        end
    end
    
    return shifts
end

"""
$(TYPEDSIGNATURES)

Generate a lifted parity-check matrix from a protograph base matrix using the PEG algorithm.
Instead of random permutations, it intelligently selects permutation edges to maximize girth.

# Arguments
* `B::Matrix{Int}`: The protograph base matrix. `B[i, j]` represents the number of parallel edges.
* `Q::Int`: The lifting factor (number of replicas).
"""
function progressive_edge_growth_protograph(B::Matrix{Int}, Q::Int)
    mb, nb = size(B)
    M = mb * Q  # Total lifted check nodes
    N = nb * Q  # Total lifted variable nodes
    
    v_adj = [Int[] for _ in 1:N]
    c_adj = [Int[] for _ in 1:M]
    
    # Preallocated BFS Trackers
    depths = fill(-1, M)
    visited_v = fill(false, N)
    q_bfs = zeros(Int, M)
    candidate_cs = Int[]
    
    for j in 1:nb
        for i in 1:mb
            # Handle parallel edges in the base matrix
            for edge_idx in 1:B[i, j]
                # We must pair the Q copies of macro-variable j perfectly 
                # with the Q copies of macro-check i.
                avail_c = collect(1:Q)
                
                for q_v in 1:Q
                    v_root = (j - 1) * Q + q_v
                    
                    # 1. Run BFS from this specific variable replica
                    fill!(depths, -1)
                    fill!(visited_v, false)
                    visited_v[v_root] = true
                    
                    q_head = 1
                    q_tail = 1
                    
                    for c in v_adj[v_root]
                        q_bfs[q_tail] = c
                        q_tail += 1
                        depths[c] = 1
                    end
                    
                    while q_head < q_tail
                        curr_c = q_bfs[q_head]
                        q_head += 1
                        curr_depth = depths[curr_c]
                        
                        for v in c_adj[curr_c]
                            if !visited_v[v]
                                visited_v[v] = true
                                for next_c in v_adj[v]
                                    if depths[next_c] == -1
                                        depths[next_c] = curr_depth + 1
                                        q_bfs[q_tail] = next_c
                                        q_tail += 1
                                    end
                                end
                            end
                        end
                    end
                    
                    # 2. Evaluate depths only for the UNUSED replicas of macro-check i
                    best_depth = -2
                    empty!(candidate_cs)
                    
                    for (idx, q_c) in enumerate(avail_c)
                        c_target = (i - 1) * Q + q_c
                        d = depths[c_target]
                        
                        if d == -1
                            d = typemax(Int) # Infinite local girth!
                        end
                        
                        if d > best_depth
                            best_depth = d
                            empty!(candidate_cs)
                            push!(candidate_cs, idx)
                        elseif d == best_depth
                            push!(candidate_cs, idx)
                        end
                    end
                    
                    # 3. Randomly break ties to avoid structured stopping sets
                    chosen_idx = rand(candidate_cs)
                    q_c_chosen = avail_c[chosen_idx]
                    
                    # Remove the chosen check replica so it cannot be used by the other variable replicas
                    deleteat!(avail_c, chosen_idx)
                    
                    # 4. Wire the edge
                    c_node = (i - 1) * Q + q_c_chosen
                    push!(v_adj[v_root], c_node)
                    push!(c_adj[c_node], v_root)
                end
            end
        end
    end
    
    return _build_sparse_from_adj(v_adj, M, N)
end

"""
$(TYPEDSIGNATURES)

Generate a random LDPC parity-check matrix using the MacKay-Neal algorithm.
Enforces a strict "no 4-cycles" (girth >= 6) rule while attempting to balance row weights.

# Arguments
* `n::Int`: Number of variable nodes.
* `m::Int`: Number of check nodes.
* `deg_v::Vector{Int}`: Target degree for each variable node.
* `max_retries::Int`: Maximum number of times to backtrack on a column before throwing an error.
"""
function _generate_mackay_neal(n::Int, m::Int, deg_v::Vector{Int}; max_retries::Int=50)
    length(deg_v) == n || throw(ArgumentError("Length of deg_v must match n"))
    
    v_adj = [Int[] for _ in 1:n]
    c_adj = [Int[] for _ in 1:m]
    c_degrees = zeros(Int, m)
    
    forbidden = fill(false, m)
    candidate_cs = Int[]
    
    for j in 1:n
        retries = 0
        while retries < max_retries
            success = true
            
            for k in 1:deg_v[j]
                fill!(forbidden, false)
                
                # Fast 2-cycle and 4-cycle detection
                for c in v_adj[j]
                    forbidden[c] = true # Prevents double-edges (2-cycles)
                    for v_other in c_adj[c]
                        for c_other in v_adj[v_other]
                            forbidden[c_other] = true # Prevents 4-cycles
                        end
                    end
                end
                
                # Find available check nodes with the lowest current degree
                empty!(candidate_cs)
                min_deg = typemax(Int)
                
                for c in 1:m
                    if !forbidden[c]
                        if c_degrees[c] < min_deg
                            min_deg = c_degrees[c]
                            empty!(candidate_cs)
                            push!(candidate_cs, c)
                        elseif c_degrees[c] == min_deg
                            push!(candidate_cs, c)
                        end
                    end
                end
                
                # If we get boxed in, fail the column and trigger a retry
                if isempty(candidate_cs)
                    success = false
                    break
                end
                
                # Add the chosen edge
                c_opt = rand(candidate_cs)
                push!(v_adj[j], c_opt)
                push!(c_adj[c_opt], j)
                c_degrees[c_opt] += 1
            end
            
            if success
                break
            else
                # Rollback column j and retry
                for c in v_adj[j]
                    pop!(c_adj[c]) 
                    c_degrees[c] -= 1
                end
                empty!(v_adj[j])
                retries += 1
            end
        end
        
        if retries == max_retries
            error("MacKay-Neal algorithm got stuck on column $j. Try increasing m or reducing the column degrees.")
        end
    end
    
    return _build_sparse_from_adj(v_adj, m, n)
end

function Mackay_Neal(n::Int, m::Int, deg_v::Vector{Int}; max_retries::Int=50)
    return LDPCCode(_generate_mackay_neal(n, m, deg_v; max_retries=max_retries))
end

"""
Internal helper to convert adjacency lists into a Julia SparseMatrixCSC.
"""
function _build_sparse_from_adj(v_adj::Vector{Vector{Int}}, m::Int, n::Int)
    total_edges = sum(length, v_adj)
    rows = Int[]
    cols = Int[]
    sizehint!(rows, total_edges)
    sizehint!(cols, total_edges)
    
    for j in 1:n
        for c in v_adj[j]
            push!(rows, c)
            push!(cols, j)
        end
    end
    
    return sparse(rows, cols, fill(1, total_edges), m, n)
end

"""
$(TYPEDSIGNATURES)

Generate a Spatially Coupled LDPC (SC-LDPC) base matrix by "braiding" component 
matrices along a diagonal window. 

# Arguments
* `B_components::Vector{Matrix{Int}}`: A sequence of component matrices `[B_0, B_1, ..., B_w]` 
  that sum to the target uncoupled base matrix.
* `L::Int`: The coupling length (how many times to repeat the sequence down the diagonal).

# Returns
* A terminated SC-LDPC matrix of size `((L + w) * m_b) × (L * n_b)`.
"""
function _generate_spatially_coupled(B_components::Vector{Matrix{Int}}, L::Int)
    length(B_components) > 1 || throw(ArgumentError("Must provide at least 2 component matrices to couple."))
    
    w = length(B_components) - 1 # The memory (coupling width)
    mb, nb = size(B_components[1])
    
    # Ensure all components match in size
    for B in B_components
        size(B) == (mb, nb) || throw(ArgumentError("All component matrices must be the exact same size."))
    end
    
    # Terminated SC-LDPC matrices are rectangular (slightly lower rate, but extreme performance)
    M = (L + w) * mb
    N = L * nb
    
    H_SC = zeros(Int, M, N)
    
    # Stamp the component matrices down the diagonal
    for t in 1:L
        for i in 0:w
            # Calculate the block indices for this stamp
            row_start = (t - 1 + i) * mb + 1
            row_end   = (t + i) * mb
            col_start = (t - 1) * nb + 1
            col_end   = t * nb
            
            # Place the component matrix B_i
            H_SC[row_start:row_end, col_start:col_end] = B_components[i + 1]
        end
    end
    
    return H_SC
end

function SC_LDPCCode(B_components::Vector{Matrix{Int}}, L::Int)
    return LDPCCode(_generate_spatially_coupled(B_components, L))
end

"""
$(TYPEDSIGNATURES)

Generate a Euclidean Geometry EG(2, p) LDPC incidence matrix over a prime field.
This construction mathematically guarantees a girth of exactly 6 (zero 4-cycles) 
and completely avoids all trapping sets.

# Arguments
* `p::Int`: A prime number defining the finite field size (e.g., 7, 11, 31).

# Notes
* Columns represent 0-flats (Points in 2D space).
* Rows represent 1-flats (Lines in 2D space).
"""
function _generate_eg2(p::Int)
    # 1. Enumerate all Points (0-flats). 
    # In a 2D grid mod p, there are exactly p^2 points.
    points = Tuple{Int, Int}[]
    for x in 0:p-1
        for y in 0:p-1
            push!(points, (x, y))
        end
    end
    
    # 2. Enumerate all Lines (1-flats).
    # There are two types of lines in a 2D grid:
    lines = []
    
    # Type A: y = mx + b (mod p). There are p slopes and p intercepts.
    for m in 0:p-1
        for b in 0:p-1
            push!(lines, (m, b, :slope))
        end
    end
    
    # Type B: x = c (mod p). Vertical lines have infinite slope, defined by constant x.
    for c in 0:p-1
        push!(lines, (c, 0, :vertical))
    end
    
    # 3. Build the Incidence Matrix
    M = length(lines)  # p^2 + p rows
    N = length(points) # p^2 columns
    
    # Using SparseArrays coordinate format for efficiency
    rows = Int[]
    cols = Int[]
    
    for (i, line) in enumerate(lines)
        for (j, point) in enumerate(points)
            x, y = point
            
            # If the point satisfies the line equation, they intersect!
            if line[3] == :slope
                m, b = line[1], line[2]
                if y == mod(m * x + b, p)
                    push!(rows, i)
                    push!(cols, j)
                end
            elseif line[3] == :vertical
                c = line[1]
                if x == c
                    push!(rows, i)
                    push!(cols, j)
                end
            end
        end
    end
    
    return sparse(rows, cols, fill(1, length(rows)), M, N)
end

"""
$(TYPEDSIGNATURES)

Generate a generalized Euclidean Geometry EG(m, p) LDPC incidence matrix over a prime field.
This uses 0-flats (points) as variables and 1-flats (lines) as parity checks in m-dimensional space.

# Arguments
* `m::Int`: The dimension of the Euclidean space (e.g., 2, 3, 4).
* `p::Int`: A prime number defining the finite field size (e.g., 2, 3, 7).

# Returns
* A sparse incidence matrix of size `M x N`.
  `N = p^m` (Total points)
  `M = p^(m-1) * (p^m - 1) / (p - 1)` (Total distinct lines)
"""
function _generate_eg(m::Int, p::Int)
    # 1. Enumerate all points (0-flats) in GF(p)^m
    # We use Iterators.product to generate the m-dimensional grid
    points_iter = Iterators.product(fill(0:p-1, m)...)
    points = [collect(pt) for pt in points_iter]
    
    # 2. Enumerate all distinct directions (1D subspaces)
    # A direction is valid if its first non-zero element is exactly 1 (vector normalization).
    # This prevents us from treating vector [0, 2, 4] and [0, 1, 2] as different directions.
    directions = Vector{Int}[]
    for pt in points
        if all(pt .== 0) 
            continue 
        end
        
        idx = findfirst(x -> x != 0, pt)
        if pt[idx] == 1
            push!(directions, pt)
        end
    end
    
    # 3. Generate all distinct lines (1-flats)
    # A line is uniquely defined by the set of points it contains.
    unique_lines = Set{Vector{Vector{Int}}}()
    
    for d in directions
        for base_pt in points
            # Generate the line: base_pt + t * direction (mod p)
            line = [mod.(base_pt .+ t .* d, p) for t in 0:p-1]
            
            # Sort the points so lines generated from different starting points match exactly
            sort!(line)
            push!(unique_lines, line)
        end
    end
    
    lines = collect(unique_lines)
    
    # 4. Build the Incidence Matrix
    M_rows = length(lines)
    N_cols = length(points)
    
    # Fast O(1) lookup dictionary for column indices
    point_to_idx = Dict(pt => i for (i, pt) in enumerate(points))
    
    rows = Int[]
    cols = Int[]
    
    for (i, line) in enumerate(lines)
        for pt in line
            push!(rows, i)
            push!(cols, point_to_idx[pt])
        end
    end
    
    return sparse(rows, cols, fill(1, length(rows)), M_rows, N_cols)
end

function EuclideanGeometryCode(m::Int, p::Int)
    if m < 2
        throw(ArgumentError("Euclidean Geometry construction requires m >= 2"))
    elseif m == 2
        return LDPCCode(_generate_eg2(p))
    else
        return LDPCCode(_generate_eg(m, p))
    end
end

"""
$(TYPEDSIGNATURES)

Generate a Projective Geometry PG(m, p) LDPC incidence matrix over a prime field.
This construction guarantees zero 4-cycles and creates a highly symmetric matrix 
where every pair of lines in a plane intersects.

# Arguments
* `m::Int`: The dimension of the Projective space (e.g., 2, 3).
* `p::Int`: A prime number defining the finite field size (e.g., 2, 3, 7).

# Returns
* A sparse incidence matrix.
  Columns (Points): `N = (p^(m+1) - 1) / (p - 1)`
  Rows (Lines): `M` scales symmetrically based on the projective Grassmannian.
"""
function _generate_pg(m::Int, p::Int)
    # 1. Enumerate Points in PG(m, p)
    # A point is a 1D subspace in GF(p)^{m+1}.
    # We represent it by a normalized non-zero vector (the first non-zero element must be 1).
    points = Vector{Int}[]
    
    # Iterate through the (m+1)-dimensional grid
    for pt in Iterators.product(fill(0:p-1, m + 1)...)
        if all(pt .== 0)
            continue
        end
        
        vec = collect(pt)
        idx = findfirst(x -> x != 0, vec)
        
        # Only keep the vector if it is the "canonical" normalized representative
        if vec[idx] == 1
            push!(points, vec)
        end
    end
    
    N_cols = length(points)
    
    # Fast O(1) lookup dictionary for column indices
    point_to_idx = Dict(pt => i for (i, pt) in enumerate(points))
    
    # 2. Enumerate Lines in PG(m, p)
    # A line is a 2D subspace in GF(p)^{m+1}.
    # It is spanned by any two distinct Projective points u and v.
    # A projective line contains exactly p + 1 points.
    unique_lines = Set{Vector{Int}}()
    
    for i in 1:N_cols
        for j in i+1:N_cols
            u = points[i]
            v = points[j]
            
            # The line contains point `u`...
            line_pts = Int[i] 
            
            # ...and `t * u + v` for all t in GF(p)
            for t in 0:p-1
                # 1. Take the linear combination
                combo = mod.(t .* u .+ v, p)
                
                # 2. Normalize it back into a valid Projective Point representation
                idx = findfirst(x -> x != 0, combo)
                inv_val = invmod(combo[idx], p) # Find the modular inverse
                norm_combo = mod.(combo .* inv_val, p)
                
                push!(line_pts, point_to_idx[norm_combo])
            end
            
            # Sort the indices so identical lines match the same Set hash
            sort!(line_pts)
            push!(unique_lines, line_pts)
        end
    end
    
    lines = collect(unique_lines)
    M_rows = length(lines)
    
    # 3. Build the Incidence Matrix
    rows = Int[]
    cols = Int[]
    
    for (i, line) in enumerate(lines)
        for pt_idx in line
            push!(rows, i)
            push!(cols, pt_idx)
        end
    end
    
    return sparse(rows, cols, fill(1, length(rows)), M_rows, N_cols)
end

function ProjectiveGeometryCode(m::Int, p::Int)
    if m < 2
        throw(ArgumentError("Projective Geometry construction requires m >= 2"))
    end
    return LDPCCode(_generate_pg(m, p))
end

"""
    _gallager_H(n::Int, wc::Int, wr::Int)

Generates a regular (wc, wr) LDPC parity-check matrix using Gallager's original construction.
"""
function _gallager_H(n::Int, wc::Int, wr::Int)
    (n * wc) % wr == 0 || throw(ArgumentError("n * wc must be perfectly divisible by wr"))
    
    m = div(n * wc, wr)
    block_size = div(n, wr)
    
    H = zeros(Int, m, n)
    
    # 1. Build the fundamental base block (Block 0)
    # The dot product of any two distinct rows here is strictly 0.
    for i in 1:block_size
        start_col = (i - 1) * wr + 1
        end_col = i * wr
        H[i, start_col:end_col] .= 1
    end
    
    # 2. Build the remaining (wc - 1) blocks via random column permutations
    base_block = H[1:block_size, :]
    
    for block in 1:(wc - 1)
        start_row = block * block_size + 1
        end_row = (block + 1) * block_size
        
        # Apply a random permutation to the columns of the base block
        perm = randperm(n)
        H[start_row:end_row, :] = base_block[:, perm]
    end
    
    return H
end

function GallagerCode(n::Int, wc::Int, wr::Int)
    return LDPCCode(_gallager_H(n, wc, wr))
end