# Copyright (c) 2024 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
       # Region Graphs
#############################

struct Region
    id::Vector{Int}              # The variable nodes in this region
    parents::Vector{Int}         # Indices of parent regions in the master array
    ancestors::Vector{Int}       # Indices of ancestor regions
    subregions::Vector{Int}      # Indices of child regions
    overcounting_number::Int
end

struct RegionGraph
    regions::Vector{Region}
end

id(r::Region) = r.id
parents(r::Region) = r.parents
ancestors(r::Region) = r.ancestors
subregions(r::Region) = r.subregions
overcounting_number(r::Region) = r.overcounting_number
regions(R::RegionGraph) = R.regions
base_regions(R::RegionGraph) = (r for r in R.regions if isempty(r.parents))
leaves(R::RegionGraph) = (r for r in R.regions if isempty(r.subregions))
outer_regions(R::RegionGraph) = base_regions(R)
basic_clusters(R::RegionGraph) = base_regions(R)

==(r1::Region, r2::Region) = r1.id == r2.id

function canonical_region_graph(H::CTMatrixTypes)
    num_check, num_var = size(H)
    
    # Initialize as BitSets for O(1) intersections later
    check_adj_list = [BitSet() for _ in 1:num_check]
    
    if H isa SparseMatrixCSC
        # Hyper-fast sweep for sparse matrices
        rows = rowvals(H)
        for c in 1:num_var
            for i in nzrange(H, c)
                push!(check_adj_list[rows[i]], c)
            end
        end
    else
        # Standard sweep for dense matrices
        for r in 1:num_check
            for c in 1:num_var
                if !iszero(H[r, c])
                    push!(check_adj_list[r], c)
                end
            end
        end
    end
    
    return region_graph_from_base_nodes(check_adj_list)
end

canonical_region_graph(L::AbstractLDPCCode) = canonical_region_graph(parity_check_matrix(L))

function region_graph_from_base_nodes(R::Vector{Vector{Int}})
    return region_graph_from_base_nodes([BitSet(r) for r in R])
end

function region_graph_from_base_nodes(base_nodes::Vector{BitSet})
    isempty(base_nodes) && return RegionGraph(Region[])

    # ---------------------------------------------------------
    # 1. Initialize Mutable Builder Arrays
    # ---------------------------------------------------------
    ids = copy(base_nodes)
    n_initial = length(ids)
    
    parents = [Int[] for _ in 1:n_initial]
    ancestors = [Int[] for _ in 1:n_initial]
    subregions = [Int[] for _ in 1:n_initial]
    overcounting = ones(Int, n_initial)

    # ---------------------------------------------------------
    # 2. The O(N^3) Intersection Closure Loop
    # ---------------------------------------------------------
    left = 1
    right = n_initial

    while left < right
        for r1 in left:(right - 1)
            for r2 in (left + 1):right
                # Fast BitSet intersection (zero heap allocations)
                cap = intersect(ids[r1], ids[r2])
                
                if !isempty(cap) && cap != ids[r1] && cap != ids[r2]
                    
                    # Search if this region intersection already exists
                    found_idx = -1
                    for i in 1:length(ids)
                        if cap == ids[i]
                            found_idx = i
                            break
                        end
                    end
                    
                    if found_idx != -1
                        # --- UPDATE EXISTING REGION ---
                        # 1. Remove parents that are ancestors of r1 or r2 [cite: 7, 8, 9, 10]
                        filter!(p -> (p ∉ ancestors[r1]) && (p ∉ ancestors[r2]), parents[found_idx])
                        
                        # 2. Add r1 and r2 as parents [cite: 11]
                        r1 ∉ parents[found_idx] && push!(parents[found_idx], r1)
                        r2 ∉ parents[found_idx] && push!(parents[found_idx], r2)
                        
                        # 3. Recompute ancestors [cite: 12, 13]
                        empty!(ancestors[found_idx])
                        for p in parents[found_idx]
                            append!(ancestors[found_idx], ancestors[p])
                            push!(ancestors[found_idx], p)
                        end
                        unique!(ancestors[found_idx])
                        
                        # 4. Add this region to its ancestors' subregions [cite: 14, 15]
                        for anc in ancestors[found_idx]
                            found_idx ∉ subregions[anc] && push!(subregions[anc], found_idx)
                        end
                    else
                        # --- CREATE NEW REGION ---
                        # 1. Determine Ancestors
                        new_ancs = Int[]
                        append!(new_ancs, ancestors[r1])
                        append!(new_ancs, ancestors[r2])
                        push!(new_ancs, r1)
                        push!(new_ancs, r2)
                        unique!(new_ancs)
                        
                        # 2. Calculate Overcounting Number
                        c_r = 1
                        for anc in new_ancs
                            c_r -= overcounting[anc]
                        end
                        
                        # 3. Push to Builder Arrays
                        push!(ids, cap)
                        push!(parents, [r1, r2])
                        push!(ancestors, new_ancs)
                        push!(subregions, Int[])
                        push!(overcounting, c_r)
                        
                        # 4. Update Ancestors' Subregions
                        new_idx = length(ids)
                        for anc in new_ancs
                            push!(subregions[anc], new_idx)
                        end
                    end
                end
            end
        end
        left = right + 1
        right = length(ids)
    end
    
    # ---------------------------------------------------------
    # 3. Compute Correct Overcounting Numbers
    # ---------------------------------------------------------
    # Overcounting numbers must be calculated top-down after all 
    # ancestors have been finalized to ensure validity.
    sizes = Int[length(ids[i]) for i in 1:length(ids)]
    order = sortperm(sizes, rev=true)
    for i in order
        overcounting[i] = 1
        for anc in ancestors[i]
            overcounting[i] -= overcounting[anc]
        end
    end
    
    # ---------------------------------------------------------
    # 4. Finalize and Freeze into Immutable Structs
    # ---------------------------------------------------------
    regions = Vector{Region}(undef, length(ids))
    for i in 1:length(ids)
        sorted_id = sort!(collect(ids[i]))
        regions[i] = Region(sorted_id, parents[i], ancestors[i], subregions[i], overcounting[i])
    end
    
    return RegionGraph(regions)
end

function is_valid_region_graph(R::RegionGraph)
    # 1. Ensure the sum of c_r for all regions containing a specific variable node is exactly 1.
    all_vars = BitSet()
    for r in R.regions
        union!(all_vars, r.id)
    end

    for v in all_vars
        c_r_sum = 0
        for r in R.regions
            if v in r.id
                c_r_sum += r.overcounting_number
            end
        end
        c_r_sum != 1 && return false
    end

    # 2. Ensure the sum of c_r for any region and all of its ancestors is exactly 1.
    for r in R.regions
        c_r_sum = r.overcounting_number
        for anc_idx in r.ancestors
            c_r_sum += R.regions[anc_idx].overcounting_number
        end
        c_r_sum != 1 && return false
    end

    return true
end

function remove_zero_overcounting_numbers(R::RegionGraph)
    isempty(R.regions) && return R

    # 1. Unpack into mutable arrays for the Builder pattern
    N = length(R.regions)
    active = trues(N)
    
    parents = [copy(r.parents) for r in R.regions]
    ancestors = [copy(r.ancestors) for r in R.regions]
    subregions = [copy(r.subregions) for r in R.regions]
    
    # 2. Iterative Removal Logic
    changed = true
    while changed
        changed = false
        for i in 1:N
            if active[i] && R.regions[i].overcounting_number == 0
                remove_flag = true
                
                # If it's a leaf, ensure removing it doesn't disconnect the graph
                if isempty(subregions[i])
                    remove_flag = false
                    
                    # Check if any two parents share a common ancestor
                    for p1 in parents[i]
                        for p2 in parents[i]
                            if p1 != p2
                                common_ancs = intersect(ancestors[p1], ancestors[p2])
                                if !isempty(common_ancs)
                                    remove_flag = true
                                    break
                                end
                            end
                        end
                        remove_flag && break
                    end
                end

                if remove_flag
                    active[i] = false
                    changed = true
                    
                    # A. Update Descendents
                    for child in subregions[i]
                        filter!(x -> x != i, ancestors[child])
                        
                        if i in parents[child]
                            filter!(x -> x != i, parents[child])
                            append!(parents[child], parents[i])
                            unique!(parents[child])
                        end
                    end
                    
                    # B. Update Ancestors
                    for anc in ancestors[i]
                        filter!(x -> x != i, subregions[anc])
                    end
                end
            end
        end
    end

    # 3. Build the Compaction Map (old index -> new index)
    old_to_new = zeros(Int, N)
    new_idx = 1
    for i in 1:N
        if active[i]
            old_to_new[i] = new_idx
            new_idx += 1
        end
    end
    
    # 4. Re-wire the Graph and Freeze
    new_regions = Vector{Region}(undef, new_idx - 1)
    
    for i in 1:N
        if active[i]
            mapped_idx = old_to_new[i]
            
            # Re-wire pointers, dropping any that point to deleted regions
            new_parents = Int[old_to_new[p] for p in parents[i] if active[p]]
            new_ancestors = Int[old_to_new[a] for a in ancestors[i] if active[a]]
            new_subregions = Int[old_to_new[s] for s in subregions[i] if active[s]]
            
            new_regions[mapped_idx] = Region(
                R.regions[i].id, 
                new_parents, 
                new_ancestors, 
                new_subregions, 
                R.regions[i].overcounting_number
            )
        end
    end

    return RegionGraph(new_regions)
end

function remove_generational_skips(R::RegionGraph)
    # 1. Unpack into mutable arrays
    N = length(R.regions)
    parents = [copy(r.parents) for r in R.regions]
    ancestors = [copy(r.ancestors) for r in R.regions]
    subregions = [copy(r.subregions) for r in R.regions]

    changed = false

    # 2. Sweep all regions to find generational skips
    for i in 1:N
        if isempty(parents[i])
            continue
        end
        
        valid_parents = Int[]
        
        for candidate_p in parents[i]
            is_skip = false
            
            # Check if this candidate is an ancestor of ANY OTHER parent
            for other_p in parents[i]
                if candidate_p != other_p
                    if candidate_p in ancestors[other_p]
                        is_skip = true
                        changed = true
                        break
                    end
                end
            end
            
            if !is_skip
                push!(valid_parents, candidate_p)
            end
        end
        
        # 3. If a skip was found, sever the redundant edges
        if length(valid_parents) != length(parents[i])
            # Identify which parents we are dropping
            dropped_parents = setdiff(parents[i], valid_parents)
            
            # Update child's parent list
            parents[i] = valid_parents
            
            # Update the dropped parents' subregion list
            for dp in dropped_parents
                filter!(x -> x != i, subregions[dp])
            end
        end
    end

    # 4. If nothing changed, return original. Otherwise, freeze the new topology.
    if !changed
        return R
    end

    new_regions = Vector{Region}(undef, N)
    for i in 1:N
        new_regions[i] = Region(
            R.regions[i].id,
            parents[i],
            ancestors[i],
            subregions[i],
            R.regions[i].overcounting_number
        )
    end

    return RegionGraph(new_regions)
end

function triangulate_base_regions(H::CTMatrixTypes)
    num_check, num_var = size(H)
    
    # 1. Build the Primal Graph (Variables connected if they share a check)
    adj = [BitSet() for _ in 1:num_var]
    for c in 1:num_check
        # Find all variables in this check
        vars_in_check = Int[]
        for v in 1:num_var
            if !iszero(H[c, v])
                push!(vars_in_check, v)
            end
        end
        
        # Fully connect the clique for this check
        for i in 1:length(vars_in_check)
            for j in (i+1):length(vars_in_check)
                v1, v2 = vars_in_check[i], vars_in_check[j]
                push!(adj[v1], v2)
                push!(adj[v2], v1)
            end
        end
    end
    
    # 2. Min-Degree Elimination (Chordalization)
    active = trues(num_var)
    maximal_cliques = Vector{BitSet}()
    
    for _ in 1:num_var
        # Find active node with minimum degree
        min_deg = Inf
        best_v = -1
        for v in 1:num_var
            if active[v]
                deg = count(n -> active[n], adj[v])
                if deg < min_deg
                    min_deg = deg
                    best_v = v
                end
            end
        end
        
        best_v == -1 && break
        
        # Form a clique with the node and its active neighbors
        clique = BitSet(best_v)
        neighbors = Int[]
        for n in adj[best_v]
            if active[n]
                push!(clique, n)
                push!(neighbors, n)
            end
        end
        
        push!(maximal_cliques, clique)
        
        # Add fill-in edges between all neighbors (Triangulation)
        for i in 1:length(neighbors)
            for j in (i+1):length(neighbors)
                n1, n2 = neighbors[i], neighbors[j]
                push!(adj[n1], n2)
                push!(adj[n2], n1)
            end
        end
        
        active[best_v] = false # Eliminate node
    end
    
    # 3. Filter subsets to keep ONLY maximal cliques
    sort!(maximal_cliques, by=length, rev=true)
    final_cliques = Vector{BitSet}()
    
    for clique in maximal_cliques
        is_subset = false
        for fc in final_cliques
            if issubset(clique, fc)
                is_subset = true
                break
            end
        end
        if !is_subset
            push!(final_cliques, clique)
        end
    end
    
    return final_cliques
end

function message_passing_order(R::RegionGraph)
    N = length(R.regions)
    
    # Track how many direct children each region is waiting on
    in_degree = zeros(Int, N)
    for i in 1:N
        for p in R.regions[i].parents
            in_degree[p] += 1
        end
    end
    
    # Initialize queue with strict leaves (in_degree == 0)
    queue = Int[]
    for i in 1:N
        if in_degree[i] == 0
            push!(queue, i)
        end
    end
    
    order = Int[]
    sizehint!(order, N)
    
    # Topological BFS (Upwards: Leaves -> Roots)
    head = 1
    while head <= length(queue)
        curr = queue[head]
        head += 1
        push!(order, curr)
        
        # Move up to parents
        for p in R.regions[curr].parents
            in_degree[p] -= 1
            if in_degree[p] == 0
                push!(queue, p)
            end
        end
    end
    
    # Sanity check for isolated cycles
    if length(order) != N
        error("Region Graph contains an invalid cycle. Cannot determine order.")
    end
    
    return order
end

# ==============================================================================
# DISPLAY & SHOW FUNCTIONS
# ==============================================================================

# --- Region Display ---

# Compact 1-liner (e.g., "Region(id={1, 2, 4}, c_r=1)")
function Base.show(io::IO, r::Region)
    # Convert BitSet to a sorted array for clean printing
    id_str = join(sort!(collect(r.id)), ", ")
    print(io, "Region(id={", id_str, "}, c_r=", r.overcounting_number, ")")
end


# --- RegionGraph Display ---

# Compact 1-liner
function Base.show(io::IO, R::RegionGraph)
    print(io, "RegionGraph(", length(R.regions), " regions)")
end

# Rich REPL Display
function Base.show(io::IO, ::MIME"text/plain", R::RegionGraph)
    num_regs = length(R.regions)
    print(io, "RegionGraph with ", num_regs, " region")
    num_regs != 1 && print(io, "s")
    
    if num_regs > 0
        println(io, ":")
        
        # Calculate some quick topological stats
        max_size = maximum(length(r.id) for r in R.regions)
        num_leaves = sum(1 for r in R.regions if isempty(r.subregions))
        num_roots = sum(1 for r in R.regions if isempty(r.parents))
        
        println(io, "  Max region size: ", max_size, " variables")
        println(io, "  Base regions:    ", num_roots)
        println(io, "  Leaf regions:    ", num_leaves)
        println(io)
        
        # Print a preview of the first few regions so it doesn't flood the terminal
        limit = min(num_regs, 10)
        for i in 1:limit
            print(io, "  [$i] ")
            
            # Print the region ID and overcounting number
            id_str = join(sort!(collect(R.regions[i].id)), ", ")
            print(io, "id={", id_str, "}, c_r=", R.regions[i].overcounting_number)
            
            # Add parent/child info context
            parents = R.regions[i].parents
            if !isempty(parents)
                print(io, "  (Parents: ", join(parents, ", "), ")")
            end
            
            println(io)
        end
        
        if num_regs > 10
            println(io, "  ⋮ (", num_regs - 10, " more regions)")
        end
    end
end

#############################
            # GBP
#############################

# ==============================================================================
# GENERALIZED BELIEF PROPAGATION WORKSPACE
# ==============================================================================

struct GBPWorkspace
    num_regions::Int
    num_edges::Int
    
    # Shared State (Log-Beliefs)
    log_beliefs::Vector{Float64}
    log_belief_offsets::Vector{Int}
    
    # Edge State (Messages: Parent -> Child)
    edge_to_parent::Vector{Int}
    edge_to_child::Vector{Int}
    messages::Vector{Float64}
    message_offsets::Vector{Int}
    
    # Marginalization Maps
    marg_maps::Vector{Vector{Int}}
    
    # Topological Schedule & Early Termination Buffers
    edge_update_order::Vector{Int}
    is_decimated::Vector{Bool}
    current_synd_buffer::Vector{UInt8}
end

function init_gbp_workspace(R::RegionGraph, H::CTMatrixTypes)
    num_regions = length(R.regions)
    num_check, num_var = size(H)
    
    # ---------------------------------------------------------
    # 1. Allocate Log-Belief Strides
    # ---------------------------------------------------------
    log_belief_offsets = zeros(Int, num_regions + 1)
    log_belief_offsets[1] = 1
    
    for i in 1:num_regions
        num_states = 2^(length(R.regions[i].id))
        log_belief_offsets[i+1] = log_belief_offsets[i] + num_states
    end
    
    total_belief_states = log_belief_offsets[end] - 1
    log_beliefs = zeros(Float64, total_belief_states)
    
    # ---------------------------------------------------------
    # 2. Enumerate Edges & Allocate Message Strides
    # ---------------------------------------------------------
    edge_to_parent = Int[]
    edge_to_child = Int[]
    message_offsets = Int[1]
    
    for c_idx in 1:num_regions
        for p_idx in R.regions[c_idx].parents
            push!(edge_to_parent, p_idx)
            push!(edge_to_child, c_idx)
            
            num_child_states = 2^(length(R.regions[c_idx].id))
            push!(message_offsets, message_offsets[end] + num_child_states)
        end
    end
    
    num_edges = length(edge_to_parent)
    total_message_states = message_offsets[end] - 1
    messages = zeros(Float64, total_message_states)
    
    # ---------------------------------------------------------
    # 3. Build the Bitwise Marginalization Maps
    # ---------------------------------------------------------
    marg_maps = Vector{Vector{Int}}(undef, num_edges)
    
    for e in 1:num_edges
        p_idx = edge_to_parent[e]
        c_idx = edge_to_child[e]
        
        p_vars = sort!(collect(R.regions[p_idx].id))
        c_vars = sort!(collect(R.regions[c_idx].id))
        
        p_bit_indices = Int[]
        for cv in c_vars
            idx = findfirst(==(cv), p_vars) 
            push!(p_bit_indices, idx - 1) 
        end
        
        num_p_states = 2^(length(p_vars))
        map_e = zeros(Int, num_p_states)
        
        for s_p in 0:(num_p_states - 1)
            s_c = 0
            for (c_bit_idx, p_bit_idx) in enumerate(p_bit_indices)
                bit_val = (s_p >> p_bit_idx) & 1
                s_c |= (bit_val << (c_bit_idx - 1))
            end
            map_e[s_p + 1] = s_c + 1
        end
        marg_maps[e] = map_e
    end
    
    # ---------------------------------------------------------
    # 4. Generate the Two-Way Topological Schedule
    # ---------------------------------------------------------
    region_order_up = message_passing_order(R)     # Leaves -> Roots
    region_order_down = reverse(region_order_up)   # Roots -> Leaves
    
    edge_update_order = Int[]
    sizehint!(edge_update_order, 2 * num_edges)
    
    # Downward Pass
    for r in region_order_down
        for e in 1:num_edges
            if edge_to_parent[e] == r
                push!(edge_update_order, e)
            end
        end
    end
    
    # Upward Pass
    for r in region_order_up
        for e in 1:num_edges
            if edge_to_child[e] == r
                push!(edge_update_order, e)
            end
        end
    end
    
    # ---------------------------------------------------------
    # 5. Initialize Decimation & Syndrome Buffers
    # ---------------------------------------------------------
    is_decimated = zeros(Bool, num_var)
    current_synd_buffer = zeros(UInt8, num_check)
    
    return GBPWorkspace(
        num_regions, num_edges,
        log_beliefs, log_belief_offsets,
        edge_to_parent, edge_to_child, messages, message_offsets,
        marg_maps,
        edge_update_order, is_decimated, current_synd_buffer
    )
end

# ==============================================================================
# GBP DECODING ENGINE
# ==============================================================================

"""
Numerically stable Log-Sum-Exp for combining probabilities in the log domain.
"""
@inline function _log_add_exp(x::Float64, y::Float64)
    x == -Inf && return y
    y == -Inf && return x
    max_val = max(x, y)
    return max_val + log1p(exp(-abs(x - y)))
end

# The Max-Log Approximation
@inline function _log_add_exp_approx(x::Float64, y::Float64)
    return max(x, y)
end

"""
Master GBP API Wrapper. 
Executes full Generalized Belief Propagation with early termination.
"""
function gbp_decode!(W::GBPWorkspace, R::RegionGraph, H::CTMatrixTypes, total_llrs::Vector{Float64};
                     target_syndrome::Vector{UInt8} = zeros(UInt8, size(H, 1)),
                     decimation_type::Val = Val(:none),
                     dec_thresh::Float64 = 10.0,
                     dec_rounds::Int = 3,
                     max_iter::Int = 50, 
                     damping::Float64 = 0.5)
    
    num_check, num_var = size(H)
    
    # Reset decimation flags
    fill!(W.is_decimated, false)
    
    # 1. Initialize Log-Beliefs (Assigns LLRs and Quantum Parity Constraints)
    init_region_beliefs!(W, R, H, total_llrs, target_syndrome)
    
    marg_buffer = Float64[]
    
    @inbounds begin
        for iter in 1:max_iter
            
            # ---------------------------------------------------------
            # 2. MESSAGE PASSING ENGINE (Parent -> Child Updates)
            # ---------------------------------------------------------
            for e in W.edge_update_order
                p_idx = W.edge_to_parent[e]
                c_idx = W.edge_to_child[e]
                
                p_off = W.log_belief_offsets[p_idx] - 1
                c_off = W.log_belief_offsets[c_idx] - 1
                m_off = W.message_offsets[e] - 1
                
                num_c_states = W.message_offsets[e+1] - W.message_offsets[e]
                num_p_states = W.log_belief_offsets[p_idx+1] - W.log_belief_offsets[p_idx]
                
                if length(marg_buffer) < num_c_states
                    resize!(marg_buffer, num_c_states)
                end
                for i in 1:num_c_states; marg_buffer[i] = -Inf; end
                
                map_e = W.marg_maps[e]
                for s_p in 1:num_p_states
                    s_c = map_e[s_p]
                    parent_val = W.log_beliefs[p_off + s_p]
                    marg_buffer[s_c] = _log_add_exp(marg_buffer[s_c], parent_val)
                end
                
                max_belief = -Inf
                for s_c in 1:num_c_states
                    old_msg = W.messages[m_off + s_c]
                    child_belief = W.log_beliefs[c_off + s_c]
                    
                    new_msg_raw = marg_buffer[s_c] - (child_belief - old_msg)
                    new_msg = (damping * new_msg_raw) + ((1.0 - damping) * old_msg)
                    
                    W.messages[m_off + s_c] = new_msg
                    updated_belief = child_belief + (new_msg - old_msg)
                    W.log_beliefs[c_off + s_c] = updated_belief
                    
                    if updated_belief > max_belief
                        max_belief = updated_belief
                    end
                end
                
                if max_belief > -Inf
                    for s_c in 1:num_c_states
                        W.log_beliefs[c_off + s_c] -= max_belief
                    end
                end
            end
            
            # ---------------------------------------------------------
            # 3. DECIMATION HOOK
            # ---------------------------------------------------------
            _apply_gbp_decimation!(decimation_type, W, R, num_var, iter, dec_thresh, dec_rounds)
            
            # ---------------------------------------------------------
            # 4. EARLY TERMINATION (Extract bits and check syndrome)
            # ---------------------------------------------------------
            current_bits = extract_hard_decisions(W, R, num_var)
            is_valid = true
            
            if H isa SparseMatrixCSC
                # FAST PATH: Column-major iteration for CSC Matrices
                fill!(W.current_synd_buffer, 0x00)
                rows = rowvals(H)
                
                for v in 1:num_var
                    if current_bits[v] == 0x01
                        for i in nzrange(H, v)
                            W.current_synd_buffer[rows[i]] ⊻= 0x01
                        end
                    end
                end
                
                for c in 1:num_check
                    if W.current_synd_buffer[c] != target_syndrome[c]
                        is_valid = false
                        break
                    end
                end
            else
                # SLOW PATH: Dense matrix check
                for c in 1:num_check
                    syn = target_syndrome[c]
                    for v in 1:num_var
                        if H[c, v] != 0
                            syn ⊻= current_bits[v]
                        end
                    end
                    if syn != 0x00
                        is_valid = false
                        break
                    end
                end
            end
            
            if is_valid
                return true, current_bits, iter
            end
        end
    end
    
    final_bits = extract_hard_decisions(W, R, num_var)
    return false, final_bits, max_iter
end

"""
Initializes the Log-Beliefs for all regions.
Assigns channel LLRs and parity checks to strictly one region to prevent double-counting.
"""
function init_region_beliefs!(W::GBPWorkspace, R::RegionGraph, H::CTMatrixTypes, channel_llrs::Vector{Float64}, target_syndrome::Vector{UInt8})
    num_check, num_var = size(H)
    num_regions = length(R.regions)
    
    # 1. Map Checks to BitSets
    check_adj = [BitSet() for _ in 1:num_check]
    if H isa SparseMatrixCSC
        rows = rowvals(H)
        for c in 1:num_var
            for i in nzrange(H, c)
                push!(check_adj[rows[i]], c)
            end
        end
    else
        for c in 1:num_check
            for v in 1:num_var
                if !iszero(H[c, v])
                    push!(check_adj[c], v)
                end
            end
        end
    end
    
    # 2. Assign Variables to EXACTLY ONE Region
    var_assignment = zeros(Int, num_var)
    for v in 1:num_var
        for r_idx in 1:num_regions
            if v in R.regions[r_idx].id
                var_assignment[v] = r_idx
                break
            end
        end
        if var_assignment[v] == 0
            error("Variable $v is not contained in any region!")
        end
    end
    
    # 3. Assign Parity Checks to EXACTLY ONE Region
    check_assignment = zeros(Int, num_check)
    for c in 1:num_check
        for r_idx in 1:num_regions
            if issubset(check_adj[c], R.regions[r_idx].id)
                check_assignment[c] = r_idx
                break
            end
        end
        if check_assignment[c] == 0
            error("Parity check $c is not fully contained in any region!")
        end
    end
    
    # 4. Initialize the Log-Beliefs Array
    fill!(W.log_beliefs, 0.0)
    fill!(W.messages, 0.0) # Reset messages to 0 (log(1))
    
    @inbounds for r_idx in 1:num_regions
        r_vars = sort!(collect(R.regions[r_idx].id))
        num_vars = length(r_vars)
        num_states = 2^num_vars
        offset = W.log_belief_offsets[r_idx] - 1
        
        # Find which variables and checks were assigned to this specific region
        assigned_vars = [i for i in 1:num_vars if var_assignment[r_vars[i]] == r_idx]
        assigned_checks = [c for c in 1:num_check if check_assignment[c] == r_idx]
        
        for s in 0:(num_states - 1)
            is_valid = true
            
            # A. Enforce Parity Checks assigned to this region
            for c in assigned_checks
                parity = 0
                for (bit_idx, v) in enumerate(r_vars)
                    if v in check_adj[c]
                        bit_val = (s >> (bit_idx - 1)) & 1
                        parity ⊻= bit_val
                    end
                end
                
                # Check against the target syndrome instead of hardcoded 0
                if parity != target_syndrome[c]
                    is_valid = false
                    break
                end
            end
            
            # Invalid states get -Inf log-probability (0 probability)
            if !is_valid
                W.log_beliefs[offset + s + 1] = -Inf
                continue
            end
            
            # B. Inject Channel Information for variables assigned to this region
            # If bit_val == 1, we penalize it by -LLR. 
            state_belief = 0.0
            for bit_idx in assigned_vars
                bit_val = (s >> (bit_idx - 1)) & 1
                if bit_val == 1
                    state_belief -= channel_llrs[r_vars[bit_idx]]
                end
            end
            
            W.log_beliefs[offset + s + 1] = state_belief
        end
    end
end

"""
Extracts the final hard decisions by marginalizing the converged Log-Beliefs.
"""
function extract_hard_decisions(W::GBPWorkspace, R::RegionGraph, num_var::Int)
    hard_decisions = zeros(UInt8, num_var)
    
    @inbounds for v in 1:num_var
        # Find the first region containing this variable
        target_r = 0
        for r_idx in 1:length(R.regions)
            if v in R.regions[r_idx].id
                target_r = r_idx
                break
            end
        end
        
        r_vars = sort!(collect(R.regions[target_r].id))
        bit_pos = findfirst(==(v), r_vars) - 1 # 0-indexed for bit shifting
        
        num_states = 2^(length(r_vars))
        offset = W.log_belief_offsets[target_r] - 1
        
        log_prob_0 = -Inf
        log_prob_1 = -Inf
        
        # Marginalize the region's beliefs down to this specific variable
        for s in 0:(num_states - 1)
            belief = W.log_beliefs[offset + s + 1]
            belief == -Inf && continue
            
            bit_val = (s >> bit_pos) & 1
            if bit_val == 0
                log_prob_0 = _log_add_exp(log_prob_0, belief)
            else
                log_prob_1 = _log_add_exp(log_prob_1, belief)
            end
        end
        
        # Hard decision based on max log-probability
        hard_decisions[v] = log_prob_1 > log_prob_0 ? 0x01 : 0x00
    end
    
    return hard_decisions
end

# ==============================================================================
# GBP DECIMATION HOOKS
# ==============================================================================

@inline _apply_gbp_decimation!(::Val{:none}, W, R, num_var, iter, thresh, rounds) = nothing

function _apply_gbp_decimation!(::Val{:hard}, W::GBPWorkspace, R::RegionGraph, num_var::Int, iter::Int, thresh::Float64, rounds::Int)
    if iter % rounds != 0
        return
    end
    
    @inbounds for v in 1:num_var
        if !W.is_decimated[v]
            # 1. Find the first region to calculate the marginal LLR
            target_r = 0
            for r_idx in 1:length(R.regions)
                if v in R.regions[r_idx].id
                    target_r = r_idx
                    break
                end
            end
            
            r_vars = sort!(collect(R.regions[target_r].id))
            bit_pos = findfirst(==(v), r_vars) - 1 
            num_states = 2^(length(r_vars))
            offset = W.log_belief_offsets[target_r] - 1
            
            log_prob_0, log_prob_1 = -Inf, -Inf
            
            for s in 0:(num_states - 1)
                belief = W.log_beliefs[offset + s + 1]
                belief == -Inf && continue
                
                bit_val = (s >> bit_pos) & 1
                if bit_val == 0
                    log_prob_0 = _log_add_exp(log_prob_0, belief)
                else
                    log_prob_1 = _log_add_exp(log_prob_1, belief)
                end
            end
            
            # LLR = Log(P(0)) - Log(P(1))
            llr = log_prob_0 - log_prob_1
            
            # 2. Check threshold and Decimate
            if abs(llr) > thresh
                W.is_decimated[v] = true
                locked_val = llr > 0 ? 0 : 1
                
                # 3. Banish contradictory states in ALL regions containing this variable
                for r_idx in 1:length(R.regions)
                    if v in R.regions[r_idx].id
                        r_vars_all = sort!(collect(R.regions[r_idx].id))
                        local_bit_pos = findfirst(==(v), r_vars_all) - 1
                        local_num_states = 2^(length(r_vars_all))
                        local_offset = W.log_belief_offsets[r_idx] - 1
                        
                        for s in 0:(local_num_states - 1)
                            # If the state's bit doesn't match the locked value, destroy it
                            if ((s >> local_bit_pos) & 1) != locked_val
                                W.log_beliefs[local_offset + s + 1] = -Inf
                            end
                        end
                    end
                end
            end
        end
    end
end

