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

"""
    bethe_region_graph(H::CTMatrixTypes) -> RegionGraph

The Bethe region graph: one outer region per parity check holding that check's
support, and one inner region per variable.

This is the region graph on which generalized belief propagation reduces exactly
to ordinary sum-product belief propagation, so it is the control that validates a
GBP implementation. Note that it is *not* what `canonical_region_graph` returns:
that takes the intersection closure of the check supports, whose inner regions
may hold several variables and which is therefore already a strictly stronger
approximation than Bethe.

Counting numbers are `1` for each check region and `1 - deg(v)` for each variable
region, which satisfies the validity conditions of `is_valid_region_graph`.
"""
function bethe_region_graph(H::CTMatrixTypes)
    num_check, num_var = size(H)

    supports = [Int[] for _ in 1:num_check]
    if H isa SparseMatrixCSC
        rows = rowvals(H)
        for v in 1:num_var, i in nzrange(H, v)
            push!(supports[rows[i]], v)
        end
    else
        for c in 1:num_check, v in 1:num_var
            iszero(H[c, v]) || push!(supports[c], v)
        end
    end
    for c in 1:num_check
        sort!(supports[c])
    end

    # Only variables that appear in some check get a region.
    live = [v for v in 1:num_var if any(v in supports[c] for c in 1:num_check)]
    var_region = Dict(v => num_check + k for (k, v) in enumerate(live))

    checks_of = Dict(v => Int[] for v in live)
    for c in 1:num_check, v in supports[c]
        push!(checks_of[v], c)
    end

    regs = Vector{Region}(undef, num_check + length(live))
    for c in 1:num_check
        kids = [var_region[v] for v in supports[c]]
        regs[c] = Region(copy(supports[c]), Int[], Int[], kids, 1)
    end
    for v in live
        ps = sort(checks_of[v])
        regs[var_region[v]] = Region([v], ps, copy(ps), Int[], 1 - length(ps))
    end
    return RegionGraph(regs)
end

#############################
            # GBP
#############################

# ==============================================================================
# GENERALIZED BELIEF PROPAGATION
# ==============================================================================
#
# The Yedidia-Freeman-Weiss parent-to-child algorithm on a region graph.
#
# Three conventions carry the correctness of this file, and getting any of them
# wrong silently degrades GBP into something weaker:
#
#   1. Sign. A positive channel LLR favours bit 0, matching `MP_decoders.jl`.
#
#   2. Local factors. `f_r` holds EVERY factor whose support lies inside `r`:
#      the channel term of each variable in `r`, and every parity check whose
#      support is contained in `r`. It is tempting to assign each factor to just
#      one region to avoid double counting, but that is wrong here -- the
#      counting numbers already prevent double counting in the region-based free
#      energy, and a region's belief is a *local joint* that needs every factor
#      it covers. This is also exactly what makes the Bethe region graph reduce
#      to ordinary sum-product BP.
#
#   3. Beliefs. Writing `D(r)` for the strict descendants of `r`,
#
#          b_r = f_r + sum of m_e over e in E(r),
#          E(r) = { (p -> c) : c in {r} u D(r), p not in {r} u D(r) },
#
#      so a region also receives the messages that other regions send into its
#      descendants. That term is what lets information travel between outer
#      regions: without it an outer region never hears from any other, its
#      belief is frozen at `f_r`, and the iteration does nothing.
#
# The message update enforces `b_c = sum_{x_p \ x_c} b_p` at a fixed point:
#
#     m_{p->c} <- marginalise(b_p onto x_c) - (b_c - m_{p->c}).

struct GBPWorkspace
    num_regions::Int
    num_edges::Int
    num_var::Int

    # Region state. `local_factors` and `log_beliefs` share `log_belief_offsets`.
    local_factors::Vector{Float64}
    log_beliefs::Vector{Float64}
    log_belief_offsets::Vector{Int}

    # Edge state. Messages run parent -> child.
    edge_to_parent::Vector{Int}
    edge_to_child::Vector{Int}
    messages::Vector{Float64}
    new_messages::Vector{Float64}
    message_offsets::Vector{Int}

    # `marg_maps[e][s_p + 1]` is the child state index for parent state `s_p`.
    marg_maps::Vector{Vector{Int}}

    # `belief_edges[r]` is E(r). `belief_bits[r][k]` gives the bit positions
    # inside `r` of the variables of edge `belief_edges[r][k]`'s child, in that
    # child's sorted variable order.
    belief_edges::Vector{Vector{Int}}
    belief_bits::Vector{Vector{Vector{Int}}}
    # Inverse of the above, so a serial schedule can apply one message delta
    # without recomputing every belief.
    edge_targets::Vector{Vector{Tuple{Int, Int}}}

    region_vars::Vector{Vector{Int}}
    counting::Vector{Int}
    var_regions::Vector{Vector{Int}}

    edge_update_order::Vector{Int}
    is_decimated::Vector{Bool}
    current_synd_buffer::Vector{UInt8}
    marg_buffer::Vector{Float64}
end

"""
Numerically stable log-sum-exp for combining probabilities in the log domain.
"""
@inline function _log_add_exp(x::Float64, y::Float64)
    x == -Inf && return y
    y == -Inf && return x
    m = max(x, y)
    return m + log1p(exp(-abs(x - y)))
end

# The max-log approximation.
@inline _log_add_exp_approx(x::Float64, y::Float64) = max(x, y)

# Strict descendants of every region, by transitive closure of `subregions`.
function _gbp_descendants(R::RegionGraph)
    n = length(R.regions)
    desc = [Set{Int}() for _ in 1:n]
    # Process shorter regions first so children are complete before parents.
    order = sortperm([length(R.regions[i].id) for i in 1:n])
    for i in order
        for c in R.regions[i].subregions
            push!(desc[i], c)
            union!(desc[i], desc[c])
        end
    end
    return desc
end

# Bit positions inside `outer` of each variable of `inner`, in `inner` order.
function _gbp_bit_positions(outer::Vector{Int}, inner::Vector{Int})
    pos = Vector{Int}(undef, length(inner))
    for (k, v) in enumerate(inner)
        idx = findfirst(==(v), outer)
        isnothing(idx) && error("region $inner is not contained in $outer")
        pos[k] = idx - 1
    end
    return pos
end

@inline function _gbp_project(s::Int, bits::Vector{Int})
    sc = 0
    @inbounds for i in eachindex(bits)
        sc |= (((s >> bits[i]) & 1) << (i - 1))
    end
    return sc
end

function init_gbp_workspace(R::RegionGraph, H::CTMatrixTypes)
    num_regions = length(R.regions)
    num_check, num_var = size(H)

    region_vars = [sort!(collect(R.regions[i].id)) for i in 1:num_regions]
    counting = [R.regions[i].overcounting_number for i in 1:num_regions]

    log_belief_offsets = zeros(Int, num_regions + 1)
    log_belief_offsets[1] = 1
    for i in 1:num_regions
        log_belief_offsets[i + 1] = log_belief_offsets[i] + 2^length(region_vars[i])
    end
    total_states = log_belief_offsets[end] - 1

    # Edges, one per (child, parent) pair.
    edge_to_parent = Int[]
    edge_to_child = Int[]
    message_offsets = Int[1]
    for c_idx in 1:num_regions
        for p_idx in R.regions[c_idx].parents
            push!(edge_to_parent, p_idx)
            push!(edge_to_child, c_idx)
            push!(message_offsets, message_offsets[end] + 2^length(region_vars[c_idx]))
        end
    end
    num_edges = length(edge_to_parent)

    # Parent-state -> child-state maps.
    marg_maps = Vector{Vector{Int}}(undef, num_edges)
    for e in 1:num_edges
        pv = region_vars[edge_to_parent[e]]
        cv = region_vars[edge_to_child[e]]
        bits = _gbp_bit_positions(pv, cv)
        map_e = Vector{Int}(undef, 2^length(pv))
        for s in 0:(2^length(pv) - 1)
            map_e[s + 1] = _gbp_project(s, bits) + 1
        end
        marg_maps[e] = map_e
    end

    # E(r) and the bit maps needed to read each of its messages.
    desc = _gbp_descendants(R)
    belief_edges = [Int[] for _ in 1:num_regions]
    belief_bits = [Vector{Vector{Int}}() for _ in 1:num_regions]
    edge_targets = [Tuple{Int, Int}[] for _ in 1:num_edges]
    for r in 1:num_regions
        inside = union(Set(r), desc[r])
        for e in 1:num_edges
            c = edge_to_child[e]
            p = edge_to_parent[e]
            (c in inside && !(p in inside)) || continue
            push!(belief_edges[r], e)
            push!(belief_bits[r], _gbp_bit_positions(region_vars[r], region_vars[c]))
            push!(edge_targets[e], (r, length(belief_edges[r])))
        end
    end

    var_regions = [Int[] for _ in 1:num_var]
    for r in 1:num_regions, v in region_vars[r]
        push!(var_regions[v], r)
    end

    # Two-way topological sweep: roots to leaves, then leaves to roots. Both
    # passes update the same parent -> child messages, which is correct here
    # because the beliefs now carry information in both directions.
    up = message_passing_order(R)
    order = Int[]
    sizehint!(order, 2 * num_edges)
    for r in reverse(up), e in 1:num_edges
        edge_to_parent[e] == r && push!(order, e)
    end
    for r in up, e in 1:num_edges
        edge_to_child[e] == r && push!(order, e)
    end

    max_child_states = num_edges == 0 ? 1 :
        maximum(2^length(region_vars[edge_to_child[e]]) for e in 1:num_edges)

    return GBPWorkspace(
        num_regions, num_edges, num_var,
        zeros(Float64, total_states), zeros(Float64, total_states), log_belief_offsets,
        edge_to_parent, edge_to_child,
        zeros(Float64, message_offsets[end] - 1), zeros(Float64, message_offsets[end] - 1),
        message_offsets,
        marg_maps,
        belief_edges, belief_bits, edge_targets,
        region_vars, counting, var_regions,
        order, zeros(Bool, num_var), zeros(UInt8, num_check),
        zeros(Float64, max_child_states)
    )
end

@inline function _gbp_normalize!(W::GBPWorkspace, r::Int)
    off = W.log_belief_offsets[r] - 1
    ns = W.log_belief_offsets[r + 1] - W.log_belief_offsets[r]
    m = -Inf
    @inbounds for s in 1:ns
        v = W.log_beliefs[off + s]
        v > m && (m = v)
    end
    m == -Inf && return
    @inbounds for s in 1:ns
        W.log_beliefs[off + s] -= m
    end
    return
end

# b_r = f_r + sum of messages in E(r).
function _gbp_recompute_beliefs!(W::GBPWorkspace)
    copyto!(W.log_beliefs, W.local_factors)
    @inbounds for r in 1:W.num_regions
        off = W.log_belief_offsets[r] - 1
        ns = W.log_belief_offsets[r + 1] - W.log_belief_offsets[r]
        es = W.belief_edges[r]
        for k in eachindex(es)
            moff = W.message_offsets[es[k]] - 1
            bits = W.belief_bits[r][k]
            for s in 0:(ns - 1)
                W.log_beliefs[off + s + 1] == -Inf && continue
                W.log_beliefs[off + s + 1] += W.messages[moff + _gbp_project(s, bits) + 1]
            end
        end
        _gbp_normalize!(W, r)
    end
    return
end

# Apply one message's change to every belief that includes it.
@inline function _gbp_apply_delta!(W::GBPWorkspace, e::Int, delta::Vector{Float64})
    @inbounds for (r, k) in W.edge_targets[e]
        off = W.log_belief_offsets[r] - 1
        ns = W.log_belief_offsets[r + 1] - W.log_belief_offsets[r]
        bits = W.belief_bits[r][k]
        for s in 0:(ns - 1)
            W.log_beliefs[off + s + 1] == -Inf && continue
            W.log_beliefs[off + s + 1] += delta[_gbp_project(s, bits) + 1]
        end
    end
    return
end

# One message: marginalise the parent belief and divide out the child's copy.
@inline function _gbp_message!(W::GBPWorkspace, e::Int, damping::Float64,
                               dest::Vector{Float64})
    p = W.edge_to_parent[e]
    c = W.edge_to_child[e]
    p_off = W.log_belief_offsets[p] - 1
    c_off = W.log_belief_offsets[c] - 1
    m_off = W.message_offsets[e] - 1
    n_c = W.message_offsets[e + 1] - W.message_offsets[e]
    n_p = W.log_belief_offsets[p + 1] - W.log_belief_offsets[p]

    buf = W.marg_buffer
    @inbounds for i in 1:n_c
        buf[i] = -Inf
    end
    map_e = W.marg_maps[e]
    @inbounds for s_p in 1:n_p
        v = W.log_beliefs[p_off + s_p]
        v == -Inf && continue
        sc = map_e[s_p]
        buf[sc] = _log_add_exp(buf[sc], v)
    end

    # Normalise the outgoing message so repeated sweeps cannot drift.
    mx = -Inf
    @inbounds for s_c in 1:n_c
        old = W.messages[m_off + s_c]
        cav = W.log_beliefs[c_off + s_c] - old
        raw = buf[s_c] == -Inf ? -Inf : (cav == -Inf ? buf[s_c] : buf[s_c] - cav)
        new = if raw == -Inf
            -Inf
        elseif old == -Inf
            raw
        else
            damping * raw + (1.0 - damping) * old
        end
        dest[s_c] = new
        new > mx && (mx = new)
    end
    if mx != -Inf
        @inbounds for s_c in 1:n_c
            dest[s_c] != -Inf && (dest[s_c] -= mx)
        end
    end
    return n_c
end

@inline function _gbp_region_llr(W::GBPWorkspace, r::Int, v::Int)
    off = W.log_belief_offsets[r] - 1
    ns = W.log_belief_offsets[r + 1] - W.log_belief_offsets[r]
    bit = findfirst(==(v), W.region_vars[r]) - 1
    l0 = -Inf
    l1 = -Inf
    @inbounds for s in 0:(ns - 1)
        b = W.log_beliefs[off + s + 1]
        b == -Inf && continue
        if ((s >> bit) & 1) == 0
            l0 = _log_add_exp(l0, b)
        else
            l1 = _log_add_exp(l1, b)
        end
    end
    return l0 - l1
end

"""
    gbp_marginal_llrs(W::GBPWorkspace; mode::Symbol = :smallest) -> Vector{Float64}

Single-variable LLRs, positive favouring bit 0.

At a consistent fixed point every region containing a variable agrees on its
marginal, so the two modes coincide. Off the fixed point they do not:

  * `:smallest` reads the innermost region containing the variable. These are the
    beliefs whose consistency the message updates actually enforce, and on the
    Bethe region graph this is exactly the ordinary BP posterior, which is what
    makes GBP reduce to BP there.
  * `:counting` takes the counting-number weighted combination over every region
    containing the variable. Valid, but it is not BP's posterior on Bethe.
"""
function gbp_marginal_llrs(W::GBPWorkspace; mode::Symbol = :smallest)
    llrs = zeros(Float64, W.num_var)
    if mode === :smallest
        @inbounds for v in 1:W.num_var
            rs = W.var_regions[v]
            isempty(rs) && continue
            best = rs[1]
            for r in rs
                length(W.region_vars[r]) < length(W.region_vars[best]) && (best = r)
            end
            llrs[v] = _gbp_region_llr(W, best, v)
        end
    elseif mode === :counting
        @inbounds for v in 1:W.num_var
            acc = 0.0
            for r in W.var_regions[v]
                # Clamp so a region that forbids one value outright cannot turn a
                # negative counting number into a sign flip.
                acc += W.counting[r] * clamp(_gbp_region_llr(W, r, v), -1.0e3, 1.0e3)
            end
            llrs[v] = acc
        end
    else
        throw(ArgumentError("mode must be :smallest or :counting, got :$mode"))
    end
    return llrs
end

"""
    extract_hard_decisions(W::GBPWorkspace, R::RegionGraph, num_var::Int)

Hard decisions from the counting-number weighted marginals.
"""
function extract_hard_decisions(W::GBPWorkspace, R::RegionGraph, num_var::Int)
    llrs = gbp_marginal_llrs(W)
    bits = zeros(UInt8, num_var)
    @inbounds for v in 1:min(num_var, W.num_var)
        bits[v] = llrs[v] < 0 ? 0x01 : 0x00
    end
    return bits
end

"""
Initializes the local factors, and hence the beliefs, of every region.

Each region receives every factor whose support it contains: the channel term of
each of its variables and each parity check lying wholly inside it. Messages are
reset to zero.
"""
function init_region_beliefs!(W::GBPWorkspace, R::RegionGraph, H::CTMatrixTypes,
                              channel_llrs::Vector{Float64},
                              target_syndrome::Vector{UInt8})
    num_check, num_var = size(H)

    check_adj = [Int[] for _ in 1:num_check]
    if H isa SparseMatrixCSC
        rows = rowvals(H)
        for v in 1:num_var, i in nzrange(H, v)
            push!(check_adj[rows[i]], v)
        end
    else
        for c in 1:num_check, v in 1:num_var
            iszero(H[c, v]) || push!(check_adj[c], v)
        end
    end

    fill!(W.local_factors, 0.0)
    fill!(W.messages, 0.0)
    fill!(W.new_messages, 0.0)

    @inbounds for r in 1:W.num_regions
        rv = W.region_vars[r]
        rset = Set(rv)
        ns = 2^length(rv)
        off = W.log_belief_offsets[r] - 1

        inside = [c for c in 1:num_check if issubset(check_adj[c], rset)]
        # Bit positions within the region for each contained check.
        cbits = [[findfirst(==(v), rv) - 1 for v in check_adj[c]] for c in inside]

        for s in 0:(ns - 1)
            bad = false
            for (ci, c) in enumerate(inside)
                par = 0x00
                for b in cbits[ci]
                    par ⊻= UInt8((s >> b) & 1)
                end
                if par != target_syndrome[c]
                    bad = true
                    break
                end
            end
            if bad
                W.local_factors[off + s + 1] = -Inf
                continue
            end
            acc = 0.0
            for (i, v) in enumerate(rv)
                ((s >> (i - 1)) & 1) == 1 && (acc -= channel_llrs[v])
            end
            W.local_factors[off + s + 1] = acc
        end
    end

    _gbp_recompute_beliefs!(W)
    return W
end

function _gbp_syndrome_ok(W::GBPWorkspace, H::CTMatrixTypes, bits::Vector{UInt8},
                          target_syndrome::Vector{UInt8})
    num_check, num_var = size(H)
    if H isa SparseMatrixCSC
        fill!(W.current_synd_buffer, 0x00)
        rows = rowvals(H)
        @inbounds for v in 1:num_var
            bits[v] == 0x01 || continue
            for i in nzrange(H, v)
                W.current_synd_buffer[rows[i]] ⊻= 0x01
            end
        end
        @inbounds for c in 1:num_check
            W.current_synd_buffer[c] == target_syndrome[c] || return false
        end
    else
        @inbounds for c in 1:num_check
            syn = target_syndrome[c]
            for v in 1:num_var
                iszero(H[c, v]) || (syn ⊻= bits[v])
            end
            syn == 0x00 || return false
        end
    end
    return true
end

"""
Master GBP API wrapper. Executes generalized belief propagation with early
termination on syndrome validity.

`schedule` is `:flooding` (all messages from the current beliefs, then one
belief rebuild) or `:serial` (each message applied before the next is computed).
Returns `(success, hard_decisions, iterations)`.
"""
function gbp_decode!(W::GBPWorkspace, R::RegionGraph, H::CTMatrixTypes,
                     total_llrs::Vector{Float64};
                     target_syndrome::Vector{UInt8} = zeros(UInt8, size(H, 1)),
                     decimation_type::Val = Val(:none),
                     dec_thresh::Float64 = 10.0,
                     dec_rounds::Int = 3,
                     max_iter::Int = 50,
                     damping::Float64 = 0.5,
                     schedule::Symbol = :flooding)
    num_check, num_var = size(H)
    fill!(W.is_decimated, false)
    init_region_beliefs!(W, R, H, total_llrs, target_syndrome)

    bits = extract_hard_decisions(W, R, num_var)
    _gbp_syndrome_ok(W, H, bits, target_syndrome) && return true, bits, 0

    delta = Float64[]
    scratch = Float64[]

    for iter in 1:max_iter
        if schedule === :serial
            for e in W.edge_update_order
                n_c = W.message_offsets[e + 1] - W.message_offsets[e]
                length(scratch) < n_c && resize!(scratch, n_c)
                length(delta) < n_c && resize!(delta, n_c)
                _gbp_message!(W, e, damping, scratch)
                m_off = W.message_offsets[e] - 1
                @inbounds for s in 1:n_c
                    old = W.messages[m_off + s]
                    new = scratch[s]
                    delta[s] = (new == -Inf || old == -Inf) ? 0.0 : new - old
                    W.messages[m_off + s] = new
                end
                _gbp_apply_delta!(W, e, delta)
            end
        else
            for e in W.edge_update_order
                n_c = W.message_offsets[e + 1] - W.message_offsets[e]
                length(scratch) < n_c && resize!(scratch, n_c)
                _gbp_message!(W, e, damping, scratch)
                m_off = W.message_offsets[e] - 1
                @inbounds for s in 1:n_c
                    W.new_messages[m_off + s] = scratch[s]
                end
            end
            copyto!(W.messages, W.new_messages)
            _gbp_recompute_beliefs!(W)
        end

        _apply_gbp_decimation!(decimation_type, W, R, num_var, iter, dec_thresh, dec_rounds)

        bits = extract_hard_decisions(W, R, num_var)
        _gbp_syndrome_ok(W, H, bits, target_syndrome) && return true, bits, iter
    end

    return false, extract_hard_decisions(W, R, num_var), max_iter
end

# ==============================================================================
# GBP DECIMATION HOOKS
# ==============================================================================

@inline _apply_gbp_decimation!(::Val{:none}, W, R, num_var, iter, thresh, rounds) = nothing

function _apply_gbp_decimation!(::Val{:hard}, W::GBPWorkspace, R::RegionGraph,
                                num_var::Int, iter::Int, thresh::Float64, rounds::Int)
    iter % rounds == 0 || return
    llrs = gbp_marginal_llrs(W)
    changed = false

    @inbounds for v in 1:num_var
        W.is_decimated[v] && continue
        abs(llrs[v]) > thresh || continue
        W.is_decimated[v] = true
        changed = true
        locked = llrs[v] > 0 ? 0 : 1
        # Banish contradictory states in every region holding this variable.
        # This edits the local factors, not the beliefs, so the constraint
        # survives the next belief rebuild.
        for r in W.var_regions[v]
            rv = W.region_vars[r]
            bit = findfirst(==(v), rv) - 1
            ns = 2^length(rv)
            off = W.log_belief_offsets[r] - 1
            for s in 0:(ns - 1)
                ((s >> bit) & 1) != locked && (W.local_factors[off + s + 1] = -Inf)
            end
        end
    end

    changed && _gbp_recompute_beliefs!(W)
    return
end
