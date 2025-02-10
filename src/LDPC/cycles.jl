# Copyright (c) 2024 - 2025 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

#############################
      # getter functions
#############################

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

"""
    girth(C::AbstractLDPCCode; max_iter::Int = 100)

Return the girth of the Tanner graph of `C`.

An error is thrown if the maximum number of iterations is reached and
`-1` is returned to represent infinite girth.
"""
function girth(C::AbstractLDPCCode; max_iter::Int = 100)
    check_adj_list, var_adj_list = _node_adjacencies(C.H)
    girth_arr = zeros(Int, C.n)
    Threads.@threads for vn in 1:C.n
        iter = 0
        not_found = true
        to_check = [(0, [vn])]
        while not_found
            iter += 1
            to_check_next = Vector{Tuple{Int, Vector{Int}}}()
            for i in 1:length(to_check)
                for (prev, v_arr) in (to_check[i], )
                    for v in v_arr
                        for cn in var_adj_list[v]
                            if cn != prev
                                if iter != 1 && vn ∈ check_adj_list[cn]
                                    not_found = false
                                    girth_arr[vn] = iter + 1
                                    break
                                else
                                    push!(to_check_next, (cn, [v2 for v2 in check_adj_list[cn] if v2 != v]))
                                end
                            end
                        end
                        !not_found && break
                    end
                    !not_found && break
                end
                !not_found && break
            end
            !not_found && break
            iter += 1
            iter > max_iter && error("Hit the maximum number of iterations")
            isempty(to_check_next) && break
            to_check = to_check_next
        end
    end
    # println(girth_arr)
    min = minimum(girth_arr)
    iseven(min) || error("Computed girth to be an odd integer")
    min == 0 ? (C.girth = -1;) : (C.girth = min;)
    return C.girth
end

"""
    computation_graph(C::AbstractLDPCCode, lvl::Int, v::Int, v_type::Symbol = :v)

Return a figure representing the expansion of the Tanner graph of `C` to level `lvl`
for node `v`. If `v_type` is `:v`, `v` is interpreted as a variable node; otherwise,
`v_type` is `:c` and `v` is interpreted as a check node.

# Note
- Run `using Makie` to activate this extension.
"""
function computation_graph end

function remove_cycles(H::Union{Matrix{<: Integer}, CTMatrixTypes}, n_max::Int)
    n_max ≥ 4 || throw(DomainError("n_max must be an even number greater than four"))
    iseven(n_max) || (n_max -= 1;)
    size(H, 1) ≥ 1 && size(H, 2) ≥ 1 || throw(ArgumentError("Matrix cannot have a zero dimension"))

    # also used in Tanner.jl, perhaps should make into its own function
    # A = _adjacency_matrix_from_code(H)
    typeof(H) <: CTMatrixTypes ? (I = _Flint_matrix_to_Julia_int_matrix(H);) : (I = H;)
    nr_H, nc_H = size(I)
    A = vcat(hcat(zeros(Int, nc_H, nc_H), transpose(I)), hcat(I, zeros(Int, nr_H, nr_H)))
    # display(A)
    matrices = [A^i for i in 0:n_max - 1]
    # for M in matrices
    #     display(M)
    #     println(" ")
    # end
    nr, nc = size(A)
    for n in 4:2:n_max
        n2 = div(n, 2)
        n2m1 = n2 - 1
        n2m2 = n2 - 2
        nm1 = n - 1
        for i in 1:nr
            for j in 1:nc
                if i ≠ j
                    if matrices[n2 + 1][i, j] ≥ 2 && iszero(matrices[n2m2 + 1][i, j])
                        # println("$i, $j")
                        # println(matrices[n2 + 1][i, j])
                        temp = Int[]
                        for k in 1:nc
                            # should be exactly two of these
                            if matrices[n2m1 + 1][i, k] > 0 && matrices[1 + 1][j, k] == 1
                                append!(temp, k)
                            end
                        end
                        
                        if !isempty(temp)
                            k = rand(temp)
                            C_e = Int[]
                            for x in 1:nc
                                if iszero(matrices[nm1 + 1][j, x]) && iszero(matrices[nm1 + 1][k, x])
                                    append!(C_e, x)
                                end
                            end

                            E_e = Vector{Tuple{Int, Int}}()
                            for v1 in C_e
                                for v2 in C_e
                                    if matrices[1 + 1][v1, v2] > 0
                                        push!(E_e, (v1, v2))
                                    end
                                end
                            end

                            # this handles the case that C_e is empty
                            if !isempty(E_e)
                                (l, m) = rand(E_e)
                                # println("here")
                                # remove two old edges
                                matrices[1 + 1][j, k] -= 1
                                matrices[1 + 1][k, j] -= 1
                                matrices[1 + 1][l, m] -= 1
                                matrices[1 + 1][m, l] -= 1

                                # add two new edges
                                matrices[1 + 1][j, m] += 1
                                matrices[1 + 1][m, j] += 1
                                matrices[1 + 1][l, k] += 1
                                matrices[1 + 1][k, l] += 1

                                # update matrices in indices j, k, l, m
                                for indx in 2:n_max - 1
                                    for c in 1:nc
                                        matrices[indx + 1][j, c] = dot(view(matrices[indx - 1 + 1], j, :), view(matrices[1 + 1], :, c))
                                        matrices[indx + 1][k, c] = dot(view(matrices[indx - 1 + 1], k, :), view(matrices[1 + 1], :, c))
                                        matrices[indx + 1][l, c] = dot(view(matrices[indx - 1 + 1], l, :), view(matrices[1 + 1], :, c))
                                        matrices[indx + 1][m, c] = dot(view(matrices[indx - 1 + 1], m, :), view(matrices[1 + 1], :, c))
                                    end

                                    for r in 1:nc
                                        matrices[indx + 1][r, j] = dot(view(matrices[indx - 1 + 1], r, :), view(matrices[1 + 1], :, j))
                                        matrices[indx + 1][r, k] = dot(view(matrices[indx - 1 + 1], r, :), view(matrices[1 + 1], :, k))
                                        matrices[indx + 1][r, l] = dot(view(matrices[indx - 1 + 1], r, :), view(matrices[1 + 1], :, l))
                                        matrices[indx + 1][r, m] = dot(view(matrices[indx - 1 + 1], r, :), view(matrices[1 + 1], :, m))
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return matrices[1 + 1][nc_H + 1:end, 1:nc_H]
end
remove_cycles(L::LDPCCode, n_max::Int) = LDPCCode(remove_cycles(parity_check_matrix(L), n_max))
# remove_cycles(L::LDPCCode, n_max::Int) = remove_cycles(parity_check_matrix(L), n_max)

#############################
       # simple cycles
#############################

# modified from implementation in Graphs.jl
function _modified_hawick_james(g::DiGraph{Int}, num_var_nodes::Int,
    max_len::Int = 16)

    # since this is bipartite, the input parameter tells you which are var nodes and which are
    # check nodes
    nvg = Graphs.nv(g)
    local_cycles = [Vector{Vector{Int}}() for _ in 1:num_var_nodes]
    Threads.@threads for i in 1:num_var_nodes
        B = [Vector{Int}() for _ in Graphs.vertices(g)]
        blocked = zeros(Bool, nvg)
        stack = Vector{Int}()
        keys_Dict = Dict{Vector{Int}, Bool}()
        _circuit_recursive!(g, Graphs.vertices(g)[i], Graphs.vertices(g)[i], blocked, B, stack,
            local_cycles[i], keys_Dict, max_len)
    end
    return reduce(vcat, local_cycles)
end

function _circuit_recursive!(g::DiGraph{Int}, v1::Int, v2::Int, blocked::Vector{Bool},
    B::Vector{Vector{Int}}, stack::Vector{Int}, cycles::Vector{Vector{Int}},
    keys_Dict::Dict{Vector{Int}, Bool}, max_len::Int)

    if length(stack) + 1 <= max_len
        flag = false
        push!(stack, v2)
        blocked[v2] = true

        # just put this entire thing in the if statement
        Av = Graphs.outneighbors(g, v2)
        for w in Av
            (w < v1) && continue
            if w == v1
                if length(stack) > 2
                    cycle = copy(stack)
                    # println("cycle: $stack")
                    # println("also checking: ", [cycle[1]; reverse(cycle[2:end])])
                    if !haskey(keys_Dict, cycle) && !haskey(keys_Dict, [cycle[1]; reverse(
                        cycle[2:end])])

                        push!(cycles, cycle)
                        keys_Dict[cycle] = true
                        # println("is new")
                    # else
                    #     println("is old")
                    end
                end
                flag = true
            elseif !blocked[w]
                # bit-wise or
                flag |= _circuit_recursive!(g, v1, w, blocked, B, stack, cycles, keys_Dict, max_len)
            end
        end

        if flag
            _unblock!(v2, blocked, B)
        else
            for w in Av
                (w < v1) && continue
                if !(v2 in B[w])
                    push!(B[w], v2)
                end
            end
        end
        
        pop!(stack)
    else
        flag = true
    end
    return flag
end

function _unblock!(v::Int, blocked::Vector{Bool}, B::Vector{Vector{Int}})
    blocked[v] = false
    wPos = 1
    Bv = B[v]
    while wPos <= length(Bv)
        w = Bv[wPos]
        old_length = length(Bv)
        filter!(v -> v != w, Bv)
        wPos += 1 - (old_length - length(Bv))
        if blocked[w]
            _unblock!(w, blocked, B)
        end
    end
    return nothing
end

#############################
       # simple cycles
#############################

function _enumerate_cycles(L::AbstractLDPCCode, len::Int)
    g, _, _ = Tanner_graph(L)
    d = DiGraph(g)
    cycles = _modified_hawick_james(d, L.n, len)

    # no way to get around this if you want it parallel
    unique_cycles = Set{Vector{Int}}()
    final_cycles = Vector{Vector{Int}}()
    for cycle in cycles
        if length(cycle) > 2
            temp = sort(cycle)
            found_flag = false
            for key in keys(unique_cycles.dict)
                if temp == key
                    found_flag = true
                    break
                end
            end
            if !found_flag
                push!(unique_cycles, temp)
                push!(final_cycles, cycle)
            end
        end
    end
    
    if !isempty(final_cycles)
        L.max_cyc_len = len
        L.simple_cycles = final_cycles
        lens = length.(final_cycles)
        girth = minimum(lens)
        if ismissing(L.girth)
            L.girth = girth
        else
            if L.girth != girth
                @warn "Known girth, $(L.girth), does not match just computed girth, $girth"
            end
        end
    end
    return final_cycles
end

"""
    enumerate_simple_cycles(L::AbstractLDPCCode; len::Int = 16)

Return the unique simple cycles up to length `len` of the Tanner graph of `L`. If
`len` is `-1`, then all simple cycles will be enumerated. An empty `Vector{Vector{Int}}`
is returned when there is no cycles.

# Note
- Simple cycles do not contain the same vertex twice.
- Cycles are returned as a vector of vertex indices, where the vertices are ordered left-to-right by
  columns of `parity_check_matrix(L)` then top-to-bottom by rows.
"""
function enumerate_simple_cycles(L::AbstractLDPCCode; len::Int = 16)
    ispositive(len) || throw(DomainError("Cycle length parameter must be a positive integer"))

    # TODO add this variable to struct
    if len > L.max_cyc_len
        return _enumerate_cycles(L, len)
    elseif len == L.max_cyc_len
        return L.simple_cycles
    else
        return filter(x -> length(x) ≤ len, L.simple_cycles)
    end
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
simple_cycle_length_distribution(L::AbstractLDPCCode; len::Int = 16) = StatsBase.countmap(length.(
    enumerate_simple_cycles(L, len = len)))

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
average_simple_cycle_length(L::AbstractLDPCCode; len::Int = 16) = mean(length.(
    enumerate_simple_cycles(L, len = len)))

"""
    median_simple_cycle_length(L::AbstractLDPCCode; len::Int = 16)

Return the median cycle length of unique simple cycles up to length `len` of the Tanner graph of
`L`. If `len` is `-1`, then all simple cycles will be enumerated.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
median_simple_cycle_length(L::AbstractLDPCCode; len::Int = 16) = median(length.(
    enumerate_simple_cycles(L, len = len)))

"""
    mode_simple_cycle_length(L::AbstractLDPCCode; len::Int = 16)

Return the most common cycle length of unique simple cycles up to length `len` of the Tanner graph
of `L`. If `len` is `-1`, then all simple cycles will be enumerated.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
mode_simple_cycle_length(L::AbstractLDPCCode; len::Int = 16) = StatsBase.mode(length.(
    enumerate_simple_cycles(L, len = len)))

# TODO: there are ways to compute this without enumerating all of the cycles, but they are complicated and perhaps not worth the effort unless explicitly requested by a user
"""
    count_simple_cycles(L::AbstractLDPCCode; len::Int = 16)

Return the total number of unique simple cycles up to length `len` of the Tanner graph
of `L`. If `len` is `-1`, then all simple cycles will be enumerated.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
count_simple_cycles(L::AbstractLDPCCode; len::Int = 16) = length(enumerate_simple_cycles(L,
    len = len))

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
simple_cycle_distribution_by_variable_node(L::AbstractLDPCCode; len::Int = 16) =
    StatsBase.countmap(filter(x -> x ≤ L.n, reduce(vcat, enumerate_simple_cycles(L, len = len))))
# this works by removing the check nodes and then counting how many times each var node appears
# not the most efficient since it makes several lists repeatedly but fine for now

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

# TODO: time this versus calling girth and then doing the upper bound function
"""
    enumerate_short_cycles(L::AbstractLDPCCode; len::Int = 16)

Return the unique short cycles up to length `len` of the Tanner graph of `L`. If
`len` is `-1`, then all short cycles will be enumerated. An empty `Vector{Vector{Int}}`
is returned when there is no short cycles.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- Cycles are returned as a vector of vertex indices, where the vertices are ordered left-to-right by
  columns of `parity_check_matrix(L)` then top-to-bottom by rows.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
enumerate_short_cycles(L::AbstractLDPCCode; len::Int = 16) = filter(x -> length(x) <=
    2 * girth(L) - 2, _enumerate_cycles(L, len))

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
short_cycle_length_distribution(L::AbstractLDPCCode; len::Int = 16) = StatsBase.countmap(length.(
    enumerate_short_cycles(L, len = len)))

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
average_short_cycle_length(L::AbstractLDPCCode; len::Int = 16) = mean(length.(
    enumerate_short_cycles(L, len = len)))

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
median_short_cycle_length(L::AbstractLDPCCode; len::Int = 16) = median(length.(
    enumerate_short_cycles(L, len = len)))

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
mode_short_cycle_length(L::AbstractLDPCCode; len::Int = 16) = StatsBase.mode(length.(
    enumerate_short_cycles(L, len = len)))

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
count_short_cycles(L::AbstractLDPCCode; len::Int = 16) = length(enumerate_short_cycles(L,
    len = len))

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
short_cycle_distribution_by_variable_node(L::AbstractLDPCCode; len::Int = 16) =
    StatsBase.countmap(filter(x -> x ≤ L.n, reduce(vcat, enumerate_short_cycles(L, len = len))))

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
        # lollipops
#############################


# one for short cycles which cuts off stems at the correct lengths
# one for doing all of them





#############################
           # ACE
#############################

mutable struct _ACEVarNode
    id::Int
    parent_id::Int
    lvl::Int
    cum_ACE::Int
    local_ACE::Int
end

mutable struct _ACECheckNode
    id::Int
    parent_id::Int
    lvl::Int
    cum_ACE::Int
end

# TODO: degree 1 nodes
# why did I make this note? is ACE defined for them differently?
"""
    shortest_cycle_ACE(C::AbstractLDPCCode, v::Int)
    shortest_cycle_ACE(C::AbstractLDPCCode, vs::Vector{Int})
    shortest_cycle_ACE(C::AbstractLDPCCode)

Return a cycle of minimum length and minimum ACE in the Tanner graph of `C`
for the vertex `v` or vertices `vs`, in the order (ACEs, cycles). If no vertices
are given, all vertices are computed by default. The cycle `v1 -- c1 -- ... -- 
cn -- vn` is returned in the format `[(v1, c1), (c1, v2), ..., (cn, vn)]`.
"""
function shortest_cycle_ACE(C::AbstractLDPCCode, vs::Vector{Int})
    isempty(vs) && throw(ArgumentError("Input variable node list cannot be empty"))
    all(x -> 1 <= x <= C.n, vs) || throw(DomainError("Variable node indices must be between 1 and length(C)"))

    # might not be efficient to have the or here
    vs_to_do = [x for x in vs if isempty(C.ACEs_per_var_node[x]) || isempty(C.shortest_cycles[x])]
    processed = false
    if !isempty(vs_to_do)
        processed = true
        check_adj_list, var_adj_list = _node_adjacencies(C.H)
        
        Threads.@threads for i in 1:length(vs_to_do)
            # moving this inside allocates more but allows for multi-threading
            check_nodes = [_ACECheckNode(j, -1, -1, -1) for j in 1:length(check_adj_list)]
            var_nodes = [_ACEVarNode(j, -1, -1, -1, length(var_adj_list[j]) - 2) for j in 1:C.n]

            ACEs = Vector{Int}()
            cycle_lens = Vector{Int}()
            cycles = Vector{Vector{Tuple{Int, Int}}}()
            not_emptied = true
            root = var_nodes[vs_to_do[i]]
            root.lvl = 0
            root.cum_ACE = root.local_ACE
            queue = Deque{Union{_ACECheckNode, _ACEVarNode}}()
            push!(queue, root)
            while length(queue) > 0
                curr = first(queue)
                if isa(curr, _ACEVarNode)
                    for cn in var_adj_list[curr.id]
                        # can't pass messages back to the same node
                        if cn != curr.parent_id
                            cn_node = check_nodes[cn]
                            if cn_node.lvl != -1
                                # have seen before
                                push!(ACEs, curr.cum_ACE + cn_node.cum_ACE - root.local_ACE)
                                push!(cycle_lens, curr.lvl + cn_node.lvl + 1)

                                # trace the cycle from curr to root and cn_node to root
                                temp = Vector{Tuple{Int, Int}}()
                                node = cn_node
                                while node.lvl != 0
                                    push!(temp, (node.parent_id, node.id))
                                    if isodd(node.lvl)
                                        node = var_nodes[node.parent_id]
                                    else
                                        node = check_nodes[node.parent_id]
                                    end
                                end
                                reverse!(temp)
                                push!(temp, (cn_node.id, curr.id))
                                node = curr
                                while node.lvl != 0
                                    push!(temp, (node.id, node.parent_id))
                                    if isodd(node.lvl)
                                        node = var_nodes[node.parent_id]
                                    else
                                        node = check_nodes[node.parent_id]
                                    end
                                end
                                push!(cycles, temp)

                                # finish this level off but don't go deeper so remove children at lower level
                                if not_emptied
                                    while length(queue) > 0
                                        back = last(queue)
                                        if back.lvl != curr.lvl
                                            pop!(queue)
                                        else
                                            break
                                        end
                                    end
                                    not_emptied = false
                                end
                            elseif not_emptied
                                cn_node.lvl = curr.lvl + 1
                                cn_node.parent_id = curr.id
                                cn_node.cum_ACE = curr.cum_ACE
                                push!(queue, cn_node)
                            end
                        end
                    end
                else
                    for vn in check_adj_list[curr.id]
                        # can't pass messages back to the same node
                        if vn != curr.parent_id
                            vn_node = var_nodes[vn]
                            if vn_node.lvl != -1
                                # have seen before
                                push!(ACEs, curr.cum_ACE + vn_node.cum_ACE - root.local_ACE)
                                push!(cycle_lens, curr.lvl + vn_node.lvl + 1)

                                # trace the cycle from curr to root and cn_node to root
                                temp = Vector{Tuple{Int, Int}}()
                                node = vn_node
                                while node.lvl != 0
                                    push!(temp, (node.parent_id, node.id))
                                    if isodd(node.lvl)
                                        node = var_nodes[node.parent_id]
                                    else
                                        node = check_nodes[node.parent_id]
                                    end
                                end
                                reverse!(temp)
                                push!(temp, (vn_node.id, curr.id))
                                node = curr
                                while node.lvl != 0
                                    push!(temp, (node.id, node.parent_id))
                                    if isodd(node.lvl)
                                        node = var_nodes[node.parent_id]
                                    else
                                        node = check_nodes[node.parent_id]
                                    end
                                end
                                push!(cycles, temp)

                                # finish this level off but don't go deeper so remove children at lower level
                                if not_emptied
                                    while length(queue) > 0
                                        back = last(queue)
                                        if back.lvl != curr.lvl
                                            pop!(queue)
                                        else
                                            break
                                        end
                                    end
                                    not_emptied = false
                                end
                            elseif not_emptied
                                vn_node.lvl = curr.lvl + 1
                                vn_node.parent_id = curr.id
                                vn_node.cum_ACE = curr.cum_ACE + vn_node.local_ACE
                                push!(queue, vn_node)
                            end
                        end
                    end
                end
                popfirst!(queue)
            end
            # println("variable node $i, cycles: $cycles")
            C.ACEs_per_var_node[vs_to_do[i]] = ACEs
            C.cycle_lens[vs_to_do[i]] = cycle_lens
            C.shortest_cycles[vs_to_do[i]] = cycles
        end
    end

    vs_ACE = zeros(Int, length(vs))
    cycles_vs = [Vector{Tuple{Int, Int}}() for _ in 1:length(vs)]
    for i in 1:length(vs)
        min, index = findmin(C.ACEs_per_var_node[vs[i]])
        vs_ACE[i] = min
        cycles_vs[i] = C.shortest_cycles[vs[i]][index]
    end

    if processed
        if all(!isempty, C.cycle_lens)
            girth = minimum([minimum(C.cycle_lens[i]) for i in 1:C.n])
            if ismissing(C.girth)
                C.girth = girth
            else
                if C.girth != girth
                    @warn "Known girth, $(C.girth), does not match just computed girth, $girth"
                end
            end
        end
    end

    return vs_ACE, cycles_vs
end
shortest_cycle_ACE(C::AbstractLDPCCode, v::Int) = shortest_cycle_ACE(C, [v])[1]
shortest_cycle_ACE(C::AbstractLDPCCode) = shortest_cycle_ACE(C, collect(1:C.n))

"""
    shortest_cycles(C::AbstractLDPCCode, v::Int)
    shortest_cycles(C::AbstractLDPCCode, vs::Vector{Int})
    shortest_cycles(C::AbstractLDPCCode)

Return all the cycles of shortest length in the Tanner graph of `C` for the vertex `v` or
vertices `vs`. If no vertices are given, all vertices are computed by default.

# Note
- The length of the shortest cycle is not necessarily the same for each vertex.
- To reduce computational complexity, the same cycle may appear under each vertex in the cycle.
"""
function shortest_cycles(C::AbstractLDPCCode, vs::Vector{Int})
    shortest_cycle_ACE(C, vs)
    return C.shortest_cycles[vs]
end
shortest_cycles(C::AbstractLDPCCode, v::Int) = shortest_cycles(C, [v])[1]
shortest_cycles(C::AbstractLDPCCode) = shortest_cycles(C, collect(1:C.n))

"""
    ACE_distribution(C::AbstractLDPCCode, v::Int)
    ACE_distribution(C::AbstractLDPCCode, vs::Vector{Int})
    ACE_distribution(C::AbstractLDPCCode)

Return the ACEs and cycle lengths for vertex `v` or vertices `vs` of the Tanner graph
of `C`. If no vertices are given, all vertices are computed by default.
"""
function ACE_distribution(C::AbstractLDPCCode, vs::Vector{Int})
    # using the original DFS approach constructs a significantly larger tree than this truncated BFS approach

    isempty(vs) && throw(ArgumentError("Input node list cannot be empty"))
    all(x -> 1 <= x <= C.n, vs) || throw(DomainError("Variable node index must be between 1 and length(C)"))

    vs_to_do = [x for x in vs if isempty(C.ACEs_per_var_node[x])]
    processed = false
    if !isempty(vs_to_do)
        processed = true
        check_adj_list, var_adj_list = _node_adjacencies(C.H)
    
        Threads.@threads for i in 1:length(vs_to_do)
            # moving this inside allocates more but allows for multi-threading
            check_nodes = [_ACECheckNode(i, -1, -1, -1) for i in 1:length(check_adj_list)]
            var_nodes = [_ACEVarNode(i, -1, -1, -1, length(var_adj_list[i]) - 2) for i in 1:C.n]

            ACEs = Vector{Int}()
            cycle_lens = Vector{Int}()
            root = var_nodes[vs[i]]
            root.lvl = 0
            root.cum_ACE = root.local_ACE
            queue = Queue{Union{_ACECheckNode, _ACEVarNode}}()
            enqueue!(queue, root)
            while length(queue) > 0
                curr = first(queue)
                if isa(curr, _ACEVarNode)
                    for cn in var_adj_list[curr.id]
                        # can't pass messages back to the same node
                        if cn != curr.parent_id
                            cn_node = check_nodes[cn]
                            if cn_node.lvl != -1
                                # have seen before
                                push!(ACEs, curr.cum_ACE + cn_node.cum_ACE - root.local_ACE)
                                push!(cycle_lens, curr.lvl + cn_node.lvl + 1)
                            else
                                cn_node.lvl = curr.lvl + 1
                                cn_node.parent_id = curr.id
                                cn_node.cum_ACE = curr.cum_ACE
                                enqueue!(queue, cn_node)
                            end
                        end
                    end
                else
                    for vn in check_adj_list[curr.id]
                        # can't pass messages back to the same node
                        if vn != curr.parent_id
                            vn_node = var_nodes[vn]
                            if vn_node.lvl != -1
                                # have seen before
                                push!(ACEs, curr.cum_ACE + vn_node.cum_ACE - root.local_ACE)
                                push!(cycle_lens, curr.lvl + vn_node.lvl + 1)
                            else
                                vn_node.lvl = curr.lvl + 1
                                vn_node.parent_id = curr.id
                                vn_node.cum_ACE = curr.cum_ACE + vn_node.local_ACE
                                enqueue!(queue, vn_node)
                            end
                        end
                    end
                end
                dequeue!(queue)
            end
            C.ACEs_per_var_node[vs_to_do[i]] = ACEs
            # C.cycle_lens[vs_to_do[i]] = cycle_lens
        end
    end

    vs_ACEs = [C.ACEs_per_var_node[i] for i in vs]
    # lengths = [C.cycle_lens[i] for i in vs]

    # if processed
    #     if all(!isempty, C.cycle_lens)
    #         girth = minimum([minimum(C.cycle_lens[i]) for i in 1:C.n])
    #         if ismissing(C.girth)
    #             C.girth = girth
    #         else
    #             if C.girth != girth
    #                 @warn "Known girth, $(C.girth), does not match just computed girth, $girth"
    #             end
    #         end
    #     end
    # end

    return vs_ACEs #, lengths
end
function ACE_distribution(C::AbstractLDPCCode, v::Int)
    vs_ACE, lengths = ACE_distribution(C, [v])
    return vs_ACE[1], lengths[1]
end

# TODO: plots
ACE_distribution(C::AbstractLDPCCode) = ACE_distribution(C, collect(1:C.n))

"""
    average_ACE_distribution(C::AbstractLDPCCode, v::Int)
    average_ACE_distribution(C::AbstractLDPCCode, vs::Vector{Int})
    average_ACE_distribution(C::AbstractLDPCCode)

Return the average ACE of the vertex `v` or vertices `vs` of the Tanner graph of `C`. If no
vertices are given, all vertices are computed (individually) by default.
"""
function average_ACE_distribution(C::AbstractLDPCCode, vs::Vector{Int})
    vs_to_do = [x for x in vs if isempty(C.ACEs_per_var_node[x])]
    isempty(vs_to_do) || ACE_distribution(C, vs_to_do)
    return [mean(C.ACEs_per_var_node[v]) for v in vs]
end
average_ACE_distribution(C::AbstractLDPCCode, v::Int) = average_ACE_distribution(C, [v])[1]
average_ACE_distribution(C::AbstractLDPCCode) = average_ACE_distribution(C, collect(1:C.n))

"""
    median_ACE_distribution(C::AbstractLDPCCode, v::Int)
    median_ACE_distribution(C::AbstractLDPCCode, vs::Vector{Int})
    median_ACE_distribution(C::AbstractLDPCCode)

Return the median ACE of the vertex `v` or vertices `vs` of the Tanner graph of `C`. If no
vertices are given, all vertices are computed (individually) by default.
"""
function median_ACE_distribution(C::AbstractLDPCCode, vs::Vector{Int})
    vs_to_do = [x for x in vs if isempty(C.ACEs_per_var_node[x])]
    isempty(vs_to_do) || ACE_distribution(C, vs_to_do)
    return [median(C.ACEs_per_var_node[v]) for v in vs]
end
median_ACE_distribution(C::AbstractLDPCCode, v::Int) = median_ACE_distribution(C, [v])[1]
median_ACE_distribution(C::AbstractLDPCCode) = median_ACE_distribution(C, collect(1:C.n))

"""
    mode_ACE_distribution(C::AbstractLDPCCode, v::Int)
    mode_ACE_distribution(C::AbstractLDPCCode, vs::Vector{Int})
    mode_ACE_distribution(C::AbstractLDPCCode)

Return the mode ACE of the vertex `v` or vertices `vs` of the Tanner graph of `C`. If no
vertices are given, all vertices are computed (individually) by default.

# Note
- In case of ties, the smallest tied value is returned.
"""
function mode_ACE_distribution(C::AbstractLDPCCode, vs::Vector{Int})
    vs_to_do = [x for x in vs if isempty(C.ACEs_per_var_node[x])]
    isempty(vs_to_do) || ACE_distribution(C, vs_to_do)
    return [StatsBase.mode(sort(C.ACEs_per_var_node[v])) for v in vs]
end
mode_ACE_distribution(C::AbstractLDPCCode, v::Int) = mode_ACE_distribution(C, [v])[1]
mode_ACE_distribution(C::AbstractLDPCCode) = mode_ACE_distribution(C, collect(1:C.n))

"""
    ACE_spectrum(C::AbstractLDPCCode)

Return the ACE spectrum of the Tanner graph of `C`.
"""
function ACE_spectrum(C::AbstractLDPCCode)
    # vs_ACEs, lengths = ACE_distribution(C, collect(1:C.n))
    vs_ACEs = ACE_distribution(C, collect(1:C.n))
    # (false) spectrum: how many nodes have that ACE for that length
    # (true) spectrum: for a given length 4 <= l <= maximum(variabledegreedistribution(C)),
    # how many var nodes have shortest cycle that ACE

    shortest_lens = [minimum(i) for i in lengths]
    girth = minimum(shortest_lens)
    if ismissing(C.girth)
        C.girth = girth
    else
        if C.girth != girth
            @warn "Known girth, $(C.girth), does not match just computed girth, $girth"
        end
    end

    counts = [Dict{Int, Int}() for _ in 1:length(girth:2:2 * girth - 2)]
    for (k, l) in enumerate(girth:2:2 * girth - 2)
        for i in 1:length(shortest_lens)
            if shortest_lens[i] == l
                for j in 1:length(lengths[i])
                    if lengths[i][j] == l
                        if vs_ACEs[i][j] ∈ keys(counts[k])
                            counts[k][vs_ACEs[i][j]] += 1
                        else
                            counts[k][vs_ACEs[i][j]] = 1
                        end
                    end
                end
            end
        end
    end
    return counts
end

"""
    ACE_spectrum_plot(C::AbstractLDPCCode)

Return an interactive figure and data for the ACE spectrum of the Tanner graph of `C`.

# Note
- Run `using Makie` to activate this extension.
"""
function ACE_spectrum_plot end
