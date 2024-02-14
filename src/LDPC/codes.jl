# Copyright (c) 2022, 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    LDPCCode(H::fq_nmod_mat)

Return the LDPC code defined by the parity-check matrix `H`.
"""
function LDPCCode(H::CTMatrixTypes)
    # TODO: remove empties
    nnz, den = _density(H)
    # den <= 0.01 || (@warn "LDPC codes typically expect a density of less than 1%.";)

    cols, rows = _degree_distribution(H)
    is_reg = true
    c1 = cols[1]
    for i in 2:length(cols)
        c1 == cols[i] || (is_reg = false; break;)
    end
    if is_reg
        r1 = rows[1]
        for i in 2:length(rows)
            r1 == rows[i] || (is_reg = false; break;)
        end
    end
    c, r = maximum(cols), maximum(rows)

    R, x = PolynomialRing(Nemo.QQ, :x)
    col_poly = R(0)
    for i in cols
        col_poly += i * x^(i - 1)
    end
    col_poly = divexact(col_poly, nnz)
    row_poly = R(0)
    for i in rows
        row_poly += i * x^(i - 1)
    end
    row_poly = divexact(row_poly, nnz)

    C = LinearCode(H, true)
    return LDPCCode(base_ring(H), C.n, C.k, C.d, C.l_bound, C.u_bound, H, nnz,
        cols, rows, c, r, maximum([c, r]), den, is_reg, missing, col_poly,
        row_poly, missing, [Vector{Int}() for _ in 1:C.n], [Vector{Int}() for _ in 1:C.n],
        [Vector{Tuple{Int, Int}}() for _ in 1:C.n],
        Dict{Int, Int}(), Dict{Int, Int}())
end

"""
    LDPCCode(C::AbstractLinearCode)

Return the LDPC code given by the parity-check matrix of `C`.
"""
LDPCCode(C::AbstractLinearCode) = LDPCCode(parity_check_matrix(C))

"""
    regular_LDPC_code(q::Int, n::Int, l::Int, r::Int [; seed=nothing])

Return a random regular LDPC code over `GF(q)` of length `n` with column degree `l`
and row degree `r`.

If a seed is given, i.e. `regulardLDPCCode(4, 1200, 3, 6, seed=123)`, the
results are reproducible.
"""
function regular_LDPC_code(q::Int, n::Int, l::Int, r::Int; seed::Union{Nothing, Int} = nothing)
    Random.seed!(seed)
    m = divexact(n * l, r)
    F = if isprime(q)
        GF(q)
    else
        factors = Nemo.factor(q)
        length(factors) == 1 || throw(DomainError("There is no finite field of order $q"))
        (p, t), = factors
        GF(p, t, :α)
    end
    elems = collect(F)[2:end]
    H = zero_matrix(F, m, n)
    col_sums = zeros(Int, n)
    for i in axes(H, 1)
        ind = reduce(vcat, shuffle(filter(k -> col_sums[k] == s, 1:n)) for s in 0:l - 1)[1:r]
        for j in ind
            H[i, j] = rand(elems)
        end
        col_sums[ind] .+= 1
    end
    @assert all(count(.!iszero.(H[:, j])) == l for j in axes(H, 2))
    @assert all(count(.!iszero.(H[i, :])) == r for i in axes(H, 1))

    R, x = PolynomialRing(Nemo.QQ, :x)
    C = LinearCode(H, true)
    return LDPCCode(C.F, C.n, C.k, C.d, C.l_bound, C.u_bound, H, n * l, l * ones(Int, n),
        r * ones(Int, m), l, r, max(l, r), r / n, true, missing, (1 // l) * x^l,
        (1 // r) * x^r, missing, [Vector{Int}() for _ in 1:C.n], [Vector{Int}() for _ in 1:C.n],
        [Vector{Tuple{Int, Int}}() for _ in 1:C.n], Dict{Int, Int}(), Dict{Int, Int}())
end
regular_LDPC_code(n::Int, l::Int, r::Int; seed::Union{Nothing, Int} = nothing) =
    regular_LDPC_code(2, n, l, r, seed = seed)

#############################
      # getter functions
#############################

"""
    variable_degree_distribution(C::AbstractLDPCCode)

Return the variable node degree distribution of `C`.
"""
variable_degree_distribution(C::AbstractLDPCCode) = C.var_degs

"""
    check_degree_distribution(C::AbstractLDPCCode)

Return the check node degree distribution of `C`.
"""
check_degree_distribution(C::AbstractLDPCCode) = C.check_degs

"""
    degree_distributions(C::AbstractLDPCCode)

Return the variable and check node degree distributions of `C`.
"""
degree_distributions(C::AbstractLDPCCode) = C.var_degs, C.check_degs

"""
    column_bound(C::AbstractLDPCCode)

Return the column bound `c` of the `(c, r)`-LDPC code `C`.
"""
column_bound(C::AbstractLDPCCode) = C.col_bound

"""
    row_bound(C::AbstractLDPCCode)

Return the row bound `r` of the `(c, r)`-LDPC code `C`.
"""
row_bound(C::AbstractLDPCCode) = C.row_bound

"""
    column_row_bounds(C::AbstractLDPCCode)

Return the column and row bounds `c, r` of the `(c, r)`-LDPC code `C`.
"""
column_row_bounds(C::AbstractLDPCCode) = C.col_bound, C.row_bound

"""
    limited(C::AbstractLDPCCode)

Return the maximum of the row and column bounds for `C`.
"""
limited(C::AbstractLDPCCode) = C.limited

"""
    density(C::AbstractLDPCCode)

Return the density of the parity-check matrix of `C`.
"""
density(C::AbstractLDPCCode) = C.density

"""
    is_regular(C::AbstractLDPCCode)

Return `true` if the `C` is a regular LDPC code.

# Notes
- An LDPC is regular if all the column degrees and equal and all the row degrees are equal.
"""
is_regular(C::AbstractLDPCCode) = C.is_reg

"""
    variable_degree_polynomial(C::AbstractLDPCCode)

Return the variable degree polynomial of `C`.
"""
variable_degree_polynomial(C::AbstractLDPCCode) = C.λ

"""
    check_degree_polynomial(C::AbstractLDPCCode)

Return the check degree polynomial of `C`.
"""
check_degree_polynomial(C::AbstractLDPCCode) = C.ρ

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

function _degree_distribution(H::Union{CTMatrixTypes,
    MatElem{AbstractAlgebra.Generic.ResidueRingElem{fpPolyRingElem}}})

    nr, nc = size(H)
    cols = zeros(Int, 1, nc)
    @inbounds @views @simd for i in 1:nc
        # cols[i] = wt(H[:,  i])
        cols[i] = count(x -> !iszero(x), H[:, i])
    end
    rows = zeros(Int, 1, nr)
    @inbounds @views @simd for i in 1:nr
        # rows[i] = wt(H[i,  :])
        rows[i] = count(x -> !iszero(x), H[i, :])
    end
    return vec(cols), vec(rows)
end

function _density(H::CTMatrixTypes)
    count = 0
    nr, nc = size(H)
    for c in 1:nc
        for r in 1:nr
            !iszero(H[r, c]) && (count += 1;)
        end
    end
    return count, count / (nr * nc)
end

# TODO: make uniform with others
function show(io::IO, C::AbstractLDPCCode)
    if ismissing(C.d)
        if C.is_reg
            println(io, "[$(C.n), $(C.k)]_$(order(C.F)) regular ($(C.col_bound), $(C.row_bound))-LDPC code with density $(C.density).")
        else
            println(io, "[$(C.n), $(C.k)]_$(order(C.F)) irregular $(C.limited)-limited LDPC code with density $(C.density).")
        end
    else
        if C.is_reg
            println(io, "[$(C.n), $(C.k), $(C.d)]_$(order(C.F)) regular ($(C.col_bound), $(C.row_bound))-LDPC code with density $(C.density).")
        else
            println(io, "[$(C.n), $(C.k), $(C.d)]_$(order(C.F)) irregular $(C.limited)-limited LDPC code with density $(C.density).")
        end
    end
    if get(io, :compact, true)
        println(io, "\nVariable degree polynomial:")
        println(io, "\t", C.λ)
        println(io, "Check degree polynomial:")
        println(io, "\t", C.ρ)
        if C.n <= 30
            # was using Horig here, which is probably what I want
            nr, nc = size(C.H)
            println(io, "Parity-check matrix: $nr × $nc")
            for i in 1:nr
                print(io, "\t")
                for j in 1:nc
                    if j != nc
                        print(io, "$(C.H[i, j]) ")
                    elseif j == nc && i != nr
                        println(io, "$(C.H[i, j])")
                    else
                        print(io, "$(C.H[i, j])")
                    end
                end
            end
        # if !ismissing(C.weightenum)
        #     println(io, "\nComplete weight enumerator:")
        #     println(io, "\t", C.weightenum.polynomial)
        # end
        end
    end
end

"""
    degree_distributions_plot(C::AbstractLDPCCode)

Return a bar plot of the column and row degree distributions of `C`.
"""
function degree_distributions_plot(C::AbstractLDPCCode)
    cols, rows = degree_distributions(C)

    occurs_cols = [(i, count(==(i), cols)) for i in unique(cols)]
    cols_x_data = [x for (x, _) in occurs_cols]
    cols_y_data = [y for (_, y) in occurs_cols]
    cols_title="Variable Nodes"
    f1 = bar(cols_x_data, cols_y_data, bar_width=1, xticks=cols_x_data, yticks=cols_y_data,
        legend=false, xlabel="Degree", ylabel="Occurrences", title=cols_title)

    occurs_rows = [(i, count(==(i), rows)) for i in unique(rows)]
    rows_x_data = [x for (x, _) in occurs_rows]
    rows_y_data = [y for (_, y) in occurs_rows]
    rows_title="Check Nodes"
    f2 = bar(rows_x_data, rows_y_data, bar_width=1, xticks=rows_x_data, yticks=rows_y_data,
        legend=false, xlabel="Degree", ylabel="Occurrences", title=rows_title)
    f = Plots.plot(f1, f2, layout=(1, 2))
    display(f)
    return f
end

"""
    girth(C::LDPCCode)

Return the girth of the Tanner graph of `C`.

An error is thrown if the maximum number of iterations is reached and
``-1`` is returned to represent infinite girth.
"""
function girth(C::LDPCCode, max_iter::Int=100)
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
    shortest_cycle_ACE(C::LDPCCode, v::Int)
    shortest_cycle_ACE(C::LDPCCode, vs::Vector{Int})
    shortest_cycle_ACE(C::LDPCCode)

Return a cycle of minimum length and minimum ACE in the Tanner graph of `C`
for the vertex `v` or vertices `vs`, in the order (ACEs, cycles). If no vertices
are given, all vertices are computed by default. The cycle `v1 -- c1 -- ... -- 
cn -- vn` is returned in the format `[(v1, c1), (c1, v2), ..., (cn, vn)]`.
"""
function shortest_cycle_ACE(C::LDPCCode, vs::Vector{Int})
    isempty(vs) && throw(ArgumentError("Input variable node list cannot be empty"))
    all(x->1 <= x <= C.n, vs) || throw(DomainError("Variable node indices must be between 1 and length(C)"))

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
shortest_cycle_ACE(C::LDPCCode, v::Int) = shortest_cycle_ACE(C, [v])[1]
shortest_cycle_ACE(C::LDPCCode) = shortest_cycle_ACE(C, collect(1:C.n))

"""
    shortest_cycles(C::LDPCCode, v::Int)
    shortest_cycles(C::LDPCCode, vs::Vector{Int})
    shortest_cycles(C::LDPCCode)

Return all the cycles of shortest length in the Tanner graph of `C` for the vertex `v` or
vertices `vs`. If no vertices are given, all vertices are computed by default.

# Note
- The length of the shortest cycle is not necessarily the same for each vertex.
- To reduce computational complexity, the same cycle may appear under each vertex in the cycle.
"""
function shortest_cycles(C::LDPCCode, vs::Vector{Int})
    # display(vs)
    shortest_cycle_ACE(C, vs)
    return C.shortest_cycles[vs]
    # return [C.shortest_cycles[v] for v in vs]
    # cycles_vs = [Vector{Tuple{Int, Int}}() for _ in 1:length(vs)]
    # for (i, v) in enuemrate(vs)
    #     cycles_vs[i] = C.shortest_cycles[v]
    # end
    # isempty(vs) && throw(ArgumentError("Input variable node list cannot be empty"))
    # all(x->1 <= x <= C.n, vs) || throw(DomainError("Variable node indices must be between 1 and length(C)"))

    # check_adj_list, var_adj_list = _node_adjacencies(C.H)
    # cycles_vs = [Vector{Vector{Tuple{Int, Int}}}() for _ in 1:length(vs)]

    # Threads.@threads for i in 1:length(vs)
    #     # moving this inside allocates more but allows for multi-threading
    #     check_nodes = [_ACECheckNode(i, -1, -1, -1) for i in 1:length(check_adj_list)]
    #     var_nodes = [_ACEVarNode(i, -1, -1, -1, length(var_adj_list[i]) - 2) for i in 1:C.n]
    #     cycles = Vector{Vector{Tuple{Int, Int}}}()
    #     not_emptied = true
    #     root = var_nodes[vs[i]]
    #     root.lvl = 0
    #     queue = Deque{Union{_ACECheckNode, _ACEVarNode}}()
    #     push!(queue, root)
    #     while length(queue) > 0
    #         curr = first(queue)
    #         if isa(curr, _ACEVarNode)
    #             for cn in var_adj_list[curr.id]
    #                 # can't pass messages back to the same node
    #                 if cn != curr.parent_id
    #                     cn_node = check_nodes[cn]
    #                     if cn_node.lvl != -1
    #                         # have seen before
    #                         # trace the cycle from curr to root and cn_node to root
    #                         temp = Vector{Tuple{Int, Int}}()
    #                         node = cn_node
    #                         while node.lvl != 0
    #                             push!(temp, (node.parent_id, node.id))
    #                             if isodd(node.lvl)
    #                                 node = var_nodes[node.parent_id]
    #                             else
    #                                 node = check_nodes[node.parent_id]
    #                             end
    #                         end
    #                         reverse!(temp)
    #                         push!(temp, (cn_node.id, curr.id))
    #                         node = curr
    #                         while node.lvl != 0
    #                             push!(temp, (node.id, node.parent_id))
    #                             if isodd(node.lvl)
    #                                 node = var_nodes[node.parent_id]
    #                             else
    #                                 node = check_nodes[node.parent_id]
    #                             end
    #                         end
    #                         push!(cycles, temp)

    #                         # finish this level off but don't go deeper so remove children at lower level
    #                         if not_emptied
    #                             while length(queue) > 0
    #                                 back = last(queue)
    #                                 if back.lvl != curr.lvl
    #                                     pop!(queue)
    #                                 else
    #                                     break
    #                                 end
    #                             end
    #                             not_emptied = false
    #                         end
    #                     elseif not_emptied
    #                         cn_node.lvl = curr.lvl + 1
    #                         cn_node.parent_id = curr.id
    #                         push!(queue, cn_node)
    #                     end
    #                 end
    #             end
    #         else
    #             for vn in check_adj_list[curr.id]
    #                  # can't pass messages back to the same node
    #                 if vn != curr.parent_id
    #                     vn_node = var_nodes[vn]
    #                     if vn_node.lvl != -1
    #                         # have seen before
    #                         temp = Vector{Tuple{Int, Int}}()
    #                         node = vn_node
    #                         while node.lvl != 0
    #                             push!(temp, (node.parent_id, node.id))
    #                             if isodd(node.lvl)
    #                                 node = var_nodes[node.parent_id]
    #                             else
    #                                 node = check_nodes[node.parent_id]
    #                             end
    #                         end
    #                         reverse!(temp)
    #                         push!(temp, (vn_node.id, curr.id))
    #                         node = curr
    #                         while node.lvl != 0
    #                             push!(temp, (node.id, node.parent_id))
    #                             if isodd(node.lvl)
    #                                 node = var_nodes[node.parent_id]
    #                             else
    #                                 node = check_nodes[node.parent_id]
    #                             end
    #                         end
    #                         push!(cycles, temp)

    #                         # finish this level off but don't go deeper so remove children at lower level
    #                         if not_emptied
    #                             while length(queue) > 0
    #                                 back = last(queue)
    #                                 if back.lvl != curr.lvl
    #                                     pop!(queue)
    #                                 else
    #                                     break
    #                                 end
    #                             end
    #                             not_emptied = false
    #                         end
    #                     elseif not_emptied
    #                         vn_node.lvl = curr.lvl + 1
    #                         vn_node.parent_id = curr.id
    #                         push!(queue, vn_node)
    #                     end
    #                 end
    #             end
    #         end
    #         popfirst!(queue)
    #     end
    #     C.shortest_cycles[vs_to_do[i]] = cycles
    #     cycles_vs[i] = cycles
    # end
    # return cycles_vs
end
shortest_cycles(C::LDPCCode, v::Int) = shortest_cycles(C, [v])[1]
shortest_cycles(C::LDPCCode) = shortest_cycles(C, collect(1:C.n))
# function shortest_cycles(C::LDPCCode)
    # cycles = shortest_cycles(C, collect(1:C.n))
    # girth = minimum([minimum([length(cycle) for cycle in cycles[i]]) for i in 1:C.n])
    # if ismissing(C.girth)
    #     C.girth = girth
    # else
    #     if C.girth != girth
    #         @warn "Known girth, $(C.girth), does not match just computed girth, $girth"
    #     end
    # end
    # C.shortest_cycles = filter.(x -> length(x) < 2 * girth - 2, C.shortest_cycles)
    # return cycles
# end

function _progressive_node_adjacencies(H::CTMatrixTypes, vs::Vector{Int}, v_type::Symbol)
    check_adj_list, var_adj_list = _node_adjacencies(H)
    unique!(sort!(vs))
    len = length(vs)
    check_adj_lists = [deepcopy(check_adj_list) for _ in 1:len]
    var_adj_lists = [deepcopy(var_adj_list) for _ in 1:len]

    for i in 2:len
        prev = vs[1:i - 1]
        if v_type == :v
            for (j, x) in enumerate(check_adj_lists[i])
                check_adj_lists[i][j] = setdiff(x, prev)
            end
        else
            for (j, x) in enumerate(var_adj_lists[i])
                var_adj_lists[i][j] = setdiff(x, prev)
            end
        end
    end
    return check_adj_lists, var_adj_lists
end

function _count_cycles(C::LDPCCode)
    check_adj_lists, var_adj_lists = _progressive_node_adjacencies(C.H, collect(1:C.n), :v)
    lengths = [Vector{Int}() for _ in 1:C.n]
    Threads.@threads for i in 1:C.n
        check_nodes = [_ACECheckNode(i, -1, -1, -1) for i in 1:length(check_adj_lists[i])]
        var_nodes = [_ACEVarNode(i, -1, -1, -1, length(var_adj_lists[i][i]) - 2) for i in 1:C.n]

        cycle_lens = Vector{Int}()
        root = var_nodes[i]
        root.lvl = 0
        queue = Queue{Union{_ACECheckNode,_ACEVarNode}}()
        enqueue!(queue, root)
        while length(queue) > 0
            curr = first(queue)
            if isa(curr, _ACEVarNode)
                for cn in var_adj_lists[i][curr.id]
                    # can't pass messages back to the same node
                    if cn != curr.parent_id
                        cn_node = check_nodes[cn]
                        if cn_node.lvl != -1
                            # have seen before
                            push!(cycle_lens, curr.lvl + cn_node.lvl + 1)
                        else
                            cn_node.lvl = curr.lvl + 1
                            cn_node.parent_id = curr.id
                            enqueue!(queue, cn_node)
                        end
                    end
                end
            else
                for vn in check_adj_lists[i][curr.id]
                    # can't pass messages back to the same node
                    if vn != curr.parent_id
                        vn_node = var_nodes[vn]
                        if vn_node.lvl != -1
                            # have seen before
                            push!(cycle_lens, curr.lvl + vn_node.lvl + 1)
                        else
                            vn_node.lvl = curr.lvl + 1
                            vn_node.parent_id = curr.id
                            enqueue!(queue, vn_node)
                        end
                    end
                end
            end
            dequeue!(queue)
        end
        lengths[i] = cycle_lens
    end

    counts = Dict{Int, Int}()
    lens = unique!(reduce(vcat, lengths))
    for i in lens
        for j in 1:C.n
            if i ∈ keys(counts)
                counts[i] += count(x -> x == i, lengths[j])
            else
                counts[i] = count(x -> x == i, lengths[j])
            end
        end
    end
    C.elementary_cycle_counts = counts

    girth = minimum([isempty(lengths[i]) ? 9999999 : minimum(lengths[i]) for i in 1:C.n])
    girth == 9999999 && (girth = -1)
    if ismissing(C.girth)
        C.girth = girth
    else
        if C.girth != girth
            @warn "Known girth, $(C.girth), does not match just computed girth, $girth"
        end
    end

    counts = Dict{Int, Int}()
    for i in girth:2:2 * girth - 2
        for j in 1:C.n
            if i ∈ keys(counts)
                counts[i] += count(x -> x == i, lengths[j])
            else
                counts[i] = count(x -> x == i, lengths[j])
            end
        end
    end
    C.short_cycle_counts = counts
    return nothing
end

"""
    count_short_cycles(C::LDPCCode)

Return a bar graph and a dictionary of (length, count) pairs for unique short
cycles in the Tanner graph of `C`. An empty graph and dictionary are returned
when there are no cycles.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
"""
function count_short_cycles(C::LDPCCode)
    if isempty(C.short_cycle_counts) || isempty(C.elementary_cycle_counts)
        _count_cycles(C)
    end
    
    len = length(C.short_cycle_counts)
    x_data = [0 for _ in 1:len]
    y_data = [0 for _ in 1:len]
    index = 1
    for (i, j) in C.short_cycle_counts
        x_data[index] = i
        y_data[index] = j
        index += 1
    end

    fig = Plots.bar(x_data, y_data, bar_width=1, xticks=x_data, yticks=y_data,
        legend=false, xlabel="Cycle Length", ylabel="Occurrences", title="Short Cycle Counts")
    display(fig)
    return fig, C.short_cycle_counts
end

"""
    count_elementary_cycles(C::LDPCCode)

Return a bar graph and a dictionary of (length, count) pairs for unique elementary
cycles in the Tanner graph of `C`. An empty graph and dictionary are returned
when there are no cycles.

# Note
- Elementary cycles do not contain the same vertex twice and are unable to be
  decomposed into a sequence of shorter cycles.
"""
function count_elementary_cycles(C::LDPCCode)
    if isempty(C.short_cycle_counts) || isempty(C.elementary_cycle_counts)
        _count_cycles(C)
    end

    len = length(C.elementary_cycle_counts)
    x_data = [0 for _ in 1:len]
    y_data = [0 for _ in 1:len]
    index = 1
    for (i, j) in C.elementary_cycle_counts
        x_data[index] = i
        y_data[index] = j
        index += 1
    end

    fig = Plots.bar(x_data, y_data, bar_width=1, xticks=x_data, yticks=y_data,
        legend=false, xlabel="Cycle Length", ylabel="Occurrences", title="Elementary Cycle Counts")
    display(fig)
    return fig, C.elementary_cycle_counts
end

"""
    ACE_distribution(C::LDPCCode, v::Int)
    ACE_distribution(C::LDPCCode, vs::Vector{Int})
    ACE_distribution(C::LDPCCode)

Return the ACEs and cycle lengths for vertex `v` or vertices `vs` of the Tanner graph
of `C`. If no vertices are given, all vertices are computed by default.
"""
function ACE_distribution(C::LDPCCode, vs::Vector{Int})
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
            C.cycle_lens[vs_to_do[i]] = cycle_lens
        end
    end

    vs_ACEs = [C.ACEs_per_var_node[i] for i in vs]
    lengths = [C.cycle_lens[i] for i in vs]

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

    return vs_ACEs, lengths
end
function ACE_distribution(C::LDPCCode, v::Int)
    vs_ACE, lengths = ACE_distribution(C, [v])
    return vs_ACE[1], lengths[1]
end

# TODO: plots
ACE_distribution(C::LDPCCode) = ACE_distribution(C, collect(1:C.n))

"""
    average_ACE_distribution(C::LDPCCode, v::Int)
    average_ACE_distribution(C::LDPCCode, vs::Vector{Int})
    average_ACE_distribution(C::LDPCCode)

Return the average ACE of the vertex `v` or vertices `vs` of the Tanner graph of `C`. If no
vertices are given, all vertices are computed (individually) by default.
"""
function average_ACE_distribution(C::LDPCCode, vs::Vector{Int})
    vs_to_do = [x for x in vs if isempty(C.ACEs_per_var_node[x])]
    isempty(vs_to_do) || ACE_distribution(C, vs_to_do)
    return [mean(C.ACEs_per_var_node[v]) for v in vs]
end
average_ACE_distribution(C::LDPCCode, v::Int) = average_ACE_distribution(C, [v])[1]
average_ACE_distribution(C::LDPCCode) = average_ACE_distribution(C, collect(1:C.n))

"""
    median_ACE_distribution(C::LDPCCode, v::Int)
    median_ACE_distribution(C::LDPCCode, vs::Vector{Int})
    median_ACE_distribution(C::LDPCCode)

Return the median ACE of the vertex `v` or vertices `vs` of the Tanner graph of `C`. If no
vertices are given, all vertices are computed (individually) by default.
"""
function median_ACE_distribution(C::LDPCCode, vs::Vector{Int})
    vs_to_do = [x for x in vs if isempty(C.ACEs_per_var_node[x])]
    isempty(vs_to_do) || ACE_distribution(C, vs_to_do)
    return [median(C.ACEs_per_var_node[v]) for v in vs]
end
median_ACE_distribution(C::LDPCCode, v::Int) = median_ACE_distribution(C, [v])[1]
median_ACE_distribution(C::LDPCCode) = median_ACE_distribution(C, collect(1:C.n))

"""
    mode_ACE_distribution(C::LDPCCode, v::Int)
    mode_ACE_distribution(C::LDPCCode, vs::Vector{Int})
    mode_ACE_distribution(C::LDPCCode)

Return the mode ACE of the vertex `v` or vertices `vs` of the Tanner graph of `C`. If no
vertices are given, all vertices are computed (individually) by default.

# Note
- In case of ties, the smallest tied value is returned.
"""
function mode_ACE_distribution(C::LDPCCode, vs::Vector{Int})
    vs_to_do = [x for x in vs if isempty(C.ACEs_per_var_node[x])]
    isempty(vs_to_do) || ACE_distribution(C, vs_to_do)
    return [StatsBase.mode(sort(C.ACEs_per_var_node[v])) for v in vs]
end
mode_ACE_distribution(C::LDPCCode, v::Int) = mode_ACE_distribution(C, [v])[1]
mode_ACE_distribution(C::LDPCCode) = mode_ACE_distribution(C, collect(1:C.n))

"""
    ACE_spectrum(C::LDPCCode)

Return an interactive figure and data for the ACE spectrum of the Tanner graph of `C`.
"""
function ACE_spectrum(C::LDPCCode)
    vs_ACEs, lengths = ACE_distribution(C, collect(1:C.n))
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
    
    # TODO: remove WGLMakie as a default use and only use for interactive plots
    fig = Figure();
    ax = Axis(fig[1, 1], xlabel="ACE", ylabel="Occurrences", title="ACE Spectrum")
    sg = SliderGrid(fig[2, 1], (label = "Cycle Length", range=girth:2:2 * girth - 2, startvalue=4))

    x_max = maximum(reduce(vcat, vs_ACEs))
    y_max = 0
    counts = [Dict{Int, Int}() for _ in 1:length(girth:2:2 * girth - 2)]
    X_data = Observable(Vector{Int}())
    Y_data = Observable(Vector{Int}())
    barplot!(ax, X_data, Y_data, bar_width=1, xticks=X_data, yticks=Y_data)
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
        ys = collect(values(counts[k]))
        y_max = maximum([y_max; ys])

        on(sg.sliders[1].value) do val
            if to_value(val) == l
                X_data.val = collect(keys(counts[k]))
                Y_data.val = ys
                notify(X_data)
                notify(Y_data)
            end
        end
    end

    GLMakie.limits!(0, x_max + 1, 0, y_max + 2)
    display(fig)
    return fig, counts
end

mutable struct _ComputationGraphNode
    id::Int
    parent_id::Int
    lvl::Int
    vertex_number::Int
    type::Symbol
end

# doesn't seem to be a point in making this dynamic with a slider, as it simply
# continues in the same tree shape and no useful information is gained from watching it
"""
    computation_graph(C::LDPCCode, lvl::Int, v::Int, v_type::Symbol=:v)

Return a figure representing the expansion of the Tanner graph of `C` to level `lvl`
for node `v`. If `v_type` is `:v`, `v` is interpreted as a variable node; otherwise,
`v_type` is `:c` and `v` is interpreted as a check node.
"""
function computation_graph(C::LDPCCode, lvl::Int, v::Int, v_type::Symbol=:v)
    v_type ∈ (:v, :c) || throw(ArgumentError("Unknown argument for v_type"))
    if v_type == :v
        1 <= v <= C.n || throw(DomainError("Variable node index must be between 1 and length(C)"))
    else
        1 <= v <= nrows(C.H) || throw(DomainError("Check node index must be between 1 and the number of rows of H"))
    end
    lvl > 0 || throw(DomainError("Graph recipe requires at least one level"))

    check_adj_list, var_adj_list = _node_adjacencies(C.H)
    G = SimpleDiGraph()
    labels = Vector{String}()
    colors = Vector{Symbol}()
    markers = Vector{Symbol}()

    if v_type == :v
        root = _ComputationGraphNode(v, -1, 0, 1, :v)
        Grphs.add_vertex!(G)
        push!(labels, L"v_{%$v}")
        push!(colors, :black)
        push!(markers, :circle)
    else
        root = _ComputationGraphNode(v, -1, 0, 1, :c)
        Grphs.add_vertex!(G)
        push!(labels, L"c_{%$v}")
        push!(colors, :red)
        push!(markers, :rect)
    end
    
    queue = Queue{_ComputationGraphNode}()
    enqueue!(queue, root)
    while length(queue) > 0
        curr = first(queue)
        curr.lvl == lvl && break
        new_lvl = curr.lvl + 1
        if curr.type != :c
            for cn in var_adj_list[curr.id]
                if cn != curr.parent_id
                    Grphs.add_vertex!(G)
                    cn_node = _ComputationGraphNode(cn, curr.id, new_lvl, Grphs.nv(G), :c)
                    Grphs.add_edge!(G, curr.vertex_number, cn_node.vertex_number)
                    enqueue!(queue, cn_node)
                    push!(labels, L"c_{%$cn}")
                    push!(colors, :red)
                    push!(markers, :rect)
                end
            end
        else
            for vn in check_adj_list[curr.id]
                if vn != curr.parent_id
                    Grphs.add_vertex!(G)
                    vn_node = _ComputationGraphNode(vn, curr.id, new_lvl, Grphs.nv(G), :v)
                    Grphs.add_edge!(G, curr.vertex_number, vn_node.vertex_number)
                    enqueue!(queue, vn_node)
                    push!(labels, L"v_{%$vn}")
                    push!(colors, :black)
                    push!(markers, :circle)
                end
            end
        end
        dequeue!(queue)
    end

    # TODO: fix plot - count number of added nodes in level and manually pass in a calculated image size
    f, ax, p = graphplot(G, layout=Buchheim(),
        nlabels=labels,
        node_marker=markers,
        node_color=colors,
        nlabels_textsize=10,
        nlabels_align=(:left, :center),
        nlabels_distance=7);
    hidedecorations!(ax)
    hidespines!(ax)
    display(f)
    # TODO: what do we want to return here and make uniform with the doc string
    return f, ax, p
end
