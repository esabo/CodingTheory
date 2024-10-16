# Copyright (c) 2023 - 2024 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
     # general functions
#############################

"""
$TYPEDSIGNATURES

Return a bar plot of the column and row degree distributions of `C`.

# Note
- Run `using Makie` to activate this extension.
"""
function CodingTheory.degree_distributions_plot(C::AbstractLDPCCode)
    cols, rows = degree_distributions(C)

    occurs_cols = [(i, count(==(i), cols)) for i in unique(cols)]
    cols_x_data = [x for (x, _) in occurs_cols]
    cols_y_data = [y for (_, y) in occurs_cols]
    # f1 = bar(cols_x_data, cols_y_data, bar_width = 1, xticks = cols_x_data, yticks = cols_y_data,
    #     legend = false, xlabel = "Degree", ylabel = "Occurrences", title = "Variable Nodes")

    occurs_rows = [(i, count(==(i), rows)) for i in unique(rows)]
    rows_x_data = [x for (x, _) in occurs_rows]
    rows_y_data = [y for (_, y) in occurs_rows]
    # f2 = bar(rows_x_data, rows_y_data, bar_width = 1, xticks = rows_x_data, yticks = rows_y_data,
    #     legend = false, xlabel = "Degree", ylabel = "Occurrences", title = "Check Nodes")
    # f = Plots.plot(f1, f2, layout = (1, 2))

    fig = Figure()
    ax1 = Axis(fig[1, 1], xlabel = "Degree", ylabel = "Occurrences", title = "Variable Nodes")
    ax2 = Axis(fig[1, 1], xlabel = "Degree", ylabel = "Occurrences", title = "Check Nodes")
    barplot!(ax1, cols_x_data, cols_y_data, bar_width = 1, xticks = cols_x_data,
        yticks = cols_y_data)
    barplot!(ax2, rows_x_data, rows_y_data, bar_width = 1, xticks = rows_x_data,
        yticks = rows_y_data)
    display(f)
    return f
end

"""
$TYPEDSIGNATURES

Return a bar graph and a dictionary of (length, count) pairs for unique short
cycles in the Tanner graph of `C`. An empty graph and dictionary are returned
when there are no cycles.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- Run `using Makie` to activate this extension.
"""
function CodingTheory.count_short_cycles_plot(C::AbstractLDPCCode)
    if isempty(C.short_cycle_counts) || isempty(C.elementary_cycle_counts)
        CodingTheory._count_cycles(C)
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

    fig = Figure();
    ax = Axis(fig[1, 1], xlabel = "Cycle Length", ylabel = "Occurrences",
        title = "Short Cycle Counts")
    barplot!(ax, x_data, y_data, bar_width = 1, xticks = x_data, yticks = y_data)
    # fig = Plots.bar(x_data, y_data, bar_width = 1, xticks = x_data, yticks = y_data,
    #     legend = false, xlabel = "Cycle Length", ylabel = "Occurrences",
    #     title = "Short Cycle Counts")
    display(fig)
    return fig, C.short_cycle_counts
end

"""
$TYPEDSIGNATURES

Return a bar graph and a dictionary of (length, count) pairs for unique elementary
cycles in the Tanner graph of `C`. An empty graph and dictionary are returned
when there are no cycles.

# Note
- Elementary cycles do not contain the same vertex twice and are unable to be
  decomposed into a sequence of shorter cycles.
- Run `using Makie` to activate this extension.
"""
function CodingTheory.count_elementary_cycles_plot(C::AbstractLDPCCode)
    if isempty(C.short_cycle_counts) || isempty(C.elementary_cycle_counts)
        CodingTheory._count_cycles(C)
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

    fig = Figure();
    ax = Axis(fig[1, 1], xlabel = "Cycle Length", ylabel = "Occurrences",
        title = "Elementary Cycle Counts")
    barplot!(ax, x_data, y_data, bar_width = 1, xticks = x_data, yticks = y_data)
    # fig = Plots.bar(x_data, y_data, bar_width = 1, xticks = x_data, yticks = y_data,
    #     legend = false, xlabel = "Cycle Length", ylabel = "Occurrences",
    #     title = "Elementary Cycle Counts")
    display(fig)
    return fig, C.elementary_cycle_counts
end

"""
$TYPEDSIGNATURES

Return an interactive figure and data for the ACE spectrum of the Tanner graph of `C`.

# Note
- Run `using Makie` to activate this extension.
"""
function CodingTheory.ACE_spectrum_plot(C::AbstractLDPCCode) 
    # TODO: remove WGLMakie as a default use and only use for interactive plots
    fig = Figure();
    ax = Axis(fig[1, 1], xlabel = "ACE", ylabel = "Occurrences", title = "ACE Spectrum")
    sg = SliderGrid(fig[2, 1], (label = "Cycle Length", range = girth:2:2 * girth - 2,
        startvalue = 4))

    x_max = maximum(reduce(vcat, vs_ACEs))
    y_max = 0
    counts = ACE_spectrum(C)
    X_data = Observable(Vector{Int}())
    Y_data = Observable(Vector{Int}())
    barplot!(ax, X_data, Y_data, bar_width = 1, xticks = X_data, yticks = Y_data)
    for (k, l) in enumerate(girth:2:2 * girth - 2)
        
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
$TYPEDSIGNATURES

Return a figure representing the expansion of the Tanner graph of `C` to level `lvl`
for node `v`. If `v_type` is `:v`, `v` is interpreted as a variable node; otherwise,
`v_type` is `:c` and `v` is interpreted as a check node.

# Note
- Run `using Makie` to activate this extension.
"""
function CodingTheory.computation_graph(C::AbstractLDPCCode, lvl::Int, v::Int, v_type::Symbol = :v)
    v_type ∈ (:v, :c) || throw(ArgumentError("Unknown argument for v_type"))
    if v_type == :v
        1 <= v <= C.n || throw(DomainError("Variable node index must be between 1 and length(C)"))
    else
        1 <= v <= nrows(C.H) || throw(DomainError("Check node index must be between 1 and the number of rows of H"))
    end
    lvl > 0 || throw(DomainError("Graph recipe requires at least one level"))

    check_adj_list, var_adj_list = CodingTheory._node_adjacencies(C.H)
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
    f, ax, p = graphplot(G, layout = Buchheim(),
        nlabels = labels,
        node_marker = markers,
        node_color = colors,
        nlabels_textsize = 10,
        nlabels_align = (:left, :center),
        nlabels_distance = 7);
    hidedecorations!(ax)
    hidespines!(ax)
    display(f)
    # TODO: what do we want to return here and make uniform with the doc string
    return f, ax, p
end
