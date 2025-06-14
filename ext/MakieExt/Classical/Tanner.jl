# Copyright (c) 2022 - 2024 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
$TYPEDSIGNATURES

Return the Tanner graph of the matrix `H` as a `Figure` object.

# Note
- Run `using Makie` to activate this extension.
"""
function CodingTheory.Tanner_graph_plot(H::Union{CodingTheory.CTMatrixTypes,Matrix{Int}})
    # convert H to A
    M = CodingTheory._Flint_matrix_to_Julia_int_matrix(H)
    nr, nc = size(M)
    A = zeros(Int, nr + nc, nr + nc)
    # put in top right corner in order to get parents, children working
    A[1:nc, (nc+1):end] = transpose(M)

    fig = CairoMakie.Figure();
    ax = CairoMakie.Axis(
        fig[1, 1],
        yreversed = true,
        xautolimitmargin = (0.15, 0.20),
        yautolimitmargin = (0.15, 0.20),
    )
    CairoMakie.hidespines!(ax)
    CairoMakie.hidedecorations!(ax)

    left_x, left_y = zeros(nc), 1.0:nc
    right_x, right_y = ones(nr) * nr, range(1, nc, nr)
    x = vcat(left_x, right_x)
    y = vcat(left_y, right_y)
    points = CairoMakie.Point2f.(zip(x, y))
    cols = (:aqua, :red, :orange, :green, :blue, :purple)

    G = Grphs.SimpleDiGraph(A)
    parents = [Grphs.inneighbors(G, i) for i in Grphs.vertices(G)]
    children = findall(x -> length(x) > 0, parents)

    for (i, v) in enumerate(children)
        for node in parents[v]
            CairoMakie.lines!(
                [CairoMakie.Point2f(x[[node, v]])...],
                [CairoMakie.Point2f(y[[node, v]])...],
                color = cols[i%6+1],
                linewidth = 5,
            )
        end
        CairoMakie.text!(points[v], text = L"c_{%$i}", offset = (20, -15))
    end

    for (i, point) in enumerate(points[1:nc])
        CairoMakie.scatter!(point, color = :black, marker = :circle, markersize = 25)
        CairoMakie.text!(point, text = L"v_{%$i}", offset = (-30, -10))
    end

    for (i, point) in enumerate(points[(nc+1):end])
        CairoMakie.scatter!(point, color = :black, marker = :rect, markersize = 25)
    end
    display(fig)
    return fig
end
