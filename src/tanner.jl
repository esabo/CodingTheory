# Copyright (c) 2022, Eric Sabo, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
    Tannergraph(H::fq_nmod_mat)

Return the Tanner graph of the matrix `H` as a `Figure` object.
"""
function Tannergraph(H::fq_nmod_mat)
    # convert H to A
    M = FpmattoJulia(H)
    nr, nc = size(M)
    A = zeros(Int, nr + nc, nr + nc)
    # bottom left corner
    # need to threshold any nonzero to a 1
    # A[nc + 1:end, 1:nc] = M
    # no, put in top right corner in order to get parents, childs working
    A[1:nc, nc + 1:end] = transpose(M)

    f = Figure();
    ax = Axis(f[1, 1], yreversed = true, xautolimitmargin = (0.15, 0.20),
        yautolimitmargin = (0.15, 0.20))
    hidespines!(ax)
    hidedecorations!(ax)

    leftx, lefty = zeros(nc), 1.:nc
    rightx, righty = ones(nr) * nr, range(1, nc, nr)
    x = vcat(leftx, rightx)
    y = vcat(lefty, righty)
    points = Point.(zip(x, y))
    cols = (:aqua, :red, :orange, :green, :blue, :purple)

    G = SimpleDiGraph(A)
    parents = [inneighbors(G, i) for i in Graphs.vertices(G)]
    childs = findall(x -> length(x) > 0, parents)
    # println(parents)
    # println(childs)

    for (i, v) in enumerate(childs)
        for node in parents[v]
            lines!(Point(x[[node, v]]...), Point(y[[node, v]]...),
                   color=cols[i % 6 + 1], linewidth=5)
        end
        text!(points[v], text=L"h_%$i", offset=(20, -15))
    end

    for (i, point) in enumerate(points[1:nc])
        CairoMakie.scatter!(point, color=:black, marker=:circle, markersize=25)
        text!(point, text=L"v_%$i", offset=(-30, -10))
    end

    for (i, point) in enumerate(points[nc + 1:end])
        CairoMakie.scatter!(point, color=:black, marker=:rect, markersize=25)
    end
    f
    return f
    # save("test.png", f)
end

"""
    Tannercode(G::Union{Graph, SparseMatrixCSC{Int, Int}}, C::AbstractLinearCode)

Return the Tanner code constructed from the graph defined by `G` and `C`.

An argument of type `SparseMatrixCSC` is assumed to be a valid adjacency matrix of a graph.
"""
function Tannercode(G::Union{SimpleGraph{Int}, SparseMatrixCSC{Int, Int}}, C::AbstractLinearCode)    
    typeof(G) <: SimpleGraph{Int} ? (A = adjacency_matrix(G);) : (A = G;)
    C.n == length(A[1, :].nzind) || throw(ArgumentError("Local code length must equal graph valency."))

    H = paritycheckmatrix(C)
    nrH = nrows(H)
    nrA = A.m # number rows
    ncA = A.n # number columns
    m = nrA * (C.n - C.k)

    F, _ = FiniteField(2, 1, "α")
    Fzero = F(0)
    Fone = F(1)
    pcm = zero_matrix(F, m, ncA)
    l = 1
    for i in 1:nrA
        # allocates
        edges = A[i, :].nzind
        for r in 1:nrH
            check = vec(H[r, :] .!= Fzero)
            supp = getindex(edges, check)
            for idx in supp
                pcm[l, idx] = Fone
            end
            l += 1
        end
    end
    return LinearCode(pcm, true)
end
