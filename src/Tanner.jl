# Copyright (c) 2022, 2023 Eric Sabo, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
    Tannergraphplot(H::Union{fq_nmod_mat, Matrix{Int}})

Return the Tanner graph of the matrix `H` as a `Figure` object.
"""
function Tannergraphplot(H::Union{fq_nmod_mat, Matrix{Int}})
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
    Tannergraph(H::Union{fq_nmod_mat, Matrix{Int}})

Return the `SimpleGraph` object repesenting the Tanner graph of the parity-check
matrix `H` along with the indices of the left and right vertices representing
the bits and parity-checks, respectively.
"""
function Tannergraph(H::Union{fq_nmod_mat, Matrix{Int}})
    typeof(H) <: fq_nmod_mat ? (I = FpmattoJulia(H);) : (I = H;)
    Itr = transpose(I)
    B = vcat(hcat(zeros(I), Itr), hcat(I, zeros(Itr)))
    G = SimpleGraph(B)
    nr, nc = size(H)
    # lhs - bits
    # rhs - parity-checks
    return G, collect(1:nc), collect(nc + 1:nc + nr)
end

"""
    Tannergraph(C::AbstractLinearCode)

Return the `SimpleGraph` object repesenting the Tanner graph of `C` along with
the indices of the left and right vertices representing the bits and parity-checks,
respectively.
"""
Tannergraph(C::AbstractLinearCode) = Tannergraph(paritycheckmatrix(C))

# compressed sparse column (CSC) format used here so data is
# colptr, nzvals, rowval
# nzvals - stores all the nonzero values of the matrix
# rowval - for every nonzero value stores the rows at which they occur
# colptr - stores at which element in nzvals does the next column start
# ex: [a 0 b; c d e; 0 0 f] gives
# nzvals = [a c d b e f]
# rowval = [1 2 2 1 2 3]
# colptr = [1 3 4]
# additionally, it is common to put the number of nonzeros + 1 at the end of colptr
# colptr = [1 3 4 7]

# TODO: branch for small and large outputs
# TODO: multi-thread the outer for loop
# TODO: check if should not sure value in Hlocind and should instead just access Hloc directly
# TODO: make checks for a binary code and if so just loop and set to 1
"""
    Tannercode(EVI::SparseMatrixCSC{Int, Int}, C::AbstractLinearCode)

Return the Tanner code obtained by applying the local code `C` to the edges of the graph with
edge-vertex incidence matrix `EVI`.
"""
function Tannercode(EVI::SparseMatrixCSC{Int, Int}, C::AbstractLinearCode)
    numE = EVI.m # rows
    numV = EVI.n # columns
    numE > numV || throw(ArgumentError("The number of rows (edges) must be larger than the number of columns (vertices)."))
    nnz(EVI) % numE == 0 || throw(ArgumentError("The number of vertices does not divide the number of non-zero entries, cannot be regular."))
    nnz(EVI) % (C.n - C.k) == 0 || throw(ArgumentError("The dimension of the local code does not divide the number of non-zero entries."))
    
    # pre-store all the information about Hloc
    # this could be a bit redundant if Hloc is sparse
    Hloc = FpmattoJulia(paritycheckmatrix(C))
    nrHloc, ncHloc = size(Hloc)
    Hlocind = Vector{Vector{Tuple{Int, Int}}}()
    for r in 1:nrHloc
        temp = Vector{Tuple{Int, Int}}()
        for c in 1:ncHloc
            Hloc[r, c] != 0 && (push!(temp, (Hloc[r, c], c)))
        end
        push!(Hlocind, temp)
    end
    Hrowsitrloc = ones(Int, 1, nrHloc)
    Hlocindlens = [length(Hlocind[i]) for i in 1:nrHloc]

    currrow = 0
    H = zeros(Int, numV * nrHloc, numE)
    # look at every edge attached to a vertex, so check every row for a fixed column
    for c in 1:numV
        count = 0
        # since these graphs are regular this is always the same gap and could be exploited to save a few clock cycles
        for r in EVI.colptr[c]:EVI.colptr[c + 1] - 1
            count += 1
            # this is the count-th edge, is the count-th entry of any row of Hloc
            # this loop handles all rows of Hloc in a single pass
            # doesn't actually do anything here because there's the if statement I think
            @simd for i in 1:nrHloc
                # instead of looping through a col of Hloc every time
                # Hlocind stores the next column index and since this is sorted
                # if the 2nd element of Hlocind at this index is count then
                # there is a 1 there, otherwise a zero since it wasn't stored
                if Hlocindlens[i] >= Hrowsitrloc[i] && Hlocind[i][Hrowsitrloc[i]][2] == count
                    H[currrow + i, EVI.rowval[r]] = Hlocind[i][Hrowsitrloc[i]][1]
                    Hrowsitrloc[i] += 1
                end
            end
        end
        currrow += nrHloc
        Hrowsitrloc[:] .= 1
    end
    return H
end 

"""
    Tannercode(G::SimpleGraph{Int}, C::AbstractLinearCode)

Return the Tanner code obtained by applying the local code `C` to the edges of `G`.
"""
function Tannercode(G::SimpleGraph{Int}, C::AbstractLinearCode)
    isregular(G) || throw(ArgumentError("Graph must be regular."))
    length(G.fadjlist[1]) == C.n || throw(ArgumentError("The degree of the verties must be equal to the length of the local code."))
    # can use G.fadjlist directly?
    # would need to make sure we don't use the same edge twice
    return Tannercode(sparse(transpose(incidence_matrix(G))), C)
end

"""
    Tannercode(G::SimpleGraph{Int}, left::Vector{Int}, right::Vector{Int}, C::AbstractLinearCode)

Return the Tanner code obtained by applying the local code `C` to the vertices `right` in the
bipartition of `G` and treating the vertices of `left` as bits.
"""
function Tannercode(G::SimpleGraph{Int}, left::Vector{Int}, right::Vector{Int}, C::AbstractLinearCode)
    # remove this to allow for overcomplete matrices like quasi-cyclic codes?
    # length(left) > length(right) || throw(ArgumentError("The size of `left` (bits) must be greater than the size of `right` (parity checks)."))
    isvalidbipartition(G, left, right) || throw(ArgumentError("The input vectors are not a valid partition for the graph."))
    
    # pre-store all the information about Hloc
    # this could be a bit redundant if Hloc is sparse
    Hloc = FpmattoJulia(paritycheckmatrix(C))
    nrHloc, ncHloc = size(Hloc)
    Hlocind = Vector{Vector{Tuple{Int, Int}}}()
    for r in 1:nrHloc
        temp = Vector{Tuple{Int, Int}}()
        for c in 1:ncHloc
            Hloc[r, c] != 0 && (push!(temp, (Hloc[r, c], c)))
        end
        push!(Hlocind, temp)
    end

    # make dictionary here mapping left to 1:|E|
    # this should use sizehint now so no longer slower than the standard loop
    edgemap = Dict(lv => i for (i, lv) in enumerate(left))
    currrow = 0
    H = zeros(Int, length(right) * nrHloc, length(left))
    for rv in right
        for r in 1:nrHloc
            @simd for c in Hlocind[r]
                H[currrow + r, edgemap[G.fadjlist[rv][c[2]]]] = c[1]
            end
        end
        currrow += nrHloc
    end
    return H
end

# TODO: currently untested
"""
    Tannercode(G::SimpleGraph{Int}, left::Vector{Int}, right1::Vector{Int}, right2::Vector{Int}, C1::AbstractLinearCode, C2::AbstractLinearCode)

Return the Tanner code obtained by applying the local code `C1` to the vertices `right1` and the local code `C2` to the vertices
`right2` in the bipartition of `G` and treating the vertices of `left` as bits.
"""
function Tannercode(G::SimpleGraph{Int}, left::Vector{Int}, right1::Vector{Int}, right2::Vector{Int}, C1::AbstractLinearCode, C2::AbstractLinearCode)
    # remove this to allow for overcomplete matrices like quasi-cyclic codes?
    # length(left) > length(right) || throw(ArgumentError("The size of `left` (bits) must be greater than the size of `right` (parity checks)."))
    isvalidbipartition(G, left, right1 ∪ right2) || throw(ArgumentError("The input vectors are not a valid partition for the graph."))
    
    # pre-store all the information about Hloc
    # this could be a bit redundant if Hloc is sparse
    Hloc = FpmattoJulia(paritycheckmatrix(C1))
    nrHloc, ncHloc = size(Hloc)
    Hlocind = Vector{Vector{Tuple{Int, Int}}}()
    for r in 1:nrHloc
        temp = Vector{Tuple{Int, Int}}()
        for c in 1:ncHloc
            Hloc[r, c] != 0 && (push!(temp, (Hloc[r, c], c)))
        end
        push!(Hlocind, temp)
    end

    # make dictionary here mapping left to 1:|E|
    # this should use sizehint now so no longer slower than the standard loop
    edgemap = Dict(lv => i for (i, lv) in enumerate(left))
    currrow = 0
    H = zeros(Int, (length(right1) + length(right2)) * nrHloc, length(left))
    for rv in right1
        for r in 1:nrHloc
            @simd for c in Hlocind[r]
                H[currrow + r, edgemap[G.fadjlist[rv][c[2]]]] = c[1]
            end
        end
        currrow += nrHloc
    end

    # do second code
    Hloc = FpmattoJulia(paritycheckmatrix(C2))
    nrHloc, ncHloc = size(Hloc)
    Hlocind = Vector{Vector{Tuple{Int, Int}}}()
    for r in 1:nrHloc
        temp = Vector{Tuple{Int, Int}}()
        for c in 1:ncHloc
            Hloc[r, c] != 0 && (push!(temp, (Hloc[r, c], c)))
        end
        push!(Hlocind, temp)
    end

    for rv in right2
        for r in 1:nrHloc
            @simd for c in Hlocind[r]
                H[currrow + r, edgemap[G.fadjlist[rv][c[2]]]] = c[1]
            end
        end
        currrow += nrHloc
    end

    return H
end
