# Copyright (c) 2022, 2023 Eric Sabo, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

# TODO: make these functions accessible for LDPC code objects (in LDPC/codes.jl)
"""
    Tanner_graph_plot(H::Union{CTMatrixTypes, Matrix{Int}})

Return the Tanner graph of the matrix `H` as a Makie `Figure` object.

# Note
- Run `using Makie` to activate this extension.
"""
function Tanner_graph_plot end

"""
    Tanner_graph(H::Union{CTMatrixTypes, Matrix{Int}})

Return the `SimpleGraph` object repesenting the Tanner graph of the parity-check
matrix `H` along with the indices of the left and right vertices representing
the bits and parity checks, respectively.
"""
function Tanner_graph(H::Union{CTMatrixTypes, Matrix{Int}})
    typeof(H) <: CTMatrixTypes ? (I = _Flint_matrix_to_Julia_int_matrix(H);) : (I = H;)
    I_tr = transpose(I)
    # TODO: fix B - no zeros for this type
    B = vcat(hcat(zeros(Int, size(I_tr)), I), hcat(I_tr, zeros(Int, size(I))))
    G = SimpleGraph(B)
    nr, nc = size(H)
    # lhs - bits
    # rhs - parity checks
    return G, collect(1:nr), collect(nr + 1:nr + nc)
end

"""
    Tanner_graph(C::AbstractLinearCode)

Return the `SimpleGraph` object repesenting the Tanner graph of `C` along with
the indices of the left and right vertices representing the bits and parity checks,
respectively.
"""
Tanner_graph(C::AbstractLinearCode) = Tanner_graph(parity_check_matrix(C))

"""
    Tanner_graph(C::AbstractLDPCCode)

Return the Tanner graph of `C` as a `Figure` object.
"""
Tanner_graph(C::AbstractLDPCCode) = ismissing(C.tangr) ? (return Tanner_graph(C.H);) : (return C.tangr;)

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
# TODO: check if should not sure value in H_loc_ind and should instead just access H_loc directly
# TODO: make checks for a binary code and if so just loop and set to 1
"""
    Tanner_code(EVI::SparseMatrixCSC{Int, Int}, C::AbstractLinearCode)

Return the Tanner code obtained by applying the local code `C` to the edges of the graph with
edge-vertex incidence matrix `EVI`.
"""
function Tanner_code(EVI::SparseMatrixCSC{Int, Int}, C::AbstractLinearCode)
    num_E = EVI.m # rows
    num_V = EVI.n # columns
    num_E > num_V || throw(ArgumentError("The number of rows (edges) must be larger than the number of columns (vertices)."))
    nnz(EVI) % num_E == 0 || throw(ArgumentError("The number of vertices does not divide the number of non-zero entries, cannot be regular."))
    nnz(EVI) % (C.n - C.k) == 0 || throw(ArgumentError("The dimension of the local code does not divide the number of non-zero entries."))
    
    # pre-store all the information about H_loc
    # this could be a bit redundant if H_loc is sparse
    # H_loc = _Flint_matrix_to_Julia_int_matrix(parity_check_matrix(C))
    H_loc = parity_check_matrix(C)
    nr_H_loc, nc_H_loc = size(H_loc)
    H_loc_ind = Vector{Vector{Tuple{typeof(H_loc[1, 1]), Int}}}()
    for r in 1:nr_H_loc
        temp = Vector{Tuple{typeof(H_loc[1, 1]), Int}}()
        for c in 1:nc_H_loc
            H_loc[r, c] != 0 && (push!(temp, (H_loc[r, c], c)))
        end
        push!(H_loc_ind, temp)
    end
    H_rows_I_tr_loc = ones(Int, 1, nr_H_loc)
    H_loc_ind_lens = [length(H_loc_ind[i]) for i in 1:nr_H_loc]

    curr_row = 0
    # H = zeros(Int, num_V * nr_H_loc, num_E)
    H = zero_matrix(C.F, num_V * nr_H_loc, num_E)
    # look at every edge attached to a vertex, so check every row for a fixed column
    for c in 1:num_V
        count = 0
        # since these graphs are regular this is always the same gap and could be exploited to save a few clock cycles
        for r in EVI.colptr[c]:EVI.colptr[c + 1] - 1
            count += 1
            # this is the count-th edge, is the count-th entry of any row of H_loc
            # this loop handles all rows of H_loc in a single pass
            # doesn't actually do anything here because there's the if statement I think
            @simd for i in 1:nr_H_loc
                # instead of looping through a col of H_loc every time
                # H_loc_ind stores the next column index and since this is sorted
                # if the 2nd element of H_loc_ind at this index is count then
                # there is a 1 there, otherwise a zero since it wasn't stored
                if H_loc_ind_lens[i] >= H_rows_I_tr_loc[i] && H_loc_ind[i][H_rows_I_tr_loc[i]][2] == count
                    H[curr_row + i, EVI.rowval[r]] = H_loc_ind[i][H_rows_I_tr_loc[i]][1]
                    H_rows_I_tr_loc[i] += 1
                end
            end
        end
        curr_row += nr_H_loc
        H_rows_I_tr_loc[:] .= 1
    end
    return LinearCode(H, true)
end 

"""
    Tanner_code(G::SimpleGraph{Int}, C::AbstractLinearCode)

Return the Tanner code obtained by applying the local code `C` to the edges of `G`.
"""
function Tanner_code(G::SimpleGraph{Int}, C::AbstractLinearCode)
    isregular(G) || throw(ArgumentError("Graph must be regular."))
    length(G.fadjlist[1]) == C.n || throw(ArgumentError("The degree of the verties must be equal to the length of the local code."))
    # can use G.fadjlist directly?
    # would need to make sure we don't use the same edge twice
    return Tanner_code(sparse(transpose(incidence_matrix(G))), C)
end

"""
    Tanner_code(G::SimpleGraph{Int}, left::Vector{Int}, right::Vector{Int}, C::AbstractLinearCode)

Return the Tanner code obtained by applying the local code `C` to the vertices `right` in the
bipartition of `G` and treating the vertices of `left` as bits.
"""
function Tanner_code(G::SimpleGraph{Int}, left::Vector{Int}, right::Vector{Int}, C::AbstractLinearCode)
    # remove this to allow for overcomplete matrices like quasi-cyclic codes?
    # length(left) > length(right) || throw(ArgumentError("The size of `left` (bits) must be greater than the size of `right` (parity checks)."))
    is_valid_bipartition(G, left, right) || throw(ArgumentError("The input vectors are not a valid partition for the graph."))
    
    # pre-store all the information about H_loc
    # this could be a bit redundant if H_loc is sparse
    H_loc = parity_check_matrix(C)
    nr_H_loc, nc_H_loc = size(H_loc)
    H_loc_ind = Vector{Vector{Tuple{typeof(H_loc[1, 1]), Int}}}()
    for r in 1:nr_H_loc
        temp = Vector{Tuple{typeof(H_loc[1, 1]), Int}}()
        for c in 1:nc_H_loc
            H_loc[r, c] != 0 && (push!(temp, (H_loc[r, c], c)))
        end
        push!(H_loc_ind, temp)
    end

    # make dictionary here mapping left to 1:|E|
    # this should use sizehint now so no longer slower than the standard loop
    edge_map = Dict(lv => i for (i, lv) in enumerate(left))
    curr_row = 0
    H = zero_matrix(C.F, length(right) * nr_H_loc, length(left))
    for rv in right
        for r in 1:nr_H_loc
            @simd for c in H_loc_ind[r]
                H[curr_row + r, edge_map[G.fadjlist[rv][c[2]]]] = c[1]
            end
        end
        curr_row += nr_H_loc
    end
    return LinearCode(H, true)
end

# TODO: currently untested
"""
    Tanner_code(G::SimpleGraph{Int}, left::Vector{Int}, right1::Vector{Int}, right2::Vector{Int}, C1::AbstractLinearCode, C2::AbstractLinearCode)

Return the Tanner code obtained by applying the local code `C1` to the vertices `right1` and the local code `C2` to the vertices
`right2` in the bipartition of `G` and treating the vertices of `left` as bits.
"""
function Tanner_code(G::SimpleGraph{Int}, left::Vector{Int}, right1::Vector{Int}, right2::Vector{Int}, C1::AbstractLinearCode, C2::AbstractLinearCode)
    # remove this to allow for overcomplete matrices like quasi-cyclic codes?
    # length(left) > length(right) || throw(ArgumentError("The size of `left` (bits) must be greater than the size of `right` (parity checks)."))
    is_valid_bipartition(G, left, right1 ∪ right2) || throw(ArgumentError("The input vectors are not a valid partition for the graph."))
    
    # pre-store all the information about H_loc
    # this could be a bit redundant if H_loc is sparse
    H_loc = parity_check_matrix(C1)
    nr_H_loc, nc_H_loc = size(H_loc)
    H_loc_ind = Vector{Vector{Tuple{typeof(H_loc[1, 1]), Int}}}()
    for r in 1:nr_H_loc
        temp = Vector{Tuple{typeof(H_loc[1, 1]), Int}}()
        for c in 1:nc_H_loc
            H_loc[r, c] != 0 && (push!(temp, (H_loc[r, c], c)))
        end
        push!(H_loc_ind, temp)
    end

    # make dictionary here mapping left to 1:|E|
    # this should use sizehint now so no longer slower than the standard loop
    edge_map = Dict(lv => i for (i, lv) in enumerate(left))
    curr_row = 0
    H = zero_matrix(C.F, (length(right1) + length(right2)) * nr_H_loc, length(left))
    for rv in right1
        for r in 1:nr_H_loc
            @simd for c in H_loc_ind[r]
                H[curr_row + r, edge_map[G.fadjlist[rv][c[2]]]] = c[1]
            end
        end
        curr_row += nr_H_loc
    end

    # do second code
    H_loc = parity_check_matrix(C2)
    nr_H_loc, nc_H_loc = size(H_loc)
    H_loc_ind = Vector{Vector{Tuple{typeof(H_loc[1, 1]), Int}}}()
    for r in 1:nr_H_loc
        temp = Vector{Tuple{typeof(H_loc[1, 1]), Int}}()
        for c in 1:nc_H_loc
            H_loc[r, c] != 0 && (push!(temp, (H_loc[r, c], c)))
        end
        push!(H_loc_ind, temp)
    end

    for rv in right2
        for r in 1:nr_H_loc
            @simd for c in H_loc_ind[r]
                H[curr_row + r, edge_map[G.fadjlist[rv][c[2]]]] = c[1]
            end
        end
        curr_row += nr_H_loc
    end

    return LinearCode(H, true)
end
