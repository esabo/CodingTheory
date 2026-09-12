# Copyright (c) 2022 - 2026 Eric Sabo, Michael Vasmer
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
$(TYPEDSIGNATURES)

Return the `SimpleGraph` object repesenting the Tanner graph of the parity-check
matrix `H` along with the indices of the left and right vertices representing
the bits and parity checks, respectively.
"""
function Tanner_graph(H::Union{CTMatrixTypes, Matrix{Int}})
    typeof(H) <: CTMatrixTypes ? (I = _Flint_matrix_to_Julia_int_matrix(H);) : (I = H;)
    nr, nc = size(I)
    B = vcat(hcat(zeros(Int, nc, nc), transpose(I)), hcat(I, zeros(Int, nr, nr)))
    G = SimpleGraph(B)
    # lhs - bits
    # rhs - parity checks
    return G, collect(1:nr), collect(nr + 1:nr + nc)
end

"""
$(TYPEDSIGNATURES)

Return the `SimpleGraph` object repesenting the Tanner graph of `C` along with
the indices of the left and right vertices representing the bits and parity checks,
respectively.
"""
Tanner_graph(C::AbstractLinearCode) = Tanner_graph(parity_check_matrix(C))

# """
#     Tanner_graph(C::AbstractLDPCCode)

# Return the Tanner graph of `C` as a `Figure` object.
# """
# Tanner_graph(C::AbstractLDPCCode) = Tanner_graph(C.H)

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
function parity_check_matrix(C::TannerCode, stand_form::Bool = false)
    cache = getfield(C, :cache)
    if !haskey(cache, :H)
        build_method = cache[:build_method]
        
        # ---------------------------------------------------------
        # RECIPE 1: Edge-Vertex Incidence Matrix
        # ---------------------------------------------------------
        if build_method == :EVI
            EVI = cache[:EVI]
            C_loc = cache[:C_local]
            H_loc = parity_check_matrix(C_loc)
            nr_H_loc, nc_H_loc = size(H_loc)
            
            num_E = EVI.m
            num_V = EVI.n
            
            H_loc_ind = [[(H_loc[r, c], c) for c in 1:nc_H_loc if H_loc[r, c] != 0] for r in 1:nr_H_loc]
            H_rows_I_tr_loc = ones(Int, 1, nr_H_loc)
            H_loc_ind_lens = length.(H_loc_ind)

            curr_row = 0
            H = zero_matrix(C.F, num_V * nr_H_loc, num_E)
            
            for c in 1:num_V
                count = 0
                for r in EVI.colptr[c]:(EVI.colptr[c + 1] - 1)
                    count += 1
                    @simd for i in 1:nr_H_loc
                        if H_loc_ind_lens[i] >= H_rows_I_tr_loc[i] && H_loc_ind[i][H_rows_I_tr_loc[i]][2] == count
                            H[curr_row + i, EVI.rowval[r]] = H_loc_ind[i][H_rows_I_tr_loc[i]][1]
                            H_rows_I_tr_loc[i] += 1
                        end
                    end
                end
                curr_row += nr_H_loc
                H_rows_I_tr_loc[:] .= 1
            end
            cache[:H] = H
            
        # ---------------------------------------------------------
        # RECIPE 2: Single Bipartition
        # ---------------------------------------------------------
        elseif build_method == :bipartition
            G = cache[:G]
            left = cache[:left]
            right = cache[:right]
            C_loc = cache[:C_local]
            
            H_loc = parity_check_matrix(C_loc)
            nr_H_loc, nc_H_loc = size(H_loc)
            H_loc_ind = [[(H_loc[r, c], c) for c in 1:nc_H_loc if H_loc[r, c] != 0] for r in 1:nr_H_loc]

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
            cache[:H] = H
            
        # ---------------------------------------------------------
        # RECIPE 3: Dual Bipartition (Two Codes)
        # ---------------------------------------------------------
        elseif build_method == :bipartition_dual
            G = cache[:G]
            left = cache[:left]
            right1 = cache[:right1]
            right2 = cache[:right2]
            C1 = cache[:C1]
            C2 = cache[:C2]
            
            edge_map = Dict(lv => i for (i, lv) in enumerate(left))
            curr_row = 0
            H = zero_matrix(C.F, (length(right1) + length(right2)) * nrows(parity_check_matrix(C1)), length(left))
            
            # Sub-routine to apply a local code to a specific right partition
            function _apply_local_code!(H_out, C_local, right_part, row_offset)
                H_l = parity_check_matrix(C_local)
                nr_H, nc_H = size(H_l)
                H_ind = [[(H_l[r, c], c) for c in 1:nc_H if H_l[r, c] != 0] for r in 1:nr_H]
                
                for rv in right_part
                    for r in 1:nr_H
                        @simd for c in H_ind[r]
                            H_out[row_offset + r, edge_map[G.fadjlist[rv][c[2]]]] = c[1]
                        end
                    end
                    row_offset += nr_H
                end
                return row_offset
            end
            
            curr_row = _apply_local_code!(H, C1, right1, curr_row)
            _apply_local_code!(H, C2, right2, curr_row)
            
            cache[:H] = H
        end
    end
    
    if stand_form
        generator_matrix(C, true)
        return cache[:H_stand]
    end
    return cache[:H]
end

function generator_matrix(C::TannerCode, stand_form::Bool = false)
    cache = getfield(C, :cache)
    if !haskey(cache, :G)
        H_lift = parity_check_matrix(C)
        _, G = right_kernel(H_lift)
        cache[:G] = transpose(G)
    end
    
    if stand_form
        generator_matrix(C, true)
        return cache[:G_stand]
    end
    return cache[:G]
end

#############################
        # constructors
#############################

"""
$(TYPEDSIGNATURES)

Return the Tanner code obtained by applying the local code `C` to the edges of the graph with
edge-vertex incidence matrix `EVI`.
"""
function TannerCode(EVI::SparseMatrixCSC{Int, Int}, C::AbstractLinearCode)
    num_E = EVI.m 
    num_V = EVI.n 
    num_E > num_V || throw(ArgumentError("The number of edges must be larger than the number of vertices."))
    nnz(EVI) % num_E == 0 || throw(ArgumentError("Graph cannot be regular."))
    nnz(EVI) % (C.n - C.k) == 0 || throw(ArgumentError("Dimension of local code does not divide non-zero entries."))
    
    n_new = num_E
    k_bound = max(1, num_E - (num_V * (C.n - C.k)))
    
    cache = Dict{Symbol, Any}(:build_method => :EVI, :EVI => EVI, :C_local => C)
    return TannerCode(C.F, n_new, k_bound, missing, 1, n_new, cache)
end

"""
$(TYPEDSIGNATURES)

Return the Tanner code obtained by applying the local code `C` to the edges of `G`.
"""
function TannerCode(G::SimpleGraph{Int}, C::AbstractLinearCode)
    isregular(G) || throw(ArgumentError("Graph must be regular."))
    length(G.fadjlist[1]) == C.n || throw(ArgumentError("The degree of the vertices must be equal to the length of the local code."))
    
    return TannerCode(sparse(transpose(incidence_matrix(G))), C)
end

"""
$(TYPEDSIGNATURES)

Return the Tanner code obtained by applying the local code `C` to the vertices `right` in the
bipartition of `G` and treating the vertices of `left` as bits.
"""
function TannerCode(G::SimpleGraph{Int}, left::Vector{Int}, right::Vector{Int}, C::AbstractLinearCode)
    # is_valid_bipartition(G, left, right) || throw(ArgumentError("The input vectors are not a valid partition for the graph."))
    
    n_new = length(left)
    k_bound = max(1, length(left) - (length(right) * (C.n - C.k)))
    
    cache = Dict{Symbol, Any}(
        :build_method => :bipartition, 
        :G => G, :left => left, :right => right, :C_local => C
    )
    return TannerCode(C.F, n_new, k_bound, missing, 1, n_new, cache)
end

"""
$(TYPEDSIGNATURES)

Return the Tanner code obtained by applying the local codes `C1` and `C2` to the vertices `right1` and
`right2` in the bipartition of `G` and treating the vertices of `left` as bits.
"""
function TannerCode(G::SimpleGraph{Int}, left::Vector{Int}, right1::Vector{Int}, right2::Vector{Int}, C1::AbstractLinearCode, C2::AbstractLinearCode)
    # is_valid_bipartition(G, left, right1 ∪ right2) || throw(ArgumentError("The input vectors are not a valid partition for the graph."))
    
    n_new = length(left)
    k_bound = max(1, length(left) - (length(right1) * (C1.n - C1.k)) - (length(right2) * (C2.n - C2.k)))
    
    cache = Dict{Symbol, Any}(
        :build_method => :bipartition_dual, 
        :G => G, :left => left, :right1 => right1, :right2 => right2, :C1 => C1, :C2 => C2
    )
    return TannerCode(C1.F, n_new, k_bound, missing, 1, n_new, cache)
end

"""
$(TYPEDSIGNATURES)

Return the eigenvalues of the base graph used to construct the Tanner code.
Computes lazily and caches the result.
"""
function graph_eigenvalues(C::TannerCode)
    cache = getfield(C, :cache)
    if !haskey(cache, :eigenvalues)
        build_method = cache[:build_method]
        
        if build_method == :EVI
            # For EVI, the underlying graph's adjacency matrix can be extracted from the incidence matrix
            EVI = cache[:EVI]
            A = EVI' * EVI
            # Subtract diagonal (degrees) to get the true adjacency matrix
            A[diagind(A)] .= 0 
            cache[:eigenvalues] = eigvals(Symmetric(Matrix(A)))
        else
            # For bipartition methods, we extract the adjacency matrix directly from the SimpleGraph
            G = cache[:G]
            A = adjacency_matrix(G)
            cache[:eigenvalues] = eigvals(Symmetric(Matrix(A)))
        end
    end
    return cache[:eigenvalues]
end

"""
$(TYPEDSIGNATURES)

Return the second-largest eigenvalue (in absolute value) of the base graph, often denoted `λ`.
"""
function spectral_gap(C::TannerCode)
    evals = graph_eigenvalues(C)
    # The eigenvalues are sorted in ascending order. 
    # For bipartite/regular graphs, the largest in absolute value are at the extremes.
    # We strip the largest eigenvalue (degree d) and find the maximum of the absolute values of the rest.
    max_eval = maximum(abs.(evals))
    
    # Filter out the trivial largest eigenvalues (which equal d or -d for bipartite graphs)
    non_trivial_evals = filter(x -> abs(x) < max_eval - 1e-7, evals)
    
    isempty(non_trivial_evals) && return 0.0 # Fully disconnected or trivial graph
    return maximum(abs.(non_trivial_evals))
end

"""
$(TYPEDSIGNATURES)

Return the Sipser-Spielman spectral lower bound on the minimum distance of the Tanner code.
"""
function Sipser_Spielman_bound(C::TannerCode)
    λ = spectral_gap(C)
    
    cache = getfield(C, :cache)
    build_method = cache[:build_method]
    
    # We need the local code's minimum distance and the graph's regularity degree
    if build_method == :EVI || build_method == :bipartition
        C_loc = cache[:C_local]
        d_0 = minimum_distance(C_loc)
        
        if build_method == :EVI
            # For EVI, left degree is the column weight
            EVI = cache[:EVI]
            d = length(EVI.colptr[1]:EVI.colptr[2]-1) 
        else
            G = cache[:G]
            left = cache[:left]
            d = length(G.fadjlist[left[1]])
        end
        
        # Sipser-Spielman Expansion Bound
        if λ < d_0
            bound = C.n * ((d_0^2 - λ * d_0) / (d^2 - λ * d_0))
            return max(1, floor(Int, bound))
        end
    end
    
    return 1 # Fallback trivial bound if spectral conditions aren't met or dual-bipartition is used
end
