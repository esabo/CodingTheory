# Copyright (c) 2024 Eric Sabo, Antonio De Martí I Olius
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

function ordered_Tanner_forest(C::AbstractStabilizerCode, BP_alg::Symbol, max_iters::Int, )

    BP_alg ∈ (:SP, :MS) || throw(ArgumentError("`BP_alg` should be `:SP` or `:MS`"))

    # 1. Run BP
    # 2. Sort soft info from output
    # 3. Make hypergraph (implicit step)
    # 4. start with least reliable node and record edges to neighbors
    # 5. Move to next least reliable node, store same making sure no loop is added
    # 6. Repeat until done
    # 7. For hyperedges that are never added, set vertex to erased in original Tanner graph
    # 8. Run BP

    if BP_alg == :SP
        flag, e, _, posteriors = sum_product(...)
    elseif BP_alg == :MS
        flag, e, _, posteriors = min_sum(...)
    end

    if !flag
        # initial BP did not converge
        ordered_indices = sortperm(posteriors, rev = true)
        erased_columns = _select_erased_columns(H, posteriors, ordered_indices)
        for i in erased_columns
            posteriors[i] = 0.0
        end

        if BP_alg == :SP
            flag, e, _, posteriors = sum_product(...)
        elseif BP_alg == :MS
            flag, e, _, posteriors = min_sum(...)
        end
    else
        # initial BP converged

    end

    return
end

# Define the decode method
function decode(H::Matrix{Bool}, priors::Vector{Float64})



    
    # Here we should add the decoding process.
    
    otf_probabilities = copy(priors)
    # This should be a very low value.
    otf_probabilities[indices .== true] .= 1e-10
    # This is where you'd return the actual decoded data
    return otf_probabilities
end

function _select_erased_columns(H::Matrix{Int}, ordered_indices::Vector{Int}, var_adj_list::Vector{Vector{Int}})

    # this is using the disjoint-set data structure/union find algorithm for merging them
    nr = size(H, 1)
    parents = collect(1:nr)
    depths = ones(Int, nr)
    output_indices = falses(size(H, 2))
    seen_roots_list = [[-1 for _ in 1:length(var_adj_list[c])] for c in 1:length(var_adj_list)]
    for col in ordered_indices
        # println("col $col")
        count = 0
        for row in var_adj_list[col]
            row_root = _find_root(parents, row)
            flag = _check_roots_list!(seen_roots_list, col, row_root)
            # println(seen_roots_list[col])
            if flag
                # cycle
                # println("loop at row $row with root $row_root")
                output_indices[col] = true
                break
            elseif count ≥ 1
                _union_by_rank!(parents, depths, seen_roots_list[col][count], seen_roots_list[col][count + 1])
                count += 1
                # println(parents)
            else
                count += 1
            end
        end
    end
    return output_indices
end

# find w/ path compression
function _find_root(parents::Vector{Int}, i::Int)
    if parents[i] ≠ i
        parents[i] = _find_root(parents, parents[i])
    end
    return parents[i]
end

# union
function _check_roots_list!(seen_roots_list::Vector{Vector{Int}}, c::Int, new_root::Int)
    i = 1
    len = length(seen_roots_list[c])
    while i ≤ len && seen_roots_list[c][i] ≠ -1
        seen_roots_list[c][i] == new_root && return true
        i += 1
    end
    # hasn't found it and not yet at the end of the list
    i ≤ len && (seen_roots_list[c][i] = new_root; return false)
    return false
end

# union by rank
function _union_by_rank!(parents::Vector{Int}, depths::Vector{Int}, i::Int, j::Int)
    if depths[i] > depths[j]
        parents[j] == i
    elseif depths[j] > depths[i]
        parents[i] = j
    else
        parents[j] = i
        depths[i] += 1
    end
end

# Example usage:
priors = [0.1, 0.3, 0.2]                      # Priors array (length 3)
# Create an instance of OTF_decoder with a 3x3 H matrix and a priors vector of length 3
H = [true true false; 
    false true true; 
    true false true]   # H matrix (3x3)
# Run the decode method
decoded_result = decode(H, priors)
# Print the decoded result
println("Decoded Result: ", decoded_result)


H_Int = CodingTheory._Flint_matrix_to_Julia_int_matrix(parity_check_matrix(HammingCode(2, 3)));
ordering = collect(1:7);
# check_adj_list = [Int[] for _ in 1:3];
var_adj_list = [Int[] for _ in 1:7];
for r in 1:3
    for c in 1:7
        if !iszero(H_Int[r, c])
            # push!(check_adj_list[r], c)
            push!(var_adj_list[c], r)
        end
    end
end
_select_erased_columns(H_Int, ordering, var_adj_list)
