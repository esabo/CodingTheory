# Copyright (c) 2024 Eric Sabo, Antonio De Martí I Olius
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

function ordered_Tanner_forest(H::T, v::T, chn::AbstractClassicalNoiseChannel, BP_alg::Symbol;
    max_iter::Int = 100, chn_inits::Union{Missing, Vector{Float64}} = missing, schedule::Symbol =
    :parallel, rand_sched::Bool = false, erasures::Vector{Int} = Int[]) where T <: CTMatrixTypes

    # TODO add MS_C
    # TODO decimation
    BP_alg ∈ (:SP, :MS) || throw(ArgumentError("`BP_alg` should be `:SP` or `:MS`"))
    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    2 <= max_iter || throw(DomainError("Number of maximum iterations must be at least two"))
    schedule ∈ (:flooding, :parallel, :serial, :layered, :semiserial) || 
        throw(ArgumentError("Unknown schedule algorithm"))
    schedule == :flooding && (schedule = :parallel;)
    schedule == :semiserial && (schedule = :layered;)

    # TODO search for wt 1 columns and add new row

    if BP_alg == :SP
        H_Int, v_Int, syndrome_based, check_adj_list, check_to_var_messages, var_to_check_messages,
            current_bits, syn = _message_passing_init_fast(H, v, chn, :BP_alg, chn_inits, :serial,
            erasures)

        if schedule == :layered
            layers = layered_schedule(H, schedule = schedule, random = rand_sched)
            flag, e, _, posteriors = _message_passing_fast_layered(H_Int, v_Int, syndrome_based,
                check_adj_list, check_to_var_messages, var_to_check_messages, current_bits, syn,
                ϕ_test, ϕ_test, max_iter, layers)
        else
            flag, e, _, posteriors = _message_passing_fast(H_Int, v_Int, syndrome_based,
                check_adj_list, check_to_var_messages, var_to_check_messages, current_bits, syn,
                ϕ_test, ϕ_test, max_iter)
        end
    elseif BP_alg == :MS
        H_Int, _, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages,
            var_to_check_messages, current_bits, totals, syn = _message_passing_init(H, v, chn,
            :MS, chn_inits, schedule, erasures)
        layers = layered_schedule(H, schedule = schedule, random = rand_sched)
        flag, e, _, posteriors = _message_passing_layered(H_Int, missing, chn_inits_2,
            _MS_check_node_message, var_adj_list, check_adj_list, max_iter, schedule,
            current_bits, totals, syn, check_to_var_messages, var_to_check_messages,
            attenuation, layers)
    end

    if !flag
        # initial BP did not converge
        ordered_indices = sortperm(posteriors, rev = true)
        erased_columns = _select_erased_columns(H, posteriors, ordered_indices)
        for i in erased_columns
            posteriors[i] = 1e-10
        end

        if BP_alg == :SP
            H_Int, v_Int, syndrome_based, check_adj_list, check_to_var_messages, var_to_check_messages,
                current_bits, syn = _message_passing_init_fast(H, v, chn, :BP_alg, chn_inits, :serial,
                erasures)
    
            if schedule == :layered
                layers = layered_schedule(H, schedule = schedule, random = rand_sched)
                flag, e, _, posteriors = _message_passing_fast_layered(H_Int, v_Int, syndrome_based,
                    check_adj_list, check_to_var_messages, var_to_check_messages, current_bits, syn,
                    ϕ_test, ϕ_test, max_iter, layers)
            else
                flag, e, _, posteriors = _message_passing_fast(H_Int, v_Int, syndrome_based,
                    check_adj_list, check_to_var_messages, var_to_check_messages, current_bits, syn,
                    ϕ_test, ϕ_test, max_iter)
            end
        elseif BP_alg == :MS
            H_Int, _, var_adj_list, check_adj_list, chn_inits_2, check_to_var_messages,
                var_to_check_messages, current_bits, totals, syn = _message_passing_init(H, v, chn,
                :MS, chn_inits, schedule, erasures)
            layers = layered_schedule(H, schedule = schedule, random = rand_sched)
            flag, e, _, posteriors = _message_passing_layered(H_Int, missing, chn_inits_2,
                _MS_check_node_message, var_adj_list, check_adj_list, max_iter, schedule,
                current_bits, totals, syn, check_to_var_messages, var_to_check_messages,
                attenuation, layers)
        end
    end

    return flag, e
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
                # println("here on $count")
                _union_by_rank!(parents, depths, seen_roots_list[col][count], seen_roots_list[col][count + 1])
                # println(parents)
            end
            count += 1
        end
        # println("parents: $parents")
        # println(depths)
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
    # println("union on $i and $j")
    # println("parents before $parents")
    if depths[i] > depths[j]
        # println("case 1")
        parents[j] = i
    elseif depths[j] > depths[i]
        parents[i] = j
        # println("case 2")
    else
        parents[j] = i
        depths[i] += 1
        # println("case 3")
    end
end

# # Example usage:
# priors = [0.1, 0.3, 0.2]                      # Priors array (length 3)
# # Create an instance of OTF_decoder with a 3x3 H matrix and a priors vector of length 3
# H = [true true false; 
#     false true true; 
#     true false true]   # H matrix (3x3)
# # Run the decode method
# decoded_result = decode(H, priors)
# # Print the decoded result
# println("Decoded Result: ", decoded_result)


# H_Int = CodingTheory._Flint_matrix_to_Julia_int_matrix(parity_check_matrix(HammingCode(2, 3)));
# ordering = collect(1:7);
# # check_adj_list = [Int[] for _ in 1:3];
# var_adj_list = [Int[] for _ in 1:7];
# for r in 1:3
#     for c in 1:7
#         if !iszero(H_Int[r, c])
#             # push!(check_adj_list[r], c)
#             push!(var_adj_list[c], r)
#         end
#     end
# end
# _select_erased_columns(H_Int, ordering, var_adj_list)
