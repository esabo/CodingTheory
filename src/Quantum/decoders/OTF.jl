# Copyright (c) 2024 Eric Sabo, Antonio De Martí I Olius
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

function ordered_Tanner_forest(
    H::T,
    v::T,
    chn::AbstractClassicalNoiseChannel;
    BP_alg::Symbol = :SP,
    max_iter::Int = 100,
    chn_inits::Union{Missing,Vector{Float64}} = missing,
    schedule::Symbol = :parallel,
    rand_sched::Bool = false,
    erasures::Vector{Int} = Int[],
) where {T<:CTMatrixTypes}

    # TODO add MS_C
    # TODO decimation
    BP_alg ∈ (:SP, :MS) || throw(ArgumentError("`BP_alg` should be `:SP` or `:MS`"))
    Int(order(base_ring(H))) == 2 ||
        throw(ArgumentError("Currently only implemented for binary codes"))
    nr, nc = size(H)
    (nr ≥ 0 && nc ≥ 0) || throw(ArgumentError("H cannot have a zero dimension"))
    2 <= max_iter || throw(DomainError("Number of maximum iterations must be at least two"))
    schedule ∈ (:flooding, :parallel, :serial, :layered, :semiserial) ||
        throw(ArgumentError("Unknown schedule algorithm"))
    schedule == :flooding && (schedule = :parallel;)
    schedule == :semiserial && (schedule = :layered;)

    # TODO search for wt 1 columns and add new row
    H_Int,
    v_Int,
    syndrome_based,
    check_adj_list,
    check_to_var_messages,
    var_to_check_messages,
    current_bits,
    syn = _message_passing_init_fast(H, v, chn, BP_alg, chn_inits, :serial, erasures)

    if BP_alg == :SP
        # H_Int, v_Int, syndrome_based, check_adj_list, check_to_var_messages, var_to_check_messages,
        #     current_bits, syn = _message_passing_init_fast(H, v, chn, BP_alg, chn_inits, :serial,
        #     erasures)

        if schedule == :layered
            layers = layered_schedule(H, schedule = schedule, random = rand_sched)
            flag, e, _, posteriors = _message_passing_fast_layered(
                H_Int,
                v_Int,
                syndrome_based,
                check_adj_list,
                check_to_var_messages,
                var_to_check_messages,
                current_bits,
                syn,
                ϕ_test,
                ϕ_test,
                max_iter,
                layers,
            )
        else
            flag, e, _, posteriors = _message_passing_fast(
                H_Int,
                v_Int,
                syndrome_based,
                check_adj_list,
                check_to_var_messages,
                var_to_check_messages,
                current_bits,
                syn,
                ϕ_test,
                ϕ_test,
                max_iter,
            )
        end
    elseif BP_alg == :MS
        H_Int,
        _,
        var_adj_list,
        check_adj_list,
        chn_inits_2,
        check_to_var_messages,
        var_to_check_messages,
        current_bits,
        totals,
        syn = _message_passing_init(H, v, chn, :MS, chn_inits, schedule, erasures)
        layers = layered_schedule(H, schedule = schedule, random = rand_sched)
        flag, e, _, posteriors = _message_passing_layered(
            H_Int,
            missing,
            chn_inits_2,
            _MS_check_node_message,
            var_adj_list,
            check_adj_list,
            max_iter,
            schedule,
            current_bits,
            totals,
            syn,
            check_to_var_messages,
            var_to_check_messages,
            attenuation,
            layers,
        )
    end

    if !flag
        # initial BP did not converge
        var_adj_list = [Int[] for _ = 1:size(H_Int, 2)];
        for r = 1:size(H_Int, 1)
            for c = 1:size(H_Int, 2)
                if !iszero(H_Int[r, c])
                    push!(var_adj_list[c], r)
                end
            end
        end

        OTF_type = :OSD
        if syndrome_based
            if OTF_type == :OSD
                # sort LLRs from greatest to least (most positive to most negative)
                # since negative implies an error is more likely here, the selection process will allow BP to run on columns which are still positive while fixing columns which are more likely to have an error to have an error
                # this is similar to selecting a test pattern in OSD
                ordered_indices = sortperm(posteriors, rev = true)
                erased_columns =
                    _select_erased_columns(H_Int, ordered_indices, var_adj_list)
                for i in erased_columns
                    posteriors[i] = -20.0
                end
            else
                ordered_indices = sortperm(posteriors, by = abs)
                erased_columns =
                    _select_erased_columns(H_Int, ordered_indices, var_adj_list)
                for i in erased_columns
                    if posteriors[i] < 0
                        posteriors[i] = -20.0
                    else
                        posteriors[i] = 20.0
                    end
                end
            end
        else
            ordered_indices = sortperm(posteriors, by = abs)
            erased_columns = _select_erased_columns(H_Int, ordered_indices, var_adj_list)
            for i in erased_columns
                if posteriors[i] < 0
                    posteriors[i] = -20.0
                else
                    posteriors[i] = 20.0
                end
            end
        end
        println(erased_columns)
        println(" ")
        println(posteriors)

        if BP_alg == :SP
            H_Int,
            v_Int,
            syndrome_based,
            check_adj_list,
            check_to_var_messages,
            var_to_check_messages,
            current_bits,
            syn =
                _message_passing_init_fast(H, v, chn, BP_alg, chn_inits, :serial, erasures)

            if schedule == :layered
                layers = layered_schedule(H, schedule = schedule, random = rand_sched)
                flag, e, _, posteriors = _message_passing_fast_layered(
                    H_Int,
                    v_Int,
                    syndrome_based,
                    check_adj_list,
                    check_to_var_messages,
                    var_to_check_messages,
                    current_bits,
                    syn,
                    ϕ_test,
                    ϕ_test,
                    max_iter,
                    layers,
                )
            else
                flag, e, _, posteriors = _message_passing_fast(
                    H_Int,
                    v_Int,
                    syndrome_based,
                    check_adj_list,
                    check_to_var_messages,
                    var_to_check_messages,
                    current_bits,
                    syn,
                    ϕ_test,
                    ϕ_test,
                    max_iter,
                )
            end
        elseif BP_alg == :MS
            H_Int,
            _,
            var_adj_list,
            check_adj_list,
            chn_inits_2,
            check_to_var_messages,
            var_to_check_messages,
            current_bits,
            totals,
            syn = _message_passing_init(H, v, chn, :MS, chn_inits, schedule, erasures)
            layers = layered_schedule(H, schedule = schedule, random = rand_sched)
            flag, e, _, posteriors = _message_passing_layered(
                H_Int,
                missing,
                chn_inits_2,
                _MS_check_node_message,
                var_adj_list,
                check_adj_list,
                max_iter,
                schedule,
                current_bits,
                totals,
                syn,
                check_to_var_messages,
                var_to_check_messages,
                attenuation,
                layers,
            )
        end
    end

    return flag, e
end

function _select_erased_columns(
    H::Matrix{UInt8},
    ordered_indices::Vector{Int},
    var_adj_list::Vector{Vector{Int}},
)

    # this is using the disjoint-set data structure/union find algorithm for merging them
    nr = size(H, 1)
    parents = collect(1:nr)
    depths = ones(Int, nr)
    output_indices = Vector{Int}()
    seen_roots_list =
        [[-1 for _ = 1:length(var_adj_list[c])] for c = 1:length(var_adj_list)]
    flag = false
    for col in ordered_indices
        # count = 0
        for row in var_adj_list[col]
            row_root = _find_root(parents, row)
            flag = _check_roots_list!(seen_roots_list, col, row_root)
            if flag
                # cycle
                push!(output_indices, col)
                break
            end
        end

        if !flag
            # Find maximum of depths[see_roots_list[col][count]], also, see if there are two or more roots of maximum depth. In which case we should grow.
            merged_root = -1
            merged_depth = 0
            growth = false
            for row_root in seen_roots_list[col]
                if depths[row_root] > merged_depth
                    merged_root = row_root
                    merged_depth = depths[row_root]
                    growth = false
                elseif depths[row_root] == merged_depth
                    growth = true
                end
            end

            for row_root in seen_roots_list[col]
                parents[row_root] = merged_root
            end

            if growth
                depths[merged_root] += 1
            end

        end
        # println(parents)
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
    if depths[i] > depths[j]
        parents[j] = i
    elseif depths[j] > depths[i]
        parents[i] = j
    else
        parents[j] = i
        depths[i] += 1
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


# if syndrome based
#     OSD picture:
#         sort LLRs from greatest to least (most positive to most negative)
#         since negative implies an error is more likely here, the selection process will allow BP to run on columns which are still positive while fixing columns which are more likely to have an error to have an error
#         this is similar to selecting a test pattern in OSD
#         ideally we can fix these, determine a syndrome for this, then run BP with this new syndrome
#         the final result comes from combing the correction and fixed error
#         I think it may suffice to skip this step and simply continue BP with the new fixed LLRs
#         it is possible a column all the way on the left of the sort (with clearly no error here) can produce a loop and would therefore be fixed to an error
#         this may be okay in the quantum picture due to degeneracy
#         this would dramatically kick BP out of its local minimum
#         we hope this doesn't kick us out so far to land in another logical coset (can we control this?)
#         codes with a ton of loops (such as QCCs) may fix wayy too many columns to be errors
#         as long as this doesn't cause BP to fail, this may be find for quantum codes due to degeneracy, although we increase the risk we jumped into another logical coset

#     BP picture:
#         the bits which are most reliably correct are those with high LLR values
#         so sort by absolute value
#         keep solving on the least reliable nodes while fixing bits we are more sure of
#         this can have the same issues as above:
#             assigning a definite value to an index which BP isn't yet sure about
#             fixing too many bits
#         as above, this relies on degeneracy
#         to mimic the above approach, we still set all loopy columns to 1, even if it is strongly known to be not in error
#         however, it may be advantageous to "listen" to BP better by fixing indices based on their sign, keeping nonerrors BP is sure about
#         this can help by reducing the fixed error wt, helping keep us in the logical coset
# else
#     we are unable to assign an interpretation to 0's and 1's in the received vector
#     we could do standard reliability decoding where we take a parameter and look for test patterns inside the k least-reliable bits
#     or we could proceed as above with the BP picture
