# Copyright (c) 2023, 2024 Eric Sabo, Benjamin Ide, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

function _reduce_rows(H::Union{CTMatrixTypes, MatElem{<: ResElem}}, rows::AbstractVector{Int},
    compressed::Bool = false, permute::Bool = false, target::Int = 3)

    # TODO: for compressed, should use the target weight as well so it can be even more compressed if possible.
    compressed && target > 3 && @warn("Compressed weight-reduction for weights greater than 3 not yet implemented.")

    isempty(rows) && return H
    nr, nc = size(H)
    if compressed
        temp = sum(count(!iszero, H[r, :]) - 3 for r in rows)
        H_new = zero_matrix(base_ring(H), nr + temp, nc + temp)
        t_C(n) = matrix(base_ring(H), diagm(n, n - 1, 0 => ones(Int, n - 1), -1 => ones(Int, n - 1)))
        i = 1
        j = nc + 1
        for r in 1:nr
            if r in rows
                wt = count(!iszero, H[r, :])
                for (k, p) in enumerate(permute ? shuffle([[i, i + wt - 3]; i:i + wt - 3]) : [[i, i + wt - 3]; i:i + wt - 3])
                    columns = getindex.(findall(!iszero, H[r, :]), 2)
                    H_new[p, columns[k]] = H[r, columns[k]]
                end
                H_new[i:i + wt - 3, j:j + wt - 4] = t_C(wt - 2)
                i += wt - 2
                j += wt - 3
            else
                H_new[i, 1:nc] = H[r, 1:nc]
                i += 1
            end
        end
    else
        temp = sum(count(!iszero, H[r, :]) - 1 for r in rows)
        H_new = zero_matrix(base_ring(H), nr + temp, nc + temp)
        t(n) = matrix(base_ring(H), diagm(n, n - 1, 0 => ones(Int, n - 1), -1 => ones(Int, n - 1)))
        i = 1
        j = nc + 1
        for r in 1:nr
            if r in rows
                wt = count(!iszero, H[r, :])
                for (k, p) in enumerate(permute ? shuffle(i:i + wt - 1) : i:i + wt - 1)
                    columns = getindex.(findall(!iszero, H[r, :]), 2)
                    H_new[p, columns[k]] = H[r, columns[k]]
                end
                H_new[i:i + wt - 1, j:j + wt - 2] = t(wt)
                i += wt
                j += wt - 1
            else
                H_new[i, 1:nc] = H[r, 1:nc]
                i += 1
            end
        end
    end
    return H_new
end

"""
weight_reduction(H::Union{CTMatrixTypes, MatElem{<: ResElem}}; rows::Bool = true,
    row_indices::AbstractVector{Int} = Int[], permute_rows::Bool = true, row_target::Int = 3,
    columns::Bool = true, column_indices::AbstractVector{Int} = Int[], permute_columns::Bool = true,
    column_target::Int = 3, compressed::Bool = false, seed::Union{Nothing, Int} = nothing)

Return the weight-reduced parity-check matrix of `H` with the given arguments.
"""
function weight_reduction(H::Union{CTMatrixTypes, MatElem{<: ResElem}}; rows::Bool = true,
    row_indices::AbstractVector{Int} = Int[], permute_rows::Bool = true, row_target::Int = 3,
    columns::Bool = true, column_indices::AbstractVector{Int} = Int[], permute_columns::Bool = true,
    column_target::Int = 3, compressed::Bool = false, seed::Union{Nothing, Int} = nothing)

    # method ∈ (:compressed, :identity) || throw(ArgumentError("Unknown reduction method"))
    all(1 <= r <= nrows(H) for r in row_indices) || throw(ArgumentError("Invalid row indices"))
    all(1 <= r <= ncols(H) for r in column_indices) || throw(ArgumentError("Invalid col indices"))
    Random.seed!(seed)

    if rows
        if isempty(row_indices)
            _, row_wts = _degree_distribution(H)
            row_indices = findall(x -> x > row_target, row_wts)
        end
        H_reduced = _reduce_rows(H, row_indices, compressed, permute_rows, row_target)
    else
        H_reduced = H
    end

    if columns
        if isempty(column_indices)
            col_wts, _ = _degree_distribution(H_reduced)
            column_indices = findall(x -> x > column_target, col_wts)
        end
        H_reduced = transpose(_reduce_rows(transpose(H_reduced), column_indices, compressed, permute_columns, column_target))
    end
    return H_reduced
end

"""
weight_reduction(C::Union{AbstractLinearCode, AbstractLDPCCode}; rows::Bool = true,
    row_indices::AbstractVector{Int} = Int[], permute_rows::Bool = true, row_target::Int = 3,
    columns::Bool = true, column_indices::AbstractVector{Int} = Int[], permute_columns::Bool = true,
    column_target::Int = 3, compressed::Bool = false, seed::Union{Nothing, Int} = nothing)

Return the code with the weight-reduced parity-check matrix of `C` with the given arguments.
"""
function weight_reduction(C::Union{AbstractLinearCode, AbstractLDPCCode}; rows::Bool = true,
    row_indices::AbstractVector{Int} = Int[], permute_rows::Bool = true, row_target::Int = 3,
    columns::Bool = true, column_indices::AbstractVector{Int} = Int[], permute_columns::Bool = true,
    column_target::Int = 3, compressed::Bool = false, seed::Union{Nothing, Int} = nothing)

    H_reduced = weight_reduction(parity_check_matrix(C), rows = rows, row_indices = row_indices,
        permute_rows = permute_rows, row_target = row_target, columns = columns,
        column_indices = column_indices, permute_columns = permute_columns,
        column_target = column_target, compressed = compressed, seed = seed)
    
    if C isa QuasiCyclicCode
        return QuasiCyclicCode(H_reduced, true)
    elseif C isa AbstractLinearCode
        return LinearCode(H_reduced, true)
    else
        return LDPCCode(H_reduced)
    end
end
