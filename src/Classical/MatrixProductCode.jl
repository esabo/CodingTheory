# Copyright (c) 2022, 2023 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
# constructors
#############################

# TODO: need to make sure codes are not overcomplete
"""
    MatrixProductCode(C::Vector{AbstractLinearCode}, A::CTMatrixTypes)

Return the matrix product code defined by the vector of linear codes `C` and matrix `A`.
"""
function MatrixProductCode(C::Vector{AbstractLinearCode}, A::CTMatrixTypes)
    isempty(C) && throw(ArgumentError("Vector of linear codes cannot be empty."))
    iszero(A) && throw(ArgumentError("Matrix A cannot be zero."))
    s, l = size(A)
    s == length(C) ||
        throw(ArgumentError("Number of rows of A must be equal to the number of codes."))
    F = C[1].F
    n = C[1].n
    for i = 2:s
        F == C[i].F || throw(ArgumentError("All codes must have the same base ring."))
        n == C[i].n || throw(ArgumentError("All codes must have the same length."))
    end
    F == base_ring(A) ||
        throw(ArgumentError("Codes and matrix must have the same base ring."))

    # H formula holds in special case only
    G = zero_matrix(F, sum([C[i].k for i = 1:s], l * n))
    # H = zero_matrix(F, sum([nrows(C[i].H) for i in 1:s], l * n))
    curr = 1
    # currH = 1
    # need to do in this row/column order
    for r = 1:s
        for c = 1:l
            G[curr:(curr+C[r].k), (1+(r-1)*n):(r*n)] =
                A[r, c] * generator_matrix(C[r], true)
            # H[currH:currH + nrows(C[r].H), 1 + (r - 1) * n:r * n] = A[r, c] * paritycheckmatrix(C[r])
        end
        curr += C[r].k
        # currH += nrows(C[r].H)
    end

    G_stand, H_stand, P, k = _standard_form(G)
    if ismissing(P)
        _, H = right_kernel(G)
        # note the H here is transpose of the standard definition
        H = _remove_empty(transpose(H), :rows)
    else
        H = H_stand * P
    end

    ub1, _ = _min_wt_row(G)
    ub2, _ = _min_wt_row(G_stand)
    ub = minimum([ub1, ub2])
    return MatrixProductCode(
        F,
        n,
        k,
        1,
        ub,
        missing,
        G,
        H,
        G_stand,
        H_stand,
        P,
        missing,
        C,
        A,
    )
end

#############################
# getter functions
#############################

# TODO: why are there no getter functions for this? C, A inputs?

#############################
# setter functions
#############################

#############################
# general functions
#############################
