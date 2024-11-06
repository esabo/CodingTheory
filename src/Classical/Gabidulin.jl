# Copyright (c) 2024 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

# TODO currently untested and have no unit tests

"""
    GeneralizedGabidulinCode(F::CTFieldTypes, eval_pts::Vector{CTFieldElem}, k::Int, s::Int; parity::Bool = false)

Return the vector representation of the dimension `k` generalized Gabidulin code given the
evaluation points `eval_pts` with respect to the subfield `F` of the base ring of `eval_pts`.

# Notes
- Distances are reported with respect to the Hamming metric.
"""
function GeneralizedGabidulinCode(F::CTFieldTypes, eval_pts::Vector{CTFieldElem}, k::Int, s::Int; parity::Bool = false)

    is_empty(eval_pts) && throw(ArgumentError("The input vector `eval_pts` cannot be empty."))
    is_positive(s) || throw(DomainError(s, "The parameter `s` must be positive."))
    E = parent(eval_pts[1])
    all(parent(pt) == E for pt in eval_pts) || throw(ArgumentError("All evaluation points must be over the same base ring."))
    flag, m = is_extension_field(E, F)
    flag || throw(ArgumentError("The input field is not a subfield of the base ring of the evaluation points."))
    n = length(eval_pts)
    n ≤ m || throw(ArgumentError("The number of evaluation points must be less than or equal to the degree of the field extension."))
    0 < k ≤ n || throw(DomainError(k, "The code dimenion must be `0 < k ≤ n`."))
    q = Int(order(F))

    B = zero_matrix(E, n, n)
    for r in 1:n
        for c in 1:n
            B[r, c] = eval_pts[c]^(q^(n - 1))
        end
    end
    iszero(det(B)) || throw(ArgumentError("The evaluation points must be linearly independent."))

    G = zero_matrix(E, k, n)
    for r in 1:k
        for c in 1:n
            G[r, c] = eval_pts[c]^(q^(s * (r - 1)))
        end
    end

    return LinearCode(G, parity)
    # reaches Singleton bound for the rank metric: d_R = n - k + 1
end

"""
    GabidulinCode(F::CTFieldTypes, eval_pts::Vector{CTFieldElem}, k::Int; parity::Bool = false)

Return the vector representation of the dimension `k` Gabidulin code given the evaluation points
`eval_pts` with respect to the subfield `F` of the base ring of `eval_pts`.

# Notes
- Distances are reported with respect to the Hamming metric.
"""
GabidulinCode(F::CTFieldTypes, eval_pts::Vector{CTFieldElem}, k::Int; parity::Bool = false) =
    GeneralizedGabidulinCode(F, eval_pts, k, 1, parity = parity)

# TODO the dual is also a Gabidulin code, need to figure out those eval points based on these
# but that would require making a type for this named code
# TODO are these MDS or just in rank metric?
