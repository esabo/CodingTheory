# Copyright (c) 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    u_u_plus_v(C1::AbstractLinearCode, C2::AbstractLinearCode)
    Plotkin_construction(C1::AbstractLinearCode, C2::AbstractLinearCode)

Return the Plotkin (u | u + v)-construction with `u ∈ C1` and `v ∈ C2`.
"""
function u_u_plus_v(C1::AbstractLinearCode, C2::AbstractLinearCode)
    C1.F == C2.F || throw(ArgumentError("Base field must be the same in the Plotkin (u|u + v)-construction."))
    C1.n == C2.n || throw(ArgumentError("Both codes must have the same length in the Plotkin (u|u + v)-construction."))
    (iszero(C1.G) || iszero(C2.G)) && throw(ArgumentError("`u_u_plus_v` not supported for zero codes."))

    G1 = generator_matrix(C1)
    G2 = generator_matrix(C2)
    G = vcat(hcat(G1, G1), hcat(parent(G2)(0), G2))
    H1 = parity_check_matrix(C1)
    H2 = parity_check_matrix(C2)
    H = vcat(hcat(H1, parent(H1)(0)), hcat(-H2, H2))
    G_stand, H_stand, P, k = _standard_form(G)
    k == C1.k + C2.k || error("Something went wrong in the Plotkin (u|u + v)-construction;
        dimension ($k) and expected dimension ($(C1.k + C2.k)) are not equal.")

    if ismissing(C1.d) || ismissing(C2.d)
        lb = min(2 * C1.l_bound, C2.l_bound)
        ub1, _ = _min_wt_row(G)
        ub2, _ = _min_wt_row(G_stand)
        ub = min(ub1, ub2, 2 * C1.u_bound, C2.u_bound)
        return LinearCode(C1.F, 2 * C1.n, k, missing, lb, ub, G, H, G_stand, H_stand, P, missing)
    else
        d = min(2 * C1.d, C2.d)
        return LinearCode(C1.F, 2 * C1.n, k, d, d, d, G, H, G_stand, H_stand, P, missing)
    end
end
Plotkin_construction(C1::AbstractLinearCode, C2::AbstractLinearCode) = u_u_plus_v(C1, C2)

"""
    u_plus_w_v_plus_w_u_plus_v_plus_w(C1::AbstractLinearCode, C2::AbstractLinearCode)

Return the code generated by the (u + w | v + w | u + v + w)-construction.
"""
function u_plus_w_v_plus_w_u_plus_v_plus_w(C1::AbstractLinearCode, C2::AbstractLinearCode)
    C1.F == C2.F || throw(ArgumentError("All codes must be over the same base ring in the (u + w | v + w | u + v + w)-construction."))
    C1.n == C2.n || throw(ArgumentError("All codes must be the same length in the (u + w | v + w | u + v + w)-construction."))

    G1 = generator_matrix(C1)
    G2 = generator_matrix(C2)
    Z = zero_matrix(C1.F, C1.k, C1.n)
    # could do some verification steps on parameters
    return LinearCode(vcat(hcat(G1, Z,  G1),
                           hcat(G2, G2, G2),
                           hcat(Z,  G1, G1)))
end

# TODO: add construction A, B, Y, B2

"""
    construction_X(C1::AbstractLinearCode, C2::AbstractLinearCode, C3::AbstractLinearCode)

Return the code generated by the construction X procedure.
"""
function construction_X(C1::AbstractLinearCode, C2::AbstractLinearCode, C3::AbstractLinearCode)
    C1 ⊆ C2 || throw(ArgumentError("The first code must be a subcode of the second in construction X."))
    are_equivalent(C1, C2) && throw(ArgumentError("The first code must be a proper subcode of the second in construction X."))
    C1.F == C2.F == C3.F || throw(ArgumentError("All codes must be over the same base ring in construction X."))
    C2.k == C1.k + dimension(C3) ||
        throw(ArgumentError("The dimension of the second code must be the sum of the dimensions of the first and third codes."))

    # could do some verification steps on parameters
    C = LinearCode(vcat(hcat(generator_matrix(C1 / C2), generator_matrix(C3)),
        hcat(C1.G, zero_matrix(C1.F, C1.k, length(C3)))))
    C.n == C1.n + C3.n || error("Something went wrong in construction X. Expected length
        $(C1.n + C3.n) but obtained length $(C.n)")

    if !ismissing(C1.d) && !ismissing(C2.d) && !ismissing(C3.d)
        C.d = minimum([C1.d, C2.d + C3.d])
    elseif ismissing(C.d)
        # pretty sure this holds
        set_distance_lower_bound!(C, min(C1.l_bound, C2.l_bound + C3.l_bound))
    end
    return C
end

"""
    construction_X3(C1::AbstractLinearCode, C2::AbstractLinearCode, C3::AbstractLinearCode,
        C4::AbstractLinearCode, C5::AbstractLinearCode))

Return the code generated by the construction X3 procedure.
"""
function construction_X3(C1::AbstractLinearCode, C2::AbstractLinearCode, C3::AbstractLinearCode,
    C4::AbstractLinearCode, C5::AbstractLinearCode)

    C1 ⊆ C2 || throw(ArgumentError("The first code must be a subcode of the second in construction X3."))
    are_equivalent(C1, C2) && throw(ArgumentError("The first code must be a proper subcode of the second in construction X3."))
    C2 ⊆ C3 || throw(ArgumentError("The second code must be a subcode of the third in construction X3."))
    are_equivalent(C2, C3) && throw(ArgumentError("The second code must be a proper subcode of the third in construction X3."))
    # the above lines check C1.F == C2.F == field(C3)
    C1.F == C4.F  == C5.F || throw(ArgumentError("All codes must be over the same base ring in construction X3."))
    C3.k == C2.k + C4.k ||
        throw(ArgumentError("The dimension of the third code must be the sum of the dimensions of the second and fourth codes in construction X3."))
    C2.k == C1.k + C5.k ||
        throw(ArgumentError("The dimension of the second code must be the sum of the dimensions of the first and fifth codes in construction X3."))

    C2_mod_C1 = C2 / C1
    C3_mod_C2 = C3 / C2
    F = C1.F

    # could do some verification steps on parameters
    G = vcat(hcat(generator_matrix(C1), zero_matrix(F, C1.k, C4.n), zero_matrix(F, C1.k, C5.n)),
    hcat(C2_mod_C1.G, C4.G, zero_matrix(F, C4.k, C5.n)),
    hcat(C3_mod_C2.G, zero_matrix(F, C5.k, C4.n), generator_matrix(C5)))
    C = LinearCode(G)
    if !ismissing(C1.d) && !ismissing(C2.d) && !ismissing(C3.d) && !ismissing(C4.d) && !ismissing(C5.d) && ismissing(C.d)
        setlowerdistancebound!(C, minimum([C1.d, C2.d + C4.d, C3.d + C5.d]))
    end
    return C
end

"""
    ⊕(C1::AbstractLinearCode, C2::AbstractLinearCode)
    direct_sum(C1::AbstractLinearCode, C2::AbstractLinearCode)

Return the direct sum code of `C1` and `C2`.
"""
function ⊕(C1::AbstractLinearCode, C2::AbstractLinearCode)
    C1.F == C2.F || throw(ArgumentError("Codes must be over the same field."))

    G1 = generator_matrix(C1)
    G2 = generator_matrix(C2)
    G = direct_sum(G1, G2)
    H1 = parity_check_matrix(C1)
    H2 = parity_check_matrix(C2)
    H = direct_sum(H1, H2)
    # ordering is tougher for standard form, easiest to just recompute:
    G_stand, H_stand, P, k = _standard_form(G)
    k == C1.k + C2.k || error("Unexpected dimension in direct sum output.")

    if !ismissing(C1.d) && !ismissing(C2.d)
        d = minimum([C1.d, C2.d])
        return LinearCode(C1.F, C1.n, k, d, d, d, G, H, G_stand, H_stand, P, missing)
    else
        lb = minimum([C1.l_bound, C2.l_bound])
        ub = minimum([C1.u_bound, C2.u_bound])
        return LinearCode(C1.F, C1.n, k, missing, lb, ub, G, H, G_stand, H_stand, P, missing)
    end
end
direct_sum(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊕ C2

"""
    ×(C1::AbstractLinearCode, C2::AbstractLinearCode)
    direct_product(C1::AbstractLinearCode, C2::AbstractLinearCode)
    product_code(C1::AbstractLinearCode, C2::AbstractLinearCode)

Return the (direct) product code of `C1` and `C2`.
"""
function ×(C1::AbstractLinearCode, C2::AbstractLinearCode)
    C1.F == C2.F || throw(ArgumentError("Codes must be over the same field."))

    G = generator_matrix(C1) ⊗ generator_matrix(C2)
    G_stand, H_stand, P, k = _standard_form(G)
    k == C1.k * C2.k || error("Unexpected dimension in direct product output.")
    H = ismissing(P) ? deepcopy(H_stand) : H_stand * P

    if !ismissing(C1.d) && !ismissing(C2.d)
        d = C1.d * C2.d
        return LinearCode(C1.F, C1.n * C2.n, k, d, d, d, G, H, G_stand, H_stand, P, missing)
    else
        return LinearCode(C1.F, C1.n * C2.n, k, missing, C1.l_bound * C2.l_bound,
            C1.u_bound * C2.u_bound, G, H, G_stand, H_stand, P, missing)

    end
end
direct_product(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 × C2
product_code(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 × C2

"""
    ⊗(C1::AbstractLinearCode, C2::AbstractLinearCode)
    kron(C1::AbstractLinearCode, C2::AbstractLinearCode)
    tensor_product(C1::AbstractLinearCode, C2::AbstractLinearCode)
    
Return the tensor product code of `C1` and `C2`.
"""
⊗(C1::AbstractLinearCode, C2::AbstractLinearCode) = LinearCode(parity_check_matrix(C1) ⊗ parity_check_matrix(C2), true)
kron(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊗ C2
tensor_product(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊗ C2

# R.Pellikaan, On decoding by error location and dependent sets of error
# positions, Discrete Mathematics, 106107 (1992), 369-381.
# the Schur product of vector spaces is highly basis dependent and is often the
# full vector space (an [n, n, 1] code)
# might be a cyclic code special case in
# "On the Schur Product of Vector Spaces over Finite Fields"
# Christiaan Koster
"""
    entrywise_product_code(C::AbstractLinearCode, D::AbstractLinearCode)
    *(C::AbstractLinearCode, D::AbstractLinearCode)
    Schur_product_code(C::AbstractLinearCode, D::AbstractLinearCode)
    Hadamard_product_code(C::AbstractLinearCode, D::AbstractLinearCode)
    componentwise_product_code(C::AbstractLinearCode, D::AbstractLinearCode)

Return the entrywise product of `C` and `D`.
"""
function entrywise_product_code(C::AbstractLinearCode, D::AbstractLinearCode)
    C.F == D.F || throw(ArgumentError("Codes must be over the same field in the Schur product."))
    C.n == D.n || throw(ArgumentError("Codes must have the same length in the Schur product."))
    if isa(C, ReedMullerCode)

        r = C.r + D.r
        if r <= C.n
            return ReedMullerCode(r, C.m)
        else
            return ReedMullerCode(C.m, C.m)
        end
    else
        G_C = generator_matrix(C)
        G_D = generator_matrix(D)
        nr_C = nrows(G_C)
        nr_D = nrows(G_D)

        # TODO: Oscar doesn't work well with dot operators
        # TODO: I think this choice of indices only works for C == D
        indices = Vector{Tuple{Int, Int}}()
        for i in 1:nr_C
            for j in 1:nr_D
                i <= j && push!(indices, (i, j))
            end
        end
        return LinearCode(matrix(C.F, reduce(vcat, G_C[i, :] .* G_D[j, :] for (i, j) in indices)))
    end
    # TODO: verify C ⊂ it?
end
*(C::AbstractLinearCode, D::AbstractLinearCode) = entrywise_product_code(C, D)
Schur_product_code(C::AbstractLinearCode, D::AbstractLinearCode) = entrywise_product_code(C, D)
Hadamard_product_code(C::AbstractLinearCode, D::AbstractLinearCode) = entrywise_product_code(C, D)
componentwise_product_code(C::AbstractLinearCode, D::AbstractLinearCode) = entrywise_product_code(C, D)

"""
    /(C1::AbstractLinearCode, C2::AbstractLinearCode)
    quo(C1::AbstractLinearCode, C2::AbstractLinearCode)
    quotient(C1::AbstractLinearCode, C2::AbstractLinearCode)
    code_complement(C2::AbstractLinearCode, C1::AbstractLinearCode)

Return the code `C2 / C1` given `C1 ⊆ C2`.
"""
function /(C2::AbstractLinearCode, C1::AbstractLinearCode)
    C1 ⊆ C2 || throw(ArgumentError("C1 ⊈ C2"))
    F = C1.F
    G1 = generator_matrix(C1)
    G2 = generator_matrix(C2)
    V = vector_space(F, C1.n)
    U, U_to_V = sub(V, [V(G1[i, :]) for i in 1:nrows(G1)])
    W, W_to_V = sub(V, [V(G2[i, :]) for i in 1:nrows(G2)])
    gens_of_U_in_W = [preimage(W_to_V, U_to_V(g)) for g in gens(U)]
    U_in_W, _ = sub(W, gens_of_U_in_W)
    Q, W_to_Q = quo(W, U_in_W)
    C2_mod_C1_basis = [W_to_V(x) for x in [preimage(W_to_Q, g) for g in gens(Q)]]
    F_basis = [[F(C2_mod_C1_basis[j][i]) for i in 1:dim(parent(C2_mod_C1_basis[1]))] for j in 1:length(C2_mod_C1_basis)]
    G = matrix(F, length(F_basis), length(F_basis[1]), reduce(vcat, F_basis))
    for r in 1:length(F_basis)
        v = G[r, :]
        (v ∈ C2 && v ∉ C1) || error("Error in creation of basis for C2 / C1.")
    end
    return LinearCode(G)
end
quo(C1::AbstractLinearCode, C2::AbstractLinearCode) = /(C1, C2)
quotient(C1::AbstractLinearCode, C2::AbstractLinearCode) = /(C1, C2)
code_complement(C2::AbstractLinearCode, C1::AbstractLinearCode) = /(C1, C2)

"""
    juxtaposition(C1::AbstractLinearCode, C2::AbstractLinearCode)

Return the code generated by the horizontal concatenation of the generator
matrices of `C1` then `C2`.
"""
function juxtaposition(C1::AbstractLinearCode, C2::AbstractLinearCode)
    C1.F == C2.F || throw(ArgumentError("Cannot juxtapose two codes over different fields."))
    G1 = generator_matrix(C1)
    G2 = generator_matrix(C2)
    nrows(G1) == nrows(G2) || throw(ArgumentError("Cannot juxtapose codes with generator matrices with a different number of rows."))

    return LinearCode(hcat(G1, G2))
end

#############################
      # getter functions
#############################

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

# keep not as a one-linear for a moment to add new properties later
"""

Return the transposed code of `C`.
"""
transpose(C::AbstractLinearCode) = LinearCode(transpose(parity_check_matrix(C)), true)

"""
    permute_code(C::AbstractLinearCode, σ::Union{PermGroupElem, Perm{Int}, Vector{Int}})

Return the code whose generator matrix is `C`'s with the columns permuted by `σ`.
"""
function permute_code(C::AbstractLinearCode, σ::Union{PermGroupElem, Perm{Int}, Vector{Int}})
    if isa(C, QuasiCyclicCode)
        P = transpose(permutation_matrix(C.F, typeof(σ) <: Perm ? σ.d : σ))
        size(P, 1) == C.n || throw(ArgumentError("Incorrect number of digits in permutation."))
        G = generator_matrix(C) * P
        H = parity_check_matrix(C) * P
        C2 = LinearCode(G, H)
        ismissing(C2.d) && !ismissing(C.d) && set_minimum_distance!(C2, C.d)
        return C2
    else
        C2 = copy(C)
        P = transpose(permutation_matrix(C.F, typeof(σ) <: Perm ? σ.d : σ))
        size(P, 1) == C.n || throw(ArgumentError("Incorrect number of digits in permutation."))
        C2.G = C2.G * P
        C2.H = C2.H * P
        C2.P_stand = ismissing(C2.P_stand) ? P : C2.P_stand * P
        return C2
    end
end

"""
    extend(C::AbstractLinearCode, a::CTMatrixTypes, c::Integer)
    extend(C::AbstractLinearCode, c::Integer)
    extend(C::AbstractLinearCode, a::CTMatrixTypes)
    extend(C::AbstractLinearCode)
    even_extension(C::AbstractLinearCode)

Return the extended code of `C` extending on column `c`. For each row
`g` of the generator matrix for `C`, a digit `-a ⋅ g` is inserted in
the `c`th position. If `c` isn't given, it is appended. If `a` isn't
given, then the all 1's vector is used giving an even extension.
"""
function extend(C::AbstractLinearCode, a::CTMatrixTypes, c::Integer)
    1 <= c <= C.n + 1 || throw(ArgumentError("The code has length $(C.n), so `c` must be between 1 and $(C.n + 1)."))
    length(a) == C.n || throw(ArgumentError("The vector `a` should have length $(C.n)."))
    b = ncols(a) == 1 ? transpose(a) : a
    nrows(b) == 1 || throw(ArgumentError("The argument `a` should be a vector."))

    F = base_ring(C.G)
    new_col = zero_matrix(F, nrows(C.G), 1)
    is_binary = order(F) == 2
    for i in axes(new_col, 1)
        new_col[i, 1] = dot(b, view(C.G, i:i, :))
        is_binary || (new_col[i, 1] *= F(-1);)
    end
    G_new = hcat(view(C.G, :, 1:c - 1), new_col, view(C.G, :, c:C.n))
    new_row = hcat(view(b, 1:1, 1:c - 1), matrix(F, 1, 1, [1]), view(b, 1:1, c:C.n))
    H_new = vcat(new_row, hcat(view(C.H, :, 1:c - 1), zero_matrix(F, nrows(C.H), 1), view(C.H, :, c:C.n)))
    C_new = LinearCode(G_new, H_new)
    if !ismissing(C.d) && ismissing(C_new.d)
        if is_binary && all(isone(x) for x in a)
            C_new.d = iseven(C.d) ? C.d : C.d + 1
            C_new.l_bound = C_new.u_bound = C_new.d
        else
            C_new.l_bound = C.d
            C_new.u_bound = C.d + 1
        end
    end

    return C_new
end
extend(C::AbstractLinearCode, c::Integer) = extend(C, matrix(base_ring(C.G), 1, C.n, ones(Int, C.n)), c)
extend(C::AbstractLinearCode, a::CTMatrixTypes) = extend(C, a, C.n + 1)
extend(C::AbstractLinearCode) = extend(C, matrix(base_ring(C.G), 1, C.n, ones(Int, C.n)), C.n + 1)
even_extension(C::AbstractLinearCode) = extend(C)

"""
    puncture(C::AbstractLinearCode, cols::Vector{<:Integer})
    puncture(C::AbstractLinearCode, cols::Integer)

Return the code of `C` punctured at the columns in `cols`.
"""
function puncture(C::AbstractLinearCode, cols::Vector{<:Integer})
    isempty(cols) && return C
    allunique(cols) || throw(ArgumentError("Columns to puncture are not unique."))
    cols ⊆ 1:C.n || throw(ArgumentError("Columns to puncture are not a subset of the index set."))
    length(cols) == C.n && throw(ArgumentError("Cannot puncture all columns of a generator matrix."))

    G = generator_matrix(C)[:, setdiff(1:C.n, cols)]
    G = _remove_empty(G, :rows)
    G_stand, H_stand, P, k = _standard_form(G)
    H = ismissing(P) ? deepcopy(H_stand) : H_stand * P

    ub1, _ = _min_wt_row(G)
    ub2, _ = _min_wt_row(G_stand)
    ub = C.l_bound > 1 ? min(ub1, ub2, C.u_bound) : min(ub1, ub2)
    lb = max(1, C.l_bound - 1)
    return LinearCode(C.F, ncols(G), k, missing, lb, ub, G, H, G_stand, H_stand, P, missing)
end
puncture(C::AbstractLinearCode, cols::Integer) = puncture(C, [cols])

"""
    expurgate(C::AbstractLinearCode, rows::Vector{<:Integer})
    expurgate(C::AbstractLinearCode, rows::Integer)

Return the code of `C` expuragated at the rows in `rows`.
"""
function expurgate(C::AbstractLinearCode, rows::Vector{<:Integer})
    isempty(rows) && return C
    allunique(rows) || throw(ArgumentError("Rows to expurgate are not unique."))
    G = generator_matrix(C)
    nr = nrows(G)
    rows ⊆ 1:nr || throw(ArgumentError("Rows to expurgate are not a subset of the index set."))
    length(rows) == nr && throw(ArgumentError("Cannot expurgate all rows of a generator matrix."))

    G = G[setdiff(1:nr, rows), :]
    G_stand, H_stand, P, k = _standard_form(G)
    H = ismissing(P) ? deepcopy(H_stand) : H_stand * P

    ub1, _ = _min_wt_row(G)
    ub2, _ = _min_wt_row(G_stand)
    ub = min(ub1, ub2)
    return LinearCode(C.F, C.n, k, missing, C.l_bound, ub, G, H, G_stand, H_stand, P, missing)
end
expurgate(C::AbstractLinearCode, rows::Integer) = expurgate(C, [rows])

"""
    shorten(C::AbstractLinearCode, L::Vector{<:Integer})
    shorten(C::AbstractLinearCode, L::Integer)

Return the code of `C` shortened on the indices `L`.
"""
shorten(C::AbstractLinearCode, L::Vector{<:Integer}) = isempty(L) ? (return C;) : (return dual(puncture(dual(C), L));)
shorten(C::AbstractLinearCode, L::Integer) = shorten(C, [L])

"""
    augment(C::AbstractLinearCode, M::CTMatrixTypes)

Return the code of `C` whose generator matrix is augmented with `M`.
"""
function augment(C::AbstractLinearCode, M::CTMatrixTypes)
    iszero(M) && throw(ArgumentError("Zero matrix passed to augment."))
    C.n == ncols(M) || throw(ArgumentError("Rows to augment must have the same number of columns as the generator matrix."))
    C.F == base_ring(M) || throw(ArgumentError("Rows to augment must have the same base field as the code."))

    M = _remove_empty(M, :rows)
    G = vcat(generator_matrix(C), M)
    G_stand, H_stand, P, k = _standard_form(G)
    H = ismissing(P) ? deepcopy(H_stand) : H_stand * P

    ub1, _ = _min_wt_row(G)
    ub2, _ = _min_wt_row(G_stand)
    ub = min(ub1, ub2, C.u_bound)
    return LinearCode(C.F, C.n, k, missing, 1, ub, G, H, G_stand, H_stand, P, missing)
end
augment(C::AbstractLinearCode, M::Matrix{<:Integer}) = augment(C, matrix(C.F, M))

"""
    lengthen(C::AbstractLinearCode)

Return the lengthened code of `C`.
"""
lengthen(C::AbstractLinearCode) = extend(augment(C, matrix(C.F, ones(Int, 1, C.n))))

"""
    subcode(C::AbstractLinearCode, k::Int)

Return a `k`-dimensional subcode of `C`.
"""
function subcode(C::AbstractLinearCode, k::Int)
    k >= 1 && k < C.k || throw(ArgumentError("Cannot construct a $k-dimensional subcode of an $(C.k)-dimensional code."))

    k != C.k || return C
    return if ismissing(C.P_stand)
        LinearCode(view(generator_matrix(C, true), 1:k, :))
    else
        LinearCode(view(generator_matrix(C, true), 1:k, :) * C.P_stand)
    end
end

"""
    subcode(C::AbstractLinearCode, rows::Vector{Int})

Return a subcode of `C` using the rows of the generator matrix of `C` listed in
`rows`.
"""
function subcode(C::AbstractLinearCode, rows::Vector{Int})
    isempty(rows) && throw(ArgumentError("Row index set empty in subcode."))
    allunique(rows) || throw(ArgumentError("Rows are not unique."))
    rows ⊆ 1:C.k || throw(ArgumentError("Rows are not a subset of the index set."))

    length(rows) == C.k && return C
    return LinearCode(generator_matrix(C)[setdiff(1:C.k, rows), :])
end

"""
    subcode_of_dimension_between_codes(C1::AbstractLinearCode, C2::AbstractLinearCode, k::Int)

Return a subcode of dimenion `k` between `C1` and `C2`.
"""
function subcode_of_dimension_between_codes(C1::AbstractLinearCode, C2::AbstractLinearCode, k::Int)
    C2 ⊆ C1 || throw(ArgumentError("C2 must be a subcode of C1"))
    C2.k <= k <= C1.k || throw(ArgumentError("The dimension must be between that of C1 and C2."))
    
    k == C2.k && return C2
    k == C1.k && return C1
    C = C1 / C2
    return if ismissing(C.P_stand)
        augment(C2, generator_matrix(C, true)[1:k - C2.k, :])
    else
        augment(C2, generator_matrix(C, true)[1:k - C2.k, :] * C.P_stand)
    end
end

"""
    expanded_code(C::AbstractLinearCode, K::CTFieldTypes, β::Vector{FqFieldElem})

Return the expanded code of `C` constructed by exapnding the generator matrix
to the subfield `K` using the basis `β` for `field(C)` over `K`.
"""
function expanded_code(C::AbstractLinearCode, K::CTFieldTypes, β::Vector{FqFieldElem})
    flag, λ = is_basis(C.F, K, β)
    flag || throw(ArgumentError("β is not a basis for the extension"))

    # D = _expansion_dict(C.F, K, λ)
    # m = div(degree(C.F), degree(K))

    G = generator_matrix(C)
    G_new = zero_matrix(C.F, nrows(G) * length(β), ncols(G))
    curr_row = 1
    for r in 1:nrows(G)
        for βi in β
            G_new[curr_row, :] = βi * view(G, r:r, :)
            curr_row += 1
        end
    end
    # G_exp = _expandmatrix(G_new, D, m)
    G_exp = expand_matrix(G_new, K, β)

    H = parity_check_matrix(C)
    H_new = zero_matrix(C.F, nrows(H) * length(β), ncols(H))
    curr_row = 1
    for r in 1:nrows(H)
        for βi in β
            H_new[curr_row, :] = βi * view(H, r:r, :)
            curr_row += 1
        end
    end
    H_exp = expand_matrix(H_new, K, λ)
    # H_exp = _expandmatrix(H_new, D, m)
    return G_exp

    C_new = LinearCode(G_exp, H_exp)
    C_new.G = change_base_ring(K, C_new.G)
    C_new.H = change_base_ring(K, C_new.H)
    C_new.G_stand = change_base_ring(K, C_new.G_stand)
    C_new.H_stand = change_base_ring(K, C_new.H_stand)
    ismissing(C_new.P_stand) || (C_new.P_stand = change_base_ring(K, new.P_stand);)
    C_new.F = K
    C_new.l_bound = C.l_bound
    C_new.u_bound, _ = _min_wt_row(C_new.G)
    return C_new
end

# """
#     subfield_subcode(C::AbstractLinearCode, K::fqPolyRepField)
#
# Return the subfield subcode code of `C` over `K` using Delsarte's theorem.
#
# Use this method if you are unsure of the dual basis to the basis you which
# to expand with.
# """
# function subfield_subcode(C::AbstractLinearCode, K::fqPolyRepField)
#     return dual(trace_code(dual(C), K))
# end

"""
    subfield_subcode(C::AbstractLinearCode, K::CTFieldTypes, basis::Vector{<:CTFieldElem})

Return the subfield subcode code of `C` over `K` using the provided dual `basis`
for the field of `C` over `K`.
"""
function subfield_subcode(C::AbstractLinearCode, K::CTFieldTypes, basis::Vector{<:CTFieldElem})
    C_new = LinearCode(transpose(expand_matrix(transpose(parity_check_matrix(C)), K, basis)), true)
    C_new.G = change_base_ring(K, C_new.G)
    C_new.H = change_base_ring(K, C_new.H)
    C_new.G_stand = change_base_ring(K, C_new.G_stand)
    C_new.H_stand = change_base_ring(K, C_new.H_stand)
    ismissing(C_new.P_stand) || (C_new.P_stand = change_base_ring(K, new.P_stand);)
    C_new.F = K
    return C_new
end

"""
    trace_code(C::AbstractLinearCode, K::CTFieldTypes, basis::Vector{<:CTFieldElem})

Return the trace code of `C` over `K` using the provided dual `basis`
for the field of `C` over `K` using Delsarte's theorem.
"""
trace_code(C::AbstractLinearCode, K::CTFieldTypes, basis::Vector{<:CTFieldElem}) = dual(subfield_subcode(dual(C), K, basis))

# needs significant testing, works so far
"""
    even_subcode(C::AbstractLinearCode)

Return the even subcode of `C`.
"""
function even_subcode(C::AbstractLinearCode)
    F = C.F
    Int(order(F)) == 2 || throw(ArgumentError("Even-ness is only defined for binary codes."))

    V_C, ψ = vector_space(C)
    G_F_VS = vector_space(F, 1)
    homo1_quad = ModuleHomomorphism(V_C, G_F_VS, matrix(F, dim(V_C), 1,
        reduce(vcat, [F(weight(ψ(g).v) % 2) for g in gens(V_C)])))
    even_sub, ϕ1 = kernel(homo1_quad)
    iszero(dim(even_sub)) ? (return missing;) : (return LinearCode(reduce(vcat, [ψ(ϕ1(g)).v for g in gens(even_sub)]));)
end

# needs significant testing, works so far
"""
    doubly_even_subcode(C::AbstractLinearCode)

Return the doubly-even subcode of `C`.
"""
function doubly_even_subcode(C::AbstractLinearCode)
    F = C.F
    Int(order(C.F)) == 2 || throw(ArgumentError("Even-ness is only defined for binary codes."))

    V_C, ψ = vector_space(C)
    G_F_VS = vector_space(F, 1)

    # first get the even subspace
    homo1_quad = ModuleHomomorphism(V_C, G_F_VS, matrix(F, dim(V_C), 1,
        reduce(vcat, [F(weight(ψ(g).v) % 2) for g in gens(V_C)])))
    even_sub, ϕ1 = kernel(homo1_quad)

    if !iszero(dim(even_sub))
        # now control the overlap (Ward's divisibility theorem)
        homo2_bi = ModuleHomomorphism(even_sub, even_sub, matrix(F, dim(even_sub),
            dim(even_sub),
            reduce(vcat, [F(weight(matrix(F, 1, C.n, ψ(ϕ1(gens(even_sub)[i])).v .*
                ψ(ϕ1(gens(even_sub)[j])).v)) % 2)
            for i in 1:dim(even_sub), j in 1:dim(even_sub)])))
        even_sub_w_overlap, μ1 = kernel(homo2_bi)

        if !iszero(dim(even_sub_w_overlap))
            # now apply the weight four condition
            homo2_quad = ModuleHomomorphism(even_sub_w_overlap, G_F_VS, matrix(F,
                dim(even_sub_w_overlap), 1, reduce(vcat, [F(div(weight(ψ(ϕ1(μ1(g))).v),
                2) % 2) for g in gens(even_sub_w_overlap)])))
            four_sub, ϕ2 = kernel(homo2_quad)

            if !iszero(dim(four_sub))
                return LinearCode(reduce(vcat, [ψ(ϕ1(μ1(ϕ2(g)))).v for g in gens(four_sub)]))
            end
        end
    end
    return missing
end

# currently incomplete
# function triply_even_subcode(C::AbstractLinearCode)
#     F = C.F
#     V_C, ψ = vector_space(C)
#     G_F_VS = vector_space(F, 1)
#
#     # first get the even subspace
#     homo1_quad = ModuleHomomorphism(V_C, G_F_VS, matrix(F, dim(V_C), 1, vcat([F(weight(ψ(g).v) % 2) for g in gens(V_C)]...)))
#     even_sub, ϕ1 = kernel(homo1_quad)
#
#     # now control the overlap (Ward's divisibility theorem)
#     homo2_bi = ModuleHomomorphism(even_sub, even_sub, matrix(F, dim(even_sub), dim(even_sub),
#         vcat([F(weight(matrix(F, 1, C.n, ψ(ϕ1(gens(even_sub)[i])).v .* ψ(ϕ1(gens(even_sub)[j])).v)) % 2)
#         for i in 1:dim(even_sub), j in 1:dim(even_sub)]...)))
#     even_sub_w_overlap, μ1 = kernel(homo2_bi)
#
#     # now apply the weight four condition
#     homo2_quad = ModuleHomomorphism(even_sub_w_overlap, G_F_VS, matrix(F, dim(even_sub_w_overlap), 1,
#         vcat([F(div(weight(ψ(ϕ1(μ1(g))).v), 2) % 2) for g in gens(even_sub_w_overlap)]...)))
#     four_sub, ϕ2 = kernel(homo2_quad)
#
#     # now control the triple overlap
#     ###########
#     # how do we actually represent this?
#
#
#     # now apply the weight eight condition
#     homo3 = ModuleHomomorphism(four_sub, G_F_VS, matrix(F, dim(four_sub), 1, vcat([F(div(weight(ψ(ϕ1(ϕ2(g))).v), 4) % 2) for g in gens(four_sub)]...)))
#     eight_sub, ϕ3 = kernel(homo3)
#     # # println(dim(eight_sub))
#     # # println(vcat([ψ(ϕ1(ϕ2(ϕ3(g)))).v for g in gens(eight_sub)]...))
#     #
#     # q, γ = quo(four_sub, eight_sub)
#     # println(dim(q))
#     # # println([preimage(γ, g).v for g in gens(q)])
#     # return LinearCode(vcat([ψ(ϕ1(ϕ2(preimage(γ, g)))).v for g in gens(q)]...))
#     # # return LinearCode(vcat([ψ(ϕ1(ϕ2(ϕ3(g)))).v for g in gens(eight_sub)]...))
#
#     return LinearCode(vcat([ψ(ϕ1(μ1(ϕ2(g)))).v for g in gens(four_sub)]...))
# end
