# Copyright (c) 2021 - 2024 Eric Sabo, Benjamin Ide, David Marquis
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
  # Generic Helper Functions
#############################

"""
    copy(C::T) where T <: AbstractCode

Returns a copy of the code `C`.
"""
function copy(C::T) where T <: AbstractCode
    C2 = deepcopy(C)
    hasfield(T, :F) && (C2.F = C.F;)
    hasfield(T, :E) && (C2.E = C.E;)
    hasfield(T, :R) && (C2.R = C.R;)
    hasfield(T, :C) && (C2.C = copy(C.C);)
    if hasfield(T, :Cvec)
        for i in eachindex(C.Cvec)
            C2.C[i] = copy(C.C[i])
        end
    end
    return C2
end

_has_equivalent_row_spaces(A::CTMatrixTypes, B::CTMatrixTypes) =
    return _remove_empty(rref(deepcopy(A))[2], :rows) == _remove_empty(rref(deepcopy(B))[2], :rows)

# """
#     reverse(v::CTMatrixTypes)
#     reverse!(v::CTMatrixTypes)

# Return the reverse of the vector `v`.
# """
# function reverse(v::CTMatrixTypes)
#     nr, nc = size(v)
#     if nr == 1
#         return v[:, end:-1:1]
#     elseif nc == 1
#         return v[end:-1:1, :]
#     else
#         throw(ArgumentError("Matrix must be a vector."))
#     end
# end

# function reverse!(v::CTMatrixTypes)
#     nr, nc = size(v)
#     if nr == 1
#         v[:, 1:end] = v[:, end:-1:1]
#     elseif nc == 1
#         v[1:end, :] = v[end:-1:1, :]
#     else
#         throw(ArgumentError("Matrix must be a vector."))
#     end
# end

# """
#     circshift(v::CTMatrixTypes, l::Int)

# Return the circular shift of the vector `v` by `l` bits.

# This is an overload of Base.circshift for type `CTMatrixTypes`.
# Either the number of rows or the number of columns must have dimension one.
# """
# function circshift(v::CTMatrixTypes, l::Int)
#     nr, nc = size(v)
#     if nr == 1
#         l = l % nc
#         l < 0 && (l = nc + l;)
#         vshift = zero_matrix(base_ring(v), 1, nc)
#         vshift[1, 1:l] = v[1, nc - l + 1:nc]
#         vshift[1, l + 1:end] = v[1, 1:nc - l]
#     elseif nc == 1
#         l = l % nr
#         l < 0 && (l = nr + l;)
#         vshift = zero_matrix(base_ring(v), nr, 1)
#         vshift[1:l, 1] = v[nr - l + 1:nr, 1]
#         vshift[l + 1:end, 1] = v[1:nr - l, 1]
#     else
#         throw(ArgumentError("Input matrix must be a vector."))
#     end
#     return vshift
# end

"""
    ⊕(A::CTMatrixTypes, B::CTMatrixTypes)
    direct_sum(A::CTMatrixTypes, B::CTMatrixTypes)

Return the direct sum of the two matrices `A` and `B`.
"""
function ⊕(A::T, B::T) where T <: CTMatrixTypes
    base_ring(A) == base_ring(B) || throw(ArgumentError("Matrices must be over the same base ring in direct_sum."))

    return vcat(hcat(A, zero_matrix(base_ring(B), nrows(A), ncols(B))),
        hcat(zero_matrix(base_ring(A), nrows(B), ncols(A)), B))
end
direct_sum(A::T, B::T) where T <: CTMatrixTypes = A ⊕ B

"""
    ⊗(A::CTMatrixTypes, B::CTMatrixTypes)
    kron(A::CTMatrixTypes, B::CTMatrixTypes)
    tensor_product(A::CTMatrixTypes, B::CTMatrixTypes)
    kronecker_product(A::CTMatrixTypes, B::CTMatrixTypes)

Return the Kronecker product of the two matrices `A` and `B`.
"""
⊗(A::Union{CTMatrixTypes, MatElem{<: ResElem}, MatElem{<: CTGroupAlgebra}}, B::Union{CTMatrixTypes, MatElem{<: ResElem}, MatElem{<: CTGroupAlgebra}}) = kronecker_product(A, B)
kron(A::Union{CTMatrixTypes, MatElem{<: ResElem}, MatElem{<: CTGroupAlgebra}}, B::Union{CTMatrixTypes, MatElem{<: ResElem}, MatElem{<: CTGroupAlgebra}}) = kronecker_product(A, B)
tensor_product(A::Union{CTMatrixTypes, MatElem{<: ResElem}, MatElem{<: CTGroupAlgebra}}, B::Union{CTMatrixTypes, MatElem{<: ResElem}, MatElem{<: CTGroupAlgebra}}) = kronecker_product(A, B)

# I think we should avoid length checking here and return it for entire matrix if given
# Hamming_weight(v::T) where T <: Union{CTMatrixTypes, gfp_mat, Vector{S}} where S <: Integer = count(i->(i != 0), v)
"""
    Hamming_weight(v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer
    weight(v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer
    wt(v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer

Return the Hamming weight of `v`.
"""
function Hamming_weight(v::T) where T <: Union{CTMatrixTypes, Vector{<:CTFieldElem}, Vector{S}, Adjoint{S, Vector{S}}, AbstractMatrix{S}} where S <: Integer
    count(x -> !iszero(x), v)
end
weight(v::T) where T <: Union{CTMatrixTypes, Vector{<:CTFieldElem}, Vector{S}, Adjoint{S, Vector{S}}, AbstractMatrix{S}} where S <: Integer = Hamming_weight(v)
wt(v::T) where T <: Union{CTMatrixTypes, Vector{<:CTFieldElem}, Vector{S}, Adjoint{S, Vector{S}}, AbstractMatrix{S}} where S <: Integer = Hamming_weight(v)

# TODO polish and export?
function row_wts_symplectic(A::Union{CTMatrixTypes, Matrix{<: Integer}, Matrix{Bool}})
    nc = size(A, 2)
    iseven(nc) || throw(ArgumentError("Input does not have even length"))
    n = div(nc, 2)
    return [count(!iszero(A[i, j]) || !iszero(A[i, j + n]) for j in 1:n) for i in axes(A, 1)]
end

# Hamming_weight(v::Matrix{Int}) = return sum(v)

# wt(v::Matrix{Int}) = count(x -> !iszero(x), v)

"""
    wt(f::CTPolyRingElem)

Return the number of nonzero coefficients of the polynomial `f`.
"""
wt(f::CTPolyRingElem) = Hamming_weight(collect(coefficients(f)))

# """
#     _min_wt_row(A::CTMatrixTypes)
#
# Return the minimum weight and corresponding index of the rows of `A`.
# """
function _min_wt_row(A::Union{CTMatrixTypes, Matrix{S}, LinearAlgebra.Adjoint{S, Matrix{S}}}) where S <: Integer
    w = size(A, 2) + 1
    i = 0
    for r in axes(A, 1)
        w_loc = 0
        for c in axes(A, 2)
            iszero(A[r,c]) || (w_loc += 1)
        end
        w_loc < w && (w = w_loc; i = r;)
    end
    return w, i
end

function _min_wt_col(A::Union{CTMatrixTypes, Matrix{S}, LinearAlgebra.Adjoint{S, Matrix{S}}}) where S <: Integer
    w = size(A, 1) + 1
    i = 0
    for c in axes(A, 2)
        w_loc = 0
        for r in axes(A, 1)
            iszero(A[r,c]) || (w_loc += 1)
        end
        w_loc < w && (w = w_loc; i = c;)
    end
    return w, i
end

"""
    Hamming_distance(u::T, v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer
    distance(u::T, v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer
    dist(u::T, v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer

Return the Hamming distance between `u` and `v`.
"""
Hamming_distance(u::T, v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer = Hamming_weight(u - v)
distance(u::T, v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer = Hamming_weight(u - v)
dist(u::T, v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer = Hamming_weight(u - v)

"""
    symplectic_inner_product(u::CTMatrixTypes, v::CTMatrixTypes)

Return the symplectic inner product of `u` and `v`.
"""
function symplectic_inner_product(u::CTMatrixTypes, v::CTMatrixTypes)
    (nrows(u) == 1 || ncols(u) == 1) || throw(ArgumentError("First argument of symplectic inner product is not a vector: dims = $(size(u, 1))."))
    (nrows(v) == 1 || ncols(v) == 1) || throw(ArgumentError("Second argument of symplectic inner product is not a vector: dims = $(size(v, 1))."))
    length(u) == length(v) || throw(ArgumentError("Vectors must be the same length in symplectic inner product."))
    iseven(length(u)) || throw(ArgumentError("Vectors must have even length in symplectic inner product."))
    base_ring(u) == base_ring(v) || throw(ArgumentError("Vectors must be over the same field in symplectic inner product."))
    
    ncols = div(length(u), 2)
    return sum(u[i + ncols] * v[i] - v[i + ncols] * u[i] for i in 1:ncols)
end

"""
    are_symplectic_orthogonal(A::CTMatrixTypes, B::CTMatrixTypes)

Return `true` if the rows of the matrices `A` and `B` are symplectic orthogonal.
"""
function are_symplectic_orthogonal(A::CTMatrixTypes, B::CTMatrixTypes)
    base_ring(A) == base_ring(B) || throw(ArgumentError("Matrices in product must both be over the same base ring."))
    
    return iszero(hcat(A[:, div(ncols(A), 2) + 1:end], -A[:, 1:div(ncols(A), 2)]) * transpose(B))
end

# function traceinnerproduct(u::CTMatrixTypes, v::CTMatrixTypes)
#
# end

"""
    Hermitian_inner_product(u::CTMatrixTypes, v::CTMatrixTypes)

Return the Hermitian inner product of `u` and `v`.
"""
function Hermitian_inner_product(u::CTMatrixTypes, v::CTMatrixTypes)
    (nrows(u) == 1 || ncols(u) == 1) || throw(ArgumentError("First argument of Hermitian inner product is not a vector: dims = $(size(u, 1))."))
    (nrows(v) == 1 || ncols(v) == 1) || throw(ArgumentError("Second argument of Hermitian inner product is not a vector: dims = $(size(v, 1))."))
    length(u) == length(v) || throw(ArgumentError("Vectors must be the same length in Hermitian inner product."))
    base_ring(u) == base_ring(v) || throw(ArgumentError("Vectors must be over the same field in Hermitian inner product."))
    q2 = order(base_ring(u))
    is_square(q2) || throw(ArgumentError("The Hermitian inner product is only defined over quadratic field extensions."))
    
    q = Int(sqrt(q2, check = false))
    return sum(u[i] * v[i]^q for i in 1:length(u))
end

"""
    Hermitian_conjugate_matrix(A::CTMatrixTypes)

Return the Hermitian conjugate of the matrix `A`.
"""
function Hermitian_conjugate_matrix(A::CTMatrixTypes)
    R = base_ring(A)
    q2 = order(R)
    is_square(q2) || throw(ArgumentError("The Hermitian conjugate is only defined over quadratic field extensions."))

    q = Int(sqrt(q2, check = false))
    # return matrix(R, A .^ q)
    return A .^ q
end

# TODO: entropy function is incomplete
# """
#     entropy(x::Real)

# Return the entropy of the real number `x`.
# """
# function entropy(x::Real)
#     x != 0 || return 0
#     (0 < x <= 1 - 1 / q) || error("Number should be in the range [0, 1 - 1/order(field)].")
#     F = parent(x)
#     q = order(F)
#     return x * (log(q, q - 1) - log(q, x)) - (1 - x) * log(q, 1 - x)
# end

_Flint_matrix_element_to_Julia_int(x::fpMatrix, i::Int, j::Int) = ccall((:nmod_mat_get_entry,
    Oscar.Nemo.libflint), Int, (Ref{fpMatrix}, Int, Int), x, i - 1 , j - 1)

_Flint_matrix_to_Julia_int_matrix(A) = [ _Flint_matrix_element_to_Julia_int(A, i, j) for i in
    1:nrows(A), j in 1:ncols(A)]

_Flint_matrix_to_Julia_int_vector(A) = vec(_Flint_matrix_to_Julia_int_matrix(A))

function _Flint_matrix_to_Julia_bit_matrix(A::CTMatrixTypes)
    order(base_ring(A)) == 2 || throw(DomainError(A, "Only works for binary matrices"))
    BitMatrix(_Flint_matrix_to_Julia_int_matrix(A))
end

function _Flint_matrix_to_Julia_bool_matrix(A::CTMatrixTypes)
    order(base_ring(A)) == 2 || throw(DomainError(A, "Only works for binary matrices"))
    Matrix{Bool}(_Flint_matrix_to_Julia_int_matrix(A))
end

function _Flint_matrix_to_Julia_T_matrix(A::CTMatrixTypes, ::Type{T}) where T <: Number
    Matrix{T}(_Flint_matrix_to_Julia_int_matrix(A))
end

"""
Assumes the input is in rref form and returns the indexs of the columns that do not contain a pivot entry.
Note that rref form here requires pivot entries have been normalized to 1.
"""
function _rref_non_pivot_cols(A::CTMatrixTypes, type::Symbol = :nsp)
    type ∈ (:sp, :nsp) || throw(DomainError(type, "Parameter should be `:sp` (sparse) or `:nsp` (not sparse)."))

    if type == :sp
        return setdiff(collect(1:ncols(A)), [x.pos[1] for x in A])
    else
        nonpivots = Vector{Int}()
        i = 1
        j = 1
        nr, nc = size(A)
        while i <= nr && j <= nc
            if isone(A[i, j])
                i += 1
                j += 1
            else
                push!(nonpivots, j)
                j += 1
            end
        end
        append!(nonpivots, j:nc)
        return nonpivots
    end
end

function _quotient_space(big::T, small::T, alg::Symbol=:sys_eqs) where T <: CTMatrixTypes
    alg ∈ [:VS, :sys_eqs] || throw(ArgumentError("Unknown algorithm type"))

    if alg == :VS
        F = base_ring(big)
        V = vector_space(F, ncols(big))
        U, U_to_V = sub(V, [V(small[i, :]) for i in 1:nrows(small)])
        W, W_to_V = sub(V, [V(big[i, :]) for i in 1:nrows(big)])
        gens_of_U_in_W = Vector{typeof(gens(U)[1])}(undef, length(gens(U)))
        # gens_of_U_in_W = [preimage(W_to_V, U_to_V(g)) for g in gens(U)]
        Threads.@threads for i in 1:length(gens(U))
            gens_of_U_in_W[i] = preimage(W_to_V, U_to_V(gens(U)[i]))
        end
        U_in_W, _ = sub(W, gens_of_U_in_W)
        Q, W_to_Q = quo(W, U_in_W)
        iszero(dim(Q)) && (return zero_matrix(F, 1, ncols(big));)
        C2_mod_C1_basis = [W_to_V(x) for x in [preimage(W_to_Q, g) for g in gens(Q)]]
        F_basis = [[F(C2_mod_C1_basis[j][i]) for i in 1:AbstractAlgebra.dim(parent(C2_mod_C1_basis[1]))] for j in 1:length(C2_mod_C1_basis)]
        return matrix(F, length(F_basis), length(F_basis[1]), reduce(vcat, F_basis))
    else
        # solve the system x big = small
        # sol contains the way to write the rows of small in terms of the rows of big
        # if big is of the form (big = small ∪ (big / small)), then this will have zeros for the rows
        # corresponding to the basis of the quotient
        # in the general case, anything without a pivot in the row reduction is not required to make
        # the elements of small and therefore lie in the quotient space
        flag, sol = can_solve_with_solution(big, small, side=:left)
        !flag && error("Cannot solve system for quotient")
        _, rref_sol = rref(sol)
        if typeof(rref_sol) <: SMat{W, Vector{W}} where W <: CTFieldElem
            nonpivots = _rref_non_pivot_cols(rref_sol, :sp)
            return reduce(vcat, [big[r, :] for r in nonpivots])
        else
            return big[_rref_non_pivot_cols(rref_sol, :nsp), :]
        end
    end
end

# NOTE: This code works for sorted vectors with unique elements, but can be improved a bit in that case. It does not work otherwise, e.g.:
#   largestconsecrun([1,1,1,4]) == 4
# function largestconsecrun(arr::Vector{Int})
#     n = length(arr)
#     maxlen = 1
#     for i in 1:n
#         mn = arr[i]
#         mx = arr[i]
#         for j in (i + 1):n
#             mn = min(mn, arr[j])
#             mx = max(mx, arr[j])
#             if (mx - mn) == (j - i)
#                 maxlen = max(maxlen, mx - mn + 1)
#             end
#         end
#     end
#     return maxlen
# end

function _remove_empty(A::Union{CTMatrixTypes, Matrix{<: Number}, BitMatrix, Matrix{Bool}},
    type::Symbol)
    
    type ∈ (:rows, :cols) || throw(ArgumentError("Unknown type in _remove_empty"))
    
    del = Vector{Int}()
    if type == :rows
        for r in axes(A, 1)
            flag = true
            for c in axes(A, 2)
                if !iszero(A[r, c])
                    flag = false
                    break
                end
            end
            flag && append!(del, r)
        end
        return isempty(del) ? A : A[setdiff(axes(A, 1), del), :]
    elseif type == :cols
        for c in axes(A, 2)
            flag = true
            for r in axes(A, 1)
                if !iszero(A[r, c])
                    flag = false
                    break
                end
            end
            flag && append!(del, c)
        end
        return isempty(del) ? A : A[:, setdiff(axes(A, 2), del)]
    end
end

function _rref_no_col_swap(M::CTMatrixTypes, row_range::AbstractUnitRange{Int} = axes(M, 1),
    col_range::AbstractUnitRange{Int} = axes(M, 2))

    A = deepcopy(M)
    _rref_no_col_swap!(A, row_range, col_range)
    return A
end

function _rref_no_col_swap_binary(A::Union{BitMatrix, Matrix{Bool}, Matrix{<: Integer}},
    row_range::AbstractUnitRange{Int} = axes(A, 1), col_range::AbstractUnitRange{Int} = axes(A, 2))

    B = deepcopy(A)
    _rref_no_col_swap_binary!(B, row_range, col_range)
    return B
end

function _rref_no_col_swap!(A::CTMatrixTypes, row_range::AbstractUnitRange{Int} = axes(A, 1),
    col_range::AbstractUnitRange{Int} = axes(A, 2))

    # don't do anything to A if the range is empty
    isempty(row_range) && return nothing
    isempty(col_range) && return nothing

    i = first(row_range)
    j = first(col_range)
    nr = last(row_range)
    nc = last(col_range)
    if Int(order(base_ring(A))) != 2
        while i <= nr && j <= nc
            # find first pivot
            ind = 0
            for k in i:nr
                if !iszero(A[k, j])
                    ind = k
                    break
                end
            end

            if !iszero(ind)
                # normalize pivot
                if !isone(A[ind, j])
                    A[ind, :] *= inv(A[ind, j])
                end

                # swap to put the pivot in the next row
                if ind != i
                    A[i, :], A[ind, :] = A[ind, :], A[i, :]
                end

                # eliminate
                for k in row_range
                    if k != i
                        # do a manual loop here to reduce allocations
                        d = A[k, j]
                        @simd for l in axes(A, 2)
                            A[k, l] = (A[k, l] - d * A[i, l])
                        end
                    end
                end
                i += 1
            end
            j += 1
        end
    else
        while i <= nr && j <= nc
            # find first pivot
            ind = 0
            for k in i:nr
                if !iszero(A[k, j])
                    ind = k
                    break
                end
            end
    
            if !iszero(ind)
                # swap to put the pivot in the next row
                if ind != i
                    A[i, :], A[ind, :] = A[ind, :], A[i, :]
                end
    
                # eliminate
                for k in row_range
                    if k != i
                        if isone(A[k, j])
                            # do a manual loop here to reduce allocations
                            @simd for l in axes(A, 2)
                                A[k, l] += A[i, l]
                            end
                        end
                    end
                end
                i += 1
            end
            j += 1
        end
    end
    return nothing
end

function _rref_no_col_swap_binary!(A::Union{BitMatrix, Matrix{Bool}, Matrix{<: Integer}},
    row_range::AbstractUnitRange{Int} = axes(A, 1), col_range::AbstractUnitRange{Int} = axes(A, 2))

    isempty(row_range) && return nothing
    isempty(col_range) && return nothing
    i = first(row_range)
    j = first(col_range)
    nr = last(row_range)
    nc = last(col_range)
    while i <= nr && j <= nc
        # find first pivot
        ind = 0
        for k in i:nr
            if !iszero(A[k, j])
                ind = k
                break
            end
        end

        if !iszero(ind)
            # swap to put the pivot in the next row
            if ind != i
                A[i, :], A[ind, :] = A[ind, :], A[i, :]
            end

            # eliminate
            for k in row_range
                if k != i
                    if !iszero(A[k, j])
                        # do a manual loop here to reduce allocations
                        @simd for l in axes(A, 2)
                            A[k, l] ⊻= A[i, l]
                        end
                    end
                end
            end
            i += 1
        end
        j += 1
    end
    return nothing
end

function _rref_col_swap(M::CTMatrixTypes, row_range::AbstractUnitRange{Int} = axes(M, 1),
    col_range::AbstractUnitRange{Int} = axes(M, 2))

    A = deepcopy(M)
    rnk, P = _rref_col_swap!(A, row_range, col_range)
    return rnk, A, P
end

function _rref_col_swap!(A::CTMatrixTypes, row_range::AbstractUnitRange{Int} = axes(A, 1),
    col_range::AbstractUnitRange{Int} = axes(A, 2))

    # don't do anything to A if the range is empty, return rank 0 and missing permutation matrix
    isempty(row_range) && return 0, missing
    isempty(col_range) && return 0, missing

    nc_A = ncols(A)

    rnk = 0
    i = first(row_range)
    j = first(col_range)
    nr = last(row_range)
    nc = last(col_range)
    # permutation matrix required to return to rowspace if column swap done
    P = identity_matrix(base_ring(A), nc_A)
    if Int(order(base_ring(A))) != 2
        while i <= nr && j <= nc
            # find first pivot
            ind = 0
            for k in i:nr
                if !iszero(A[k, j])
                    ind = k
                    break
                end
            end

            # need to column swap
            if iszero(ind)
                for k in j + 1:nc
                    for l in i:nr
                        if !iszero(A[l, k])
                            ismissing(P) && (P = identity_matrix(base_ring(A), nc_A);)
                            swap_cols!(A, k, j)
                            swap_rows!(P, k, j)
                            ind = l
                            break
                        end
                    end
                end
            end

            # if true, the rest of the submatrix is zero
            if iszero(ind)
                return rnk, P
            else
                # normalize pivot
                if !isone(A[ind, j])
                    A[ind:ind, :] *= inv(A[ind, j])
                end

                # swap to put the pivot in the next row
                ind != i && swap_rows!(A, ind, i)

                # eliminate
                for k = first(row_range):nr
                    if k != i
                        # do a manual loop here to reduce allocations
                        d = A[k, j]
                        @simd for l = 1:nc_A
                            A[k, l] = (A[k, l] - d * A[i, l])
                        end
                    end
                end
            end
            i += 1
            j += 1
            rnk += 1
        end
    else
        while i <= nr && j <= nc
            # find first pivot
            ind = 0
            for k in i:nr
                if !iszero(A[k, j])
                    ind = k
                    break
                end
            end
    
            # need to column swap
            if iszero(ind)
                for k in j + 1:nc
                    for l in i:nr
                        if !iszero(A[l, k])
                            ismissing(P) && (P = identity_matrix(base_ring(A), nc_A);)
                            swap_cols!(A, k, j)
                            swap_rows!(P, k, j)
                            ind = l
                            break
                        end
                    end
                end
            end
    
            # if true, the rest of the submatrix is zero
            if iszero(ind)
                return rnk, P
            else    
                # swap to put the pivot in the next row
                ind != i && swap_rows!(A, ind, i)
    
                # eliminate
                for k = first(row_range):nr
                    if k != i
                        if isone(A[k, j])
                            # do a manual loop here to reduce allocations
                            @simd for l = 1:nc_A
                                A[k, l] += A[i, l]
                            end
                        end
                    end
                end
            end
            i += 1
            j += 1
            rnk += 1
        end
    end
    return rnk, P
end

function _rref_col_swap_perm(M::CTMatrixTypes, row_range::AbstractUnitRange{Int} = axes(M, 1), col_range::AbstractUnitRange{Int} = axes(M, 2))

    A = deepcopy(M)
    rnk, P = _rref_col_swap_perm!(A, row_range, col_range)
    return rnk, A, P
end

function _rref_col_swap_perm!(A::CTMatrixTypes, row_range::AbstractUnitRange{Int} = axes(A, 1),
    col_range::AbstractUnitRange{Int} = axes(A, 2))

    # don't do anything to A if the range is empty, return rank 0 and missing permutation matrix
    isempty(row_range) && return 0, missing
    isempty(col_range) && return 0, missing

    nc_A = ncols(A)

    rnk = 0
    i = first(row_range)
    j = first(col_range)
    nr = last(row_range)
    nc = last(col_range)
    sym_group = symmetric_group(nc)
    # permutation to return to rowspace if column swap done
    P = cperm(sym_group) 
    if Int(order(base_ring(A))) != 2
        while i <= nr && j <= nc
            # find first pivot
            ind = 0
            for k in i:nr
                if !iszero(A[k, j])
                    ind = k
                    break
                end
            end

            # need to column swap
            if iszero(ind)
                for k in j + 1:nc
                    for l in i:nr
                        if !iszero(A[l, k])
                            swap_cols!(A, k, j)
                            P = P * cperm(sym_group, [k,j])
                            ind = l
                            break
                        end
                    end
                end
            end

            # if true, the rest of the submatrix is zero
            if iszero(ind)
                return rnk, P
            else
                # normalize pivot
                if !isone(A[ind, j])
                    A[ind:ind, :] *= inv(A[ind, j])
                end

                # swap to put the pivot in the next row
                ind != i && swap_rows!(A, ind, i)

                # eliminate
                for k = first(row_range):nr
                    if k != i
                        # do a manual loop here to reduce allocations
                        d = A[k, j]
                        @simd for l = 1:nc_A
                            A[k, l] = (A[k, l] - d * A[i, l])
                        end
                    end
                end
            end
            i += 1
            j += 1
            rnk += 1
        end
    else
        while i <= nr && j <= nc
            # find first pivot
            ind = 0
            for k in i:nr
                if !iszero(A[k, j])
                    ind = k
                    break
                end
            end
    
            # need to column swap
            if iszero(ind)
                for k in j + 1:nc
                    for l in i:nr
                        if !iszero(A[l, k])
                            swap_cols!(A, k, j)
                            P = P * cperm(sym_group, [k,j])
                            ind = l
                            break
                        end
                    end
                end
            end
    
            # if true, the rest of the submatrix is zero
            if iszero(ind)
                return rnk, P
            else    
                # swap to put the pivot in the next row
                ind != i && swap_rows!(A, ind, i)
    
                # eliminate
                for k = first(row_range):nr
                    if k != i
                        if isone(A[k, j])
                            # do a manual loop here to reduce allocations
                            @simd for l = 1:nc_A
                                A[k, l] += A[i, l]
                            end
                        end
                    end
                end
            end
            i += 1
            j += 1
            rnk += 1
        end
    end
    return rnk, inv(P)
end

function _rref_symp_col_swap(A::CTMatrixTypes, row_range::AbstractUnitRange{Int} = axes(A, 1),
    col_range::AbstractUnitRange{Int} = axes(A, 2))
    B = deepcopy(A)
    _rref_symp_col_swap!(B, row_range, col_range)
    return B
end

function _rref_symp_col_swap!(A::CTMatrixTypes, row_range::AbstractUnitRange{Int} = axes(A, 1),
    col_range::AbstractUnitRange{Int} = axes(A, 2))

    # don't do anything to A if the range is empty, return rank 0 and missing permutation matrix
    isempty(row_range) && return 0, missing
    isempty(col_range) && return 0, missing

    # permutation matrix required to return to rowspace if column swap done
    P = missing
    nc_A = ncols(A)

    rnk = 0
    i = first(row_range)
    j = first(col_range)
    nr = last(row_range)
    nc = last(col_range)
    if Int(order(base_ring(A))) != 2
        while i <= nr && j <= nc
            # find first pivot
            ind = 0
            for k in i:nr
                if !iszero(A[k, j])
                    ind = k
                    break
                end
            end

            # need to column swap
            if iszero(ind)
                for k in j + 1:nc
                    for l in i:nr
                        if !iszero(A[l, k])
                            ismissing(P) && (P = identity_matrix(base_ring(A), nc_A);)
                            k_symp = mod1(k + div(nc_A, 2), nc_A)
                            j_symp = mod1(j + div(nc_A, 2), nc_A)
                            swap_cols!(A, k, j)
                            swap_rows!(P, k, j)
                            swap_cols!(A, k_symp, j_symp)
                            swap_rows!(P, k_symp, j_symp)
                            ind = l
                            break
                        end
                    end
                end
            end

            # if true, the rest of the submatrix is zero
            if iszero(ind)
                return rnk, P
            else
                # normalize pivot
                if !isone(A[ind, j])
                    A[ind, :] *= inv(A[ind, j])
                end

                # swap to put the pivot in the next row
                ind != i && swap_rows!(A, ind, i)

                # eliminate
                for k = first(row_range):nr
                    if k != i
                        # do a manual loop here to reduce allocations
                        d = A[k, j]
                        @simd for l = 1:nc_A
                            A[k, l] = (A[k, l] - d * A[i, l])
                        end
                    end
                end
            end
            i += 1
            j += 1
            rnk += 1
        end
    else
        while i <= nr && j <= nc
            # find first pivot
            ind = 0
            for k in i:nr
                if !iszero(A[k, j])
                    ind = k
                    break
                end
            end

            # need to column swap
            if iszero(ind)
                for k in j + 1:nc
                    for l in i:nr
                        if !iszero(A[l, k])
                            ismissing(P) && (P = identity_matrix(base_ring(A), nc_A);)
                            k_symp = mod1(k + div(nc_A, 2), nc_A)
                            j_symp = mod1(j + div(nc_A, 2), nc_A)
                            swap_cols!(A, k, j)
                            swap_rows!(P, k, j)
                            swap_cols!(A, k_symp, j_symp)
                            swap_rows!(P, k_symp, j_symp)
                            ind = l
                            break
                        end
                    end
                end
            end

            # if true, the rest of the submatrix is zero
            if iszero(ind)
                return rnk, P
            else
                # swap to put the pivot in the next row
                ind != i && swap_rows!(A, ind, i)

                # eliminate
                for k = first(row_range):nr
                    if k != i
                        if isone(A[k, j])
                            # do a manual loop here to reduce allocations
                            @simd for l = 1:nc_A
                                A[k, l] += A[i, l]
                            end
                        end
                    end
                end
            end
            i += 1
            j += 1
            rnk += 1
        end
    end
    return rnk, P
end

function _col_permutation!(X::Matrix{T}, A::Matrix{T}, p::AbstractVector{Int}) where T
    length(p) == size(A, 2) || throw(ArgumentError("`p` should have length `size(A, 2)`."))
    size(X) == size(A) || throw(ArgumentError("`X` and `A` should have the same shape."))
    for j in axes(X, 2)
        for i in axes(X, 1)
            X[i, j] = A[i, p[j]]
        end
    end
    return nothing
end

function _col_permutation_symp!(X::Matrix{T}, A::Matrix{T}, p::AbstractVector{Int}) where T
    n = length(p)
    # 2n == size(A, 2) || throw(ArgumentError("`p` should have length `size(A, 2)/2`."))
    # size(X) == size(A) || throw(ArgumentError("`X` and `A` should have the same shape."))
    for j in 1:n
        for i in axes(X, 1)
            X[i, j] = A[i, p[j]]
        end
    end
    for j in 1:n
        for i in axes(X, 1)
            X[i, j + n] = A[i, p[j] + n]
        end
    end
    return nothing
end

function digits_to_int(x::Vector{Int}, base::Int=2)
    res = 0
    for digit in x
        res = digit + base * res
    end
    return res
end

function _CT_adjoint(A::MatElem{T}) where T <: ResElem
    R = parent(A[1, 1])
    S = base_ring(A[1, 1])
    f = modulus(R)
    l = degree(f)
    Fz = base_ring(S)(0)

    A_tr = transpose(A)
    nr, nc = size(A_tr)
    for c in 1:nc
        for r in 1:nr
            h_coeffs = collect(coefficients(Nemo.lift(A_tr[r, c])))
            for _ in 1:l - length(h_coeffs)
                push!(h_coeffs, Fz)
            end
            h_coeffs[2:end] = reverse(h_coeffs[2:end])
            A_tr[r, c] = R(S(h_coeffs))
        end
    end
    return A_tr
end

function _CT_adjoint(A::MatElem{T}) where T <: CTGroupAlgebra
    FG = parent(A[1, 1])
    f = hom(FG, FG, basis_matrix([inv(g) for g in basis(FG)]))
    # for each element of A, map the sum \sum_i c_i g_i to \sum_i c_i (g_i)^-1
    return map(f, transpose(A))
end

"""
    residue_polynomial_to_circulant_matrix(f::ResElem)

    Return the circulant matrix whose first row or column is the coefficients of `f` if `type` is `:row` or `:col`, respectively.
"""
function residue_polynomial_to_circulant_matrix(f::ResElem, type::Symbol=:col)
    type ∈ (:col, :row) || throw(ArgumentError("Unknown type"))

    R = parent(f)
    S = base_ring(R)
    F = base_ring(S)
    g = modulus(R)
    l = degree(g)
    g == gen(S)^l - 1 || throw(ArgumentError("Residue ring not of the form x^l - 1."))
    # gcd(l, Int(characteristic(F))) == 1 || throw(ArgumentError("Residue ring over F_q[x] must have modulus x^l - 1 with gcd(l, q) = 1."))

    A = zero_matrix(F, l, l)
    if type == :col
        F_coeffs = zero_matrix(F, l, 1)
        temp = collect(coefficients(Nemo.lift(f)))
        F_coeffs[1:length(temp), 1] = temp
        A[:, 1] = F_coeffs
        for c in 2:l
            # A[:, c] = circshift(F_coeffs, c - 1)
            A[1:c - 1, c] = F_coeffs[l - (c - 1) + 1:l, 1]
            A[c:end, c] = F_coeffs[1:l - (c - 1), 1]
        end
    elseif type == :row
        F_coeffs = zero_matrix(F, 1, l)
        temp = collect(coefficients(Nemo.lift(f)))
        F_coeffs[1, 1:length(temp)] = temp
        A[1, :] = F_coeffs
        for c in 2:l
            A[c, 1:c - 1] = F_coeffs[1, l - (c - 1) + 1:l]
            A[c, c:end] = F_coeffs[1, 1:l - (c - 1)]
        end
    end
    return A
end

"""
    group_algebra_element_to_circulant_matrix(x::CTGroupAlgebra; type::Symbol=:col)

Return the circulant matrix whose first row or column is the coefficients of `x` if `type` is `:row` or `:col`, respectively.
"""
function group_algebra_element_to_circulant_matrix(x::CTGroupAlgebra, type::Symbol=:col)
    type ∈ (:col, :row) || throw(ArgumentError("Unknown type"))
    
    F = base_ring(parent(x))
    F_coeffs = coefficients(x)
    l = length(F_coeffs)
    A = zero_matrix(F, l, l)
    if type == :col
        F_coeffs = matrix(F, l, 1, F_coeffs)
        A[:, 1] = F_coeffs
        for c in 2:l
            A[1:c - 1, c] = F_coeffs[l - (c - 1) + 1:l, 1]
            A[c:end, c] = F_coeffs[1:l - (c - 1), 1]
        end
    elseif type == :row
        F_coeffs = matrix(F, 1, l, F_coeffs)
        A[1, :] = F_coeffs
        for c in 2:l
            A[c, 1:c - 1] = F_coeffs[1, l - (c - 1) + 1:l]
            A[c, c:end] = F_coeffs[1, 1:l - (c - 1)]
        end
    end
    return A
end

"""
    lift(A::MatElem{T}, type::Symbol=:col) where T <: ResElem

Return the matrix whose residue polynomial elements are converted to circulant matrices
over the base field.
"""
function lift(A::MatElem{T}, type::Symbol=:col) where T <: ResElem
    type ∈ (:col, :row) || throw(ArgumentError("Unknown type"))

    R = parent(A[1, 1])
    S = base_ring(R)
    F = base_ring(S)
    g = modulus(R)
    l = degree(g)
    g == gen(S)^l - 1 || throw(ArgumentError("Residue ring not of the form x^l - 1."))
    # gcd(l, Int(characteristic(F))) == 1 || throw(ArgumentError("Residue ring over F_q[x] must have modulus x^l - 1 with gcd(l, q) = 1."))

    nr, nc = size(A)
    A_lift = zero_matrix(F, nr * l, nc * l)
    for c in axes(A, 2)
        for r in axes(A, 1)
            if !iszero(A[r, c])
                A_lift[(r - 1) * l + 1:r * l, (c - 1) * l + 1:c * l] =
                    residue_polynomial_to_circulant_matrix(A[r, c], type)
            end
        end
    end
    return A_lift
end

"""
    lift(A::MatElem{T}, type::Symbol=:col) where T <: CTGroupAlgebra

Return the matrix whose group algebra elements are converted to circulant matrices
over the base field.
"""
function lift(A::MatElem{T}, type::Symbol=:col) where T <: CTGroupAlgebra
    type ∈ (:col, :row) || throw(ArgumentError("Unknown type"))

    F = base_ring(parent(A[1, 1]))
    l = length(coefficients(A[1, 1]))
    nr, nc = size(A)
    A_lift = zero_matrix(F, nr * l, nc * l)
    for c in axes(A, 2)
        for r in axes(A, 1)
            if !iszero(A[r, c])
                A_lift[(r - 1) * l + 1:r * l, (c - 1) * l + 1:c * l] =
                    group_algebra_element_to_circulant_matrix(A[r, c], type)
            end
        end
    end
    return A_lift
end

# Creates a matrix with copies of `M` at every nonzero entry of `locations`.
function _concat(locations::Union{CTMatrixTypes, Matrix}, M::CTMatrixTypes)
    nr_M, nc_M = size(M)
    nr_L, ncL = size(locations)
    output = zero_matrix(base_ring(M), nr_M * nr_L, nc_M * ncL)
    for j_outer in 1:ncL
        for i_outer in 1:nr_L
            if !iszero(locations[i_outer, j_outer])
                for j in 1:nc_M
                    for i in 1:nr_M
                        row = i + nr_M * (i_outer - 1)
                        col = j + nc_M * (j_outer - 1)
                        output[row, col] = M[i, j]
                    end
                end
            end
        end
    end
    return output
end

""""
    row_supports(M::CTMatrixTypes)

Returns a vector where the ith entry lists the indices of the nonzero
entries of `M[i, :]`
"""
function row_supports(M::Union{CTMatrixTypes,
    MatElem{EuclideanRingResidueRingElem{fpPolyRingElem}}})

    output = [Int[] for _ in axes(M, 1)]
    for j in axes(M, 2)
        for i in axes(M, 1)
            iszero(M[i, j]) || push!(output[i], j)
        end
    end
    return output
end

""""
    row_supports_symplectic(M::CTMatrixTypes)

Returns a vector where the ith entry is a 2-tuple of lists with the
indices of the nonzero X and Z entries of `M[i, :]`
"""
function row_supports_symplectic(M::CTMatrixTypes)
    iseven(ncols(M)) || throw(ArgumentError("Matrix should have an even number of cols"))
    n = div(ncols(M), 2)
    X = row_supports(view(M, :, 1:n))
    Z = row_supports(view(M, :, 1 + n:2n))
    collect(zip(X, Z))
end

function _node_adjacencies(H::CTMatrixTypes)
    check_adj_list = [Int[] for _ in 1:nrows(H)]
    var_adj_list = [Int[] for _ in 1:ncols(H)]
    for r in 1:nrows(H)
        for c in 1:ncols(H)
            if !iszero(H[r, c])
                push!(check_adj_list[r], c)
                push!(var_adj_list[c], r)
            end
        end
    end
    return check_adj_list, var_adj_list
end

"""
    strongly_lower_triangular_reduction(A::CTMatrixTypes)

Return a strongly lower triangular basis for the kernel of `A` and
a unit vector basis for the complement of the image of `transpose(A)`.

* Note
- This implements Algorithm 1 from https://doi.org/10.48550/arXiv.2204.10812
"""
function strongly_lower_triangular_reduction(A::CTMatrixTypes)
    B = deepcopy(A)
    F = base_ring(B)
    nr, nc = size(B)
    id_mat = identity_matrix(F, nc)
    κ = deepcopy(id_mat)
    π = collect(1:nc)
    for j in 1:nc
        i = 1
        while i < nr && !isone(B[i, j])
            i += 1
        end

        if isone(B[i, j])
            # more natural and probably faster to push pivots to a list
            π = setdiff(π, [j])
            for l in j+1:nc
                if isone(B[i, l])
                    B[:, l] += B[:, j]
                    κ[:, l] += κ[:, j]
                end
            end
        end
    end

    ker = zero_matrix(F, nc, length(π))
    im = deepcopy(ker)
    i = 1
    for j in π
        ker[:, i] = κ[:, j]
        im[:, i] = id_mat[:, j]
        i += 1
    end
    return ker, im
end

"""
    load_alist(file::String)

Return a `Matrix{Int}` object from the matrix stored in the alist file format in `file`.
"""
function load_alist(file::String)
    contents = split.(readlines(file))
    length(contents[1]) == 2 || throw(ArgumentError("Not a valid alist file. First line wrong."))
    length(contents[2]) == 2 || throw(ArgumentError("Not a valid alist file. Second line wrong."))
    M = parse(Int, contents[1][1])
    N = parse(Int, contents[1][2])
    length(contents) == M + N + 4 || throw(ArgumentError("Not a valid alist file. Wrong number of lines."))
    row_wts = parse.(Int, contents[3])
    col_wts = parse.(Int, contents[4])
    length(row_wts) == M || throw(ArgumentError("Not a valid alist file. Wrong number of row weights."))
    length(col_wts) == N || throw(ArgumentError("Not a valid alist file. Wrong number of column weights."))
    maximum(row_wts) == parse(Int, contents[2][1]) || throw(ArgumentError("Not a valid alist file. Row weight mismatch."))
    maximum(col_wts) == parse(Int, contents[2][2]) || throw(ArgumentError("Not a valid alist file. Col weight mismatch."))
    all(allunique, contents[5:end]) || throw(ArgumentError("Not a valid alist file. Indices are repeated."))
    all(length(contents[j]) == row_wts[i] for (i, j) in enumerate(5:M + 4)) || throw(ArgumentError("Not a valid alist file. Row weights don't match."))
    all(length(contents[j]) == col_wts[i] for (i, j) in enumerate(M + 5:M + N + 4)) || throw(ArgumentError("Not a valid alist file. Column weights don't match."))

    mat = zeros(Int, M, N)
    for (col, i) in enumerate(M + 5:M + N + 4)
        rows = parse.(Int, contents[i])
        for row in rows
            mat[row, col] = 1
        end
    end

    for (row, i) in enumerate(5:M + 4)
        cols = parse.(Int, contents[i])
        for col in 1:N
            mat[row, col] == (col in cols) || throw(ArgumentError("Not a valid alist file. The two matrix representations don't match."))
        end
    end
    return mat
end

# TODO polish this and add/export
# function save_mtx(h, path, info)
#     m, n = size(h)
#     num_nonzero = sum(count(!=(0), h, dims=2))
#     open(path, "w") do io
#         write(io, "%%MatrixMarket matrix coordinate integer general\n")
#         write(io, "%%$(info)")
#         write(io, "\n$m $n $(num_nonzero)\n")
#         for i in 1:m
#             for j in 1:n
#                 if h[i, j] == 1
#                     write(io, "$i $j 1\n")
#                 end
#             end
#         end
#     end
# end

#############################
  # Quantum Helper Functions
#############################

"""
    is_triorthogonal(G::CTMatrixTypes, verbose::Bool=false)
    is_triorthogonal(G::Matrix{Int}, verbose::Bool=false)

Return `true` if the binary matrix `G` is triorthogonal.

# Notes
* If the optional parameter `verbos` is set to `true`, the first pair or triple of
  non-orthogonal rows will be identified on the console.
"""
function is_triorthogonal(G::CTMatrixTypes, verbose::Bool=false)
    Int(order(base_ring(G))) == 2 || throw(ArgumentError("Triothogonality is only defined over 𝔽₂."))
    nr, nc = size(G)
    for r1 in 1:nr
        for r2 in (r1 + 1):nr
            if !iszero(sum(G[r1, i] * G[r2, i] for i in 1:nc))
                verbose && println("Rows $r1 and $r2 are not orthogonal.")
                return false
            end
            for r3 in (r2 + 1):nr
                if !iszero(sum(G[r1, i] * G[r2, i] * G[r3, i] for i in 1:nc))
                    verbose && println("Rows $r1, $r2, and $r3 are not orthogonal.")
                    return false
                end
            end
        end
    end
    return true
end

function is_triorthogonal(G::Matrix{Int}, verbose::Bool=false)
    nr, nc = size(G)
    for r1 in 1:nr
        for r2 in (r1 + 1):nr
            if !iszero(sum(G[r1, i] * G[r2, i] for i in 1:nc) % 2)
                verbose && println("Rows $r1 and $r2 are not orthogonal.")
                return false
            end
            for r3 in (r2 + 1):nr
                if !iszero(sum(G[r1, i] * G[r2, i] * G[r3, i] for i in 1:nc) % 2)
                    verbose && println("Rows $r1, $r2, and $r3 are not orthogonal.")
                    return false
                end
            end
        end
    end
    return true
end

function print_string_array(A::Vector{String}, without_Is=false)
    for a in A
        if !withoutIs
            println(a)
        else
            for i in a
                if i == 'I'
                    print(' ')
                else
                    print(i)
                end
            end
            print('\n')
        end
    end
end
# BUG: do these set functions exist anymore?
print_char_array(A::Vector{Vector{Char}}, without_Is=false) = print_string_array(set_char_to_string_array(A), without_Is)
printsymplecticarray(A::Vector{Vector{T}}, without_Is=false) where T <: Int = print_string_array(set_symplectic_to_string_array(A), without_Is)

"""
    pseudoinverse(M::CTMatrixTypes)

Return the pseudoinverse of a stabilizer matrix `M` over a quadratic extension.

# Notes
* This is not the Penrose-Moore pseudoinverse.
"""
function pseudoinverse(M::CTMatrixTypes)
    # let this fail elsewhere if not actually over a quadratic extension
    if degree(base_ring(M)) != 1
        # TODO: quadratic_to_symplectic is no longer defined, this will need changed
        M = transpose(quadratic_to_symplectic(M))
    else
        M = transpose(M)
    end

    nr, nc = size(M)
    _, E = rref(hcat(M, identity_matrix(base_ring(M), nr)))
    E = E[:, (nc + 1):end]
    p_inv = E[1:nc, :]
    dual = E[nc + 1:nr, :]

    # verify
    _, M_rref = rref(M)
    E * M == M_rref || error("Pseudoinverse calculation failed (transformation incorrect).")
    M_rref[1:nc, 1:nc] == identity_matrix(base_ring(M), nc) || error("Pseudoinverse calculation failed (failed to get I).")
    iszero(M_rref[nc + 1:nr, :]) || error("Pseudoinverse calculation failed (failed to get zero).")
    p_inv * M == identity_matrix(base_ring(M), nc) || error("Pseudoinverse calculation failed (eq 1).")
    transpose(M) * transpose(p_inv) == identity_matrix(base_ring(M), nc) || error("Pseudoinverse calculation failed (eq 2).")
    iszero(transpose(M) * transpose(dual)) || error("Failed to correctly compute dual (rhs).")
    iszero(dual * M) || error("Failed to correctly compute dual (lhs).")
    return p_inv
end

# """
#     quadratic_to_symplectic(M::CTMatrixTypes)

# Return the matrix `M` converted from the quadratic to the symplectic form.
# """
# function quadratic_to_symplectic(M::CTMatrixTypes)
#     E = base_ring(M)
#     iseven(degree(E)) || error("The base ring of the given matrix is not a quadratic extension.")
#     F = GF(Int(characteristic(E)), div(degree(E), 2), :ω)
#     nr = nrows(M)
#     nc = ncols(M)
#     Msym = zero_matrix(F, nr, 2 * nc)
#     for c in 1:nc
#         for r in 1:nr
#             # TODO: benchmark this without the branching
#             if !iszero(M[r, c])
#                 Msym[r, c] = F(coeff(M[r, c], 0))
#                 Msym[r, c + nc] = F(coeff(M[r, c], 1))
#             end
#         end
#     end
#     return Msym
# end

# """
#     symplectictoquadratic(M::CTMatrixTypes)

# Return the matrix `M` converted from the symplectic to the quadratic form.
# """
# function symplectictoquadratic(M::CTMatrixTypes)
#     iseven(ncols(M)) || error("Input to symplectictoquadratic is not of even length.")
#     nr = nrows(M)
#     nc = div(ncols(M), 2)
#     F = base_ring(M)
#     E = GF(Int(characteristic(F)), 2 * degree(F), :ω)
#     ω = gen(E)
#     ϕ = embed(F, E)
#     Mquad = zero_matrix(E, nr, nc)
#     for c in 1:nc
#         for r in 1:nr
#             Mquad[r, c] = ϕ(M[r, c]) + ϕ(M[r, c + nc]) * ω
#         end
#     end
#     return Mquad
# end

function _Pauli_string_to_symplectic(str::T) where T <: Union{String, Vector{Char}}
    n = length(str)
    # F = GF(2, 1, :ω)
    F = Oscar.Nemo.Native.GF(2)
    sym = zero_matrix(F, 1, 2 * n)
    for (i, c) in enumerate(str)
        if c == 'X'
            sym[1, i] = F(1)
        elseif c == 'Z'
            sym[1, i + n] = F(1)
        elseif c == 'Y'
            sym[1, i] = 1
            sym[1, i + n] = F(1)
        elseif c != 'I'
            error("Encountered non-{I, X, Y, Z} character in Pauli string. This function is only defined for binary strings.")
        end
    end
    return sym
end
_Pauli_string_to_symplectic(A::Vector{T}) where T <: Union{String, Vector{Char}} = reduce(vcat, [_Pauli_string_to_symplectic(s) for s in A])
# need symplectictoPaulistring

# charvec::Union{Vector{zzModRingElem}, Missing}=missing)
function _process_strings(SPauli::Vector{T}) where T <: Union{String, Vector{Char}}
    # Paulisigns = Vector{Int}()
    S_tr_Pauli_stripped = Vector{String}()
    for (i, s) in enumerate(SPauli)
        if s[1] ∈ ['I', 'X', 'Y', 'Z']
            # append!(Paulisigns, 1)
            push!(S_tr_Pauli_stripped, s)
        elseif s[1] == '+'
            # append!(Paulisigns, 1)
            push!(S_tr_Pauli_stripped, s[2:end])
        elseif s[1] == '-'
            # append!(Paulisigns, -1)
            push!(S_tr_Pauli_stripped, s[2:end])
        else
            error("The first element of Pauli string $i is neither a Pauli character or +/-: $s.")
        end
    end

    n = length(S_tr_Pauli_stripped[1])
    for s in S_tr_Pauli_stripped
        for i in s
            i ∈ ['I', 'X', 'Y', 'Z'] || error("Element of provided Pauli string is not a Pauli character: $s.")
        end
        length(s) == n || error("Not all Pauli strings are the same length.")
    end

    # if !ismissing(charvec)
    #     2 * n == length(charvec) || error("The characteristic value is of incorrect length.")
    #     R = ResidueRing(Nemo.ZZ, 4)
    #     for s in charvec
    #         modulus(s) == modulus(R) || error("Phases are not in the correct ring.")
    #     end
    # else
    #     R = ResidueRing(Nemo.ZZ, 4)
    #     charvec = [R(0) for _ in 1:2 * n]
    # end
    return S_tr_Pauli_stripped#, charvec
end

#############################
        # Finite Fields
#############################

"""
    tr(x::fqPolyRepFieldElem, K::fqPolyRepField, verify::Bool=false)

Return the relative trace of `x` from its base field to the field `K`.

# Notes
* If the optional parameter `verify` is set to `true`, the two fields are checked
  for compatibility.
"""
function tr(x::CTFieldElem, K::CTFieldTypes; verify::Bool = false)
    L = parent(x)
    q = order(K)
    if verify
        # # shouldn't need Int casting here but just in case...
        # Int(characteristic(L)) == Int(characteristic(K)) || error("The given field is not a subfield of the base ring of the element.")
        # degree(L) % degree(K) == 0 || error("The given field is not a subfield of the base ring of the element.")
        flag, m = is_extension(L, K)
        flag || throw(ArgumentError("The given field is not a subfield of the base ring of the matrix."))
    end
    n = div(degree(L), degree(K))
    return sum([x^(q^i) for i in 0:(n - 1)])
end

# function _expandelement(x::CTFieldElem, K::CTFieldTypes, basis::Vector{<:CTFieldElem}, verify::Bool=false)
#     return [tr(x * i) for i in basis] #, K, verify
# end

# function _expandrow(row::CTMatrixTypes, K::CTFieldTypes, basis::Vector{<:CTFieldElem}, verify::Bool=false)
#     # new_row = _expandelement(row[1], K, basis, verify)
#     new_row = [tr(row[1] * β) for β in basis]
#     for i in 2:ncols(row)
#         # new_row = vcat(new_row, _expandelement(row[i], K, basis, verify))
#         new_row = vcat(new_row, [tr(row[i] * β) for β in basis])
#     end
#     return matrix(K, 1, length(new_row), new_row)
# end

function _expansion_dict(L::CTFieldTypes, K::CTFieldTypes, λ::Vector{<:CTFieldElem})
    m = div(degree(L), degree(K))
    L_elms = collect(L)
    D = Dict{FqFieldElem, FqMatrix}()
    for x in L_elms
        D[x] = matrix(L, 1, m, [CodingTheory.tr(x * λi, K) for λi in λ])
    end
    return D
end

# BUG this is building the expanded matrix in the wrong ring, added change_base_ring below
function _expand_matrix(M::CTMatrixTypes, D::Dict{FqFieldElem, FqMatrix}, m::Int)
    m > 0 || throw(DomainError("Expansion factor must be positive"))

    M_exp = zero_matrix(base_ring(M), nrows(M), ncols(M) * m)
    for r in 1:nrows(M)
        for c in 1:ncols(M)
            M_exp[r, (c - 1) * m + 1:c * m] = D[M[r, c]]
        end
    end
    return M_exp
end

"""
    expand_matrix(M::CTMatrixTypes, K::fqPolyRepField, β::Vector{fqPolyRepFieldElem})

Return the matrix constructed by expanding the elements of `M` to the subfield
`K` using the basis `β` for the base ring of `M` over `K`.
"""
function expand_matrix(M::CTMatrixTypes, K::CTFieldTypes, β::Vector{<:CTFieldElem})
    L = base_ring(M)
    L == K && return M
    flag, m = is_extension(L, K)
    flag || throw(ArgumentError("The given field is not a subfield of the base ring of the matrix."))
    m == length(β) || throw(ArgumentError("Basis does not have length degree of the extension."))
    flag, λ = _is_basis(L, β, Int(order(K)))
    flag || throw(ArgumentError("The provided vector is not a basis for the extension."))

    # λ = dual_basis(L, K, β)
    D = _expansion_dict(L, K, λ)
    return change_base_ring(K, _expand_matrix(M, D, m))
end

"""
    quadratic_residues(q::Int, n::Int)

Return the sets of quadratic resides and quadratic non-residues of `q` and `n`.
"""
function quadratic_residues(q::Int, n::Int)
    isodd(n) && is_prime(n) || throw(ArgumentError("n must be an odd prime in quadratic residues"))
    q^div(n - 1, 2) % n == 1 || throw(ArgumentError("q^(n - 1)/2 ≅ 1 mod n in quadratic residues"))

    # F = GF(n, 1, :α)
    # elms = collect(F)
    # # skip 0
    # q_res = unique!([i^2 for i in elms[2:end]])
    # not_q_res = setdiff!(elms[2:end], q_res)
    # return q_res, not_q_res

    # don't want this returning in the field
    q_res = sort!(unique!([i^2 % n for i in 1:n - 1]))
    not_q_res = setdiff(1:(n - 1), q_res)
    return q_res, not_q_res
end

function _is_basis(E::CTFieldTypes, basis::Vector{<:CTFieldElem}, q::Int)
    m = length(basis)
    B = zero_matrix(E, m, m)
    for r in 1:m
        for c in 1:m
            B[r, c] = basis[r]^(q^(c - 1))
        end
    end
    iszero(det(B)) && return false, missing
    
    try
        B_inv = inv(B)
        λ = [B_inv[1, i] for i in 1:m]
        return true, λ
    catch
        return false, missing
    end
end

"""
    is_extension(E::CTFieldTypes, F::CTFieldTypes)

Return `true` if `E/F` is a valid field extension and the degree of the extension; otherwise return
`false, -1`.
"""
function is_extension(E::CTFieldTypes, F::CTFieldTypes)
    p = Int(characteristic(E))
    Int(characteristic(F)) == p || return false, -1
    deg_E = degree(E)
    deg_F = degree(F)
    deg_E % deg_F == 0 || return false, -1
    return true, div(deg_E, deg_F)
    # the below allows you to embed GF(2) into GF(5) without error but is not an extension
    # try
    #     embed(F, E)
    #     return true, div(degree(E), degree(F))
    # catch
    #     return false, missing
    # end
end

"""
    is_subfield(F::CTFieldTypes, E::CTFieldTypes)

Return `true` if `E/F` is a valid field extension and the degree of the extension; otherwise return
`false, -1`.
"""
is_subfield(F::CTFieldTypes, E::CTFieldTypes) = is_extension(E, F)

"""
    is_basis(E::fqPolyRepField, F::fqPolyRepField, basis::Vector{fqPolyRepFieldElem})

Return `true` and the dual (complementary) basis if `basis` is a basis for `E/F`,
otherwise return `false, missing`.
"""
function is_basis(E::CTFieldTypes, F::CTFieldTypes, basis::Vector{<:CTFieldElem})
    flag, m = is_extension(E, F)
    flag || throw(ArgumentError("Second field is not a subfield of the first."))
    length(basis) == m || throw(ArgumentError("Basis does not have length degree of the extension."))
    for i in 1:m
        parent(basis[i]) == E || throw(ArgumentError("The basis must be elements of the extension field."))
    end

    return _is_basis(E, basis, Int(order(F)))
end

"""
    primitive_basis(E::fqPolyRepField, F::fqPolyRepField)

Return a primitive basis for `E/F` and its dual (complementary) basis.
"""
function primitive_basis(E::CTFieldTypes, F::CTFieldTypes)
    flag, m = is_extension(E, F)
    flag || throw(ArgumentError("Second field is not a subfield of the first."))
    α = gen(E)
    basis = [α^i for i in 0:m - 1]
    flag, λ = _is_basis(E, basis, Int(order(F)))
    return basis, λ
end
# these are slightly different
# polynomialbasis(E::fqPolyRepField, F::fqPolyRepField) = primitive_basis(E, F)
# monomialbasis(E::fqPolyRepField, F::fqPolyRepField) = primitive_basis(E, F)

"""
    normal_basis(E::fqPolyRepField, F::fqPolyRepField)

Return a normal basis for `E/F` and its dual (complementary) basis.
"""
# "Normal Bases over Finite Fields" by Shuhong Gao has algorithms for this but they are
# complicated for the field sizes intended in this work
function normal_basis(E::CTFieldTypes, F::CTFieldTypes)
    flag, m = is_extension(E, F)
    flag || throw(ArgumentError("Second field is not a subfield of the first."))

    q = Int(order(F))
    elms = collect(E)
    for e in elms
        basis = [e^(q^i) for i in 0:m - 1]
        flag, dual_basis = _is_basis(E, basis, q)
        flag && return basis, dual_basis
    end
    error("Somehow failed to final a normal element for the extension.")
end

"""
    dual_basis(E::fqPolyRepField, F::fqPolyRepField, basis::Vector{fqPolyRepFieldElem})
    complementary_basis(E::fqPolyRepField, F::fqPolyRepField, basis::Vector{fqPolyRepFieldElem})

Return the dual (complentary) basis of `basis` for the extension `E/F`.
"""
function dual_basis(E::CTFieldTypes, F::CTFieldTypes, basis::Vector{<:CTFieldElem})
    flag, λ = is_basis(E, F, basis)
    flag || throw(ArgumentError("The provided vector is not a basis for the extension."))
    return λ
end
complementary_basis(E::CTFieldTypes, F::CTFieldTypes, basis::Vector{<:CTFieldElem}) = dual_basis(E, F, basis)

"""
    verify_dual_basis(E::fqPolyRepField, F::fqPolyRepField, basis::Vector{fqPolyRepFieldElem}, dual_basis::Vector{fqPolyRepFieldElem})
    verify_complementary_basis(E::fqPolyRepField, F::fqPolyRepField, basis::Vector{fqPolyRepFieldElem}, dual_basis::Vector{fqPolyRepFieldElem})

Return `true` if `basis` is the dual of `dual_basis` for `E/F`, otherwise return `false`.
"""
function verify_dual_basis(E::CTFieldTypes, F::CTFieldTypes, basis::Vector{<:CTFieldElem}, dual_basis::Vector{<:CTFieldElem})
    flag, m = is_extension(E, F)
    flag || throw(ArgumentError("Second field is not a subfield of the first."))

    m = length(basis)
    length(dual_basis) == m || throw(ArgumentError("The basis and dual basis must have the same length."))
    E = parent(basis[1])
    for i in 1:m
        parent(basis[i]) == E || throw(ArgumentError("Elements must be over the same field."))
        parent(dual_basis[i]) == E || throw(ArgumentError("Elements must be over the same field."))
    end

    q = Int(order(F))
    B = zero_matrix(E, m, m)
    for r in 1:m
        for c in 1:m
            B[r, c] = basis[r]^(q^(c - 1))
        end
    end
    B_inv = zero_matrix(E, m, m)
    for r in 1:m
        for c in 1:m
            B_inv[r, c] = dual_basis[c]^(q^(r - 1))
        end
    end
    return B * B_inv == identity_matrix(E, m)
end
verify_complementary_basis(E::CTFieldTypes, F::CTFieldTypes, basis::Vector{<:CTFieldElem}, dual_basis::Vector{<:CTFieldElem}) = verify_dual_basis(E, F, basis, dual_basis)

"""
    are_equivalent_basis(basis::Vector{fqPolyRepFieldElem}, basis2::Vector{fqPolyRepFieldElem})

Return `true` if `basis` is a scalar multiple of `basis2`.
"""
function are_equivalent_basis(basis::Vector{<:CTFieldElem}, basis2::Vector{<:CTFieldElem})
    m = length(basis)
    length(basis2) == m || throw(ArgumentError("The two vectors must have the same length."))
    c = basis[1] * basis2[1]^-1
    E = parent(basis[1])
    for i in 1:m
        parent(basis[i]) == E || throw(ArgumentError("Elements must be over the same field."))
        parent(basis2[i]) == E || throw(ArgumentError("Elements must be over the same field."))
        # for logical consistency should probably do this in a separate loop
        basis[i] == c * basis2[i] || return false
    end
    return true
end

"""
    is_self_dual_basis(E::fqPolyRepField, F::fqPolyRepField, basis::Vector{fqPolyRepFieldElem})

Return `true` if `basis` is equal to its dual.
"""
function is_self_dual_basis(E::CTFieldTypes, F::CTFieldTypes, basis::Vector{<:CTFieldElem})
    flag, λ = is_basis(E, F, basis)
    flag || throw(ArgumentError("The provided vector is not a basis for the extension."))
    return basis == λ
end

"""
    is_primitive_basis(E::fqPolyRepField, F::fqPolyRepField, basis::Vector{fqPolyRepFieldElem})

Return `true` if `basis` is a primitive basis for `E/F`.
"""
function is_primitive_basis(E::CTFieldTypes, F::CTFieldTypes, basis::Vector{<:CTFieldElem})
    flag, _ = is_basis(E, F, basis)
    flag || throw(ArgumentError("The provided vector is not a basis for the extension."))
    isone(basis[1]) ? (x = basis[2];) : (x = basis[1];)
    m = length(basis)
    for i in 0:m - 1
        x^i ∈ basis || return false
    end
    return true
end

"""
    is_normal_basis(E::fqPolyRepField, F::fqPolyRepField, basis::Vector{fqPolyRepFieldElem})

Return `true` if `basis` is a normal basis for `E/F`.
"""
function is_normal_basis(E::CTFieldTypes, F::CTFieldTypes, basis::Vector{<:CTFieldElem})
    flag, _ = is_basis(E, F, basis)
    flag || throw(ArgumentError("The provided vector is not a basis for the extension."))
    isone(basis[1]) ? (return false;) : (x = basis[1];)
    m = length(basis)
    q = Int(order(F))
    for i in 0:m - 1
        x^(q^i) ∈ basis || return false
    end
    return true
end

# The finite field F_q^n has a pair of self-dual bases for the following parameters.
# (1) q is an even prime power.
# (2) q is an odd prime power and n = 2k + 1.
# (Imamura 1983) The finite field F_q^n has no self-dual power bases.

#############################
          # Graphs
#############################

"""
    is_regular(G::SimpleGraph{Int})

Return `true` if `G` is regular.
"""
function is_regular(G::SimpleGraph{Int})
    deg = length(G.fadjlist[1])
    all(length(v) == deg for v in G.fadjlist)
end

"""
    edge_vertex_incidence_matrix(G::SimpleGraph{Int})

Return the edge-vertex incidence matrix of `G` along with the vertex incides of the left
and right bipartition.
"""
function edge_vertex_incidence_matrix(G::SimpleGraph{Int})
    I = Array(Grphs.incidence_matrix(G))
    nr, nc = size(I)
    I_tr = transpose(I)
    B = vcat(hcat(zeros(Int, nc, nc), I_tr), hcat(I, zeros(Int, nr, nr)))
    return B, collect(1:nc), collect(nc + 1:nr + nc)
end

"""
    edge_vertex_incidence_graph(G::SimpleGraph{Int})

Return the edge-vertex incidence graph of `G` along with the vertex incides of the left
and right bipartition.
"""
function edge_vertex_incidence_graph(G::SimpleGraph{Int})
    B, left, right = edge_vertex_incidence_matrix(G)
    return SimpleGraph(B), left, right
end

"""
    is_valid_bipartition(G::SimpleGraph{Int}, left::Vector{Int}, right::Vector{Int})

Return `true` if the vertices indexed by `left` and `right` form a valid bipartition for `G`.
"""
function is_valid_bipartition(G::SimpleGraph{Int}, left::Vector{Int}, right::Vector{Int})
    Grphs.nv(G) == length(left) + length(right) || throw(ArgumentError("Too few vertices in lists."))
    l = sort(left)
    r = sort(right)
    temp = l ∩ r # can manually do using sorted knowledge if this is slow
    isempty(temp) || throw(ArgumentError("Bipartition must be disjoint."))
    # could consider getting rid of the above in exchange for changing the elseif to an if

    for (i, v) in enumerate(G.fadjlist)
        if insorted(i, l)
            for e in v
                insorted(e, r) || return false
            end
        elseif insorted(i, r)
            for e in v
                insorted(e, l) || return false
            end
        else
            return false
        end
    end
    return true
end

"""
    extract_bipartition(G::SimpleGraph{Int})

Return two vectors representing the vertex indices of each side of the bipartition.
"""
function extract_bipartition(G::SimpleGraph{Int})
    temp = bipartite_map(G)
    # this is the definition of the function is_bipartite in Graphs.jl
    length(temp) == nv(G) || throw(ArgumentError("Input graph is not bipartite."))
    left = Vector{Int}()
    right = Vector{Int}()
    for i in temp
        i == 1 ? (push!(i, left);) : (push!(i, right);)
    end
    return left, right
end

function _rand_invertible_matrix(F::CTFieldTypes, n::Integer)
    n > 0 || throw(DomainError(n, "The dimension `n` must be positive."))

    # start with random, nonzero, 1×1 matrix
    A = matrix(F, 1, 1, [rand(collect(F)[2:end])])

    # extend from (k-1)×(k-1) to k×k, repeat up to n×n
    for k in 2:n

        # pick a random nonzero vector of length n
        v = zero_matrix(F, 1, k)
        while iszero(v)
            v = matrix(F, 1, k, rand(F, k))
        end

        # pick random nonzero index of v
        r = rand(findall(!iszero, v))[2]

        # create identity matrix with r-th row replaced by v (note, this is invertible)
        I_v = identity_matrix(F, k)
        I_v[r, :] = v

        # copy data from A into the right places
        B = zero_matrix(F, k, k)
        B[1, r] = one(F)
        B[2:end, 1:r - 1] = A[:, 1:r - 1]
        B[2:end, r + 1:end] = A[:, r:end]

        A = B * I_v
    end

    return A
end

function extended_binomial(x::UInt, y::UInt)
  z = UInt(0)
  if y <= x
    z = binomial(big(x), big(y))
  end
  return z
end

# #=
# Example of using the repeated iterator inside of product.
#
# It turns out that this is faster than the Nemo iterator and doesn't allocate.
#
# julia> @benchmark for i in Base.Iterators.product(Base.Iterators.repeated(0:1, 10)...) i end
# BenchmarkTools.Trial: 10000 samples with 137 evaluations.
#  Range (min … max):  713.022 ns …  1.064 μs  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     755.949 ns              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   760.380 ns ± 24.121 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%
#
#  Memory estimate: 0 bytes, allocs estimate: 0.
#
# julia> @benchmark for i in Nemo.AbstractAlgebra.ProductIterator([0:1 for _ in 1:10]) i end
# BenchmarkTools.Trial: 10000 samples with 1 evaluation.
#  Range (min … max):  34.064 μs …   2.604 ms  ┊ GC (min … max):  0.00% … 97.51%
#  Time  (median):     36.970 μs               ┊ GC (median):     0.00%
#  Time  (mean ± σ):   46.342 μs ± 124.916 μs  ┊ GC (mean ± σ):  16.57% ±  6.04%
#
#  Memory estimate: 176.50 KiB, allocs estimate: 2051.
#
# julia> @benchmark for i in Base.Iterators.product([0:1 for _ in 1:10]...) i end
# BenchmarkTools.Trial: 10000 samples with 1 evaluation.
#  Range (min … max):  53.741 μs …   1.465 ms  ┊ GC (min … max):  0.00% … 87.86%
#  Time  (median):     63.790 μs               ┊ GC (median):     0.00%
#  Time  (mean ± σ):   76.919 μs ± 104.655 μs  ┊ GC (mean ± σ):  12.40% ±  8.59%
#
#  Memory estimate: 432.88 KiB, allocs estimate: 2061.
# =#
#
#
# # Gray code iterator, naive formula, gives Ints instead of vectors
#
# struct GrayCodeNaive
#     n::Int
# end
#
# Base.iterate(G::GrayCodeNaive) = G.n < 64 ? (0, 1) : error("Don't handle cases this large")
#
# function Base.iterate(G::GrayCodeNaive, k)
#     k == 2^G.n && return nothing
#     return (k ⊻ (k >> 1), k + 1)
# end
#
# Base.length(G::GrayCodeNaive) = 2^G.n
#
# #=
# Benchmark result:
#
# julia> @benchmark for g in GrayCodeNaive(25) g end
# BenchmarkTools.Trial: 25 samples with 1 evaluation.
#  Range (min … max):  200.014 ms … 202.460 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     200.305 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   200.626 ms ± 659.369 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
#
#  Memory estimate: 0 bytes, allocs estimate: 0.
# =#
#
#
#
# # Gray code iterator, chooses next based on previous, gives Ints instead of vectors
#

# #=
# Benchmark result:
#
# julia> @benchmark for g in GrayCode(25) g end
# BenchmarkTools.Trial: 15 samples with 1 evaluation.
#  Range (min … max):  348.661 ms … 349.411 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     349.050 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   349.048 ms ± 225.710 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
#
#  Memory estimate: 0 bytes, allocs estimate: 0.
# =#
