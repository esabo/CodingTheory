# Copyright (c) 2021, 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
  # Generic Helper Functions
#############################

_hasequivalentrowspaces(A::T, B::T) where T <: CTMatrixTypes =
    return _removeempty(rref(deepcopy(A))[2], :rows) == _removeempty(rref(deepcopy(B))[2], :rows)

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
    directsum(A::CTMatrixTypes, B::CTMatrixTypes)

Return the direct sum of the two matrices `A` and `B`.
"""
function ⊕(A::T, B::T) where T <: CTMatrixTypes
    base_ring(A) == base_ring(B) || throw(ArgumentError("Matrices must be over the same base ring in directsum."))

    return vcat(hcat(A, zero_matrix(base_ring(B), nrows(A), ncols(B))),
        hcat(zero_matrix(base_ring(A), nrows(B), ncols(A)), B))
end
directsum(A::T, B::T) where T <: CTMatrixTypes = A ⊕ B

"""
    ⊗(A::CTMatrixTypes, B::CTMatrixTypes)
    kron(A::CTMatrixTypes, B::CTMatrixTypes)
    tensorproduct(A::CTMatrixTypes, B::CTMatrixTypes)
    kroneckerproduct(A::CTMatrixTypes, B::CTMatrixTypes)

Return the Kronecker product of the two matrices `A` and `B`.
"""
⊗(A::CTMatrixTypes, B::CTMatrixTypes) = kronecker_product(A, B)
kron(A::CTMatrixTypes, B::CTMatrixTypes) = kronecker_product(A, B)
tensorproduct(A::CTMatrixTypes, B::CTMatrixTypes) = kronecker_product(A, B)
kroneckerproduct(A::CTMatrixTypes, B::CTMatrixTypes) = kronecker_product(A, B)

# I think we should avoid length checking here and return it for entire matrix if given
# Hammingweight(v::T) where T <: Union{CTMatrixTypes, gfp_mat, Vector{S}} where S <: Integer = count(i->(i != 0), v)
"""
    Hammingweight(v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer
    weight(v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer
    wt(v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer

Return the Hamming weight of `v`.
"""
function Hammingweight(v::T) where T <: Union{CTMatrixTypes, Vector{<:CTFieldElem}, Vector{S}, Adjoint{S, Vector{S}}} where S <: Integer
    count(x -> !iszero(x), v)
end
weight(v::T) where T <: Union{CTMatrixTypes, Vector{<:CTFieldElem}, Vector{S}, Adjoint{S, Vector{S}}} where S <: Integer = Hammingweight(v)
wt(v::T) where T <: Union{CTMatrixTypes, Vector{<:CTFieldElem}, Vector{S}, Adjoint{S, Vector{S}}} where S <: Integer = Hammingweight(v)

# Hammingweight(v::Matrix{Int}) = return sum(v)

# wt(v::Matrix{Int}) = count(x -> !iszero(x), v)

"""
    wt(f::CTPolyRingElem)

Return the number of nonzero coefficients of the polynomial `f`.
"""
wt(f::CTPolyRingElem) = Hammingweight(collect(coefficients(f)))

# """
#     _minwtrow(A::CTMatrixTypes)
#
# Return the minimum weight and corresponding index of the rows of `A`.
# """
function _minwtrow(A::Union{CTMatrixTypes, Matrix{S}, LinearAlgebra.Adjoint{S, Matrix{S}}}) where S <: Integer
    w = size(A, 2) + 1
    i = 0
    for r in axes(A, 1)
        wloc = 0
        for c in axes(A, 2)
            iszero(A[r,c]) || (wloc += 1)
        end
        wloc < w && (w = wloc; i = r;)
    end
    return w, i
end

function _minwtcol(A::Union{CTMatrixTypes, Matrix{S}, LinearAlgebra.Adjoint{S, Matrix{S}}}) where S <: Integer
    w = size(A, 1) + 1
    i = 0
    for c in axes(A, 2)
        wloc = 0
        for r in axes(A, 1)
            iszero(A[r,c]) || (wloc += 1)
        end
        wloc < w && (w = wloc; i = c;)
    end
    return w, i
end

"""
    Hammingdistance(u::T, v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer
    distance(u::T, v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer
    dist(u::T, v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer

Return the Hamming distance between `u` and `v`.
"""
Hammingdistance(u::T, v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer = Hammingweight(u - v)
distance(u::T, v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer = Hammingweight(u - v)
dist(u::T, v::T) where T <: Union{CTMatrixTypes, Vector{S}} where S <: Integer = Hammingweight(u - v)

"""
    symplecticinnerproduct(u::CTMatrixTypes, v::CTMatrixTypes)

Return the symplectic inner product of `u` and `v`.
"""
function symplecticinnerproduct(u::CTMatrixTypes, v::CTMatrixTypes)
    (nrows(u) == 1 || ncols(u) == 1) || throw(ArgumentError("First argument of symplectic inner product is not a vector: dims = $(size(u, 1))."))
    (nrows(v) == 1 || ncols(v) == 1) || throw(ArgumentError("Second argument of symplectic inner product is not a vector: dims = $(size(v, 1))."))
    length(u) == length(v) || throw(ArgumentError("Vectors must be the same length in symplectic inner product."))
    iseven(length(u)) || throw(ArgumentError("Vectors must have even length in symplectic inner product."))
    base_ring(u) == base_ring(v) || throw(ArgumentError("Vectors must be over the same field in symplectic inner product."))
    
    ncols = div(length(u), 2)
    return sum(u[i + ncols] * v[i] - v[i + ncols] * u[i] for i in 1:ncols)
end

"""
    aresymplecticorthogonal(A::CTMatrixTypes, B::CTMatrixTypes)

Return `true` if the rows of the matrices `A` and `B` are symplectic orthogonal.
"""
function aresymplecticorthogonal(A::CTMatrixTypes, B::CTMatrixTypes)
    base_ring(A) == base_ring(B) || throw(ArgumentError("Matrices in product must both be over the same base ring."))
    
    return iszero(hcat(A[:, div(ncols(A), 2) + 1:end], -A[:, 1:div(ncols(A), 2)]) * transpose(B))
end

# function traceinnerproduct(u::CTMatrixTypes, v::CTMatrixTypes)
#
# end

"""
    Hermitianinnerproduct(u::CTMatrixTypes, v::CTMatrixTypes)

Return the Hermitian inner product of `u` and `v`.
"""
function Hermitianinnerproduct(u::CTMatrixTypes, v::CTMatrixTypes)
    (nrows(u) == 1 || ncols(u) == 1) || throw(ArgumentError("First argument of Hermitian inner product is not a vector: dims = $(size(u, 1))."))
    (nrows(v) == 1 || ncols(v) == 1) || throw(ArgumentError("Second argument of Hermitian inner product is not a vector: dims = $(size(v, 1))."))
    length(u) == length(v) || throw(ArgumentError("Vectors must be the same length in Hermitian inner product."))
    base_ring(u) == base_ring(v) || throw(ArgumentError("Vectors must be over the same field in Hermitian inner product."))
    q2 = order(base_ring(u))
    issquare(q2) || throw(ArgumentError("The Hermitian inner product is only defined over quadratic field extensions."))
    
    q = Int(sqrt(q2, check = false))
    return sum(u[i] * v[i]^q for i in 1:length(u))
end

"""
    Hermitianconjugatematrix(A::CTMatrixTypes)

Return the Hermitian conjugate of the matrix `A`.
"""
function Hermitianconjugatematrix(A::CTMatrixTypes)
    q2 = order(base_ring(A))
    issquare(q2) || throw(ArgumentError("The Hermitian conjugate is only defined over quadratic field extensions."))

    q = Int(sqrt(q2, check = false))
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

# """
#     FpmattoJulia(M::CTMatrixTypes)

# Return the `CTMatrixTypes` matrix `M` as a Julia Int matrix.
# """
# TODO: want to remove and cease use of FpmattoJulia
function FpmattoJulia(M::CTMatrixTypes)
    degree(base_ring(M)) == 1 || throw(ArgumentError("Cannot promote higher order elements to the Ints."))

    A = zeros(Int, size(M))
    for c in 1:ncols(M)
        for r in 1:nrows(M)
            A[r, c] = coeff(M[r, c], 0)
        end
    end
    return A
end
FpmattoJulia(M::fpMatrix) = data.(M)

function _nonpivotcols(A::CTMatrixTypes, type::Symbol=:nsp)
    type ∈ [:sp, :nsp]
    if type == :sp
        return setdiff(collect(1:ncols(A)), [x.pos[1] for x in A])
    else #if type == :nsp - not sparse
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

function _quotientspace(big::T, small::T, alg::Symbol=:syseqs) where T <: CTMatrixTypes
    alg ∈ [:VS, :syseqs] || throw(ArgumentError("Unknown algorithm type"))

    if alg == :VS
        F = base_ring(big)
        V = VectorSpace(F, ncols(big))
        U, UtoV = sub(V, [V(small[i, :]) for i in 1:nrows(small)])
        W, WtoV = sub(V, [V(big[i, :]) for i in 1:nrows(big)])
        gensofUinW = Vector{typeof(gens(U)[1])}(undef, length(gens(U)))
        # gensofUinW = [preimage(WtoV, UtoV(g)) for g in gens(U)]
        Threads.@threads for i in 1:length(gens(U))
            gensofUinW[i] = preimage(WtoV, UtoV(gens(U)[i]))
        end
        UinW, _ = sub(W, gensofUinW)
        Q, WtoQ = quo(W, UinW)
        iszero(dim(Q)) && (return zero_matrix(F, 1, ncols(big));)
        C2modC1basis = [WtoV(x) for x in [preimage(WtoQ, g) for g in gens(Q)]]
        Fbasis = [[F(C2modC1basis[j][i]) for i in 1:AbstractAlgebra.dim(parent(C2modC1basis[1]))] for j in 1:length(C2modC1basis)]
        return matrix(F, length(Fbasis), length(Fbasis[1]), reduce(vcat, Fbasis))
    else
        # solve the system x big = small
        # sol contains the way to write the rows of small in terms of the rows of big
        # if big is of the form (big = small ∪ (big / small)), then this will have zeros for the rows
        # corresponding to the basis of the quotient
        # in the general case, anything without a pivot in the row reduction is not required to make
        # the elements of small and therefore lie in the quotient space
        flag, sol = can_solve_with_solution(big, small, side=:left)
        !flag && error("Cannot solve system for quotient")
        _, rrefsol = rref(sol)
        if typeof(rrefsol) <: SMat{W, Vector{W}} where W <: CTFieldElem
            nonpivots = _nonpivotcols(rrefsol, :sp)
            return reduce(vcat, [big[r, :] for r in nonpivots])
        else
            return big[_nonpivotcols(rrefsol, :nsp), :]
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

function _removeempty(A::CTMatrixTypes, type::Symbol)
    type ∈ [:rows, :cols] || throw(ArgumentError("Unknown type in _removeempty"))
    
    del = Vector{Int}()
    if type == :rows
        for r in axes(A, 1)
            if iszero(A[r, :])
                append!(del, r)
            end
        end
        return isempty(del) ? A : A[setdiff(1:nrows(A), del), :]
    elseif type == :cols
        for c in axes(A, 2)
            if iszero(A[:, c])
                append!(del, c)
            end
        end
        return isempty(del) ? A : A[:, setdiff(1:ncols(A), del)]
    end
end

function _rref_no_col_swap(M::CTMatrixTypes, rowrange::UnitRange{Int}, colrange::UnitRange{Int})
    A = deepcopy(M)
    _rref_no_col_swap!(A, rowrange, colrange)
    return A
end
_rref_no_col_swap(M::CTMatrixTypes, rowrange::Base.OneTo{Int}, colrange::Base.OneTo{Int}) = _rref_no_col_swap(M, 1:rowrange.stop, 1:colrange.stop)
_rref_no_col_swap(M::CTMatrixTypes) = _rref_no_col_swap(M, axes(M, 1), axes(M, 2))

function _rref_no_col_swap!(A::CTMatrixTypes, rowrange::UnitRange{Int}, colrange::UnitRange{Int})
    # don't do anything to A if the range is empty
    isempty(rowrange) && return nothing
    isempty(colrange) && return nothing

    i = rowrange.start
    j = colrange.start
    nr = rowrange.stop
    nc = colrange.stop
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
                for k in rowrange
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
                for k in rowrange
                    if k != i
                        # do a manual loop here to reduce allocations
                        @simd for l in axes(A, 2)
                            A[k, l] = (A[k, l] - A[i, l])
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

function _rref_col_swap(M::CTMatrixTypes, rowrange::UnitRange{Int}, colrange::UnitRange{Int})
    A = deepcopy(M)
    rnk, P = _rref_col_swap!(A, rowrange, colrange)
    return rnk, A, P
end
_rref_col_swap(M::CTMatrixTypes, rowrange::Base.OneTo{Int}, colrange::Base.OneTo{Int}) = _rref_col_swap(M, 1:rowrange.stop, 1:colrange.stop)
_rref_col_swap(M::CTMatrixTypes) = _rref_col_swap(M, axes(M, 1), axes(M, 2))

function _rref_col_swap!(A::CTMatrixTypes, rowrange::UnitRange{Int}, colrange::UnitRange{Int})
    # don't do anything to A if the range is empty, return rank 0 and missing permutation matrix
    isempty(rowrange) && return 0, missing
    isempty(colrange) && return 0, missing

    # permutation matrix required to return to rowspace if column swap done
    P = missing
    ncA = ncols(A)

    rnk = 0
    i = rowrange.start
    j = colrange.start
    nr = rowrange.stop
    nc = colrange.stop
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
                            ismissing(P) && (P = identity_matrix(base_ring(A), ncA);)
                            swap_cols!(A, k, j)
                            swap_cols!(P, k, j)
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
                for k = rowrange.start:nr
                    if k != i
                        # do a manual loop here to reduce allocations
                        d = A[k, j]
                        @simd for l = 1:ncA
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
                            ismissing(P) && (P = identity_matrix(base_ring(A), ncA);)
                            swap_cols!(A, k, j)
                            swap_cols!(P, k, j)
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
                for k = rowrange.start:nr
                    if k != i
                        # do a manual loop here to reduce allocations
                        @simd for l = 1:ncA
                            A[k, l] = (A[k, l] - A[i, l])
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

function _rref_symp_col_swap!(A::CTMatrixTypes, rowrange::UnitRange{Int}, colrange::UnitRange{Int})
    # don't do anything to A if the range is empty, return rank 0 and missing permutation matrix
    isempty(rowrange) && return 0, missing
    isempty(colrange) && return 0, missing

    # permutation matrix required to return to rowspace if column swap done
    P = missing
    ncA = ncols(A)

    rnk = 0
    i = rowrange.start
    j = colrange.start
    nr = rowrange.stop
    nc = colrange.stop
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
                            ismissing(P) && (P = identity_matrix(base_ring(A), ncA);)
                            k_symp = mod1(k + div(ncA, 2), ncA)
                            j_symp = mod1(j + div(ncA, 2), ncA)
                            swap_cols!(A, k, j)
                            swap_cols!(P, k, j)
                            swap_cols!(A, k_symp, j_symp)
                            swap_cols!(P, k_symp, j_symp)
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
                for k = rowrange.start:nr
                    if k != i
                        # do a manual loop here to reduce allocations
                        d = A[k, j]
                        @simd for l = 1:ncA
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
                            ismissing(P) && (P = identity_matrix(base_ring(A), ncA);)
                            k_symp = mod1(k + div(ncA, 2), ncA)
                            j_symp = mod1(j + div(ncA, 2), ncA)
                            swap_cols!(A, k, j)
                            swap_cols!(P, k, j)
                            swap_cols!(A, k_symp, j_symp)
                            swap_cols!(P, k_symp, j_symp)
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
                for k = rowrange.start:nr
                    if k != i
                        # do a manual loop here to reduce allocations
                        @simd for l = 1:ncA
                            A[k, l] = (A[k, l] - A[i, l])
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

function digitstoint(x::Vector{Int}, base::Int=2)
    res = 0
    for digit in x
        res = digit + base * res
    end
    return res
end

"""
    polytocircmatrix(f::AbstractAlgebra.Generic.Res{T}) where T <: CTPolyRingElem

Return the circulant matrix whose first column is the coefficients of `f`.
"""
function polytocircmatrix(f::AbstractAlgebra.Generic.Res{T}) where T <: CTPolyRingElem
    R = parent(f)
    S = base_ring(R)
    F = base_ring(S)
    g = modulus(R)
    l = degree(g)
    g == gen(S)^l - 1 || throw(ArgumentError("Residue ring not of the form x^l - 1."))
    # gcd(l, Int(characteristic(F))) == 1 || throw(ArgumentError("Residue ring over F_q[x] must have modulus x^l - 1 with gcd(l, q) = 1."))

    A = zero_matrix(F, l, l)
    fcoeffs = zero_matrix(F, l, 1)
    temp = collect(coefficients(Nemo.lift(f)))
    fcoeffs[1:length(temp), 1] = temp
    A[:, 1] = fcoeffs
    for c in 2:l
        # A[:, c] = circshift(fcoeffs, c - 1)
        A[1:c-1, c] = fcoeffs[l - (c - 1) + 1:l, 1]
        A[c:end, c] = fcoeffs[1:l - (c - 1), 1]
    end
    return A
end

"""
    lift(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{T}}) where T <: CTPolyRingElem

Return the matrix whose polynomial elements are converted to circulant matrices over the base field.
"""
function lift(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{T}}) where T <: CTPolyRingElem
    R = parent(A[1, 1])
    S = base_ring(R)
    F = base_ring(S)
    g = modulus(R)
    l = degree(g)
    g == gen(S)^l - 1 || throw(ArgumentError("Residue ring not of the form x^l - 1."))
    # gcd(l, Int(characteristic(F))) == 1 || throw(ArgumentError("Residue ring over F_q[x] must have modulus x^l - 1 with gcd(l, q) = 1."))

    nr, nc = size(A)
    Alift = zero_matrix(F, nr * l, nc * l)
    for c in axes(A, 2)
        for r in axes(A, 1)
            if !iszero(A[r, c])
                Alift[(r - 1) * l + 1:r * l, (c - 1) * l + 1:c * l] = polytocircmatrix(A[r, c])
            end
        end
    end
    return Alift
end

# Creates a matrix with copies of `M` at every nonzero entry of `locations`.
function _concat(locations::Union{CTMatrixTypes, Matrix}, M::CTMatrixTypes)
    nrM, ncM = size(M)
    nrL, ncL = size(locations)
    output = zero_matrix(base_ring(M), nrM * nrL, ncM * ncL)
    for jouter in 1:ncL
        for iouter in 1:nrL
            if !iszero(locations[iouter, jouter])
                for j in 1:ncM
                    for i in 1:nrM
                        row = i + nrM * (iouter - 1)
                        col = j + ncM * (jouter - 1)
                        output[row, col] = M[i, j]
                    end
                end
            end
        end
    end
    return output
end

#############################
  # Quantum Helper Functions
#############################

"""
    istriorthogonal(G::CTMatrixTypes, verbose::Bool=false)
    istriorthogonal(G::Matrix{Int}, verbose::Bool=false)

Return `true` if the binary matrix `G` is triorthogonal.

# Notes
* If the optional parameter `verbos` is set to `true`, the first pair or triple of
  non-orthogonal rows will be identified on the console.
"""
function istriorthogonal(G::CTMatrixTypes, verbose::Bool=false)
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

function istriorthogonal(G::Matrix{Int}, verbose::Bool=false)
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

function printstringarray(A::Vector{String}, withoutIs=false)
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
printchararray(A::Vector{Vector{Char}}, withoutIs=false) = printstringarray(setchartostringarray(A), withoutIs)
printsymplecticarray(A::Vector{Vector{T}}, withoutIs=false) where T <: Int = printstringarray(setsymplectictostringarray(A), withoutIs)

"""
    pseudoinverse(M::CTMatrixTypes)

Return the pseudoinverse of a stabilizer matrix `M` over a quadratic extension.

# Notes
* This is not the Penrose-Moore pseudoinverse.
"""
function pseudoinverse(M::CTMatrixTypes)
    # let this fail elsewhere if not actually over a quadratic extension
    if degree(base_ring(M)) != 1
        # TODO: quadratictosymplectic is no longer defined, this will need changed
        M = transpose(quadratictosymplectic(M))
    else
        M = transpose(M)
    end

    nr, nc = size(M)
    MS = MatrixSpace(base_ring(M), nr, nr)
    _, E = rref(hcat(M, MS(1)))
    E = E[:, (nc + 1):end]
    pinv = E[1:nc, :]
    dual = E[nc + 1:nr, :]

    # verify
    _, Mrref = rref(M)
    MScols = MatrixSpace(base_ring(M), nc, nc)
    E * M == Mrref || error("Pseudoinverse calculation failed (transformation incorrect).")
    Mrref[1:nc, 1:nc] == MScols(1) || error("Pseudoinverse calculation failed (failed to get I).")
    iszero(Mrref[nc + 1:nr, :]) || error("Pseudoinverse calculation failed (failed to get zero).")
    pinv * M == MScols(1) || error("Pseudoinverse calculation failed (eq 1).")
    transpose(M) * transpose(pinv) == MScols(1) || error("Pseudoinverse calculation failed (eq 2).")
    iszero(transpose(M) * transpose(dual)) || error("Failed to correctly compute dual (rhs).")
    iszero(dual * M) || error("Failed to correctly compute dual (lhs).")
    return pinv
end

# """
#     quadratictosymplectic(M::CTMatrixTypes)

# Return the matrix `M` converted from the quadratic to the symplectic form.
# """
# function quadratictosymplectic(M::CTMatrixTypes)
#     E = base_ring(M)
#     iseven(degree(E)) || error("The base ring of the given matrix is not a quadratic extension.")
#     F, _ = FiniteField(Int(characteristic(E)), div(degree(E), 2), "ω")
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
#     E, ω = FiniteField(Int(characteristic(F)), 2 * degree(F), "ω")
#     ϕ = embed(F, E)
#     Mquad = zero_matrix(E, nr, nc)
#     for c in 1:nc
#         for r in 1:nr
#             Mquad[r, c] = ϕ(M[r, c]) + ϕ(M[r, c + nc]) * ω
#         end
#     end
#     return Mquad
# end

function _Paulistringtosymplectic(str::T) where T <: Union{String, Vector{Char}}
    n = length(str)
    F, _ = FiniteField(2, 1, "ω")
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
_Paulistringtosymplectic(A::Vector{T}) where T <: Union{String, Vector{Char}} = reduce(vcat, [_Paulistringtosymplectic(s) for s in A])
# need symplectictoPaulistring

# charvec::Union{Vector{nmod}, Missing}=missing)
function _processstrings(SPauli::Vector{T}) where T <: Union{String, Vector{Char}}
    # Paulisigns = Vector{Int}()
    StrPaulistripped = Vector{String}()
    for (i, s) in enumerate(SPauli)
        if s[1] ∈ ['I', 'X', 'Y', 'Z']
            # append!(Paulisigns, 1)
            push!(StrPaulistripped, s)
        elseif s[1] == '+'
            # append!(Paulisigns, 1)
            push!(StrPaulistripped, s[2:end])
        elseif s[1] == '-'
            # append!(Paulisigns, -1)
            push!(StrPaulistripped, s[2:end])
        else
            error("The first element of Pauli string $i is neither a Pauli character or +/-: $s.")
        end
    end

    n = length(StrPaulistripped[1])
    for s in StrPaulistripped
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
    return StrPaulistripped#, charvec
end

#############################
        # Finite Fields
#############################

"""
    tr(x::fq_nmod, K::FqNmodFiniteField, verify::Bool=false)

Return the relative trace of `x` from its base field to the field `K`.

# Notes
* If the optional parameter `verify` is set to `true`, the two fields are checked
  for compatibility.
"""
function tr(x::CTFieldElem, K::CTFieldTypes, verify::Bool=false)
    L = parent(x)
    q = order(K)
    if verify
        # # shouldn't need Int casting here but just in case...
        # Int(characteristic(L)) == Int(characteristic(K)) || error("The given field is not a subfield of the base ring of the element.")
        # degree(L) % degree(K) == 0 || error("The given field is not a subfield of the base ring of the element.")
        flag, m = isextension(L, K)
        flag || throw(ArgumentError("The given field is not a subfield of the base ring of the matrix."))
    end
    n = div(degree(L), degree(K))
    return sum([x^(q^i) for i in 0:(n - 1)])
end

function _expandelement(x::CTFieldElem, K::CTFieldTypes, basis::Vector{<:CTFieldElem}, verify::Bool=false)
    return [tr(x * i) for i in basis] #, K, verify
end

function _expandrow(row::CTMatrixTypes, K::CTFieldTypes, basis::Vector{<:CTFieldElem}, verify::Bool=false)
    new_row = _expandelement(row[1], K, basis, verify)
    for i in 2:ncols(row)
        new_row = vcat(new_row, _expandelement(row[i], K, basis, verify))
    end
    return matrix(K, 1, length(new_row), new_row)
end

"""
    expandmatrix(M::CTMatrixTypes, K::FqNmodFiniteField, basis::Vector{fq_nmod})

Return the matrix constructed by expanding the elements of `M` to the subfield
`K` using the provided `basis` for the base ring of `M` over `K`.
"""
function expandmatrix(M::CTMatrixTypes, K::CTFieldTypes, basis::Vector{<:CTFieldElem})
    L = base_ring(M)
    L == K && return M
    Int(characteristic(L)) == Int(characteristic(K)) || throw(ArgumentError("The given field is not a subfield of the base ring of the element."))
    degree(L) % degree(K) == 0 || throw(ArgumentError("The given field is not a subfield of the base ring of the element."))
    n = div(degree(L), degree(K))
    n == length(basis) || throw(ArgumentError("Provided basis is of incorrect size for the given field and subfield."))
    # should really check if it is a basis
    flag, m = isextension(L, K)
    flag || throw(ArgumentError("The given field is not a subfield of the base ring of the matrix."))
    m == length(basis) || throw(ArgumentError("Basis does not have length degree of the extension."))
    flag, _ = _isbasis(L, basis, Int(order(K)))
    flag || throw(ArgumentError("The provided vector is not a basis for the extension."))

    return reduce(vcat, [_expandrow(M[r, :], K, basis) for r in 1:nrows(M)])
end

"""
    quadraticresidues(q::Int, n::Int)

Return the sets of quadratic resides and quadratic non-residues of `q` and `n`.
"""
function quadraticresidues(q::Int, n::Int)
    isodd(n) && isprime(n) || throw(ArgumentError("n must be an odd prime in quadratic residues"))
    q^div(n - 1, 2) % n == 1 || throw(ArgumentError("q^(n - 1)/2 ≅ 1 mod n in quadratic residues"))

    # F, _ = FiniteField(n, 1, "α")
    # elms = collect(F)
    # # skip 0
    # qres = unique!([i^2 for i in elms[2:end]])
    # nqres = setdiff!(elms[2:end], qres)
    # return qres, nqres

    # don't want this returning in the field
    qres = sort!(unique!([i^2 % n for i in 1:n - 1]))
    nqres = setdiff(1:(n - 1), qres)
    return qres, nqres
end

function _isbasis(E::CTFieldTypes, basis::Vector{<:CTFieldElem}, q::Int)
    m = length(basis)
    B = zero_matrix(E, m, m)
    for r in 1:m
        for c in 1:m
            B[r, c] = basis[r]^(q^(c - 1))
        end
    end
    iszero(det(B)) && return false, missing
    
    try
        Binv = inv(B)
        λ = [Binv[1, i] for i in 1:m]
        return true, λ
    catch
        return false, missing
    end
end

"""
    isextension(E::FqNmodFiniteField, F::FqNmodFiniteField)

Return `true` if `E/F` is a valid field extension.
"""
function isextension(E::CTFieldTypes, F::CTFieldTypes)
    p = Int(characteristic(E))
    Int(characteristic(F)) == p || return false, missing
    degE = degree(E)
    degF = degree(F)
    degE % degF == 0 || return false, missing
    return true, div(degE, degF)
    # the below allows you to embed GF(2) into GF(5) without error but is not an extension
    # try
    #     embed(F, E)
    #     return true, div(degree(E), degree(F))
    # catch
    #     return false, missing
    # end
end

"""
    isbasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod})

Return `true` and the dual (complementary) basis if `basis` is a basis for `E/F`,
otherwise return `false, missing`.
"""
function isbasis(E::CTFieldTypes, F::CTFieldTypes, basis::Vector{<:CTFieldElem})
    flag, m = isextension(E, F)
    flag || throw(ArgumentError("Second field is not a subfield of the first."))
    length(basis) == m || throw(ArgumentError("Basis does not have length degree of the extension."))
    for i in 1:m
        parent(basis[i]) == E || throw(ArgumentError("The basis must be elements of the extension field."))
    end

    return _isbasis(E, basis, Int(order(F)))
end

"""
    primitivebasis(E::FqNmodFiniteField, F::FqNmodFiniteField)

Return a primitive basis for `E/F` and its dual (complementary) basis.
"""
function primitivebasis(E::CTFieldTypes, F::CTFieldTypes)
    flag, m = isextension(E, F)
    flag || throw(ArgumentError("Second field is not a subfield of the first."))
    α = gen(E)
    basis = [α^i for i in 0:m - 1]
    flag, λ = _isbasis(E, basis, Int(order(F)))
    return basis, λ
end
# these are slightly different
# polynomialbasis(E::FqNmodFiniteField, F::FqNmodFiniteField) = primitivebasis(E, F)
# monomialbasis(E::FqNmodFiniteField, F::FqNmodFiniteField) = primitivebasis(E, F)

"""
    normalbasis(E::FqNmodFiniteField, F::FqNmodFiniteField)

Return a normal basis for `E/F` and its dual (complementary) basis.
"""
# "Normal Bases over Finite Fields" by Shuhong Gao has algorithms for this but they are
# complicated for the field sizes intended in this work
function normalbasis(E::CTFieldTypes, F::CTFieldTypes)
    flag, m = isextension(E, F)
    flag || throw(ArgumentError("Second field is not a subfield of the first."))

    q = Int(order(F))
    elms = collect(E)
    for e in elms
        basis = [e^(q^i) for i in 0:m - 1]
        flag, dualbasis = _isbasis(E, basis, q)
        flag && return basis, dualbasis
    end
    error("Somehow failed to final a normal element for the extension.")
end

"""
    dualbasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod})
    complementarybasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod})

Return the dual (complentary) basis of `basis` for the extension `E/F`.
"""
function dualbasis(E::CTFieldTypes, F::CTFieldTypes, basis::Vector{<:CTFieldElem})
    flag, λ = isbasis(E, F, basis)
    flag || throw(ArgumentError("The provided vector is not a basis for the extension."))
    return λ
end
complementarybasis(E::CTFieldTypes, F::CTFieldTypes, basis::Vector{<:CTFieldElem}) = dualbasis(E, F, basis)

"""
    verifydualbasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod}, dualbasis::Vector{fq_nmod})
    verifycomplementarybasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod}, dualbasis::Vector{fq_nmod})

Return `true` if `basis` is the dual of `dualbasis` for `E/F`, otherwise return `false`.
"""
function verifydualbasis(E::CTFieldTypes, F::CTFieldTypes, basis::Vector{<:CTFieldElem}, dualbasis::Vector{<:CTFieldElem})
    flag, m = isextension(E, F)
    flag || throw(ArgumentError("Second field is not a subfield of the first."))

    m = length(basis)
    length(dualbasis) == m || throw(ArgumentError("The basis and dual basis must have the same length."))
    E = parent(basis[1])
    for i in 1:m
        parent(basis[i]) == E || throw(ArgumentError("Elements must be over the same field."))
        parent(dualbasis[i]) == E || throw(ArgumentError("Elements must be over the same field."))
    end

    q = Int(order(F))
    M = MatrixSpace(E, m, m)
    B = M(0)
    for r in 1:m
        for c in 1:m
            B[r, c] = basis[r]^(q^(c - 1))
        end
    end
    Binv = M(0)
    for r in 1:m
        for c in 1:m
            Binv[r, c] = dualbasis[c]^(q^(r - 1))
        end
    end
    return B * Binv == M(1)
end
verifycomplementarybasis(E::CTFieldTypes, F::CTFieldTypes, basis::Vector{<:CTFieldElem}, dualbasis::Vector{<:CTFieldElem}) = verifydualbasis(E, F, basis, dualbasis)

"""
    areequivalentbasis(basis::Vector{fq_nmod}, basis2::Vector{fq_nmod})

Return `true` if `basis` is a scalar multiple of `basis2`.
"""
function areequivalentbasis(basis::Vector{<:CTFieldElem}, basis2::Vector{<:CTFieldElem})
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
    isselfdualbasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod})

Return `true` if `basis` is equal to its dual.
"""
function isselfdualbasis(E::CTFieldTypes, F::CTFieldTypes, basis::Vector{<:CTFieldElem})
    flag, λ = isbasis(E, F, basis)
    flag || throw(ArgumentError("The provided vector is not a basis for the extension."))
    return basis == λ
end

"""
    isprimitivebasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod})

Return `true` if `basis` is a primitive basis for `E/F`.
"""
function isprimitivebasis(E::CTFieldTypes, F::CTFieldTypes, basis::Vector{<:CTFieldElem})
    flag, _ = isbasis(E, F, basis)
    flag || throw(ArgumentError("The provided vector is not a basis for the extension."))
    isone(basis[1]) ? (x = basis[2];) : (x = basis[1];)
    m = length(basis)
    for i in 0:m - 1
        x^i ∈ basis || return false
    end
    return true
end

"""
    isnormalbasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod})

Return `true` if `basis` is a normal basis for `E/F`.
"""
function isnormalbasis(E::CTFieldTypes, F::CTFieldTypes, basis::Vector{<:CTFieldElem})
    flag, _ = isbasis(E, F, basis)
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
    isregular(G::SimpleGraph{Int})

Return `true` if `G` is regular.
"""
function isregular(G::SimpleGraph{Int})
    deg = length(G.fadjlist[1])
    all(length(v) == deg for v in G.fadjlist)
end

"""
    edgevertexincidencematrix(G::SimpleGraph{Int})

Return the edge-vertex incidence matrix of `G` along with the vertex incides of the left
and right bipartition.
"""
function edgevertexincidencematrix(G::SimpleGraph{Int})
    I = Array(incidence_matrix(G))
    nr, nc = size(I)
    Itr = transpose(I)
    B = vcat(hcat(zeros(Int, nc, nc), Itr), hcat(I, zeros(Int, nr, nr)))
    return B, collect(1:nc), collect(nc + 1:nr + nc)
end

"""
    edgevertexincidencegraph(G::SimpleGraph{Int})

Return the edge-vertex incidence graph of `G` along with the vertex incides of the left
and right bipartition.
"""
function edgevertexincidencegraph(G::SimpleGraph{Int})
    B, left, right = edgevertexincidencematrix(G)
    return SimpleGraph(B), left, right
end

"""
    isvalidbipartition(G::SimpleGraph{Int}, left::Vector{Int}, right::Vector{Int})

Return `true` if the vertices indexed by `left` and `right` form a valid bipartition for `G`.
"""
function isvalidbipartition(G::SimpleGraph{Int}, left::Vector{Int}, right::Vector{Int})
    nv(G) == length(left) + length(right) || throw(ArgumentError("Too few vertices in lists."))
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
    extractbipartition(G::SimpleGraph{Int})

Return two vectors representing the vertex indices of each side of the bipartition.
"""
function extractbipartition(G::SimpleGraph{Int})
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
