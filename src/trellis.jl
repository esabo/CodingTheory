# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
         # General
#############################

mutable struct Vertex
    label::BigInt
    prev::Int64
    next::Int64
    value::Float64
    edge_loc::Int64
    polynomial::Union{fmpz_mpoly, AbstractAlgebra.Generic.MPoly{nf_elem}, Missing}
end

# could redo this for classical to remove sign
# for a trellis with 10 million edges, this would save 10 MB
# to do this, make EdgeC and EdgeQ and change all following functions to Union both
mutable struct Edge
    label::Union{fq_nmod, fq_nmod_mat}
    weight::Float64
    out_vertex::Int64
    sign::Union{Missing, nmod}
end

mutable struct Trellis
    vertices::Vector{Vector{Vertex}}
    edges::Vector{Vector{Vector{Edge}}}
    code::T where T <: AbstractCode
    # complete weight enumerator
    CWE::Union{WeightEnumerator, Missing}
    # shifted::Bool
    shift::fq_nmod_mat
end

"""
    vertices(T::Trellis)

Return the set of vertices of the trellis `T`.
"""
vertices(T::Trellis) = T.vertices

"""
    edges(T::Trellis)

Return the set of edges of the trellis `T`.
"""
edges(T::Trellis) = T.edges

"""
    isshifted(T::Trellis)

Return `true` if the trellis now represents a shifted version of the original
code.
"""
is_shifted(T::Trellis) = !iszero(T.shift)

function ==(a::Vertex, b::Vertex)
    return a.label == b.label && a.prev == b.prev && a.next == b.next &&
        a.value == b.value && a.edge_loc == b.edge_loc && a.polynomial == b.polynomial
end

function ==(a::Vector{Vertex}, b::Vector{Vertex})
    if length(a) != length(b)
        return false
    end

    for i in 1:length(a)
        if a[i] != b[i]
            return false
        end
    end
    return true
end

function ==(a::Edge, b::Edge)
    return a.label == b.label && a.weight == b.weight && a.out_vertex == b.out_vertex
end

function ==(a::Vector{Edge}, b::Vector{Edge})
    if length(a) != length(b)
        return false
    end

    for i in 1:length(a)
        if a[i] != b[i]
            return false
        end
    end

    return true
end

function ==(a::Vector{Vector{Edge}}, b::Vector{Vector{Edge}})
    if length(a) != length(b)
        return false
    end

    for i in 1:length(a)
        if a[i] != b[i]
            return false
        end
    end

    return true
end

function ==(a::Vector{Vector{Vector{Edge}}}, b::Vector{Vector{Vector{Edge}}})
    if length(a) != length(b)
        return false
    end

    for i in 1:length(a)
        if a[i] != b[i]
            return false
        end
    end

    return true
end

function is_isomorphic(V1::Vector{Vector{Vertex}}, V2::Vector{Vector{Vertex}})
    if length(V1) != length(V2)
        return false
    else
        for i in 1:length(V1)
            if length(V1[i]) != length(V2[i])
                return false
            end
        end
        return true
    end
end

function is_isomorphic(E1::Vector{Vector{Vector{Edge}}}, E2::Vector{Vector{Vector{Edge}}})
    if length(E1) != length(E2)
        return false
    else
        for i in 1:length(E1)
            if length(E1[i]) != length(E2[i]) || length(E1[i][1]) != length(E2[i][1])
                return false
            end
        end
        return true
    end
end

function is_isomorphic(V1::Vector{Vector{Vertex}}, E1::Vector{Vector{Vector{Edge}}},
    V2::Vector{Vector{Vertex}}, E2::Vector{Vector{Vector{Edge}}})

    iso_V = is_isomorphic(V1, V2)
    if !iso_V
        return false
    end
    return is_isomorphic(E1, E2)
end

function is_isomorphic(V1::Vector{Vector{Vertex}}, V2::Vector{Vector{Vertex}},
    E1::Vector{Vector{Vector{Edge}}}, E2::Vector{Vector{Vector{Edge}}})

    return is_isomorphic(V1, E1, V2, E2)
end

# provide "sorted" option
function is_equal(V1::Vector{Vector{Vertex}}, V2::Vector{Vector{Vertex}}, verbose=false)
    is_iso = is_isomorphic(V1, V2)
    verbose && println("V1 and V2 are isomorphic.")

    if !is_iso
        verbose && println("V1 and V2 are not isomorphic.")
        return false
    end

    for (i, V1i) in enumerate(V1)
        verbose && println("Checking V1[", i, "]")

        for (j, v1) in enumerate(V1i)
            found_v = false
            for (k, v2) in enumerate(V2[i])
                if v1.label == v2.label
                    found_v = true
                    verbose && println("V1[", i, "][", j, "] = V2[", i, "][", k, "]")
                    break
                end
            end

            if !found_v
                verbose && println("Vertex V1[", i, "][", j, "] not found in V2[", i, "]")
                return false
            end
        end
    end

    return true
end

function ==(V1::Vector{Vector{Vertex}}, V2::Vector{Vector{Vertex}})
    return is_equal(V1, V2, false)
end

function is_equal(V1::Vector{Vector{Vertex}}, E1::Vector{Vector{Vector{Edge}}},
    V2::Vector{Vector{Vertex}}, E2::Vector{Vector{Vector{Edge}}},
    verbose=false, rhs_sorted=false)

    is_iso = is_isomorphic(V1, E1, V2, E2)
    verbose && println("V1 and V2 and E1 and E2 are isomorphic, respectively.")
    if !is_iso
        verbose && println("Not  isomorphic.")
        return false
    end

    for (i, V1i) in enumerate(V1)
        len_right = length(V1i)
        verbose && println("Checking V1[", i, "]")

        for (j, v1) in enumerate(V1i)
            found_v = false
            v_loc = 0
            if rhs_sorted
                bin_search_left = 1
                bin_search_right = len_right
                while bin_search_left <= bin_search_right
                    mid = fld((bin_search_left + bin_search_right), 2)
                    if V2[i][mid].label < v1.label
                        bin_search_left = mid + 1
                    elseif V2[i][mid].label > v1.label
                        bin_search_right = mid - 1
                    else
                        found_v = true
                        v_loc = mid
                        verbose && println("V1[", i, "][", j, "] = V2[", i, "][", mid, "]")
                        break
                    end
                end
            else
                for (k, v2) in enumerate(V2[i])
                    if v1.label == v2.label
                        found_v = true
                        v_loc = k
                        verbose && println("V1[", i, "][", j, "] = V2[", i, "][", k, "]")
                        break
                    end
                end
            end

            if !found_v
                verbose && println("Vertex V1[", i, "][", j, "] not found in V2[", i, "]")
                return false
            else
                if i != 1
                    verbose && println("Checking for equality in corresponding edges")
                    for (k, e1) in enumerate(E1[i - 1][j])
                        found_e = false
                        for (l, e2) in enumerate(E2[i - 1][v_loc])
                            if e1.label == e2.label && e1.weight == e2.weight && V1[i - 1][e1.out_vertex].label == V2[i - 1][e2.out_vertex].label
                                found_e = true
                                verbose && println("E1[", i - 1, "][", j, "][", k, "] = E2[", i - 1, "][", v_loc, "][", l, "]")
                                break
                            end
                        end

                        if !found_e
                            verbose && println("Edge E1[", i - 1, "][", j, "][", k, "] not found in E2[", i - 1, "][", v_loc, "]")
                            return false
                        end
                    end
                end
            end
        end
    end
    return true
end

function is_equal(V1::Vector{Vector{Vertex}}, V2::Vector{Vector{Vertex}},
    E1::Vector{Vector{Vector{Edge}}}, E2::Vector{Vector{Vector{Edge}}},
    verbose=false, rhs_sorted=false)

    return is_equal(V1, E1, V2, E2, verbose, rhs_sorted)
end

function _sort_by_left_index(A::fq_nmod_mat)
    nr, nc = size(A)
    arr = []
    for r in 1:nr
        for c in 1:nc
            if !iszero(A[r, c])
                push!(arr, [c, A[r, :]])
                break
            end
        end
        # push!(arr, [numcols, A[r, :]])
    end
    # println("sort1: $arr")
    sort!(arr, by = x -> x[1])
    # println("sort2: $arr")
    # vcat is faster and cheaper than the following
    # return A[[arr2[i][2] for i in 1:num_rows], :]
    return vcat([arr[i][2] for i in 1:nr]...)
end

function _sort_by_right_index(A::fq_nmod_mat)
    nr, nc = size(A)
    arr = []
    for r in 1:nr
        for c in nc:-1:1
            if !iszero(A[r, c])
                push!(arr, [c, A[r, :]])
                break
            end
        end
        # push!(arr, [1, A[r, :]])
    end
    sort!(arr, by = x -> x[1]) #, rev=true
    return vcat([arr[i][2] for i in 1:nr]...)
end

function _left_right_indices(A::fq_nmod_mat)
    # iseven(size(A, 2)) || error("Vectors in leftrightindices must have even length.")
    # n = div(size(A, 2), 2)

#take quadratic only here

    nr, nc = size(A)
    left = Vector{Int64}()
    right = Vector{Int64}()
    for r in 1:nr
        for c in 1:nc
            if !iszero(A[r, c])
                push!(left, c)
                break
            end
        end
    end

    for r in 1:nr
        for c in nc:-1:1
            if !iszero(A[r, c])
                push!(right, c)
                break
            end
        end
    end
    return left, right
end

function _find_active(A::fq_nmod_mat, edges::Bool=false)
    # need to make sure quadratic extension
    nr, nc = size(A)
    left_right = [[0, 0] for _ in 1:nr]
    for r in 1:nr
        for c in 1:nc
            if !iszero(A[r, c])
                left_right[r][c] = c
                break
            end
        end
    end

    for r in 1:nr
        for c in nc:-1:1
            if !iszero(A[r, c])
                left_right[i][2] = c
                break
            end
        end
    end
    # display(A)
    # println(left_right)

    active = Vector{Vector{Int64}}()
    if edges
        for c in 1:nc
            arr_c = Vector{Int64}()
            for r in 1:nr
                # c == 5 && println(c, ", ", r, ", ", left_right[r])
                if left_right[r][1] <= c <= left_right[r][2]
                    # c == 5 && println("added")
                    push!(arr_c, c)
                end
            end
            if !isempty(arr_c)
                push!(active, arr_c)
            end
        end
    else
        for c in 1:nc
            arr_i = Vector{Int64}()
            for r in 1:nr
                if left_right[r][1] <= c < left_right[r][2]
                    push!(arr_c, r)
                end
            end
            if !isempty(arr_c)
                push!(active, arr_c)
            end
        end
    end
    # println(active)
    return active
end

# finds all elements of B which have zero symplectic inner product with all of A
function _kernel_inner_prod(A::fq_nmod_mat, B::fq_nmod_mat,
    inner_prod::String="Euclidean", return_ker::Bool=false)

    inner_prod ∈ ["Euclidean", "symplectic"] || error("Unsupported inner product type in _kernel_inner_prod; expected: `Euclidean`, `symplectic`, received: $inner_prod.")
    size(A, 2) == size(B, 2) || error("Length mismatch in _kernel_inner_prod.")
    if inner_prod == "symplectic"
        iseven(size(A, 2)) || error("Vectors in symplectic_kernel must have even length.")
    end

    # unclear if it will compile away the choice if we put the check just around
    # the two lines of changed code inside the loops or if it will check it
    # every time
    # can check compiled code later
    nr_A, nc_A = size(A)
    half_nc_A = div(nc_A, 2)
    nr_B = nrows(B)
    if inner_prod == "symplectic"
        if return_ker
            ker = []
            for rb in 1:nr_B
                for ra in 1:nr_A
                    @views if !iszero(symplectic_inner_product(A[ra, :], B[rb, :]))
                        push!(ker, B[rb, :])
                        break
                    end
                end
            end
            return length(ker), reduce(vcat, ker)
        else
            A_Euc = hcat(A[:, half_nc_A + 1:end], -A[:, 1:half_nc_A])
            prod = A_Euc * transpose(B)
            iszero(prod) && return 0
            nc_prod = ncols(prod)
            count = 0
            for i in 1:nc_prod
                if !iszero(prod[:, i])
                    count += 1
                end
            end
            return count
        end
    else
        if return_ker
            ker = []
            for rb in 1:nr_B
                for ra in 1:nr_A
                    @views if !iszero(sum([A[ra, i] * B[rb, i] for i in 1:ncols(A)]))
                    # @views if !iszero(A[ra, :] ⋅ B[rb, :])
                        push!(ker, B[rb, :])
                        break
                    end
                end
            end
            return length(ker), reduce(vcat, ker)
        else
            prod = A * transpose(B)
            iszero(prod) && return 0
            nc_prod = ncols(prod)
            count = 0
            for i in 1:nc_prod
                if !iszero(prod[:, i])
                    count += 1
                end
            end
            return count
        end
    end
end

function _past_future(A::fq_nmod_mat)
    nr, nc = size(A)
    past = zeros(Int64, nc + 1)
    future = zeros(Int64, nc + 1)
    left, right = _left_right_indices(A)
    past[1] = 0
    future[1] = nr
    for i in 2:nc + 1
        past[i] = length(right[right .<= i - 1])
        future[i] = length(left[left .> i - 1])
    end
    return past, future
end

function load_balanced_code(profile::Vector{Int64})
    left_sum = 0
    left_loc = 1
    right_sum = 0
    right_loc = length(profile)
    while right_loc - left_loc > 2
        if left_sum <= right_sum
            left_sum += profile[left_loc]
            left_loc += 1
        else
            right_sum += profile[right_loc]
            right_loc -= 1
        end
    end
    # println(left_sum, ", ", right_sum)
    return left_loc, right_loc
end

function load_balanced_code(profiles::Vector{Vector{Int64}})
    length(profiles) == 4 || error("Expected a length 4 profile vector. Pass in all or just the edges.")
    return load_balanced_code(profiles[2])
end

"""
    trellis_oriented_form_linear(A::fq_nmod_mat)

Return the trellis oriented form of the matrix `A` assuming the row space is
linear.
"""
function trellis_oriented_form_linear(A::fq_nmod_mat)
    A = _sort_by_left_index(A)
    nc = ncols(A)
    for c in 1:nc
        left, right = _left_right_indices(A)
        rows = findall(x -> x == c, left)
        if length(rows) == 1
            @views A[rows[1], :] *= inv(A[rows[1], c])
        elseif length(rows) > 1
            # edges = []
            # for row in rows
            #     push!(edges, (row, A[row, c]))
            # end
            #
            # pivot = false
            # if !isempty(edges)
            #     pivot = true
            #     row, X = edges[1]
            #     @views A[row, :] *= inv(X)
            #     for i in 2:length(edges)
            #         @views A[edges[i][1], :] -= edges[i][2] * A[row, :]
            #     end
            # end

            #take first row, normalize to 1
            @views A[rows[1], :] *= inv(A[rows[1], c])
            # for rest of them, remove
            for i in 2:length(rows)
                @views A[rows[i], :] -= A[rows[i], c] * A[rows[1], :]
            end
        end
    end
    A = _sort_by_left_index(A)
    left, right = _left_right_indices(A)

    for c in nc:-1:1
        rows = findall(x -> x == c, right)
        if length(rows) == 1
            @views A[rows[1], :] *= inv(A[rows[1], c])
        elseif length(rows) > 1
            # edges = []
            # for row in rows
            #     push!(edges, (row, A[row, c]))
            # end
            #
            # pivot = false
            # if !isempty(edges)
            #     pivot = true
            #     row, X = edges[end]
            #     @views A[row, :] *= inv(X)
            #     # no problems here if this is only length 1
            #     for i in length(edges) - 1:-1:1
            #         @views A[edges[i][1], :] -= edges[i][2] * A[row, :]
            #     end
            # end

            #take last row, normalize to 1
            @views A[rows[end], :] *= inv(A[rows[end], c])
            # for rest of them, remove
            for i in length(rows) - 1:-1:1
                @views A[rows[i], :] -= A[rows[i], c] * A[rows[end], :]
            end
        end
        left, right = _left_right_indices(A)
    end
    return A
end

# seeks to eliminate edges of the form a + bω in a quadratic extension
# can't check that it's quadratic directly so note this will fail if a higher degree
# should work for degree 1 extensions given the way AbstractAlgebra stores coefficints?
"""
    trellis_oriented_form_additive(A::fq_nmod_mat)

Return the trellis oriented form of the matrix `A` assuming the row space is
additive.

# Note
* So far this is only implemented for quadratic extensions over a prime subfield, i.e., `GF(p^2)`.
"""
function trellis_oriented_form_additive(A::fq_nmod_mat)
    E = base_ring(A)
    degree(E) == 2 || error("So far this is only implemented for quadratic extensions over a prime subfield.")

    A = _sort_by_left_index(A)
    nc = ncols(A)
    # display(A)
    # println(" ")
    @views for c in 1:nc
        left, right = _left_right_indices(A)
        rows = findall(x -> x == c, left)
        if length(rows) == 1
            if !iszero(coeff(A[rows[1], c], 0)) && !iszero(coeff(A[rows[1], c], 1))
                A[rows[1], :] *= inv(E(coeff(A[rows[1], c], 0)))
            elseif !iszero(coeff(A[rows[1], c], 0))
                A[rows[1], :] *= inv(E(coeff(A[rows[1], c], 0)))
            else
                A[rows[1], :] *= inv(E(coeff(A[rows[1], c], 1)))
            end
        elseif length(rows) > 1
            X_edges = []
            Z_edges = []
            mixed_edges = []
            for row in rows
                if !iszero(coeff(A[row, c], 0)) && !iszero(coeff(A[row, c], 1))
                    push!(mixed_edges, (row, A[row, c]))
                elseif !iszero(coeff(A[row, c], 0)) && iszero(coeff(A[row, c], 1))
                    push!(X_edges, (row, A[row, c]))
                elseif iszero(coeff(A[row, c], 0)) && !iszero(coeff(A[row, c], 1))
                    push!(Z_edges, (row, A[row, c]))
                end
            end
            # println(X_edges)
            # println(Z_edges)
            # println(mixed_edges)

            if length(X_edges) <= 1 && length(Z_edges) <= 1 && length(mixed_edges) <= 1 && (length(X_edges) + length(Z_edges) + length(mixed_edges)) <= 2
                continue
            else
                X_pivot = false
                Z_pivot = false
                if !isempty(X_edges)
                    X_pivot = true
                    row, X = X_edges[1]
                    A[row, :] *= inv(E(coeff(X, 0)))
                    # println("X")
                    # display(A)
                    # println(" ")
                    # no problems here if this is only length 1
                    for i in 2:length(X_edges)
                        A[X_edges[i][1], :] -= E(coeff(X_edges[i][2], 0)) * A[row, :]
                    end
                end
                if !isempty(Z_edges)
                    Z_pivot = true
                    row, Z = Z_edges[1]
                    A[row, :] *= inv(E(coeff(Z, 1)))
                    # println("Z")
                    # display(A)
                    # println(" ")
                    # no problems here if this is only length 1
                    for i in 2:length(Z_edges)
                        A[Z_edges[i][1], :] -= E(coeff(Z_edges[i][2], 1)) * A[row, :]
                    end
                end
                if !isempty(mixed_edges)
                    if X_pivot && Z_pivot
                        for i in 1:length(mixed_edges)
                            A[mixed_edges[i][1], :] -= E(coeff(mixed_edges[i][2], 0)) * A[X_edges[1][1], :] + E(coeff(mixed_edges[i][2], 1)) * A[Z_edges[1][1], :]
                        end
                    elseif X_pivot
                        A[mixed_edges[1][1], :] -= E(coeff(mixed_edges[1][2], 0)) * A[X_edges[1][1], :]
                        # no problems here if this is only length 1
                        for i in 2:length(mixed_edges)
                            A[mixed_edges[i][1], :] -= E(coeff(mixed_edges[i][2], 0)) * A[X_edges[1][1], :] + E(coeff(mixed_edges[i][2], 1)) * A[mixed_edges[1][1], :]
                        end
                    elseif Z_pivot
                        A[mixed_edges[1][1], :] -= E(coeff(mixed_edges[1][2], 1)) * A[Z_edges[1][1], :]
                        # no problems here if this is only length 1
                        for i in 2:length(mixed_edges)
                            A[mixed_edges[i][1], :] -= E(coeff(mixed_edges[i][2], 0)) * A[mixed_edges[1][1], :] + E(coeff(mixed_edges[i][2], 1)) * A[Z_edges[1][1], :]
                        end
                    else
                        A[mixed_edges[1][1], :] *= inv(E(coeff(mixed_edges[1][2], 0)))
                        if length(mixed_edges) > 1
                            A[mixed_edges[2][1], :] -= E(coeff(mixed_edges[2][2], 0)) * A[mixed_edges[1][1], :]
                            A[mixed_edges[2][1], :] *= inv(E(coeff(A[mixed_edges[2][1], c], 1)))
                            if length(mixed_edges) > 2
                                A[mixed_edges[3][1], :] -= E(coeff(mixed_edges[3][2], 1)) * A[mixed_edges[2][1], :]
                                A[mixed_edges[3][1], :] *= inv(E(coeff(A[mixed_edges[3][1], c], 0)))
                                # no problems here if this is only length 3
                                for i in 3:length(mixed_edges)
                                    A[mixed_edges[i][1], :] -= E(coeff(mixed_edges[i][2], 0)) * A[mixed_edges[3][1], :] + E(coeff(mixed_edges[i][2], 1)) * A[mixed_edges[2][1], :]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    # println("after1")
    # display(A)
    # println(" ")
    A = _sort_by_left_index(A)
    left, right = _left_right_indices(A)
    # println("after2")
    # display(A)
    # println(" ")

    @views for c in nc:-1:1
        rows = findall(x -> x == c, right)
        if length(rows) == 1
            if !iszero(coeff(A[rows[1], c], 0)) && !iszero(coeff(A[rows[1], c], 1))
                A[rows[1], :] *= inv(E(coeff(A[rows[1], c], 0)))
            elseif !iszero(coeff(A[rows[1], c], 0))
                A[rows[1], :] *= inv(E(coeff(A[rows[1], c], 0)))
            else
                A[rows[1], :] *= inv(E(coeff(A[rows[1], c], 1)))
            end
        elseif length(rows) > 1
            X_edges = []
            Z_edges = []
            mixed_edges = []
            for row in rows
                if !iszero(coeff(A[row, c], 0)) && !iszero(coeff(A[row, c], 1))
                    push!(mixed_edges, (row, A[row, c]))
                elseif !iszero(coeff(A[row, c], 0)) && iszero(coeff(A[row, c], 1))
                    push!(X_edges, (row, A[row, c]))
                elseif iszero(coeff(A[row, c], 0)) && !iszero(coeff(A[row, c], 1))
                    push!(Z_edges, (row, A[row, c]))
                end
            end

            if length(X_edges) <= 1 && length(Z_edges) <= 1 && length(mixed_edges) <= 1 && (length(X_edges) + length(Z_edges) + length(mixed_edges)) <= 2
                continue
            else
                # need to determine if X and/or Z pivots are below Y and if not, use Y to
                # create one of these pivots and reset the X/Z pivots to a potentially lower row
                if !isempty(mixed_edges)
                    # can only set one coefficient of a + bω to 1 in additive, do a
                    row, Y = mixed_edges[end]
                    A[row, :] *= inv(E(coeff(Y, 0)))
                    for i in length(mixed_edges) - 1:-1:1
                        # can't use a + bω to eliminate a c + dω, so first use 1 + b'ω to eliminate c
                        A[mixed_edges[i][1], :] -= E(coeff(mixed_edges[i][2], 0)) * A[row, :]
                        # now turn d to 1 - these are all pure Z's now and some could be pivots
                        if !iszero(coeff(A[mixed_edges[i][1], c], 1))
                            A[mixed_edges[i][1], :] *= inv(E(coeff(A[mixed_edges[i][1], c], 1)))
                            append!(Z_edges, (mixed_edges[i][1], A[mixed_edges[i][1], c]))
                        end
                    end
                    mixed_edges = [mixed_edges[end]]
                    sort!(Z_edges, by=x->x[1])
                end

                if !isempty(X_edges)
                    row, X = X_edges[end]
                    A[row, :] *= inv(E(coeff(X, 0)))
                    # no problems here if this is only length 1
                    for i in length(X_edges) - 1:-1:1
                        A[X_edges[i][1], :] -= E(coeff(X_edges[i][2], 0)) * A[row, :]
                    end
                    X_edges = [X_edges[end]]
                end
                if !isempty(Z_edges)
                    row, Z = Z_edges[end]
                    A[row, :] *= inv(E(coeff(Z, 1)))
                    # no problems here if this is only length 1
                    for i in length(Z_edges) - 1:-1:1
                        A[Z_edges[i][1], :] -= E(coeff(Z_edges[i][2], 1)) * A[row, :]
                    end
                    Z_edges = [Z_edges[end]]
                end

                # now could have three pivots, remove the top one using the bottom two
                if !isempty(mixed_edges) && !isempty(X_edges) && !isempty(Z_edges)
                    _, loc = findmin([mixed_edges[1], X_edges[1], Z_edges[1]])
                    if loc == 1
                        # Y edge is top, apply both X and Z pivots to remove it
                        # pivot to be removed is of the form 1 + bω
                        A[mixed_edges[1][1], :] -= A[X_edges[1][1], :]
                        A[mixed_edges[1][1], :] -= E(coeff(A[mixed_edges[1][1], c], 1)) * A[Z_edges[1][1], :]
                    elseif loc == 2
                        # X edge is top, apply both Y and Z pivots to remove it
                        A[X_edges[1][1], :] -= A[mixed_edges[1][1], :]
                        # pivot to be remove is now of the form bω coming from the Y
                        A[X_edges[1][1], :] -= E(coeff(A[X_edges[1][1], c], 1)) * A[Z_edges[1][1], :]
                    else
                        # Z edge is top, apply both Y and X pivots to remove it
                        A[Z_edges[1][1], :] -= inv(E(coeff(mixed_edges[1][2], 1))) * A[mixed_edges[1][1], :]
                        # pivot to be remove is now of the form b^{-1}
                        A[Z_edges[1][1], :] -= E(coeff(A[Z_edges[1][1], c], 0)) * A[X_edges[1][1], :]
                    end
                end
            end
        end
        left, right = _left_right_indices(A)
    end
    return A
end

# need to brainstorem parallel edges
function trellis_profiles(wrt_V::fq_nmod_mat, wrt_E::fq_nmod_mat,
    boundaries::Union{Vector{Int64}, Missing}, inner_prod::String="Euclidean")

    inner_prod ∈ ["Euclidean", "symplectic"] || error("Unsupported inner product type in _trellisprofiles; expected: `Euclidean`, `symplectic`, received: $inner_prod.")

    if ismissing(boundaries)
        bds = [1:size(wrt_V, 2) + 1...]
    else
        bds = deepcopy(boundaries)
        if 0 ∈ bds
            bds .+= 1
        end
    end

    n = length(bds) - 1
    state_profile = zeros(Int64, n + 1)
    branch_profile = zeros(Int64, n)
    in_degrees = zeros(Int64, n)
    out_degrees = zeros(Int64, n)
    past, future = _past_future(wrt_V)
    past = past[bds]
    future = future[bds]

    # if inner_prod == "Euclidean"
    #     dim_ker = _kernel_inner_prod(wrt_V, wrt_E, inner_prod)
    # else
    #     dim_ker = _kernel_inner_prod(quadratic_to_symplectic(wrt_V),
    #         quadratic_to_symplectic(wrt_E), inner_prod)
    # end
    dim_ker = 0

    p = Int64(characteristic(base_ring(wrt_V)))
    for i in 1:n + 1
        state_profile[i] = p^(size(wrt_V, 1) - past[i] - future[i] - dim_ker)
    end

    left, right = _left_right_indices(wrt_E)
    past, future = _past_future(wrt_E)
    past = past[bds]
    future = future[bds]
    for i in 1:n
        dim_parallel = 0
        branch_profile[i] = p^(size(wrt_E, 1) - past[i] - future[i + 1] - dim_ker - dim_parallel)
        in_degrees[i] = div(branch_profile[i], state_profile[i + 1])
        out_degrees[i] = div(branch_profile[i], state_profile[i])
    end
    return [state_profile, branch_profile, in_degrees, out_degrees]
end

#############################
        # Classical
#############################

# TODO: handle lookup table better - temporarily skipping
# TODO: remove dictionaries, iterate once to find left, once for right
# keep first bipartite structure and shift it as a coset to find the next ones - fix for sectionalization
function syndrome_trellis(C::AbstractCode, type::String="primal", sect::Bool=true,
    verbose::Bool=false)

    (typeof(C) <: AbstractLinearCode || typeof(C) <: AbstractStabilizerCode) ||
        error("Syndrome trellises are so far only implemented for linear and stabilizer codes.")

    if typeof(C) <: AbstractLinearCode
        wrt_V = trellis_oriented_form_linear(parity_check_matrix(C))
        wrt_E = trellis_oriented_form_linear(generator_matrix(C))
        if sect
            boundaries, num_E_sect = optimal_sectionalization_C(wrt_V, wrt_E)
            profiles = trellis_profiles(wrt_V, wrt_E, boundaries, "Euclidean")
            profiles_no_sect = trellis_profiles(wrt_V, wrt_E, missing, "Euclidean")
            if verbose
                num_E_no_sect = sum(profiles_no_sect[2])
                println("|E| original: $num_E_no_sect, |E| sectionalized: $num_E_sect")
            end
            if length(boundaries) == 2
                if verbose
                    println("Sectionalized to lookup table, skipping")
                end
                boundaries = missing
                profiles = profiles_no_sect
            end
        else
            boundaries = missing
            profiles = trellis_profiles(wrt_V, wrt_E, missing, "Euclidean")
            if verbose
                boundaries2, num_E_sect = optimal_sectionalization_C(wrt_V, wrt_E)
                profiles_sect = trellis_profiles(wrt_V, wrt_E, boundaries2, "Euclidean")
                num_E_no_sect = sum(profiles[2])
                println("|E| original: $num_E_no_sect, |E| sectionalized: $num_E_sect")
            end
        end
    else
        if type == "primal"
            wrt_V = trellis_oriented_form_additive(stabilizers(C))
            # BUG: normalizermatrix function no longer exists
            wrt_E = trellis_oriented_form_additive(normalizermatrix(C))
        else
            wrt_V = trellis_oriented_form_additive(normalizermatrix(C))
            wrt_E = trellis_oriented_form_additive(stabilizers(C))
        end
        if sect
            boundaries, num_E_sect = optimal_sectionalization_Q(wrt_V, wrt_E)
            profiles = trellis_profiles(wrt_V, wrt_E, boundaries, "symplectic")
            profiles_no_sect = trellis_profiles(wrt_V, wrt_E, missing, "symplectic")
            if verbose
                num_E_no_sect = sum(profiles_no_sect[2])
                println("|E| original: $num_E_no_sect, |E| sectionalized: $num_E_sect")
            end
            if length(boundaries) == 2
                if verbose
                    println("Sectionalized to lookup table, skipping")
                end
                boundaries = missing
                profiles = profiles_no_sect
            end
        else
            boundaries = missing
            profiles = trellis_profiles(wrt_V, wrt_E, missing, "symplectic")
            if verbose
                boundaries2, num_E_sect = optimal_sectionalization_Q(wrt_V, wrt_E)
                profiles_sect = trellis_profiles(wrt_V, wrt_E, boundaries2, "symplectic")
                num_E_no_sect = sum(profiles[2])
                println("|E| original: $num_E_no_sect, |E| sectionalized: $num_E_sect")
            end
        end
    end

    if ismissing(boundaries)
        bds = [0:length(C)...]
    else
        bds = deepcopy(boundaries)
    end

    if typeof(C) <: AbstractLinearCode
        K = field(C)
    else
        # BUG: no longer exists
        K = quadratic_field(C)
        R = parent(C.char_vec[1])
    end
    p = Int64(characteristic(K))
    n = length(C)
    V = Vector{Vertex}[Vertex[] for _ in 1:length(bds)]
    Threads.@threads for i in 1:length(profiles[1])
        V[i] = [Vertex(999, 0, 0, 0.0, 0, missing) for _ in 1:profiles[1][i]]
    end
    V[1] = [Vertex(0, 0, 0, 0.0, 0, missing)]
    V[end] = [Vertex(0, 0, 0, 0.0, 0, missing)]
    verbose && println("Vertex preallocation completed.")

    E = Vector{Vector{Edge}}[[Edge[]] for _ in 1:length(profiles[3])]
    Threads.@threads for i in 1:length(profiles[3])
        # the j-th element of Ei is going to be all of the edges going into Vi[j]
        E[i] = [[Edge(K(0), 0.0, 0, missing) for j in 1:profiles[3][i]] for _ in 1:profiles[1][i + 1]]
    end
    verbose && println("Edge preallocation completed.")

    biz = BigInt(0)
    bio = BigInt(1)
    syn_len = size(wrt_V, 1)
    active = _find_active(wrt_V)
    active = active[bds[2:end - 1]]
    Threads.@threads for i in 2:length(bds) - 1
        Vi_size = profiles[1][i]
        for num in 0:Vi_size - 1
            bin = reverse(digits(num, base=2, pad=length(active[i - 1])))
            temp_label = zeros(Int64, syn_len)
            loc = 1
            for j in active[i - 1]
                temp_label[j] = bin[loc]
                loc += 1
            end

            cur_label = biz
            for (shift, val) in enumerate(reverse(temp_label, dims=1))
                if val == 1
                    cur_label += bio << (shift - 1)
                end
            end
            V[i][num + 1].label = cur_label
        end
    end
    verbose && println("Vertex construction completed.")

    # there has to be a one-liner for the below
    left, right = _left_right_indices(wrt_E)
    active = _find_active(wrt_E, true)
    if ismissing(boundaries)
        active_temp = active
    else
        active_temp = Vector{Vector{Int64}}()
        for i in 1:length(bds) - 1
            # temp = Vector{Int64}()
            # for j in bds[i] + 1:bds[i + 1]
            #     append!(temp, active[j])
            # end
            # push!(active_temp, sort!(unique!(temp)))
            push!(active_temp, sort!(unique!(vcat([active[j] for j in bds[i] + 1:bds[i + 1]]...))))
        end
    end

    if typeof(C) <: AbstractLinearCode
        H = FpmattoJulia(wrt_V)
    else
        sym_wrt_V = quadratic_to_symplectic(wrt_V)
        H = FpmattoJulia(hcat(sym_wrt_V[:, n + 1:end], -sym_wrt_V[:, 1:n]))
    end
    # Threads.@threads 
    for i in length(bds) - 1:-1:1
        verbose && println("Starting E[$i]")
        # seclen = bds[i + 1] - bds[i]
        valid_edges = Vector{fq_nmod_mat}()
        edge_contrib = Dict{fq_nmod_mat, Vector{Int64}}()
        contrib_edge = Dict{Vector{Int64}, fq_nmod_mat}()

        for a in active_temp[i]
            temp = wrt_E[a, bds[i] + 1:bds[i + 1]]
            if !iszero(temp)
                push!(valid_edges, temp)
            end
        end
        unique!(valid_edges)
        # i == 1 && display(valid_edges)

        # may need to use Oscar here now that the import is switched
        for iter in Nemo.AbstractAlgebra.ProductIterator(collect(0:p - 1), length(valid_edges))
            e = K(iter[1]) * valid_edges[1]
            for r in 2:length(valid_edges)
                if !iszero(iter[r])
                    e += K(iter[r]) * valid_edges[r]
                end
            end

            if typeof(C) <: AbstractLinearCode
                P = zeros(Int64, n)
                for (j, k) in enumerate(e)
                    P[bds[i] + j] = coeff(k, 0)
                end
            else
                P = zeros(Int64, 2 * n)
                for (j, k) in enumerate(e)
                    P[bds[i] + j] = coeff(k, 0)
                    P[bds[i] + j + n] = coeff(k, 1)
                end
            end
            syn = H * P .% p
            edge_contrib[e] = syn
            contrib_edge[syn] = e
            # when it's so small it sectionalizes to a lookup table, there is
            # only one syndrome = 0 at V_right = V_n, so contrib_edge only has
            # one value. to fix, can only handle this case entirely manually
        end
        verbose && println("Edges dictionaries completed for E[$i].")
        # i == 1 && display(edge_contrib)
        # display(edge_contrib)
        # display(contrib_edge)

        V_left = V[i]
        V_right = V[i + 1]
        len_left = length(V_left)
        len_right = length(V_right)
        V_right_locs = trues(len_right)
        starting_right_index = 1
        # keep below here instead of preallocating above or else the answer comes out wrong
        blank = K(0)
        # println(len_right)

        while starting_right_index <= len_right
            starting_right_v = V_right[starting_right_index].label
            left_vertices = Vector{Tuple{Int64, BigInt, Vector{Int64}}}()
            right_vertices = Vector{Tuple{Int64, BigInt, Vector{Int64}}}()
            sizehint!(left_vertices, profiles[4][i])
            sizehint!(right_vertices, profiles[3][i])

            starting_right_v_syn = reverse(digits(starting_right_v, base=2, pad=syn_len))
            # println("i: $i, syn_len: $syn_len, srsyn: $starting_right_v_syn, srv: $starting_right_v")
            push!(right_vertices, (starting_right_index, starting_right_v, starting_right_v_syn))
            connecting_starts = blank
            starting_left_v = biz

            # start with a fixed right vertex and find all left vertices
            for lab in keys(edge_contrib)
                temp = (starting_right_v_syn .- edge_contrib[lab])
                for t in 1:length(temp)
                    if temp[t] < 0
                        temp[t] = p + temp[t]
                    end
                end
                temp = temp .% p

                left_label = biz
                for (shift, val) in enumerate(reverse(temp, dims=1))
                    if val == 1
                        left_label += bio << (shift - 1)
                    end
                end

                bin_search_left = 1
                bin_search_right = len_left
                while bin_search_left <= bin_search_right
                    mid = fld((bin_search_left + bin_search_right), 2)
                    if V_left[mid].label < left_label
                        bin_search_left = mid + 1
                    elseif V_left[mid].label > left_label
                        bin_search_right = mid - 1
                    else
                        push!(left_vertices, (mid, left_label, temp))
                        if connecting_starts == blank
                            connecting_starts = lab
                            starting_left_v = temp
                        end
                        break
                    end
                end

                if length(left_vertices) == profiles[3][i]
                    break
                end
            end

            # start with first left vertex and find all right vertices
            if length(right_vertices) != profiles[4][i]
                for lab in keys(edge_contrib)
                    if lab != connecting_starts
                        temp = (starting_left_v .+ edge_contrib[lab]) .% p
                        right_label = biz
                        for (shift, val) in enumerate(reverse(temp, dims=1))
                            if val == 1
                                right_label += bio << (shift - 1)
                            end
                        end

                        bin_search_left = 1
                        bin_search_right = len_right
                        while bin_search_left <= bin_search_right
                            mid = fld((bin_search_left + bin_search_right), 2)
                            if V_right[mid].label < right_label
                                bin_search_left = mid + 1
                            elseif V_right[mid].label > right_label
                                bin_search_right = mid - 1
                            else
                                push!(right_vertices, (mid, right_label, temp))
                                break
                            end
                        end
                    end

                    if length(right_vertices) == profiles[4][i]
                        break
                    end
                end
            end

            # can probably skip this recalculation of temp by immediately storing
            # instead of building right and left vertex lists
            # should now have all vertices
            for (right_index, _, right_syn) in right_vertices
                count = 1
                for (left_index, _, left_syn) in left_vertices
                    temp = right_syn .- left_syn
                    for t in 1:length(temp)
                        if temp[t] < 0
                            temp[t] = p + temp[t]
                        end
                    end
                    temp = temp .% p
                    lab = contrib_edge[temp]

                    E[i][right_index][count].label = lab
                    E[i][right_index][count].out_vertex = left_index
                    if typeof(C) <: AbstractStabilizerCode
                        sign = R(0)
                        for (j, k) in enumerate(lab)
                            if !iszero(coeff(k, 0))
                                sign += character_vector(C)[bds[i] + j]
                            end
                            if !iszero(coeff(k, 1))
                                sign += character_vector(C)[bds[i] + j + n]
                            end
                        end
                        E[i][right_index][count].sign = sign
                    end

                    count += 1
                end
                V_right_locs[right_index] = false
            end

            # should have ==
            while starting_right_index <= len_right
                starting_right_index += 1
                if starting_right_index <= len_right && !V_right_locs[starting_right_index]
                    starting_right_index += 1
                else
                    break
                end
            end

        end
        verbose && println("E[$i] complete")
    end

    if typeof(C) <: AbstractLinearCode
        return Trellis(V, E, C, missing, zero_matrix(K, 1, n))
    else
        return Trellis(V, E, C, missing, zero_matrix(K, 1, 2 * n))
    end
end

# should probably return missing to unify if statements in the function below
function trellis_profiles(C::AbstractLinearCode, sect::Bool=false)
    G_TOF = trellis_oriented_form_linear(generator_matrix(C))
    H_TOF = trellis_oriented_form_linear(parity_check_matrix(C))
    if sect
        opt, _ = optimal_sectionalization_C(H_TOF, G_TOF)
        return trellis_profiles(H_TOF, G_TOF, opt, "Euclidean"), opt
    end
    return trellis_profiles(H_TOF, G_TOF, missing, "Euclidean")
end

# # recopy above so I don't call it and redo the TOF multiple times
# function syndrome_trellis(C::AbstractLinearCode, sect::Bool=false)
#     G_TOF = trellis_oriented_form_linear(generator_matrix(C))
#     H_TOF = trellis_oriented_form_linear(parity_check_matrix(C))
#     if sect
#         opt, _ = optimal_sectionalization_C(H_TOF, G_TOF)
#         profiles = trellis_profiles(H_TOF, G_TOF, opt, "Euclidean")
#         return _syndrometrellisC(profiles, opt, H_TOF, G_TOF, false)
#     else
#         profiles = trellis_profiles(H_TOF, G_TOF, missing, "Euclidean")
#         return _syndrometrellisC(profiles, missing, H_TOF, G_TOF, false)
#         # return _syndrometrellisC(profiles, missing, H_TOF, G_TOF, [1 for _ in 1:length(C)], ' ', false)
#     end
# end

#############################
         # Quantum
#############################

# only valid for quantum codes
function optimal_sectionalization_Q(wrt_V::fq_nmod_mat, wrt_E::fq_nmod_mat)
    K = base_ring(wrt_E)
    base_ring(wrt_V) == K || error("Vertices and edges must have the same base ring.")

    n = size(wrt_V, 2)
    p = Int64(characteristic(K))
    V = [Vertex(i, 0, 0, 0.0, 0, missing) for i in 1:n + 1]
    E = [[Edge(K(0), 0.0, 0, missing) for i in 1:j] for j in n:-1:1]

    sym_V = quadratic_to_symplectic(wrt_V)
    sym_E = quadratic_to_symplectic(wrt_E)
    # dim_ker = _kernel_inner_prod(sym_V, sym_E, "symplectic")
    dim_ker = 0
    # println(dim_ker)
    past, future = _past_future(wrt_E)
    for i in 1:n
        for j in i:n
            # arbitrary size cutoff
            if size(wrt_E, 1) - past[i] - future[j + 1] - dim_ker > 50
                E[i][j - i + 1].weight = Inf
            else
                E[i][j - i + 1].weight = p^(size(wrt_E, 1) - past[i] - future[j + 1] - dim_ker)
            end
        end
    end

    for i in 1:n
        left_b = 1
        right_b = i
        arr = [E[left_b][right_b].weight + V[left_b].value]

        while left_b < i
            left_b += 1
            right_b -= 1
            append!(arr, E[left_b][right_b].weight + V[left_b].value)
        end

        m, arg = findmin(arr)
        V[i + 1].prev = arg
        V[i + 1].value = m
    end

    sect_boundaries = [n]
    next = V[end].prev
    # val = V[end].value
    while next > 0
        append!(sect_boundaries, next - 1)
        # val += V[next].value
        next = V[next].prev
    end
    return reverse(sect_boundaries), Int(V[end].value) #Int(val)
end

# TODO: remove dictionaries, iterate once to find left, once for right
# keep first bipartite structure and shift it as a coset to find the next ones - fix for sectionalization
# function _syndrometrellisQ(profiles::Vector{Vector{T}}, boundaries::Union{Vector{Int64}, Missing},
#     wrt_V::fq_nmod_mat, wrt_E::fq_nmod_mat, char_vec::Vector{Int64}, Pauli::Char=' ',
#     verbose=false) where T <: Integer
#
#     Pauli ∈ [' ', 'X', 'Z'] || error("Pauli parameter needs to be ' ', 'X', or 'Z'; received $Pauli.")
#     K = base_ring(wrt_E)
#     base_ring(wrt_V) == K || error("Vertices and edges must have the same base ring.")
#     if ismissing(boundaries)
#         bds = [0:size(wrt_V, 2)...]
#     else
#         bds = deepcopy(boundaries)
#     end
#
#     ω = gen(K)
#     p = Int64(characteristic(K))
#     n = length(profiles[1]) - 1
#     symsize = size(wrt_V, 2)
#     V = Vector{Vertex}[Vertex[] for i = 1:n + 1]
#     Threads.@threads for i = 1:n + 1
#         V[i] = [Vertex(999, 0, 0, 0.0, 0, missing) for j = 1:profiles[1][i]]
#     end
#     V[1] = [Vertex(0, 0, 0, 0.0, 0, missing)]
#     V[end] = [Vertex(0, 0, 0, 0.0, 0, missing)]
#     verbose && println("Vertex preallocation completed.")
#
#     E = Vector{Vector{Edge}}[[Edge[]] for i = 1:n]
#     Threads.@threads for i in 1:n
#         # the j-th element of Ei is going to be all of the edges going into Vi[j]
#         E[i] = [[Edge(K(0), 0.0, 0, 1) for j = 1:profiles[3][i]] for k = 1:profiles[1][i + 1]]
#     end
#     verbose && println("Edge preallocation completed.")
#
#     biz = BigInt(0)
#     bio = BigInt(1)
#     syn_len = size(wrt_V, 1)
#     active = _find_active(wrt_V)
#     active = active[bds[2:end - 1]]
#     Threads.@threads for i = 2:n
#         Vi_size = profiles[1][i]
#         for num in 0:Vi_size - 1
#             bin = reverse(digits(num, base=2, pad=length(active[i - 1])))
#             temp_label = zeros(Int64, syn_len)
#             loc = 1
#             for j in active[i - 1]
#                 temp_label[j] = bin[loc]
#                 loc += 1
#             end
#
#             cur_label = biz
#             for (shift, val) in enumerate(reverse(temp_label, dims=1))
#                 if val == 1
#                     cur_label += bio << (shift - 1)
#                 end
#             end
#             V[i][num + 1].label = cur_label
#         end
#     end
#     verbose && println("Vertex construction completed.")
#
#     # there has to be a one-liner for the below
#     left, right = _left_right_indices(wrt_E)
#     active = _find_active(wrt_E, true)
#     if ismissing(boundaries)
#         active_temp = active
#     else
#         active_temp = Vector{Vector{Int64}}()
#         for i in 1:length(bds) - 1
#             temp = Vector{Int64}()
#             for j in bds[i] + 1:bds[i + 1]
#                 if !(bds[i] <= left[j] && right[j] <= bds[i + 1])
#                     append!(temp, active[j])
#                 end
#             end
#             push!(active_temp, sort!(unique!(temp)))
#         end
#     end
#
#     sym_wrt_V = quadratic_to_symplectic(wrt_V)
#     G = FpmattoJulia(hcat(sym_wrt_V[:, symsize + 1:end], -sym_wrt_V[:, 1:symsize]))
#     Threads.@threads for i = n:-1:1
#         verbose && println("Starting E[$i]")
#         seclen = bds[i + 1] - bds[i]
#         valid_edges = Vector{fq_nmod_mat}()
#         edge_contrib = Dict{fq_nmod_mat, Vector{Int64}}()
#         contrib_edge = Dict{Vector{Int64}, fq_nmod_mat}()
#
#         for a in active_temp[i]
#             temp = wrt_E[a, bds[i] + 1:bds[i + 1]]
#             if !iszero(temp)
#                 push!(valid_edges, temp)
#             end
#         end
#         unique!(valid_edges)
#
#         for iter in Nemo.AbstractAlgebra.ProductIterator(collect(0:p - 1), length(valid_edges))
#             e = K(iter[1]) * valid_edges[1]
#             for r in 2:length(valid_edges)
#                 if !iszero(iter[r])
#                     e += K(iter[r]) * valid_edges[r]
#                 end
#             end
#
#             P = zeros(Int64, 2 * symsize)
#             for (j, k) in enumerate(e)
#                 P[bds[i] + j] = coeff(k, 0)
#                 P[bds[i] + j + symsize] = coeff(k, 1)
#             end
#             syn = G * P .% p
#             edge_contrib[e] = syn
#             contrib_edge[syn] = e
#         end
#         verbose && println("Edges dictionaries completed for E[$i].")
#
#         V_left = V[i]
#         V_right = V[i + 1]
#         len_left = length(V_left)
#         len_right = length(V_right)
#         V_right_locs = trues(len_right)
#         starting_right_index = 1
#         # keep below here instead of preallocating above or else the answer comes out wrong
#         blank = K(0)
#
#         while starting_right_index <= len_right
#             starting_right_v = V_right[starting_right_index].label
#             left_vertices = Vector{Tuple{Int64, BigInt, Vector{Int64}}}()
#             right_vertices = Vector{Tuple{Int64, BigInt, Vector{Int64}}}()
#             sizehint!(left_vertices, profiles[4][i])
#             sizehint!(right_vertices, profiles[3][i])
#
#             starting_right_v_syn = reverse(digits(starting_right_v, base=2, pad=syn_len))
#             push!(right_vertices, (starting_right_index, starting_right_v, starting_right_v_syn))
#             connecting_starts = blank
#             starting_left_v = biz
#
#             # start with a fixed right vertex and find all left vertices
#             for lab in keys(edge_contrib)
#                 temp = (starting_right_v_syn .- edge_contrib[lab])
#                 for t in 1:length(temp)
#                     if temp[t] < 0
#                         temp[t] = p + temp[t]
#                     end
#                 end
#                 temp = temp .% p
#
#                 left_label = biz
#                 for (shift, val) in enumerate(reverse(temp, dims=1))
#                     if val == 1
#                         left_label += bio << (shift - 1)
#                     end
#                 end
#
#                 bin_search_left = 1
#                 bin_search_right = len_left
#                 while bin_search_left <= bin_search_right
#                     mid = fld((bin_search_left + bin_search_right), 2)
#                     if V_left[mid].label < left_label
#                         bin_search_left = mid + 1
#                     elseif V_left[mid].label > left_label
#                         bin_search_right = mid - 1
#                     else
#                         push!(left_vertices, (mid, left_label, temp))
#                         if connecting_starts == blank
#                             connecting_starts = lab
#                             starting_left_v = temp
#                         end
#                         break
#                     end
#                 end
#
#                 if length(left_vertices) == profiles[3][i]
#                     break
#                 end
#             end
#
#             # start with first left vertex and find all right vertices
#             if length(right_vertices) != profiles[4][i]
#                 for lab in keys(edge_contrib)
#                     if lab != connecting_starts
#                         temp = (starting_left_v .+ edge_contrib[lab]) .% p
#                         right_label = biz
#                         for (shift, val) in enumerate(reverse(temp, dims=1))
#                             if val == 1
#                                 right_label += bio << (shift - 1)
#                             end
#                         end
#
#                         bin_search_left = 1
#                         bin_search_right = len_right
#                         while bin_search_left <= bin_search_right
#                             mid = fld((bin_search_left + bin_search_right), 2)
#                             if V_right[mid].label < right_label
#                                 bin_search_left = mid + 1
#                             elseif V_right[mid].label > right_label
#                                 bin_search_right = mid - 1
#                             else
#                                 push!(right_vertices, (mid, right_label, temp))
#                                 break
#                             end
#                         end
#                     end
#
#                     if length(right_vertices) == profiles[4][i]
#                         break
#                     end
#                 end
#             end
#
#             # can probably skip this recalculation of temp by immediately storing
#             # instead of building right and left vertex lists
#             # should now have all vertices
#             for (right_index, right_label, right_syn) in right_vertices
#                 count = 1
#                 for (left_index, left_label, left_syn) in left_vertices
#                     temp = right_syn .- left_syn
#                     for t in 1:length(temp)
#                         if temp[t] < 0
#                             temp[t] = p + temp[t]
#                         end
#                     end
#                     temp = temp .% p
#                     lab = contrib_edge[temp]
#                     sign = 1 # should be K(1) when implementing as roots of unity or in \C?
#                     for (j, k) in enumerate(lab)
#                         if !iszero(coeff(k, 0))
#                             sign *= char_vec[bds[i] + j]
#                         end
#                         if !iszero(coeff(k, 1))
#                             sign *= char_vec[bds[i] + j + symsize]
#                         end
#                     end
#
#                     E[i][right_index][count].label = lab
#                     E[i][right_index][count].out_vertex = left_index
#                     E[i][right_index][count].sign = sign
#                     count += 1
#                 end
#                 V_right_locs[right_index] = false
#             end
#
#             # should have ==
#             while starting_right_index <= len_right
#                 starting_right_index += 1
#                 if starting_right_index <= len_right && !V_right_locs[starting_right_index]
#                     starting_right_index += 1
#                 else
#                     break
#                 end
#             end
#
#         end
#         verbose && println("E[$i] complete")
#     end
#     return Trellis(V, E)
# end

# error models need to take CSS combinations into account
# Pauli == 'X'
# I -> I + X
# Z -> Z + Y
# Pauli == 'Z'
# I -> I + Z
# X -> X + Y
function weight_Q!(T::Trellis, Ps::fq_nmod_mat, err_models::Vector{Dict{String, Float64}},
    weight_type::String="additive")

    weight_type ∈ ["additive", "multiplicative"] || error("Weight type needs to be 'additive' or 'multiplicative'.")

    V = vertices(T)
    E = edges(T)
    for i in 1:length(E)
        model = err_models[i]
        for j in 1:length(V[i + 1])
            Threads.@threads for e in E[i][j]
                if weight_type == "additive"
                    weight = 0.0
                    for k in 1:length(e.label)
                        weight += model[e.label[k]]
                    end
                else
                    weight = 1.0
                    for k in 1:length(e.label)
                        weight *= model[e.label[k]]
                    end
                end
                e.weight = weight
            end
        end
    end
end

# error models need to take CSS combinations into account
# Pauli == 'X'
# I -> I + X
# Z -> Z + Y
# Pauli == 'Z'
# I -> I + Z
# X -> X + Y
# remove wrt_V here
function shift_and_weight_Q!(T::Trellis, Ps::fq_nmod_mat, boundaries::Union{Vector{Int64}, Missing},
    err_models::Vector{Dict{String, Float64}}, char_vec::Vector{Int64}, Pauli::Char=' ',
    weight_type::String="additive")

    Pauli ∈ [' ', 'X', 'Z'] || error("Pauli parameter needs to be ' ', 'X', or 'Z'; received $Pauli.")
    weight_type ∈ ["additive", "multiplicative"] || error("Weight type needs to be 'additive' or 'multiplicative'.")
    length(char_vec) == 2 * length(err_models) || error("Lengths of character vector and error models are not consistent.")

    K = base_ring(Ps)
    V = vertices(T)
    E = edges(T)
    code_n = length(err_models)
    if ismissing(boundaries)
        bds = [0:code_n...]
    else
        bds = deepcopy(boundaries)
    end

    for i in 1:length(E)
        model = err_models[i]
        for j in 1:length(V[i + 1])
            Threads.@threads for e in E[i][j]
                e.label += Ps[bds[i] + 1:bds[i + 1]]
                if Pauli == 'X'
                    for k in e.label
                        k -= K(coeff(k, 0))
                    end
                elseif Pauli == 'Z'
                    for k in e.label
                        k -= K(coeff(k, 1))
                    end
                end

                sign = 1 # should be K(1) when implementing as roots of unity or in \C?
                if weight_type == "additive"
                    weight = 0.0
                else
                    weight = 1.0
                end
                for (l, k) in enumerate(e.label)
                    if !iszero(coeff(k, 0))
                        sign *= char_vec[bds[i] + l]
                    end
                    if !iszero(coeff(k, 1))
                        sign *= char_vec[bds[i] + l + code_n]
                    end
                    if weight_type == "additive"
                        weight += model[k]
                    else
                        weight *= model[k]
                    end
                end
                e.sign = sign
                e.weight = weight
            end
        end
    end
end

# do I actually care about updating the signs here?
function shift_and_decode_Q!!(T::Trellis, Ps::fq_nmod_mat, boundaries::Union{Vector{Int64}, Missing},
    err_models::Vector{Dict{fq_nmod, Float64}}, char_vec::Vector{Int64}, Pauli::Char=' ',
    weight_type::String="additive")

    # Pauli ∈ [' ', 'X', 'Z'] || error("Pauli parameter needs to be ' ', 'X', or 'Z'; received $Pauli.")
    # weight_type ∈ ["additive", "multiplicative"] || error("Weight type needs to be 'additive' or 'multiplicative'.")
    # length(char_vec) == 2 * length(err_models) || error("Lengths of character vector and error models are not consistent.")

    K = base_ring(Ps)
    V = vertices(T)
    E = edges(T)
    code_n = length(err_models)
    if ismissing(boundaries)
        bds = [0:code_n...]
    else
        bds = deepcopy(boundaries)
    end

    for i in 1:length(E)
        model = err_models[i]
        for (j, v) in enumerate(V[i + 1])
            # don't Threads.@threads the below loop or you get a >100x slow down due to locking
            for e in E[i][j]
                e.label += Ps[1, bds[i] + 1:bds[i + 1]]
                if Pauli == 'X'
                    for k in e.label
                        k -= K(coeff(k, 0))
                    end
                elseif Pauli == 'Z'
                    for k in e.label
                        k -= K(coeff(k, 1))
                    end
                end

                # sign = 1 # should be K(1) when implementing as roots of unity or in \C?
                if weight_type == "additive"
                    weight = 0.0
                else
                    weight = 1.0
                end
                for (j, k) in enumerate(e.label)
                    # if !iszero(coeff(k, 0))
                    #     sign *= char_vec[bds[i] + j]
                    # end
                    # if !iszero(coeff(k, 1))
                    #     sign *= char_vec[bds[i] + j + code_n]
                    # end
                    if weight_type == "additive"
                        weight += model[k]
                    else
                        weight *= model[k]
                    end
                end
                # e.sign = sign
                e.weight = weight
            end

            # ignoring random tie breaker, breaks ties here by order of entered into trellis
            _, loc = findmin([e.weight + V[i][e.out_vertex].value for e in E[i][j]])
            v.value = E[i][j][loc].weight + V[i][E[i][j][loc].out_vertex].value
            v.prev = E[i][j][loc].out_vertex
            v.edge_loc = loc
        end
    end

    # redo here for the two-sided approach
    path = zero_matrix(base_ring(Ps), 1, code_n)
    curr = code_n
    prev = 1
    for i in length(E) + 1:-1:2
        e_label = E[i - 1][prev][V[i][prev].edge_loc].label
        path[1, (curr - length(e_label) + 1):curr] = e_label
        prev = V[i][prev].prev
        curr -= length(e_label)
    end
    return path
end

function shift!(T::Trellis, Ps::fq_nmod_mat, boundaries::Union{Vector{Int64}, Missing},
    err_models::Vector{Dict{fq_nmod, Float64}}, char_vec::Vector{Int64}, Pauli::Char=' ',
    weight_type::String="additive")

    # Pauli ∈ [' ', 'X', 'Z'] || error("Pauli parameter needs to be ' ', 'X', or 'Z'; received $Pauli.")
    # weight_type ∈ ["additive", "multiplicative"] || error("Weight type needs to be 'additive' or 'multiplicative'.")
    # length(char_vec) == 2 * length(err_models) || error("Lengths of character vector and error models are not consistent.")

    K = base_ring(Ps)
    V = vertices(T)
    E = edges(T)
    code_n = length(err_models)
    if ismissing(boundaries)
        bds = [0:code_n...]
    else
        bds = deepcopy(boundaries)
    end

    for i in 1:length(E)
        model = err_models[i]
        for j in 1:length(V[i + 1])
            Threads.@threads for e in E[i][j]
                e.label += Ps[1, bds[i] + 1:bds[i + 1]]
                if Pauli == 'X'
                    for k in e.label
                        k -= K(coeff(k, 0))
                    end
                elseif Pauli == 'Z'
                    for k in e.label
                        k -= K(coeff(k, 1))
                    end
                end

                # sign = 1 # should be K(1) when implementing as roots of unity or in \C?
                if weight_type == "additive"
                    weight = 0.0
                else
                    weight = 1.0
                end
                for (j, k) in enumerate(e.label)
                    # if !iszero(coeff(k, 0))
                    #     sign *= char_vec[bds[i] + j]
                    # end
                    # if !iszero(coeff(k, 1))
                    #     sign *= char_vec[bds[i] + j + code_n]
                    # end
                    if weight_type == "additive"
                        weight += model[k]
                    else
                        weight *= model[k]
                    end
                end
                # e.sign = sign
                e.weight = weight
            end
        end
    end
end

# think of more scenarios
# could allow general trellises given partial stabilizers for use in trellis product
function trellis_profiles(Q::AbstractStabilizerCode, type::String="weight", Pauli::Char=' ',
    sect::Bool=false)

    type ∈ ["weight", "decoding"] || error("Unknown type parameter in trellis_profiles.")
    # (Pauli != ' ' && typeof(Q) <: CSSCode) && error("Pauli parameter is non-empty but the code is not CSS.")
    Pauli ∈ [' ', 'X', 'Z'] || error("Unknown Pauli parameter $Pauli; must be ' ', 'X', or 'Z'.")

    if type == "weight"
        if Pauli == ' '
            S_TOF = trellis_oriented_form_additive(stabilizers(Q))
            # BUG: normalizermatrix no longer exists
            n_TOF = trellis_oriented_form_additive(normalizermatrix(Q))
            if sect
                opt, _ = optimal_sectionalization_Q(n_TOF, S_TOF)
                return trellis_profiles(n_TOF, S_TOF, opt, "symplectic"), opt
            end
            return trellis_profiles(n_TOF, S_TOF, missing, "symplectic")
        elseif Pauli == 'X'
            _, _, Z_perp, _, _, _ = split_symplectic_stabilizers(quadratic_to_symplectic(normalizermatrix(Q)), ones(Int64, size(normalizermatrix(Q), 1)))
            X = hcat(X_stabilizers(Q), zero_matrix(field(Q), size(X_stabilizers(Q), 1), size(X_stabilizers(Q), 2)))
            Z_perp_TOF = trellis_oriented_form_additive(symplectic_to_quadratic(Z_perp))
            X_TOF = trellis_oriented_form_additive(symplectic_to_quadratic(X))
            if sect
                opt, _ = optimal_sectionalization_Q(Z_perp_TOF, X_TOF)
                return trellis_profiles(Z_perp_TOF, X_TOF, opt, "symplectic"), opt
            end
            return trellis_profiles(Z_perp_TOF, X_TOF, missing, "symplectic")
        else
            X_perp, _, _, _, _, _ = split_symplectic_stabilizers(quadratic_to_symplectic(normalizermatrix(Q)), ones(Int64, size(normalizermatrix(Q), 1)))
            Z = hcat(zero_matrix(field(Q), size(Z_stabilizers(Q), 1), size(Z_stabilizers(Q), 2)), Z_stabilizers(Q))
            X_perp_TOF = trellis_oriented_form_additive(symplectic_to_quadratic(X_perp))
            Z_TOF = trellis_oriented_form_additive(symplectic_to_quadratic(Z))
            if sect
                opt, _ = optimal_sectionalization_Q(X_perp_TOF, Z_TOF)
                return trellis_profiles(X_perp_TOF, Z_TOF, opt, "symplectic"), opt
            end
            return trellis_profiles(X_perp_TOF, Z_TOF, missing, "symplectic")
        end
    else
        if Pauli == ' '
            S_TOF = trellis_oriented_form_additive(stabilizers(Q))
            n_TOF = trellis_oriented_form_additive(normalizermatrix(Q))
            if sect
                opt, _ = optimal_sectionalization_Q(S_TOF, n_TOF)
                return trellis_profiles(S_TOF, n_TOF, opt, "symplectic"), opt
            end
            return trellis_profiles(S_TOF, n_TOF, missing, "symplectic")
        elseif Pauli == 'X'
            _, _, Z_perp, _, _, _ = split_symplectic_stabilizers(quadratic_to_symplectic(normalizermatrix(Q)), ones(Int64, size(normalizermatrix(Q), 1)))
            X = hcat(X_stabilizers(Q), zero_matrix(field(Q), size(X_stabilizers(Q), 1), size(X_stabilizers(Q), 2)))
            Z_perp_TOF = trellis_oriented_form_additive(symplectic_to_quadratic(Z_perp))
            X_TOF = trellis_oriented_form_additive(symplectic_to_quadratic(X))
            if sect
                opt, _ = optimal_sectionalization_Q(X_TOF, Z_perp_TOF)
                return trellis_profiles(X_TOF, Z_perp_TOF, opt, "symplectic"), opt
            end
            return trellis_profiles(X_TOF, Z_perp_TOF, missing, "symplectic")
        else
            X_perp, _, _, _, _, _ = split_symplectic_stabilizers(quadratic_to_symplectic(normalizermatrix(Q)), ones(Int64, size(normalizermatrix(Q), 1)))
            Z = hcat(zero_matrix(field(Q), size(Z_stabilizers(Q), 1), size(Z_stabilizers(Q), 2)), Z_stabilizers(Q))
            X_perp_TOF = trellis_oriented_form_additive(symplectic_to_quadratic(X_perp))
            Z_TOF = trellis_oriented_form_additive(symplectic_to_quadratic(Z))
            if sect
                opt, _ = optimal_sectionalization_Q(Z_TOF, X_perp_TOF)
                return trellis_profiles(Z_TOF, X_perp_TOF, opt, "symplectic"), opt
            end
            return trellis_profiles(Z_TOF, X_perp_TOF, missing, "symplectic")
        end
    end
end

# function syndrome_trellis(Q::AbstractStabilizerCode, type::String="weight", Pauli::Char=' ',
#     sect::Bool=false)
#
#     type ∈ ["weight", "decoding"] || error("Unknown type parameter in syndrome_trellis.")
#     (Pauli != ' ' && typeof(Q) <: CSSCode) && error("Pauli parameter is non-empty but the code is not CSS.")
#     Pauli ∈ [' ', 'X', 'Z'] || error("Unknown Pauli parameter $Pauli; must be ' ', 'X', or 'Z'.")
#
#     if type == "weight"
#         if Pauli == ' '
#             S_TOF = trellis_oriented_form_additive(stabilizers(Q))
#             n_TOF = trellis_oriented_form_additive(normalizermatrix(Q))
#             if sect
#                 profiles, opt = trellis_profiles(Q, type, Pauli, sect)
#                 return _syndrometrellisQ(profiles, opt, n_TOF, S_TOF, character_vector(Q), Pauli, false)
#             else
#                 profiles = trellis_profiles(Q, type, Pauli)
#                 return _syndrometrellisQ(profiles, missing, n_TOF, S_TOF, character_vector(Q), Pauli, false)
#             end
#         elseif Pauli == 'X'
#             _, _, Z_perp, _, _, _ = split_symplectic_stabilizers(quadratic_to_symplectic(normalizermatrix(Q)), ones(Int64, size(normalizermatrix(Q), 1)))
#             X = hcat(X_stabilizers(Q), zero_matrix(field(Q), size(X_stabilizers(Q), 1), size(X_stabilizers(Q), 2)))
#             Z_perp_TOF = trellis_oriented_form_additive(symplectic_to_quadratic(Z_perp))
#             X_TOF = trellis_oriented_form_additive(symplectic_to_quadratic(X))
#             if sect
#                 profiles, opt = trellis_profiles(Q, type, Pauli, sect)
#                 return _syndrometrellisQ(profiles, opt, Z_perp_TOF, X_TOF, character_vector(Q), 'Z', false)
#             else
#                 profiles = trellis_profiles(Q, type, Pauli)
#                 return _syndrometrellisQ(profiles, missing, Z_perp_TOF, X_TOF, character_vector(Q), 'Z', false)
#             end
#         else
#             X_perp, _, _, _, _, _ = split_symplectic_stabilizers(quadratic_to_symplectic(normalizermatrix(Q)), ones(Int64, size(normalizermatrix(Q), 1)))
#             Z = hcat(zero_matrix(field(Q), size(Z_stabilizers(Q), 1), size(Z_stabilizers(Q), 2)), Z_stabilizers(Q))
#             X_perp_TOF = trellis_oriented_form_additive(symplectic_to_quadratic(X_perp))
#             Z_TOF = trellis_oriented_form_additive(symplectic_to_quadratic(Z))
#             if sect
#                 profiles, opt = trellis_profiles(Q, type, Pauli, sect)
#                 return _syndrometrellisQ(profiles, opt, X_perp_TOF, Z_TOF, character_vector(Q), 'X', false)
#             else
#                 profiles = trellis_profiles(Q, type, Pauli)
#                 return _syndrometrellisQ(profiles, missing, X_perp_TOF, Z_TOF, character_vector(Q), 'X', false)
#             end
#         end
#     else
#         if Pauli == ' '
#             S_TOF = trellis_oriented_form_additive(stabilizers(Q))
#             n_TOF = trellis_oriented_form_additive(normalizermatrix(Q))
#             if sect
#                 profiles, opt = trellis_profiles(Q, type, Pauli, sect)
#                 return _syndrometrellisQ(profiles, opt, S_TOF, n_TOF, character_vector(Q), Pauli, false)
#             else
#                 profiles = trellis_profiles(Q, type, Pauli)
#                 return _syndrometrellisQ(profiles, missing, S_TOF, n_TOF, character_vector(Q), Pauli, false)
#             end
#         elseif Pauli == 'X'
#             _, _, Z_perp, _, _, _ = split_symplectic_stabilizers(quadratic_to_symplectic(normalizermatrix(Q)), ones(Int64, size(normalizermatrix(Q), 1)))
#             X = hcat(X_stabilizers(Q), zero_matrix(field(Q), size(X_stabilizers(Q), 1), size(X_stabilizers(Q), 2)))
#             Z_perp_TOF = trellis_oriented_form_additive(symplectic_to_quadratic(Z_perp))
#             X_TOF = trellis_oriented_form_additive(symplectic_to_quadratic(X))
#             if sect
#                 profiles, opt = trellis_profiles(Q, type, Pauli, sect)
#                 return _syndrometrellisQ(profiles, opt, X_TOF, Z_perp_TOF, character_vector(Q), Pauli, false)
#             else
#                 profiles = trellis_profiles(Q, type, Pauli)
#                 return _syndrometrellisQ(profiles, missing, X_TOF, Z_perp_TOF, character_vector(Q), Pauli, false)
#             end
#         else
#             X_perp, _, _, _, _, _ = split_symplectic_stabilizers(quadratic_to_symplectic(normalizermatrix(Q)), ones(Int64, size(normalizermatrix(Q), 1)))
#             Z = hcat(zero_matrix(field(Q), size(Z_stabilizers(Q), 1), size(Z_stabilizers(Q), 2)), Z_stabilizers(Q))
#             X_perp_TOF = trellis_oriented_form_additive(symplectic_to_quadratic(X_perp))
#             Z_TOF = trellis_oriented_form_additive(symplectic_to_quadratic(Z))
#             if sect
#                 profiles, opt = trellis_profiles(Q, type, Pauli, sect)
#                 return _syndrometrellisQ(profiles, opt, Z_TOF, X_perp_TOF, character_vector(Q), Pauli, false)
#             else
#                 profiles = trellis_profiles(Q, type, Pauli)
#                 return _syndrometrellisQ(profiles, missing, Z_TOF, X_perp_TOF, character_vector(Q), Pauli, false)
#             end
#         end
#     end
# end



function sect(C::AbstractCode, type::String="primal", sect::Bool=true, verbose::Bool=false)

    (typeof(C) <: AbstractLinearCode || typeof(C) <: AbstractStabilizerCode) ||
        error("Syndrome trellises are so far only implemented for linear and stabilizer codes.")

    if typeof(C) <: AbstractLinearCode
        wrt_V = trellis_oriented_form_linear(parity_check_matrix(C))
        wrt_E = trellis_oriented_form_linear(generator_matrix(C))
        if sect
            boundaries, num_E_sect = optimal_sectionalization_C(wrt_V, wrt_E)
            profiles = trellis_profiles(wrt_V, wrt_E, boundaries, "Euclidean")
            profiles_no_sect = trellis_profiles(wrt_V, wrt_E, missing, "Euclidean")
            if verbose
                num_E_no_sect = sum(profiles_no_sect[2])
                println("|E| original: $num_E_no_sect, |E| sectionalized: $num_E_sect")
            end
            if length(boundaries) == 2
                if verbose
                    println("Sectionalized to lookup table, skipping")
                end
                boundaries = missing
                profiles = profiles_no_sect
            end
        else
            boundaries = missing
            profiles = trellis_profiles(wrt_V, wrt_E, missing, "Euclidean")
            if verbose
                boundaries2, num_E_sect = optimal_sectionalization_C(wrt_V, wrt_E)
                profiles_sect = trellis_profiles(wrt_V, wrt_E, boundaries2, "Euclidean")
                num_E_no_sect = sum(profiles[2])
                println("|E| original: $num_E_no_sect, |E| sectionalized: $num_E_sect")
            end
        end
    else
        if type == "primal"
            wrt_V = trellis_oriented_form_additive(stabilizers(C))
            wrt_E = trellis_oriented_form_additive(normalizermatrix(C))           
        else
            wrt_V = trellis_oriented_form_additive(normalizermatrix(C))
            wrt_E = trellis_oriented_form_additive(stabilizers(C))
        end
        if sect
            boundaries, num_E_sect = optimal_sectionalization_Q(wrt_V, wrt_E)
            profiles = trellis_profiles(wrt_V, wrt_E, boundaries, "symplectic")
            profiles_no_sect = trellis_profiles(wrt_V, wrt_E, missing, "symplectic")
            if verbose
                num_E_no_sect = sum(profiles_no_sect[2])
                println("|E| original: $num_E_no_sect, |E| sectionalized: $num_E_sect")
            end
            if length(boundaries) == 2
                if verbose
                    println("Sectionalized to lookup table, skipping")
                end
                boundaries = missing
                profiles = profiles_no_sect
            end
        else
            boundaries = missing
            profiles = trellis_profiles(wrt_V, wrt_E, missing, "symplectic")
            if verbose
                boundaries2, num_E_sect = optimal_sectionalization_Q(wrt_V, wrt_E)
                profiles_sect = trellis_profiles(wrt_V, wrt_E, boundaries2, "symplectic")
                num_E_no_sect = sum(profiles[2])
                println("|E| original: $num_E_no_sect, |E| sectionalized: $num_E_sect")
            end
        end
    end

    if ismissing(boundaries)
        bds = [0:C.n...]
    else
        bds = deepcopy(boundaries)
    end
    # println(bds)
    # println(profiles)
    # return

    if typeof(C) <: AbstractLinearCode
        K = C.F
    else
        K = C.E
        R = parent(C.char_vec[1])
    end
    p = Int64(characteristic(K))
    n = C.n
    V = Vector{Vertex}[Vertex[] for _ in 1:length(bds)]
    Threads.@threads for i in 1:length(profiles[1])
        V[i] = [Vertex(999, 0, 0, 0.0, 0, missing) for _ in 1:profiles[1][i]]
    end
    V[1] = [Vertex(0, 0, 0, 0.0, 0, missing)]
    V[end] = [Vertex(0, 0, 0, 0.0, 0, missing)]
    verbose && println("Vertex preallocation completed.")

    E = Vector{Vector{Edge}}[[Edge[]] for _ in 1:length(profiles[3])]
    Threads.@threads for i in 1:length(profiles[3])
        # the j-th element of Ei is going to be all of the edges going into Vi[j]
        E[i] = [[Edge(K(0), 0.0, 0, missing) for _ in 1:profiles[3][i]] for _ in 1:profiles[1][i + 1]]
    end
    verbose && println("Edge preallocation completed.")

    bio = BigInt(1)
    syn_len = nrows(wrt_V)
    active_Vs = _find_active(wrt_V)
    active_Vs = active_Vs[bds[2:end - 1]]
    Threads.@threads for i in 2:length(bds) - 1
        Vi_size = profiles[1][i]
        len_act = length(active_Vs[i - 1])
        # TODO: can I get away with not reversing throughout
        # TODO: do I gain anything from making the V.label correct?
        for num in 0:Vi_size - 1
            # int to small active digits array
            bin = reverse(digits(num, base=p, pad=len_act))
            # to full syn length size
            temp_label = zeros(Int64, syn_len)
            loc = 1
            for j in active_Vs[i - 1]
                temp_label[j] = bin[loc]
                loc += 1
            end
            # i == 3 && println("i = 3: $temp_label")
            # i == 2 && println("i = 2: $temp_label")
            # back to int
            V[i][num + 1].label = digitstoint(reverse(temp_label, dims=1), p)
        end
    end
    verbose && println("Vertex construction completed.")
    # display(wrt_V)
    # display(V)
    # return

    left, right = _left_right_indices(wrt_E)
    # println("left: $left")
    # println("right $right")
    active = _find_active(wrt_E, true)
    # println("active: $active")
    if ismissing(boundaries)
        active_temp = active
        parallel = missing
    else
        active_temp = Vector{Vector{Int64}}()
        parallel = Vector{Vector{Int64}}()
        for i in 1:length(bds) - 1
            temp = sort!(unique!(vcat([active[j] for j in bds[i] + 1:bds[i + 1]]...)))
            # println("temp: $temp")
            act = Vector{Int64}()
            par = Vector{Int64}()
            for a in temp
                # i == 1 && println(a)
                # i == 1 && println(bds[i] + 1, ", ", left[a], ", ", right[a], ", ", bds[i + 1])
                # i == 1 && println(bds[i] + 1 <= left[a] && right[a] <= bds[i + 1])
                (bds[i] + 1 <= left[a] && right[a] <= bds[i + 1]) ? append!(par, a) : append!(act, a)
            end
            push!(active_temp, act)
            push!(parallel, par)
        end
    end
    # println("par: $parallel")

    if typeof(C) <: AbstractLinearCode
        H = FpmattoJulia(wrt_V)
    else
        sym_wrt_V = quadratic_to_symplectic(wrt_V)
        H = FpmattoJulia(hcat(sym_wrt_V[:, n + 1:end], -sym_wrt_V[:, 1:n]))
    end
    
    # Threads.@threads 
    for i in length(bds) - 1:-1:1
        verbose && println("Starting E[$i]")
        
        valid_edges = Vector{fq_nmod_mat}()
        edge_contrib = Dict{fq_nmod_mat, Vector{Int64}}()
        contrib_edge = Dict{Vector{Int64}, fq_nmod_mat}()

        par_flag = false
        if !ismissing(parallel) && !isempty(parallel[i])
            par_flag = true
            p_edges = Vector{fq_nmod_mat}()
            for a in parallel[i]
                # should never be zero because the entire row is between this
                # same argument says it's always unique
                push!(p_edges, wrt_E[a, bds[i] + 1:bds[i + 1]])
            end

            parallel_edges = Vector{fq_nmod_mat}()
            for iter in Nemo.AbstractAlgebra.ProductIterator(collect(0:p - 1), length(p_edges))
                e = K(iter[1]) * p_edges[1]
                for r in 2:length(p_edges)
                    if !iszero(iter[r])
                        e += K(iter[r]) * p_edges[r]
                    end
                end
                !iszero(e) && push!(parallel_edges, e)
            end
            
            if length(parallel_edges) > 1
                p_e_mat_sym = quadratic_to_symplectic(reduce(vcat, parallel_edges))
                temp = symplectic_to_quadratic(_remove_empty(_rref_no_col_swap(p_e_mat_sym, 1:nrows(p_e_mat_sym), 1:ncols(p_e_mat_sym)), :rows))
                parallel_edges = [temp[i, :] for i in 1:nrows(temp)]
            else
                p_e_mat_sym = quadratic_to_symplectic(reduce(vcat, parallel_edges))
            end
        end
        println("i = $i")
        display(parallel_edges)
        # i == 2 && return

        for a in active_temp[i]
            temp = wrt_E[a, bds[i] + 1:bds[i + 1]]
            if !iszero(temp)
                push!(valid_edges, temp)
            end
        end
        # unique!(valid_edges)
        if !isempty(parallel[i])
            v_e_mat_sym = quadratic_to_symplectic(reduce(vcat, valid_edges))
            F = base_ring(v_e_mat_sym)
            VS = VectorSpace(F, ncols(v_e_mat_sym))
            U, U_to_VS = sub(VS, [VS(p_e_mat_sym[i, :]) for i in 1:nrows(p_e_mat_sym)])
            W, W_to_VS = sub(VS, [VS(v_e_mat_sym[i, :]) for i in 1:nrows(v_e_mat_sym)])
            I, _ = intersect(U, W)
            if !iszero(AbstractAlgebra.dim(I))
                println("i = $i, here quo")
                gens_of_U_in_W = [preimage(W_to_VS, U_to_VS(g)) for g in gens(U)]
                U_in_W, _ = sub(W, gens_of_U_in_W)
                Q, W_to_Q = quo(W, U_in_W)
                C2_mod_C1_basis = [W_to_VS(x) for x in [preimage(W_to_Q, g) for g in gens(Q)]]
                F_basis = [[F(C2_mod_C1_basis[j][i]) for i in 1:AbstractAlgebra.dim(parent(C2_mod_C1_basis[1]))] for j in 1:length(C2_mod_C1_basis)]
                temp = symplectic_to_quadratic(matrix(F, length(F_basis), length(F_basis[1]), vcat(F_basis...)))
                valid_edges = [temp[i, :] for i in 1:nrows(temp)]
            else
                temp = symplectic_to_quadratic(_remove_empty(_rref_no_col_swap(v_e_mat_sym, 1:nrows(v_e_mat_sym), 1:ncols(v_e_mat_sym)), :rows))
                valid_edges = [temp[i, :] for i in 1:nrows(temp)]
            end
        else
            temp = symplectic_to_quadratic(_remove_empty(_rref_no_col_swap(v_e_mat_sym, 1:nrows(v_e_mat_sym), 1:ncols(v_e_mat_sym)), :rows))
            valid_edges = [temp[i, :] for i in 1:nrows(temp)]
        end
        println("i = $i")
        display(valid_edges)
        # return
        # i == 2 && return

        for iter in Nemo.AbstractAlgebra.ProductIterator(collect(0:p - 1), length(valid_edges))
            e = K(iter[1]) * valid_edges[1]
            for r in 2:length(valid_edges)
                if !iszero(iter[r])
                    e += K(iter[r]) * valid_edges[r]
                end
            end

            if typeof(C) <: AbstractLinearCode
                P = zeros(Int64, n)
                for (j, k) in enumerate(e)
                    P[bds[i] + j] = coeff(k, 0)
                end
            else
                P = zeros(Int64, 2 * n)
                for (j, k) in enumerate(e)
                    P[bds[i] + j] = coeff(k, 0)
                    P[bds[i] + j + n] = coeff(k, 1)
                end
            end
            syn = H * P .% p
            # if !iszero(syn)
                edge_contrib[e] = syn
                contrib_edge[syn] = e
            # end
        end
        # i == 2 && display(edge_contrib)
        # i == 2 && display(contrib_edge)
        verbose && println("Edges dictionaries completed for E[$i].")

        Vl_len = profiles[1][i]
        Vr_len = profiles[1][i + 1]
        left_vertices = Vector{Tuple{BigInt, Vector{Int}}}()
        right_vertices = Vector{Tuple{BigInt, Vector{Int}}}()
        sizehint!(left_vertices, profiles[4][i])
        sizehint!(right_vertices, profiles[3][i])

        # find fundamental edge configuration
        # find all v-e-0
        left_syn = right_syn = zeros(Int64, syn_len)
        fundamental = Vector{Tuple{BigInt, Vector{Int}, Vector{fq_nmod_mat}, Vector{Int}, BigInt}}()
        for lab in keys(edge_contrib)
            left_syn = (right_syn .- edge_contrib[lab] .+ p) .% p
            if i == 1 && iszero(left_syn)
                edgs = [lab]
                if par_flag
                    for a in parallel_edges
                        push!(edgs, a + lab)
                    end
                end
                push!(fundamental, (bio, left_syn, edgs, right_syn, bio))
                push!(left_vertices, (bio, left_syn))
            elseif i != 1 && iszero(left_syn[setdiff(1:syn_len, active_Vs[i - 1])])
                left_loc = BigInt(digitstoint(reverse(left_syn[active_Vs[i - 1]], dims=1), p)) + 1
                if left_loc <= Vl_len
                    edgs = [lab]
                    if par_flag
                        for a in parallel_edges
                            push!(edgs, a + lab) # this idea is mentally incorrect
                        end
                    end
                    push!(fundamental, (left_loc, left_syn, edgs, right_syn, bio))
                    push!(left_vertices, (left_loc, left_syn))
                end
            end
        end

        # find all 0-e-v
        left_syn = zeros(Int64, syn_len)
        for lab in keys(edge_contrib)
            right_syn = edge_contrib[lab]
            if i != length(bds) - 1 && iszero(right_syn[setdiff(1:syn_len, active_Vs[i])])
                right_loc = BigInt(digitstoint(reverse(right_syn[active_Vs[i]], dims=1), p)) + 1
                if !isone(right_loc)
                    edgs = [lab]
                    if par_flag
                        for a in parallel_edges
                            push!(edgs, a + lab)
                        end
                    end
                    push!(fundamental, (bio, left_syn, edgs, right_syn, right_loc))
                    push!(right_vertices, (right_loc, right_syn))
                end
            end
        end

        # use the above v-e-0 and 0-e-v' to find all v-e-v'
        for (ll, lv) in left_vertices
            for (rl, rv) in right_vertices
                temp = (rv .- lv .+ p) .% p
                # if !iszero(temp)
                    lab = contrib_edge[temp]
                    edgs = [lab]
                    if par_flag
                        for a in parallel_edges
                            push!(edgs, a + lab)
                        end
                    end
                    tup = (ll, lv, edgs, rv, rl)
                    if tup ∉ fundamental
                        push!(fundamental, tup)
                    end
                # end
            end
        end

        # record fundamental in E[i1]
        sort!(fundamental, by=last)
        # i == 2 && 
        println("i = $i")
        display(fundamental)
        # return fundamental
        # println(profiles)
        count = 1
        cur = fundamental[1][5]
        V_right_locs = trues(Vr_len) # TODO: switch to profiles reference
        for (ll, _, edgs, _, rl) in fundamental
            rl == cur || (count = 1; cur = rl;)
            for e in edgs
                # i == 2 && println("here: $i, $cur, $rl, $count")
                E[i][rl][count].label = e
                E[i][rl][count].out_vertex = ll
                if typeof(C) <: AbstractStabilizerCode
                    sign = R(0)
                    for (j, k) in enumerate(e)
                        if !iszero(coeff(k, 0))
                            sign += character_vector(C)[bds[i] + j]
                        end
                        if !iszero(coeff(k, 1))
                            sign += character_vector(C)[bds[i] + j + n]
                        end
                    end                    
                    E[i][rl][count].sign = sign
                end
                count += 1
            end
            V_right_locs[rl] = false
        end

        # i == 2 && display(edge_contrib)
        # i == 2 && return
        # there's nothing but the fundamental in E[1] and E[end]
        if i != 1 && i != length(bds) - 1
            # shift fundamental edge configuration
            r_loc = findfirst(x -> x == true, V_right_locs)
            len_act = length(active_Vs[i])
            while !isnothing(r_loc)
                for edg in keys(edge_contrib)
                    error_syn = edge_contrib[edg]
                    # int to small active digits array
                    bin = reverse(digits(r_loc - 1, base=p, pad=len_act))
                    # to full syn length size
                    right_syn = zeros(Int64, syn_len)
                    loc = 1
                    for j in active_Vs[i]
                        right_syn[j] = bin[loc]
                        loc += 1
                    end
                    # i == 2 && println("right: $right_syn")
                    left_syn = (right_syn .- error_syn .+ p) .% p
                    # i == 2 && println("left: $left_syn")
                    # check if this exists and only shift if it does
                    if iszero(left_syn[setdiff(1:syn_len, active_Vs[i - 1])])
                        # now have v-e-v' not in the fundamental edge configuration
                        # use it to shift
                        # i == 2 && println("in")
                        count = 1
                        cur = fundamental[1][5]
                        for (_, lv, edgs, rv, rl) in fundamental
                            rl == cur || (count = 1; cur = rl;)
                            right_v = (rv .+ right_syn .+ p) .% p
                            # okay to reuse variable here
                            rl = BigInt(digitstoint(reverse(right_v[active_Vs[i]], dims=1), p)) + 1
                            left_v = (lv .+ left_syn .+ p) .% p
                            E[i][rl][count].out_vertex = BigInt(digitstoint(reverse(left_v[active_Vs[i - 1]], dims=1), p)) + 1
                            for e in edgs
                                newe = e + edg
                                E[i][rl][count].label = newe
                                if typeof(C) <: AbstractStabilizerCode
                                    sign = R(0)
                                    for (j, k) in enumerate(newe)
                                        if !iszero(coeff(k, 0))
                                            sign += character_vector(C)[bds[i] + j]
                                        end
                                        if !iszero(coeff(k, 1))
                                            sign += character_vector(C)[bds[i] + j + n]
                                        end
                                    end
                                    E[i][rl][count].sign = sign
                                end
                                count += 1
                            end
                            V_right_locs[r_loc] = false
                        end
                    end
                end
                r_loc = findfirst(x -> x == true, V_right_locs)
                # i == 2 && println(r_loc)
                # i == 2 && println(V_right_locs)
                # i == 2 && display(E[3])
                # i == 2 && return
            end
        end
        verbose && println("E[$i] complete")
    end

    if typeof(C) <: AbstractLinearCode
        return Trellis(V, E, C, missing, zero_matrix(K, 1, n))
    else
        return Trellis(V, E, C, missing, zero_matrix(K, 1, 2 * n))
    end
end
