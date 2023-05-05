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
    edgeloc::Int64
    polynomial::Union{fmpz_mpoly, AbstractAlgebra.Generic.MPoly{nf_elem}, Missing}
end

# could redo this for classical to remove sign
# for a trellis with 10 million edges, this would save 10 MB
# to do this, make EdgeC and EdgeQ and change all following functions to Union both
mutable struct Edge
    label::Union{fq_nmod, fq_nmod_mat}
    weight::Float64
    outvertex::Int64
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
isshifted(T::Trellis) = !iszero(T.shift)

function ==(a::Vertex, b::Vertex)
    return a.label == b.label && a.prev == b.prev && a.next == b.next &&
        a.value == b.value && a.edgeloc == b.edgeloc && a.polynomial == b.polynomial
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
    return a.label == b.label && a.weight == b.weight && a.outvertex == b.outvertex
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

function isisomorphic(V1::Vector{Vector{Vertex}}, V2::Vector{Vector{Vertex}})
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

function isisomorphic(E1::Vector{Vector{Vector{Edge}}}, E2::Vector{Vector{Vector{Edge}}})
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

function isisomorphic(V1::Vector{Vector{Vertex}}, E1::Vector{Vector{Vector{Edge}}},
    V2::Vector{Vector{Vertex}}, E2::Vector{Vector{Vector{Edge}}})

    isoV = isisomorphic(V1, V2)
    if !isoV
        return false
    end
    return isisomorphic(E1, E2)
end

function isisomorphic(V1::Vector{Vector{Vertex}}, V2::Vector{Vector{Vertex}},
    E1::Vector{Vector{Vector{Edge}}}, E2::Vector{Vector{Vector{Edge}}})

    return isisomorphic(V1, E1, V2, E2)
end

# provide "sorted" option
function isequal(V1::Vector{Vector{Vertex}}, V2::Vector{Vector{Vertex}}, verbose=false)
    isiso = isisomorphic(V1, V2)
    verbose && println("V1 and V2 are isomorphic.")

    if !isiso
        verbose && println("V1 and V2 are not isomorphic.")
        return false
    end

    for (i, V1i) in enumerate(V1)
        verbose && println("Checking V1[", i, "]")

        for (j, v1) in enumerate(V1i)
            foundv = false
            for (k, v2) in enumerate(V2[i])
                if v1.label == v2.label
                    foundv = true
                    verbose && println("V1[", i, "][", j, "] = V2[", i, "][", k, "]")
                    break
                end
            end

            if !foundv
                verbose && println("Vertex V1[", i, "][", j, "] not found in V2[", i, "]")
                return false
            end
        end
    end

    return true
end

function ==(V1::Vector{Vector{Vertex}}, V2::Vector{Vector{Vertex}})
    return isequal(V1, V2, false)
end

function isequal(V1::Vector{Vector{Vertex}}, E1::Vector{Vector{Vector{Edge}}},
    V2::Vector{Vector{Vertex}}, E2::Vector{Vector{Vector{Edge}}},
    verbose=false, rhssorted=false)

    isiso = isisomorphic(V1, E1, V2, E2)
    verbose && println("V1 and V2 and E1 and E2 are isomorphic, respectively.")
    if !isiso
        verbose && println("Not  isomorphic.")
        return false
    end

    for (i, V1i) in enumerate(V1)
        lenright = length(V1i)
        verbose && println("Checking V1[", i, "]")

        for (j, v1) in enumerate(V1i)
            foundv = false
            vloc = 0
            if rhssorted
                binsearchleft = 1
                binsearchright = lenright
                while binsearchleft <= binsearchright
                    mid = fld((binsearchleft + binsearchright), 2)
                    if V2[i][mid].label < v1.label
                        binsearchleft = mid + 1
                    elseif V2[i][mid].label > v1.label
                        binsearchright = mid - 1
                    else
                        foundv = true
                        vloc = mid
                        verbose && println("V1[", i, "][", j, "] = V2[", i, "][", mid, "]")
                        break
                    end
                end
            else
                for (k, v2) in enumerate(V2[i])
                    if v1.label == v2.label
                        foundv = true
                        vloc = k
                        verbose && println("V1[", i, "][", j, "] = V2[", i, "][", k, "]")
                        break
                    end
                end
            end

            if !foundv
                verbose && println("Vertex V1[", i, "][", j, "] not found in V2[", i, "]")
                return false
            else
                if i != 1
                    verbose && println("Checking for equality in corresponding edges")
                    for (k, e1) in enumerate(E1[i - 1][j])
                        founde = false
                        for (l, e2) in enumerate(E2[i - 1][vloc])
                            if e1.label == e2.label && e1.weight == e2.weight && V1[i - 1][e1.outvertex].label == V2[i - 1][e2.outvertex].label
                                founde = true
                                verbose && println("E1[", i - 1, "][", j, "][", k, "] = E2[", i - 1, "][", vloc, "][", l, "]")
                                break
                            end
                        end

                        if !founde
                            verbose && println("Edge E1[", i - 1, "][", j, "][", k, "] not found in E2[", i - 1, "][", vloc, "]")
                            return false
                        end
                    end
                end
            end
        end
    end
    return true
end

function isequal(V1::Vector{Vector{Vertex}}, V2::Vector{Vector{Vertex}},
    E1::Vector{Vector{Vector{Edge}}}, E2::Vector{Vector{Vector{Edge}}},
    verbose=false, rhssorted=false)

    return isequal(V1, E1, V2, E2, verbose, rhssorted)
end

function _sortbyleftindex(A::fq_nmod_mat)
    numrows = size(A, 1)
    numcols = size(A, 2)
    arr = []
    for r = 1:numrows
        for c in 1:numcols
            if !iszero(A[r, c])
                push!(arr, [c, A[r, :]])
                break
            end
        end
        # push!(arr, [numcols, A[r, :]])
    end
    # println("sort1: $arr")
    sort!(arr, by=x->x[1])
    # println("sort2: $arr")
    # vcat is faster and cheaper than the following
    # return A[[arr2[i][2] for i in 1:numrows], :]
    return vcat([arr[i][2] for i in 1:numrows]...)
end

function _sortbyrightindex(A::fq_nmod_mat)
    numrows = size(A, 1)
    numcols = size(A, 2)
    arr = []
    for r = 1:numrows
        for c in numcols:-1:1
            if !iszero(A[r, c])
                push!(arr, [c, A[r, :]])
                break
            end
        end
        # push!(arr, [1, A[r, :]])
    end
    sort!(arr, by=x->x[1]) #, rev=true
    return vcat([arr[i][2] for i in 1:numrows]...)
end

function _leftrightindices(A::fq_nmod_mat)
    # iseven(size(A, 2)) || error("Vectors in leftrightindices must have even length.")
    # n = div(size(A, 2), 2)

#take quadratic only here


    left = Vector{Int64}()
    right = Vector{Int64}()
    for i in 1:size(A, 1)
        for j in 1:size(A, 2)
            if !iszero(A[i, j])
                push!(left, j)
                break
            end
        end
    end
    for i in 1:size(A, 1)
        for j in size(A, 2):-1:1
            if !iszero(A[i, j])
                push!(right, j)
                break
            end
        end
    end
    return left, right
end

function _findactive(A::fq_nmod_mat, edges::Bool=false)
    # need to make sure quadratic extension
    k, n = size(A)
    leftright = [[0, 0] for i in 1:k]
    for i in 1:k
        for j in 1:n
            if !iszero(A[i, j])
                leftright[i][1] = j
                break
            end
        end
    end
    for i in 1:k
        for j in n:-1:1
            if !iszero(A[i, j])
                leftright[i][2] = j
                break
            end
        end
    end
    # display(A)
    # println(leftright)

    active = Vector{Vector{Int64}}()
    if edges
        for i in 1:n
            arri = Vector{Int64}()
            for j in 1:k
                # i == 5 && println(i, ", ", j, ", ", leftright[j])
                if leftright[j][1] <= i <= leftright[j][2]
                    # i == 5 && println("added")
                    push!(arri, j)
                end
            end
            if !isempty(arri)
                push!(active, arri)
            end
        end
    else
        for i in 1:n
            arri = Vector{Int64}()
            for j in 1:k
                if leftright[j][1] <= i < leftright[j][2]
                    push!(arri, j)
                end
            end
            if !isempty(arri)
                push!(active, arri)
            end
        end
    end
    # println(active)
    return active
end

# finds all elements of B which have zero symplectic inner product with all of A
function _kernelinnerprod(A::fq_nmod_mat, B::fq_nmod_mat,
    innerprod::String="Euclidean", returnker::Bool=false)

    innerprod ∈ ["Euclidean", "symplectic"] || error("Unsupported inner product type in _kernelinnerprod; expected: `Euclidean`, `symplectic`, received: $innerprod.")
    size(A, 2) == size(B, 2) || error("Length mismatch in _kernelinnerprod.")
    if innerprod == "symplectic"
        iseven(size(A, 2)) || error("Vectors in symplectickernel must have even length.")
    end

    # unclear if it will compile away the choice if we put the check just around
    # the two lines of changed code inside the loops or if it will check it
    # every time
    # can check compiled code later
    if innerprod == "symplectic"
        if returnker
            ker = []
            for rb in 1:size(B, 1)
                for ra in 1:size(A, 1)
                    @views if !iszero(symplecticinnerproduct(A[ra, :], B[rb, :]))
                        push!(ker, B[rb, :])
                        break
                    end
                end
            end
            return length(ker), vcat(ker...)
        else
            AEuc = hcat(A[:, div(size(A, 2), 2) + 1:end], -A[:, 1:div(size(A, 2), 2)])
            prod = AEuc * transpose(B)
            iszero(prod) && return 0
            count = 0
            for i in 1:size(prod, 2)
                if !iszero(prod[:, i])
                    count += 1
                end
            end
            return count
        end
    else
        if returnker
            ker = []
            for rb in 1:size(B, 1)
                for ra in 1:size(A, 1)
                    @views if !iszero(sum([A[ra, i] * B[rb, i] for i in 1:ncols(A)]))
                    # @views if !iszero(A[ra, :] ⋅ B[rb, :])
                        push!(ker, B[rb, :])
                        break
                    end
                end
            end
            return length(ker), vcat(ker...)
        else
            prod = A * transpose(B)
            iszero(prod) && return 0
            count = 0
            for i in 1:size(prod, 2)
                if !iszero(prod[:, i])
                    count += 1
                end
            end
            return count
        end
    end
end

function _pastfuture(A::fq_nmod_mat)
    past = zeros(Int64, size(A, 2) + 1)
    future = zeros(Int64, size(A, 2) + 1)
    left, right = _leftrightindices(A)
    past[1] = 0
    future[1] = size(A, 1)
    for i in 2:size(A, 2) + 1
        past[i] = length(right[right .<= i - 1])
        future[i] = length(left[left .> i - 1])
    end
    return past, future
end

function loadbalancedecode(profile::Vector{Int64})
    leftsum = 0
    leftloc = 1
    rightsum = 0
    rightloc = length(profile)
    while rightloc - leftloc > 2
        if leftsum <= rightsum
            leftsum += profile[leftloc]
            leftloc += 1
        else
            rightsum += profile[rightloc]
            rightloc -= 1
        end
    end
    # println(leftsum, ", ", rightsum)
    return leftloc, rightloc
end

function loadbalancedcode(profiles::Vector{Vector{Int64}})
    length(profiles) == 4 || error("Expected a length 4 profile vector. Pass in all or just the edges.")
    return loadbalancedcode(profiles[2])
end

"""
    trellisorientedformlinear(A::fq_nmod_mat)

Return the trellis oriented form of the matrix `A` assuming the row space is
linear.
"""
function trellisorientedformlinear(A::fq_nmod_mat)
    A = _sortbyleftindex(A)
    for c in 1:size(A, 2)
        left, right = _leftrightindices(A)
        rows = findall(x->x==c, left)
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
    A = _sortbyleftindex(A)
    left, right = _leftrightindices(A)

    for c in size(A, 2):-1:1
        rows = findall(x->x==c, right)
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
        left, right = _leftrightindices(A)
    end
    return A
end

# seeks to eliminate edges of the form a + bω in a quadratic extension
# can't check that it's quadratic directly so note this will fail if a higher degree
# should work for degree 1 extensions given the way AbstractAlgebra stores coefficints?
"""
    trellisorientedformadditive(A::fq_nmod_mat)

Return the trellis oriented form of the matrix `A` assuming the row space is
additive.

# Note
* So far this is only implemented for quadratic extensions over a prime subfield,
 i.e., `GF(p^2)`.
"""
function trellisorientedformadditive(A::fq_nmod_mat)
    E = base_ring(A)
    degree(E) == 2 || error("So far this is only implemented for quadratic extensions over a prime subfield.")

    A = _sortbyleftindex(A)
    # display(A)
    # println(" ")
    @views for c in 1:size(A, 2)
        left, right = _leftrightindices(A)
        rows = findall(x->x==c, left)
        if length(rows) == 1
            if !iszero(coeff(A[rows[1], c], 0)) && !iszero(coeff(A[rows[1], c], 1))
                A[rows[1], :] *= inv(E(coeff(A[rows[1], c], 0)))
            elseif !iszero(coeff(A[rows[1], c], 0))
                A[rows[1], :] *= inv(E(coeff(A[rows[1], c], 0)))
            else
                A[rows[1], :] *= inv(E(coeff(A[rows[1], c], 1)))
            end
        elseif length(rows) > 1
            Xedges = []
            Zedges = []
            mixededges = []
            for row in rows
                if !iszero(coeff(A[row, c], 0)) && !iszero(coeff(A[row, c], 1))
                    push!(mixededges, (row, A[row, c]))
                elseif !iszero(coeff(A[row, c], 0)) && iszero(coeff(A[row, c], 1))
                    push!(Xedges, (row, A[row, c]))
                elseif iszero(coeff(A[row, c], 0)) && !iszero(coeff(A[row, c], 1))
                    push!(Zedges, (row, A[row, c]))
                end
            end
            # println(Xedges)
            # println(Zedges)
            # println(mixededges)

            if length(Xedges) <= 1 && length(Zedges) <= 1 && length(mixededges) <= 1 && (length(Xedges) + length(Zedges) + length(mixededges)) <= 2
                continue
            else
                Xpivot = false
                Zpivot = false
                if !isempty(Xedges)
                    Xpivot = true
                    row, X = Xedges[1]
                    A[row, :] *= inv(E(coeff(X, 0)))
                    # println("X")
                    # display(A)
                    # println(" ")
                    # no problems here if this is only length 1
                    for i in 2:length(Xedges)
                        A[Xedges[i][1], :] -= E(coeff(Xedges[i][2], 0)) * A[row, :]
                    end
                end
                if !isempty(Zedges)
                    Zpivot = true
                    row, Z = Zedges[1]
                    A[row, :] *= inv(E(coeff(Z, 1)))
                    # println("Z")
                    # display(A)
                    # println(" ")
                    # no problems here if this is only length 1
                    for i in 2:length(Zedges)
                        A[Zedges[i][1], :] -= E(coeff(Zedges[i][2], 1)) * A[row, :]
                    end
                end
                if !isempty(mixededges)
                    if Xpivot && Zpivot
                        for i in 1:length(mixededges)
                            A[mixededges[i][1], :] -= E(coeff(mixededges[i][2], 0)) * A[Xedges[1][1], :] + E(coeff(mixededges[i][2], 1)) * A[Zedges[1][1], :]
                        end
                    elseif Xpivot
                        A[mixededges[1][1], :] -= E(coeff(mixededges[1][2], 0)) * A[Xedges[1][1], :]
                        # no problems here if this is only length 1
                        for i in 2:length(mixededges)
                            A[mixededges[i][1], :] -= E(coeff(mixededges[i][2], 0)) * A[Xedges[1][1], :] + E(coeff(mixededges[i][2], 1)) * A[mixededges[1][1], :]
                        end
                    elseif Zpivot
                        A[mixededges[1][1], :] -= E(coeff(mixededges[1][2], 1)) * A[Zedges[1][1], :]
                        # no problems here if this is only length 1
                        for i in 2:length(mixededges)
                            A[mixededges[i][1], :] -= E(coeff(mixededges[i][2], 0)) * A[mixededges[1][1], :] + E(coeff(mixededges[i][2], 1)) * A[Zedges[1][1], :]
                        end
                    else
                        A[mixededges[1][1], :] *= inv(E(coeff(mixededges[1][2], 0)))
                        if length(mixededges) > 1
                            A[mixededges[2][1], :] -= E(coeff(mixededges[2][2], 0)) * A[mixededges[1][1], :]
                            A[mixededges[2][1], :] *= inv(E(coeff(A[mixededges[2][1], c], 1)))
                            if length(mixededges) > 2
                                A[mixededges[3][1], :] -= E(coeff(mixededges[3][2], 1)) * A[mixededges[2][1], :]
                                A[mixededges[3][1], :] *= inv(E(coeff(A[mixededges[3][1], c], 0)))
                                # no problems here if this is only length 3
                                for i in 3:length(mixededges)
                                    A[mixededges[i][1], :] -= E(coeff(mixededges[i][2], 0)) * A[mixededges[3][1], :] + E(coeff(mixededges[i][2], 1)) * A[mixededges[2][1], :]
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
    A = _sortbyleftindex(A)
    left, right = _leftrightindices(A)
    # println("after2")
    # display(A)
    # println(" ")

    @views for c in size(A, 2):-1:1
        rows = findall(x->x==c, right)
        if length(rows) == 1
            if !iszero(coeff(A[rows[1], c], 0)) && !iszero(coeff(A[rows[1], c], 1))
                A[rows[1], :] *= inv(E(coeff(A[rows[1], c], 0)))
            elseif !iszero(coeff(A[rows[1], c], 0))
                A[rows[1], :] *= inv(E(coeff(A[rows[1], c], 0)))
            else
                A[rows[1], :] *= inv(E(coeff(A[rows[1], c], 1)))
            end
        elseif length(rows) > 1
            Xedges = []
            Zedges = []
            mixededges = []
            for row in rows
                if !iszero(coeff(A[row, c], 0)) && !iszero(coeff(A[row, c], 1))
                    push!(mixededges, (row, A[row, c]))
                elseif !iszero(coeff(A[row, c], 0)) && iszero(coeff(A[row, c], 1))
                    push!(Xedges, (row, A[row, c]))
                elseif iszero(coeff(A[row, c], 0)) && !iszero(coeff(A[row, c], 1))
                    push!(Zedges, (row, A[row, c]))
                end
            end

            if length(Xedges) <= 1 && length(Zedges) <= 1 && length(mixededges) <= 1 && (length(Xedges) + length(Zedges) + length(mixededges)) <= 2
                continue
            else
                # need to determine if X and/or Z pivots are below Y and if not, use Y to
                # create one of these pivots and reset the X/Z pivots to a potentially lower row
                if !isempty(mixededges)
                    # can only set one coefficient of a + bω to 1 in additive, do a
                    row, Y = mixededges[end]
                    A[row, :] *= inv(E(coeff(Y, 0)))
                    for i in length(mixededges) - 1:-1:1
                        # can't use a + bω to eliminate a c + dω, so first use 1 + b'ω to eliminate c
                        A[mixededges[i][1], :] -= E(coeff(mixededges[i][2], 0)) * A[row, :]
                        # now turn d to 1 - these are all pure Z's now and some could be pivots
                        if !iszero(coeff(A[mixededges[i][1], c], 1))
                            A[mixededges[i][1], :] *= inv(E(coeff(A[mixededges[i][1], c], 1)))
                            append!(Zedges, (mixededges[i][1], A[mixededges[i][1], c]))
                        end
                    end
                    mixededges = [mixededges[end]]
                    sort!(Zedges, by=x->x[1])
                end

                if !isempty(Xedges)
                    row, X = Xedges[end]
                    A[row, :] *= inv(E(coeff(X, 0)))
                    # no problems here if this is only length 1
                    for i in length(Xedges) - 1:-1:1
                        A[Xedges[i][1], :] -= E(coeff(Xedges[i][2], 0)) * A[row, :]
                    end
                    Xedges = [Xedges[end]]
                end
                if !isempty(Zedges)
                    row, Z = Zedges[end]
                    A[row, :] *= inv(E(coeff(Z, 1)))
                    # no problems here if this is only length 1
                    for i in length(Zedges) - 1:-1:1
                        A[Zedges[i][1], :] -= E(coeff(Zedges[i][2], 1)) * A[row, :]
                    end
                    Zedges = [Zedges[end]]
                end

                # now could have three pivots, remove the top one using the bottom two
                if !isempty(mixededges) && !isempty(Xedges) && !isempty(Zedges)
                    _, loc = findmin([mixededges[1], Xedges[1], Zedges[1]])
                    if loc == 1
                        # Y edge is top, apply both X and Z pivots to remove it
                        # pivot to be removed is of the form 1 + bω
                        A[mixededges[1][1], :] -= A[Xedges[1][1], :]
                        A[mixededges[1][1], :] -= E(coeff(A[mixededges[1][1], c], 1)) * A[Zedges[1][1], :]
                    elseif loc == 2
                        # X edge is top, apply both Y and Z pivots to remove it
                        A[Xedges[1][1], :] -= A[mixededges[1][1], :]
                        # pivot to be remove is now of the form bω coming from the Y
                        A[Xedges[1][1], :] -= E(coeff(A[Xedges[1][1], c], 1)) * A[Zedges[1][1], :]
                    else
                        # Z edge is top, apply both Y and X pivots to remove it
                        A[Zedges[1][1], :] -= inv(E(coeff(mixededges[1][2], 1))) * A[mixededges[1][1], :]
                        # pivot to be remove is now of the form b^{-1}
                        A[Zedges[1][1], :] -= E(coeff(A[Zedges[1][1], c], 0)) * A[Xedges[1][1], :]
                    end
                end
            end
        end
        left, right = _leftrightindices(A)
    end
    return A
end

# need to brainstorem parallel edges
function trellisprofiles(wrtV::fq_nmod_mat, wrtE::fq_nmod_mat,
    boundaries::Union{Vector{Int64}, Missing}, innerprod::String="Euclidean")

    innerprod ∈ ["Euclidean", "symplectic"] || error("Unsupported inner product type in _trellisprofiles; expected: `Euclidean`, `symplectic`, received: $innerprod.")

    if ismissing(boundaries)
        bds = [1:size(wrtV, 2) + 1...]
    else
        bds = deepcopy(boundaries)
        if 0 ∈ bds
            bds .+= 1
        end
    end

    n = length(bds) - 1
    stateprofile = zeros(Int64, n + 1)
    branchprofile = zeros(Int64, n)
    indegrees = zeros(Int64, n)
    outdegrees = zeros(Int64, n)
    past, future = _pastfuture(wrtV)
    past = past[bds]
    future = future[bds]

    # if innerprod == "Euclidean"
    #     dimker = _kernelinnerprod(wrtV, wrtE, innerprod)
    # else
    #     dimker = _kernelinnerprod(quadratictosymplectic(wrtV),
    #         quadratictosymplectic(wrtE), innerprod)
    # end
    dimker = 0

    p = Int64(characteristic(base_ring(wrtV)))
    for i in 1:n + 1
        stateprofile[i] = p^(size(wrtV, 1) - past[i] - future[i] - dimker)
    end

    left, right = _leftrightindices(wrtE)
    past, future = _pastfuture(wrtE)
    past = past[bds]
    future = future[bds]
    for i in 1:n
        dimparallel = 0
        branchprofile[i] = p^(size(wrtE, 1) - past[i] - future[i + 1] - dimker - dimparallel)
        indegrees[i] = div(branchprofile[i], stateprofile[i + 1])
        outdegrees[i] = div(branchprofile[i], stateprofile[i])
    end
    return [stateprofile, branchprofile, indegrees, outdegrees]
end

#############################
        # Classical
#############################

# TODO: handle lookup table better - temporarily skipping
# TODO: remove dictionaries, iterate once to find left, once for right
# keep first bipartite structure and shift it as a coset to find the next ones - fix for sectionalization
function syndrometrellis(C::AbstractCode, type::String="primal", sect::Bool=true,
    verbose::Bool=false)

    (typeof(C) <: AbstractLinearCode || typeof(C) <: AbstractStabilizerCode) ||
        error("Syndrome trellises are so far only implemented for linear and stabilizer codes.")

    if typeof(C) <: AbstractLinearCode
        wrtV = trellisorientedformlinear(paritycheckmatrix(C))
        wrtE = trellisorientedformlinear(generatormatrix(C))
        if sect
            boundaries, numEsect = optimalsectionalizationC(wrtV, wrtE)
            profiles = trellisprofiles(wrtV, wrtE, boundaries, "Euclidean")
            profilesnosect = trellisprofiles(wrtV, wrtE, missing, "Euclidean")
            if verbose
                numEnosect = sum(profilesnosect[2])
                println("|E| original: $numEnosect, |E| sectionalized: $numEsect")
            end
            if length(boundaries) == 2
                if verbose
                    println("Sectionalized to lookup table, skipping")
                end
                boundaries = missing
                profiles = profilesnosect
            end
        else
            boundaries = missing
            profiles = trellisprofiles(wrtV, wrtE, missing, "Euclidean")
            if verbose
                boundaries2, numEsect = optimalsectionalizationC(wrtV, wrtE)
                profilessect = trellisprofiles(wrtV, wrtE, boundaries2, "Euclidean")
                numEnosect = sum(profiles[2])
                println("|E| original: $numEnosect, |E| sectionalized: $numEsect")
            end
        end
    else
        if type == "primal"
            wrtV = trellisorientedformadditive(stabilizers(C))
            wrtE = trellisorientedformadditive(normalizermatrix(C))           
        else
            wrtV = trellisorientedformadditive(normalizermatrix(C))
            wrtE = trellisorientedformadditive(stabilizers(C))
        end
        if sect
            boundaries, numEsect = optimalsectionalizationQ(wrtV, wrtE)
            profiles = trellisprofiles(wrtV, wrtE, boundaries, "symplectic")
            profilesnosect = trellisprofiles(wrtV, wrtE, missing, "symplectic")
            if verbose
                numEnosect = sum(profilesnosect[2])
                println("|E| original: $numEnosect, |E| sectionalized: $numEsect")
            end
            if length(boundaries) == 2
                if verbose
                    println("Sectionalized to lookup table, skipping")
                end
                boundaries = missing
                profiles = profilesnosect
            end
        else
            boundaries = missing
            profiles = trellisprofiles(wrtV, wrtE, missing, "symplectic")
            if verbose
                boundaries2, numEsect = optimalsectionalizationQ(wrtV, wrtE)
                profilessect = trellisprofiles(wrtV, wrtE, boundaries2, "symplectic")
                numEnosect = sum(profiles[2])
                println("|E| original: $numEnosect, |E| sectionalized: $numEsect")
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
        K = quadraticfield(C)
        R = parent(C.charvec[1])
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
        E[i] = [[Edge(K(0), 0.0, 0, missing) for j = 1:profiles[3][i]] for _ in 1:profiles[1][i + 1]]
    end
    verbose && println("Edge preallocation completed.")

    biz = BigInt(0)
    bio = BigInt(1)
    synlen = size(wrtV, 1)
    active = _findactive(wrtV)
    active = active[bds[2:end - 1]]
    Threads.@threads for i in 2:length(bds) - 1
        Visize = profiles[1][i]
        for num in 0:Visize - 1
            bin = reverse(digits(num, base=2, pad=length(active[i - 1])))
            templabel = zeros(Int64, synlen)
            loc = 1
            for j in active[i - 1]
                templabel[j] = bin[loc]
                loc += 1
            end

            curlabel = biz
            for (shift, val) in enumerate(reverse(templabel, dims=1))
                if val == 1
                    curlabel += bio << (shift - 1)
                end
            end
            V[i][num + 1].label = curlabel
        end
    end
    verbose && println("Vertex construction completed.")

    # there has to be a one-liner for the below
    left, right = _leftrightindices(wrtE)
    active = _findactive(wrtE, true)
    if ismissing(boundaries)
        activetemp = active
    else
        activetemp = Vector{Vector{Int64}}()
        for i in 1:length(bds) - 1
            # temp = Vector{Int64}()
            # for j in bds[i] + 1:bds[i + 1]
            #     append!(temp, active[j])
            # end
            # push!(activetemp, sort!(unique!(temp)))
            push!(activetemp, sort!(unique!(vcat([active[j] for j in bds[i] + 1:bds[i + 1]]...))))
        end
    end

    if typeof(C) <: AbstractLinearCode
        H = FpmattoJulia(wrtV)
    else
        symwrtV = quadratictosymplectic(wrtV)
        H = FpmattoJulia(hcat(symwrtV[:, n + 1:end], -symwrtV[:, 1:n]))
    end
    # Threads.@threads 
    for i in length(bds) - 1:-1:1
        verbose && println("Starting E[$i]")
        # seclen = bds[i + 1] - bds[i]
        validedges = Vector{fq_nmod_mat}()
        edgecontrib = Dict{fq_nmod_mat, Vector{Int64}}()
        contribedge = Dict{Vector{Int64}, fq_nmod_mat}()

        for a in activetemp[i]
            temp = wrtE[a, bds[i] + 1:bds[i + 1]]
            if !iszero(temp)
                push!(validedges, temp)
            end
        end
        unique!(validedges)
        # i == 1 && display(validedges)

        for iter in Nemo.AbstractAlgebra.ProductIterator(collect(0:p - 1), length(validedges))
            e = K(iter[1]) * validedges[1]
            for r in 2:length(validedges)
                if !iszero(iter[r])
                    e += K(iter[r]) * validedges[r]
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
            edgecontrib[e] = syn
            contribedge[syn] = e
            # when it's so small it sectionalizes to a lookup table, there is
            # only one syndrome = 0 at Vright = V_n, so contribedge only has
            # one value. to fix, can only handle this case entirely manually
        end
        verbose && println("Edges dictionaries completed for E[$i].")
        # i == 1 && display(edgecontrib)
        # display(edgecontrib)
        # display(contribedge)

        Vleft = V[i]
        Vright = V[i + 1]
        lenleft = length(Vleft)
        lenright = length(Vright)
        Vrightlocs = trues(lenright)
        startingrightindex = 1
        # keep below here instead of preallocating above or else the answer comes out wrong
        blank = K(0)
        # println(lenright)

        while startingrightindex <= lenright
            startingrightv = Vright[startingrightindex].label
            leftvertices = Vector{Tuple{Int64, BigInt, Vector{Int64}}}()
            rightvertices = Vector{Tuple{Int64, BigInt, Vector{Int64}}}()
            sizehint!(leftvertices, profiles[4][i])
            sizehint!(rightvertices, profiles[3][i])

            startingrightvsyn = reverse(digits(startingrightv, base=2, pad=synlen))
            # println("i: $i, synlen: $synlen, srsyn: $startingrightvsyn, srv: $startingrightv")
            push!(rightvertices, (startingrightindex, startingrightv, startingrightvsyn))
            connectingstarts = blank
            startingleftv = biz

            # start with a fixed right vertex and find all left vertices
            for lab in keys(edgecontrib)
                temp = (startingrightvsyn .- edgecontrib[lab])
                for t in 1:length(temp)
                    if temp[t] < 0
                        temp[t] = p + temp[t]
                    end
                end
                temp = temp .% p

                leftlabel = biz
                for (shift, val) in enumerate(reverse(temp, dims=1))
                    if val == 1
                        leftlabel += bio << (shift - 1)
                    end
                end

                binsearchleft = 1
                binsearchright = lenleft
                while binsearchleft <= binsearchright
                    mid = fld((binsearchleft + binsearchright), 2)
                    if Vleft[mid].label < leftlabel
                        binsearchleft = mid + 1
                    elseif Vleft[mid].label > leftlabel
                        binsearchright = mid - 1
                    else
                        push!(leftvertices, (mid, leftlabel, temp))
                        if connectingstarts == blank
                            connectingstarts = lab
                            startingleftv = temp
                        end
                        break
                    end
                end

                if length(leftvertices) == profiles[3][i]
                    break
                end
            end

            # start with first left vertex and find all right vertices
            if length(rightvertices) != profiles[4][i]
                for lab in keys(edgecontrib)
                    if lab != connectingstarts
                        temp = (startingleftv .+ edgecontrib[lab]) .% p
                        rightlabel = biz
                        for (shift, val) in enumerate(reverse(temp, dims=1))
                            if val == 1
                                rightlabel += bio << (shift - 1)
                            end
                        end

                        binsearchleft = 1
                        binsearchright = lenright
                        while binsearchleft <= binsearchright
                            mid = fld((binsearchleft + binsearchright), 2)
                            if Vright[mid].label < rightlabel
                                binsearchleft = mid + 1
                            elseif Vright[mid].label > rightlabel
                                binsearchright = mid - 1
                            else
                                push!(rightvertices, (mid, rightlabel, temp))
                                break
                            end
                        end
                    end

                    if length(rightvertices) == profiles[4][i]
                        break
                    end
                end
            end

            # can probably skip this recalculation of temp by immediately storing
            # instead of building right and left vertex lists
            # should now have all vertices
            for (rightindex, _, rightsyn) in rightvertices
                count = 1
                for (leftindex, _, leftsyn) in leftvertices
                    temp = rightsyn .- leftsyn
                    for t in 1:length(temp)
                        if temp[t] < 0
                            temp[t] = p + temp[t]
                        end
                    end
                    temp = temp .% p
                    lab = contribedge[temp]

                    E[i][rightindex][count].label = lab
                    E[i][rightindex][count].outvertex = leftindex
                    if typeof(C) <: AbstractStabilizerCode
                        sign = R(0)
                        for (j, k) in enumerate(lab)
                            if !iszero(coeff(k, 0))
                                sign += charactervector(C)[bds[i] + j]
                            end
                            if !iszero(coeff(k, 1))
                                sign += charactervector(C)[bds[i] + j + n]
                            end
                        end
                        E[i][rightindex][count].sign = sign
                    end

                    count += 1
                end
                Vrightlocs[rightindex] = false
            end

            # should have ==
            while startingrightindex <= lenright
                startingrightindex += 1
                if startingrightindex <= lenright && !Vrightlocs[startingrightindex]
                    startingrightindex += 1
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
function trellisprofiles(C::AbstractLinearCode, sect::Bool=false)
    GTOF = trellisorientedformlinear(generatormatrix(C))
    HTOF = trellisorientedformlinear(paritycheckmatrix(C))
    if sect
        opt, _ = optimalsectionalizationC(HTOF, GTOF)
        return trellisprofiles(HTOF, GTOF, opt, "Euclidean"), opt
    end
    return trellisprofiles(HTOF, GTOF, missing, "Euclidean")
end

# # recopy above so I don't call it and redo the TOF multiple times
# function syndrometrellis(C::AbstractLinearCode, sect::Bool=false)
#     GTOF = trellisorientedformlinear(generatormatrix(C))
#     HTOF = trellisorientedformlinear(paritycheckmatrix(C))
#     if sect
#         opt, _ = optimalsectionalizationC(HTOF, GTOF)
#         profiles = trellisprofiles(HTOF, GTOF, opt, "Euclidean")
#         return _syndrometrellisC(profiles, opt, HTOF, GTOF, false)
#     else
#         profiles = trellisprofiles(HTOF, GTOF, missing, "Euclidean")
#         return _syndrometrellisC(profiles, missing, HTOF, GTOF, false)
#         # return _syndrometrellisC(profiles, missing, HTOF, GTOF, [1 for _ in 1:length(C)], ' ', false)
#     end
# end

#############################
         # Quantum
#############################

# only valid for quantum codes
function optimalsectionalizationQ(wrtV::fq_nmod_mat, wrtE::fq_nmod_mat)
    K = base_ring(wrtE)
    base_ring(wrtV) == K || error("Vertices and edges must have the same base ring.")

    n = size(wrtV, 2)
    p = Int64(characteristic(K))
    V = [Vertex(i, 0, 0, 0.0, 0, missing) for i in 1:n + 1]
    E = [[Edge(K(0), 0.0, 0, missing) for i in 1:j] for j in n:-1:1]

    symV = quadratictosymplectic(wrtV)
    symE = quadratictosymplectic(wrtE)
    # dimker = _kernelinnerprod(symV, symE, "symplectic")
    dimker = 0
    # println(dimker)
    past, future = _pastfuture(wrtE)
    for i in 1:n
        for j in i:n
            # arbitrary size cutoff
            if size(wrtE, 1) - past[i] - future[j + 1] - dimker > 50
                E[i][j - i + 1].weight = Inf
            else
                E[i][j - i + 1].weight = p^(size(wrtE, 1) - past[i] - future[j + 1] - dimker)
            end
        end
    end

    for i in 1:n
        leftb = 1
        rightb = i
        arr = [E[leftb][rightb].weight + V[leftb].value]

        while leftb < i
            leftb += 1
            rightb -= 1
            append!(arr, E[leftb][rightb].weight + V[leftb].value)
        end

        m, arg = findmin(arr)
        V[i + 1].prev = arg
        V[i + 1].value = m
    end

    sectboundaries = [n]
    next = V[end].prev
    # val = V[end].value
    while next > 0
        append!(sectboundaries, next - 1)
        # val += V[next].value
        next = V[next].prev
    end
    return reverse(sectboundaries), Int(V[end].value) #Int(val)
end

# TODO: remove dictionaries, iterate once to find left, once for right
# keep first bipartite structure and shift it as a coset to find the next ones - fix for sectionalization
# function _syndrometrellisQ(profiles::Vector{Vector{T}}, boundaries::Union{Vector{Int64}, Missing},
#     wrtV::fq_nmod_mat, wrtE::fq_nmod_mat, charvec::Vector{Int64}, Pauli::Char=' ',
#     verbose=false) where T <: Integer
#
#     Pauli ∈ [' ', 'X', 'Z'] || error("Pauli parameter needs to be ' ', 'X', or 'Z'; received $Pauli.")
#     K = base_ring(wrtE)
#     base_ring(wrtV) == K || error("Vertices and edges must have the same base ring.")
#     if ismissing(boundaries)
#         bds = [0:size(wrtV, 2)...]
#     else
#         bds = deepcopy(boundaries)
#     end
#
#     ω = gen(K)
#     p = Int64(characteristic(K))
#     n = length(profiles[1]) - 1
#     symsize = size(wrtV, 2)
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
#     synlen = size(wrtV, 1)
#     active = _findactive(wrtV)
#     active = active[bds[2:end - 1]]
#     Threads.@threads for i = 2:n
#         Visize = profiles[1][i]
#         for num in 0:Visize - 1
#             bin = reverse(digits(num, base=2, pad=length(active[i - 1])))
#             templabel = zeros(Int64, synlen)
#             loc = 1
#             for j in active[i - 1]
#                 templabel[j] = bin[loc]
#                 loc += 1
#             end
#
#             curlabel = biz
#             for (shift, val) in enumerate(reverse(templabel, dims=1))
#                 if val == 1
#                     curlabel += bio << (shift - 1)
#                 end
#             end
#             V[i][num + 1].label = curlabel
#         end
#     end
#     verbose && println("Vertex construction completed.")
#
#     # there has to be a one-liner for the below
#     left, right = _leftrightindices(wrtE)
#     active = _findactive(wrtE, true)
#     if ismissing(boundaries)
#         activetemp = active
#     else
#         activetemp = Vector{Vector{Int64}}()
#         for i in 1:length(bds) - 1
#             temp = Vector{Int64}()
#             for j in bds[i] + 1:bds[i + 1]
#                 if !(bds[i] <= left[j] && right[j] <= bds[i + 1])
#                     append!(temp, active[j])
#                 end
#             end
#             push!(activetemp, sort!(unique!(temp)))
#         end
#     end
#
#     symwrtV = quadratictosymplectic(wrtV)
#     G = FpmattoJulia(hcat(symwrtV[:, symsize + 1:end], -symwrtV[:, 1:symsize]))
#     Threads.@threads for i = n:-1:1
#         verbose && println("Starting E[$i]")
#         seclen = bds[i + 1] - bds[i]
#         validedges = Vector{fq_nmod_mat}()
#         edgecontrib = Dict{fq_nmod_mat, Vector{Int64}}()
#         contribedge = Dict{Vector{Int64}, fq_nmod_mat}()
#
#         for a in activetemp[i]
#             temp = wrtE[a, bds[i] + 1:bds[i + 1]]
#             if !iszero(temp)
#                 push!(validedges, temp)
#             end
#         end
#         unique!(validedges)
#
#         for iter in Nemo.AbstractAlgebra.ProductIterator(collect(0:p - 1), length(validedges))
#             e = K(iter[1]) * validedges[1]
#             for r in 2:length(validedges)
#                 if !iszero(iter[r])
#                     e += K(iter[r]) * validedges[r]
#                 end
#             end
#
#             P = zeros(Int64, 2 * symsize)
#             for (j, k) in enumerate(e)
#                 P[bds[i] + j] = coeff(k, 0)
#                 P[bds[i] + j + symsize] = coeff(k, 1)
#             end
#             syn = G * P .% p
#             edgecontrib[e] = syn
#             contribedge[syn] = e
#         end
#         verbose && println("Edges dictionaries completed for E[$i].")
#
#         Vleft = V[i]
#         Vright = V[i + 1]
#         lenleft = length(Vleft)
#         lenright = length(Vright)
#         Vrightlocs = trues(lenright)
#         startingrightindex = 1
#         # keep below here instead of preallocating above or else the answer comes out wrong
#         blank = K(0)
#
#         while startingrightindex <= lenright
#             startingrightv = Vright[startingrightindex].label
#             leftvertices = Vector{Tuple{Int64, BigInt, Vector{Int64}}}()
#             rightvertices = Vector{Tuple{Int64, BigInt, Vector{Int64}}}()
#             sizehint!(leftvertices, profiles[4][i])
#             sizehint!(rightvertices, profiles[3][i])
#
#             startingrightvsyn = reverse(digits(startingrightv, base=2, pad=synlen))
#             push!(rightvertices, (startingrightindex, startingrightv, startingrightvsyn))
#             connectingstarts = blank
#             startingleftv = biz
#
#             # start with a fixed right vertex and find all left vertices
#             for lab in keys(edgecontrib)
#                 temp = (startingrightvsyn .- edgecontrib[lab])
#                 for t in 1:length(temp)
#                     if temp[t] < 0
#                         temp[t] = p + temp[t]
#                     end
#                 end
#                 temp = temp .% p
#
#                 leftlabel = biz
#                 for (shift, val) in enumerate(reverse(temp, dims=1))
#                     if val == 1
#                         leftlabel += bio << (shift - 1)
#                     end
#                 end
#
#                 binsearchleft = 1
#                 binsearchright = lenleft
#                 while binsearchleft <= binsearchright
#                     mid = fld((binsearchleft + binsearchright), 2)
#                     if Vleft[mid].label < leftlabel
#                         binsearchleft = mid + 1
#                     elseif Vleft[mid].label > leftlabel
#                         binsearchright = mid - 1
#                     else
#                         push!(leftvertices, (mid, leftlabel, temp))
#                         if connectingstarts == blank
#                             connectingstarts = lab
#                             startingleftv = temp
#                         end
#                         break
#                     end
#                 end
#
#                 if length(leftvertices) == profiles[3][i]
#                     break
#                 end
#             end
#
#             # start with first left vertex and find all right vertices
#             if length(rightvertices) != profiles[4][i]
#                 for lab in keys(edgecontrib)
#                     if lab != connectingstarts
#                         temp = (startingleftv .+ edgecontrib[lab]) .% p
#                         rightlabel = biz
#                         for (shift, val) in enumerate(reverse(temp, dims=1))
#                             if val == 1
#                                 rightlabel += bio << (shift - 1)
#                             end
#                         end
#
#                         binsearchleft = 1
#                         binsearchright = lenright
#                         while binsearchleft <= binsearchright
#                             mid = fld((binsearchleft + binsearchright), 2)
#                             if Vright[mid].label < rightlabel
#                                 binsearchleft = mid + 1
#                             elseif Vright[mid].label > rightlabel
#                                 binsearchright = mid - 1
#                             else
#                                 push!(rightvertices, (mid, rightlabel, temp))
#                                 break
#                             end
#                         end
#                     end
#
#                     if length(rightvertices) == profiles[4][i]
#                         break
#                     end
#                 end
#             end
#
#             # can probably skip this recalculation of temp by immediately storing
#             # instead of building right and left vertex lists
#             # should now have all vertices
#             for (rightindex, rightlabel, rightsyn) in rightvertices
#                 count = 1
#                 for (leftindex, leftlabel, leftsyn) in leftvertices
#                     temp = rightsyn .- leftsyn
#                     for t in 1:length(temp)
#                         if temp[t] < 0
#                             temp[t] = p + temp[t]
#                         end
#                     end
#                     temp = temp .% p
#                     lab = contribedge[temp]
#                     sign = 1 # should be K(1) when implementing as roots of unity or in \C?
#                     for (j, k) in enumerate(lab)
#                         if !iszero(coeff(k, 0))
#                             sign *= charvec[bds[i] + j]
#                         end
#                         if !iszero(coeff(k, 1))
#                             sign *= charvec[bds[i] + j + symsize]
#                         end
#                     end
#
#                     E[i][rightindex][count].label = lab
#                     E[i][rightindex][count].outvertex = leftindex
#                     E[i][rightindex][count].sign = sign
#                     count += 1
#                 end
#                 Vrightlocs[rightindex] = false
#             end
#
#             # should have ==
#             while startingrightindex <= lenright
#                 startingrightindex += 1
#                 if startingrightindex <= lenright && !Vrightlocs[startingrightindex]
#                     startingrightindex += 1
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
function weightQ!(T::Trellis, Ps::fq_nmod_mat, err_models::Vector{Dict{String, Float64}},
    weighttype::String="additive")

    weighttype ∈ ["additive", "multiplicative"] || error("Weight type needs to be 'additive' or 'multiplicative'.")

    V = vertices(T)
    E = edges(T)
    for i in 1:length(E)
        model = err_models[i]
        for (j, v) in enumerate(V[i + 1])
            Threads.@threads for e in E[i][j]
                if weighttype == "additive"
                    weight = 0.0
                    for (j, k) in enumerate(e.label)
                        weight += model[j]
                    end
                else
                    weight = 1.0
                    for (j, k) in enumerate(e.label)
                        weight += model[j]
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
# remove wrtV here
function shiftandweightQ!(T::Trellis, Ps::fq_nmod_mat, boundaries::Union{Vector{Int64}, Missing},
    err_models::Vector{Dict{String, Float64}}, charvec::Vector{Int64}, Pauli::Char=' ',
    weighttype::String="additive")

    Pauli ∈ [' ', 'X', 'Z'] || error("Pauli parameter needs to be ' ', 'X', or 'Z'; received $Pauli.")
    weighttype ∈ ["additive", "multiplicative"] || error("Weight type needs to be 'additive' or 'multiplicative'.")
    length(charvec) == 2 * length(err_models) || error("Lengths of character vector and error models are not consistent.")

    K = base_ring(Ps)
    V = vertices(T)
    E = edges(T)
    coden = length(err_models)
    if ismissing(boundaries)
        bds = [0:coden...]
    else
        bds = deepcopy(boundaries)
    end

    for i in 1:length(E)
        model = err_models[i]
        for (j, v) in enumerate(V[i + 1])
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
                if weighttype == "additive"
                    weight = 0.0
                else
                    weight = 1.0
                end
                for (j, k) in enumerate(e.label)
                    if !iszero(coeff(k, 0))
                        sign *= charvec[bds[i] + j]
                    end
                    if !iszero(coeff(k, 1))
                        sign *= charvec[bds[i] + j + coden]
                    end
                    if weighttype == "additive"
                        weight += model[j]
                    else
                        weight *= model[j]
                    end
                end
                e.sign = sign
                e.weight = weight
            end
        end
    end
end

# do I actually care about updating the signs here?
function shiftanddecodeQ!(T::Trellis, Ps::fq_nmod_mat, boundaries::Union{Vector{Int64}, Missing},
    err_models::Vector{Dict{fq_nmod, Float64}}, charvec::Vector{Int64}, Pauli::Char=' ',
    weighttype::String="additive")

    # Pauli ∈ [' ', 'X', 'Z'] || error("Pauli parameter needs to be ' ', 'X', or 'Z'; received $Pauli.")
    # weighttype ∈ ["additive", "multiplicative"] || error("Weight type needs to be 'additive' or 'multiplicative'.")
    # length(charvec) == 2 * length(err_models) || error("Lengths of character vector and error models are not consistent.")

    K = base_ring(Ps)
    V = vertices(T)
    E = edges(T)
    coden = length(err_models)
    if ismissing(boundaries)
        bds = [0:coden...]
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
                if weighttype == "additive"
                    weight = 0.0
                else
                    weight = 1.0
                end
                for (j, k) in enumerate(e.label)
                    # if !iszero(coeff(k, 0))
                    #     sign *= charvec[bds[i] + j]
                    # end
                    # if !iszero(coeff(k, 1))
                    #     sign *= charvec[bds[i] + j + coden]
                    # end
                    if weighttype == "additive"
                        weight += model[k]
                    else
                        weight *= model[k]
                    end
                end
                # e.sign = sign
                e.weight = weight
            end

            # ignoring random tie breaker, breaks ties here by order of entered into trellis
            w, loc = findmin([e.weight + V[i][e.outvertex].value for e in E[i][j]])
            v.value = E[i][j][loc].weight + V[i][E[i][j][loc].outvertex].value
            v.prev = E[i][j][loc].outvertex
            v.edgeloc = loc
        end
    end

    # redo here for the two-sided approach
    path = zero_matrix(base_ring(Ps), 1, coden)
    curr = coden
    prev = 1
    for i in length(E) + 1:-1:2
        elabel = E[i - 1][prev][V[i][prev].edgeloc].label
        path[1, (curr - length(elabel) + 1):curr] = elabel
        prev = V[i][prev].prev
        curr -= length(elabel)
    end
    return path
end

function shift!(T::Trellis, Ps::fq_nmod_mat, boundaries::Union{Vector{Int64}, Missing},
    err_models::Vector{Dict{fq_nmod, Float64}}, charvec::Vector{Int64}, Pauli::Char=' ',
    weighttype::String="additive")

    # Pauli ∈ [' ', 'X', 'Z'] || error("Pauli parameter needs to be ' ', 'X', or 'Z'; received $Pauli.")
    # weighttype ∈ ["additive", "multiplicative"] || error("Weight type needs to be 'additive' or 'multiplicative'.")
    # length(charvec) == 2 * length(err_models) || error("Lengths of character vector and error models are not consistent.")

    K = base_ring(Ps)
    V = vertices(T)
    E = edges(T)
    coden = length(err_models)
    if ismissing(boundaries)
        bds = [0:coden...]
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
                if weighttype == "additive"
                    weight = 0.0
                else
                    weight = 1.0
                end
                for (j, k) in enumerate(e.label)
                    # if !iszero(coeff(k, 0))
                    #     sign *= charvec[bds[i] + j]
                    # end
                    # if !iszero(coeff(k, 1))
                    #     sign *= charvec[bds[i] + j + coden]
                    # end
                    if weighttype == "additive"
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
function trellisprofiles(Q::AbstractStabilizerCode, type::String="weight", Pauli::Char=' ',
    sect::Bool=false)

    type ∈ ["weight", "decoding"] || error("Unknown type parameter in trellisprofiles.")
    # (Pauli != ' ' && typeof(Q) <: CSSCode) && error("Pauli parameter is non-empty but the code is not CSS.")
    Pauli ∈ [' ', 'X', 'Z'] || error("Unknown Pauli parameter $Pauli; must be ' ', 'X', or 'Z'.")

    if type == "weight"
        if Pauli == ' '
            STOF = trellisorientedformadditive(stabilizers(Q))
            nTOF = trellisorientedformadditive(normalizermatrix(Q))
            if sect
                opt, _ = optimalsectionalizationQ(nTOF, STOF)
                return trellisprofiles(nTOF, STOF, opt, "symplectic"), opt
            end
            return trellisprofiles(nTOF, STOF, missing, "symplectic")
        elseif Pauli == 'X'
            _, _, Zperp, _, _, _ = splitsymplecticstabilizers(quadratictosymplectic(normalizermatrix(Q)), ones(Int64, size(normalizermatrix(Q), 1)))
            X = hcat(Xstabilizers(Q), zero_matrix(field(Q), size(Xstabilizers(Q), 1), size(Xstabilizers(Q), 2)))
            ZperpTOF = trellisorientedformadditive(symplectictoquadratic(Zperp))
            XTOF = trellisorientedformadditive(symplectictoquadratic(X))
            if sect
                opt, _ = optimalsectionalizationQ(ZperpTOF, XTOF)
                return trellisprofiles(ZperpTOF, XTOF, opt, "symplectic"), opt
            end
            return trellisprofiles(ZperpTOF, XTOF, missing, "symplectic")
        else
            Xperp, _, _, _, _, _ = splitsymplecticstabilizers(quadratictosymplectic(normalizermatrix(Q)), ones(Int64, size(normalizermatrix(Q), 1)))
            Z = hcat(zero_matrix(field(Q), size(Zstabilizers(Q), 1), size(Zstabilizers(Q), 2)), Zstabilizers(Q))
            XperpTOF = trellisorientedformadditive(symplectictoquadratic(Xperp))
            ZTOF = trellisorientedformadditive(symplectictoquadratic(Z))
            if sect
                opt, _ = optimalsectionalizationQ(XperpTOF, ZTOF)
                return trellisprofiles(XperpTOF, ZTOF, opt, "symplectic"), opt
            end
            return trellisprofiles(XperpTOF, ZTOF, missing, "symplectic")
        end
    else
        if Pauli == ' '
            STOF = trellisorientedformadditive(stabilizers(Q))
            nTOF = trellisorientedformadditive(normalizermatrix(Q))
            if sect
                opt, _ = optimalsectionalizationQ(STOF, nTOF)
                return trellisprofiles(STOF, nTOF, opt, "symplectic"), opt
            end
            return trellisprofiles(STOF, nTOF, missing, "symplectic")
        elseif Pauli == 'X'
            _, _, Zperp, _, _, _ = splitsymplecticstabilizers(quadratictosymplectic(normalizermatrix(Q)), ones(Int64, size(normalizermatrix(Q), 1)))
            X = hcat(Xstabilizers(Q), zero_matrix(field(Q), size(Xstabilizers(Q), 1), size(Xstabilizers(Q), 2)))
            ZperpTOF = trellisorientedformadditive(symplectictoquadratic(Zperp))
            XTOF = trellisorientedformadditive(symplectictoquadratic(X))
            if sect
                opt, _ = optimalsectionalizationQ(XTOF, ZperpTOF)
                return trellisprofiles(XTOF, ZperpTOF, opt, "symplectic"), opt
            end
            return trellisprofiles(XTOF, ZperpTOF, missing, "symplectic")
        else
            Xperp, _, _, _, _, _ = splitsymplecticstabilizers(quadratictosymplectic(normalizermatrix(Q)), ones(Int64, size(normalizermatrix(Q), 1)))
            Z = hcat(zero_matrix(field(Q), size(Zstabilizers(Q), 1), size(Zstabilizers(Q), 2)), Zstabilizers(Q))
            XperpTOF = trellisorientedformadditive(symplectictoquadratic(Xperp))
            ZTOF = trellisorientedformadditive(symplectictoquadratic(Z))
            if sect
                opt, _ = optimalsectionalizationQ(ZTOF, XperpTOF)
                return trellisprofiles(ZTOF, XperpTOF, opt, "symplectic"), opt
            end
            return trellisprofiles(ZTOF, XperpTOF, missing, "symplectic")
        end
    end
end

# function syndrometrellis(Q::AbstractStabilizerCode, type::String="weight", Pauli::Char=' ',
#     sect::Bool=false)
#
#     type ∈ ["weight", "decoding"] || error("Unknown type parameter in syndrometrellis.")
#     (Pauli != ' ' && typeof(Q) <: CSSCode) && error("Pauli parameter is non-empty but the code is not CSS.")
#     Pauli ∈ [' ', 'X', 'Z'] || error("Unknown Pauli parameter $Pauli; must be ' ', 'X', or 'Z'.")
#
#     if type == "weight"
#         if Pauli == ' '
#             STOF = trellisorientedformadditive(stabilizers(Q))
#             nTOF = trellisorientedformadditive(normalizermatrix(Q))
#             if sect
#                 profiles, opt = trellisprofiles(Q, type, Pauli, sect)
#                 return _syndrometrellisQ(profiles, opt, nTOF, STOF, charactervector(Q), Pauli, false)
#             else
#                 profiles = trellisprofiles(Q, type, Pauli)
#                 return _syndrometrellisQ(profiles, missing, nTOF, STOF, charactervector(Q), Pauli, false)
#             end
#         elseif Pauli == 'X'
#             _, _, Zperp, _, _, _ = splitsymplecticstabilizers(quadratictosymplectic(normalizermatrix(Q)), ones(Int64, size(normalizermatrix(Q), 1)))
#             X = hcat(Xstabilizers(Q), zero_matrix(field(Q), size(Xstabilizers(Q), 1), size(Xstabilizers(Q), 2)))
#             ZperpTOF = trellisorientedformadditive(symplectictoquadratic(Zperp))
#             XTOF = trellisorientedformadditive(symplectictoquadratic(X))
#             if sect
#                 profiles, opt = trellisprofiles(Q, type, Pauli, sect)
#                 return _syndrometrellisQ(profiles, opt, ZperpTOF, XTOF, charactervector(Q), 'Z', false)
#             else
#                 profiles = trellisprofiles(Q, type, Pauli)
#                 return _syndrometrellisQ(profiles, missing, ZperpTOF, XTOF, charactervector(Q), 'Z', false)
#             end
#         else
#             Xperp, _, _, _, _, _ = splitsymplecticstabilizers(quadratictosymplectic(normalizermatrix(Q)), ones(Int64, size(normalizermatrix(Q), 1)))
#             Z = hcat(zero_matrix(field(Q), size(Zstabilizers(Q), 1), size(Zstabilizers(Q), 2)), Zstabilizers(Q))
#             XperpTOF = trellisorientedformadditive(symplectictoquadratic(Xperp))
#             ZTOF = trellisorientedformadditive(symplectictoquadratic(Z))
#             if sect
#                 profiles, opt = trellisprofiles(Q, type, Pauli, sect)
#                 return _syndrometrellisQ(profiles, opt, XperpTOF, ZTOF, charactervector(Q), 'X', false)
#             else
#                 profiles = trellisprofiles(Q, type, Pauli)
#                 return _syndrometrellisQ(profiles, missing, XperpTOF, ZTOF, charactervector(Q), 'X', false)
#             end
#         end
#     else
#         if Pauli == ' '
#             STOF = trellisorientedformadditive(stabilizers(Q))
#             nTOF = trellisorientedformadditive(normalizermatrix(Q))
#             if sect
#                 profiles, opt = trellisprofiles(Q, type, Pauli, sect)
#                 return _syndrometrellisQ(profiles, opt, STOF, nTOF, charactervector(Q), Pauli, false)
#             else
#                 profiles = trellisprofiles(Q, type, Pauli)
#                 return _syndrometrellisQ(profiles, missing, STOF, nTOF, charactervector(Q), Pauli, false)
#             end
#         elseif Pauli == 'X'
#             _, _, Zperp, _, _, _ = splitsymplecticstabilizers(quadratictosymplectic(normalizermatrix(Q)), ones(Int64, size(normalizermatrix(Q), 1)))
#             X = hcat(Xstabilizers(Q), zero_matrix(field(Q), size(Xstabilizers(Q), 1), size(Xstabilizers(Q), 2)))
#             ZperpTOF = trellisorientedformadditive(symplectictoquadratic(Zperp))
#             XTOF = trellisorientedformadditive(symplectictoquadratic(X))
#             if sect
#                 profiles, opt = trellisprofiles(Q, type, Pauli, sect)
#                 return _syndrometrellisQ(profiles, opt, XTOF, ZperpTOF, charactervector(Q), Pauli, false)
#             else
#                 profiles = trellisprofiles(Q, type, Pauli)
#                 return _syndrometrellisQ(profiles, missing, XTOF, ZperpTOF, charactervector(Q), Pauli, false)
#             end
#         else
#             Xperp, _, _, _, _, _ = splitsymplecticstabilizers(quadratictosymplectic(normalizermatrix(Q)), ones(Int64, size(normalizermatrix(Q), 1)))
#             Z = hcat(zero_matrix(field(Q), size(Zstabilizers(Q), 1), size(Zstabilizers(Q), 2)), Zstabilizers(Q))
#             XperpTOF = trellisorientedformadditive(symplectictoquadratic(Xperp))
#             ZTOF = trellisorientedformadditive(symplectictoquadratic(Z))
#             if sect
#                 profiles, opt = trellisprofiles(Q, type, Pauli, sect)
#                 return _syndrometrellisQ(profiles, opt, ZTOF, XperpTOF, charactervector(Q), Pauli, false)
#             else
#                 profiles = trellisprofiles(Q, type, Pauli)
#                 return _syndrometrellisQ(profiles, missing, ZTOF, XperpTOF, charactervector(Q), Pauli, false)
#             end
#         end
#     end
# end



function sect(C::AbstractCode, type::String="primal", sect::Bool=true, verbose::Bool=false)

    (typeof(C) <: AbstractLinearCode || typeof(C) <: AbstractStabilizerCode) ||
        error("Syndrome trellises are so far only implemented for linear and stabilizer codes.")

    if typeof(C) <: AbstractLinearCode
        wrtV = trellisorientedformlinear(paritycheckmatrix(C))
        wrtE = trellisorientedformlinear(generatormatrix(C))
        if sect
            boundaries, numEsect = optimalsectionalizationC(wrtV, wrtE)
            profiles = trellisprofiles(wrtV, wrtE, boundaries, "Euclidean")
            profilesnosect = trellisprofiles(wrtV, wrtE, missing, "Euclidean")
            if verbose
                numEnosect = sum(profilesnosect[2])
                println("|E| original: $numEnosect, |E| sectionalized: $numEsect")
            end
            if length(boundaries) == 2
                if verbose
                    println("Sectionalized to lookup table, skipping")
                end
                boundaries = missing
                profiles = profilesnosect
            end
        else
            boundaries = missing
            profiles = trellisprofiles(wrtV, wrtE, missing, "Euclidean")
            if verbose
                boundaries2, numEsect = optimalsectionalizationC(wrtV, wrtE)
                profilessect = trellisprofiles(wrtV, wrtE, boundaries2, "Euclidean")
                numEnosect = sum(profiles[2])
                println("|E| original: $numEnosect, |E| sectionalized: $numEsect")
            end
        end
    else
        if type == "primal"
            wrtV = trellisorientedformadditive(stabilizers(C))
            wrtE = trellisorientedformadditive(normalizermatrix(C))           
        else
            wrtV = trellisorientedformadditive(normalizermatrix(C))
            wrtE = trellisorientedformadditive(stabilizers(C))
        end
        if sect
            boundaries, numEsect = optimalsectionalizationQ(wrtV, wrtE)
            profiles = trellisprofiles(wrtV, wrtE, boundaries, "symplectic")
            profilesnosect = trellisprofiles(wrtV, wrtE, missing, "symplectic")
            if verbose
                numEnosect = sum(profilesnosect[2])
                println("|E| original: $numEnosect, |E| sectionalized: $numEsect")
            end
            if length(boundaries) == 2
                if verbose
                    println("Sectionalized to lookup table, skipping")
                end
                boundaries = missing
                profiles = profilesnosect
            end
        else
            boundaries = missing
            profiles = trellisprofiles(wrtV, wrtE, missing, "symplectic")
            if verbose
                boundaries2, numEsect = optimalsectionalizationQ(wrtV, wrtE)
                profilessect = trellisprofiles(wrtV, wrtE, boundaries2, "symplectic")
                numEnosect = sum(profiles[2])
                println("|E| original: $numEnosect, |E| sectionalized: $numEsect")
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
        R = parent(C.charvec[1])
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
        E[i] = [[Edge(K(0), 0.0, 0, missing) for j = 1:profiles[3][i]] for _ in 1:profiles[1][i + 1]]
    end
    verbose && println("Edge preallocation completed.")

    bio = BigInt(1)
    synlen = nrows(wrtV)
    activeVs = _findactive(wrtV)
    activeVs = activeVs[bds[2:end - 1]]
    Threads.@threads for i in 2:length(bds) - 1
        Visize = profiles[1][i]
        lenact = length(activeVs[i - 1])
        # TODO: can I get away with not reversing throughout
        # TODO: do I gain anything from making the V.label correct?
        for num in 0:Visize - 1
            # int to small active digits array
            bin = reverse(digits(num, base=p, pad=lenact))
            # to full syn length size
            templabel = zeros(Int64, synlen)
            loc = 1
            for j in activeVs[i - 1]
                templabel[j] = bin[loc]
                loc += 1
            end
            # i == 3 && println("i = 3: $templabel")
            # i == 2 && println("i = 2: $templabel")
            # back to int
            V[i][num + 1].label = digitstoint(reverse(templabel, dims=1), p)
        end
    end
    verbose && println("Vertex construction completed.")
    # display(wrtV)
    # display(V)
    # return

    left, right = _leftrightindices(wrtE)
    # println("left: $left")
    # println("right $right")
    active = _findactive(wrtE, true)
    # println("active: $active")
    if ismissing(boundaries)
        activetemp = active
        parallel = missing
    else
        activetemp = Vector{Vector{Int64}}()
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
            push!(activetemp, act)
            push!(parallel, par)
        end
    end
    # println("par: $parallel")

    if typeof(C) <: AbstractLinearCode
        H = FpmattoJulia(wrtV)
    else
        symwrtV = quadratictosymplectic(wrtV)
        H = FpmattoJulia(hcat(symwrtV[:, n + 1:end], -symwrtV[:, 1:n]))
    end
    
    # Threads.@threads 
    for i in length(bds) - 1:-1:1
        verbose && println("Starting E[$i]")
        
        validedges = Vector{fq_nmod_mat}()
        edgecontrib = Dict{fq_nmod_mat, Vector{Int64}}()
        contribedge = Dict{Vector{Int64}, fq_nmod_mat}()

        parflag = false
        if !ismissing(parallel) && !isempty(parallel[i])
            parflag = true
            pedges = Vector{fq_nmod_mat}()
            for a in parallel[i]
                # should never be zero because the entire row is between this
                # same argument says it's always unique
                push!(pedges, wrtE[a, bds[i] + 1:bds[i + 1]])
            end

            paralleledges = Vector{fq_nmod_mat}()
            for iter in Nemo.AbstractAlgebra.ProductIterator(collect(0:p - 1), length(pedges))
                e = K(iter[1]) * pedges[1]
                for r in 2:length(pedges)
                    if !iszero(iter[r])
                        e += K(iter[r]) * pedges[r]
                    end
                end
                !iszero(e) && push!(paralleledges, e)
            end
            
            if length(paralleledges) > 1
                pemathsym = quadratictosymplectic(vcat(paralleledges...))
                temp = symplectictoquadratic(_removeempty(_rref_no_col_swap(pemathsym, 1:nrows(pemathsym), 1:ncols(pemathsym)), :rows))
                paralleledges = [temp[i, :] for i in 1:nrows(temp)]
            else
                pemathsym = quadratictosymplectic(vcat(paralleledges...))
            end
        end
        println("i = $i")
        display(paralleledges)
        # i == 2 && return

        for a in activetemp[i]
            temp = wrtE[a, bds[i] + 1:bds[i + 1]]
            if !iszero(temp)
                push!(validedges, temp)
            end
        end
        # unique!(validedges)
        if !isempty(parallel[i])
            vematsym = quadratictosymplectic(vcat(validedges...))
            F = base_ring(vematsym)
            VS = VectorSpace(F, ncols(vematsym))
            U, UtoVS = sub(VS, [VS(pemathsym[i, :]) for i in 1:nrows(pemathsym)])
            W, WtoVS = sub(VS, [VS(vematsym[i, :]) for i in 1:nrows(vematsym)])
            I, _ = intersect(U, W)
            if !iszero(AbstractAlgebra.dim(I))
                println("i = $i, here quo")
                gensofUinW = [preimage(WtoVS, UtoVS(g)) for g in gens(U)]
                UinW, _ = sub(W, gensofUinW)
                Q, WtoQ = quo(W, UinW)
                C2modC1basis = [WtoVS(x) for x in [preimage(WtoQ, g) for g in gens(Q)]]
                Fbasis = [[F(C2modC1basis[j][i]) for i in 1:AbstractAlgebra.dim(parent(C2modC1basis[1]))] for j in 1:length(C2modC1basis)]
                temp = symplectictoquadratic(matrix(F, length(Fbasis), length(Fbasis[1]), vcat(Fbasis...)))
                validedges = [temp[i, :] for i in 1:nrows(temp)]
            else
                temp = symplectictoquadratic(_removeempty(_rref_no_col_swap(vematsym, 1:nrows(vematsym), 1:ncols(vematsym)), :rows))
                validedges = [temp[i, :] for i in 1:nrows(temp)]
            end
        else
            temp = symplectictoquadratic(_removeempty(_rref_no_col_swap(vematsym, 1:nrows(vematsym), 1:ncols(vematsym)), :rows))
            validedges = [temp[i, :] for i in 1:nrows(temp)]
        end
        println("i = $i")
        display(validedges)
        # return
        # i == 2 && return

        for iter in Nemo.AbstractAlgebra.ProductIterator(collect(0:p - 1), length(validedges))
            e = K(iter[1]) * validedges[1]
            for r in 2:length(validedges)
                if !iszero(iter[r])
                    e += K(iter[r]) * validedges[r]
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
                edgecontrib[e] = syn
                contribedge[syn] = e
            # end
        end
        # i == 2 && display(edgecontrib)
        # i == 2 && display(contribedge)
        verbose && println("Edges dictionaries completed for E[$i].")

        Vllen = profiles[1][i]
        Vrlen = profiles[1][i + 1]
        leftvertices = Vector{Tuple{BigInt, Vector{Int}}}()
        rightvertices = Vector{Tuple{BigInt, Vector{Int}}}()
        sizehint!(leftvertices, profiles[4][i])
        sizehint!(rightvertices, profiles[3][i])

        # find fundamental edge configuration
        # find all v-e-0
        leftsyn = rightsyn = zeros(Int64, synlen)
        fundamental = Vector{Tuple{BigInt, Vector{Int}, Vector{fq_nmod_mat}, Vector{Int}, BigInt}}()
        for lab in keys(edgecontrib)
            leftsyn = (rightsyn .- edgecontrib[lab] .+ p) .% p
            if i == 1 && iszero(leftsyn)
                edgs = [lab]
                if parflag
                    for a in paralleledges
                        push!(edgs, a + lab)
                    end
                end
                push!(fundamental, (bio, leftsyn, edgs, rightsyn, bio))
                push!(leftvertices, (bio, leftsyn))
            elseif i != 1 && iszero(leftsyn[setdiff(1:synlen, activeVs[i - 1])])
                leftloc = BigInt(digitstoint(reverse(leftsyn[activeVs[i - 1]], dims=1), p)) + 1
                if leftloc <= Vllen
                    edgs = [lab]
                    if parflag
                        for a in paralleledges
                            push!(edgs, a + lab) # this idea is mentally incorrect
                        end
                    end
                    push!(fundamental, (leftloc, leftsyn, edgs, rightsyn, bio))
                    push!(leftvertices, (leftloc, leftsyn))
                end
            end
        end

        # find all 0-e-v
        leftsyn = zeros(Int64, synlen)
        for lab in keys(edgecontrib)
            rightsyn = edgecontrib[lab]
            if i != length(bds) - 1 && iszero(rightsyn[setdiff(1:synlen, activeVs[i])])
                rightloc = BigInt(digitstoint(reverse(rightsyn[activeVs[i]], dims=1), p)) + 1
                if !isone(rightloc)
                    edgs = [lab]
                    if parflag
                        for a in paralleledges
                            push!(edgs, a + lab)
                        end
                    end
                    push!(fundamental, (bio, leftsyn, edgs, rightsyn, rightloc))
                    push!(rightvertices, (rightloc, rightsyn))
                end
            end
        end

        # use the above v-e-0 and 0-e-v' to find all v-e-v'
        for (ll, lv) in leftvertices
            for (rl, rv) in rightvertices
                temp = (rv .- lv .+ p) .% p
                # if !iszero(temp)
                    lab = contribedge[temp]
                    edgs = [lab]
                    if parflag
                        for a in paralleledges
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
        Vrightlocs = trues(Vrlen) # TODO: switch to profiles reference
        for (ll, _, edgs, _, rl) in fundamental
            rl == cur || (count = 1; cur = rl;)
            for e in edgs
                # i == 2 && println("here: $i, $cur, $rl, $count")
                E[i][rl][count].label = e
                E[i][rl][count].outvertex = ll
                if typeof(C) <: AbstractStabilizerCode
                    sign = R(0)
                    for (j, k) in enumerate(e)
                        if !iszero(coeff(k, 0))
                            sign += charactervector(C)[bds[i] + j]
                        end
                        if !iszero(coeff(k, 1))
                            sign += charactervector(C)[bds[i] + j + n]
                        end
                    end                    
                    E[i][rl][count].sign = sign
                end
                count += 1
            end
            Vrightlocs[rl] = false
        end

        # i == 2 && display(edgecontrib)
        # i == 2 && return
        # there's nothing but the fundamental in E[1] and E[end]
        if i != 1 && i != length(bds) - 1
            # shift fundamental edge configuration
            rloc = findfirst(x->x==true, Vrightlocs)
            lenact = length(activeVs[i])
            while !isnothing(rloc)
                for edg in keys(edgecontrib)
                    errorsyn = edgecontrib[edg]
                    # int to small active digits array
                    bin = reverse(digits(rloc - 1, base=p, pad=lenact))
                    # to full syn length size
                    rightsyn = zeros(Int64, synlen)
                    loc = 1
                    for j in activeVs[i]
                        rightsyn[j] = bin[loc]
                        loc += 1
                    end
                    # i == 2 && println("right: $rightsyn")
                    leftsyn = (rightsyn .- errorsyn .+ p) .% p
                    # i == 2 && println("left: $leftsyn")
                    # check if this exists and only shift if it does
                    if iszero(leftsyn[setdiff(1:synlen, activeVs[i - 1])])
                        # now have v-e-v' not in the fundamental edge configuration
                        # use it to shift
                        # i == 2 && println("in")
                        count = 1
                        cur = fundamental[1][5]
                        for (_, lv, edgs, rv, rl) in fundamental
                            rl == cur || (count = 1; cur = rl;)
                            rightv = (rv .+ rightsyn .+ p) .% p
                            # okay to reuse variable here
                            rl = BigInt(digitstoint(reverse(rightv[activeVs[i]], dims=1), p)) + 1
                            leftv = (lv .+ leftsyn .+ p) .% p
                            E[i][rl][count].outvertex = BigInt(digitstoint(reverse(leftv[activeVs[i - 1]], dims=1), p)) + 1
                            for e in edgs
                                newe = e + edg
                                E[i][rl][count].label = newe
                                if typeof(C) <: AbstractStabilizerCode
                                    sign = R(0)
                                    for (j, k) in enumerate(newe)
                                        if !iszero(coeff(k, 0))
                                            sign += charactervector(C)[bds[i] + j]
                                        end
                                        if !iszero(coeff(k, 1))
                                            sign += charactervector(C)[bds[i] + j + n]
                                        end
                                    end
                                    E[i][rl][count].sign = sign
                                end
                                count += 1
                            end
                            Vrightlocs[rloc] = false
                        end
                    end
                end
                rloc = findfirst(x->x==true, Vrightlocs)
                # i == 2 && println(rloc)
                # i == 2 && println(Vrightlocs)
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
