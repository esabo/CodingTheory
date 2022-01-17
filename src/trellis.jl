# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

import Base: ==, show

mutable struct Vertex
    label::BigInt
    prev::Int64
    next::Int64
    value::Float64
    edgeloc::Int64
    polynomial::Vector{Vector{Int64}}
end

mutable struct Edge
    label::Union{fq_nmod, fq_nmod_mat}
    weight::Float64
    outvertex::Int64
    sign::Int8
end

mutable struct Trellis
    vertices::Vector{Vector{Vertex}}
    edges::Vector{Vector{Vector{Edge}}}
end

vertices(T::Trellis) = T.vertices
edges(T::Trellis) = T.edges

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
    end
    sort!(arr, by=x->x[1])
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
    n = size(A, 2)
    k = size(A, 1)
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

    active = Vector{Vector{Int64}}()
    if edges
        for i in 1:n
            arri = Vector{Int64}()
            for j in 1:k
                if leftright[j][1] <= i <= leftright[j][2]
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
    return active
end

# finds all elements of B which have zero symplectic inner product with all of A
function _symplectickernel(A::fq_nmod_mat, B::fq_nmod_mat)
    size(A, 2) == size(B, 2) || error("Length mismatch in computedimkernel.")
    iseven(size(A, 2)) || error("Vectors in symplectickernel must have even length.")

    n = div(size(A, 2), 2)
    ker = []
    for rb in 1:size(B, 1)
        flag = false
        for ra in 1:size(A, 1)
            @views if !iszero(symplecticinnerproduct(A[ra, :], B[rb, :]))
                flag = true
                break
            end
        end
        if flag
            push!(ker, B[rb, :])
        end
    end
    return length(ker), vcat(ker...)
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

# seeks to eliminate edges of the form a + bω in a quadratic extension
# can't check that it's quadratic directly so note this will fail if a higher degree
# should work for degree 1 extensions given the way AbstractAlgebra stores coefficints?
# currently untested for characteristic > 2
function trellisorientedform(A::fq_nmod_mat)
    # check for quadratic field extension
    E = base_ring(A)
    A = _sortbyleftindex(A)
    for c in 1:size(A, 2)
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

            if length(Xedges) <= 1 && length(Zedges) <= 1 && length(mixededges) <= 1 && (length(Xedges) + length(Zedges) + length(mixededges)) <= 2
                continue
            else
                Xpivot = false
                Zpivot = false
                if !isempty(Xedges)
                    Xpivot = true
                    row, X = Xedges[1]
                    A[row, :] *= inv(E(coeff(X, 0)))
                    # no problems here if this is only length 1
                    for i in 2:length(Xedges)
                        A[Xedges[i][1], :] -= E(coeff(Xedges[i][2], 0)) * A[row, :]
                    end
                end
                if !isempty(Zedges)
                    Zpivot = true
                    row, Z = Zedges[1]
                    A[row, :] *= inv(E(coeff(Z, 1)))
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
    A = _sortbyleftindex(A)
    left, right = _leftrightindices(A)

    for c in size(A, 2):-1:1
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
                Xpivot = false
                Zpivot = false
                if !isempty(Xedges)
                    Xpivot = true
                    row, X = Xedges[end]
                    A[row, :] *= inv(E(coeff(X, 0)))
                    # no problems here if this is only length 1
                    for i in length(Xedges) - 1:-1:1
                        A[Xedges[i][1], :] -= E(coeff(Xedges[i][2], 0)) * A[row, :]
                    end
                end
                if !isempty(Zedges)
                    Zpivot = true
                    row, Z = Zedges[end]
                    A[row, :] *= inv(E(coeff(Z, 1)))
                    # no problems here if this is only length 1
                    for i in length(Zedges) - 1:-1:1
                        A[Zedges[i][1], :] -= E(coeff(Zedges[i][2], 1)) * A[row, :]
                    end
                end
                if !isempty(mixededges)
                    if Xpivot && Zpivot
                        for i in 1:length(mixededges)
                            A[mixededges[i][1], :] -= E(coeff(mixededges[i][2], 0)) * A[Xedges[end][1], :] + E(coeff(mixededges[i][2], 1)) * A[Zedges[end][1], :]
                        end
                    elseif Xpivot
                        A[mixededges[end][1], :] -= E(coeff(mixededges[1][2], 0)) * A[Xedges[end][1], :]
                        # no problems here if this is only length 1
                        for i in length(mixededges) - 1:-1:1
                            A[mixededges[i][1], :] -= E(coeff(mixededges[i][2], 0)) * A[Xedges[end][1], :] + E(coeff(mixededges[i][2], 1)) * A[mixededges[end][1], :]
                        end
                    elseif Zpivot
                        A[mixededges[end][1], :] -= E(coeff(mixededges[1][2], 1)) * A[Zedges[end][1], :]
                        # no problems here if this is only length 1
                        for i in length(mixededges) - 1:-1:1
                            A[mixededges[i][1], :] -= E(coeff(mixededges[i][2], 0)) * A[mixededges[end][1], :] + E(coeff(mixededges[i][2], 1)) * A[Zedges[end][1], :]
                        end
                    else
                        A[mixededges[end][1], :] *= inv(E(coeff(mixededges[end][2], 0)))
                        if length(mixededges) > 1
                            A[mixededges[end - 1][1], :] -= E(coeff(mixededges[end - 1][2], 0)) * A[mixededges[end][1], :]
                            A[mixededges[end - 1][1], :] *= inv(E(coeff(A[mixededges[end - 1][1], c], 1)))
                            if length(mixededges) > 2
                                A[mixededges[end - 2][1], :] -= E(coeff(mixededges[end - 2][2], 1)) * A[mixededges[end - 1][1], :]
                                A[mixededges[end - 2][1], :] *= inv(E(coeff(A[mixededges[end - 2][1], c], 0)))
                                # no problems here if this is only length 3
                                for i in length(mixededges) - 1:-1:1
                                    A[mixededges[i][1], :] -= E(coeff(mixededges[i][2], 0)) * A[mixededges[end - 2][1], :] + E(coeff(mixededges[i][2], 1)) * A[mixededges[end - 1][1], :]
                                end
                            end
                        end
                    end
                end
            end
        end
        left, right = _leftrightindices(A)
    end
    return A
end

# only valid for quantum codes
function optimalsectionalization(wrtV::fq_nmod_mat, wrtE::fq_nmod_mat)
    K = base_ring(wrtE)
    base_ring(wrtV) == K || error("Vertices and edges must have the same base ring.")

    n = size(wrtV, 2)
    p = Int64(characteristic(K))
    V = [Vertex(i, 0, 0, 0.0, 0, [[0, 0, 0, 0]]) for i in 1:n + 1]
    E = [[Edge(K(0), 0.0, 0, 1) for i in 1:j] for j in n:-1:1]

    symV = quadratictosymplectic(wrtV)
    symE = quadratictosymplectic(wrtE)
    # dimker, _ = _symplectickernel(symV, symE)
    dimker = 0
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
    val = V[end].value
    while next > 0
        append!(sectboundaries, next - 1)
        val += V[next].value
        next = V[next].prev
    end
    return reverse(sectboundaries), Int(val)#Int(V[end].value)
end

function _trellisprofiles(wrtV::fq_nmod_mat, wrtE::fq_nmod_mat, boundaries::Union{Vector{Int64}, Missing})
    if ismissing(boundaries)
        bds = [1:size(wrtV, 2) + 1...]
    else
        bds = deepcopy(boundaries)
        if 0 ∈ bds
            bds .+= 1
        end
    end

    symV = quadratictosymplectic(wrtV)
    symE = quadratictosymplectic(wrtE)
    # dimker, _ = _symplectickernel(symV, symE)
    dimker = 0

    n = length(bds) - 1
    stateprofile = zeros(Int64, n + 1)
    branchprofile = zeros(Int64, n)
    indegrees = zeros(Int64, n)
    outdegrees = zeros(Int64, n)
    past, future = _pastfuture(wrtV)
    past = past[bds]
    future = future[bds]

    p = Int64(characteristic(base_ring(wrtV)))
    for i in 1:n + 1
    #     # same formula works wrtE but still needs a separate loop for the n + 1
    #     # if only doing wrtE need to think about role of dimker
        stateprofile[i] = p^(size(wrtV, 1) - past[i] - future[i] - dimker)
    end

    left, right = _leftrightindices(wrtE)
    past, future = _pastfuture(wrtE)
    past = past[bds]
    future = future[bds]
    for i in 1:n
        dimparallel = 0
        for k in 1:length(left)
            # maybw subtract 1 from all of these?
            if bds[i] <= left[k] && right[k] <= bds[i + 1]
                # println("Going from V[$(boundaries[i])] to V[$(boundaries[i + 1])] and claiming stabilizer $k has left index $(left[k]) and right index $(right[k]).")
                dimparallel += 1
            end
        end
        # stateprofile[i] = p^(size(wrtE, 1) - past[i] - future[i] - dimker)
        branchprofile[i] = p^(size(wrtE, 1) - past[i] - future[i + 1] - dimker - dimparallel)
        indegrees[i] = div(branchprofile[i], stateprofile[i + 1])
        outdegrees[i] = div(branchprofile[i], stateprofile[i])
    end
    # stateprofile[n + 1] = p^(size(wrtE, 1) - past[n + 1] - future[n + 1] - dimker)
    return [stateprofile, branchprofile, indegrees, outdegrees]
end

# TODO: remove dictionaries, iterate once to find left, once for right
# keep first bipartite structure and shift it as a coset to find the next ones - fix for sectionalization
function _syndrometrellisquantum(profiles::Vector{Vector{T}}, boundaries::Union{Vector{Int64}, Missing},
    wrtV::fq_nmod_mat, wrtE::fq_nmod_mat, charvec::Vector{Int64}, Pauli::Char=' ',
    verbose=false) where T <: Integer

    Pauli ∈ [' ', 'X', 'Z'] || error("Pauli parameter needs to be ' ', 'X', or 'Z'; received $Pauli.")
    K = base_ring(wrtE)
    base_ring(wrtV) == K || error("Vertices and edges must have the same base ring.")
    if ismissing(boundaries)
        bds = [0:size(wrtV, 2)...]
    else
        bds = deepcopy(boundaries)
    end

    ω = gen(K)
    p = Int64(characteristic(K))
    n = length(profiles[1]) - 1
    symsize = size(wrtV, 2)
    V = Vector{Vertex}[Vertex[] for i = 1:n + 1]
    Threads.@threads for i = 1:n + 1
        V[i] = [Vertex(999, 0, 0, 0.0, 0, [[0, 0, 0, 0]]) for j = 1:profiles[1][i]]
    end
    V[1] = [Vertex(0, 0, 0, 0.0, 0, [[1, 0, 0, 0]])]
    V[end] = [Vertex(0, 0, 0, 0.0, 0, [[0, 0, 0, 0]])]
    verbose && println("Vertex preallocation completed.")

    E = Vector{Vector{Edge}}[[Edge[]] for i = 1:n]
    Threads.@threads for i in 1:n
        # the j-th element of Ei is going to be all of the edges going into Vi[j]
        E[i] = [[Edge(K(0), 0.0, 0, 1) for j = 1:profiles[3][i]] for k = 1:profiles[1][i + 1]]
    end
    verbose && println("Edge preallocation completed.")

    biz = BigInt(0)
    bio = BigInt(1)
    synlen = size(wrtV, 1)
    active = _findactive(wrtV)
    active = active[bds[2:end - 1]]
    Threads.@threads for i = 2:n
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
            temp = Vector{Int64}()
            for j in bds[i] + 1:bds[i + 1]
                if !(bds[i] <= left[j] && right[j] <= bds[i + 1])
                    append!(temp, active[j])
                end
            end
            push!(activetemp, sort!(unique!(temp)))
        end
    end

    symwrtV = quadratictosymplectic(wrtV)
    G = FpmattoJulia(hcat(symwrtV[:, symsize + 1:end], -symwrtV[:, 1:symsize]))
    Threads.@threads for i = n:-1:1
        verbose && println("Starting E[$i]")
        seclen = bds[i + 1] - bds[i]
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

        for iter in Nemo.AbstractAlgebra.ProductIterator(collect(0:Int64(characteristic(K)) - 1), length(validedges))
            e = K(iter[1]) * validedges[1]
            for r in 2:length(validedges)
                if !iszero(iter[r])
                    e += K(iter[r]) * validedges[r]
                end
            end

            P = zeros(Int64, 2 * symsize)
            for (j, k) in enumerate(e)
                P[bds[i] + j] = coeff(k, 0)
                P[bds[i] + j + symsize] = coeff(k, 1)
            end
            syn = G * P .% p
            edgecontrib[e] = syn
            contribedge[syn] = e
        end
        verbose && println("Edges dictionaries completed for E[$i].")

        Vleft = V[i]
        Vright = V[i + 1]
        lenleft = length(Vleft)
        lenright = length(Vright)
        Vrightlocs = trues(lenright)
        startingrightindex = 1
        # keep below here instead of preallocating above or else the answer comes out wrong
        blank = K(0)

        while startingrightindex <= lenright
            startingrightv = Vright[startingrightindex].label
            leftvertices = Vector{Tuple{Int64, BigInt, Vector{Int64}}}()
            rightvertices = Vector{Tuple{Int64, BigInt, Vector{Int64}}}()
            sizehint!(leftvertices, profiles[4][i])
            sizehint!(rightvertices, profiles[3][i])

            startingrightvsyn = reverse(digits(startingrightv, base=2, pad=synlen))
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
            for (rightindex, rightlabel, rightsyn) in rightvertices
                count = 1
                for (leftindex, leftlabel, leftsyn) in leftvertices
                    temp = rightsyn .- leftsyn
                    for t in 1:length(temp)
                        if temp[t] < 0
                            temp[t] = p + temp[t]
                        end
                    end
                    temp = temp .% p
                    lab = contribedge[temp]
                    sign = 1 # should be K(1) when implementing as roots of unity or in \C?
                    for (j, k) in enumerate(lab)
                        if !iszero(coeff(k, 0))
                            sign *= charvec[bds[i] + j]
                        end
                        if !iszero(coeff(k, 1))
                            sign *= charvec[bds[i] + j + symsize]
                        end
                    end

                    E[i][rightindex][count].label = lab
                    E[i][rightindex][count].outvertex = leftindex
                    E[i][rightindex][count].sign = sign
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
    return Trellis(V, E)
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
    println(leftsum, ", ", rightsum)
    return leftloc, rightloc
end

function loadbalancedcode(profiles::Vector{Vector{Int64}})
    length(profiles) == 4 || error("Expected a length 4 profile vector. Pass in all or just the edges.")
    return loadbalancedcode(profiles[2])
end

# error models need to take CSS combinations into account
# Pauli == 'X'
# I -> I + X
# Z -> Z + Y
# Pauli == 'Z'
# I -> I + Z
# X -> X + Y
function weight!(T::Trellis, Ps::fq_nmod_mat, err_models::Vector{Dict{String, Float64}},
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
function shiftandweight!(T::Trellis, Ps::fq_nmod_mat, boundaries::Union{Vector{Int64}, Missing},
    err_models::Vector{Dict{String, Float64}}, charvec::Vector{Int64}, Pauli::Char=' ',
    weighttype::String="additive")

    Pauli ∈ [' ', 'X', 'Z'] || error("Pauli parameter needs to be ' ', 'X', or 'Z'; received $Pauli.")
    weighttype ∈ ["additive", "multiplicative"] || error("Weight type needs to be 'additive' or 'multiplicative'.")
    length(charvec) == 2 * length(err_models) || error("Lengths of character vector and error models are not consistent.")

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
                        coeff(k, 0) # = 0
                    end
                elseif Pauli == 'Z'
                    for k in e.label
                        coeff(k, 1) # = 0
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
function shiftanddecode!(T::Trellis, Ps::fq_nmod_mat, boundaries::Union{Vector{Int64}, Missing},
    err_models::Vector{Dict{fq_nmod, Float64}}, charvec::Vector{Int64}, Pauli::Char=' ',
    weighttype::String="additive")

    # Pauli ∈ [' ', 'X', 'Z'] || error("Pauli parameter needs to be ' ', 'X', or 'Z'; received $Pauli.")
    # weighttype ∈ ["additive", "multiplicative"] || error("Weight type needs to be 'additive' or 'multiplicative'.")
    # length(charvec) == 2 * length(err_models) || error("Lengths of character vector and error models are not consistent.")

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
                e.label += Ps[1, bds[i] + 1:bds[i + 1]]
                if Pauli == 'X'
                    for k in e.label
                        coeff(k, 0) # = 0
                    end
                elseif Pauli == 'Z'
                    for k in e.label
                        coeff(k, 1) # = 0
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
                # if Pauli == 'X'
                #     for k in e.label
                #         coeff(k, 0) # = 0
                #     end
                # elseif Pauli == 'Z'
                #     for k in e.label
                #         coeff(k, 1) # = 0
                #     end
                # end

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
