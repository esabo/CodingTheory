# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

include("quantumcode.jl")

import Base: ==, show

mutable struct Vertex
    label::BigInt
    prev::Int64
    value::Float64
    edgeloc::Int64
    polynomial::Vector{Vector{Int64}}
end

mutable struct Edge
    label::String
    weight::Float64
    outvertex::Int64
end

mutable struct Trellis
    vertices::Vector{Vector{Vertex}}
    edges::Vector{Vector{Vector{Edge}}}
end

struct WeightEnumerator
    polynomial::Vector{Vector{Int64}}
end

function islessPoly(a::Vector{Int64}, b::Vector{Int64})
    if isless(a[2:4], b[2:4])
        return true
    elseif a[2:4] == b[2:4]
        if isless(a[1], b[1])
            return true
        else
            return false
        end
    else
        return false
    end
end

function show(io::IO, W::WeightEnumerator)
    # sort lexicographically
    for (i, term) in enumerate(W.polynomial)
        if term == [1, 0, 0, 0]
            print("1")
        end
        if term[1] != 1
            print(io, term[1])
        end
        if term[2] != 0
            if term[2] > 1
                print(io, "x^", term[2])
            else
                print(io, "x")
            end
        end
        if term[3] != 0
            if term[3] > 1
                print(io, "y^", term[3])
            else
                print(io, "y")
            end
        end
        if term[4] != 0
            if term[4] > 1
                print(io, "z^", term[4])
            else
                print(io, "z")
            end
        end
        if i != length(W.polynomial)
            print(io, " + ")
        else
            print(io, "\n")
        end
    end
end

vertices(T::Trellis) = T.vertices
edges(T::Trellis) = T.edges

function ==(a::Vertex, b::Vertex)
    return a.label == b.label && a.prev == b.prev && a.value == b.value && a.edgeloc == b.edgeloc
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

function _findactive(A::fq_nmod_mat)
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
            if !iszero(symplecticinnerproduct(A[ra, :], B[rb, :]))
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

function trellisorientedform(A::fq_nmod_mat)
    E = base_ring(A)
    Int64(characteristic(E)) == 2 || error("Have not yet programmed the elimination of the TOF for p ≥ 3.")
    A = _sortbyleftindex(A)
    numrows = size(A, 1)
    numcols = size(A, 2)

    row = 1
    col = 1
    while row < numrows
        firstpivot = [row, A[row, col]]
        secondpivot = [0, E(0)]
        others = []
        while row <= numrows && !iszero(A[row, col])
            if secondpivot == [0, E(0)] && A[row, col] != firstpivot[2]
                secondpivot = [row, A[row, col]]
            elseif [row, A[row, col]] != firstpivot
                push!(others, [row, A[row, col]])
            end
            row += 1
        end

        if secondpivot == [0, E(0)]
            row = firstpivot[1] + 1
        else
            row = firstpivot[1] + 2
        end

        if !isempty(others)
            Afspiv = 0
            if !iszero(secondpivot[1])
                Afspiv = A[firstpivot[1], :] + A[secondpivot[1], :]
            end
            for o in others
                if o[2] == firstpivot[2]
                    A[o[1], :] = A[firstpivot[1], :] + A[o[1], :]
                elseif o[2] == secondpivot[2]
                    A[o[1], :] = A[secondpivot[1], :] + A[o[1], :]
                elseif !iszero(secondpivot[1])
                    A[o[1], :] = Afspiv + A[o[1], :]
                end
            end
        end
        # println(A)
        A = _sortbyleftindex(A)
        col += 1
    end

    A = _sortbyrightindex(A)
    row = numrows
    col = numcols
    while row > 1
        others = []
        temprow = row
        while temprow >= 1 && !iszero(A[temprow, col])
            # here must choose pivots by left index
            for i in 1:numcols
                if !iszero(A[temprow, i])
                    push!(others, [temprow, A[temprow, col], i])
                    break
                end
            end
            temprow -= 1
        end

        if !isempty(others)
            sort!(others, by=x->x[3])
            firstpivot = others[end]
            pop!(others)
            secondpivot = [0, E(0)]
            for (i, o) in enumerate(reverse!(others))
                if o[2] != firstpivot[2]
                    secondpivot = o
                    deleteat!(others, i)
                    break
                end
            end

            if secondpivot == [0, E(0)]
                row -= 1
            else
                row -= 2
            end

            Afspiv = 0
            if !iszero(secondpivot[1])
                Afspiv = A[firstpivot[1], :] + A[secondpivot[1], :]
            end
            for o in others
                if o[2] == firstpivot[2]
                    A[o[1], :] = A[firstpivot[1], :] + A[o[1], :]
                elseif o[2] == secondpivot[2]
                    A[o[1], :] = A[secondpivot[1], :] + A[o[1], :]
                elseif !iszero(secondpivot[1])
                    A[o[1], :] = Afspiv + A[o[1], :]
                end
            end
            A = _sortbyrightindex(A)
        end
        col -= 1
    end
    A = _sortbyleftindex(A)
    return A
end

# need to include sectionalization - check other codebase for this
function _trellisprofiles(wrtV::fq_nmod_mat, wrtE::fq_nmod_mat)
    symV = quadratictosymplectic(wrtV)
    symE = quadratictosymplectic(wrtE)
    dimker, _ = _symplectickernel(symV, symE)

    n = size(wrtV, 2)
    stateprofile = zeros(Int64, n + 1)
    branchprofile = zeros(Int64, n)
    indegrees = zeros(Int64, n)
    outdegrees = zeros(Int64, n)
    past, future = _pastfuture(wrtV)

    p = Int64(characteristic(base_ring(wrtV)))
    for i in 1:n + 1
        # same formula works wrtE but still needs a separate loop for the n + 1
        # if only doing wrtE need to think about role of dimker
        stateprofile[i] = p^(size(wrtV, 1) - past[i] - future[i] - dimker)
    end

    past, future = _pastfuture(wrtE)
    for i in 1:n
        branchprofile[i] = p^(size(wrtE, 1) - past[i] - future[i + 1] - dimker)
        indegrees[i] = div(branchprofile[i], stateprofile[i + 1])
        outdegrees[i] = div(branchprofile[i], stateprofile[i])
    end
    return [stateprofile, branchprofile, indegrees, outdegrees]
end

# think of more scenarios
# could allow general trellises given partial stabilizers for use in trellis product
function trellisprofiles(Q::AbstractQuantumCode, type::String="weight")
    type ∈ ["weight", "decoding"] || error("Unknown type parameter in trellisprofiles.")
    if type == "weight"
        return _trellisprofiles(Q.dualgens, stabilizers(Q))
    else
        return _trellisprofiles(stabilizers(Q), Q.dualgens)
    end
end

# need to fix all of these stabilizer grabs
# I need to grab the stabs and make them quadratic before passing them
function trellisprofilesCSS(Q::CSSCode, type::String="weight", Pauli::Char=' ')
    type ∈ ["weight", "decoding"] || error("Unknown type parameter in trellisprofilesCSS.")
    Pauli ∈ ['X', 'Z'] || error("Pauli parameter in trellisprofilesCSS must be either 'X' or 'Z'.")
    if type == "weight"
        if Pauli == 'X'
            # this should be wrtV = SperpZ, wrtE = SX
            return _trellisprofiles(Q.dualgens, stabilizers(Q))
        else
            # this should be wrtV = SperpX, wrtE = SZ
            return _trellisprofiles(Q.dualgens, stabilizers(Q))
        end
    else
        if Pauli == 'X'
            # this should be wrtV = SX, wrtE = SperpZ
            return _trellisprofiles(stabilizers(Q), Q.dualgens)
        else
            # this should be wrtV = SZ, wrtE = SperpX
            return _trellisprofiles(stabilizers(Q), Q.dualgens)
        end
    end
end

# ignoring signs for now
# straight copy for now, redo signature later
# assumes all stabilizers are given in positive powers
function _syndrometrellis(profiles::Vector{Vector{T}}, wrtV::fq_nmod_mat, wrtE::fq_nmod_mat,
    CSS::Bool=false, type::Char=' ', verbose=false) where T <: Integer

    if CSS && type ∉ ['X', 'Z']
        error("CSS flag detected but type flag is neither X or Z.")
    end

    q = Int64(characteristic(base_ring(wrtV)))
    n = size(wrtV, 2)
    V = Vector{Vertex}[Vertex[] for i = 1:n + 1]
    Threads.@threads for i = 1:length(profiles[1])
        V[i] = [Vertex(999, 0, 0.0, 0, [[0, 0, 0, 0]]) for j = 1:profiles[1][i]]
    end
    V[1] = [Vertex(0, 0, 0.0, 0, [[1, 0, 0, 0]])]
    V[end] = [Vertex(0, 0, 0.0, 0, [[0, 0, 0, 0]])]
    verbose && println("Vertex preallocation completed.")

    # this code and the edge code below it has large problems with GC due to the lack of proper preallocation
    # however, using, for example, the code below to fix this switchs from 60% GC time to 99% compile time
    # rerunning it then works with almost no GC but is slower than doing the above code in parallel with
    # tons of GC
    # V = [if i == 1 || i == C.n + 1
    #         SVector{1, Vertex}(Vertex(0, 0x000000, 0.0, 0x000000))
    #     else
    #         SVector{Int(profiles[1][i]), Vertex}(Vertex(999, 0x000000, 0.0, 0x000000) for j = 1:profiles[1][i])
    #     end
    #     for i = 1:C.n + 1]

    E = Vector{Vector{Edge}}[[Edge[]] for i = 1:n]
    Threads.@threads for i in 1:n
        # the j-th element of Ei is going to be all of the edges going into Vi[j]
        E[i] = [[Edge("e", 0.0, 0) for j = 1:profiles[3][i]] for k = 1:profiles[1][i + 1]]
    end
    verbose && println("Edge preallocation completed.")

    biz = BigInt(0)
    bio = BigInt(1)
    # blank = [0xff, 0xff]
    synlen = size(wrtV, 1)
    # M = MatrixSpace(base_ring(wrtV), 1, synlen)
    active = _findactive(wrtV)
    Threads.@threads for i = 2:n
        Visize = profiles[1][i]
        for num in 0:Visize - 1
            bin = digits(num, base=2, pad=length(active[i - 1])) |> reverse
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
    # println(V)

    valid = []
    if !CSS
        valid = ["I", "X", "Y", "Z"]
    else
        if type == 'X'
            valid = ["I", "Z"]
        else
            valid = ["I", "X"]
        end
    end

    symwrtV = quadratictosymplectic(wrtV)
    G = FpmattoJulia(hcat(symwrtV[:, n + 1:end], -symwrtV[:, 1:n]))
    # make this n?
    # Threads.@threads
    for i = length(E):-1:1
        verbose && println("Starting E[$i]")
        edgecontrib = Dict{String, Vector{Int64}}()
        contribedge = Dict{Vector{Int64}, String}()
        for lab in valid
            if lab == "I"
                syn = zeros(Int64, size(wrtV, 1))
                edgecontrib[lab] = syn
                contribedge[syn] = lab
            else
                # I can deterine this value without building this by something like
                # syn1 = [wrtV[j][i] for j in 1:length(wrtV)]
                PL = zeros(Int64, 2 * n)
                if lab == "X"
                    PL[i] = 1
                elseif lab == "Y"
                    PL[i] = 1
                    PL[i + n] = 1
                elseif lab == "Z"
                    PL[i + n] = 1
                else
                    error("Invalid edge label ($lab) in _syndrometrellis.")
                end
                syn = G * PL
                edgecontrib[lab] = syn
                contribedge[syn] = lab
            end
        end
        # println(edgecontrib)
        # println(contribedge)

        Vleft = V[i]
        Vright = V[i + 1]
        lenleft = length(Vleft)
        lenright = length(Vright)
        Vrightlocs = trues(lenright)
        startingrightindex = 1
        # keep below here instead of preallocating above or else the answer comes out wrong
        blank = "e"

        while startingrightindex <= lenright
            startingrightv = Vright[startingrightindex].label
            leftvertices = Vector{Tuple{Int64, BigInt, Vector{Int64}}}()
            rightvertices = Vector{Tuple{Int64, BigInt, Vector{Int64}}}()
            sizehint!(leftvertices, profiles[4][i])
            sizehint!(rightvertices, profiles[3][i])

            startingrightvsyn = digits(startingrightv, base=2, pad=synlen) |> reverse
            push!(rightvertices, (startingrightindex, startingrightv, startingrightvsyn))
            connectingstarts = blank
            startingleftv = biz

            # start with a fixed right vertex and find all left vertices
            for lab in valid
                temp = (startingrightvsyn .- edgecontrib[lab])
                for t in 1:length(temp)
                    if temp[t] < 0
                        temp[t] = q + temp[t]
                    end
                end
                temp = temp .% q

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
                for lab in valid
                    if lab != connectingstarts
                        temp = (startingleftv .+ edgecontrib[lab]) .% q
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
                    # else
                    #     println("skipping $lab here")
                    end

                    if length(rightvertices) == profiles[4][i]
                        break
                    end
                end
            end

            # if i == 1
            #     println(rightvertices)
            # end
            # should now have all vertices
            for (rightindex, rightlabel, rightsyn) in rightvertices
                count = 1
                for (leftindex, leftlabel, leftsyn) in leftvertices
                    temp = rightsyn .- leftsyn
                    for t in 1:length(temp)
                        if temp[t] < 0
                            temp[t] = q + temp[t]
                        end
                    end
                    temp = temp .% q
                    E[i][rightindex][count].label = contribedge[temp]
                    E[i][rightindex][count].outvertex = leftindex
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


############# need TOF here
function syndrometrellis(Q::AbstractQuantumCode, type::String="weight")
    type ∈ ["weight", "decoding"] || error("Unknown type parameter in syndrometrellis.")
    profiles = trellisprofiles(Q, type)
    if type == "weight"
        return _syndrometrellis(Q.dualgens, stabilizers(Q), profiles)
    else
        return _syndrometrellis(stabilizers(Q), Q.dualgens, profiles)
    end
end

# need to fix all of these stabilizer grabs
# I need to grab the stabs and make them quadratic before passing them
# how to pass CSS into this function? should separate?
function syndrometrellisCSS(Q::CSSCode, type::String="weight", Pauli::Char=' ')
    type ∈ ["weight", "decoding"] || error("Unknown type parameter in syndrometrellisCSS.")
    Pauli ∈ ['X', 'Z'] || error("Pauli parameter in syndrometrellisCSS must be either 'X' or 'Z'.")
    profiles = trellisprofilesCSS(Q, type, Pauli)

    if type == "weight"
        if Pauli == 'X'
            # this should be wrtV = SperpZ, wrtE = SX
            return _syndrometrellis(Q.dualgens, stabilizers(Q), profiles)
        else
            # this should be wrtV = SperpX, wrtE = SZ
            return _syndrometrellis(Q.dualgens, stabilizers(Q), profiles)
        end
    else
        if Pauli == 'X'
            # this should be wrtV = SX, wrtE = SperpZ
            return _syndrometrellis(stabilizers(Q), Q.dualgens, profiles)
        else
            # this should be wrtV = SZ, wrtE = SperpX
            return _syndrometrellis(stabilizers(Q), Q.dualgens, profiles)
        end
    end
end

function _reducepoly(poly::Vector{Vector{Int64}})
    reducedpoly = Vector{Vector{Int64}}()
    processed = trues(length(poly))
    # println("received: $poly")
    # println("processed: $processed")
    for (i, term) in enumerate(poly)
        if processed[i]
            for j in (i + 1):length(poly)
                if processed[j]
                    # println("comparing $term and $(poly[j])")
                    if term[2] == poly[j][2] && term[3] == poly[j][3] && term[4] == poly[j][4]
                        term[1] += poly[j][1]
                        processed[j] = false
                    end
                end
            end
            push!(reducedpoly, term)
            processed[i] = false
        end
    end
    return reducedpoly
end

function weightenumerator(T::Trellis, cleanV::Bool=true)
    V = vertices(T)
    E = edges(T)
    for i in 2:length(V)
        for (j, v) in enumerate(V[i])
            outer = Vector{Vector{Int64}}()
            for e in E[i - 1][j]
                inner = deepcopy(V[i - 1][e.outvertex].polynomial)
                # println("edge: $(e.label), inner: $inner")
                # if e.label == "I"
                #     continue
                if e.label != "I"
                    for term in inner
                        term[2] += 1
                    end
                else
                    for term in inner
                        term[3] += 1
                    end
                end
                append!(outer, inner)
            end
            # println(outer)
            v.polynomial = _reducepoly(outer)
            # println(v.polynomial)
        end
    end
    W = WeightEnumerator(sort!(V[end][1].polynomial, lt=islessPoly))

    if cleanV
        for i in 2:length(V)
            for v in V[i]
                v.polynomial = [[0, 0, 0, 0]]
            end
        end
    end
    return W
end

function Pauliweightenumerator(T::Trellis, cleanV::Bool=true)
    V = vertices(T)
    E = edges(T)
    for i in 2:length(V)
        for (j, v) in enumerate(V[i])
            outer = Vector{Vector{Int64}}()
            for e in E[i - 1][j]
                inner = deepcopy(V[i - 1][e.outvertex].polynomial)
                # println("edge: $(e.label), inner: $inner")
                # if e.label == "I"
                #     continue
                if e.label == "X"
                    for term in inner
                        if term[2] < 0 || term[3] < 0 || term[4] < 0
                            term[2] -= 1
                        else
                            term[2] += 1
                        end
                    end
                elseif e.label == "Y"
                    for term in inner
                        if term[2] < 0 || term[3] < 0 || term[4] < 0
                            term[3] -= 1
                        else
                            term[3] += 1
                        end
                    end
                elseif e.label == "Z"
                    for term in inner
                        if term[2] < 0 || term[3] < 0 || term[4] < 0
                            term[4] -= 1
                        else
                            term[4] += 1
                        end
                    end
                elseif e.label == "-X"
                    for term in inner
                        term[2] *= -1
                        term[3] *= -1
                        term[4] *= -1
                        if term[2] < 0 || term[3] < 0 || term[4] < 0
                            term[2] += 1
                        else
                            term[2] += 1
                        end
                    end
                elseif e.label == "-Y"
                    for term in inner
                        term[2] *= -1
                        term[3] *= -1
                        term[4] *= -1
                        if term[2] < 0 || term[3] < 0 || term[4] < 0
                            term[3] += 1
                        else
                            term[3] += 1
                        end
                    end
                elseif e.label == "-Z"
                    for term in inner
                        term[2] *= -1
                        term[3] *= -1
                        term[4] *= -1
                        if term[2] < 0 || term[3] < 0 || term[4] < 0
                            term[4] += 1
                        else
                            term[4] += 1
                        end
                    end
                end
                append!(outer, inner)
            end
            # println(outer)
            v.polynomial = _reducepoly(outer)
            # println(v.polynomial)
        end
    end
    W = WeightEnumerator(sort!(V[end][1].polynomial, lt=islessPoly))

    if cleanV
        for i in 2:length(V)
            for v in V[i]
                v.polynomial = [[0, 0, 0, 0]]
            end
        end
    end
    return W
end

function testtrellisSteane()
    Q = QuantumCode(["XXIXXII", "IXXIXXI", "IIIXXXX", "ZZIZZII", "IZZIZZI", "IIIZZZZ"], [1, 1, 1, 1, 1, 1]);
    F = field(Q)
    test = symplectictoquadratic(matrix(F,
        [1 1 1 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 1 1 1 0 0 0 0;
         0 1 0 1 0 1 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 1 0 1 0 1 0;
         0 0 1 1 1 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 1 1 1 0 0;
         0 0 0 1 1 1 1 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 1 1 1 1]))
    profiles = _trellisprofiles(stabilizers(Q), test)
    # println(profiles)
    T = _syndrometrellis(profiles, stabilizers(Q), test, false, ' ', false)
    x = Pauliweightenumerator(T)
    return x
end

function testtrellis422()
    Q = QuantumCode(["XXXX", "-ZZII", "-IIZZ"], [1, -1, -1]);
    F = field(Q)
    test = symplectictoquadratic(matrix(F,
        [1 1 0 0 0 0 0 0;
         0 0 1 1 0 0 0 0;
         0 0 0 0 1 1 0 0;
         0 0 0 0 0 1 1 0;
         0 0 0 0 0 0 1 1]))
    profiles = _trellisprofiles(stabilizers(Q), test)
    # println(profiles)
    T = _syndrometrellis(profiles, stabilizers(Q), test, false, ' ', false)
    # x = Pauliweightenumerator(T)
    x = weightenumerator(T)
    return x
end
