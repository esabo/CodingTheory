# Copyright (c) 2022, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.


#############################
         # General
#############################

# sort print on type here
function show(io::IO, W::WeightEnumerator)
    if length(W.polynomial[1]) == 3
        for (i, term) in enumerate(W.polynomial)
            if term == [1, 0, 0]
                print("1")
            end
            if term[1] != 1
                print(io, term[1])
            end
            if term[2] != 0
                if term[2] != 1
                    print(io, "x^", term[2])
                else
                    print(io, "x")
                end
            end
            if term[3] != 0
                if term[3] != 1
                    print(io, "y^", term[3])
                else
                    print(io, "y")
                end
            end
            if i != length(W.polynomial)
                print(io, " + ")
            # else
            #     print(io, "\n")
            end
        end
    else
        for (i, term) in enumerate(W.polynomial)
            if term == [1, 0, 0, 0]
                print("1")
            end
            if term[1] != 1
                print(io, term[1])
            end
            if term[2] != 0
                if term[2] != 1
                    print(io, "x^", term[2])
                else
                    print(io, "x")
                end
            end
            if term[3] != 0
                if term[3] != 1
                    print(io, "y^", term[3])
                else
                    print(io, "y")
                end
            end
            if term[4] != 0
                if term[4] != 1
                    print(io, "z^", term[4])
                else
                    print(io, "z")
                end
            end
            if i != length(W.polynomial)
                print(io, " + ")
            # else
            #     print(io, "\n")
            end
        end
    end
end

function _weightenumerator(G::fq_nmod_mat)
    E = base_ring(G)
    n = ncols(G)
    weightdist = zeros(Int, n + 1)
    for iter in Base.Iterators.product([0:(Int64(characteristic(E)) - 1) for _ in 1:nrows(G)]...)
        row = E(iter[1]) * G[1, :]
        for r in 2:nrows(G)
            if !iszero(iter[r])
                row += E(iter[r]) * G[r, :]
            end
        end
        wt = weight(row)
        weightdist[wt + 1] += 1
    end

    # weight distribution to weight enumartor
    poly = Vector{Vector{Int64}}()
    for i in 1:n + 1
        if !iszero(weightdist[i])
            push!(poly, [weightdist[i], i - 1, n - i + 1])
        end
    end
    return WeightEnumerator(poly, "Hamming")
end

#############################
        # Classical
#############################

#############################
    # Weight Enumerators
#############################

# speed boost by reducing operation on poly[j][4] at the cost of compile size
function _reducepolyC(poly::Vector{Vector{Int64}})
    reducedpoly = Vector{Vector{Int64}}()
    processed = trues(length(poly))
    for (i, term) in enumerate(poly)
        if processed[i]
            for j in (i + 1):length(poly)
                if processed[j]
                    if term[2] == poly[j][2] && term[3] == poly[j][3]
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

# for some reason this is broken?
# sorting somehow changes all the elements
function _islessPolyC(a::Vector{Int64}, b::Vector{Int64})
    if isless(a[2:3], b[2:3])
        return true
    elseif a[2:3] == b[2:3]
        if isless(a[1], b[1])
            return true
        else
            return false
        end
    else
        return false
    end
end

# set to struct
function weightenumeratorC(T::Trellis, cleanV::Bool=true)
    V = vertices(T)
    E = edges(T)
    for i in 2:length(V)
        for (j, v) in enumerate(V[i])
            outer = Vector{Vector{Int64}}()
            for e in E[i - 1][j]
                inner = deepcopy(V[i - 1][e.outvertex].polynomial)
                for k in e.label
                    if iszero(k)
                        for term in inner
                            term[3] += 1
                        end
                    else
                        for term in inner
                            term[2] += 1
                        end
                    end
                end
                append!(outer, inner)
            end
            v.polynomial = _reducepolyC(outer)
        end
    end
    W = WeightEnumerator(sort!(V[end][1].polynomial, lt=_islessPolyC), "Hamming")

    if cleanV
        for i in 2:length(V)
            for v in V[i]
                v.polynomial = [[0, 0, 0]]
            end
        end
    end
    return W
end

function weightenumerator(C::AbstractLinearCode, alg::String="trellis", sect::Bool=false)
    !ismissing(C.weightenum) && return C.weightenum
    alg ∈ ["trellis", "bruteforce"] || error("Algorithm `$alg` is not implemented in weightdistribution.")
    if alg == "trellis"
        T = syndrometrellis(C, sect)
        C.weightenum = weightenumeratorC(T, true)
    elseif alg == "bruteforce"
        C.weightenum = _weightenumerator(generatormatrix(C))
    end
    return C.weightenum
end

# MAGMA returns
# [ <0, 1>, <4, 105>, <6, 280>, <8, 435>, <10, 168>, <12, 35> ]
function weightdistribution(C::AbstractLinearCode, alg::String="trellis", sect::Bool=false)
    ismissing(C.weightenum) && weightenumerator(C, alg, sect)
    temp = zeros(Int64, 1, length(C) + 1)
    for term in C.weightenum.polynomial
        temp[term[2] + 1] = term[1]
    end
    return temp
end

#############################
     # Minimum Distance
#############################

"""
    minimumdistance(C::AbstractLinearCode, alg::String="trellis", sect::Bool=false)

Return the minimum distance of the linear code if known, otherwise computes it
using the algorithm of `alg`. If `alg = "trellis"`, the sectionalization flag
`sect` can be set to true to further compactify the reprsentation.
"""
function minimumdistance(C::AbstractLinearCode, alg::String="trellis", sect::Bool=false)
    !ismissing(C.d) && return C.d
    weightenumerator(C, alg, sect)
    # first element is always the identity y^n, want power of first x term
    C.d = C.weightenum.polynomial[2][2]
    return C.d
end

#############################
         # Quantum
#############################

#############################
    # Weight Enumerators
#############################

function _reducepolyQ(poly::Vector{Vector{Int64}})
    reducedpoly = Vector{Vector{Int64}}()
    processed = trues(length(poly))
    for (i, term) in enumerate(poly)
        if processed[i]
            for j in (i + 1):length(poly)
                if processed[j]
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

# for some reason this is broken?
# sorting somehow changes all the elements
function _islessPolyQ(a::Vector{Int64}, b::Vector{Int64})
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

# right now this only makes sense for qubit codes
# for nonqubit, can add coeff(l, ⋅) instead of 1 but then the concept of Y
# is itself fuzzy
# can just define a qudit version with X, Z tracking only

# does this actually modify T in a way I don't want to keep? should always cleanV?
function Pauliweightenumerator(T::Trellis, cleanV::Bool=true)
    V = vertices(T)
    E = edges(T)
    Int(characteristic(base_ring(E[1][1][1].label))) > 2 &&
        error("The PWE is currently only valid for qubit codes. Run the standard Hamming weight enumerator for qudit codes.")
    for i in 2:length(V)
        for (j, v) in enumerate(V[i])
            outer = Vector{Vector{Int64}}()
            for e in E[i - 1][j]
                inner = deepcopy(V[i - 1][e.outvertex].polynomial)
                for k in e.label
                    if coeff(k, 0) == 1 && iszero(coeff(k, 1))
                        if e.sign == 1
                            for term in inner
                                if term[2] < 0 || term[3] < 0 || term[4] < 0
                                    term[2] -= 1
                                else
                                    term[2] += 1
                                end
                            end
                        else
                            for term in inner
                                term[2] *= -1
                                term[3] *= -1
                                term[4] *= -1
                                if (term[2] < 0 || term[3] < 0 || term[4] < 0) || (term[2] == 0 && term[3] == 0 && term[4] == 0)
                                    term[2] -= 1
                                else
                                    term[2] += 1
                                end
                            end
                        end
                    elseif coeff(k, 0) == 1 && coeff(k, 1) == 1
                        if e.sign == 1
                            for term in inner
                                if term[2] < 0 || term[3] < 0 || term[4] < 0
                                    term[3] -= 1
                                else
                                    term[3] += 1
                                end
                            end
                        else
                            for term in inner
                                term[2] *= -1
                                term[3] *= -1
                                term[4] *= -1
                                if (term[2] < 0 || term[3] < 0 || term[4] < 0) || (term[2] == 0 && term[3] == 0 && term[4] == 0)
                                    term[3] -= 1
                                else
                                    term[3] += 1
                                end
                            end
                        end
                    elseif iszero(coeff(k, 0)) && coeff(k, 1) == 1
                        if e.sign == 1
                            for term in inner
                                if term[2] < 0 || term[3] < 0 || term[4] < 0
                                    term[4] -= 1
                                else
                                    term[4] += 1
                                end
                            end
                        else
                            for term in inner
                                term[2] *= -1
                                term[3] *= -1
                                term[4] *= -1
                                if (term[2] < 0 || term[3] < 0 || term[4] < 0) || (term[2] == 0 && term[3] == 0 && term[4] == 0)
                                    term[4] -= 1
                                else
                                    term[4] += 1
                                end
                            end
                        end
                    else
                        continue
                    end
                end
                append!(outer, inner)
            end
            v.polynomial = _reducepolyQ(outer)
        end
    end
    # sort is somehow broken
    # println(V[end][1].polynomial)
    W = WeightEnumerator(sort!(V[end][1].polynomial, lt=_islessPolyQ), "Pauli")
    # W = WeightEnumerator(V[end][1].polynomial, "Pauli")
    # println(W)

    if cleanV
        for i in 2:length(V)
            for v in V[i]
                v.polynomial = [[0, 0, 0, 0]]
            end
        end
    end
    return W
end

# function Pauliweightenumerator(Q::AbstractStabilizerCode, Pauli::Char=' ',
#     keeptrellis::Bool=true, cleanV::Bool=true, sect::Bool=false)
function Pauliweightenumerator(Q::AbstractStabilizerCode, Pauli::Char=' ',
    sect::Bool=false)

    !ismissing(Q.Pauliweightenum) && return Q.Pauliweightenum
    T = syndrometrellis(Q, "weight", Pauli, sect)
    # if keeptrellis
    #     return Pauliweightenumerator(T, cleanV), T
    # end
    # return PauliweightenumeratorQ(T, cleanV)
    return Pauliweightenumerator(T, false)
end

function PWEtoHWE(PWE::WeightEnumerator)
    poly = deepcopy(PWE.polynomial)
    n = 7
    for term in poly
        tot = abs(term[2] + term[3] + term[4])
        term[2] = tot
        term[3] = 7 - tot
        term[4] = 0
    end
    return WeightEnumerator(_reducepolyQ(poly), "Hamming")
end

function PWEtoXWE(PWE::WeightEnumerator)
    poly = deepcopy(PWE.polynomial)
    for term in poly
        tot = abs(term[2] + term[3])
        term[2] = tot
        term[3] = length(S) - tot
        term[4] = 0
    end
    return WeightEnumerator(_reducepoly(poly), "X")
end

function PWEtoZWE(PWE::WeightEnumerator)
    poly = deepcopy(PWE.polynomial)
    for term in poly
        tot = abs(term[3] + term[4])
        term[2] = length(S) - tot
        term[3] = 0
        term[4] = tot
    end
    W = WeightEnumerator(_reducepoly(poly), "Z")
    if ismissing(S.Zwtenum)
        S.Zwtenum = W
    end
    return W
end

function HammingweightenumeratorQ(T::Trellis, cleanV::Bool=true)
    V = vertices(T)
    E = edges(T)
    for i in 2:length(V)
        for (j, v) in enumerate(V[i])
            outer = Vector{Vector{Int64}}()
            for e in E[i - 1][j]
                inner = deepcopy(V[i - 1][e.outvertex].polynomial)
                for k in e.label
                    if iszero(k)
                        for term in inner
                            term[3] += 1
                        end
                    else
                        for term in inner
                            term[2] += 1
                        end
                    end
                end
                append!(outer, inner)
            end
            # println(outer)
            v.polynomial = _reducepolyQ(outer)
            # println(v.polynomial)
        end
    end
    W = WeightEnumerator(sort!(V[end][1].polynomial, lt=_islessPolyQ), "Hamming")

    if cleanV
        for i in 2:length(V)
            for v in V[i]
                v.polynomial = [[0, 0, 0, 0]]
            end
        end
    end
    return W
end

# function HammingweightenumeratorQ(Q::AbstractStabilizerCode, Paui::Char=' ',
#     keeptrellis::Bool=true, cleanV::Bool=true, sect::Bool=false)
function Hammingweightenumerator(Q::AbstractStabilizerCode, Paui::Char=' ',
    sect::Bool=false)

    !ismissing(Q.Pauliweightenum) && return PWEtoHWE(Q.Pauliweightenum)
    T = syndrometrellis(Q, "weight", Pauli, sect)
    # if keeptrellis
    #     return HammingweightenumeratorQ(T, cleanV), T
    # end
    # return HammingweightenumeratorQ(T, cleanV)
    return HammingweightenumeratorQ(T, false)
end

# function weightenumerator(Q::AbstractStabilizerCode, alg::String="trellis", sect::Bool=false,
#     Pauli::Char=' ', keeptrellis::Bool=true, cleanV::Bool=true, verbose::Bool=false)


# we can organize the S⟂ \ S any way we want as long as we have these elements on the
function weightenumerator(Q::AbstractStabilizerCode, alg::String="trellis", sect::Bool=false,
    Pauli::Char=' ', verbose::Bool=true)

    alg ∈ ["trellis"] || error("Algorithm `$alg` is not implemented in weightdistribution.")

    # S⟂ - this is not a stabilizer code
    # if ismissing(Q.dualweightenum)
    #     Q.dualweightenum = weightenumerator(LinearCode(normalizermatrix(Q)))
    # end
    verbose && println("S⟂: ", Q.dualweightenum)

    # S⟂ \ S - this is not a stabilizer code
    # k is always small, so just brute force it
    if ismissing(Q.logsweightenum)
        logspace = LinearCode(logicalspace(Q))
        Q.logsweightenum = _weightenumerator(logspace)
        Q.d = Q.logsweightenum.polynomial[2][2]
    end
    verbose && println("S⟂ / S: ", Q.logsweightenum)

    # S
    if ismissing(Q.Pauliweightenum)
        # this is not valid for qudit codes
        Q.Pauliweightenum = Pauliweightenumerator(Q, Pauli, sect)
    end
    verbose && println("S: ", Q.Pauliweightenum)

    return [Q.dualweightenum, Q.logsweightenum, Q.Pauliweightenum]
end

function weightdistribution(Q::AbstractStabilizerCode, alg::String="trellis", sect::Bool=false,
    Pauli::Char=' ', keeptrellis::Bool=true, cleanV::Bool=true, verbose::Bool=true)

    weightenumerator(Q, alg, sect, Pauli, false)
    # S⟂
    temp1 = zeros(Int64, 1, length(Q) + 1)
    for term in Q.dualweightenum
        temp1[term[2] + 1] = term[1]
    end
    verbose && println("S⟂: ", temp1)

    # S⟂ \  S
    temp2 = zeros(Int64, 1, length(Q) + 1)
    for term in Q.logsweightenum
        temp2[term[2] + 1] = term[1]
    end
    verbose && println("S⟂ / S: ", temp2)

    # S
    temp3 = zeros(Int64, 1, length(Q) + 1)
    for term in Q.Pauliweightenum
        temp3[term[2] + 1] = term[1]
    end
    verbose && println("S: ", temp3)

    return [temp1, temp2, temp3]
end

# function weightdistribution(Q::AbstractStabilizerCode, alg::String="trellis", sect::Bool=false)
#     alg ∈ ["trellis"] || error("Algorithm `$alg` is not implemented in weightdistribution.")
#     T = syndrometrellis(Q, "weight", ' ', sect)
#     W = PauliweightenumeratorQ(T, true)
#     H = PWEtoHWE(W)
#     temp = zeros(Int64, 1, length(Q))
#     for term in H
#         temp[term[2]] = term[1]
#     end
#
#     if typeof(Q) <: CSSCode
#         X = PWEtoXWE(W)
#         Z = PWEtoZWE(W)
#         dist = Vector{Vector{Int64}}()
#         push!(dist, temp)
#
#         temp = zeros(Int64, 1, length(Q))
#         for term in X
#             temp[term[2] + 1] = term[1]
#         end
#         push!(dist, temp)
#
#         temp = zeros(Int64, 1, length(Q))
#         for term in Z
#             temp[term[4] + 1] = term[1]
#         end
#         push!(dist, temp)
#         return dist
#     end
#     return temp
# end

#############################
     # Minimum Distance
#############################

"""
    minimumdistance(Q::AbstractStabilizerCode, alg::String="trellis", sect::Bool=false)

Return the minimum distance of the stabilizer code if known, otherwise computes it.

"""
function minimumdistance(Q::AbstractStabilizerCode)
    !ismissing(Q.d) && return Q.d
    if ismissing(Q.logsweightenum)
        logspace = logicalspace(Q)
        Q.logsweightenum = _weightenumerator(logspace)
    end
    Q.d = Q.logsweightenum.polynomial[2][2]
    return Q.d
end

# function minimumdistance(Q::AbstractCSSCode, alg::String="trellis", sect::Bool=false)
#
#
#
#
#
#
#     !ismissing(Q.d) && return Q.d
#     if ismissing(Q.logsweightenum)
#         Q.logsweightenum = weightenumerator(LinearCode(logicalspace(Q)))
#     end
#     # first element is always the identity y^n, want power of first x term
#     Q.d = Q.logsweightenum.polynomial[2][2]
#     return Q.d
# end

minimumdistance(S::CSSCode) = S.d
minimumdistanceX(S::CSSCode) = S.dx
minimumdistanceZ(S::CSSCode) = S.dz
