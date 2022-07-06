# Copyright (c) 2022, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
         # General
#############################

function _islessLex(a::Vector{Int64}, b::Vector{Int64})
    if a[2:end] == b[2:end] && a[1] < b[1]
        return true
    elseif a[2:end] == b[2:end] && a[1] > b[1]
        return false
    end
    for i = 2:length(a)
        if a[i] < b[i]
            return true
        elseif a[i] > b[i]
            return false
        end
    end
    return false
end

function _reducepoly(poly::Vector{Vector{Int64}})
    reducedpoly = Vector{Vector{Int64}}()
    processed = trues(length(poly))
    for (i, term) in enumerate(poly)
        if processed[i]
            for j in (i + 1):length(poly)
                if processed[j]
                    if term[2:end] == poly[j][2:end]
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

# things should enter this reduced and sorted, so no need to do this here? -
function -(W1::WeightEnumerator, W2::WeightEnumerator)
    W1.type == W2.type || error("Cannot subtract two weight enumerators of different types.")
    polyW1 = W1.polynomial
    n = length(W1.polynomial[1])
    n == length(W2.polynomial[1]) || error("Cannot subtract two weight enumerators with a different number of variables.")

    n -= 1
    newpoly = Vector{Vector{Int}}()
    for termW1 in polyW1
        flag = false
        for termW2 in W2.polynomial
            if termW1[2:n] == termW2[2:n]
                push!(newpoly, [termW1[1] - termW2[1]; termW1[2:n]])
                flag = true
                break
            end
        end
        if !flag
            push!(newpoly, termW1)
        end
    end
    return WeightEnumerator(newpoly, W1.type)
end

function ==(W1::WeightEnumerator, W2::WeightEnumerator)
    W1.type == W2.type || return false
    a = WeightEnumerator(sort!(_reducepoly(W1.polynomial), lt=_islessLex), W1.type)
    b = WeightEnumerator(sort!(_reducepoly(W2.polynomial), lt=_islessLex), W2.type)
    length(a.polynomial) == length(b.polynomial) || return false
    # for i in 1:length(a.polynomial)
    #     a.polynomial[i] == b.polynomial[i] || return false
    # end
    # return true
    a.polynomial == b.polynomial && return true
    return false
end

"""
prints in order of collect(F) where F is the field of the code

"""
function show(io::IO, W::WeightEnumerator)
    poly = W.polynomial
    if W.type == "complete"
        # println(io, "A complete weight enumerator for a classical code.")
        len = length(poly[1])
        for (i, term) in enumerate(poly)
            if !iszero(term[1]) && iszero(term[2:end])
                print(io, term[1])
            else
                if term[1] != 1
                    print(io, term[1])
                end
                for j in 2:len
                    if !iszero(term[j])
                        # TODO: fix this
                        # print(io, "x$(Base.REPLCompletions.latex_symbols["\\_$(j - 1)"])",
                        #     "$(Base.REPLCompletions.latex_symbols["\\^$term[j]"])")
                        print(io, "x_$(j - 1)^$(term[j])")
                    end
                end
            end
            if i != length(poly)
                print(io, " + ")
            end
        end
    elseif W.type == "Hamming"
        # println(io, "A Hamming weight enumerator for a classical code.")
        for (i, term) in enumerate(poly)
            if !iszero(term[1]) && iszero(term[2:end])
                print(io, term[1])
            else
                if term[1] != 1
                    print(io, term[1])
                end
                if !iszero(term[2])
                    print(io, "x^", term[2])
                end
                if !iszero(term[3])
                    print(io, "y^", term[3])
                end
            end
            if i != length(poly)
                print(io, " + ")
            end
        end
    end
end

# TODO: test with other iterator
function _weightenumeratorBF(G::fq_nmod_mat)
    E = base_ring(G)
    ordE = Int(order(E))
    lookup = Dict(value => key for (key, value) in enumerate(collect(E)))
    poly = Vector{Vector{Int}}()
    # Nemo.AbstractAlgebra.ProductIterator
    for iter in Base.Iterators.product([0:(Int64(characteristic(E)) - 1) for _ in 1:nrows(G)]...)
        row = E(iter[1]) * G[1, :]
        for r in 2:nrows(G)
            if !iszero(iter[r])
                row += E(iter[r]) * G[r, :]
            end
        end

        term = zeros(Int, 1, ordE + 1)
        for x in row
            term[lookup[x]] += 1
        end

        # this is slow but keeps the memory down of reducing afterwares
        flag = true
        for x in poly
            if term == x[2:end]
                x[1] += 1
                flag = false
                break
            end
        end
        flag && push!(poly, vcat([1], vec(term)))
    end

    return WeightEnumerator(sort!(_reducepoly(poly), lt=_islessLex), "complete")
end

#############################
        # Classical
#############################

#############################
    # Weight Enumerators
#############################

# make private?
function CWEtoHWE(CWE::WeightEnumerator)
    poly = Vector{Vector{Int64}}()
    # is homogeneous
    n = sum(CWE.polynomial[1][2:end])
    for term in CWE.polynomial
        # absolute value for quantum, signed CWEs?
        # tot = sum(term[3:end])
        # CWE.n >= tot || error("Total degree of term in complete weight enumerator is larger than the code length.")
        # push!(poly, [term[1], tot, n - tot])

        # actually can just do this since it's homogeneous
        push!(poly, [term[1], n - term[2], term[2]])
    end
    return WeightEnumerator(sort!(_reducepoly(poly), lt=_islessLex), "Hamming")
end

function weightenumeratorC(T::Trellis, type::String="complete")
    type ∈ ["complete", "Hamming"] || error("Unsupported weight enumerator type '$type'. Expected 'complete' or 'Hamming'.")

    if type == "complete" && !ismissing(T.CWE)
        return T.CWE
    elseif type == "Hamming" && !ismissing(T.CWE)
        return CWEtoHWE(T.CWE)
    end

    # if this ever changes or permutes will have to store with T
    # TODO: fix this
    # elms = collect(field(T.code))
    elms = collect(T.code.F)
    V = vertices(T)
    E = edges(T)
    V[1][1].polynomial = [zeros(Int, length(elms) + 1)]
    V[1][1].polynomial[1][1] = 1
    for i in 2:length(V)
        for (j, v) in enumerate(V[i])
            outer = Vector{Vector{Int64}}()
            for e in E[i - 1][j]
                inner = deepcopy(V[i - 1][e.outvertex].polynomial)
                for k in e.label
                    # probably a better way to do this
                    # TODO: dictionary it for faster lookup?
                    for term in inner
                        term[findfirst(x->x==k, elms) + 1] += 1
                    end
                end
                append!(outer, inner)
            end
            v.polynomial = _reducepoly(outer)
        end
    end
    T.CWE = WeightEnumerator(sort!(V[end][1].polynomial, lt=_islessLex),
        "complete")

    # currently Missing is not an option but how to implement dual trellis
    if !isshifted(T) && !ismissing(T.code)
        T.code.weightenum = T.CWE
    end

    # clean up vertices
    for i in 1:length(V)
        for v in V[i]
            v.polynomial = missing
        end
    end

    if type == "Hamming"
        return CWEtoHWE(T.CWE)
    end
    return T.CWE
end

# TODO: MacWilliams identities, check size of dual, if small just bruteforce
function weightenumerator(C::AbstractLinearCode, type::String="complete",
    alg::String="auto")

    type ∈ ["complete", "Hamming"] || error("Unsupported weight enumerator type '$type'. Expected 'complete' or 'Hamming'.")
    alg ∈ ["auto", "trellis", "bruteforce"] || error("Algorithm `$alg` is not implemented in weightenumerator.")

    if type == "complete" && !ismissing(C.weightenum)
        return C.weightenum
    elseif type == "Hamming" && !ismissing(C.weightenum)
        return CWEtoHWE(C.weightenum)
    end

    if alg == "auto"
        if BigInt(characteristic(LinearCodeMod.field(C)))^LinearCodeMod.dimension(C) <= 1e6 # random cutoff
            C.weightenum = _weightenumeratorBF(generatormatrix(C))
            type == "Hamming" && return CWEtoHWE(C.weightenum)
            return C.weightenum
        else
            return weightenumeratorC(syndrometrellis(C, "primal", false), type)
        end
    elseif alg == "trellis"
        return weightenumeratorC(syndrometrellis(C, "primal", false), type)
    elseif alg == "bruteforce"
        C.weightenum = _weightenumeratorBF(generatormatrix(C))
        type == "Hamming" && return CWEtoHWE(C.weightenum)
        return C.weightenum
    end
end

# MAGMA returns this format
# [ <0, 1>, <4, 105>, <6, 280>, <8, 435>, <10, 168>, <12, 35> ]
function weightdistribution(C::AbstractLinearCode, alg::String="auto", format::String="full")
    alg ∈ ["auto", "trellis", "bruteforce"] || error("Algorithm `$alg` is not implemented in weightenumerator.")
    format ∈ ["full", "compact"] || error("Unknown value for parameter format: $format; expected `full` or `compact`.")

    ismissing(C.weightenum) && weightenumerator(C, "trellis")
    wtdist = zeros(Int64, 1, length(C) + 1)
    for term in C.weightenum.polynomial
        wtdist[term[2] + 1] = term[1]
    end

    if format == "compact"
        compact = Vector{Tuple}()
        for (i, x) in enumerate(wtdist)
            !iszero(x) && push!(compact, (i, x))
        end
        return compact
    end
    return wtdist
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
function minimumdistance(C::AbstractLinearCode, alg::String="trellis")
    !ismissing(C.d) && return C.d
    weightenumerator(C, alg)
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
            v.polynomial = _reducepoly(outer)
        end
    end
    # sort is somehow broken
    # println(V[end][1].polynomial)
    W = WeightEnumerator(sort!(V[end][1].polynomial, lt=_islessLex), "Pauli")
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
    return WeightEnumerator(_reducepoly(poly), "Hamming")
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
            v.polynomial = _reducepoly(outer)
            # println(v.polynomial)
        end
    end
    W = WeightEnumerator(sort!(V[end][1].polynomial, lt=_islessLex), "Hamming")

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
        Q.logsweightenum = _weightenumeratorBF(logspace)
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
        Q.logsweightenum = _weightenumeratorBF(logspace)
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
