# Copyright (c) 2022, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
         # General
#############################

"""
    polynomial(W::WeightEnumerator)

Returns the polynomial of the weight enumerator `W`.
"""
polynomial(W::WeightEnumerator) = W.polynomial

"""
    type(W::WeightEnumerator)

Returns the type of the weight enumerator `W`.
"""
type(W::WeightEnumerator) = W.type

function ==(W1::WeightEnumerator, W2::WeightEnumerator)
    W1.type == W2.type || return false
    W1.polynomial == W2.polynomial && return true
    return false
end

# TODO: test with other iterator
function _weightenumeratorBF(G::fq_nmod_mat)
    E = base_ring(G)
    ordE = Int(order(E))
    lookup = Dict(value => key for (key, value) in enumerate(collect(E)))
    R, vars = PolynomialRing(Nemo.ZZ, ordE)
    poly = R(0)
    nr = nrows(G)

    # Nemo.AbstractAlgebra.ProductIterator
    for iter in Base.Iterators.product([0:(Int64(characteristic(E)) - 1) for _ in 1:nr]...)
        row = E(iter[1]) * G[1, :]
        for r in 2:nr
            if !iszero(iter[r])
                row += E(iter[r]) * G[r, :]
            end
        end

        # TODO: can do this in one step, but is it faster?
        term = zeros(Int, 1, ordE)
        for x in row
            term[lookup[x]] += 1
        end
        termpoly = R(1)
        for i in 1:ordE
            termpoly *= vars[i]^term[i]
        end
        poly += termpoly
    end
    return WeightEnumerator(poly, "complete")
end

#############################
        # Classical
#############################

#############################
    # Weight Enumerators
#############################

# make private?
function CWEtoHWE(CWE::WeightEnumerator)
    R, (x, y) = PolynomialRing(Nemo.ZZ, ["x", "y"])
    poly = R(0)
    for i in 1:length(CWE.polynomial)
        exps = exponent_vector(CWE.polynomial, i)
        poly += coeff(CWE.polynomial, i) * x^sum(exps[2:end]) * y^exps[1]
    end
    return WeightEnumerator(poly, "Hamming")
end

function weightenumeratorC(T::Trellis, type::String="complete")
    type ∈ ["complete", "Hamming"] || error("Unsupported weight enumerator type '$type'. Expected 'complete' or 'Hamming'.")

    if type == "complete" && !ismissing(T.CWE)
        return T.CWE
    elseif type == "Hamming" && !ismissing(T.CWE)
        return CWEtoHWE(T.CWE)
    end

    # if this ever changes or permutes will have to store with T
    elms = collect(field(T.code))
    lookup = Dict(value => key for (key, value) in enumerate(elms))
    R, vars = PolynomialRing(Nemo.ZZ, length(elms))

    V = T.vertices
    E = T.edges
    V[1][1].polynomial = R(1)
    # V[1][1].polynomial[1][1] = 1
    for i in 2:length(V)
        for (j, v) in enumerate(V[i])
            outer = R(0)
            for e in E[i - 1][j]
                inner = deepcopy(V[i - 1][e.outvertex].polynomial)
                for k in e.label
                    inner *= vars[lookup[k]]
                end
                outer += inner
            end
            v.polynomial = outer
        end
    end
    T.CWE = WeightEnumerator(V[end][1].polynomial, "complete")

    # currently Missing is not an option but how to implement dual trellis
    if !isshifted(T) && !ismissing(T.code)
        T.code.weightenum = T.CWE
        HWE = CWEtoHWE(T.CWE)
        T.code.d = minimum([collect(exponent_vectors(polynomial(HWE)))[i][1]
            for i in 1:length(polynomial(HWE))])
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

# ah, there's nothing about this function which is correct. I've been changing the
# exponents while I need to evaluate the functions then expand using the same exponents
function MacWilliamsIdentity(C::AbstractLinearCode, W::WeightEnumerator, dual::String="Euclidean")
    dual ∈ ["Euclidean", "Hermitian"] ||
        error("The MacWilliams identities are only programmed for the Euclidean and Hermitian duals.")
    (dual == "Hermitian" && Int(order(C.F)) != 4) &&
        error("The MacWilliams identity for the Hermitian dual is only programmed for GF(4).")

    if W.type == "Hamming"
        # (1/|C|)W(y - x, y + (q - 1)x)
        R = parent(W.polynomial)
        vars = gens(R)
        return WeightEnumerator(divexact(W.polynomial(vars[2] - vars[1], vars[2] +
            (Int(order(C.F)) - 1) * vars[1]), cardinality(C)), "Hamming")
    end

    # complete weight enumerators
    if Int(order(C.F)) == 2
        # the complete and Hamming weight enumerators are the same in binary
        # (1/|C|)W(x_0 + (q - 1)x_1, x_0 - x_1)
        R = parent(W.polynomial)
        vars = gens(R)
        return WeightEnumerator(divexact(W.polynomial(vars[1] +
            (Int(order(C.F)) - 1)*vars[2], vars[1] - vars[2]),
            cardinality(C)), "complete")
    elseif Int(order(C.F)) == 3
        # (1/|C|)W(x_0 + x_1 + x_2, x_0 + ω x_1 + ω^2 x_2, x_0 + ω^2 x_1 + ω x_2)
        K, ζ = CyclotomicField(3, "ζ")
        R, vars = PolynomialRing(K, 3)
        poly = divexact(W.polynomial(
            vars[1] + vars[2] + vars[3],
            vars[1] + ζ*vars[2] + ζ^2*vars[3],
            vars[1] + ζ^2*vars[2] + ζ*vars[3]), cardinality(C))
        # works so far but now needs to recast down to the integer ring
        return WeightEnumerator(map_coefficients(c -> Nemo.ZZ(coeff(c, 0)), poly,
            parent=parent(W.polynomial)), "complete")
    elseif Int(order(C.F)) == 4
        if dual == "Euclidean"
            # for Euclidean dual
            # (1/|C|)W(x_0 + x_1 + x_2 + x_3, x_0 + x_1 - x_2 - x_3, x_0 - x_1 - x_2 + x_3, x_0 - x_1 + x_2 - x_3)
            R = parent(W.polynomial)
            vars = gens(R)
            return WeightEnumerator(divexact(W.polynomial(
                vars[1] + vars[2] + vars[3] + vars[4],
                vars[1] + vars[2] - vars[3] - vars[4],
                vars[1] - vars[2] - vars[3] + vars[4],
                vars[1] - vars[2] + vars[3] - vars[4]), cardinality(C)),
                "complete")
        else
            # for Hermitian dual
            # (1/|C|)W(x_0 + x_1 + x_2 + x_3, x_0 + x_1 - x_2 - x_3, x_0 - x_1 + x_2 - x_3, x_0 - x_1 - x_2 + x_3)
            R = parent(W.polynomial)
            vars = gens(R)
            return WeightEnumerator(divexact(W.polynomial(
                vars[1] + vars[2] + vars[3] + vars[4],
                vars[1] + vars[2] - vars[3] - vars[4],
                vars[1] - vars[2] + vars[3] - vars[4],
                vars[1] - vars[2] - vars[3] + vars[4]), cardinality(C)),
                "complete")
        end
    else
        # do the full manual thing here
        error("The MacWilliams identities for complete weight enumerators have
            not yet been implemented over higher fields. See MachWilliams and
            Sloane for the proper implementation.")
    end
end

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
        if cardinality(C) <= 1e6 # random cutoff
            C.weightenum = _weightenumeratorBF(C.G)
            HWE = CWEtoHWE(C.weightenum)
            C.d = minimum(filter(x->x!=0, [collect(exponent_vectors(HWE.polynomial))[i][1]
                for i in 1:length(HWE.polynomial)]))
            type == "Hamming" && return HWE
            return C.weightenum
        elseif rate(C) > 0.5
            D = dual(C)
            if cardinality(D) <= 1e6 # random cutoff
                D.weightenum = _weightenumeratorBF(D.G)
            else
                weightenumeratorC(syndrometrellis(D, "primal", false), type)
            end
            C.weightenum = MacWilliamsIdentity(D, D.weightenum)
            HWE = CWEtoHWE(C.weightenum)
            C.d = minimum(filter(x->x!=0, [collect(exponent_vectors(HWE.polynomial))[i][1]
                for i in 1:length(HWE.polynomial)]))
            type == "Hamming" && return HWE
            return C.weightenum
        else
            return weightenumeratorC(syndrometrellis(C, "primal", false), type)
        end
    elseif alg == "trellis"
        return weightenumeratorC(syndrometrellis(C, "primal", false), type)
    elseif alg == "bruteforce"
        C.weightenum = _weightenumeratorBF(C.G)
        HWE = CWEtoHWE(C.weightenum)
        C.d = minimum(filter(x->x!=0, [collect(exponent_vectors(HWE.polynomial))[i][1]
            for i in 1:length(HWE.polynomial)]))
        type == "Hamming" && return HWE
        return C.weightenum
    end
end

# MAGMA returns this format
# [ <0, 1>, <4, 105>, <6, 280>, <8, 435>, <10, 168>, <12, 35> ]
function weightdistribution(C::AbstractLinearCode, alg::String="auto", format::String="full")
    alg ∈ ["auto", "trellis", "bruteforce"] || error("Algorithm `$alg` is not implemented in weightenumerator.")
    format ∈ ["full", "compact"] || error("Unknown value for parameter format: $format; expected `full` or `compact`.")

    ismissing(C.weightenum) && weightenumerator(C, "complete", alg)
    HWE = CWEtoHWE(C.weightenum)

    if format == "compact"
        wtdist = Vector{Tuple}()
        for i in 1:length(HWE.polynomial)
            push!(wtdist, (exponent_vector(HWE.polynomial, i)[1],
                coeff(HWE.polynomial, i)))
        end
    else
        wtdist = zeros(Int, 1, C.n + 1)
        for i in 1:length(HWE.polynomial)
            wtdist[exponent_vector(HWE.polynomial, i)[1]] = coeff(HWE.polynomial, i)
        end
    end
    return wtdist
end

"""
    support(C::AbstractLinearCode)

Returns the support of `C`.

The support of `C` is the collection of nonzero exponents of the Hamming
weight enumerator of `C`.
"""
support(C::AbstractLinearCode) = [i for (i, _) in weightdistribution(C, "auto", "compact")]

#############################
     # Minimum Distance
#############################

"""
    minimumdistance(C::AbstractLinearCode, alg::String="trellis", sect::Bool=false)

Return the minimum distance of the linear code if known, otherwise computes it
using the algorithm of `alg`. If `alg = "trellis"`, the sectionalization flag
`sect` can be set to true to further compactify the reprsentation.
"""
function minimumdistance(C::AbstractLinearCode, alg::String="auto")
    !ismissing(C.d) && return C.d
    HWE = weightenumerator(C, "Hamming", alg)
    !ismissing(C.d) && return C.d
    # this line should only be needed to be run if the weight enumerator is known
    # but the minimum distance is intentionally set to missing
    # ordering here can be a bit weird
    C.d = minimum(filter(x->x!=0, [collect(exponent_vectors(HWE.polynomial))[i][1]
        for i in 1:length(HWE.polynomial)]))
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
    V = T.vertices
    E = T.edges
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
    V = T.vertices
    E = T.edges
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
    temp1 = zeros(Int64, 1, Q.n + 1)
    for term in Q.dualweightenum
        temp1[term[2] + 1] = term[1]
    end
    verbose && println("S⟂: ", temp1)

    # S⟂ \  S
    temp2 = zeros(Int64, 1, Q.n + 1)
    for term in Q.logsweightenum
        temp2[term[2] + 1] = term[1]
    end
    verbose && println("S⟂ / S: ", temp2)

    # S
    temp3 = zeros(Int64, 1, Q.n + 1)
    for term in Q.Pauliweightenum
        temp3[term[2] + 1] = term[1]
    end
    verbose && println("S: ", temp3)

    return [temp1, temp2, temp3]
end

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

minimumdistance(S::CSSCode) = S.d
minimumdistanceX(S::CSSCode) = S.dx
minimumdistanceZ(S::CSSCode) = S.dz
