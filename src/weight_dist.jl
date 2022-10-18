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

# make private?
function CWEtoHWE(CWE::WeightEnumerator)
    R, (x, y) = PolynomialRing(base_ring(CWE.polynomial), ["x", "y"])
    poly = R(0)
    for i in 1:length(CWE.polynomial)
        exps = exponent_vector(CWE.polynomial, i)
        poly += coeff(CWE.polynomial, i) * x^sum(exps[2:end]) * y^exps[1]
    end
    return WeightEnumerator(poly, "Hamming")
end

#############################
        # Classical
#############################

#############################
    # Weight Enumerators
#############################

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
            Threads.@threads for e in E[i - 1][j]
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
        # might have to switch this here
        poly = divexact(W.polynomial(
            vars[1] + vars[2] + vars[3],
            vars[1] + ζ*vars[2] + ζ^2*vars[3],
            vars[1] + ζ^2*vars[2] + ζ*vars[3]), cardinality(C))
        # works so far but now needs to recast down to the integer ring
        return WeightEnumerator(map_coefficients(c -> Nemo.ZZ(coeff(c, 0)), poly,
            parent=parent(W.polynomial)), "complete")
    elseif Int(order(C.F)) == 4
        # these order 4 formulas are from "Self-Dual Codes" by Rains and Sloane without proof
        # the differ in order from the formula in MacWilliams and Sloane used in the general
        # case below:
        #    x1 + x2 + x3 + x4
        #    x1 - x2 + x3 - x4
        #    x1 + x2 - x3 - x4
        #    x1 - x2 - x3 + x4
        # But that formula should depend on the chosen basis and character so I assume it's okay
        if dual == "Euclidean"
            # for Euclidean dual
            # (1/|C|)W(x_0 + x_1 + x_2 + x_3, x_0 + x_1 - x_2 - x_3, x_0 - x_1 - x_2 + x_3, x_0 - x_1 + x_2 - x_3)
            R = parent(W.polynomial)
            vars = gens(R)
            # switched lines 2 and 3 from Rains & Sloane (Huffman & Press) formula because it
            # appears to implicitly assuming a primitive basis and here we permute for our basis
            return WeightEnumerator(divexact(W.polynomial(
                vars[1] + vars[2] + vars[3] + vars[4],
                vars[1] - vars[2] - vars[3] + vars[4],
                vars[1] + vars[2] - vars[3] - vars[4],
                vars[1] - vars[2] + vars[3] - vars[4]), cardinality(C)),
                "complete")
        else
            # for Hermitian dual
            # (1/|C|)W(x_0 + x_1 + x_2 + x_3, x_0 + x_1 - x_2 - x_3, x_0 - x_1 + x_2 - x_3, x_0 - x_1 - x_2 + x_3)
            R = parent(W.polynomial)
            vars = gens(R)
            # switched lines 2 and 3 from Rains & Sloane (Huffman & Press) formula because it
            # appears to implicitly assuming a primitive basis and here we permute for our basis
            return WeightEnumerator(divexact(W.polynomial(
                vars[1] + vars[2] + vars[3] + vars[4],
                vars[1] - vars[2] + vars[3] - vars[4],
                vars[1] + vars[2] - vars[3] - vars[4],
                vars[1] - vars[2] - vars[3] + vars[4]), cardinality(C)),
                "complete")
        end
    else
        q = Int(order(C.F))
        if isprime(q)
            K, ω = CyclotomicField(Int(characteristic(C.F)), "ω")
            R, vars = PolynomialRing(K, q)
            elms = collect(C.F)
            funcargs = []
            for i in 1:q
                innersum = R(0)
                for j in 1:q
                    innersum += ω^coeff(elms[i] * elms[j], 0) * vars[j]
                end
                append!(funcargs, innersum)
            end
            return WeightEnumerator(divexact(W.polynomial(funcargs), cardinality(C)),
                "complete")
        else
            K, ω = CyclotomicField(Int(characteristic(C.F)), "ω")
            R, vars = PolynomialRing(K, q)
            primefield, _ = FiniteField(Int(characteristic(C.F)), 1, "α2")
            _, λ = primitivebasis(C.F, primefield)
            elms = collect(C.F)
            funcargs = []
            for i in 1:q
                innersum = R(0)
                for j in 1:q
                    β = elms[i] * elms[j]
                    βexp = _expandelement(β, primefield, λ, false)
                    innersum += ω^coeff(βexp[1], 0) * vars[j]
                end
                push!(funcargs, innersum)
            end
            display(funcargs)
            return WeightEnumerator(divexact(W.polynomial(funcargs...), cardinality(C)),
                "complete")
        end
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
            wtdist[exponent_vector(HWE.polynomial, i)[1] + 1] = coeff(HWE.polynomial, i)
        end
    end
    return wtdist
end

"""
    weighthplot(C::AbstractLinearCode, alg::String="auto")

Return a bar plot of the weight distribution of `C`.
"""
function weightplot(C::AbstractLinearCode, alg::String="auto")
    wtdist = weightdistribution(C, alg, "full")
    xticks = findall(x->x>0, vec(wtdist)) .- 1
    yticks = [wtdist[i] for i in 1:length(wtdist) if !iszero(wtdist[i])]
    ismissing(C.d) ? (title="Weight Distribution - [$(C.n), $(C.k)]";) :
        title="Weight Distribution - [$(C.n), $(C.k), $(C.d)]"
    f = bar(0:C.n, wtdist', bar_width=1, xticks=xticks, yticks=yticks,
        legend=false, xlabel="Weight", ylabel="Number of Terms", title=title)
    show(f)
    return f
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

# TODO: update doc strings for these and this whole file in general
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

# TODO: test with other iterator
function _weightenumeratorBFQ(G::fq_nmod_mat, charvec::Vector{nmod},
    R::Union{AbstractAlgebra.Generic.MPolyRing{nf_elem}, Missing})
    # this should be the quadratic extension field
    E = base_ring(G)
    iseven(Int(degree(E))) || error("Matrix passed to weight enumerator does not appear to be over the quadratic extension.")
    ordE = Int(order(E))
    lookup = Dict(value => key for (key, value) in enumerate(collect(E)))
    
    p = Int(characteristic(E))
    iseven(p) ? nth = 2 * p : nth = p
    if ismissing(R)
        K, ω = CyclotomicField(nth, "ω")
        R, vars = PolynomialRing(K, ordE)
    else
        ω = gen(base_ring(R))
        vars = gens(R)
    end
    poly = R(0)
    nr, nc = size(G)

    # Nemo.AbstractAlgebra.ProductIterator
    for iter in Base.Iterators.product([0:(p - 1) for _ in 1:nr]...)
        row = E(iter[1]) * G[1, :]
        for r in 2:nr
            if !iszero(iter[r])
                row += E(iter[r]) * G[r, :]
            end
        end
        rowsym = quadratictosymplectic(row)

        # to do process signs here
        parity = 0
        for c in 1:2 * nc
            iszero(rowsym[c]) || (parity += data(charvec[c]);)
        end

        # TODO: can do this in one step, but is it faster?
        term = zeros(Int, 1, ordE)
        for x in row
            term[lookup[x]] += 1
        end
        # println(term, ", ", typeof(term))
        termpoly = ω^parity
        for i in 1:ordE
            termpoly *= vars[i]^term[i]
        end
        poly += termpoly
    end
    # display(poly)
    return WeightEnumerator(poly, "complete")
    # return poly
end

# formulas from
# "Weight enumerators for nonbinary asymmetric quantum codes and their applications"
# by Chuangqiang Hu, Shudi Yang, Stephen S.-T.Yau
function MacWilliamsIdentity(S::AbstractStabilizerCode, W::WeightEnumerator, dual::Bool=false)
    dual ? (card = BigInt(characteristic(S.F))^(S.n + S.k);) : (card = cardinality(S);)
    if W.type == "Hamming"
        # (1/(q^n|S|))W(y - x, y + (q^2 - 1)x)
        R = parent(W.polynomial)
        vars = gens(R)
        q = Int(order(S.F))
        return WeightEnumerator(divexact(W.polynomial(vars[2] - vars[1], vars[2] +
            (q^2 - 1) * vars[1]), card), "Hamming")
        # could probably put the /q under each variable and remove the q^n
    end

    # complete weight enumerators
    if Int(order(S.E)) == 4
        # 1/|S| W((x - y - z + w)/2, (-x + y - z + w)/2, (-x - y + z + w)/2, (x + y + z + w)/2)
        R = parent(W.polynomial)
        vars = gens(R)
        # this is the same as the classical Hermitian dual formula
        # switched lines 2 and 3 from citation for our basis
        return WeightEnumerator(divexact(W.polynomial(
            vars[1] + vars[2] + vars[3] + vars[4],
            vars[1] + vars[2] - vars[3] - vars[4],
            vars[1] - vars[2] + vars[3] - vars[4],
            vars[1] - vars[2] - vars[3] + vars[4]),
            card), "complete") # need the /2 to connect to the original Shor-Laflamme def
    else
        error("The quantum MacWilliams identity for higher fields has a bug and is currently unavailable.")
        # BUG: in the below it's unclear what the proper permutation is given the paper
        # the various combinations I've tried always fix one but break the dual
        # need to set ω ↦ ω^2 and then match the equations above (try Q15RM())
        # want perm = [1, 3, 2, 4]
        # R = parent(W.polynomial)
        # vars = gens(R)
        # ω = gen(base_ring(R)) # if Int(order(S.F)) == 2, ω ↦ ω^2 in below
        # elms = collect(S.E)
        # q = Int(order(S.E))
        # α = gen(S.E)
        # basis = [S.E(0); [α^i for i in 1:q - 1]]
        # # perm = [findfirst(x->x==b, elms) for b in basis]
        # perm = [findfirst(x->x==b, basis) for b in elms]
        # funcargs = []
        # for i in 1:q
        #     innersum = R(0)
        #     for j in 1:q
        #         # innersum += ω^(2*tr(coeff(elms[i], 0) * coeff(elms[j], 1) -
        #         #     coeff(elms[j], 0) * coeff(elms[i], 1))) * vars[j]
        #         innersum += ω^(2*tr(coeff(basis[i], 0) * coeff(basis[j], 1) -
        #             coeff(basis[j], 0) * coeff(basis[i], 1))) * vars[perm[j]]
        #     end
        #     push!(funcargs, innersum) # /q for Shor-Laflamme
        # end
        # println(basis)
        # println(elms)
        # println(perm)
        # display(funcargs)
        # display(funcargs[perm])
        # return WeightEnumerator(divexact(W.polynomial(funcargs[perm]...), card),
        #     "complete")
    end
end

function weightenumerator(S::AbstractStabilizerCode, type::String="complete",
    alg::String="auto", set::String="all")

    type ∈ ["complete", "Hamming"] || error("Unsupported weight enumerator type '$type'. Expected 'complete' or 'Hamming'.")
    alg ∈ ["auto", "trellis", "bruteforce"] || error("Algorithm `$alg` is not implemented in weightenumerator.")
    set ∈ ["all", "stabilizers", "logicals", "quotient"] || throw(ArgumentError("Unsupported set type '$set'. Expected 'all', 'stabilizers', 'logicals', 'quotient'."))

    if set ∈ ["all", "logicals"] && ismissing(S.sCWElogs)
        logsmat = logicalsmatrix(S)
        S.sCWElogs = _weightenumeratorBFQ(logsmat, S.charvec, missing)
    end

    if set != "logicals" && ismissing(S.sCWEstabs)
        if alg == "bruteforce" || cardinality(S) <= 1e6
            S.sCWEstabs = _weightenumeratorBFQ(S.stabs, S.charvec, parent(S.sCWElogs.polynomial))
        else
            # trellis solution here
        end
    end

    if set ∈ ["all", "quotient"] && ismissing(S.sCWEdual)
        if alg == "bruteforce" || BigInt(characteristic(S.F))^(S.n + S.k) <= 3e6
            S.sCWEdual = _weightenumeratorBFQ(vcat(S.stabs, logicalsmatrix(S)), S.charvec, parent(S.sCWElogs.polynomial))
        else
            # trellis solution here
        end
    end
    
    if !ismissing(S.sCWEstabs) && !ismissing(S.sCWEdual)
        # compute minimum distance here
        poly = WeightEnumerator(S.sCWEdual.polynomial - S.sCWEstabs.polynomial, "complete")
        HWE = CWEtoHWE(poly)
        S.d = minimum(filter(x->x!=0, [collect(exponent_vectors(HWE.polynomial))[i][1]
            for i in 1:length(HWE.polynomial)]))
    end

    if type == "complete"
        set == "all" && return S.sCWEstabs, S.sCWEdual, S.sCWElogs, poly
        set == "stabilizers" && return S.sCWEstabs
        set == "logicals" && return S.sCWElogs
        return poly
    else
        set == "all" && return CWEtoHWE(S.sCWEstabs), CWEtoHWE(S.sCWEdual), CWEtoHWE(S.sCWElogs), CWEtoHWE(poly)
        set == "stabilizers" && return CWEtoHWE(S.sCWEstabs)
        set == "logicals" && return CWEtoHWE(S.sCWElogs)
        return HWE
    end
end

# MAGMA returns this format
# [ <0, 1>, <4, 105>, <6, 280>, <8, 435>, <10, 168>, <12, 35> ]
function weightdistribution(S::AbstractStabilizerCode, alg::String="auto", format::String="full",
    set::String="all")

    alg ∈ ["auto", "trellis", "bruteforce"] || error("Algorithm `$alg` is not implemented in weightenumerator.")
    format ∈ ["full", "compact"] || error("Unknown value for parameter format: $format; expected `full` or `compact`.")
    set ∈ ["all", "stabilizers", "logicals", "quotient"] || throw(ArgumentError("Unsupported set type '$set'. Expected 'all', 'stabilizers', 'logicals', 'quotient'."))

    wtenums = weightenumerator(S, "Hamming", alg, set)

    if format == "compact"
        if length(wtenums) == 1
            wtdist = Vector{Tuple}()
            for i in 1:length(wtenums.polynomial)
                push!(wtdist, (exponent_vector(wtenums.polynomial, i)[1],
                    coeff(wtenums.polynomial, i)))
            end
        else
            wtdist = Vector{Vector{Tuple}}()
            for wtenum in wtenums
                wtdistinner = Vector{Tuple}()
                for i in 1:length(wtenum.polynomial)
                    push!(wtdistinner, (exponent_vector(wtenum.polynomial, i)[1],
                        coeff(wtenum.polynomial, i)))
                end
                push!(wtdist, wtdistinner)
            end
        end
    else
        if length(wtenums) == 1
            K = base_ring(wtenums.polynomial)
            wtdist = zero_matrix(K, 1, S.n + 1)
            for i in 1:length(wtenums.polynomial)
                wtdist[1, exponent_vector(wtenums.polynomial, i)[1] + 1] = coeff(wtenums.polynomial, i)
            end
        else
            K = base_ring(wtenums[1].polynomial)
            wtdist = [] #Vector{Vector{K}}()
            for wtenum in wtenums
                wtdistinner = zero_matrix(K, 1, S.n + 1)
                for i in 1:length(wtenum.polynomial)
                    # println(coeff(wtenum.polynomial, i))
                    wtdistinner[1, exponent_vector(wtenum.polynomial, i)[1] + 1] = coeff(wtenum.polynomial, i)
                end
                push!(wtdist, wtdistinner)
            end
        end
    end
    return wtdist
end

function weightenumeratorQ(T::Trellis, type::String="complete")
    type ∈ ["complete", "Hamming"] || error("Unsupported weight enumerator type '$type'. Expected 'complete' or 'Hamming'.")

    if type == "complete" && !ismissing(T.CWE)
        return T.CWE
    elseif type == "Hamming" && !ismissing(T.CWE)
        return CWEtoHWE(T.CWE)
    end

    # if this ever changes or permutes will have to store with T
    elms = collect(T.code.E)
    lookup = Dict(value => key for (key, value) in enumerate(elms))

    p = Int(characteristic(T.code.E))
    iseven(p) ? nth = 2 * p : nth = p
    K, ω = CyclotomicField(nth, "ω")
    R, vars = PolynomialRing(K, length(elms))

    n = T.code.n
    charvec = T.code.charvec
    V = T.vertices
    E = T.edges
    V[1][1].polynomial = R(1)
    bit = 1
    for i in 2:length(V)
        # for (j, v) in enumerate(V[i])
        Threads.@threads for j in 1:length(V[i])
            v = V[i][j]
            outer = R(0)
            for e in E[i - 1][j]
                innerbit = deepcopy(bit)
                parity = 0
                for k in e.label
                    if !iszero(coeff(k, 0))
                        parity += data(charvec[innerbit])
                    end
                    if !iszero(coeff(k, 1))
                        parity += data(charvec[innerbit + n])
                    end
                    innerbit += 1
                end

                inner = deepcopy(V[i - 1][e.outvertex].polynomial)
                for k in e.label
                    inner *= ω^parity * vars[lookup[k]]
                end
                outer += inner # this prevents the e loop from being thread-safe
            end
            v.polynomial = outer
        end
        bit += length(E[i - 1][1][1].label)
    end
    T.CWE = WeightEnumerator(V[end][1].polynomial, "complete")

    # # currently Missing is not an option but how to implement dual trellis
    if !isshifted(T) && !ismissing(T.code)
        T.code.sCWEstabs = T.CWE
    end

    # # clean up vertices
    # for i in 1:length(V)
    #     for v in V[i]
    #         v.polynomial = missing
    #     end
    # end

    # display(T.CWE.polynomial)
    if type == "Hamming"
        return CWEtoHWE(T.CWE)
    end
    return T.CWE
end

"""
    weightplot(S::AbstractStabilizerCode, alg::String="auto", type::String="stabilizer")

Return a bar plot of the weight distribution related to `S`.

If `type` is `stabilizer`, the weight distribution of the stabilizers are computed.
If `type` is `normalizer`, the weight distrbution of the normalizer of the stabilizers
are computed. If `type` is `quotient`, the weight distrbution of the normalizer mod the
stabilizers (logical representatives only) is computed.
"""
function weightplot(S::AbstractStabilizerCode, alg::String="auto", type::String="stabilizer")
    type ∈ ["stabilizer", "normalizer", "quotient"] || throw(ArgumentError("Unknown value $type for parameter type."))

    wtdist = weightdistribution(S, alg, type, "full")
    xticks = findall(x->x>0, vec(wtdist)) .- 1
    yticks = [wtdist[i] for i in 1:length(wtdist) if !iszero(wtdist[i])]
    if type == "stabilizer"
        titlestr = "Stabilizer Weight Distribution"
    elseif type == "normalizer"
        titlestr = "Normalizer Weight Distribution"
    else
        titlestr = "Quotient Weight Distribution"
    end
    ismissing(S.d) ? (title="$titlestr - [$(S.n), $(S.k)]";) :
        title="$titlestr - [$(S.n), $(S.k), $(S.d)]"
    f = bar(0:S.n, wtdist', bar_width=1, xticks=xticks, yticks=yticks,
        legend=false, xlabel="Weight", ylabel="Number of Terms", title=title)
    show(f)
    return f
end

"""
    weightplotCSSX(S::AbstractCSSCode, alg::String="auto")

Return a bar plot of the weight distribution of the `X` stabilizers.
"""
function weightplotCSSX(S::AbstractCSSCode, alg::String="auto")
    C = LinearCode(S.Xstabs)
    wtdist = weightdistribution(C, alg, "full")
    xticks = findall(x->x>0, vec(wtdist)) .- 1
    yticks = [wtdist[i] for i in 1:length(wtdist) if !iszero(wtdist[i])]
    f = bar(0:C.n, wtdist', bar_width=1, xticks=xticks, yticks=yticks,
        legend=false, xlabel="Weight", ylabel="Number of Terms",
        title="X-Weight Distribution")
    show(f)
    return f
end

"""
    weightplotCSSZ(S::AbstractCSSCode, alg::String="auto")

Return a bar plot of the weight distribution of the `Z` stabilizers.
"""
function weightplotCSSZ(S::AbstractCSSCode, alg::String="auto")
    C = LinearCode(S.Zstabs)
    wtdist = weightdistribution(C, alg, "full")
    xticks = findall(x->x>0, vec(wtdist)) .- 1
    yticks = [wtdist[i] for i in 1:length(wtdist) if !iszero(wtdist[i])]
    f = bar(0:C.n, wtdist', bar_width=1, xticks=xticks, yticks=yticks,
        legend=false, xlabel="Weight", ylabel="Number of Terms",
        title="Z-Weight Distribution")
    show(f)
    return f

end

"""
    weightplotCSS(S::AbstractCSSCode, alg::String="auto")

Return bar plots of the weight distribution of the both the
`X` and 'Z' stabilizers, separately.
"""
function weightplotCSS(S::AbstractCSSCode, alg::String="auto")
    C = LinearCode(S.Xstabs)
    wtdist = weightdistribution(C, alg, "full")
    xticks = findall(x->x>0, vec(wtdist)) .- 1
    yticks = [wtdist[i] for i in 1:length(wtdist) if !iszero(wtdist[i])]
    fX = bar(0:C.n, wtdist', bar_width=1, xticks=xticks, yticks=yticks,
        legend=false, xlabel="Weight", ylabel="Number of Terms",
        title="X-Weight Distribution")

    # okay to overwrite
    C = LinearCode(S.Zstabs)
    wtdist = weightdistribution(C, alg, "full")
    xticks = findall(x->x>0, vec(wtdist)) .- 1
    yticks = [wtdist[i] for i in 1:length(wtdist) if !iszero(wtdist[i])]
    fZ = bar(0:C.n, wtdist', bar_width=1, xticks=xticks, yticks=yticks,
        legend=false, xlabel="Weight", ylabel="Number of Terms",
        title="Z-Weight Distribution")
    
    f = Plots.plot(fX, fZ, layout=(1, 2))
    show(f)
    return f
end

"""
    support(S::AbstractStabilizerCode, alg::String="auto", type::String="stabilizer")

Returns the support related to `S`.

The support is the collection of nonzero exponents of the Hamming
weight enumerator. If `type` is `stabilizer`, the support of the stabilizers are computed.
If `type` is `normalizer`, the support of the normalizer of the stabilizers
are computed. If `type` is `quotient`, the support of the normalizer mod the
stabilizers (logical representatives only) is computed.
"""
support(S::AbstractStabilizerCode, alg::String="auto", type::String="stabilizer") =
    [i for (i, _) in weightdistribution(S, alg, type, "compact")]

#############################
     # Minimum Distance
#############################

"""
    minimumdistance(Q::AbstractStabilizerCode, alg::String="trellis", sect::Bool=false)

Return the minimum distance of the stabilizer code if known, otherwise computes it.

"""
function minimumdistance(S::AbstractStabilizerCode, alg::String="auto", verbose::Bool=false)
    !ismissing(S.d) && return S.d

    # these should be different? weight? auto? BZ?
    alg ∈ ["auto", "trellis", "bruteforce"] || error("Algorithm `$alg` is not implemented in weightenumerator.")

    if iszero(S.k)
        # "Quantum Error Correction Via Codes Over GF(4)"
        # the distance of an [𝑛,0] code is defined as the smallest non-zero weight of any stabilizer in the code
    else
        # something like this
        if alg == "auto"
            weightenumerator(S, "Hamming", "auto", "quotient")
        elseif alg == "trellis"
            TOFstabs = trellisorientedformadditive(S.stabs)
            TOFnorm = trellisorientedformadditive(S.dualgens)
            boundaries, numEsectprimal = optimalsectionalizationQ(TOFstabs, TOFnorm)
            verbose && println("Primal edges: $numEsectprimal")
            profilesprimal = trellisprofiles(TOFstabs, TOFnorm, boundaries, "symplectic")
            boundaries, numEsectdual = optimalsectionalizationQ(TOFnorm, TOFstabs)
            verbose && println("Dual edges: $numEsectdual")
            profilesdual = trellisprofiles(TOFnorm, TOFstabs, boundaries, "symplectic")
            if sum(profilesprimal[2]) <= sum(profilesdual[2])
                Tprimal = sect(S, "primal", true, false)
                TprimalHWE = weightenumeratorQ(Tprimal, "complete")
                TdualHWE = MacWilliamsIdentity(S, TprimalHWE, true)
                poly = TdualHWE.polynomial - TprimalHWE.polynomial
                S.d = minimum(filter(x->x!=0, [collect(exponent_vectors(poly))[i][1]
                    for i in 1:length(poly)]))
            else
                Tdual = sect(S, "dual", true, false)
                TdualHWE = weightenumeratorQ(Tdual, "Hamming")
                TprimalHWE = MacWilliamsIdentity(S, TdualHWE)
                poly = TdualHWE.polynomial - TprimalHWE.polynomial
                S.d = minimum(filter(x->x!=0, [collect(exponent_vectors(poly))[i][1]
                    for i in 1:length(poly)]))
            end

            # Tdual = syndrometrellis(S, "primal", true, true)
            # TdualHWE = weightenumeratorQ(Tdual, "Hamming")
            # Tdual = missing
            # println("Primal trellis complete")
            # Tstabs = syndrometrellis(S, "dual", true, true)
            # THWE = weightenumeratorQ(Tstabs, "Hamming")
            # Tstabs = missing
            # poly = TdualHWE.polynomial - THWE.polynomial
            # S.d = minimum(filter(x->x!=0, [collect(exponent_vectors(poly))[i][1]
            #     for i in 1:length(poly)]))
        else
            # brute force solution here
        end
         #TODO: purity - 
    end
   
    return S.d
end

function minimumdistanceXZ(S::AbstractCSSCode)
    (!ismissing(S.dz) && !ismissing(S.dx)) && return S.dz, S.dx

    # dz = min(CX^⟂ \ CZ)
    # dx = min(CZ^⟂ \ CX)

    # need to make these if they are missing
    if !ismissing(S.ZorigCode)
        C1 = S.ZorigCode
        C2 = S.Xorigcode
    else
        C1 = LinearCode(S.Zstabs)
        C2 = LinearCode(S.Xstabs)
    end
    C1wtenum = weightenumerator(C1, "Hamming")
    C2wtenum = weightenumerator(C2, "Hamming")
    C1dualwtenum = MacWilliamsIdentity(C1, C1wtenum)
    C2dualwtenum = MacWilliamsIdentity(C2, C2wtenum)
    C1setdiffC2wtenum = C1dualwtenum.polynomial - C2dualwtenum.polynomial
    C2dualsetdiffC1dualwtenum = C2dualwtenum.polynomial - C1dualwtenum.polynomial
    S.dz = minimum(filter(x->x!=0, [collect(exponent_vectors(C1setdiffC2wtenum))[i][1]
        for i in 1:length(C1setdiffC2wtenum)]))
    S.dx = minimum(filter(x->x!=0, [collect(exponent_vectors(C2dualsetdiffC1dualwtenum))[i][1]
        for i in 1:length(C2dualsetdiffC1dualwtenum)]))
    # the above commands will set Ci.d
    (S.dx == C2.d && S.dz == C1.d) ? (S.pure = true;) : (S.pure = false;)
    return S.dz, S.dx

    # some other paper has this as the formula
    # expsC1setdiffC2 = filter(x->x!=0, [collect(exponent_vectors(C1setdiffC2wtenum))[i][1]
    #     for i in 1:length(C1setdiffC2wtenum)])
    # expsC2dualsetdiffC2dual = filter(x->x!=0, [collect(exponent_vectors(C2dualsetdiffC1dualwtenum))[i][1]
    #     for i in 1:length(C2dualsetdiffC1dualwtenum)])
    # exps = vcat(expsC1setdiffC2, expsC2dualsetdiffC2dual)
    # S.dx = minimum(exps)
    # S.dz = maximum(exps)
end

function minimumdistanceX(S::AbstractCSSCode)
    ismissing(S.dx) || return S.dx
    
     # need to make these if they are missing
     if !ismissing(S.ZorigCode)
        C1 = S.ZorigCode
        C2 = S.Xorigcode
    else
        C1 = LinearCode(S.Zstabs)
        C2 = LinearCode(S.Xstabs)
    end
    C1wtenum = weightenumerator(C1, "Hamming")
    C2wtenum = weightenumerator(C2, "Hamming")
    C1dualwtenum = MacWilliamsIdentity(C1, C1wtenum)
    C2dualwtenum = MacWilliamsIdentity(C2, C2wtenum)
    C2dualsetdiffC1dualwtenum = C2dualwtenum.polynomial - C1dualwtenum.polynomial
    S.dx = minimum(filter(x->x!=0, [collect(exponent_vectors(C2dualsetdiffC1dualwtenum))[i][1]
        for i in 1:length(C2dualsetdiffC1dualwtenum)]))
    return S.dx
end

function minimumdistanceZ(S::AbstractCSSCode)
    ismissing(S.dz) || return S.dz

    # need to make these if they are missing
    if !ismissing(S.ZorigCode)
        C1 = S.ZorigCode
        C2 = S.Xorigcode
    else
        C1 = LinearCode(S.Zstabs)
        C2 = LinearCode(S.Xstabs)
    end
    C1wtenum = weightenumerator(C1, "Hamming")
    C2wtenum = weightenumerator(C2, "Hamming")
    C1dualwtenum = MacWilliamsIdentity(C1, C1wtenum)
    C2dualwtenum = MacWilliamsIdentity(C2, C2wtenum)
    C1setdiffC2wtenum = C1dualwtenum.polynomial - C2dualwtenum.polynomial
    S.dz = minimum(filter(x->x!=0, [collect(exponent_vectors(C1setdiffC2wtenum))[i][1]
        for i in 1:length(C1setdiffC2wtenum)]))
    return S.dz
end

function ispure(S::AbstractStabilizerCode)
    ismissing(S.pure) || return S.pure
    minimumdistance(S) # this needs to force the weight enumerator approach
    return S.pure
end

function ispure(S::AbstractCSSCode)
    ismissing(S.pure) || return S.pure
    minimumdistanceXZ(S)
    return S.pure
end
