# Copyright (c) 2022, Eric Sabo, Benjamin Ide, Michael Vasmer
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
    return WeightEnumerator(poly, :complete)
end

"""
    CWEtoHWE(CWE::WeightEnumerator)

Return the Hamming weight enumerator associated with the complete weight enumerator `CWE`.
"""
function CWEtoHWE(CWE::WeightEnumerator)
    CWE.type == :complete || throw(ArgumentError("Not a complete weight enumerator"))

    R, (x, y) = PolynomialRing(base_ring(CWE.polynomial), ["x", "y"])
    poly = R(0)
    for i in 1:length(CWE.polynomial)
        exps = exponent_vector(CWE.polynomial, i)
        poly += coeff(CWE.polynomial, i) * x^sum(exps[2:end]) * y^exps[1]
    end
    return WeightEnumerator(poly, :Hamming)
end

#############################
        # GrayCode
#############################

# will be sent to a new package by Ben
struct GrayCode
    n::Int # length of codewords
    k::Int # weight of codewords
    ns::Int # n for the subcode
    ks::Int # k for the subcode
    prefix::Vector{Int}
    prefix_length::Int
    mutate::Bool
end

GrayCode(n::Int, k::Int; mutate=false) = GrayCode(n, k, Int[], mutate=mutate)

function GrayCode(n::Int, k::Int, prefix::Vector{Int}; mutate=false)
    GrayCode(n, k, n - length(prefix), k - count(prefix .!= 0), prefix,
        length(prefix), mutate)
end

Base.IteratorEltype(::GrayCode) = Base.HasEltype()
Base.eltype(::GrayCode) = Array{Int, 1}
Base.IteratorSize(::GrayCode) = Base.HasLength()
@inline function Base.length(G::GrayCode)
    if 0 <= G.ks <= G.ns
        factorial(big(G.ns)) ÷ (factorial(big(G.ks)) * factorial(big(G.ns - G.ks)))
    else
        0
    end
end
Base.in(v::Vector{Int}, G::GrayCode) = length(v) == G.n && count(v .!= 0) == G.k && view(v, 1:G.prefix_length) == G.prefix

@inline function Base.iterate(G::GrayCode)
    0 <= G.ks <= G.ns || return nothing

    g = [i <= G.ks ? 1 : 0 for i = 1:G.ns + 1]
    τ = collect(2:G.ns + 2)
    τ[1] = G.ks + 1
    # to force stopping with returning the only valid vector when ks == 0 and ns > 0
    iszero(G.ks) && (τ[1] = G.ns + 1;)
    v = [G.prefix; g[end - 1:-1:1]]
    ((G.mutate ? v : copy(v);), (g, τ, G.ks, v))
end

@inline function Base.iterate(G::GrayCode, state)
    g, τ, t, v = state
    @inbounds begin
        i = τ[1]
        i < G.ns + 1 || return nothing

        τ[1] = τ[i]
        τ[i] = i + 1
        if g[i] == 1
            if t != 0
                g[t] = g[t] == 0 ? 1 : 0
                if t < G.ns + 1
                    v[G.prefix_length + G.ns + 1 - t] = g[t]
                end
            else
                g[i - 1] = g[i - 1] == 0 ? 1 : 0
                if i - 1 < G.ns + 1
                    v[G.prefix_length + G.ns + 2 - i] = g[i - 1]
                end
            end
            t = t + 1
        else
            if t != 1
                g[t - 1] = g[t - 1] == 0 ? 1 : 0
                if t - 1 < G.ns + 1
                    v[G.prefix_length + G.ns + 2 - t] = g[t - 1]
                end
            else
                g[i - 1] = g[i - 1] == 0 ? 1 : 0
                if i - 1 < G.ns + 1
                    v[G.prefix_length + G.ns + 2 - i] = g[i - 1]
                end
            end
            t = t - 1
        end

        g[i] = g[i] == 0 ? 1 : 0
        if i < G.ns + 1
            v[G.prefix_length + G.ns + 1 - i] = g[i]
        end

        if t == i - 1 || t == 0
            t = t + 1
        else
            t = t - g[i - 1]
            τ[i - 1] = τ[1]
            if t == 0
                τ[1] = i - 1
            else
                τ[1] = t + 1
            end
        end
    end

    ((G.mutate ? v : copy(v);), (g, τ, t, v))
end

#############################
        # Classical
#############################

"""
    Sternsattack(C::AbstractLinearCode, w::Int, p::Int, l::Int)

Search for codewords of `C` of weight `w` using Stern's attack and return any found.
"""
function Sternsattack(C::AbstractLinearCode, w::Int, p::Int, l::Int, numfind::Int=2, maxitrs::Int=50000)
    # requires 2 * x = 0
    Int(order(C.F)) == 2 || throw(ArgumentError("Only valid for binary codes."))
    if !ismissing(C.d)
        C.d <= w <= C.n || throw(ArgumentError("Target weight must be between the minimum distance and code length."))
    else
        1 <= w <= C.n || throw(ArgumentError("Target weight must be positive and no more than the code length."))
    end
    1 <= p <= ceil(C.k / 2) || throw(ArgumentError("p must be between 1 and k/2"))

    # G = generatormatrix(C, true)
    H = paritycheckmatrix(C, true)
    nr, nc = size(H)
    # nr, nc = size(G)
    1 <= l <= nr || throw(ArgumentError("l must be between 1 and n - k"))
    H2 = deepcopy(H)
    # G2 = deepcopy(G)
    # lenwvecs = Vector{fq_nmod_mat}()
    lenwvecs = Set{fq_nmod_mat}()
    count = 0
    found = 0
    while count < maxitrs
        # choose n - k random columns and row reduce these to the identity
        # going to use a modified, partial Fisher-Yates algorithm for the column indices
        # however, not all choices are going to give an information set, so we are going
        # to do the randomness as needed
        colindices = collect(1:nc)
        rowpivots = falses(nr)
        rowpivotsloc = zeros(Int, nr)
        have = 0
        loc = nc
        while have < nr
            i = rand(1:loc)
            nonzeros = []
            for j in 1:nr
                if !iszero(H2[j, i])
                    append!(nonzeros, j)
                end
            end
            used = false
            for nz in nonzeros
                # nonzero rows of columns may already be pivots and should be ignored
                if rowpivots[nz] == false
                    # eliminate here because eliminating later can choose pivots which get removed
                    # this probably only works for binary, so likely unnecessary
                    # if !iszero(H2[nz, colindices[i]]) && !isone(H2[nz, colindices[i]])
                    #     H2[nz, :] = H2[nz, colindices[i]]^-1 * H2[nz, :]
                    # end

                    for j in 1:nr
                        # go through all rows and eliminate all nonzeros in column using the chosen row
                        if j != nz && !iszero(H2[j, colindices[i]])
                            if !isone(H2[j, colindices[i]])
                                H2[j, :] = H2[j, :] - H2[j, :]^-1 * H2[nz, :]
                            else
                                H2[j, :] = H2[j, :] - H2[nz, :]
                            end
                        end
                    end

                    rowpivots[nz] = true
                    rowpivotsloc[nz] = i
                    if i != loc
                        temp = colindices[loc]
                        colindices[loc] = colindices[i]
                        colindices[i] = temp
                    end
                    have += 1
                    loc -= 1
                    used = true
                    break
                end
            end

            if !used
                if i != loc
                    temp = colindices[loc]
                    colindices[loc] = colindices[i]
                    colindices[i] = temp
                end
                loc -= 1
            end

            # if we need more than are available, restart
            if nr - have > loc
                colindices = collect(1:nc)
                rowpivots = falses(nr)
                rowpivotsloc = zeros(Int, nr)
                have = 0
                loc = nc
                count += 1
                # H2 = deepcopy(H)
            end
        end

        # randomly split the remaning column indices
        X = Vector{Int}()
        Y = Vector{Int}()
        for i in colindices[1:nc - nr]
            rand(Float16) <= 0.5 ? (append!(X, i);) : (append!(Y, i);)
        end

        # choose a random size-l subset Z of rows, again using Fisher-Yates
        Z = collect(1:nr)
        loc = nr
        num = nr - l
        while loc > num
            i = rand(1:loc)
            if i != loc
                temp = Z[loc]
                Z[loc] = Z[i]
                Z[i] = temp
            end
            loc -= 1
        end

        # search for codewords that have:
        # exactly p nonzero bits in X
        # exactly p nonzero bits in Y
        # 0 nonzero bits in Z
        # and exactly w − 2p nonzero bits in the remaining columns

        # storing this is a terrible idea but avoids repeated calculations of πB for each A
        # will have to check later if that's fast enough in parallel to do
        left = nr - l
        AπAs = Vector{Tuple{Vector{Int}, Vector{fq_nmod}}}()
        # for every size-p subset A of X
        for A in powerset(X, p, p)
            # compute the sum of the columns in A for each of those l rows
            # this represents the contribution from these columns to the dot product
            πA = [sum(H2[Z[left + 1], A])]
            for i in 2:l
                push!(πA, sum(H2[Z[left + i], A]))
            end
            push!(AπAs, (A, πA))
        end

        # compute π(B) for every size-p subset B of Y
        BπBs = Vector{Tuple{Vector{Int}, Vector{fq_nmod}}}()
        for B in powerset(Y, p, p)
            πB = [sum(H2[Z[left + 1], B])]
            for i in 2:l
                push!(πB, sum(H2[Z[left + i], B]))
            end
            push!(BπBs, (B, πB))
        end

        left = nc - nr
        for i in AπAs
            for j in BπBs
                # for each collision π(A) = π(B)
                # if they are the same, then the sum of them will have 0 dot product
                if i[2] == j[2]
                    AB = i[1] ∪ j[1]
                    # compute the sum of the 2p columns in A ∪ B over all rows
                    πAB = [sum(H2[1, AB])]
                    for k in 2:nr
                        push!(πAB, sum(H2[k, AB]))
                    end

                    # if the sum has weight w − 2p
                    if wt(πAB) == w - 2 * p
                        # add the corresponding w − 2p columns in the (n − k) × (n − k) submatrix
                        # these, together with A and B, form a codeword of weight w
                        # so any time there is a 1 at index i, put a 1 in the i-th position of the end of colindices
                        # that will give w - 2p 1's at locations outside of A ∪ B
                        veclenw = zeros(C.F, 1, nc)
                        for k in 1:nr
                            if isone(πAB[k])
                                veclenw[1, rowpivotsloc[k]] = C.F(1)
                            end
                        end
                        # now add back in the 2p 1's at locations in A ∪ B to get a vector of weight w
                        for k in AB
                            veclenw[1, k] = C.F(1)
                        end
                        # println(veclenw)
                        test = matrix(C.F, 1, nc, [coeff(x, 0) for x in veclenw])
                        # println(iszero(H * transpose(test)))
                        if iszero(H * transpose(test))
                            println(test)
                            push!(lenwvecs, test)
                            found += 1
                            if found == numfind
                                return lenwvecs
                            end
                        end
                    end
                end
            end
        end
        count += 1
    end
    return lenwvecs
end

#############################
  # Enumeration Based Algs
#############################

# TODO:
# 2. extend GrayCode to non-binary
# 3. reduce ops by going to A_i instead of G_i
# 4. unit tests
# 7. decoding attacks from thesis

# TODO: this does not produce the optimal set of matrices
# see section 7.3 of White's thesis for comments on this
function _getinformationsets(G::fq_nmod_mat, alg::String)
    # alg ∈ ["Brouwer", "Zimmermann", "White", "Chen"] || throw(ArgumentError("Expected 'Brouwer', 'Zimmermann', or 'White'."))
    nr, nc = size(G)
    genmats = []
    # Gi = deepcopy(G)
    rnk = nr
    if alg == "Brouwer"
        for i in 0:Int(floor(nc / nr)) - 1
            rnk, Gi, Pi = CodingTheory._rref_col_swap(G, 1:nr, i * rnk + 1:nc)
            push!(genmats, (rnk, Gi, Pi))
            # Ai = Gi[:, setdiff(1:nc, i * nr + 1:(i + 1) * nr)]
            # push!(genmats, Ai)
        end
    elseif alg == "Zimmermann"
        for i in 0:Int(floor(nc / nr))
            rnk, Gi, Pi = CodingTheory._rref_col_swap(G, 1:nr, i * rnk + 1:nc)
            push!(genmats, (rnk, Gi, Pi))
            # display(Gi)
            # println(rnk)
            # Ai = Gi[:, setdiff(1:nc, i * nr + 1:i * nr + rnk)]
            # display(Ai)
            # println(" ")
            # push!(genmats, (rnk, Ai))
        end
    elseif alg == "White"
        # the expansion factor of the code
        for i in 0:div(nc, nr) - 1
            # could use Gi here instead of G
            rnk, Gi, Pi = CodingTheory._rref_col_swap(G, 1:nr, i * nr + 1:(i + 1) * nr)
            # display(Gi)
            # println(rnk)
            push!(genmats, (rnk, Gi, Pi))
        end
    elseif alg == "Chen"
        # should remove this in the future as this should simply be (k, Gstand, P)
        rnk, Gi, Pi = CodingTheory._rref_col_swap(G, 1:nr, 1:nc)
        push!(genmats, (rnk, Gi, Pi))
    end
    return genmats
end

@inline function _lowerbounds(r::Int, n::Int, k::Int, l::Int, rankdefs::Vector{Int}, type::String)
    if type == "BZ"
        h = length(rankdefs)
        lower = 0
        for i in 1:h
            lower += maximum([0, r - rankdefs[i]])
        end
    elseif type == "Chen"
        lower = Int(ceil(n * r / k))
    elseif type == "White"
        lower = 0
        for i in 1:l
	        lower += Int(ceil(n * maximum([0, r - rankdefs[i]]) / (l * (k + rankdefs[i]))))
        end
    end
    return lower
end

"""
    Graycodemindist(C::AbstractLinearCode, verbose::Bool=false)

Return the minimum distance of `C` using a deterministic algorithm based on enumerating
constant weight codewords of the binary reflected Gray code. If a word of minimum weight
is found before the lower and upper bounds cross, it is returned; otherwise, `missing`
is returned.
"""
# In the thesis, the h and GrayCode for loops are switched. This seems inefficient
# but does allow one to update the lower bound using a single matrix after every
# GrayCode loop has finished. I believe only enumerating the very expensive GrayCode
# once will overcome this. Could be tested though.
function Graycodemindist(C::AbstractLinearCode, verbose::Bool=false)
    ordF = Int(order(C.F))
    ordF == 2 || throw(ArgumentError("Currently only implemented for binary codes."))

    p = Int(characteristic(C.F))
    found = missing
    perm = missing

    G = generatormatrix(C)
    if typeof(C) <: AbstractCyclicCode
        verbose && println("Detected a cyclic code, using Chen's adaption.")
        genmats = _getinformationsets(G, "Chen")
    elseif typeof(C) <: AbstractQuasiCyclicCode
        verbose && println("Detected a quasi-cyclic code, using White's adaption.")
        genmats = _getinformationsets(G, "White")
    else
        genmats = _getinformationsets(G, "Zimmermann")
    end
    genmatsJulia = [deepcopy(FpmattoJulia(x[2])') for x in genmats]
    h = length(genmatsJulia)
    rankdefs = zeros(Int, h)
    if verbose
        print("Generated $h information sets with ranks: ")
        for i in 1:h
            i == h ? (println(genmats[i][1]);) : (print("$(genmats[i][1]), "))
            # will only be using the rank deficits here
            # at the moment, the information sets are always disjoint so the relative
            # rank is zero
            rankdefs[i] = C.k - genmats[i][1]
        end
    end
    
    evenflag = false
    devenflag = false
    tevenflag = false
    ordF == 2 && (evenflag = iseven(C);)
    evenflag && (devenflag = isdoublyeven(C);)
    devenflag && (tevenflag = istriplyeven(C);)
    if verbose
        tevenflag && println("Detected a triply even code.")
        (!tevenflag && devenflag) && println("Detected a doubly even code.")
        (!tevenflag && !devenflag && evenflag) && println("Detected an even code.")
    end
    
    flag = true
    while flag
        for (j, g) in enumerate(genmatsJulia)
            # redo j = 1 to grab index and perm
            w, i = _minwtcol(g)
            if w <= C.ubound
                found = g[:, i]
                C.ubound = w
                perm = j
            end
        end
        # if restarted, found is missing
        ismissing(found) ? (C.ubound = C.n;) : (flag = false;)
    end
    verbose && println("Upper bound after row analysis: $(C.ubound)")
    verbose && !ismissing(found) && println("Found element matching upper bound.")

    numthrds = Threads.nthreads()
    verbose && println("Detected $numthrds threads.")
    power = 0
    for i in 1:20
        if 2^i > numthrds
            power = i - 1
            break
        end
    end

    for r in 1:C.k
        if typeof(C) <: AbstractCyclicCode
            newlower = _lowerbounds(r, C.n, C.k, 0, [0], "Chen")
            C.lbound < newlower && (C.lbound = newlower;)
        elseif typeof(C) <: AbstractQuasiCyclicCode
            newlower = _lowerbounds(r, C.n, C.k, C.l, rankdefs, "White")
            C.lbound < newlower && (C.lbound = newlower;)
        else
            newlower = _lowerbounds(r, C.n, C.k, 0, rankdefs, "BZ")
            C.lbound < newlower && (C.lbound = newlower;)
        end
        # an even code can't have have an odd minimum weight
        (!tevenflag && !devenflag && evenflag) && (C.lbound += C.lbound % 2;)
        (!tevenflag && devenflag) && (C.lbound += 4 - C.lbound % 4;)
        tevenflag && (C.lbound += 8 - C.lbound % 8;)

        verbose && println("r: $r")
        verbose && println("Lower bound: $(C.lbound)")
        verbose && println("Upper bound: $(C.ubound)")
        if C.lbound >= C.ubound
            C.d = C.ubound
            if !ismissing(found)
                y = matrix(C.F, 1, C.n, found)
                ismissing(genmats[perm][3]) || (y = y * genmats[perm][3];)
                iszero(C.H * transpose(y)) ? (return C.ubound, y;) : (return C.ubound, missing;)
            else
                return C.ubound, found
            end
        end

        flag = Threads.Atomic{Bool}(true)
        # numthrds = 1
        # power = 0
        # total = Threads.Atomic{Int}(0)
        uppers = [C.ubound for _ in 1:numthrds]
        founds = [found for _ in 1:numthrds]
        perms = [perm for _ in 1:numthrds]
        Threads.@threads for m in 1:numthrds
            c = zeros(Int, C.n)
            prefix = digits(m - 1, base=2, pad=power)
            # count = Threads.Atomic{Int}(0)
            for u in GrayCode(C.k, r, prefix, mutate=true)
                # count += 1
                # Threads.atomic_add!(count, 1)
                if flag[]
                    for i in 1:h
                        if r - rankdefs[i] > 0
                            LinearAlgebra.mul!(c, genmatsJulia[i], u)
                            w = 0
                            @inbounds for j in 1:C.n
                                c[j] % p != 0 && (w += 1;)
                            end

                            if uppers[m] > w
                                uppers[m] = w
                                founds[m] = c
                                perms[m] = i
                                # verbose && println("Adjusting upper bound: $upperlocal")
                                if C.lbound == uppers[m]
                                    Threads.atomic_cas!(flag, true, false)
                                    break
                                end
                            end
                        end
                    end
                end
            end
            # println(count)
            # Threads.atomic_add!(total, count[])
        end
        loc = argmin(uppers)
        C.ubound = uppers[loc]
        found = founds[loc]
        perm = perms[loc]
        # println("total: $total")
        if !flag[]
            C.d = C.ubound
            if !ismissing(found)
                y = matrix(C.F, 1, C.n, found)
                ismissing(genmats[perm][3]) || (y = y * genmats[perm][3];)
                iszero(C.H * transpose(y)) ? (return C.ubound, y;) : (return C.ubound, missing;)
            else
                return C.ubound, found
            end
        end
    end
end

"""
    wordsofweight(C::AbstractLinearCode, lbound::Int, ubound::Int, verbose::Bool=false)

Return all the codewords of `C` of Hamming weight in the range `[lbound, ubound]`.
"""
function wordsofweight(C::AbstractLinearCode, lbound::Int, ubound::Int, verbose::Bool=false)
    ordF = Int(order(C.F))
    ordF == 2 || throw(ArgumentError("Currently only implemented for binary codes."))

    1 <= lbound <= ubound <= C.n || throw(ArgumentError("Expected 1 <= lbound <= ubound <= C.n"))
    if lbound < C.n / 2 && ordF == 2
        # faster to enumerate backwards, but only in binary
        return _wordsofweighthigh(C, lbound, ubound, verbose)
    end

    p = Int(characteristic(C.F))
    G = generatormatrix(C)
    if typeof(C) <: AbstractCyclicCode
        verbose && println("Detected a cyclic code, using Chen's adaption.")
        genmats = _getinformationsets(G, "Chen")
    elseif typeof(C) <: AbstractQuasiCyclicCode
        verbose && println("Detected a quasi-cyclic code, using White's adaption.")
        genmats = _getinformationsets(G, "White")
    else
        genmats = _getinformationsets(G, "Zimmermann")
    end
    genmatsJulia = [FpmattoJulia(x[2])' for x in genmats]
    h = length(genmatsJulia)
    rankdefs = zeros(Int, h)
    if verbose
        print("Generated $h information sets with ranks: ")
        for i in 1:h
            i == h ? (println(genmats[i][1]);) : (print("$(genmats[i][1]), "))
            # will only be using the rank deficits here
            # at the moment, the information sets are always disjoint so the relative
            # rank is zero
            rankdefs[i] = C.k - genmats[i][1]
        end
    end

    evenflag = false
    devenflag = false
    tevenflag = false
    ordF == 2 && (evenflag = iseven(C);)
    evenflag && (devenflag = isdoublyeven(C);)
    devenflag && (tevenflag = istriplyeven(C);)
    if verbose
        tevenflag && println("Detected a triply even code.")
        (!tevenflag && devenflag) && println("Detected a doubly even code.")
        (!tevenflag && !devenflag && evenflag) && println("Detected an even code.")
    end

    numthrds = Threads.nthreads()
    verbose && println("Detected $numthrds threads.")
    power = 0
    for i in 1:20
        if 2^i > numthrds
            power = i - 1
            break
        end
    end

    W = Set{fq_nmod_mat}()
    for r in 1:C.k
        if typeof(C) <: AbstractCyclicCode
            lower = _lowerbounds(r, C.n, C.k, 0, [0], "Chen")
        elseif typeof(C) <: AbstractQuasiCyclicCode
            lower = _lowerbounds(r, C.n, C.k, C.l, rankdefs, "White")
        else
            lower = _lowerbounds(r, C.n, C.k, 0, rankdefs, "BZ")
        end
        # an even code can't have have an odd minimum weight
        (!tevenflag && !devenflag && evenflag) && (lower += lower % 2;)
        (!tevenflag && devenflag) && (lower += 4 - lower % 4;)
        tevenflag && (lower += 8 - lower % 8;)

        verbose && println("r: $r")
        verbose && println("Lower bound: $lower")
        if lower >= ubound
            return W
        end

        Ws = [Set{fq_nmod_mat}() for _ in 1:numthrds]
        Threads.@threads for m in 1:numthrds
            c = zeros(Int, C.n)
            prefix = digits(m - 1, base=2, pad=power)
            for u in GrayCode(C.k, r, prefix, mutate=true)
                for i in 1:h
                    LinearAlgebra.mul!(c, genmatsJulia[i], u)
                    w = 0
                    @inbounds for j in 1:C.n
                        c[j] % p != 0 && (w += 1;)
                    end

                    if lbound <= w <= ubound
                        c2 = matrix(C.F, 1, C.n, c)
                        ismissing(genmats[i][3]) || (c2 = c2 * genmats[i][3];)
                        push!(Ws[m], c2)
                    end
                end
            end
        end
        for m in 1:numthrds
            union!(W, Ws[m])
        end
    end
end

"""
    wordsofweight(C::AbstractLinearCode, bound::Int, verbose::Bool=false)

Return all the codewords of `C` of Hamming weight `bound`.
"""
wordsofweight(C::AbstractLinearCode, bound::Int, verbose::Bool=false) = wordsofweight(C, bound, bound, verbose)

# untested
# TODO: figure out if even weight upper weight needs to subtract or not
function _wordsofweighthigh(C::AbstractLinearCode, lbound::Int, ubound::Int, verbose::Bool=false)
    p = Int(characteristic(C.F))
    G = generatormatrix(C)
    if typeof(C) <: AbstractCyclicCode
        verbose && println("Detected a cyclic code, using Chen's adaption.")
        genmats = _getinformationsets(G, "Chen")
    elseif typeof(C) <: AbstractQuasiCyclicCode
        verbose && println("Detected a quasi-cyclic code, using White's adaption.")
        genmats = _getinformationsets(G, "White")
    else
        genmats = _getinformationsets(G, "Zimmermann")
    end
    genmatsJulia = [FpmattoJulia(x[2])' for x in genmats]
    h = length(genmatsJulia)
    rankdefs = zeros(Int, h)
    if verbose
        print("Generated $h information sets with ranks: ")
        for i in 1:h
            i == h ? (println(genmats[i][1]);) : (print("$(genmats[i][1]), "))
            # will only be using the rank deficits here
            # at the moment, the information sets are always disjoint so the relative
            # rank is zero
            rankdefs[i] = C.k - genmats[i][1]
        end
    end

    evenflag = false
    devenflag = false
    tevenflag = false
    ordF == 2 && (evenflag = iseven(C);)
    evenflag && (devenflag = isdoublyeven(C);)
    devenflag && (tevenflag = istriplyeven(C);)
    if verbose
        tevenflag && println("Detected a triply even code.")
        (!tevenflag && devenflag) && println("Detected a doubly even code.")
        (!tevenflag && !devenflag && evenflag) && println("Detected an even code.")
    end

    numthrds = Threads.nthreads()
    verbose && println("Detected $numthrds threads.")
    power = 0
    for i in 1:20
        if 2^i > numthrds
            power = i - 1
            break
        end
    end

    W = Set{fq_nmod_mat}()
    for r in C.k:-1:1
        if typeof(C) <: AbstractCyclicCode
            upper = _lowerbounds(r, C.n, C.k, 0, [0], "Chen")
        elseif typeof(C) <: AbstractQuasiCyclicCode
            upper = _lowerbounds(r, C.n, C.k, C.l, rankdefs, "White")
        else
            upper = _lowerbounds(r, C.n, C.k, 0, rankdefs, "BZ")
        end
        # # an even code can't have have an odd minimum weight
        # (!tevenflag && !devenflag && evenflag) && (upper += upper % 2;)
        # (!tevenflag && devenflag) && (upper += 4 - upper % 4;)
        # tevenflag && (upper += 8 - upper % 8;)

        verbose && println("r: $r")
        verbose && println("Upper bound: $upper")
        if upper < lbound
            return W
        end

        Ws = [Set{fq_nmod_mat}() for _ in 1:numthrds]
        Threads.@threads for m in 1:numthrds
            c = zeros(Int, C.n)
            prefix = digits(m - 1, base=2, pad=power)
            for u in GrayCode(C.k, r, prefix, mutate=true)
                for i in 1:h
                    LinearAlgebra.mul!(c, genmatsJulia[i], u)
                    w = 0
                    @inbounds for j in 1:C.n
                        c[j] % p != 0 && (w += 1;)
                    end

                    if lbound <= w <= ubound
                        c2 = matrix(C.F, 1, C.n, c)
                        ismissing(genmats[i][3]) || (c2 = c2 * genmats[i][3];)
                        push!(Ws[m], c2)
                    end
                end
            end
        end
        for m in 1:numthrds
            union!(W, Ws[m])
        end
    end
end

#TODO: change above to BigInts
"""
    partialweightdistribution(C::AbstractLinearCode, bound::Int, compact::Bool=false)

Return the partial weight distribution of `C` up to weight `bound`. If `compact` is false,
the result will be a `Vector{BigInt}` of length `length(C) + 1` whose `i`th entry is the
number of codewords of `C` of Hamming weight `i - 1`. Otherwise, the result is a
`Vector{Tuple{Int, BigInt}}` whose entries specify the nonzero indices and values of the
above.
"""
# TODO: type stability of this
function partialweightdistribution(C::AbstractLinearCode, bound::Int, compact::Bool=false)
	1 <= bound <= C.n || throw(ArgumentError("Bound must be between 1 and n."))

	W = wordsofweight(C, 1, bound)
    wtdist = zeros(BigInt, bound + 1)
    bio = BigInt(1)
    wtdist[1] = bio    
	for c in W
		wtdist[wt(c) + 1] += bio
	end
    
    if compact
        wtdistcomp = Vector{Tuple{Int, BigInt}}()
        for (i, x) in enumerate(wtdist)
            !iszero(x) && (push!(wtdistcomp, (i - 1, x)))
        end
        return wtdistcomp
    else
        return wtdist
    end
end

"""
	minimumwords(C::AbstractLinearCode)

Return the set of codewords of `C` with weight equal to the minimum distance.

This algorithm simultaneously computes the minimum distance and stores the words of
this weight that it finds, removing the repeated work of calling
`w = minimumdistance(C); W = wordsofweight(C, w);`
"""
function minimumwords(C::AbstractLinearCode)
    ordF = Int(order(C.F))
    ordF == 2 || throw(ArgumentError("Currently only implemented for binary codes."))

    p = Int(characteristic(C.F))
    found = missing
    perm = missing

    G = generatormatrix(C)
    if typeof(C) <: AbstractCyclicCode
        verbose && println("Detected a cyclic code, using Chen's adaption.")
        genmats = _getinformationsets(G, "Chen")
    elseif typeof(C) <: AbstractQuasiCyclicCode
        verbose && println("Detected a quasi-cyclic code, using White's adaption.")
        genmats = _getinformationsets(G, "White")
    else
        genmats = _getinformationsets(G, "Zimmermann")
    end
    genmatsJulia = [FpmattoJulia(x[2])' for x in genmats] # deepcopy doesn't save anythign here
    # genmatsJulia = [FpmattoJulia(x[2]) for x in genmats]
    h = length(genmatsJulia)
    rankdefs = zeros(Int, h)
    if verbose
        print("Generated $h information sets with ranks: ")
        for i in 1:h
            i == h ? (println(genmats[i][1]);) : (print("$(genmats[i][1]), "))
            # will only be using the rank deficits here
            # at the moment, the information sets are always disjoint so the relative
            # rank is zero
            rankdefs[i] = C.k - genmats[i][1]
        end
    end
    
    evenflag = false
    devenflag = false
    tevenflag = false
    ordF == 2 && (evenflag = iseven(C);)
    evenflag && (devenflag = isdoublyeven(C);)
    devenflag && (tevenflag = istriplyeven(C);)
    if verbose
        tevenflag && println("Detected a triply even code.")
        (!tevenflag && devenflag) && println("Detected a doubly even code.")
        (!tevenflag && !devenflag && evenflag) && println("Detected an even code.")
    end
    
    upper = C.n - C.k + 1
    verbose && println("Singleton upper bound: $upper")
    for (j, g) in enumerate(genmatsJulia)
        for i in 1:C.k
            w = wt(g[:, i]) # for transposed matrix
            if w < upper
                found = g[:, i]
                upper = w
                perm = j
            end
        end
    end
    verbose && println("Upper bound after row analysis: $upper")
    verbose && !ismissing(found) && println("Found vector for upper bound.")

    numthrds = Threads.nthreads()
    verbose && println("Detected $numthrds threads.")
    power = 0
    for i in 1:20
        if 2^i > numthrds
            power = i - 1
            break
        end
    end

    W = Set{fq_nmod_mat}()
    for r in 1:C.k
        if typeof(C) <: AbstractCyclicCode
            lower = _lowerbounds(r, C.n, C.k, 0, [0], "Chen")
        elseif typeof(C) <: AbstractQuasiCyclicCode
            lower = _lowerbounds(r, C.n, C.k, C.l, rankdefs, "White")
        else
            lower = _lowerbounds(r, C.n, C.k, 0, rankdefs, "BZ")
        end
        # an even code can't have have an odd minimum weight
        (!tevenflag && !devenflag && evenflag) && (lower += lower % 2;)
        (!tevenflag && devenflag) && (lower += 4 - lower % 4;)
        tevenflag && (lower += 8 - lower % 8;)

        verbose && println("r: $r")
        verbose && println("Lower bound: $lower")
        verbose && println("Upper bound: $upper")
        if lower >= upper
            C.d = upper
            return upper, W
        end

        uppers = [upper for _ in 1:numthrds]
        Ws = [Set{fq_nmod_mat}() for _ in 1:numthrds]
        Threads.@threads for m in 1:numthrds
            c = zeros(Int, C.n)
            prefix = digits(m - 1, base=2, pad=power)
            for u in GrayCode(C.k, r, prefix, mutate=true)
                for i in 1:h
                    LinearAlgebra.mul!(c, genmatsJulia[i], u)
                    w = 0
                    @inbounds for j in 1:C.n
                        c[j] % p != 0 && (w += 1;)
                    end

                    if w <= uppers[m]
                        if w < uppers[m]
                            uppers[m] = w
                            verbose && println("Adjusting upper bound: $upper")
                            Ws[m] = Set{fq_nmod_mat}()
                        end
                        # TODO: this is very expensive just to get erased
                        # maybe keep perms[m] and adjust outside loop
                        # allocations are more expensive inside threads
                        c2 = matrix(C.F, 1, C.n, c)
                        ismissing(genmats[i][3]) || (c2 = c2 * genmats[i][3];)
                        push!(Ws[m], c2)
                    end
                end
            end
        end
        loc = argmin(uppers)
        if upper > uppers[loc]
            upper = uppers[loc]
            W = Set{fq_nmod_mat}()
        end
        for m in 1:numthrds
            if uppers[m] == upper
                union!(W, Ws[m])
            end
        end
    end
end

#############################
    # Weight Enumerators
#############################

function weightenumeratorC(T::Trellis, type::Symbol=:complete)
    type ∈ [:complete, :Hamming] ||
        throw(ArgumentError("Unsupported weight enumerator type '$type'. Expected ':complete' or ':Hamming'."))

    if type == :complete && !ismissing(T.CWE)
        return T.CWE
    elseif type == :Hamming && !ismissing(T.CWE)
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
    T.CWE = WeightEnumerator(V[end][1].polynomial, :complete)

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

    type == :Hamming && return CWEtoHWE(T.CWE)
    return T.CWE
end

"""
    MacWilliamsIdentity(C::AbstractLinearCode, W::WeightEnumerator, dual::Symbol=:Euclidean)

Return the weight enumerator of the dual (`:Euclidean` or `:Hermitian`) of `C` obtained
by applying the MacWilliams identities to `W`.
"""
# TODO: remove C from this, store in WE struct
function MacWilliamsIdentity(C::AbstractLinearCode, W::WeightEnumerator, dual::Symbol=:Euclidean)
    dual ∈ [:Euclidean, :Hermitian] ||
        throw(ArgumentError("The MacWilliams identities are only programmed for the Euclidean and Hermitian duals."))
    (dual == :Hermitian && Int(order(C.F)) != 4) &&
        throw(ArgumentError("The MacWilliams identity for the Hermitian dual is only programmed for GF(4)."))

    if W.type == :Hamming
        # (1/|C|)W(y - x, y + (q - 1)x)
        R = parent(W.polynomial)
        vars = gens(R)
        return WeightEnumerator(divexact(W.polynomial(vars[2] - vars[1], vars[2] +
            (Int(order(C.F)) - 1) * vars[1]), cardinality(C)), :Hamming)
    end

    # complete weight enumerators
    if Int(order(C.F)) == 2
        # the complete and Hamming weight enumerators are the same in binary
        # (1/|C|)W(x_0 + (q - 1)x_1, x_0 - x_1)
        R = parent(W.polynomial)
        vars = gens(R)
        return WeightEnumerator(divexact(W.polynomial(vars[1] +
            (Int(order(C.F)) - 1) * vars[2], vars[1] - vars[2]),
            cardinality(C)), :complete)
    elseif Int(order(C.F)) == 3
        # (1/|C|)W(x_0 + x_1 + x_2, x_0 + ω x_1 + ω^2 x_2, x_0 + ω^2 x_1 + ω x_2)
        K, ζ = CyclotomicField(3, "ζ")
        R, vars = PolynomialRing(K, 3)
        # might have to switch this here
        poly = divexact(W.polynomial(
            vars[1] + vars[2] + vars[3],
            vars[1] + ζ * vars[2] + ζ^2 * vars[3],
            vars[1] + ζ^2 * vars[2] + ζ * vars[3]), cardinality(C))
        # works so far but now needs to recast down to the integer ring
        return WeightEnumerator(map_coefficients(c -> Nemo.ZZ(coeff(c, 0)), poly,
            parent=parent(W.polynomial)), :complete)
    elseif Int(order(C.F)) == 4
        # these order 4 formulas are from "Self-Dual Codes" by Rains and Sloane without proof
        # the differ in order from the formula in MacWilliams and Sloane used in the general
        # case below:
        #    x1 + x2 + x3 + x4
        #    x1 - x2 + x3 - x4
        #    x1 + x2 - x3 - x4
        #    x1 - x2 - x3 + x4
        # But that formula should depend on the chosen basis and character so I assume it's okay
        if dual == :Euclidean
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
                :complete)
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
                :complete)
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
                :complete)
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
            return WeightEnumerator(divexact(W.polynomial(funcargs...), cardinality(C)), :complete)
        end
    end
end

"""
    weightenumerator(C::AbstractLinearCode, type::Symbol=:complete, alg::String="auto")

Return either the `:complete` or `:Hamming` weight enumerator of `C` using the algorithm `alg`.
"""
function weightenumerator(C::AbstractLinearCode, type::Symbol=:complete, alg::String="auto")
    type ∈ [:complete, :Hamming] ||
        throw(ArgumentError("Unsupported weight enumerator type '$type'. Expected ':complete' or ':Hamming'."))
    alg ∈ ["auto", "trellis", "bruteforce"] ||
        throw(ArgumentError("Algorithm `$alg` is not implemented in weightenumerator."))

    if type == :complete && !ismissing(C.weightenum)
        return C.weightenum
    elseif type == :Hamming && !ismissing(C.weightenum)
        return CWEtoHWE(C.weightenum)
    end

    if alg == "auto"
        if cardinality(C) <= 1e6 # random cutoff
            C.weightenum = _weightenumeratorBF(C.G)
            HWE = CWEtoHWE(C.weightenum)
            C.d = minimum(filter(x->x!=0, [collect(exponent_vectors(HWE.polynomial))[i][1]
                for i in 1:length(HWE.polynomial)]))
            type == :Hamming && return HWE
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
            type == :Hamming && return HWE
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
        type == :Hamming && return HWE
        return C.weightenum
    end
end

"""
    weightdistribution(C::AbstractLinearCode, alg::String="auto", compact::Bool=true)

Return the weight distribution of `C` using the algorithm `alg`. If `compact` is false,
the result will be a `Vector{BigInt}` of length `length(C) + 1` whose `i`th entry is the
number of codewords of `C` of Hamming weight `i - 1`. Otherwise, the result is a
`Vector{Tuple{Int, BigInt}}` whose entries specify the nonzero indices and values of the
above.
"""
function weightdistribution(C::AbstractLinearCode, alg::String="auto", compact::Bool=true)
    alg ∈ ["auto", "trellis", "bruteforce"] ||
        throw(ArgumentError("Algorithm `$alg` is not implemented in weightenumerator."))

    ismissing(C.weightenum) && weightenumerator(C, :complete, alg)
    HWE = CWEtoHWE(C.weightenum)

    if compact
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
    weightplot(C::AbstractLinearCode, alg::String="auto")

Return a bar plot of the weight distribution of `C`.
"""
function weightplot(C::AbstractLinearCode, alg::String="auto")
    wtdist = weightdistribution(C, alg, true)
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
support(C::AbstractLinearCode) = [i for (i, _) in weightdistribution(C, "auto", true)]

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
function minimumdistance(C::AbstractLinearCode, alg::String="auto", sect::Bool=false, verbose::Bool=false)
    !ismissing(C.d) && return C.d

    alg ∈ ["auto", "Gray", "trellis", "Leon", "bruteforce", "wtdist"] ||
        throw(ArgumentError("Unexpected algorithm '$alg'."))
    
    if alg == "auto"
        D = dual(C)
        if cardinality(C) <= 1e6 # random cutoff
            C.weightenum = _weightenumeratorBF(C.G)
            HWE = CWEtoHWE(C.weightenum)
            C.d = minimum(filter(x->x!=0, [collect(exponent_vectors(HWE.polynomial))[i][1]
                for i in 1:length(HWE.polynomial)]))
            return C.d
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
            return C.d
        else
            return Graycodemindist(C, verbose)
        end
    elseif alg == "Gray"
        return Graycodemindist(C, verbose)
    elseif alg == "trellis"
        weightenumeratorC(syndrometrellis(C, "primal", false), type)
        return C.d
    elseif alg == "bruteforce"
        C.weightenum = _weightenumeratorBF(C.G)
        HWE = CWEtoHWE(C.weightenum)
        C.d = minimum(filter(x->x!=0, [collect(exponent_vectors(HWE.polynomial))[i][1]
            for i in 1:length(HWE.polynomial)]))
        return C.d
    elseif alg == "wtdist"
        HWE = weightenumerator(C, :Hamming, alg)
        !ismissing(C.d) && return C.d
        # this line should only be needed to be run if the weight enumerator is known
        # but the minimum distance is intentionally set to missing
        # ordering here can be a bit weird
        C.d = minimum(filter(x->x!=0, [collect(exponent_vectors(HWE.polynomial))[i][1]
            for i in 1:length(HWE.polynomial)]))
        return C.d
    # elseif alg == "Leon"
    #     Leon(C)
    end
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
    weightplotCSSX(S::AbstractStabilizerCodeCSS, alg::String="auto")

Return a bar plot of the weight distribution of the `X` stabilizers.
"""
function weightplotCSSX(S::AbstractStabilizerCodeCSS, alg::String="auto")
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
    weightplotCSSZ(S::AbstractStabilizerCodeCSS, alg::String="auto")

Return a bar plot of the weight distribution of the `Z` stabilizers.
"""
function weightplotCSSZ(S::AbstractStabilizerCodeCSS, alg::String="auto")
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
    weightplotCSS(S::AbstractStabilizerCodeCSS, alg::String="auto")

Return bar plots of the weight distribution of the both the
`X` and 'Z' stabilizers, separately.
"""
function weightplotCSS(S::AbstractStabilizerCodeCSS, alg::String="auto")
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

function minimumdistanceXZ(S::AbstractStabilizerCodeCSS)
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

function minimumdistanceX(S::AbstractStabilizerCodeCSS)
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

function minimumdistanceZ(S::AbstractStabilizerCodeCSS)
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

function ispure(S::AbstractStabilizerCodeCSS)
    ismissing(S.pure) || return S.pure
    minimumdistanceXZ(S)
    return S.pure
end

# TODO: pure for subsystem if no weight of gauge group is less than min dist



"""
    distrandCSS(hx::Matrix{Int}, hz::Matrix{Int}, num::Int, mindist::Int=0, debug::Int=0, field::GapObj=GAP.Globals.GF(2), maxav=Nothing)
Wrapper for the QDistRnd function DistRandCSS.
## QDistRnd documentation
- `num`: number of information sets to construct (should be large).
- `mindist`: the algorithm stops when distance equal or below `mindist` 
    is found and returns the result with negative sign. Set 
    `mindist` to 0 if you want the actual distance.
- `debug`: optional integer argument containing debug bitmap (default: `0`).
    - 1 (0s  bit set): print 1st of the vectors found.
    - 2 (1st bit set): check orthogonality of matrices and of the final vector.
    - 4 (2nd bit set): show occasional progress update.
    - 8 (3rd bit set): maintain cw count and estimate the success probability.
- `field` (Options stack): Galois field, default: GF(2).
- `maxav` (Options stack): if set, terminate when `<n>` greater than `maxav`, 
    see Section Emprirical. Not set by default.
"""
# Michael Vasmer
function QDistRndCSS(hx::Matrix{Int}, hz::Matrix{Int}, num::Int, mindist::Int,
    debug::Int=0, field::GapObj=GAP.Globals.GF(2), maxav=missing)

    GAP.Packages.load("QDistRnd");
    # Convert input matrices to GAP matrices over the given field
    e = GAP.Globals.One(field)
    hxg = GapObj(hx) * e
    hzg = GapObj(hz) * e
    if ismissing(maxav)
        dist = GAP.Globals.DistRandCSS(hxg, hzg, num, mindist, debug, field)
    else
        dist = GAP.Globals.DistRandCSS(hxg, hzg, num, mindist, debug, field, maxav)
    end
    return dist
end
