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

function _weight_enumerator_BF(G::CTMatrixTypes)
    E = base_ring(G)
    ord_E = Int(order(E))
    R, vars = PolynomialRing(Nemo.ZZ, ord_E)

    # See if we can return immediately
    iszero(G) && return WeightEnumerator(vars[1]^ncols(G), :complete)

    poly = R(0)
    nr = nrows(G)
    lookup = Dict(value => key for (key, value) in enumerate(collect(E)))

    # for iter in Iterators.product(Iterators.repeated(E, nr)...)
    for iter in Nemo.AbstractAlgebra.ProductIterator([E for _ in 1:nr], inplace = true)
        row = iter[1] * view(G, 1:1, :)
        for r in 2:nr
            if !iszero(iter[r])
                row += iter[r] * view(G, r:r, :)
            end
        end

        # TODO: can do this in one step, but is it faster?
        term = zeros(Int, 1, ord_E)
        for x in row
            term[lookup[x]] += 1
        end
        term_poly = R(1)
        for i in 1:ord_E
            term_poly *= vars[i]^term[i]
        end
        poly += term_poly
    end
    return WeightEnumerator(poly, :complete)
end

"""
    CWE_to_HWE(CWE::WeightEnumerator)

Return the Hamming weight enumerator associated with the complete weight enumerator `CWE`.
"""
function CWE_to_HWE(CWE::WeightEnumerator)
    CWE.type == :complete || throw(ArgumentError("Not a complete weight enumerator"))

    R, (x, y) = PolynomialRing(base_ring(CWE.polynomial), (:x, :y))
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

# TODO: fix doc string
"""
    Sterns_attack(C::AbstractLinearCode, w::Int, p::Int, l::Int)

Search for codewords of `C` of weight `w` using Stern's attack and return any found.
"""
function Sterns_attack(C::AbstractLinearCode, w::Int, p::Int, l::Int, num_find::Int=2, max_iters::Int=50000)
    # requires 2 * x = 0
    Int(order(C.F)) == 2 || throw(ArgumentError("Only valid for binary codes."))
    if !ismissing(C.d)
        C.d <= w <= C.n || throw(ArgumentError("Target weight must be between the minimum distance and code length."))
    else
        1 <= w <= C.n || throw(ArgumentError("Target weight must be positive and no more than the code length."))
    end
    1 <= p <= ceil(C.k / 2) || throw(ArgumentError("p must be between 1 and k/2"))

    # G = generator_matrix(C, true)
    H = paritycheckmatrix(C, true)
    nr, nc = size(H)
    # nr, nc = size(G)
    1 <= l <= nr || throw(ArgumentError("l must be between 1 and n - k"))
    H2 = deepcopy(H)
    # G2 = deepcopy(G)
    # len_w_vecs = Vector{fq_nmod_mat}()
    len_w_vecs = Set{fq_nmod_mat}()
    count = 0
    found = 0
    while count < max_iters
        # choose n - k random columns and row reduce these to the identity
        # going to use a modified, partial Fisher-Yates algorithm for the column indices
        # however, not all choices are going to give an information set, so we are going
        # to do the randomness as needed
        col_indices = collect(1:nc)
        row_pivots = falses(nr)
        row_pivots_loc = zeros(Int, nr)
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
                if row_pivots[nz] == false
                    # eliminate here because eliminating later can choose pivots which get removed
                    # this probably only works for binary, so likely unnecessary
                    # if !iszero(H2[nz, col_indices[i]]) && !isone(H2[nz, col_indices[i]])
                    #     H2[nz, :] = H2[nz, col_indices[i]]^-1 * H2[nz, :]
                    # end

                    for j in 1:nr
                        # go through all rows and eliminate all nonzeros in column using the chosen row
                        if j != nz && !iszero(H2[j, col_indices[i]])
                            if !isone(H2[j, col_indices[i]])
                                H2[j, :] = H2[j, :] - H2[j, :]^-1 * H2[nz, :]
                            else
                                H2[j, :] = H2[j, :] - H2[nz, :]
                            end
                        end
                    end

                    row_pivots[nz] = true
                    row_pivots_loc[nz] = i
                    if i != loc
                        temp = col_indices[loc]
                        col_indices[loc] = col_indices[i]
                        col_indices[i] = temp
                    end
                    have += 1
                    loc -= 1
                    used = true
                    break
                end
            end

            if !used
                if i != loc
                    temp = col_indices[loc]
                    col_indices[loc] = col_indices[i]
                    col_indices[i] = temp
                end
                loc -= 1
            end

            # if we need more than are available, restart
            if nr - have > loc
                col_indices = collect(1:nc)
                row_pivots = falses(nr)
                row_pivots_loc = zeros(Int, nr)
                have = 0
                loc = nc
                count += 1
                # H2 = deepcopy(H)
            end
        end

        # randomly split the remaning column indices
        X = Vector{Int}()
        Y = Vector{Int}()
        for i in col_indices[1:nc - nr]
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
                        # so any time there is a 1 at index i, put a 1 in the i-th position of the end of col_indices
                        # that will give w - 2p 1's at locations outside of A ∪ B
                        vec_len_w = zeros(C.F, 1, nc)
                        for k in 1:nr
                            if isone(πAB[k])
                                vec_len_w[1, row_pivots_loc[k]] = C.F(1)
                            end
                        end
                        # now add back in the 2p 1's at locations in A ∪ B to get a vector of weight w
                        for k in AB
                            vec_len_w[1, k] = C.F(1)
                        end
                        # println(vec_len_w)
                        test = matrix(C.F, 1, nc, [coeff(x, 0) for x in vec_len_w])
                        # println(iszero(H * transpose(test)))
                        if iszero(H * transpose(test))
                            println(test)
                            push!(len_w_vecs, test)
                            found += 1
                            if found == num_find
                                return len_w_vecs
                            end
                        end
                    end
                end
            end
        end
        count += 1
    end
    return len_w_vecs
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
# TODO: switch to symbols instead of strings
# see section 7.3 of White's thesis for comments on this
function _get_information_sets(G::CTMatrixTypes, alg::String)
    # alg ∈ ["Brouwer", "Zimmermann", "White", "Chen"] || throw(ArgumentError("Expected 'Brouwer', 'Zimmermann', or 'White'."))
    nr, nc = size(G)
    gen_mats = []
    # Gi = deepcopy(G)
    rnk = nr
    if alg == "Brouwer"
        for i in 0:Int(floor(nc / nr)) - 1
            rnk, Gi, Pi = CodingTheory._rref_col_swap(G, 1:nr, i * rnk + 1:nc)
            push!(gen_mats, (rnk, Gi, Pi))
            # Ai = Gi[:, setdiff(1:nc, i * nr + 1:(i + 1) * nr)]
            # push!(gen_mats, Ai)
        end
    elseif alg == "Zimmermann"
        for i in 0:Int(floor(nc / nr))
            rnk, Gi, Pi = CodingTheory._rref_col_swap(G, 1:nr, i * rnk + 1:nc)
            push!(gen_mats, (rnk, Gi, Pi))
            # display(Gi)
            # println(rnk)
            # Ai = Gi[:, setdiff(1:nc, i * nr + 1:i * nr + rnk)]
            # display(Ai)
            # println(" ")
            # push!(gen_mats, (rnk, Ai))
        end
    elseif alg == "White"
        # the expansion factor of the code
        for i in 0:div(nc, nr) - 1
            # could use Gi here instead of G
            rnk, Gi, Pi = CodingTheory._rref_col_swap(G, 1:nr, i * nr + 1:(i + 1) * nr)
            # display(Gi)
            # println(rnk)
            push!(gen_mats, (rnk, Gi, Pi))
        end
    elseif alg == "Chen"
        # should remove this in the future as this should simply be (k, Gstand, P)
        rnk, Gi, Pi = CodingTheory._rref_col_swap(G, 1:nr, 1:nc)
        push!(gen_mats, (rnk, Gi, Pi))
    end
    return gen_mats
end

# TODO: switch to symbols instead of strings
@inline function _lower_bounds(r::Int, n::Int, k::Int, l::Int, rank_defs::Vector{Int}, type::String)
    if type == "BZ"
        h = length(rank_defs)
        lower = 0
        for i in 1:h
            lower += maximum([0, r - rank_defs[i]])
        end
    elseif type == "Chen"
        lower = Int(ceil(n * r / k))
    elseif type == "White"
        lower = 0
        for i in 1:l
	        lower += Int(ceil(n * maximum([0, r - rank_defs[i]]) / (l * (k + rank_defs[i]))))
        end
    end
    return lower
end

"""
    Gray_code_min_dist(C::AbstractLinearCode, verbose::Bool=false)

Return the minimum distance of `C` using a deterministic algorithm based on enumerating
constant weight codewords of the binary reflected Gray code. If a word of minimum weight
is found before the lower and upper bounds cross, it is returned; otherwise, `missing`
is returned.
"""
# In the thesis, the h and GrayCode for loops are switched. This seems inefficient
# but does allow one to update the lower bound using a single matrix after every
# GrayCode loop has finished. I believe only enumerating the very expensive GrayCode
# once will overcome this. Could be tested though.
function Gray_code_min_dist(C::AbstractLinearCode, verbose::Bool=false)
    ord_F = Int(order(C.F))
    ord_F == 2 || throw(ArgumentError("Currently only implemented for binary codes."))

    p = Int(characteristic(C.F))
    found = missing
    perm = missing

    G = generator_matrix(C)
    if typeof(C) <: AbstractCyclicCode
        verbose && println("Detected a cyclic code, using Chen's adaption.")
        gen_mats = _get_information_sets(G, "Chen")
    elseif typeof(C) <: AbstractQuasiCyclicCode
        verbose && println("Detected a quasi-cyclic code, using White's adaption.")
        gen_mats = _get_information_sets(G, "White")
    else
        gen_mats = _get_information_sets(G, "Zimmermann")
    end
    gen_mats_Julia = [deepcopy(FpmattoJulia(x[2])') for x in gen_mats]
    h = length(gen_mats_Julia)
    rank_defs = zeros(Int, h)
    if verbose
        print("Generated $h information sets with ranks: ")
        for i in 1:h
            i == h ? (println(gen_mats[i][1]);) : (print("$(gen_mats[i][1]), "))
            # will only be using the rank deficits here
            # at the moment, the information sets are always disjoint so the relative
            # rank is zero
            rank_defs[i] = C.k - gen_mats[i][1]
        end
    end
    
    even_flag = false
    doubly_even_flag = false
    triply_even_flag = false
    ord_F == 2 && (even_flag = is_even(C);)
    even_flag && (doubly_even_flag = is_doubly_even(C);)
    doubly_even_flag && (triply_even_flag = is_triply_even(C);)
    if verbose
        triply_even_flag && println("Detected a triply even code.")
        (!triply_even_flag && doubly_even_flag) && println("Detected a doubly even code.")
        (!triply_even_flag && !doubly_even_flag && even_flag) && println("Detected an even code.")
    end
    
    flag = true
    while flag
        for (j, g) in enumerate(gen_mats_Julia)
            # redo j = 1 to grab index and perm
            w, i = _min_wt_col(g)
            if w <= C.u_bound
                found = g[:, i]
                C.u_bound = w
                perm = j
            end
        end
        # if restarted, found is missing
        ismissing(found) ? (C.u_bound = C.n;) : (flag = false;)
    end
    verbose && println("Upper bound after row analysis: $(C.u_bound)")
    verbose && !ismissing(found) && println("Found element matching upper bound.")

    num_thrds = Threads.nthreads()
    verbose && println("Detected $num_thrds threads.")
    power = 0
    for i in 1:20
        if 2^i > num_thrds
            power = i - 1
            break
        end
    end

    for r in 1:C.k
        if typeof(C) <: AbstractCyclicCode
            new_lower = _lower_bounds(r, C.n, C.k, 0, [0], "Chen")
            C.l_bound < new_lower && (C.l_bound = new_lower;)
        elseif typeof(C) <: AbstractQuasiCyclicCode
            new_lower = _lower_bounds(r, C.n, C.k, C.l, rank_defs, "White")
            C.l_bound < new_lower && (C.l_bound = new_lower;)
        else
            new_lower = _lower_bounds(r, C.n, C.k, 0, rank_defs, "BZ")
            C.l_bound < new_lower && (C.l_bound = new_lower;)
        end
        # an even code can't have have an odd minimum weight
        (!triply_even_flag && !doubly_even_flag && even_flag) && (C.l_bound += C.l_bound % 2;)
        (!triply_even_flag && doubly_even_flag) && (C.l_bound += 4 - C.l_bound % 4;)
        triply_even_flag && (C.l_bound += 8 - C.l_bound % 8;)

        verbose && println("r: $r")
        verbose && println("Lower bound: $(C.l_bound)")
        verbose && println("Upper bound: $(C.u_bound)")
        if C.l_bound >= C.u_bound
            C.d = C.u_bound
            if !ismissing(found)
                y = matrix(C.F, 1, C.n, found)
                ismissing(gen_mats[perm][3]) || (y = y * gen_mats[perm][3];)
                iszero(C.H * transpose(y)) ? (return C.u_bound, y;) : (return C.u_bound, missing;)
            else
                return C.u_bound, found
            end
        end

        flag = Threads.Atomic{Bool}(true)
        # num_thrds = 1
        # power = 0
        # total = Threads.Atomic{Int}(0)
        uppers = [C.u_bound for _ in 1:num_thrds]
        founds = [found for _ in 1:num_thrds]
        perms = [perm for _ in 1:num_thrds]
        Threads.@threads for m in 1:num_thrds
            c = zeros(Int, C.n)
            prefix = digits(m - 1, base=2, pad=power)
            # count = Threads.Atomic{Int}(0)
            for u in GrayCode(C.k, r, prefix, mutate=true)
                # count += 1
                # Threads.atomic_add!(count, 1)
                if flag[]
                    for i in 1:h
                        if r - rank_defs[i] > 0
                            LinearAlgebra.mul!(c, gen_mats_Julia[i], u)
                            w = 0
                            @inbounds for j in 1:C.n
                                c[j] % p != 0 && (w += 1;)
                            end

                            if uppers[m] > w
                                uppers[m] = w
                                founds[m] = c
                                perms[m] = i
                                # verbose && println("Adjusting upper bound: $upperlocal")
                                if C.l_bound == uppers[m]
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
        C.u_bound = uppers[loc]
        found = founds[loc]
        perm = perms[loc]
        # println("total: $total")
        if !flag[]
            C.d = C.u_bound
            if !ismissing(found)
                y = matrix(C.F, 1, C.n, found)
                ismissing(gen_mats[perm][3]) || (y = y * gen_mats[perm][3];)
                iszero(C.H * transpose(y)) ? (return C.u_bound, y;) : (return C.u_bound, missing;)
            else
                return C.u_bound, found
            end
        end
    end
end

"""
    words_of_weight(C::AbstractLinearCode, l_bound::Int, u_bound::Int, verbose::Bool=false)

Return all the codewords of `C` of Hamming weight in the range `[l_bound, u_bound]`.
"""
function words_of_weight(C::AbstractLinearCode, l_bound::Int, u_bound::Int, verbose::Bool=false)
    ord_F = Int(order(C.F))
    ord_F == 2 || throw(ArgumentError("Currently only implemented for binary codes."))

    1 <= l_bound <= u_bound <= C.n || throw(ArgumentError("Expected 1 <= l_bound <= u_bound <= C.n"))
    if l_bound < C.n / 2 && ord_F == 2
        # faster to enumerate backwards, but only in binary
        return _words_of_weight_high(C, l_bound, u_bound, verbose)
    end

    p = Int(characteristic(C.F))
    G = generator_matrix(C)
    if typeof(C) <: AbstractCyclicCode
        verbose && println("Detected a cyclic code, using Chen's adaption.")
        gen_mats = _get_information_sets(G, "Chen")
    elseif typeof(C) <: AbstractQuasiCyclicCode
        verbose && println("Detected a quasi-cyclic code, using White's adaption.")
        gen_mats = _get_information_sets(G, "White")
    else
        gen_mats = _get_information_sets(G, "Zimmermann")
    end
    gen_mats_Julia = [FpmattoJulia(x[2])' for x in gen_mats]
    h = length(gen_mats_Julia)
    rank_defs = zeros(Int, h)
    if verbose
        print("Generated $h information sets with ranks: ")
        for i in 1:h
            i == h ? (println(gen_mats[i][1]);) : (print("$(gen_mats[i][1]), "))
            # will only be using the rank deficits here
            # at the moment, the information sets are always disjoint so the relative
            # rank is zero
            rank_defs[i] = C.k - gen_mats[i][1]
        end
    end

    even_flag = false
    doubly_even_flag = false
    triply_even_flag = false
    ord_F == 2 && (even_flag = is_even(C);)
    even_flag && (doubly_even_flag = is_doubly_even(C);)
    doubly_even_flag && (triply_even_flag = is_triply_even(C);)
    if verbose
        triply_even_flag && println("Detected a triply even code.")
        (!triply_even_flag && doubly_even_flag) && println("Detected a doubly even code.")
        (!triply_even_flag && !doubly_even_flag && even_flag) && println("Detected an even code.")
    end

    num_thrds = Threads.nthreads()
    verbose && println("Detected $num_thrds threads.")
    power = 0
    for i in 1:20
        if 2^i > num_thrds
            power = i - 1
            break
        end
    end

    W = Set{typeof(G)}()
    for r in 1:C.k
        if typeof(C) <: AbstractCyclicCode
            lower = _lower_bounds(r, C.n, C.k, 0, [0], "Chen")
        elseif typeof(C) <: AbstractQuasiCyclicCode
            lower = _lower_bounds(r, C.n, C.k, C.l, rank_defs, "White")
        else
            lower = _lower_bounds(r, C.n, C.k, 0, rank_defs, "BZ")
        end
        # an even code can't have have an odd minimum weight
        (!triply_even_flag && !doubly_even_flag && even_flag) && (lower += lower % 2;)
        (!triply_even_flag && doubly_even_flag) && (lower += 4 - lower % 4;)
        triply_even_flag && (lower += 8 - lower % 8;)

        verbose && println("r: $r")
        verbose && println("Lower bound: $lower")
        if lower >= u_bound
            return W
        end

        Ws = [Set{typeof(G)}() for _ in 1:num_thrds]
        Threads.@threads for m in 1:num_thrds
            c = zeros(Int, C.n)
            prefix = digits(m - 1, base=2, pad=power)
            for u in GrayCode(C.k, r, prefix, mutate=true)
                for i in 1:h
                    LinearAlgebra.mul!(c, gen_mats_Julia[i], u)
                    w = 0
                    @inbounds for j in 1:C.n
                        c[j] % p != 0 && (w += 1;)
                    end

                    if l_bound <= w <= u_bound
                        c2 = matrix(C.F, 1, C.n, c)
                        ismissing(gen_mats[i][3]) || (c2 = c2 * gen_mats[i][3];)
                        push!(Ws[m], c2)
                    end
                end
            end
        end
        for m in 1:num_thrds
            union!(W, Ws[m])
        end
    end
end

"""
    words_of_weight(C::AbstractLinearCode, bound::Int, verbose::Bool=false)

Return all the codewords of `C` of Hamming weight `bound`.
"""
words_of_weight(C::AbstractLinearCode, bound::Int, verbose::Bool=false) = words_of_weight(C, bound, bound, verbose)

# untested
# TODO: figure out if even weight upper weight needs to subtract or not
function _words_of_weight_high(C::AbstractLinearCode, l_bound::Int, u_bound::Int, verbose::Bool=false)
    p = Int(characteristic(C.F))
    G = generator_matrix(C)
    if typeof(C) <: AbstractCyclicCode
        verbose && println("Detected a cyclic code, using Chen's adaption.")
        gen_mats = _get_information_sets(G, "Chen")
    elseif typeof(C) <: AbstractQuasiCyclicCode
        verbose && println("Detected a quasi-cyclic code, using White's adaption.")
        gen_mats = _get_information_sets(G, "White")
    else
        gen_mats = _get_information_sets(G, "Zimmermann")
    end
    gen_mats_Julia = [FpmattoJulia(x[2])' for x in gen_mats]
    h = length(gen_mats_Julia)
    rank_defs = zeros(Int, h)
    if verbose
        print("Generated $h information sets with ranks: ")
        for i in 1:h
            i == h ? (println(gen_mats[i][1]);) : (print("$(gen_mats[i][1]), "))
            # will only be using the rank deficits here
            # at the moment, the information sets are always disjoint so the relative
            # rank is zero
            rank_defs[i] = C.k - gen_mats[i][1]
        end
    end

    even_flag = false
    doubly_even_flag = false
    triply_even_flag = false
    ord_F == 2 && (even_flag = is_even(C);)
    even_flag && (doubly_even_flag = is_doubly_even(C);)
    doubly_even_flag && (triply_even_flag = is_triply_even(C);)
    if verbose
        triply_even_flag && println("Detected a triply even code.")
        (!triply_even_flag && doubly_even_flag) && println("Detected a doubly even code.")
        (!triply_even_flag && !doubly_even_flag && even_flag) && println("Detected an even code.")
    end

    num_thrds = Threads.nthreads()
    verbose && println("Detected $num_thrds threads.")
    power = 0
    for i in 1:20
        if 2^i > num_thrds
            power = i - 1
            break
        end
    end

    W = Set{typeof(G)}()
    for r in C.k:-1:1
        if typeof(C) <: AbstractCyclicCode
            upper = _lower_bounds(r, C.n, C.k, 0, [0], "Chen")
        elseif typeof(C) <: AbstractQuasiCyclicCode
            upper = _lower_bounds(r, C.n, C.k, C.l, rank_defs, "White")
        else
            upper = _lower_bounds(r, C.n, C.k, 0, rank_defs, "BZ")
        end
        # # an even code can't have have an odd minimum weight
        # (!triply_even_flag && !doubly_even_flag && even_flag) && (upper += upper % 2;)
        # (!triply_even_flag && doubly_even_flag) && (upper += 4 - upper % 4;)
        # triply_even_flag && (upper += 8 - upper % 8;)

        verbose && println("r: $r")
        verbose && println("Upper bound: $upper")
        if upper < l_bound
            return W
        end

        Ws = [Set{typeof(G)}() for _ in 1:num_thrds]
        Threads.@threads for m in 1:num_thrds
            c = zeros(Int, C.n)
            prefix = digits(m - 1, base=2, pad=power)
            for u in GrayCode(C.k, r, prefix, mutate=true)
                for i in 1:h
                    LinearAlgebra.mul!(c, gen_mats_Julia[i], u)
                    w = 0
                    @inbounds for j in 1:C.n
                        c[j] % p != 0 && (w += 1;)
                    end

                    if l_bound <= w <= u_bound
                        c2 = matrix(C.F, 1, C.n, c)
                        ismissing(gen_mats[i][3]) || (c2 = c2 * gen_mats[i][3];)
                        push!(Ws[m], c2)
                    end
                end
            end
        end
        for m in 1:num_thrds
            union!(W, Ws[m])
        end
    end
end

#TODO: change above to BigInts
"""
    partial_weight_distribution(C::AbstractLinearCode, bound::Int, compact::Bool=false)

Return the partial weight distribution of `C` up to weight `bound`. If `compact` is false,
the result will be a `Vector{BigInt}` of length `length(C) + 1` whose `i`th entry is the
number of codewords of `C` of Hamming weight `i - 1`. Otherwise, the result is a
`Vector{Tuple{Int, BigInt}}` whose entries specify the nonzero indices and values of the
above.
"""
# TODO: type stability of this
function partial_weight_distribution(C::AbstractLinearCode, bound::Int, compact::Bool=false)
	1 <= bound <= C.n || throw(ArgumentError("Bound must be between 1 and n."))

	W = words_of_weight(C, 1, bound)
    wt_dist = zeros(BigInt, bound + 1)
    bio = BigInt(1)
    wt_dist[1] = bio    
	for c in W
		wt_dist[wt(c) + 1] += bio
	end
    
    if compact
        wt_dist_comp = Vector{Tuple{Int, BigInt}}()
        for (i, x) in enumerate(wt_dist)
            !iszero(x) && (push!(wt_dist_comp, (i - 1, x)))
        end
        return wt_dist_comp
    else
        return wt_dist
    end
end

"""
	minimum_words(C::AbstractLinearCode)

Return the set of codewords of `C` with weight equal to the minimum distance.

# Notes
* This algorithm simultaneously computes the minimum distance and stores the words of
  this weight that it finds, removing the repeated work of calling
  `w = minimum_distance(C); W = words_of_weight(C, w);`
"""
function minimum_words(C::AbstractLinearCode)
    ord_F = Int(order(C.F))
    ord_F == 2 || throw(ArgumentError("Currently only implemented for binary codes."))

    p = Int(characteristic(C.F))
    found = missing
    perm = missing

    G = generator_matrix(C)
    if typeof(C) <: AbstractCyclicCode
        verbose && println("Detected a cyclic code, using Chen's adaption.")
        gen_mats = _get_information_sets(G, "Chen")
    elseif typeof(C) <: AbstractQuasiCyclicCode
        verbose && println("Detected a quasi-cyclic code, using White's adaption.")
        gen_mats = _get_information_sets(G, "White")
    else
        gen_mats = _get_information_sets(G, "Zimmermann")
    end
    gen_mats_Julia = [FpmattoJulia(x[2])' for x in gen_mats] # deepcopy doesn't save anythign here
    # gen_mats_Julia = [FpmattoJulia(x[2]) for x in gen_mats]
    h = length(gen_mats_Julia)
    rank_defs = zeros(Int, h)
    if verbose
        print("Generated $h information sets with ranks: ")
        for i in 1:h
            i == h ? (println(gen_mats[i][1]);) : (print("$(gen_mats[i][1]), "))
            # will only be using the rank deficits here
            # at the moment, the information sets are always disjoint so the relative
            # rank is zero
            rank_defs[i] = C.k - gen_mats[i][1]
        end
    end
    
    even_flag = false
    doubly_even_flag = false
    triply_even_flag = false
    ord_F == 2 && (even_flag = is_even(C);)
    even_flag && (doubly_even_flag = is_doubly_even(C);)
    doubly_even_flag && (triply_even_flag = is_triply_even(C);)
    if verbose
        triply_even_flag && println("Detected a triply even code.")
        (!triply_even_flag && doubly_even_flag) && println("Detected a doubly even code.")
        (!triply_even_flag && !doubly_even_flag && even_flag) && println("Detected an even code.")
    end
    
    upper = C.n - C.k + 1
    verbose && println("Singleton upper bound: $upper")
    for (j, g) in enumerate(gen_mats_Julia)
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

    num_thrds = Threads.nthreads()
    verbose && println("Detected $num_thrds threads.")
    power = 0
    for i in 1:20
        if 2^i > num_thrds
            power = i - 1
            break
        end
    end

    W = Set{typeof(G)}()
    for r in 1:C.k
        if typeof(C) <: AbstractCyclicCode
            lower = _lower_bounds(r, C.n, C.k, 0, [0], "Chen")
        elseif typeof(C) <: AbstractQuasiCyclicCode
            lower = _lower_bounds(r, C.n, C.k, C.l, rank_defs, "White")
        else
            lower = _lower_bounds(r, C.n, C.k, 0, rank_defs, "BZ")
        end
        # an even code can't have have an odd minimum weight
        (!triply_even_flag && !doubly_even_flag && even_flag) && (lower += lower % 2;)
        (!triply_even_flag && doubly_even_flag) && (lower += 4 - lower % 4;)
        triply_even_flag && (lower += 8 - lower % 8;)

        verbose && println("r: $r")
        verbose && println("Lower bound: $lower")
        verbose && println("Upper bound: $upper")
        if lower >= upper
            C.d = upper
            return upper, W
        end

        uppers = [upper for _ in 1:num_thrds]
        Ws = [Set{typeof(G)}() for _ in 1:num_thrds]
        Threads.@threads for m in 1:num_thrds
            c = zeros(Int, C.n)
            prefix = digits(m - 1, base=2, pad=power)
            for u in GrayCode(C.k, r, prefix, mutate=true)
                for i in 1:h
                    LinearAlgebra.mul!(c, gen_mats_Julia[i], u)
                    w = 0
                    @inbounds for j in 1:C.n
                        c[j] % p != 0 && (w += 1;)
                    end

                    if w <= uppers[m]
                        if w < uppers[m]
                            uppers[m] = w
                            verbose && println("Adjusting upper bound: $upper")
                            Ws[m] = Set{typeof(G)}()
                        end
                        # TODO: this is very expensive just to get erased
                        # maybe keep perms[m] and adjust outside loop
                        # allocations are more expensive inside threads
                        c2 = matrix(C.F, 1, C.n, c)
                        ismissing(gen_mats[i][3]) || (c2 = c2 * gen_mats[i][3];)
                        push!(Ws[m], c2)
                    end
                end
            end
        end
        loc = argmin(uppers)
        if upper > uppers[loc]
            upper = uppers[loc]
            W = Set{typeof(G)}()
        end
        for m in 1:num_thrds
            if uppers[m] == upper
                union!(W, Ws[m])
            end
        end
    end
end

#############################
    # Weight Enumerators
#############################

# TODO: name here, C is a placeholder for classical
function weight_enumerator_C(T::Trellis, type::Symbol=:complete)
    type ∈ [:complete, :Hamming] ||
        throw(ArgumentError("Unsupported weight enumerator type '$type'. Expected ':complete' or ':Hamming'."))

    if type == :complete && !ismissing(T.CWE)
        return T.CWE
    elseif type == :Hamming && !ismissing(T.CWE)
        return CWE_to_HWE(T.CWE)
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
        T.code.weight_enum = T.CWE
        HWE = CWE_to_HWE(T.CWE)
        T.code.d = minimum([collect(exponent_vectors(polynomial(HWE)))[i][1]
            for i in 1:length(polynomial(HWE))])
    end

    # clean up vertices
    for i in 1:length(V)
        for v in V[i]
            v.polynomial = missing
        end
    end

    type == :Hamming && return CWE_to_HWE(T.CWE)
    return T.CWE
end

"""
    MacWilliams_identity(C::AbstractLinearCode, W::WeightEnumerator, dual::Symbol=:Euclidean)

Return the weight enumerator of the dual (`:Euclidean` or `:Hermitian`) of `C` obtained
by applying the MacWilliams identities to `W`.
"""
# TODO: remove C from this, store in WE struct
function MacWilliams_identity(C::AbstractLinearCode, W::WeightEnumerator, dual::Symbol=:Euclidean)
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
        K, ζ = CyclotomicField(3, :ζ)
        R, vars = PolynomialRing(K, 3)
        # might have to switch this here
        poly = divexact(W.polynomial(
            vars[1] + vars[2] + vars[3],
            vars[1] + ζ * vars[2] + ζ^2 * vars[3],
            vars[1] + ζ^2 * vars[2] + ζ * vars[3]), cardinality(C))
        # works so far but now needs to recast down to the integer ring
        return WeightEnumerator(Oscar.map_coefficients(c -> Nemo.ZZ(coeff(c, 0)), poly,
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
            K, ω = CyclotomicField(Int(characteristic(C.F)), :ω)
            R, vars = PolynomialRing(K, q)
            elms = collect(C.F)
            func_args = []
            for i in 1:q
                inner_sum = R(0)
                for j in 1:q
                    inner_sum += ω^coeff(elms[i] * elms[j], 0) * vars[j]
                end
                append!(func_args, inner_sum)
            end
            return WeightEnumerator(divexact(W.polynomial(func_args), cardinality(C)),
                :complete)
        else
            K, ω = CyclotomicField(Int(characteristic(C.F)), :ω)
            R, vars = PolynomialRing(K, q)
            prime_field = GF(Int(characteristic(C.F)))
            _, λ = primitivebasis(C.F, prime_field)
            elms = collect(C.F)
            func_args = []
            for i in 1:q
                inner_sum = R(0)
                for j in 1:q
                    β = elms[i] * elms[j]
                    β_exp = _expandelement(β, prime_field, λ, false)
                    inner_sum += ω^coeff(β_exp[1], 0) * vars[j]
                end
                push!(func_args, inner_sum)
            end
            display(func_args)
            return WeightEnumerator(divexact(W.polynomial(func_args...), cardinality(C)), :complete)
        end
    end
end

"""
    weight_enumerator(C::AbstractLinearCode, type::Symbol=:complete, alg::String="auto")

Return either the `:complete` or `:Hamming` weight enumerator of `C` using the algorithm `alg`.
"""
function weight_enumerator(C::AbstractLinearCode, type::Symbol=:complete, alg::String="auto")
    type ∈ [:complete, :Hamming] ||
        throw(ArgumentError("Unsupported weight enumerator type '$type'. Expected ':complete' or ':Hamming'."))
    alg ∈ ["auto", "trellis", "bruteforce"] ||
        throw(ArgumentError("Algorithm `$alg` is not implemented in weight_enumerator."))

    if type == :complete && !ismissing(C.weight_enum)
        return C.weight_enum
    elseif type == :Hamming && !ismissing(C.weight_enum)
        return CWE_to_HWE(C.weight_enum)
    end

    if alg == "auto"
        if cardinality(C) <= 1e6 # random cutoff
            C.weight_enum = _weight_enumerator_BF(C.G)
            HWE = CWE_to_HWE(C.weight_enum)
            C.d = minimum(filter(x->x!=0, [collect(exponent_vectors(HWE.polynomial))[i][1]
                for i in 1:length(HWE.polynomial)]))
            type == :Hamming && return HWE
            return C.weight_enum
        elseif rate(C) > 0.5
            D = dual(C)
            if cardinality(D) <= 1e6 # random cutoff
                D.weight_enum = _weight_enumerator_BF(D.G)
            else
                weight_enumerator_C(syndrome_trellis(D, "primal", false), type)
            end
            C.weight_enum = MacWilliams_identity(D, D.weight_enum)
            HWE = CWE_to_HWE(C.weight_enum)
            C.d = minimum(filter(x->x!=0, [collect(exponent_vectors(HWE.polynomial))[i][1]
                for i in 1:length(HWE.polynomial)]))
            type == :Hamming && return HWE
            return C.weight_enum
        else
            return weight_enumerator_C(syndrome_trellis(C, "primal", false), type)
        end
    elseif alg == "trellis"
        return weight_enumerator_C(syndrome_trellis(C, "primal", false), type)
    elseif alg == "bruteforce"
        C.weight_enum = _weight_enumerator_BF(C.G)
        HWE = CWE_to_HWE(C.weight_enum)
        C.d = minimum(filter(x->x!=0, [collect(exponent_vectors(HWE.polynomial))[i][1]
            for i in 1:length(HWE.polynomial)]))
        type == :Hamming && return HWE
        return C.weight_enum
    end
end

"""
    weight_distribution(C::AbstractLinearCode, alg::String="auto", compact::Bool=true)

Return the weight distribution of `C` using the algorithm `alg`. If `compact` is false,
the result will be a `Vector{BigInt}` of length `length(C) + 1` whose `i`th entry is the
number of codewords of `C` of Hamming weight `i - 1`. Otherwise, the result is a
`Vector{Tuple{Int, BigInt}}` whose entries specify the nonzero indices and values of the
above.
"""
function weight_distribution(C::AbstractLinearCode, alg::String="auto", compact::Bool=true)
    alg ∈ ["auto", "trellis", "bruteforce"] ||
        throw(ArgumentError("Algorithm `$alg` is not implemented in weight_enumerator."))

    ismissing(C.weight_enum) && weight_enumerator(C, :complete, alg)
    HWE = CWE_to_HWE(C.weight_enum)

    if compact
        wt_dist = Vector{Tuple}()
        for i in 1:length(HWE.polynomial)
            push!(wt_dist, (exponent_vector(HWE.polynomial, i)[1],
                coeff(HWE.polynomial, i)))
        end
    else
        wt_dist = zeros(Int, 1, C.n + 1)
        for i in 1:length(HWE.polynomial)
            wt_dist[exponent_vector(HWE.polynomial, i)[1] + 1] = coeff(HWE.polynomial, i)
        end
    end
    return wt_dist
end

"""
    weight_plot(C::AbstractLinearCode, alg::String="auto")

Return a bar plot of the weight distribution of `C`.
"""
function weight_plot(C::AbstractLinearCode, alg::String="auto")
    wt_dist = weight_distribution(C, alg, true)
    x_ticks = findall(x->x>0, vec(wt_dist)) .- 1
    y_ticks = [wt_dist[i] for i in 1:length(wt_dist) if !iszero(wt_dist[i])]
    ismissing(C.d) ? (title="Weight Distribution - [$(C.n), $(C.k)]";) :
        title="Weight Distribution - [$(C.n), $(C.k), $(C.d)]"
    f = bar(0:C.n, wt_dist', bar_width=1, xticks=x_ticks, yticks=y_ticks,
        legend=false, xlabel="Weight", ylabel="Number of Terms", title=title)
    display(f)
    return f
end

"""
    support(C::AbstractLinearCode)

Returns the support of `C`.

# Notes
* The support of `C` is the collection of nonzero exponents of the Hamming weight enumerator of `C`.
"""
support(C::AbstractLinearCode) = [i for (i, _) in weight_distribution(C, "auto", true)]

#############################
     # Minimum Distance
#############################

# TODO: update doc strings for these and this whole file in general
"""
    minimum_distance(C::AbstractLinearCode, alg::String="trellis", sect::Bool=false)

Return the minimum distance of the linear code if known, otherwise computes it
using the algorithm of `alg`. If `alg = "trellis"`, the sectionalization flag
`sect` can be set to true to further compactify the reprsentation.
"""
function minimum_distance(C::AbstractLinearCode, alg::String="auto", sect::Bool=false, verbose::Bool=false)
    !ismissing(C.d) && return C.d

    alg ∈ ["auto", "Gray", "trellis", "Leon", "bruteforce", "wt_dist"] ||
        throw(ArgumentError("Unexpected algorithm '$alg'."))
    
    if alg == "auto"
        D = dual(C)
        if cardinality(C) <= 1e6 # random cutoff
            C.weight_enum = _weight_enumerator_BF(C.G)
            HWE = CWE_to_HWE(C.weight_enum)
            C.d = minimum(filter(x->x!=0, [collect(exponent_vectors(HWE.polynomial))[i][1]
                for i in 1:length(HWE.polynomial)]))
            return C.d
        elseif rate(C) > 0.5
            D = dual(C)
            if cardinality(D) <= 1e6 # random cutoff
                D.weight_enum = _weight_enumerator_BF(D.G)
            else
                weight_enumerator_C(syndrome_trellis(D, "primal", false), type)
            end
            C.weight_enum = MacWilliams_identity(D, D.weight_enum)
            HWE = CWE_to_HWE(C.weight_enum)
            C.d = minimum(filter(x->x!=0, [collect(exponent_vectors(HWE.polynomial))[i][1]
                for i in 1:length(HWE.polynomial)]))
            return C.d
        else
            return Gray_code_min_dist(C, verbose)
        end
    elseif alg == "Gray"
        return Gray_code_min_dist(C, verbose)
    elseif alg == "trellis"
        weight_enumerator_C(syndrome_trellis(C, "primal", false), type)
        return C.d
    elseif alg == "bruteforce"
        C.weight_enum = _weight_enumerator_BF(C.G)
        HWE = CWE_to_HWE(C.weight_enum)
        C.d = minimum(filter(x->x!=0, [collect(exponent_vectors(HWE.polynomial))[i][1]
            for i in 1:length(HWE.polynomial)]))
        return C.d
    elseif alg == "wt_dist"
        HWE = weight_enumerator(C, :Hamming, alg)
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
# TODO: remove quadratic extension
function _weight_enumerator_BF_Q(G::CTMatrixTypes, char_vec::Vector{nmod},
    R::Union{AbstractAlgebra.Generic.MPolyRing{nf_elem}, Missing})
    # this should be the quadratic extension field
    E = base_ring(G)
    is_even(Int(degree(E))) || error("Matrix passed to weight enumerator does not appear to be over the quadratic extension.")
    ord_E = Int(order(E))
    lookup = Dict(value => key for (key, value) in enumerate(collect(E)))
    
    p = Int(characteristic(E))
    is_even(p) ? nth = 2 * p : nth = p
    if ismissing(R)
        K, ω = CyclotomicField(nth, :ω)
        R, vars = PolynomialRing(K, ord_E)
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
        row_sym = quadratic_to_symplectic(row)

        # to do process signs here
        parity = 0
        for c in 1:2 * nc
            iszero(row_sym[c]) || (parity += data(char_vec[c]);)
        end

        # TODO: can do this in one step, but is it faster?
        term = zeros(Int, 1, ord_E)
        for x in row
            term[lookup[x]] += 1
        end
        # println(term, ", ", typeof(term))
        term_poly = ω^parity
        for i in 1:ord_E
            term_poly *= vars[i]^term[i]
        end
        poly += term_poly
    end
    # display(poly)
    return WeightEnumerator(poly, :complete)
    # return poly
end

# formulas from
# "Weight enumerators for nonbinary asymmetric quantum codes and their applications"
# by Chuangqiang Hu, Shudi Yang, Stephen S.-T.Yau
function MacWilliams_identity(S::AbstractStabilizerCode, W::WeightEnumerator, dual::Bool=false)
    dual ? (card = BigInt(characteristic(S.F))^(S.n + S.k);) : (card = cardinality(S);)
    if W.type == :Hamming
        # (1/(q^n|S|))W(y - x, y + (q^2 - 1)x)
        R = parent(W.polynomial)
        vars = gens(R)
        q = Int(order(S.F))
        return WeightEnumerator(divexact(W.polynomial(vars[2] - vars[1], vars[2] +
            (q^2 - 1) * vars[1]), card), :Hamming)
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
            card), :complete) # need the /2 to connect to the original Shor-Laflamme def
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
        # func_args = []
        # for i in 1:q
        #     inner_sum = R(0)
        #     for j in 1:q
        #         # inner_sum += ω^(2*tr(coeff(elms[i], 0) * coeff(elms[j], 1) -
        #         #     coeff(elms[j], 0) * coeff(elms[i], 1))) * vars[j]
        #         inner_sum += ω^(2*tr(coeff(basis[i], 0) * coeff(basis[j], 1) -
        #             coeff(basis[j], 0) * coeff(basis[i], 1))) * vars[perm[j]]
        #     end
        #     push!(func_args, inner_sum) # /q for Shor-Laflamme
        # end
        # println(basis)
        # println(elms)
        # println(perm)
        # display(func_args)
        # display(func_args[perm])
        # return WeightEnumerator(divexact(W.polynomial(func_args[perm]...), card),
        #     "complete")
    end
end

# TODO: switch to symbols
function weight_enumerator(S::AbstractStabilizerCode, type::Symbol=:complete,
    alg::String="auto", set::String="all")

    type ∈ [:complete, :Hamming] || error("Unsupported weight enumerator type '$type'. Expected ':complete' or ':Hamming'.")
    alg ∈ ["auto", "trellis", "bruteforce"] || error("Algorithm `$alg` is not implemented in weight_enumerator.")
    set ∈ ["all", "stabilizers", "logicals", "quotient"] || throw(ArgumentError("Unsupported set type '$set'. Expected 'all', 'stabilizers', 'logicals', 'quotient'."))

    if set ∈ ["all", "logicals"] && ismissing(S.sgn_CWE_logs)
        logs_mat = logicals_matrix(S)
        S.sgn_CWE_logs = _weight_enumerator_BF_Q(logs_mat, S.char_vec, missing)
    end

    if set != "logicals" && ismissing(S.sgn_CWE_stabs)
        if alg == "bruteforce" || cardinality(S) <= 1e6
            S.sgn_CWE_stabs = _weight_enumerator_BF_Q(S.stabs, S.char_vec, parent(S.sgn_CWE_logs.polynomial))
        else
            # trellis solution here
        end
    end

    if set ∈ ["all", "quotient"] && ismissing(S.sgn_CWE_dual)
        if alg == "bruteforce" || BigInt(characteristic(S.F))^(S.n + S.k) <= 3e6
            S.sgn_CWE_dual = _weight_enumerator_BF_Q(vcat(S.stabs, logicals_matrix(S)), S.char_vec, parent(S.sgn_CWE_logs.polynomial))
        else
            # trellis solution here
        end
    end
    
    if !ismissing(S.sgn_CWE_stabs) && !ismissing(S.sgn_CWE_dual)
        # compute minimum distance here
        poly = WeightEnumerator(S.sgn_CWE_dual.polynomial - S.sgn_CWE_stabs.polynomial, :complete)
        HWE = CWE_to_HWE(poly)
        S.d = minimum(filter(x->x!=0, [collect(exponent_vectors(HWE.polynomial))[i][1]
            for i in 1:length(HWE.polynomial)]))
    end

    if type == :complete
        set == "all" && return S.sgn_CWE_stabs, S.sgn_CWE_dual, S.sgn_CWE_logs, poly
        set == "stabilizers" && return S.sgn_CWE_stabs
        set == "logicals" && return S.sgn_CWE_logs
        return poly
    else
        set == "all" && return CWE_to_HWE(S.sgn_CWE_stabs), CWE_to_HWE(S.sgn_CWE_dual), CWE_to_HWE(S.sgn_CWE_logs), CWE_to_HWE(poly)
        set == "stabilizers" && return CWE_to_HWE(S.sgn_CWE_stabs)
        set == "logicals" && return CWE_to_HWE(S.sgn_CWE_logs)
        return HWE
    end
end

# MAGMA returns this format
# [ <0, 1>, <4, 105>, <6, 280>, <8, 435>, <10, 168>, <12, 35> ]
function weight_distribution(S::AbstractStabilizerCode, alg::String="auto", compact::Bool=true, set::String="all")
    alg ∈ ["auto", "trellis", "bruteforce"] || error("Algorithm `$alg` is not implemented in weight_enumerator.")
    set ∈ ["all", "stabilizers", "logicals", "quotient"] || throw(ArgumentError("Unsupported set type '$set'. Expected 'all', 'stabilizers', 'logicals', 'quotient'."))

    wt_enums = weight_enumerator(S, :Hamming, alg, set)

    if compact
        if length(wt_enums) == 1
            wt_dist = Vector{Tuple}()
            for i in 1:length(wt_enums.polynomial)
                push!(wt_dist, (exponent_vector(wt_enums.polynomial, i)[1],
                    coeff(wt_enums.polynomial, i)))
            end
        else
            wt_dist = Vector{Vector{Tuple}}()
            for wt_enum in wt_enums
                wt_dist_inner = Vector{Tuple}()
                for i in 1:length(wt_enum.polynomial)
                    push!(wt_dist_inner, (exponent_vector(wt_enum.polynomial, i)[1],
                        coeff(wt_enum.polynomial, i)))
                end
                push!(wt_dist, wt_dist_inner)
            end
        end
    else
        if length(wt_enums) == 1
            K = base_ring(wt_enums.polynomial)
            wt_dist = zero_matrix(K, 1, S.n + 1)
            for i in 1:length(wt_enums.polynomial)
                wt_dist[1, exponent_vector(wt_enums.polynomial, i)[1] + 1] = coeff(wt_enums.polynomial, i)
            end
        else
            K = base_ring(wt_enums[1].polynomial)
            wt_dist = [] #Vector{Vector{K}}()
            for wt_enum in wt_enums
                wt_dist_inner = zero_matrix(K, 1, S.n + 1)
                for i in 1:length(wt_enum.polynomial)
                    # println(coeff(wt_enum.polynomial, i))
                    wt_dist_inner[1, exponent_vector(wt_enum.polynomial, i)[1] + 1] = coeff(wt_enum.polynomial, i)
                end
                push!(wt_dist, wt_dist_inner)
            end
        end
    end
    return wt_dist
end

function weight_enumerator_Q(T::Trellis, type::Symbol=:complete)
    type ∈ [:complete, :Hamming] || error("Unsupported weight enumerator type '$type'. Expected ':complete' or ':Hamming'.")

    if type == :complete && !ismissing(T.CWE)
        return T.CWE
    elseif type == :Hamming && !ismissing(T.CWE)
        return CWE_to_HWE(T.CWE)
    end

    # if this ever changes or permutes will have to store with T
    elms = collect(T.code.E)
    lookup = Dict(value => key for (key, value) in enumerate(elms))

    p = Int(characteristic(T.code.E))
    is_even(p) ? nth = 2 * p : nth = p
    K, ω = CyclotomicField(nth, :ω)
    R, vars = PolynomialRing(K, length(elms))

    n = T.code.n
    char_vec = T.code.char_vec
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
                inner_bit = deepcopy(bit)
                parity = 0
                for k in e.label
                    if !iszero(coeff(k, 0))
                        parity += data(char_vec[inner_bit])
                    end
                    if !iszero(coeff(k, 1))
                        parity += data(char_vec[inner_bit + n])
                    end
                    inner_bit += 1
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
    T.CWE = WeightEnumerator(V[end][1].polynomial, :complete)

    # # currently Missing is not an option but how to implement dual trellis
    if !isshifted(T) && !ismissing(T.code)
        T.code.sgn_CWE_stabs = T.CWE
    end

    # # clean up vertices
    # for i in 1:length(V)
    #     for v in V[i]
    #         v.polynomial = missing
    #     end
    # end

    # display(T.CWE.polynomial)
    if type == :Hamming
        return CWE_to_HWE(T.CWE)
    end
    return T.CWE
end

"""
    weight_plot(S::AbstractStabilizerCode, alg::String="auto", type::String="stabilizer")

Return a bar plot of the weight distribution related to `S`.

If `type` is `stabilizer`, the weight distribution of the stabilizers are computed.
If `type` is `normalizer`, the weight distrbution of the normalizer of the stabilizers
are computed. If `type` is `quotient`, the weight distrbution of the normalizer mod the
stabilizers (logical representatives only) is computed.
"""
function weight_plot(S::AbstractStabilizerCode, alg::String="auto", type::String="stabilizer")
    type ∈ ["stabilizer", "normalizer", "quotient"] || throw(ArgumentError("Unknown value $type for parameter type."))

    wt_dist = weight_distribution(S, alg, type, false)
    x_ticks = findall(x -> x > 0, vec(wt_dist)) .- 1
    y_ticks = [wt_dist[i] for i in 1:length(wt_dist) if !iszero(wt_dist[i])]
    if type == "stabilizer"
        title_str = "Stabilizer Weight Distribution"
    elseif type == "normalizer"
        title_str = "Normalizer Weight Distribution"
    else
        title_str = "Quotient Weight Distribution"
    end
    ismissing(S.d) ? (title="$title_str - [$(S.n), $(S.k)]";) :
        title="$title_str - [$(S.n), $(S.k), $(S.d)]"
    f = bar(0:S.n, wt_dist', bar_width=1, xticks=x_ticks, yticks=y_ticks,
        legend=false, xlabel="Weight", ylabel="Number of Terms", title=title)
    display(f)
    return f
end

"""
    weight_plot_CSS_X(S::AbstractStabilizerCodeCSS, alg::String="auto")

Return a bar plot of the weight distribution of the `X` stabilizers.
"""
function weight_plot_CSS_X(S::AbstractStabilizerCodeCSS, alg::String="auto")
    C = LinearCode(S.X_stabs)
    wt_dist = weight_distribution(C, alg, false)
    x_ticks = findall(x->x>0, vec(wt_dist)) .- 1
    y_ticks = [wt_dist[i] for i in 1:length(wt_dist) if !iszero(wt_dist[i])]
    f = bar(0:C.n, wt_dist', bar_width=1, xticks=x_ticks, yticks=y_ticks,
        legend=false, xlabel="Weight", ylabel="Number of Terms",
        title="X-Weight Distribution")
    display(f)
    return f
end

"""
    weight_plot_CSS_Z(S::AbstractStabilizerCodeCSS, alg::String="auto")

Return a bar plot of the weight distribution of the `Z` stabilizers.
"""
function weight_plot_CSS_Z(S::AbstractStabilizerCodeCSS, alg::String="auto")
    C = LinearCode(S.Z_stabs)
    wt_dist = weight_distribution(C, alg, false)
    x_ticks = findall(x->x>0, vec(wt_dist)) .- 1
    y_ticks = [wt_dist[i] for i in 1:length(wt_dist) if !iszero(wt_dist[i])]
    f = bar(0:C.n, wt_dist', bar_width=1, xticks=x_ticks, yticks=y_ticks,
        legend=false, xlabel="Weight", ylabel="Number of Terms",
        title="Z-Weight Distribution")
    display(f)
    return f

end

"""
    weight_plot_CSS(S::AbstractStabilizerCodeCSS, alg::String="auto")

Return bar plots of the weight distribution of the both the
`X` and 'Z' stabilizers, separately.
"""
function weight_plot_CSS(S::AbstractStabilizerCodeCSS, alg::String="auto")
    C = LinearCode(S.X_stabs)
    wt_dist = weight_distribution(C, alg, false)
    x_ticks = findall(x->x>0, vec(wt_dist)) .- 1
    y_ticks = [wt_dist[i] for i in 1:length(wt_dist) if !iszero(wt_dist[i])]
    f_X = bar(0:C.n, wt_dist', bar_width=1, xticks=x_ticks, yticks=y_ticks,
        legend=false, xlabel="Weight", ylabel="Number of Terms",
        title="X-Weight Distribution")

    # okay to overwrite
    C = LinearCode(S.Z_stabs)
    wt_dist = weight_distribution(C, alg, false)
    x_ticks = findall(x->x>0, vec(wt_dist)) .- 1
    y_ticks = [wt_dist[i] for i in 1:length(wt_dist) if !iszero(wt_dist[i])]
    f_Z = bar(0:C.n, wt_dist', bar_width=1, xticks=x_ticks, yticks=y_ticks,
        legend=false, xlabel="Weight", ylabel="Number of Terms",
        title="Z-Weight Distribution")
    
    f = Plots.plot(f_X, f_Z, layout=(1, 2))
    display(f)
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
    [i for (i, _) in weight_distribution(S, alg, type, true)]

#############################
     # Minimum Distance
#############################

"""
    minimum_distance(Q::AbstractStabilizerCode, alg::String="trellis", sect::Bool=false)

Return the minimum distance of the stabilizer code if known, otherwise computes it.

"""
function minimum_distance(S::AbstractStabilizerCode, alg::String="auto", verbose::Bool=false)
    !ismissing(S.d) && return S.d

    # these should be different? weight? auto? BZ?
    alg ∈ ["auto", "trellis", "bruteforce"] || error("Algorithm `$alg` is not implemented in weight_enumerator.")

    if iszero(S.k)
        # "Quantum Error Correction Via Codes Over GF(4)"
        # the distance of an [𝑛,0] code is defined as the smallest non-zero weight of any stabilizer in the code
    else
        # something like this
        if alg == "auto"
            weight_enumerator(S, :Hamming, "auto", "quotient")
        elseif alg == "trellis"
            TOF_stabs = trellis_oriented_form_additive(S.stabs)
            TOF_norm = trellis_oriented_form_additive(S.dualgens)
            boundaries, num_E_sect_primal = optimal_sectionalization_Q(TOF_stabs, TOF_norm)
            verbose && println("Primal edges: $num_E_sect_primal")
            profiles_primal = trellis_profiles(TOF_stabs, TOF_norm, boundaries, "symplectic")
            boundaries, num_E_sect_dual = optimal_sectionalization_Q(TOF_norm, TOF_stabs)
            verbose && println("Dual edges: $num_E_sect_dual")
            profiles_dual = trellis_profiles(TOF_norm, TOF_stabs, boundaries, "symplectic")
            if sum(profiles_primal[2]) <= sum(profiles_dual[2])
                T_primal = sect(S, "primal", true, false)
                T_primal_HWE = weight_enumerator_Q(T_primal, "complete")
                T_dual_HWE = MacWilliams_identity(S, T_primal_HWE, true)
                poly = T_dual_HWE.polynomial - T_primal_HWE.polynomial
                S.d = minimum(filter(x->x!=0, [collect(exponent_vectors(poly))[i][1]
                    for i in 1:length(poly)]))
            else
                T_dual = sect(S, "dual", true, false)
                T_dual_HWE = weight_enumerator_Q(T_dual, :Hamming)
                T_primal_HWE = MacWilliams_identity(S, T_dual_HWE)
                poly = T_dual_HWE.polynomial - T_primal_HWE.polynomial
                S.d = minimum(filter(x->x!=0, [collect(exponent_vectors(poly))[i][1]
                    for i in 1:length(poly)]))
            end

            # T_dual = syndrome_trellis(S, "primal", true, true)
            # T_dual_HWE = weight_enumerator_Q(T_dual, "Hamming")
            # T_dual = missing
            # println("Primal trellis complete")
            # Tstabs = syndrome_trellis(S, "dual", true, true)
            # THWE = weight_enumerator_Q(Tstabs, "Hamming")
            # Tstabs = missing
            # poly = T_dual_HWE.polynomial - THWE.polynomial
            # S.d = minimum(filter(x->x!=0, [collect(exponent_vectors(poly))[i][1]
            #     for i in 1:length(poly)]))
        else
            # brute force solution here
        end
         #TODO: purity - 
    end
   
    return S.d
end

function minimum_distance_X_Z(S::AbstractStabilizerCodeCSS)
    (!ismissing(S.d_z) && !ismissing(S.d_x)) && return S.d_z, S.d_x

    # d_z = min(CX^⟂ \ CZ)
    # d_x = min(CZ^⟂ \ CX)

    # need to make these if they are missing
    if !ismissing(S.Z_orig_code)
        C1 = S.Z_orig_code
        C2 = S.X_orig_code
    else
        C1 = LinearCode(S.Z_stabs)
        C2 = LinearCode(S.X_stabs)
    end
    C1_wt_enum = weight_enumerator(C1, :Hamming)
    C2_wt_enum = weight_enumerator(C2, :Hamming)
    C1_dual_wt_enum = MacWilliams_identity(C1, C1_wt_enum)
    C2_dual_wt_enum = MacWilliams_identity(C2, C2_wt_enum)
    C1_set_diff_C2_wt_enum = C1_dual_wt_enum.polynomial - C2_dual_wt_enum.polynomial
    C2_dual_set_diff_C1_dual_wt_enum = C2_dual_wt_enum.polynomial - C1_dual_wt_enum.polynomial
    S.d_z = minimum(filter(x->x!=0, [collect(exponent_vectors(C1_set_diff_C2_wt_enum))[i][1]
        for i in 1:length(C1_set_diff_C2_wt_enum)]))
    S.d_x = minimum(filter(x->x!=0, [collect(exponent_vectors(C2_dual_set_diff_C1_dual_wt_enum))[i][1]
        for i in 1:length(C2_dual_set_diff_C1_dual_wt_enum)]))
    # the above commands will set Ci.d
    (S.d_x == C2.d && S.d_z == C1.d) ? (S.pure = true;) : (S.pure = false;)
    return S.d_z, S.d_x

    # some other paper has this as the formula
    # expsC1setdiffC2 = filter(x->x!=0, [collect(exponent_vectors(C1_set_diff_C2_wt_enum))[i][1]
    #     for i in 1:length(C1_set_diff_C2_wt_enum)])
    # expsC2dualsetdiffC2dual = filter(x->x!=0, [collect(exponent_vectors(C2_dual_set_diff_C1_dual_wt_enum))[i][1]
    #     for i in 1:length(C2_dual_set_diff_C1_dual_wt_enum)])
    # exps = vcat(expsC1setdiffC2, expsC2dualsetdiffC2dual)
    # S.d_x = minimum(exps)
    # S.d_z = maximum(exps)
end

function minimum_distance_X(S::AbstractStabilizerCodeCSS)
    ismissing(S.d_x) || return S.d_x
    
     # need to make these if they are missing
     if !ismissing(S.Z_orig_code)
        C1 = S.Z_orig_code
        C2 = S.X_orig_code
    else
        C1 = LinearCode(S.Z_stabs)
        C2 = LinearCode(S.X_stabs)
    end
    C1_wt_enum = weight_enumerator(C1, :Hamming)
    C2_wt_enum = weight_enumerator(C2, :Hamming)
    C1_dual_wt_enum = MacWilliams_identity(C1, C1_wt_enum)
    C2_dual_wt_enum = MacWilliams_identity(C2, C2_wt_enum)
    C2_dual_set_diff_C1_dual_wt_enum = C2_dual_wt_enum.polynomial - C1_dual_wt_enum.polynomial
    S.d_x = minimum(filter(x->x!=0, [collect(exponent_vectors(C2_dual_set_diff_C1_dual_wt_enum))[i][1]
        for i in 1:length(C2_dual_set_diff_C1_dual_wt_enum)]))
    return S.d_x
end

function minimum_distance_Z(S::AbstractStabilizerCodeCSS)
    ismissing(S.d_z) || return S.d_z

    # need to make these if they are missing
    if !ismissing(S.Z_orig_code)
        C1 = S.Z_orig_code
        C2 = S.X_orig_code
    else
        C1 = LinearCode(S.Z_stabs)
        C2 = LinearCode(S.X_stabs)
    end
    C1_wt_enum = weight_enumerator(C1, :Hamming)
    C2_wt_enum = weight_enumerator(C2, :Hamming)
    C1_dual_wt_enum = MacWilliams_identity(C1, C1_wt_enum)
    C2_dual_wt_enum = MacWilliams_identity(C2, C2_wt_enum)
    C1_set_diff_C2_wt_enum = C1_dual_wt_enum.polynomial - C2_dual_wt_enum.polynomial
    S.d_z = minimum(filter(x->x!=0, [collect(exponent_vectors(C1_set_diff_C2_wt_enum))[i][1]
        for i in 1:length(C1_set_diff_C2_wt_enum)]))
    return S.d_z
end

function is_pure(S::AbstractStabilizerCode)
    ismissing(S.pure) || return S.pure
    minimum_distance(S) # this needs to force the weight enumerator approach
    return S.pure
end

function is_pure(S::AbstractStabilizerCodeCSS)
    ismissing(S.pure) || return S.pure
    minimum_distance_X_Z(S)
    return S.pure
end

# TODO: pure for subsystem if no weight of gauge group is less than min dist



"""
    distrandCSS(H_X::Matrix{Int}, H_Z::Matrix{Int}, num::Int, min_dist::Int=0, debug::Int=0, field::GapObj=GAP.Globals.GF(2), max_av=Nothing)

Wrapper for the QDistRnd function DistRandCSS.
## QDistRnd documentation
- `num`: number of information sets to construct (should be large).
- `min_dist`: the algorithm stops when distance equal or below `min_dist` 
    is found and returns the result with negative sign. Set 
    `min_dist` to 0 if you want the actual distance.
- `debug`: optional integer argument containing debug bitmap (default: `0`).
    - 1 (0s  bit set): print 1st of the vectors found.
    - 2 (1st bit set): check orthogonality of matrices and of the final vector.
    - 4 (2nd bit set): show occasional progress update.
    - 8 (3rd bit set): maintain cw count and estimate the success probability.
- `field` (Options stack): Galois field, default: GF(2).
- `max_av` (Options stack): if set, terminate when `<n>` greater than `max_av`, 
    see Section Emprirical. Not set by default.
"""
# Michael Vasmer
function QDistRndCSS(H_X::Matrix{Int}, H_Z::Matrix{Int}, num::Int, min_dist::Int,
    debug::Int=0, field::GapObj=GAP.Globals.GF(2), max_av=missing)

    GAP.Packages.load("QDistRnd");
    # Convert input matrices to GAP matrices over the given field
    e = GAP.Globals.One(field)
    H_X_GAP = GapObj(H_X) * e
    H_Z_GAP = GapObj(H_Z) * e
    if ismissing(max_av)
        dist = GAP.Globals.DistRandCSS(H_X_GAP, H_Z_GAP, num, min_dist, debug, field)
    else
        dist = GAP.Globals.DistRandCSS(H_X_GAP, H_Z_GAP, num, min_dist, debug, field, max_av)
    end
    return dist
end
