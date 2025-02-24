# Copyright (c) 2022 Eric Sabo, Benjamin Ide, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
         # General
#############################
# using AllocCheck
# using Profile
using ProgressMeter
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
    R, vars = polynomial_ring(Nemo.ZZ, ord_E)

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

    R, (x, y) = polynomial_ring(base_ring(CWE.polynomial), [:x, :y])
    poly = R(0)
    for i in 1:length(CWE.polynomial)
        exps = exponent_vector(CWE.polynomial, i)
        poly += coeff(CWE.polynomial, i) * x^sum(exps[2:end]) * y^exps[1]
    end
    return WeightEnumerator(poly, :Hamming)
end

# TODO: fix doc string
"""
    Sterns_attack(C::AbstractLinearCode, w::Int, p::Int, l::Int; num_find::Int =  2, max_iters::Int = 50000)

Search for codewords of `C` of weight `w` using Stern's attack and return any found.
"""
function Sterns_attack(C::AbstractLinearCode, w::Int, p::Int, l::Int; num_find::Int =  2, max_iters::Int = 50000)
    # requires 2 * x = 0
    Int(order(C.F)) == 2 || throw(ArgumentError("Only valid for binary codes."))
    if !ismissing(C.d)
        C.d <= w <= C.n || throw(ArgumentError("Target weight must be between the minimum distance and code length."))
    else
        1 <= w <= C.n || throw(ArgumentError("Target weight must be positive and no more than the code length."))
    end
    1 <= p <= ceil(C.k / 2) || throw(ArgumentError("p must be between 1 and k/2"))

    # G = generator_matrix(C, true)
    H = parity_check_matrix(C, true)
    nr, nc = size(H)
    # nr, nc = size(G)
    1 <= l <= nr || throw(ArgumentError("l must be between 1 and n - k"))
    H2 = deepcopy(H)
    # G2 = deepcopy(G)
    # len_w_vecs = Vector{fqPolyRepMatrix}()
    len_w_vecs = Set{fqPolyRepMatrix}()
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
        AπAs = Vector{Tuple{Vector{Int}, Vector{fqPolyRepFieldElem}}}()
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
        BπBs = Vector{Tuple{Vector{Int}, Vector{fqPolyRepFieldElem}}}()
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

# TODO: this does not produce the optimal set of matrices
# see section 7.3 of White's thesis for comments on this
function information_sets(G::CTMatrixTypes, alg::Symbol = :Edmonds; permute::Bool = false, only_A::Bool = false)

    alg ∈ (:Brouwer, :Zimmermann, :White, :Chen, :Bouyuklieva, :Edmonds) || throw(ArgumentError("Unknown information set algorithm. Expected `:Brouwer`, `:Zimmermann`, `:White`, `:Chen`, `:Bouyuklieva`, or `:Edmonds`."))
    # TODO should rref to begin with and remove empty rows?
    nr, nc = size(G)
    gen_mats = Vector{}()
    perms = Vector{}()
    rnks = Vector{Int}()
    
    rnk = nr
    #TODO Brouw and Zimm should use the same code for rref then discard the last matrix
    if alg == :Brouwer
        for i in 0:Int(floor(nc / nr)) - 1
            start_ind = i * nr + 1
            rnk, Gi, Pi = _rref_col_swap(G, 1:nr, start_ind:nc)
            if rnk < nr # for Brouwer the Gi must all have full rank
                break
            end

            if only_A
                Ai = Gi[:, setdiff(1:nc, start_ind:(i + 1) * nr)]
                push!(gen_mats, Ai)
                push!(perms, Pi)
                push!(rnks, rnk)
            else
                if permute
                    # permute identities to the front
                    pivots = collect(start_ind:(i + 1) * nr)
                    σ = [pivots; setdiff(1:nc, pivots)]
                    Gi = Gi[:, σ]
                    Pi = Pi[:, σ]
                end
                push!(gen_mats, Gi)
                push!(perms, Pi)
                push!(rnks, rnk)
            end
        end
    elseif alg == :Zimmermann
        for i in 0:Int(floor(nc / nr))
            rnk, Gi, Pi = _rref_col_swap(G, 1:nr, i * rnk + 1:nc)
            if only_A
                Ai = Gi[:, setdiff(1:nc, i * nr + 1:i * nr + rnk)]
                push!(gen_mats, Ai)
                push!(perms, Pi)
                push!(rnks, rnk)
            else
                if permute
                    # permute identities to the front
                    pivots = collect(i * nr + 1:(i + 1) * nr)
                    σ = [pivots; setdiff(1:nc, pivots)]
                    Gi = Gi[:, σ]
                    Pi = Pi[:, σ]
                end
                push!(gen_mats, Gi)
                push!(perms, Pi)
                push!(rnks, rnk)
            end
        end
    elseif alg == :White
        # TODO: this is not true when the parity-check matrix is true
        # the expansion factor of the code
        for i in 0:div(nc, nr) - 1
            # could use Gi here instead of G
            rnk, Gi, Pi = _rref_col_swap(G, 1:nr, i * nr + 1:(i + 1) * nr)
            # display(Gi)
            # println(rnk)
            push!(gen_mats, Gi)
            push!(perms, Pi)
            push!(rnks, rnk)
        end
    elseif alg == :Chen
        Gi, _, Pi, rnk = _standard_form(G)
        if only_A
            Ai = Gi[:, rnk + 1:nc]
            push!(gen_mats, Ai)
            push!(perms, Pi)
            push!(rnks, rnk)
        else
            push!(gen_mats, Gi)
            push!(perms, Pi)
            push!(rnks, rnk)
        end
    end
    return gen_mats, perms, rnks
end

function information_set_lower_bound(r::Int, n::Int, k::Int, l::Int, rank_defs::Vector{Int},
    info_set_alg::Symbol; even::Bool = false, doubly_even::Bool = false, triply_even::Bool = false)

    info_set_alg ∈ (:auto, :Brouwer, :Zimmermann, :White, :Chen, :Bouyuklieva, :Edmonds) || throw(ArgumentError("Unknown information set algorithm. Expected `:auto`, `:Brouwer`, `:Zimmermann`, `:White`, `:Chen`, `:Bouyuklieva`, or `:Edmonds`."))

    lower = 0
    if info_set_alg == :Brouwer
        lower = r * length(rank_defs)
    elseif info_set_alg == :Zimmermann
        h = length(rank_defs)
        lower = count(x -> x != 0, rank_defs) 
        for i in 1:h
            lower += maximum([0, r - rank_defs[i]]) 
        end
    elseif info_set_alg == :Chen
        lower = Int(ceil(n * r / k))
    elseif info_set_alg == :White
        lower = 0
        for i in 1:l
	        lower += Int(ceil(n * maximum([0, r - rank_defs[i]]) / (l * (k + rank_defs[i]))))
        end
    # elseif info_set_alg == :Bouyuklieva
    #     continue
    # elseif info_set_alg == :Edmonds
    #     continue
    end
    if lower == 0
        println("initial lower bound is set to 0")
    end

    (!triply_even && !doubly_even && even) && (lower += lower % 2;)
    (!triply_even && doubly_even) && (lower += 4 - lower % 4;)
    triply_even && (lower += 8 - lower % 8;)
    return lower
end

"""
    minimum_distance_Gray(C::AbstractLinearCode; alg::Symbol = :Zimmermann, v::Bool = false)
Return the minimum distance of `C` using a deterministic algorithm based on enumerating
constant weight codewords of the binary reflected Gray code. If a word of minimum weight
is found before the lower and upper bounds cross, it is returned; otherwise, the zero vector 
is returned.

show_progress will display a progress meter for each iteration of a weight that takes longer than 
    a second
"""
function minimum_distance_Gray(C::AbstractLinearCode; alg::Symbol = :auto, v::Bool = false, 
    show_progress = true)

    ord_F = Int(order(C.F))
    ord_F == 2 || throw(ArgumentError("Currently only implemented for binary codes."))
    #TODO add :Bouyuklieva, and :Edmonds

    # We implement the following algorithms described in White's thesis:
    # :Brouwer      Algo 2.2
    # :Zimmermann   Algo 2.4 
    # :Chen         Algo 2.6 
    # :White        Algo 3.1 
    alg ∈ (:auto, :Brouwer, :Zimmermann, :White, :Chen) || throw(ArgumentError("Unknown information set algorithm. Expected `:auto`, `:Brouwer`, `:Zimmermann`, `:White`, `:Chen`"))

    if alg == :auto
        if typeof(C) <: AbstractCyclicCode
            v && println("Detected a cyclic code, using Chen's adaption.")
            alg = :Chen
            # TODO: fix this case
        elseif typeof(C) <: AbstractQuasiCyclicCode
            v && println("Detected a quasi-cyclic code, using White's adaption.")
            alg = :White
        else
            v && println("Using Zimmermann's algorithm.")
            alg = :Zimmermann
        end
    end
    alg == :auto && throw(ErrorException("Could not determine minimum distance algo automatically"))

    if alg in (:Brouwer, :Zimmermann) 
        return _minimum_distance_BZ(C::AbstractLinearCode; info_set_alg = alg, verbose = v, show_progress = show_progress)
    end
    println("Warning: old enumeration algorithm selected. Performance will be slow") # TODO remove when all code updated
    return _minimum_distance_enumeration_with_matrix_multiply(C::AbstractLinearCode; info_set_alg = alg)
end

function _minimum_distance_BZ(C::AbstractLinearCode; info_set_alg::Symbol = :Zimmermann,
    verbose::Bool = false, dbg = Dict(), show_progress=false)
    dbg_key_exit_r = "exit_r"

    ord_F = Int(order(C.F))
    ord_F == 2 || throw(ArgumentError("Currently only implemented for binary codes."))
    C.k < 2^16 || throw(DomainError("The given linear code has length k >= 2^16 which is not supported"))
    info_set_alg ∈ (:Brouwer, :Zimmermann) || throw(ArgumentError("Unknown information set algorithm. Expected `:Brouwer`, `:Zimmermann`"))

    generator_matrix(C, true) # ensure G_stand exists
    if _has_empty_vec(C.G, :cols) 
        #TODO err string can instruct the user to construct a new code without 0 cols and tell them the function for that
        throw(ArgumentError("Codes with standard form of generator matrix having 0 columns not supported")) 
    end
    # generate if not pre-stored
    parity_check_matrix(C)

    A_mats, perms_mats, rnks = information_sets(C.G, info_set_alg, permute = true, only_A = false)

    A_mats = [deepcopy(_Flint_matrix_to_Julia_T_matrix(Ai, UInt16)') for Ai in A_mats]
    # A_mats_trunc = () 
    perms_mats = [deepcopy(_Flint_matrix_to_Julia_T_matrix(Pi, UInt16)') for Pi in perms_mats]
    h = length(A_mats)
    # println("Starting loop to refine upper bound. Initial upper bound ", C.u_bound, " num of mats is ", length(A_mats), " dimension ", size(A_mats[1]))
    rank_defs = zeros(Int, h)

    if length(keys(dbg)) > 0
        println("Debug mode ON")
    end

    if haskey(dbg, dbg_key_exit_r)
        verbose && println("dbg Dict: largest message weight searched stored @key=$dbg_key_exit_r")
        dbg[dbg_key_exit_r] = -1
    end

    k, n = size(C.G)
    A_mats_trunc = [Matrix{UInt16}(undef, k, n-k) for _ in 1:length(A_mats)]
    for i in 1:size(A_mats, 1) 
       A_mats_trunc[i] = deepcopy(A_mats[i][k+1 : n, :])
    end

    if info_set_alg == :Brouwer && rnks[h] != k
        println("Rank of last matrix too small")
        return
    end
    if verbose
        print("Generated $h information sets with ranks: ")
        for i in 1:h
            i == h ? (println(rnks[i]);) : (print("$(rnks[i]), ");)
            # will only be using the rank deficits here
            # at the moment, the information sets are always disjoint so the relative
            # rank is zero
            # TODO huh? check this comment and setup properly
            rank_defs[i] = C.k - rnks[i]
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

    # initial_perm_ind will match the permutation we use for the 'found' vector if the found vector is nonzero. To simplify the code below we're going to choose an initial permutation arbitrarily.  
    initial_perm_ind = 1 
    # following loop is the r=1 case of the enumeration. We do this case here because we want to make a good guess at the terminating r before we start multiple threads
    for (j, g) in enumerate(A_mats) # loop over the A_mats rather than the original G because it would add another case to deal with later 
        # can make this faster with dots and views
        w, i = _min_wt_col(g)
        if w <= C.u_bound
            found = g[:, i] 
            C.u_bound = w
            y = perms_mats[j] * found 
        end
    end

    verbose && println("Current upper bound: $(C.u_bound)")
    found = A_mats[1][:, 1]

    l = 0
    if verbose 
        _, _, b_rnks = information_sets(C.G, :Brouwer, permute = true, only_A = false)
        b_h = length(b_rnks)
        b_lower_bounds = [information_set_lower_bound(r+1, n, k, l, [k - 0 for i in 1:b_h], :Brouwer, even = even_flag, doubly_even = doubly_even_flag, triply_even = triply_even_flag) for r in 1:k-1]
        b_r_term = findfirst(x -> x ≥ C.u_bound, b_lower_bounds)

        # _, _, z_rnks = information_sets(G, :Zimmermann, permute = true, only_A = false)
        # z_h = length(b_rnks)
        # z_lower_bounds = [information_set_lower_bound(r+1, n, k, l, [k - z_rnks[i] for i in 1:z_h], :Zimmermann, even = even_flag, doubly_even = doubly_even_flag, triply_even = triply_even_flag) for r in 1:k-1]
        # z_r_term = findfirst(x -> x ≥ C.u_bound, z_lower_bounds)
        # verbose && println("ranks: Brouwer $b_rnks Zimm $z_rnks")
        # verbose && println("Predicted termination weight based on current upper bound: Brouwer $b_r_term Zimm $z_r_term")
    end

    #Note the r+1 here. 
    lower_bounds_for_prediction = [information_set_lower_bound(r+1, n, k, l, rank_defs, info_set_alg, even = even_flag, doubly_even = doubly_even_flag, triply_even = triply_even_flag) for r in 1:k-1]
    r_term = findfirst(x -> x ≥ C.u_bound, lower_bounds_for_prediction)
    if isnothing(r_term)
        raise(DomainError("invalid termination r")) 
    end
    verbose && println("Predicted termination weight based on current upper bound: $r_term")

    #In the main loop we check if lower bound > upper bound before we enumerate and so the lower bounds for the loop use r not r+1
    lower_bounds = [information_set_lower_bound(r, n, k, l, rank_defs, info_set_alg, even = even_flag, doubly_even = doubly_even_flag, triply_even = triply_even_flag) for r in 1:k-1]

    predicted_work_factor = fld(n, k) * sum([binomial(k, i) for i in 1:r_term])
    verbose && println("Predicted work factor: $predicted_work_factor")
    if show_progress 
        prog_bar = Progress(predicted_work_factor, dt=1.0, showspeed=true) # updates no faster than once every 1s
    end
    weight_sum_bound = min(2 * C.u_bound + 5, n-k)
    verbose && println("Codeword weights initially checked on first $weight_sum_bound entries")

    num_thrds = Threads.nthreads()
    verbose && println("Number of threads ", num_thrds)
    for r in 2:k
        if r > 2^16
            verbose && println("Warning: Reached an r larger than 2^16") 
        end
        C.l_bound < lower_bounds[r] && (C.l_bound = lower_bounds[r];)
        # an even code can't have have an odd minimum weight
        # (!triply_even_flag && !doubly_even_flag && even_flag) && (C.l_bound += C.l_bound % 2;)
        # (!triply_even_flag && doubly_even_flag) && (C.l_bound += 4 - C.l_bound % 4;)
        # triply_even_flag && (C.l_bound += 8 - C.l_bound % 8;)
        if C.l_bound >= C.u_bound
            dbg[dbg_key_exit_r] = r-1
            break
        end
        verbose && println("r: $r")
        verbose && println("Lower bound: $(C.l_bound)")
        verbose && println("Upper bound: $(C.u_bound)")

        if verbose
            i_count = 0
            for i in 1:h
                r - rank_defs[i] ≤ 0 && (i_count += 1;)
            end
            i_count > 0 && println("$i_count of the original $h information sets no longer contribute to the lower bound")
        end
        p = Int(characteristic(C.F))

        uppers = [C.u_bound for _ in 1:num_thrds]
        founds = [found for _ in 1:num_thrds]
        exit_thread_indicator_vec = [initial_perm_ind for _ in 1:num_thrds]
        keep_going = Threads.Atomic{Bool}(true)

        bin = extended_binomial(C.k, r)

        thrd_stop_msg = "Stopping current thread, main loop finished"

        Threads.@threads for ind in 1:num_thrds 
            len = (ind == num_thrds) ? bin - (num_thrds - 1) * fld(bin, num_thrds) : fld(bin, num_thrds)

            # iteration begins with a single matrix multiplication of the generator matrix by first_vec
            init_rank = 1 + (ind - 1) * fld(bin, num_thrds)
            first_vec = zeros(Int, k)
            if init_rank == 1
                for i in 1:r
                    first_vec[i] = 1
                end
            else
                CodingTheory._subset_unrank_to_vec!(init_rank, UInt64(r), first_vec)
            end
 
            # as in White Algo 7.1 we loop over matrices first 
            for i in 1:h
                if keep_going[] == false
                    verbose && println(thrd_stop_msg)
                    break
                end

                c_itr = zeros(UInt16, C.n - C.k) 
                is_first = true
                curr_mat = A_mats_trunc[i]
                count = UInt128(0)

                for u in SubsetGrayCode(k, r, len, init_rank)
                    if keep_going[] == false
                        println(thrd_stop_msg)
                        break
                    end
                    show_progress && ProgressMeter.next!(prog_bar) 
                    if r - rank_defs[i] > 0
                        if is_first 
                            LinearAlgebra.mul!(c_itr, curr_mat, first_vec)
                            @inbounds @simd for j in eachindex(c_itr) 
                                c_itr[j] %= p
                            end
                            is_first = false
                        else
                            for ci in u 
                                if ci != -1
                                    @simd for i in eachindex(c_itr)
                                        @inbounds c_itr[i] = xor(c_itr[i], curr_mat[i, ci])
                                    end
                                end
                            end
                        end

                        partial_weight = r + sum(view(c_itr, 1:weight_sum_bound))

                        if uppers[ind] > partial_weight
                            w = r + sum(c_itr) 
                            verbose && @assert w != 0
                            if uppers[ind] > w 
                                subset_vec_full = zeros(Int, k)
                                CodingTheory._subset_unrank_to_vec!(UInt128(init_rank + count), UInt64(r), subset_vec_full)

                                uppers[ind] = w 
                                founds[ind] = vcat(subset_vec_full, c_itr)
                                verbose && @assert size(founds[ind], 1) == C.n "found vector has length $(size(founds[ind], 1)) but should be n=$(C.n)"
                                exit_thread_indicator_vec[ind] = i 

                                println("Adjusting (local) upper bound: $w for c_itr=$(Int.(c_itr))")
                                if C.l_bound == uppers[ind]
                                    println("early exit")
                                    Threads.atomic_cas!(keep_going, true, false)
                                else
                                    r_term = findfirst(x -> x ≥ C.u_bound, lower_bounds)
                                    isnothing(r_term) && (r_term = k;)
                                    verbose && println("Updated termination weight: $r_term")
                                end
                            end
                        end
                    end
                end
            end
            loc = argmin(uppers) 
            C.u_bound = uppers[loc]
            found = founds[loc]
            initial_perm_ind = exit_thread_indicator_vec[loc]
        end
    end

    C.d = C.u_bound
    y = matrix(C.F, 1, n, perms_mats[initial_perm_ind] * found) # weight(y) >= C.d, with equality not being the typical case
    verbose && @assert iszero(C.H * transpose(y))
    if dbg[dbg_key_exit_r] == -1
        dbg[dbg_key_exit_r] = r
    end 
    show_progress && ProgressMeter.finish!(prog_bar)
    verbose && println("Computation complete")
    return C.u_bound, y
end

function _minimum_distance_enumeration_with_matrix_multiply(C::AbstractLinearCode; info_set_alg::Symbol = :auto,
    verbose::Bool = false, dbg = Dict())

    ord_F = Int(order(C.F))
    ord_F == 2 || throw(ArgumentError("Currently only implemented for binary codes."))
    info_set_alg ∈ (:auto, :Brouwer, :Zimmermann, :White, :Chen, :Bouyuklieva, :Edmonds) || throw(ArgumentError("Unknown information set algorithm. Expected `:auto`, `:Brouwer`, `:Zimmermann`, `:White`, `:Chen`, `:Bouyuklieva`, or `:Edmonds`."))

    generator_matrix(C, true) # ensure G_stand exists
    if _has_empty_vec(C.G, :cols) 
        #TODO err string can instruct the user to construct a new code without 0 cols and tell them the function for that
        throw(ArgumentError("Codes with standard form of generator matrix having 0 columns not supported")) 
    end
    G = C.G 
    # generate if not pre-stored
    parity_check_matrix(C)

    if info_set_alg == :auto
        if typeof(C) <: AbstractCyclicCode
            verbose && println("Detected a cyclic code, using Chen's adaption.")
            info_set_alg = :Chen
            # TODO: fix this case
        elseif typeof(C) <: AbstractQuasiCyclicCode
            verbose && println("Detected a quasi-cyclic code, using White's adaption.")
            info_set_alg = :White
        else
            info_set_alg = :Edmonds
        end
    end
    # should have no need to permute to find better ranks because of Edmond's?
    A_mats, perms_mats, rnks = information_sets(G, info_set_alg, permute = true, only_A = false)
    A_mats = [deepcopy(_Flint_matrix_to_Julia_int_matrix(Ai)') for Ai in A_mats]
    perms_mats = [deepcopy(_Flint_matrix_to_Julia_int_matrix(Pi)') for Pi in perms_mats]
    h = length(A_mats)
    rank_defs = zeros(Int, h)

    k, n = size(G)
    if verbose
        print("Generated $h information sets with ranks: ")
        for i in 1:h
            i == h ? (println(rnks[i]);) : (print("$(rnks[i]), ");)
            # will only be using the rank deficits here
            # at the moment, the information sets are always disjoint so the relative
            # rank is zero
            # TODO huh? check this comment and setup properly
            rank_defs[i] = C.k - rnks[i]
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

    # initial_perm_ind will match the permutation we use for the 'found' vector if the found vector is nonzero. To simplify the code below we're going to choose an initial permutation arbitrarily.  
    initial_perm_ind = 1 
    found = A_mats[1][:, 1] 
    for (j, g) in enumerate(A_mats) # loop over the A_mats rather than the original G because it would add another case to deal with later 
        # can make this faster with dots and views
        w, i = _min_wt_col(g)
        if w <= C.u_bound
            found = g[:, i] 
            C.u_bound = w
            initial_perm_ind = j
            y = perms_mats[initial_perm_ind] * found 
        end
    end

    verbose && println("Current upper bound: $(C.u_bound)")
    verbose && !iszero(found) && println("Found element matching upper bound.")

    info_set_alg == :Chen ? (l = C.l;) : (l = 0;)
    lower_bounds = [information_set_lower_bound(r, n, k, l, rank_defs, info_set_alg, even = even_flag, doubly_even = doubly_even_flag, triply_even = triply_even_flag) for r in 1:k]
    r_term = findfirst(x -> x ≥ C.u_bound, lower_bounds)
    isnothing(r_term) && (r_term = k;)
    verbose && println("Predicted termination weight based on current upper bound: $r_term")

    num_thrds = Threads.nthreads()
    verbose && println("Detected $num_thrds threads.")
    power = Int(floor(log(2, num_thrds)))

    for r in 1:k
        C.l_bound < lower_bounds[r] && (C.l_bound = lower_bounds[r];)
        # an even code can't have have an odd minimum weight
        # (!triply_even_flag && !doubly_even_flag && even_flag) && (C.l_bound += C.l_bound % 2;)
        # (!triply_even_flag && doubly_even_flag) && (C.l_bound += 4 - C.l_bound % 4;)
        # triply_even_flag && (C.l_bound += 8 - C.l_bound % 8;)
        verbose && println("r: $r")
        verbose && println("Lower bound: $(C.l_bound)")
        verbose && println("Upper bound: $(C.u_bound)")
        if C.l_bound >= C.u_bound
            C.d = C.u_bound
            y = perms_mats[initial_perm_ind] * found 
            return C.u_bound, (C.F).(y) 
        end

        if verbose
            count = 0
            for i in 1:h
                r - rank_defs[i] ≤ 0 && (count += 1;)
            end
            count > 0 && println("$count of the original $h information sets no longer contribute to the lower bound")
        end

        flag = Threads.Atomic{Bool}(true)
        verbose && println("Detected $num_thrds threads.")

        p = Int(characteristic(C.F))
        uppers = [C.u_bound for _ in 1:num_thrds]
        founds = [found for _ in 1:num_thrds]
        exit_thread_indicator_vec = [initial_perm_ind for _ in 1:num_thrds]

        c = zeros(Int, C.n)
        itrs = [GrayCode(C.k, r, digits(m - 1, base = 2, pad = power), mutate = true) for m in 1:num_thrds]
        Threads.@threads for m in 1:num_thrds 
            itr = itrs[m]
            for u in itr 
                if flag[]
                    for i in 1:h 
                        vec = u
                        if r - rank_defs[i] > 0 
                            LinearAlgebra.mul!(c, A_mats[i], vec)

                            w = r
                            @inbounds for j in 1:n 
                                c[j] % p != 0 && (w += 1;)
                            end

                            if uppers[m] > w
                                uppers[m] = w
                                founds[m] .= c 
                                exit_thread_indicator_vec[m] = i 
                                verbose && println("Adjusting (local) upper bound: $w")
                                if C.l_bound == uppers[m]
                                    Threads.atomic_cas!(flag, true, false)
                                    break
                                else
                                    r_term = findfirst(x -> x ≥ C.u_bound, lower_bounds)
                                    isnothing(r_term) && (r_term = k;)
                                    verbose && println("Updated termination weight: $r_term")
                                end
                            end
                        end
                    end
                else
                    break
                end
            end
        end
        loc = argmin(uppers)
        C.u_bound = uppers[loc]
        found = founds[loc]
        initial_perm_ind = exit_thread_indicator_vec[loc]
    end

    # at this point we are guaranteed to have found the answer
    C.d = C.u_bound
    y = perms_mats[initial_perm_ind] * found 
    # iszero(C.H * transpose(y))
    return C.u_bound, (C.F).(y)
end

"""
    words_of_weight(C::AbstractLinearCode, l_bound::Int, u_bound::Int; verbose::Bool = false)

Return all the codewords of `C` of Hamming weight in the range `[l_bound, u_bound]`.
"""
function words_of_weight(C::AbstractLinearCode, l_bound::Int, u_bound::Int; verbose::Bool = false)
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
        gen_mats = information_sets(G, :Chen)
    elseif typeof(C) <: AbstractQuasiCyclicCode
        verbose && println("Detected a quasi-cyclic code, using White's adaption.")
        gen_mats = information_sets(G, :White)
    else
        gen_mats = information_sets(G, :Zimmermann)
    end
    gen_mats_Julia = [_Flint_matrix_to_Julia_int_matrix(x[2])' for x in gen_mats]
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
            lower = _lower_bounds(r, C.n, C.k, 0, [0], :Chen)
        elseif typeof(C) <: AbstractQuasiCyclicCode
            lower = _lower_bounds(r, C.n, C.k, C.l, rank_defs, :White)
        else
            lower = _lower_bounds(r, C.n, C.k, 0, rank_defs, :BZ)
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
    words_of_weight(C::AbstractLinearCode, bound::Int; verbose::Bool = false)

Return all the codewords of `C` of Hamming weight `bound`.
"""
words_of_weight(C::AbstractLinearCode, bound::Int; verbose::Bool = false) = words_of_weight(C, bound, bound, verbose = verbose)

# untested
# TODO: figure out if even weight upper weight needs to subtract or not
function _words_of_weight_high(C::AbstractLinearCode, l_bound::Int, u_bound::Int;
    verbose::Bool = false)

    p = Int(characteristic(C.F))
    G = generator_matrix(C)
    if typeof(C) <: AbstractCyclicCode
        verbose && println("Detected a cyclic code, using Chen's adaption.")
        gen_mats = information_sets(G, :Chen)
    elseif typeof(C) <: AbstractQuasiCyclicCode
        verbose && println("Detected a quasi-cyclic code, using White's adaption.")
        gen_mats = information_sets(G, :White)
    else
        gen_mats = information_sets(G, :Zimmermann)
    end
    gen_mats_Julia = [_Flint_matrix_to_Julia_int_matrix(x[2])' for x in gen_mats]
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
            upper = _lower_bounds(r, C.n, C.k, 0, [0], :Chen)
        elseif typeof(C) <: AbstractQuasiCyclicCode
            upper = _lower_bounds(r, C.n, C.k, C.l, rank_defs, :White)
        else
            upper = _lower_bounds(r, C.n, C.k, 0, rank_defs, :BZ)
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

# TODO: change above to BigInts
# TODO: type stability of this
"""
    partial_weight_distribution(C::AbstractLinearCode, bound::Int; compact::Bool = false)

Return the partial weight distribution of `C` up to weight `bound`. If `compact` is false,
the result will be a `Vector{BigInt}` of length `length(C) + 1` whose `i`th entry is the
number of codewords of `C` of Hamming weight `i - 1`. Otherwise, the result is a
`Vector{Tuple{Int, BigInt}}` whose entries specify the nonzero indices and values of the
above.
"""
function partial_weight_distribution(C::AbstractLinearCode, bound::Int; compact::Bool = false)
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
- This algorithm simultaneously computes the minimum distance and stores the words of
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
        gen_mats = information_sets(G, :Chen)
    elseif typeof(C) <: AbstractQuasiCyclicCode
        verbose && println("Detected a quasi-cyclic code, using White's adaption.")
        gen_mats = information_sets(G, :White)
    else
        gen_mats = information_sets(G, :Zimmermann)
    end
    gen_mats_Julia = [_Flint_matrix_to_Julia_int_matrix(x[2])' for x in gen_mats] # deepcopy doesn't save anythign here
    # gen_mats_Julia = [_Flint_matrix_to_Julia_int_matrix(x[2]) for x in gen_mats]
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
            lower = _lower_bounds(r, C.n, C.k, 0, [0], :Chen)
        elseif typeof(C) <: AbstractQuasiCyclicCode
            lower = _lower_bounds(r, C.n, C.k, C.l, rank_defs, :White)
        else
            lower = _lower_bounds(r, C.n, C.k, 0, rank_defs, :BZ)
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
            prefix = digits(m - 1, base = 2, pad = power)
            for u in GrayCode(C.k, r, prefix, mutate = true)
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

"""
    random_information_set_minimum_distance_bound!(C::AbstractLinearCode; max_iters::Int = 10000, verbose::Bool = false)

Return an upper bound on the minimum distance of `C` and a codeword of that weight using
`max_iters` random information sets.
"""
function random_information_set_minimum_distance_bound!(C::AbstractLinearCode; max_iters::Int = 10000, verbose::Bool = false)

    order(field(C)) == 2 || throw(DomainError(C, "Currently only implemented for binary codes."))
    is_positive(max_iters) || throw(DomainError(max_iters, "The number of iterations must be a positive integer."))

    !ismissing(C.d) && (println("Distance already known"); return C.d;)
    verbose && println("Bounding the distance")
    G = _Flint_matrix_to_Julia_T_matrix(_rref_no_col_swap(C.G), UInt8)
    G = _remove_empty(G, :rows)
    upper, ind  = _min_wt_row(G)
    found = G[ind, :]
    # TODO this can contradict C.u_bound cause that requires a different found
    verbose && println("Starting lower bound: $(C.l_bound)")
    verbose && println("Starting upper bound: $upper")

    uppers, founds = _RIS_bound_loop(G, C.l_bound, upper, found, max_iters, C.n, verbose)
    loc = argmin(uppers)
    C.u_bound = uppers[loc]
    verbose && println("Ending $max_iters iterations with an upper bound of $(uppers[loc])")
    return uppers[loc], matrix(field(C), permutedims(founds[loc]))
end

function _RIS_bound_loop(operators_to_reduce::Matrix{T}, curr_l_bound::Int, curr_u_bound::Int,
    found::Vector{T}, max_iters::Int, n::Int, verbose::Bool) where T <: Integer

    num_thrds = Threads.nthreads()
    verbose && println("Detected $num_thrds threads.")

    flag = Threads.Atomic{Bool}(true)
    uppers = [curr_u_bound for _ in 1:num_thrds]
    founds = [found for _ in 1:num_thrds]
    thread_load = Int(floor(max_iters / num_thrds))
    remaining = max_iters - thread_load * num_thrds
    verbose && (prog_meter = Progress(max_iters);)
    Threads.@threads for t in 1:num_thrds
        orig_ops = deepcopy(operators_to_reduce)
        perm_ops = similar(orig_ops)
        ops = similar(orig_ops)
        perm = collect(1:n)
        for _ in 1:(thread_load + (t <= remaining ? 1 : 0))
            if flag[]
                shuffle!(perm)
                _col_permutation!(perm_ops, orig_ops, perm)
                _rref_no_col_swap_binary!(perm_ops)
                _col_permutation!(ops, perm_ops, invperm(perm))

                for i in axes(perm_ops, 1)
                    w = 0
                    @inbounds for j in 1:n
                        isodd(ops[i, j]) && (w += 1;)
                    end

                    if uppers[t] > w
                        uppers[t] = w
                        founds[t] .= ops[i, :]
                        verbose && println("Adjusting (thread's local) upper bound: $w")
                        if curr_l_bound == w
                            verbose && println("Found a logical that matched the lower bound of $curr_l_bound")
                            Threads.atomic_cas!(flag, true, false)
                            break
                        end
                    end
                end
            end
            verbose && next!(prog_meter)
        end
    end
    verbose && finish!(prog_meter)

    return uppers, founds
end

#############################
    # Weight Enumerators
#############################
# TODO: doc string?
function weight_enumerator_classical(T::Trellis; type::Symbol = :complete)
    type ∈ (:complete, :Hamming) ||
        throw(ArgumentError("Unsupported weight enumerator type '$type'. Expected ':complete' or ':Hamming'."))

    if type == :complete && !ismissing(T.CWE)
        return T.CWE
    elseif type == :Hamming && !ismissing(T.CWE)
        return CWE_to_HWE(T.CWE)
    end

    # if this ever changes or permutes will have to store with T
    elms = collect(field(T.code))
    lookup = Dict(value => key for (key, value) in enumerate(elms))
    R, vars = polynomial_ring(Nemo.ZZ, length(elms))

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

# TODO: remove C from this, store in WE struct
"""
    MacWilliams_identity(C::AbstractLinearCode, W::WeightEnumerator; dual::Symbol = :Euclidean)

Return the weight enumerator of the dual (`:Euclidean` or `:Hermitian`) of `C` obtained
by applying the MacWilliams identities to `W`.
"""
function MacWilliams_identity(C::AbstractLinearCode, W::WeightEnumerator; dual::Symbol = :Euclidean)
    dual ∈ (:Euclidean, :Hermitian) ||
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
        K, ζ = cyclotomic_field(3, :ζ)
        R, vars = polynomial_ring(K, 3)
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
        if is_prime(q)
            K, ω = cyclotomic_field(Int(characteristic(C.F)), :ω)
            R, vars = polynomial_ring(K, q)
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
            K, ω = cyclotomic_field(Int(characteristic(C.F)), :ω)
            R, vars = polynomial_ring(K, q)
            prime_field = GF(Int(characteristic(C.F)))
            _, λ = primitive_basis(C.F, prime_field)
            elms = collect(C.F)
            func_args = []
            for i in 1:q
                inner_sum = R(0)
                for j in 1:q
                    β = elms[i] * elms[j]
                    β_exp = _expand_element(β, prime_field, λ, false)
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
    weight_enumerator(C::AbstractLinearCode; type::Symbol = :complete, alg::Symbol = :auto)

Return either the `:complete` or `:Hamming` weight enumerator of `C` using the algorithm `alg`.
"""
function weight_enumerator(C::AbstractLinearCode; type::Symbol = :complete, alg::Symbol = :auto)
    type ∈ (:complete, :Hamming) ||
        throw(ArgumentError("Unsupported weight enumerator type '$type'. Expected ':complete' or ':Hamming'."))
    alg ∈ (:auto, :trellis, :bruteforce) ||
        throw(ArgumentError("Algorithm `$alg` is not implemented in weight_enumerator."))

    if type == :complete && !ismissing(C.weight_enum)
        return C.weight_enum
    elseif type == :Hamming && !ismissing(C.weight_enum)
        return CWE_to_HWE(C.weight_enum)
    end

    if alg == :auto
        if cardinality(C) <= 1e6 # random cutoff
            C.weight_enum = _weight_enumerator_BF(C.G)
            HWE = CWE_to_HWE(C.weight_enum)
            C.d = minimum(filter(x -> x != 0, [collect(exponent_vectors(HWE.polynomial))[i][1]
                for i in 1:length(HWE.polynomial)]))
            type == :Hamming && return HWE
            return C.weight_enum
        elseif rate(C) > 0.5
            D = dual(C)
            if cardinality(D) <= 1e6 # random cutoff
                D.weight_enum = _weight_enumerator_BF(D.G)
            else
                weight_enumerator_classical(syndrome_trellis(D, "primal", false), type = type)
            end
            C.weight_enum = MacWilliams_identity(D, D.weight_enum)
            HWE = CWE_to_HWE(C.weight_enum)
            C.d = minimum(filter(x -> x != 0, [collect(exponent_vectors(HWE.polynomial))[i][1]
                for i in 1:length(HWE.polynomial)]))
            type == :Hamming && return HWE
            return C.weight_enum
        else
            return weight_enumerator_classical(syndrome_trellis(C, "primal", false), type = type)
        end
    elseif alg == :trellis
        return weight_enumerator_classical(syndrome_trellis(C, "primal", false), type = type)
    elseif alg == :bruteforce
        C.weight_enum = _weight_enumerator_BF(C.G)
        HWE = CWE_to_HWE(C.weight_enum)
        C.d = minimum(filter(x -> x != 0, [collect(exponent_vectors(HWE.polynomial))[i][1]
            for i in 1:length(HWE.polynomial)]))
        type == :Hamming && return HWE
        return C.weight_enum
    end
end

"""
    weight_distribution(C::AbstractLinearCode; alg::Symbol = :auto, compact::Bool = true)

Return the weight distribution of `C` using the algorithm `alg`. If `compact` is false,
the result will be a `Vector{BigInt}` of length `length(C) + 1` whose `i`th entry is the
number of codewords of `C` of Hamming weight `i - 1`. Otherwise, the result is a
`Vector{Tuple{Int, BigInt}}` whose entries specify the nonzero indices and values of the
above.
"""
function weight_distribution(C::AbstractLinearCode; alg::Symbol = :auto, compact::Bool = true)
    alg ∈ (:auto, :trellis, :bruteforce) ||
        throw(ArgumentError("Algorithm `$alg` is not implemented in weight_enumerator."))

    ismissing(C.weight_enum) && weight_enumerator(C, type = :complete, alg = alg)
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
    weight_plot(C::AbstractLinearCode; alg::Symbol = :auto)

Return a bar graph of the weight distribution of `C`.

# Note
- Run `using Makie` to activate this extension.
"""
function weight_plot end

"""
    support(C::AbstractLinearCode)

Returns the support of `C`.

# Notes
- The support of `C` is the collection of nonzero exponents of the Hamming weight enumerator of `C`.
"""
support(C::AbstractLinearCode) = [i for (i, _) in weight_distribution(C, alg = :auto,
    compact = true)]

#############################
     # Minimum Distance
#############################

# TODO: update doc strings for these and this whole file in general
"""
    minimum_distance(C::AbstractLinearCode; alg::Symbol = :trellis, sect::Bool = false, verbose::Bool = false)

Return the minimum distance of the linear code if known, otherwise computes it
using the algorithm of `alg`. If `alg = "trellis"`, the sectionalization flag
`sect` can be set to true to further compactify the reprsentation.
"""
function minimum_distance(C::AbstractLinearCode; alg::Symbol = :trellis, sect::Bool = false,
    verbose::Bool = false)

    !ismissing(C.d) && return C.d

    alg ∈ (:auto, :Gray, :trellis, :Leon, :bruteforce, :wt_dist) ||
        throw(ArgumentError("Unexpected algorithm '$alg'."))
    
    if alg == :auto
        D = dual(C)
        if cardinality(C) <= 1e6 # random cutoff
            C.weight_enum = _weight_enumerator_BF(C.G)
            HWE = CWE_to_HWE(C.weight_enum)
            C.d = minimum(filter(x -> x != 0, [collect(exponent_vectors(HWE.polynomial))[i][1]
                for i in 1:length(HWE.polynomial)]))
            return C.d
        elseif rate(C) > 0.5
            D = dual(C)
            if cardinality(D) <= 1e6 # random cutoff
                D.weight_enum = _weight_enumerator_BF(D.G)
            else
                weight_enumerator_classical(syndrome_trellis(D, "primal", false), type = type)
            end
            C.weight_enum = MacWilliams_identity(D, D.weight_enum)
            HWE = CWE_to_HWE(C.weight_enum)
            C.d = minimum(filter(x -> x != 0, [collect(exponent_vectors(HWE.polynomial))[i][1]
                for i in 1:length(HWE.polynomial)]))
            return C.d
        else
            return minimum_distance_Gray(C, verbose = verbose)
        end
    elseif alg == :Gray
        return minimum_distance_Gray(C, verbose = verbose)
    elseif alg == :trellis
        weight_enumerator_classical(syndrome_trellis(C, "primal", false), type = type)
        return C.d
    elseif alg == :bruteforce
        C.weight_enum = _weight_enumerator_BF(C.G)
        HWE = CWE_to_HWE(C.weight_enum)
        C.d = minimum(filter(x -> x != 0, [collect(exponent_vectors(HWE.polynomial))[i][1]
            for i in 1:length(HWE.polynomial)]))
        return C.d
    elseif alg == :wt_dist
        HWE = weight_enumerator(C, type = :Hamming, alg = alg)
        !ismissing(C.d) && return C.d
        # this line should only be needed to be run if the weight enumerator is known
        # but the minimum distance is intentionally set to missing
        # ordering here can be a bit weird
        C.d = minimum(filter(x -> x != 0, [collect(exponent_vectors(HWE.polynomial))[i][1]
            for i in 1:length(HWE.polynomial)]))
        return C.d
    # elseif alg == "Leon"
    #     Leon(C)
    end
end
