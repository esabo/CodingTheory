#############################
     # Weight Enumerators
#############################
using ProgressMeter
"""
    polynomial(W::WeightEnumerator)

Returns the polynomial of the weight enumerator `W`.
"""
function polynomial(W::WeightEnumerator)
    W.type == :complete && (return _cwe_to_polynomial(W.data);)
    return _hwe_to_polynomial(W.data) #TODO
end

"""
    _cwe_to_polynomial(cwe_dict::Dict{NTuple{Q, Int}, BigInt}) where Q

Converts a Complete Weight Enumerator dictionary into an Oscar.jl multivariate polynomial.
The variables are named x0, x1, ..., x_{q-1} corresponding to the field elements.
"""
function _cwe_to_polynomial(cwe_dict::Dict{NTuple{T, Int}, BigInt}) where T
    num = length(keys(cwe_dict)[1])
    # Create a multivariate polynomial ring over the Integers (ZZ)
    var_names = ["x$(i-1)" for i in 1:num]
    R, vars = polynomial_ring(ZZ, var_names)
    
    poly = R(0)
    for (counts, coeff) in cwe_dict
        term = R(coeff)
        for i in 1:num
            term *= vars[i]^counts[i]
        end
        poly += term
    end
    
    return poly
end

"""
    _cwe_to_hwe_polynomial(cwe_dict::Dict{NTuple{T, Int}, BigInt}) where T

Reduces a Complete Weight Enumerator dictionary down to a Homogeneous Weight Enumerator 
(Hamming weight distribution) and returns it as an Oscar.jl bivariate polynomial W(x, y).
`x` represents the non-zero elements (weight) and `y` represents the zero elements.
"""
function _cwe_to_hwe_polynomial(cwe_dict::Dict{NTuple{T, Int}, BigInt}) where T
    # Create a bivariate polynomial ring over the Integers (ZZ)
    R, (x, y) = polynomial_ring(ZZ, [:x, :y])
    poly = R(0)
    # We assume the first element of the tuple counts the '0' field element,
    # which is standard for Oscar's collect(F).
    for (counts, coeff) in cwe_dict
        zero_count = counts[1]
        hamming_weight = sum(counts) - zero_count
        poly += coeff * (x^hamming_weight) * (y^zero_count)
    end

    return poly
end

"""
    type(W::WeightEnumerator)

Returns the type of the weight enumerator `W`.
"""
type(W::WeightEnumerator) = W.type

function ==(W1::WeightEnumerator, W2::WeightEnumerator)
    W1.type == W2.type || return false
    # Get the union of all unique tuple configurations present in either CWE
    all_terms = union(keys(W1.data), keys(W2.data))
    
    for term in all_terms
        # Fetch the count, defaulting to BigInt(0) if the term is absent
        # Robustly handles cases where one dictionary implicitly omits a target weight 
        # while the other explicitly stores it with a coefficient of 0.   
        count1 = get(W1.data, term, big(0))
        count2 = get(W2.data, term, big(0))
        if count1 != count2
            return false
        end
    end
    
    return true
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
    