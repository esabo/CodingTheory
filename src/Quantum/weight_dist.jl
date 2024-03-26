# Copyright (c) 2022 - 2024 Eric Sabo, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

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
function MacWilliams_identity(S::AbstractStabilizerCode, W::WeightEnumerator; dual::Bool = false)
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

function weight_enumerator(S::AbstractStabilizerCode; type::Symbol = :complete,
    alg::Symbol = :auto, set::Symbol = :all)

    type ∈ (:complete, :Hamming) || throw(ArgumentError("Unsupported weight enumerator type '$type'. Expected ':complete' or ':Hamming'."))
    alg ∈ (:auto, :trellis, :bruteforce) || throw(ArgumentError("Algorithm `$alg` is not implemented in weight_enumerator."))
    set ∈ (:all, :stabilizers, :logicals, :quotient) || throw(ArgumentError("Unsupported set type '$set'. Expected ':all', ':stabilizers', ':logicals', ':quotient'."))

    if set ∈ (:all, :logicals) && ismissing(S.sgn_CWE_logs)
        logs_mat = logicals_matrix(S)
        S.sgn_CWE_logs = _weight_enumerator_BF_Q(logs_mat, S.char_vec, missing)
    end

    if set != :logicals && ismissing(S.sgn_CWE_stabs)
        if alg == :bruteforce || cardinality(S) <= 1e6
            S.sgn_CWE_stabs = _weight_enumerator_BF_Q(S.stabs, S.char_vec,
                parent(S.sgn_CWE_logs.polynomial))
        else
            # trellis solution here
        end
    end

    if set ∈ (:all, :quotient) && ismissing(S.sgn_CWE_dual)
        if alg == :bruteforce || BigInt(characteristic(S.F))^(S.n + S.k) <= 3e6
            S.sgn_CWE_dual = _weight_enumerator_BF_Q(vcat(S.stabs, logicals_matrix(S)), S.char_vec,
                parent(S.sgn_CWE_logs.polynomial))
        else
            # trellis solution here
        end
    end
    
    if !ismissing(S.sgn_CWE_stabs) && !ismissing(S.sgn_CWE_dual)
        # compute minimum distance here
        poly = WeightEnumerator(S.sgn_CWE_dual.polynomial - S.sgn_CWE_stabs.polynomial, :complete)
        HWE = CWE_to_HWE(poly)
        S.d = minimum(filter(x -> x != 0, [collect(exponent_vectors(HWE.polynomial))[i][1]
            for i in 1:length(HWE.polynomial)]))
    end

    if type == :complete
        set == :all && return S.sgn_CWE_stabs, S.sgn_CWE_dual, S.sgn_CWE_logs, poly
        set == :stabilizers && return S.sgn_CWE_stabs
        set == :logicals && return S.sgn_CWE_logs
        return poly
    else
        set == :all && return CWE_to_HWE(S.sgn_CWE_stabs), CWE_to_HWE(S.sgn_CWE_dual), CWE_to_HWE(S.sgn_CWE_logs), CWE_to_HWE(poly)
        set == :stabilizers && return CWE_to_HWE(S.sgn_CWE_stabs)
        set == :logicals && return CWE_to_HWE(S.sgn_CWE_logs)
        return HWE
    end
end

# MAGMA returns this format
# [ <0, 1>, <4, 105>, <6, 280>, <8, 435>, <10, 168>, <12, 35> ]
function weight_distribution(S::AbstractStabilizerCode; alg::Symbol = :auto, compact::Bool = true, set::Symbol = :all)

    alg ∈ (:auto, :trellis, :bruteforce) || throw(ArgumentError("Algorithm `$alg` is not implemented in weight_enumerator."))
    set ∈ (:all, :stabilizers, :logicals, :quotient) || throw(ArgumentError("Unsupported set type '$set'. Expected ':all', ':stabilizers', ':logicals', ':quotient'."))

    wt_enums = weight_enumerator(S, type = :Hamming, alg = alg, set = set)

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

function weight_enumerator_quantum(T::Trellis; type::Symbol = :complete)
    type ∈ (:complete, :Hamming) || throw(ArgumentError("Unsupported weight enumerator type '$type'. Expected ':complete' or ':Hamming'."))

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
    weight_plot(S::AbstractStabilizerCode; alg::Symbol = :auto, type::Symbol = :stabilizer)

Return a bar graph of the weight distribution related to `S`.

If `type` is `:stabilizer`, the weight distribution of the stabilizers are computed.
If `type` is `:normalizer`, the weight distrbution of the normalizer of the stabilizers
are computed. If `type` is `:quotient`, the weight distrbution of the normalizer mod the
stabilizers is computed.

# Note
- Run `using Makie` to activate this extension.
"""
function weight_plot end

"""
    weight_plot_CSS_X(S::AbstractStabilizerCodeCSS; alg::Symbol = :auto)

Return a bar graph of the weight distribution of the `X` stabilizers.

# Note
- Run `using Makie` to activate this extension.
"""
function weight_plot_CSS_X end

"""
    weight_plot_CSS_Z(S::AbstractStabilizerCodeCSS; alg::Symbol = :auto)

Return a bar graph of the weight distribution of the `Z` stabilizers.

# Note
- Run `using Makie` to activate this extension.
"""
function weight_plot_CSS_Z end

"""
    weight_plot_CSS(S::AbstractStabilizerCodeCSS; alg::Symbol = :auto)

Return bar plots of the weight distribution of the both the
`X` and 'Z' stabilizers, separately.

# Note
- Run `using Makie` to activate this extension.
"""
function weight_plot_CSS end

# TODO: standardize set and type throughout file
"""
    support(S::AbstractStabilizerCode; alg::Symbol = :auto, type::Symbol = :stabilizer)

Returns the support related to `S`.

The support is the collection of nonzero exponents of the Hamming
weight enumerator. If `type` is `stabilizer`, the support of the stabilizers are computed.
If `type` is `normalizer`, the support of the normalizer of the stabilizers
are computed. If `type` is `quotient`, the support of the normalizer mod the
stabilizers (logical representatives only) is computed.
"""
support(S::AbstractStabilizerCode; alg::Symbol = :auto, type::Symbol = :stabilizer) =
    [i for (i, _) in weight_distribution(S, alg = alg, set = type, compact = true)]

#############################
     # Minimum Distance
#############################

"""
    minimum_distance(Q::AbstractStabilizerCode; alg::Symbol = :auto, sect::Bool=false)

Return the minimum distance of the stabilizer code if known, otherwise computes it.

"""
function minimum_distance(S::AbstractStabilizerCode; alg::Symbol = :auto, verbose::Bool = false)
    !ismissing(S.d) && return S.d

    # these should be different? weight? auto? BZ?
    alg ∈ (:auto, :trellis, :bruteforce) || throw(ArgumentError("Algorithm `$alg` is not implemented in weight_enumerator."))

    if iszero(S.k)
        # "Quantum Error Correction Via Codes Over GF(4)"
        # the distance of an [𝑛,0] code is defined as the smallest non-zero weight of any stabilizer in the code
    else
        # something like this
        if alg == :auto
            weight_enumerator(S, type = :Hamming, alg = auto, set = :quotient)
        elseif alg == :trellis
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
                T_primal_HWE = weight_enumerator_quantum(T_primal, type = :complete)
                T_dual_HWE = MacWilliams_identity(S, T_primal_HWE, dual = true)
                poly = T_dual_HWE.polynomial - T_primal_HWE.polynomial
                S.d = minimum(filter(x -> x != 0, [collect(exponent_vectors(poly))[i][1]
                    for i in 1:length(poly)]))
            else
                T_dual = sect(S, "dual", true, false)
                T_dual_HWE = weight_enumerator_quantum(T_dual, type = :Hamming)
                T_primal_HWE = MacWilliams_identity(S, T_dual_HWE)
                poly = T_dual_HWE.polynomial - T_primal_HWE.polynomial
                S.d = minimum(filter(x -> x != 0, [collect(exponent_vectors(poly))[i][1]
                    for i in 1:length(poly)]))
            end

            # T_dual = syndrome_trellis(S, "primal", true, true)
            # T_dual_HWE = weight_enumerator_quantum(T_dual, "Hamming")
            # T_dual = missing
            # println("Primal trellis complete")
            # Tstabs = syndrome_trellis(S, "dual", true, true)
            # THWE = weight_enumerator_quantum(Tstabs, "Hamming")
            # Tstabs = missing
            # poly = T_dual_HWE.polynomial - THWE.polynomial
            # S.d = minimum(filter(x -> x != 0, [collect(exponent_vectors(poly))[i][1]
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
    C1_wt_enum = weight_enumerator(C1, type = :Hamming)
    C2_wt_enum = weight_enumerator(C2, type = :Hamming)
    C1_dual_wt_enum = MacWilliams_identity(C1, C1_wt_enum)
    C2_dual_wt_enum = MacWilliams_identity(C2, C2_wt_enum)
    C1_set_diff_C2_wt_enum = C1_dual_wt_enum.polynomial - C2_dual_wt_enum.polynomial
    C2_dual_set_diff_C1_dual_wt_enum = C2_dual_wt_enum.polynomial - C1_dual_wt_enum.polynomial
    S.d_z = minimum(filter(x -> x != 0, [collect(exponent_vectors(C1_set_diff_C2_wt_enum))[i][1]
        for i in 1:length(C1_set_diff_C2_wt_enum)]))
    S.d_x = minimum(filter(x -> x != 0, [collect(exponent_vectors(
        C2_dual_set_diff_C1_dual_wt_enum))[i][1] for i in eachindex(
        C2_dual_set_diff_C1_dual_wt_enum)]))
    # the above commands will set Ci.d
    (S.d_x == C2.d && S.d_z == C1.d) ? (S.pure = true;) : (S.pure = false;)
    return S.d_z, S.d_x

    # some other paper has this as the formula
    # expsC1setdiffC2 = filter(x -> x != 0, [collect(exponent_vectors(C1_set_diff_C2_wt_enum))[i][1]
    #     for i in 1:length(C1_set_diff_C2_wt_enum)])
    # expsC2dualsetdiffC2dual = filter(x -> x != 0, [collect(exponent_vectors(C2_dual_set_diff_C1_dual_wt_enum))[i][1]
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
    C1_wt_enum = weight_enumerator(C1, type = :Hamming)
    C2_wt_enum = weight_enumerator(C2, type = :Hamming)
    C1_dual_wt_enum = MacWilliams_identity(C1, C1_wt_enum)
    C2_dual_wt_enum = MacWilliams_identity(C2, C2_wt_enum)
    C2_dual_set_diff_C1_dual_wt_enum = C2_dual_wt_enum.polynomial - C1_dual_wt_enum.polynomial
    S.d_x = minimum(filter(x -> x != 0, [collect(exponent_vectors(
        C2_dual_set_diff_C1_dual_wt_enum))[i][1] for i in eachindex(
        C2_dual_set_diff_C1_dual_wt_enum)]))
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
    C1_wt_enum = weight_enumerator(C1, type = :Hamming)
    C2_wt_enum = weight_enumerator(C2, type = :Hamming)
    C1_dual_wt_enum = MacWilliams_identity(C1, C1_wt_enum)
    C2_dual_wt_enum = MacWilliams_identity(C2, C2_wt_enum)
    C1_set_diff_C2_wt_enum = C1_dual_wt_enum.polynomial - C2_dual_wt_enum.polynomial
    S.d_z = minimum(filter(x -> x != 0, [collect(exponent_vectors(C1_set_diff_C2_wt_enum))[i][1]
        for i in eachindex(C1_set_diff_C2_wt_enum)]))
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

# TODO: max_av is what type?
"""
    distrandCSS(H_X::Matrix{Int}, H_Z::Matrix{Int}, num::Int; min_dist::Int = 0, debug::Int = 0, field::GapObj = GAP.Globals.GF(2), max_av = Nothing)

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
function QDistRndCSS(H_X::Matrix{Int}, H_Z::Matrix{Int}, num::Int; min_dist::Int = 0,
    debug::Int = 0, field::GapObj = GAP.Globals.GF(2), max_av = missing)

    Packages.load("QDistRnd");
    # Convert input matrices to GAP matrices over the given field
    e = Globals.One(field)
    H_X_GAP = GapObj(H_X) * e
    H_Z_GAP = GapObj(H_Z) * e
    if ismissing(max_av)
        dist = Globals.DistRandCSS(H_X_GAP, H_Z_GAP, num, min_dist, debug, field)
    else
        dist = Globals.DistRandCSS(H_X_GAP, H_Z_GAP, num, min_dist, debug, field, max_av)
    end
    return dist
end
