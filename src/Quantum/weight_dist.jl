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
function _weight_enumerator_BF_Q(G::CTMatrixTypes, char_vec::Vector{zzModRingElem},
    R::Union{AbsSimpleNumFieldElem, Missing})
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

# TODO
# use new system of bounds
# functions for bare and dressed

# TODO add graph states
function minimum_distance_upper_bound!(S::AbstractSubsystemCode)
    # subsystem code
    if GaugeTrait(typeof(S)) == HasGauges()
        if CSSTrait(typeof(S)) == IsCSS()
            # bare
            _, mat = rref(vcat(S.X_stabs, reduce(vcat, [log[1] for log in S.logicals])))
            u_bound_dx_bare, _ = _min_wt_row(mat)
            _, mat = rref(vcat(S.Z_stabs, reduce(vcat, [log[2] for log in S.logicals])))
            u_bound_dz_bare, _ = _min_wt_row(mat)

            S.u_bound_dx_bare = u_bound_dx_bare
            S.u_bound_dz_bare = u_bound_dz_bare
            S.u_bound_bare = min(u_bound_dx_bare, u_bound_dz_bare)

            # dressed
            dressed_logs = dressed_logicals(S)
            _, mat = rref(vcat(S.X_stabs, reduce(vcat, [log[1] for log in dressed_logs])))
            u_bound_dx_dressed, _ = _min_wt_row(mat)
            _, mat = rref(vcat(S.Z_stabs, reduce(vcat, [log[2] for log in dressed_logs])))
            u_bound_dz_dressed, _ = _min_wt_row(mat)

            S.u_bound_dx_dressed = u_bound_dx_dressed
            S.u_bound_dz_dressed = u_bound_dz_dressed
            S.u_bound_dressed = minimum([u_bound_dx_dressed, u_bound_dz_dressed, S.u_bound_bare])
        else
            # bare
            _, mat = _rref_symp_col_swap(vcat(S.stabs, S.logs_mat))
            u_bound_bare, _ = _min_wt_row(mat)
            S.u_bound_bare = u_bound_bare

            # dressed
            _, mat = _rref_symp_col_swap(vcat(S.stabs, S.logs_mat, S.g_ops_mat))
            u_bound_dressed, _ = _min_wt_row(mat)
            S.u_bound_dressed = u_bound_dressed
        end
    # stabilizer code
    else
        # is a CSS code
        if CSSTrait(typeof(S)) == IsCSS()
            _, mat = rref(vcat(S.X_stabs, reduce(vcat, [log[1] for log in S.logicals])))
            u_bound_dx, _ = _min_wt_row(mat)

            _, mat = rref(vcat(S.Z_stabs, reduce(vcat, [log[2] for log in S.logicals])))
            u_bound_dz, _ = _min_wt_row(mat)

            S.u_bound_dx = u_bound_dx
            S.u_bound_dz = u_bound_dz
            S.u_bound = min(u_bound_dx, u_bound_dz)
        # is not a CSS code
        else
            _, mat = _rref_symp_col_swap(vcat(S.stabs, S.logs_mat))
            u_bound, _ = _min_wt_row(mat)

            S.u_bound = u_bound
        end
    end
    return nothing
end

# TODO: need to add subsystem support throughout here
"""
    minimum_distance(Q::AbstractStabilizerCode; alg::Symbol = :auto, sect::Bool=false)

Return the minimum distance of the stabilizer code if known, otherwise computes it.

"""
function minimum_distance_Gray(S::AbstractStabilizerCode; alg::Symbol = :auto, verbose::Bool = false)
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

function XZ_minimum_distance(S::AbstractStabilizerCodeCSS)
    (!ismissing(S.dz) && !ismissing(S.dx)) && return S.d_z, S.d_x

    # dz = min(CX^⟂ \ CZ)
    # dx = min(CZ^⟂ \ CX)

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
    S.dz = minimum(filter(x -> x != 0, [collect(exponent_vectors(C1_set_diff_C2_wt_enum))[i][1]
        for i in 1:length(C1_set_diff_C2_wt_enum)]))
    S.dx = minimum(filter(x -> x != 0, [collect(exponent_vectors(
        C2_dual_set_diff_C1_dual_wt_enum))[i][1] for i in eachindex(
        C2_dual_set_diff_C1_dual_wt_enum)]))
    # the above commands will set Ci.d
    (S.dx == C2.d && S.d_z == C1.d) ? (S.pure = true;) : (S.pure = false;)
    return S.dz, S.dx

    # some other paper has this as the formula
    # expsC1setdiffC2 = filter(x -> x != 0, [collect(exponent_vectors(C1_set_diff_C2_wt_enum))[i][1]
    #     for i in 1:length(C1_set_diff_C2_wt_enum)])
    # expsC2dualsetdiffC2dual = filter(x -> x != 0, [collect(exponent_vectors(C2_dual_set_diff_C1_dual_wt_enum))[i][1]
    #     for i in 1:length(C2_dual_set_diff_C1_dual_wt_enum)])
    # exps = vcat(expsC1setdiffC2, expsC2dualsetdiffC2dual)
    # S.dx = minimum(exps)
    # S.dz = maximum(exps)
end

function X_minimum_distance(S::AbstractStabilizerCodeCSS)
    ismissing(S.dx) || return S.dx
    
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
    S.dx = minimum(filter(x -> x != 0, [collect(exponent_vectors(
        C2_dual_set_diff_C1_dual_wt_enum))[i][1] for i in eachindex(
        C2_dual_set_diff_C1_dual_wt_enum)]))
    return S.dx
end

function Z_minimum_distance(S::AbstractStabilizerCodeCSS)
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
    minimum_distance_Gray(S) # this needs to force the weight enumerator approach
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
    _QDistRndCSS_GAP(H_X::Matrix{Int}, H_Z::Matrix{Int}, num::Int; min_dist::Int = 0, debug::Int = 0, field::GapObj = GAP.Globals.GF(2), max_av = Nothing)

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
function _QDistRndCSS_GAP(H_X::Matrix{Int}, H_Z::Matrix{Int}, num::Int; min_dist::Int = 0,
    debug::Int = 0, field::GapObj = GAP.Globals.GF(2), max_av = missing)

    # this requires a check on the install and load flags but since this is being moved to private
    # we will ignore it for now
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

"""
    random_information_set_minimum_distance_bound!(S::AbstractSubsystemCode, which::Symbol = :full; dressed::Bool = true, max_iters::Int = 10000, verbose::Bool = false)
    QDistRnd!(S::AbstractSubsystemCode, which::Symbol = :full; dressed::Bool = true, max_iters::Int = 10000, verbose::Bool = false)

Return an upper bound on the minimum distance of `S` and a logical of that weight using
`max_iters` random information sets. Terminates early if a logical of weight equal to the code's
currently stored lower bound is found. The current bounds on the distance are updated if the
algorithm finishes normally.

# Arguments
- `which` - either `:full`, `:X`, or `:Z` selecting which distance to bound; this parameter is
  ignored for non-CSS codes
- `dressed` - set to `true` to bound the dressed distance and `false` to bound the bare distance; this parameter is ignored for stabilizer codes
- `max_iters` - the number of random iterations
"""
function random_information_set_minimum_distance_bound!(S::T, which::Symbol = :full;
        dressed::Bool = true, max_iters::Int = 10000, verbose::Bool = false) where T <: AbstractSubsystemCode

    which ∈ (:full, :X, :Z) || throw(DomainError(which, "Must choose `:full`, `:X` or `:Z`."))
    order(field(S)) == 2 || throw(DomainError(S, "Currently only implemented for binary codes."))
    is_positive(max_iters) || throw(DomainError(max_iters, "The number of iterations must be a positive integer."))

    return random_information_set_minimum_distance_bound!(GaugeTrait(T), CSSTrait(T),
        LogicalTrait(T), S, which, dressed, max_iters, verbose)
end
QDistRnd!(S::T, which::Symbol = :full; dressed::Bool = true, max_iters::Int = 10000,
    verbose::Bool = false) where T <: AbstractSubsystemCode =
    random_information_set_minimum_distance_bound!(S, which; dressed = dressed, max_iters =
    max_iters, verbose = verbose)

function random_information_set_minimum_distance_bound!(::HasGauges, ::IsCSS, ::HasLogicals,
    S::AbstractSubsystemCode, which::Symbol, dressed::Bool, max_iters::Int, verbose::Bool)
    # CSS subsystem code

    n = S.n
    if which == :full
        if dressed && !ismissing(S.d_dressed)
            println("Dressed distance already known")
            return S.d_dressed
        elseif !dressed && !ismissing(S.d_bare)
            println("Bare distance already known")
            return S.d_bare
        end
        
        stabs = _Flint_matrix_to_Julia_T_matrix(stabilizers(S), UInt8)
        if dressed
            verbose && println("Bounding the full dressed distance")
            gauges = _Flint_matrix_to_Julia_T_matrix(gauges_matrix(S), UInt8)
            stabs = vcat(stabs, gauges)
            curr_l_bound = S.l_bound_dressed
        else
            verbose && println("Bounding the full bare distance")
            curr_l_bound = S.l_bound_bare
        end
        verbose && println("Starting lower bound: $curr_l_bound")

        _rref_no_col_swap_binary!(stabs)
        stabs = _remove_empty(stabs, :rows)
        logs = _Flint_matrix_to_Julia_T_matrix(logicals_matrix(S), UInt8)
        operators_to_reduce = vcat(stabs, logs)
        check_against = permutedims(logs[:, [n + 1:2n; 1:n]])

        # this is done in the constructor but the logical is not stored at the time
        # so must redo here
        mat = _rref_no_col_swap_binary(operators_to_reduce)
        anti = mat * check_against
        curr_u_bound, index = findmin(row_wts_symplectic(mat[findall(!iszero(anti[i, :]) for i in axes(anti, 1)), :]))
        found = operators_to_reduce[index, :]
        verbose && println("Starting upper bound: $curr_u_bound")
    else
        if which == :X
            if dressed && !ismissing(S.dx_dressed)
                println("Dressed X-distance already known")
                return S.dx_dressed
            elseif !dressed && !ismissing(S.dx_bare)
                println("Bare X-distance already known")
                return S.dx_bare
            end
        else
            if dressed && !ismissing(S.dz_dressed)
                println("Dressed Z-distance already known")
                return S.dz_dressed
            elseif !dressed && !ismissing(S.dz_bare)
                println("Bare Z-distance already known")
                return S.dz_bare
            end
        end

        stabs = _Flint_matrix_to_Julia_T_matrix(stabilizers(S)[:, (which == :X ? (1:n) : (n + 1:2n))], UInt8)
        if dressed
            gauges = _Flint_matrix_to_Julia_T_matrix(gauge_operators_matrix(S)[:, (which == :X ? (1:n) : (n + 1:2n))], UInt8)
            stabs = vcat(stabs, gauges)
        end
        _rref_no_col_swap_binary!(stabs)
        stabs = _remove_empty(stabs, :rows)
        logs = _Flint_matrix_to_Julia_T_matrix(logicals_matrix(S)[:, (which == :X ? (1:n) : (n + 1:2n))], UInt8)
        logs = _remove_empty(logs, :rows)
        operators_to_reduce = vcat(stabs, logs)
        check_against = _Flint_matrix_to_Julia_T_matrix(logicals_matrix(S)[:, (which == :X ? (n + 1:2n) : (1:n))], UInt8)
        check_against = permutedims(_remove_empty(check_against, :rows))
        curr_l_bound = if dressed
            which == :X ? S.l_bound_dx_dressed : S.l_bound_dz_dressed
        else
            which == :X ? S.l_bound_dx_bare : S.l_bound_dz_bare
        end
        verbose && println("Starting lower bound: $curr_l_bound")
        curr_u_bound, index = findmin(count(!iszero, logs[i, :]) for i in 1:size(logs, 1))
        found = logs[index, :]
        verbose && println("Starting upper bound: $curr_u_bound")
    end

    uppers, founds = _RIS_bound_loop!(operators_to_reduce, check_against, curr_l_bound,
        curr_u_bound, found, max_iters, n, verbose)
    loc = argmin(uppers)
    verbose && println("Ending $max_iters iterations with an upper bound of $(uppers[loc])")
    if dressed
        if which == :full
            S.u_bound_dressed = uppers[loc]
            flint_mat_found = matrix(field(S), permutedims(founds[loc]))
        elseif which == :X
            S.u_bound_dx_dressed = uppers[loc]
            flint_mat_found = matrix(field(S), [permutedims(founds[loc]) zeros(Int, 1, n)])
        else
            S.u_bound_dz_dressed = uppers[loc]
            flint_mat_found = matrix(field(S), [zeros(Int, 1, n) permutedims(founds[loc])])
        end
    else
        if which == :full
            S.u_bound_bare = uppers[loc]
            flint_mat_found = matrix(field(S), permutedims(founds[loc]))
        elseif which == :X
            S.u_bound_dx_bare = uppers[loc]
            flint_mat_found = matrix(field(S), [permutedims(founds[loc]) zeros(Int, 1, n)])
        else
            S.u_bound_dz_bare = uppers[loc]
            flint_mat_found = matrix(field(S), [zeros(Int, 1, n) permutedims(founds[loc])])
        end
    end
    
    return uppers[loc], flint_mat_found
end

function random_information_set_minimum_distance_bound!(::HasNoGauges, ::IsCSS, ::HasLogicals,
    S::AbstractSubsystemCode, which::Symbol, dressed::Bool, max_iters::Int, verbose::Bool)
    # CSS stabilizer code

    n = S.n
    if which == :full
        !ismissing(S.d) && (println("Distance already known"); return S.d;)

        stabs = _Flint_matrix_to_Julia_T_matrix(stabilizers(S), UInt8)
        _rref_no_col_swap_binary!(stabs)
        stabs = _remove_empty(stabs, :rows)
        logs = _Flint_matrix_to_Julia_T_matrix(logicals_matrix(S), UInt8)
        operators_to_reduce = vcat(stabs, logs)
        check_against = permutedims(logs)
        curr_l_bound = S.l_bound
        verbose && println("Starting lower bound: $curr_l_bound")
        curr_u_bound, index = findmin(count(!iszero, logs[i, :]) for i in 1:size(logs, 1))
        found = logs[index, :]
        verbose && println("Starting upper bound: $curr_u_bound")
    else
        if verbose && which == :X
            verbose && println("Bounding the X-distance")
        elseif verbose && which == :Z
            verbose && println("Bounding the Z-distance")
        end
        stabs = _Flint_matrix_to_Julia_T_matrix(stabilizers(S)[:, (which == :X ? (1:n) : (n + 1:2n))], UInt8)
        _rref_no_col_swap_binary!(stabs)
        stabs = _remove_empty(stabs, :rows)
        logs = _Flint_matrix_to_Julia_T_matrix(logicals_matrix(S)[:, (which == :X ? (1:n) : (n + 1:2n))], UInt8)
        logs = _remove_empty(logs, :rows)
        operators_to_reduce = vcat(stabs, logs)
        check_against = _Flint_matrix_to_Julia_T_matrix(logicals_matrix(S)[:, (which == :X ? (n + 1:2n) : (1:n))], UInt8)
        check_against = permutedims(_remove_empty(check_against, :rows))
        which == :X ? (curr_l_bound = S.l_bound_dx;) : (curr_l_bound = S.l_bound_dz;)
        verbose && println("Starting lower bound: $curr_l_bound")
        curr_u_bound, index = findmin(count(!iszero, logs[i, :]) for i in 1:size(logs, 1))
        found = logs[index, :]
        verbose && println("Starting upper bound: $curr_u_bound")
    end

    uppers, founds = if which == :full
        _RIS_bound_loop_symp!(operators_to_reduce, check_against, curr_l_bound,
            curr_u_bound, found, max_iters, n, verbose)
    else
        _RIS_bound_loop!(operators_to_reduce, check_against, curr_l_bound,
            curr_u_bound, found, max_iters, n, verbose)
    end
    loc = argmin(uppers)
    verbose && println("Ending $max_iters iterations with an upper bound of $(uppers[loc])")
    if which == :full
        S.u_bound = uppers[loc]
        flint_mat_found = matrix(field(S), permutedims(founds[loc]))
    elseif which == :X
        S.u_bound_dx = uppers[loc]
        S.u_bound = min(S.u_bound_dx, S.u_bound_dz)
        flint_mat_found = matrix(field(S), [permutedims(founds[loc]) zeros(Int, 1, n)])
    else
        S.u_bound_dz = uppers[loc]
        S.u_bound = min(S.u_bound_dx, S.u_bound_dz)
        flint_mat_found = matrix(field(S), [zeros(Int, 1, n) permutedims(founds[loc])])
    end

    return uppers[loc], flint_mat_found
end

function random_information_set_minimum_distance_bound!(::HasGauges, ::IsNotCSS, ::HasLogicals,
    S::AbstractSubsystemCode, which::Symbol, dressed::Bool, max_iters::Int, verbose::Bool)
    # non-CSS subsystem code

    which == :full || throw(ArguementError(which, "Parameter is not valid for non-CSS codes."))

    n = S.n
    if dressed && !ismissing(S.d_dressed)
        println("Dressed distance already known")
        return S.d_dressed
    elseif !dressed && !ismissing(S.d_bare)
        println("Bare distance already known")
        return S.d_bare
    end

    stabs = _Flint_matrix_to_Julia_T_matrix(stabilizers(S), UInt8)
    _rref_no_col_swap_binary!(stabs)
    stabs = _remove_empty(stabs, :rows)
    if dressed
        verbose && println("Bounding the full dressed distance")
        gauges = _Flint_matrix_to_Julia_T_matrix(gauge_operators_matrix(S), UInt8)
        stabs = vcat(stabs, gauges)
        curr_l_bound = S.l_bound_dressed
    else
        verbose && println("Bounding the full bare distance")
        curr_l_bound = S.l_bound_bare
    end
    verbose && println("Starting lower bound: $curr_l_bound")
    logs = _Flint_matrix_to_Julia_T_matrix(logicals_matrix(S), UInt8)
    operators_to_reduce = vcat(stabs, logs)
    check_against = permutedims(logs[:, [n + 1:2n; 1:n]])

    # this is done in the constructor but the logical is not stored at the time
    # so must redo here
    mat = _rref_no_col_swap_binary(operators_to_reduce)
    anti = mat * check_against
    curr_u_bound, index = findmin(row_wts_symplectic(mat[findall(!iszero(anti[i, :]) for i in axes(anti, 1)), :]))
    found = operators_to_reduce[index, :]
    verbose && println("Starting upper bound: $curr_u_bound")

    uppers, founds = _RIS_bound_loop_symp!(operators_to_reduce, check_against, curr_l_bound,
        curr_u_bound, found, max_iters, n, verbose)
    loc = argmin(uppers)
    if dressed
        S.u_bound_dressed = uppers[loc]
    else
        S.u_bound_bare = uppers[loc]
    end
    verbose && println("Ending $max_iters iterations with an upper bound of $(uppers[loc])")
    return uppers[loc], matrix(field(S), permutedims(founds[loc]))
end

function random_information_set_minimum_distance_bound!(::HasNoGauges, ::IsNotCSS, ::HasLogicals,
    S::AbstractSubsystemCode, which::Symbol, dressed::Bool, max_iters::Int, verbose::Bool)
    # non-CSS stabilizer code

    which == :full || throw(ArguementError(which, "Parameter is not valid for non-CSS codes."))

    n = S.n
    !ismissing(S.d) && (println("Distance already known"); return S.d;)
    verbose && println("Bounding the full distance")
    stabs = _Flint_matrix_to_Julia_T_matrix(stabilizers(S), UInt8)
    _rref_no_col_swap_binary!(stabs)
    # display(stabs)
    stabs = _remove_empty(stabs, :rows)
    logs = _Flint_matrix_to_Julia_T_matrix(logicals_matrix(S), UInt8)
    operators_to_reduce = vcat(stabs, logs)
    # println(" ")
    # display(operators_to_reduce)
    # println(" ")
    check_against = permutedims(logs[:, [n + 1:2n; 1:n]])
    curr_l_bound = S.l_bound
    verbose && println("Starting lower bound: $curr_l_bound")

    # this is done in the constructor but the logical is not stored at the time
    # so must redo here
    mat = _rref_no_col_swap_binary(operators_to_reduce)
    anti = mat * check_against
    curr_u_bound, index = findmin(row_wts_symplectic(mat[findall(!iszero(anti[i, :]) for i in axes(anti, 1)), :]))
    found = operators_to_reduce[index, :]
    verbose && println("Starting upper bound: $curr_u_bound")

    uppers, founds = _RIS_bound_loop_symp!(operators_to_reduce, check_against, curr_l_bound,
        curr_u_bound, found, max_iters, n, verbose)
    loc = argmin(uppers)
    S.u_bound = uppers[loc]
    verbose && println("Ending $max_iters iterations with an upper bound of $(uppers[loc])")
    return uppers[loc], matrix(field(S), permutedims(founds[loc]))
end

function random_information_set_minimum_distance_bound(::HasGauges, ::IsNotCSS, ::HasNoLogicals, S::AbstractSubsystemCode, which::Symbol, dressed::Bool, max_iters::Int, verbose::Bool)
    # non-CSS subsystem graph state

    which == :full || throw(ArguementError(which, "Parameter is not valid for non-CSS codes."))

    n = S.n
    stabs = _Flint_matrix_to_Julia_T_matrix(stabilizers(S), UInt8)
    _rref_no_col_swap_binary!(stabs)
    stabs = _remove_empty(stabs, :rows)
    if dressed
        !ismissing(S.d_dressed) && (println("Dressed distance already known"); return S.d_dressed;)
        verbose && println("Bounding the full dressed distance")
        gauges = _Flint_matrix_to_Julia_T_matrix(gauge_operators_matrix(S), UInt8)
        stabs = vcat(stabs, gauges)
    else
        verbose && println("Bounding the full bare distance")
    end
    curr_l_bound = if dressed
        S.l_bound_dressed
    else
        S.l_bound_bare
    end
    verbose && println("Starting lower bound: $curr_l_bound")

    # this is done in the constructor but the logical is not stored at the time
    # so must redo here
    curr_u_bound, index = findmin(row_wts_symplectic(stabs))
    found = stabs[index, :]
    verbose && println("Starting upper bound: $curr_u_bound")

    uppers, founds = _RIS_bound_loop_symp!(stabs, curr_l_bound, curr_u_bound, found, max_iters, n,
        verbose)
    loc = argmin(uppers)
    if dressed
        S.u_bound_dressed = uppers[loc]
    else
        S.u_bound_bare = uppers[loc]
    end
    verbose && println("Ending $max_iters iterations with an upper bound of $(uppers[loc])")
    return uppers[loc], matrix(field(S), permutedims(founds[loc]))
end

function random_information_set_minimum_distance_bound(::HasGauges, ::IsCSS, ::HasNoLogicals, S::AbstractSubsystemCode, which::Symbol, dressed::Bool, max_iters::Int, verbose::Bool)
    # CSS subsystem graph state

    n = S.n
    if which == :full
        !ismissing(S.d_dressed) && (println("Dressed distance already known"); return S.d_dressed;)
        verbose && println("Bounding the full dressed distance")
        stabs = _Flint_matrix_to_Julia_T_matrix(stabilizers(S), UInt8)
        if dressed
            gauges = _Flint_matrix_to_Julia_T_matrix(gauges_matrix(S), UInt8)
            stabs = vcat(stabs, gauges)
        end
        _rref_no_col_swap_binary!(stabs)
        stabs = _remove_empty(stabs, :rows)
        curr_l_bound = if dressed
            S.l_bound_dressed
        else
            S.l_bound_bare
        end
        verbose && println("Starting lower bound: $curr_l_bound")

        # this is done in the constructor but the logical is not stored at the time
        # so must redo here
        curr_u_bound, index = findmin(row_wts_symplectic(stabs))
        found = stabs[index, :]
        verbose && println("Starting upper bound: $curr_u_bound")
    else
        stabs = _Flint_matrix_to_Julia_T_matrix(stabilizers(S)[:, (which == :X ? (1:n) : (n + 1:2n))], UInt8)
        if dressed
            gauges = _Flint_matrix_to_Julia_T_matrix(gauge_operators_matrix(S)[:, (which == :X ? (1:n) : (n + 1:2n))], UInt8)
            stabs = vcat(stabs, gauges)
        end
        _rref_no_col_swap_binary!(stabs)
        stabs = _remove_empty(stabs, :rows)
        curr_l_bound = if dressed
            which == :X ? S.l_bound_dx_dressed : S.l_bound_dz_dressed
        else
            which == :X ? S.l_bound_dx_bare : S.l_bound_dz_bare
        end
        verbose && println("Starting lower bound: $curr_l_bound")
        curr_u_bound, index = _min_wt_row(stabs)
        found = stabs[index, :]
        verbose && println("Starting upper bound: $curr_u_bound")
    end

    uppers, founds = _RIS_bound_loop!(stabs, curr_l_bound, curr_u_bound, found, max_iters, n,
        verbose)
    loc = argmin(uppers)
    verbose && println("Ending $max_iters iterations with an upper bound of $(uppers[loc])")
    if dressed
        if which == :full
            S.u_bound_dressed = uppers[loc]
            flint_mat_found = matrix(field(S), permutedims(founds[loc]))
        elseif which == :X
            S.u_bound_dx_dressed = uppers[loc]
            flint_mat_found = matrix(field(S), [permutedims(founds[loc]) zeros(Int, 1, n)])
        else
            S.u_bound_dz_dressed = uppers[loc]
            flint_mat_found = matrix(field(S), [zeros(Int, 1, n) permutedims(founds[loc])])
        end
    else
        if which == :full
            S.u_bound_bare = uppers[loc]
            flint_mat_found = matrix(field(S), permutedims(founds[loc]))
        elseif which == :X
            S.u_bound_dx_bare = uppers[loc]
            flint_mat_found = matrix(field(S), [permutedims(founds[loc]) zeros(Int, 1, n)])
        else
            S.u_bound_dz_bare = uppers[loc]
            flint_mat_found = matrix(field(S), [zeros(Int, 1, n) permutedims(founds[loc])])
        end
    end
    
    return uppers[loc], flint_mat_found
end

function random_information_set_minimum_distance_bound(::HasNoGauges, ::IsNotCSS, ::HasNoLogicals, S::AbstractSubsystemCode, which::Symbol, dressed::Bool, max_iters::Int, verbose::Bool)
    # non-CSS stabilizer graph state

    which == :full || throw(ArguementError(which, "Parameter is not valid for non-CSS codes."))

    n = S.n
    !ismissing(S.d) && (println("Distance already known"); return S.d;)
    verbose && println("Bounding the full distance")
    curr_l_bound = S.l_bound
    verbose && println("Starting lower bound: $curr_l_bound")

    # this is done in the constructor but the logical is not stored at the time
    # so must redo here
    stabs = _Flint_matrix_to_Julia_T_matrix(stabilizers(S), UInt8)
    _rref_no_col_swap_binary!(stabs)
    stabs = _remove_empty(stabs, :rows)
    curr_u_bound, index = findmin(row_wts_symplectic(stabs))
    found = stabs[index, :]
    verbose && println("Starting upper bound: $curr_u_bound")

    # TODO write
    uppers, founds = _RIS_bound_loop_symp!(stabs, curr_l_bound, curr_u_bound, found, max_iters, n, 
        verbose)
    loc = argmin(uppers)
    S.u_bound = uppers[loc]
    verbose && println("Ending $max_iters iterations with an upper bound of $(uppers[loc])")
    return uppers[loc], matrix(field(S), permutedims(founds[loc]))
end

function random_information_set_minimum_distance_bound(::HasNoGauges, ::IsCSS, ::HasNoLogicals, S::AbstractSubsystemCode, which::Symbol, dressed::Bool, max_iters::Int, verbose::Bool)
    # CSS stabilizer graph state

    n = S.n
    if which == :full
        !ismissing(S.d) && (println("Distance already known"); return S.d;)
        verbose && println("Bounding the full distance")
        curr_l_bound = S.l_bound
        verbose && println("Starting lower bound: $curr_l_bound")

        # this is done in the constructor but the logical is not stored at the time
        # so must redo here
        stabs = _Flint_matrix_to_Julia_T_matrix(stabilizers(S), UInt8)
        _rref_no_col_swap_binary!(stabs)
        stabs = _remove_empty(stabs, :rows)
        curr_u_bound, index = findmin(row_wts_symplectic(stabs))
        found = stabs[index, :]
        verbose && println("Starting upper bound: $curr_u_bound")
    else
        if verbose && which == :X
            verbose && println("Bounding the X-distance")
        elseif verbose && which == :Z
            verbose && println("Bounding the Z-distance")
        end
        stabs = _Flint_matrix_to_Julia_T_matrix(stabilizers(S)[:, (which == :X ? (1:n) : (n + 1:2n))], UInt8)
        _rref_no_col_swap_binary!(stabs)
        stabs = _remove_empty(stabs, :rows)
        which == :X ? (curr_l_bound = S.l_bound_dx;) : (curr_l_bound = S.l_bound_dz;)
        verbose && println("Starting lower bound: $curr_l_bound")
        curr_u_bound, index = _min_wt_row(stabs)
        found = stabs[index, :]
        verbose && println("Starting upper bound: $curr_u_bound")
    end

    uppers, founds = if which == :full
        _RIS_bound_loop_symp!(stabs, curr_l_bound, curr_u_bound, found, max_iters, n, verbose)
    else
        _RIS_bound_loop!(stabs, curr_l_bound, curr_u_bound, found, max_iters, n, verbose)
    end
    loc = argmin(uppers)
    verbose && println("Ending $max_iters iterations with an upper bound of $(uppers[loc])")
    if which == :full
        S.u_bound = uppers[loc]
        flint_mat_found = matrix(field(S), permutedims(founds[loc]))
    elseif which == :X
        S.u_bound_dx = uppers[loc]
        S.u_bound = min(S.u_bound_dx, S.u_bound_dz)
        flint_mat_found = matrix(field(S), [permutedims(founds[loc]) zeros(Int, 1, n)])
    else
        S.u_bound_dz = uppers[loc]
        S.u_bound = min(S.u_bound_dx, S.u_bound_dz)
        flint_mat_found = matrix(field(S), [zeros(Int, 1, n) permutedims(founds[loc])])
    end

    return uppers[loc], flint_mat_found
end

function _RIS_bound_loop_symp!(operators_to_reduce::Matrix{T}, check_against::Matrix{T},
    curr_l_bound::Int, curr_u_bound::Int, found::Vector{T}, max_iters::Int, n::Int,
    verbose::Bool) where T <: Integer

    num_thrds = Threads.nthreads()
    verbose && println("Detected $num_thrds threads.")

    flag = Threads.Atomic{Bool}(true)
    uppers = [curr_u_bound for _ in 1:num_thrds]
    founds = [found for _ in 1:num_thrds]
    thread_load = Int(floor(max_iters / num_thrds))
    remaining = max_iters - thread_load * num_thrds
    verbose && (prog_meter = Progress(max_iters);)
    # Threads.@threads for t in 1:num_thrds
    for t in 1:num_thrds
        orig_ops = deepcopy(operators_to_reduce)
        log_test = zeros(Int, size(orig_ops, 1), size(check_against, 2))
        perm_ops = similar(orig_ops)
        ops = similar(orig_ops)
        perm = collect(1:n)
        for _ in 1:(thread_load + (t <= remaining ? 1 : 0))
            if flag[]
                shuffle!(perm)
                # println("original")
                # display(orig_ops)
                # display(perm)
                _col_permutation_symp!(perm_ops, orig_ops, perm)
                # println("permuted")
                # display(perm_ops)
                # return
                _rref_no_col_swap_binary!(perm_ops)
                _col_permutation_symp!(ops, perm_ops, invperm(perm))
                LinearAlgebra.mul!(log_test, ops, check_against)
                for i in axes(log_test, 1)
                    # then ops[i, :] is a logical
                    if any(isodd, log_test[i, :])
                        w = 0
                        @inbounds for j in 1:n
                            (isodd(ops[i, j]) || isodd(ops[i, j + n])) && (w += 1;)
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
            end
            verbose && next!(prog_meter)
        end
    end
    verbose && finish!(prog_meter)

    return uppers, founds
end

function _RIS_bound_loop!(operators_to_reduce::Matrix{T}, check_against::Matrix{T},
    curr_l_bound::Int, curr_u_bound::Int, found::Vector{T}, max_iters::Int, n::Int,
    verbose::Bool) where T <: Integer

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
        log_test = zeros(Int, size(orig_ops, 1), size(check_against, 2))
        perm_ops = similar(orig_ops)
        ops = similar(orig_ops)
        perm = collect(1:n)
        for _ in 1:(thread_load + (t <= remaining ? 1 : 0))
            if flag[]
                shuffle!(perm)
                _col_permutation!(perm_ops, orig_ops, perm)
                # modifying this in place is not thread safe (apparently)
                _rref_no_col_swap_binary!(perm_ops)
                _col_permutation!(ops, perm_ops, invperm(perm))
                LinearAlgebra.mul!(log_test, ops, check_against)
                for i in axes(log_test, 1)
                    # then ops[i, :] is a logical
                    if any(isodd, log_test[i, :])
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
            end
            verbose && next!(prog_meter)
        end
    end
    verbose && finish!(prog_meter)

    return uppers, founds
end
