# Copyright (c) 2021, 2022, 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    LinearCode(G::CTMatrixTypes, parity::Bool=false, brute_force_WE::Bool=true)

Return the linear code constructed with generator matrix `G`. If the optional paramater `parity` is
set to `true`, a linear code is built with `G` as the parity-check matrix. If the optional parameter
`brute_force_WE` is `true`, the weight enumerator and (and therefore the distance) is calculated when
there are fewer than 1.5e5 codewords.
"""
function LinearCode(G::CTMatrixTypes, parity::Bool = false, brute_force_WE::Bool = true)
    iszero(G) && return parity ? IdentityCode(base_ring(G), ncols(G)) : ZeroCode(base_ring(G), ncols(G))

    G_new = deepcopy(G)
    G_new = _remove_empty(G_new, :rows)

    C = if parity
            H = kernel(G_new, side = :right)
            rnk_H = rank(H)
            # println(rnk_H)
            # display(H)
        if ncols(H) == rnk_H
            H_tr = transpose(H)
        else
            # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
            nr = nrows(H)
            H_tr = zero_matrix(base_ring(H), rnk_H, nr)
            for r in 1:nr
                for c in 1:rnk_H
                    !iszero(H[r, c]) && (H_tr[c, r] = H[r, c];)
                end
            end
        end
        H = H_tr
        H_stand, G_stand, P, k = _standard_form(H)
        k == ncols(H) && return IdentityCode(base_ring(G_new), ncols(H))
        ub1, _ = _min_wt_row(H)
        ub2, _ = _min_wt_row(H_stand)
        ub = min(ub1, ub2)
        # treat G as the parity-check matrix H
        LinearCode(base_ring(G_new), ncols(H), nrows(H_stand), missing, 1, ub, H, G_new, H_stand, G_stand, P, missing)
    else
        G_stand, H_stand, P, k = _standard_form(G_new)
        if k == ncols(G)
            C = IdentityCode(base_ring(G_new), ncols(G))
            C.G = G_new
            return C
        end
        H = ismissing(P) ? H_stand : H_stand * P
        ub1, _ = _min_wt_row(G_new)
        ub2, _ = _min_wt_row(G_stand)
        ub = min(ub1, ub2)
        LinearCode(base_ring(G_new), ncols(G_new), k, missing, 1, ub, G_new, H, G_stand, H_stand, P, missing)
    end

    if k == 0
        set_minimum_distance!(C, C.n)
    else
        if brute_force_WE && BigInt(order(base_ring(G)))^min(k, ncols(G) - k) <= 1.5e5
            # BUG only adding this new case because MacWilliams has a new Oscar problem
            if BigInt(order(base_ring(G)))^k ≤ 1.5e5
                C.weight_enum = _weight_enumerator_BF(C.G_stand)
            else
                C.weight_enum = if 2k <= ncols(G)
                    _weight_enumerator_BF(C.G_stand)
                else
                    MacWilliams_identity(dual(C), _weight_enumerator_BF(C.H_stand))
                end
            end
            d = minimum(filter(is_positive, first.(exponent_vectors(CWE_to_HWE(C.weight_enum).polynomial))))
            set_minimum_distance!(C, d)
        end
    end

    return C
end

# TODO: add doc strings
function LinearCode(G::T, H::T, brute_force_WE::Bool = true) where T <: CTMatrixTypes
    ncols(G) == ncols(H) ||
        throw(ArgumentError("The number of columns of G and H should be the same (received ncols(G) = $(ncols(G)), ncols(H) = $(ncols(H)))"))
    base_ring(G) == base_ring(H) || throw(ArgumentError("G and H are not over the same field"))
    G_new = _remove_empty(G, :rows)
    H_new = _remove_empty(H, :rows)
    iszero(G_new * transpose(H_new)) || throw(ArgumentError("H isn't orthogonal to G"))
    G_stand, H_stand, P, k = _standard_form(G_new)
    rank(H) == ncols(G) - k || throw(ArgumentError("The given matrix H is not a parity check matrix for G"))

    ub1, _ = _min_wt_row(G_new)
    ub2, _ = _min_wt_row(G_stand)
    ub = min(ub1, ub2)
    C = LinearCode(base_ring(G_new), ncols(G_new), k, missing, 1, ub, G_new, H_new, G_stand, H_stand, P, missing)

    if brute_force_WE && BigInt(order(base_ring(G)))^min(k, ncols(G) - k) <= 1.5e5
        C.weight_enum = if 2k <= ncols(G)
            _weight_enumerator_BF(C.G_stand)
        else
            MacWilliams_identity(dual(C), _weight_enumerator_BF(C.H_stand))
        end
        d = minimum(filter(is_positive, first.(exponent_vectors(CWE_to_HWE(C.weight_enum).polynomial))))
        set_minimum_distance!(C, d)
    end

    return C
end

function LinearCode(G::Matrix{Int}, q::Int, parity::Bool = false)
    factors = Nemo.factor(q)
    (length(factors) == 1 && q > 1) || throw(ArgumentError("There is no finite field of order $q."))

    p, m = first(factors)
    F = m == 1 ? GF(p) : GF(p, m, :ω)
    G1 = matrix(F, G)
    rref!(G1)
    return LinearCode(G1, parity)
end

function LinearCode(Gs::Vector{<:CTMatrixTypes})
    s = size(Gs[1])
    all(s == size(Gs[i]) for i in 2:length(Gs)) || throw(ArgumentError("Not all vectors in `Gs` were the same size."))

    G = reduce(vcat, Gs)
    rref!(G)
    return LinearCode(G)
end

function LinearCode(Gs::Vector{Vector{Int}}, q::Int, parity::Bool = false)
    s = size(Gs[1])
    all(s == size(Gs[i]) for i in 2:length(Gs)) || throw(ArgumentError("Not all vectors in `Gs` were the same size."))

    return LinearCode(reduce(vcat, Gs), q, parity)
end

#############################
      # getter functions
#############################

"""
    field(C::AbstractLinearCode)

Return the base ring of the generator matrix.
"""
field(C::AbstractLinearCode) = C.F

"""
    length(C::AbstractLinearCode)

Return the length of the code.
"""
length(C::AbstractLinearCode) = C.n

"""
    dimension(C::AbstractLinearCode)

Return the dimension of the code.
"""
dimension(C::AbstractLinearCode) = C.k

"""
    cardinality(C::AbstractLinearCode)

Return the cardinality of the code.
"""
cardinality(C::AbstractLinearCode) = BigInt(order(C.F))^C.k

"""
    rate(C::AbstractLinearCode)

Return the rate, ``R = k/n``, of the code.
"""
rate(C::AbstractLinearCode) = C.k / C.n

"""
    generator_matrix(C::AbstractLinearCode, stand_form::Bool=false)

Return the generator matrix of the code. If the optional parameter `stand_form`
is set to `true`, the standard form of the generator matrix is returned instead.
"""
function generator_matrix(C::AbstractLinearCode, stand_form::Bool=false)
    if isa(C, QuasiCyclicCode)
        stand_form && !ismissing(C.G_stand) && (return C.G_stand;)
        if ismissing(C.G)
            if C.A_type == :G
                G = lift(C.A)
                C.G = G
            else
                G = kernel(lift(C.A), side = :right)
                rnk_G = rank(G)
                # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
                nr = nrows(G)
                G_tr = zero_matrix(base_ring(G), rnk_G, nr)
                for r in 1:nr
                    for c in 1:rnk_G
                        !iszero(G[r, c]) && (G_tr[c, r] = G[r, c];)
                    end
                end
                C.G = G_tr
            end
        end
        if stand_form
            C.G_stand, C.H_stand, C.P_stand, _ = _standard_form(C.G)
            return C.G_stand
        end
        return C.G
    else
        stand_form ? (return C.G_stand;) : (return C.G;)
    end
end

"""
    parity_check_matrix(C::AbstractLinearCode, stand_form::Bool=false)

Return the parity-check matrix of the code. If the optional parameter
`stand_form` is set to `true`, the standard form of the parity-check matrix
is returned instead.
"""
function parity_check_matrix(C::AbstractLinearCode, stand_form::Bool=false)
    if isa(C, QuasiCyclicCode)
        if stand_form
            ismissing(C.H_stand) || (return C.H_stand;)
            if ismissing(C.G)
                if C.A_type == :G
                    G = lift(C.A)
                    C.G = G
                else
                    G = kernel(lift(C.A), side = :right)
                    rnk_G = rank(G)
                    # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
                    nr = nrows(G)
                    G_tr = zero_matrix(base_ring(G), rnk_G, nr)
                    for r in 1:nr
                        for c in 1:rnk_G
                            !iszero(G[r, c]) && (G_tr[c, r] = G[r, c];)
                        end
                    end
                    C.G = G_tr
                end
            end
            C.G_stand, C.H_stand, C.P_stand, _ = _standard_form(C.G)
            return C.H_stand
        elseif ismissing(C.H)
            C.H = C.A_type == :H ? lift(C.A) : transpose(kernel(lift(C.A), side = :right))
        end
        return C.H
    elseif isa(C, LDPCCode)
        # TODO: fix this since LDPCCode no longer has a standard form or a C object
        stand_form ? (return C.C.H_stand;) : (return C.H;)
    else
        stand_form ? (return C.H_stand;) : (return C.H;)
    end
end

"""
    standard_form_permutation(C::AbstractLinearCode)

Return the permutation matrix required to permute the columns of the code mA_trices to have the same
row space as the mA_trices in standard form. Returns `missing` is no such permutation is required.
"""
standard_form_permutation(C::AbstractLinearCode) = C.P_stand

"""
    relative_distance(C::AbstractLinearCode)

Return the relative minimum distance, ``\\delta = d / n`` of the code if ``d`` is known;
otherwise return `missing`.
"""
relative_distance(C::AbstractLinearCode) = C.d / C.n

"""
    genus(C::AbstractLinearCode)

Return the genus, ``n + 1 - k - d``, of the code.
"""
genus(C::AbstractLinearCode) = C.n + 1 - C.k - C.d

"""
    minimum_distance_lower_bound(C::AbstractLinearCode)

Return the current lower bound on the minimum distance of `C`.
"""
minimum_distance_lower_bound(C::AbstractLinearCode) = C.l_bound

"""
    minimum_distance_upper_bound(C::AbstractLinearCode)

Return the current upper bound on the minimum distance of `C`.
"""
minimum_distance_upper_bound(C::AbstractLinearCode) = C.u_bound

"""
    is_MDS(C::AbstractLinearCode)

Return `true` if code is maximum distance separable (MDS).
"""
function is_MDS(C::AbstractLinearCode)
    ismissing(C.d) && minimum_distance(C)
    return C.d == Singleton_bound(C.n, C.k)
end

"""
    number_correctable_errors(C::AbstractLinearCode)

Return the number of correctable errors for the code.

# Notes
* The number of correctable errors is ``t = \\floor{(d - 1) / 2}``.
"""
number_correctable_errors(C::AbstractLinearCode) = ismissing(C.d) ? missing : Int(fld(C.d - 1, 2))

"""
    is_overcomplete(C::AbstractLinearCode, which::Symbol=:G)

Return `true` if the generator matrix is over complete, or if the optional parameter is
set to :H and the parity-check matrix is over complete.
"""
function is_overcomplete(C::AbstractLinearCode, which::Symbol=:G)
    if which == :G
        return nrows(C.G) > C.k
    elseif which == :H
        return nrows(C.H) > C.n - C.k
    end
    throw(ArgumentError("Received symbol $which, expected :G or :H"))
end

#############################
      # setter functions
#############################

"""
    set_distance_lower_bound!(C::AbstractLinearCode, l::Int)

Set the lower bound on the minimum distance of `C`, if `l` is better than the current bound.
"""
function set_distance_lower_bound!(C::AbstractLinearCode, l::Int)
    1 <= l <= C.u_bound || throw(DomainError(l, "The lower bound must be between 1 and the upper bound."))
    C.l_bound < l && (C.l_bound = l;)
    if C.l_bound == C.u_bound
        @info "The new lower bound is equal to the upper bound; setting the minimum distance."
        C.d = C.l_bound
    end
end

"""
    set_distance_upper_bound!!(C::AbstractLinearCode, u::Int)

Set the upper bound on the minimum distance of `C`, if `u` is better than the current bound.
"""
function set_distance_upper_bound!!(C::AbstractLinearCode, u::Int)
    C.l_bound <= u <= C.n || throw(DomainError(u, "The upper bound must be between the lower bound and the code length."))
    u < C.u_bound && (C.u_bound = u;)
    if C.l_bound == C.u_bound
        @info "The new upper bound is equal to the lower bound; setting the minimum distance."
        C.d = C.l_bound
    end
end

"""
    set_minimum_distance!(C::AbstractLinearCode, d::Int)

Set the minimum distance of the code to `d`.

# Notes
* The only check done on the value of `d` is that ``1 \\leq d \\leq n``.
"""
function set_minimum_distance!(C::AbstractLinearCode, d::Int)
    0 < d <= C.n || throw(DomainError("The minimum distance of a code must be ≥ 1; received: d = $d."))
    C.d = d
    C.l_bound = d
    C.u_bound = d
end

function change_field!(C::T, F::CTFieldTypes) where T <: AbstractLinearCode
    # TODO: this is doesn't work for cyclic codes yet
    T <: AbstractCyclicCode && @error "Not implemented for cyclic codes yet"
    C.G = change_base_ring(F, C.G)
    C.H = change_base_ring(F, C.H)
    C.G_stand = change_base_ring(F, C.G_stand)
    C.H_stand = change_base_ring(F, C.H_stand)
    ismissing(C.P_stand) || (C.P_stand = change_base_ring(F, C.P_stand);)

    if order(F) != order(C.F)
        # TODO: should be able to compute the new weight_enum easily as well
        C.weight_enum = missing
    end

    C.F = F
    
    return nothing
end

function change_field(C::AbstractLinearCode, F::CTFieldTypes)
    C2 = deepcopy(C)
    change_field!(C2, F)
    return C2
end

#############################
     # general functions
#############################
"""
    information_set(C::AbstractLinearCode)

Return an information set of C defined by the pivot column indexs of C.G. 

# Arguments
- `C` - a linear code 
"""
function information_set(C::AbstractLinearCode)
    _, _, perm = _rref_col_swap_perm(C.G)
    return on_sets(collect(1 : C.k), inv(perm)) 
end

"""
    random_information_set(C::AbstractLinearCode; rng::AbstractRNG = Random.seed!())

Return a random information set of C defined by pivot column indexs of C.G. 

# Arguments
- `C` - a linear code
- `rng` - random number generator
"""
function random_information_set(C::AbstractLinearCode; rng::AbstractRNG = Random.seed!())
    shuffle_perm_julia = shuffle(rng, collect(1 : C.n)) 
    shuffle_perm = perm(shuffle_perm_julia)
    # apply shuffle permutation
    permuted_mat = C.G[: , invperm(shuffle_perm_julia)] 
    # transform to rref. The first k columns of the rref matrix are pivot columns
    _, _, rref_perm = CodingTheory._rref_col_swap_perm(permuted_mat)
    inv_perm = inv(shuffle_perm * rref_perm)
    # apply the inverse permutation to the pivots
    return on_sets(collect(1 : C.k), inv_perm) 
end

"""
    random_linear_code(F::CTFieldTypes, n::Int, k::Int; rng::AbstractRNG = Random.seed!())

Return a random [n, k] linear code over F.

# Arguments
- `F` - finite field
- `n` - length of the code
- `k` - dimension of the code
- `rng` - random number generator
"""
function random_linear_code(F::CTFieldTypes, n::Int, k::Int; rng::AbstractRNG = Random.seed!())
    rand_mat = zero_matrix(F, k, n - k)
    for r in 1:nrows(rand_mat) 
        for c in 1:ncols(rand_mat) 
            rand_mat[r, c] = rand(rng, F) 
        end
    end
    full_mat = hcat(identity_matrix(F, k), rand_mat)
    return LinearCode(full_mat)
end 

"""
    random_linear_code(q::Int, n::Int, k::Int; rng::AbstractRNG = Random.seed!())

Return a random [n, k] linear code over GF(q).

# Arguments
- `q` - a prime power 
- `n` - length of the code
- `k` - dimension of the code
- `rng` - random number generator
"""
function random_linear_code(q::Int, n::Int, k::Int; rng::AbstractRNG = Random.seed!())
    is_pp, e, p = is_prime_power_with_data(q)
    if e == 1 
        field = Oscar.Nemo.Native.GF(p)
    elseif e > 1
        field = GF(p, e, :x)
    else
        throw(ArgumentError("q must be a prime power"))
    end
    return random_linear_code(field, n, k, rng = rng)
end 

function _standard_form(G::CTMatrixTypes)
    rnk, G_stand, P = _rref_col_swap(G, 1:nrows(G), 1:ncols(G))
    F = base_ring(G_stand)
    nrows(G_stand) > rnk && (G_stand = _remove_empty(G_stand, :rows);)
    A_tr = transpose(view(G_stand, :, (nrows(G_stand) + 1):ncols(G_stand)))
    H_stand = hcat(order(F) == 2 ? A_tr : -A_tr, identity_matrix(F, ncols(G_stand) - nrows(G_stand)))
    return G_stand, H_stand, P, rnk
end

# a bit odd to handle all subtypes in the supertype but saves hundreds
# of repeated lines and keeps uniform
function show(io::IO, C::AbstractLinearCode)
    print(io, "[$(C.n), $(C.k)")
    !ismissing(C.d) && print(io, ", $(C.d)")
    typeof(C) <: AbstractBCHCode && print(io, "; $(C.b)")
    print(io, "]_$(order(C.F)) ")
    if isa(C, ReedSolomonCode)
        println(io, "Reed-Solomon code")
    elseif isa(C, BCHCode)
        println(io, "BCH code")
    elseif isa(C, CyclicCode)
        println(io, "cyclic code")
    elseif isa(C, ReedMullerCode)
        println(io, "Reed-Muller code RM($(C.r), $(C.m))")
    elseif isa(C, QuasiCyclicCode)
        println(io, "quasi-cyclic code of index $(C.l)")
    elseif isa(C, GeneralizedReedSolomonCode)
        println(io, "generalized Reed-Solomon code")
    elseif isa(C, GoppaCode)
        println(io, "Goppa code")
    else
        println(io, "linear code")
    end

    if get(io, :compact, true)
        if typeof(C) <: AbstractCyclicCode
            println(io, "$(order(C.F))-Cyclotomic cosets: ")
            len = length(qcosets_reps(C))
            if len == 1
                println("\tC_$(qcosets_reps(C)[1])")
            else
                for (i, x) in enumerate(qcosets_reps(C))
                    if i == 1
                        print(io, "\tC_$x ∪ ")
                    elseif i == 1 && i == len
                        println(io, "\tC_$x")
                    elseif i != len
                        print(io, "C_$x ∪ ")
                    else
                        println(io, "C_$x")
                    end
                end
            end
            println(io, "Generator polynomial:")
            println(io, "\t", generator_polynomial(C))
        elseif isa(C, GoppaCode)
            len = length(C.L)
            if len ≤ 20
                println(io, "Evaluated at:")
                print(io, "\t[")
                for i in 1:len
                    if i ≠ len
                        print(C.L[i], ", ")
                    else
                        println(C.L[i], "]")
                    end
                end
            end
            println(io, "Goppa polynomial:")
            println(io, "\t", C.g)
        end

        if C.n ≤ 30 && C.k ≠ 0
            if isa(C, QuasiCyclicCode)
                if C.A_type == :G
                    M = generator_matrix(C)
                    nr, nc = size(M)
                    println(io, "Generator matrix: $(nr) × $(nc)")
                else
                    M = parity_check_matrix(C)
                    nr, nc = size(M)
                    println(io, "Parity-check matrix: $(nr) × $(nc)")
                end
                for i in 1:nr
                    print(io, "\t")
                    for j in 1:nc
                        if j != nc
                            print(io, "$(M[i, j]) ")
                        elseif j == nc && i != nr
                            println(io, "$(M[i, j])")
                        else
                            print(io, "$(M[i, j])")
                        end
                    end
                end
            else
                G = generator_matrix(C)
                nr, nc = size(G)
                println(io, "Generator matrix: $nr × $nc")
                println(io, "\t" * replace(repr(MIME("text/plain"), G), r"\n" => "\n\t"))
            end
        end
        # if !ismissing(C.weight_enum)
        #     println(io, "\nComplete weight enumerator:")
        #     print(io, "\t", polynomial(C.weight_enum))
        # end
    end
end

"""
    Singleton_bound(n::Int, a::Int)

Return the Singleton bound ``d \\leq n - k + 1`` or ``k \\leq n - d + 1`` depending on the interpretation of `a`.
"""
Singleton_bound(n::Int, a::Int) = 0 <= a <= n ? (return n - a + 1) : 
    error("Invalid parameters for the Singleton bound. Received n = $n, k/d = $a")

"""
    Singleton_bound(C::AbstractLinearCode)

Return the Singleton bound on the minimum distance of the code (``d \\leq n - k + 1``).
"""
Singleton_bound(C::AbstractLinearCode) = Singleton_bound(C.n, C.k)

"""
    encode(C::AbstractLinearCode, v::Union{CTMatrixTypes, Vector{Int}})

Return the encoding of `v` into `C`
"""
function encode(C::AbstractLinearCode, v::Union{CTMatrixTypes, Vector{Int}})
    w = isa(v, Vector{Int}) ? matrix(C.F, 1, length(v), v) : v
    G = generator_matrix(C)
    nr = nrows(G)
    (size(w) != (1, nr) && size(w) != (nr, 1)) &&
        throw(ArgumentError("Vector has incorrect dimension; expected length $nr, received: $(size(v))."))
    base_ring(w) == C.F || throw(ArgumentError("Vector must have the same base ring as the generator matrix."))
    nrows(w) != 1 || return w * G
    return transpose(w) * G
end

"""
    syndrome(C::AbstractLinearCode, v::Union{CTMatrixTypes, Vector{Int}})

Return the syndrome of `v` with respect to `C`.
"""
function syndrome(C::AbstractLinearCode, v::Union{CTMatrixTypes, Vector{Int}})
    w = isa(v, Vector{Int}) ? matrix(C.F, length(v), 1, v) : v
    H = parity_check_matrix(C)
    nc = ncols(H)
    (size(w) != (nc, 1) && size(w) != (1, nc)) &&
        throw(ArgumentError("Vector has incorrect dimension; expected length $nc, received: $(size(v))."))
    if base_ring(w) != C.F
        if order(base_ring(w)) == order(C.F)
            @warn "Fields are of different types, but have the same order."
        else
            throw(ArgumentError("Vector must have the same base ring as the parity-check matrix."))
        end
    end
    return nrows(w) == 1 ? H * transpose(w) : H * w
end

"""
    in(v::Union{CTMatrixTypes, Vector{Int}}, C::AbstractLinearCode)
    ∈(v::Union{CTMatrixTypes, Vector{Int}}, C::AbstractLinearCode)

Return whether or not `v` is a codeword of `C`.
"""
in(v::Union{CTMatrixTypes, Vector{Int}}, C::AbstractLinearCode) = iszero(syndrome(C, v))

"""
    ⊆(C1::AbstractLinearCode, C2::AbstractLinearCode)
    ⊂(C1::AbstractLinearCode, C2::AbstractLinearCode)
    is_subcode(C1::AbstractLinearCode, C2::AbstractLinearCode)

Return whether or not `C1` is a subcode of `C2`.
"""
function ⊆(C1::AbstractLinearCode, C2::AbstractLinearCode)
    C1.F == C2.F || (order(C1.F) == order(C2.F) ? (@warn "Fields are of different types, but have the same order.") : (return false;))
    (C1.n == C2.n && C1.k <= C2.k) || return false

    G1 = generator_matrix(C1)
    return all(view(G1, r:r, 1:C1.n) ∈ C2 for r in axes(G1, 1))
end
⊂(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊆ C2
is_subcode(C1::AbstractLinearCode, C2::AbstractLinearCode) = C1 ⊆ C2

"""
    dual(C::AbstractLinearCode)
    Euclidean_dual(C::AbstractLinearCode)

Return the (Euclidean) dual of the code `C`.
"""
function dual(C::AbstractLinearCode)
    if typeof(C) <: AbstractCyclicCode
        return CyclicCode(Int(order(C.F)), C.n, dual_qcosets(Int(order(C.F)), C.n, C.qcosets))
    elseif isa(C, GeneralizedReedSolomonCode)
        d = C.k + 1
        return GeneralizedReedSolomonCode(C.F, C.n, C.n - C.k, d, d, d,
            deepcopy(C.dual_scalars), deepcopy(C.scalars), deepcopy(C.eval_pts),
            deepcopy(C.H), deepcopy(C.G), deepcopy(C.H_stand),
            deepcopy(C.G_stand), deepcopy(C.P_stand), missing)
    elseif isa(C, MatrixProductCode)
        nr, nc = size(C.A)
        # probably not going to work
        nr == nc || return LinearCode.dual(C)
        D = Vector{LinearCode}()
        for i in 1:length(C.C)
            push!(D, dual(C.C[i]))
        end
        
        try
            A_inv = inv(C.A)
        catch
            return LinearCode.dual(C)
        end
        return MatrixProductCode(D, transpose(A_inv))
    elseif isa(C, ReedMullerCode)
        d = 2^(C.r + 1)
        return ReedMullerCode(C.F, C.n, C.n - C.k, d, d, d, C.m - C.r - 1, C.m, C.H, C.G,
            C.H_stand, C.G_stand, C.P_stand, missing)
    else
        G = deepcopy(generator_matrix(C))
        H = deepcopy(parity_check_matrix(C))

        H_stand, G_stand, P, _ = _standard_form(C.H)
        P_stand_old = ismissing(C.P_stand) ? identity_matrix(C.F, C.n) : C.P_stand
        P_stand = ismissing(P) ? P_stand_old : P_stand_old * P

        # Possible alternative to above, might be faster:
        # G_standold = generator_matrix(C, true)
        # G_stand = hcat(view(G_standold, :, C.k + 1:C.n), view(G_standold, :, 1:C.k))
        # H_standold = parity_check_matrix(C, true)
        # H_stand = hcat(view(H_standold, :, C.k + 1:C.n), view(H_standold, :, 1:C.k))
        # P_stand_old = ismissing(C.P_stand) ? identity_matrix(C.F, C.n) : C.P_stand
        # P_stand = vcat(view(C.P_stand, C.k + 1:C.n, :), view(C.P_stand, 1:C.k, :))

        if !ismissing(C.weight_enum)
            dual_wt_enum = MacWilliams_identity(C, C.weight_enum)
            dual_HWE_poly = CWE_to_HWE(dual_wt_enum).polynomial
            d = minimum(filter(>(0), first.(exponent_vectors(dual_HWE_poly))))
            return LinearCode(C.F, C.n, C.n - C.k, d, d, d, H, G,
                H_stand, G_stand, P_stand, dual_wt_enum)
        else
            ub1, _ = _min_wt_row(H)
            ub2, _ = _min_wt_row(H_stand)
            ub = min(ub1, ub2)
            return LinearCode(C.F, C.n, C.n - C.k, missing, 1, ub, H,
                            G, H_stand, G_stand, P_stand, missing)
        end
    end
end
Euclidean_dual(C::AbstractLinearCode) = dual(C)

"""
    Hermitian_dual(C::AbstractLinearCode)

Return the Hermitian dual of a code defined over a quadratic extension.
"""
function Hermitian_dual(C::AbstractLinearCode)
    if isa(C, MatrixProductCode)
        # the inner functions here should complain if not quadratic, so don't have to check here
        nr, nc = size(C.A)
        # probably not going to work
        nr == nc || return LinearCode.Hermitian_dual(C)
        D = Vector{LinearCode}()
        for i in 1:length(C.C)
            push!(D, Hermitian_dual(C.C[i]))
        end

        try
            A_inv = inv(Hermitian_conjugate_matrix(C.A))
        catch
            return LinearCode.dual(C)
        end
        return MatrixProductCode(D, transpose(A_inv))
    else
        return dual(LinearCode(Hermitian_conjugate_matrix(C.G)))
    end
end

"""
    are_equivalent(C1::AbstractLinearCode, C2::AbstractLinearCode)

Return `true` if `C1 ⊆ C2` and `C2 ⊆ C1`.
"""
are_equivalent(C1::AbstractLinearCode, C2::AbstractLinearCode) = (C1 ⊆ C2) && (C2 ⊆ C1)

"""
    are_permutation_equivalent(C1::AbstractLinearCode, C2::AbstractLinearCode)

Returns a `Tuple{Bool, Union{Missing, Perm}}`. The first entry is
`true` if `C1` and `C2` are permutation equivalent codes. The second
entry gives the permutation that, when applied to the columns of the
generator matrix of C1, would give C2.
"""
function are_permutation_equivalent(C1::AbstractLinearCode, C2::AbstractLinearCode)
    G1 = generator_matrix(C1, true)
    G2 = generator_matrix(C2, true)
    if G1 == G2
        if ismissing(C1.P_stand) && ismissing(C2.P_stand)
            (true, Perm(1:C1.n))
        elseif ismissing(C1.P_stand)
            (true, Perm(Int.(data.(transpose(C2.P_stand))) * collect(1:C1.n)))
        elseif ismissing(C2.P_stand)
            (true, Perm(Int.(data.(C1.P_stand)) * collect(1:C1.n)))
        else
            (true, Perm(Int.(data.(transpose(C2.P_stand) * C1.P_stand)) * collect(1:C1.n)))
        end
    else
        (false, missing)
    end
end

"""
Checks permutation equivalence by brute force for the purpose of doing small tests. 
"""
function _are_perm_equivalent_exhaustive_search(C1::AbstractLinearCode, C2::AbstractLinearCode)
    G1 = C1.G 
    G2 = C2.G 
    if G1 == G2 
        return are_permutation_equivalent(C1, C2) 
    end 
    nc = ncols(G1)
    sym = symmetric_group(nc)
    for e in collect(sym) 
        P = permutation_matrix(C1.F, e)
        if G1 * P == G2 
            return (true, P)
        end
    end
    return (false, missing)
end

"""
    is_self_dual(C::AbstractLinearCode)

Return `true` if `are_equivalent(C, dual(C))`.
"""
is_self_dual(C::AbstractLinearCode) = are_equivalent(C, dual(C))

"""
    is_self_orthogonal(C::AbstractLinearCode)
    is_weakly_self_dual(C::AbstractLinearCode)
    is_Euclidean_self_orthogonal(C::AbstractLinearCode)

Return `true` if `C ⊆ dual(C)`.
"""
is_self_orthogonal(C::AbstractLinearCode) = C ⊆ dual(C)
is_weakly_self_dual(C::AbstractLinearCode) = is_self_orthogonal(C)
is_Euclidean_self_orthogonal(C::AbstractLinearCode) = is_self_orthogonal(C)

"""
    is_dual_containing(C::AbstractLinearCode)
    is_Euclidean_dual_containing(C::AbstractLinearCode)

Return `true` if `dual(C) ⊆ C`.
"""
is_dual_containing(C::AbstractLinearCode) = dual(C) ⊆ C
is_Euclidean_dual_containing(C::AbstractLinearCode) = is_dual_containing(C)

"""
    is_Hermitian_self_dual(C::AbstractLinearCode)

Return `true` if `are_equivalent(C, Hermitian_dual(C))`.
"""
is_Hermitian_self_dual(C::AbstractLinearCode) = are_equivalent(C, Hermitian_dual(C))

"""
    is_Hermitian_self_orthogonal(C::AbstractLinearCode)
    is_Hermitian_weakly_self_dual(C::AbstractLinearCode)

Return `true` if `C ⊆ Hermitian_dual(C)`.
"""
is_Hermitian_self_orthogonal(C::AbstractLinearCode) = C ⊆ Hermitian_dual(C)
is_Hermitian_weakly_self_dual(C::AbstractLinearCode) = is_Hermitian_self_orthogonal(C)

"""
    is_Hermitian_dual_containing(C::AbstractLinearCode)

Return `true` if `Hermitian_dual(C) ⊆ C`.
"""
is_Hermitian_dual_containing(C::AbstractLinearCode) = Hermitian_dual(C) ⊆ C

# TODO: add l-Galois dual and self-dual/orthogonal based functions
# LGaloisdual(C::AbstractLinearCode) = dual(LinearCode(LGaloisconjugatematrix(C.G)))
# need LGaloisconjugatematrix(C.G)

"""
    characteristic_polynomial(C::AbstractLinearCode)

Return the characteristic polynomial of `C`.

# Notes
* The characteristic polynomial is defined in [Lin1999]_
"""
function characteristic_polynomial(C::AbstractLinearCode)
    _, x = polynomial_ring(Nemo.QQ, :x)
    D = dual(C)
    supD = support(D)
    q = Int(order(C.F))
    return q^(n - k) * prod(1 - x / j for j in supD if j > 0)
end

"""
    vector_space(C::AbstractLinearCode)

Return the code `C` as a vector space object.
"""
function vector_space(C::AbstractLinearCode)
    V = vector_space(C.F, C.n)
    G = generator_matrix(C)
    return sub(V, [V(view(G, i:i, :)) for i in 1:nrows(G)])
end
# vector_space(C::AbstractLinearCode) = vector_space(C)

"""
    is_even(C::AbstractLinearCode)

Return `true` if `C` is even.
"""
function is_even(C::AbstractLinearCode)
    Int(order(C.F)) == 2 || throw(ArgumentError("Even-ness is only defined for binary codes."))
    
    # A binary code generated by G is even if and only if each row of G has
    # even weight.
    G = generator_matrix(C)
    return all(wt(view(G, r:r, :)) % 2 == 0 for r in 1:nrows(G))
end

"""
    is_doubly_even(C::AbstractLinearCode)

Return `true` if `C` is doubly-even.
"""
function is_doubly_even(C::AbstractLinearCode)
    Int(order(C.F)) == 2 || throw(ArgumentError("Even-ness is only defined for binary codes."))

    # A binary code generated by G is doubly-even if and only if each row of G
    # has weight divisible by 4 and the sum of any two rows of G has weight
    # divisible by 4.
    G = generator_matrix(C)
    nr = nrows(G)
    all(wt(view(G, r:r, :)) % 4 == 0 for r in 1:nr) || (return false;)
    all(wt(view(G, r1:r1, :) + view(G, r2:r2, :)) % 4 == 0
        for r1 in 1:nr, r2 in 1:nr) || (return false;)
    return true
end

"""
    is_triply_even(C::AbstractLinearCode)

Return `true` if `C` is triply-even.
"""
function is_triply_even(C::AbstractLinearCode)
    Int(order(C.F)) == 2 || throw(ArgumentError("Even-ness is only defined for binary codes."))

    # following Ward's divisibility theorem
    G = _Flint_matrix_to_Julia_int_matrix(generator_matrix(C))
    nr, _ = size(G)
    all(wt(view(G, r:r, :)) % 8 == 0
        for r in 1:nr) || (return false;)
    all(wt(view(G, r1:r1, :) .* view(G, r2:r2, :)) % 4 == 0
        for r1 in 1:nr, r2 in 1:nr) || (return false;)
    all(wt(view(G, r1:r1, :) .* view(G, r2:r2, :) .* view(G, r3:r3, :)) % 2 == 0
        for r1 in 1:nr, r2 in 1:nr, r3 in 1:nr) || (return false;)
    return true
end

"""
    words(C::AbstractLinearCode, only_print::Bool=false)
    codewords(C::AbstractLinearCode, only_print::Bool=false)
    elements(C::AbstractLinearCode, only_print::Bool=false)

Return the elements of `C`.

# Notes
* If `only_print` is `true`, the elements are only printed to the console and not
  returned.
"""
function words(C::AbstractLinearCode, only_print::Bool=false)
    words = only_print ? nothing : Vector{typeof(C.G)}()
    G = ismissing(C.P_stand) ? generator_matrix(C, true) : generator_matrix(C, true) * C.P_stand
    E = base_ring(G)

    if iszero(G)
        row = zero_matrix(E, 1, C.n)
        only_print ? println(row) : push!(words, row)
        # TODO: looks pretty returned to me
        return words
    end

    # TODO new bug?
    # for iter in Iterators.product(Iterators.repeated(E, nrows(G))...)
    for iter in Nemo.AbstractAlgebra.ProductIterator([E for _ in 1:nrows(G)], inplace = true)
        row = iter[1] * view(G, 1:1, :)
        for r in 2:nrows(G)
            if !iszero(iter[r])
                row += iter[r] * view(G, r:r, :)
            end
        end
        only_print ? println(row) : push!(words, row)
    end
    return words
end
codewords(C::AbstractLinearCode, only_print::Bool = false) = words(C, only_print)
elements(C::AbstractLinearCode, only_print::Bool = false) = words(C, only_print)

"""
    hull(C::AbstractLinearCode)
    Euclidean_hull(C::AbstractLinearCode)

Return the (Euclidean) hull of `C` and its dimension.

# Notes
* The hull of a code is the intersection of it and its dual.
"""
function hull(C::AbstractLinearCode)
    G = generator_matrix(C)
    H = parity_check_matrix(C)
    F = field(C)
    VS = vector_space(F, C.n)
    U, U_to_VS = sub(VS, [VS(G[i, :]) for i in 1:nrows(G)])
    W, _ = sub(VS, [VS(H[i, :]) for i in 1:nrows(H)])
    I, I_to_W = intersect(U, W)
    if !iszero(AbstractAlgebra.dim(I))
        I_basis = [U_to_VS(I_to_W(g)) for g in gens(I)]
        G_I = reduce(vcat, I_basis)
        F_basis = [[F(G_I[j][i]) for i in 1:C.n] for j in 1:AbstractAlgebra.dim(I)]
        G_hull = matrix(F, length(F_basis), length(F_basis[1]), reduce(vcat, F_basis))
        return LinearCode(G_hull), AbstractAlgebra.dim(I)
    else
        return missing, 0 # is this the behavior I want?
    end
end
Euclidean_hull(C::AbstractLinearCode) = hull(C)

"""
    Hermitian_hull::AbstractLinearCode)
    
Return the Hermitian hull of `C` and its dimension.

# Notes
* The Hermitian hull of a code is the intersection of it and its Hermitian dual.
"""
function Hermitian_hull(C::AbstractLinearCode)
    D = Hermitian_dual(C)
    G = generator_matrix(C)
    H = parity_check_matrix(D)
    F = field(C)
    VS = vector_space(F, C.n)
    U, U_to_VS = sub(VS, [VS(G[i, :]) for i in 1:nrows(G)])
    W, _ = sub(VS, [VS(H[i, :]) for i in 1:nrows(H)])
    I, I_to_W = intersect(U, W)
    if !iszero(AbstractAlgebra.dim(I))
        I_basis = [U_to_VS(I_to_W(g)) for g in gens(I)]
        G_I = reduce(vcat, I_basis)
        F_basis = [[F(G_I[j][i]) for i in 1:C.n] for j in 1:AbstractAlgebra.dim(I)]
        G_hull = matrix(F, length(F_basis), length(F_basis[1]), reduce(vcat, F_basis))
        return LinearCode(G_hull), AbstractAlgebra.dim(I)
    else
        return missing, 0 # is this the behavior I want?
    end
end

"""
    is_LCD(C::AbstractLinearCode)

Return `true` if `C` is linear complementary dual.

# Notes
* A code is linear complementary dual if the dimension of `hull(C)` is zero.
"""
function is_LCD(C::AbstractLinearCode)
    # return true if dim(hull) = 0
    _, dim_hull_C = hull(C)
    dim_hull_C == 0 ? (return true;) : (return false;)
end

"""
    is_Hermitian_LCD(C::AbstractLinearCode)

Return `true` if `C` is linear complementary Hermitian dual.

# Notes
* A code is linear complementary Hermitian dual if the dimension of `Hermitian_hull(C)` is zero.
"""
function is_Hermitian_LCD(C::AbstractLinearCode)
    # return true if dim(hull) = 0
    _, dim_hull_C = Hermitian_hull(C)
    dim_hull_C == 0 ? (return true;) : (return false;)
end

# TODO: add l-Galois hull and is_LCD functions for that dual

"""
    contains_self_dual_subcode(C::AbstractLinearCode)

Return `true` if `C` contains a self-dual subcode.
"""
function contains_self_dual_subcode(C::AbstractLinearCode)
    # conditions in https://arxiv.org/abs/2310.16504v1
    q = Int(order(C.F))
    factors = Nemo.factor(q)
    (p, t), = factors

    if p == 2 # also works for powers of 2
        return iseven(C.n) && is_self_orthogonal(dual(C))
    elseif q % 4 == 1 # is it necessary to require k >= n / 2 here?
        return iseven(C.n) && is_self_orthogonal(dual(C))
    elseif q % 4 == 3 # is it necessary to require k >= n / 2 here?
        return C.n % 4 == 0 && is_self_orthogonal(dual(C))
    else
        error("Unknown case for finite field of order $p^$t")
    end
end
