# Copyright (c) 2021, 2023 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
     # general functions
#############################

"""
    ord(n::Int, q::Int)

Return the order of `n` mod `q`.
"""
function ord(n::Int, q::Int)
    (q <= 0 || n <= 0) && 
        throw(DomainError("q and n both need to be positive. Passed: q = $q, n = $n"))

    # finite stop instead of while
    for i in 1:200
        if mod(BigInt(q)^i, n) == 1
            return i
        end
    end
    error("Unable to compute ord($n, $q).")
end

"""
    cyclotomic_coset(x::Int, q::Int, n::Int; to_sort::Bool=true, verbose::Bool=false)

Return the `q`-cyclotomic coset of `x` modulo `n`.

# Notes
* If the optional parameter `to_sort` is set to `false`, the result will not be
sorted. If the optional parameter `verbose` is set to `true`, the result will
pretty print.
"""
function cyclotomic_coset(x::Int, q::Int, n::Int; to_sort::Bool = true,
    verbose::Bool = false)

    temp = [mod(x, n)]
    for i in 0:(n - 1)
        y = mod(temp[end] * q, n)
        if y ∉ temp
            append!(temp, y)
        else
            break
        end
    end
    if to_sort
        sort!(temp)
    end

    if verbose
        len = length(temp)
        print("C_$x = {")
        for (i, y) in enumerate(temp)
            if i != len
                print("$y, ")
            else
                println("$y}")
            end
        end
    end
    return temp
end

"""
    all_cyclotomic_cosets(q::Int, n::Int; to_sort::Bool=true, verbose::Bool=false)

Return all `q`-cyclotomic cosets modulo `n`.

# Notes
* If the optional parameter `to_sort` is set to `false`, the result will not be
sorted. If the optional parameter `verbose` is set to `true`, the result will
pretty print.
"""
function all_cyclotomic_cosets(q::Int, n::Int; to_sort::Bool = true,
    verbose::Bool = false)

    n % q == 0 && throw(DomainError("Cyclotomic coset requires gcd(n, q) = 1"))

    arr = [[0]]
    for x in 1:(n - 1)
        found = false
        for a in arr
            if x ∈ a
                found = true
                break
            end
        end

        if !found
            Cx = cyclotomic_coset(x, q, n, to_sort = to_sort, verbose = false)
            push!(arr, Cx)
        end
    end

    # sort!(arr, by=x->x[1])

    if verbose
        for Cx in arr
            len = length(Cx)
            print("C_$(Cx[1]) = {")
            for (i, y) in enumerate(Cx)
                if i != len
                    print("$y, ")
                else
                    println("$y}")
                end
            end
        end
    end

    if sort!(reduce(vcat, arr)) != [i for i in 0:(n - 1)]
        error("Missed some")
    end
    return arr
end

"""
    complement_qcosets(q::Int, n::Int, qcosets::Vector{Vector{Int64}})

Return the complement of the `q`-cyclotomic cosets modulo `n` of `qcosets`.
"""
function complement_qcosets(q::Int, n::Int, qcosets::Vector{Vector{Int64}})
    all = all_cyclotomic_cosets(q, n)
    comp_cosets = Vector{Vector{Int64}}()
    for a in all
        # if a != [0]
            found = false
            for b in qcosets
                if a[1] == b[1]
                    found = true
                    break
                end
            end
            found || (push!(comp_cosets, a);)
        # end
    end
    return comp_cosets
end

"""
    qcoset_pairings(arr::Vector{Vector{Int64}}, n::Int)

Return the `q`-cyclotomic cosets modulo `n` collected into complementary pairs.
"""
function qcoset_pairings(arr::Vector{Vector{Int64}}, n::Int)
    coset_rep_list = Vector{Tuple{Int64, Int64}}()
    coset_pair_list = Vector{Tuple{Vector{Int64}, Vector{Int64}}}()
    for a in arr
        found = false
        for pair in coset_rep_list
            if a[1] == pair[1] || a[1] == pair[2]
                found = true
            end
        end

        if !found
            neg = sort!([mod(n - i, n) for i in a])
            if neg == a
                push!(coset_rep_list, (a[1], a[1]))
                push!(coset_pair_list, (a, a))
            else
                for b in arr
                    if neg == b
                        push!(coset_rep_list, (a[1], b[1]))
                        push!(coset_pair_list, (a, b))
                        break
                    end
                end
            end
        end
    end
    return coset_pair_list, coset_rep_list
end
qcoset_pairings(q::Int, n::Int) = qcoset_pairings(all_cyclotomic_cosets(q, n, to_sort = false), n)

# TODO: redo this with an abstract range
"""
    qcoset_table(a::Int, b::Int, q::Int)

Print all `q`-cyclotomic cosets modulo `n` for `n` between `a` and `b`.
"""
function qcoset_table(a::Int, b::Int, q::Int)
    for n in a:b
        if n % q != 0
            println("n = $n")
            all_cyclotomic_cosets(q, n, to_sort = true, verbose = true)
            println(" ")
        end
    end
end

"""
    dual_qcosets(q::Int, n::Int, qcosets::Vector{Vector{Int64}})

Return the dual of the `q`-cyclotomic cosets modulo `n` of `qcosets`.
"""
function dual_qcosets(q::Int, n::Int, qcosets::Vector{Vector{Int64}})
    comp_cosets = complement_qcosets(q, n, qcosets)
    for a in comp_cosets
        for (i, x) in enumerate(a)
            a[i] = mod(n - x, n)
        end
        sort!(a)
    end
    return comp_cosets
end

function _coerce_Kx_to_Ky_x(p::PolyRingElem{T}, A) where {T <: RingElement}
    # coerce p(x) in K[x] into A = Ky[x] (where Ky = K[y]) by coefficient embedding.
    Kx = parent(p)
    xA = gen(A)
    coeffs = Oscar.coefficients(p)  
    q = zero(A)
    for i in 0:length(coeffs)
        q += A(coeffs[i]) * xA^i
    end
    return q
end

function _composed_product(
    f::PolyRingElem{T},
    g::PolyRingElem{T}
) where {T <: RingElement}
    # computes the composed-product of two univariate polynomials using the resultant 

    # the multivariate resultant called here only accepts AbstractAlgebra.Generic.Poly type as input
    Kx = parent(f)
    parent(g) === Kx || throw(ArgumentError("f and g must have the same parent K[x]."))
    K = base_ring(Kx)
    d = degree(g)
    d < 0 && throw(ArgumentError("g must be nonzero."))
    Kx, x = polynomial_ring(K, :x)
    A,  y = polynomial_ring(Kx, :y)  # resultant(f, g) now eliminates y
    fA = _coerce_Kx_to_Ky_x(f, A)
    gA = _coerce_Kx_to_Ky_x(g, A)
    fy = evaluate(fA, y)

    fy = (fy + zero(A)) # converts the type of fy to AbstractAlgebra.Generic.Poly
    # gy_scaled = x^d * g(y/x) = sum_{i=0}^d a_i * y^i * x^(d-i)
    gy_scaled = zero(A)
    for i in 0:d
        ai = coeff(gA, i)           
        gy_scaled += ai * y^i * x^(d - i)
    end
    res = resultant(fy, gy_scaled)
    return res
end

"""
    construct_field(f::PolyRingElem{T}, g::PolyRingElem{T}) where {T <: RingElement}

constructs the smallest finite field extension where f and g split simultaneously
"""
function construct_field(f::PolyRingElem{T}, g::PolyRingElem{T}) where {T <: RingElement}
    xfacs = [x[1] for x in factor(f)]
    yfacs = [y[1] for y in factor(g)]
    comp_prods = [_composed_product(f1,f2) for f1 in xfacs, f2 in yfacs]
    comp_prods = [begin
        facs = collect(factor(c))                
        facs[argmax(degree(f[1]) for f in facs)][1]
    end for c in comp_prods] # select irreducible factor of highest degree
    dgs = unique([degree(c) for c in comp_prods])

    m = 1 # lcm of the degrees
    for dg in dgs
        m = lcm(m, dg)
    end
    # |K| = 2^k where k is the smallest integer with 2^k=1 (mod lcm(d_i))
    ZZm, _ = residue_ring(Nemo.ZZ, m) 
    two = ZZm(2)
    one = ZZm(1)
    i = 0
    for i in 0:m
        if two^i == one
            break
        end
    end
    if i == m
        throw(Error("failed to construct finite field"))
    end
    return Oscar.GF(2, m, :α), m, comp_prods
end