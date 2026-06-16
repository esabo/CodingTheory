# Copyright (c) 2021, 2023 - 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.
#############################
     # general functions
#############################

"""
$(TYPEDSIGNATURES)

Return the multiplicative order of `q` mod `n`.
"""
function ord(n::Int, q::Int)
    (q <= 0 || n <= 0) && throw(DomainError((q, n), "q and n both need to be positive."))
    gcd(n, q) == 1 || throw(ArgumentError("n and q must be coprime to compute multiplicative order (gcd($n, $q) != 1)."))

    val = mod(q, n)
    t = 1
    # Euler's Totient Theorem guarantees this will terminate in <= n steps
    while val != 1
        val = mod(val * q, n)
        t += 1
    end
    return t
end

"""
$(TYPEDSIGNATURES)

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
$(TYPEDSIGNATURES)

Return all `q`-cyclotomic cosets modulo `n`.
# Notes
* If the optional parameter `to_sort` is set to `false`, the result will not be
sorted.
If the optional parameter `verbose` is set to `true`, the result will
pretty print.
"""
function all_cyclotomic_cosets(q::Int, n::Int; to_sort::Bool = true,
    verbose::Bool = false)

    gcd(n, q) == 1 || throw(DomainError((n, q), "Cyclotomic cosets require gcd(n, q) = 1"))

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
$(TYPEDSIGNATURES)

Return the complement of the `q`-cyclotomic cosets modulo `n` of `qcosets`.
"""
function complement_qcosets(q::Int, n::Int, qcosets::Vector{Vector{Int64}})
    all = all_cyclotomic_cosets(q, n)
    comp_cosets = Vector{Vector{Int64}}()
    for a in all
        found = false
        for b in qcosets
            if a[1] == b[1]
                found = true
                break
            end
        end
        found || (push!(comp_cosets, a);)
    end
    return comp_cosets
end

"""
$(TYPEDSIGNATURES)

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
                break
            end
        end

        if !found
            neg_sorted = sort([mod(n - i, n) for i in a])
            a_sorted = sort(copy(a))
            
            # Safely compare set equality regardless of generation order
            if neg_sorted == a_sorted
                push!(coset_rep_list, (a[1], a[1]))
                push!(coset_pair_list, (a, a))
            else
                for b in arr
                    b_sorted = sort(copy(b))
                    if neg_sorted == b_sorted
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

"""
$(TYPEDSIGNATURES)

Print all `q`-cyclotomic cosets modulo `n` for `n` between `a` and `b`.
"""
function qcoset_table(a::Int, b::Int, q::Int)
    for n in a:b
        if gcd(n, q) == 1
            println("n = $n")
            all_cyclotomic_cosets(q, n, to_sort = true, verbose = true)
            println(" ")
        end
    end
end

"""
$(TYPEDSIGNATURES)

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

"""
$(TYPEDSIGNATURES)

Return the minimal polynomial of `α` defined by the `q`-cyclotomic coset `coset`.
# Notes
* The minimal polynomial is computed over the parent field of `α`, but mathematically 
  its coefficients are guaranteed to lie in the base field `GF(q)`.
"""
function minimal_polynomial(coset::Vector{Int}, α::CTFieldElem)
    E = parent(α)
    R, z = polynomial_ring(E, "z")
    
    M = one(R)
    for j in coset
        M *= (z - α^j)
    end
    
    return M
end

"""
$(TYPEDSIGNATURES)

Return `true` if `x` and `y` are conjugates over the subfield of order `q`.
"""
function are_conjugates(x::CTFieldElem, y::CTFieldElem, q::Int)
    parent(x) == parent(y) || return false
    
    # x and y are conjugates if y = x^(q^i) for some integer i
    E = parent(x)
    n_max = degree(E) # They must map to each other within the extension degree
    
    curr = x
    for _ in 1:n_max
        curr == y && return true
        curr = curr^q
    end
    return false
end