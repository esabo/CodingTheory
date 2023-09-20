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
    for i = 1:200
        if mod(BigInt(q)^i, n) == 1
            return i
        end
    end
    error("Unable to compute ord($n, $q).")
end

"""
    cyclotomic_coset(x::Int, q::Int, n::Int, to_sort::Bool=true, verbose::Bool=false)

Return the `q`-cyclotomic coset of `x` modulo `n`.

# Notes
* If the optional parameter `to_sort` is set to `false`, the result will not be
sorted. If the optional parameter `verbose` is set to `true`, the result will
pretty print.
"""
function cyclotomic_coset(x::Int, q::Int, n::Int, to_sort::Bool=true,
    verbose::Bool=false)

    temp = [mod(x, n)]
    for i = 0:(n - 1)
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
    all_cyclotomic_cosets(q::Int, n::Int, to_sort::Bool=true, verbose::Bool=false)

Return all `q`-cyclotomic cosets modulo `n`.

# Notes
* If the optional parameter `to_sort` is set to `false`, the result will not be
sorted. If the optional parameter `verbose` is set to `true`, the result will
pretty print.
"""
function all_cyclotomic_cosets(q::Int, n::Int, to_sort::Bool=true,
    verbose::Bool=false)

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
            Cx = cyclotomic_coset(x, q, n, to_sort, false)
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

qcoset_pairings(q::Int, n::Int) = qcoset_pairings(all_cyclotomic_cosets(q, n, false), n)

function qcoset_table(a::Int, b::Int, q::Int)
    for n in a:b
        if n % q != 0
            println("n = $n")
            all_cyclotomic_cosets(q, n, true)
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
            a[i] = mod(n - a[i], n)
        end
        sort!(a)
    end
    return comp_cosets
end
