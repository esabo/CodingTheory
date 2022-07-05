# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
    ord(n::Integer, q::Integer)

Return the order of `n` mod `q`.
"""
function ord(n::Integer, q::Integer)
    if q <= 0 || n <= 0
        error("q and n both need to be positive. Passed: q = $q, n = $n")
    end

    # finite stop instead of while
    for i = 1:200
        if mod(BigInt(q)^i, n) == 1
            return i
        end
    end

    error("Unable to compute ord($n, $q).")
end

"""
    cyclotomiccoset(x::Integer, q::Integer, n::Integer, tosort::Bool=true, verbose::Bool=false)

Return the `q`-cyclotomic coset of `x` modulo `n`.

If the optional parameter `tosort` is set to `false`, the result will not be
sorted. If the optional parameter `verbose` is set to `true`, the result will
pretty print.
"""
function cyclotomiccoset(x::Integer, q::Integer, n::Integer, tosort::Bool=true,
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
    if tosort
        sort!(temp)
    end
    len = length(temp)

    if verbose
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
    allcyclotomiccosets(q::Integer, n::Integer, tosort::Bool=true, verbose::Bool=false)

Return all `q`-cyclotomic cosets modulo `n`.

If the optional parameter `tosort` is set to `false`, the result will not be
sorted. If the optional parameter `verbose` is set to `true`, the result will
pretty print.
"""
function allcyclotomiccosets(q::Integer, n::Integer, tosort::Bool=true,
    verbose::Bool=false)

    if n % q == 0
        error("Cyclotomic coset requires gcd(n, q) = 1")
    end

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
            Cx = cyclotomiccoset(x, q, n, tosort, false)
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

    if sort!(vcat(arr...)) != [i for i in 0:(n - 1)]
        error("Missed some")
    end

    return arr
end

"""
    complementqcosets(q::Integer, n::Integer, qcosets::Vector{Vector{Int64}})

Return the complement of the `q`-cyclotomic cosets modulo `n` of `qcosets`.
"""
function complementqcosets(q::Integer, n::Integer, qcosets::Vector{Vector{Int64}})
    all = allcyclotomiccosets(q, n)
    compcosets = Vector{Vector{Int64}}()
    for a in all
        # if a != [0]
            found = false
            for b in qcosets
                if a[1] == b[1]
                    found = true
                    break
                end
            end

            if !found
                push!(compcosets, a)
            end
        # end
    end

    return compcosets
end

"""
    qcosetpairings(arr::Vector{Vector{Int64}}, n::Integer)

Return the .
"""
function qcosetpairings(arr::Vector{Vector{Int64}}, n::Integer)
    cosetreplist = Vector{Tuple{Int64, Int64}}()
    cosetpairlist = Vector{Tuple{Vector{Int64}, Vector{Int64}}}()
    for a in arr
        found = false
        for pair in cosetreplist
            if a[1] == pair[1] || a[1] == pair[2]
                found = true
            end
        end

        if !found
            neg = sort!([mod(n - i, n) for i in a])
            if neg == a
                push!(cosetreplist, (a[1], a[1]))
                push!(cosetpairlist, (a, a))
            else
                for b in arr
                    if neg == b
                        push!(cosetreplist, (a[1], b[1]))
                        push!(cosetpairlist, (a, b))
                        break
                    end
                end
            end
        end
    end

    return cosetpairlist, cosetreplist
end

function qcosetpairings(q::Integer, n::Integer)
    arr = allcyclotomiccosets(q, n, false)
    return qcosetpairings(arr, n)
end

function qcosettable(a::Integer, b::Integer, q::Integer)
    for n in a:b
        if n % q != 0
            println("n = $n")
            allcyclotomiccosets(q, n, true)
            println(" ")
        end
    end
end

"""
    dualqcosets(q::Integer, n::Integer, qcosets::Vector{Vector{Int64}})

Return the dual of the `q`-cyclotomic cosets modulo `n` of `qcosets`.
"""
function dualqcosets(q::Integer, n::Integer, qcosets::Vector{Vector{Int64}})
    compcosets = complementqcosets(q, n, qcosets)
    for a in compcosets
        for (i, x) in enumerate(a)
            a[i] = mod(n - a[i], n)
        end
        sort!(a)
    end

    return compcosets
end
