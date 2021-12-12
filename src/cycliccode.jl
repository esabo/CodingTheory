# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

import Base: ==, ∩, +

include("cyclotomic.jl")
include("linearcode.jl")

abstract type AbstractCyclicCode <: AbstractLinearCode end
abstract type AbstractBCHCode <: AbstractCyclicCode end
abstract type AbstractReedSolomonCode <: AbstractBCHCode end

# need to write == functions
struct CyclicCode <: AbstractCyclicCode
    F::FqNmodFiniteField # base field
    E::FqNmodFiniteField # splitting field
    R::FqNmodPolyRing # polynomial ring of generator polynomial
    β::fq_nmod # n-th root of primitive element of splitting field
    n::Integer
    k::Integer
    d::Union{Integer, Missing}
    b::Integer
    δ::Integer
    qcosets::Vector{Vector{Int64}}
    qcosetsreps::Vector{Int64}
    defset::Vector{Int64}
    g::fq_nmod_poly
    h::fq_nmod_poly
    e::fq_nmod_poly
    G::fq_nmod_mat
    Gorig::Union{fq_nmod_mat, Missing}
    H::fq_nmod_mat
    Horig::Union{fq_nmod_mat, Missing}
    Gstand::fq_nmod_mat
    Hstand::fq_nmod_mat
end

struct BCHCode <: AbstractBCHCode
    F::FqNmodFiniteField # base field
    E::FqNmodFiniteField # splitting field
    R::FqNmodPolyRing # polynomial ring of generator polynomial
    β::fq_nmod # n-th root of primitive element of splitting field
    n::Integer
    k::Integer
    d::Union{Integer, Missing}
    b::Integer
    δ::Integer
    qcosets::Vector{Vector{Int64}}
    qcosetsreps::Vector{Int64}
    defset::Vector{Int64}
    g::fq_nmod_poly
    h::fq_nmod_poly
    e::fq_nmod_poly
    G::fq_nmod_mat
    Gorig::Union{fq_nmod_mat, Missing}
    H::fq_nmod_mat
    Horig::Union{fq_nmod_mat, Missing}
    Gstand::fq_nmod_mat
    Hstand::fq_nmod_mat
end

struct ReedSolomonCode <: AbstractReedSolomonCode
    F::FqNmodFiniteField # base field
    E::FqNmodFiniteField # splitting field
    R::FqNmodPolyRing # polynomial ring of generator polynomial
    β::fq_nmod # n-th root of primitive element of splitting field
    n::Integer
    k::Integer
    d::Union{Integer, Missing}
    b::Integer
    δ::Integer
    qcosets::Vector{Vector{Int64}}
    qcosetsreps::Vector{Int64}
    defset::Vector{Int64}
    g::fq_nmod_poly
    h::fq_nmod_poly
    e::fq_nmod_poly
    G::fq_nmod_mat
    Gorig::Union{fq_nmod_mat, Missing}
    H::fq_nmod_mat
    Horig::Union{fq_nmod_mat, Missing}
    Gstand::fq_nmod_mat
    Hstand::fq_nmod_mat
end

function _generatorpolynomial(R::FqNmodPolyRing, β::fq_nmod, Z::Vector{Int64})
    # from_roots(R, [β^i for i in Z]) - R has wrong type for this
    g = one(R)
    for i in Z
        g *= (gen(R) - β^i)
    end

    return g
end
_generatorpolynomial(R::FqNmodPolyRing, β::fq_nmod, qcosets::Vector{Vector{Int64}}) = _generatorpolynomial(R, β, vcat(qcosets...))

function _generatormatrix(F::FqNmodFiniteField, n::Integer, k::Integer, g::fq_nmod_poly)
    # if g = x^10 + α^2*x^9 + x^8 + α*x^7 + x^3 + α^2*x^2 + x + α
    # g.coeffs = [α  1  α^2  1  0  0  0  α  1  α^2  1]
    coeffs = collect(coefficients(g))
    len = length(coeffs)
    k + len - 1 <= n || error("Too many coefficients for $k shifts in _generatormatrix.")

    G = zero_matrix(F, k, n)
    for i in 1:k
        G[i, i:i + len - 1] = coeffs
    end
    return G
end

function definingset(nums::Vector{Int64}, q::Integer, n::Integer, flat::Bool=true)
    arr = Vector{Vector{Int64}}()
    arrflat = Vector{Int64}()
    for x in nums
        Cx = cyclotomiccoset(x, q, n)
        if Cx[1] ∉ arrflat
            arrflat = [arrflat; Cx]
            push!(arr, Cx)
        end
    end

    !flat || return sort!(vcat(arr...))
    return arr
end

function _idempotent(g::fq_nmod_poly, h::fq_nmod_poly, n::Integer)
    # solve 1 = a(x) g(x) + b(x) h(x) for a(x) then e(x) = a(x) g(x) mod x^n - 1
    d, a, b = gcdx(g, h)
    return mod(g * a, gen(parent(g))^n - 1)
end

# MattsonSolomontransform(f, n)
# inverseMattsonSolomontransform


basefield(C::AbstractCyclicCode) = C.F
splittingfield(C::AbstractCyclicCode) = C.E
polynomialring(C::AbstractCyclicCode) = C.R
primitiveroot(C::AbstractCyclicCode) = C.β
offset(C::AbstractCyclicCode) = C.b
designdistance(C::AbstractCyclicCode) = C.δ
qcosets(C::AbstractCyclicCode) = C.qcosets
qcosetsreps(C::AbstractCyclicCode) = C.qcosetsreps
definingset(C::AbstractCyclicCode) = C.defset
generatorpolynomial(C::AbstractCyclicCode) = C.g
paritycheckpolynomial(C::AbstractCyclicCode) = C.h
idempotent(C::AbstractCyclicCode) = C.e
isprimitive(C::AbstractCyclicCode) = length(n) == order(basefield(C)) - 1
isnarrowsense(C::AbstractCyclicCode) = iszero(offset(C)) # should we define this as b = 1 instead?
isreversible(C::AbstractCyclicCode) = [length(C) - i for i in defset] ⊆ defset

function Base.show(io::IO, C::AbstractCyclicCode)
    if get(io, :compact, false)
        # to use type "show(IOContext(stdout, :compact=>true), C)" instead
        if typeof(C) <: ReedSolomonCode
            println(io, "[$(length(C)), $(dimension(C)), ≥$(designdistance(C)); $(offset(C))]_$(order(field(C))) Reed Solomon code.")
        elseif typeof(C) <: BCHCode
            println(io, "[$(length(C)), $(dimension(C)), ≥$(designdistance(C)); $(offset(C))]_$(order(field(C))) BCH code over splitting field GF($(order(splittingfield(C)))).")
        else
            println(io, "[$(length(C)), $(dimension(C)), ≥$(designdistance(C)); $(offset(C))]_$(order(field(C))) cyclic code over splitting field GF($(order(splittingfield(C)))).")
        end
    else
        if typeof(C) <: ReedSolomonCode
            println(io, "[$(length(C)), $(dimension(C)), ≥$(designdistance(C)); $(offset(C))]_$(order(field(C))) Reed Solomon code.")
        elseif typeof(C) <: BCHCode
            println(io, "[$(length(C)), $(dimension(C)), ≥$(designdistance(C)); $(offset(C))]_$(order(field(C))) BCH code over splitting field GF($(order(splittingfield(C)))).")
        else
            println(io, "[$(length(C)), $(dimension(C)), ≥$(designdistance(C)); $(offset(C))]_$(order(field(C))) cyclic code over splitting field GF($(order(splittingfield(C)))).")
        end
        println(io, "$(order(field(C)))-Cyclotomic cosets: ")
        for (i, x) in enumerate(qcosetsreps(C))
            if i == 1
                print(io, "\tC_$x ∪ ")
            elseif i == 1 && i == length(qcosetsreps(C))
                println(io, "\tC_$x")
            elseif i != length(qcosetsreps(C))
                print(io, "C_$x ∪ ")
            else
                println(io, "C_$x")
            end
        end
        println(io, "Generator polynomial:")
        println(io, "\t", generatorpolynomial(C))
        println(io, "Generator matrix: $(dimension(C)) × $(length(C))")
        for i in 1:dimension(C)
            print(io, "\t")
            for j in 1:length(C)
                if j != length(C)
                    print(io, "$(C.G[i, j]) ")
                elseif j == length(C) && i != dimension(C)
                    println(io, "$(C.G[i, j])")
                else
                    print(io, "$(C.G[i, j])")
                end
            end
        end
    end
end

function finddelta(n::Integer, cosets::Vector{Vector{Int64}})
    # println(cosets)
    defset = sort!(vcat(cosets...))
    runs = Vector{Vector{Int64}}()
    for x in defset
        useddefset = Vector{Int64}()
        reps = Vector{Int64}()
        cosetnum = 0
        for i in 1:length(cosets)
            if x ∈ cosets[i]
                cosetnum = i
                append!(useddefset, cosets[i])
                append!(reps, x)
                break
            end
        end

        y = x + 1
        while y ∈ defset
            if y ∈ useddefset
                append!(reps, y)
            else
                cosetnum = 0
                for i in 1:length(cosets)
                    if y ∈ cosets[i]
                        cosetnum = i
                        append!(useddefset, cosets[i])
                        append!(reps, y)
                        break
                    end
                end
            end
            y += 1
        end
        push!(runs, reps)
    end

    runlens = [length(i) for i in runs]
    (consec, ind) = findmax(runlens)
    # there are δ - 1 consecutive numbers for designed distance δ
    δ = consec + 1
    # start of run
    offset = runs[ind][1]
    # BCH Bound is thus d≥δ

    # moving to Hartmann-Tzeng Bound refinement
    currbound = δ
    if consec > 1
        for A in runs
            if length(A) == consec
                for b in 1:(n - 1)
                    if gcd(b, n) ≤ δ
                        for s in 0:(δ - 2)
                            B = [mod(j * b, n) for j in 0:s]
                            AB = [x + y for x in A for y in B]
                            if AB ⊆ defset
                                if currbound < δ + s
                                    currbound = δ + s
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return δ, offset, currbound
end

function largestconsecrun(arr::Vector{Int64})
    n = length(arr)
    maxlen = 1
    for i = 1:n
        mn = arr[i]
        mx = arr[i]

        for j = (i + 1):n
            mn = min(mn, arr[j])
            mx = max(mx, arr[j])

            if (mx - mn) == (j - i)
                maxlen = max(maxlen, mx - mn + 1)
            end
        end
    end

    return maxlen
end

function dualdefiningset(arr::Vector{Int64}, n::Integer)
    full = [i for i in 0:(n - 1)]
    temp = Vector{Int64}()
    for i in full
        if i ∉ arr
            append!(temp, i)
        end
    end

    return sort!([mod(n - i, n) for i in temp])
end

function CyclicCode(q::Integer, n::Integer, cosets::Vector{Vector{Int64}}, verify::Bool=true)
    !(q <= 1 || n <= 1) ||error("Invalid parameters past to CyclicCode constructor: q = $q, n = $n.")

    if !isprime(q)
        factors = factor(q)
        if length(factors) != 1
            error("There is no finite field of order $(prod(factors)).")
        end
        (p, t), = factors
    else
        p = q
        t = 1
    end

    F, _ = FiniteField(p, t, "α") # changed to keep ReedSolomonCodes printing α's'
    deg = ord(n, q)
    E, α = FiniteField(p, t * deg, "α")
    R, _ = PolynomialRing(E, "x")
    β = α^(div(q^deg - 1, n))

    defset = sort!(vcat(cosets...))
    k = n - length(defset)
    comcosets = complementqcosets(q, n, cosets)
    g = _generatorpolynomial(R, β, defset)
    h = _generatorpolynomial(R, β, vcat(comcosets...))
    e = _idempotent(g, h, n)
    G = _generatormatrix(F, n, k, g)
    H = _generatormatrix(F, n, n - k, reverse(h))
    Gstand, Hstand = _standardform(G)
    δ, b, HT = finddelta(n, cosets)

    if verify
        flag, htest = divides(gen(R)^n - 1, g)
        flag || error("Incorrect generator polynomial, does not divide x^$n - 1.")
        htest == h || error("Division of x^$n - 1 by the generator polynomial does not yield the constructed parity check polynomial.")
        if size(H) == (n - k, k)
            H = deepcopy(H')
        end
        !(!iszero(G * H') || !iszero(H * G')) || error("Generator and parity check matrices are not transpose orthogonal.")
        for r in 1:size(Gstand, 1)
            iszero(Gstand[r, :] * H') || error("Column swap appeared in _standardform.")
        end

        # check e=e^2
    end

    if δ >= 2 && defset == definingset([i for i = b:(b + δ - 2)], q, n, true)
        if deg == 1 && n == q - 1
            return ReedSolomonCode(F, E, R, β, n, k, n - k + 1, b, HT, cosets, sort!([arr[1] for arr in cosets]),
                defset, g, h, e, G, missing, H, missing, Gstand, Hstand)
        end

        return BCHCode(F, E, R, β, n, k, missing, b, HT, cosets, sort!([arr[1] for arr in cosets]),
            defset, g, h, e, G, missing, H, missing, Gstand, Hstand)
    end

    return CyclicCode(F, E, R, β, n, k, missing, b, HT, cosets, sort!([arr[1] for arr in cosets]),
        defset, g, h, e, G, missing, H, missing, Gstand, Hstand)
end

# currently untested - not fully fixed yet
function CyclicCode(q::Integer, n::Integer, g::fq_nmod_poly, verify::Bool=true)
    flag, htest = divides(gen(R)^n - 1, g)
    flag || error("Given polynomial does not divide x^$n - 1.")

    if !isprime(q)
        factors = factor(q)
        if length(factors) != 1
            error("There is no finite field of order $(prod(factors)).")
        end
        (q, n), = factors
    end

    R = parent(g)
    F = base_ring(R)
    α = gen(F)
    t = ord(n, q)
    β = α^(div(q^t - 1, n))

    dic = Dict{fq_nmod, Int64}()
    for i in 1:n
        dic[β^i] = i
    end
    rts = roots(g)
    defset = sort!([dic[rt] for rt in rts])
    qcosets = definingset(defset, q, n, false)
    k = n - length(defset)
    comcosets = complementqcosets(q, n, qcosets)
    e = _idempotent(g, h)
    G = _generatormatrix(n, k, g)
    H = _generatormatrix(n, n - k, reverse(htest))
    Gstand, Hstand = _standardform(G)
    _, _, HT = finddelta(n, qcosets)

    if verify
        h, _, _, _ = _generatorpolynomial(q, n, vcat(comcosets...))
        htest == h || error("Division of x^$n - 1 by the generator polynomial does not yield the constructed parity check polynomial.")
        if size(H) == (n - k, k)
            H = deepcopy(H')
        end
        !(!iszero(G * H') || !iszero(H * G')) || error("Generator and parity check matrices are not transpose orthogonal.")
        for r in 1:size(Gstand, 1)
            iszero(Gstand[r, :] * H') || error("Column swap appeared in _standardform.")
        end
    end

    return CyclicCode(F, E, R, β, n, k, missing, b, HT, qcosets, sort!([arr[1] for arr in qcosets]),
        defset, g, h, e, G, H, Gstand, Hstand)
end

# self orthogonal cyclic codes are even-like
# does this require them too have even minimum distance?
# self orthogonal code must contain all of its self orthogonal q-cosets and at least one of every q-coset pair
function BCHCode(q::Integer, n::Integer, δ::Integer, b::Integer=0, verify::Bool=true)
    δ >= 2 || error("BCH codes require δ ≥ 2 but the constructor was given δ = $δ.")
    !(q <= 1 || n <= 1) || error("Invalid parameters past to BCHCode constructor: q = $q, n = $n.")

    if !isprime(q)
        factors = factor(q)
        if length(factors) != 1
            error("There is no finite field of order $(prod(factors)).")
        end
        (p, t), = factors
    else
        p = q
        t = 1
    end

    F, _ = FiniteField(p, t, "α") # changed to keep ReedSolomonCodes printing α's'
    deg = ord(n, q)
    E, α = FiniteField(p, t * deg, "α")
    R, _ = PolynomialRing(E, "x")
    β = α^(div(q^deg - 1, n))

    cosets = definingset([i for i = b:(b + δ - 2)], q, n, false)
    defset = sort!(vcat(cosets...))
    k = n - length(defset)
    comcosets = complementqcosets(q, n, cosets)
    g = _generatorpolynomial(R, β, defset)
    h = _generatorpolynomial(R, β, vcat(comcosets...))
    e = _idempotent(g, h, n)
    G = _generatormatrix(F, n, k, g)
    H = _generatormatrix(F, n, n - k, reverse(h))
    Gstand, Hstand = _standardform(G)
    δ, b, HT = finddelta(n, cosets)

    if verify
        flag, htest = divides(gen(R)^n - 1, g)
        flag || error("Incorrect generator polynomial, does not divide x^$n - 1.")
        htest == h || error("Division of x^$n - 1 by the generator polynomial does not yield the constructed parity check polynomial.")
        if size(H) == (n - k, k)
            H = deepcopy(H')
        end
        !(!iszero(G * H') || !iszero(H * G')) || error("Generator and parity check matrices are not transpose orthogonal.")
        for r in 1:size(Gstand, 1)
            iszero(Gstand[r, :] * H') || error("Column swap appeared in _standardform.")
        end
    end

    if deg == 1 && n == q - 1
        return ReedSolomonCode(F, E, R, β, n, k, n - k + 1, b, HT, cosets, sort!([arr[1] for arr in cosets]),
            defset, g, h, e, G, missing, H, missing, Gstand, Hstand)
    end

    return BCHCode(F, E, R, β, n, k, n - k + 1, b, HT, cosets, sort!([arr[1] for arr in cosets]),
        defset, g, h, e, G, missing, H, missing, Gstand, Hstand)
end

function ReedSolomonCode(q::Integer, δ::Integer, b::Integer=0, verify::Bool=true)
    δ >= 2 || error("Reed Solomon codes require δ ≥ 2 but the constructor was given δ = $δ.")
    q > 4 || error("Invalid or too small parameters past to ReedSolomonCode constructor: q = $q.")

    # n = q - 1
    # if ord(n, q) != 1
    #     error("Reed Solomon codes require n = q - 1.")
    # end

    if !isprime(q)
        factors = factor(q)
        if length(factors) != 1
            error("There is no finite field of order $(prod(factors)).")
        end
        (p, t), = factors
    else
        p = q
        t = 1
    end

    F, α = FiniteField(p, t, "α") # changed to keep ReedSolomonCodes printing α's'
    R, _ = PolynomialRing(F, "x")

    n = q - 1
    cosets = definingset([i for i = b:(b + δ - 2)], q, n, false)
    defset = sort!(vcat(cosets...))
    k = n - length(defset)
    comcosets = complementqcosets(q, n, cosets)
    g = _generatorpolynomial(R, α, defset)
    h = _generatorpolynomial(R, α, vcat(comcosets...))
    e = _idempotent(g, h, n)
    G = _generatormatrix(F, n, k, g)
    H = _generatormatrix(F, n, n - k, reverse(h))
    Gstand, Hstand = _standardform(G)
    δ, b, HT = finddelta(n, cosets)

    if verify
        flag, htest = divides(gen(R)^n - 1, g)
        flag || error("Incorrect generator polynomial, does not divide x^$n - 1.")
        htest == h || error("Division of x^$n - 1 by the generator polynomial does not yield the constructed parity check polynomial.")
        if size(H) == (n - k, k)
            H = deepcopy(H')
        end
        !(!iszero(G * H') || !iszero(H * G')) || error("Generator and parity check matrices are not transpose orthogonal.")
        for r in 1:size(Gstand, 1)
            iszero(Gstand[r, :] * H) || error("Column swap appeared in _standardform.")
        end
    end

    return ReedSolomonCode(F, F, R, α, n, k, n - k + 1, b, HT, cosets, sort!([arr[1] for arr in cosets]),
        defset, g, h, e, G, missing, H, missing, Gstand, Hstand)
end

function complement(C::AbstractCyclicCode, verify::Bool=true)
    D = CyclicCode(Int64(order(field(C))), length(C),
        complementqcosets(Int64(order(field(C))), length(C), qcosets(C)))
    if verify
        if paritycheckpolynomial(C) != generatorpolynomial(D) || idempotent(D) != (1 - idempotent(C))
            error("Error constructing the complement cyclic code.")
        end
    end

    return D
end

# C1 ⊆ C2 iff g_2(x) | g_1(x) iff T_2 ⊆ T_1
⊆(C1::AbstractCyclicCode, C2::AbstractCyclicCode) = C2.defset ⊆ C1.defset
issubcode(C1::AbstractCyclicCode, C2::AbstractCyclicCode) = C1 ⊆ C2

function ==(C1::AbstractCyclicCode, C2::AbstractCyclicCode)
    # should also check primitive root but so far the user is not given a choice here
    return field(C1) == field(C2) && length(C1) == length(C2) && definingset(C1) == definingset(C2)
end

function dual(C::AbstractCyclicCode)
    # one is even-like and the other is odd-like
    return CyclicCode(Int64(order(field(C))), length(C),
        dualqcosets(Int64(order(field(C))), length(C), qcosets(C)))
end

# this checks def set, need to rewrite == for linear first
isselfdual(C::AbstractCyclicCode) = C == dual(C)

# don't think this is necessary in order to invoke the ⊆ for CyclicCode
# function isselforthogonal(C::AbstractCyclicCode)
#     # A code is self-orthogonal if it is a subcode of its dual.
#     return C ⊆ dual(C)
# end

# function μa(C::CyclicCode)
#     # check gcd(a, n) = 1
#     # technically changes g(x) and e(x) but the q-cosets are the same?
# end

function ∩(C1::AbstractCyclicCode, C2::AbstractCyclicCode)
    # has generator polynomial lcm(g_1(x), g_2(x))
    # has generator idempotent e_1(x) e_2(x)

    if field(C1) == field(C2) && length(C1) == length(C2)
        return CyclicCode(Int64(order(field(C1))), length(C1),
            definingset(definingset(C1) ∪ definingset(C2), Int64(order(field(C1))), length(C1), false))
    else
        error("Cannot intersect two codes over different base fields or lengths.")
    end
end

function +(C1::AbstractCyclicCode, C2::AbstractCyclicCode)
    # has generator polynomial gcd(g_1(x), g_2(x))
    # has generator idempotent e_1(x) + e_2(x) - e_1(x) e_2(x)
    if field(C1) == field(C2) && length(C1) == length(C2)
        defset = definingset(C1) ∩ definingset(C2)
        if length(defset) != 0
            return CyclicCode(Int64(order(field(C1))), length(C1),
                definingset(defset, Int64(order(field(C1))), length(C1), false))
        else
            error("Addition of codes has empty defining set.")
        end
    else
        error("Cannot add two codes over different base fields or lengths.")
    end
end
