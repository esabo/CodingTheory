# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

mutable struct CyclicCode <: AbstractCyclicCode
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
    weightenum::Union{WeightEnumerator, Missing}
end

mutable struct BCHCode <: AbstractBCHCode
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
    weightenum::Union{WeightEnumerator, Missing}
end

mutable struct ReedSolomonCode <: AbstractReedSolomonCode
    F::FqNmodFiniteField # base field
    E::FqNmodFiniteField # splitting field
    R::FqNmodPolyRing # polynomial ring of generator polynomial
    β::fq_nmod # n-th root of primitive element of splitting field
    n::Integer
    k::Integer
    d::Union{Integer, Missing}
    b::Integer
    δ::Integer # BCH bound
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
    weightenum::Union{WeightEnumerator, Missing}
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

"""
    definingset(nums::Vector{Int64}, q::Integer, n::Integer, flat::Bool=true)

Returns the set of `q`-cyclotomic cosets of the numbers in `nums` modulo
`n`.

If `flat` is set to true, the result will be a single flattened and sorted
array.
"""
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

    flat && return sort!(vcat(arr...))
    return arr
end

function _idempotent(g::fq_nmod_poly, h::fq_nmod_poly, n::Integer)
    # solve 1 = a(x) g(x) + b(x) h(x) for a(x) then e(x) = a(x) g(x) mod x^n - 1
    d, a, b = gcdx(g, h)
    return mod(g * a, gen(parent(g))^n - 1)
end

# MattsonSolomontransform(f, n)
# inverseMattsonSolomontransform

"""
    field(C::AbstractCyclicCode)

Return the base field of the generator matrix as a Nemo object.
"""
field(C::AbstractCyclicCode) = C.F

"""
    splittingfield(C::AbstractCyclicCode)

Return the splitting field of the generator polynomial as a Nemo object.
"""
splittingfield(C::AbstractCyclicCode) = C.E

"""
    polynomialring(C::AbstractCyclicCode)

Return the polynomial ring of the generator polynomial as a Nemo object.
"""
polynomialring(C::AbstractCyclicCode) = C.R

"""
    primitiveroot(C::AbstractCyclicCode)

Return the primitive root of the splitting field as a Nemo object.
"""
primitiveroot(C::AbstractCyclicCode) = C.β

"""
    offset(C::AbstractBCHCode)

Return the offset of the BCH code.
"""
offset(C::AbstractBCHCode) = C.b

"""
    designdistance(C::AbstractBCHCode)

Return the design distance of the BCH code.
"""
designdistance(C::AbstractBCHCode) = C.δ

"""
    mindistlowerbound(C::AbstractCyclicCode)

Return a lower bound on the minimum distance of the code.

At the moment, this is only the BCH bound with the Hartmann-Tzeng Bound
refinement. The minimum distance is returned if known.
"""
mindistlowerbound(C::AbstractCyclicCode) = C.δ

"""
    qcosets(C::AbstractCyclicCode)

Return the q-cyclotomic cosets of the cyclic code.
"""
qcosets(C::AbstractCyclicCode) = C.qcosets

"""
    qcosetsreps(C::AbstractCyclicCode)

Return the set of representatives for the q-cyclotomic cosets of the cyclic code.
"""
qcosetsreps(C::AbstractCyclicCode) = C.qcosetsreps

"""
    definingset(C::AbstractCyclicCode)

Return the defining set of the cyclic code.
"""
definingset(C::AbstractCyclicCode) = C.defset

"""
    zeros(C::AbstractCyclicCode)

Return the zeros of `C`.
"""
zeros(C::AbstractCyclicCode) = [primitiveroot(C)^i for i in definingset(C)]

"""
    nonzeros(C::AbstractCyclicCode)

Return the nonzeros of `C`.
"""
nonzeros(C::AbstractCyclicCode) = [primitiveroot(C)^i for i in setdiff(0:length(C)
    - 1, definingset(C))]

"""
    generatorpolynomial(C::AbstractCyclicCode)

Return the generator polynomial of the cyclic code as a Nemo object.
"""
generatorpolynomial(C::AbstractCyclicCode) = C.g

"""
    paritycheckpolynomial(C::AbstractCyclicCode)

Return the parity-check polynomial of the cyclic code as a Nemo object.
"""
paritycheckpolynomial(C::AbstractCyclicCode) = C.h

"""
    idempotent(C::AbstractCyclicCode)

Return the idempotent (polynomial) of the cyclic code as a Nemo object.
"""
idempotent(C::AbstractCyclicCode) = C.e

"""
    isprimitive(C::AbstractBCHCode)

Return `true` if the BCH code is primitive.
"""
isprimitive(C::AbstractBCHCode) = length(n) == order(basefield(C)) - 1

"""
    isnarrowsense(C::AbstractBCHCode)

Return `true` if the BCH code is narrowsense.
"""
isnarrowsense(C::AbstractBCHCode) = iszero(offset(C)) # should we define this as b = 1 instead?

"""
    isreversible(C::AbstractCyclicCode)

Return `true` if the cyclic code is reversible.
"""
isreversible(C::AbstractCyclicCode) = [length(C) - i for i in definingset(C)] ⊆ definingset(C)


"""
    isdegenerate(C::AbstractCyclicCode)

Return `true` if the cyclic code is degenerate.

A cyclic code is degenerate if the parity-check polynomial divides `x^r - 1` for
some `r` less than the length of the code.
"""
function isdegenerate(C::AbstractCyclicCode)
    for r in 1:length(C) - 1
        flag, _ = divides(gen(polynomialring(C))^n - 1, paritycheckpolynomial(C))
        flag && return true
    end
    return false
end

# TODO: why doesn't this check C.d like the other shows?
function show(io::IO, C::AbstractCyclicCode)
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
        if !ismissing(C.weightenum)
            println(io, "\nComplete weight enumerator:")
            println(io, "\t", polynomial(C.weightenum))
        end
    end
end

"""
    finddelta(n::Integer, cosets::Vector{Vector{Int64}})

Return the number of consecutive elements of `cosets`, the offset for this, and
a lower bound on the distance of the code defined with length `n` and
cyclotomic cosets `cosets`.

The lower bound is determined by applying the Hartmann-Tzeng Bound refinement to
the BCH bound.
"""
function finddelta(n::Integer, cosets::Vector{Vector{Int64}})
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
    # BCH Bound is thus d ≥ δ

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

"""
    dualdefiningset(defset::Vector{Int64}, n::Integer)

Return the defining set of the dual code of length `n` and defining set `defset`.
"""
dualdefiningset(defset::Vector{Int64}, n::Integer) = sort!([mod(n - i, n) for i in setdiff(0:n - 1, defset)])

"""
    CyclicCode(q::Integer, n::Integer, cosets::Vector{Vector{Int64}})

Return the CyclicCode of length `n` over `GF(q)` with `q`-cyclotomic cosets `cosets`.

This function will auto determine if the constructed code is BCH or Reed-Solomon
and call the appropriate constructor.

# Examples
```julia
julia> q = 2; n = 15; b = 3; δ = 4;
julia> cosets = definingset([i for i = b:(b + δ - 2)], q, n, false);
julia> C = CyclicCode(q, n, cosets)
```
"""
function CyclicCode(q::Integer, n::Integer, cosets::Vector{Vector{Int64}})
    !(q <= 1 || n <= 1) ||error("Invalid parameters past to CyclicCode constructor: q = $q, n = $n.")

    if !isprime(q)
        factors = factor(q)
        length(factors) == 1 || error("There is no finite field of order $q.")
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

    # verify
    flag, htest = divides(gen(R)^n - 1, g)
    flag || error("Incorrect generator polynomial, does not divide x^$n - 1.")
    htest == h || error("Division of x^$n - 1 by the generator polynomial does not yield the constructed parity check polynomial.")
    if size(H) == (n - k, k)
        H = transpose(H)
    end
    !(!iszero(G * transpose(H)) || !iszero(H * transpose(G))) || error("Generator and parity check matrices are not transpose orthogonal.")
    for r in 1:nrows(Gstand)
        iszero(Gstand[r, :] * transpose(H)) || error("Column swap appeared in _standardform.")
    end

    # check e=e^2

    if δ >= 2 && defset == definingset([i for i = b:(b + δ - 2)], q, n, true)
        if deg == 1 && n == q - 1
            return ReedSolomonCode(F, E, R, β, n, k, n - k + 1, b, HT, cosets,
                sort!([arr[1] for arr in cosets]), defset, g, h, e, G, missing,
                H, missing, Gstand, Hstand, missing)
        end

        return BCHCode(F, E, R, β, n, k, missing, b, HT, cosets,
            sort!([arr[1] for arr in cosets]), defset, g, h, e, G, missing, H,
            missing, Gstand, Hstand, missing)
    end

    return CyclicCode(F, E, R, β, n, k, missing, b, HT, cosets,
        sort!([arr[1] for arr in cosets]), defset, g, h, e, G, missing, H, missing,
        Gstand, Hstand, missing)
end

# currently untested - not fully fixed yet
function CyclicCode(q::Integer, n::Integer, g::fq_nmod_poly)
    flag, htest = divides(gen(R)^n - 1, g)
    flag || error("Given polynomial does not divide x^$n - 1.")

    if !isprime(q)
        factors = factor(q)
        length(factors) == 1 || error("There is no finite field of order $q.")
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

    # verify
    h, _, _, _ = _generatorpolynomial(q, n, vcat(comcosets...))
    htest == h || error("Division of x^$n - 1 by the generator polynomial does not yield the constructed parity check polynomial.")
    if size(H) == (n - k, k)
        H = transpose(H)
    end
    !(!iszero(G * transpose(H)) || !iszero(H * transpose(G))) || error("Generator and parity check matrices are not transpose orthogonal.")
    for r in 1:nrows(Gstand)
        iszero(Gstand[r, :] * transpose(H)) || error("Column swap appeared in _standardform.")
    end

    return CyclicCode(F, E, R, β, n, k, missing, b, HT, qcosets,
        sort!([arr[1] for arr in qcosets]), defset, g, h, e, G, H, Gstand, Hstand,
        missing)
end

# self orthogonal cyclic codes are even-like
# does this require them too have even minimum distance?
# self orthogonal code must contain all of its self orthogonal q-cosets and at least one of every q-coset pair
"""
    BCHCode(q::Integer, n::Integer, δ::Integer, b::Integer=0)

Return the BCHCode of length `n` over `GF(q)` with design distance `δ` and offset
`b`.

This function will auto determine if the constructed code is Reed-Solomon
and call the appropriate constructor.

# Examples
```julia
julia> q = 2; n = 15; b = 3; δ = 4;
julia> B = BCHCode(q, n, δ, b)
[15, 5, ≥7; 1]_2 BCH code over splitting field GF(16).
2-Cyclotomic cosets:
        C_1 ∪ C_3 ∪ C_5
Generator polynomial:
        x^10 + x^8 + x^5 + x^4 + x^2 + x + 1
Generator matrix: 5 × 15
        1 1 1 0 1 1 0 0 1 0 1 0 0 0 0
        0 1 1 1 0 1 1 0 0 1 0 1 0 0 0

        0 0 1 1 1 0 1 1 0 0 1 0 1 0 0
        0 0 0 1 1 1 0 1 1 0 0 1 0 1 0
        0 0 0 0 1 1 1 0 1 1 0 0 1 0 1
```
"""
function BCHCode(q::Integer, n::Integer, δ::Integer, b::Integer=0)
    δ >= 2 || error("BCH codes require δ ≥ 2 but the constructor was given δ = $δ.")
    !(q <= 1 || n <= 1) || error("Invalid parameters past to BCHCode constructor: q = $q, n = $n.")

    if !isprime(q)
        factors = factor(q)
        length(factors) == 1 || error("There is no finite field of order $q.")
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

    # verify
    flag, htest = divides(gen(R)^n - 1, g)
    flag || error("Incorrect generator polynomial, does not divide x^$n - 1.")
    htest == h || error("Division of x^$n - 1 by the generator polynomial does not yield the constructed parity check polynomial.")
    if size(H) == (n - k, k)
        H = transpose(H)
    end
    !(!iszero(G * transpose(H)) || !iszero(H * transpose(G))) || error("Generator and parity check matrices are not transpose orthogonal.")
    for r in 1:nrows(Gstand)
        iszero(Gstand[r, :] * transpose(H)) || error("Column swap appeared in _standardform.")
    end

    if deg == 1 && n == q - 1
        return ReedSolomonCode(F, E, R, β, n, k, n - k + 1, b, HT, cosets,
            sort!([arr[1] for arr in cosets]), defset, g, h, e, G, missing, H,
            missing, Gstand, Hstand, missing)
    end

    return BCHCode(F, E, R, β, n, k, missing, b, HT, cosets,
        sort!([arr[1] for arr in cosets]), defset, g, h, e, G, missing, H,
        missing, Gstand, Hstand, missing)
end

"""
    ReedSolomonCode(q::Integer, δ::Integer, b::Integer=0)

Return the ReedSolomonCode over `GF(q)` with distance `d` and offset `b`.

# Examples
```julia
julia> ReedSolomonCode(8, 3, 0)
[7, 5, ≥3; 0]_8 Reed Solomon code.
8-Cyclotomic cosets:
        C_0 ∪ C_1
Generator polynomial:
        x^2 + (α + 1)*x + α
Generator matrix: 5 × 7
        α α + 1 1 0 0 0 0
        0 α α + 1 1 0 0 0
        0 0 α α + 1 1 0 0
        0 0 0 α α + 1 1 0
        0 0 0 0 α α + 1 1

julia> ReedSolomonCode(13, 5, 1)
[12, 8, ≥5; 1]_13 Reed Solomon code.
13-Cyclotomic cosets:
        C_1 ∪ C_2 ∪ C_3 ∪ C_4
Generator polynomial:
        x^4 + 9*x^3 + 7*x^2 + 2*x + 10
Generator matrix: 8 × 12
        10 2 7 9 1 0 0 0 0 0 0 0
        0 10 2 7 9 1 0 0 0 0 0 0
        0 0 10 2 7 9 1 0 0 0 0 0
        0 0 0 10 2 7 9 1 0 0 0 0
        0 0 0 0 10 2 7 9 1 0 0 0
        0 0 0 0 0 10 2 7 9 1 0 0
        0 0 0 0 0 0 10 2 7 9 1 0
        0 0 0 0 0 0 0 10 2 7 9 1
```
"""
function ReedSolomonCode(q::Integer, d::Integer, b::Integer=0)
    d >= 2 || error("Reed Solomon codes require δ ≥ 2 but the constructor was given d = $d.")
    q > 4 || error("Invalid or too small parameters past to ReedSolomonCode constructor: q = $q.")

    # n = q - 1
    # if ord(n, q) != 1
    #     error("Reed Solomon codes require n = q - 1.")
    # end

    if !isprime(q)
        factors = factor(q)
        length(factors) == 1 || error("There is no finite field of order $q.")
        (p, t), = factors
    else
        p = q
        t = 1
    end

    F, α = FiniteField(p, t, "α")
    R, _ = PolynomialRing(F, "x")

    n = q - 1
    cosets = definingset([i for i = b:(b + d - 2)], q, n, false)
    defset = sort!(vcat(cosets...))
    k = n - length(defset)
    comcosets = complementqcosets(q, n, cosets)
    g = _generatorpolynomial(R, α, defset)
    h = _generatorpolynomial(R, α, vcat(comcosets...))
    e = _idempotent(g, h, n)
    G = _generatormatrix(F, n, k, g)
    H = _generatormatrix(F, n, n - k, reverse(h))
    Gstand, Hstand = _standardform(G)

    # verify
    flag, htest = divides(gen(R)^n - 1, g)
    flag || error("Incorrect generator polynomial, does not divide x^$n - 1.")
    htest == h || error("Division of x^$n - 1 by the generator polynomial does not yield the constructed parity check polynomial.")
    if size(H) != (n - k, k)
        H = transpose(H)
    end
    !(!iszero(G * H) || !iszero(transpose(H) * transpose(G))) || error("Generator and parity check matrices are not transpose orthogonal.")
    for r in 1:nrows(Gstand)
        iszero(Gstand[r, :] * H) || error("Column swap appeared in _standardform.")
    end

    return ReedSolomonCode(F, F, R, α, n, k, n - k + 1, b, d, cosets,
        sort!([arr[1] for arr in cosets]), defset, g, h, e, G, missing, H,
        missing, Gstand, Hstand, missing)
end

"""
    complement(C::AbstractCyclicCode)

Return the cyclic code whose cyclotomic cosets are the completement of `C`'s.
"""
function complement(C::AbstractCyclicCode)
    D = CyclicCode(Int64(order(field(C))), length(C),
        complementqcosets(Int64(order(field(C))), length(C), qcosets(C)))
    (paritycheckpolynomial(C) != generatorpolynomial(D) || idempotent(D) != (1 -
        idempotent(C))) && error("Error constructing the complement cyclic code.")
    return D
end

# C1 ⊆ C2 iff g_2(x) | g_1(x) iff T_2 ⊆ T_1
"""
    ⊆(C1::AbstractCyclicCode, C2::AbstractCyclicCode)
    ⊂(C1::AbstractCyclicCode, C2::AbstractCyclicCode)
    issubcode(C1::AbstractCyclicCode, C2::AbstractCyclicCode)

Return whether or not `C1` is a subcode of `C2`.
"""
⊆(C1::AbstractCyclicCode, C2::AbstractCyclicCode) = C2.defset ⊆ C1.defset
⊂(C1::AbstractCyclicCode, C2::AbstractCyclicCode) = C1 ⊆ C2
issubcode(C1::AbstractCyclicCode, C2::AbstractCyclicCode) = C1 ⊆ C2

"""
    ==(C1::AbstractCyclicCode, C2::AbstractCyclicCode)

Return whether or not `C1` and `C2` have the same fields, lengths, and defining sets.
"""
function ==(C1::AbstractCyclicCode, C2::AbstractCyclicCode)
    # should also check primitive root but so far the user is not given a choice here
    return field(C1) == field(C2) && length(C1) == length(C2) && definingset(C1) == definingset(C2)
end

"""
    dual(C::AbstractCyclicCode)

Return the dual of the cyclic code `C`.

Unlike with `LinearCode`, everything is recomputed here so the proper
polynomials and cyclotomic cosets are stored.
"""
function dual(C::AbstractCyclicCode)
    # one is even-like and the other is odd-like
    return CyclicCode(Int64(order(field(C))), length(C),
        dualqcosets(Int64(order(field(C))), length(C), qcosets(C)))
end

# this checks def set, need to rewrite == for linear first
"""
    isselfdual(C::AbstractCyclicCode)

Return whether or not `C == dual(C)`.
"""
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

"""
    ∩(C1::AbstractCyclicCode, C2::AbstractCyclicCode)

Return the intersection code of `C1` and `C2`.
"""
function ∩(C1::AbstractCyclicCode, C2::AbstractCyclicCode)
    # has generator polynomial lcm(g_1(x), g_2(x))
    # has generator idempotent e_1(x) e_2(x)

    if field(C1) == field(C2) && length(C1) == length(C2)
        return CyclicCode(Int64(order(field(C1))), length(C1),
            definingset(definingset(C1) ∪ definingset(C2), Int64(order(field(C1))),
            length(C1), false))
    else
        error("Cannot intersect two codes over different base fields or lengths.")
    end
end

"""
    +(C1::AbstractCyclicCode, C2::AbstractCyclicCode)

Return the addition code of `C1` and `C2`.
"""
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

# "Schur products of linear codes: a study of parameters"
# Diego Mirandola
# """
#     entrywiseproductcode(C::AbstractCyclicCode)
#     *(C::AbstractCyclicCode)
#     Schurproductcode(C::AbstractCyclicCode)
#     Hadamardproductcode(C::AbstractCyclicCode)
#     componentwiseproductcode(C::AbstractCyclicCode)
#
# Return the entrywise product of `C` with itself, which is also a cyclic code.
#
# Note that this is known to often be the full ambient space.
# """
# function entrywiseproductcode(C::AbstractCyclicCode)
#     # generator polynomial is gcd(g*g, g*g*x, g*g*x^{k - 1})
#     R = parent(g)
#     g = generatorpolynomial(C)
#     coefsg = collect(coefficients(g))
#     n = length(coefsg)
#     cur = R([coefsg[i] * coefsg[i] for i in 1:n])
#     for i in 1:dimension(C) - 1
#         coefsgx = collect(coefficents(g * x^i))
#         cur = gcd(cur, R([coefsg[i] * coefsgx[i] for i in 1:n]))
#     end
#     return CyclicCode(cur)
# end
# *(C::AbstractCyclicCode) = entrywiseproductcode(C)
# Schurproductcode(C::AbstractCyclicCode) = entrywiseproductcode(C)
# Hadamardproductcode(C::AbstractCyclicCode) = entrywiseproductcode(C)
# componentwiseproductcode(C::AbstractCyclicCode) = entrywiseproductcode(C)
