# Copyright (c) 2022, 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    LDPCCode(H::fq_nmod_mat)

Return the LDPC code defined by the parity-check matrix `H`.
"""
function LDPCCode(H::CTMatrixTypes)
    # TODO: remove empties
    nnz, den = _density(H)
    # den <= 0.01 || (@warn "LDPC codes typically expect a density of less than 1%.";)

    cols, rows = _degreedistribution(H)
    isreg = true
    c1 = cols[1]
    for i in 2:length(cols)
        c1 == cols[i] || (isreg = false; break;)
    end
    if isreg
        r1 = rows[1]
        for i in 2:length(rows)
            r1 == rows[i] || (isreg = false; break;)
        end
    end
    c, r = maximum(cols), maximum(rows)

    R, x = PolynomialRing(Nemo.QQ, :x)
    colpoly = R(0)
    for i in cols
        colpoly += i * x^(i - 1)
    end
    colpoly = divexact(colpoly, nnz)
    rowpoly = R(0)
    for i in rows
        rowpoly += i * x^(i - 1)
    end
    rowpoly = divexact(rowpoly, nnz)

    C = LinearCode(H, true)
    return LDPCCode(base_ring(H), C.n, C.k, C.d, C.lbound, C.ubound, H, nnz,
        cols, rows, c, r, maximum([c, r]), den, isreg, missing, colpoly,
        rowpoly, missing)
end

"""
    LDPCCode(C::AbstractLinearCode)

Return the LDPC code given by the parity-check matrix of `C`.
"""
LDPCCode(C::AbstractLinearCode) = LDPCCode(paritycheckmatrix(C))

"""
    regularLDPCCode(q::Int, n::Int, l::Int, r::Int [; seed=nothing])

Return a random regular LDPC code over `GF(q)` of length `n` with column degree `l`
and row degree `r`.

If a seed is given, i.e. `regulardLDPCCode(4, 1200, 3, 6, seed=123)`, the
results are reproducible.
"""
function regularLDPCCode(q::Int, n::Int, l::Int, r::Int; seed::Union{Nothing, Int}=nothing)
    Random.seed!(seed)
    m = divexact(n * l, r)
    F = if isprime(q)
        GF(q)
    else
        factors = Nemo.factor(q)
        length(factors) == 1 || throw(DomainError("There is no finite field of order $q"))
        (p, t), = factors
        GF(p, t, :α)
    end
    elems = collect(F)[2:end]
    H = zero_matrix(F, m, n);
    colsums = zeros(Int, n)
    for i in axes(H, 1)
        ind = reduce(vcat, shuffle(filter(k -> colsums[k] == s, 1:n)) for s in 0:l-1)[1:r]
        for j in ind
            H[i, j] = rand(elems)
        end
        colsums[ind] .+= 1
    end
    @assert all(count(.!iszero.(H[:, j])) == l for j in axes(H, 2))
    @assert all(count(.!iszero.(H[i, :])) == r for i in axes(H, 1))

    R, x = PolynomialRing(Nemo.QQ, :x)
    C = LinearCode(H, true)
    return LDPCCode(C.F, C.n, C.k, C.d, C.lbound, C.ubound, H, n * l, l * ones(Int, n),
        r * ones(Int, m), l, r, max(l, r), r / n, true, missing, (1 // l) * x^l,
        (1 // r) * x^r, missing)
end
regularLDPCCode(n::Int, l::Int, r::Int; seed::Union{Nothing, Int} = nothing) =
    regularLDPCCode(2, n, l, r, seed = seed)

#############################
      # getter functions
#############################

"""
    variabledegreedistribution(C::AbstractLDPCCode)

Return the variable node degree distribution of `C`.
"""
variabledegreedistribution(C::AbstractLDPCCode) = C.vardegs

"""
    checkdegreedistribution(C::AbstractLDPCCode)

Return the check node degree distribution of `C`.
"""
checkdegreedistribution(C::AbstractLDPCCode) = C.checkdegs

"""
    degreedistributions(C::AbstractLDPCCode)

Return the variable and check node degree distributions of `C`.
"""
degreedistributions(C::AbstractLDPCCode) = C.vardegs, C.checkdegs

"""
    columnbound(C::AbstractLDPCCode)

Return the column bound `c` of the `(c, r)`-LDPC code `C`.
"""
columnbound(C::AbstractLDPCCode) = C.colbound

"""
    rowbound(C::AbstractLDPCCode)

Return the row bound `r` of the `(c, r)`-LDPC code `C`.
"""
rowbound(C::AbstractLDPCCode) = C.rowbound

"""
    columnrowbounds(C::AbstractLDPCCode)

Return the column and row bounds `c, r` of the `(c, r)`-LDPC code `C`.
"""
columnrowbounds(C::AbstractLDPCCode) = C.colbound, C.rowbound

"""
    limited(C::AbstractLDPCCode)

Return the maximum of the row and column bounds for `C`.
"""
limited(C::AbstractLDPCCode) = C.limited

"""
    density(C::AbstractLDPCCode)

Return the density of the parity-check matrix of `C`.
"""
density(C::AbstractLDPCCode) = C.density

"""
    isregular(C::AbstractLDPCCode)

Return `true` if the `C` is a regular LDPC code.

# Notes
* An LDPC is regular if all the column degrees and equal and all the row degrees are equal.
"""
isregular(C::AbstractLDPCCode) = C.isreg

"""
    variabledegreepolynomial(C::AbstractLDPCCode)

Return the variable degree polynomial of `C`.
"""
variabledegreepolynomial(C::AbstractLDPCCode) = C.λ

"""
    checkdegreepolynomial(C::AbstractLDPCCode)

Return the check degree polynomial of `C`.
"""
checkdegreepolynomial(C::AbstractLDPCCode) = C.ρ

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

function _degreedistribution(H::CTMatrixTypes)
    nr, nc = size(H)
    cols = zeros(Int, 1, nc)
    @inbounds @views @simd for i in 1:nc
        cols[i] = wt(H[:,  i])
    end
    rows = zeros(Int, 1, nr)
    @inbounds @views @simd for i in 1:nr
        rows[i] = wt(H[i,  :])
    end
    return vec(cols), vec(rows)
end

function _density(H::CTMatrixTypes)
    count = 0
    nr, nc = size(H)
    for c in 1:nc
        for r in 1:nr
            !iszero(H[r, c]) && (count += 1;)
        end
    end
    return count, count / (nr * nc)
end

# TODO: make uniform with others
function show(io::IO, C::AbstractLDPCCode)
    if ismissing(C.d)
        if C.isreg
            println(io, "[$(C.n), $(C.k)]_$(order(C.F)) regular ($(C.colbound), $(C.rowbound))-LDPC code with density $(C.density).")
        else
            println(io, "[$(C.n), $(C.k)]_$(order(C.F)) irregular $(C.limited)-limited LDPC code with density $(C.density).")
        end
    else
        if C.isreg
            println(io, "[$(C.n), $(C.k), $(C.d)]_$(order(C.F)) regular ($(C.colbound), $(C.rowbound))-LDPC code with density $(C.density).")
        else
            println(io, "[$(C.n), $(C.k), $(C.d)]_$(order(C.F)) irregular $(C.limited)-limited LDPC code with density $(C.density).")
        end
    end
    if get(io, :compact, true)
        println(io, "\nVariable degree polynomial:")
        println(io, "\t", C.λ)
        println(io, "Check degree polynomial:")
        println(io, "\t", C.ρ)
        if C.n <= 30
            # was using Horig here, which is probably what I want
            nr, nc = size(C.H)
            println(io, "Parity-check matrix: $nr × $nc")
            for i in 1:nr
                print(io, "\t")
                for j in 1:nc
                    if j != nc
                        print(io, "$(C.H[i, j]) ")
                    elseif j == nc && i != nr
                        println(io, "$(C.H[i, j])")
                    else
                        print(io, "$(C.H[i, j])")
                    end
                end
            end
        # if !ismissing(C.weightenum)
        #     println(io, "\nComplete weight enumerator:")
        #     println(io, "\t", C.weightenum.polynomial)
        # end
        end
    end
end

"""
    degreedistributionsplot(C::AbstractLDPCCode)

Return a bar plot of the column and row degree distributions of `C`.
"""
function degreedistributionsplot(C::AbstractLDPCCode)
    cols, rows = degreedistributions(C)

    occurscols = [(i, count(==(i), cols)) for i in unique(cols)]
    colsxdata = [x for (x, _) in occurscols]
    colsydata = [y for (_, y) in occurscols]
    colstitle="Variable Nodes"
    f1 = bar(colsxdata, colsydata, bar_width=1, xticks=colsxdata, yticks=colsydata,
        legend=false, xlabel="Degree", ylabel="Occurrences", title=colstitle)

    occursrows = [(i, count(==(i), rows)) for i in unique(rows)]
    rowsxdata = [x for (x, _) in occursrows]
    rowsydata = [y for (_, y) in occursrows]
    rowstitle="Check Nodes"
    f2 = bar(rowsxdata, rowsydata, bar_width=1, xticks=rowsxdata, yticks=rowsydata,
        legend=false, xlabel="Degree", ylabel="Occurrences", title=rowstitle)
    f = Plots.plot(f1, f2, layout=(1, 2))
    show(f)
    return f
end

"""
    girth(C::LDPCCode)

Return the girth of the Tanner graph of `C`.

An error is thrown if the maximum number of iterations is reached and
``-1`` is returned to represent infinite girth.
"""
function girth(C::LDPCCode, max_iter::Int=100)
    checkadlist, varadlist = _nodeadjacencies(C.H)
    girtharr = zeros(Int, C.n)
    Threads.@threads for vn in 1:C.n
        iter = 0
        notfound = true
        tocheck = [vn]
        while notfound
            iter += 1
            tochecknext = Vector{Int}()
            for v in tocheck
                for cn in varadlist[v]
                    if iter != 1 && vn ∈ checkadlist[cn]
                        notfound = false
                        girtharr[vn] = iter + 1
                        break
                    else
                        for v2 in checkadlist[cn]
                            v != v2 && append!(tochecknext, v2)
                        end
                    end
                end
                iter += 1
                iter > max_iter && error("Hit the maximum number of iterations")
                !notfound && break
                unique!(sort!(tochecknext))
            end
        end
    end
    min = minimum(girtharr)
    iseven(min) || error("Computed girth to be an odd integer")
    min == 0 ? (C.girth = -1;) : (C.girth = min;)
    return C.girth
end

mutable struct _ACEvarnode
    id::Int
    parentid::Int
    lvl::Int
    cumACE::Int
    localACE::Int
    # children::Vector{_ACEchecknode}
end

mutable struct _ACEchecknode
    id::Int
    parentid::Int
    lvl::Int
    cumACE::Int
    # children::Vector{_ACEvarnode}
end

# TODO: girth, degree 1 nodes
"""
    shortestcycleACE(C::LDPCCode, v::Int)
    shortestcycleACE(C::LDPCCode, vs::Vector{Int})
    shortestcycleACE(C::LDPCCode)

Return a cycle of minimum length and minimum ACE in the Tanner graph of `C`
for the vertex `v` or vertices `vs`, in the order (ACEs, cycles). If no vertices
are given, all vertices are computed by default. The cycle `v1 -- c1 -- ... -- 
cn -- vn` are returned in the format `[(v1, c1), (c1, v2), ..., (cn, vn)]`.
"""
function shortestcycleACE(C::LDPCCode, vs::Vector{Int})
    isempty(vs) && throw(ArgumentError("Input variable node list cannot be empty"))
    all(x->1 <= x <= C.n, vs) || throw(DomainError("Variable node indices must be between 1 and length(C)"))

    checkadlist, varadlist = _nodeadjacencies(C.H)
    vsACE = zeros(Int, length(vs))
    cyclesvs = [Vector{Tuple{Int, Int}}() for _ in 1:length(vs)]

    Threads.@threads for i in 1:length(vs)
        # moving this inside allocates more but allows for multi-threading
        checknodes = [_ACEchecknode(i, -1, -1, -1) for i in 1:length(checkadlist)]
        varnodes = [_ACEvarnode(i, -1, -1, -1, length(varadlist[i]) - 2) for i in 1:C.n]

        ACEs = Vector{Int}()
        cycles = Vector{Vector{Tuple{Int, Int}}}()
        notemptied = true
        root = varnodes[vs[i]]
        root.lvl = 0
        root.cumACE = root.localACE
        queue = Deque{Union{_ACEchecknode, _ACEvarnode}}()
        push!(queue, root)
        while length(queue) > 0
            curr = first(queue)
            if isa(curr, _ACEvarnode)
                for cn in varadlist[curr.id]
                    # can't pass messages back to the same node
                    if cn != curr.parentid
                        cnnode = checknodes[cn]
                        if cnnode.lvl != -1
                            # have seen before
                            push!(ACEs, curr.cumACE + cnnode.cumACE - root.localACE)

                            # trace the cycle from curr to root and cnnode to root
                            temp = Vector{Tuple{Int, Int}}()
                            node = cnnode
                            while node.lvl != 0
                                push!(temp, (node.parentid, node.id))
                                if isodd(node.lvl)
                                    node = varnodes[node.parentid]
                                else
                                    node = checknodes[node.parentid]
                                end
                            end
                            reverse!(temp)
                            push!(temp, (cnnode.id, curr.id))
                            node = curr
                            while node.lvl != 0
                                push!(temp, (node.id, node.parentid))
                                if isodd(node.lvl)
                                    node = varnodes[node.parentid]
                                else
                                    node = checknodes[node.parentid]
                                end
                            end
                            push!(cycles, temp)

                            # finish this level off but don't go deeper so remove children at lower level
                            if notemptied
                                while length(queue) > 0
                                    back = last(queue)
                                    if back.lvl != curr.lvl
                                        pop!(queue)
                                    else
                                        break
                                    end
                                end
                                notemptied = false
                            end
                        elseif notemptied
                            cnnode.lvl = curr.lvl + 1
                            cnnode.parentid = curr.id
                            cnnode.cumACE = curr.cumACE
                            push!(queue, cnnode)
                        end
                    end
                end
            else
                for vn in checkadlist[curr.id]
                     # can't pass messages back to the same node
                    if vn != curr.parentid
                        vnnode = varnodes[vn]
                        if vnnode.lvl != -1
                            # have seen before
                            push!(ACEs, curr.cumACE + vnnode.cumACE - root.localACE)

                             # trace the cycle from curr to root and cnnode to root
                            temp = Vector{Tuple{Int, Int}}()
                            node = vnnode
                            while node.lvl != 0
                                push!(temp, (node.parentid, node.id))
                                if isodd(node.lvl)
                                    node = varnodes[node.parentid]
                                else
                                    node = checknodes[node.parentid]
                                end
                            end
                            reverse!(temp)
                            push!(temp, (vnnode.id, curr.id))
                            node = curr
                            while node.lvl != 0
                                push!(temp, (node.id, node.parentid))
                                if isodd(node.lvl)
                                    node = varnodes[node.parentid]
                                else
                                    node = checknodes[node.parentid]
                                end
                            end
                            push!(cycles, temp)

                            # finish this level off but don't go deeper so remove children at lower level
                            if notemptied
                                while length(queue) > 0
                                    back = last(queue)
                                    if back.lvl != curr.lvl
                                        pop!(queue)
                                    else
                                        break
                                    end
                                end
                                notemptied = false
                            end
                        elseif notemptied
                            vnnode.lvl = curr.lvl + 1
                            vnnode.parentid = curr.id
                            vnnode.cumACE = curr.cumACE + vnnode.localACE
                            push!(queue, vnnode)
                        end
                    end
                end
            end
            popfirst!(queue)
        end

        min, index = findmin(ACEs)
        vsACE[i] = min
        cyclesvs[i] = cycles[index]
    end
    return vsACE, cyclesvs
end
shortestcycleACE(C::LDPCCode, v::Int) = shortestcycleACE(C, [v])[1]
shortestcycleACE(C::LDPCCode) = shortestcycleACE(C, collect(1:C.n))

"""
    shortestcycles(C::LDPCCode, v::Int)
    shortestcycles(C::LDPCCode, vs::Vector{Int})
    shortestcycles(C::LDPCCode)

Return all the cycles of shortest length in the Tanner graph of `C` for the vertex `v` or
vertices `vs`. If no vertices are given, all vertices are computed by default.

# Note
* The length of the shortest cycle is not necessarily the same for each vertex.
* To reduce computational complexity, the same cycle may appear under each vertex in the cycle.
"""
function shortestcycles(C::LDPCCode, vs::Vector{Int})
    isempty(vs) && throw(ArgumentError("Input variable node list cannot be empty"))
    all(x->1 <= x <= C.n, vs) || throw(DomainError("Variable node indices must be between 1 and length(C)"))

    checkadlist, varadlist = _nodeadjacencies(C.H)
    cyclesvs = [Vector{Vector{Tuple{Int, Int}}}() for _ in 1:length(vs)]

    Threads.@threads for i in 1:length(vs)
        # moving this inside allocates more but allows for multi-threading
        checknodes = [_ACEchecknode(i, -1, -1, -1) for i in 1:length(checkadlist)]
        varnodes = [_ACEvarnode(i, -1, -1, -1, length(varadlist[i]) - 2) for i in 1:C.n]
        cycles = Vector{Vector{Tuple{Int, Int}}}()
        notemptied = true
        root = varnodes[vs[i]]
        root.lvl = 0
        queue = Deque{Union{_ACEchecknode, _ACEvarnode}}()
        push!(queue, root)
        while length(queue) > 0
            curr = first(queue)
            if isa(curr, _ACEvarnode)
                for cn in varadlist[curr.id]
                    # can't pass messages back to the same node
                    if cn != curr.parentid
                        cnnode = checknodes[cn]
                        if cnnode.lvl != -1
                            # have seen before
                            # trace the cycle from curr to root and cnnode to root
                            temp = Vector{Tuple{Int, Int}}()
                            node = cnnode
                            while node.lvl != 0
                                push!(temp, (node.parentid, node.id))
                                if isodd(node.lvl)
                                    node = varnodes[node.parentid]
                                else
                                    node = checknodes[node.parentid]
                                end
                            end
                            reverse!(temp)
                            push!(temp, (cnnode.id, curr.id))
                            node = curr
                            while node.lvl != 0
                                push!(temp, (node.id, node.parentid))
                                if isodd(node.lvl)
                                    node = varnodes[node.parentid]
                                else
                                    node = checknodes[node.parentid]
                                end
                            end
                            push!(cycles, temp)

                            # finish this level off but don't go deeper so remove children at lower level
                            if notemptied
                                while length(queue) > 0
                                    back = last(queue)
                                    if back.lvl != curr.lvl
                                        pop!(queue)
                                    else
                                        break
                                    end
                                end
                                notemptied = false
                            end
                        elseif notemptied
                            cnnode.lvl = curr.lvl + 1
                            cnnode.parentid = curr.id
                            push!(queue, cnnode)
                        end
                    end
                end
            else
                for vn in checkadlist[curr.id]
                     # can't pass messages back to the same node
                    if vn != curr.parentid
                        vnnode = varnodes[vn]
                        if vnnode.lvl != -1
                            # have seen before
                            temp = Vector{Tuple{Int, Int}}()
                            node = vnnode
                            while node.lvl != 0
                                push!(temp, (node.parentid, node.id))
                                if isodd(node.lvl)
                                    node = varnodes[node.parentid]
                                else
                                    node = checknodes[node.parentid]
                                end
                            end
                            reverse!(temp)
                            push!(temp, (vnnode.id, curr.id))
                            node = curr
                            while node.lvl != 0
                                push!(temp, (node.id, node.parentid))
                                if isodd(node.lvl)
                                    node = varnodes[node.parentid]
                                else
                                    node = checknodes[node.parentid]
                                end
                            end
                            push!(cycles, temp)

                            # finish this level off but don't go deeper so remove children at lower level
                            if notemptied
                                while length(queue) > 0
                                    back = last(queue)
                                    if back.lvl != curr.lvl
                                        pop!(queue)
                                    else
                                        break
                                    end
                                end
                                notemptied = false
                            end
                        elseif notemptied
                            vnnode.lvl = curr.lvl + 1
                            vnnode.parentid = curr.id
                            push!(queue, vnnode)
                        end
                    end
                end
            end
            popfirst!(queue)
        end
        cyclesvs[i] = cycles
    end
    return cyclesvs
end
shortestcycles(C::LDPCCode, v::Int) = shortestcycles(C, [v])[1]
function shortestcycles(C::LDPCCode)
    cycles = shortestcycles(C, collect(1:C.n))
    girth = minimum([minimum([length(cycle) for cycle in cycles[i]]) for i in 1:C.n])
    if ismissing(C.girth)
        C.girth = girth
    else
        if C.girth != girth
            @warn "Known girth, $(C.girth), does not match just computed girth, $girth"
        end
    end
    return cycles
end

function _progressivenodeadjacencies(H::CTMatrixTypes, vs::Vector{Int}, vtype::Symbol)
    checkadlist, varadlist = _nodeadjacencies(H)
    unique!(sort!(vs))
    len = length(vs)
    checkadlists = [deepcopy(checkadlist) for _ in 1:len]
    varadlists = [deepcopy(varadlist) for _ in 1:len]

    for i in 2:len
        prev = vs[1:i - 1]
        if vtype == :v
            for (j, x) in enumerate(checkadlists[i])
                checkadlists[i][j] = setdiff(x, prev)
            end
        else
            for (j, x) in enumerate(varadlists[i])
                varadlists[i][j] = setdiff(x, prev)
            end
        end
    end
    return checkadlists, varadlists
end

"""
    countshortcycles(C::LDPCCode)

Return a bar graph and a dictionary of (length, count) pairs for unique short
cycles in the Tanner graph of `C`.

# Note
* Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
"""
function countshortcycles(C::LDPCCode)
    checkadlists, varadlists =  _progressivenodeadjacencies(C.H, collect(1:C.n), :v)
    lengths = [Vector{Int}() for _ in 1:C.n]
    Threads.@threads for i in 1:C.n
        checknodes = [_ACEchecknode(i, -1, -1, -1) for i in 1:length(checkadlists[i])]
        varnodes = [_ACEvarnode(i, -1, -1, -1, length(varadlists[i][i]) - 2) for i in 1:C.n]

        cyclens = Vector{Int}()
        root = varnodes[i]
        root.lvl = 0
        queue = Queue{Union{_ACEchecknode, _ACEvarnode}}()
        enqueue!(queue, root)
        while length(queue) > 0
            curr = first(queue)
            if isa(curr, _ACEvarnode)
                for cn in varadlists[i][curr.id]
                    # can't pass messages back to the same node
                    if cn != curr.parentid
                        cnnode = checknodes[cn]
                        if cnnode.lvl != -1
                            # have seen before
                            push!(cyclens, curr.lvl + cnnode.lvl + 1)
                        else
                            cnnode.lvl = curr.lvl + 1
                            cnnode.parentid = curr.id
                            enqueue!(queue, cnnode)
                        end
                    end
                end
            else
                for vn in checkadlists[i][curr.id]
                     # can't pass messages back to the same node
                    if vn != curr.parentid
                        vnnode = varnodes[vn]
                        if vnnode.lvl != -1
                            # have seen before
                            push!(cyclens, curr.lvl + vnnode.lvl + 1)
                        else
                            vnnode.lvl = curr.lvl + 1
                            vnnode.parentid = curr.id
                            enqueue!(queue, vnnode)
                        end
                    end
                end
            end
            dequeue!(queue)
        end
        lengths[i] = cyclens
    end

    girth = minimum([isempty(lengths[i]) ? 9999999 : minimum(lengths[i]) for i in 1:C.n])
    if ismissing(C.girth)
        C.girth = girth
    else
        if C.girth != girth
            @warn "Known girth, $(C.girth), does not match just computed girth, $girth"
        end
    end

    counts = Dict{Int, Int}()
    for i in girth:2:2 * girth - 2
        for j in 1:C.n
            if i ∈ keys(counts)
                counts[i] += count(x -> x == i, lengths[j])
            else
                counts[i] = count(x -> x == i, lengths[j])
            end
        end
    end

    len = length(counts)
    x_data = [0 for _ in 1:len]
    y_data = [0 for _ in 1:len]
    index = 1
    for (i, j) in counts
        x_data[index] = i
        y_data[index] = j
        index += 1
    end

    fig = bar(x_data, y_data, bar_width=1, xticks=x_data, yticks=y_data,
        legend=false, xlabel="Cycle Length", ylabel="Occurrences", title="Short Cycle Counts")
    show(fig)
    return fig, counts
end

"""
    ACEspectrumofnode(C::LDPCCode, v::Int)
    ACEspectrumofnode(C::LDPCCode, vs::Vector{Int})
    ACEspectrumofnode(C::LDPCCode)

Return the ACEs, average ACE, median ACE, mode ACE, cycle lengths, average length, median
length, and mode length of the cycles for vertex `v` or vertices `vs` of the Tanner graph
of `C`. If no vertices are given, all vertices are computed by default.

# Note
* In case of ties, the smallest tied value is returned.
"""
function ACEspectrumofnode(C::LDPCCode, vs::Vector{Int})
    # using the original DFS approach constructs a significantly larger tree than this truncated BFS approach

    isempty(vs) && throw(ArgumentError("Input node list cannot be empty"))
    all(x -> 1 <= x <= C.n, vs) || throw(DomainError("Variable node index must be between 1 and length(C)"))

    checkadlist, varadlist = _nodeadjacencies(C.H)

    vsACEs = [Vector{Int}() for _ in 1:length(vs)]
    avgACEs = zeros(Float64, length(vs))
    medACEs = zeros(Float64, length(vs))
    modeACEs = zeros(Int, length(vs))
    lengths = [Vector{Int}() for _ in 1:length(vs)]
    avglens = zeros(Float64, length(vs))
    medlens = zeros(Float64, length(vs))
    modelens = zeros(Int, length(vs))
 
    Threads.@threads for i in 1:length(vs)
        checknodes = [_ACEchecknode(i, -1, -1, -1) for i in 1:length(checkadlist)]
        varnodes = [_ACEvarnode(i, -1, -1, -1, length(varadlist[i]) - 2) for i in 1:C.n]

        ACEs = Vector{Int}()
        cyclens = Vector{Int}()
        root = varnodes[vs[i]]
        root.lvl = 0
        root.cumACE = root.localACE
        queue = Queue{Union{_ACEchecknode, _ACEvarnode}}()
        enqueue!(queue, root)
        while length(queue) > 0
            curr = first(queue)
            if isa(curr, _ACEvarnode)
                for cn in varadlist[curr.id]
                    # can't pass messages back to the same node
                    if cn != curr.parentid
                        cnnode = checknodes[cn]
                        if cnnode.lvl != -1
                            # have seen before
                            push!(ACEs, curr.cumACE + cnnode.cumACE - root.localACE)
                            push!(cyclens, curr.lvl + cnnode.lvl + 1)
                        else
                            cnnode.lvl = curr.lvl + 1
                            cnnode.parentid = curr.id
                            cnnode.cumACE = curr.cumACE
                            enqueue!(queue, cnnode)
                        end
                    end
                end
            else
                for vn in checkadlist[curr.id]
                     # can't pass messages back to the same node
                    if vn != curr.parentid
                        vnnode = varnodes[vn]
                        if vnnode.lvl != -1
                            # have seen before
                            push!(ACEs, curr.cumACE + vnnode.cumACE - root.localACE)
                            push!(cyclens, curr.lvl + vnnode.lvl + 1)
                        else
                            vnnode.lvl = curr.lvl + 1
                            vnnode.parentid = curr.id
                            vnnode.cumACE = curr.cumACE + vnnode.localACE
                            enqueue!(queue, vnnode)
                        end
                    end
                end
            end
            dequeue!(queue)
        end

        vsACEs[i] = ACEs
        avgACEs[i] = mean(ACEs)
        medACEs[i] = median(ACEs)
        modeACEs[i] = StatsBase.mode(sort(ACEs))
        lengths[i] = cyclens
        avglens[i] = mean(cyclens)
        medlens[i] = median(cyclens)
        modelens[i] = StatsBase.mode(sort(cyclens))
    end
    return vsACEs, avgACEs, medACEs, modeACEs, lengths, avglens, medlens, modelens
end
ACEspectrumofnode(C::LDPCCode, v::Int) = ACEspectrumofnode(C, [v]) # TODO: return slightly different than above
ACEspectrumofnode(C::LDPCCode) = ACEspectrumofnode(C, collect(1:C.n))

"""
    ACEspectrum(C::LDPCCode)

Return an interactive figure and data for the ACE spectrum of the Tanner graph of `C`.
"""
function ACEspectrum(C::LDPCCode)
    vsACEs, _, _, _, lengths, _, _, _ = ACEspectrumofnode(C, collect(1:C.n))
    # (false) spectrum: how many nodes have that ACE for that length
    # (true) spectrum: for a given length 4 <= l <= maximum(variabledegreedistribution(C)),
    # how many var nodes have shortest cycle that ACE

    shortestlens = [minimum(i) for i in lengths]
    girth = minimum(shortestlens)
    if ismissing(C.girth)
        C.girth = girth
    else
        if C.girth != girth
            @warn "Known girth, $(C.girth), does not match just computed girth, $girth"
        end
    end
    
    # TODO: remove WGLMakie as a default use and only use for interactive plots
    fig = Figure();
    ax = Axis(fig[1, 1], xlabel="ACE", ylabel="Occurrences", title="ACE Spectrum")
    sg = SliderGrid(fig[2, 1], (label = "Cycle Length", range=girth:2:2 * girth - 2, startvalue=4))

    xmax = maximum(reduce(vcat, vsACEs))
    ymax = 0
    counts = [Dict{Int, Int}() for _ in 1:length(girth:2:2 * girth - 2)]
    Xdata = Observable(Vector{Int}())
    Ydata = Observable(Vector{Int}())
    barplot!(ax, Xdata, Ydata, bar_width=1, xticks=Xdata, yticks=Ydata)
    for (k, l) in enumerate(girth:2:2 * girth - 2)
        for i in 1:length(shortestlens)
            if shortestlens[i] == l
                for j in 1:length(lengths[i])
                    if lengths[i][j] == l
                        if vsACEs[i][j] ∈ keys(counts[k])
                            counts[k][vsACEs[i][j]] += 1
                        else
                            counts[k][vsACEs[i][j]] = 1
                        end
                    end
                end
            end
        end
        ys = collect(values(counts[k]))
        ymax = maximum([ymax; ys])

        on(sg.sliders[1].value) do val
            if to_value(val) == l
                Xdata.val = collect(keys(counts[k]))
                Ydata.val = ys
                notify(Xdata)
                notify(Ydata)
            end
        end
    end

    GLMakie.limits!(0, xmax + 1, 0, ymax + 2)
    show(fig)
    return fig, counts
end

mutable struct _compgraphnode
    id::Int
    parentid::Int
    lvl::Int
    vertexnumber::Int
    type::Symbol
end

# doesn't seem to be a point in making this dynamic with a slider, as it simply
# continues in the same tree shape and no useful information is gained from watching it
"""
    computationgraph(C::LDPCCode, lvl::Int, v::Int, vtype::Symbol=:v)

Return a figure representing the expansion of the Tanner graph of `C` to level `lvl`
for node `v`. If `vtype` is `:v`, `v` is interpreted as a variable node; otherwise,
`vtype` is `:c` and `v` is interpreted as a check node.
"""
function computationgraph(C::LDPCCode, lvl::Int, v::Int, vtype::Symbol=:v)
    vtype ∈ (:v, :c) || throw(ArgumentError("Unknown argument for vtype"))
    if vtype == :v
        1 <= v <= C.n || throw(DomainError("Variable node index must be between 1 and length(C)"))
    else
        1 <= v <= nrows(C.H) || throw(DomainError("Check node index must be between 1 and the number of rows of H"))
    end
    lvl > 0 || throw(DomainError("Graph recipe requires at least one level"))

    checkadlist, varadlist = _nodeadjacencies(C.H)
    G = SimpleDiGraph()
    labels = Vector{String}()
    colors = Vector{Symbol}()
    markers = Vector{Symbol}()

    if vtype == :v
        root = _compgraphnode(v, -1, 0, 1, :v)
        Grphs.add_vertex!(G)
        push!(labels, L"v_{%$v}")
        push!(colors, :black)
        push!(markers, :circle)
    else
        root = _compgraphnode(v, -1, 0, 1, :c)
        Grphs.add_vertex!(G)
        push!(labels, L"c_{%$v}")
        push!(colors, :red)
        push!(markers, :rect)
    end
    
    queue = Queue{_compgraphnode}()
    enqueue!(queue, root)
    while length(queue) > 0
        curr = first(queue)
        curr.lvl == lvl && break
        newlvl = curr.lvl + 1
        if curr.type != :c
            for cn in varadlist[curr.id]
                if cn != curr.parentid
                    Grphs.add_vertex!(G)
                    cnnode = _compgraphnode(cn, curr.id, newlvl, Grphs.nv(G), :c)
                    Grphs.add_edge!(G, curr.vertexnumber, cnnode.vertexnumber)
                    enqueue!(queue, cnnode)
                    push!(labels, L"c_{%$cn}")
                    push!(colors, :red)
                    push!(markers, :rect)
                end
            end
        else
            for vn in checkadlist[curr.id]
                if vn != curr.parentid
                    Grphs.add_vertex!(G)
                    vnnode = _compgraphnode(vn, curr.id, newlvl, Grphs.nv(G), :v)
                    Grphs.add_edge!(G, curr.vertexnumber, vnnode.vertexnumber)
                    enqueue!(queue, vnnode)
                    push!(labels, L"v_{%$vn}")
                    push!(colors, :black)
                    push!(markers, :circle)
                end
            end
        end
        dequeue!(queue)
    end

    # TODO: fix plot - count number of added nodes in level and manually pass in a calculated image size
    f, ax, p = graphplot(G, layout=Buchheim(),
        nlabels=labels,
        node_marker=markers,
        node_color=colors,
        nlabels_textsize=10,
        nlabels_align=(:left, :center),
        nlabels_distance=7);
    hidedecorations!(ax)
    hidespines!(ax)
    show(f)
    # TODO: what do we want to return here and make uniform with the doc string
    return f, ax, p
end

# TODO: ACE distribution function
# store results
