# Copyright (c) 2022, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

mutable struct LDPCCode <: AbstractLDPCCode
    F::Union{FqNmodFiniteField}
    n::Integer
    k::Integer
    d::Union{Integer, Missing}
    G::fq_nmod_mat
    Gorig::Union{fq_nmod_mat, Missing}
    H::fq_nmod_mat
    Horig::Union{fq_nmod_mat, Missing}
    Gstand::fq_nmod_mat
    Hstand::fq_nmod_mat
    weightenum::Union{WeightEnumerator, Missing}
    numedges::Int
    vardegs::Vector{Int}
    checkdegs::Vector{Int}
    colbound::Int
    rowbound::Int
    density::Float64
    isreg::Bool
    tangr::Union{Figure, Missing}
    λ::fmpq_poly
    ρ::fmpq_poly
end

"""
    Tannergraph(H::fq_nmod_mat)

Return the Tanner graph of the matrix `H` as a `Figure` object.
"""
function Tannergraph(H::fq_nmod_mat)
    # convert H to A
    M = FpmattoJulia(H)
    nr, nc = size(M)
    A = zeros(Int, nr + nc, nr + nc)
    # bottom left corner
    # need to threshold any nonzero to a 1
    # A[nc + 1:end, 1:nc] = M
    # no, put in top right corner in order to get parents, childs working
    A[1:nc, nc + 1:end] = transpose(M)

    f = Figure();
    ax = Axis(f[1, 1], yreversed = true, xautolimitmargin = (0.15, 0.20),
        yautolimitmargin = (0.15, 0.20))
    hidespines!(ax)
    hidedecorations!(ax)

    leftx, lefty = zeros(nc), 1.:nc
    rightx, righty = ones(nr) * nr, range(1, nc, nr)
    x = vcat(leftx, rightx)
    y = vcat(lefty, righty)
    points = Point.(zip(x, y))
    cols = (:aqua, :red, :orange, :green, :blue, :purple)

    G = SimpleDiGraph(A)
    parents = [inneighbors(G, i) for i in Graphs.vertices(G)]
    childs = findall(x -> length(x) > 0, parents)
    # println(parents)
    # println(childs)

    for (i, v) in enumerate(childs)
        for node in parents[v]
            lines!(Point(x[[node, v]]...), Point(y[[node, v]]...),
                   color=cols[i % 6 + 1], linewidth=5)
        end
        text!(points[v], text=L"h_%$i", offset=(20, -15))
    end

    for (i, point) in enumerate(points[1:nc])
        CairoMakie.scatter!(point, color=:black, marker=:circle, markersize=25)
        text!(point, text=L"v_%$i", offset=(-30, -10))
    end

    for (i, point) in enumerate(points[nc + 1:end])
        CairoMakie.scatter!(point, color=:black, marker=:rect, markersize=25)
    end
    f
    return f
    # save("test.png", f)
end

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
    density(C::AbstractLDPCCode)

Return the density of the parity-check matrix of `C`.
"""
density(C::AbstractLDPCCode) = C.density

"""
    isregular(C::AbstractLDPCCode)

Return `true` if the `C` is a regular LDPC code.

An LDPC is regular if all the column degrees and equal and all the row degrees
are equal.
"""
isregular(C::AbstractLDPCCode) = C.isreg

"""
    Tannergraph(C::AbstractLDPCCode)

Return the Tanner graph of `C` as a `Figure` object.
"""
Tannergraph(C::AbstractLDPCCode) = ismissing(C.tangr) ? (return Tannergraph(C.H);) : (return C.tangr;)

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

function _degreedistribution(H::fq_nmod_mat)
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

function _density(H::fq_nmod_mat)
    count = 0
    nr, nc = size(H)
    for c in 1:nc
        for r in 1:nr
            !iszero(H[r, c]) && (count += 1;)
        end
    end
    return count, count / (nr * nc)
end

"""
    LDPCCode(H::fq_nmod_mat)

Return the LDPC code defined by the parity-check matrix `H`.

LDPC codes are typically required to have a matrix density of less than 1%.
"""
function LDPCCode(H::fq_nmod_mat)
    nnz, den = _density(H)
    den <= 0.01 || (@warn "LDPC codes (generally) require a density of less than 1%.";)

    C = LinearCode(H, true)
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

    R, x = PolynomialRing(Nemo.QQ, "x")
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

    return LDPCCode(C.F, C.n, C.k, C.d, C.G, C.Gorig, C.H, C.Horig, C.Gstand,
        C.Hstand, C.weightenum, nnz, cols, rows, c, r, den, isreg, missing,
        colpoly, rowpoly)
end

"""
    LDPCCode(C::AbstractLinearCode)

Return the LDPC code given by `C`.

LDPC codes are typically required to have a matrix density of less than 1%.
"""
function LDPCCode(C::AbstractLinearCode)
    nnz, den = _density(C.H)
    den <= 0.01 || (@warn "LDPC codes (generally) require a density of less than 1%.";)

    cols, rows = _degreedistribution(C.H)
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

    R, x = PolynomialRing(Nemo.QQ, "x")
    colpoly = R(0)
    for i in cols
        colpoly += x^i
    end
    colpoly = divexact(colpoly, nnz)
    rowpoly = R(0)
    for i in rows
        rowpoly += x^i
    end
    rowpoly = divexact(rowpoly, nnz)

    return LDPCCode(C.F, C.n, C.k, C.d, C.G, C.Gorig, C.H, C.Horig, C.Gstand,
        C.Hstand, C.weightenum, nnz, cols, rows, c, r, den, isreg, missing,
        colpoly, rowpoly)
end

function show(io::IO, C::AbstractLDPCCode)
    if get(io, :compact, false)
        if ismissing(C.d)
            if C.isreg
                println(io, "[$(C.n), $(C.k)]_$(order(C.F)) regular ($(C.colbound), $(C.rowbound))-LDPC code with density $(C.density).")
            else
                println(io, "[$(C.n), $(C.k)]_$(order(C.F)) irregular ($(C.colbound), $(C.rowbound))-LDPC code with density $(C.density).")
            end
        else
            if C.isreg
                println(io, "[$(C.n), $(C.k), $(C.d)]_$(order(C.F)) regular ($(C.colbound), $(C.rowbound))-LDPC code with density $(C.density).")
            else
                println(io, "[$(C.n), $(C.k), $(C.d)]_$(order(C.F)) irregular ($(C.colbound), $(C.rowbound))-LDPC code with density $(C.density).")
            end
        end
    else
        if ismissing(C.d)
            if C.isreg
                println(io, "[$(C.n), $(C.k)]_$(order(C.F)) regular ($(C.colbound), $(C.rowbound))-LDPC code with density $(C.density).")
            else
                println(io, "[$(C.n), $(C.k)]_$(order(C.F)) irregular ($(C.colbound), $(C.rowbound))-LDPC code with density $(C.density).")
            end
        else
            if C.isreg
                println(io, "[$(C.n), $(C.k), $(C.d)]_$(order(C.F)) regular ($(C.colbound), $(C.rowbound))-LDPC code with density $(C.density).")
            else
                println(io, "[$(C.n), $(C.k), $(C.d)]_$(order(C.F)) irregular ($(C.colbound), $(C.rowbound))-LDPC code with density $(C.density).")
            end
        end
        nr, nc = size(C.Horig)
        println(io, "Parity-check matrix: $nr × $nc")
        for i in 1:nr
            print(io, "\t")
            for j in 1:nc
                if j != nc
                    print(io, "$(C.Horig[i, j]) ")
                elseif j == nc && i != nr
                    println(io, "$(C.Horig[i, j])")
                else
                    print(io, "$(C.Horig[i, j])")
                end
            end
        end
        if !ismissing(C.weightenum)
            println(io, "\nComplete weight enumerator:")
            println(io, "\t", C.weightenum.polynomial)
        end
        println(io, "\nVariable degree polynomial:")
        println(io, "\t", C.λ)
        println(io, "Check degree polynomial:")
        print(io, "\t", C.ρ)
    end
end

"""
    degreedistributionssplot(C::AbstractLDPCCode)

Return a bar plot of the column and row degree distributions of `C`.
"""
function degreedistributionssplot(C::AbstractLDPCCode)
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
