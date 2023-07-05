# Copyright (c) 2022, 2023 Eric Sabo
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

# Notes
* LDPC codes are typically required to have a matrix density of less than 1%.
"""
function LDPCCode(H::CTMatrixTypes)
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

    return LDPCCode(C, nnz, cols, rows, c, r, maximum([c, r]), den, isreg,
        missing, colpoly, rowpoly)
end

"""
    LDPCCode(C::AbstractLinearCode)

Return the LDPC code given by `C`.

# Notes
* LDPC codes are typically required to have a matrix density of less than 1%.
"""
function LDPCCode(C::AbstractLinearCode)
    H = paritycheckmatrix(C)
    nnz, den = _density(H)
    den <= 0.01 || (@warn "LDPC codes (generally) require a density of less than 1%.";)

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
        colpoly += x^i
    end
    colpoly = divexact(colpoly, nnz)
    rowpoly = R(0)
    for i in rows
        rowpoly += x^i
    end
    rowpoly = divexact(rowpoly, nnz)

    return LDPCCode(C, nnz, cols, rows, c, r, maximum([c, r]), den, isreg, missing,
        colpoly, rowpoly)
end

"""
    regularLDPCCode(q::Int, n::Int, l::Int, r::Int [; seed = nothing])

Return a random regular LDPC code over GF(q) of length `n` with parameters `l` and `r`.

If a seed is given, i.e. `regulardLDPCCode(4, 1200, 3, 6, seed = 123)`, the
results are reproducible.
"""
function regularLDPCcode(q::Int, n::Int, l::Int, r::Int; seed::Union{Nothing, Int} = nothing)
    q == 2 && (return regularLDPCcode(n, l, r, seed = seed);)
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
    return LDPCCode(LinearCode(H, true), n * l, l * ones(Int, n), r * ones(Int, m),
                    l, r, max(l, r), r / n, true, missing, (1 // l) * x^l, (1 // r) * x^r)
end

"""
    regularLDPCCode(n::Int, l::Int, r::Int [; seed = nothing])

Return a random regular LDPC code over GF(2) of length `n` with parameters `l` and `r`.

If a seed is given, i.e. `regulardLDPCCode(1200, 3, 6, seed = 123)`, the
results are reproducible.
"""
function regularLDPCcode(n::Int, l::Int, r::Int; seed::Union{Nothing, Int} = nothing)
    Random.seed!(seed)
    m = divexact(n * l, r)
    F = GF(2)
    Fone = F(1)
    H = zero_matrix(F, m, n);
    colsums = zeros(Int, n)
    for i in axes(H, 1)
        ind = reduce(vcat, shuffle(filter(k -> colsums[k] == s, 1:n)) for s in 0:l-1)[1:r]
        for j in ind
            H[i, j] = Fone
        end
        colsums[ind] .+= 1
    end
    @assert all(count(.!iszero.(H[:, j])) == l for j in axes(H, 2))
    @assert all(count(.!iszero.(H[i, :])) == r for i in axes(H, 1))

    R, x = PolynomialRing(Nemo.QQ, :x)
    return LDPCCode(LinearCode(H, true), n * l, l * ones(Int, n), r * ones(Int, m),
                    l, r, max(l, r), r / n, true, missing, (1 // l) * x^l, (1 // r) * x^r)
end

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
    if ismissing(C.C.d)
        if C.isreg
            println(io, "[$(C.C.n), $(C.C.k)]_$(order(C.C.F)) regular ($(C.colbound), $(C.rowbound))-LDPC code with density $(C.density).")
        else
            println(io, "[$(C.C.n), $(C.C.k)]_$(order(C.C.F)) irregular $(C.limited)-limited LDPC code with density $(C.density).")
        end
    else
        if C.isreg
            println(io, "[$(C.C.n), $(C.C.k), $(C.C.d)]_$(order(C.C.F)) regular ($(C.colbound), $(C.rowbound))-LDPC code with density $(C.density).")
        else
            println(io, "[$(C.C.n), $(C.C.k), $(C.C.d)]_$(order(C.C.F)) irregular $(C.limited)-limited LDPC code with density $(C.density).")
        end
    end
    if get(io, :compact, true)
        println(io, "\nVariable degree polynomial:")
        println(io, "\t", C.λ)
        println(io, "Check degree polynomial:")
        print(io, "\t", C.ρ)
        if C.C.n <= 30
            # was using Horig here, which is probably what I want
            H = paritycheckmatrix(C.C)
            nr, nc = size(H)
            println(io, "Parity-check matrix: $nr × $nc")
            for i in 1:nr
                print(io, "\t")
                for j in 1:nc
                    if j != nc
                        print(io, "$(H[i, j]) ")
                    elseif j == nc && i != nr
                        println(io, "$(H[i, j])")
                    else
                        print(io, "$(H[i, j])")
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
