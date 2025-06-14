# Copyright (c) 2022, 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
# constructors
#############################

"""
    LDPCCode(H::fqPolyRepMatrix)

Return the LDPC code defined by the parity-check matrix `H`.
"""
function LDPCCode(H::CTMatrixTypes)
    # TODO: remove empties
    nnz, den = _density(H)
    # den <= 0.01 || (@warn "LDPC codes typically expect a density of less than 1%.";)

    cols, rows = _degree_distribution(H)
    is_reg = true
    c1 = cols[1]
    for i = 2:length(cols)
        c1 == cols[i] || (is_reg = false; break;)
    end
    if is_reg
        r1 = rows[1]
        for i = 2:length(rows)
            r1 == rows[i] || (is_reg = false; break;)
        end
    end
    c, r = maximum(cols), maximum(rows)

    R, x = polynomial_ring(Nemo.QQ, :x)
    col_poly = R(0)
    for i in cols
        col_poly += i * x^(i - 1)
    end
    col_poly = divexact(col_poly, nnz)
    row_poly = R(0)
    for i in rows
        row_poly += i * x^(i - 1)
    end
    row_poly = divexact(row_poly, nnz)

    C = LinearCode(H, true)
    return LDPCCode(
        base_ring(H),
        C.n,
        C.k,
        C.d,
        C.l_bound,
        C.u_bound,
        H,
        nnz,
        cols,
        rows,
        c,
        r,
        maximum([c, r]),
        den,
        is_reg,
        col_poly,
        row_poly,
        missing,
        [Vector{Int}() for _ = 1:C.n],
        Vector{Vector{Int}}(),
        0,
    )
end

"""
    LDPCCode(C::AbstractLinearCode)

Return the LDPC code given by the parity-check matrix of `C`.
"""
LDPCCode(C::AbstractLinearCode) = LDPCCode(parity_check_matrix(C))

"""
    regular_LDPC_code(q::Int, n::Int, l::Int, r::Int [; seed=nothing])

Return a random regular LDPC code over `GF(q)` of length `n` with column degree `l`
and row degree `r`.

If a seed is given, i.e. `regular_LDPC_Code(4, 1200, 3, 6, seed=123)`, the
results are reproducible.
"""
function regular_LDPC_code(
    q::Int,
    n::Int,
    l::Int,
    r::Int;
    seed::Union{Nothing,Int} = nothing,
)
    Random.seed!(seed)
    m = divexact(n * l, r)
    F = if is_prime(q)
        Oscar.Nemo.Native.GF(q)
    else
        factors = Nemo.factor(q)
        length(factors) == 1 ||
            throw(DomainError("There is no finite field of order $q"))
        (p, t), = factors
        GF(p, t, :α)
    end
    elems = collect(F)[2:end]
    H = zero_matrix(F, m, n)
    col_sums = zeros(Int, n)
    for i in axes(H, 1)
        ind = reduce(vcat, shuffle(filter(k -> col_sums[k] == s, 1:n)) for s = 0:(l-1))[1:r]
        for j in ind
            H[i, j] = rand(elems)
        end
        col_sums[ind] .+= 1
    end
    @assert all(count(.!iszero.(H[:, j])) == l for j in axes(H, 2))
    @assert all(count(.!iszero.(H[i, :])) == r for i in axes(H, 1))

    R, x = polynomial_ring(Nemo.QQ, :x)
    C = LinearCode(H, true)
    return LDPCCode(
        C.F,
        C.n,
        C.k,
        C.d,
        C.l_bound,
        C.u_bound,
        H,
        n * l,
        l * ones(Int, n),
        r * ones(Int, m),
        l,
        r,
        max(l, r),
        r / n,
        true,
        (1 // l) * x^l,
        (1 // r) * x^r,
        missing,
        [Vector{Int}() for _ = 1:C.n],
        Vector{Vector{Int}}(),
        0,
    )
end
regular_LDPC_code(n::Int, l::Int, r::Int; seed::Union{Nothing,Int} = nothing) =
    regular_LDPC_code(2, n, l, r, seed = seed)

#############################
# getter functions
#############################

"""
    variable_degree_distribution(C::AbstractLDPCCode)

Return the variable node degree distribution of `C`.
"""
variable_degree_distribution(C::AbstractLDPCCode) = C.var_degs

"""
    check_degree_distribution(C::AbstractLDPCCode)

Return the check node degree distribution of `C`.
"""
check_degree_distribution(C::AbstractLDPCCode) = C.check_degs

"""
    degree_distributions(C::AbstractLDPCCode)

Return the variable and check node degree distributions of `C`.
"""
degree_distributions(C::AbstractLDPCCode) = C.var_degs, C.check_degs

"""
    column_bound(C::AbstractLDPCCode)

Return the column bound `c` of the `(c, r)`-LDPC code `C`.
"""
column_bound(C::AbstractLDPCCode) = C.col_bound

"""
    row_bound(C::AbstractLDPCCode)

Return the row bound `r` of the `(c, r)`-LDPC code `C`.
"""
row_bound(C::AbstractLDPCCode) = C.row_bound

"""
    column_row_bounds(C::AbstractLDPCCode)

Return the column and row bounds `c, r` of the `(c, r)`-LDPC code `C`.
"""
column_row_bounds(C::AbstractLDPCCode) = C.col_bound, C.row_bound

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
    is_regular(C::AbstractLDPCCode)

Return `true` if the `C` is a regular LDPC code.

# Notes
- An LDPC is regular if all the column degrees and equal and all the row degrees are equal.
"""
is_regular(C::AbstractLDPCCode) = C.is_reg

"""
    variable_degree_polynomial(C::AbstractLDPCCode)

Return the variable degree polynomial of `C`.
"""
variable_degree_polynomial(C::AbstractLDPCCode) = C.λ

"""
    check_degree_polynomial(C::AbstractLDPCCode)

Return the check degree polynomial of `C`.
"""
check_degree_polynomial(C::AbstractLDPCCode) = C.ρ

#############################
# setter functions
#############################

#############################
# general functions
#############################

function _degree_distribution(
    H::Union{CTMatrixTypes,MatElem{EuclideanRingResidueRingElem{fpPolyRingElem}}},
)

    nr, nc = size(H)
    cols = zeros(Int, 1, nc)
    @inbounds @views @simd for i = 1:nc
        # cols[i] = wt(H[:,  i])
        cols[i] = count(x -> !iszero(x), H[:, i])
    end
    rows = zeros(Int, 1, nr)
    @inbounds @views @simd for i = 1:nr
        # rows[i] = wt(H[i,  :])
        rows[i] = count(x -> !iszero(x), H[i, :])
    end
    return vec(cols), vec(rows)
end

function _density(H::CTMatrixTypes)
    count = 0
    nr, nc = size(H)
    for c = 1:nc
        for r = 1:nr
            !iszero(H[r, c]) && (count += 1;)
        end
    end
    return count, count / (nr * nc)
end

# TODO: make uniform with others
function show(io::IO, C::AbstractLDPCCode)
    if ismissing(C.d)
        if C.is_reg
            println(
                io,
                "[$(C.n), $(C.k)]_$(order(C.F)) regular ($(C.col_bound), $(C.row_bound))-LDPC code with density $(C.density).",
            )
        else
            println(
                io,
                "[$(C.n), $(C.k)]_$(order(C.F)) irregular $(C.limited)-limited LDPC code with density $(C.density).",
            )
        end
    else
        if C.is_reg
            println(
                io,
                "[$(C.n), $(C.k), $(C.d)]_$(order(C.F)) regular ($(C.col_bound), $(C.row_bound))-LDPC code with density $(C.density).",
            )
        else
            println(
                io,
                "[$(C.n), $(C.k), $(C.d)]_$(order(C.F)) irregular $(C.limited)-limited LDPC code with density $(C.density).",
            )
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
            for i = 1:nr
                print(io, "\t")
                for j = 1:nc
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
    degree_distributions_plot(C::AbstractLDPCCode)

Return a bar plot of the column and row degree distributions of `C`.

# Note
- Run `using Makie` to activate this extension.
"""
function degree_distributions_plot end
