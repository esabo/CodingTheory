# Copyright (c) 2022 - 2026 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
$(TYPEDSIGNATURES)

Return the LDPC code defined by the parity-check matrix `H`.
"""
function LDPCCode(H::CTMatrixTypes)
    nr, nc = size(H)
    
    # Uses the O(|E|) sparse-aware functions we just wrote!
    cols, rows = _degree_distribution(H)
    num_edges, den = _density(H)
    
    c_min, c_max = minimum(cols), maximum(cols)
    r_min, r_max = minimum(rows), maximum(rows)
    is_reg = (c_min == c_max) && (r_min == r_max)
    
    R, x = polynomial_ring(Nemo.QQ, :x)
    col_poly = divexact(sum(i * x^(i - 1) for i in cols), num_edges)
    row_poly = divexact(sum(i * x^(i - 1) for i in rows), num_edges)
    
    k_design = max(1, nc - nr)
    
    cache = Dict{Symbol, Any}(
        :col_degs => cols,
        :row_degs => rows,
        :c_bound => c_max,
        :r_bound => r_max,
        :density => den,
        :num_edges => num_edges
    )
    
    return LDPCCode(base_ring(H), nc, missing, k_design, missing, 1, nc, H, col_poly, row_poly, is_reg, cache)
end

"""
$(TYPEDSIGNATURES)

Return the LDPC code given by the parity-check matrix of `C`.
"""
LDPCCode(C::AbstractLinearCode) = LDPCCode(parity_check_matrix(C))

"""
$(TYPEDSIGNATURES)

Return a random regular LDPC code over `GF(q)` of length `n` with column degree `l`
and row degree `r`.
"""
function regular_LDPC_code(q::Int, n::Int, l::Int, r::Int; seed::Union{Nothing, Int} = nothing)
    !isnothing(seed) && Random.seed!(seed)
    
    m = divexact(n * l, r)
    F = if is_prime(q)
        Oscar.Nemo.Native.GF(q)
    else
        factors = Nemo.factor(q)
        length(factors) == 1 || throw(DomainError("There is no finite field of order $q"))
        (p, t), = factors
        GF(p, t, :α)
    end
    
    elems = collect(F)[2:end]
    
    # Coordinate accumulators for ultra-fast sparse matrix construction
    I_idx = Int[]
    J_idx = Int[]
    V_val = typeof(F(1))[]
    
    col_sums = zeros(Int, n)
    for i in 1:m
        ind = reduce(vcat, shuffle(filter(k -> col_sums[k] == s, 1:n)) for s in 0:l - 1)[1:r]
        for j in ind
            push!(I_idx, i)
            push!(J_idx, j)
            push!(V_val, rand(elems))
        end
        col_sums[ind] .+= 1
    end
    
    # Construct the sparse matrix natively (Oscar's native sparse constructor)
    H = sparse_matrix(F, m, n, I_idx, J_idx, V_val)
    
    R, x = polynomial_ring(Nemo.QQ, :x)
    k_design = max(1, n - m)
    den = r / n
    num_edges = n * l
    
    cache = Dict{Symbol, Any}(
        :col_degs => fill(l, n),
        :row_degs => fill(r, m),
        :c_bound => l,
        :r_bound => r,
        :density => den,
        :num_edges => num_edges
    )
    
    return LDPCCode(F, n, missing, k_design, missing, 1, n, H, (1 // l) * x^l, (1 // r) * x^r, true, cache)
end

#############################
      # getter functions
#############################

"""
    variable_degree_distribution(C::AbstractLDPCCode)

Return the variable node degree distribution of `C`.
"""
variable_degree_distribution(C::LDPCCode) = C.cache[:col_degs]
"""
    check_degree_distribution(C::AbstractLDPCCode)

Return the check node degree distribution of `C`.
"""
check_degree_distribution(C::LDPCCode) = C.cache[:row_degs]

"""
    degree_distributions(C::AbstractLDPCCode)

Return the variable and check node degree distributions of `C`.
"""
degree_distributions(C::LDPCCode) = (C.cache[:col_degs], C.cache[:row_degs])

"""
    column_bound(C::AbstractLDPCCode)

Return the column bound `c` of the `(c, r)`-LDPC code `C`.
"""
column_bound(C::LDPCCode) = C.cache[:c_bound]

"""
    row_bound(C::AbstractLDPCCode)

Return the row bound `r` of the `(c, r)`-LDPC code `C`.
"""
row_bound(C::LDPCCode) = C.cache[:r_bound]

"""
    column_row_bounds(C::AbstractLDPCCode)

Return the column and row bounds `c, r` of the `(c, r)`-LDPC code `C`.
"""
column_row_bounds(C::LDPCCode) = (C.cache[:c_bound], C.cache[:r_bound])

"""
    limited(C::AbstractLDPCCode)

Return the maximum of the row and column bounds for `C`.
"""
limited(C::LDPCCode) = max(C.cache[:c_bound], C.cache[:r_bound])

"""
    density(C::AbstractLDPCCode)

Return the density of the parity-check matrix of `C`.
"""
density(C::LDPCCode) = C.cache[:density]

"""
    is_regular(C::AbstractLDPCCode)

Return `true` if the `C` is a regular LDPC code.

# Notes
- An LDPC is regular if all the column degrees and equal and all the row degrees are equal.
"""
is_regular(C::LDPCCode) = C.is_reg

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

"""
$(TYPEDSIGNATURES)

Return the exact rank-adjusted dimension `k` of the LDPC code `C`.

# Notes
* If not previously computed, this runs the `O(n^3)` rank computation, caches the 
  result in `C.k`, and returns it.
"""
function dimension(C::LDPCCode)
    if ismissing(C.k)
        C.k = C.n - rank(C.H)
    end
    return C.k
end

"""
$(TYPEDSIGNATURES)

Return the design dimension of the LDPC code, computed instantly as `n - m`.
"""
design_dimension(C::LDPCCode) = C.k_design

design_rate(C::LDPCCode) = C.k_design / C.n
rate(C::LDPCCode) = dimension(C) / C.n

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

function _degree_distribution(H::CTMatrixTypes)
    nr, nc = size(H)
    cols = zeros(Int, nc)
    rows = zeros(Int, nr)
    
    # Check if the matrix is a sparse Oscar matrix (SMatElem)
    if typeof(H) <: SMatElem
        for (r, row) in enumerate(H)
            for (c, val) in row
                if !iszero(val)
                    cols[c] += 1
                    rows[r] += 1
                end
            end
        end
    else
        # Fallback for dense matrices
        for c in 1:nc
            for r in 1:nr
                if !iszero(H[r, c])
                    cols[c] += 1
                    rows[r] += 1
                end
            end
        end
    end
    return cols, rows
end

function _density(H::CTMatrixTypes)
    nr, nc = size(H)
    
    if typeof(H) <: SMatElem
        # Most sparse matrix implementations track the number of non-zeros intrinsically
        count = nnz(H)
        return count, count / (nr * nc)
    else
        count = 0
        for c in 1:nc
            for r in 1:nr
                if !iszero(H[r, c])
                    count += 1
                end
            end
        end
        return count, count / (nr * nc)
    end
end

function Base.show(io::IO, C::AbstractLDPCCode)
    # Safely extract via generic getters
    den = density(C)
    cb, rb = column_row_bounds(C)
    lim = limited(C)
    
    # Safely handle the dimension string depending on what the struct exposes
    if hasproperty(C, :k) && !ismissing(C.k)
        k_str = "$(C.k)"
    elseif hasproperty(C, :k_design)
        k_str = "$(C.k_design) (design)"
    else
        k_str = "?"
    end
    
    d_str = (hasproperty(C, :d) && !ismissing(C.d)) ? ", $(C.d)" : ""
    F_order = hasproperty(C, :F) ? order(C.F) : "?"
    
    if is_regular(C)
        println(io, "[$(C.n), $k_str$d_str]_$F_order regular ($cb, $rb)-LDPC code with density $den.")
    else
        println(io, "[$(C.n), $k_str$d_str]_$F_order irregular $lim-limited LDPC code with density $den.")
    end
    
    if get(io, :compact, true)
        println(io, "\nVariable degree polynomial:")
        println(io, "\t", variable_degree_polynomial(C))
        println(io, "Check degree polynomial:")
        println(io, "\t", check_degree_polynomial(C))
        
        if C.n <= 30
            # Route through the universal getter
            H = parity_check_matrix(C)
            nr, nc = size(H)
            println(io, "Parity-check matrix: $nr × $nc")
            
            for i in 1:nr
                print(io, "\t")
                for j in 1:nc
                    print(io, "$(H[i, j])")
                    j != nc && print(io, " ")
                end
                println(io)
            end
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
