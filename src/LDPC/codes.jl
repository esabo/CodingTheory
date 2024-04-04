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

    cols, rows = _degree_distribution(H)
    is_reg = true
    c1 = cols[1]
    for i in 2:length(cols)
        c1 == cols[i] || (is_reg = false; break;)
    end
    if is_reg
        r1 = rows[1]
        for i in 2:length(rows)
            r1 == rows[i] || (is_reg = false; break;)
        end
    end
    c, r = maximum(cols), maximum(rows)

    R, x = PolynomialRing(Nemo.QQ, :x)
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
    return LDPCCode(base_ring(H), C.n, C.k, C.d, C.l_bound, C.u_bound, H, nnz,
        cols, rows, c, r, maximum([c, r]), den, is_reg, col_poly,
        row_poly, missing, [Vector{Int}() for _ in 1:C.n], Vector{Vector{Int}}())
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

If a seed is given, i.e. `regulardLDPCCode(4, 1200, 3, 6, seed=123)`, the
results are reproducible.
"""
function regular_LDPC_code(q::Int, n::Int, l::Int, r::Int; seed::Union{Nothing, Int} = nothing)
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
    H = zero_matrix(F, m, n)
    col_sums = zeros(Int, n)
    for i in axes(H, 1)
        ind = reduce(vcat, shuffle(filter(k -> col_sums[k] == s, 1:n)) for s in 0:l - 1)[1:r]
        for j in ind
            H[i, j] = rand(elems)
        end
        col_sums[ind] .+= 1
    end
    @assert all(count(.!iszero.(H[:, j])) == l for j in axes(H, 2))
    @assert all(count(.!iszero.(H[i, :])) == r for i in axes(H, 1))

    R, x = PolynomialRing(Nemo.QQ, :x)
    C = LinearCode(H, true)
    return LDPCCode(C.F, C.n, C.k, C.d, C.l_bound, C.u_bound, H, n * l, l * ones(Int, n),
        r * ones(Int, m), l, r, max(l, r), r / n, true, (1 // l) * x^l, (1 // r) * x^r, missing,
        [Vector{Int}() for _ in 1:C.n], Vector{Vector{Int}}())
end
regular_LDPC_code(n::Int, l::Int, r::Int; seed::Union{Nothing, Int} = nothing) =
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

function _degree_distribution(H::Union{CTMatrixTypes,
    MatElem{AbstractAlgebra.Generic.ResidueRingElem{fpPolyRingElem}}})

    nr, nc = size(H)
    cols = zeros(Int, 1, nc)
    @inbounds @views @simd for i in 1:nc
        # cols[i] = wt(H[:,  i])
        cols[i] = count(x -> !iszero(x), H[:, i])
    end
    rows = zeros(Int, 1, nr)
    @inbounds @views @simd for i in 1:nr
        # rows[i] = wt(H[i,  :])
        rows[i] = count(x -> !iszero(x), H[i, :])
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
        if C.is_reg
            println(io, "[$(C.n), $(C.k)]_$(order(C.F)) regular ($(C.col_bound), $(C.row_bound))-LDPC code with density $(C.density).")
        else
            println(io, "[$(C.n), $(C.k)]_$(order(C.F)) irregular $(C.limited)-limited LDPC code with density $(C.density).")
        end
    else
        if C.is_reg
            println(io, "[$(C.n), $(C.k), $(C.d)]_$(order(C.F)) regular ($(C.col_bound), $(C.row_bound))-LDPC code with density $(C.density).")
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
    degree_distributions_plot(C::AbstractLDPCCode)

Return a bar plot of the column and row degree distributions of `C`.

# Note
- Run `using Makie` to activate this extension.
"""
function degree_distributions_plot end

"""
    girth(C::AbstractLDPCCode; max_iter::Int = 100)

Return the girth of the Tanner graph of `C`.

An error is thrown if the maximum number of iterations is reached and
`-1` is returned to represent infinite girth.
"""
function girth(C::AbstractLDPCCode; max_iter::Int = 100)
    check_adj_list, var_adj_list = _node_adjacencies(C.H)
    girth_arr = zeros(Int, C.n)
    Threads.@threads for vn in 1:C.n
        iter = 0
        not_found = true
        to_check = [(0, [vn])]
        while not_found
            iter += 1
            to_check_next = Vector{Tuple{Int, Vector{Int}}}()
            for i in 1:length(to_check)
                for (prev, v_arr) in (to_check[i], )
                    for v in v_arr
                        for cn in var_adj_list[v]
                            if cn != prev
                                if iter != 1 && vn ∈ check_adj_list[cn]
                                    not_found = false
                                    girth_arr[vn] = iter + 1
                                    break
                                else
                                    push!(to_check_next, (cn, [v2 for v2 in check_adj_list[cn] if v2 != v]))
                                end
                            end
                        end
                        !not_found && break
                    end
                    !not_found && break
                end
                !not_found && break
            end
            !not_found && break
            iter += 1
            iter > max_iter && error("Hit the maximum number of iterations")
            isempty(to_check_next) && break
            to_check = to_check_next
        end
    end
    # println(girth_arr)
    min = minimum(girth_arr)
    iseven(min) || error("Computed girth to be an odd integer")
    min == 0 ? (C.girth = -1;) : (C.girth = min;)
    return C.girth
end

"""
    computation_graph(C::AbstractLDPCCode, lvl::Int, v::Int, v_type::Symbol = :v)

Return a figure representing the expansion of the Tanner graph of `C` to level `lvl`
for node `v`. If `v_type` is `:v`, `v` is interpreted as a variable node; otherwise,
`v_type` is `:c` and `v` is interpreted as a check node.

# Note
- Run `using Makie` to activate this extension.
"""
function computation_graph end

mutable struct _ACEVarNode
    id::Int
    parent_id::Int
    lvl::Int
    cum_ACE::Int
    local_ACE::Int
end

mutable struct _ACECheckNode
    id::Int
    parent_id::Int
    lvl::Int
    cum_ACE::Int
end

#############################
      # simple cycles
#############################

function _enumerate_cycles(L::AbstractLDPCCode, len::Int)
    g, _, _ = Tanner_graph(L)
    d = DiGraph(g)
    if len == -1
        cycles = simplecycles_hawick_james(d)
    else
        cycles = simplecycles_limited_length(d, len)
    end

    unique_cycles = Set{Vector{Int}}()
    final_cycles = Vector{Vector{Int}}()
    for cycle in cycles
        if length(cycle) > 2
            temp = sort(cycle)
            found_flag = false
            for key in keys(unique_cycles.dict)
                if temp == key
                    found_flag = true
                    break
                end
            end
            if !found_flag
                push!(unique_cycles, temp)
                push!(final_cycles, cycle)
            end
        end
    end
    
    (len == -1) && (L.simple_cycles = final_cycles;)
    lens = length.(final_cycles)
    girth = minimum(lens)
    girth == 9999999 && (girth = -1)
    if ismissing(L.girth)
        L.girth = girth
    else
        if L.girth != girth
            @warn "Known girth, $(L.girth), does not match just computed girth, $girth"
        end
    end
    return final_cycles
end

"""
    enumerate_simple_cycles(L::AbstractLDPCCode; len::Int = -1)

Return the unique simple cycles up to length `len` of the Tanner graph of `L`. If
`len` is `-1`, then all simple cycles will be enumerated. An empty `Vector{Vector{Int}}`
is returned when there is no cycles.

# Note
- Simple cycles do not contain the same vertex twice.
- Cycles are returned as a vector of vertex indices, where the vertices are ordered left-to-right by
  columns of `parity_check_matrix(L)` then top-to-bottom by rows.
"""
function enumerate_simple_cycles(L::AbstractLDPCCode; len::Int = -1)
    (len == -1 || ispositive(len)) ||
        throw(DomainError("Cycle length parameter must be `-1` or a positive integer"))

    if len == -1 && !isempty(L.simple_cycles)
        return L.simple_cycles
    elseif len != -1 && !isempty(L.simple_cycles)
        return filter(x -> length(x) <= len, L.simple_cycles)
    else
        return _enumerate_cycles(L, len)
    end
end

"""
    simple_cycle_length_distribution(L::AbstractLDPCCode; len::Int = -1)

Return a dictionary of (length, count) pairs for the unique simple cycles up to length `len` of the
Tanner graph of `L`. If `len` is `-1`, then all simple cycles will be enumerated. An empty
dictionary is returned when there are no cycles.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
simple_cycle_length_distribution(L::AbstractLDPCCode; len::Int = -1) = StatsBase.countmap(length.(
    enumerate_simple_cycles(L, len = len)))

"""
    simple_cycle_length_distribution_plot(L::AbstractLDPCCode; len::Int = -1)

Return a bar graph and dictionary of (length, count) pairs for the unique simple cycles up to
length `len` of the Tanner graph of `L`. If `len` is `-1`, then all simple cycles will be
enumerated. An empty figure and dictionary are returned when there are no cycles.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
- Run `using Makie` to activate this extension.
"""
function simple_cycle_length_distribution_plot end

"""
    average_simple_cycle_length(L::AbstractLDPCCode; len::Int = -1)

Return the average cycle length of unique simple cycles up to length `len` of the Tanner graph of
`L`. If `len` is `-1`, then all simple cycles will be enumerated.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
average_simple_cycle_length(L::AbstractLDPCCode; len::Int = -1) = mean(length.(
    enumerate_simple_cycles(L, len = len)))

"""
    median_simple_cycle_length(L::AbstractLDPCCode; len::Int = -1)

Return the median cycle length of unique simple cycles up to length `len` of the Tanner graph of
`L`. If `len` is `-1`, then all simple cycles will be enumerated.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
median_simple_cycle_length(L::AbstractLDPCCode; len::Int = -1) = median(length.(
    enumerate_simple_cycles(L, len = len)))

"""
    mode_simple_cycle_length(L::AbstractLDPCCode; len::Int = -1)

Return the most common cycle length of unique simple cycles up to length `len` of the Tanner graph
of `L`. If `len` is `-1`, then all simple cycles will be enumerated.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
mode_simple_cycle_length(L::AbstractLDPCCode; len::Int = -1) = StatsBase.mode(length.(
    enumerate_simple_cycles(L, len = len)))

# TODO: there are ways to compute this without enumerating all of the cycles, but they are complicated and perhaps not worth the effort unless explicitly requested by a user
"""
    count_simple_cycles(L::AbstractLDPCCode; len::Int = -1)

Return the total number of unique simple cycles up to length `len` of the Tanner graph
of `L`. If `len` is `-1`, then all simple cycles will be enumerated.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
count_simple_cycles(L::AbstractLDPCCode; len::Int = -1) = length(enumerate_simple_cycles(L,
    len = len))

"""
    simple_cycle_distribution_by_variable_node(L::AbstractLDPCCode; len::Int = -1)

Return a dictionary of (node, count) pairs for the unique simple cycles up to length `len` of the
Tanner graph of `L`. If `len` is `-1`, then all simple cycles will be enumerated. An empty
dictionary is returned when there are no cycles.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
simple_cycle_distribution_by_variable_node(L::AbstractLDPCCode; len::Int = -1) =
    StatsBase.countmap(reduce(vcat, enumerate_simple_cycles(L, len = len)))

# TODO: write this function in MakieExt
"""
    simple_cycle_distribution_by_variable_node_plot(L::AbstractLDPCCode; len::Int = -1)

Return bar graph and a dictionary of (node, count) pairs for the unique simple cycles up to length
`len` of the Tanner graph of `L`. If `len` is `-1`, then all simple cycles will be enumerated. An
empty figure and dictionary are returned when there are no cycles.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
    already cached.
- Run `using Makie` to activate this extension.
"""
function simple_cycle_distribution_by_variable_node_plot end

#############################
        # short cycles
#############################

# TODO: time this versus calling girth and then doing the upper bound function
"""
    enumerate_short_cycles(L::AbstractLDPCCode; len::Int = -1)

Return the unique short cycles up to length `len` of the Tanner graph of `L`. If
`len` is `-1`, then all short cycles will be enumerated. An empty `Vector{Vector{Int}}`
is returned when there is no short cycles.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- Cycles are returned as a vector of vertex indices, where the vertices are ordered left-to-right by
  columns of `parity_check_matrix(L)` then top-to-bottom by rows.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
enumerate_short_cycles(L::AbstractLDPCCode; len::Int = -1) = filter(x -> length(x) <=
    2 * girth(L) - 2, _enumerate_cycles(L, len))

"""
    short_cycle_length_distribution(L::AbstractLDPCCode; len::Int = -1)

Return a dictionary of (length, count) pairs for the unique short cycles up to length `len` of the
Tanner graph of `L`. If `len` is `-1`, then all short cycles will be enumerated. An empty
dictionary is returned when there are no cycles.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
short_cycle_length_distribution(L::AbstractLDPCCode; len::Int = -1) = StatsBase.countmap(length.(
    enumerate_short_cycles(L, len = len)))

"""
    short_cycle_length_distribution_plot(L::AbstractLDPCCode; len::Int = -1)

Return a bar graph and dictionary of (length, count) pairs for the unique short cycles up to
length `len` of the Tanner graph of `L`. If `len` is `-1`, then all short cycles will be
enumerated. An empty figure and dictionary are returned when there are no cycles.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
- Run `using Makie` to activate this extension.
"""
function short_cycle_length_distribution_plot end

# TODO: axis and move to MakieExt
function short_cycle_length_distribution_plot(L::AbstractLDPCCode; len::Int = -1)
    dist = short_cycle_length_distribution(L)
    x_data = collect(keys(dist))
    y_data = collect(values(dist))
    fig = Figure()
    ax = Axis(fig[1, 2])
    barplot!(ax, x_data, y_data, bar_width = 1, x_ticks = collect(girth(L):2 * girth(L) - 2),
        y_ticks = y_data)
    display(fig)
    return fig
end

"""
    average_short_cycle_length(L::AbstractLDPCCode; len::Int = -1)

Return the average cycle length of unique short cycles up to length `len` of the Tanner graph of
`L`. If `len` is `-1`, then all short cycles will be enumerated.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
average_short_cycle_length(L::AbstractLDPCCode; len::Int = -1) = mean(length.(
    enumerate_short_cycles(L, len = len)))

"""
    median_short_cycle_length(L::AbstractLDPCCode; len::Int = -1)

Return the median cycle length of unique short cycles up to length `len` of the Tanner graph of
`L`. If `len` is `-1`, then all short cycles will be enumerated.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
median_short_cycle_length(L::AbstractLDPCCode; len::Int = -1) = median(length.(
    enumerate_short_cycles(L, len = len)))

"""
    mode_short_cycle_length(L::AbstractLDPCCode; len::Int = -1)

Return the most common cycle length of unique short cycles up to length `len` of the Tanner graph
of `L`. If `len` is `-1`, then all short cycles will be enumerated.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
mode_short_cycle_length(L::AbstractLDPCCode; len::Int = -1) = StatsBase.mode(length.(
    enumerate_short_cycles(L, len = len)))

 """
    count_short_cycles(L::AbstractLDPCCode; len::Int = -1)

Return the total number of unique short cycles up to length `len` of the Tanner graph
of `L`. If `len` is `-1`, then all short cycles will be enumerated.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
count_short_cycles(L::AbstractLDPCCode; len::Int = -1) = length(enumerate_short_cycles(L,
    len = len))

"""
    short_cycle_distribution_by_variable_node(L::AbstractLDPCCode; len::Int = -1)

Return a dictionary of (node, count) pairs for the unique short cycles up to length `len` of the
Tanner graph of `L`. If `len` is `-1`, then all short cycles will be enumerated. An empty
dictionary is returned when there are no cycles.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
"""
short_cycle_distribution_by_variable_node(L::AbstractLDPCCode; len::Int = -1) =
    StatsBase.countmap(reduce(vcat, enumerate_short_cycles(L, len = len)))

# TODO: write this function in MakieExt
"""
    short_cycle_distribution_by_variable_node_plot(L::AbstractLDPCCode; len::Int = -1)

Return bar graph and a dictionary of (node, count) pairs for the unique short cycles up to length
`len` of the Tanner graph of `L`. If `len` is `-1`, then all short cycles will be enumerated. An
empty figure and dictionary are returned when there are no cycles.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
    already cached.
- Run `using Makie` to activate this extension.
"""
function short_cycle_distribution_by_variable_node_plot end

#############################
        # lollipops
#############################








#############################
           # ACE
#############################

# TODO: degree 1 nodes
# why did I make this note? is ACE defined for them differently?
"""
    shortest_cycle_ACE(C::AbstractLDPCCode, v::Int)
    shortest_cycle_ACE(C::AbstractLDPCCode, vs::Vector{Int})
    shortest_cycle_ACE(C::AbstractLDPCCode)

Return a cycle of minimum length and minimum ACE in the Tanner graph of `C`
for the vertex `v` or vertices `vs`, in the order (ACEs, cycles). If no vertices
are given, all vertices are computed by default. The cycle `v1 -- c1 -- ... -- 
cn -- vn` is returned in the format `[(v1, c1), (c1, v2), ..., (cn, vn)]`.
"""
function shortest_cycle_ACE(C::AbstractLDPCCode, vs::Vector{Int})
    isempty(vs) && throw(ArgumentError("Input variable node list cannot be empty"))
    all(x -> 1 <= x <= C.n, vs) || throw(DomainError("Variable node indices must be between 1 and length(C)"))

    # might not be efficient to have the or here
    vs_to_do = [x for x in vs if isempty(C.ACEs_per_var_node[x]) || isempty(C.shortest_cycles[x])]
    processed = false
    if !isempty(vs_to_do)
        processed = true
        check_adj_list, var_adj_list = _node_adjacencies(C.H)
        
        Threads.@threads for i in 1:length(vs_to_do)
            # moving this inside allocates more but allows for multi-threading
            check_nodes = [_ACECheckNode(j, -1, -1, -1) for j in 1:length(check_adj_list)]
            var_nodes = [_ACEVarNode(j, -1, -1, -1, length(var_adj_list[j]) - 2) for j in 1:C.n]

            ACEs = Vector{Int}()
            cycle_lens = Vector{Int}()
            cycles = Vector{Vector{Tuple{Int, Int}}}()
            not_emptied = true
            root = var_nodes[vs_to_do[i]]
            root.lvl = 0
            root.cum_ACE = root.local_ACE
            queue = Deque{Union{_ACECheckNode, _ACEVarNode}}()
            push!(queue, root)
            while length(queue) > 0
                curr = first(queue)
                if isa(curr, _ACEVarNode)
                    for cn in var_adj_list[curr.id]
                        # can't pass messages back to the same node
                        if cn != curr.parent_id
                            cn_node = check_nodes[cn]
                            if cn_node.lvl != -1
                                # have seen before
                                push!(ACEs, curr.cum_ACE + cn_node.cum_ACE - root.local_ACE)
                                push!(cycle_lens, curr.lvl + cn_node.lvl + 1)

                                # trace the cycle from curr to root and cn_node to root
                                temp = Vector{Tuple{Int, Int}}()
                                node = cn_node
                                while node.lvl != 0
                                    push!(temp, (node.parent_id, node.id))
                                    if isodd(node.lvl)
                                        node = var_nodes[node.parent_id]
                                    else
                                        node = check_nodes[node.parent_id]
                                    end
                                end
                                reverse!(temp)
                                push!(temp, (cn_node.id, curr.id))
                                node = curr
                                while node.lvl != 0
                                    push!(temp, (node.id, node.parent_id))
                                    if isodd(node.lvl)
                                        node = var_nodes[node.parent_id]
                                    else
                                        node = check_nodes[node.parent_id]
                                    end
                                end
                                push!(cycles, temp)

                                # finish this level off but don't go deeper so remove children at lower level
                                if not_emptied
                                    while length(queue) > 0
                                        back = last(queue)
                                        if back.lvl != curr.lvl
                                            pop!(queue)
                                        else
                                            break
                                        end
                                    end
                                    not_emptied = false
                                end
                            elseif not_emptied
                                cn_node.lvl = curr.lvl + 1
                                cn_node.parent_id = curr.id
                                cn_node.cum_ACE = curr.cum_ACE
                                push!(queue, cn_node)
                            end
                        end
                    end
                else
                    for vn in check_adj_list[curr.id]
                        # can't pass messages back to the same node
                        if vn != curr.parent_id
                            vn_node = var_nodes[vn]
                            if vn_node.lvl != -1
                                # have seen before
                                push!(ACEs, curr.cum_ACE + vn_node.cum_ACE - root.local_ACE)
                                push!(cycle_lens, curr.lvl + vn_node.lvl + 1)

                                # trace the cycle from curr to root and cn_node to root
                                temp = Vector{Tuple{Int, Int}}()
                                node = vn_node
                                while node.lvl != 0
                                    push!(temp, (node.parent_id, node.id))
                                    if isodd(node.lvl)
                                        node = var_nodes[node.parent_id]
                                    else
                                        node = check_nodes[node.parent_id]
                                    end
                                end
                                reverse!(temp)
                                push!(temp, (vn_node.id, curr.id))
                                node = curr
                                while node.lvl != 0
                                    push!(temp, (node.id, node.parent_id))
                                    if isodd(node.lvl)
                                        node = var_nodes[node.parent_id]
                                    else
                                        node = check_nodes[node.parent_id]
                                    end
                                end
                                push!(cycles, temp)

                                # finish this level off but don't go deeper so remove children at lower level
                                if not_emptied
                                    while length(queue) > 0
                                        back = last(queue)
                                        if back.lvl != curr.lvl
                                            pop!(queue)
                                        else
                                            break
                                        end
                                    end
                                    not_emptied = false
                                end
                            elseif not_emptied
                                vn_node.lvl = curr.lvl + 1
                                vn_node.parent_id = curr.id
                                vn_node.cum_ACE = curr.cum_ACE + vn_node.local_ACE
                                push!(queue, vn_node)
                            end
                        end
                    end
                end
                popfirst!(queue)
            end
            # println("variable node $i, cycles: $cycles")
            C.ACEs_per_var_node[vs_to_do[i]] = ACEs
            C.cycle_lens[vs_to_do[i]] = cycle_lens
            C.shortest_cycles[vs_to_do[i]] = cycles
        end
    end

    vs_ACE = zeros(Int, length(vs))
    cycles_vs = [Vector{Tuple{Int, Int}}() for _ in 1:length(vs)]
    for i in 1:length(vs)
        min, index = findmin(C.ACEs_per_var_node[vs[i]])
        vs_ACE[i] = min
        cycles_vs[i] = C.shortest_cycles[vs[i]][index]
    end

    if processed
        if all(!isempty, C.cycle_lens)
            girth = minimum([minimum(C.cycle_lens[i]) for i in 1:C.n])
            if ismissing(C.girth)
                C.girth = girth
            else
                if C.girth != girth
                    @warn "Known girth, $(C.girth), does not match just computed girth, $girth"
                end
            end
        end
    end

    return vs_ACE, cycles_vs
end
shortest_cycle_ACE(C::AbstractLDPCCode, v::Int) = shortest_cycle_ACE(C, [v])[1]
shortest_cycle_ACE(C::AbstractLDPCCode) = shortest_cycle_ACE(C, collect(1:C.n))

"""
    shortest_cycles(C::AbstractLDPCCode, v::Int)
    shortest_cycles(C::AbstractLDPCCode, vs::Vector{Int})
    shortest_cycles(C::AbstractLDPCCode)

Return all the cycles of shortest length in the Tanner graph of `C` for the vertex `v` or
vertices `vs`. If no vertices are given, all vertices are computed by default.

# Note
- The length of the shortest cycle is not necessarily the same for each vertex.
- To reduce computational complexity, the same cycle may appear under each vertex in the cycle.
"""
function shortest_cycles(C::AbstractLDPCCode, vs::Vector{Int})
    shortest_cycle_ACE(C, vs)
    return C.shortest_cycles[vs]
end
shortest_cycles(C::AbstractLDPCCode, v::Int) = shortest_cycles(C, [v])[1]
shortest_cycles(C::AbstractLDPCCode) = shortest_cycles(C, collect(1:C.n))

"""
    ACE_distribution(C::AbstractLDPCCode, v::Int)
    ACE_distribution(C::AbstractLDPCCode, vs::Vector{Int})
    ACE_distribution(C::AbstractLDPCCode)

Return the ACEs and cycle lengths for vertex `v` or vertices `vs` of the Tanner graph
of `C`. If no vertices are given, all vertices are computed by default.
"""
function ACE_distribution(C::AbstractLDPCCode, vs::Vector{Int})
    # using the original DFS approach constructs a significantly larger tree than this truncated BFS approach

    isempty(vs) && throw(ArgumentError("Input node list cannot be empty"))
    all(x -> 1 <= x <= C.n, vs) || throw(DomainError("Variable node index must be between 1 and length(C)"))

    vs_to_do = [x for x in vs if isempty(C.ACEs_per_var_node[x])]
    processed = false
    if !isempty(vs_to_do)
        processed = true
        check_adj_list, var_adj_list = _node_adjacencies(C.H)
    
        Threads.@threads for i in 1:length(vs_to_do)
            # moving this inside allocates more but allows for multi-threading
            check_nodes = [_ACECheckNode(i, -1, -1, -1) for i in 1:length(check_adj_list)]
            var_nodes = [_ACEVarNode(i, -1, -1, -1, length(var_adj_list[i]) - 2) for i in 1:C.n]

            ACEs = Vector{Int}()
            cycle_lens = Vector{Int}()
            root = var_nodes[vs[i]]
            root.lvl = 0
            root.cum_ACE = root.local_ACE
            queue = Queue{Union{_ACECheckNode, _ACEVarNode}}()
            enqueue!(queue, root)
            while length(queue) > 0
                curr = first(queue)
                if isa(curr, _ACEVarNode)
                    for cn in var_adj_list[curr.id]
                        # can't pass messages back to the same node
                        if cn != curr.parent_id
                            cn_node = check_nodes[cn]
                            if cn_node.lvl != -1
                                # have seen before
                                push!(ACEs, curr.cum_ACE + cn_node.cum_ACE - root.local_ACE)
                                push!(cycle_lens, curr.lvl + cn_node.lvl + 1)
                            else
                                cn_node.lvl = curr.lvl + 1
                                cn_node.parent_id = curr.id
                                cn_node.cum_ACE = curr.cum_ACE
                                enqueue!(queue, cn_node)
                            end
                        end
                    end
                else
                    for vn in check_adj_list[curr.id]
                        # can't pass messages back to the same node
                        if vn != curr.parent_id
                            vn_node = var_nodes[vn]
                            if vn_node.lvl != -1
                                # have seen before
                                push!(ACEs, curr.cum_ACE + vn_node.cum_ACE - root.local_ACE)
                                push!(cycle_lens, curr.lvl + vn_node.lvl + 1)
                            else
                                vn_node.lvl = curr.lvl + 1
                                vn_node.parent_id = curr.id
                                vn_node.cum_ACE = curr.cum_ACE + vn_node.local_ACE
                                enqueue!(queue, vn_node)
                            end
                        end
                    end
                end
                dequeue!(queue)
            end
            C.ACEs_per_var_node[vs_to_do[i]] = ACEs
            C.cycle_lens[vs_to_do[i]] = cycle_lens
        end
    end

    vs_ACEs = [C.ACEs_per_var_node[i] for i in vs]
    lengths = [C.cycle_lens[i] for i in vs]

    if processed
        if all(!isempty, C.cycle_lens)
            girth = minimum([minimum(C.cycle_lens[i]) for i in 1:C.n])
            if ismissing(C.girth)
                C.girth = girth
            else
                if C.girth != girth
                    @warn "Known girth, $(C.girth), does not match just computed girth, $girth"
                end
            end
        end
    end

    return vs_ACEs, lengths
end
function ACE_distribution(C::AbstractLDPCCode, v::Int)
    vs_ACE, lengths = ACE_distribution(C, [v])
    return vs_ACE[1], lengths[1]
end

# TODO: plots
ACE_distribution(C::AbstractLDPCCode) = ACE_distribution(C, collect(1:C.n))

"""
    average_ACE_distribution(C::AbstractLDPCCode, v::Int)
    average_ACE_distribution(C::AbstractLDPCCode, vs::Vector{Int})
    average_ACE_distribution(C::AbstractLDPCCode)

Return the average ACE of the vertex `v` or vertices `vs` of the Tanner graph of `C`. If no
vertices are given, all vertices are computed (individually) by default.
"""
function average_ACE_distribution(C::AbstractLDPCCode, vs::Vector{Int})
    vs_to_do = [x for x in vs if isempty(C.ACEs_per_var_node[x])]
    isempty(vs_to_do) || ACE_distribution(C, vs_to_do)
    return [mean(C.ACEs_per_var_node[v]) for v in vs]
end
average_ACE_distribution(C::AbstractLDPCCode, v::Int) = average_ACE_distribution(C, [v])[1]
average_ACE_distribution(C::AbstractLDPCCode) = average_ACE_distribution(C, collect(1:C.n))

"""
    median_ACE_distribution(C::AbstractLDPCCode, v::Int)
    median_ACE_distribution(C::AbstractLDPCCode, vs::Vector{Int})
    median_ACE_distribution(C::AbstractLDPCCode)

Return the median ACE of the vertex `v` or vertices `vs` of the Tanner graph of `C`. If no
vertices are given, all vertices are computed (individually) by default.
"""
function median_ACE_distribution(C::AbstractLDPCCode, vs::Vector{Int})
    vs_to_do = [x for x in vs if isempty(C.ACEs_per_var_node[x])]
    isempty(vs_to_do) || ACE_distribution(C, vs_to_do)
    return [median(C.ACEs_per_var_node[v]) for v in vs]
end
median_ACE_distribution(C::AbstractLDPCCode, v::Int) = median_ACE_distribution(C, [v])[1]
median_ACE_distribution(C::AbstractLDPCCode) = median_ACE_distribution(C, collect(1:C.n))

"""
    mode_ACE_distribution(C::AbstractLDPCCode, v::Int)
    mode_ACE_distribution(C::AbstractLDPCCode, vs::Vector{Int})
    mode_ACE_distribution(C::AbstractLDPCCode)

Return the mode ACE of the vertex `v` or vertices `vs` of the Tanner graph of `C`. If no
vertices are given, all vertices are computed (individually) by default.

# Note
- In case of ties, the smallest tied value is returned.
"""
function mode_ACE_distribution(C::AbstractLDPCCode, vs::Vector{Int})
    vs_to_do = [x for x in vs if isempty(C.ACEs_per_var_node[x])]
    isempty(vs_to_do) || ACE_distribution(C, vs_to_do)
    return [StatsBase.mode(sort(C.ACEs_per_var_node[v])) for v in vs]
end
mode_ACE_distribution(C::AbstractLDPCCode, v::Int) = mode_ACE_distribution(C, [v])[1]
mode_ACE_distribution(C::AbstractLDPCCode) = mode_ACE_distribution(C, collect(1:C.n))

"""
    ACE_spectrum(C::AbstractLDPCCode)

Return the ACE spectrum of the Tanner graph of `C`.
"""
function ACE_spectrum(C::AbstractLDPCCode)
    vs_ACEs, lengths = ACE_distribution(C, collect(1:C.n))
    # (false) spectrum: how many nodes have that ACE for that length
    # (true) spectrum: for a given length 4 <= l <= maximum(variabledegreedistribution(C)),
    # how many var nodes have shortest cycle that ACE

    shortest_lens = [minimum(i) for i in lengths]
    girth = minimum(shortest_lens)
    if ismissing(C.girth)
        C.girth = girth
    else
        if C.girth != girth
            @warn "Known girth, $(C.girth), does not match just computed girth, $girth"
        end
    end

    counts = [Dict{Int, Int}() for _ in 1:length(girth:2:2 * girth - 2)]
    for (k, l) in enumerate(girth:2:2 * girth - 2)
        for i in 1:length(shortest_lens)
            if shortest_lens[i] == l
                for j in 1:length(lengths[i])
                    if lengths[i][j] == l
                        if vs_ACEs[i][j] ∈ keys(counts[k])
                            counts[k][vs_ACEs[i][j]] += 1
                        else
                            counts[k][vs_ACEs[i][j]] = 1
                        end
                    end
                end
            end
        end
    end
    return counts
end

"""
    ACE_spectrum_plot(C::AbstractLDPCCode)

Return an interactive figure and data for the ACE spectrum of the Tanner graph of `C`.

# Note
- Run `using Makie` to activate this extension.
"""
function ACE_spectrum_plot end
