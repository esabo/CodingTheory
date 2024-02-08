# Copyright (c) 2023, 2024 Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

function _copying_Hastings(H_X::CTMatrixTypes, H_Z::CTMatrixTypes)
    nr, n = size(H_X)

    # get column weight
    q_X = maximum(count(!iszero, H_X[:, j]) for j in 1:n)

    # X stabilizers
    X = zero_matrix(base_ring(H_X), nr + (q_X - 1) * n, q_X * n)

    # copied X stabilizers
    for j in 1:n
        for (k, i) in enumerate(getindex.(findall(isone, H_X[:, j]), 1))
            X[i, q_X * (j - 1) + k] = H_X[i, j]
        end
    end

    # new X stabilizers
    for j in 1:n
        for i in 1:q_X - 1
            row = nr + i + (q_X - 1) * (j - 1)
            col = q_X * (j - 1) + i
            X[row, col] = 1
            X[row, col + 1] = 1
        end
    end

    # Z stabilizers
    nr = nrows(H_Z)
    Z = zero_matrix(base_ring(H_Z), nr, q_X * n)
    for i in 1:nr
        for j in 1:n
            for k in 1:q_X
                Z[i, (j - 1) * q_X + k] = H_Z[i, j]
            end
        end
    end
    return X, Z
end

function _copying_reduced(H_X::CTMatrixTypes, H_Z::CTMatrixTypes)
    n_X, n = size(H_X)

    # get column weights
    qubit_expansion = [count(!iszero, H_X[:, j]) for j in 1:n]
    q_X = maximum(qubit_expansion)

    # X stabilizers
    X = zero_matrix(base_ring(H_X), n_X + sum(qubit_expansion .- 1), sum(qubit_expansion))

    # keep track of how many edges (in the tanner graph) remain to be used for each qubit
    qubits_available = zeros(Int, q_X, n)
    for (j, n2) in enumerate(qubit_expansion)
        qubits_available[1:n2, j] .= 1
    end

    # copy in the old X stabilizers
    for i in 1:n_X
        for j in 1:n
            iszero(H_X[i, j]) && continue
            k = findfirst(!iszero, qubits_available[:, j])
            qubits_available[k, j] -= 1
            X[i, sum(qubit_expansion[1:j - 1]) + k] = H_X[i, j]
        end
    end

    # new X stabilizers
    row = n_X + 1
    for (i, n) in enumerate(qubit_expansion)
        for j in 1:n - 1
            a = sum(qubit_expansion[1:i - 1]) + j
            X[row, a] = 1
            X[row, a + 1] = 1
            row += 1
        end
    end

    # Z stabilizers
    n_Z = nrows(H_Z)
    Z = zero_matrix(base_ring(H_Z), n_Z, sum(qubit_expansion))
    for i in 1:n_Z
        for j in 1:n
            for k in sum(qubit_expansion[1:j - 1]) + 1:sum(qubit_expansion[1:j])
                Z[i, k] = H_Z[i, j]
            end
        end
    end

    return X, Z
end

function _copying_target(H_X::CTMatrixTypes, H_Z::CTMatrixTypes, target_q_X::Int = 3)
    target_q_X >= 3 || throw(DomainError(target_q_X, "Target column weight must be at least 3."))

    n_X, n = size(H_X)
    # get column weights
    q_Xs = [count(!iszero, H_X[:, j]) for j in 1:n]

    # repetition expansion of qubits
    qubit_expansion = [q_X <= target_q_X ? 1 : 2 + cld(q_X - 2 * (target_q_X - 1), target_q_X - 2) for q_X in q_Xs]

    # X stabilizers
    X = zero_matrix(base_ring(H_X), n_X + sum(qubit_expansion .- 1), sum(qubit_expansion))

    # keep track of how many edges (in the tanner graph) remain to be used for each qubit
    qubits_available = zeros(Int, maximum(qubit_expansion), n)
    for (j, n2) in enumerate(qubit_expansion)
        if n2 == 1
            qubits_available[1, j] = target_q_X
            continue
        end
        qubits_available[unique([1, n2]), j] .= target_q_X - 1
        if n2 > 2
            qubits_available[2:n2 - 1, j] .= target_q_X - 2
        end
    end

    # copy in the old X stabilizers
    for i in 1:n_X
        for j in 1:n
            iszero(H_X[i, j]) && continue
            k = findfirst(!iszero, qubits_available[:, j])
            qubits_available[k, j] -= 1
            X[i, sum(qubit_expansion[1:j - 1]) + k] = H_X[i, j]
        end
    end

    # new X stabilizers
    row = n_X + 1
    for (i, n) in enumerate(qubit_expansion)
        for j in 1:n - 1
            a = sum(qubit_expansion[1:i - 1]) + j
            X[row, a] = 1
            X[row, a + 1] = 1
            row += 1
        end
    end

    # Z stabilizers
    n_Z = nrows(H_Z)
    Z = zero_matrix(base_ring(H_Z), n_Z, sum(qubit_expansion))
    for i in 1:n_Z
        for j in 1:n
            for k in sum(qubit_expansion[1:j - 1]) + 1:sum(qubit_expansion[1:j])
                Z[i, k] = H_Z[i, j]
            end
        end
    end
    return X, Z
end

"""
    copying(H_X::CTMatrixTypes, H_Z::CTMatrixTypes; method::Symbol = :Hastings, target_q_X::Int = 3)

Return the result of copying on `H_X` and `H_Z` using either the Hastings, reduced, or targeted
methods.
"""
function copying(H_X::CTMatrixTypes, H_Z::CTMatrixTypes; method::Symbol = :Hastings,
    target_q_X::Int = 3)

    method ∈ (:Hastings, :reduced, :target) || throw(ArgumentError("Unknown method type"))
    target_q_X >= 3 || throw(DomainError(target_q_X, "Target must be at least 3"))
    # should we check these commute or trust the user?

    if method == :Hastings
       return _copying_Hastings(H_X, H_Z)
    elseif method == :reduced
        return _copying_reduced(H_X, H_Z)
    else
        return _copying_target(H_X, H_Z, target_q_X)
    end
end

"""
    copying(S::AbstractStabilizerCode, method::Symbol = :Hastings, target_q_X::Int = 3)

Return the result of copying on `S` using either the Hastings, reduced, or targeted methods.
"""
copying(S::T; method::Symbol = :Hastings, target_q_X::Int = 3) where {T <: AbstractStabilizerCode} =
    copying(CSSTrait(T), S, method, target_q_X)
function copying(::IsCSS, S::AbstractStabilizerCode, method::Symbol, target_q_X::Int)
    method ∈ (:Hastings, :reduced, :target) || throw(ArgumentError("Unknown method type"))
    target_q_X >= 3 || throw(DomainError(target_q_X, "Target must be at least 3"))

    H_X, H_Z = copying(S.X_stabs, S.Z_stabs, method = method, target_q_X = target_q_X)
    return CSSCode(H_X, H_Z)
end
copying(::IsNotCSS, S::AbstractStabilizerCode, method::Symbol, target_q_X::Int) =
    error("Only valid for CSS codes.")

"""
    gauging(H_X::CTMatrixTypes, H_Z::CTMatrixTypes)

Return the result of gauging on `H_X` and `H_Z`.
"""
function gauging(H_X::CTMatrixTypes, H_Z::CTMatrixTypes)
    # get row weights and number of new X stabilizers
    n_X, n = size(H_X)
    n_Z = nrows(H_Z)
    weights = zeros(Int, n_X)
    num_X = 0
    num_new_qubits = 0
    for i in 1:n_X
        weights[i] = count(!iszero, H_X[i, :])
        num_X += max(1, weights[i] - 2)
        num_new_qubits += max(0, weights[i] - 3)
    end
    iszero(num_new_qubits) && (return H_X, H_Z;)

    # copy in the original Z stabilizers, new qubits will be adjusted below
    Z = hcat(H_Z, zero_matrix(base_ring(H_X), n_Z, num_new_qubits))

    # new X stabilizers
    X = zero_matrix(base_ring(H_X), num_X, n + num_new_qubits)

    i_new = 1
    new_qubit_index = 1 + n
    for (i, wt) in enumerate(weights)
        if wt < 4
            # just copy in the old stabilizer
            for j in 1:n
                X[i_new, j] = H_X[i, j]
            end
            i_new += 1
        else
            new_qubits = new_qubit_index:new_qubit_index + wt - 3

            # cut the X stabilizer up in to weight 3 stabilizers
            for (l, k) in enumerate(new_qubits)
                nonzeros = getindex.(findall(!iszero, H_X[i, :]), 2)
                if k == new_qubit_index
                    X[i_new, nonzeros[1]] = H_X[i, nonzeros[1]]
                    X[i_new, nonzeros[2]] = H_X[i, nonzeros[2]]
                    X[i_new, k] = 1
                elseif k == wt - 3 + new_qubit_index
                    X[i_new, nonzeros[end - 1]] = H_X[i, nonzeros[end - 1]]
                    X[i_new, nonzeros[end]] = H_X[i, nonzeros[end]]
                    X[i_new, k - 1] = 1
                else
                    X[i_new, nonzeros[l + 1]] = H_X[i, nonzeros[l + 1]]
                    X[i_new, k - 1] = 1
                    X[i_new, k] = 1
                end
                i_new += 1
            end

            # adjust the new qubits of the Z stabilizers so that they commute with the new X stabilizers
            for j in 1:n_Z
                Z_qubits = getindex.(findall(!iszero, Z[j, :]), 2)
                X_qubits = getindex.(findall(!iszero, H_X[i, :]), 2)
                for m in 1:wt - 3
                    if isodd(length(Z_qubits ∩ X_qubits[1:m + 1]))
                        Z[j, new_qubits[m]] = 1
                    end
                end
            end

            new_qubit_index += wt - 3
        end
    end
    return X, Z
end

"""
    gauging(S::AbstractStabilizerCode)

Return the result of gauging on `S`.
"""
gauging(S::T) where {T <: AbstractStabilizerCode} = gauging(CSSTrait(T), S)
gauging(::IsCSS, S::AbstractStabilizerCode) = CSSCode(gauging(S.X_stabs, S.Z_stabs)...)
gauging(::IsNotCSS, S::AbstractStabilizerCode) = error("Only valid for CSS codes.")

# should just call distance balancing and keep with that function
# # this is thickening on its own
# """
#     thickening(S::AbstractStabilizerCode, l::Integer)

# Return the result of thickening on `S`.
# """
# function thickening(S::AbstractStabilizerCode, l::Integer)
#     # C is just the repetition code, but thickening requires the parity check matrix to be in a particular form
#     H = matrix(GF(2), diagm(l - 1, l, 0 => ones(Int, l - 1), 1 => ones(Int, l - 1)))
#     # TODO: make a version of this function to take in two matrices
#     return distance_balancing(S, LinearCode(H, true))
# end

# for weight reduction, it's easier to choose heights with thickening
"""
    thickening_and_choose_heights(H_X::CTMatrixTypes, H_Z::CTMatrixTypes, l::Integer, heights::Vector{Int})

Return the result of thickening and choosing heights on `H_X` and `H_Z`.
"""
function thickening_and_choose_heights(H_X::CTMatrixTypes, H_Z::CTMatrixTypes, l::Integer, heights::Vector{Int})

    F = base_ring(H_X)
    n_Z = nrows(H_Z)
    n_X, n = size(H_X)
    @assert length(heights) == n_Z
    @assert all(1 <= x <= l for x in heights)
    H = matrix(F, diagm(l - 1, l, 0 => ones(Int, l - 1), 1 => ones(Int, l - 1)))
    X = hcat(identity_matrix(F, n_X) ⊗ transpose(H), H_X ⊗ identity_matrix(F, l))
    Z1 = hcat(transpose(H_X) ⊗ identity_matrix(F, l - 1), identity_matrix(F, n) ⊗ H)
    Z2 = hcat(zero_matrix(F, l * n_Z, (l - 1) * n_X), H_Z ⊗ identity_matrix(F, l))
    Z3 = Z2[[h + l * (i - 1) for (i, h) in enumerate(heights)], :]
    return X, vcat(Z1, Z3)
end

"""
    thickening_and_choose_heights(S::AbstractStabilizerCode, l::Integer, heights::Vector{Int})

Return the result of thickening and choosing heights on `S`.
"""
thickening_and_choose_heights(S::T, l::Integer, heights::Vector{Int}) where {T <:
    AbstractStabilizerCode} = thickening_and_choose_heights(CSSTrait(T), S, l, heights)
thickening_and_choose_heights(::IsCSS, S::AbstractStabilizerCode, l::Integer,
    heights::Vector{Int}) =  CSSCode(thickening_and_choose_heights(S.X_stabs, S.Z_stabs,
    l, heights)...)
thickening_and_choose_heights(::IsNotCSS, S::AbstractStabilizerCode, l::Integer,
    heights::Vector{Int}) = error("Only valid for CSS codes.")

function _cycle_basis_decongestion(_edges::Vector{Tuple{T, T}}) where T
    edges = Vector{T}[[e...] for e in _edges]
    vertices = unique(union(_edges...))
    degrees = zeros(Int, length(vertices))
    for (a, b) in _edges
        degrees[findfirst(a .== vertices)] += 1
        degrees[findfirst(b .== vertices)] += 1
    end
    cycles = Vector{T}[]
    (isempty(_edges) || maximum(degrees) <= 1) && (return cycles;)

    # I added step 0. Steps 1-3 implement the algorithm in the decongestion lemma.
    iter = 0
    max_iter = 1000
    while sum(degrees) > 0 && iter < max_iter
        iter += 1

        # step 0B: remove vertices with no edges (not necessary, will test if it's faster to include later)
        # i = findall(isequal(0), degrees)
        # if !isnothing(i)
        #     deleteat!(vertices, i)
        #     deleteat!(degrees, i)
        # end

        # step 1: remove dangling edges that can't participate in cycles
        while 1 in degrees
            i = findfirst(isequal(1), degrees)
            j = findfirst(vertices[i] == first(e) || vertices[i] == last(e) for e in edges)
            k = findfirst(isequal(setdiff([edges[j][1], edges[j][end]], vertices[i])[1]), vertices)
            degrees[k] -= 1
            degrees[i] -= 1
            deleteat!(edges, j)
        end

        # step 2A: find all self-edges, save the cycles and remove from graph
        # note that this is modified from 2A in the paper, but it doesn't harm their proof.
        j = findall(e -> first(e) == last(e), edges)
        for k in j
            push!(cycles, edges[k][1:end - 1])
            degrees[findfirst(isequal(first(edges[k])), vertices)] -= 2
        end
        isnothing(j) || deleteat!(edges, j)

        # step 2B: if there is a degree 2 vertex that isn't from a self-edge (because all self-edges were removed in 2A), compress the two edges into a single edge
        if 2 ∈ degrees && 1 ∉ degrees
            i = rand(findall(isequal(2), degrees)) # randomize which vertex we choose
            j = findall(vertices[i] == first(e) || vertices[i] == last(e) for e in edges)
            last(edges[j[1]]) == vertices[i] || reverse!(edges[j[1]])
            first(edges[j[2]]) == vertices[i] || reverse!(edges[j[2]])
            append!(edges[j[1]], edges[j[2]][2:end])
            deleteat!(edges, j[2])
            degrees[i] -= 2

        # step 3: if all vertices have degree > 2 then find and save a short cycle, remove an edge of that cycle from the graph
        elseif maximum(degrees) > 2 && 1 ∉ degrees && 2 ∉ degrees
            # At this point, the graph does not have self-edges because I modified step 2A.

            # detect if there are any double-edges, if there are pick randomly:
            indices = nothing
            for e in shuffle(edges)
                j = findall(sort([first(e), last(e)]) == sort([first(e2), last(e2)]) for e2 in edges)
                if length(j) > 1
                    indices = j[randperm(length(j))[1:2]]
                    break
                end
            end

            if !isnothing(indices) # there is a double-edge, so add that cycle and remove one of the edges
                last(edges[indices[1]]) == first(edges[indices[2]]) || reverse!(edges[indices[2]])
                push!(cycles, [edges[indices[1]]; edges[indices[2]][2:end - 1]])
                degrees[findfirst(isequal(first(edges[indices[1]])), vertices)] -= 1
                degrees[findfirst(isequal(last(edges[indices[1]])), vertices)] -= 1
                deleteat!(edges, indices[2])

            else # graph is not a multigraph, must find a short cycle more traditionally
                temp = [(first(e), last(e)) for e in edges]
                tempcycles = Grphs.cycle_basis(Grphs.SimpleGraph(Grphs.SimpleEdge.(temp)))
                _, i = findmin(length, tempcycles)
                c = tempcycles[i]
                length(c) > 2log2(count(!iszero, degrees)) && @warn "cycle $(length(cycles) + 1) is longer than necessary"
                cycle = [c[1]]
                for i in 2:length(c)
                    j = findfirst(sort([first(e), last(e)]) == sort([last(cycle), c[i]]) for e in edges)
                    append!(cycle, first(edges[j]) == last(cycle) ? edges[j][2:end] : edges[j][end - 1:-1:1])
                end
                j = findfirst(sort([first(e), last(e)]) == sort([last(cycle), c[1]]) for e in edges)
                if length(edges[j]) > 2
                    append!(cycle, first(edges[j]) == last(cycle) ? edges[j][2:end - 1] : edges[j][end - 1:-1:2])
                end
                push!(cycles, cycle)
                degrees[findfirst(isequal(first(edges[j])), vertices)] -= 1
                degrees[findfirst(isequal(last(edges[j])), vertices)] -= 1
                deleteat!(edges, j)
            end
        end
    end

    iter == max_iter && @warn "reached maximum iterations, probably not a full cycle basis"
    return cycles
end

"""
    coning(H_X::T, H_Z::T, whichZ::AbstractVector{Int}; l::Int = 0, target_q_X::Int = 3) where T <: CTMatrixTypes

Return the result of coning on `H_X` and `H_Z` by reducing the `Z` stabilizers in
`whichZ` and using the optional arguments `l` and `target_q_X` for an optional round of
thickening and choosing heights.
"""
function coning(H_X::T, H_Z::T, whichZ::AbstractVector{Int}; l::Int = 0, target_q_X::Int = 3) where T <: CTMatrixTypes
    
    F = base_ring(H_X)
    n_X, n = size(H_X)
    n_Z = size(H_Z, 1)

    issubset(whichZ, 1:n_Z) || throw(DomainError(whichZ, "Choice of Z stabilizers is out of bounds."))
    target_q_X >= 3 || throw(DomainError(target_q_X, "Can't target a q_X below 3."))

    # Create the spaces (B_i)_1. Each space is the support of a Z stabilizer
    Bi1 = [getindex.(findall(!iszero, H_Z[i, :]), 2) for i in whichZ]

    # Create the spaces (B_i)_0
    # TODO: This assumes overlap of Bi1[i] and any X stabilizer is 0 or 2 qubits. This is true in Hastings algorithm but won't be true in general. Need to extend to work in the general case.
    Xsupports = [getindex.(findall(!iszero, H_X[i, :]), 2) for i in 1:n_X]
    Bi0 = [Tuple{Int, Tuple{Int, Int}}[(i, (sort(X ∩ Q)...,)) for (i, X) in enumerate(Xsupports) if !iszero(length(X ∩ Q))] for Q in Bi1]

    # Create the spaces (B_i)_{-1} and append some edges to (B_i)_0
    Bim1 = [Vector{Int}[] for _ in eachindex(Bi0)]
    for i in eachindex(Bim1)
        # extract edges from Bi0[i] and find a cycle basis c for the graph they define
        edges = getindex.(Bi0[i], 2)
        isempty(edges) && (edges = Tuple{Int, Int}[];)
        unique_edges = unique(edges)

        # c = Grphs.cycle_basis(Grphs.SimpleGraph(Grphs.SimpleEdge.(unique_edges)))
        c = _cycle_basis_decongestion(unique_edges)

        # Bim1[i] = c

        # cellulate big cycles
        for j in eachindex(c)
            if length(c[j]) > 4 # split it up
                a = 1
                b = length(c[j])
                while b - a - 1 >= 2 # at least 2 points between a and b
                    push!(Bim1[i], c[j][[a, a + 1, b - 1, b]])
                    a += 1
                    b -= 1
                    if b - a > 1 # keep track of new edges that need added
                        new_edge = (-1, (sort(c[j][[a, b]])...,))
                        new_edge ∈ Bi0[i] || push!(Bi0[i], new_edge)
                    end
                end
                if b - a - 1 > 0 # if there are any points left between a and b
                    push!(Bim1[i], c[j][[a; a + 1:b - 1; b]])
                end
            else # keep the cycle as is
                push!(Bim1[i], c[j])
            end
        end

        # add extra cycles in case this is actually a multigraph
        for (j, e) in enumerate(edges)
            if e ∈ edges[1:j - 1]
                push!(Bim1[i], [e...])
            end
        end
    end

    # Build ∂1
    ∂i1 = [zero_matrix(F, length(Bi0[i]), length(Bi1[i])) for i in eachindex(Bi1)]
    for i in eachindex(∂i1)
        for j in axes(∂i1[i], 1)
            for k in axes(∂i1[i], 2)
                if Bi1[i][k] ∈ Bi0[i][j][2] # && Bi0[i][j][1] != -1
                    ∂i1[i][j, k] = 1
                end
            end
        end
    end
    ∂1 = reduce(direct_sum, ∂i1)

    # Build ∂0
    ∂i0 = [zero_matrix(F, length(Bim1[i]), length(Bi0[i])) for i in eachindex(Bi1)]
    for i in eachindex(∂i0)
        cycle_edges = [[(sort(Bim1[i][j][k == 1 ? [1, length(Bim1[i][j])] : [k - 1, k]])...,) for k in 1:length(Bim1[i][j])] for j in eachindex(Bim1[i])]
        for j in axes(∂i0[i], 1)
            if length(cycle_edges[j]) == 2
                # this cycle was added to take care of multigraphs, so special case
                indices = findall(isequal(cycle_edges[j][1]), getindex.(Bi0[i], 2))
                num = count(isequal(cycle_edges[j]), cycle_edges[1:j])
                ∂i0[i][j, indices[num]] = 1
                ∂i0[i][j, indices[num + 1]] = 1
            else
                # normal cycle from the (non-multi) graph cycle decomposition
                for k in axes(∂i0[i], 2)
                    if Bi0[i][k][2] ∈ cycle_edges[j] && Bi0[i][k][2] ∉ getindex.(Bi0[i][1:k - 1], 2)
                        ∂i0[i][j, k] = 1
                    end
                end
            end
        end
    end
    ∂0 = reduce(direct_sum, ∂i0)

    offset = size(∂1, 1)
    ∂1 = direct_sum(∂1, transpose(H_Z[setdiff(1:size(H_Z, 1), whichZ), :]))
    j = 0
    for Q in Bi1
        for q in Q
            j += 1
            ∂1[q + offset, j] = 1
        end
    end

    offset = size(∂0, 1)
    ∂0 = direct_sum(∂0, H_X)
    j = 0
    for S in Bi0
        for (X, _) in S
            j += 1
            if X > 0
                ∂0[X + offset, j] = 1
            end
        end
    end

    # We have the cone code, but it might need thickened "dually" (swap X and Z, thicken, swap X and Z back) because we could have high column weight in X.
    q_X = maximum(count(!iszero, ∂0[:, j]) for j in 1:size(∂0, 2))
    ∂1, ∂0 = if q_X > target_q_X && l > 1
        n_Z = ncols(∂1)
        n_X, n = size(∂0)
        H = matrix(F, diagm(l - 1, l, 0 => ones(Int, l - 1), 1 => ones(Int, l - 1)))
        Z = hcat(identity_matrix(F, n_Z) ⊗ transpose(H), transpose(∂1) ⊗ identity_matrix(F, l))
        X1 = hcat(∂1 ⊗ identity_matrix(F, l - 1), identity_matrix(F, n) ⊗ H)
        X2 = hcat(zero_matrix(F, l * n_X, (l - 1) * n_Z), ∂0 ⊗ identity_matrix(F, l))

        # TODO: figure out a good way to compute heights
        # choosing heights here is not done in the same way as before. We only choose heights for part of it.
        num_heights = n_X - size(H_X, 1)
        heights = rand(1:l, num_heights)

        X3 = X2[[h + l * (i - 1) for (i, h) in enumerate(heights)] ∪ (num_heights * l + 1:size(X2, 1)), :]
        transpose(Z), vcat(X1, X3)
    else
        ∂1, ∂0
    end
    return ∂0, transpose(∂1) # this is the new H_X, H_Z
end

"""
    coning(S::AbstractStabilizerCode, whichZ::AbstractVector{Int}; l::Int = 0, target_q_X::Int = 3) where T <: CTMatrixTypes

Return the result of coning on `S` by reducing the `Z` stabilizers in `whichZ` and using the
optional arguments `l` and `target_q_X` for an optional round of thickening and choosing heights.
"""
coning(S::T, whichZ::AbstractVector{Int}; l::Int, target_q_X::Int = 3) where {T <:
    AbstractStabilizerCode} = coning(CSSTrait(T), S, whichZ, l = l, target_q_X = target_q_X)
coning(::IsCSS, S::AbstractStabilizerCode, whichZ::AbstractVector{Int}; l::Int,
    target_q_X::Int) = CSSCode(coning(S.X_stabs, S.Z_stabs, whichZ, l = l, target_q_X = target_q_X)...)
coning(::IsNotCSS, S::AbstractStabilizerCode, whichZ::AbstractVector{Int}; l::Int,
    target_q_X::Int) = error("Only valid for CSS codes.")

"""
    weight_reduction(S::AbstractStabilizerCode, copying_type::Symbol=:Hastings, copying_target::Int = 3, l1::Int, heights::Vector{Int}, l2::Int = 1, desired_q_X::Int = 3, seed::Union{Nothing, Int} = nothing)
    quantum_weight_reduction(S::AbstractStabilizerCode, copying_type::Symbol=:Hastings, copying_target::Int = 3, l1::Int, heights::Vector{Int}, l2::Int = 1, desired_q_X::Int = 3, seed::Union{Nothing, Int} = nothing)

Return the weight-reduced CSS code of `S`.
"""
quantum_weight_reduction(S::T, l1::Int, heights::Vector{Int}; copying_type::Symbol = :Hastings,
    copying_target::Int = 3, l2::Int = 1, target_q_X::Int = 3,
    seed::Union{Nothing, Int} = nothing) where {T <: AbstractStabilizerCode} =
    quantum_weight_reduction(CSSTrait(T), S, l1, heights, copying_type = copying_type,
        copying_target = copying_target, l2 = l2, target_q_X = target_q_X, seed = seed)
function quantum_weight_reduction(::IsCSS, S::AbstractStabilizerCode, l1::Int, heights::Vector{Int};
    copying_type::Symbol, copying_target::Int, l2::Int, target_q_X::Int,
    seed::Union{Nothing, Int} = nothing)

    Random.seed!(seed)

    copying_type ∈ (:Hastings, :reduced, :target) || throw(ArgumentError("Unknown copying method"))
    # check copying target
    # check coning target
    # check second optional t&ch
    # check everything
    H_X, H_Z = copying(S.X_stabs, S.Z_stabs, method = copying_type, target_q_X = copying_target)
    H_X, H_Z = gauging(H_X, H_Z)
    a = nrows(H_Z)
    H_X, H_Z = thickening_and_choose_heights(H_X, H_Z, l1, heights)
    b = nrows(H_Z)
    whichZ = b - a + 1:b
    H_X, H_Z = coning(H_X, H_Z, whichZ, l = l2, target_q_X = target_q_X)
    return CSSCode(H_X, H_Z)
end
quantum_weight_reduction(::IsNotCSS, S::AbstractStabilizerCode, l1::Int, heights::Vector{Int};
    copying_type::Symbol = :Hastings, copying_target::Int = 3, l2::Int = 1, target_q_X::Int = 3,
    seed::Union{Nothing, Int} = nothing) =
    error("Only valid for CSS codes.")

weight_reduction(S::AbstractStabilizerCode, l1::Int, heights::Vector{Int};
    copying_type::Symbol = :Hastings, copying_target::Int = 3, l2::Int = 1, target_q_X::Int = 3,
    seed::Union{Nothing, Int} = nothing) =
    quantum_weight_reduction(S, l1, heights, copying_type = copying_type,
        copying_target = copying_target, l2 = l2, target_q_X = target_q_X, seed = seed)
