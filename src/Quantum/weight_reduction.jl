# Copyright (c) 2023, 2024 Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
# Copying
#############################

function _copying_Hastings(H_X::CTMatrixTypes, H_Z::CTMatrixTypes)
    nr, n = size(H_X)

    # get column weight
    q_X = maximum(count(!iszero, H_X[:, j]) for j = 1:n)

    # X stabilizers
    X = zero_matrix(base_ring(H_X), nr + (q_X - 1) * n, q_X * n)

    # copied X stabilizers
    for j = 1:n
        for (k, i) in enumerate(getindex.(findall(isone, H_X[:, j]), 1))
            X[i, q_X*(j-1)+k] = H_X[i, j]
        end
    end

    # new X stabilizers
    for j = 1:n
        for i = 1:(q_X-1)
            row = nr + i + (q_X - 1) * (j - 1)
            col = q_X * (j - 1) + i
            X[row, col] = 1
            X[row, col+1] = 1
        end
    end

    # Z stabilizers
    nr = nrows(H_Z)
    Z = zero_matrix(base_ring(H_Z), nr, q_X * n)
    for i = 1:nr
        for j = 1:n
            for k = 1:q_X
                Z[i, (j-1)*q_X+k] = H_Z[i, j]
            end
        end
    end
    return X, Z
end

function _copying_reduced(H_X::CTMatrixTypes, H_Z::CTMatrixTypes)
    n_X, n = size(H_X)

    # get column weights
    qubit_expansion = [count(!iszero, H_X[:, j]) for j = 1:n]
    q_X = maximum(qubit_expansion)

    # X stabilizers
    X = zero_matrix(base_ring(H_X), n_X + sum(qubit_expansion .- 1), sum(qubit_expansion))

    # keep track of how many edges (in the tanner graph) remain to be used for each qubit
    qubits_available = zeros(Int, q_X, n)
    for (j, n2) in enumerate(qubit_expansion)
        qubits_available[1:n2, j] .= 1
    end

    # copy in the old X stabilizers
    for i = 1:n_X
        for j = 1:n
            iszero(H_X[i, j]) && continue
            k = findfirst(!iszero, qubits_available[:, j])
            qubits_available[k, j] -= 1
            X[i, sum(qubit_expansion[1:(j-1)])+k] = H_X[i, j]
        end
    end

    # new X stabilizers
    row = n_X + 1
    for (i, n) in enumerate(qubit_expansion)
        for j = 1:(n-1)
            a = sum(qubit_expansion[1:(i-1)]) + j
            X[row, a] = 1
            X[row, a+1] = 1
            row += 1
        end
    end

    # Z stabilizers
    n_Z = nrows(H_Z)
    Z = zero_matrix(base_ring(H_Z), n_Z, sum(qubit_expansion))
    for i = 1:n_Z
        for j = 1:n
            for k = (sum(qubit_expansion[1:(j-1)])+1):sum(qubit_expansion[1:j])
                Z[i, k] = H_Z[i, j]
            end
        end
    end

    return X, Z
end

function _copying_target(H_X::CTMatrixTypes, H_Z::CTMatrixTypes, target_q_X::Int = 3)
    target_q_X >= 3 ||
        throw(DomainError(target_q_X, "Target column weight must be at least 3."))

    n_X, n = size(H_X)
    # get column weights
    q_Xs = [count(!iszero, H_X[:, j]) for j = 1:n]

    # repetition expansion of qubits
    qubit_expansion = [
        q_X <= target_q_X ? 1 : 2 + cld(q_X - 2 * (target_q_X - 1), target_q_X - 2) for
        q_X in q_Xs
    ]

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
            qubits_available[2:(n2-1), j] .= target_q_X - 2
        end
    end

    # copy in the old X stabilizers
    for i = 1:n_X
        for j = 1:n
            iszero(H_X[i, j]) && continue
            k = findfirst(!iszero, qubits_available[:, j])
            qubits_available[k, j] -= 1
            X[i, sum(qubit_expansion[1:(j-1)])+k] = H_X[i, j]
        end
    end

    # new X stabilizers
    row = n_X + 1
    for (i, n) in enumerate(qubit_expansion)
        for j = 1:(n-1)
            a = sum(qubit_expansion[1:(i-1)]) + j
            X[row, a] = 1
            X[row, a+1] = 1
            row += 1
        end
    end

    # Z stabilizers
    n_Z = nrows(H_Z)
    Z = zero_matrix(base_ring(H_Z), n_Z, sum(qubit_expansion))
    for i = 1:n_Z
        for j = 1:n
            for k = (sum(qubit_expansion[1:(j-1)])+1):sum(qubit_expansion[1:j])
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
function copying(
    H_X::CTMatrixTypes,
    H_Z::CTMatrixTypes;
    method::Symbol = :Hastings,
    target_q_X::Int = 3,
)

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
copying(
    S::T;
    method::Symbol = :Hastings,
    target_q_X::Int = 3,
) where {T<:AbstractStabilizerCode} = copying(CSSTrait(T), S, method, target_q_X)
function copying(::IsCSS, S::AbstractStabilizerCode, method::Symbol, target_q_X::Int)
    method ∈ (:Hastings, :reduced, :target) || throw(ArgumentError("Unknown method type"))
    target_q_X >= 3 || throw(DomainError(target_q_X, "Target must be at least 3"))

    H_X, H_Z = copying(S.X_stabs, S.Z_stabs, method = method, target_q_X = target_q_X)
    return CSSCode(H_X, H_Z)
end
copying(::IsNotCSS, S::AbstractStabilizerCode, method::Symbol, target_q_X::Int) =
    error("Only valid for CSS codes.")

function _copying_as_coning_Hastings(
    H_X::CTMatrixTypes,
    H_Z::CTMatrixTypes;
    permute = false,
    rng::AbstractRNG = Random.seed!(),
)
    q_X = maximum(count(!iszero, H_X[:, j]) for j = 1:ncols(H_X))
    q_X == 1 && return H_X, H_Z
    n_X = nrows(H_X)
    F = base_ring(H_X)
    for i = 1:ncols(H_X)
        H = matrix(F, diagm(q_X - 1, q_X, 0 => ones(Int, q_X - 1), 1 => ones(Int, q_X - 1)))
        f_1 = zero_matrix(F, q_X, nrows(H_X))
        for j in (permute ? shuffle(rng, 1:n_X) : 1:n_X)
            if H_X[j, 1] == 1
                f_1[count(!iszero, f_1)+1, j] = 1
            end
        end
        H_X = H_X[:, 2:end]
        H_Z = H_Z[:, 2:end]
        flag, f_2 = can_solve_with_solution(
            hcat(f_1, transpose(H)),
            hcat(H_Z * transpose(H_X), zero_matrix(GF(2), nrows(H_Z), nrows(H))),
            side = :left,
        )
        flag || error("there was no solution for f_2")
        left_kernel(hcat(f_1, transpose(H)))[1] == 0 ||
            @warn "there was more than one possible f_2, column $i"
        H_X = hcat(
            vcat(H_X, zero_matrix(GF(2), nrows(H), ncols(H_X))),
            vcat(transpose(f_1), H),
        )
        H_Z = hcat(H_Z, f_2)
    end
    return H_X, H_Z
end

function _copying_as_coning_reduced(
    H_X::CTMatrixTypes,
    H_Z::CTMatrixTypes;
    permute = false,
    rng::AbstractRNG = Random.seed!(),
)
    F = base_ring(H_X)
    n_X = nrows(H_X)
    for i = 1:ncols(H_X)
        q = count(!iszero, H_X[:, 1])
        if q <= 1
            H_X = hcat(H_X[:, 2:end], H_X[:, 1:1])
            H_Z = hcat(H_Z[:, 2:end], H_Z[:, 1:1])
            continue
        end
        H = matrix(F, diagm(q - 1, q, 0 => ones(Int, q - 1), 1 => ones(Int, q - 1)))
        f_1 = zero_matrix(F, q, nrows(H_X))
        for j in (permute ? shuffle(rng, 1:n_X) : 1:n_X)
            if H_X[j, 1] == 1
                f_1[count(!iszero, f_1)+1, j] = 1
            end
        end
        H_X = H_X[:, 2:end]
        H_Z = H_Z[:, 2:end]
        flag, f_2 = can_solve_with_solution(
            hcat(f_1, transpose(H)),
            hcat(H_Z * transpose(H_X), zero_matrix(GF(2), nrows(H_Z), nrows(H))),
            side = :left,
        )
        flag || error("there was no solution for f_2")
        left_kernel(hcat(f_1, transpose(H)))[1] == 0 ||
            @warn "there was more than one possible f_2"
        H_X = hcat(
            vcat(H_X, zero_matrix(GF(2), nrows(H), ncols(H_X))),
            vcat(transpose(f_1), H),
        )
        H_Z = hcat(H_Z, f_2)
    end
    return H_X, H_Z
end

function _copying_as_coning_target(
    H_X::CTMatrixTypes,
    H_Z::CTMatrixTypes,
    target_q_X::Int = 3;
    permute = false,
    rng::AbstractRNG = Random.seed!(),
)

    target_q_X < 3 && throw(DomainError(target_q_X, "Must be at least 3"))
    F = base_ring(H_X)
    n_X = nrows(H_X)
    for i = 1:ncols(H_X)
        q = count(!iszero, H_X[:, 1])
        if q <= target_q_X
            H_X = hcat(H_X[:, 2:end], H_X[:, 1:1])
            H_Z = hcat(H_Z[:, 2:end], H_Z[:, 1:1])
            continue
        end
        H = matrix(
            F,
            diagm(
                q - target_q_X,
                q - target_q_X + 1,
                0 => ones(Int, q - target_q_X),
                1 => ones(Int, q - target_q_X),
            ),
        )
        f_1 = zero_matrix(F, q - target_q_X + 1, nrows(H_X))
        for j in (permute ? shuffle(rng, 1:n_X) : 1:n_X)
            if H_X[j, 1] == 1
                k = 1
                while count(!iszero, f_1[k, :]) ==
                      target_q_X - 2 + isone(k) + (k == ncols(H))
                    k += 1
                end
                f_1[k, j] = 1
            end
        end
        H_X = H_X[:, 2:end]
        H_Z = H_Z[:, 2:end]
        flag, f_2 = can_solve_with_solution(
            hcat(f_1, transpose(H)),
            hcat(H_Z * transpose(H_X), zero_matrix(GF(2), nrows(H_Z), nrows(H))),
            side = :left,
        )
        flag || error("there was no solution for f_2")
        left_kernel(hcat(f_1, transpose(H)))[1] == 0 ||
            @warn "there was more than one possible solution for f_2"
        H_X = hcat(
            vcat(H_X, zero_matrix(GF(2), nrows(H), ncols(H_X))),
            vcat(transpose(f_1), H),
        )
        H_Z = hcat(H_Z, f_2)
    end
    return H_X, H_Z
end

"""
    copying_as_coning(H_X::CTMatrixTypes, H_Z::CTMatrixTypes; method::Symbol = :Hastings, target_q_X::Int = 3)

Return the result of copying on `H_X` and `H_Z` using either the Hastings, reduced, or targeted
methods by using the mapping cone.
"""
function copying_as_coning(
    H_X::CTMatrixTypes,
    H_Z::CTMatrixTypes;
    method::Symbol = :Hastings,
    target_q_X::Int = 3,
    rng::AbstractRNG = Random.seed!(),
)

    method ∈ (:Hastings, :reduced, :target) || throw(ArgumentError("Unknown method type"))
    target_q_X >= 3 || throw(DomainError(target_q_X, "Target must be at least 3"))
    # should we check these commute or trust the user?

    if method == :Hastings
        return _copying_as_coning_Hastings(H_X, H_Z, rng = rng)
    elseif method == :reduced
        return _copying_as_coning_reduced(H_X, H_Z, rng = rng)
    else
        return _copying_as_coning_target(H_X, H_Z, target_q_X, rng = rng)
    end
end

"""
    copying_as_coning(S::AbstractStabilizerCode, method::Symbol = :Hastings, target_q_X::Int = 3)

Return the result of copying on `S` using either the Hastings, reduced, or targeted methods
by using the mapping cone.
"""
copying_as_coning(
    S::T;
    method::Symbol = :Hastings,
    target_q_X::Int = 3,
    rng::AbstractRNG = Random.seed!(),
) where {T<:AbstractStabilizerCode} =
    copying_as_coning(CSSTrait(T), S, method, target_q_X, rng)
function copying_as_coning(
    ::IsCSS,
    S::AbstractStabilizerCode,
    method::Symbol,
    target_q_X::Int,
    rng::AbstractRNG = Random.seed!(),
)
    method ∈ (:Hastings, :reduced, :target) || throw(ArgumentError("Unknown method type"))
    target_q_X >= 3 || throw(DomainError(target_q_X, "Target must be at least 3"))

    H_X, H_Z = copying_as_coning(
        S.X_stabs,
        S.Z_stabs,
        method = method,
        target_q_X = target_q_X,
        rng = rng,
    )
    return CSSCode(H_X, H_Z)
end
copying_as_coning(
    ::IsNotCSS,
    S::AbstractStabilizerCode,
    method::Symbol,
    target_q_X::Int,
    rng::AbstractRNG = Random.seed!(),
) = error("Only valid for CSS codes.")

#############################
# Gauging
#############################

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
    for i = 1:n_X
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
            for j = 1:n
                X[i_new, j] = H_X[i, j]
            end
            i_new += 1
        else
            new_qubits = new_qubit_index:(new_qubit_index+wt-3)

            # cut the X stabilizer up in to weight 3 stabilizers
            for (l, k) in enumerate(new_qubits)
                nonzeros = getindex.(findall(!iszero, H_X[i, :]), 2)
                if k == new_qubit_index
                    X[i_new, nonzeros[1]] = H_X[i, nonzeros[1]]
                    X[i_new, nonzeros[2]] = H_X[i, nonzeros[2]]
                    X[i_new, k] = 1
                elseif k == wt - 3 + new_qubit_index
                    X[i_new, nonzeros[end-1]] = H_X[i, nonzeros[end-1]]
                    X[i_new, nonzeros[end]] = H_X[i, nonzeros[end]]
                    X[i_new, k-1] = 1
                else
                    X[i_new, nonzeros[l+1]] = H_X[i, nonzeros[l+1]]
                    X[i_new, k-1] = 1
                    X[i_new, k] = 1
                end
                i_new += 1
            end

            # adjust the new qubits of the Z stabilizers so that they commute with the new X stabilizers
            for j = 1:n_Z
                Z_qubits = getindex.(findall(!iszero, Z[j, :]), 2)
                X_qubits = getindex.(findall(!iszero, H_X[i, :]), 2)
                for m = 1:(wt-3)
                    if isodd(length(Z_qubits ∩ X_qubits[1:(m+1)]))
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
gauging(S::T) where {T<:AbstractStabilizerCode} = gauging(CSSTrait(T), S)
gauging(::IsCSS, S::AbstractStabilizerCode) = CSSCode(gauging(S.X_stabs, S.Z_stabs)...)
gauging(::IsNotCSS, S::AbstractStabilizerCode) = error("Only valid for CSS codes.")

# have not yet introduced the generalization for these other parameters here
"""
    gauging_as_coning(H_X::CTMatrixTypes, H_Z::CTMatrixTypes)

Return the result of gauging on `H_X` and `H_Z` by using the mapping cone.
"""
function gauging_as_coning(
    H_X::CTMatrixTypes,
    H_Z::CTMatrixTypes;
    target_w_X::Int = 3,
    permute = false,
    rng::AbstractRNG = Random.seed!(),
)

    target_w_X < 3 && throw(DomainError(target_w_X, "Must be at least 3"))
    F = base_ring(H_X)
    n = ncols(H_X)
    for i = 1:nrows(H_X)
        w = count(!iszero, H_X[1, :])
        if w <= target_w_X
            H_X = vcat(H_X[2:end, :], H_X[1:1, :])
            continue
        end
        H = matrix(
            F,
            diagm(
                w - target_w_X,
                w - target_w_X + 1,
                0 => ones(Int, w - target_w_X),
                1 => ones(Int, w - target_w_X),
            ),
        )
        f_1 = zero_matrix(F, w - target_w_X + 1, ncols(H_X))
        for j in (permute ? shuffle(rng, 1:n) : 1:n)
            if isone(H_X[1, j])
                k = 1
                while count(!iszero, f_1[k, :]) ==
                      target_w_X - 2 + isone(k) + (k == ncols(H))
                    k += 1
                end
                f_1[k, j] = 1
            end
        end
        flag, f_2 = can_solve_with_solution(transpose(H), f_1 * transpose(H_Z))
        flag || error("there was no solution for f_2")
        H_X = vcat(
            hcat(H_X[2:end, :], zero_matrix(F, nrows(H_X) - 1, nrows(H))),
            hcat(f_1, transpose(H)),
        )
        H_Z = hcat(H_Z, transpose(f_2))
    end
    return H_X, H_Z
end

"""
    gauging_as_coning(S::AbstractStabilizerCode)

Return the result of gauging on `S` by using the mapping cone.
"""
gauging_as_coning(
    S::T;
    rng::AbstractRNG = Random.seed!(),
) where {T<:AbstractStabilizerCode} = gauging_as_coning(CSSTrait(T), S, rng = rng)
gauging_as_coning(::IsCSS, S::AbstractStabilizerCode; rng::AbstractRNG = Random.seed!()) =
    CSSCode(gauging_as_coning(S.X_stabs, S.Z_stabs, rng = rng)...)
gauging_as_coning(
    ::IsNotCSS,
    S::AbstractStabilizerCode;
    rng::AbstractRNG = Random.seed!(),
) = error("Only valid for CSS codes.")

#############################
# Thickening And Choosing Heights
#############################

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
function thickening_and_choose_heights(
    H_X::CTMatrixTypes,
    H_Z::CTMatrixTypes,
    l::Integer,
    heights::Vector{Int},
)

    F = base_ring(H_X)
    n_Z = nrows(H_Z)
    n_X, n = size(H_X)
    @assert length(heights) == n_Z
    @assert all(1 <= x <= l for x in heights)
    H = matrix(F, diagm(l - 1, l, 0 => ones(Int, l - 1), 1 => ones(Int, l - 1)))
    X = hcat(H_X ⊗ identity_matrix(F, l), identity_matrix(F, n_X) ⊗ transpose(H))
    Z1 = hcat(identity_matrix(F, n) ⊗ H, transpose(H_X) ⊗ identity_matrix(F, l - 1))
    Z2 = hcat(H_Z ⊗ identity_matrix(F, l), zero_matrix(F, l * n_Z, (l - 1) * n_X))
    Z3 = Z2[[h + l * (i - 1) for (i, h) in enumerate(heights)], :]
    return X, vcat(Z3, Z1)
end

"""
    thickening_and_choose_heights(S::AbstractStabilizerCode, l::Integer, heights::Vector{Int})

Return the result of thickening and choosing heights on `S`.
"""
thickening_and_choose_heights(
    S::T,
    l::Integer,
    heights::Vector{Int},
) where {T<:AbstractStabilizerCode} =
    thickening_and_choose_heights(CSSTrait(T), S, l, heights)
thickening_and_choose_heights(
    ::IsCSS,
    S::AbstractStabilizerCode,
    l::Integer,
    heights::Vector{Int},
) = CSSCode(thickening_and_choose_heights(S.X_stabs, S.Z_stabs, l, heights)...)
thickening_and_choose_heights(
    ::IsNotCSS,
    S::AbstractStabilizerCode,
    l::Integer,
    heights::Vector{Int},
) = error("Only valid for CSS codes.")

#############################
# Coning
#############################

function _cycle_basis_decongestion(
    _edges::Vector{Tuple{T,T}};
    rng::AbstractRNG = Random.seed!(),
) where {T}
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
            k = findfirst(
                isequal(setdiff([edges[j][1], edges[j][end]], vertices[i])[1]),
                vertices,
            )
            degrees[k] -= 1
            degrees[i] -= 1
            deleteat!(edges, j)
        end

        # step 2A: find all self-edges, save the cycles and remove from graph
        # note that this is modified from 2A in the paper, but it doesn't harm their proof.
        j = findall(e -> first(e) == last(e), edges)
        for k in j
            push!(cycles, edges[k][1:(end-1)])
            degrees[findfirst(isequal(first(edges[k])), vertices)] -= 2
        end
        isnothing(j) || deleteat!(edges, j)

        # step 2B: if there is a degree 2 vertex that isn't from a self-edge (because all self-edges were removed in 2A), compress the two edges into a single edge
        if 2 ∈ degrees && 1 ∉ degrees
            i = rand(rng, findall(isequal(2), degrees)) # randomize which vertex we choose
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
            for e in shuffle(rng, edges)
                j = findall(
                    sort([first(e), last(e)]) == sort([first(e2), last(e2)]) for
                    e2 in edges
                )
                if length(j) > 1
                    indices = j[randperm(length(j))[1:2]]
                    break
                end
            end

            if !isnothing(indices) # there is a double-edge, so add that cycle and remove one of the edges
                last(edges[indices[1]]) == first(edges[indices[2]]) ||
                    reverse!(edges[indices[2]])
                push!(cycles, [edges[indices[1]]; edges[indices[2]][2:(end-1)]])
                degrees[findfirst(isequal(first(edges[indices[1]])), vertices)] -= 1
                degrees[findfirst(isequal(last(edges[indices[1]])), vertices)] -= 1
                deleteat!(edges, indices[2])

            else # graph is not a multigraph, must find a short cycle more traditionally
                temp = [(first(e), last(e)) for e in edges]
                tempcycles = Grphs.cycle_basis(Grphs.SimpleGraph(Grphs.SimpleEdge.(temp)))
                _, i = findmin(length, tempcycles)
                c = tempcycles[i]
                length(c) > 2log2(count(!iszero, degrees)) &&
                    @warn "cycle $(length(cycles) + 1) is longer than necessary"
                cycle = [c[1]]
                for i = 2:length(c)
                    j = findfirst(
                        sort([first(e), last(e)]) == sort([last(cycle), c[i]]) for
                        e in edges
                    )
                    append!(
                        cycle,
                        first(edges[j]) == last(cycle) ? edges[j][2:end] :
                        edges[j][(end-1):-1:1],
                    )
                end
                j = findfirst(
                    sort([first(e), last(e)]) == sort([last(cycle), c[1]]) for e in edges
                )
                if length(edges[j]) > 2
                    append!(
                        cycle,
                        first(edges[j]) == last(cycle) ? edges[j][2:(end-1)] :
                        edges[j][(end-1):-1:2],
                    )
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
    coning(H_X::T, H_Z::T, row_indices::AbstractVector{Int}; rng::AbstractRNG = Random.seed!()) where T <: CTMatrixTypes

Return the result of coning on `H_X` and `H_Z` by reducing the `Z` stabilizers in `row_indices`. The
optional argument `rng` can be used to make the output of this function reproducible.
"""
function coning(
    H_X::T,
    H_Z::T,
    row_indices::AbstractVector{Int};
    rng::AbstractRNG = Random.seed!(),
) where {T<:CTMatrixTypes}

    F = base_ring(H_X)
    n_X, n = size(H_X)
    n_Z = size(H_Z, 1)

    issubset(row_indices, 1:n_Z) ||
        throw(DomainError(row_indices, "Choice of Z stabilizers is out of bounds."))

    hx = _Flint_matrix_to_Julia_int_matrix(H_X)
    hz = _Flint_matrix_to_Julia_int_matrix(H_Z)
    f1_all = Matrix{Int}[]
    f0_all = Matrix{Int}[]
    p1_all = Matrix{Int}[]
    p0_all = Matrix{Int}[]

    for r in row_indices
        Q = findall(isone, hz[r, :])

        f1 = zeros(Int, n, length(Q))
        for (i, q) in enumerate(Q)
            f1[q, i] = 1
        end

        p1 = hx[:, Q]
        zero_rows = findall(iszero, p1[i, :] for i = 1:n_X)
        nonzero_rows = setdiff(1:n_X, zero_rows)
        p1 = p1[nonzero_rows, :]
        all(2 == count(isone, p1[i, :]) for i = 1:size(p1, 1)) || throw(
            DomainError(
                row_indices,
                "Z stabilizers chosen for reducing must have overlap of at most 2 with all X stabilizers.",
            ),
        )

        f0 = zeros(Int, n_X, length(nonzero_rows))
        for (j, i) in enumerate(nonzero_rows)
            f0[i, j] = 1
        end

        p0 = zeros(Int, 0, length(nonzero_rows))
        edges = [tuple(findall(isone, p1[i, :])...) for i in axes(p1, 1)]
        unique_edges = unique(edges)
        cycles = _cycle_basis_decongestion(unique_edges, rng = rng)
        # cycles = Grphs.cycle_basis(Grphs.SimpleGraph(Grphs.SimpleEdge.(unique_edges)))

        # take care of length 2 cycles
        for (i, e) in enumerate(unique_edges)
            if count(isequal(e), edges) > 1
                repeats = findall(1 == p1[j, e[1]] == p1[j, e[2]] for j = 1:size(p1, 1))
                nr, nc = size(p0)
                p0 = vcat(p0, zeros(Int, length(repeats) - 1, nc))
                for j = 1:(length(repeats)-1)
                    p0[nr+j, repeats[j:(j+1)]] .= 1
                end
            end
        end

        # take care of length 3 and longer cycles
        for c in cycles
            num_cells = max(1, div(length(c) - 1, 2))
            rows, cols = size(p0)
            p0 = vcat(p0, zeros(Int, num_cells, cols))
            if num_cells == 1 # i.e., if length(c) <= 4, don't cellulate
                if length(c) == 3
                    e1 = rand(
                        rng,
                        findall(1 == p1[j, c[1]] == p1[j, c[2]] for j = 1:size(p1, 1)),
                    )
                    e2 = rand(
                        rng,
                        findall(1 == p1[j, c[2]] == p1[j, c[3]] for j = 1:size(p1, 1)),
                    )
                    e3 = rand(
                        rng,
                        findall(1 == p1[j, c[3]] == p1[j, c[1]] for j = 1:size(p1, 1)),
                    )
                    p0[end, [e1, e2, e3]] .= 1
                else
                    e1 = rand(
                        rng,
                        findall(1 == p1[j, c[1]] == p1[j, c[2]] for j = 1:size(p1, 1)),
                    )
                    e2 = rand(
                        rng,
                        findall(1 == p1[j, c[2]] == p1[j, c[3]] for j = 1:size(p1, 1)),
                    )
                    e3 = rand(
                        rng,
                        findall(1 == p1[j, c[3]] == p1[j, c[4]] for j = 1:size(p1, 1)),
                    )
                    e4 = rand(
                        rng,
                        findall(1 == p1[j, c[4]] == p1[j, c[1]] for j = 1:size(p1, 1)),
                    )
                    p0[end, [e1, e2, e3, e4]] .= 1
                end
            else # cellulate
                p0 = hcat(p0, zeros(Int, rows + num_cells, num_cells - 1))
                p1rows = size(p1, 1)
                p1 = vcat(p1, zeros(Int, num_cells - 1, size(p1, 2)))
                f0 = hcat(f0, zeros(Int, size(f0, 1), num_cells - 1))
                for i = 1:(num_cells-1)
                    p0[(rows+i):(rows+i+1), cols+i] .= 1
                    p1[p1rows+i, c[[1 + i, end - i]]] .= 1
                end
                for i = 1:num_cells
                    e1 = rand(
                        rng,
                        findall(1 == p1[j, c[i]] == p1[j, c[i+1]] for j = 1:size(p1, 1)),
                    )
                    e2 = rand(
                        rng,
                        findall(
                            1 == p1[j, c[end-i+1]] == p1[j, c[end-i]] for j = 1:size(p1, 1)
                        ),
                    )
                    p0[rows+i, [e1, e2]] .= 1
                    if i == 1
                        e3 = rand(
                            rng,
                            findall(
                                1 == p1[j, c[1]] == p1[j, c[end]] for j = 1:size(p1, 1)
                            ),
                        )
                        p0[rows+i, e3] = 1
                    elseif i == num_cells && iseven(length(c))
                        e3 = rand(
                            rng,
                            findall(
                                1 == p1[j, c[end-i]] == p1[j, c[i+1]] for j = 1:size(p1, 1)
                            ),
                        )
                        p0[rows+i, e3] = 1
                    end
                end
            end
        end

        push!(f1_all, f1)
        push!(f0_all, f0)
        push!(p1_all, p1)
        push!(p0_all, p0)
    end

    p1 = mapreduce(x -> matrix(F, x), direct_sum, p1_all)
    p0 = mapreduce(x -> matrix(F, x), direct_sum, p0_all)
    f1 = mapreduce(x -> matrix(F, x), hcat, f1_all)
    f0 = mapreduce(x -> matrix(F, x), hcat, f0_all)

    H_X_new = hcat(vcat(p0, f0), vcat(zero_matrix(F, size(p0, 1), n), H_X))
    H_Z_remaining_tr = transpose(H_Z[setdiff(1:n_Z, row_indices), :])
    H_Z_new_tr = hcat(
        vcat(p1, f1),
        vcat(zero_matrix(F, size(p1, 1), size(H_Z_remaining_tr, 2)), H_Z_remaining_tr),
    )

    return H_X_new, transpose(H_Z_new_tr)
end


"""
    coning(S::AbstractStabilizerCode, row_indices::AbstractVector{Int}; rng::AbstractRNG = Random.seed!()) where T <: CTMatrixTypes

Return the result of coning on `S` by reducing the `Z` stabilizers in `row_indices`. The optional
argument `rng` can be used to make the output of this function reproducible.
"""
coning(
    S::T,
    row_indices::AbstractVector{Int};
    rng::AbstractRNG = Random.seed!(),
) where {T<:AbstractStabilizerCode} = coning(CSSTrait(T), S, row_indices, rng = rng)
coning(
    ::IsCSS,
    S::AbstractStabilizerCode,
    row_indices::AbstractVector{Int};
    rng::AbstractRNG = Random.seed!(),
) = CSSCode(coning(S.X_stabs, S.Z_stabs, row_indices, rng = rng)...)
coning(
    ::IsNotCSS,
    S::AbstractStabilizerCode,
    row_indices::AbstractVector{Int};
    rng::AbstractRNG = Random.seed!(),
) = error("Only valid for CSS codes.")


#############################
# All
#############################

"""
    weight_reduction(S::AbstractStabilizerCode, copying_type::Symbol=:Hastings, copying_target::Int = 3, l1::Int, heights::Vector{Int}, l2::Int = 1, rng::AbstractRNG = Random.seed!())
    quantum_weight_reduction(S::AbstractStabilizerCode, copying_type::Symbol=:Hastings, copying_target::Int = 3, l1::Int, heights::Vector{Int}, l2::Int = 1, rng::AbstractRNG = Random.seed!())

Return the weight-reduced CSS code of `S`.
"""
quantum_weight_reduction(
    S::T,
    l1::Int,
    heights::Vector{Int};
    copying_type::Symbol = :Hastings,
    copying_target::Int = 3,
    l2::Int = 1,
    rng::AbstractRNG = Random.seed!(),
) where {T<:AbstractStabilizerCode} = quantum_weight_reduction(
    CSSTrait(T),
    S,
    l1,
    heights,
    copying_type = copying_type,
    copying_target = copying_target,
    l2 = l2,
    rng = rng,
)
function quantum_weight_reduction(
    ::IsCSS,
    S::AbstractStabilizerCode,
    l1::Int,
    heights::Vector{Int};
    copying_type::Symbol,
    copying_target::Int,
    l2::Int,
    rng::AbstractRNG = Random.seed!(),
)

    copying_type ∈ (:Hastings, :reduced, :target) ||
        throw(ArgumentError("Unknown copying method"))
    H_X, H_Z =
        copying(S.X_stabs, S.Z_stabs, method = copying_type, target_q_X = copying_target)
    H_X, H_Z = gauging(H_X, H_Z)
    a = nrows(H_Z)
    H_X, H_Z = thickening_and_choose_heights(H_X, H_Z, l1, heights)
    row_indices = 1:a
    n_X_precone = size(H_X, 1)
    H_X, H_Z = coning(H_X, H_Z, row_indices, rng = rng)
    H_X, H_Z = if l2 > 1
        n_Z = size(H_Z, 1)
        n_X, n = size(H_X)
        F = base_ring(H_X)
        H = matrix(F, diagm(l2 - 1, l2, 0 => ones(Int, l2 - 1), 1 => ones(Int, l2 - 1)))
        Z = hcat(identity_matrix(F, n_Z) ⊗ transpose(H), H_Z ⊗ identity_matrix(F, l2))
        X1 =
            hcat(transpose(H_Z) ⊗ identity_matrix(F, l2 - 1), identity_matrix(F, n) ⊗ H)
        X2 =
            hcat(zero_matrix(F, l2 * n_X, (l2 - 1) * n_Z), H_X ⊗ identity_matrix(F, l2))

        # TODO: figure out a good way to compute heights
        # choosing heights here is not done in the same way as before. We only choose heights for part of it.
        num_heights = n_X - n_X_precone
        heights = rand(1:l2, num_heights)

        X3 = X2[
            [h+l2*(i-1) for (i, h) in enumerate(heights)]∪((num_heights*l2+1):size(X2, 1)),
            :,
        ]
        vcat(X1, X3), Z
    else
        H_X, H_Z
    end
    return CSSCode(H_X, H_Z)
end
quantum_weight_reduction(
    ::IsNotCSS,
    S::AbstractStabilizerCode,
    l1::Int,
    heights::Vector{Int};
    copying_type::Symbol = :Hastings,
    copying_target::Int = 3,
    l2::Int = 1,
    rng::AbstractRNG = Random.seed!(),
) = error("Only valid for CSS codes.")

weight_reduction(
    S::AbstractStabilizerCode,
    l1::Int,
    heights::Vector{Int};
    copying_type::Symbol = :Hastings,
    copying_target::Int = 3,
    l2::Int = 1,
    rng::AbstractRNG = Random.seed!(),
) = quantum_weight_reduction(
    S,
    l1,
    heights,
    copying_type = copying_type,
    copying_target = copying_target,
    l2 = l2,
    rng = rng,
)
