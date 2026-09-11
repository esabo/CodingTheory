# Copyright (c) 2021 - 2024, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
      # Subsystem codes
#############################

function GaugedShorCode()
    # Poulin, "Stabilizer Formalism for Operator Quantum Error Correction", (2008)
    # [[9, 1, 4, 3]] gauged Shor code
    S = ["XXXXXXIII", "XXXIIIXXX", "ZZIZZIZZI","IZZIZZIZZ"]
    # these are the {X, Z} pairings
    Gops = ["IZZIIIIII", "IIXIIIIIX", "IIIIZZIII", "IIIIIXIIX", "ZZIIIIIII", "XIIIIIXII", "IIIZZIIII", "IIIXIIXII"]
    # G = S ∪ Gops
    L = ["ZZZZZZZZZ", "XXXXXXXXX"]
    Q = SubsystemCode(S, L, Gops)
    set_dressed_minimum_distance!(Q, 3)
    return Q
end
Q9143() = GaugedShorCode()

"""
    BaconShorCode(m::Int, n::Int)

Return the Bacon-Shor subsystem code on a `m x n` lattice.
"""
function BaconShorCode(m::Int, n::Int)
    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)
    num_qubits = m * n

    # X stabilizers: X[r, :] = X[r + 1, :] = 1
    # X gauge: X[r, c] = X[r + 1, c] = 1
    X_stabs = zero_matrix(F, m - 1, num_qubits)
    X_gauges = zero_matrix(F, (m - 1) * n, num_qubits)
    curr_row = 1
    curr_row_G = 1
    for r in 1:m - 1
        for c in 1:n
            X_stabs[curr_row, (r - 1) * n + c] = F_one
            X_stabs[curr_row, r * n + c] = F_one
            X_gauges[curr_row_G, (r - 1) * n + c] = F_one
            X_gauges[curr_row_G, r * n + c] = F_one
            curr_row_G += 1
        end
        curr_row += 1
    end

    # Z stabilizers: Z[:, c] = Z[:, c + 1] = 1
    # Z gauge: Z[r, c] = Z[r, c + 1] = 1
    Z_stabs = zero_matrix(F, n - 1, num_qubits)
    Z_gauges = zero_matrix(F, (n - 1) * m, num_qubits)
    curr_row = 1
    curr_row_G = 1
    for c in 1:n - 1
        for r in 1:m
            Z_stabs[curr_row, (r - 1) * n + c] = F_one
            Z_stabs[curr_row, (r - 1) * n + c + 1] = F_one
            Z_gauges[curr_row_G, (r - 1) * n + c] = F_one
            Z_gauges[curr_row_G, (r - 1) * n + c + 1] = F_one
            curr_row_G += 1
        end
        curr_row += 1
    end

    # X logical: X[1, :] = 1
    X_logical = zero_matrix(F, 1, num_qubits)
    for c in 1:n
        X_logical[1, c] = F_one
    end

    # Z logical: Z[:, 1] = 1
    Z_logical = zero_matrix(F, 1, num_qubits)
    for r in 1:m
        Z_logical[1, (r - 1) * n + 1] = F_one
    end
    
    stabs = X_stabs ⊕ Z_stabs
    logs = X_logical ⊕ Z_logical
    gauges = X_gauges ⊕ Z_gauges

    S = SubsystemCode(gauges)
    set_stabilizers!(S, stabs)
    
    # Update CSS properties in the struct
    S.X_stabs = X_stabs
    S.Z_stabs = Z_stabs
    
    set_logicals!(S, logs)
    
    # Cache routing for distances
    S.cache[:dx_dressed] = n
    S.cache[:dz_dressed] = m
    return S
end

"""
    BaconShorCode(d::Int)

Return the Bacon-Shor subsystem code on a `d x d` lattice.
"""
BaconShorCode(d::Int) = BaconShorCode(d, d)

"""
    BravyiBaconShorCode(A::fqPolyRepMatrix)
    GeneralizedBaconShorCode(A::fqPolyRepMatrix)

Return the generalied Bacon-Shor code defined by Bravyi in "Subsystem Codes With Spatially Local
Generators", (2011).
"""
function BravyiBaconShorCode(A::CTMatrixTypes)
    iszero(A) && throw(ArgumentError("The input matrix cannot be zero."))
    F = base_ring(A)
    Int(order(F)) == 2 || throw(ArgumentError("Construction is only valid for binary martices."))

    n = 0
    nr, nc = size(A)
    row_wts = zeros(Int, 1, nr)
    col_wts = zeros(Int, 1, nc)
    linear_index = Dict{Tuple{Int, Int}, Int}()
    for r in 1:nr
        for c in 1:nc
            if !iszero(A[r, c])
                n += 1
                row_wts[r] += 1
                col_wts[c] += 1
                linear_index[(r, c)] = n
            end
        end
    end

    for i in 1:nr
        row_wts[i] == 1 && throw(ArgumentError("The input matrix cannot have a row of weight one."))
    end
    for i in 1:nc
        col_wts[i] == 1 && throw(ArgumentError("The input matrix cannot have a column of weight one."))
    end

    # the original paper appears to switch X and Z compared to Bacon-Shor
    # but this is corrected here to match
    # X - consequetive column pairs
    X_gauges = zero_matrix(F, sum([col_wts[i] - 1 for i in 1:length(col_wts)]), n)
    curr_row = 1
    F_one = F(1)
    for c in 1:nc
        r1 = 1
        while r1 < nr
            if !iszero(A[r1, c])
                at_end = true
                for r2 in r1 + 1:nr
                    if !iszero(A[r2, c])
                        X_gauges[curr_row, linear_index[r1, c]] = F_one
                        X_gauges[curr_row, linear_index[r2, c]] = F_one
                        curr_row += 1
                        r1 = r2
                        at_end = false
                        break
                    end
                end
                at_end && (r1 = nr;)
            else
                r1 += 1
            end
        end
    end

    # Z - consequetive row pairs
    Z_gauges = zero_matrix(F, sum([row_wts[i] - 1 for i in 1:length(row_wts)]), n)
    curr_row = 1
    for r in 1:nr
        c1 = 1
        while c1 < nc
            if !iszero(A[r, c1])
                at_end = true
                for c2 in c1 + 1:nc
                    if !iszero(A[r, c2])
                        Z_gauges[curr_row, linear_index[r, c1]] = F_one
                        Z_gauges[curr_row, linear_index[r, c2]] = F_one
                        curr_row += 1
                        c1 = c2
                        at_end = false
                        break
                    end
                end
                at_end && (c1 = nc;)
            else
                c1 += 1
            end
        end
    end

    S = SubsystemCode(X_gauges ⊕ Z_gauges)
    set_dressed_X_minimum_distance!(S, minimum(row_wts))
    set_dressed_Z_minimum_distance!(S, minimum(col_wts))
    return S
end
GeneralizedBaconShorCode(A::CTMatrixTypes) = BravyiBaconShorCode(A)

function LocalBravyiBaconShorCode(A::CTMatrixTypes)
    iszero(A) && throw(ArgumentError("The input matrix cannot be zero."))
    F = base_ring(A)
    Int(order(F)) == 2 || throw(ArgumentError("Construction is only valid for binary martices."))

    nr, nc = size(A)
    row_firsts = zeros(Int, 1, nr)
    row_lasts = zeros(Int, 1, nr)
    for r in 1:nr
        for c in 1:nc
            if !iszero(A[r, c])
                row_firsts[r] = c
                break
            end
        end
    end
    for r in 1:nr
        for c in nc:-1:1
            if !iszero(A[r, c])
                row_lasts[r] = c
                break
            end
        end
    end

    for i in 1:nr
        row_firsts[i] == row_lasts[i] && throw(ArgumentError("The input matrix cannot have a row of weight one."))
    end

    col_firsts = zeros(Int, 1, nc)
    col_lasts = zeros(Int, 1, nc)
    for c in 1:nc
        for r in 1:nr
            if !iszero(A[r, c])
                col_firsts[c] = r
                break
            end
        end
    end
    for c in 1:nc
        for r in nr:-1:1
            if !iszero(A[r, c])
                col_lasts[c] = r
                break
            end
        end
    end

    for i in 1:nc
        col_firsts[i] == col_lasts[i] && throw(ArgumentError("The input matrix cannot have a column of weight one."))
    end

    n = 0
    extra_Xs = 0
    extra_Zs = 0
    linear_index = Dict{Tuple{Int, Int}, Int}()
    for r in 1:nr
        for c in 1:nc
            if row_firsts[r] <= c <= row_lasts[r] || col_firsts[c] <= r <= col_lasts[c]
                n += 1
                linear_index[(r, c)] = n
            end
            if iszero(A[r, c])
                row_firsts[r] <= c <= row_lasts[r] && (extra_Xs += 1;)
                col_firsts[c] <= r <= col_lasts[c] && (extra_Zs += 1;)
            end
        end
    end

    # X - consequetive column pairs + single qubit operators for each 0 in a row
    X_gauges = zero_matrix(F, sum([row_lasts[i] - row_firsts[i] for i in 1:nr]) + extra_Xs, n)
    # Z - consequetive row pairs + single qubit operators for each 0 in a column
    Z_gauges = zero_matrix(F, sum([col_lasts[i] - col_firsts[i] for i in 1:nc]) + extra_Zs, n)
    curr_row_X = 1
    curr_row_Z = 1

    # the original paper appears to switch X and Z compared to Bacon-Shor
    # but this is corrected here to match
    F_one = F(1)
    for c in 1:nc
        for r in col_firsts[c]:col_lasts[c]
            if r != col_lasts[c]
                X_gauges[curr_row_X, linear_index[r, c]] = F_one
                X_gauges[curr_row_X, linear_index[r + 1, c]] = F_one
                curr_row_X += 1
            end
            if iszero(A[r, c])
                Z_gauges[curr_row_Z, linear_index[r, c]] = F_one
                curr_row_Z += 1
            end
        end
    end

    for r in 1:nr
        for c in row_firsts[r]:row_lasts[r]
            if c != row_lasts[r]
                Z_gauges[curr_row_Z, linear_index[r, c]] = F_one
                Z_gauges[curr_row_Z, linear_index[r, c + 1]] = F_one
                curr_row_Z += 1
            end
            if iszero(A[r, c])
                X_gauges[curr_row_X, linear_index[r, c]] = F_one
                curr_row_X += 1
            end
        end
    end

    S = SubsystemCode(X_gauges ⊕ Z_gauges, logs_alg = :sys_eqs)
    return S
end
AugmentedBravyiBaconShorCode(A::CTMatrixTypes) = LocalBravyiBaconShorCode(A)

# subsystem codes were described by Bacon and Casaccino in [19]. The construction of [19] starts from a pair of classical linear codes C1 = [n1, k1, d1] and C2 = [n2, k2, d2]. A quantum subsystem code is then defined by placing a physical qubit at every cell of a ma- trix A of size n1 × n2. The X-part of the gauge group is defined by replicating the parity checks of C1 in every col- umn of A (in the X-basis). Similarly, the Z-part of the gauge group is defined by replicating the parity checks of C2 in every row of A (in the Z-basis). The resulting subsystem code has parameters [n1n2, k1k2, min (d1, d2)].

# Napp & Preskill, "Optimal Bacon-Shor Codes", (2012)
"""
    NappPreskill3DCode(m::Int, n::Int, k::Int)

Return the Napp and Preskill 3D, modifed Bacon-Shor code.
"""
function NappPreskill3DCode(m::Int, n::Int, k::Int)
    (2 <= m && 2 <= n && 2 <= k) || throw(DomainError("Lattice dimensions must be at least two"))

    linear_index = Dict{Tuple{Int, Int, Int}, Int}()
    for (i, tup) in enumerate(Base.Iterators.product(1:m, 1:n, 1:k))
        linear_index[tup] = i
    end

    len = m * n * k
    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)
    gauges = zero_matrix(F, (m - 1) * n * k + m * (n - 1) * k + m * n * (k - 1), 2 * len)
    curr_row = 1
    # X's
    for l in 1:k
        for i in 1:m
            for j in 1:n
                if j != n
                    gauges[curr_row, linear_index[i, j, l]] = gauges[curr_row, linear_index[i, j + 1, l]] = F_one
                    curr_row += 1
                end
                if i != m
                    gauges[curr_row, linear_index[i, j, l]] = gauges[curr_row, linear_index[i + 1, j, l]] = F_one
                    curr_row += 1
                end
            end
        end
    end

    # Z's
    for i in 1:m
        for j in 1:n
            for l in 1:k - 1
                gauges[curr_row, linear_index[i, j, l] + len] = gauges[curr_row, linear_index[i, j, l + 1] + len] = F_one
                curr_row += 1
            end
        end
    end

    S = SubsystemCode(gauges, logs_alg = :sys_eqs)
    set_dressed_X_minimum_distance!(S, k)
    set_dressed_Z_minimum_distance!(S, m * n)
    return S
end

# Napp & Preskill, "Optimal Bacon-Shor Codes", (2012)
"""
    NappPreskill4DCode(x::Int, y::Int, z::Int, w::Int)

Return the Napp and Preskill 4D, modifed Bacon-Shor code.
"""
function NappPreskill4DCode(x::Int, y::Int, z::Int, w::Int)
    (2 <= x && 2 <= y && 2 <= z && 2 <= w) || throw(DomainError("Lattice dimensions must be at least two"))

    linear_index = Dict{Tuple{Int, Int, Int, Int}, Int}()
    for (i, tup) in enumerate(Base.Iterators.product(1:x, 1:y, 1:z, 1:w))
        linear_index[tup] = i
    end

    len = x * y * z * w
    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)
    gauges = zero_matrix(F, (x - 1) * y * z * w + x * (y - 1) * z * w + x * y * (z - 1) * w + x * y * z * (w - 1), 2 * len)
    curr_row = 1
    ## XX acting on each neighboring qubits in each xy-plane with z, w fixed
    ## ZZ in zw-plane with x, y fixed

    # X's
    for k in 1:z
        for l in 1:w
            for i in 1:x
                for j in 1:y
                    if i != x
                        gauges[curr_row, linear_index[i, j, k, l]] = gauges[curr_row, linear_index[i + 1, j, k, l]] = F_one
                        curr_row += 1
                    end
                    if j != y
                        gauges[curr_row, linear_index[i, j, k, l]] = gauges[curr_row, linear_index[i, j + 1, k, l]] = F_one
                        curr_row += 1
                    end
                end
            end
        end
    end

    # Z's
    for i in 1:x
        for j in 1:y
            for k in 1:z
                for l in 1:w
                    if k != z
                        gauges[curr_row, linear_index[i, j, k, l] + len] = gauges[curr_row, linear_index[i, j, k + 1, l] + len] = F_one
                        curr_row += 1
                    end
                    if l != w
                        gauges[curr_row, linear_index[i, j, k, l] + len] = gauges[curr_row, linear_index[i, j, k, l + 1] + len] = F_one
                        curr_row += 1
                    end
                end
            end
        end
    end

    return SubsystemCode(gauges, logs_alg = :sys_eqs)
end

# Bravyi et al, "Subsystem surface codes with three-qubit check operators", (2013)
"""
    SubsystemToricCode(m::Int, n::Int)

Return the subsystem toric code on a rectangular lattice with `m` rows and `n`
columns of squares.
"""
function SubsystemToricCode(m::Int, n::Int)
    (2 <= m && 2 <= n) || throw(DomainError("Lattice dimensions must be at least two"))

    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)
    len = 3 * m * n
    gauges = zero_matrix(F, 4 * m * n, 2 * len)
    stabs = zero_matrix(F, 2 * m * n, 2 * len)
    curr_row_stab = 1
    curr_row_gauge = 1
    for r in 1:m
        top_left = 3 * n * (r - 1) + 1
        row_right = top_left + 2 * n - 1
        for c in 1:n
            # top left - Z
            gauges[curr_row_gauge, top_left + len] = F_one
            gauges[curr_row_gauge, top_left + 1 + len] = F_one
            gauges[curr_row_gauge, row_right + c + len] = F_one
            curr_row_gauge += 1

            # top right - X
            gauges[curr_row_gauge, top_left + 1] = F_one
            if c != n
                gauges[curr_row_gauge, top_left + 2] = F_one
                gauges[curr_row_gauge, row_right + c + 1] = F_one
            else
                gauges[curr_row_gauge, row_right - 2 * n + 1] = F_one
                gauges[curr_row_gauge, row_right + 1] = F_one
            end
            curr_row_gauge += 1

            # bottom left - X
            gauges[curr_row_gauge, row_right + c] = F_one
            if r != m
                gauges[curr_row_gauge, row_right + n + 2 * (c - 1) + 1] = F_one
                gauges[curr_row_gauge, row_right + n + 2 * (c - 1) + 2] = F_one
            else
                gauges[curr_row_gauge, 2 * (c - 1) + 1] = F_one
                gauges[curr_row_gauge, 2 * (c - 1) + 2] = F_one
            end
            curr_row_gauge += 1

            # bottom right - Z
            if r == m && c == n
                gauges[curr_row_gauge, 1 + len] = F_one
                gauges[curr_row_gauge, 2 * n + len] = F_one
                gauges[curr_row_gauge, row_right + 1 + len] = F_one
            elseif r == m
                gauges[curr_row_gauge, row_right + c + 1 + len] = F_one
                gauges[curr_row_gauge, 2 * (c - 1) + 2 + len] = F_one
                gauges[curr_row_gauge, 2 * (c - 1) + 3 + len] = F_one
            elseif c != n
                gauges[curr_row_gauge, row_right + c + 1 + len] = F_one
                gauges[curr_row_gauge, row_right + n + 2 * (c - 1) + 2 + len] = F_one
                gauges[curr_row_gauge, row_right + n + 2 * (c - 1) + 3 + len] = F_one
            else
                gauges[curr_row_gauge, row_right + 1 + len] = F_one
                gauges[curr_row_gauge, row_right + n + 1 + len] = F_one
                gauges[curr_row_gauge, row_right + n + 2 * (c - 1) + 2 + len] = F_one
            end
            curr_row_gauge += 1

            # X
            stabs[curr_row_stab, top_left + 1] = F_one
            if c != n
                stabs[curr_row_stab, top_left + 2] = F_one
                stabs[curr_row_stab, row_right + c + 1] = F_one
            else
                stabs[curr_row_stab, row_right - 2 * n + 1] = F_one
                stabs[curr_row_stab, row_right + 1] = F_one
            end
            stabs[curr_row_stab, row_right + c] = F_one
            if r != m
                stabs[curr_row_stab, row_right + n + 2 * (c - 1) + 1] = F_one
                stabs[curr_row_stab, row_right + n + 2 * (c - 1) + 2] = F_one
            else
                stabs[curr_row_stab, 2 * (c - 1) + 1] = F_one
                stabs[curr_row_stab, 2 * (c - 1) + 2] = F_one
            end
            curr_row_stab += 1

            # Z
            stabs[curr_row_stab, top_left + len] = F_one
            stabs[curr_row_stab, top_left + 1 + len] = F_one
            stabs[curr_row_stab, row_right + c + len] = F_one
            if r == m && c == n
                stabs[curr_row_stab, 1 + len] = F_one
                stabs[curr_row_stab, 2 * n + len] = F_one
                stabs[curr_row_stab, row_right + 1 + len] = F_one
            elseif r == m
                stabs[curr_row_stab, row_right + c + 1 + len] = F_one
                stabs[curr_row_stab, 2 * (c - 1) + 2 + len] = F_one
                stabs[curr_row_stab, 2 * (c - 1) + 3 + len] = F_one
            elseif c != n
                stabs[curr_row_stab, row_right + c + 1 + len] = F_one
                stabs[curr_row_stab, row_right + n + 2 * (c - 1) + 2 + len] = F_one
                stabs[curr_row_stab, row_right + n + 2 * (c - 1) + 3 + len] = F_one
            else
                stabs[curr_row_stab, row_right + 1 + len] = F_one
                stabs[curr_row_stab, row_right + n + 1 + len] = F_one
                stabs[curr_row_stab, row_right + n + 2 * (c - 1) + 2 + len] = F_one
            end
            curr_row_stab += 1
            top_left += 2
        end
    end
    
    logs = zero_matrix(F, 4, 2 * len)
    # top row is an X for pair one and a Z for pair two
    for c in 1:2 * n
        logs[1, c] = F_one
        logs[4, c + len] = F_one
    end
    # left column is a Z for pair one and an X for pair two
    for r in 1:m
        top_left = 3 * n * (r - 1) + 1
        logs[2, top_left + len] = F_one
        logs[2, top_left + 2 * n + len] = F_one
        logs[3, top_left] = F_one
        logs[3, top_left + 2 * n] = F_one
    end
    
    S = SubsystemCode(gauges, logs_alg = :sys_eqs)
    S.k == 2 || error("Got wrong dimension for periodic case.")
    set_logicals!(S, logs)

    dist = min(m, n)
    S.cache[:d_dressed] = dist
    S.cache[:l_bound_dressed] = dist
    S.cache[:u_bound_dressed] = dist
    
    S.cache[:dx_dressed] = dist
    S.cache[:dz_dressed] = dist
    return S
end

"""
    SubsystemToricCode(d::Int)

Return the subsystem toric code on a square lattice with `d` rows and `d`
columns of squares.
"""
SubsystemToricCode(d::Int) = SubsystemToricCode(d, d)

# Bravyi et al, "Subsystem surface codes with three-qubit check operators", (2013)
"""
    SubsystemSurfaceCode(m::Int, n::Int)

Return the subsystem surface code on a rectangular lattice with `m` rows and `n`
columns of squares.
"""
function SubsystemSurfaceCode(m::Int, n::Int)
    (2 <= m && 2 <= n) || throw(DomainError("Lattice dimensions must be at least two"))

    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)
    len = (3n + 2) * m + 2 * n + 1
    gauges = zero_matrix(F, 4 * m * n + 2 * n + 2 * m, 2 * len)
    stabs = zero_matrix(F, 2 * m * n + 2 * n + 2 * m, 2 * len)
    curr_row_stab = 1
    curr_row_gauge = 1
    for r in 1:m
        top_left = (3 * n + 2) * (r - 1) + 1
        row_right = top_left + 2 * n
        for c in 1:n
            # top left - Z
            gauges[curr_row_gauge, top_left + len] = F_one
            gauges[curr_row_gauge, top_left + 1 + len] = F_one
            gauges[curr_row_gauge, row_right + c + len] = F_one
            curr_row_gauge += 1

            # top right - X
            gauges[curr_row_gauge, top_left + 1] = F_one
            gauges[curr_row_gauge, top_left + 2] = F_one
            gauges[curr_row_gauge, row_right + c + 1] = F_one
            curr_row_gauge += 1

            # bottom left - X
            gauges[curr_row_gauge, row_right + c] = F_one
            gauges[curr_row_gauge, row_right + n + 2 * (c - 1) + 2] = F_one
            gauges[curr_row_gauge, row_right + n + 2 * (c - 1) + 3] = F_one
            curr_row_gauge += 1

            # bottom right - Z
            gauges[curr_row_gauge, row_right + c + 1 + len] = F_one
            gauges[curr_row_gauge, row_right + n + 2 * (c - 1) + 3 + len] = F_one
            gauges[curr_row_gauge, row_right + n + 2 * (c - 1) + 4 + len] = F_one
            curr_row_gauge += 1

            # X
            stabs[curr_row_stab, top_left + 1] = F_one
            stabs[curr_row_stab, top_left + 2] = F_one
            stabs[curr_row_stab, row_right + c + 1] = F_one
            stabs[curr_row_stab, row_right + c] = F_one
            stabs[curr_row_stab, row_right + n + 2 * (c - 1) + 2] = F_one
            stabs[curr_row_stab, row_right + n + 2 * (c - 1) + 3] = F_one
            curr_row_stab += 1

            # Z
            stabs[curr_row_stab, top_left + len] = F_one
            stabs[curr_row_stab, top_left + 1 + len] = F_one
            stabs[curr_row_stab, row_right + c + len] = F_one
            stabs[curr_row_stab, row_right + c + 1 + len] = F_one
            stabs[curr_row_stab, row_right + n + 2 * (c - 1) + 3 + len] = F_one
            stabs[curr_row_stab, row_right + n + 2 * (c - 1) + 4 + len] = F_one
            curr_row_stab += 1
            top_left += 2
        end
    end

    # left - X
    for r in 1:m
        top_left = (3 * n + 2) * (r - 1) + 1
        shift = 2 * n + 1
        stabs[curr_row_stab, top_left] = F_one
        stabs[curr_row_stab, top_left + shift] = F_one
        curr_row_stab += 1
        gauges[curr_row_gauge, top_left] = F_one
        gauges[curr_row_gauge, top_left + shift] = F_one
        curr_row_gauge += 1
    end

    # right - X
    for r in 1:m
        top_left = (3 * n + 2) * (r - 1) + 1
        row_right = top_left + 2 * n
        stabs[curr_row_stab, row_right + n + 1] = F_one
        stabs[curr_row_stab, row_right + 3 * n + 2] = F_one
        curr_row_stab += 1
        gauges[curr_row_gauge, row_right + n + 1] = F_one
        gauges[curr_row_gauge, row_right + 3 * n + 2] = F_one
        curr_row_gauge += 1
    end

    # top - Z
    top_left = 1
    for c in 1:n
        stabs[curr_row_stab, top_left + 1 + len] = F_one
        stabs[curr_row_stab, top_left + 2 + len] = F_one
        curr_row_stab += 1
        gauges[curr_row_gauge, top_left + 1 + len] = F_one
        gauges[curr_row_gauge, top_left + 2 + len] = F_one
        curr_row_gauge += 1
        top_left += 2
    end

    # bottom - Z
    bottom_left = len - 2 * n
    for c in 1:n
        stabs[curr_row_stab, bottom_left + len] = F_one
        stabs[curr_row_stab, bottom_left + 1 + len] = F_one
        curr_row_stab += 1
        gauges[curr_row_gauge, bottom_left + len] = F_one
        gauges[curr_row_gauge, bottom_left + 1 + len] = F_one
        curr_row_gauge += 1
        bottom_left += 2
    end

    logs = zero_matrix(F, 2, 2 * len)
    # top row is a logical X
    for c in 1:2 * n + 1
        logs[1, c] = F_one
    end
    # left column is a logical Z 
    for r in 1:m
        logs[2, (3 * n + 2) * (r - 1) + 1 + len] = F_one
        logs[2, (3 * n + 2) * (r - 1) + 2 * n + 2 + len] = F_one
    end
    logs[2, 2 * len -  2 * n] = F_one
    
    S = SubsystemCode(gauges, logs_alg = :sys_eqs)
    S.k == 1 || error("Got wrong dimension for non-periodic case.")
    set_logicals!(S, logs)
    set_minimum_distance!(S, minimum([m, n]))
    
    # Route dressed bounds directly into the cache, bypassing direct struct fields
    S.cache[:dx_dressed] = 2 * n + 1
    S.cache[:dz_dressed] = 2 * m + 1
    return S
end

"""
    SubsystemSurfaceCode(d::Int)

Return the subsystem surface code on a square lattice with `d` rows and `d`
columns of squares.
"""
SubsystemSurfaceCode(d::Int) = SubsystemSurfaceCode(d, d)

# TODO: charvec in all of these

#############################
      # Stabilizer codes
#############################

#############################
        # Misc codes
#############################

"""
    QC6()

Return the `C_6` stabilizer code defined by Knill.
"""
function QC6()
    S = StabilizerCode(["XIIXXX", "XXXIIX", "ZIIZZZ", "ZZZIIZ"])
    set_minimum_distance!(S, 2)
    return S
end


"""
    FiveQubitCode()
    Q513()

Return the `[[5, 1, 3]]` perfect qubit stabilizer code.
"""
function FiveQubitCode()
    S = StabilizerCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"])
    _set_SL_weight_enumerator!(
        S, Dict(0 => BigInt(1), 4 => BigInt(15)))
    set_minimum_distance!(S, 3)
    return S
end
Q513() = FiveQubitCode()
# is a perfect code

# should also do a test for other CSS construction via Hamming code and actually make that one default
"""
    Steanecode()
    Q713()

Return the `[[7, 1, 3]]` Steane code.
"""
function SteaneCode()
    S = CSSCode(["XXXXIII", "XXIIXXI", "XIXIXIX", "ZZZZIII", "ZZIIZZI", "ZIZIZIZ"])
    _set_SL_weight_enumerator!(
        S, Dict(0 => BigInt(1), 4 => BigInt(21), 6 => BigInt(42)))
    set_minimum_distance!(S, 3)
    return S
end
Q713() = SteaneCode()
# _SteaneCodeTrellis() = CSSCode(["XXIXXII", "IXXIXXI", "IIIXXXX", "ZZIZZII", "IZZIZZI", "IIIZZZZ"])
# also ZZIZZII, ZIZZIZI, IZZZIIZ, XXIXXII, XIXXIXI, IXXXIIX

"""
    Shorcode()
    Q913()

Return the `[[9, 1, 3]]` Shor code.
"""
function ShorCode()
    S = CSSCode(["ZZIIIIIII", "IZZIIIIII", "IIIZZIIII", "IIIIZZIII", "IIIIIIZZI", "IIIIIIIZZ",
    "XXXXXXIII", "IIIXXXXXX"])
    set_minimum_distance!(S, 3)
    return S
end
Q913() = ShorCode()

function Q412()
    S = CSSCode(["XXXX", "ZZII", "IIZZ"])
    set_minimum_distance!(S, 2)
    return S
end

"""
    QC4()
    Q422()

Return the `[[4, 2, 2]]`` stabilizer code `C_4`` defined by Knill.
"""
function Q422()
    S = CSSCode(["XXXX", "ZZZZ"])
    set_minimum_distance!(S, 2)
    return S
end
QC4() = Q422()

function Q511()
    S = StabilizerCode(["ZXIII", "XZXII", "IXZXI", "IIXZX"])
    set_minimum_distance!(S, 1)
    return S
end

function Q823()
    F = Oscar.Nemo.Native.GF(2)
    stabs = matrix(F, [1 0 0 0 1 0 0 0 1 1 1 1 0 0 0 0;
    0 0 0 1 0 1 0 0 1 0 0 0 0 1 0 0;
    0 1 0 0 1 1 1 0 0 0 1 1 1 0 1 0;
    0 0 1 0 1 1 1 0 0 1 1 0 1 1 0 0;
    0 0 1 1 1 0 1 0 0 0 0 1 0 1 1 1;
    0 0 0 0 0 0 1 1 0 0 1 0 0 0 1 0]);
    S = StabilizerCode(stabs)
    set_minimum_distance!(S, 3)
    return S
end

function Q832()
    S = StabilizerCode(["XXXXXXXX", "ZZZZIIII", "ZZIIZZII", "ZIZIZIZI", "ZZZZZZZZ"])
    set_minimum_distance!(S, 2)
    return S
end
SmallestInterestingColorCode() = Q832()

"""
    Q15RM()
    Q1513()

Return the `[[15, 1, 3]]` quantum Reed-Muller code.
"""
function Q15RM()
    S = StabilizerCode(["ZIZIZIZIZIZIZIZ", "IZZIIZZIIZZIIZZ", "IIIZZZZIIIIZZZZ", "IIIIIIIZZZZZZZZ",
    "IIZIIIZIIIZIIIZ", "IIIIZIZIIIIIZIZ", "IIIIIZZIIIIIIZZ", "IIIIIIIIIZZIIZZ", "IIIIIIIIIIIZZZZ",
    "IIIIIIIIZIZIZIZ", "XIXIXIXIXIXIXIX", "IXXIIXXIIXXIIXX", "IIIXXXXIIIIXXXX", "IIIIIIIXXXXXXXX"])
    set_minimum_distance!(S, 3)
    return S
end
Q1513() = Q15RM()

"""
    Q1573()

Return the `[[15, 7, 3]]` quantum Hamming code.
"""
function Q1573()
    S = StabilizerCode(["IIIIIIIXXXXXXXX", "IIIXXXXIIIIXXXX", "IXXIIXXIIXXIIXX", "XIXIXIXIXIXIXIX",
    "IIIIIIIZZZZZZZZ", "IIIZZZZIIIIZZZZ", "IZZIIZZIIZZIIZZ", "ZIZIZIZIZIZIZIZ"])
    # one can use a basis for this such that the first logical pair is transversal X, Z
    set_minimum_distance!(S, 3)
    return S
end

"""
    GrossCode()

Return the `[[144, 12, 12]]` gross code.
"""
function GrossCode()
    S, (x, y) = polynomial_ring(Oscar.Nemo.Native.GF(2), [:x, :y])
    l = 12
    m = 6
    R, _ = quo(S, ideal(S, [x^l - 1, y^m - 1]))
    a = R(x^3 + y + y^2)
    b = R(y^3 + x + x^2)
    code = BBCode(a, b)
    set_minimum_distance!(code, 12)
    return code
end

"""
    QuantumRepetitionCode(d::Int; error_type::Symbol = :phase_flip)

Return the `[[d, 1, d]]` quantum repetition code. 
`error_type` can be `:phase_flip` (X stabilizers) or `:bit_flip` (Z stabilizers).
"""
function QuantumRepetitionCode(d::Int; error_type::Symbol = :phase_flip)
    (d >= 3) || throw(DomainError("Distance must be at least 3."))
    
    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)
    stabs = zero_matrix(F, d - 1, 2 * d)
    
    for i in 1:d-1
        if error_type == :phase_flip
            stabs[i, i] = F_one
            stabs[i, i + 1] = F_one
        elseif error_type == :bit_flip
            stabs[i, i + d] = F_one
            stabs[i, i + 1 + d] = F_one
        else
            throw(ArgumentError("error_type must be :phase_flip or :bit_flip"))
        end
    end
    
    S = StabilizerCode(stabs)
    
    if error_type == :phase_flip
        S.cache[:dx] = d
        S.cache[:dz] = 1
    else
        S.cache[:dx] = 1
        S.cache[:dz] = d
    end
    S.cache[:d] = 1
    
    return S
end

"""
    QuantumGolayCode()

Return the `[[23, 1, 7]]` quantum Golay CSS code.
"""
function QuantumGolayCode()
    F = Oscar.Nemo.Native.GF(2)
    S, x = polynomial_ring(F, :x)
    R, _ = residue_ring(S, x^23 - 1)
    
    # Generator polynomial for the classical [23, 12, 7] Golay code
    g = R(x^11 + x^9 + x^7 + x^6 + x^5 + x + 1)
    H_circ = residue_polynomial_to_circulant_matrix(g)
    
    # The circulant matrix is 23x23 but rank 11. We reduce it to get the clean parity check matrix.
    _, H_rref = rref(H_circ)
    H_clean = _remove_empty(H_rref, :rows)
    
    # CSS construction requires H_X * transpose(H_Z) = 0. 
    # Because G_23^⟂ ⊆ G_23, H_clean * transpose(H_clean) = 0 is satisfied automatically.
    S = CSSCode(H_clean, H_clean)
    
    # Set explicit boundaries
    S.cache[:d] = 7
    S.cache[:l_bound] = 7
    S.cache[:u_bound] = 7
    
    S.cache[:dx] = 7
    S.cache[:l_bound_dx] = 7
    S.cache[:u_bound_dx] = 7
    
    S.cache[:dz] = 7
    S.cache[:l_bound_dz] = 7
    S.cache[:u_bound_dz] = 7
    
    return S
end

#############################
 # Triangular Surface Codes
#############################

function _triangular_lattice(L::Int)
    # 0 - vertical
    # 1 - horizontal
    # 2 - diagonal
    numbering = zeros(Int, L, L, 3)
    num = 1
    for i in 1:L
        for j in 1:L
            for k in 1:3
                numbering[i, j, k] = num
                num += 1
            end
        end
    end
    return numbering
end

function _triangular_lattice_X_stabilizers(L::Int, numbering::Array{Int64, 3})
    F = Oscar.Nemo.Native.GF(2)
    stabilizers = zero_matrix(F, L^2, 3 * L^2)
    r = 1
    for i in 1:L
        for j in 1:L
            stabilizers[r, numbering[i, j, 1]] = 1
            stabilizers[r, numbering[i, j, 2]] = 1
            stabilizers[r, numbering[i, j, 3]] = 1
            if i == 1
                stabilizers[r, numbering[end, j, 1]] = 1
            else
                stabilizers[r, numbering[i - 1, j, 1]] = 1
            end
            if j == 1
                stabilizers[r, numbering[i, end, 2]] = 1
            else
                stabilizers[r, numbering[i, j - 1, 2]] = 1
            end
            if i == 1 && j == 1
                stabilizers[r, numbering[end, end, 3]] = 1
            elseif i != 1 && j == 1
                stabilizers[r, numbering[i - 1, end, 3]] = 1
            elseif i == 1 && j != 1
                stabilizers[r, numbering[end, j - 1, 3]] = 1
            else
                stabilizers[r, numbering[i - 1, j - 1, 3]] = 1
            end
            r += 1
        end
    end
    return hcat(stabilizers, zero_matrix(F, L^2, 3 * L^2))
end

function _triangular_lattice_Z_stabilizers(L::Int, numbering::Array{Int64, 3})
    F = Oscar.Nemo.Native.GF(2)
    stabilizers = zero_matrix(F, 2 * L^2, 3 * L^2)
    r = 1
    for i in 1:L
        for j in 1:L
            stabilizers[r, numbering[i, j, 2]] = 1
            stabilizers[r, numbering[i, j, 3]] = 1
            if j == L
                stabilizers[r, numbering[i, 1, 1]] = 1
            else
                stabilizers[r, numbering[i, j + 1, 1]] = 1
            end

            r += 1
            stabilizers[r, numbering[i, j, 1]] = 1
            stabilizers[r, numbering[i, j, 3]] = 1
            if i == L
                stabilizers[r, numbering[1, j, 2]] = 1
            else
                stabilizers[r, numbering[i + 1, j, 2]] = 1
            end
            r += 1
        end
    end
    return hcat(zero_matrix(F, 2 * L^2, 3 * L^2), stabilizers)
end

function _triangular_lattice_X_logicals(L::Int, numbering::Array{Int64, 3})
    # should be 0110110110
    z = zeros(UInt8, 3 * L^2)
    logical1 = zeros(UInt8, 3 * L^2)
    for j in 0:L - 1
        for k in 0:2
            if k == 0
                logical1[3 * j + k + 1] = 0x01
            elseif k == 1
                logical1[3 * j + k + 1] = 0x00
            else
                logical1[3 * j + k + 1] = 0x01
            end
        end
    end
    logical1 = [logical1; z]

    logical2 = zeros(UInt8, 3 * L^2)
    for j in 1:L
        for k in 0:2
            if k == 0
                logical2[3 * j + k + 1] = 0x01
            elseif k == 1
                logical2[3 * j + k + 1] = 0x00
            else
                logical2[3 * j + k + 1] = 0x01
            end
        end
    end
    logical2 = [logical2; z]
    return [logical1, logical2]
end

function _triangular_lattice_Z_logicals(L::Int, numbering::Array{Int64, 3}, symp::Bool=true)
    # should be 1001001001
    x = zeros(UInt8, 3 * L^2)
    logical1 = zeros(UInt8, 3 * L^2)
    for j in 0:L - 1
        for k in 0:2
            if k == 0
                logical1[3 * j + k + 1] = 0x01
            elseif k == 1
                logical1[3 * j + k + 1] = 0x00
            else
                logical1[3 * j + k + 1] = 0x01
            end
        end
    end
    logical1 = [z; logical1]

    logical2 = zeros(UInt8, 3 * L^2)
    for j in 1:L
        for k in 0:2
            if k == 0
                logical2[3 * j + k + 1] = 0x01
            elseif k == 1
                logical2[3 * j + k + 1] = 0x00
            else
                logical2[3 * j + k + 1] = 0x01
            end
        end
    end
    logical2 = [z; logical2]
    return [logical1, logical2]
end

function TriangularSurfaceCode(L::Int)
    numbering = _triangular_lattice(L)
    X_stabs = _triangular_lattice_X_stabilizers(L, numbering)
    # println(rank(X_stabs))
    Z_stabs = _triangular_lattice_Z_stabilizers(L, numbering)
    # println(Z_stabs)
    # logicals = [_triangular_lattice_X_logicals(L, numbering), _triangular_lattice_Z_logicals(L, numbering)]
    # TODO distances
    return CSSCode(X_stabs[1:end - 1, :], Z_stabs[1:end - 1, :])
end

#############################
   # Rotated Surface Codes
#############################

function _R_Surf_stabs(d::Int)
    n = d^2
    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)
    S = zero_matrix(F, n - 1, 2 * n)
    row = 1

    # X stabilizers (Interior)
    for r in 1:d-1
        for c in 1:d-1
            if iseven(r + c)
                S[row, (r-1)*d + c] = F_one
                S[row, (r-1)*d + c + 1] = F_one
                S[row, r*d + c] = F_one
                S[row, r*d + c + 1] = F_one
                row += 1
            end
        end
    end

    # Top X boundary (r = 1)
    for c in 1:d-1
        if iseven(0 + c)
            S[row, c] = F_one
            S[row, c + 1] = F_one
            row += 1
        end
    end

    # Bottom X boundary (r = d)
    for c in 1:d-1
        if iseven(d + c)
            S[row, (d-1)*d + c] = F_one
            S[row, (d-1)*d + c + 1] = F_one
            row += 1
        end
    end

    # Z stabilizers (Interior)
    for r in 1:d-1
        for c in 1:d-1
            if isodd(r + c)
                S[row, n + (r-1)*d + c] = F_one
                S[row, n + (r-1)*d + c + 1] = F_one
                S[row, n + r*d + c] = F_one
                S[row, n + r*d + c + 1] = F_one
                row += 1
            end
        end
    end

    # Left Z boundary (c = 1)
    for r in 1:d-1
        if isodd(r + 0)
            S[row, n + (r-1)*d + 1] = F_one
            S[row, n + r*d + 1] = F_one
            row += 1
        end
    end

    # Right Z boundary (c = d)
    for r in 1:d-1
        if isodd(r + d)
            S[row, n + (r-1)*d + d] = F_one
            S[row, n + r*d + d] = F_one
            row += 1
        end
    end

    return S
end

function _R_Surf_logs(F::CTFieldTypes, d::Int)
    n = d^2
    F_one = F(1)
    logs = zero_matrix(F, 2, 2 * n)
    
    # Logical X: Rightmost column (connects Top and Bottom X boundaries)
    for r in 1:d
        logs[1, (r-1)*d + d] = F_one
    end
    
    # Logical Z: Top row (connects Left and Right Z boundaries)
    for c in 1:d
        logs[2, n + c] = F_one
    end

    return logs
end

"""
    RotatedSurfaceCode(d::Int)

Return the `[[d^2, 1, d]]` rotated surface code.

This is the surface-13/17 configuration found in "Low-distance surface codes under realistic quantum noise"
by Tomita and Svore. The standard planar surface code is equivalent to their surface-25 configuration, which
can be seen by viewing the stabilizers of PlanarSurfaceCode as an adjacency matrix.
"""
function RotatedSurfaceCode(d::Int)
    d >= 3 || throw(DomainError("Current implementation requires d ≥ 3."))

    stabs = _R_Surf_stabs(d)
    logs = _R_Surf_logs(base_ring(stabs), d)
    
    S = StabilizerCode(stabs)
    set_logicals!(S, logs)
    
    # Cache boundaries safely
    S.cache[:dx] = d
    S.cache[:l_bound_dx] = d
    S.cache[:u_bound_dx] = d
    
    S.cache[:dz] = d
    S.cache[:l_bound_dz] = d
    S.cache[:u_bound_dz] = d
    
    S.cache[:d] = d
    S.cache[:l_bound] = d
    S.cache[:u_bound] = d
    
    return S
end

#############################
     # XZZX Surface Codes
#############################

function _XZZX_stabs_logs(d::Int)
    n = d^2
    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)
    S = zero_matrix(F, n - 1, 2 * n)
    row = 1

    # Interior XZZX stabilizers
    for r in 1:d-1
        for c in 1:d-1
            # top-left X, top-right Z, bottom-left Z, bottom-right X
            S[row, (r-1)*d + c] = F_one
            S[row, n + (r-1)*d + c + 1] = F_one
            S[row, n + r*d + c] = F_one
            S[row, r*d + c + 1] = F_one
            row += 1
        end
    end

    # Top boundary (r = 1) -> cut bottom half -> leaves Z, X
    for c in 1:d-1
        if iseven(c)
            S[row, n + c] = F_one       # Z
            S[row, c + 1] = F_one       # X
            row += 1
        end
    end

    # Bottom boundary (r = d) -> cut top half -> leaves X, Z
    for c in 1:d-1
        if iseven(d + c)
            S[row, (d-1)*d + c] = F_one       # X
            S[row, n + (d-1)*d + c + 1] = F_one   # Z
            row += 1
        end
    end

    # Left boundary (c = 1) -> cut right half -> leaves Z, X
    for r in 1:d-1
        if isodd(r)
            S[row, n + (r-1)*d + 1] = F_one   # Z
            S[row, r*d + 1] = F_one           # X
            row += 1
        end
    end

    # Right boundary (c = d) -> cut left half -> leaves X, Z
    for r in 1:d-1
        if isodd(r + d)
            S[row, (r-1)*d + d] = F_one       # X
            S[row, n + r*d + d] = F_one       # Z
            row += 1
        end
    end

    logs = zero_matrix(F, 2, 2 * n)
    # Logical X: Alternates X, Z, X along right column
    for r in 1:d
        if isodd(r)
            logs[1, (r-1)*d + d] = F_one # X
        else
            logs[1, n + (r-1)*d + d] = F_one # Z
        end
    end

    # Logical Z: Alternates Z, X, Z along top row
    for c in 1:d
        if isodd(c)
            logs[2, n + c] = F_one # Z
        else
            logs[2, c] = F_one # X
        end
    end

    return S, logs
end

"""
    XZZXSurfaceCode(d::Int)

Return the `[[d^2, 1, d]]` XZZX surface code.
"""
function XZZXSurfaceCode(d::Int)
    d >= 3 || throw(DomainError("Current implementation requires d ≥ 3."))

    stabs, logs = _XZZX_stabs_logs(d)
    
    S = StabilizerCode(stabs)
    set_logicals!(S, logs)
    
    S.cache[:d] = d
    S.cache[:l_bound] = d
    S.cache[:u_bound] = d
    
    return S
end

################################
 # Triangular Color Codes 6.6.6
################################

"""
    TriangularColorCode666(d::Int)

Return the 6.6.6 triangular color code of distance `d`.
Generates the lattice programmatically via 2D coordinate mapping.
"""
function _TriangularColorCode666_generated(d::Int)
    (d >= 3 && isodd(d)) || throw(DomainError("Distance must be an odd integer >= 3."))
    
    # 1. Generate Qubit Coordinates
    # We use a skewed 2D coordinate system where qubits lie on a triangular grid boundary.
    qubits = Tuple{Int, Int}[]
    for y in 0:(d-1)
        # The row length shrinks as we go up the triangle
        row_width = d - y
        offset = cld(y, 2)
        for x in offset:(offset + row_width - 1)
            push!(qubits, (x, y))
        end
    end
    
    n = length(qubits)
    q_idx = Dict(q => i for (i, q) in enumerate(qubits))
    
    # 2. Generate Stabilizer Faces
    # A face in the 6.6.6 lattice is defined by the 3 adjacent qubits in a downward triangle
    # plus the 3 qubits in the upward triangle adjacent to it.
    F = Oscar.Nemo.Native.GF(2)
    stabs = zero_matrix(F, 0, 2 * n)
    
    # Sweep the grid for downward pointing triangles that form the centers of the hexagons
    for y in 0:(d-2)
        for x in 0:(d-1)
            # Center of a face is defined by a downward triangle of qubits
            q1, q2, q3 = (x, y), (x+1, y), (x, y+1)
            
            # If the core triangle exists in our qubit dictionary, it's a valid face
            if haskey(q_idx, q1) && haskey(q_idx, q2) && haskey(q_idx, q3)
                row = zero_matrix(F, 1, 2 * n)
                
                # The 6 possible vertices of the hexagon
                hex_vertices = [
                    (x, y), (x+1, y), (x, y+1),       # Core triangle
                    (x-1, y), (x+1, y-1), (x-1, y+1)  # Expanded hexagon points
                ]
                
                weight = 0
                for v in hex_vertices
                    if haskey(q_idx, v)
                        row[1, q_idx[v]] = F(1)
                        row[1, q_idx[v] + n] = F(1) # Color codes are CSS (H_X = H_Z)
                        weight += 1
                    end
                end
                
                # Valid boundary faces in 6.6.6 are weight 4 or weight 6
                if weight == 4 || weight == 6
                    stabs = vcat(stabs, row)
                end
            end
        end
    end
    
    S = StabilizerCode(stabs)
    S.cache[:d] = d
    S.cache[:l_bound] = d
    S.cache[:u_bound] = d
    return S
end

################################
 # Triangular Color Codes 4.8.8
################################

"""
    TriangularColorCode488(d::Int)

Return the 4.8.8 triangular color code of distance `d`.
Generates the lattice programmatically via 2D coordinate mapping.
"""
function _TriangularColorCode488_generated(d::Int)
    (d >= 3 && isodd(d)) || throw(DomainError("Distance must be an odd integer >= 3."))
    
    # 4.8.8 maps to an inflated square grid. 
    # We define qubits at the centers of the edges of a standard grid.
    qubits = Tuple{Int, Int, Symbol}[]
    
    for y in 0:(d-1)
        for x in 0:(d-1 - y)
            # Horizontal edges
            if x < d - 1 - y
                push!(qubits, (x, y, :H))
            end
            # Vertical edges
            if y < d - 1 - x
                push!(qubits, (x, y, :V))
            end
        end
    end
    
    n = length(qubits)
    q_idx = Dict(q => i for (i, q) in enumerate(qubits))
    
    F = Oscar.Nemo.Native.GF(2)
    stabs = zero_matrix(F, 0, 2 * n)
    
    # Generate Square Faces (Weight 4)
    for y in 0:(d-2)
        for x in 0:(d-2 - y)
            row = zero_matrix(F, 1, 2 * n)
            face_qs = [(x, y, :H), (x, y, :V), (x+1, y, :V), (x, y+1, :H)]
            for q in face_qs
                if haskey(q_idx, q)
                    row[1, q_idx[q]] = F(1)
                    row[1, q_idx[q] + n] = F(1)
                end
            end
            stabs = vcat(stabs, row)
        end
    end
    
    # Generate Octagon Faces (Weight 8)
    for y in 0:(d-2)
        for x in 0:(d-2 - y)
            # Octagons sit between the squares
            if iseven(x + y) 
                row = zero_matrix(F, 1, 2 * n)
                # 8 surrounding edges
                oct_qs = [
                    (x, y, :H), (x+1, y, :H), 
                    (x, y, :V), (x, y+1, :V),
                    (x-1, y+1, :H), (x, y+1, :H),
                    (x+1, y-1, :V), (x+1, y, :V)
                ]
                
                weight = 0
                for q in oct_qs
                    if haskey(q_idx, q)
                        row[1, q_idx[q]] = F(1)
                        row[1, q_idx[q] + n] = F(1)
                        weight += 1
                    end
                end
                
                # Valid boundary octagons are weight 4, 6, or 8
                if weight >= 4
                    stabs = vcat(stabs, row)
                end
            end
        end
    end
    
    S = StabilizerCode(stabs)
    S.cache[:d] = d
    S.cache[:l_bound] = d
    S.cache[:u_bound] = d
    return S
end

################################
         # Toric Codes
################################

"""
    ToricCode(d::Int)

Return the `[[2d^2, 2, d]]` toric code.

The lattice orientation used here follows the picture at https://errorcorrectionzoo.org/c/surface.
"""
function ToricCode(d::Int)
    2 <= d || throw(DomainError("Distance must be at least two."))

    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)
    A = zero_matrix(F, d^2, 2 * d^2) # stars, X stabilizers
    B = zero_matrix(F, d^2, 2 * d^2) # faces, Z stabilizers
    qubit = 1
    row_A = 1
    row_B = 1
    for r in 1:2 * d
        if isodd(r)
            for c in 1:d
                # println("r = $r, c = $c, row_A = $row_A")
                if r != 2 * d - 1 && c != d
                    A[row_A, qubit] = A[row_A, qubit + d] = A[row_A, qubit + d + 1] = A[row_A, qubit + 2 * d] = F_one
                elseif r == 2 * d - 1 && c != d
                    A[row_A, qubit] = A[row_A, qubit + d] = A[row_A, qubit + d + 1] = A[row_A, c] = F_one
                elseif r != 2 * d - 1 && c == d
                    A[row_A, qubit] = A[row_A, qubit + d] = A[row_A, qubit + 1] = A[row_A, qubit + 2 * d] = F_one
                elseif r == 2 * d - 1 && c == d
                    A[row_A, qubit] = A[row_A, qubit + d] = A[row_A, qubit + 1] = A[row_A, c] = F_one
                else
                    error("Ran into unaccounted for case in creating the toric code lattice.")
                end
                row_A += 1
                qubit += 1
            end
        else
            for c in 1:d
                # println("r = $r, c = $c, row_B = $row_B")
                if r != 2 * d && c == 1
                    B[row_B, qubit] = B[row_B, qubit + d] = B[row_B, qubit + 2 * d] = B[row_B, qubit + 2 * d - 1] = F_one
                elseif r != 2 * d && c != 1
                    B[row_B, qubit] = B[row_B, qubit + d - 1] = B[row_B, qubit + d] = B[row_B, qubit + 2 * d] = F_one
                elseif r == 2 * d && c == 1
                    B[row_B, qubit] = B[row_B, d] = B[row_B, d + 1] = B[row_B, 1] = F_one
                elseif r == 2 * d && c != 1
                    B[row_B, qubit] = B[row_B, c - 1] = B[row_B, c] = B[row_B, c + d] = F_one
                else
                    error("Ran into unaccounted for case in creating the toric code lattice.")
                end
                row_B += 1
                qubit += 1
            end
        end
    end
    S = CSSCode(A, B)

    Z1 = zero_matrix(S.F, 1, 4 * d^2)
    for c in 1:d
        Z1[1, c + d + S.n] = F_one
    end
    X1 = zero_matrix(S.F, 1, 4 * d^2)
    for r in 1:2:2 * d
        X1[1, r * d + 1] = F_one
    end
    Z2 = zero_matrix(S.F, 1, 4 * d^2)
    for r in 1:2:2 * d
        Z2[1, (r - 1) * d + 1 + S.n] = F_one
    end
    X2 = zero_matrix(S.F, 1, 4 * d^2)
    for c in 1:d
        X2[1, c] = F_one
    end
    set_logicals!(S, vcat(X1, Z1, X2, Z2))
    set_X_minimum_distance!(S, d)
    set_Z_minimum_distance!(S, d)
    return S
end

################################
     # Planar Surface Codes
################################

"""
    PlanarSurfaceCode(d_x::Int, d_z::Int)
    PlanarSurfaceCode(d::Int)

Return the `[[d_x * d_z + (d_x - 1) * (d_z - 1), 1, d_x/d_z]]` planar surface code.

The top and bottom boundaries are "smooth" (`Z`) and the left and right are "rough" (`X`).
"""
function PlanarSurfaceCode(d_x::Int, d_z::Int)
    (2 <= d_x && 2 <= d_z) || throw(DomainError("Distances must be at least two."))

    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)
    num_V = d_x * d_z + (d_x - 1) * (d_z - 1)
    A = zero_matrix(F, d_x * (d_z - 1) + 1, num_V) # stars, X stabilizers
    B = zero_matrix(F, d_z * (d_x - 1), num_V) # faces, Z stabilizers
    qubit = 1
    row_A = 1
    row_B = 1
    for r in 1:d_z
        for c in 1:d_x
            if r != d_z
                if c == 1
                    B[row_B, qubit] = B[row_B, qubit + d_x] = B[row_B, qubit + 2 * d_x - 1] = F_one
                    row_B += 1
                elseif c == d_x
                    B[row_B, qubit] = B[row_B, qubit + d_x - 1] = B[row_B, qubit + 2 * d_x - 1] = F_one
                    row_B += 1
                else
                    B[row_B, qubit] = B[row_B, qubit + d_x - 1] = B[row_B, qubit + d_x] = B[row_B, qubit + 2 * d_x - 1] = F_one
                    row_B += 1
                end
            end

            if c != d_x
                if r == 1
                    A[row_A, qubit] = A[row_A, qubit + 1] = A[row_A, qubit + d_x] = F_one
                    row_A += 1
                elseif r == d_z
                    A[row_A, qubit] = A[row_A, qubit + 1] = A[row_A, qubit - d_x + 1] = F_one
                    row_A += 1
                else
                    A[row_A, qubit] = A[row_A, qubit + 1] = A[row_A, qubit + d_x] = A[row_A, qubit - d_x + 1] = F_one
                    row_A += 1
                end
            end
            qubit += 1
        end
        qubit += d_x - 1
    end
    S = CSSCode(A, B)

    X1 = zero_matrix(S.F, 1, 2 * S.n)
    for r in 1:2:d_x
        X1[1, d_z * (r - 1) + (d_z - 1) * (r - 1) + 1] = F_one
    end
    Z1 = zero_matrix(S.F, 1, 2 * S.n)
    for c in 1:d_z
        Z1[1, c + S.n] = F_one
    end
    set_logicals!(S, vcat(X1, Z1))
    # Planar surface code is a Stabilizer Code, so use regular X/Z distance bounds
    set_X_minimum_distance!(S, d_x)
    set_Z_minimum_distance!(S, d_z)
    return S
end
PlanarSurfaceCode(d::Int) = PlanarSurfaceCode(d, d)

################################
     # 3D PlanarSurfaceCode
################################

"""
    PlanarSurfaceCode3D(d::Int)

Return the 3D planar surface code of distance `d`.
"""
PlanarSurfaceCode3D(d::Int) = _build_3D_surface_code(d, periodic=false)

################################
       # XY Surface Codes
################################

# TODO remove quadratic
"""
    XYSurfaceCode(d_x::Int, d_z::Int)
    XYSurfaceCode(d::Int)

Return the `[[d_x * d_y + (d_x - 1) * (d_y - 1), 1, d_x/d_y]]` XY surface code of
"Ultrahigh Error Threshold for Surface Codes with Biased Noise" by Tuckett, Bartlett, and Flammia.

The top and bottom boundaries are "smooth" (`Y`) and the left and right are "rough" (`X`).
"""
function XYSurfaceCode(d_x::Int, d_y::Int)
    (2 <= d_x && 2 <= d_y) || throw(DomainError("Distances must be at least two."))

    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)
    num_V = d_x * d_y + (d_x - 1) * (d_y - 1)
    M = zero_matrix(F, num_V - 1, 2 * num_V)
    qubit = 1
    row = 1
    for r in 1:d_y
        for c in 1:d_x
            if r != d_z
                if c == 1
                    M[row, qubit] = M[row, qubit + d_x] = M[row, qubit + 2 * d_x - 1] = F_one
                    M[row, qubit + num_V] = M[row, qubit + d_x + num_V] = M[row, qubit + 2 * d_x - 1 + num_V] = F_one
                    row += 1
                elseif c == d_x
                    M[row, qubit] = M[row, qubit + d_x - 1] = M[row, qubit + 2 * d_x - 1] = F_one
                    M[row, qubit + num_V] = M[row, qubit + d_x - 1 + num_V] = M[row, qubit + 2 * d_x - 1 + num_V] = F_one
                    row += 1
                else
                    M[row, qubit] = M[row, qubit + d_x - 1] = M[row, qubit + d_x] = M[row, qubit + 2 * d_x - 1] = F_one
                    M[row, qubit + num_V] = M[row, qubit + d_x - 1 + num_V] = M[row, qubit + d_x + num_V] = M[row, qubit + 2 * d_x - 1 + num_V] = F_one
                    row += 1
                end
            end

            if c != d_x
                if r == 1
                    M[row, qubit] = M[row, qubit + 1] = M[row, qubit + d_x] = F_one
                    row += 1
                elseif r == d_z
                    M[row, qubit] = M[row, qubit + 1] = M[row, qubit - d_x + 1] = F_one
                    row += 1
                else
                    M[row, qubit] = M[row, qubit + 1] = M[row, qubit + d_x] = M[row, qubit - d_x + 1] = F_one
                    row += 1
                end
            end
            qubit += 1
        end
        qubit += d_x - 1
    end
    S = StabilizerCode(M)
    # Eone = S.E(1)
    # ω = gen(S.E)
    # X1 = zero_matrix(S.E, 1, num_V)
    # for r in 1:2:d_x
    #     X1[1, d_z * (r - 1) + (d_z - 1) * (r - 1) + 1] = Eone
    # end
    # Z1 = zero_matrix(S.E, 1, num_V)
    # for c in 1:d_z
    #     Z1[1, c] = ω
    # end
    # set_logicals!(S, vcat(X1, Z1))
    set_X_minimum_distance!(S, d_x)
    set_Z_minimum_distance!(S, d_y)
    return S
end
XYSurfaceCode(d::Int) = XYSurfaceCode(d, d)

# ################################
#          # XYZ^2 Codes
# ################################

# """
#     XYZ2Code(d::Int)

# Return the `[[2d^2, 1, d]]` XYZ^2 (XYZXYZ) code of "The XYZ^2 hexagonal stabilizer code"
# by Srivastava, Kockum, and Granath.
# """
# function XYZ2Code(d::Int)
#     3 <= d && isodd(d) || throw(DomainError("The distance must be an odd, positive integer."))

#     E = GF(2, 2, :ω)
#     ω = gen(E)
#     Eone = E(1)
#     M = zero_matrix(E, 2 * d^2 - 1, 2 * d^2) # stars, X stabilizers
#     qubit = d + 1
#     row = 1
#     # comments refer to rotating Figure 1 of paper to the right by 45 degrees such that it's a rectangle
#     for r in 2:2 * d - 1
#         for c in 1:d
#             if isodd(r)
#                 # weight-3 stabilizers on bottom
#                 if r == 2 * d - 1 && isodd(c) && c != d
#                     M[row, qubit] = M[row, qubit + d] = M[row, qubit + d + 1] = Eone
#                     row += 1
#                 # weight-3 stabilizers on left
#                 elseif c == 1 && (r + 1) % 4 == 0 # r != 2 * d - 1 && - never need this since restricting to d odd
#                     M[row, qubit] = M[row, qubit - d] = M[row, qubit + d] = Eone
#                     row += 1
#                 end
#             else
#                 # full hex
#                 if c != d
#                     M[row, qubit] = M[row, qubit - d] = M[row, qubit + 1] =  M[row, qubit + d] = M[row, qubit + d + 1] = M[row, qubit + 2 * d + 1] = Eone
#                     row += 1
#                 end
#                 # weight-3 stabilizers on top
#                 if r == 2 && isodd(c) && c != 1
#                     M[row, qubit] = M[row, qubit - d] = M[row, qubit - d - 1] = Eone
#                     row += 1
#                 # weight-3 stabilizers on right
#                 elseif r != 2 && c == d && r % 4 == 0
#                     M[row, qubit] = M[row, qubit - d] = M[row, qubit + d] = Eone
#                     row += 1
#                 end
#             end
#             qubit += 1
#         end
#     end
#     display(M)
#     return
#     # Y distance is also 2 * d^2
    # set_dressed_X_minimum_distance!(S, d)
    # set_dressed_Z_minimum_distance!(S, 2 * d^2)
# end

################################
           # H Codes
################################

"""
    HCode(k::Int)

Return the `[[k + 4, k, 2]]` H code from `https://errorcorrectionzoo.org/c/quantum_h`.
"""
function HCode(k::Int)
    (2 <= k && iseven(k)) || throw(DomainError("Input must be >= 2 and even.")) 

    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)
    X = zero_matrix(F, 2, k + 4)
    Z = zero_matrix(F, 2, k + 4)
    X[1, 1] = X[1, 2] = X[1, 3] = X[1, 4] = F_one
    Z[1, 1] = Z[1, 2] = Z[1, 3] = Z[1, 4] = F_one
    X[2, 1] = X[2, 2] = F_one
    Z[2, 1] = Z[2, 2] = F_one
    for c in 5:k + 3
        X[2, c] = X[2, c + 1] = F_one
        Z[2, c] = Z[2, c + 1] = F_one
    end
    S = CSSCode(X, Z)
    set_minimum_distance!(S, 2)
    return S
end

#################################
        # 3D Toric codes
#################################

"""
    ToricCode3D(d::Int)

Return the 3D toric code of distance `d`.
"""
ToricCode3D(d::Int) = _build_3D_surface_code(d, periodic=true)

#################################
        # 4D Toric codes
#################################

@auto_hash_equals struct _Vertex
    coordinates::Vector{Int}
end

@auto_hash_equals struct _Cell
    vertices::Set{_Vertex}
end

"""
    Compute the n-cells of a d-dimensional hypercubic lattice with linear size l.
    Supports both periodic (closed) and open (planar) boundaries.
"""
function _compute_cells(l::Int, n::Int, d::Int = 4; periodic::Bool = true)
    cells = Set{_Cell}()    
    coords_to_change = collect(combinations(1:d, n))
    max_coord = periodic ? l - 1 : l
    
    for coord in Iterators.product([0:max_coord for i in 1:d]...)
        for directions in coords_to_change
            coord_vec = collect(coord)
            new_coords = Vector([copy(coord_vec) for _ in 1:2^n - 1])
            valid = true
            
            for i in 1:2^n - 1, j in 1:n
                if ((i >> (j - 1)) & 1) == 1
                    new_coords[i][directions[j]] += 1
                end
            end
            
            vertices = Set{_Vertex}()
            push!(vertices, _Vertex(coord_vec))
            for new_coord in new_coords
                if periodic
                    l > 2 ? new_coord .%= l : nothing
                else
                    if any(x -> x > l, new_coord)
                        valid = false
                        break
                    end
                end
                push!(vertices, _Vertex(new_coord))
            end

            if valid
                push!(cells, _Cell(vertices))
            end
        end
    end
    return cells
end

function _build_3D_surface_code(l::Int; periodic::Bool=true)
    edges = _compute_cells(l, 1, 3, periodic=periodic)
    faces = _compute_cells(l, 2, 3, periodic=periodic)
    volumes = _compute_cells(l, 3, 3, periodic=periodic)

    # In 3D, we map Qubits to Faces (2-cells)
    q_dict = _build_q_dict(faces)
    edge_dict = _build_q_dict(edges)
    volume_dict = _build_q_dict(volumes)

    # Z Stabilizers: Volumes (3-cells) acting on their boundary Faces
    Z_dict = _contains(volumes, faces, 2, l)
    
    # X Stabilizers: Edges (1-cells) acting on incident Faces 
    # (Inverse of Faces containing Edges)
    X_dict = _inverse_dict(_contains(faces, edges, 1, l))

    X_stabs = spzeros(Bool, length(edge_dict), length(q_dict))
    Z_stabs = spzeros(Bool, length(volume_dict), length(q_dict))

    for (e, f_set) in X_dict
        for f in f_set
            X_stabs[edge_dict[e], q_dict[f]] = true
        end
    end

    for (v, f_set) in Z_dict
        for f in f_set
            Z_stabs[volume_dict[v], q_dict[f]] = true
        end
    end

    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)

    X = zero_matrix(F, size(X_stabs)[1], size(X_stabs)[2])
    Z = zero_matrix(F, size(Z_stabs)[1], size(Z_stabs)[2])

    I_x, J_x, _ = findnz(X_stabs)
    for i in 1:length(I_x)
        X[I_x[i], J_x[i]] = F_one
    end
    I_z, J_z, _ = findnz(Z_stabs)
    for i in 1:length(I_z)
        Z[I_z[i], J_z[i]] = F_one
    end

    S = CSSCode(X, Z)
    S.cache[:d] = l
    S.cache[:l_bound] = l
    S.cache[:u_bound] = l
    return S
end

"""
    Computes a dictionary associating with a cell a set of all the subcells it contains,
    where n is the dimension of the sub_cells, and l the linear size of the system. 
"""
function _contains(cells::Set{_Cell}, sub_cells::Set{_Cell}, n::Int, l::Int)
    cell_dict = Dict{_Cell, Set{_Cell}}()
    for cell in cells
        for sub_cell in combinations(collect(cell.vertices), 2^n)
            l == 2 ? new_cell = _Cell(Set(_identify_cells(sub_cell))) : new_cell = _Cell(Set(sub_cell))
            !(new_cell in sub_cells) && continue
            haskey(cell_dict, cell) ? push!(cell_dict[cell], new_cell) : cell_dict[cell] = Set([new_cell])
        end
    end
    return cell_dict
end

"""
    Inverts a dictionary.
"""
function _inverse_dict(d::Dict)
    d_inv = Dict{_Cell, Set{_Cell}}()
    for (k, v) in d, c in v
        haskey(d_inv, c) ? push!(d_inv[c], k) : d_inv[c] = Set([k])
    end
    return d_inv
end

"""
    Identify the cells on the boundaries in the l = 2 case.
"""
function _identify_cells(c::Vector{_Vertex})
    new_c = deepcopy(c)
    l = length(c[1].coordinates)
    for i in 1:l
        fix_boundary::Bool = true
        for v in c
            v.coordinates[i] < 2 ? fix_boundary = false : nothing
        end
        if fix_boundary
            for v in new_c
                v.coordinates[i] = 0
            end
        end
    end
    return new_c
end

"""
    Create and fills check_matrices. 
"""
function _compute_check_matrices(X_dict, Z_dict, q_dict, edge_dict, volume_dict)
    X_stabs = spzeros(Bool, length(X_dict), length(q_dict) )
    Z_stabs = spzeros(Bool, length(Z_dict), length(q_dict) )

    for (e, v) in X_dict
        for c in v
            X_stabs[edge_dict[e], q_dict[c]] = true
        end
    end

    for (v_, v) in Z_dict
        for c in v
            Z_stabs[volume_dict[v_], q_dict[c]] = true
        end
    end

    return X_stabs, Z_stabs
end

"""
    Create vectors of logicals.
"""
function _compute_logical_vectors(X_dict, Z_dict, q_dict)
    X_logicals, Z_logicals = Vector{Vector{Int}}(), Vector{Vector{Int}}()
    n = 1
    for (_, v) in X_dict
        logical = []
        for c in v
            push!(logical, q_dict[c])
        end
        push!(X_logicals, copy(logical))
    end

    for (_, v) in Z_dict
        logical = []
        for c in v
            push!(logical, q_dict[c])
        end
        push!(Z_logicals, copy(logical))
    end

    return X_logicals, Z_logicals
end

"""
    Build qubit dictionary.
"""
function _build_q_dict(faces::Set{_Cell})
    q_dict = Dict{_Cell, Int}()
    n = 0
    for f in faces
        n += 1
        q_dict[f] = n
    end 
    return q_dict
end

"""
    Compute logical operators.
"""
function _compute_logicals(l::Int, n::Int = 2, d::Int = 4)
    Z_logicals = Dict{Vector{Int}, Set{_Cell}}()
    X_logicals = Dict{Vector{Int}, Set{_Cell}}()
    dirs = collect(1:d)
    coords_to_change = collect(combinations(1:d, n))
    for directions in coords_to_change
        original_face = []
        for c in Iterators.product([0:l - 1 for i in 1:2]...)
            z_coord = [0, 0, 0, 0]
            for i in 1:length(directions)
                z_coord[directions[i]] = c[i]
            end
            other_directions = setdiff(dirs, directions)
            z_new_coords = Vector([copy(z_coord) for _ in 1:2^n - 1])
            for i in 1:2^n - 1, j in 1:n
                (i >> (j - 1)) & 1 == 1 ? z_new_coords[i][directions[j]] += 1 : nothing
            end
            vertices = Set{_Vertex}()
            push!(vertices, _Vertex(z_coord))
            for new_coord in z_new_coords
                l > 2 ? new_coord .%= l : nothing
                push!(vertices, _Vertex(new_coord))
            end
            haskey(Z_logicals, directions) ? push!(Z_logicals[directions], _Cell(vertices)) : Z_logicals[directions] = Set([_Cell(vertices)])
            if isempty(original_face)
                original_face = copy(z_new_coords)
                push!(original_face, z_coord)
            end
            vertices = Set{_Vertex}()
            new_face = deepcopy(original_face)
            for pt in new_face, i in 1:length(other_directions)
                pt[other_directions[i]] = c[i]
            end
            vertices = Set{_Vertex}()
            for new_coord in new_face
                l > 2 ? new_coord .%= l : nothing
                push!(vertices, _Vertex(new_coord)) 
            end
            haskey(X_logicals, directions) ? push!(X_logicals[directions], _Cell(vertices)) : X_logicals[directions] = Set([_Cell(vertices)])
        end
    end
    return X_logicals, Z_logicals
end

"""
    Get the indices of the redundant stabilizers.
"""
function _compute_redundant(redundancy::Dict{_Cell, Set{_Cell}}, re_dict::Dict{_Cell, Int})
    redundant = Vector{Vector{Int}}()
    for (_, stabs) in redundancy
        new_v = Vector{Int}()
        for stab in stabs
            push!(new_v, re_dict[stab])
        end
        push!(redundant, new_v)
    end
    return redundant
end

"""
    Function returning the stabilizers and logicals of periodic 4d surface codes of linear size l >= 2. Constructions taken from J. Math. Phys. 43, 4452-4505 (2002).
"""
function ToricCode4D(l::Int)
    l < 2 && throw(DomainError("Input must be >= 2."))

    # computing the logicals and stabilizers of the code
    vertices = _compute_cells(l, 0, periodic = true)
    edges = _compute_cells(l, 1, periodic = true)
    faces = _compute_cells(l, 2, periodic = true)
    volumes = _compute_cells(l, 3, periodic = true)
    hyper_volumes = _compute_cells(l, 4, periodic = true)

    Z_dict = _contains(volumes, faces, 2, l)
    X_dict = _inverse_dict(_contains(faces, edges, 1, l))

    Z_redundancy = _contains(hyper_volumes, volumes, 3, l)
    X_redundancy = _inverse_dict(_contains(edges, vertices, 0, l))

    q_dict = _build_q_dict(faces)
    volume_dict = _build_q_dict(volumes)
    edge_dict = _build_q_dict(edges)

    X_stabs, Z_stabs = _compute_check_matrices(X_dict, Z_dict, q_dict, edge_dict, volume_dict)
    Z_redundant = _compute_redundant(Z_redundancy, volume_dict)
    X_redundant = _compute_redundant(X_redundancy, edge_dict)

    # X_logicals, Z_logicals = _compute_logicals(l)
    # X_logical, Z_logical = _compute_logical_vectors(X_logicals, Z_logicals, q_dict)

    # defining the code objects
    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)

    X = zero_matrix(F, size(X_stabs)[1], size(X_stabs)[2])
    Z = zero_matrix(F, size(Z_stabs)[1], size(Z_stabs)[2])

    I, J, _ = findnz(X_stabs)
    for i in 1:length(I)
        X[I[i], J[i]] = F_one
    end
    I, J, _ = findnz(Z_stabs)
    for i in 1:length(I)
        Z[I[i], J[i]] = F_one
    end

    S = CSSCode(X, Z)
    set_X_metacheck!(S, X_redundant)
    set_Z_metacheck!(S, Z_redundant)
    return S
end

#################################
        # 3D Fracton Codes
#################################

"""
    XCubeModel(L::Int)

Return the X-Cube model on a cubic lattice of linear size `L` with periodic boundaries.
The X-Cube model is a 3D fracton stabilizer code with parameters `[[3L^3, 6L - 3, L]]`.
"""
function XCubeModel(L::Int)
    L >= 2 || throw(DomainError("Lattice size must be >= 2."))
    
    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)

    # 3 L^3 qubits (one per edge: x-edge, y-edge, z-edge)
    # Edges are indexed by (x, y, z, dir) where dir in 1 (x-edge), 2 (y-edge), 3 (z-edge)
    function edge_idx(x, y, z, dir)
        return (mod(x, L) + mod(y, L) * L + mod(z, L) * L^2) * 3 + dir
    end

    n = 3 * L^3
    num_X = L^3
    num_Z = 3 * L^3

    H_X = zero_matrix(F, num_X, n)
    H_Z = zero_matrix(F, num_Z, n)

    row_X = 1
    row_Z = 1

    for z in 0:L-1, y in 0:L-1, x in 0:L-1
        # X stabilizer: Cube at (x,y,z) (12 edges)
        edges_X = [
            edge_idx(x, y, z, 1), edge_idx(x, y+1, z, 1), edge_idx(x, y, z+1, 1), edge_idx(x, y+1, z+1, 1),
            edge_idx(x, y, z, 2), edge_idx(x+1, y, z, 2), edge_idx(x, y, z+1, 2), edge_idx(x+1, y, z+1, 2),
            edge_idx(x, y, z, 3), edge_idx(x+1, y, z, 3), edge_idx(x, y+1, z, 3), edge_idx(x+1, y+1, z, 3)
        ]
        for e in edges_X
            H_X[row_X, e] = F_one
        end
        row_X += 1

        # Z stabilizers: Crosses at vertex (x,y,z)
        
        # xy-plane cross
        H_Z[row_Z, edge_idx(x, y, z, 1)] = F_one
        H_Z[row_Z, edge_idx(x-1, y, z, 1)] = F_one
        H_Z[row_Z, edge_idx(x, y, z, 2)] = F_one
        H_Z[row_Z, edge_idx(x, y-1, z, 2)] = F_one
        row_Z += 1

        # yz-plane cross
        H_Z[row_Z, edge_idx(x, y, z, 2)] = F_one
        H_Z[row_Z, edge_idx(x, y-1, z, 2)] = F_one
        H_Z[row_Z, edge_idx(x, y, z, 3)] = F_one
        H_Z[row_Z, edge_idx(x, y, z-1, 3)] = F_one
        row_Z += 1

        # zx-plane cross
        H_Z[row_Z, edge_idx(x, y, z, 1)] = F_one
        H_Z[row_Z, edge_idx(x-1, y, z, 1)] = F_one
        H_Z[row_Z, edge_idx(x, y, z, 3)] = F_one
        H_Z[row_Z, edge_idx(x, y, z-1, 3)] = F_one
        row_Z += 1
    end

    # Pass directly into the new lazy cache constructor!
    S = CSSCode(H_X, H_Z)
    
    # Distance is strictly L
    S.cache[:d] = L
    S.cache[:l_bound] = L
    S.cache[:u_bound] = L
    S.cache[:dx] = L
    S.cache[:dz] = L
    
    return S
end

"""
    HaahsCubicCode(L::Int)

Return Haah's Cubic Code (Type 15) on a 3D periodic lattice of size `L`.
This is a CSS fracton model with 2 qubits per vertex. The dimension `k` 
fluctuates based on the number-theoretic properties of `L`.
"""
function HaahsCubicCode(L::Int)
    L >= 2 || throw(DomainError("Lattice size must be >= 2."))
    
    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)

    # 2 qubits per vertex. Indexing: q1 = 1, q2 = 2.
    function q_idx(x, y, z, q)
        return (mod(x, L) + mod(y, L) * L + mod(z, L) * L^2) * 2 + q
    end

    n = 2 * L^3
    num_stabs = L^3

    H_X = zero_matrix(F, num_stabs, n)
    H_Z = zero_matrix(F, num_stabs, n)

    row = 1
    for z in 0:L-1, y in 0:L-1, x in 0:L-1
        # X Stabilizer at (x,y,z)
        # Derived from polynomials h1 = 1 + x + y + z, h2 = 1 + xy + yz + zx
        X_q1 = [(x,y,z), (x+1,y,z), (x,y+1,z), (x,y,z+1)]
        X_q2 = [(x,y,z), (x+1,y+1,z), (x,y+1,z+1), (x+1,y,z+1)]
        
        for (vx, vy, vz) in X_q1
            H_X[row, q_idx(vx, vy, vz, 1)] = F_one
        end
        for (vx, vy, vz) in X_q2
            H_X[row, q_idx(vx, vy, vz, 2)] = F_one
        end

        # Z Stabilizer at (x,y,z)
        # Derived from dual polynomials h1*, h2* (inverted coordinates)
        Z_q1 = [(x,y,z), (x-1,y,z), (x,y-1,z), (x,y,z-1)]
        Z_q2 = [(x,y,z), (x-1,y-1,z), (x,y-1,z-1), (x-1,y,z-1)]
        
        for (vx, vy, vz) in Z_q1
            H_Z[row, q_idx(vx, vy, vz, 1)] = F_one
        end
        for (vx, vy, vz) in Z_q2
            H_Z[row, q_idx(vx, vy, vz, 2)] = F_one
        end
        
        row += 1
    end

    # Automatic rank, dimension, and logical computation via the robust CSS pipeline
    return CSSCode(H_X, H_Z)
end

"""
    ToricColorCode666(L::Int)

Return the Toric 6.6.6 Color Code on a periodic lattice of size `L`.
`L` must be a multiple of 3 for the lattice to be 3-colorable.
Generates a `[[2L^2, 4, d]]` CSS code.
"""
function ToricColorCode666(L::Int)
    (L >= 3 && L % 3 == 0) || throw(DomainError("Lattice size L must be a multiple of 3 for periodic 3-coloring."))
    
    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)

    # Qubits are the two triangular faces (Upper and Lower) per square in the grid
    function q_idx(x, y, triangle_type)
        return (mod(x, L) + mod(y, L) * L) * 2 + triangle_type
    end

    n = 2 * L^2
    num_stabs = L^2
    H = zero_matrix(F, num_stabs, n)

    row = 1
    for y in 0:L-1, x in 0:L-1
        # Each vertex (x,y) is a stabilizer of weight 6.
        # It touches 6 surrounding triangles in the triangulated square grid.
        
        # 1: Upper triangles (x,y to x+1,y+1)
        # 2: Lower triangles (x,y to x+1,y+1)
        
        faces = [
            (x, y, 1),       # Upper right
            (x, y, 2),       # Lower right
            (x-1, y, 1),     # Upper left
            (x, y-1, 2),     # Lower right (from cell below)
            (x-1, y-1, 1),   # Upper left (from cell below left)
            (x-1, y-1, 2)    # Lower left
        ]
        
        for (fx, fy, t) in faces
            H[row, q_idx(fx, fy, t)] = F_one
        end
        
        row += 1
    end

    # Color codes are symmetric CSS codes (H_X = H_Z)
    S = CSSCode(H, H)
    
    # Distance is bounded by L
    S.cache[:u_bound] = L
    S.cache[:u_bound_dx] = L
    S.cache[:u_bound_dz] = L
    
    return S
end

"""
    CleveGottesmanCode()

Return the `[[8, 3, 3]]` Cleve-Gottesman (Eight-qubit) code.
This is a small non-degenerate stabilizer code.
"""
function CleveGottesmanCode()
    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)
    
    # We define the 5 mutually commuting generators for the [[8, 3, 3]] code
    # Format: [X_block | Z_block]
    stabs = zero_matrix(F, 5, 16)
    
    # Generator 1: X X X X X X X X
    for i in 1:8
        stabs[1, i] = F_one
    end
    
    # Generator 2: Z Z Z Z Z Z Z Z
    for i in 1:8
        stabs[2, 8 + i] = F_one
    end
    
    # Generator 3: X X Y Y Z Z I I 
    # (Y has both X and Z components)
    stabs[3, 1] = F_one; stabs[3, 2] = F_one; # X
    stabs[3, 3] = F_one; stabs[3, 4] = F_one; stabs[3, 8+3] = F_one; stabs[3, 8+4] = F_one; # Y
    stabs[3, 8+5] = F_one; stabs[3, 8+6] = F_one; # Z
    
    # Generator 4: X Y Z I X Y Z I
    stabs[4, 1] = F_one; stabs[4, 5] = F_one; # X
    stabs[4, 2] = F_one; stabs[4, 8+2] = F_one; stabs[4, 6] = F_one; stabs[4, 8+6] = F_one; # Y
    stabs[4, 8+3] = F_one; stabs[4, 8+7] = F_one; # Z
    
    # Generator 5: Y Z I X I X Y Z
    stabs[5, 1] = F_one; stabs[5, 8+1] = F_one; stabs[5, 7] = F_one; stabs[5, 8+7] = F_one; # Y
    stabs[5, 8+2] = F_one; stabs[5, 8+8] = F_one; # Z
    stabs[5, 4] = F_one; stabs[5, 6] = F_one; # X

    S = StabilizerCode(stabs)
    
    S.cache[:d] = 3
    S.cache[:l_bound] = 3
    S.cache[:u_bound] = 3
    
    return S
end

"""
    TwistDefectSurfaceCode(d::Int)

Return a distance `d` planar surface code containing a single twist defect in the bulk.
The defect increases the logical dimension `k` by merging an X and Z stabilizer.
"""
function TwistDefectSurfaceCode(d::Int)
    d >= 3 || throw(DomainError("Distance must be >= 3."))
    
    # Generate the standard rotated surface code stabilizers
    base_stabs = _R_Surf_stabs(d)
    F = base_ring(base_stabs)
    n = d^2
    num_stabs = nrows(base_stabs)
    
    # To create a defect, we find an adjacent X and Z stabilizer in the center
    # and combine them into a single Y stabilizer, reducing the total rank by 1.
    stabs = zero_matrix(F, num_stabs - 1, 2 * n)
    
    # The X stabilizers sit on even coordinate sums, Z on odd.
    # In the _R_Surf_stabs matrix, row 1 is a central-ish X stab, 
    # and the first Z stab is at row index: 1 + number of X stabs.
    num_X = div(d^2 - 1, 2)
    center_X_idx = div(num_X, 2)
    center_Z_idx = num_X + div(num_X, 2)
    
    row_out = 1
    for r in 1:num_stabs
        if r == center_X_idx
            # Merge the X stabilizer and the Z stabilizer
            for c in 1:(2*n)
                stabs[row_out, c] = base_stabs[center_X_idx, c] + base_stabs[center_Z_idx, c]
            end
            row_out += 1
        elseif r == center_Z_idx
            # Skip this row entirely, as it was merged into the X row above
            continue
        else
            # Copy standard stabilizer
            for c in 1:(2*n)
                stabs[row_out, c] = base_stabs[r, c]
            end
            row_out += 1
        end
    end
    
    # Returns a StabilizerCode (it is no longer CSS due to the Y defect!)
    return StabilizerCode(stabs)
end

"""
    HeavyHexCode(d::Int)

Return a Heavy-Hex subsystem code on a periodic grid of size `d`.
Data qubits sit on the vertices and edges of a hexagonal lattice.
"""
function HeavyHexCode(d::Int)
    d >= 2 || throw(DomainError("Grid size must be >= 2."))
    
    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)
    
    # Generate Hexagonal Grid Coordinates (Brick-wall lattice)
    vertices = Tuple{Int, Int}[]
    for y in 0:d-1, x in 0:d-1
        push!(vertices, (x, y))
    end
    
    # Heavy Hex introduces data qubits on the vertical and horizontal edges
    edges = Tuple{Int, Int, Symbol}[]
    for (x, y) in vertices
        push!(edges, (x, y, :H)) # Horizontal edge extending right
        if iseven(x + y)
            push!(edges, (x, y, :V)) # Vertical edge extending down
        end
    end
    
    qubits = vcat(vertices, edges)
    n = length(qubits)
    q_idx = Dict(q => i for (i, q) in enumerate(qubits))
    
    gauge_ops = zero_matrix(F, 0, 2 * n)
    
    # Weight-2 Gauge Operators: An X and Z gauge exist on every physical link 
    # connecting a vertex to its adjacent edge-qubit.
    for (x, y) in vertices
        # Horizontal link to right edge
        h_edge = (x, y, :H)
        if haskey(q_idx, h_edge)
            row_x = zero_matrix(F, 1, 2 * n)
            row_z = zero_matrix(F, 1, 2 * n)
            row_x[1, q_idx[(x, y)]] = F_one; row_x[1, q_idx[h_edge]] = F_one
            row_z[1, n + q_idx[(x, y)]] = F_one; row_z[1, n + q_idx[h_edge]] = F_one
            gauge_ops = vcat(gauge_ops, row_x, row_z)
        end
        
        # Vertical link to down edge
        if iseven(x + y)
            v_edge = (x, y, :V)
            if haskey(q_idx, v_edge)
                row_x = zero_matrix(F, 1, 2 * n)
                row_z = zero_matrix(F, 1, 2 * n)
                row_x[1, q_idx[(x, y)]] = F_one; row_x[1, q_idx[v_edge]] = F_one
                row_z[1, n + q_idx[(x, y)]] = F_one; row_z[1, n + q_idx[v_edge]] = F_one
                gauge_ops = vcat(gauge_ops, row_x, row_z)
            end
        end
    end
    
    # Passing the gauge operators dynamically computes the subsystem structure
    return SubsystemCode(gauge_ops)
end

"""
    HeavySquareCode(d::Int)

Return a Heavy-Square subsystem code on a periodic grid of size `d`.
Data qubits sit on the vertices and edges of a square lattice.
"""
function HeavySquareCode(d::Int)
    d >= 2 || throw(DomainError("Grid size must be >= 2."))
    
    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)
    
    vertices = Tuple{Int, Int}[]
    for y in 0:d-1, x in 0:d-1
        push!(vertices, (x, y))
    end
    
    # Heavy Square has edge qubits on ALL vertical and horizontal edges
    edges = Tuple{Int, Int, Symbol}[]
    for (x, y) in vertices
        push!(edges, (x, y, :H)) # Horizontal edge
        push!(edges, (x, y, :V)) # Vertical edge
    end
    
    qubits = vcat(vertices, edges)
    n = length(qubits)
    q_idx = Dict(q => i for (i, q) in enumerate(qubits))
    
    gauge_ops = zero_matrix(F, 0, 2 * n)
    
    # Create gauge pairs connecting vertices to their local edges
    for (x, y) in vertices
        neighbors = [(x, y, :H), (x, y, :V)]
        for edge in neighbors
            if haskey(q_idx, edge)
                row_x = zero_matrix(F, 1, 2 * n)
                row_z = zero_matrix(F, 1, 2 * n)
                row_x[1, q_idx[(x, y)]] = F_one; row_x[1, q_idx[edge]] = F_one
                row_z[1, n + q_idx[(x, y)]] = F_one; row_z[1, n + q_idx[edge]] = F_one
                gauge_ops = vcat(gauge_ops, row_x, row_z)
            end
        end
    end
    
    return SubsystemCode(gauge_ops) 
end

"""
    ColorCode4612(d::Int)

Return the 4.6.12 Archimedean Color Code on a grid of size `d`.
Stabilizers sit on the square (4), hexagonal (6), and dodecagonal (12) faces.
"""
function ColorCode4612(d::Int)
    d >= 2 || throw(DomainError("Grid size must be >= 2."))
    
    F = Oscar.Nemo.Native.GF(2)
    F_one = F(1)
    
    # A 4.6.12 lattice can be constructed by truncating a hexagonal grid.
    # Qubits sit in clusters of 6 around the central dodecagons.
    qubits = Tuple{Int, Int, Int}[]
    for y in 0:d-1, x in 0:d-1
        # Each coordinate (x,y) contains a hexagonal cluster of 6 qubits
        for i in 1:6
            push!(qubits, (x, y, i))
        end
    end
    
    n = length(qubits)
    q_idx = Dict(q => i for (i, q) in enumerate(qubits))
    
    stabs = zero_matrix(F, 0, 2 * n)
    
    for y in 0:d-1, x in 0:d-1
        # 1. Hexagonal Face (Weight 6) - Internal to the cluster
        hex_row = zero_matrix(F, 1, 2 * n)
        for i in 1:6
            hex_row[1, q_idx[(x, y, i)]] = F_one
            hex_row[1, n + q_idx[(x, y, i)]] = F_one
        end
        stabs = vcat(stabs, hex_row)
        
        # 2. Square Faces (Weight 4) - Connecting adjacent clusters
        if x < d-1
            sq_row = zero_matrix(F, 1, 2 * n)
            qs = [(x, y, 2), (x, y, 3), (x+1, y, 5), (x+1, y, 6)]
            for q in qs
                sq_row[1, q_idx[q]] = F_one; sq_row[1, n + q_idx[q]] = F_one
            end
            stabs = vcat(stabs, sq_row)
        end
        
        if y < d-1
            sq_row = zero_matrix(F, 1, 2 * n)
            qs = [(x, y, 4), (x, y, 5), (x, y+1, 1), (x, y+1, 2)]
            for q in qs
                sq_row[1, q_idx[q]] = F_one; sq_row[1, n + q_idx[q]] = F_one
            end
            stabs = vcat(stabs, sq_row)
        end
    end
    
    return StabilizerCode(stabs)
end


