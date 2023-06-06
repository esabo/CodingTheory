# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

# include("tricolorcodes488trellis.jl")
# include("tricolorcodes666trellis.jl")

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
    return SubsystemCode(S, L, Gops)
end
Q9143() = GaugedShorCode()

# BUG: doesn't work for m ≂̸ n
"""
    BaconShorCode(m::Int, n::Int)

Return the Bacon-Shor subsystem code on a `m x n` lattice.
"""
function BaconShorCode(m::Int, n::Int)
    F = GF(2)
    Fone = F(1)
    numqubits = m * n

    # X stabilizers: X[r, :] = X[r + 1, :] = 1
    # X gauge: X[r, c] = X[r + 1, c] = 1
    Xstabs = zero_matrix(F, m - 1, numqubits)
    Xgauges = zero_matrix(F, (m - 1) * n, numqubits)
    currrow = 1
    currrowG = 1
    for r in 1:m - 1
        for c in 1:n
            Xstabs[currrow, (r - 1) * m + c] = Fone
            Xstabs[currrow, r * m + c] = Fone
            Xgauges[currrowG, (r - 1) * m + c] = Fone
            Xgauges[currrowG, r * m + c] = Fone
            currrowG += 1
        end
        currrow += 1
    end

    # Z stabilizers: Z[:, c] = Z[:, c + 1] = 1
    # Z gauge: Z[r, c] = Z[r, c + 1] = 1
    Zstabs = zero_matrix(F, n - 1, numqubits)
    Zgauges = zero_matrix(F, (n - 1) * m, numqubits)
    currrow = 1
    currrowG = 1
    for c in 1:n - 1
        for r in 1:m
            Zstabs[currrow, (r - 1) * m + c] = Fone
            Zstabs[currrow, (r - 1) * m + c + 1] = Fone
            Zgauges[currrowG, (r - 1) * m + c] = Fone
            Zgauges[currrowG, (r - 1) * m + c + 1] = Fone
            currrowG += 1
        end
        currrow += 1
    end

    # X logical: X[1, :] = 1
    Xlogical = zero_matrix(F, 1, numqubits)
    # TODO: consider @simd or @unroll here
    for c in 1:n
        Xlogical[1, c] = Fone
    end

    # Z logical: Z[:, 1] = 1
    Zlogical = zero_matrix(F, 1, numqubits)
    # TODO: consider @simd or @unroll here
    for r in 1:m
        Zlogical[1, (r - 1) * m + 1] = Fone
    end
    
    stabs = Xstabs ⊕ Zstabs
    logs = Xlogical ⊕ Zlogical
    gauges = Xgauges ⊕ Zgauges
    # S = SubsystemCodeCSS(Xstabs, Zstabs, (Xlogical, Zlogical), {})
    # S = SubsystemCode(stabs, logs, gauges, true)
    S = SubsystemCode(gauges)
    setstabilizers!(S, stabs)
    S.Xstabs = Xstabs
    S.Zstabs = Zstabs
    # CSS Xsigns and Zsigns don't need to be updated, should be same length and still chi(0)
    setlogicals!(S, logs)
    m == n && setminimumdistance!(S, m)
    # Z distance is m
    # X distance is n
    return S
end

"""
    BaconShorCode(d::Int)

Return the Bacon-Shor subsystem code on a `d x d` lattice.
"""
BaconShorCode(d::Int) = BaconShorCode(d, d)

"""
    BravyiSubsystemCode(A::fq_nmod_mat)
    GeneralizedBaconShorCode(A::fq_nmod_mat)

Return the generalied Bacon-Shor code defined by Bravyi in "Subsystem Codes With Spatially Local
Generators", (2011).
"""
function BravyiSubsystemCode(A::CTMatrixTypes)
    iszero(A) && throw(ArgumentError("The input matrix cannot be zero."))
    F = base_ring(A)
    Int(order(F)) == 2 || throw(ArgumentError("Construction is only valid for binary martices."))

    n = 0
    nr, nc = size(A)
    rowwts = zeros(Int, 1, nr)
    colwts = zeros(Int, 1, nc)
    linearindex = Dict{Tuple{Int, Int}, Int}()
    for r in 1:nr
        for c in 1:nc
            if !iszero(A[r, c])
                n += 1
                rowwts[r] += 1
                colwts[c] += 1
                linearindex[(r, c)] = n
            end
        end
    end

    for i in 1:nr
        rowwts[i] == 1 && throw(ArgumentError("The input matrix cannot have a row of weight one."))
    end
    for i in 1:nc
        colwts[i] == 1 && throw(ArgumentError("The input matrix cannot have a column of weight one."))
    end

    # really only need to generate every consequetive pair
    Xgauges = zero_matrix(F, sum([colwts[i] - 1 for i in 1:length(colwts)]), n)
    currrow = 1
    Fone = F(1)
    for c in 1:nc
        r1 = 1
        while r1 < nr
            if !iszero(A[r1, c])
                atend = true
                for r2 in r1 + 1:nr
                    if !iszero(A[r2, c])
                        Xgauges[currrow, linearindex[r1, c]] = Fone
                        Xgauges[currrow, linearindex[r2, c]] = Fone
                        currrow += 1
                        r1 = r2
                        atend = false
                        break
                    end
                end
                atend && (r1 = nr;)
            else
                r1 += 1
            end
        end
    end

    # really only need to generate every consequetive pair
    Zgauges = zero_matrix(F, sum([rowwts[i] - 1 for i in 1:length(rowwts)]), n)
    currrow = 1
    for r in 1:nr
        c1 = 1
        while c1 < nc
            if !iszero(A[r, c1])
                atend = true
                for c2 in c1 + 1:nc
                    if !iszero(A[r, c2])
                        Zgauges[currrow, linearindex[r, c1]] = Fone
                        Zgauges[currrow, linearindex[r, c2]] = Fone
                        currrow += 1
                        c1 = c2
                        atend = false
                        break
                    end
                end
                atend && (c1 = nc;)
            else
                c1 += 1
            end
        end
    end

    # # every pair of 1's in row of A gets an Z gauge operator
    # Zgauges = zero_matrix(F, totZgauges, n)
    # currrow = 1
    # Fone = F(1)
    # for r in 1:nr
    #     for c1 in 1:nc - 1
    #         if !iszero(A[r, c1])
    #             for c2 in c1 + 1:nc
    #                 if !iszero(A[r, c2])
    #                     Zgauges[currrow, linearindex[r, c1]] = Fone
    #                     Zgauges[currrow, linearindex[r, c2]] = Fone
    #                     currrow += 1
    #                 end
    #             end
    #         end
    #     end
    # end

    # # every pair of 1's in col of A gets a X gauge operator
    # Xgauges = zero_matrix(F, totXgauges, n)
    # currrow = 1
    # for c in 1:nc
    #     for r1 in 1:nr - 1
    #         if !iszero(A[r1, c])
    #             for r2 in r1 + 1:nr
    #                 if !iszero(A[r2, c])
    #                     Xgauges[currrow, linearindex[r1, c]] = Fone
    #                     Xgauges[currrow, linearindex[r2, c]] = Fone
    #                     currrow += 1
    #                 end
    #             end
    #         end
    #     end
    # end

    S = SubsystemCode(Xgauges ⊕ Zgauges)
    minrowwt = minimum(rowwts)
    mincolwt = minimum(colwts)
    setminimumdistance!(S, minimum([minrowwt, mincolwt]))
    # TODO: also set dx and dz
    return S
end
GeneralizedBaconShorCode(A::CTMatrixTypes) = BravyiSubsystemCode(A)

# subsystem codes were described by Bacon and Casaccino in [19]. The construction of [19] starts from a pair of classical linear codes C1 = [n1, k1, d1] and C2 = [n2, k2, d2]. A quantum subsystem code is then defined by placing a physical qubit at every cell of a ma- trix A of size n1 × n2. The X-part of the gauge group is defined by replicating the parity checks of C1 in every col- umn of A (in the X-basis). Similarly, the Z-part of the gauge group is defined by replicating the parity checks of C2 in every row of A (in the Z-basis). The resulting subsystem code has parameters [n1n2, k1k2, min (d1, d2)].

# Napp & Preskill, "Optimal Bacon-Shor Codes", (2012)
"""
    NappPreskill3DCode(m::Int, n::Int, k::Int)

Return the Napp and Preskill 3D, modifed Bacon-Shor code.
"""
function NappPreskill3DCode(m::Int, n::Int, k::Int)
    (2 <= m && 2 <= n && 2 <= k) || throw(DomainError("Lattice dimensions must be at least two"))

    linearindex = Dict{Tuple{Int, Int, Int}, Int}()
    for (i, tup) in enumerate(Base.Iterators.product(1:m, 1:n, 1:k))
        linearindex[tup] = i
    end

    len = m * n * k
    F = GF(2)
    Fone = F(1)
    gauges = zero_matrix(F, (m - 1) * n * k + m * (n - 1) * k + m * n * (k - 1), 2 * len)
    currrow = 1
    # X's
    for l in 1:k
        for i in 1:m
            for j in 1:n
                if j != n
                    gauges[currrow, linearindex[i, j, l]] = gauges[currrow, linearindex[i, j + 1, l]] = Fone
                    currrow += 1
                end
                if i != m
                    gauges[currrow, linearindex[i, j, l]] = gauges[currrow, linearindex[i + 1, j, l]] = Fone
                    currrow += 1
                end
            end
        end
    end

    # Z's
    for i in 1:m
        for j in 1:n
            for l in 1:k - 1
                gauges[currrow, linearindex[i, j, l] + len] = gauges[currrow, linearindex[i, j, l + 1] + len] = Fone
                currrow += 1
            end
        end
    end

    S = SubsystemCode(gauges, missing, :VS)
    S.dx = k
    S.dz = m * n
    S.d = minimum([k, m * n])
    return S
end

# Napp & Preskill, "Optimal Bacon-Shor Codes", (2012)
"""
    NappPreskill4DCode(x::Int, y::Int, z::Int, w::Int)

Return the Napp and Preskill 4D, modifed Bacon-Shor code.
"""
function NappPreskill4DCode(x::Int, y::Int, z::Int, w::Int)
    (2 <= x && 2 <= y && 2 <= z && 2 <= w) || throw(DomainError("Lattice dimensions must be at least two"))

    linearindex = Dict{Tuple{Int, Int, Int, Int}, Int}()
    for (i, tup) in enumerate(Base.Iterators.product(1:x, 1:y, 1:z, 1:w))
        linearindex[tup] = i
    end

    len = x * y * z * w
    F = GF(2)
    Fone = F(1)
    gauges = zero_matrix(F, (x - 1) * y * z * w + x * (y - 1) * z * w + x * y * (z - 1) * w + x * y * z * (w - 1), 2 * len)
    currrow = 1
    ## XX acting on each neighboring qubits in each xy-plane with z, w fixed
    ## ZZ in zw-plane with x, y fixed

    # X's
    for k in 1:z
        for l in 1:w
            for i in 1:x
                for j in 1:y
                    if i != x
                        gauges[currrow, linearindex[i, j, k, l]] = gauges[currrow, linearindex[i + 1, j, k, l]] = Fone
                        currrow += 1
                    end
                    if j != y
                        gauges[currrow, linearindex[i, j, k, l]] = gauges[currrow, linearindex[i, j + 1, k, l]] = Fone
                        currrow += 1
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
                        gauges[currrow, linearindex[i, j, k, l] + len] = gauges[currrow, linearindex[i, j, k + 1, l] + len] = Fone
                        currrow += 1
                    end
                    if l != w
                        gauges[currrow, linearindex[i, j, k, l] + len] = gauges[currrow, linearindex[i, j, k, l + 1] + len] = Fone
                        currrow += 1
                    end
                end
            end
        end
    end

    return SubsystemCode(gauges, missing, :VS)
end

# Bravyi et al, "Subsystem surface codes with three-qubit check operators", (2013)
"""
    SubsystemToricCode(m::Int, n::Int)

Return the subsystem toric code on a rectangular lattice with `m` rows and `n`
columns of squares.
"""
function SubsystemToricCode(m::Int, n::Int)
    (2 <= m && 2 <= n) || throw(DomainError("Lattice dimensions must be at least two"))

    F = GF(2)
    Fone = F(1)
    len = 3 * m * n
    gauges = zero_matrix(F, 4 * m * n, 2 * len)
    stabs = zero_matrix(F, 2 * m * n, 2 * len)
    currrowstab = 1
    currrowgauge = 1
    for r in 1:m
        topleft = 3 * n * (r - 1) + 1
        rowright = topleft + 2 * n - 1
        for c in 1:n
            # top left - Z
            gauges[currrowgauge, topleft + len] = Fone
            gauges[currrowgauge, topleft + 1 + len] = Fone
            gauges[currrowgauge, rowright + c + len] = Fone
            currrowgauge += 1

            # top right - X
            gauges[currrowgauge, topleft + 1] = Fone
            if c != n
                gauges[currrowgauge, topleft + 2] = Fone
                gauges[currrowgauge, rowright + c + 1] = Fone
            else
                gauges[currrowgauge, rowright - 2 * n + 1] = Fone
                gauges[currrowgauge, rowright + 1] = Fone
            end
            currrowgauge += 1

            # bottom left - X
            gauges[currrowgauge, rowright + c] = Fone
            if r != m
                gauges[currrowgauge, rowright + n + 2 * (c - 1) + 1] = Fone
                gauges[currrowgauge, rowright + n + 2 * (c - 1) + 2] = Fone
            else
                gauges[currrowgauge, 2 * (c - 1) + 1] = Fone
                gauges[currrowgauge, 2 * (c - 1) + 2] = Fone
            end
            currrowgauge += 1

            # bottom right - Z
            if r == m && c == n
                gauges[currrowgauge, 1 + len] = Fone
                gauges[currrowgauge, 2 * n + len] = Fone
                gauges[currrowgauge, rowright + 1 + len] = Fone
            elseif r == m
                gauges[currrowgauge, rowright + c + 1 + len] = Fone
                gauges[currrowgauge, 2 * (c - 1) + 2 + len] = Fone
                gauges[currrowgauge, 2 * (c - 1) + 3 + len] = Fone
            elseif c != n
                gauges[currrowgauge, rowright + c + 1 + len] = Fone
                gauges[currrowgauge, rowright + n + 2 * (c - 1) + 2 + len] = Fone
                gauges[currrowgauge, rowright + n + 2 * (c - 1) + 3 + len] = Fone
            else
                gauges[currrowgauge, rowright + 1 + len] = Fone
                gauges[currrowgauge, rowright + n + 1 + len] = Fone
                gauges[currrowgauge, rowright + n + 2 * (c - 1) + 2 + len] = Fone
            end
            currrowgauge += 1

            # X
            stabs[currrowstab, topleft + 1] = Fone
            if c != n
                stabs[currrowstab, topleft + 2] = Fone
                stabs[currrowstab, rowright + c + 1] = Fone
            else
                stabs[currrowstab, rowright - 2 * n + 1] = Fone
                stabs[currrowstab, rowright + 1] = Fone
            end
            stabs[currrowstab, rowright + c] = Fone
            if r != m
                stabs[currrowstab, rowright + n + 2 * (c - 1) + 1] = Fone
                stabs[currrowstab, rowright + n + 2 * (c - 1) + 2] = Fone
            else
                stabs[currrowstab, 2 * (c - 1) + 1] = Fone
                stabs[currrowstab, 2 * (c - 1) + 2] = Fone
            end
            currrowstab += 1

            # Z
            stabs[currrowstab, topleft + len] = Fone
            stabs[currrowstab, topleft + 1 + len] = Fone
            stabs[currrowstab, rowright + c + len] = Fone
            if r == m && c == n
                stabs[currrowstab, 1 + len] = Fone
                stabs[currrowstab, 2 * n + len] = Fone
                stabs[currrowstab, rowright + 1 + len] = Fone
            elseif r == m
                stabs[currrowstab, rowright + c + 1 + len] = Fone
                stabs[currrowstab, 2 * (c - 1) + 2 + len] = Fone
                stabs[currrowstab, 2 * (c - 1) + 3 + len] = Fone
            elseif c != n
                stabs[currrowstab, rowright + c + 1 + len] = Fone
                stabs[currrowstab, rowright + n + 2 * (c - 1) + 2 + len] = Fone
                stabs[currrowstab, rowright + n + 2 * (c - 1) + 3 + len] = Fone
            else
                stabs[currrowstab, rowright + 1 + len] = Fone
                stabs[currrowstab, rowright + n + 1 + len] = Fone
                stabs[currrowstab, rowright + n + 2 * (c - 1) + 2 + len] = Fone
            end
            currrowstab += 1
            topleft += 2
        end
    end
    
    logs = zero_matrix(F, 4, 2 * len)
    # top row is an X for pair one and a Z for pair two
    for c in 1:2 * n
        logs[1, c] = Fone
        logs[4, c + len] = Fone
    end
    # left column is a Z for pair one and an X for pair two
    for r in 1:m
        topleft = 3 * n * (r - 1) + 1
        logs[2, topleft + len] = Fone
        logs[2, topleft + 2 * n + len] = Fone
        logs[3, topleft] = Fone
        logs[3, topleft + 2 * n] = Fone
    end
    
    S = SubsystemCode(gauges, missing, :VS)
    S.k == 2 || error("Got wrong dimension for periodic case.")
    setstabilizers!(S, stabs)
    setlogicals!(S, logs)
    setminimumdistance!(S, minimum([m, n]))
    S.dx = S.d
    S.dx = S.d
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

    F = GF(2)
    Fone = F(1)
    len = (3n + 2) * m + 2 * n + 1
    gauges = zero_matrix(F, 4 * m * n + 2 * n + 2 * m, 2 * len)
    stabs = zero_matrix(F, 2 * m * n + 2 * n + 2 * m, 2 * len)
    currrowstab = 1
    currrowgauge = 1
    for r in 1:m
        topleft = (3 * n + 2) * (r - 1) + 1
        rowright = topleft + 2 * n
        for c in 1:n
            # top left - Z
            gauges[currrowgauge, topleft + len] = Fone
            gauges[currrowgauge, topleft + 1 + len] = Fone
            gauges[currrowgauge, rowright + c + len] = Fone
            currrowgauge += 1

            # top right - X
            gauges[currrowgauge, topleft + 1] = Fone
            gauges[currrowgauge, topleft + 2] = Fone
            gauges[currrowgauge, rowright + c + 1] = Fone
            currrowgauge += 1

            # bottom left - X
            gauges[currrowgauge, rowright + c] = Fone
            gauges[currrowgauge, rowright + n + 2 * (c - 1) + 2] = Fone
            gauges[currrowgauge, rowright + n + 2 * (c - 1) + 3] = Fone
            currrowgauge += 1

            # bottom right - Z
            gauges[currrowgauge, rowright + c + 1 + len] = Fone
            gauges[currrowgauge, rowright + n + 2 * (c - 1) + 3 + len] = Fone
            gauges[currrowgauge, rowright + n + 2 * (c - 1) + 4 + len] = Fone
            currrowgauge += 1

            # X
            stabs[currrowstab, topleft + 1] = Fone
            stabs[currrowstab, topleft + 2] = Fone
            stabs[currrowstab, rowright + c + 1] = Fone
            stabs[currrowstab, rowright + c] = Fone
            stabs[currrowstab, rowright + n + 2 * (c - 1) + 2] = Fone
            stabs[currrowstab, rowright + n + 2 * (c - 1) + 3] = Fone
            currrowstab += 1

            # Z
            stabs[currrowstab, topleft + len] = Fone
            stabs[currrowstab, topleft + 1 + len] = Fone
            stabs[currrowstab, rowright + c + len] = Fone
            stabs[currrowstab, rowright + c + 1 + len] = Fone
            stabs[currrowstab, rowright + n + 2 * (c - 1) + 3 + len] = Fone
            stabs[currrowstab, rowright + n + 2 * (c - 1) + 4 + len] = Fone
            currrowstab += 1
            topleft += 2
        end
    end

    # left - X
    for r in 1:m
        topleft = (3 * n + 2) * (r - 1) + 1
        shift = 2 * n + 1
        stabs[currrowstab, topleft] = Fone
        stabs[currrowstab, topleft + shift] = Fone
        currrowstab += 1
        gauges[currrowgauge, topleft] = Fone
        gauges[currrowgauge, topleft + shift] = Fone
        currrowgauge += 1
    end

    # right - X
    for r in 1:m
        topleft = (3 * n + 2) * (r - 1) + 1
        rowright = topleft + 2 * n
        stabs[currrowstab, rowright + n + 1] = Fone
        stabs[currrowstab, rowright + 3 * n + 2] = Fone
        currrowstab += 1
        gauges[currrowgauge, rowright + n + 1] = Fone
        gauges[currrowgauge, rowright + 3 * n + 2] = Fone
        currrowgauge += 1
    end

    # top - Z
    topleft = 1
    for c in 1:n
        stabs[currrowstab, topleft + 1 + len] = Fone
        stabs[currrowstab, topleft + 2 + len] = Fone
        currrowstab += 1
        gauges[currrowgauge, topleft + 1 + len] = Fone
        gauges[currrowgauge, topleft + 2 + len] = Fone
        currrowgauge += 1
        topleft += 2
    end

    # bottom - Z
    bottomleft = len - 2 * n
    for c in 1:n
        stabs[currrowstab, bottomleft + len] = Fone
        stabs[currrowstab, bottomleft + 1 + len] = Fone
        currrowstab += 1
        gauges[currrowgauge, bottomleft + len] = Fone
        gauges[currrowgauge, bottomleft + 1 + len] = Fone
        currrowgauge += 1
        bottomleft += 2
    end

    logs = zero_matrix(F, 2, 2 * len)
    # top row is a logical X
    for c in 1:2 * n + 1
        logs[1, c] = Fone
    end
    # left column is a logical Z 
    for r in 1:m
        logs[2, (3 * n + 2) * (r - 1) + 1 + len] = Fone
        logs[2, (3 * n + 2) * (r - 1) + 2 * n + 2 + len] = Fone
    end
    logs[2, 2 * len -  2 * n] = Fone
    
    S = SubsystemCode(gauges, missing, :VS)
    S.k == 1 || error("Got wrong dimension for non-periodic case.")
    setstabilizers!(S, stabs)
    setlogicals!(S, logs)
    setminimumdistance!(S, minimum([m, n]))
    S.dx = 2 * n + 1
    S.dz = 2 * m + 1
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
    FiveQubitCode()
    Q513()

Return the `[[5, 1, 3]]` perfect qubit stabilizer code.
"""
# is a perfect code
FiveQubitCode() = StabilizerCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"])
Q513() = FiveQubitCode()

# should also do a test for other CSS construction via Hamming code and actually make that one default
"""
    Steanecode()
    Q713()

Return the `[[7, 1, 3]]` Steane code.
"""
SteaneCode() = CSSCode(["XXXXIII", "XXIIXXI", "XIXIXIX", "ZZZZIII", "ZZIIZZI", "ZIZIZIZ"])
Q713() = SteaneCode()
# _SteaneCodeTrellis() = CSSCode(["XXIXXII", "IXXIXXI", "IIIXXXX", "ZZIZZII", "IZZIZZI", "IIIZZZZ"])
# also ZZIZZII, ZIZZIZI, IZZZIIZ, XXIXXII, XIXXIXI, IXXXIIX

"""
    Shorcode()
    Q913()

Return the `[[9, 1, 3]]` Shor code.
"""
ShorCode() = CSSCode(["ZZIIIIIII", "IZZIIIIII", "IIIZZIIII", "IIIIZZIII", "IIIIIIZZI", "IIIIIIIZZ",
    "XXXXXXIII", "IIIXXXXXX"])
Q913() = ShorCode()

Q412() = CSSCode(["XXXX", "ZZII", "IIZZ"])
Q422() = CSSCode(["XXXX", "ZZZZ"])
Q511() = StabilizerCode(["ZXIII", "XZXII", "IXZXI", "IIXZX"])

function Q823()
    F = GF(2)
    stabs = matrix(F, [1 0 0 0 1 0 0 0 1 1 1 1 0 0 0 0;
    0 0 0 1 0 1 0 0 1 0 0 0 0 1 0 0;
    0 1 0 0 1 1 1 0 0 0 1 1 1 0 1 0;
    0 0 1 0 1 1 1 0 0 1 1 0 1 1 0 0;
    0 0 1 1 1 0 1 0 0 0 0 1 0 1 1 1;
    0 0 0 0 0 0 1 1 0 0 1 0 0 0 1 0]);
    return StabilizerCode(stabs)
end

"""
    Q15RM()
    Q1513()

Return the `[[15, 1, 3]]` quantum Reed-Muller code.
"""
Q15RM() = StabilizerCode(["ZIZIZIZIZIZIZIZ", "IZZIIZZIIZZIIZZ", "IIIZZZZIIIIZZZZ", "IIIIIIIZZZZZZZZ",
    "IIZIIIZIIIZIIIZ", "IIIIZIZIIIIIZIZ", "IIIIIZZIIIIIIZZ", "IIIIIIIIIZZIIZZ", "IIIIIIIIIIIZZZZ",
    "IIIIIIIIZIZIZIZ", "XIXIXIXIXIXIXIX", "IXXIIXXIIXXIIXX", "IIIXXXXIIIIXXXX", "IIIIIIIXXXXXXXX"])
Q1513() = Q15RM()

"""
    Q1573()

Return the `[[15, 7, 3]]` quantum Hamming code.
"""
Q1573() = StabilizerCode(["IIIIIIIXXXXXXXX", "IIIXXXXIIIIXXXX", "IXXIIXXIIXXIIXX", "XIXIXIXIXIXIXIX",
    "IIIIIIIZZZZZZZZ", "IIIZZZZIIIIZZZZ", "IZZIIZZIIZZIIZZ", "ZIZIZIZIZIZIZIZ"])
    # one can use a basis for this such that the first logical pair is transversal X, Z

#############################
 # Triangular Surface Codes
#############################

function _triangularlattice(L::Int)
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

function _triangularlatticeXstabilizers(L::Int, numbering::Array{Int64, 3})
    F = GF(2)
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

function _triangularlatticeZstabilizers(L::Int, numbering::Array{Int64, 3})
    F = GF(2)
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

function _triangularlatticeXlogicals(L::Int, numbering::Array{Int64, 3})
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

function _triangularlatticeZlogicals(L::Int, numbering::Array{Int64, 3}, symp::Bool=true)
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
    numbering = _triangularlattice(L)
    Xstabs = _triangularlatticeXstabilizers(L, numbering)
    # println(rank(Xstabs))
    Zstabs = _triangularlatticeZstabilizers(L, numbering)
    # println(Zstabs)
    # logicals = [triangularlatticeXlogicals(L, numbering), triangularlatticeZlogicals(L, numbering)]
    return CSSCode(Xstabs[1:end - 1, :], Zstabs[1:end - 1, :])
end

#############################
   # Rotated Surface Codes
#############################

function _RSurfstabs(d::Int)
    n = d^2
    F = GF(2)
    Fone = F(1)
    S = zero_matrix(F, n - 1, 2 * n)
    row = 1

    # X's
    i = 1
    while i <= n - d
        S[row, i] = Fone
        S[row, i + 1] = Fone
        S[row, i + d] = Fone
        S[row, i + d + 1] = Fone
        row += 1
        if (i + 2) % d == 0
            i += 4
        else
            i += 2
        end
    end

    # top row X's
    i = 2
    while i <= d - 1
        S[row, i] = Fone
        S[row, i + 1] = Fone
        row += 1
        i += 2
    end

    # bottom row X's
    i = d * (d - 1) + 1
    while i <= d * d - 2
        S[row, i] = Fone
        S[row, i + 1] = Fone
        row += 1
        i += 2
    end

    # Z's
    i = 2
    while i < n - d
        S[row, i + n] = Fone
        S[row, i + 1 + n] = Fone
        S[row, i + d + n] = Fone
        S[row, i + d + 1 + n] = Fone
        row += 1
        if (i + 2) % d == 0
            i += 4
        else
            i += 2
        end
    end

    # left Z's
    i = 1
    while i < d * (d - 1)
        S[row, i + n] = Fone
        S[row, i + d + n] = Fone
        row += 1
        i += 2 * d
    end

    # right Z's
    i = 2 * d
    while i < d * d
        S[row, i + n] = Fone
        S[row, i + d + n] = Fone
        row += 1
        i += 2 * d
    end

    return S
end

function _RSurflogs(F::CTFieldTypes, d::Int)
    n = d^2
    Fone = F(1)
    logs = zero_matrix(F, 2, 2 * n)
    i = d
    while i <= d * d
        logs[1, i] = Fone
        i += d
    end
    i = 1
    while i <= d
        logs[2, i + n] = Fone
        i += 1
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
# BUG: doesn't work for even distances
function RotatedSurfaceCode(d::Int)
    d >= 3 || throw(DomainError("Current implementation requires d ≥ 3."))

    stabs = _RSurfstabs(d)
    S = StabilizerCode(stabs)
    d <= 10 && setlogicals!(S, _RSurflogs(base_ring(stabs), d))
    return S
end

#############################
     # XZZX Surface Codes
#############################

function _XZZXstabslogs(d::Int)
    n = d^2
    F = GF(2)
    S = zero_matrix(F, n - 1, 2 * n)
    row = 1
    Fone = F(1)

    i = 1
    for i in 1:n - d
        if i % d != 0
            S[row, i] = Fone
            S[row, i + 1 + n] = Fone
            S[row, i + d + n] = Fone
            S[row, i + d + 1] = Fone
            row += 1;
        end
    end

    # top row ZX's
    i = 2
    while i <= d - 1
        S[row, i + n] = Fone
        S[row, i + 1] = Fone
        row += 1
        i += 2
    end

    # bottom row XZ's
    i = d * (d - 1) + 1
    while i <= d * d - 2
        S[row, i] = Fone
        S[row, i + 1 + n] = Fone
        row += 1
        i += 2
    end

    # left ZX's
    i = 1
    while i < d * (d - 1)
        S[row, i + n] = Fone
        S[row, i + d] = Fone
        row += 1
        i += 2 * d
    end

    # right XZ's
    i = 2 * d
    while i < d * d
        S[row, i] = Fone
        S[row, i + d + n] = Fone
        row += 1
        i += 2 * d
    end

    logs = zero_matrix(F, 2, 2 * n)
    i = d
    count = 1
    while i <= d * d
        if count % 2 == 1
            logs[1, i] = Fone
        else
            logs[1, i + n] = Fone
        end
        i += d
        count += 1
    end
    i = 1
    count = 1
    while i <= d
        if count % 2 == 1
            logs[2, i + n] = Fone
        else
            logs[2, i] = Fone
        end
        i += 1
        count += 1
    end

    return S, logs
end

"""
    XZZXSurfaceCode(d::Int)

Return the `[[d^2, 1, d]]` XZZX surface code.
"""
function XZZXSurfaceCode(d::Int)
    d >= 3 || throw(DomainError("Current implementation requires d ≥ 3."))

    stabs, logs = _XZZXstabslogs(d)
    S = StabilizerCode(stabs)
    setlogicals!(S, logs)
    return S
end

################################
 # Triangular Color Codes 4.8.8
################################

"""
    TriangularColorCode488(d::Int)

Return the 4.8.8 triangular color code of distance `d` with trellis numbering.
"""
function TriangularColorCode488(d::Int)
    3 <= d <= 19 || throw(DomainError("Current implementation requires 3 ≤ d ≤ 19."))

    if d == 3
        # S, logs = _488d3trellis()
        @load "../data/488d3stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    elseif d == 5
        # S, logs = _488d5trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/488d5stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    elseif d == 7
        # S, logs = _488d7trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/488d7stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    elseif d == 9
        # S, logs = _488d9trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/488d9stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    elseif d == 11
        # S, logs = _488d11trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/488d11stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    elseif d == 13
        # S, logs = _488d13trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/488d13stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    elseif d == 15
        # S, logs = _488d15trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/488d15stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    elseif d == 17
        # S, logs = _488d17trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/488d17stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    elseif d == 19
        # S, logs = _488d19trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/488d19stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    # elseif d == 21
    #     # S, logs = _488d21trellis()
    #     # S = StabilizerCode(stabs)
    #     # setlogicals!(Q, logs)
    #     @load "../data/488d21stabslogs_trellis.jld2" S l
    #     F = GF(2)
    #     stabs = matrix(F, S)
    #     S = StabilizerCode(stabs)
    #     l = matrix(F, l)
    #     # setlogicals!(S, [(l[1, :], l[2, :])])
    #     setlogicals!(S, l)
    #     return S
    end
end

################################
 # Triangular Color Codes 6.6.6
################################

"""
    TriangularColorCode666(d::Int)

Return the 6.6.6 triangular color code of distance `d` with trellis numbering.
"""
function TriangularColorCode666(d::Int)
    3 <= d <= 21 || throw(DomainError("Current implementation requires 3 ≤ d ≤ 21."))

    if d == 3
        # same as 4.8.8
        # S, logs = _488d3trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/488d3stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    elseif d == 5
        # S, logs = _666d5trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/666d5stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    elseif d == 7
        # S, logs = _666d7trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/666d7stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    elseif d == 9
        # S, logs = _666d9trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/666d9stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    elseif d == 11
        # S, logs = _666d11trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/666d11stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    elseif d == 13
        # S, logs = _666d13trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/666d13stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    elseif d == 15
        # S, logs = _666d15trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/666d15stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    elseif d == 17
        # S, logs = _666d17trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/666d17stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    elseif d == 19
        # S, logs = _666d19trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/666d19stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    elseif d == 21
        # S, logs = _666d21trellis()
        # S = StabilizerCode(stabs)
        # setlogicals!(Q, logs)
        @load "../data/666d21stabslogs_trellis.jld2" S l
        F = GF(2)
        stabs = matrix(F, S)
        S = StabilizerCode(stabs)
        l = matrix(F, l)
        # setlogicals!(S, [(l[1, :], l[2, :])])
        setlogicals!(S, l)
        return S
    end
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

    F = GF(2)
    Fone = F(1)
    A = zero_matrix(F, d^2, 2 * d^2) # stars, X stabilizers
    B = zero_matrix(F, d^2, 2 * d^2) # faces, Z stabilizers
    qubit = 1
    rowA = 1
    rowB = 1
    for r in 1:2 * d
        if isodd(r)
            for c in 1:d
                # println("r = $r, c = $c, rowA = $rowA")
                if r != 2 * d - 1 && c != d
                    A[rowA, qubit] = A[rowA, qubit + d] = A[rowA, qubit + d + 1] = A[rowA, qubit + 2 * d] = Fone
                elseif r == 2 * d - 1 && c != d
                    A[rowA, qubit] = A[rowA, qubit + d] = A[rowA, qubit + d + 1] = A[rowA, c] = Fone
                elseif r != 2 * d - 1 && c == d
                    A[rowA, qubit] = A[rowA, qubit + d] = A[rowA, qubit + 1] = A[rowA, qubit + 2 * d] = Fone
                elseif r == 2 * d - 1 && c == d
                    A[rowA, qubit] = A[rowA, qubit + d] = A[rowA, qubit + 1] = A[rowA, c] = Fone
                else
                    error("Ran into unaccounted for case in creating the toric code lattice.")
                end
                rowA += 1
                qubit += 1
            end
        else
            for c in 1:d
                # println("r = $r, c = $c, rowB = $rowB")
                if r != 2 * d && c == 1
                    B[rowB, qubit] = B[rowB, qubit + d] = B[rowB, qubit + 2 * d] = B[rowB, qubit + 2 * d - 1] = Fone
                elseif r != 2 * d && c != 1
                    B[rowB, qubit] = B[rowB, qubit + d - 1] = B[rowB, qubit + d] = B[rowB, qubit + 2 * d] = Fone
                elseif r == 2 * d && c == 1
                    B[rowB, qubit] = B[rowB, d] = B[rowB, d + 1] = B[rowB, 1] = Fone
                elseif r == 2 * d && c != 1
                    B[rowB, qubit] = B[rowB, c - 1] = B[rowB, c] = B[rowB, c + d] = Fone
                else
                    println("here")
                    error("Ran into unaccounted for case in creating the toric code lattice.")
                end
                rowB += 1
                qubit += 1
            end
        end
    end
    S = CSSCode(A, B)
    X1 = zero_matrix(S.F, 1, 4 * d^2)
    for c in 1:d
        X1[1, c + d] = Fone
    end
    Z1 = zero_matrix(S.F, 1, 4 * d^2)
    for r in 1:2:2 * d
        Z1[1, r * d + 1 + S.n] = Fone
    end
    X2 = zero_matrix(S.F, 1, 4 * d^2)
    for r in 1:2:2 * d
        X2[1, (r - 1) * d + 1] = Fone
    end
    Z2 = zero_matrix(S.F, 1, 4 * d^2)
    for c in 1:d
        Z2[1, c + S.n] = Fone
    end
    S.logicals = [(X1, Z1), (X2, Z2)]
    S.dx = d
    S.dz = d
    S.d = d
    return S
end

################################
     # Planar Surface Codes
################################

"""
    PlanarSurfaceCode(dx::Int, dz::Int)
    PlanarSurfaceCode(d::Int)

Return the `[[dx * dz + (dx - 1) * (dz - 1), 1, dx/dz]]` planar surface code.

The top and bottom boundaries are "smooth" (`Z`) and the left and right are "rough" (`X`).
"""
function PlanarSurfaceCode(dx::Int, dz::Int)
    (2 <= dx && 2 <= dz) || throw(DomainError("Distances must be at least two."))

    F = GF(2)
    Fone = F(1)
    numV = dx * dz + (dx - 1) * (dz - 1)
    A = zero_matrix(F, dx * (dz - 1) + 1, numV) # stars, X stabilizers
    B = zero_matrix(F, dz * (dx - 1), numV) # faces, Z stabilizers
    qubit = 1
    rowA = 1
    rowB = 1
    for r in 1:dz
        for c in 1:dx
            if r != dz
                if c == 1
                    B[rowB, qubit] = B[rowB, qubit + dx] = B[rowB, qubit + 2 * dx - 1] = Fone
                    rowB += 1
                elseif c == dx
                    B[rowB, qubit] = B[rowB, qubit + dx - 1] = B[rowB, qubit + 2 * dx - 1] = Fone
                    rowB += 1
                else
                    B[rowB, qubit] = B[rowB, qubit + dx - 1] = B[rowB, qubit + dx] = B[rowB, qubit + 2 * dx - 1] = Fone
                    rowB += 1
                end
            end

            if c != dx
                if r == 1
                    A[rowA, qubit] = A[rowA, qubit + 1] = A[rowA, qubit + dx] = Fone
                    rowA += 1
                elseif r == dz
                    A[rowA, qubit] = A[rowA, qubit + 1] = A[rowA, qubit - dx + 1] = Fone
                    rowA += 1
                else
                    A[rowA, qubit] = A[rowA, qubit + 1] = A[rowA, qubit + dx] = A[rowA, qubit - dx + 1] = Fone
                    rowA += 1
                end
            end
            qubit += 1
        end
        qubit += dx - 1
    end

    S = CSSCode(A, B)
    X1 = zero_matrix(S.F, 1, 2 * S.n)
    for r in 1:2:dx
        X1[1, dz * (r - 1) + (dz - 1) * (r - 1) + 1] = Fone
    end
    Z1 = zero_matrix(S.F, 1, 2 * S.n)
    for c in 1:dz
        Z1[1, c + S.n] = Fone
    end
    S.logicals = [(X1, Z1)]
    S.dx = dx
    S.dz = dz
    S.d = minimum([dx, dz])
    return S
end
PlanarSurfaceCode(d::Int) = PlanarSurfaceCode(d, d)

################################
       # XY Surface Codes
################################

"""
    XYSurfaceCode(dx::Int, dz::Int)
    XYSurfaceCode(d::Int)

Return the `[[dx * dy + (dx - 1) * (dy - 1), 1, dx/dy]]` XY surface code of
"Ultrahigh Error Threshold for Surface Codes with Biased Noise" by Tuckett, Bartlett, and Flammia.

The top and bottom boundaries are "smooth" (`Y`) and the left and right are "rough" (`X`).
"""
# TODO: remove quadratic
function XYSurfaceCode(dx::Int, dy::Int)
    (2 <= dx && 2 <= dy) || throw(DomainError("Distances must be at least two."))

    F = GF(2)
    Fone = F(1)
    numV = dx * dy + (dx - 1) * (dy - 1)
    M = zero_matrix(F, numV - 1, 2 * numV)
    qubit = 1
    row = 1
    for r in 1:dy
        for c in 1:dx
            if r != dz
                if c == 1
                    M[row, qubit] = M[row, qubit + dx] = M[row, qubit + 2 * dx - 1] = Fone
                    M[row, qubit + numV] = M[row, qubit + dx + numV] = M[row, qubit + 2 * dx - 1 + numV] = Fone
                    row += 1
                elseif c == dx
                    M[row, qubit] = M[row, qubit + dx - 1] = M[row, qubit + 2 * dx - 1] = Fone
                    M[row, qubit + numV] = M[row, qubit + dx - 1 + numV] = M[row, qubit + 2 * dx - 1 + numV] = Fone
                    row += 1
                else
                    M[row, qubit] = M[row, qubit + dx - 1] = M[row, qubit + dx] = M[row, qubit + 2 * dx - 1] = Fone
                    M[row, qubit + numV] = M[row, qubit + dx - 1 + numV] = M[row, qubit + dx + numV] = M[row, qubit + 2 * dx - 1 + numV] = Fone
                    row += 1
                end
            end

            if c != dx
                if r == 1
                    M[row, qubit] = M[row, qubit + 1] = M[row, qubit + dx] = Fone
                    row += 1
                elseif r == dz
                    M[row, qubit] = M[row, qubit + 1] = M[row, qubit - dx + 1] = Fone
                    row += 1
                else
                    M[row, qubit] = M[row, qubit + 1] = M[row, qubit + dx] = M[row, qubit - dx + 1] = Fone
                    row += 1
                end
            end
            qubit += 1
        end
        qubit += dx - 1
    end
    S = StabilizerCode(M)
    # Eone = S.E(1)
    # ω = gen(S.E)
    # X1 = zero_matrix(S.E, 1, numV)
    # for r in 1:2:dx
    #     X1[1, dz * (r - 1) + (dz - 1) * (r - 1) + 1] = Eone
    # end
    # Z1 = zero_matrix(S.E, 1, numV)
    # for c in 1:dz
    #     Z1[1, c] = ω
    # end
    # S.logicals = [(X1, Z1)]
    # S.dx = dx
    # S.dz = dz
    # S.d = minimum([dx, dz])
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
#     S.d = d
#     S.dx = d
#     S.dz = 2 * d^2
#     # Y distance is also 2 * d^2
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

    F = GF(2)
    Fone = F(1)
    X = zero_matrix(F, 2, k + 4)
    Z = zero_matrix(F, 2, k + 4)
    X[1, 1] = X[1, 2] = X[1, 3] = X[1, 4] = Fone
    Z[1, 1] = Z[1, 2] = Z[1, 3] = Z[1, 4] = Fone
    X[2, 1] = X[2, 2] = Fone
    Z[2, 1] = Z[2, 2] = Fone
    for c in 5:k + 3
        X[2, c] = X[2, c + 1] = Fone
        Z[2, c] = Z[2, c + 1] = Fone
    end
    return CSSCode(X, Z)
end

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
    Compute the n-cells of a periodic d-dimensional hypercubic lattice with linear size l.
"""
function _compute_cells_periodic(l::Int, n::Int, d::Int = 4)
    cells = Set{_Cell}()    
    coords_to_change = collect(combinations(1:d, n))
    for coord in Iterators.product([0:l-1 for i=1:d]...), directions in coords_to_change
        coord = collect(coord)
        new_coords = Vector([copy(coord) for _ in 1:2^n-1])
        for i in 1:2^n-1, j in 1:n
            (i >> (j-1)) & 1 == 1 ? new_coords[i][directions[j]] += 1 : nothing
        end
        vertices = Set{_Vertex}()
        push!(vertices, _Vertex(coord))
        for new_coord in new_coords
            l > 2 ? new_coord .%= l : nothing
            push!(vertices, _Vertex(new_coord))
        end

        push!(cells, _Cell(vertices) )
    end
    return cells
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
    x_check_mat = spzeros(Bool, length(X_dict), length(q_dict) )
    z_check_mat = spzeros(Bool, length(Z_dict), length(q_dict) )

    for (e, v) in X_dict
        for c in v
            x_check_mat[edge_dict[e], q_dict[c]] = true
        end
    end

    for (v_, v) in Z_dict
        for c in v
            z_check_mat[volume_dict[v_], q_dict[c]] = true
        end
    end

    return x_check_mat, z_check_mat
end

"""
    Create vectors of logicals.
"""
function _compute_logical_vectors(X_dict, Z_dict, q_dict)
    x_logicals, z_logicals = Vector{Vector{Int}}(), Vector{Vector{Int}}()
    n = 1
    for (_, v) in X_dict
        logical = []
        for c in v
            push!(logical, q_dict[c])
        end
        push!(x_logicals, copy(logical))
    end

    for (_, v) in Z_dict
        logical = []
        for c in v
            push!(logical, q_dict[c])
        end
        push!(z_logicals, copy(logical))
    end

    return x_logicals, z_logicals
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
        for c in Iterators.product([0:l-1 for i=1:2]...)
            z_coord = [0, 0, 0, 0]
            for i in 1:length(directions)
                z_coord[directions[i]] = c[i]
            end
            other_directions = setdiff(dirs, directions)
            z_new_coords = Vector([copy(z_coord) for _ in 1:2^n-1])
            for i in 1:2^n-1, j in 1:n
                (i >> (j-1)) & 1 == 1 ? z_new_coords[i][directions[j]] += 1 : nothing
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
    vertices = _compute_cells_periodic(l, 0)
    edges = _compute_cells_periodic(l, 1)
    faces = _compute_cells_periodic(l, 2)
    volumes = _compute_cells_periodic(l, 3)
    hyper_volumes = _compute_cells_periodic(l, 4)

    Z_dict = _contains(volumes, faces, 2, l)
    X_dict = _inverse_dict(_contains(faces, edges, 1, l))

    Z_redundancy = _contains(hyper_volumes, volumes, 3, l)
    X_redundancy = _inverse_dict(_contains(edges, vertices, 0, l))

    q_dict = _build_q_dict(faces)
    volume_dict = _build_q_dict(volumes)
    edge_dict = _build_q_dict(edges)

    x_check_matrix, z_check_matrix = _compute_check_matrices(X_dict, Z_dict, q_dict, edge_dict, volume_dict)
    z_redundant = _compute_redundant(Z_redundancy, volume_dict)
    x_redundant = _compute_redundant(X_redundancy, edge_dict)

    X_logicals, Z_logicals = _compute_logicals(l)
    x_logical, z_logical = _compute_logical_vectors(X_logicals, Z_logicals, q_dict)

    # defining the code objects
    F = GF(2)
    Fone = F(1)

    X = zero_matrix(F, size(x_check_matrix)[1], size(x_check_matrix)[2])
    Z = zero_matrix(F, size(z_check_matrix)[1], size(z_check_matrix)[2])

    I, J, _ = findnz(x_check_matrix)
    for i in 1:length(I)
        X[I[i], J[i]] = Fone
    end
    I, J, _ = findnz(z_check_matrix)
    for i in 1:length(I)
        Z[I[i], J[i]] = Fone
    end

    return CSSCode(X, Z)
end