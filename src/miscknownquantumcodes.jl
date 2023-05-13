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

"""
    BaconShorCode(m::Int, n::Int)

Return the Bacon-Shor subsystem code on a `m x n` lattice.
"""
function BaconShorCode(m::Int, n::Int)
    # F, _ = FiniteField(2, 1, "α")
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
# Bravyi, "Subsystem Codes With Spatially Local Generators", (2011)
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

    totXgauges = 0

    # TODO: really only need to generate every consequetive pair

    # every pair of 1's in row of A gets an X gauge operator
    Xgauges = zero_matrix(F, totXgauges, n)
    currrow = 1
    Fone = F(1)
    for r in 1:nr
        for c1 in 1:nc - 1
            if !iszero(A[r, c1])
                for c2 in c1 + 1:nc
                    if !iszero(A[r, c2])
                        Xgauges[currrow, linearindex[r, c1]] = Fone
                        Xgauges[currrow, linearindex[r, c2]] = Fone
                        currrow += 1
                    end
                end
            end
        end
    end

    # every pair of 1's in col of A gets a Z gauge operator
    Zgauges = zero_matrix(F, totZgauges, n)
    currrow = 1
    for c in 1:nc
        for r1 in 1:nr - 1
            if !iszero(A[r1, c])
                for r2 in r1 + 1:nr
                    if !iszero(A[r2, c])
                        Zgauges[currrow, linearindex[r1, c]] = Fone
                        Zgauges[currrow, linearindex[r2, c]] = Fone
                        currrow += 1
                    end
                end
            end
        end
    end
    S = SubsystemCode(Xgauges ⊕ Zgauges)
    minrowwt = minimum(rowwts)
    mincolwt = minimum(colwts)
    setminimumdistance!(S, minimum([minrowwt, mincolwt]))
    # TODO: also set dx and dz
    return S
end
GeneralizedBaconShorCode(A::CTMatrixTypes) = BravyiSubsystemCode(A)

# subsystem codes was described by Bacon and Casaccino in [19]. The construction of [19] starts from a pair of classical linear codes C1 = [n1, k1, d1] and C2 = [n2, k2, d2]. A quantum subsystem code is then defined by placing a physical qubit at every cell of a ma- trix A of size n1 × n2. The X-part of the gauge group is defined by replicating the parity checks of C1 in every col- umn of A (in the X-basis). Similarly, the Z-part of the gauge group is defined by replicating the parity checks of C2 in every row of A (in the Z-basis). The resulting subsystem code has parameters [n1n2, k1k2, min (d1, d2)].

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
    # F, _ = FiniteField(2, 1, "α")
    F = GF(2)
    S = matrix(F, [1 0 0 0 1 0 0 0 1 1 1 1 0 0 0 0;
    0 0 0 1 0 1 0 0 1 0 0 0 0 1 0 0;
    0 1 0 0 1 1 1 0 0 0 1 1 1 0 1 0;
    0 0 1 0 1 1 1 0 0 1 1 0 1 1 0 0;
    0 0 1 1 1 0 1 0 0 0 0 1 0 1 1 1;
    0 0 0 0 0 0 1 1 0 0 1 0 0 0 1 0]);
    return StabilizerCode(S)
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
    # F, _ = FiniteField(2, 1, "α")
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
    # F, _ = FiniteField(2, 1, "α")
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

function _RSurfstabslogs(d::Int)
    n = d^2
    # F, _ = FiniteField(2, 1, "ω")
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

    return S, logs
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

    stabs, logs = _RSurfstabslogs(d)
    S = StabilizerCode(stabs)
    setlogicals!(S, logs)
    return S
end

#############################
     # XZZX Surface Codes
#############################

function _XZZXstabslogs(d::Int)
    n = d^2
    # F, _ = FiniteField(2, 1, "ω")
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
        # F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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
    #     F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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
        # F, _ = FiniteField(2, 1, "α")
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

    # F, _ = FiniteField(2, 1, "α")
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

    # F, _ = FiniteField(2, 1, "α")
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

    # F, _ = FiniteField(2, 1, "ω")
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

#     E, ω = FiniteField(2, 2, "ω")
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

    # F, _ = FiniteField(2, 1, "α")
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
