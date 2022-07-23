# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

include("tricolorcodes488trellis.jl")
include("tricolorcodes666trellis.jl")

#############################
        # Misc codes
#############################

"""
    fivequbitcode()
    Q513()

Return the `[[5, 1, 3]]` perfect qubit stabilizer code.
"""
# is a perfect code
fivequbitcode() = return QuantumCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"])
Q513() = fivequbitcode()

# should also do a test for other CSS construction via Hamming code and actually make that one default
"""
    Steanecode()
    Q713()

Return the `[[7, 1, 3]]` Steane code with stabilizers in standard ordering.
"""
SteaneCode() = return CSSCode(["XXXXIII", "XXIIXXI", "XIXIXIX", "ZZZZIII", "ZZIIZZI", "ZIZIZIZ"])
Q713() = SteaneCode()
_SteaneCodetrellis() = return CSSCode(["XXIXXII", "IXXIXXI", "IIIXXXX", "ZZIZZII", "IZZIZZI", "IIIZZZZ"])
# also ZZIZZII, ZIZZIZI, IZZZIIZ, XXIXXII, XIXXIXI, IXXXIIX

"""
    Shorcode()
    Q913()

Return the `[[9, 1, 3]]` Shor code.
"""
ShorCode() = return CSSCode(["ZZIIIIIII", "IZZIIIIII", "IIIZZIIII", "IIIIZZIII", "IIIIIIZZI", "IIIIIIIZZ", "XXXXXXIII", "IIIXXXXXX"])
Q913() = ShorCode()

Q412() = return CSSCode(["XXXX", "ZZII", "IIZZ"])
Q422() = return CSSCode(["XXXX", "ZZZZ"])
Q511() = return QuantumCode(["ZXIII", "XZXII", "IXZXI", "IIXZX"])

function Q823()
    F, _ = FiniteField(2, 1, "α")
    S = matrix(F, [1 0 0 0 1 0 0 0 1 1 1 1 0 0 0 0;
    0 0 0 1 0 1 0 0 1 0 0 0 0 1 0 0;
    0 1 0 0 1 1 1 0 0 0 1 1 1 0 1 0;
    0 0 1 0 1 1 1 0 0 1 1 0 1 1 0 0;
    0 0 1 1 1 0 1 0 0 0 0 1 0 1 1 1;
    0 0 0 0 0 0 1 1 0 0 1 0 0 0 1 0]);
    return QuantumCode(S, true)
end

"""
    Q15RM()
    Q1513()

Return the `[[15, 1, 3]]` quantum Reed-Muller code with stabilizers in standard
ordering.
"""
Q15RM() = return QuantumCode(["ZIZIZIZIZIZIZIZ", "IZZIIZZIIZZIIZZ", "IIIZZZZIIIIZZZZ",
    "IIIIIIIZZZZZZZZ", "IIZIIIZIIIZIIIZ", "IIIIZIZIIIIIZIZ", "IIIIIZZIIIIIIZZ",
    "IIIIIIIIIZZIIZZ", "IIIIIIIIIIIZZZZ", "IIIIIIIIZIZIZIZ", "XIXIXIXIXIXIXIX",
    "IXXIIXXIIXXIIXX", "IIIXXXXIIIIXXXX", "IIIIIIIXXXXXXXX"])
Q1513() = Q15RM()

"""
    Q1573()

Return the `[[15, 7, 3]]` quantum Hamming code.
"""
Q1573() = return QuantumCode(["IIIIIIIXXXXXXXX", "IIIXXXXIIIIXXXX", "IXXIIXXIIXXIIXX",
    "XIXIXIXIXIXIXIX", "IIIIIIIZZZZZZZZ", "IIIZZZZIIIIZZZZ", "IZZIIZZIIZZIIZZ",
    "ZIZIZIZIZIZIZIZ"])
    # one can use a basis for this such that the first logical pair is transversal X, Z

#############################
#  Triangular Surface Codes #
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

function _triangularlatticeXstabilizers(L::Int, numbering::Array{Int64, 3}, symp::Bool=true)
    F, _ = FiniteField(2, 1, "α")
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
    if symp
        return hcat(stabilizers, zero_matrix(F, L^2, 3 * L^2))
    end
    return stabilizers
end

function _triangularlatticeZstabilizers(L::Int, numbering::Array{Int64, 3}, symp::Bool=true)
    F, _ = FiniteField(2, 1, "α")
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
    if symp
        return hcat(zero_matrix(F, 2 * L^2, 3 * L^2), stabilizers)
    end
    return stabilizers
end

function _triangularlatticeXlogicals(L::Int, numbering::Array{Int64, 3}, symp::Bool=true)
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
    if symp
        logical1 = [logical1; z]
    end

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
    if symp
        logical2 = [logical2; z]
    end

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
    if symp
        logical1 = [z; logical1]
    end

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
    if symp
        logical2 = [z; logical2]
    end

  return [logical1, logical2]
end

function triangularsurfacecode(L::Int)
    numbering = _triangularlattice(L)
    Xstabs = _triangularlatticeXstabilizers(L, numbering, false)
    # println(rank(Xstabs))
    Zstabs = _triangularlatticeZstabilizers(L, numbering, false)
    # println(Zstabs)
    # logicals = [triangularlatticeXlogicals(L, numbering), triangularlatticeZlogicals(L, numbering)]
    return CSSCode(Xstabs[1:end - 1, :], Zstabs[1:end - 1, :])
end

#############################
#  Rotated Surface Codes #
#############################

function _RSurfstabslogs(d::Int)
    n = d^2
    E, ω = FiniteField(2, 2, "ω")
    X = E(1)
    S = zero_matrix(E, n - 1, n)
    row = 1

    # X's
    i = 1
    while i <= n - d
        S[row, i] = X
        S[row, i + 1] = X
        S[row, i + d] = X
        S[row, i + d + 1] = X
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
        S[row, i] = X
        S[row, i + 1] = X
        row += 1
        i += 2
    end

    # bottom row X's
    i = d * (d - 1) + 1
    while i <= d * d - 2
        S[row, i] = X
        S[row, i + 1] = X
        row += 1
        i += 2
    end

    # Z's
    i = 2
    while i < n - d
        S[row, i] = ω
        S[row, i + 1] = ω
        S[row, i + d] = ω
        S[row, i + d + 1] = ω
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
        S[row, i] = ω
        S[row, i + d] = ω
        row += 1
        i += 2 * d
    end

    # right Z's
    i = 2 * d
    while i < d * d
        S[row, i] = ω
        S[row, i + d] = ω
        row += 1
        i += 2 * d
    end

    logs = zero_matrix(E, 2, n)
    i = d
    while i <= d * d
        logs[1, i] = 0x01
        i += d
    end
    i = 1
    while i <= d
        logs[2, i] = ω
        i += 1
    end

    return S, logs
end

# flint is being ridiculous here
# julia> @time Q = rotatedsurfacecode(21);
# 226.128107 seconds (1.73 G allocations: 80.659 GiB, 18.33% gc time, 0.09% compilation time)
# julia> Base.summarysize(Q)
# 16872
"""
    rotatedsurfacecode(d::Int)

Return the `[[d^2, 1, d]]` rotated surface code.
"""
function rotatedsurfacecode(d::Int)
    d >= 3 || error("Current implementation requires d ≥ 3.")

    S, logs = _RSurfstabslogs(d)
    Q = QuantumCode(S)
    setlogicals!(Q, logs)
    return Q
end

#############################
#     XZZX Surface Codes    #
#############################

function _XZZXstabslogs(d::Int)
    n = d^2
    E, ω = FiniteField(2, 2, "ω")
    S = zero_matrix(E, n - 1, n)
    row = 1
    X = E(1)

    i = 1
    for i in 1:n - d
        if i % d != 0
            S[row, i] = X
            S[row, i + 1] = ω
            S[row, i + d] = ω
            S[row, i + d + 1] = X
            row += 1;
        end
    end

    # top row ZX's
    i = 2
    while i <= d - 1
        S[row, i] = ω
        S[row, i + 1] = X
        row += 1
        i += 2
    end

    # bottom row XZ's
    i = d * (d - 1) + 1
    while i <= d * d - 2
        S[row, i] = X
        S[row, i + 1] = ω
        row += 1
        i += 2
    end

    # left ZX's
    i = 1
    while i < d * (d - 1)
        S[row, i] = ω
        S[row, i + d] = X
        row += 1
        i += 2 * d
    end

    # right XZ's
    i = 2 * d
    while i < d * d
        S[row, i] = X
        S[row, i + d] = ω
        row += 1
        i += 2 * d
    end

    logs = zero_matrix(E, 2, n)
    i = d
    count = 1
    while i <= d * d
        if count % 2 == 1
            logs[1, i] = X
        else
            logs[1, i] = ω
        end
        i += d
        count += 1
    end
    i = 1
    count = 1
    while i <= d
        if count % 2 == 1
            logs[2, i] = ω
        else
            logs[2, i] = X
        end
        i += 1
        count += 1
    end

    return S, logs
end

"""
    XZZXsurfacecode(d::Int)

Return the `[[d^2, 1, d]]` XZZX surface code.
"""
function XZZXsurfacecode(d::Int)
    d >= 3 || error("Current implementation requires d ≥ 3.")

    S, logs = _XZZXstabslogs(d)
    Q = QuantumCode(S)
    setlogicals!(Q, logs)
    return Q
end

################################
# Triangular Color Codes 4.8.8 #
################################

"""
    tricolorcode488(d::Int)

Return the 4.8.8 triangular color code of distance `d` with trellis numbering.
"""
function tricolorcode488(d::Int)
    3 <= d <= 21 || error("Current implementation requires 3 ≤ d ≤ 21.")

    if d == 3
        S, logs = _488d3trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 5
        S, logs = _488d5trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 7
        S, logs = _488d7trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 9
        S, logs = _488d9trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 11
        S, logs = _488d11trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 13
        S, logs = _488d13trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 15
        S, logs = _488d15trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 17
        S, logs = _488d17trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 19
        S, logs = _488d19trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 21
        S, logs = _488d21trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    end
end

################################
# Triangular Color Codes 6.6.6 #
################################

"""
    tricolorcode666(d::Int)

Return the 6.6.6 triangular color code of distance `d` with trellis numbering.
"""
function tricolorcode666(d::Int)
    3 <= d <= 21 || error("Current implementation requires 3 ≤ d ≤ 21.")

    if d == 3
        # same as 4.8.8
        S, logs = _488d3trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 5
        S, logs = _666d5trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 7
        S, logs = _666d7trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 9
        S, logs = _666d9trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 11
        S, logs = _666d11trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 13
        S, logs = _666d13trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 15
        S, logs = _666d15trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 17
        S, logs = _666d17trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 19
        S, logs = _666d19trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    elseif d == 21
        S, logs = _666d21trellis()
        Q = QuantumCode(S)
        setlogicals!(Q, logs)
        return Q
    end
end
