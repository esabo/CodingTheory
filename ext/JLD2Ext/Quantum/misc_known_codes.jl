# Copyright (c) 2021 - 2024, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

################################
 # Triangular Color Codes 4.8.8
################################

"""
    TriangularColorCode488(d::Int)

Return the 4.8.8 triangular color code of distance `d` with trellis numbering.
"""
function CodingTheory.TriangularColorCode488(d::Int)
    3 <= d <= 19 || throw(DomainError("Current implementation requires 3 ≤ d ≤ 19."))

    path = joinpath(@__DIR__, "../../../data/color_codes")

    if d == 3
        @load joinpath(path, "488d3stabslogs_trellis.jld2") S l
    elseif d == 5
        @load joinpath(path, "488d5stabslogs_trellis.jld2") S l
    elseif d == 7
        @load joinpath(path, "488d7stabslogs_trellis.jld2") S l
    elseif d == 9
        @load joinpath(path, "488d9stabslogs_trellis.jld2") S l
    elseif d == 11
        @load joinpath(path, "488d11stabslogs_trellis.jld2") S l
    elseif d == 13
        @load joinpath(path, "488d13stabslogs_trellis.jld2") S l
    elseif d == 15
        @load joinpath(path, "488d15stabslogs_trellis.jld2") S l
    elseif d == 17
        @load joinpath(path, "488d17stabslogs_trellis.jld2") S l
    elseif d == 19
        @load joinpath(path, "488d19stabslogs_trellis.jld2") S l
    end
    F = CodingTheory.Oscar.GF(2)
    Q = StabilizerCode(CodingTheory.Oscar.matrix(F, S))
    set_logicals!(Q, CodingTheory.Oscar.matrix(F, l))
    return Q
end

################################
 # Triangular Color Codes 6.6.6
################################

"""
    TriangularColorCode666(d::Int)

Return the 6.6.6 triangular color code of distance `d` with trellis numbering.
"""
function CodingTheory.TriangularColorCode666(d::Int)
    3 <= d <= 21 || throw(DomainError("Current implementation requires 3 ≤ d ≤ 21."))

    path = joinpath(@__DIR__, "../../../data/color_codes")

    if d == 3
        # same as 4.8.8 d= 3
        @load joinpath(path, "488d3stabslogs_trellis.jld2") S l
    elseif d == 5
        @load joinpath(path, "666d5stabslogs_trellis.jld2") S l
    elseif d == 7
        @load joinpath(path, "666d7stabslogs_trellis.jld2") S l
    elseif d == 9
        @load joinpath(path, "666d9stabslogs_trellis.jld2") S l
    elseif d == 11
        @load joinpath(path, "666d11stabslogs_trellis.jld2") S l
    elseif d == 13
        @load joinpath(path, "666d13stabslogs_trellis.jld2") S l
    elseif d == 15
        @load joinpath(path, "666d15stabslogs_trellis.jld2") S l
    elseif d == 17
        @load joinpath(path, "666d17stabslogs_trellis.jld2") S l
    elseif d == 19
        @load joinpath(path, "666d19stabslogs_trellis.jld2") S l
    elseif d == 21
        @load joinpath(path, "666d21stabslogs_trellis.jld2") S l    
    end
    F = CodingTheory.Oscar.GF(2)
    Q = StabilizerCode(CodingTheory.Oscar.matrix(F, S))
    set_logicals!(Q, CodingTheory.Oscar.matrix(F, l))
    return Q
end

################################
     # 3D PlanarSurfaceCode
################################

# TODO missing Z matrices
"""
    PlanarSurfaceCode3D(d::Int)

Return the 3D planar surface code of distance `d`.

# Note
- Run `using JLD2` to activate this extension.
- For the moment, these are not computed but loaded from file (from MikeVasmer) and are limited to
  `3 ≤ d ≤ 9`.
"""
function CodingTheory.PlanarSurfaceCode3D_X(d::Int)
    3 ≤ d ≤ 9 || throw(DomainError("Current implementation requires 3 ≤ d ≤ 9."))

    path = joinpath(@__DIR__, "../../../data/surface3D")
    if d == 3
        X_stabs = load_alist(joinpath(path, "surface3D_3_hx.alist"))
        X_logs = load_alist(joinpath(path, "surface3D_3_lx.alist"))
        X_meta = load_alist(joinpath(path, "surface3D_3_mx.alist"))
    elseif d == 4
        X_stabs = load_alist(joinpath(path, "surface3D_4_hx.alist"))
        X_logs = load_alist(joinpath(path, "surface3D_4_lx.alist"))
        X_meta = load_alist(joinpath(path, "surface3D_4_mx.alist"))
    elseif d == 5
        X_stabs = load_alist(joinpath(path, "surface3D_5_hx.alist"))
        X_logs = load_alist(joinpath(path, "surface3D_5_lx.alist"))
        X_meta = load_alist(joinpath(path, "surface3D_5_mx.alist"))
    elseif d == 6
        X_stabs = load_alist(joinpath(path, "surface3D_6_hx.alist"))
        X_logs = load_alist(joinpath(path, "surface3D_6_lx.alist"))
        X_meta = load_alist(joinpath(path, "surface3D_6_mx.alist"))
    elseif d == 7
        X_stabs = load_alist(joinpath(path, "surface3D_7_hx.alist"))
        X_logs = load_alist(joinpath(path, "surface3D_7_lx.alist"))
        X_meta = load_alist(joinpath(path, "surface3D_7_mx.alist"))
    elseif d == 8
        X_stabs = load_alist(joinpath(path, "surface3D_8_hx.alist"))
        X_logs = load_alist(joinpath(path, "surface3D_8_lx.alist"))
        X_meta = load_alist(joinpath(path, "surface3D_8_mx.alist"))
    elseif d == 9
        X_stabs = load_alist(joinpath(path, "surface3D_9_hx.alist"))
        X_logs = load_alist(joinpath(path, "surface3D_9_lx.alist"))
        X_meta = load_alist(joinpath(path, "surface3D_9_mx.alist"))
    end
    F = CodingTheory.Oscar.GF(2)
    return CodingTheory.Oscar.matrix(F, X_stabs), CodingTheory.Oscar.matrix(F, X_logs),
        CodingTheory.Oscar.matrix(F, X_meta)
end

#################################
        # 3D Toric codes
#################################

# TODO missing Z matrices
"""
    ToricCode3D(d::Int)

Return the 3D toric code of distance `d`.

# Note
- Run `using JLD2` to activate this extension.
- For the moment, these are not computed but loaded from file (from MikeVasmer) and are limited to
  `2 ≤ d ≤ 13`.
"""
function CodingTheory.ToricCode3D_X(d::Int)
    2 ≤ d ≤ 13 || throw(DomainError("Current implementation requires 2 ≤ d ≤ 13."))

    path = joinpath(@__DIR__, "../../../data/toric3D")
    if d == 2
        X_stabs = load_alist(joinpath(path, "toric3D_2_hx.alist"))
        X_logs = load_alist(joinpath(path, "toric3D_2_lx.alist"))
        X_meta = load_alist(joinpath(path, "toric3D_2_mx.alist"))
    elseif d == 3
        X_stabs = load_alist(joinpath(path, "toric3D_3_hx.alist"))
        X_logs = load_alist(joinpath(path, "toric3D_3_lx.alist"))
        X_meta = load_alist(joinpath(path, "toric3D_3_mx.alist"))
    elseif d == 4
        X_stabs = load_alist(joinpath(path, "toric3D_4_hx.alist"))
        X_logs = load_alist(joinpath(path, "toric3D_4_lx.alist"))
        X_meta = load_alist(joinpath(path, "toric3D_4_mx.alist"))
    elseif d == 5
        X_stabs = load_alist(joinpath(path, "toric3D_5_hx.alist"))
        X_logs = load_alist(joinpath(path, "toric3D_5_lx.alist"))
        X_meta = load_alist(joinpath(path, "toric3D_5_mx.alist"))
    elseif d == 6
        X_stabs = load_alist(joinpath(path, "toric3D_6_hx.alist"))
        X_logs = load_alist(joinpath(path, "toric3D_6_lx.alist"))
        X_meta = load_alist(joinpath(path, "toric3D_6_mx.alist"))
    elseif d == 7
        X_stabs = load_alist(joinpath(path, "toric3D_7_hx.alist"))
        X_logs = load_alist(joinpath(path, "toric3D_7_lx.alist"))
        X_meta = load_alist(joinpath(path, "toric3D_7_mx.alist"))
    elseif d == 8
        X_stabs = load_alist(joinpath(path, "toric3D_8_hx.alist"))
        X_logs = load_alist(joinpath(path, "toric3D_8_lx.alist"))
        X_meta = load_alist(joinpath(path, "toric3D_8_mx.alist"))
    elseif d == 9
        X_stabs = load_alist(joinpath(path, "toric3D_9_hx.alist"))
        X_logs = load_alist(joinpath(path, "toric3D_9_lx.alist"))
        X_meta = load_alist(joinpath(path, "toric3D_9_mx.alist"))
    elseif d == 10
        X_stabs = load_alist(joinpath(path, "toric3D_10_hx.alist"))
        X_logs = load_alist(joinpath(path, "toric3D_10_lx.alist"))
        X_meta = load_alist(joinpath(path, "toric3D_10_mx.alist"))
    elseif d == 11
        X_stabs = load_alist(joinpath(path, "toric3D_11_hx.alist"))
        X_logs = load_alist(joinpath(path, "toric3D_11_lx.alist"))
        X_meta = load_alist(joinpath(path, "toric3D_11_mx.alist"))
    elseif d == 12
        X_stabs = load_alist(joinpath(path, "toric3D_12_hx.alist"))
        X_logs = load_alist(joinpath(path, "toric3D_12_lx.alist"))
        X_meta = load_alist(joinpath(path, "toric3D_12_mx.alist"))
    elseif d == 13
        X_stabs = load_alist(joinpath(path, "toric3D_13_hx.alist"))
        X_logs = load_alist(joinpath(path, "toric3D_13_lx.alist"))
        X_meta = load_alist(joinpath(path, "toric3D_13_mx.alist"))
    end
    F = CodingTheory.Oscar.GF(2)
    return CodingTheory.Oscar.matrix(F, X_stabs), CodingTheory.Oscar.matrix(F, X_logs),
        CodingTheory.Oscar.matrix(F, X_meta)
end
