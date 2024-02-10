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

    if d == 3
        # S, logs = _488d3trellis()
        # print(pwd())
        #TODO: these links are relative to the current path and not the file structure
        # @load "../../../data/488d3stabslogs_trellis.jld2" S l
        @load "data/488d3stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    elseif d == 5
        # S, logs = _488d5trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/488d5stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    elseif d == 7
        # S, logs = _488d7trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/488d7stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    elseif d == 9
        # S, logs = _488d9trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/488d9stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    elseif d == 11
        # S, logs = _488d11trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/488d11stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    elseif d == 13
        # S, logs = _488d13trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/488d13stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    elseif d == 15
        # S, logs = _488d15trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/488d15stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    elseif d == 17
        # S, logs = _488d17trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/488d17stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    elseif d == 19
        # S, logs = _488d19trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/488d19stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    # elseif d == 21
    #     # S, logs = _488d21trellis()
    #     # S = StabilizerCode(stabs)
    #     # set_logicals!(Q, logs)
    #     @load "../data/488d21stabslogs_trellis.jld2" S l
    #     F = CodingTheory.Oscar.GF(2)
    #     stabs = CodingTheory.Oscar.matrix(F, S)
    #     S = StabilizerCode(stabs)
    #     l = CodingTheory.Oscar.matrix(F, l)
    #     # set_logicals!(S, [(l[1, :], l[2, :])])
    #     set_logicals!(S, l)
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
function CodingTheory.TriangularColorCode666(d::Int)
    3 <= d <= 21 || throw(DomainError("Current implementation requires 3 ≤ d ≤ 21."))

    if d == 3
        # same as 4.8.8
        # S, logs = _488d3trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/488d3stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    elseif d == 5
        # S, logs = _666d5trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/666d5stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    elseif d == 7
        # S, logs = _666d7trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/666d7stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    elseif d == 9
        # S, logs = _666d9trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/666d9stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    elseif d == 11
        # S, logs = _666d11trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/666d11stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    elseif d == 13
        # S, logs = _666d13trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/666d13stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    elseif d == 15
        # S, logs = _666d15trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/666d15stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    elseif d == 17
        # S, logs = _666d17trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/666d17stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    elseif d == 19
        # S, logs = _666d19trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/666d19stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    elseif d == 21
        # S, logs = _666d21trellis()
        # S = StabilizerCode(stabs)
        # set_logicals!(Q, logs)
        @load "../../../data/666d21stabslogs_trellis.jld2" S l
        F = CodingTheory.Oscar.GF(2)
        stabs = CodingTheory.Oscar.matrix(F, S)
        S = StabilizerCode(stabs)
        l = CodingTheory.Oscar.matrix(F, l)
        # set_logicals!(S, [(l[1, :], l[2, :])])
        set_logicals!(S, l)
        return S
    end
end
