# Copyright (c) 2022, 2023 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
    GraphState(G::SimpleGraph{Int64})

Return the graph state defined by the graph `G`.
"""
function GraphState(G::SimpleGraph{Int})
    # probably need some checks on G here but maybe the function args are good enough
    A = adjacency_matrix(G)
    _, nc = size(A)
    for i = 1:nc
        iszero(A[i, i]) || throw(ArgumentError("Graph cannot have self-loops."))
    end
    # are there non-binary graph states?
    F = Oscar.Nemo.Native.GF(2)
    Fone = F(1)
    sym_stabs = zero_matrix(F, nc, 2 * nc)
    for r = 1:nc
        sym_stabs[r, r] = Fone
        for c = 1:nc
            isone(A[r, c]) && (sym_stabs[r, c+nc] = Fone;)
        end
    end
    # this should automatically compute everything for the GraphState constructor
    return StabilizerCode(sym_stabs)
end

"""
    ClusterState(w::Int, h::Int)

Return the cluster state (graph state) on the rectangular lattice with width `w` and height `h`.
"""
function ClusterState(w::Int, h::Int)
    (0 <= w && 0 <= h) || throw(ArgumentError("Rectangle dimensions must be positive."))

    # BUG this function has never been updated from the quadratic form
    E = GF(2, 2, :ω)
    ω = gen(E)
    Eone = E(1)
    A = zero_matrix(E, w * h, w * h)
    curr = 1
    for r = 1:h
        for c = 1:w
            A[curr, curr] = Eone
            c != 1 && (A[curr, curr-1] = ω;)
            c != w && (A[curr, curr+1] = ω;)
            r != 1 && (A[curr, curr-w] = ω;)
            r != h && (A[curr, curr+w] = ω;)
            curr += 1
        end
    end
    return StabilizerCode(A)
end
