# Copyright (c) 2022 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

mutable struct GraphState <: AbstractGraphState
    F::FqNmodFiniteField # base field (symplectic)
    E::FqNmodFiniteField # additive field
    n::Int
    k::Int
    d::Union{Int, Missing}
    stabs::fq_nmod_mat
    charvec::Vector{nmod}
    signs::Vector{nmod}
    wtenum::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    overcomplete::Bool
end

mutable struct GraphStateCSS <: AbstractGraphStateCSS
    F::FqNmodFiniteField # base field (symplectic)
    E::FqNmodFiniteField # additive field
    n::Int
    k::Int
    d::Union{Int, Missing}
    dx::Union{Int, Missing}
    dz::Union{Int, Missing}
    stabs::fq_nmod_mat
    Xstabs::fq_nmod_mat
    Zstabs::fq_nmod_mat
    Xorigcode::Union{LinearCode, Missing}
    ZorigCode::Union{LinearCode, Missing}
    signs::Vector{nmod}
    Xsigns::Vector{nmod}
    Zsigns::Vector{nmod}
    charvec::Vector{nmod}
    wtenum::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    overcomplete::Bool
end

function graphstate(G::SimpleGraph{Int64})
    # probably need some checks on G here but maybe the function args are good enough
    A = adjacency_matrix(G)
    _, nc = size(A)
    for i in 1:nc
        iszero(A[i, i]) || error("Graph cannot have self-loops.")
    end
    # are there non-binary graph states?
    F, _ = FiniteField(2, 1, "α")
    fone = F(1)
    symstabs = zero_matrix(F, nc, 2 * nc)
    for r in 1:nc
        symstabs[r, r] = fone
        for c in 1:nc
            isone(A[r, c]) && (symstabs[r, c + nc] = fone;)
        end
    end
    # this should automatically compute everything for the GraphState constructor
    return QuantumCode(symstabs, true, missing)
end

"""
    ClusterState(w::Int, h::Int)

Return the cluster state (graph state) on the rectangular lattice with width `w` and height `h`.
"""
function ClusterState(w::Int, h::Int)
    (0 <= w && 0 <= h) || throw(ArgumentError("Rectangle dimensions must be positive."))

    E, ω = FiniteField(2, 2, "ω")
    Eone = E(1)
    A = zero_matrix(E, w * h, w * h)
    curr = 1
    for r in 1:h
        for c in 1:w
            A[curr, curr] = Eone
            c != 1 && (A[curr, curr - 1] = ω;)
            c != w && (A[curr, curr + 1] = ω;)
            r != 1 && (A[curr, curr - w] = ω;)
            r != h && (A[curr, curr + w] = ω;)
            curr += 1
        end
    end
    return QuantumCode(A, false, missing)
end

function show(io::IO, S::Union{GraphState, GraphStateCSS})
    if typeof(S) == GraphStateCSS
        if ismissing(S.d)
            println(io, "[[$(S.n), 0]]_$(order(S.F)) CSS graph state.")
        else
            println(io, "[[$(S.n), 0, $(S.d)]]_$(order(S.F)) CSS graph state.")
        end
    else
        if ismissing(S.d)
            println(io, "[[$(S.n), 0]]_$(order(S.F)) graph state.")
        else
            println(io, "[[$(S.n), 0, $(S.d)]]_$(order(S.F)) graph state.")
        end
    end
    if !get(io, :compact, false)
        if typeof(S) == GraphStateCSS
            if S.overcomplete
                println(io, "X-stabilizer matrix (overcomplete): $(numXstabs(S)) × $(S.n)")
            else
                println(io, "X-stabilizer matrix: $(numXstabs(S)) × $(S.n)")
            end
            for i in 1:numXstabs(S)
                print(io, "\t chi($(S.Xsigns[i])) ")
                for j in 1:S.n
                    if j != S.n
                        print(io, "$(S.Xstabs[i, j]) ")
                    elseif j == S.n && i != S.n
                        println(io, "$(S.Xstabs[i, j])")
                    else
                        print(io, "$(S.Xstabs[i, j])")
                    end
                end
            end
            if isovercomplete(S)
                println(io, "Z-stabilizer matrix (overcomplete): $(numZstabs(S)) × $(S.n)")
            else
                println(io, "Z-stabilizer matrix: $(numZstabs(S)) × $(S.n)")
            end
            for i in 1:numZstabs(S)
                print(io, "\t chi($(S.Zsigns[i])) ")
                for j in 1:S.n
                    if j != S.n
                        print(io, "$(S.Zstabs[i, j]) ")
                    elseif j == S.n && i != S.n
                        println(io, "$(S.Zstabs[i, j])")
                    else
                        print(io, "$(S.Zstabs[i, j])")
                    end
                end
            end
        else
            if isovercomplete(S)
                println(io, "Stabilizer matrix (overcomplete): $(nrows(S.stabs)) × $(S.n)")
            else
                println(io, "Stabilizer matrix: $(nrows(S.stabs)) × $(S.n)")
            end
            for i in 1:nrows(S.stabs)
                print(io, "\t chi($(S.signs[i])) ")
                for j in 1:S.n
                    if j != S.n
                        print(io, "$(S.stabs[i, j]) ")
                    elseif j == S.n && i != S.n
                        println(io, "$(S.stabs[i, j])")
                    else
                        print(io, "$(S.stabs[i, j])")
                    end
                end
            end
        end
    end
end
