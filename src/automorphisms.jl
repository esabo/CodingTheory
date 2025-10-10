# Copyright (c) 2025 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
       # Tanner Graphs
#############################

function Tanner_graph_automorphism_group(H::CTMatrixTypes)
    nr, nc = size(H)
    A = vcat(hcat(zero_matrix(base_ring(H), nc, nc), transpose(H)), hcat(H, zero_matrix(base_ring(H), nr, nr)))
    G = graph_from_adjacency_matrix(Undirected, A)
    Aut_gens = automorphism_group_generators(G)

    # with the convention above, 1:nc are the var nodes, nc + 1:end are the check nodes

    return Aut_gens
end

Tanner_graph_automorphism_group(L::AbstractLDPCCode) = Tanner_graph_automorphism_group(L.H)

# I think this needs to be different
Tanner_graph_automorphism_group(S::AbstractStabilizerCode) = Tanner_graph_automorphism_group(S.stabs)

# these need to pull out the zero's in the symplectic format
Tanner_graph_automorphism_group_X(S::AbstractStabilizerCSSCode) = Tanner_graph_automorphism_group(S.X_stabs)

Tanner_graph_automorphism_group_Z(S::AbstractStabilizerCSSCode) = Tanner_graph_automorphism_group(S.Z_stabs)



_foiliated_stabilizer_matrix(H::CodingTheory.CTMatrixTypes, d::Int) = 
    hcat(identity_matrix(base_ring(H), d) ⊗ H, CodingTheory._rep_pcm_tr(base_ring(H), d) ⊗
    identity_matrix(base_ring(H), Oscar.nrows(H)))
