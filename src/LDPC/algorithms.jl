# # Copyright (c) 2022 Eric Sabo
# # All rights reserved.
# #
# # This source code is licensed under the BSD-style license found in the
# # LICENSE file in the root directory of this source tree.

# function PEGalg(numvar::Int, numchk::Int, vardegdist::Vector{Int})
#     numvar == length(vardegdist) || throw(ArgumentError("Vector length must match number of variable nodes."))




#     create empty PCM H
#     also store graph?

#     what data structure is required to efficently build the subgraph
#     need to be able to immediately know which vertices connect to which
#     if we have to have some overhead to get that format, so be it
#     therefore, we will use the adjacency list format
#     we don't need the edges to be sorted, so we will implement ourselves instead of usin Graphs.jl
#     the variable node lists can be preallocated since we know the degree distribution
#     but the check node lists can be higher since the number of these is less than the number of variables
#     need to come up with an estimated average check node degree to preallocate


#     don't want the subgraph to be stored in an array because it will constantly grow and shrink
#     linked-list based structure most appropriate then
#     start with a vertex and create leaves for everything in adj_list[v]
#     simultaneously create a RB tree out of these nodes
#     check node: Int ID, Int dist, Int currdeg, Vector{Int} vertices
#         can most likely skip the latter and simply point to location in other array given ID
#     loop over this next layer and add them as children and make a RB tree out of them
#     arrange RB trees in each layer into a single array of RB tree pointers

# end
