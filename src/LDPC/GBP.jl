# Copyright (c) 2024 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
       # Region Graphs
#############################

mutable struct Region
    # maybe rename to label
    id::Vector{Int}
    # edges
    parents::Vector{Region}
    ancestors::Vector{Region}
    subregions::Vector{Region}
    overcounting_number::Int
    # size of parents, element is indices of parent marginalized to get child (this) 
    missing_parent_indices::Vector{Vector{Int}}
    # size of parents, element is numerator for each edge in (region index, parent index) pairs
    message_numerators::Vector{Vector{Tuple{Int, Int}}}
    # size of parents, element is denominator for each edge in (region index, parent index) pairs
    message_denominators::Vector{Vector{Tuple{Int, Int}}}
end

Region(id::Vector{Int}) = Region(id, Vector{Region}(), Vector{Region}(), Vector{Region}(), 1,
    Vector{Vector{Int}}(), Vector{Vector{Tuple{Int, Int}}}(), Vector{Vector{Tuple{Int, Int}}}())
Region(id::Vector{Int}, c_r::Int) = Region(id, Vector{Region}(), Vector{Region}(),
    Vector{Region}(), c_r, Vector{Vector{Int}}(), Vector{Vector{Tuple{Int, Int}}}(),
    Vector{Vector{Tuple{Int, Int}}}())

id(r::Region) = r.id
label(r::Region) = id(r)
parents(r::Region) = r.parents
ancestors(r::Region) = r.ancestors
subregions(r::Region) = r.subregions
descendents(r::Region) = subregions(r)
overcounting_number(r::Region) = r.overcounting_number
counting_number(r::Region) = overcounting_number(r)
missing_parent_indices(r::Region) = r.missing_parent_indices
message_numerators(r::Region) = r.message_numerators
message_numerators(r::Region, i::Int) = r.message_numerators[i]
message_denominators(r::Region) = r.message_denominators
message_denominators(r::Region, i::Int) = r.message_denominators[i]

function ==(r1::Region, r2::Region)
# do this entirely in terms of id's at every comparison

end

mutable struct RegionGraph
    regions::Vector{Region}
end

# iterate(R::RegionGraph, state = 1) = state > length(R.regions) ? nothing : (R.regions[state], state + 1)

regions(R::RegionGraph) = R.regions
base_regions(R::RegionGraph) = [r for r in R.regions if isempty(r.parents)]
outer_regions(R::RegionGraph) = base_regions(R)
basic_clusters(R::RegionGraph) = base_regions(R)
leaves(R::RegionGraph) = [r for r in R.regions if isempty(r.subregions)]

function canonical_region_graph(H::CTMatrixTypes)
    num_check, num_var = size(H)
    check_adj_list = [Int[] for _ in 1:num_check]
    for r in 1:num_check
        for c in 1:num_var
            iszero(H[r, c]) || push!(check_adj_list[r], c)
        end
    end
    return region_graph_from_base_nodes(check_adj_list)
end
canonical_region_graph(L::AbstractLDPCCode) = canonical_region_graph(parity_check_matrix(L))

function region_graph_from_base_nodes(regions::Vector{Region})
    isempty(regions) && return regions

    left = 1
    right = length(regions)
    while left < right
        for r1 in left:right - 1
            for r2 in left + 1:right
                if r1 ≠ r2
                    cap = regions[r1].id ∩ regions[r2].id
                    if !isempty(cap) && cap ≠ regions[r1].id && cap ≠ regions[r2].id
                        found = false
                        for r in regions
                            if cap == r.id
                                found = true
                                for par in r.parents
                                    if par ∈ regions[r1].ancestors
                                        r.parents = [r3 for r3 in r.parents if r3 ≠ par]
                                    end

                                    if par ∈ regions[r2].ancestors
                                        r.parents = [r3 for r3 in r.parents if r3 ≠ par]
                                    end
                                end

                                regions[r1] ∉ r.parents && push!(r.parents, regions[r1])
                                regions[r2] ∉ r.parents && push!(r.parents, regions[r2])
                                #### TODO here find missing labels between parents and child


                                
                                r.ancestors = unique!(reduce(vcat, [[r3.ancestors for r3 in r.parents]; r.parents]))
                                # do I want to blank this out?
                                # r.subregions = Vector{Region}()

                                for r3 in r.ancestors
                                    r ∉ r3.subregions && push!(r3.subregions, r)
                                end
                                break
                            end
                        end

                        if !found
                            ancestors = unique!([regions[r1].ancestors; regions[r2].ancestors; [regions[r1], regions[r2]]])

                            c_r = 1
                            for r3 in ancestors
                                c_r -= r3.overcounting_number
                            end

                            push!(regions, Region(cap, [regions[r1], regions[r2]], ancestors, Vector{Region}(), c_r))
                            
                            # can we combine this with the above loop?
                            for r3 in ancestors
                                push!(r3.subregions, regions[end])
                            end

                        end
                    end
                end
            end
        end
        left = right + 1
        right = length(regions)
    end
    return RegionGraph(regions)
end

function region_graph_from_base_nodes(R::Vector{Vector{Int}})
    regions = Vector{Region}()
    for r in R
        push!(regions, Region(r, Vector{Region}(), Vector{Region}(), Vector{Region}(), 1))
    end
    return region_graph_from_base_nodes(regions)
end

function is_valid_region_graph(R::RegionGraph)
    ids = unique!([r.id for r in R.regions])
    for id in ids
        c_r = 0
        for r in R.regions
            if id ∈ r.ids
                c_r += r.overcounting_number
            end
        end
        c_r ≠ 1 && return false
    end

    for r in R.regions
        c_r = r.overcounting_number
        for r2 in r.ancestors
            c_r += r2.overcounting_number
        end
        c_r ≠ 1 && return false
    end

    return true
end

function remove_zero_overcounting_numbers(R::RegionGraph)
    isempty(R.regions) && return R

    regions = copy(R.regions)
    i = 1
    while i ≤ length(regions)
        r = regions[i]
        if iszero(r.overcounting_number)
            remove_flag = true
            if isempty(r.subregions)
                # make sure removing this does not make the graph disconnected
                for r2 in r.parents
                    remove_flag = false
                    for r3 in r.parents
                        if r2 !== r3
                            # check if parents have a common ancestor
                            for anc2 in r2.ancestors
                                for anc3 in r3.ancestors
                                    if anc2 === anc3
                                        remove_flag = true
                                        break
                                    end
                                end
                                remove_flag && break
                            end
                        end
                        remove_flag && break
                    end
                end
            end

            if remove_flag
                for r2 in r.subregions
                    # remove r from all ancestors of subregion
                    r2.ancestors = [r3 for r3 in r2.ancestors if r3 !== r]
                    # if r is a parent of descendent, then remove and set parents of r
                    # as parents of descendent
                    if r ∈ r2.parents
                        r2.parents = hcat([r3 for r3 in r2.parents if r3 !== r], r.parents)
                    end
                end

                # remove r from all subregions of ancestors
                for r2 in r.ancestors
                    r2.subregions = [r3 for r3 in r2.subregions if r3 !== r]
                end

                # remove from R
                regions = [r2 for r2 in regions if r2 !== r]
            else
                i += 1
            end
        else
            i += 1
        end
    end
    return RegionGraph(regions)
end

# TODO: print function?
# TODO: sort by id length?
# Add node option?
# Triangulate/optimize

function remove_diamonds(R::RegionGraph)

end

function remove_generational_skips(R::RegionGraph)
    # look at parents of parents and see if they are also your parents, if so, remove
end

function reduce_memory_footprint(R::RegionGraph)
    for r in R.regions
        # what exactly do I not need anymore once it's all created?
        r.ancestors = Vector{Region}()

    end
end

function message_passing_order(R::RegionGraph)
    # node update order is a BFS up from each leaf
end

# R = CodingTheory.region_graph_from_base_nodes([[1, 2, 4, 5], [2, 3, 5, 6], [4, 5, 7, 8], [5, 6, 8, 9]]);
# for r in R
#     println("id = ", r.id)
#     for p in r.parents
#         println("parent: ", p.id)
#     end
#     for p in r.subregions
#        println("sub: ", p.id)
#     end
#     for p in r.ancestors
#        println("anc: ", p.id)
#     end
# end

#############################
            # GBP
#############################

