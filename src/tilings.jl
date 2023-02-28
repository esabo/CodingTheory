# Copyright (c) 2022, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # Hyperbolic
#############################

struct ReflectionGroup
    group::GapObj
    generators::Vector{GapObj}
    orders::Vector{Int}
    dimension::Int
end

"""
    trianglegroup(l::Int, m::Int, n::Int)

Return the (`l`,`m`,`n`) triangle group.
"""
function trianglegroup(l::Int, m::Int, n::Int)
    l >= 0 && m >= 0 && n >= 0 || throw(ArgumentError("Arguments must be non-negative.")) 
    # f = GAP.Globals.FreeGroup(g"a", g"b", g"c")
    # f = GAP.Globals.FreeGroup(GAP.Obj.(["a","b","c"]))
    # f = free_group(["a", "b", "c"]).X
    f = GAP.Globals.FreeGroup(GAP.GapObj(["a", "b", "c"]; recursive=true))
    g = f / GapObj([f.:1^2, f.:2^2, f.:3^2, (f.:1 * f.:2)^l, (f.:1 * f.:3)^m, (f.:2 * f.:3)^n])
    return ReflectionGroup(g, [g.:1, g.:2, g.:3], [l, m, n], 3)
end

"""
    rsgroup(r::Int, s::Int)

Return the Coxeter group corresponding to Schläfli symbol {`r`,`s`}.

Corresponding Coxeter diagram:
```
o---o---o
  r   s
```
"""
rsgroup(r::Int, s::Int) = trianglegroup(r, 2, s)

"""
    tetrahedrongroup(orders::Vector{Int})

Return the tetrahedron group with relations given by `orders`.
"""
function tetrahedrongroup(orders::Vector{Int})
    all(>=(0), orders) || throw(ArgumentError("Arguments must be non-negative."))
    # f = GAP.Globals.FreeGroup(g"a", g"b", g"c", g"d")
    # f = GAP.Globals.FreeGroup(GAP.Obj.(["a","b","c", "d"]))
    # f = free_group(["a", "b", "c", "d"]).X
    f = GAP.Globals.FreeGroup(GAP.GapObj(["a", "b", "c", "d"]; recursive=true))
    g = f / GapObj([f.:1^2, f.:2^2, f.:3^2, f.:4^2,
        (f.:1 * f.:2)^orders[1], (f.:1 * f.:3)^orders[2], (f.:1 * f.:4)^orders[3],
        (f.:2 * f.:3)^orders[4], (f.:2 * f.:4)^orders[5], (f.:3 * f.:4)^orders[6]])
    return ReflectionGroup(g, [g.:1, g.:2, g.:3, g.:4], orders, 4)
end

"""
    qrsgroup(q::Int, r::Int, s::Int)

Return the Coxeter group corresponding to Schläfli symbol {`q`,`r`,`s`}.

Corresponding Coxeter diagram:
```
o---o---o---o
  q   r   s
```
"""
qrsgroup(q::Int, r::Int, s::Int) = tetrahedrongroup([q, 2, 2, r, 2, s])

"""
    startetrahedrongroup(q::Int, r::Int, s::Int)

Return the "star" Coxeter group with higher-order (>2) relations given by `q`,`r`,`s`.

Corresponding Coxeter diagram:
```
      o
     / r
o---o
  q  \\ s
      o
```
"""
startetrahedrongroup(q::Int, r::Int, s::Int) = tetrahedrongroup([q, r, s, 2, 2, 2])

"""
    cycletetrahedrongroup(q::Int, r::Int, s::Int, t::Int)

Return the "cycle" Coxeter group with high-order (>2) relations given by `q`,`r`,`s`,`t`.

Corresponding Coxeter diagram:
```
   q
 o---o
t|   |r
 o---o
   s
```
"""
cycletetrahedrongroup(q::Int, r::Int, s::Int, t::Int) = tetrahedrongroup([q, 2, t, r, 2, s])

"""
    normalsubgroups(g::ReflectionGroup, maxindex::Int)

Return all normal subgroups of `g` with index up to `maxindex`.
"""
function normalsubgroups(g::ReflectionGroup, maxindex::Int)
    gr = GAP.Globals.LowIndexNormalSubgroupsSearchForAll(g.group, maxindex)
    lns = GAP.Globals.List(gr)
    sbgrps = Vector{GapObj}()
    len = GAP.Globals.Length(lns)
    for i = 1:len
        push!(sbgrps, GAP.Globals.Grp(lns[i]))
    end
    return sbgrps
end

"""
    fixedpointfree(subgroup::GapObj, g::ReflectionGroup)

Return `true` if the `subgroup` of `g` is fixed-point free; otherwise `false`.
"""
function fixedpointfree(subgroup::GapObj, g::ReflectionGroup)
    hom = GAP.Globals.NaturalHomomorphismByNormalSubgroup(g.group, subgroup)
    fpf = true
    G = GapObj(())
    for i in 1:g.dimension
        fpf = fpf && GAP.Globals.Image(hom, g.generators[i]) != G
    end
    i = 1
    for pair in combinations(1:g.dimension, 2)
        for j in 1:g.orders[i] - 1
            fpf = fpf && GAP.Globals.Image(hom, (g.generators[pair[1]] * g.generators[pair[2]])^j) != G
        end
        i += 1
    end
    return fpf
end

"""
    orientable(subgroup::GapObj, g::ReflectionGroup)

Return `true' if the `subgroup` of `g` is orientable; otherwise `false`.
"""
function orientable(sbgrps::GapObj, g::ReflectionGroup)
    gens = Vector{GapObj}()
    for pair in combinations(1:g.dimension, 2)
        push!(gens, g.generators[pair[1]] * g.generators[pair[2]])
    end
    gens = GapObj(gens)
    gplus = GAP.Globals.Subgroup(g.group, gens)
    return GAP.Globals.IsSubgroup(gplus, sbgrps)
end

"""
    kcolorable(k::Int, genidx::Vector{GapObj}, translations::Vector{GapObj}, subgroup::GapObj, g::ReflectionGroup)

Return `true` if the group elements corresponding to `genidx` in `g/subgroup` are
`k`-colorable; otherwise `false`.
"""
function kcolorable(k::Int, genidx::Vector{Int}, translations::Vector{GapObj}, subgroup::GapObj, g::ReflectionGroup)
    gens = GAP.Globals.List(GAP.Globals.GeneratorsOfGroup(subgroup))
    GAP.Globals.Append(gens, GapObj(getindex(g.generators, genidx)))
    GAP.Globals.Append(gens, GapObj(translations))
    subgroupT = GAP.Globals.GroupByGenerators(gens)
    return GAP.Globals.Length(GAP.Globals.RightCosets(g.group, subgroupT)) == k
end

"""
    cosetintersection(genidxA::Vector{Int}, genidxB::Vector{Int}, subgroup::GapObj, g::ReflectionGroup)

Return the intersection of the cosets of `g/subgroup` wrt `genidxA` and wrt `genidxB`.

Outputs a sparse matrix with rows indexing the `genidxA` cosets and columns indexing the `genidxB` cosets.
"""
function cosetintersection(genidxA::Vector{Int}, genidxB::Vector{Int}, subgroup::GapObj, g::ReflectionGroup)
    gens = GAP.Globals.List(GAP.Globals.GeneratorsOfGroup(subgroup))
    gensA = deepcopy(gens)
    gensB = deepcopy(gens)
    GAP.Globals.Append(gensA, GapObj(getindex(g.generators, genidxA)))
    GAP.Globals.Append(gensB, GapObj(getindex(g.generators, genidxB)))
    subgroupA = GAP.Globals.Subgroup(g.group, gensA)
    subgroupB = GAP.Globals.Subgroup(g.group, gensB)
    transversalB = GAP.Globals.RightTransversal(g.group, subgroupB)
    intersectionAB = GAP.Globals.Intersection(subgroupA, subgroupB)
    i = 1
    I = Vector{Int}()
    J = Vector{Int}()
    for cosetA in GAP.Globals.RightCosets(g.group, subgroupA)
        repA = GAP.Globals.Representative(cosetA)
        for cosetAB in GAP.Globals.RightCosets(subgroupA, intersectionAB)
            repAB = GAP.Globals.Representative(cosetAB)
            push!(I, i)
            push!(J, GAP.Globals.PositionCanonical(transversalB, repAB * repA))
        end
        i += 1
    end
    return sparse(I, J, ones(Int, size(I)))
end
