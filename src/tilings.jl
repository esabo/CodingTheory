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

__gapobj(G::ReflectionGroup) = G.group
gens(G::ReflectionGroup) = G.generators
orders(G::ReflectionGroup) = G.orders
dim(G::ReflectionGroup) = G.dimension # length(gens(G)) ???

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
    normalsubgroups(G::ReflectionGroup, maxindex::Int)

Return all normal subgroups of `G` with index up to `maxindex`.
"""
function normalsubgroups(G::ReflectionGroup, maxindex::Integer)
    lins_search = GAP.Globals.LowIndexNormalSubgroupsSearchForAll(__gapobj(G), maxindex)
    sbgrps = GapObj[H for H in GAP.Globals.List(lins_search)]
    return sbgrps
end 

"""
    fixedpointfree(subgroup::GapObj, G::ReflectionGroup)

Return `true` if the `subgroup` of `G` is fixed-point free; otherwise `false`.
"""
function fixedpointfree(subgroup::GapObj, G::ReflectionGroup)
    hom = GAP.Globals.NaturalHomomorphismByNormalSubgroup(__gapobj(G), subgroup)
    h(g) = GAP.Globals.Image(hom, g)
    fixedpoint = any(isone ∘ h, gens(G))
    fixedpoint && return false
    
    for g1 in gens(G)
        for (g2, ord) in zip(gens(G), orders(G))
            fixedpoint = any(isone ∘ h, (g1*g2^j for j in 1:ord - 1))
            fixedpoint && return false
        end
    end
    return true
end 

"""
    orientable(subgroup::GapObj, G::ReflectionGroup)

Return `true' if the `subgroup` of `G` is orientable; otherwise `false`.
"""
function orientable(sbgrps::GapObj, G::ReflectionGroup)
    S = GapObj[a*b for (a,b) in combinations(gens(G), 2)]
    G⁺ = GAP.Globals.Subgroup(__gapobj(G), GapObj(S))
    return GAP.Globals.IsSubgroup(G⁺, sbgrps)
end

"""
    kcolorable(k, genidx, translations, subgroup::GapObj, g::ReflectionGroup)

Return `true` if the group elements corresponding to `genidx` in `G/subgroup` are
`k`-colorable; otherwise `false`.
"""
function kcolorable(
    k::Integer,
    genidx::AbstractVector{<:Integer},
    translations::AbstractVector{<:GapObj},
    subgroup::GapObj,
    G::ReflectionGroup,
)   
    Tgens = GAP.Globals.List(GAP.Globals.GeneratorsOfGroup(subgroup))
    GAP.Globals.Append(Tgens, GapObj(gens(G)[genidx]))
    GAP.Globals.Append(Tgens, GapObj(translations))
    subgroupT = GAP.Globals.GroupByGenerators(Tgens)
    return GAP.Globals.Index(__gapobj(G), subgroupT) ≤ k
end 

"""
    cosetintersection(genidxA, genidxB, subgroup::GapObj, G::ReflectionGroup)

Return the intersection of the cosets of `G/subgroup` wrt `genidxA` and wrt `genidxB`.

Outputs a sparse matrix with rows indexing the `genidxA` cosets and columns indexing the `genidxB` cosets.
"""
function cosetintersection(
    genidxA::AbstractVector{<:Integer},
    genidxB::AbstractVector{<:Integer},
    subgroup::GapObj,
    G::ReflectionGroup,
)   
    A, B = let S = GAP.Globals.List(GAP.Globals.GeneratorsOfGroup(subgroup))
        map((genidxA, genidxB)) do idx
            S = deepcopy(S)
            GAP.Globals.Append(S, GapObj(gens(G)[idx]))
            GAP.Globals.Subgroup(__gapobj(G), S)
        end
    end
    
    transversalB = GAP.Globals.RightTransversal(__gapobj(G), B)
    AB = GAP.Globals.Intersection(A, B)
    I = Vector{Int}()
    J = Vector{Int}()
    for (i, Ag) in enumerate(GAP.Globals.RightCosets(__gapobj(G), A))
        g = GAP.Globals.Representative(Ag)
        for ABh in GAP.Globals.RightCosets(A, AB)
            h = GAP.Globals.Representative(ABh)
            push!(I, i)
            push!(J, GAP.Globals.PositionCanonical(transversalB, h * g))
        end
    end
    return sparse(I, J, ones(Int, length(I)))
end

