# Copyright (c) 2022, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # Hyperbolic
#############################

struct CoxeterMatrix <: AbstractMatrix{Int}
    n::Int
    vec::Vector{Int}
end
Base.size(m::CoxeterMatrix) = (m.n, m.n)

function Base.getindex(coxmat::CoxeterMatrix, i::Int, j::Int)
    smallidx = min(i, j)
    bigidx = max(i, j)
    n = size(coxmat)[1]
    linearidx = (n * (n - 1) - (n - smallidx) * (n - smallidx + 1)) / 2 + bigidx
    coxmat.vec[convert(Int, linearidx)]
end

struct ReflectionGroup
    group::GapObj
    generators::GapObj
    coxmat::CoxeterMatrix
end

__gapobj(G::ReflectionGroup) = G.group
gens(G::ReflectionGroup) = G.generators
coxmat(G::ReflectionGroup) = G.coxmat
dim(G::ReflectionGroup) = length(gens(G))

"""
    simplexgroup(coxmat::CoxeterMatrix)

Return the reflection group with fundamental simplex specified by `coxmat`.
"""
function simplexgroup(coxmat::CoxeterMatrix)
    n = size(coxmat)[1]
    f = GAP.Globals.FreeGroup(n)
    gens = GAP.Globals.GeneratorsOfGroup(f)
    relations = []
    for (i, g1) in enumerate(gens)
        for (j, g2) in enumerate(gens)
            if j >= i append!(relations, [(g1*g2)^coxmat[i,j]]) end
        end
    end
    g = f / GapObj(relations)
    return ReflectionGroup(g, GAP.Globals.GeneratorsOfGroup(g), coxmat)
end


"""
    trianglegroup(l::Int, m::Int, n::Int)

Return the (`l`,`m`,`n`) triangle group.
"""
trianglegroup(l::Int, m::Int, n::Int) = simplexgroup(CoxeterMatrix(3, [1, l, m, 1, n, 1]))

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
    qrsgroup(q::Int, r::Int, s::Int)

Return the Coxeter group corresponding to Schläfli symbol {`q`,`r`,`s`}.

Corresponding Coxeter diagram:
```
o---o---o---o
  q   r   s
```
"""
qrsgroup(q::Int, r::Int, s::Int) = simplexgroup(CoxeterMatrix(4, [1, q, 2, 2, 1, r, 2, 1, s, 1]))

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
startetrahedrongroup(q::Int, r::Int, s::Int) = simplexgroup(CoxeterMatrix(4, [1, q, r, s, 1, 2, 2, 1, 2, 1]))

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
cycletetrahedrongroup(q::Int, r::Int, s::Int, t::Int) = simplexgroup(CoxeterMatrix(4, [1, q, 2, t, 1, r, 2, 1, s, 1]))

"""
    normalsubgroups(G::ReflectionGroup, maxindex::Int)

Return all normal subgroups of `G` with index up to `maxindex`.
"""
function normalsubgroups(G::ReflectionGroup, maxindex::Integer)
    lins_search = GAP.Globals.LowIndexNormalSubgroupsSearchForAll(__gapobj(G), maxindex)
    sbgrps = GapObj[GAP.Globals.Grp(H) for H in GAP.Globals.List(lins_search)]
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
    
    for (i, g1) in enumerate(gens(G))
        for (j, g2) in enumerate(gens(G))
            if j > i
                k = coxmat(G)[i,j]
                fixedpoint = any(isone ∘ h, ((g1*g2)^m for m in 1:k - 1))
                fixedpoint && return false
            end
        end
    end
    return true
end 

"""
    orientable(subgroup::GapObj, G::ReflectionGroup)

Return `true' if the `subgroup` of `G` is orientable; otherwise `false`.
"""
function orientable(sbgrp::GapObj, G::ReflectionGroup)
    S = GapObj[a*b for (a,b) in combinations(gens(G), 2)]
    G⁺ = GAP.Globals.Subgroup(__gapobj(G), GapObj(S))
    return GAP.Globals.IsSubgroup(G⁺, sbgrp)
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
    GAP.Globals.Append(Tgens, gens(G)[genidx])
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
    transversal::Union{GapObj,Nothing}=nothing,
)   
    A, B = let S = GAP.Globals.List(GAP.Globals.GeneratorsOfGroup(subgroup))
        map((genidxA, genidxB)) do idx
            Scopy = deepcopy(S)
            GAP.Globals.Append(Scopy, gens(G)[idx])
            GAP.Globals.Subgroup(__gapobj(G), Scopy)
        end
    end
    if isnothing(transversal)
        transversalB = GAP.Globals.RightTransversal(__gapobj(G), B)
    else
        transversalB = transversal
    end
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

