# Copyright (c) 2022, 2023 Michael Vasmer, Eric Sabo
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
    simplex_group(coxmat::CoxeterMatrix)

Return the reflection group with fundamental simplex specified by `coxmat`.
"""
function simplex_group(coxmat::CoxeterMatrix)
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
    triangle_group(l::Int, m::Int, n::Int)

Return the (`l`, `m`, `n`) triangle group.
"""
triangle_group(l::Int, m::Int, n::Int) = simplex_group(CoxeterMatrix(3, [1, l, m, 1, n, 1]))
# function triangle_group(l::Int, m::Int, n::Int)
#     l >= 0 && m >= 0 && n >= 0 || throw(ArgumentError("Arguments must be non-negative.")) 
#     # f = GAP.Globals.FreeGroup(g"a", g"b", g"c")
#     # f = GAP.Globals.FreeGroup(GAP.Obj.(["a","b","c"]))
#     # f = free_group(["a", "b", "c"]).X
#     f = GAP.Globals.FreeGroup(GAP.GapObj(["a", "b", "c"]; recursive=true))
#     g = f / GapObj([f.:1^2, f.:2^2, f.:3^2, (f.:1 * f.:2)^l, (f.:1 * f.:3)^m, (f.:2 * f.:3)^n])
#     return ReflectionGroup(g, [g.:1, g.:2, g.:3], [l, m, n], 3)
# end

"""
    r_s_group(r::Int, s::Int)

Return the Coxeter group corresponding to Schläfli symbol `{r, s}`.

# Corresponding Coxeter diagram:
```
o---o---o
  r   s
```
"""
r_s_group(r::Int, s::Int) = triangle_group(r, 2, s)

"""
    tetrahedron_group(orders::Vector{Int})

Return the tetrahedron group with relations given by `orders`.
"""
function tetrahedron_group(orders::Vector{Int})
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
    q_r_s_group(q::Int, r::Int, s::Int)

Return the Coxeter group corresponding to Schläfli symbol {`q`, `r`, `s`}.

# Corresponding Coxeter diagram:
```
o---o---o---o
  q   r   s
```
"""
q_r_s_group(q::Int, r::Int, s::Int) = simplex_group(CoxeterMatrix(4, [1, q, 2, 2, 1, r, 2, 1, s, 1]))
# qr_s_group(q::Int, r::Int, s::Int) = tetrahedron_group([q, 2, 2, r, 2, s])

"""
    star_tetrahedron_group(q::Int, r::Int, s::Int)

Return the "star" Coxeter group with higher-order (>2) relations given by `q`, `r`, and `s`.

# Corresponding Coxeter diagram:
```
      o
     / r
o---o
  q  \\ s
      o
```
"""
star_tetrahedron_group(q::Int, r::Int, s::Int) = simplex_group(CoxeterMatrix(4, [1, q, r, s, 1, 2, 2, 1, 2, 1]))
# star_tetrahedron_group(q::Int, r::Int, s::Int) = tetrahedron_group([q, r, s, 2, 2, 2])

"""
    cycle_tetrahedron_group(q::Int, r::Int, s::Int, t::Int)

Return the "cycle" Coxeter group with high-order (>2) relations given by `q`, `r`, `s`, and `t`.

# Corresponding Coxeter diagram:
```
   q
 o---o
t|   |r
 o---o
   s
```
"""
cycle_tetrahedron_group(q::Int, r::Int, s::Int, t::Int) = simplex_group(CoxeterMatrix(4, [1, q, 2, t, 1, r, 2, 1, s, 1]))
# cycle_tetrahedron_group(q::Int, r::Int, s::Int, t::Int) = tetrahedron_group([q, 2, r, s, 2, t])

"""
    normal_subgroups(g::ReflectionGroup, max_index::Int)

Return all normal subgroups of `g` with index up to `max_index`.
"""
function normal_subgroups(G::ReflectionGroup, max_index::Integer)
    lins_search = GAP.Globals.LowIndexNormalSubgroupsSearchForAll(__gapobj(G), max_index)
    sbgrps = GapObj[GAP.Globals.Grp(H) for H in GAP.Globals.List(lins_search)]
    return sbgrps
end 
# function normal_subgroups(g::ReflectionGroup, max_index::Int)
#     gr = GAP.Globals.LowIndexNormalSubgroupsSearchForAll(g.group, max_index)
#     lns = GAP.Globals.List(gr)
#     subgroups = Vector{GapObj}()
#     len = GAP.Globals.Length(lns)
#     for i in 1:len
#         push!(subgroups, GAP.Globals.Grp(lns[i]))
#     end
#     return subgroups
# end

"""
    is_fixed_point_free(subgroup::GapObj, g::ReflectionGroup)

Return `true` if the `subgroup` of `g` is fixed-point free; otherwise `false`.
"""
function is_fixed_point_free(subgroup::GapObj, G::ReflectionGroup)
    hom = GAP.Globals.NaturalHomomorphismByNormalSubgroup(__gapobj(G), subgroup)
    h(g) = GAP.Globals.Image(hom, g)
    fixed_point = any(isone ∘ h, gens(G))
    fixed_point && return false
    
    for (i, g1) in enumerate(gens(G))
        for (j, g2) in enumerate(gens(G))
            if j > i
                k = coxmat(G)[i,j]
                fixed_point = any(isone ∘ h, ((g1 * g2)^m for m in 1:k - 1))
                fixed_point && return false
            end
        end
    end
    return true
end
# function is_fixed_point_free(subgroup::GapObj, g::ReflectionGroup)
#     hom = GAP.Globals.NaturalHomomorphismByNormalSubgroup(g.group, subgroup)
#     fpf = true
#     G = GapObj(())
#     for i in 1:g.dimension
#         fpf = fpf && GAP.Globals.Image(hom, g.generators[i]) != G
#     end
#     i = 1
#     for pair in combinations(1:g.dimension, 2)
#         for j in 1:g.orders[i] - 1
#             fpf = fpf && GAP.Globals.Image(hom, (g.generators[pair[1]] * g.generators[pair[2]])^j) != G
#         end
#     end
#     return true
# end

"""
    is_orientable(subgroup::GapObj, F::ReflectionGroup)

Return `true` if the `subgroup` of `F` is is_orientable; otherwise `false`.
"""
function is_orientable(subgroup::GapObj, G::ReflectionGroup)
    S = GapObj[a*b for (a,b) in combinations(gens(G), 2)]
    G⁺ = GAP.Globals.Subgroup(__gapobj(G), GapObj(S))
    return GAP.Globals.IsSubgroup(G⁺, subgroup)
end
# function is_orientable(subgroups::GapObj, g::ReflectionGroup)
#     gens = Vector{GapObj}()
#     for pair in combinations(1:g.dimension, 2)
#         push!(gens, g.generators[pair[1]] * g.generators[pair[2]])
#     end
#     gens = GapObj(gens)
#     g_plus = GAP.Globals.Subgroup(g.group, gens)
#     return GAP.Globals.IsSubgroup(g_plus, subgroups)
# end

"""
    is_k_colorable(k::Int, gen_idx::Vector{GapObj}, translations::Vector{GapObj}, subgroup::GapObj, g::ReflectionGroup)

Return `true` if the group elements corresponding to `gen_idx` in `g/subgroup` are
`k`-colorable; otherwise `false`.
"""
function is_k_colorable(k::Int, gen_idx::AbstractVector{<: Int},
    translations::AbstractVector{<:GapObj}, subgroup::GapObj, G::ReflectionGroup)   

    T_gens = GAP.Globals.List(GAP.Globals.GeneratorsOfGroup(subgroup))
    GAP.Globals.Append(T_gens, gens(G)[gen_idx])
    GAP.Globals.Append(T_gens, GapObj(translations))
    subgroup_T = GAP.Globals.GroupByGenerators(T_gens)
    return GAP.Globals.Index(__gapobj(G), subgroup_T) ≤ k
end 
# function is_k_colorable(k::Int, gen_idx::Vector{Int}, translations::Vector{GapObj}, subgroup::GapObj, g::ReflectionGroup)
#     gens = GAP.Globals.List(GAP.Globals.GeneratorsOfGroup(subgroup))
#     GAP.Globals.Append(gens, GapObj(getindex(g.generators, gen_idx)))
#     GAP.Globals.Append(gens, GapObj(translations))
#     subgroup_T = GAP.Globals.GroupByGenerators(gens)
#     return GAP.Globals.Length(GAP.Globals.RightCosets(g.group, subgroup_T)) == k
# end

"""
    coset_intersection(gen_idx_A::Vector{Int}, gen_idx_B::Vector{Int}, subgroup::GapObj, g::ReflectionGroup)

Return the intersection of the cosets of `g/subgroup` wrt `gen_idx_A` and wrt `gen_idx_B`.

# Notes
* This outputs a sparse matrix with rows indexing the `gen_idx_A` cosets and columns indexing the `gen_idx_B` cosets.
"""
function coset_intersection(gen_idx_A::Vector{Int}, gen_idx_B::Vector{Int}, subgroup::GapObj, G::ReflectionGroup, transversal::Union{GapObj, Nothing} = nothing)

    A, B = let S = GAP.Globals.List(GAP.Globals.GeneratorsOfGroup(subgroup))
        map((gen_idx_A, gen_idx_B)) do idx
            S_copy = deepcopy(S)
            GAP.Globals.Append(S_copy, gens(G)[idx])
            GAP.Globals.Subgroup(__gapobj(G), S_copy)
        end
    end
    if isnothing(transversal)
        transversal_B = GAP.Globals.RightTransversal(__gapobj(G), B)
    else
        transversal_B = transversal
    end
    AB = GAP.Globals.Intersection(A, B)
    I = Vector{Int}()
    J = Vector{Int}()
    for (i, Ag) in enumerate(GAP.Globals.RightCosets(__gapobj(G), A))
        g = GAP.Globals.Representative(Ag)
        for ABh in GAP.Globals.RightCosets(A, AB)
            h = GAP.Globals.Representative(ABh)
            push!(I, i)
            push!(J, GAP.Globals.PositionCanonical(transversal_B, h * g))
        end
    end
    return sparse(I, J, ones(Int, length(I)))
end

# function coset_intersection(gen_idx_A::Vector{Int}, gen_idx_B::Vector{Int}, subgroup::GapObj, g::ReflectionGroup)
#     gens = GAP.Globals.List(GAP.Globals.GeneratorsOfGroup(subgroup))
#     gens_A = deepcopy(gens)
#     gens_B = deepcopy(gens)
#     GAP.Globals.Append(gens_A, GapObj(getindex(g.generators, gen_idx_A)))
#     GAP.Globals.Append(gens_B, GapObj(getindex(g.generators, gen_idx_B)))
#     subgroup_A = GAP.Globals.Subgroup(g.group, gens_A)
#     subgroup_B = GAP.Globals.Subgroup(g.group, gens_B)
#     transversal_B = GAP.Globals.RightTransversal(g.group, subgroup_B)
#     intersection_AB = GAP.Globals.Intersection(subgroup_A, subgroup_B)
#     i = 1
#     I = Vector{Int}()
#     J = Vector{Int}()
#     for coset_A in GAP.Globals.RightCosets(g.group, subgroup_A)
#         rep_A = GAP.Globals.Representative(coset_A)
#         for coset_AB in GAP.Globals.RightCosets(subgroup_A, intersection_AB)
#             rep_AB = GAP.Globals.Representative(coset_AB)
#             push!(I, i)
#             push!(J, GAP.Globals.PositionCanonical(transversal_B, rep_AB * rep_A))
#         end
#     end
#     return sparse(I, J, ones(Int, length(I)))
# end
