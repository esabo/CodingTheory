# Copyright (c) 2022, 2023 Michael Vasmer, Eric Sabo
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
    triangle_group(l::Int, m::Int, n::Int)

Return the (`l`, `m`, `n`) triangle group.
"""
function triangle_group(l::Int, m::Int, n::Int)
    l >= 0 && m >= 0 && n >= 0 || throw(ArgumentError("Arguments must be non-negative.")) 
    # f = GAP.Globals.FreeGroup(g"a", g"b", g"c")
    # f = GAP.Globals.FreeGroup(GAP.Obj.(["a","b","c"]))
    # f = free_group(["a", "b", "c"]).X
    f = GAP.Globals.FreeGroup(GAP.GapObj(["a", "b", "c"]; recursive=true))
    g = f / GapObj([f.:1^2, f.:2^2, f.:3^2, (f.:1 * f.:2)^l, (f.:1 * f.:3)^m, (f.:2 * f.:3)^n])
    return ReflectionGroup(g, [g.:1, g.:2, g.:3], [l, m, n], 3)
end

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
    qr_s_group(q::Int, r::Int, s::Int)

Return the Coxeter group corresponding to Schläfli symbol {`q`, `r`, `s`}.

# Corresponding Coxeter diagram:
```
o---o---o---o
  q   r   s
```
"""
qr_s_group(q::Int, r::Int, s::Int) = tetrahedron_group([q, 2, 2, r, 2, s])

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
star_tetrahedron_group(q::Int, r::Int, s::Int) = tetrahedron_group([q, r, s, 2, 2, 2])

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
cycle_tetrahedron_group(q::Int, r::Int, s::Int, t::Int) = tetrahedron_group([q, 2, r, s, 2, t])

"""
    normal_subgroups(g::ReflectionGroup, max_index::Int)

Return all normal subgroups of `g` with index up to `max_index`.
"""
function normal_subgroups(g::ReflectionGroup, max_index::Int)
    gr = GAP.Globals.LowIndexNormalSubgroupsSearchForAll(g.group, max_index)
    lns = GAP.Globals.List(gr)
    subgroups = Vector{GapObj}()
    len = GAP.Globals.Length(lns)
    for i in 1:len
        push!(subgroups, GAP.Globals.Grp(lns[i]))
    end
    return subgroups
end

"""
    is_fixed_point_free(subgroup::GapObj, g::ReflectionGroup)

Return `true` if the `subgroup` of `g` is fixed-point free; otherwise `false`.
"""
function is_fixed_point_free(subgroup::GapObj, g::ReflectionGroup)
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
    is_orientable(subgroup::GapObj, g::ReflectionGroup)

Return `true` if the `subgroup` of `g` is is_orientable; otherwise `false`.
"""
function is_orientable(subgroups::GapObj, g::ReflectionGroup)
    gens = Vector{GapObj}()
    for pair in combinations(1:g.dimension, 2)
        push!(gens, g.generators[pair[1]] * g.generators[pair[2]])
    end
    gens = GapObj(gens)
    g_plus = GAP.Globals.Subgroup(g.group, gens)
    return GAP.Globals.IsSubgroup(g_plus, subgroups)
end

"""
    is_k_colorable(k::Int, gen_idx::Vector{GapObj}, translations::Vector{GapObj}, subgroup::GapObj, g::ReflectionGroup)

Return `true` if the group elements corresponding to `gen_idx` in `g/subgroup` are
`k`-colorable; otherwise `false`.
"""
function is_k_colorable(k::Int, gen_idx::Vector{Int}, translations::Vector{GapObj}, subgroup::GapObj, g::ReflectionGroup)
    gens = GAP.Globals.List(GAP.Globals.GeneratorsOfGroup(subgroup))
    GAP.Globals.Append(gens, GapObj(getindex(g.generators, gen_idx)))
    GAP.Globals.Append(gens, GapObj(translations))
    subgroup_T = GAP.Globals.GroupByGenerators(gens)
    return GAP.Globals.Length(GAP.Globals.RightCosets(g.group, subgroup_T)) == k
end

"""
    coset_intersection(gen_idx_A::Vector{Int}, gen_idx_B::Vector{Int}, subgroup::GapObj, g::ReflectionGroup)

Return the intersection of the cosets of `g/subgroup` wrt `gen_idx_A` and wrt `gen_idx_B`.

# Notes
* This outputs a sparse matrix with rows indexing the `gen_idx_A` cosets and columns indexing the `gen_idx_B` cosets.
"""
function coset_intersection(gen_idx_A::Vector{Int}, gen_idx_B::Vector{Int}, subgroup::GapObj, g::ReflectionGroup)
    gens = GAP.Globals.List(GAP.Globals.GeneratorsOfGroup(subgroup))
    gens_A = deepcopy(gens)
    gens_B = deepcopy(gens)
    GAP.Globals.Append(gens_A, GapObj(getindex(g.generators, gen_idx_A)))
    GAP.Globals.Append(gens_B, GapObj(getindex(g.generators, gen_idx_B)))
    subgroup_A = GAP.Globals.Subgroup(g.group, gens_A)
    subgroup_B = GAP.Globals.Subgroup(g.group, gens_B)
    transversal_B = GAP.Globals.RightTransversal(g.group, subgroup_B)
    intersection_AB = GAP.Globals.Intersection(subgroup_A, subgroup_B)
    i = 1
    I = Vector{Int}()
    J = Vector{Int}()
    for coset_A in GAP.Globals.RightCosets(g.group, subgroup_A)
        rep_A = GAP.Globals.Representative(coset_A)
        for coset_AB in GAP.Globals.RightCosets(subgroup_A, intersection_AB)
            rep_AB = GAP.Globals.Representative(coset_AB)
            push!(I, i)
            push!(J, GAP.Globals.PositionCanonical(transversal_B, rep_AB * rep_A))
        end
        i += 1
    end
    return sparse(I, J, ones(Int, size(I)))
end
