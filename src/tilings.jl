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

function Base.getindex(cox_mat::CoxeterMatrix, i::Int, j::Int)
    small_idx = min(i, j)
    big_idx = max(i, j)
    n = size(cox_mat)[1]
    linear_idx = (n * (n - 1) - (n - small_idx) * (n - small_idx + 1)) / 2 + big_idx
    cox_mat.vec[convert(Int, linear_idx)]
end

struct ReflectionGroup
    group::GapObj
    generators::GapObj
    cox_mat::CoxeterMatrix
end

group(G::ReflectionGroup) = G.group
gens(G::ReflectionGroup) = G.generators
cox_mat(G::ReflectionGroup) = G.cox_mat
dim(G::ReflectionGroup) = length(gens(G))

"""
    simplex_group(cox_mat::CoxeterMatrix)

Return the reflection group with fundamental simplex specified by `cox_mat`.
"""
function simplex_group(cox_mat::CoxeterMatrix)
    n = size(cox_mat)[1]
    f = Globals.FreeGroup(n)
    gen = Globals.GeneratorsOfGroup(f)
    relations = []
    # doing it the Oscar way causes severe type problems passing into LINS
    # significantly nicer code just ignoring their wrapper
    # f = free_group(n)
    # gen = gens(f)
    # relations = elem_type(f)[]
    for (i, g1) in enumerate(gen)
        for (j, g2) in enumerate(gen)
            if j >= i append!(relations, [(g1 * g2)^cox_mat[i, j]]) end
        end
    end
    g = f / GapObj(relations)
    return ReflectionGroup(g, Globals.GeneratorsOfGroup(g), cox_mat)
    # Oscar way
    # g, _ = quo(f, relations)
    # return ReflectionGroup(g, gens(g), cox_mat)
end

"""
    triangle_group(l::Int, m::Int, n::Int)

Return the (`l`, `m`, `n`) triangle group.
"""
triangle_group(l::Int, m::Int, n::Int) = simplex_group(CoxeterMatrix(3, [1, l, m, 1, n, 1]))

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
    f = Globals.FreeGroup(GapObj(["a", "b", "c", "d"]; recursive=true))
    g = f / GapObj([f.:1^2, f.:2^2, f.:3^2, f.:4^2,
        (f.:1 * f.:2)^orders[1], (f.:1 * f.:3)^orders[2], (f.:1 * f.:4)^orders[3],
        (f.:2 * f.:3)^orders[4], (f.:2 * f.:4)^orders[5], (f.:3 * f.:4)^orders[6]])
    # BUG: call here is off, the dimension need not be saved but orders is not a CoxeterMatrix
    # return ReflectionGroup(g, [g.:1, g.:2, g.:3, g.:4], orders, 4)
    return ReflectionGroup(g, [g.:1, g.:2, g.:3, g.:4], CoxeterMatrix(4, orders))
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

"""
    normal_subgroups(g::ReflectionGroup, max_index::Int)

Return all normal subgroups of `g` with index up to `max_index`.
"""
function normal_subgroups(G::ReflectionGroup, max_index::Integer)
    lins_search = Globals.LowIndexNormalSubgroupsSearchForAll(group(G), max_index)
    sbgrps = GapObj[Globals.Grp(H) for H in Globals.List(lins_search)]
    return sbgrps
end 

"""
    is_fixed_point_free(subgroup::GapObj, g::ReflectionGroup)

Return `true` if the `subgroup` of `g` is fixed-point free; otherwise `false`.
"""
function is_fixed_point_free(subgroup::GapObj, G::ReflectionGroup)
    hom = Globals.NaturalHomomorphismByNormalSubgroup(group(G), subgroup)
    h(g) = Globals.Image(hom, g)
    fixed_point = any(isone ∘ h, gens(G))
    fixed_point && return false
    
    for (i, g1) in enumerate(gens(G))
        for (j, g2) in enumerate(gens(G))
            if j > i
                k = cox_mat(G)[i, j]
                fixed_point = any(isone ∘ h, ((g1 * g2)^m for m in 1:k - 1))
                fixed_point && return false
            end
        end
    end
    return true
end

"""
    is_orientable(subgroup::GapObj, F::ReflectionGroup)

Return `true` if the `subgroup` of `F` is is_orientable; otherwise `false`.
"""
function is_orientable(subgroup::GapObj, G::ReflectionGroup)
    S = GapObj[a * b for (a, b) in combinations(gens(G), 2)]
    G⁺ = Globals.Subgroup(group(G), GapObj(S))
    return Globals.IsSubgroup(G⁺, subgroup)
end

"""
    is_k_colorable(k::Int, gen_idx::Vector{GapObj}, translations::Vector{GapObj}, subgroup::GapObj, g::ReflectionGroup)

Return `true` if the group elements corresponding to `gen_idx` in `g/subgroup` are
`k`-colorable; otherwise `false`.
"""
function is_k_colorable(k::Int, gen_idx::AbstractVector{<: Int},
    translations::AbstractVector{<:GapObj}, subgroup::GapObj, G::ReflectionGroup)   

    T_gens = Globals.List(Globals.GeneratorsOfGroup(subgroup))
    Globals.Append(T_gens, gens(G)[gen_idx])
    Globals.Append(T_gens, GapObj(translations))
    subgroup_T = Globals.GroupByGenerators(T_gens)
    return Globals.Index(group(G), subgroup_T) ≤ k
end 

"""
    coset_intersection(gen_idx_A::Vector{Int}, gen_idx_B::Vector{Int}, subgroup::GapObj, g::ReflectionGroup)

Return the intersection of the cosets of `g/subgroup` wrt `gen_idx_A` and wrt `gen_idx_B`.

# Notes
* This outputs a sparse matrix with rows indexing the `gen_idx_A` cosets and columns indexing the `gen_idx_B` cosets.
"""
function coset_intersection(gen_idx_A::Vector{Int}, gen_idx_B::Vector{Int}, subgroup::GapObj, G::ReflectionGroup, transversal::Union{GapObj, Nothing} = nothing)

    A, B = let S = Globals.List(Globals.GeneratorsOfGroup(subgroup))
        map((gen_idx_A, gen_idx_B)) do idx
            S_copy = deepcopy(S)
            Globals.Append(S_copy, gens(G)[idx])
            Globals.Subgroup(group(G), S_copy)
        end
    end
    if isnothing(transversal)
        transversal_B = Globals.RightTransversal(group(G), B)
    else
        transversal_B = transversal
    end
    AB = Globals.Intersection(A, B)
    I = Vector{Int}()
    J = Vector{Int}()
    for (i, Ag) in enumerate(Globals.RightCosets(group(G), A))
        g = Globals.Representative(Ag)
        for ABh in Globals.RightCosets(A, AB)
            h = Globals.Representative(ABh)
            push!(I, i)
            push!(J, Globals.PositionCanonical(transversal_B, h * g))
        end
    end
    return sparse(I, J, ones(Int, length(I)))
end
