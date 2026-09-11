# Expansion and local testability

CodingTheory provides exact finite-size computations for binary Tanner-graph
expansion and reduced syndrome profiles, together with spectral graph
estimates. The exact routines are exponential in the requested subset or
error weight and are intended for small instances and bounded-weight checks.

## Tanner-graph expansion

For a binary matrix `H`, columns are left vertices and rows are right
vertices. The exact predicate

```julia
is_expander(H, γ, A)
```

checks

```math
|N(S)| \geq A |S|
```

for every nonempty column set of size at most
``\lfloor \gamma n\rfloor``. Use `expansion_witness` to obtain a failing
subset and `bipartite_expansion_profile` to compute the exact minimum at each
subset size. `estimated_bipartite_vertex_expansion` uses a Fiedler sweep and
returns only an upper bound on the exact minimum.

`is_left_right_expander` checks both orientations. Code-object methods accept
`:X`, `:Z`, or `:both` for CSS codes. Non-CSS codes support `:both`, where an
X or Z entry connects the stabilizer to that qubit.

## Reduced syndrome profiles

`confinement_profile(H; max_error_weight=t)` performs breadth-first search in
the syndrome space. Consequently, the reported error weight is
``\operatorname{dist}(e,\ker H)``, not the weight of an arbitrary error
representative.

`deterministic_QLTC_soundness` computes

```math
\min \frac{|He|}{\operatorname{dist}(e,\ker H)}
```

over syndrome cosets with a leader of weight at most `max_error_weight`.
`evaluate_single_shot_soundness` restricts the resulting profile by syndrome
weight.

## Spectral estimates

`algebraic_connectivity`, `fiedler_vector`, `estimated_edge_expansion`, and
`estimated_vertex_expansion` operate on ordinary `Graphs.jl` graphs.
`edge_expansion_bounds` returns the combinatorial Cheeger interval. These
graph quantities should not be confused with one-sided qubit-to-check
expansion; use the matrix routines above for that purpose.

## Restricted cosystolic expansion

`cosystolic_expansion(boundary, incoming_boundary; max_weight=t)` computes

```math
\frac{|\partial_k x|}
{\operatorname{dist}(x,\operatorname{im}\partial_{k+1})}
```

exactly for binary vectors of weight at most `t`. It verifies
``\partial_k\partial_{k+1}=0`` and enumerates the incoming image, so
`max_boundary_rank` limits the permitted image rank.

Group-ring and module-relative expansion are not currently exposed. Those
require a canonical group-algebra lifting API and block-metric coset solver;
treating physical Hamming distance as block distance would not be correct.

## API

```@docs
normalized_laplacian_matrix
algebraic_connectivity
normalized_spectral_gap
nontrivial_adjacency_spectral_radius
fiedler_vector
is_topologically_connected
estimated_edge_expansion
estimated_vertex_expansion
edge_expansion_bounds
expansion_witness
is_expander
bipartite_expansion_profile
is_bipartite_expander
is_left_right_expander
estimated_bipartite_vertex_expansion
confinement_profile
deterministic_QLTC_soundness
verify_QLTC_soundness
verify_confinement
evaluate_single_shot_soundness
evaluate_confinement
sipser_spielman_guarantees
cosystolic_expansion
```
