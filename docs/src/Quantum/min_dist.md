# [Quantum Minimum Distance](@id quantum-minimum-distance-api)

The minimum distance of a quantum code is the lowest weight of a logical
operator, that is, of an element of the normalizer that is not in the
stabilizer group. Degeneracy is what makes this harder than the classical
problem: many distinct operators represent the same logical action, so the
search is over cosets rather than over codewords.

For a CSS code the problem splits, and `X_minimum_distance` and
`Z_minimum_distance` may be computed independently; the code distance is the
smaller of the two. For a subsystem code there are two different answers and
they must not be conflated. The *bare* distance minimizes over operators that
commute with the whole gauge group, while the *dressed* distance allows
multiplication by gauge operators and so can be smaller. The bare and dressed
functions below are separate for exactly this reason.

As in the classical case, the exact solvers return certified values, while the
probabilistic and heuristic searches return upper bounds witnessed by an
operator they found. The setter functions record externally supplied bounds and
do not verify them. Distance computations are the most expensive operation in
this library; see the
[minimum distance tutorial](@ref minimum-distance-tutorial) for how to choose a
solver and how to supply known bounds to prune the search.

## Bounds and bookkeeping

```@autodocs
Modules = [CodingTheory]
Pages = ["Quantum/min_dist_bounds.jl"]
Private = false
```

## Exact algorithms

```@autodocs
Modules = [CodingTheory]
Pages = ["Quantum/min_dist_exact.jl"]
Private = false
```

## Probabilistic algorithms

```@autodocs
Modules = [CodingTheory]
Pages = ["Quantum/min_dist_probabilistic.jl"]
Private = false
```

## Heuristic algorithms

```@autodocs
Modules = [CodingTheory]
Pages = ["Quantum/min_dist_heuristics.jl"]
Private = false
```
