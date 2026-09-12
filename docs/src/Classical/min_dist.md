# [Minimum Distance](@id classical-minimum-distance-api)

Computing the minimum distance of a linear code is NP-hard, so this library
exposes three distinct kinds of routine and it is important not to confuse
them. See the [minimum distance tutorial](@ref minimum-distance-tutorial) for
guidance on choosing among them.

**Exact algorithms** return the true minimum distance and cache it on the code.
These are enumeration-based methods built on bit-packed Gray-code sweeps with
Brouwer-Zimmermann style pruning across disjoint information sets, and they are
the only functions here whose output is a proof.

**Probabilistic algorithms** are information-set decoding searches. Each one
returns the weight of the lowest-weight codeword it found, which is an upper
bound on the distance, together with a witness. They never certify a lower
bound, so a returned value equal to the true distance is not distinguishable
from an unlucky run without further information.

**Heuristic algorithms** are local searches such as genetic algorithms and ant
colony optimization. Like the probabilistic methods they yield upper bounds,
and they are intended for codes far too large for the exact methods.

Bounds discovered by any of these are recorded on the code, so a later call to
an exact algorithm can start from a better upper bound.

## Exact algorithms

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/min_dist_exact.jl"]
Private = false
```

## Probabilistic algorithms

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/min_dist_probabilistic.jl"]
Private = false
```

## Heuristic algorithms

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/min_dist_heuristics.jl"]
Private = false
```
