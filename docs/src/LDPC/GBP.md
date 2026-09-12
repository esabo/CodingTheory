# [Generalized Belief Propagation](@id ldpc-gbp-api)

Ordinary belief propagation minimizes the Bethe free energy, an approximation
that is exact only on a tree. Generalized belief propagation (GBP) instead
works on a *region graph*: vertices are regions, each a set of variable and
check nodes, and edges record containment between regions. Passing messages
between regions rather than between individual nodes accounts for the short
cycles inside each region exactly, which is precisely what ordinary belief
propagation gets wrong on an LDPC code.

Each region carries an overcounting number, the Möbius coefficient that makes
every variable and check counted exactly once across the region graph. A region
graph is valid when those coefficients sum correctly, and the helper functions
below construct candidate region graphs from base regions or clusters, repair
them, and check that condition before decoding.

The cost of GBP grows quickly with region size, so it is best used with small
regions chosen around the problematic cycles identified by the
[cycle and ACE tools](@ref ldpc-cycles-api).

```@autodocs
Modules = [CodingTheory]
Pages = ["LDPC/GBP.jl"]
Private = false
```
