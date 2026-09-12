# Designs, Self-Dual Codes, and Invariant Theory

The weight enumerator of a self-dual code is invariant under the MacWilliams
transform, so it lies in a ring of invariants of a finite group. Gleason's
theorem identifies generators of that ring, and bounding the coefficients of
the invariants that can be weight enumerators of an actual code gives upper
bounds on the minimum distance of a self-dual code. That is the thread
connecting the three groups of functions below.

The designs functions come from the Assmus-Mattson theorem: under a condition
on the weight distributions of a code and its dual, the supports of the
minimum-weight codewords form a combinatorial design. `design_strength`
reports the resulting strength and `minimum_weight_blocks` returns the blocks.

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/designs.jl"]
Private = false
```

## Gleason bounds

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/Gleason.jl"]
Private = false
```

## Invariant theory

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/invariant_theory.jl"]
Private = false
```
