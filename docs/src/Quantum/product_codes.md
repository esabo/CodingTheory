# Product and BB Codes

These constructors return ordinary quantum-code objects, so accessors such as
`X_stabilizers`, `Z_stabilizers`, `stabilizer_weights`, and
`minimum_distance` apply uniformly.

- [Hypergraph product](https://errorcorrectionzoo.org/c/hypergraph_product): [Tillich_2014](@cite)
- [Generalized Shor](https://errorcorrectionzoo.org/c/generalized_shor): [bacon2006quantum](@cite)
- Hyperbicycle: [pryadko2013quantum](@cite)
- [Generalized bicycle](https://errorcorrectionzoo.org/c/generalized_bicycle): [pryadko2013quantum](@cite), [Kovalev_2013](@cite), [panteleev2021degenerate](@cite)
- Generalized hypergraph product: [panteleev2021degenerate](@cite)
- Bias-tailored lifted product: [roffe2023bias](@cite)
- Bivariate bicycle (BB): [wang2024coprime](@cite)
- Coprime bivariate bicycle: [wang2024coprime](@cite)

`InfiniteBBCode` and `Generalized3DToricCode` retain algebraic family data
before a finite lattice is selected. `BBCode` and
`FiniteGeneralized3DToricCode` are finite CSS codes and support the standard
quantum-code accessors, distance metadata, and solvers.

For example, the coprime univariate BB form can be constructed over a
polynomial ring:

```julia
using Oscar
using CodingTheory

F = GF(2)
R, z = polynomial_ring(F, :z)
a = one(R) + z
b = one(R) + z^2
S = BBCode(a, b, 7)

(S.n, S.k)
is_CSS(S)
X_stabilizers(S)
```

Finite BB matrices are evaluated lazily and cached on first access. Twisted
two-dimensional lattices use the Laurent-polynomial constructor with two
lattice vectors. The `twisted` property records which finite presentation was
selected; use `twist_vectors(S)` to inspect those vectors.


Quantum Tanner codes are included here as well: like the product
constructions, they build a quantum code out of classical ingredients placed on
a graph, and their distance guarantees come from expansion of that graph.

```@autodocs
Modules = [CodingTheory]
Pages = [
    "Quantum/BB_codes.jl",
    "Quantum/generalized_3d_toric_codes.jl",
    "Quantum/concatenated_codes.jl",
    "Quantum/hypergraph_product_codes.jl",
    "Quantum/generalized_shor_codes.jl",
    "Quantum/bicycle_codes.jl",
    "Quantum/hyperbicycle_codes.jl",
    "Quantum/lifted_product_codes.jl",
    "Quantum/fold_product_codes.jl",
    "Quantum/homological_product_codes.jl",
    "Quantum/Tanner.jl",
]
Private = false
```
