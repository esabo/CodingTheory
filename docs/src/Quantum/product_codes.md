# Product Codes

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


```@autodocs
Modules = [CodingTheory]
Pages = [
    "BB_codes.jl",
    "generalized_3d_toric_codes.jl",
    "concatenated_codes.jl",
    "hypergraph_product_codes.jl",
    "generalized_shor_codes.jl",
    "bicycle_codes.jl",
    "hyperbicycle_codes.jl",
    "lifted_product_codes.jl",
    "fold_product_codes.jl",
    "homological_product_codes.jl",
]
Private = false
```
