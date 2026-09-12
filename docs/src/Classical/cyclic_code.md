# [Cyclic Codes](@id cyclic-codes-api)

Cyclic codes are a subtype of `LinearCode` and inherit its methods. For a
worked introduction see the [Cyclic Codes tutorial](@ref cyclic-codes-tutorial).

A cyclic code of length ``n`` over ``\mathbb{F}_q`` is determined by its
defining set, a union of ``q``-cyclotomic cosets modulo ``n``. The cyclotomic
functions below are therefore the natural way to specify and inspect these
codes, and several of them are useful on their own when hunting for codes with
prescribed parameters.

Reed-Solomon and BCH codes are the classical special cases. The generalized
Reed-Solomon view of a Reed-Solomon code, along with alternant and Goppa codes,
is documented in
[Generalized Reed-Solomon codes](@ref generalized-reed-solomon-api).

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/cyclic_code.jl", "Classical/cyclotomic.jl"]
Private = false
```

The following are not exported but may be useful.

```@docs
CodingTheory.is_degenerate
CodingTheory.is_cyclic
```
