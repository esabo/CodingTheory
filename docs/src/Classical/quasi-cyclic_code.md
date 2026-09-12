# Quasi-Cyclic Codes

Quasi-cyclic codes are a subtype of `LinearCode` and inherit its methods. They
may be viewed as a generalization of cyclic codes, but here they are treated as
an independent topic.

A quasi-cyclic code is presented by a matrix over a polynomial quotient ring,
each entry standing for a circulant block. The type parameter is either `:G` or
`:H`, recording whether that polynomial matrix represents the generator or the
parity-check matrix. The noncirculant forms below expand the blocks back out
over the base field; they are not stored at construction and are computed when
first requested.

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/quasi-cyclic_code.jl"]
Private = false
```

The following are not exported but may be useful.

```@docs
CodingTheory.index
CodingTheory.generators
```
