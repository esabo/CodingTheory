# [Trellises](@id trellises-api)

A trellis is a layered graph whose paths are exactly the codewords of a code,
built from a generator matrix in trellis-oriented (minimal-span) form. Its size
depends on the coordinate ordering, so the permutation and sectionalization
functions below search for an ordering that keeps the vertex and edge counts
small. The past and future profiles are what drive that search: they record how
many generators are active on each side of every coordinate boundary.

Because a trellis represents the whole code compactly, it yields weight
distributions and minimum distances without enumerating codewords; see
[Weight enumerators](@ref weight-enumerators-api).

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/trellis.jl"]
Private = false
```
