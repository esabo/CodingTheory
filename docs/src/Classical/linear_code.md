# [Linear Codes](@id linear-codes-api)

Generic linear codes may be constructed from a matrix or a vector-space object.
If the optional parameter `parity` is true, a matrix input is interpreted as a
parity-check matrix. Zero rows are removed automatically, while zero columns are
retained. See the [Linear Codes tutorial](@ref linear-codes-tutorial) for usage
examples.

Accessors are provided for every stored attribute, and users are strongly
encouraged to use them rather than reaching into the structs directly: many
functions rely on the stored data being in a specific order, and the accessors
compute and cache derived quantities on first use.

If a code was created from a generator or parity-check matrix, that matrix is
stored alongside the standard form. It is potentially overcomplete, meaning it
has more rows than its rank. Passing `stand_form = true` returns the standard
form instead. Some families are not built from an explicit matrix; there the
matrices start out `missing` and are computed on demand.

Putting a matrix into standard form may require column permutations. When it
does, `standard_form_permutation` returns the permutation matrix ``P`` with
``\mathrm{rowspace}(G) = \mathrm{rowspace}(G_\mathrm{stand} P)``, and `missing`
otherwise.

The minimum distance of some families is known and is set during construction,
and it is computed automatically for codes deemed small enough. Otherwise it is
`missing` and must be requested explicitly; see
[Minimum distance](@ref classical-minimum-distance-api). Functions that depend
on knowing the distance return `missing` when it is unknown. The distance and
its bounds may also be set by hand, in which case nothing is done to check the
value for correctness.

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/types.jl", "Classical/linear_code.jl"]
Private = false
```

The following are not exported but may be useful.

```@docs
CodingTheory.genus
```
