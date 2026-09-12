# CodingTheory.jl

`CodingTheory.jl` is a Julia library for classical, LDPC, and quantum
error-correcting codes. It uses [Oscar.jl](https://www.oscar-system.org/) for
exact finite-field and polynomial arithmetic and native Julia data structures
for performance-sensitive sparse and iterative algorithms.

## Installation

The package is under active development. Install the development version from
Julia's package prompt:

```julia
] add https://github.com/esabo/CodingTheory
```

Then load the package together with Oscar:

```julia
using Oscar
using CodingTheory
```

Start with [Linear Codes](@ref linear-codes-tutorial),
[Quantum Codes](@ref quantum-codes-tutorial), or
[Message-passing Decoding](@ref message-passing-tutorial). The API pages document constructors and
specialized code families after the tutorials establish the common workflow.

## A first classical code

```julia
F = GF(2)
G = matrix(F, [
    1 0 0 0 0 1 1
    0 1 0 0 1 0 1
    0 0 1 0 1 1 0
    0 0 0 1 1 1 1
])
C = LinearCode(G)

(length(C), dimension(C), minimum_distance(C)[1])
```

Code objects retain the presentation supplied by the user while caching
derived data such as standard forms, logical operators, enumerators, and
certified distance bounds. Use accessors such as `generator_matrix`,
`parity_check_matrix`, `stabilizers`, and `logicals_matrix`; do not depend on
internal struct fields.

## Conventions

- Use `GF(p)` for a prime field. Do not use `GF(p, 1)`: extension-field
  representations are substantially more expensive.
- A parity-check or stabilizer presentation may be overcomplete. Parameters
  are computed from ranks, not from the number of supplied rows.
- Many expensive quantities are cached. Use `copy(C)` when an independent code
  object is needed.
- Prefer predicates and traits such as `is_CSS` and `GaugeTrait` to exact
  `typeof` checks. Constructors may return a more specific supported subtype.
- Exact minimum-distance routines can be exponential. Consult
  [Minimum-distance Computation](@ref minimum-distance-tutorial) before
  running them on large codes.

The Oscar banner can be suppressed by starting Julia with `julia -q`.

## Contributing

Bug reports and contributions are welcome on
[GitHub](https://github.com/esabo/CodingTheory). Development discussion also
takes place in the `#codingtheory` channel of the Julia Slack workspace.
