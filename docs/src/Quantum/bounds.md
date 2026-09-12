# Bounds and Weight Enumerators

These functions answer the question of what parameters are possible, as
opposed to what a particular code achieves. They fall into three groups.

*Closed-form bounds* such as the quantum Singleton and quantum Hamming bounds
are cheap arithmetic relations among ``n``, ``k``, and ``d``, with companion
predicates that test whether a given triple satisfies them and existence
bounds of Gilbert-Varshamov type that say when a code must exist.

*Linear-programming bounds* are stronger and more expensive. The Shor-Laflamme
weight enumerator of a quantum code satisfies a set of linear constraints, so
the nonexistence of a code with given parameters can be certified by showing
the corresponding linear program is infeasible. Solving these requires the
`JuMP` extension, which loads when `JuMP` and a solver are available. Take
care with the numerics: enumerator coefficients span many orders of magnitude,
with ``A_0 = 1`` while other terms may reach ``10^{23}``, so the constraint
rows are built exactly in `BigInt` and normalized per row before being handed
to the solver.

*Check-weight bounds* restrict the stabilizer generator weights in addition to
``n``, ``k``, and ``d``, which is the regime relevant to quantum LDPC codes.
These are the bounds where a low-weight constraint genuinely changes the
answer, and the corresponding functions report both the bound and, where
applicable, whether a construction attaining it is known.

```@autodocs
Modules = [CodingTheory]
Pages = ["Quantum/bounds.jl"]
Private = false
```

## Weight enumerators

```@autodocs
Modules = [CodingTheory]
Pages = ["Quantum/weight_enumerators.jl"]
Private = false
```
