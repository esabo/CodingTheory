# New Quantum Codes From Old

These operations include quantum direct sums, puncturing, shortening,
augmentation and expurgation, local Fourier transformations, exchanging the
``X`` and ``Z`` sectors, and conversions between stabilizer and subsystem
presentations.

Modifying a quantum code is more delicate than modifying a classical one,
because the result must still be a valid code: the checks have to remain
mutually commuting, and the number of logical qubits changes with the number of
independent checks. The constructors here validate that and recompute the
dependent data rather than copying stale values, so the returned code carries
no distance information it cannot justify.

```@autodocs
Modules = [CodingTheory]
Pages = ["Quantum/new_codes_from_old.jl"]
Private = false
```
