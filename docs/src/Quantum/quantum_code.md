# [Quantum Code API](@id quantum-code-api)

Quantum code constructors accept symplectic matrices, CSS check pairs, or
Pauli strings. They validate commutation and return the most specific
supported code type, so a constructor given data with no gauge operators
returns a stabilizer code, and a stabilizer code whose checks split by type
returns a CSS code. Use the accessors rather than reaching into struct fields,
and prefer dispatching on the traits (`LogicalTrait`, `GaugeTrait`,
`CSSTrait`) over testing concrete types.

For a guided introduction, see [Quantum Codes](@ref quantum-codes-tutorial).

Because every stabilizer code is a subsystem code with no gauge qubits, the
bulk of the shared accessors are documented on the
[subsystem code page](@ref quantum-subsystem-api); this page covers the type
hierarchy, the traits, and the stabilizer-specific constructors and solvers.

```@autodocs
Modules = [CodingTheory]
Pages = ["Quantum/types.jl", "Quantum/stabilizer_code.jl"]
Private = false
```
