# Quantum-code API

Quantum code constructors accept symplectic matrices, CSS check pairs, or
binary Pauli strings. Constructors validate commutation and return the most
specific supported code type. Use the common accessors documented below
instead of relying on struct fields.

For a guided introduction, see [Quantum Codes](@ref quantum-codes-tutorial).

## Types and core operations

```@autodocs
Modules = [CodingTheory]
Pages = [
    "Quantum/types.jl",
    "Quantum/stabilizer_code.jl",
    "Quantum/subsystem_code.jl",
]
Private = false
```

## New codes from old

These operations include quantum direct sums, puncturing, shortening,
augmentation and expurgation, local Fourier transformations, and conversions
between stabilizer and subsystem presentations.

```@autodocs
Modules = [CodingTheory]
Pages = ["Quantum/new_codes_from_old.jl"]
Private = false
```

## Distance and bounds

```@autodocs
Modules = [CodingTheory]
Pages = [
    "Quantum/min_dist_bounds.jl",
    "Quantum/min_dist_exact.jl",
    "Quantum/min_dist_probabilistic.jl",
    "Quantum/min_dist_heuristics.jl",
    "Quantum/bounds.jl",
]
Private = false
```

## Weight enumerators

```@autodocs
Modules = [CodingTheory]
Pages = ["Quantum/weight_enumerators.jl"]
Private = false
```

## Input and output

```@autodocs
Modules = [CodingTheory]
Pages = ["Quantum/io.jl"]
Private = false
```
