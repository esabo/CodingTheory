# [Quantum Codes](@id quantum-codes-tutorial)

This tutorial introduces the common stabilizer and CSS-code workflow. Quantum
generators are represented in symplectic form: an ``n``-qubit Pauli operator is
a row of length ``2n`` whose first half is its ``X`` support and whose second
half is its ``Z`` support.

## Stabilizer codes

Small binary codes can be constructed from Pauli strings:

```julia
using Oscar
using CodingTheory

five_qubit = StabilizerCode([
    "XZZXI",
    "IXZZX",
    "XIXZZ",
    "ZXIXZ",
])

(five_qubit.n, five_qubit.k)
stabilizers(five_qubit)
logicals_matrix(five_qubit)
```

Signs on Pauli strings are not accepted as implicit phase data. For
phase-sensitive work, pass an explicit character vector to the matrix
constructor.

The equivalent matrix interface accepts Oscar matrices and ordinary dense or
sparse Julia matrices:

```julia
F = GF(2)
G = matrix(F, [
    1 1 0 0   0 0 0 0
    0 0 1 1   0 0 0 0
    0 0 0 0   1 1 0 0
    0 0 0 0   0 0 1 1
])
bell_pairs = StabilizerCode(G)
```

Construction verifies that the rows commute. Use
`are_symplectic_orthogonal(G, G)` to check a presentation before constructing
the code.

## CSS codes

For a CSS code, provide separate ``X``- and ``Z``-check matrices. The
constructor verifies ``H_X H_Z^T = 0``.

```julia
H = matrix(F, [
    0 0 0 1 1 1 1
    0 1 1 0 0 1 1
    1 0 1 0 1 0 1
])
steane = CSSCode(H, H)

(steane.n, steane.k)
is_CSS(steane)
X_stabilizers(steane)
Z_stabilizers(steane)
```

Use accessors instead of internal fields. Important structural queries include:

```julia
stabilizer_weights(steane)
qubit_degrees(steane)
quantum_LDPC_parameters(steane)
normalizer_matrix(steane)
```

Subsystem codes use the same conventions, with noncommuting generators passed
to `SubsystemCode` or separate gauge sectors passed to `CSSSubsystemCode`.
Their center is available from `stabilizers`, and the full gauge group from
`gauge_group`.

## Distance

Exact quantum distance is the minimum weight of a nontrivial logical operator,
not the minimum weight of a stabilizer. For binary CSS codes, solving the two
sectors separately is usually preferable:

```julia
dX, logical_X = minimum_distance(steane; which=:X, alg=:Gray)
dZ, logical_Z = minimum_distance(steane; which=:Z, alg=:Gray)
d = min(dX, dZ)
```

The witness is returned with the distance and is used to validate cached upper
bounds. Exact methods can be exponential; see
[Minimum-distance Computation](@ref minimum-distance-tutorial) for solver selection, cached bounds, and
probabilistic alternatives.

## Code families

The library includes named small codes and several construction families. For
example:

```julia
surface = ToricCode(3)
shor = ShorCode()
```

BB, hypergraph-product, bicycle, hyperbicycle, homological-product, and
concatenated constructions are documented under the Quantum API. Family
constructors return ordinary stabilizer or subsystem code objects, so the same
accessors and distance routines apply.

## Transforming and saving codes

Quantum direct sums, puncturing, shortening, local Fourier transformations,
gauge fixing, weight reduction, homological measurements, and code expansion
construct new code objects without requiring direct struct manipulation.

Portable TOML preserves the finite field, generators, character vector, sparse
preference, and validated cache entries:

```julia
save_code("steane.toml", steane; type=:toml)
restored = load_quantum_code("steane.toml")
```

CSV is intended for matrix interchange and does not contain enough metadata to
reconstruct a quantum code. Binary Pauli strings can be written and read with
`write_pauli_strings` and `read_pauli_strings`.
