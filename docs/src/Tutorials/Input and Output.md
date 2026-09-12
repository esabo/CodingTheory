# [Input and Output](@id io-tutorial)

`save_code(path, code; type=...)` is the common export entry point. The `type`
keyword is a symbol and dispatches through `Val(type)`. If it is omitted, the
format is inferred from the filename when possible.

## Matrix CSV

CSV stores one integer matrix and is intended for interchange with numerical
tools:

```julia
using Oscar
using CodingTheory

C = HammingCode(2, 3)
save_code("hamming.csv", C; type=:csv)
```

Classical codes export their generator matrix by default. Use
`code_matrix_array(C; representation=:parity_check)` when selecting a matrix
explicitly. LDPC codes default to their parity-check presentation.

For a quantum code, CSV contains its additive generator matrix:

```julia
F = GF(2)
H = matrix(F, [
    0 0 0 1 1 1 1
    0 1 1 0 0 1 1
    1 0 1 0 1 0 1
])
S = CSSCode(H, H)
save_code("steane.csv", S)
```

CSV is not a lossless code format: it does not preserve the field,
character-vector phases, sparse preference, code subtype, or cached results.

## Portable quantum TOML

TOML is the portable, versioned quantum-code format:

```julia
save_code("steane.toml", S; type=:toml)
restored = load_quantum_code("steane.toml")
```

It preserves the base field, additive generators, character vector, sparse
preference, and validated cache entries. Set `restore_cache=false` when loading
data from an untrusted or incompatible computation.

The dictionary representation is available directly:

```julia
data = quantum_code_data(S)
restored = quantum_code_from_data(data)
```

## Pauli strings

Binary phase-free stabilizer presentations can be exchanged as one Pauli
string per line:

```julia
write_pauli_strings("steane.pauli", S)
restored = read_pauli_strings("steane.pauli")
```

Use `pauli_strings(S)` when an in-memory vector is sufficient. Explicitly
phased codes are rejected rather than silently losing phase information.

## Optional binary formats

JLD2 and NumPy NPZ support are package extensions. Load the corresponding
package before requesting the format:

```julia
using JLD2
save_code("code.jld2", S; type=:jld2)
```

```julia
using NPZ
save_code("code.npz", S; type=:npz)
```

If an optional dependency is not loaded, `save_code` reports which package is
required.
