# [Input and Output](@id io-api)

`save_code` is the single entry point for writing a code to disk. The
destination format is chosen by the `type` keyword, which is dispatched on
internally, so adding a format does not change the call site. See the
[input and output tutorial](@ref io-tutorial) for worked examples of each format.

The formats differ in what they preserve, and the difference matters:

* `:csv` writes a bare matrix. It is the most portable option and the most
  lossy, since the field, the code type, and any cached metadata are gone.
* `:toml` writes a human-readable, portable description of a quantum code,
  including the field, the checks, and the parameters, along with a SHA-256
  fingerprint of the check data so that corruption is detectable on load.
* `:pauli` writes stabilizers as Pauli strings, which is the interchange format
  most other quantum software understands.
* `:jld2` and `:npz` (also accepted as `:nz`) require their respective package
  extensions to be loaded. `:jld2` round-trips native Julia objects, while
  `:npz` targets NumPy consumers; its reader accepts SciPy's zero-based CSR
  index convention.

```@autodocs
Modules = [CodingTheory]
Pages = ["src/io.jl"]
Private = false
```

## Quantum codes

```@autodocs
Modules = [CodingTheory]
Pages = ["Quantum/io.jl"]
Private = false
```
