* `char_vec`: a length `2n` vector with elements in the `Z/(2p)` if
  `characteristic(field(C1))` is 2 and `Z/(p)` otherwise. The first `n` elements
  specify the exponents of the `X` phases and second `n` the exponents of the
  `Z` phases; a missing argument will be set to the all-zero vector

# Notes
* A `+1` phase should be entered as `0` since the character vector stores the
  exponents.
* Stabilizer signs are automatically computed given the character vector.
* The orthogonality of the stabilizers are automatically checked and will error
  upon failure.

assumed to be in
symplectic form over the base field.

* This is intended to be a simple function wrapper for `typeof(S)` since the
 constructor for `SubsystemCode` automatically returns a `SubsystemCodeCSS` if possible.
 Manually changing the elements of the struct `S` without using the helper
 functions provided here is therefore not recommended.

# Stabilizer Codes

```@autodocs
Modules = [CodingTheory]
Pages = ["stabilizercode.jl"]
Private = false
```

# Subsystem Codes

```@autodocs
Modules = [CodingTheory]
Pages = ["subsystemcode.jl"]
Private = false
```

# Graph States

```@autodocs
Modules = [CodingTheory]
Pages = ["graphstate.jl"]
Private = false
```
