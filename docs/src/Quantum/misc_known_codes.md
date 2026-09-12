# Known Quantum Codes

Named constructors return ordinary stabilizer or subsystem code objects and
support the common quantum API.

```@autodocs
Modules = [CodingTheory]
Pages = ["Quantum/misc_known_codes.jl"]
Private = false
```

## Stored lattices

The following constructors read stabilizers, logicals, and metachecks from data
files shipped with the package, so their qubit numbering is fixed and, for the
color codes, chosen to give a small trellis. Run `using JLD2` to activate them.

```@docs
TriangularColorCode488
TriangularColorCode666
PlanarSurfaceCode3D_X
ToricCode3D_X
```
