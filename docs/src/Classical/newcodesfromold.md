# Modifying Codes

```@autodocs
Modules = [CodingTheory]
Pages = ["newcodesfromold.jl"]
Private = false
```

Notes:
- notice how the generalized `extend` works (explained in its docstring)
- `expurgate` and `augment` both are working with potentially overcomplete generator matrices (not directly on the standard form)
- `uuplusv` throws an error when one of the arguments is the zero code (there is a reasonable way to define it, but it's unnecessary)
- `juxtaposition` works on potentially overcomplete generator matrices (not on the standard form)
- Schur product is probably not valid except for `C * C`
