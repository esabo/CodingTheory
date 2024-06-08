# Cyclic Codes

Cyclic codes are a subtype of `LinearCode` and inherit its methods. For more information on how to use these functions, see the [Cyclic Codes Tutorial](@ref).

## Cyclotomic Cosets

The following set of functions are useful for defining cyclic codes.

```@docs
ord
```

```@docs
cyclotomic_coset
```

```@docs
all_cyclotomic_cosets
```

```@docs
complement_qcosets
```

```@docs
qcoset_pairings
```

```@docs
qcoset_table
```

```@docs
dual_qcosets
```

## Constructors

```@docs
CyclicCode
```

```@docs
BCHCode
```

```@docs
ReedSolomonCode
```

```@docs
QuadraticResidueCode
```

## Attributes

```@docs
splitting_field
```

```@docs
polynomial_ring
```

```@docs
primitive_root
```

```@docs
offset
```

```@docs
design_distance
```

```@docs
qcosets
```

```@docs
qcosets_reps
```

```@docs
defining_set
```

```@docs
zeros
```

```@docs
nonzeros
```

```@docs
generator_polynomial
```

```@docs
parity_check_polynomial
```

```@docs
idempotent
```

```@docs
BCH_bound
```

```@docs
is_narrowsense
```

```@docs
is_reversible
```

```@docs
CodingTheory.is_degenerate
```

```@docs
is_primitive
```

```@docs
is_antiprimitive
```

## Methods

```@docs
defining_set
```

```@docs
dual_defining_set
```

```@docs
CodingTheory.is_cyclic
```

```@docs
complement
```

```@docs
∩
```

```@docs
+
```
