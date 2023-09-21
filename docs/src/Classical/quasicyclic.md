# Quasi-Cyclic Codes

Quasi-cyclic codes are a subtype of `LinearCode` and inherit its methods. While quasi-cyclic codes may be seen as generalizations of cyclic codes, here they are treated as independent topics.

## Constructors

```@docs
QuasiCyclicCode
```

## Attributes

```@docs
CodingTheory.index
```

```@docs
expansion_factor
```

```@docs
is_single_generator
```

```@docs
polynomial_matrix
```

```@docs
polynomial_matrix_type
```

The type parameter is either `:G` or `:H`, specifying whether the polynomial matrix represents the generator or parity-check matrix.
```@docs
type
```

## Methods

The following are not computed and stored at the time of construction and must be computed by using these methods.

```@docs
weight_matrix
```

```@docs
noncirculant_generator_matrix
```

```@docs
noncirculant_parity_check_matrix
```

```@docs
CodingTheory.generators
```

```@docs
circulants
```
