# Linear Codes

## Constructors

Generic linear codes may be constructed in two ways: via a matrix or via a vector space object. If a vector space is used, the basis of the vector space is used as a generator matrix for the code. If the optional parameter `parity` is set to true, the input is considered a parity-check matrix instead of a generator matrix. At the moment, no convention is used for the zero code and an error is thrown for such imputs. Zero rows are automatically removed from the input but zero columns are not. See the [Linear Codes Over Finite Fields](@ref) for usage examples.
```@docs
LinearCode
```

## Attributes

Various getter/accessor functions are provided for accessing attributes about the codes. The user is strongly encouraged to use these functions and never to work with the underlying structs directly, as many functions rely on the information in the structs to be in a specific order and don't check if information has been updated.

```@docs
field
```

```@docs
length
```

```@docs
dimension
```

```@docs
cardinality
```

```@docs
rate
```

See also: encode`

If the linear code was created by passing in a generator (parity-check) matrix, then this matrix is stored in addition to the standard form. Note that this matrix is potentially over complete (has more rows than its rank). The standard form is returned when the optional parameter `stand_form` is set to true. Some code families are not constructed using these matrices. In these cases, the matrices are initially `missing` and are computed and cached when these functions are called for the first time. Direct access to the underlying structs is not recommended.
```@docs
generator_matrix
```

```@docs
parity_check_matrix
```

```@docs
is_overcomplete
```

Recall that putting the matrix into standard form may require column permutations. If this is the case, the column permutation matrix $P$ such that $\mathrm{rowspace}(G) = \mathrm{rowspace}(G_\mathrm{stand} * P)$ may be accessed using the following function. If no column permutations are required, this returns `missing`.
```@docs
standard_form_permutation
```

The minimum distance of some code families are known and are set during construction. The minimum distance is automatically computed in the constructor for codes which are deemed "small enough". Otherwise, the minimum distance is `missing`. Primitive bounds on the minimum distance are given by

```@docs
minimum_distance_upper_bound
```

If the minimum distance of the code is known, the following functions return useful properties; otherwise they return `missing`.

```@docs
relative_distance
```

```@docs
CodingTheory.genus
```

```@docs
is_MDS
```

```@docs
number_correctable_errors
```

The minimum distance and its bounds may be manually set as well. Nothing is done to check this value for correctness.
```@docs
set_distance_lower_bound!
```

```@docs
set_minimum_distance!
```

See: `minimum_distance_lower_bound`, `set_distance_upper_bound!`

## Methods

```@docs
Singleton_bound
```

```@docs
syndrome
```

```@docs
in
```

```@docs
⊆
```

```@docs
are_equivalent
```

```@docs
dual
```

```@docs
is_self_dual
```

```@docs
is_self_orthogonal
```

```@docs
is_dual_containing
```

```@docs
hull
```

```@docs
is_LCD
```

```@docs
Hermitian_dual
```

```@docs
is_Hermitian_self_dual
```

```@docs
is_Hermitian_self_orthogonal
```

```@docs
is_Hermitian_dual_containing
```

```@docs
Hermitian_hull
```

```@docs
is_Hermitian_LCD
```

```@docs
is_even
```

```@docs
is_doubly_even
```

```@docs
is_triply_even
```

```@docs
characteristic_polynomial
```

```@docs
VectorSpace
```

```@docs
words
```
