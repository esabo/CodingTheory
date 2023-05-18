# Linear Codes

## Constructors

Generic linear codes may be constructed in two ways: via a matrix or via a vector space object. If a vector space is used, the basis of the vector space is used as a generator matrix for the code. If the optional parameter `parity` is set to true, the input is considered a parity-check matrix instead of a generator matrix. At the moment, no convention is used for the zero code and an error is thrown for such imputs. Zero rows are automatically removed from the input but zero columns are not. See the [tutorials](link) for usage examples.
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

If the linear code was created by passing in a generator (parity-check) matrix, then this matrix is stored in addition to the standard form. Note that this matrix is potentially over complete (has more rows than its rank). The standard form is returned when the optional parameter `standform` is set to true. Some code families are not constructed using these matrices. In these cases, the matrices are initially `missing` and are computed and cached when these functions are called for the first time. Direct access to the underlying structs is not recommended.
```@docs
generatormatrix
```

```@docs
paritycheckmatrix
```

Recall that putting the matrix into standard form may require column permutations. If this is the case, the column permutation matrix $P$ such that $\mathrm{rowspace}(G) = \mathrm{rowspace}(G_\mathrm{stand} * P)$ may be accessed using the following function. If no column permutations are required, this returns `missing`.
```@docs
standardformpermutation
```

The minimum distance of some code families are known and are set during construction. The minimum distance is automatically computed in the constructor for codes which are deemed "small enough". Otherwise, the minimum distance is `missing`. Primitive bounds on the minimum distance are given by
```@docs
minimumdistancelowerbound
```

```@docs
minimumdistanceupperbound
```

If the minimum distance of the code is known, the following functions return useful properties; otherwise they return `missing`.

```@docs
relativedistance
```

```@docs
genus
```

```@docs
isMDS
```

```@docs
numbercorrectableerrors
```

The minimum distance and its bounds may be manually set as well. Nothing is done to check this value for correctness.
```@docs
setdistancelowerbound!
```

```@docs
setdistanceupperbound!
```

```@docs
setminimumdistance!
```

## Methods

```@docs
Singletonbound
```

```@docs
encode
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
areequivalent
```

```@docs
dual
```

```@docs
isselfdual
```

```@docs
isselforthogonal
```

```@docs
isdualcontaining
```

```@docs
hull
```

```@docs
isLCD
```

```@docs
Hermitiandual
```

```@docs
isHermitianselfdual
```

```@docs
isHermitianselforthogonal
```

```@docs
isHermitiandualcontaining
```

```@docs
Hermitianhull
```

```@docs
isHermitianLCD
```

```@docs
iseven
```

```@docs
isdoublyeven
```

```@docs
istriplyeven
```

```@docs
characteristicpolynomial
```

```@docs
VectorSpace
```

```@docs
words
```
