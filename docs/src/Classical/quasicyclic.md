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
expansionfactor
```

```@docs
issinglegenerator
```

```@docs
polynomialmatrix
```

```@docs
polynomialmatrixtype
```

The type parameter is either `:G` or `:H`, specifying whether the polynomial matrix represents the generator or parity-check matrix.
```@docs
type
```

## Methods

The following are not computed and stored at the time of construction and must be computed by using these methods.

```@docs
weightmatrix
```

```@docs
noncirculantgeneratormatrix
```

```@docs
noncirculantparitycheckmatrix
```

```@docs
CodingTheory.generators
```

```@docs
circulants
```
