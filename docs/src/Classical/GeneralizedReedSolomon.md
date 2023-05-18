# Generalized Reed-Solomon Codes

Generalized Reed-Solomon codes are a subtype of `LinearCode` and inherit its methods.

## Constructors

```@docs
GeneralizedReedSolomonCode
```

## Attributes

```@docs
scalars
```

```@docs
dualscalars
```

```@docs
evaluationpoints
```

## Methods

The dual of a generalized Reed-Solomon code is another generalized Reed-Solomon code. This function override the default `LinearCode` method to take this into account.
```@docs
dual
```
