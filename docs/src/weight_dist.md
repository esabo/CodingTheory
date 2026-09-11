# Weight Enumerators and Distributions

Classical weight routines enumerate codewords with bit-packed Gray-code
methods and cache the resulting exact `BigInt` coefficients. For some inputs,
computing the dual distribution and applying the MacWilliams transform is
cheaper.

```@docs
HammingWeightEnumerator
weight_distribution
weight_distribution_array
weight_enumerator
complete_weight_distribution
MacWilliams_transform
words_of_weight
```

`complete_weight_enumerator` constructs the corresponding multivariate
polynomial from a complete distribution.
