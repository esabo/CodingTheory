# [Weight Enumerators and Distributions](@id weight-enumerators-api)

The weight distribution of a code is the multiplicity of each Hamming weight
among its codewords; the weight enumerator is the same data as a polynomial.
Coefficients are kept as exact `BigInt` values, since they grow like the size
of the code.

Enumeration uses bit-packed Gray-code sweeps and the result is cached on the
code. For a high-rate code this is the wrong way around, because the dual has
far fewer codewords: `weight_distribution` therefore computes the dual
distribution and applies the MacWilliams transform automatically when that is
cheaper.

The complete weight distribution refines the Hamming distribution by recording
how many coordinates take each field value, so it is only interesting over
non-binary fields, and it reduces to the Hamming version by summing over the
nonzero values.

Plotting requires a Makie backend to be loaded.

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/words_of_weight.jl"]
Private = false
```

## The MacWilliams transform

The transform maps the Hamming weight distribution of a code to that of its
dual by evaluating Krawtchouk polynomials. `MacWilliams_transform` works on a
`HammingWeightEnumerator` and `MacWilliams_HWE_transform` on a raw
weight-to-multiplicity dictionary; both keep exact `BigInt` coefficients.

## Distributions from a trellis

A trellis gives a second route to the weight distribution that does not
enumerate codewords one at a time; see [Trellises](@ref trellises-api).
