# LDPC Codes

## Constructors
An LDPC code is defined by a specific choice of a parity-check matrix for a code. Different parity-check matrices for the same linear code produce different LDPC codes. As such, the `LDPCCode` constructor does not accept a code, but rather a matrix.
```@docs
LDPCCode
```

```
julia> H = matrix(GF(2), 6, 9, [
          1 0 1 0 1 0 0 0 1;
          0 1 1 0 1 1 1 0 0;
          0 0 0 1 0 1 0 0 0;
          0 0 0 1 1 0 1 1 0;
          0 1 1 1 0 1 0 0 1;
          1 1 0 0 0 0 1 1 1]);

julia> L = LDPCCode(H)
[9, 3, 3]_2 irregular 5-limited LDPC code with density 0.46296296296296297.

Variable degree polynomial:
        21//25*x^2 + 4//25*x
Check degree polynomial:
        3//5*x^4 + 8//25*x^3 + 2//25*x
Parity-check matrix: 6 × 9
        1 0 1 0 1 0 0 0 1
        0 1 1 0 1 1 1 0 0
        0 0 0 1 0 1 0 0 0
        0 0 0 1 1 0 1 1 0
        0 1 1 1 0 1 0 0 1
        1 1 0 0 0 0 1 1 1
```

Random regular LDPC codes maybe be constructed via
```@docs
regularLDPCCode
```
and irregular LDPC codes via
```@docs
irregularLDPCCode
```

## Attributes
The polynomials ``\lambda(x)`` and ``\rho(x)`` as well as the degrees of each variable and check nodes are computed upon construction.
```@docs
variabledegreepolynomial
```

```@docs
checkdegreepolynomial
```

```@docs
variabledegreedistribution
```

```@docs
checkdegreedistribution
```

```@docs
degreedistributions
```

A bar graph of the degree distributions is available
```@docs
degreedistributionssplot
```

For convenience, the maximum degrees are also stored.
```@docs
columnbound
```

```@docs
rowbound
```

```@docs
columnrowbounds
```

```@docs
limited
```

```@docs
CodingTheory.density
```

```@docs
isregular
```

The Tanner graph corresponding to the parity-matrix defining the LDPC code can be generated as a `SimpleDiGraph` and visually in a `Figure` object.
```@docs
Tannergraph
```

```@docs
Tannergraphplot
```

## Methods
Occassionally useful for small examples, the following function produces a `Figure` of the Tanner graph unrolled to a given level.
```@docs
computationgraph
```

To count or explicitly enumerate the short cycles of the Tanner graph, use
```@docs
countshortcycles
```

```@docs
shortestcycles
```

Various information about the ACE values of cycles in the Tanner graph may be computed with the following functions.
```@docs
ACEspectrumofnode
```

```@docs
shortestcycleACE
```

```@docs
ACEspectrum
```




## Greedy Construction Algorithms
