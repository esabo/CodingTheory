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
regular_LDPC_code
```
and irregular LDPC codes via
```@docs
irregular_LDPC_code
```

## Attributes
The polynomials ``\lambda(x)`` and ``\rho(x)`` as well as the degrees of each variable and check nodes are computed upon construction.
```@docs
variable_degree_polynomial
```

```@docs
check_degree_polynomial
```

```@docs
variable_degree_distribution
```

```@docs
check_degree_distribution
```

```@docs
degree_distributions
```

A bar graph of the degree distributions is available
```@docs
degree_distributions_plot
```

For convenience, the maximum degrees are also stored.
```@docs
column_bound
```

```@docs
row_bound
```

```@docs
column_row_bounds
```

```@docs
limited
```

```@docs
CodingTheory.density
```

```@docs
is_regular
```

The Tanner graph corresponding to the parity-matrix defining the LDPC code can be generated as a `SimpleDiGraph` and visually in a `Figure` object.
```@docs
Tanner_graph
```

```@docs
Tanner_graph_plot
```

## Methods
Occassionally useful for small examples, the following function produces a `Figure` of the Tanner graph unrolled to a given level.
```@docs
computation_graph
```

```@docs
girth
```

To count or explicitly enumerate the short cycles of the Tanner graph, use
```@docs
count_short_cycles
```

```@docs
shortest_cycles
```

Various information about the ACE values of cycles in the Tanner graph may be computed with the following functions.
```@docs
ACE_spectrum
```

```@docs
shortest_cycle_ACE
```

```@docs
ACE_distribution
```

```@docs
average_ACE_distribution
```

```@docs
median_ACE_distribution
```

```@docs
mode_ACE_distribution
```

## Greedy Construction Algorithms
