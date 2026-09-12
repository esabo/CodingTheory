# LDPC Codes

An LDPC code is defined by a specific choice of parity-check matrix for a code.
Different parity-check matrices for the same linear code produce different LDPC
codes, so the `LDPCCode` constructor does not accept a code but rather a
matrix.

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

The degree polynomials ``\lambda(x)`` and ``\rho(x)``, the degrees of the
individual variable and check nodes, and the maximum degrees are all computed
at construction. A bar graph of the degree distributions is available through
`degree_distributions_plot`, which requires a Makie backend to be loaded.

The Tanner graph of the defining parity-check matrix can be produced as a
`SimpleDiGraph` or drawn into a `Figure`, and `computation_graph` draws the
graph unrolled to a given level, which is occasionally useful for small
examples. Cycle structure, girth, and ACE data have their own page; see
[Cycles and ACE](@ref ldpc-cycles-api).

```@autodocs
Modules = [CodingTheory]
Pages = ["LDPC/codes.jl"]
Private = false
```

The following is not exported but may be useful.

```@docs
CodingTheory.density
```

## Construction algorithms

Beyond random regular codes, the library provides the progressive-edge-growth
family, which greedily adds edges so as to maximize the local girth, and
several named algebraic and pseudorandom families: Gallager's original
construction, MacKay-Neal codes, spatially coupled codes, and codes from
Euclidean and projective geometries.

```@autodocs
Modules = [CodingTheory]
Pages = ["LDPC/algorithms.jl"]
Private = false
```
