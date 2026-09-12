# [Minimum-distance Computation](@id minimum-distance-tutorial)

Minimum distance is expensive: exact computation is exponential in the worst
case. `CodingTheory.jl` therefore distinguishes exact solvers, witnessed upper
bounds, certified lower bounds, and heuristic searches.

## Classical codes

`minimum_distance` returns the distance and a minimum-weight codeword when a
witness is available:

```julia
using Oscar
using CodingTheory

C = HammingCode(2, 3)
d, word = minimum_distance(C)
```

For small codes, exhaustive enumeration or a weight distribution is often
adequate. For larger codes, information-set and probabilistic algorithms can
find low-weight codewords and tighten the upper bound without certifying the
distance.

The code object caches successful results. Query the current interval with
`minimum_distance_lower_bound(C)` and `minimum_distance_upper_bound(C)`.

## Quantum CSS codes

For a binary CSS stabilizer code, use:

```julia
F = GF(2)
H = matrix(F, [
    0 0 0 1 1 1 1
    0 1 1 0 0 1 1
    1 0 1 0 1 0 1
])
S = CSSCode(H, H)

dX, logical_X = minimum_distance(S; which=:X, alg=:Gray)
dZ, logical_Z = minimum_distance(S; which=:Z, alg=:Gray)
```

The available exact binary CSS methods are:

- `:Gray`: threaded, bit-packed enumeration of normalizer combinations.
- `:Wagner`: a quotient-aware meet-in-the-middle search in physical weight.
- `:ILP`: a JuMP extension using the free HiGHS optimizer.
- `:auto`: selects among the available exact implementations.

Solving the ``X`` and ``Z`` sectors separately exposes useful intermediate
results and is commonly faster than a full search. The quantum distance is
`min(dX, dZ)`.

## Bounds and witnesses

An upper bound is only accepted with a validated nontrivial logical witness:

```julia
set_minimum_distance_upper_bound!(
    S,
    dX,
    logical_X;
    which=:X,
)
```

A mathematically certified lower bound can be recorded separately:

```julia
set_minimum_distance_lower_bound!(S, 2; which=:X)
```

Exact solvers begin at the cached lower bound and search below the witnessed
incumbent. This allows results from different methods and sessions to
cooperate without treating a heuristic result as a proof.

## Probabilistic and heuristic searches

Information-set decoding finds quantum logical witnesses and updates the shared
upper bound:

```julia
d, logical = probabilistic_minimum_distance(
    S;
    which=:X,
    alg=:Stern,
    p=2,
    l=4,
    max_iters=10_000,
    seed=1,
)
```

The binary CSS variants are `:Prange`, `:LeeBrickell`, and `:Stern`.
`heuristic_minimum_distance` additionally provides `:GGAOrder`, `:NNCS`,
`:GA`, and `:ACO`. These searches can improve upper bounds but cannot certify
new lower bounds.

Reusable physical-qubit symmetries can be registered with
`set_distance_automorphisms!`. The setter verifies that every permutation
preserves both CSS stabilizer row spaces.

## Practical guidance

1. Inspect cached bounds before starting a solver.
2. Search ``X`` and ``Z`` separately for CSS codes.
3. Use heuristics to obtain a good witnessed incumbent.
4. Run an exact method only when certification is required.
5. Do not set a solver time limit if an exact answer is required: a timeout is
   inconclusive.
