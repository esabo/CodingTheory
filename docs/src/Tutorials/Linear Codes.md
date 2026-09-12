# [Linear Codes](@id linear-codes-tutorial)

A linear ``[n,k,d]_q`` code is a ``k``-dimensional subspace of
``\mathbb F_q^n``. Construct one from a generator matrix or a parity-check
matrix; `CodingTheory.jl` computes ranks and standard forms and retains the
presentation supplied by the user.

## Construction

```julia
using Oscar
using CodingTheory

F = GF(2)
G = matrix(F, [
    1 0 0 0 0 1 1
    0 1 0 0 1 0 1
    0 0 1 0 1 1 0
    0 0 0 1 1 1 1
])
C = LinearCode(G)
```

The basic parameters use Julia's standard interfaces where possible:

```julia
length(C)
dimension(C)
cardinality(C)
rate(C)
```

Constructing from a parity-check matrix requires the second argument:

```julia
H = matrix(F, [
    0 0 0 1 1 1 1
    0 1 1 0 0 1 1
    1 0 1 0 1 0 1
])
C_from_H = LinearCode(H, true)
```

Use `generator_matrix(C)` and `parity_check_matrix(C)` to recover the current
presentation. Passing `true` requests the corresponding standard-form matrix.
Overcomplete inputs are supported; code parameters come from matrix ranks, not
the number of supplied rows.

## Encoding and syndromes

Messages are row vectors over the code's field:

```julia
message = matrix(F, 1, dimension(C), [1, 0, 0, 0])
word = encode(C, message)
iszero(syndrome(C, word))
```

The dual code exchanges generator and parity-check spaces:

```julia
D = dual(C)
iszero(generator_matrix(C) * transpose(generator_matrix(D)))
```

Use `C1 ⊆ C2`, `is_self_orthogonal`, and `is_self_dual` for subspace and
duality queries. `are_equivalent(C1, C2)` tests equality of code spaces in the
current coordinate order; it does not search over column permutations.

## Distance and weight data

```julia
d, minimum_word = minimum_distance(C)
number_correctable_errors(C)
relative_distance(C)
```

Exact distance may require exponential work. A code stores certified lower
bounds and witnessed upper bounds so that multiple algorithms can cooperate.
See [Minimum-distance Computation](@ref minimum-distance-tutorial) before applying exact solvers to large
codes.

Weight enumerators and distributions are also cached:

```julia
W = weight_enumerator(C)
distribution = weight_distribution(C)
```

For some codes it is cheaper to enumerate the dual and apply a MacWilliams
transform. The library uses available structural information and cached data
when selecting a method.

## Named families

Common families have dedicated constructors:

```julia
hamming = HammingCode(2, 3)
simplex = SimplexCode(2, 3)
reed_muller = ReedMullerCode(1, 3)

are_equivalent(hamming, dual(simplex))
```

Cyclic, BCH, Reed--Solomon, generalized Reed--Solomon, Reed--Muller, Gabidulin,
Goppa, quasi-cyclic, Tanner, product, concatenated, and other constructions
have dedicated API pages.

## Extension fields

Use `GF(p)` for prime fields and `GF(p, m, :α)` only for a genuine extension:

```julia
E = GF(2, 3, :α)
α = gen(E)
extension_code = LinearCode(matrix(E, [one(E) α α + 1]))
```

Extension-field matrices are more expensive than prime-field matrices.
Functions including `primitive_basis`, `normal_basis`, `dual_basis`,
`subfield_subcode`, and `expanded_code` support movement between related
fields where mathematically valid.

## New codes from old

Standard constructions include direct sums, products, augmentation,
expurgation, puncturing, shortening, extension, subcodes, and subfield
subcodes. These return new code objects; use `copy(C)` rather than mutating
internal struct fields.

```julia
punctured = puncture(C, 7)
shortened = shorten(C, 7)
extended = extend(C)
```
