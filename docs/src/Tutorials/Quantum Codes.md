# Quantum Codes

# Background

## Stabilizer Codes
Like classical codes, a genric stabilizer code is constructed by passing the stabilizers into the constructor. For small codes, it may be convenient to use Pauli strings
```
julia> StabilizerCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"])
[[5, 1]]_2 stabilizer code.
Stabilizer matrix: 4 × 5
         chi(0) 1 0 0 1 0
         chi(0) 0 1 0 0 1
         chi(0) 1 0 1 0 0
         chi(0) 0 1 0 1 0
```
Any $\pm$ signs on the Pauli strings are ignored. More generally,
```
julia> F = GF(2)
Galois field with characteristic 2

julia> stabs = matrix(F, [1 0 0 0 1 0 0 0 1 1 1 1 0 0 0 0;
           0 0 0 1 0 1 0 0 1 0 0 0 0 1 0 0;
           0 1 0 0 1 1 1 0 0 0 1 1 1 0 1 0;
           0 0 1 0 1 1 1 0 0 1 1 0 1 1 0 0;
           0 0 1 1 1 0 1 0 0 0 0 1 0 1 1 1;
           0 0 0 0 0 0 1 1 0 0 1 0 0 0 1 0]);

julia> S = StabilizerCode(stabs)
[[8, 2]]_2 stabilizer code.
Stabilizer matrix: 6 × 8
         chi(0) 1 0 0 0 1 0 0 0
         chi(0) 0 0 0 1 0 1 0 0
         chi(0) 0 1 0 0 1 1 1 0
         chi(0) 0 0 1 0 1 1 1 0
         chi(0) 0 0 1 1 1 0 1 0
         chi(0) 0 0 0 0 0 0 1 1
```
Matrices are required to be written in symplectic form. (A pervious version of this library supported stabilizer codes in quadratic form but the use of extension fields severely limited the practical size of the codes able to be represented.) Errors are thrown for inputs which are not in a possible symplectic form. The inputs must also be symplectic orthogonal. One can use ``symplectic_inner_product`` to check two vectors but ``are_symplectic_orthogonal`` is a more efficient implementation for checking collections of vectors.
```
julia> are_symplectic_orthogonal(stabs, stabs)
true
```
As with classical codes, the ``GF(p)`` constructor is strongly preferred over ``GF(p, 1, :α)``. Over complete matrices are allowed and the code parameters will be correctly computed.

The standard form of the stabilizer matrix as well as the corresponding logical operators are automatically computed during the construction of a code object.
```
julia> logicals(S)
2-element Vector{Tuple{fpMatrix, fpMatrix}}:
 ([0 0 0 0 1 1 0 0 0 1 0 1 0 0 0 0], [0 0 0 0 0 0 0 0 1 1 1 0 1 0 0 0])
 ([0 0 0 0 0 1 1 0 0 1 1 1 0 0 0 1], [0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1])
```
Here, the logicals are stored in $X-Z$ anti-commuting pairs; each pair commuting with the stabilizers and all other logical pairs. These are also available in matrix format with rows stacked in the order $X$ then $Z$ (left to right) from top to bottom.
```
julia> L = logicals_matrix(S)
[0   0   0   0   1   1   0   0   0   1   0   1   0   0   0   0]
[0   0   0   0   0   0   0   0   1   1   1   0   1   0   0   0]
[0   0   0   0   0   1   1   0   0   1   1   1   0   0   0   1]
[0   0   0   0   0   0   0   0   0   1   1   0   0   0   1   1]

julia> are_symplectic_orthogonal(stabs, L)
true

julia> are_symplectic_orthogonal(L, L)
false
```
Several library functions assume this ordering, so these properties should never be accessed directly. If one would like to work with a specific, known form of the logical operators, one may set them using `set_logicals!`. This errors if the automatically computed set of logicals are not equivalent to the function input up to multiplication by stabilizers. One may similiarly `set_stabilizers!`.

As with linear codes, permutations may be required to compute the standard form. If this is the case, the column permutation matrix $P$ such that $\mathrm{rowspace}(stabilizers(S)) = \mathrm{rowspace}(stabilizers(S, true) * standardformpermutation(S))$ may be accessed using the following function. If no column permutations are required, this returns `missing`. The logicals derived from the standard form are always returned in the original qudit ordering of the stabilizer matrix.

Codes encoding no qudits ($k = 0$) are called graph states. Having no logical operators, related functions will not work on these codes. It is unreliable to detect whether or not a code has logical operators using ``typeof``. Instead, one can use traits.
```
julia> LogicalTrait(typeof(S))
HasLogicals()
```
Graph states have trait ``HasNoLogicals()``.

Given two linear codes $\mathcal{C}_1$ and $\mathcal{C}_2$ with $\mathcal{C}_2 ⊆ \mathcal{C}_1$ or a single-linear code $\mathcal{C}$ with $\mathcal{C} \subseteq \mathcal{C}^\perp$, the CSS construction
```
(two examples here)
```
One may also explicitly specify the $X$ and $Z$ stabilizers directly.
```
example for CSSCode(X_matrix, Z_matrix)

example using X_stabilizers, Z_stabilizers
```
Contrary to the above, it is assumed that these matrices are of length $n$ instead of $2n$ with $n$ columns of zeros. No check is done to determine if a code *can* be made CSS, only that the passed in representation *is* CSS. All constructors will automatically return a CSS code if detected. To reliably check whether a code is CSS use `isCSS`.

Stabilizer signs are automatically determined in a consistent manner by a length $2n$ character vector whose first $n$ elements specify the $X$ phase and second $n$ elements the $Z$ phase. This may be passed into any constructor and a missing argument is automatically set to the all-no-phase vector. The signs may be changed after a code has been constructed using ``set_signs!``. Signs for codes over a finite field with characteristic $p$ are $2p$-th roots of unity if $p = 2$ or $p$-th roots of unity otherwise. Since it is more difficult to represent these exactly, signs are stored by keeping track of the exponents.
```
example using the following idea
R = residue_ring(Nemo.ZZ, 4)
char_vec = [R(0) for _ in 1:n]


character_vector(S)
```

## Miscellaneous Known Stabilizer Codes


## Subsystem Codes



## Miscellaneous Known Subsystem Codes
