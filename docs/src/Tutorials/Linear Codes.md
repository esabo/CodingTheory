# Linear Codes

To create a linear code, simply pass a generator matrix into the `LinearCode` constructor.
```
julia> using Oscar, CodingTheory

julia> F = GF(2)
Galois field with characteristic 2

julia> G = matrix(F, [1 0 0 0 0 1 1;
                  0 1 0 0 1 0 1;
                  0 0 1 0 1 1 0;
                  0 0 0 1 1 1 1]);

julia> C = LinearCode(G)
[7, 4, 3]_2 linear code
Generator matrix: 4 × 7
        1 0 0 0 0 1 1
        0 1 0 0 1 0 1
        0 0 1 0 1 1 0
        0 0 0 1 1 1 1
```

We can get the basic information about the code.
```
julia> length(C)
7

julia> dimension(C)
4

julia> cardinality(C)
16

julia> rate(C)
0.5714285714285714
```
Since we passed in a full-rank matrix, the rank should equal the dimension of the code.
```
julia> rank(G) == dimension(C)
true
```
Since the minimum distance of this code is known (since it was small enough to determine in the constructor), we can also get some more information.
```
julia> minimumdistance(C)
3

julia> relativedistance(C)
0.42857142857142855

julia> CodingTheory.genus(C)
1

julia> isMDS(C)
true

julia> numbercorrectableerrors(C)
1
```
We can also manually set the minimum distance using `setminimumdistance!(C, 3)`.

From the output, we see that this is a $[n, k, d] = [7, 4, 3]$ linear code over $\mathbb{F}_2$. The parameters are correctly computed regardless of the input.
```
julia> C2 = LinearCode(vcat(G, G))
[7, 4, 3]_2 linear code
Generator matrix: 8 × 7
        1 0 0 0 0 1 1
        0 1 0 0 1 0 1
        0 0 1 0 1 1 0
        0 0 0 1 1 1 1
        1 0 0 0 0 1 1
        0 1 0 0 1 0 1
        0 0 1 0 1 1 0
        0 0 0 1 1 1 1
```
We can also specify a code by its parity-check matrix
```
julia> H = matrix(F, [0 0 0 1 1 1 1;
                0 1 1 0 0 1 1;
                1 0 1 0 1 0 1]);

julia> C3 = LinearCode(H, true)
[7, 4, 3]_2 linear code
Generator matrix: 4 × 7
        1 1 1 0 0 0 0
        1 1 0 1 0 0 1
        0 1 0 0 1 0 1
        1 0 0 0 0 1 1

julia> paritycheckmatrix(C3)
[0   0   0   1   1   1   1]
[0   1   1   0   0   1   1]
[1   0   1   0   1   0   1]
```

The standard form generator and parity-check matrices are also accessible by passing the optional parameter `true` to each method.
```
julia> generatormatrix(C)
[1   0   0   0   0   1   1]
[0   1   0   0   1   0   1]
[0   0   1   0   1   1   0]
[0   0   0   1   1   1   1]

julia> generatormatrix(C2)
[1   0   0   0   0   1   1]
[0   1   0   0   1   0   1]
[0   0   1   0   1   1   0]
[0   0   0   1   1   1   1]
[1   0   0   0   0   1   1]
[0   1   0   0   1   0   1]
[0   0   1   0   1   1   0]
[0   0   0   1   1   1   1]

julia> generatormatrix(C2, true)
[1   0   0   0   0   1   1]
[0   1   0   0   1   0   1]
[0   0   1   0   1   1   0]
[0   0   0   1   1   1   1]

julia> paritycheckmatrix(C3, true)
[1   0   0   1   1   0   1]
[0   1   0   1   1   1   0]
[0   0   1   0   1   1   1]
```
Recall that column permutations may be required to make the standard form. If this is true, the permutation matrix can be accessed via `standardformpermutation(C)` with the convention that `generatormatrix(C)` and `generatormatrix(C, true) * standardformpermutation(C)` have equivalent row spaces. If no permutation is required, this will return `missing` instead of storing a potentially large identity matrix.

As expected the basic relationship between the matrices holds.
```
julia> iszero(generatormatrix(C) * transpose(paritycheckmatrix(C)))
true

julia> iszero(paritycheckmatrix(C) * transpose(generatormatrix(C)))
true
```

The reader may recognize `C3` as the $[7, 4, 3]$ binary Hamming code.
```
julia> C4 = HammingCode(2, 3)
[7, 4, 3]_2 linear code
Generator matrix: 4 × 7
        1 1 1 0 0 0 0
        1 1 0 1 0 0 1
        0 1 0 0 1 0 1
        1 0 0 0 0 1 1

julia> areequivalent(C3, C4)
true
```
The function `areequivalent` does *not* test if two codes are equivalent up to column permutations.
```
julia> S7 = SymmetricGroup(7)
Sym( [ 1 .. 7 ] )

julia> σ = S7([3, 2, 1, 4, 5, 6, 7])
(1,3)

julia> C3perm = permutecode(C3, σ)
[7, 4]_2 linear code
Generator matrix: 4 × 7
        1 1 1 0 0 0 0
        0 1 1 1 0 0 1
        0 1 0 0 1 0 1
        0 0 1 0 0 1 1

julia> areequivalent(C3perm, C4)
false
```

Of course we know that the Hamming codes are dual to the simplex codes.
```
julia> C5 = SimplexCode(2, 3)
[7, 3, 4]_2 linear code
Generator matrix: 3 × 7
        0 0 0 1 1 1 1
        0 1 1 0 0 1 1
        1 0 1 0 1 0 1

julia> areequivalent(C4, dual(C5))
true
```

A vector $v$ is in the code $C$ if it has zero syndrome.
```
julia> iszero(syndrome(C, generatormatrix(C)[1, :]))
true
```
Similary, we can encode a vector into the codespace.
```
julia> v = encode(C, matrix(F, 1, 4, [1, 0, 0, 0]))
[1   0   0   0   0   1   1]

julia> iszero(syndrome(C, v))
true
```

A code $C_1$ is a subset of $C_2$ if every row of the generator matrix of $C_1$ is in $C_2$.
```
julia> C ⊆ C
true

julia> C ⊆ dual(C)
false
```
Two codes $C_1$ and $C_2$ are equivalent if $C_1 \subseteq C_2$ and $C_2 \subseteq C_1$. A code is self dual if it is equivalent to its dual and self orthogonal if it is a subcode of its dual.
```
julia> isselfdual(C)
false

julia> isselforthogonal(C)
false
```
These are taken with respect to the Euclidean dual/metric/inner product. Similar functions exist for the Hermitian case.
```
julia> C6 = Hexacode()
[6, 3, 4]_4 linear code
Generator matrix: 3 × 6
        1 0 0 1 ω ω
        0 1 0 ω 1 ω
        0 0 1 ω ω 1

julia> isHermitianselfdual(C6)
true
```

To create codes over higher fields, use the `GF(p, l, :α)` constructor. Do not use this when `l = 1`. Note that `α` may be replaced with any symbol.
```
julia> E = GF(2, 3, :α)
Finite field of degree 3 over F_2

julia> α = gen(E)
α

julia> G2 = matrix(E, [α α + 1 1 0 0 0 0;
               0 α α + 1 1 0 0 0;
               0 0 α α + 1 1 0 0;
               0 0 0 α α + 1 1 0;
               0 0 0 0 α α + 1 1])
[α   α + 1       1       0       0       0   0]
[0       α   α + 1       1       0       0   0]
[0       0       α   α + 1       1       0   0]
[0       0       0       α   α + 1       1   0]
[0       0       0       0       α   α + 1   1]

julia> C5 = LinearCode(G2)
[7, 5]_8 linear code
Generator matrix: 5 × 7
        α α + 1 1 0 0 0 0
        0 α α + 1 1 0 0 0
        0 0 α α + 1 1 0 0
        0 0 0 α α + 1 1 0
        0 0 0 0 α α + 1 1
```
As is apparent from the generator matrix, this code is actually cyclic.
```
julia> CodingTheory.iscyclic(C5, false)
true
```

## Reed-Muller Codes
So far, only the standard (binary) Reed-Muller codes have been implemented; the generalized (non-binary) Reed-Muller codes have *not* yet been implemented.

This library constructs Reed-Muller codes using the standard recursive definition of the generator matrices. The literature has conflicting conventions for the base case generator matrix of $\mathcal{RM}(1, 1)$. To use the convention that this should be the identity matrix, set `alt` to `true`; otherwise, $\begin{pmatrix} 1 & 1\\ 0 & 1\end{pmatrix}$ is used.
```
julia> C7 = ReedMullerCode(1, 3)
[8, 4, 4]_2 Reed-Muller code RM(1, 3)
Generator matrix: 4 × 8
        1 1 1 1 1 1 1 1
        0 1 0 1 0 1 0 1
        0 0 1 1 0 0 1 1
        0 0 0 0 1 1 1 1

julia> C8 = ReedMullerCode(1, 3, true)
[8, 4, 4]_2 Reed-Muller code RM(1, 3)
Generator matrix: 4 × 8
        1 0 1 0 1 0 1 0
        0 1 0 1 0 1 0 1
        0 0 1 1 0 0 1 1
        0 0 0 0 1 1 1 1

julia> areequivalent(C7, C8)
true

julia> isselfdual(C7)
true
```

## Modifying Codes And Building New Codes From Old Codes

## Finite Fields And Expanded Codes

At the time of development, the finite field objects used in this library through Oscar do not support any concept of ``relationships``. Some elementary functions for this are provided, although they are intended to only be used for small field sizes.

Let $E$ be an extension field of $F$.
```
julia> F = GF(2)
Galois field with characteristic 2

julia> E = GF(2, 3, :α)
Finite field of degree 3 over F_2

julia> isextension(E, F)
(true, 3)
```
The most two common types of bases for $E/F$ can be computed via
```
julia> primitivebasis(E, F)
(fqPolyRepFieldElem[1, α, α^2], fqPolyRepFieldElem[1, α^2, α])

julia> normalbasis(E, F)
(fqPolyRepFieldElem[α + 1, α^2 + 1, α^2 + α + 1], fqPolyRepFieldElem[α + 1, α^2 + 1, α^2 + α + 1])
```
which return both the basis and its dual (complementary) basis. Alternatively, one specify a basis manually and check its properties.
```
julia> α  = gen(E)
α

julia> β = [α^3, α^5, α^6]
3-element Vector{fqPolyRepFieldElem}:
 α + 1
 α^2 + α + 1
 α^2 + 1

julia> isbasis(E, F, β)
(true, fqPolyRepFieldElem[α + 1, α^2 + α + 1, α^2 + 1])

julia> isselfdualbasis(E, F, β)
true

julia> isprimitivebasis(E, F, β)
false

julia> isnormalbasis(E, F, β)
true

julia> λ = dualbasis(E, F, β)
3-element Vector{fqPolyRepFieldElem}:
 α + 1
 α^2 + α + 1
 α^2 + 1

julia> verifydualbasis(E, F, β, λ)
true

julia> β2 = α .* β
3-element Vector{fqPolyRepFieldElem}:
 α^2 + α
 α^2 + 1
 1

julia> areequivalentbasis(β, β2)
true
```

Using these tools, we can construct expanded codes such as "binary" Reed-Solomon codes.

```
julia> C9 = ReedSolomonCode(8, 3, 5)
[7, 5, 3; 5]_8 Reed-Solomon code
8-Cyclotomic cosets: 
        C_5 ∪ C_6
Generator polynomial:
        x^2 + α*x + α^2 + α
Generator matrix: 5 × 7
        α^2 + α α 1 0 0 0 0
        0 α^2 + α α 1 0 0 0
        0 0 α^2 + α α 1 0 0
        0 0 0 α^2 + α α 1 0
        0 0 0 0 α^2 + α α 1

julia> F8 = field(C9)
Finite field of degree 3 over F_2

julia> α  = gen(F8)
α

julia> β = [field(C9)(1), α, α^6]
3-element Vector{fqPolyRepFieldElem}:
 1
 α
 α^2 + 1

julia> C10 = expandedcode(C9, F, β)
[21, 15]_2 linear code
Generator matrix: 15 × 21
        1 1 1 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        1 1 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0
        1 0 0 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 1 1 1 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 1 1 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0
        0 0 0 1 0 0 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 1 1 1 1 0 1 1 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 1 1 0 0 0 1 0 1 0 0 0 0 0 0 0
        0 0 0 0 0 0 1 0 0 1 1 1 0 0 1 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 1 1 1 1 0 1 1 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 0 1 0 0 0 0
        0 0 0 0 0 0 0 0 0 1 0 0 1 1 1 0 0 1 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 1 1 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 0 1 0
        0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 1 1 0 0 1
```

This is a special example because it is known that this specific exapanded code is also cyclic.
```
julia>  C11 = BCHCode(2, 21, 3, 19)
[21, 15; 19]_2 BCH code
2-Cyclotomic cosets: 
        C_5
Generator polynomial:
        x^6 + x^4 + x^2 + x + 1
Generator matrix: 15 × 21
        1 1 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 1 1 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 1 1 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 1 1 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 1 1 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 1 1 1 0 1 0 1 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 1 1 1 0 1 0 1 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 1 1 1 0 1 0 1 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 1 1 1 0 1 0 1 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 1 1 1 0 1 0 1 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 1 1 1 0 1 0 1 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 1 0 1 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 1 0 1 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 1 0 1 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 1 0 1

julia> areequivalent(C10, C11)
true
```
# TODO: get iscyclic working

# TODO:
* go through runtests and find codes which went something is done to them they are equivalent
