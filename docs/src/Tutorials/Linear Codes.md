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
We see that this is a $[n, k, d] = [7, 4, 3]$ linear code over $\mathbb{F}_2$. The parameters are correctly computed regardless of the input.
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
