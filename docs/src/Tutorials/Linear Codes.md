# Linear Codes Over Finite Fields

## Background
The following represents the notation and conventions used for linear codes throughout the library.

A (classical) error correcting code $\mathcal{C}$ is a $k$-dimensional subspace of $\mathbb{F}^n_q$. Elements of $\mathcal{C}$ are called codewords. The number of codewords in $\mathcal{C}$ is denoted $|\mathcal{C}|$. The dimension of $\mathcal{C}$, $\mathrm{dim}(\mathcal{C})$, is defined to be the dimension of $\mathcal{C}$ as a vector space over $\mathbb{F}_q$, i.e., $\displaystyle |\mathcal{C}| = q^{\mathrm{dim}(\mathcal{C})}$. It is customary to denote $\mathrm{dim}(\mathcal{C})$ by $k$ such that $\mathcal{C}$ is an $(n, q^k)_q$ code, or an $[n, k]_q$ code. The notation of choice depends on whether or not it is easier to make an argument about $|\mathcal{C}|$ or $\mathrm{dim}(\mathcal{C})$, although here we always stick to the latter. An $[n, k]_q$ code is written $[n, k]$ when $q = 2$.

A $k \times n$ matrix $G$ is a generator matrix for $\mathcal{C}$ if $\mathcal{C}$ is the row space of $G$. An $(n - k) \times n$ parity check matrix $H$ for $\mathcal{C}$ is a generator matrix for the row space of the vector space orthogonal to $\mathcal{C}$ in $\mathbb{F}_q^n$ with respect to the standard Euclidean inner product, $\mathcal{C}^\perp$, i.e., $\mathcal{C} = \mathrm{ker} H$. This is called the dual code of $\mathcal{C}$ and the generator and parity-check matrices of $\mathcal{C}$ and $\mathcal{C}^\perp$ are switched. A code is called self-orthogonal if $\mathcal{C} \subseteq \mathcal{C}^\perp$ and self-dual if $\mathcal{C} = \mathcal{C}^\perp$. The orthogonality of $\mathcal{C}$ and $\mathcal{C}^\perp$ gives $G^T H = G H^T = 0$. The product $Hv$ is called the syndrome of $v$ and a zero syndrome implies $v \in \mathcal{C}$.

A generator matrix is said to be in standard form if $G = (I_k \mid A)$, where $I_k$ is the $k \times k$ identity matrix, and a parity-check matrix is said to be in standard form if $H = (B \mid I_{n - k})$. The relationship between $G$ and $H$ gives $B = -A^T$. By elementary row *and column* operations, any linear code is equivalent to a linear code with a generator matrix in standard form.

It is often convenient to define a code using a matrix with linearly dependent rows. In this case, we say that the matrix, or code, is over complete. The standard form matrices represent a basis for the row space and cannot be over complete. The code whose only element is the zero vector is called the zero code. The $1 \times n$ zero matrix is an over complete generator matrix for this code. Since the zero vector cannot be part of a basis by definition and the span of the empty set is zero, the standard form of this code is given by a $0 \times n$ matrix. This causes no problems with the library nor the underlying Oscar framework.

The (Hamming) weight of $x \in \mathbb{F}^n_q$, $\mathrm{wt}(x)$, is the number of nonzero components in the vector. The (Hamming) distance between $x \in \mathbb{F}^n_q$ and $y \in \mathbb{F}^n_q$, denoted by $d(x, y)$, is defined to be the number of places at which $x$ and $y$ differ, i.e.,  $d(x, y) = \mathrm{wt}(x - y)$. For a code $\mathcal{C}$ with $|\mathcal{C}| \geq 2$, the (minimum) distance of $\mathcal{C}$, denoted by $d = d(\mathcal{C})$, is

$$d(\mathcal{C}) = \min \{d(x, y) \mid x, y \in \mathcal{C}, x \neq y\} = \min \{\mathrm{wt}(c) \mid c \in \mathcal{C}\},$$

where the second equality holds only for the linear codes considered in this work. An $[n, k]_q$ code with minimum weight $d$ is denoted by $[n, k, d]_q$. The homogenous, Hamming weight enumerator of $\mathcal{C}$ is the bivariate polynomial

$$W(\mathcal{C}; x, y) = \sum_{i = 0}^n A_i x^i y^{n - i},$$

where $A_i$ is the number of elements of $\mathcal{C}$ with weight $i$. The weight distribution of $\mathcal{C}$ is the ordered sequence $\{A_i\}_{i = 0}^n$. The minimum distance is hence the smallest index $i$ such that $A_i \neq 0$. The weight enumerator of $\mathcal{C}$ and $\mathcal{C}^\perp$ are related via the MacWilliams identity

$$W(\mathcal{C}^\perp; x, y) = \frac{1}{|C|} W(\mathcal{C}; y - x, y + x).$$

The complete weight enumerator is the multivariate polynomial

$$\sum_{i = 0}^n A_i x^{\mathrm{wt}_1}_1 x^{\mathrm{wt}_2}_2 \dots x^{\mathrm{wt}_{|\mathbb{F}|}}_{|\mathbb{F}|},$$

where $A_i$ is the number of elements which have $\mathrm{wt}_j$ occurences of element $j$ in some fixed enumeration of the field. The convention used here is the enumeration defined by applying ``collect`` to the field. MacWilliams identities also exist for complete weight enumerators, but are more complicated. Given a complete weight enumerator, one may always derive the Hamming weight enumerator, but the Hamming weight enumerator is not enough to uniquely specify the complete weight enumerator.

Information is encoded in $\mathcal{C}$ via $\mathrm{enc}: \mathbb{F}^k_q \to \mathbb{F}^n_q$, $v \mapsto vG$. The parameter $d$ is related to the error-correcting process called decoding, which should not be confused with "unencoding".
!!! note "Theorem"
	An $[n, k, d]_q$ code can correct $t = \lfloor(d -1)/2\rfloor$ or fewer errors. Conversely, a code which can correct $t = \lfloor(d -1)/2\rfloor$ or fewer errors has minimum weight $d$.

The variables $n$, $k$, $d$, and $t$ will be reserved for these quantities.

The most common form of extending a code is to add an extra column to the generator matrix such that the sum of the coordinates of each row is 0. Augmenting a code adjoins rows to the generator matrix. Expurgating a code deletes rows from the generator matrix and then removes any potentially resulting zero columns. Puncturing a code deletes columns from the generator matrix and then removes any potentially resulting zero rows. Shortening is expurgating followed by puncturing. Codes with a single punctured column are often denoted by $\mathcal{C}^*$. Shortened codes are often denoted by $\overline{\mathcal{C}}$.

Let $\mathcal{C}$ be an $[n, k, d]_{p^m}$ code. Then the subfield subcode of $\mathcal{C}$ over a subfield $\mathbb{F} \leq \mathbb{F}_{p^m}$, denoted $\mathcal{C}|_\mathbb{F}$, is given by $\mathcal{C} \cap \mathbb{F}^n$, i.e., the collection of codewords of $\mathcal{C}$ whose components lie entirely in $\mathbb{F}$. The code $\mathcal{C}$ is called the supercode of $\mathcal{C}|_\mathbb{F}$. If $\mathcal{C}$ has parameters $[n, k, d]_{p^m}$, $\mathcal{C}|_\mathbb{F}$ has parameters $[n, k^\prime, \geq d]$ over $\mathbb{F}$, where $n - k \leq n - k^\prime \leq \ell (n - k)$ and $\ell = [\mathbb{F}_{p^m} : \mathbb{F}]$ (the index of $\mathbb{F}$ in $\mathbb{F}_{p^m}$). As the codewords of $\mathcal{C}|_\mathbb{F}$ are codewords of $\mathcal{C}$, it follows immediately that the minimum distance of $\mathcal{C}|_\mathbb{F}$ is at least the minimum distance of $\mathcal{C}$, and $\mathcal{C}|_\mathbb{F}$ can be decoded using the same algorithm as $\mathcal{C}$, although perhaps not efficiently as a native algorithm over $\mathbb{F}$ designed specifically for the subfield subcode.

An $[mn, mk, \geq d]_p$ code may be constructed from an $[n, k, d]_{p^m}$ code by expanding its elements using a basis of $\mathbb{F}_{p^m}/\mathbb{F}_p$. The first code is called the expanded code of the second. To see why the minimum distance of the code could increase, let $\beta = \{\beta_j\}_1^m$ be a basis of $\mathbb{F}_{p^m}/\mathbb{F}_p$ and let $c = (c_1, \dots, c_n) \in \mathbb{F}^n_{p^m}$ be a minimum weight codeword in an $[n, k, d]_{p^m}$ code. Expressing each $c_i$ with respect to $\beta$, $c_i = \sum_j c_{ij} \beta_j$, we can replace each element with its corresponding $m$-tuple, $(c_{i1}, \dots, c_{im})$. If  $c_i \neq 0$, then the Hamming weight of its expansion is at least one and therefore the Hamming weight of the expansion of $c$ is at least $d$.

Recall that the inner product over finite fields is given by the trace. In particular, if $\beta$ is a basis of $\mathbb{F}_{p^m}/\mathbb{F}_p$ such that $x = \sum_{j = 1}^m x_j \beta_j$ for $x \in \mathbb{F}_{p^m}$, then $x_j = \mathrm{Tr}_{\mathbb{F}_{p^m}/\mathbb{F}_p} (x \beta_j^\perp) \in \mathbb{F}_p$, where $\beta^\perp$ is the unique trace-orthogonal dual of $\beta$ such that $\mathrm{Tr}_{\mathbb{F}_{p^m}/\mathbb{F}_p}(x_i y_j) = \delta_{ij}$ for $x_i \in \beta$ and $y_j \in \beta^\perp$. (The dual basis always exists and is easy to compute given $\beta$.) For $c = (c_1, \dots, c_n) \in \mathbb{F}_{p^m}^n$ denote the expansion with respect to $\beta$ by the isomorphism $\phi_\beta: \mathbb{F}_{p^m}^n \to \mathbb{F}_p^{nm}$ given by

$$\begin{aligned}
	\phi_\beta(c) &= (\phi_\beta(c_1), \dots, \phi_\beta(c_n))\\
		&= (\mathrm{Tr}_{\mathbb{F}_{p^m}/\mathbb{F}_p}(c_1 \beta_1^\perp), \dots, \mathrm{Tr}_{\mathbb{F}_{p^m}/\mathbb{F}_p}(c_1 \beta_m^\perp), \mathrm{Tr}_{\mathbb{F}_{p^m}/\mathbb{F}_p}(c_2 \beta_1^\perp), \dots, \mathrm{Tr}_{\mathbb{F}_{p^m}/\mathbb{F}_p}(c_n \beta_m^\perp)).
\end{aligned}$$

Generator and parity check matrices for expanded codes are given by

$$G_\phi = \begin{pmatrix}
		\phi_\beta(\beta_1 g_1)\\
		\vdots\\
		\phi_\beta(\beta_m g_1)\\
		\phi_\beta(\beta_1 g_2)\\
		\vdots\\
		\phi_\beta(\beta_m g_k)
	\end{pmatrix}
	\qquad , \qquad
	H_\phi = \begin{pmatrix}
		\phi_\beta^\perp(\beta^\perp_1 h_1)\\
		\vdots\\
		\phi_\beta^\perp(\beta^\perp_m h_1)\\
		\phi_\beta^\perp(\beta^\perp_1 h_2)\\
		\vdots\\
		\phi_\beta^\perp(\beta^\perp_m h_{n - k})
	\end{pmatrix}.$$

In general, an expanded code loses the properties of its parent code and different bases could produce different expanded codes with different parameters and properties. It is still not yet known how to choose a basis to *a priori* maximize the minimum distance of the expanded code. One crucial property that might not be maintained by a basis expansion is orthogonality. To see this, let $\beta$ be an arbitrary basis for $\mathbb{F}_{p^m}/\mathbb{F}_p$. If $\mathcal{C}_2 \subseteq \mathcal{C}_1$ over $\mathbb{F}_{p^m}$, then $\phi_\beta(\mathcal{C}_2) \subseteq \phi_\beta(\mathcal{C}_1)$ over $\mathbb{F}_p$ trivially, since if $x \in \mathcal{C}_2$ then $x \in \mathcal{C}_1$ and $\phi_\beta(x) \in \phi_\beta(\mathcal{C}_2)$ and $\phi_\beta(x) \in \phi_\beta(\mathcal{C}_1)$. It is well-known in classical coding theory, and can be verified by direct computation, that $(\phi_\beta(\mathcal{C}))^\perp = \phi_{\beta^\perp}(\mathcal{C}^\perp)$. Now suppose $\mathcal{C} \subseteq \mathcal{C}^\perp$. Then $\phi_\beta(\mathcal{C}) \subseteq \phi_\beta(\mathcal{C}^\perp)$ and $\phi_\beta(\mathcal{C})$ is self-orthogonal if and only if $\phi_\beta(\mathcal{C}^\perp) \subseteq (\phi_\beta(\mathcal{C}))^\perp = \phi_{\beta^\perp}(\mathcal{C}^\perp)$. It is sufficient for $\beta = \beta^\perp$ but not every field extension has a self-dual basis. Even if a self-dual basis for the extension exists, it is often difficult to find. The two most common bases are the polynomial bases of the form $\{1, \alpha, \dots, \alpha^{m - 1}\}$ and the normal bases of the form $\{\alpha, \alpha^p, \alpha^{p^2}, \dots, \alpha^{p^{m - 1}}\}$. If $\alpha$ is primitive, then the polynomial basis is called a primitive (polynomial) basis.

## Basics
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
julia> minimum_distance(C)
3

julia> relative_distance(C)
0.42857142857142855

julia> CodingTheory.genus(C)
1

julia> isMDS(C)
true

julia> number_correctable_errors(C)
1
```
We can also manually set the minimum distance using `set_minimum_distance!(C, 3)`.

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

julia> parity_check_matrix(C3)
[0   0   0   1   1   1   1]
[0   1   1   0   0   1   1]
[1   0   1   0   1   0   1]
```

The standard form generator and parity-check matrices are also accessible by passing the optional parameter `true` to each method.
```
julia> generator_matrix(C)
[1   0   0   0   0   1   1]
[0   1   0   0   1   0   1]
[0   0   1   0   1   1   0]
[0   0   0   1   1   1   1]

julia> generator_matrix(C2)
[1   0   0   0   0   1   1]
[0   1   0   0   1   0   1]
[0   0   1   0   1   1   0]
[0   0   0   1   1   1   1]
[1   0   0   0   0   1   1]
[0   1   0   0   1   0   1]
[0   0   1   0   1   1   0]
[0   0   0   1   1   1   1]

julia> generator_matrix(C2, true)
[1   0   0   0   0   1   1]
[0   1   0   0   1   0   1]
[0   0   1   0   1   1   0]
[0   0   0   1   1   1   1]

julia> parity_check_matrix(C3, true)
[1   0   0   1   1   0   1]
[0   1   0   1   1   1   0]
[0   0   1   0   1   1   1]
```
Recall that column permutations may be required to make the standard form. If this is true, the permutation matrix can be accessed via `standard_form_permutation(C)` with the convention that `generator_matrix(C)` and `generator_matrix(C, true) * standard_form_permutation(C)` have equivalent row spaces. If no permutation is required, this will return `missing` instead of storing a potentially large identity matrix.

As expected the basic relationship between the matrices holds.
```
julia> iszero(generator_matrix(C) * transpose(parity_check_matrix(C)))
true

julia> iszero(parity_check_matrix(C) * transpose(generator_matrix(C)))
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

julia> are_equivalent(C3, C4)
true
```
The function `are_equivalent` does *not* test if two codes are equivalent up to column permutations.
```
julia> S7 = SymmetricGroup(7)
Sym( [ 1 .. 7 ] )

julia> σ = S7([3, 2, 1, 4, 5, 6, 7])
(1,3)

julia> C3_perm = permute_code(C3, σ)
[7, 4]_2 linear code
Generator matrix: 4 × 7
        1 1 1 0 0 0 0
        0 1 1 1 0 0 1
        0 1 0 0 1 0 1
        0 0 1 0 0 1 1

julia> are_equivalent(C3_perm, C4)
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

julia> are_equivalent(C4, dual(C5))
true
```

A vector $v$ is in the code $C$ if it has zero syndrome.
```
julia> iszero(syndrome(C, generator_matrix(C)[1, :]))
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
julia> is_self_dual(C)
false

julia> is_self_orthogonal(C)
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

julia> is_Hermitian_self_dual(C6)
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
julia> is_cyclic(C5, false)
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

julia> are_equivalent(C7, C8)
true

julia> is_self_dual(C7)
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

julia> is_extension(E, F)
(true, 3)
```
The most two common types of bases for $E/F$ can be computed via
```
julia> primitive_basis(E, F)
(fqPolyRepFieldElem[1, α, α^2], fqPolyRepFieldElem[1, α^2, α])

julia> normal_basis(E, F)
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

julia> is_basis(E, F, β)
(true, fqPolyRepFieldElem[α + 1, α^2 + α + 1, α^2 + 1])

julia> is_self_dual_basis(E, F, β)
true

julia> is_primitive_basis(E, F, β)
false

julia> is_normal_basis(E, F, β)
true

julia> λ = dual_basis(E, F, β)
3-element Vector{fqPolyRepFieldElem}:
 α + 1
 α^2 + α + 1
 α^2 + 1

julia> verify_dual_basis(E, F, β, λ)
true

julia> β2 = α .* β
3-element Vector{fqPolyRepFieldElem}:
 α^2 + α
 α^2 + 1
 1

julia> are_equivalent_basis(β, β2)
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

julia> C10 = expanded_code(C9, F, β)
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

julia> are_equivalent(C10, C11)
true
```

## Weight Enumerators, Distributions, And Minimum Distance


