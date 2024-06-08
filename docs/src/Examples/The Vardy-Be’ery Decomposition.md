# The Vardy-Be’ery Decomposition

## Background
This example draws from and uses the notation from the background sections in the [Linear Codes Over Finite Fields](@ref) and [Cyclic Codes](@ref) tutorials.

## The Almost-Block-Diagonal Form
Vardy and Be'ery showed that expanded (cyclic) Reed-Solomon codes may be seen as interleaved BCH codes plus some extra "glue".
!!! note "The Vardy-Be’ery Decomposition"
	Let $\mathcal{C}$ be a cyclic Reed-Solomon code over $\mathbb{F}_{p^m}$ for $p, m \geq 2$. Using column permutations, the expanded generator matrix can be put into the form of Figure \ref{fig:blockstruct}, where $B$ is the generator matrix of the corresponding BCH subfield subcode $\mathcal{B}$.

The proof of this is easy. Instead of expanding $\mathcal{C}$ in its entirety, first perform the expansion on $\mathcal{B} \subset \mathcal{C}$ only. As a subfield subcode, the codewords of $\mathcal{B}$ are closed in $\mathcal{C}$ under scalar multiplication in $\mathbb{F}_{p^m}$. Choose a basis, $\beta$, of $\mathbb{F}_{p^m}/\mathbb{F}_p$ and let $\mathcal{B}_i = \{ \beta_i b \, : \, b \in \mathcal{B} \}$. Clearly the $\mathcal{B}_i$ are disjoint subcodes of $\mathcal{C}$. For $b = (b_1, \dots, b_n)$, a row of the generator matrix of $\mathcal{B}$, the expansion of $\beta_i b$ gives $m$ rows of the form

$$\begin{aligned}
	&(b_1, 0, \dots, 0, b_2, 0, \dots, 0, b_n, 0, \dots, 0),\\
	&(0, b_1, 0, \dots, 0, b_2, 0, \dots, 0, b_n, 0, \dots, 0),\\
	&(0, 0, b_1, 0, \dots, 0, b_2, 0, \dots, 0, b_n, 0, \dots, 0),
\end{aligned}$$

where we have used the fact that $\mathrm{Tr}_{\mathbb{F}_{p^m}/\mathbb{F}_p}$ is $\mathbb{F}_p$-linear and $\mathrm{Tr}_{\mathbb{F}_{p^m}/\mathbb{F}_p}(b_\ell \beta_i \beta_j^\perp) = b_\ell \delta_{ij}$. There are $m - 1$ zeros (cyclically) between each non-zero element. Permuting columns put these into the form

$$\begin{aligned}
	&(b_1, \dots, b_n, 0, \dots, 0, 0, \dots 0),\\
	&(0, \dots, 0, b_1, \dots, b_n, 0, \dots 0),\\
	&(0, \dots, 0, 0, \dots 0, b_1, \dots, b_n).
\end{aligned}$$

Repeating this for all of the rows of the generator matrix of $\mathcal{B}$ then permuting rows completes the $m$ factors of $B \oplus \dots \oplus B$.

If $\mathcal{C}$ has dimension $k$ and $\mathcal{B}$, $k^\prime$, then $m(k - k^\prime)$ more "glue vectors" are required to span the expanded code, $\phi_\beta (\mathcal{C})$. These must be inside of $\phi_\beta (\mathcal{C})$ but with a nonzero component outside out $\mathcal{B}_i$ if we want the generator matrix to be full rank. The remaining vectors are therefore a basis of the row space of $\phi_\beta(\mathcal{C}) / \left(\oplus_{i = 1}^m B_i\right)$. The literature often describes the glue vectors as sums of minimum-weight coset leaders of $\mathcal{B}$ considered as polynomials also satisfying the zeros of $\mathcal{C}$ \cite{halford2005soft}; however, the standard coset leaders algorithm makes this difficult to use for even small codes. Instead, these may be computed, even for large codes, using elementary linear algebra using the same algorithm one would for computing the quotient space of two modules. If being used to, for example, connect different hardware modules (corresponding to each $\mathcal{B}$) on a quantum computer, it may be experimentally advantageous to further enforce that the glue vectors be of specific weight or as low-weight as possible. This can be done by selecting appropriate elements from $\phi_\beta(\mathcal{C}) / \left(\oplus_{i = 1}^m B_i\right)$ using `minimum_words` that together have full rank.

This method is highly constrained as the resulting code always has length $m(p^m - 1)$. In a separate 1994 paper \cite{vardy1994maximum}, Vardy and Be'ery showed that binary, primitive BCH codes and binary BCH codes of composite block length may also be put into the form above. We will also refer to this as a Vardy-Be’ery decomposition since the proper technique should be clear from context. For primitive BCH codes, they extended the code and then split the zeros into partitions satisfying certain properties. The direct-sum subcodes are then obtained by puncturing on the set complement of the indices corresponding to the defining sets of each partition. This applies directly to Reed-Muller codes. We will not use this approach here but instead consider BCH codes of composite block length. The two approaches are almost identical except that in the latter case the partitions are immediate from the structure of the code. The following applies to cyclic codes in general.

Let $\mathcal{C}$ be a BCH code over $\mathbb{F}_2$ of composite length $n = n_h n_q$ with defining set $C^n_b \cup \dots \cup C^n_{b + \delta - 2}$. Consider the sets $\mathcal{I}_i = \{1 + j + i \cdot n_h\}$ where $0 \leq j \leq n_h - 1$ for fixed $0 \leq i \leq n_q - 1$.
!!! note "Lemma"
	The code obtained from $\mathcal{C}$ punctured on the complement, $\mathcal{I}_i^c$, is a BCH code of length $n_q$ with defining set $C^{n_q}_b \cup \dots \cup C^{n_q}_{b + \delta - 2}$.

While the second paper does not cite the first paper, the first result may seen as a special case of the second, where the sets $\mathcal{I}_i$ are the non-zero locations of the matrices above.

## Explicit Example
### Example 1: First Result
For an example of the first result, consider the $[7, 4, 4; 0]_8$ Reed-Solomon code.
```
julia> C = ReedSolomonCode(8, 4, 0)
[7, 4, 4; 0]_8 Reed-Solomon code
8-Cyclotomic cosets: 
        C_0 ∪ C_1 ∪ C_2
Generator polynomial:
        x^3 + (α^2 + α + 1)*x^2 + (α^2 + 1)*x + α + 1
Generator matrix: 4 × 7
        α + 1 α^2 + 1 α^2 + α + 1 1 0 0 0
        0 α + 1 α^2 + 1 α^2 + α + 1 1 0 0
        0 0 α + 1 α^2 + 1 α^2 + α + 1 1 0
        0 0 0 α + 1 α^2 + 1 α^2 + α + 1 1
```
To expand
```
julia> F = GF(2)
Galois field with characteristic 2

julia> primitive_basis(field(C), F)
(fqPolyRepFieldElem[1, α, α^2], fqPolyRepFieldElem[1, α^2, α])

julia> β, λ = primitive_basis(field(C), F)
(fqPolyRepFieldElem[1, α, α^2], fqPolyRepFieldElem[1, α^2, α])

julia> C_exp = expanded_code(C, F, β)
[21, 12]_2 linear code
Generator matrix: 12 × 21
        1 1 0 1 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
        0 1 1 1 0 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0
        1 1 1 0 1 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0
        0 0 0 1 1 0 1 0 1 1 1 1 1 0 0 0 0 0 0 0 0
        0 0 0 0 1 1 1 0 0 1 0 1 0 1 0 0 0 0 0 0 0
        0 0 0 1 1 1 0 1 0 1 0 0 0 0 1 0 0 0 0 0 0
        0 0 0 0 0 0 1 1 0 1 0 1 1 1 1 1 0 0 0 0 0
        0 0 0 0 0 0 0 1 1 1 0 0 1 0 1 0 1 0 0 0 0
        0 0 0 0 0 0 1 1 1 0 1 0 1 0 0 0 0 1 0 0 0
        0 0 0 0 0 0 0 0 0 1 1 0 1 0 1 1 1 1 1 0 0
        0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 1 0 1 0 1 0
        0 0 0 0 0 0 0 0 0 1 1 1 0 1 0 1 0 0 0 0 1

julia> function permutation_matrix(F::CodingTheory.CTFieldTypes, n1::Int, n2::Int)
                  # usage: P = permutation_matrix(GF(2), 15, 3) for 3 modules of size 15 each
                  arr = [1 + j + i * n2 for j in 0:(n2 - 1) for i in 0:(n1 - 1)]
                  P = zero_matrix(F, n1 * n2, n1 * n2)
                  F_one = F(1)
                  for i in 1:(n1 * n2)
                      P[arr[i], i] = F_one
                  end
                  return P
              end
permutation_matrix (generic function with 3 methods)

julia> P = permutation_matrix(field(C_exp), 7, 3);

julia> C_exp_P = LinearCode(generator_matrix(C_exp) * P)
[21, 12]_2 linear code
Generator matrix: 12 × 21
        1 1 1 1 0 0 0 1 0 1 0 0 0 0 0 1 1 0 0 0 0
        0 1 1 0 0 0 0 1 0 0 1 0 0 0 1 0 1 0 0 0 0
        1 0 1 0 0 0 0 1 1 0 0 0 0 0 1 0 0 1 0 0 0
        0 1 1 1 1 0 0 0 1 0 1 0 0 0 0 0 1 1 0 0 0
        0 0 1 1 0 0 0 0 1 0 0 1 0 0 0 1 0 1 0 0 0
        0 1 0 1 0 0 0 0 1 1 0 0 0 0 0 1 0 0 1 0 0
        0 0 1 1 1 1 0 0 0 1 0 1 0 0 0 0 0 1 1 0 0
        0 0 0 1 1 0 0 0 0 1 0 0 1 0 0 0 1 0 1 0 0
        0 0 1 0 1 0 0 0 0 1 1 0 0 0 0 0 1 0 0 1 0
        0 0 0 1 1 1 1 0 0 0 1 0 1 0 0 0 0 0 1 1 0
        0 0 0 0 1 1 0 0 0 0 1 0 0 1 0 0 0 1 0 1 0
        0 0 0 1 0 1 0 0 0 0 1 1 0 0 0 0 0 1 0 0 1

julia> B = subfield_subcode(C, F, β)
[7, 3]_2 linear code
Generator matrix: 3 × 7
        1 0 1 1 1 0 0
        1 1 1 0 0 1 0
        0 1 1 1 0 0 1

julia> B_block = LinearCode(generator_matrix(B) ⊕ generator_matrix(B) ⊕ generator_matrix(B))
[21, 9]_2 linear code
Generator matrix: 9 × 21
        1 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 1 0 1 1 1 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 1 1 1 0 0 1 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 1 1 1 0 0 1 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 1 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 1 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 1

julia> Quo = C_exp_P / B_block
[21, 3]_2 linear code
Generator matrix: 3 × 21
        0 1 1 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1 0 0
        1 1 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1 1 1
        1 1 1 0 0 0 1 0 0 0 1 1 1 1 0 0 0 0 1 0 1

julia> C_full = augment(B_block, generator_matrix(Quo))
[21, 12]_2 linear code
Generator matrix: 12 × 21
        1 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 1 0 1 1 1 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 1 1 1 0 0 1 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 1 1 1 0 0 1 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 1 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 1 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 1
        0 1 1 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1 0 0
        1 1 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1 1 1
        1 1 1 0 0 0 1 0 0 0 1 1 1 1 0 0 0 0 1 0 1

julia> are_equivalent(C_exp_P, C_full)
true
```

### Example 2: Second Result
For an example of the second result, consider the BCH code with $b = 0$ and $\delta = 5$ with $n_h = 3$ and $n_q = 15$. The 2-cosets modulo 45 are
```
julia> q = 2; nh = 3; nq = 15; n = nh * nq; b = 0; δ = 5;

julia> all_cyclotomic_cosets(q, n)
8-element Vector{Vector{Int64}}:
 [0]
 [1, 2, 4, 8, 16, 17, 19, 23, 31, 32, 34, 38]
 [3, 6, 12, 24]
 [5, 10, 20, 25, 35, 40]
 [7, 11, 13, 14, 22, 26, 28, 29, 37, 41, 43, 44]
 [9, 18, 27, 36]
 [15, 30]
 [21, 33, 39, 42]

julia> cosets = defining_set([i for i = b:(b + δ - 2)], q, n, false)
3-element Vector{Vector{Int64}}:
 [0]
 [1, 2, 4, 8, 16, 17, 19, 23, 31, 32, 34, 38]
 [3, 6, 12, 24]

julia> C_big = CyclicCode(q, n, cosets)
[45, 28; 0]_2 BCH code
2-Cyclotomic cosets: 
        C_0 ∪ C_1 ∪ C_3
Generator polynomial:
        x^17 + x^16 + x^14 + x^12 + x^8 + x^7 + x^4 + x^3 + x^2 + 1
```
Alternatively, we could have called ``BCHCode(q, n, δ, b)`` directly. Then

$$\mathcal{I}_1 = \{ 1, 4, 7, \dots, 43\} \quad , \quad \mathcal{I}_2 = \{ 2, 5, 8, \dots, 44\} \quad , \quad \mathcal{I}_3 = \{3, 6, 9, \dots, 45\}.$$

The BCH subcode has the defining set $C^{15}_0 \cup \dots \cup C^{15}_3$. The 2-cosets modulo 15 are
```
julia> all_cyclotomic_cosets(q, nq)
5-element Vector{Vector{Int64}}:
 [0]
 [1, 2, 4, 8]
 [3, 6, 9, 12]
 [5, 10]
 [7, 11, 13, 14]
```
so this will be a $[15, 6]$ code.
```
julia> C_small = BCHCode(q, nq, δ, b)
[15, 6; 0]_2 BCH code
2-Cyclotomic cosets: 
        C_0 ∪ C_1 ∪ C_3
Generator polynomial:
        x^9 + x^6 + x^5 + x^4 + x + 1
Generator matrix: 6 × 15
        1 1 0 0 1 1 1 0 0 1 0 0 0 0 0
        0 1 1 0 0 1 1 1 0 0 1 0 0 0 0
        0 0 1 1 0 0 1 1 1 0 0 1 0 0 0
        0 0 0 1 1 0 0 1 1 1 0 0 1 0 0
        0 0 0 0 1 1 0 0 1 1 1 0 0 1 0
        0 0 0 0 0 1 1 0 0 1 1 1 0 0 1
```
Both codes happen to have minimum distance $6 > \delta$.
```
julia> minimum_distance(C_big)
6

julia> minimum_distance(C_small)
6
```
Permuting the indices of $\mathcal{I}_1$ to indices 1 - 15, $\mathcal{I}_2$ to 16 - 30, and $\mathcal{I}_3$ to 31 - 45 completes the direct-sum subcode.
```
julia> P = permutation_matrix(field(C_big), nq, nh);

julia> C_big_P = LinearCode(generator_matrix(C_big) * P);

julia> C_block = LinearCode(generator_matrix(C_small) ⊕ generator_matrix(C_small) ⊕ generator_matrix(C_small));

julia> QC = C_big_P / C_block;

julia> C_P_full = augment(C_block, generator_matrix(QC));

julia> are_equivalent(C_P_full, C_big_P)
true
```

It is also possible to work these theorems backwards by choosing a desired subcode then seeing if a supercode can be built to contain it.
