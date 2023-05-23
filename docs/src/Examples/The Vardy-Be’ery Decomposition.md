# The Vardy-Be’ery Decomposition

## Background
Before getting into the main application, let us briefly review the necessary background from classical coding theory.

Let $\mathcal{C}$ be an $[n, k, d]_{p^m}$ code. Then the subfield subcode of $\mathcal{C}$ over a subfield $\mathbb{F} \leq \mathbb{F}_{p^m}$, denoted $\mathcal{C}|_\mathbb{F}$, is given by $\mathcal{C} \cap \mathbb{F}^n$, i.e., the collection of codewords of $\mathcal{C}$ whose components lie entirely in $\mathbb{F}$. The code $\mathcal{C}$ is called the supercode of $\mathcal{C}|_\mathbb{F}$. If $\mathcal{C}$ has parameters $[n, k, d]_{p^m}$, $\mathcal{C}|_\mathbb{F}$ has parameters $[n, k^\prime, \geq d]$ over $\mathbb{F}$, where $n - k \leq n - k^\prime \leq \ell (n - k)$ and $\ell = [\mathbb{F}_{p^m} : \mathbb{F}]$ (the index of $\mathbb{F}$ in $\mathbb{F}_{p^m}$). As the codewords of $\mathcal{C}|_\mathbb{F}$ are codewords of $\mathcal{C}$, it follows immediately that the minimum distance of $\mathcal{C}|_\mathbb{F}$ is at least the minimum distance of $\mathcal{C}$, and $\mathcal{C}|_\mathbb{F}$ can be decoded using the same algorithm as $\mathcal{C}$, although perhaps not efficiently as a native algorithm over $\mathbb{F}$ designed specifically for the subfield subcode.

If $\mathcal{C}$ is an $[n, k, d]_{p^m}$ Reed-Solomon code, then $\mathcal{C}|_{\mathbb{F}_p}$ is the BCH code over $\mathbb{F}_p$ of length $n$ and designed distance $d$. The proof of this follows immediately from the fact that the codewords of the BCH code are elements of $\mathbb{F}_p^n$ and the zero set of the Reed-Solomon code is a subset of the zero set of the BCH code.

An $[mn, mk, \geq d]_p$ code may be constructed from an $[n, k, d]_{p^m}$ code by expanding its elements using a basis of $\mathbb{F}_{p^m}/\mathbb{F}_p$. The first code is called the expanded code of the second. To see why the minimum distance of the code could increase, let $\beta = \{\beta_j\}_1^m$ be a basis of $\mathbb{F}_{p^m}/\mathbb{F}_p$ and let $c = (c_1, \dots, c_n) \in \mathbb{F}^n_{p^m}$ be a minimum weight codeword in an $[n, k, d]_{p^m}$ code. Expressing each $c_i$ with respect to $\beta$, $c_i = \sum_j c_{ij} \beta_j$, we can replace each element with its corresponding $m$-tuple, $(c_{i1}, \dots, c_{im})$. If  $c_i \neq 0$, then the Hamming weight of its expansion is at least one and therefore the Hamming weight of the expansion of $c$ is at least $d$.

Recall that the inner product over finite fields is given by the trace. In particular, if $\beta$ is a basis of $\mathbb{F}_{p^m}/\mathbb{F}_p$ such that $x = \sum_{j = 1}^m x_j \beta_j$ for $x \in \mathbb{F}_{p^m}$, then $x_j = \mathrm{Tr}_{\mathbb{F}_{p^m}/\mathbb{F}_p} (x \beta_j^\perp) \in \mathbb{F}_p$, where $\beta^\perp$ is the unique trace-orthogonal dual of $\beta$ such that $\mathrm{Tr}_{\mathbb{F}_{p^m}/\mathbb{F}_p}(x_i y_j) = \delta_{ij}$ for $x_i \in \beta$ and $y_j \in \beta^\perp$. (The dual basis always exists and is easy to compute given $\beta$.) For $c = (c_1, \dots, c_n) \in \mathbb{F}_{p^m}^n$ denote the expansion with respect to $\beta$ by the isomorphism $\phi_\beta: \mathbb{F}_{p^m}^n \to \mathbb{F}_p^{nm}$ given by

$$\begin{aligned}
	\phi_\beta(c) &= (\phi_\beta(c_1), \dots, \phi_\beta(c_n))\\
		&= (\mathrm{Tr}_{\mathbb{F}_{p^m}/\mathbb{F}_p}(c_1 \beta_1^\perp), \dots, \mathrm{Tr}_{\mathbb{F}_{p^m}/\mathbb{F}_p}(c_1 \beta_m^\perp), \mathrm{Tr}_{\mathbb{F}_{p^m}/\mathbb{F}_p}(c_2 \beta_1^\perp), \dots, \mathrm{Tr}_{\mathbb{F}_{p^m}/\mathbb{F}_p}(c_n \beta_m^\perp)).
\end{aligned}$$

If $\mathcal{C}$ is a code, let $\phi_\beta(C)$ denote the corresponding expanded code. In general, an expanded code loses the properties of its parent code and different bases could produce different expanded codes with different parameters and properties. It is still not yet known how to choose a basis to *a priori* maximize the minimum distance of the expanded code.

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

If $\mathcal{C}$ has dimension $k$ and $\mathcal{B}$, $k^\prime$, then $m(k - k^\prime)$ more "glue vectors" are required to span the expanded code, $\phi_\beta (\mathcal{C})$. These must be inside of $\phi_\beta (\mathcal{C})$ but with a nonzero component outside out $\mathcal{B}_i$ if we want the generator matrix to be full rank. The remaining vectors are therefore a basis of the row space of $\phi_\beta(\mathcal{C}) / \left(\oplus_{i = 1}^m B_i\right)$. The literature often describes the glue vectors as sums of minimum-weight coset leaders of $\mathcal{B}$ considered as polynomials also satisfying the zeros of $\mathcal{C}$ \cite{halford2005soft}; however, the standard coset leaders algorithm makes this difficult to use for even small codes. Instead, these may be computed, even for large codes, using elementary linear algebra using the same algorithm one would for computing the quotient space of two modules. If being used to, for example, connect different hardware modules (corresponding to each $\mathcal{B}$) on a quantum computer, it may be experimentally advantageous to further enforce that the glue vectors be of specific weight or as low-weight as possible. This can be done by selecting appropriate elements from $\phi_\beta(\mathcal{C}) / \left(\oplus_{i = 1}^m B_i\right)$ using `minimumwords` that together have full rank.

This method is highly constrained as the resulting code always has length $m(p^m - 1)$. In a separate 1994 paper \cite{vardy1994maximum}, Vardy and Be'ery showed that binary, primitive BCH codes and binary BCH codes of composite block length may also be put into the form above. We will also refer to this as a Vardy-Be’ery decomposition since the proper technique should be clear from context. For primitive BCH codes, they extended the code and then split the zeros into partitions satisfying certain properties. The direct-sum subcodes are then obtained by puncturing on the set complement of the indices corresponding to the defining sets of each partition. This applies directly to Reed-Muller codes. We will not use this approach here but instead consider BCH codes of composite block length. The two approaches are almost identical except that in the latter case the partitions are immediate from the structure of the code. The following applies to cyclic codes in general.

Let $\mathcal{C}$ be a BCH code over $\mathbb{F}_2$ of composite length $n = n_h n_q$ with defining set $C^n_b \cup \dots \cup C^n_{b + \delta - 2}$. Consider the sets $\mathcal{I}_i = \{1 + j + i \cdot n_h\}$ where $0 \leq j \leq n_h - 1$ for fixed $0 \leq i \leq n_q - 1$.
!!! note "Lemma"
	The code obtained from $\mathcal{C}$ punctured on the complement, $\mathcal{I}_i^c$, is a BCH code of length $n_q$ with defining set $C^{n_q}_b \cup \dots \cup C^{n_q}_{b + \delta - 2}$.

While the second paper does not cite the first paper, the first result may seen as a special case of the second, where the sets $\mathcal{I}_i$ are the non-zero locations of the matrices above.

## Explicit Example
TODO: work everything out by hand here and come up with an example for both decomps






To demo this, we first define the permutation matrix required to put the columns into the almost-block-diagonal form above.
```
julia> function permutationmatrix(F::CodingTheory.CTFieldTypes, n1::Int, n2::Int)
           # usage: P = permutationmatrix(GF(2), 15, 3) for 3 modules of size 15 each
           arr = [1 + j + i * n2 for j in 0:(n2 - 1) for i in 0:(n1 - 1)]
           P = zero_matrix(F, n1 * n2, n1 * n2)
           Fone = F(1)
           for i in 1:(n1 * n2)
               P[arr[i], i] = Fone
           end
           return P
       end
permutationmatrix (generic function with 1 method)
```

```
q = 2; n = 45; nq = 15; nh = 3; bX = 0; δX = 7;
CBigX = BCHCode(q, n, δX, bX);
CSmX = BCHCode(q, nq, δX, bX);
PX = permutationmatrix(field(CBigX), nq, nh);
CBigXP = LinearCode(generatormatrix(CBigX) * PX);
CblockX = LinearCode(generatormatrix(CSmX) ⊕ generatormatrix(CSmX) ⊕ generatormatrix(CSmX));
QCX = CBigXP / CblockX;
CPfullX = augment(CblockX, generatormatrix(QCX));
areequivalent(CPfullX, CBigXP)

```
