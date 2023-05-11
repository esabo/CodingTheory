A common way to circumvent the Eastin-Knill theorem preventing any single code from having a universal gate set is to switch between two codes which together have a universal gate set. The canonical example of this are the Steane code, which supports the transversal Clifford group, and the 15-qubit quantum Reed-Muller code, which supports a transversal T and/or CCZ.

There are multiple ways to define a family of quantum codes from the classical Reed-Muller family, and the literature is evenly split between the various possibilities. The quantum Reed-Muller code family we are interested in here derived from shortened Reed-Muller codes and the CSS construction. We will follow the reference and denote these by QRM(m). The Steane code is QRM(3) and the [[15, 1, 3]] code is QRM(4). Let's see how we can generate these codes and study the code switching in the library.

The 15-qubit QRM code is formed from the shortened Reed-Muller codes $\overline{\mathcal{RM}}(1, 4)$ and $\overline{\mathcal{RM}}(2, 4)$ with $X$ stabilizers given by $\overline{G}(1, 4)$ and $Z$ stabilizers by $\overline{G}(2, 4)$. Recalling the relationship between shortened and punctured codes, the $X$ stabilizers are equivalent to a parity check matrix for $\mathcal{RM}^*(2, 4)$ and the $Z$ stabilizers are equivalent to a parity check matrix for $\mathcal{RM}^*(1, 4)$. It's instructive to construct these explicitly.

The easiest way to generate QRM(3) is
```
m = 3
RM13 = shorten(ReedMullerCode(2, 1, m), 1)
QRM3 = CSSCode(RM13)
```
where we have shortened the codes by deleting the first row and column of the (u | u + v) form of the generator matrix of the Reed-Muller codes RM(r, m). We can check that \overline(RM(1, 3)) satisfies the requirements of the single-code CSS construction
```
isselforthogonal(RM13)
```
We also know that the Steane code may be constructed via the [7, 4, 3] Hamming code. In the C \subseteq C^\perp versus C^\perp \subseteq C convention used in the library, the Steane code may be derived from the dual of this code, which is called the simplex code.
```
D = dual(HammingCode(2, 3))
SteaneHamming = CSSCode(D)
```
The command `SimplexCode(2, 3)` also would have worked.

Upon immediate inspection, `QRM3` and `SteaneHamming` are not equivalent. In fact, neither are equivalent to the built-in `SteaneCode()`. Let us show that these are all equivalent up to permutation of the qubits.
    (implement permutation on the quantum side and demo here)

Similar to the above, we can generate QRM(4) via
```
m = 4
RM14 = shorten(ReedMullerCode(2, 1, m), 1)
RM24 = shorten(ReedMullerCode(2, m - 2, m), 1)
```
In general, the X stabilizers of QRM(m) are given by the generators of \overline{RM(1, m)} and the Z stabilizers are given by the generators of \overline(RM(m - 2,  m)). The CSS construction uses a parity-check matrix, so
```
CSSCode(RM24, RM14)
```
is actually a [[15, 6]] code. The code we are looking for therefore requires either using the constructor where we pass in the stabilizer matrices directly or taking the dual. We can check to see if taking the dual works for the CSS construction,
```
RM14 ⊆ dual(RM24)
```
Finally,
```
QRM4 = CSSCode(dual(RM24), RM14)
```
gives us the [[15, 1, 3]] code we are looking for.

The built-in command `Q15RM()` returns a different looking stabilizer code whose form is hardcoded from the stabilizers listed in [Chamberlin].

It is easy to analyze subsystem codes which arise from classical cyclic codes. Using a well-known relationship between Reed-Muller and BCH codes from classical coding theory, we can see that the $QRM(m)$ family is indeed cyclic. Let us show this for $QRM(3)$ and $QRM(4)$. Let $n = p^m - 1$. The $p$-weight of an integer $0 \leq a \leq n$ is $\mathrm{wt}_p(a) = \sum_{j = 0}^{m - 1} a_j$, where $a = \sum_{j = 0}^{m - 1} a_j p^j$, $0 \leq a_j \leq p - 1$ is the $p$-adic expansion of $a$. Equivalently, we may interpret a vector in $\F_p^m$ as the coefficients of a $p$-adic expansion and define the $p$-weight as the sum of the elements. Consider monomials of the form $x_1^{i_1} x_2^{i_2}\hdots x_m^{i_m}$ for  $i_1 + \hdots + i_m \leq r$. Interpreting $(i_1, \hdots, i_m)$ as a $p$-adic expansion, all cyclic shifts are also valid monomials of total degree less than $r$, have constant $p$-weight, and generate the $p$-coset $C_i$ where $i = \sum_{j = 0}^{m - 1} i_j p^j$. In this way we establish a correspondence between multivariate polynomials of $\mathrm{RM}^*(r, m)$ and univariate polynomials of BCH codes.
\begin{theorem}[\cite{kasami1968new}]
	Let $\alpha$ be a primitive root of $\F_{p^m}^\times$ and define $g^*_{r, m}(x) = \prod (x - \alpha^a)$ where $0 < \mathrm{wt}_p(a) \leq m(p - 1) - r - 1$. Then $\mathcal{RM}_{p^m}^*(r, m)$ is permutation equivalent to the subfield subcode of $\mathcal{C}^*_{r, m} = (g^*_{r, m}(x))$ over $\F_p$.
\end{theorem}
\noindent Since the defining set of $\mathcal{C}^*_{r, m}$ is comprised of complete $p$-cosets, $g^*_{r, m} \in \F_p[x]$ and hence $\mathcal{C}^*_{r, m} = \mathcal{C}^*_{r, m} \cap \F_p$.

We begin with $\mathcal{RM}^*(1, 4)$. The set of all integers $a$ with Hamming weight $0 < \mathrm{wt}_2(a) \leq 2$ is
```
a = Vector{Int}()
for i in 1:15
    sum(digits(i, base=2)) <= 2 && push!(a, i)
end
b = sort(cyclotomiccoset(1, 2, 15) ∪ cyclotomiccoset(3, 2, 15) ∪ cyclotomiccoset(5, 2, 15))
a == b
```
\begin{equation*}
	\{1, 2, 3, 4, 5, 6, 8, 9, 10, 12\} = C^{15}_1 \cup C^{15}_3 \cup C^{15}_5.
\end{equation*}
The corresponding generator polynomial for the cyclic code is $g^*_{1, 4}(x) = 1 + x + x^2 + x^4 + x^5 + x^8 + x^{10}$. While this is a binary code, it must be built over the extension field $\F_{16}$ for the root of unity.
```
C1 = CyclicCode(16, 15, [[1], [2], [3], [4], [5], [6], [8], [9], [10], [12]])
```
For $\mathcal{RM}^*(2, 4)$, the set of all integers $a$ with Hamming weight $0 < \mathrm{wt}_2(a) \leq 1$ is $\{1, 2, 4, 8\} = C^{15}_1$. The corresponding generator polynomial is $g^*_{2, 4}(x) = 1 + x^2 + x^3 + x^4$.
```
C2 = CyclicCode(16, 15, [[1], [2], [4], [8]])
```
The dual codes have generator polynomials $(g^*_{1, 4})^\perp(x) = 1 + x^2 + x^4 + x^5$ and $(g^*_{2, 4})^\perp(x) = 1 + x^3 + x^4 + x^6 + x^8 + x^9 + x^{10} + x^{11}$. Generator matrices for these are
\begin{equation}
	\overline{G}(1, 4) = \begin{pmatrix}
		1 & 0 & 0 & 1 & 1 & 0 & 1 & 0 & 1 & 1 & 1 & 1 & 0 & 0 & 0\\
		0 & 1 & 0 & 0 & 1 & 1 & 0 & 1 & 0 & 1 & 1 & 1 & 1 & 0 & 0\\
		0 & 0 & 1 & 0 & 0 & 1 & 1 & 0 & 1 & 0 & 1 & 1 & 1 & 1 & 0\\
		0 & 0 & 0 & 1 & 0 & 0 & 1 & 1 & 0 & 1 & 0 & 1 & 1 & 1 & 1
	\end{pmatrix},
\end{equation}
and
\begin{equation}
	\overline{G}(2, 4) = \begin{pmatrix}
		1 & 0 & 1 & 0 & 1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
		0 & 1 & 0 & 1 & 0 & 1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
		0 & 0 & 1 & 0 & 1 & 0 & 1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
		0 & 0 & 0 & 1 & 0 & 1 & 0 & 1 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\
		0 & 0 & 0 & 0 & 1 & 0 & 1 & 0 & 1 & 1 & 0 & 0 & 0 & 0 & 0\\
		0 & 0 & 0 & 0 & 0 & 1 & 0 & 1 & 0 & 1 & 1 & 0 & 0 & 0 & 0\\
		0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 1 & 0 & 1 & 1 & 0 & 0 & 0\\
		0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 1 & 0 & 1 & 1 & 0 & 0\\
		0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 1 & 0 & 1 & 1 & 0\\
		0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 1 & 0 & 1 & 1
	\end{pmatrix}.
\end{equation}
A set of explicit stabilizers for the 15-qubit \gls{qrm} code are given in \cite{chamberland2017error} as
\begin{equation}\label{RMX}
	G_X = \begin{pmatrix}
		1 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 1\\
   		 0 & 1 & 1 & 0 & 0 & 1 & 1 & 0 & 0 & 1 & 1 & 0 & 0 & 1 & 1\\
   		 0 & 0 & 0 & 1 & 1 & 1 & 1 & 0 & 0 & 0 & 0 & 1 & 1 & 1 & 1\\
   		 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1
	\end{pmatrix},
\end{equation}
and
\begin{equation}\label{RMZ}
	G_Z = \begin{pmatrix}
		1 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 1\\
		0 & 1 & 1 & 0 & 0 & 1 & 1 & 0 & 0 & 1 & 1 & 0 & 0 & 1 & 1\\
		0 & 0 & 0 & 1 & 1 & 1 & 1 & 0 & 0 & 0 & 0 & 1 & 1 & 1 & 1\\
		0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1\\
		0 & 0 & 1 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 1\\
		0 & 0 & 0 & 0 & 1 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 1\\
		0 & 0 & 0 & 0 & 0 & 1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 1\\
		0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 0 & 0 & 1 & 1\\
		0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 1 & 1\\
		0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 1
	\end{pmatrix}.
\end{equation}
One may check that $\mathrm{rowspace}(\overline{G}(1, 4)) \cong \mathrm{rowspace}(G_X)$ and $\mathrm{rowspace}(\overline{G}(2, 4)) \cong$\\
$\mathrm{rowspace}(G_Z)$ via the permutation $(3 \,\, 5 \,\, 9 \,\, 15 \,\, 13 \,\, 14 \,\, 12 \,\, 7 \,\, 11 \,\, 8 \,\, 4)$.






first discussed by [] and then by many authors. 

[Anderson, Duclos-Cianci, Poulin, "Fault-tolerant conversion between the Steane and Reed-Muller quantum codes", (2014)] showed that this can be viewed as a different gauge fixings of a single subsystem code. Let's see how 

what are the elements in Z outside of X
[[2^m − 1, 1, 3]]
logicals are all 1's

Steane code is shortened RM(1, 3), which is [7, 4, 3] the Hamming code
the 15QRM code has two copies of the Steane code plus extra
demo this and setup the code in this fashion (hcat)

intersect stabilizers of QRM(m + 1) and QRM(m) to find common subset
need to build this in a way that the logical doesn't move







RM13 = ReedMullerCode(2, 1, 3)
RM14 = ReedMullerCode(2, 1, 4)
RM14s = shorten(C2, 1)
RM14salt = LinearCode(generatormatrix(RM14s)[2:end, 2:end])
areequivalent(RM14s, RM14salt)
RM24 = ReedMullerCode(2, 2, 4)
RM24s = shorten(RM24, 1)
RM24salt = LinearCode(generatormatrix(RM24s)[2:end, 2:end])
areequivalent(RM24s, RM24salt)
However, it is convenient to use the explicit form of the generator matrices in `RM14salt` and `RM24salt`. We can do this by either passing explicit `X`- and `Z`-stabilizer matrices into the constructor directly via `CSSCode(generatormatrix(RM24salt), generatormatrix(RM14salt))` or using the command `setstabilizers!` to change the form of the stabilizers of an already existing code. The latter automatically checks that the old stabilizers and the new stabilizers have equivalent row spaces and errors if they don't.
```
setXstabilizers!(S, generatormatrix(RM14salt))
setZstabilizers!(S, generatormatrix(RM24salt))
S
```




julia> ReedMullerCode(2, 2, 4) / ReedMullerCode(2, 1, 4)
[16, 6]_2 linear code
Generator matrix: 6 × 16
        0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1
        0 0 0 0 0 1 0 1 0 0 0 0 0 1 0 1
        1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1
        0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 0
        1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0
        1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0

