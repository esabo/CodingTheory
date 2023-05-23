# Quantum Reed-Muller Codes
A common way to circumvent the Eastin-Knill theorem preventing any single code from having a universal gate set is to switch between two codes which together have a universal gate set. The canonical example of this are the Steane code, which supports the transversal Clifford group, and the 15-qubit quantum Reed-Muller code, which supports a transversal $T$ and/or $CCZ$.

There are multiple ways to define a family of quantum codes from the classical Reed-Muller family, and the literature is evenly split between the various possibilities. The quantum Reed-Muller code family we are interested in here is derived from shortened Reed-Muller codes and the CSS construction. We will follow [1] and denote these by $QRM(m)$. The Steane code is $QRM(3)$ and the $[[15, 1, 3]]$ code is $QRM(4)$. Let's see how we can generate these codes and study the code switching in the library.

The 15-qubit QRM code is formed from the shortened Reed-Muller codes $\overline{\mathcal{RM}}(1, 4)$ and $\overline{\mathcal{RM}}(2, 4)$ with $X$ stabilizers given by $\overline{G}(1, 4)$ and $Z$ stabilizers by $\overline{G}(2, 4)$. Recalling the relationship between shortened and punctured codes, the $X$ stabilizers are equivalent to a parity check matrix for $\mathcal{RM}^*(2, 4)$ and the $Z$ stabilizers are equivalent to a parity check matrix for $\mathcal{RM}^*(1, 4)$. It's instructive to construct these explicitly.

The easiest way to generate $QRM(3)$ is
```
m = 3
RM13 = ReedMullerCode(2, 1, m)
RM13s = shorten(RM13, 1)
QRM3 = CSSCode(RM13s)
```
We can check that $\overline{\mathcal{RM}}(1, 3)$ satisfies the requirements of the single-code CSS construction
```
isselforthogonal(RM13s)
```
Theory tells us that the shortened codes are given by deleting the first row and column of the $(u \mid u + v)$ form of the generator matrix of the Reed-Muller codes. Unfortunately, the `shorten` function obscures this fact when it computes a kernel, but we can test whether or not this is true.
```
RM13salt = LinearCode(generatormatrix(RM13)[2:end, 2:end])
areequivalent(RM13s, RM13salt)
```
We want the stabilizers of our code in this alternate form, so we can either remake `QRM3` using `RM13salt`, use the other `CSSCode` constructor where we explicitly pass in the $X$ and $Z$ stabilizer matrices, or replace the stabilizers of the already constructed object `QRM3`. This last option automatically checks that the old stabilizers and the new stabilizers have equivalent row spaces and errors if they don't.
```
setXstabilizers!(QRM3, generatormatrix(RM13salt))
setZstabilizers!(QRM3, generatormatrix(RM13salt))
```

Similarly, we know the logicals of this code is the all-ones vector and can use this form if desired.
```
logicals(QRM3)
F = field(QRM3)
newlogs = zero_matrix(F, 2, 2 * length(QRM3))
for i in 1:length(QRM3)
    newlogs[1, i] = F(1)
    newlogs[2, i + length(QRM3)] = F(1)
end
setlogicals!(QRM3, newlogs)
```
As before, this will automatically check if the input is equivalent to the automatically computed logicals up to stabilizers.

In general, the $X$ stabilizers of $QRM(m)$ are given by the generators of $\overline{\mathcal{RM}}(1, m)$ and the $Z$ stabilizers are given by the generators of $\overline{\mathcal{RM}}(m - 2,  m)$, producing the parameters $$[[2^m − 1, 1, 3]]$$. For $QRM(4)$, we have,
```
m = 4
RM14 = ReedMullerCode(2, 1, m)
RM24 = ReedMullerCode(2, m - 2, m)
RM14s = shorten(RM14, 1)
RM24s = shorten(RM24, 1)
RM14salt = LinearCode(generatormatrix(RM14)[2:end, 2:end])
RM24salt = LinearCode(generatormatrix(RM24)[2:end, 2:end])
areequivalent(RM14s, RM14salt)
areequivalent(RM24s, RM24salt)
```

In this library, we choose the convention that `C2 ⊆ C1` for the `CSSCode(C1, C2)`. Thus, to make our code, we actually require `CSSCode(dual(RM24s), RM14s)`, and we can check that `RM14 ⊆ dual(RM24)`. We can do this instead with the alternative form of the generator matrix or repeat what we did above, but instead let's use the other constructor
```
QRM4 = CSSCode(generatormatrix(RM14salt), generatormatrix(RM24salt))
```
One may compare these stabilizers to the built-in commands `SteaneCode()` and `Q15RM()` and also against the explicit set of stabilizers listed in [2].

## Viewing $QRM(m)$ As Subsystem Codes
It was long known that the Steane code is contained in the 15-qubit Reed-Muller code, but [1] extended this idea to show that this is not only true for $QRM(m)$ and $QRM(m + 1)$ but also that this can be viewed as gauge fixes of a single subsystem code. To understand this, consider the generator matrices of the Reed-Muller family. They are constructed recursively via

$$G(r, m) = \begin{pmatrix}
    G(r, m - 1) & G(r, m - 1)\\
    0 & G(r - 1, m - 1)
\end{pmatrix},$$

with the base case that

$$G(1, 1) = \begin{pmatrix} 1 & 1\\ 0 & 1 \end{pmatrix},$$

$G(m, m)$ is the identity otherwise, and $G(0, m)$ is the length $2^m$ all-ones vector. Thus,

$$G(1, 2) = \begin{pmatrix}
    1 & 1 & 1 & 1\\
    0 & 1 & 0 & 1\\
    0 & 0 & 1 & 1
\end{pmatrix}$$

$$G(1, 3) = \begin{pmatrix}
    1 & 1 & 1 & 1 & 1 & 1 & 1 & 1\\
    0 & 1 & 0 & 1 & 0 & 1 & 0 & 1\\
    0 & 0 & 1 & 1 & 0 & 0 & 1 & 1\\
    0 & 0 & 0 & 0 & 1 & 1 & 1 & 1
\end{pmatrix}$$

$$G(1, 4) = \begin{pmatrix}
    1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1\\
    0 & 1 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 1\\
    0 & 0 & 1 & 1 & 0 & 0 & 1 & 1 & 0 & 0 & 1 & 1 & 0 & 0 & 1 & 1\\
    0 & 0 & 0 & 0 & 1 & 1 & 1 & 1 & 0 & 0 & 0 & 0 & 1 & 1 & 1 & 1\\
    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1
\end{pmatrix}$$

The generator matrices of the shortened codes are therefore

$$\overline{G}(1, 3) = \begin{pmatrix}
    1 & 0 & 1 & 0 & 1 & 0 & 1\\
    0 & 1 & 1 & 0 & 0 & 1 & 1\\
    0 & 0 & 0 & 1 & 1 & 1 & 1
\end{pmatrix}$$

and

$$\overline{G}(1, 4) = \begin{pmatrix}
    1 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 1\\
    0 & 1 & 1 & 0 & 0 & 1 & 1 & 0 & 0 & 1 & 1 & 0 & 0 & 1 & 1\\
    0 & 0 & 0 & 1 & 1 & 1 & 1 & 0 & 0 & 0 & 0 & 1 & 1 & 1 & 1\\
    0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1
\end{pmatrix}$$

Notice that the first three rows of $\overline{G}(1, 4)$ are of the form $(\overline{G}(1, 3) \mid 0 \mid \overline{G}(1, 3))$, where the notation $( \, \mid \, )$ denotes horizontal concatenation. In this sense we see that the 15-qubit Reed-Muller code really contains *two* copies of the Steane code. We can see this in the $X$ stabilizers of `QRM4`,
```
Xstabilizers(QRM4)[1:3, :] == hcat(generatormatrix(RM13salt), zero_matrix(F, 3, 1), generatormatrix(RM13salt))
```
It's less clear that the $Z$ stabilizers also contain the two copies of the Steane code in this sense. To see this, let's first define a new stabilizer code whose $X$ and $Z$ stabilizers are of this form.
```
test = CSSCode(Xstabilizers(QRM4)[1:3, :], Xstabilizers(QRM4)[1:3, :])
```
Now we can remove these stabilizers from `QRM4`,
```
quo1 = CodingTheory._quotientspace(stabilizers(QRM4), stabilizers(test))
```
Let's set the stabilizers of `QRM4` to make this more explicit.
```
setstabilizers!(QRM4, vcat(stabilizers(test), quo1))
```



In order for the information to not be disturbed...

```

L = vcat(hcat(logicalsmatrix(QRM3)[1, :], zero_matrix(F, 1, length(QRM4) + 1)),
	hcat(zero_matrix(F, 1, length(QRM4)), logicalsmatrix(QRM3)[1, :], zero_matrix(F, 1, 1)))
CodingTheory._quotientspace(logicalsmatrix(test), L)
```


[1]: Anderson, Duclos-Cianci, Poulin, "Fault-tolerant conversion between the Steane and Reed-Muller quantum codes", (2014)

[2]: Chamberlin paper


showed that this can be viewed as a different gauge fixings of a single subsystem code. Let's see how 

what are the elements in Z outside of X


Steane code is shortened RM(1, 3), which is [7, 4, 3] the Hamming code
the 15QRM code has two copies of the Steane code plus extra
demo this and setup the code in this fashion (hcat)

intersect stabilizers of QRM(m + 1) and QRM(m) to find common subset
need to build this in a way that the logical doesn't move







RM13 = ReedMullerCode(2, 1, 3)
RM14 = ReedMullerCode(2, 1, 4)
RM14s = shorten(C2, 1)
RM14salt = LinearCode(generatormatrix(RM14)[2:end, 2:end])
areequivalent(RM14s, RM14salt)
RM24 = ReedMullerCode(2, 2, 4)
RM24s = shorten(RM24, 1)
RM24salt = LinearCode(generatormatrix(RM24)[2:end, 2:end])
areequivalent(RM24s, RM24salt)
However, it is convenient to use the explicit form of the generator matrices in `RM14salt` and `RM24salt`. We can do this by either passing explicit `X`- and `Z`-stabilizer matrices into the constructor directly via `CSSCode(generatormatrix(RM24salt), generatormatrix(RM14salt))` or using the command `setstabilizers!` to change the form of the stabilizers of an already existing code. The latter automatically checks that 
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



## Other Representations

We also know that the Steane code may be constructed via the [7, 4, 3] Hamming code. In the C \subseteq C^\perp versus C^\perp \subseteq C convention used in the library, the Steane code may be derived from the dual of this code, which is called the simplex code.
```
D = dual(HammingCode(2, 3))
SteaneHamming = CSSCode(D)
```
The command `SimplexCode(2, 3)` also would have worked.

Upon immediate inspection, `QRM3` and `SteaneHamming` are not equivalent. In fact, neither are equivalent to the built-in `SteaneCode()`. Let us show that these are all equivalent up to permutation of the qubits.
    (implement permutation on the quantum side and demo here)





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
