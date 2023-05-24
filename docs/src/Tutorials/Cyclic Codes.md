# Cyclic Codes

## Background
The following represents the notation and conventions used for cyclic codes throughout the library.

A code is cyclic if for every $(f_0, \dots, f_{n - 1}) \in \mathcal{C}$ the vector $(f_{n - 1}, f_0, \dots, f_{n - 2})$ is also in $C$. Consider the polynomial $f = f_0 + f_1 x + \dots + f_{n - 1} x^{n - 1}$. Multiplying by $x$ and setting $x^n = 1$ gives $f_{n - 1} + f_0 x + \dots + f_{n - 2} x^{n - 1}$. Thus, elements of cyclic codes are naturally viewed as coefficient vectors of polynomials in $\mathbb{F}_p[x]/(x^n - 1)$. Cyclic codes are in bijection with ideals of this ring and hence with divisors of $x^n - 1$.

Let $x^n - 1 = g(x)h(x)$ for some $g(x), h(x) \in \mathbb{F}_p[x]$. Then $\mathcal{C} = (g(x))$ is a cyclic code with generator polynomial $g(x)$ viewed as an ideal of $\mathbb{F}_p[x]/(x^n - 1)$. Let $c(x) = a(x) g(x) \in \mathcal{C}$; then $c(x) h(x) = a(x) g(x) h(x) = a(x) (x^n - 1) \equiv 0$. In analogy with $H$, $h(x)$ is called the parity check polynomial. The code $\mathcal{C}$ has parameters $[n, n - \mathrm{deg}(g)]$. Recall that the reciprocal (reverse) of a polynomial $f(x)$ of degree $n$ is $f^r(x) = x^n f(x^{-1})$. An easy argument shows that $\mathcal{C}^\perp = (h^r(x)/h(0))$, where we have introduced a normalization factor to keep the generator polynomial monic. Note that if $h(x) \mid x^n - 1$, then $h^r(x)/h(0) \mid x^n - 1$, so the dual code of a cyclic code is cyclic. Let $\mathrm{deg}(g) = n - k$. Then the generator and parity check matrices for $\mathcal{C}$ are given by

$$G = \begin{pmatrix}
		g_0 & g_1 & \dots & \dots & \dots & \dots & g_{n - k}& & & \\
		 & g_0 & g_1 & \dots & \dots & \dots & \dots & g_{n - k} & & \\
		 & & \ddots & \ddots & \ddots & \ddots & \ddots & \ddots & \ddots & \\
		 & & & g_0 & g_1 & \dots & \dots & \dots & \dots & g_{n - k}
	\end{pmatrix},$$

and

$$H = \begin{pmatrix}
		h_k & h_{k - 1} & \dots & \dots & \dots & \dots & h_0 & & & \\
		 & h_k & h_{k - 1} & \dots & \dots & \dots & \dots & h_0 & & \\
		 &  & \ddots & \ddots & \ddots & \ddots & \ddots & \ddots & \ddots & & \\
		 & & & h_k & h_{k - 1} & \dots & \dots & \dots & \dots & h_0
	\end{pmatrix},$$

respectively.

The general idea to factoring $x^n - 1$ is always that over the splitting field, $\displaystyle x^n - 1 = \prod_{\alpha^n = 1} (x - \alpha)$, where the product is taken over all $n$-th roots of unity, not necessarily primitive. While not over the splitting field, some of these terms need to be grouped together into irreducible factors, $\displaystyle x^n - 1 = \prod \mathrm{min}_\alpha(x)$, where $\mathrm{min}_\alpha(x)$ is the minimal polynomial for $\alpha$ over the appropriate base field. It follows from the binomial theorem that for $f(x) \in \mathbb{F}_p[x]$, $f(x^p) = f(x)^p$. Hence, for $\alpha$ a root of $f(x)$ in some extension field of $\mathbb{F}_p$, $f(\alpha^p) = f(\alpha)^p = 0$, implying $\alpha, \alpha^p, \alpha^{p^2}, \dots$ are all roots of $f(x)$. This sequence stops when $\alpha^{p^r} = \alpha$ for some natural number $r$. Let $E$ be the splitting field of $x^n - 1$ with $\mathrm{gcd}(n, p) = 1$, and let $\beta$ be a primitive element of $E$. Then $\alpha = \beta^d$ is a primitive $n$-th root of unity with $d = (|E| - 1)/n$. Then $\alpha^{p^r} = \alpha \to \beta^{dp^r - d} = 1$, or $dp^r \equiv d \mod (|E| - 1)$. Note that $|E| = p^{\mathrm{ord}_n(p)}$, where $\mathrm{ord}_n(p)$ is the smallest positive integer $m$ such that $p^m \equiv 1\mod n$.

We can collect this sequence of roots into the $p$-cyclotomic cosets modulo $n$ ($p$-cosets),

$$C_s = \{s, sp, \dots, sp^{r - 1}\} \mod (n - 1),$$

such that

$$\mathrm{min}_\alpha(x) = \mathrm{min}_{\beta^d}(x) = \prod_{j \in C_d} (x - \beta^j).$$

If the minimal polynomial contained any other roots, it would also need to contain all of its $p$ powers and we could separate all of these new terms into a polynomial which divides $\mathrm{min}_\alpha(x)$, contradicting the irreducibility of the minimal polynomial. Hence, over $\mathbb{F}_p$, we have

$$x^n - 1 = \prod_{d \mid n} \mathrm{min}_{\alpha^d}(x).$$

The assumption $\mathrm{gcd}(n, p) = 1$ ensures there are no repeated roots in the factorization since $f(x^p) = f(x)^p$.

Let $\mathcal{C} = (g(x))$ be a cyclic code. Then $\mathcal{T} = \cup_s C_s$, where $C_s$ are the $p$-cosets present in the construction of $g(x)$, is called the defining set of $\mathcal{C}$. As is clear from the definition, $\mathcal{T}$ completely determines $g(x)$ and vice versa, $g(x) = \prod_{j \in \mathcal{T}} (x - \alpha^j)$. The set of powers of $\alpha$ that are roots of $g(x)$ is called the variety (zero set) of $\mathcal{C}$, $\{\alpha^j \, : \, j \in \mathcal{T} \}$, and elements of the set are called zeros of the code. For two cyclic codes $\mathcal{C}_1$ and $\mathcal{C}_2$ with defining sets $\mathcal{T}_{\mathcal{C}_1}$ and $\mathcal{T}_{\mathcal{C}_2}$, respectfully, $\mathcal{C}_1 \subseteq \mathcal{C}_2$ if and only if $\mathcal{T}_{\mathcal{C}_2} \subseteq \mathcal{T}_{\mathcal{C}_1}$ ($g_2(x) \mid g_1(x)$).

BCH codes are a way of constructing a cyclic code with high minimum distance and high dimension by choosing $\mathcal{T}$ as small as possible that is a union of cyclotomic cosets with $\delta - 1$ consecutive elements. A BCH code $\mathcal{C}$ over $\mathbb{F}_p$ of length $n$ and designed distance $2 \leq \delta < n$ is a cyclic code with defining set $\mathcal{T} = C_b \cup C_{b + 1} \cup \dots \cup C_{b + \delta - 2}$ and zeros generated by a primitive element $\alpha \in \mathbb{F}_{p^m}$, where $m = \mathrm{ord}_n(p)$. A BCH code is called narrow-sense if $b = 1$ and primitive if $n = p^m - 1$. The number $b$ is called the offset of the code. It is crucial to note that many sources define narrow-sense as $b = 0$. This definition uses the zero set $\{ \alpha^{j + b} \}$ which is a shifted version of the definition above. While less common, the emphasis on the defining set over the zero set makes the choice $b = 1$ more natural for this work. The ``Magma`` and ``Sagemath`` coding theory libraries define narrow-sense as $b = 1$ and default to this parameter, although previous versions of the latter used $b = 0$.

The BCH bound says that if the defining set of a cyclic code contains a set of $\delta - 1$ consecutive integers (modulo $n$), then the minimum distance of the code is at least $\delta$. BCH codes therefore have minimum distance at least $\delta$ by design and maximize dimension by containing no extra roots. The dual of a BCH is, in general, not a BCH code, as the remaining cyclotomic cosets giving $h(x)$ need no longer be consecutive; however, the dual of narrow-sense BCH codes are BCH codes.

Reed-Solomon codes are primitive BCH codes over $\mathbb{F}_{p^m}$ for an integer $m \geq 1$. In this case, $\mathbb{F}_{p^m}$ is the splitting field of $x^{p^m - 1} - 1$ and each element $\alpha_i$ has minimal polynomial $x - \alpha_i$ with cyclotomic cosets of cardinality one. Hence, BCH codes are related to two fields while Reed-Solomon codes are only related to one. Reed-Solomon codes have the theoretically maximum possible distance with parameters $[n, k, n - k + 1]$. If $\mathcal{C}$ is an $[n, k, d]_{p^m}$ Reed-Solomon code, then $\mathcal{C}|_{\mathbb{F}_p}$ is the BCH code over $\mathbb{F}_p$ of length $n$ and designed distance $d$. The proof of this follows immediately from the fact that the codewords of the BCH code are elements of $\mathbb{F}_p^n$ and the zero set of the Reed-Solomon code is a subset of the zero set of the BCH code. Unlike BCH codes which can have any length relatively prime with the characteristic of the field, Reed-Solomon codes over $\mathbb{F}_p$ have $n \leq p$ and therefore do not make good binary codes directly. Instead, one may construct a ``binary Reed-Solomon code" using the expansion procedure for $\mathbb{F}_{2^m}/\mathbb{F}_2$.










# Remove?
We will only be concerned with cyclic Reed-Solomon codes, but the more general, and original, definition of Reed-Solomon codes will lead us into the final family of codes we will use in this work. Let $\mathcal{P}_k(x)$ denote the set of polynomials of degree less than $k$ in $\mathbb{F}_{p^m}[x]$. The Reed-Solomon code of length $n \leq p^m$ and dimension $k < n$ is given by

$$\mathrm{RS}_{p^m}(n, k) = \{ (f(\alpha_1), \dots, f(\alpha_n)) \mid f(x) \in \mathcal{P}_k(x)\},$$

where $\alpha_i \in \mathbb{F}_{p^m}$. The most common case $n = p^m$ is the extended code of the cyclic definition, but only the case $n = p^m -1$ is, in general, cyclic. The proof of this is direct application of the Chinese Remainder Theorem.