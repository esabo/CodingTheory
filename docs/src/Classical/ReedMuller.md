# Reed-Muller Codes

Reed-Muller codes are a subtype of `LinearCode` and inherit its methods.

The binary family is generated from the recursive ``(u \mid u + v)`` form of the
generator matrix. Sources differ on the base case: if `alt` is `true` the
identity is used as the generator matrix of ``\mathcal{RM}(1, 1)``, and
otherwise ``\begin{pmatrix} 1 & 1\\ 0 & 1\end{pmatrix}`` is used.

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/ReedMuller.jl"]
Private = false
```
