# Reed-Muller Codes

Reed-Muller codes are a subtype of `LinearCode` and inherit its methods.

## Constructors
The (binary) Reed-Muller family is generated using the recursive, $(u \mid u + v)$-form of the generator matrices. Different sources use different conventions for the base case generator matrix. If `alt` is `true`, the identity is used for the generator matrix for $\mathcal{RM}(1, 1)$; otherwise, $\begin{pmatrix} 1 & 1\\ 0 & 1\end{pmatrix}$ is used.

```@docs
ReedMullerCode
```

## Attributes

```@docs
order
```

```@docs
number_of_variables
```

## Methods
