# [New Codes From Old](@id new-codes-from-old-api)

## Combining two or more codes

`u_u_plus_v` and `u_plus_w_v_plus_w_u_plus_v_plus_w` throw an error when one of
the arguments is the zero code. For the latter, let `C1` be an
``[n, k_1, d_1]`` and `C2` an ``[n, k_2, d_2]`` linear code; the construction
produces an ``[3n, 2k_1 + k_2]`` linear code, and for binary codes

```math
\mathrm{wt}(u + w \mid v + w \mid u + v + w)
    = 2\,\mathrm{wt}(u \veebar v) - \mathrm{wt}(w) + 4s,
    \qquad s = |\{i \mid u_i = v_i = 0,\, w_i = 1\}|.
```

Construction X takes an ``[n, k, d]`` code `C1`, an ``[n, k - l, d + e]`` code
`C2`, and an ``[m, l, e]`` code `C3` with `C2` a proper subcode of `C1`, and
returns an ``[n + m, k, d + e]`` code. Construction X3 takes
``[n, k_1, d_1]``, ``[n, k_2, d_2]``, ``[n, k_3, d_3]``,
``[n_4, k_2 - k_1, d_4]``, and ``[n_5, k_3 - k_2, d_5]`` codes with
`C1 ⊂ C2 ⊂ C3` and returns an ``[n + n_4 + n_5, k_3, d]`` code with
``d \geq \min\{d_1, d_2 + d_4, d_3 + d_5\}``.

The direct sum `⊕` has generator matrix `G1 ⊕ G2` and parity-check matrix
`H1 ⊕ H2`. The generator matrix of the direct product `×` is the Kronecker
product of the input generator matrices, and the parity-check matrix of the
tensor product is the Kronecker product of the input parity-check matrices.
There is some debate over how to define the entrywise (Schur) product of two
codes; the result is known to often be the full ambient space.

`juxtaposition` is representation dependent and therefore works on the
potentially overcomplete generator matrices rather than on the standard form.

## Modifying a single code

Extending adds a column to the generator matrix whose values make the row sums
zero. This even extension is the default for `extend(C)`, and the new column
may instead be inserted at any index with `extend(C, c)`. In the general case
one supplies a vector `a` and the new entries are the negated inner products of
`a` with the rows; the standard definition is the special case where `a` is the
all-ones vector.

Puncturing deletes columns from the generator matrix and removes any resulting
zero rows. Expurgating deletes rows and removes any resulting zero columns,
working directly on the potentially overcomplete generator matrix rather than
the standard form. Shortening is expurgating followed by puncturing; the
implementation uses the theorem that the code shortened on `L` is the dual of
the dual punctured on `L`, that is `dual(puncture(dual(C), L))`.

Augmentation vertically joins a matrix to the bottom of the generator matrix,
again working on the potentially overcomplete form. Lengthening augments the
all-ones row and then extends.

`subcode_of_dimension_between_codes` adds generators of `C1 / C2` to `C2` until
the desired dimension is reached. The subfield subcode is computed directly via
an expansion, whereas the trace code is computed using Delsarte's theorem.

If `C` is a quasi-cyclic code, `permute_code` returns a `LinearCode`.

!!! warning "Experimental"
    `even_subcode`, `doubly_even_subcode`, and `triply_even_subcode` need
    significantly more testing, but appear to work so far.

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/new_codes_from_old.jl"]
Private = false
```

## Weight reduction

Weight reduction rewrites a parity-check matrix so that its row and column
weights meet a target, at the cost of extra rows and columns. See the
[weight reduction tutorial](@ref weight-reduction-tutorial) for a detailed explanation.

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/weight_reduction.jl"]
Private = false
```
