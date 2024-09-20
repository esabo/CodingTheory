# Modifying Codes

## Constructors
The first two constructors throw an error when one of the arguments is the zero code.
```@docs
u_u_plus_v
```

Let `C1` be an $[n, k1, d1]$ and `C2` an $[n, k2, d2]$ linear code. This construction produces an $[3n, 2k1 + k2]$ linear code. For binary codes, $\mathrm{wt}(u + w \mid v + w \mid u + v + w) = 2 \mathrm{wt}(u \veebar v) - \mathrm{wt}(w) + 4s$, where $s = |\{i \mid u_i = v_i = 0, w_i = 1\}|$.
```@docs
u_plus_w_v_plus_w_u_plus_v_plus_w
```

Let `C1` be an $[n, k, d]$, `C2` an $[n, k - l, d + e]$, and `C3` an $[m, l, e]$ linear code with `C2 ⊂ C1` be proper. Construction X creates a $[n + m, k, d + e]$ linear code.
```@docs
construction_X
```

Let `C1` be an $[n, k1, d1]$, `C2` an $[n, k2, d2]$, `C3` an $[n, k3, d3]$, `C4` an $[n4, k2 - k1, d4]$, and `C5` an $[n5, k3 - k2, d5]$ linear code with `C1 ⊂ C2 ⊂ C3`. Construction X3 creates an $[n + n4 + n5, k3, d]$ linear code with $d \geq \min\{d1, d2 + d4, d3 + d5\}$.
```@docs
construction_X3
```

The direct sum code has generator matrix `G1 ⊕ G2` and parity-check matrix `H1 ⊕ H2`.

```@docs
⊕
```

The generator matrix of the (direct) product code is the kronecker product of the generator matrices of the inputs.

```@docs
×
```

The parity-check matrix of the tensor product code is the kronecker product of the parity-check matrices of the inputs.

There is some debate on how to define this product. This is known to often be the full ambient space.
```@docs
entrywise_product_code
```

```@docs
/
```

`juxtaposition` is representation dependent and therefore works on the potentially over-complete generator matrices, not on the standard form.
```@docs
juxtaposition
```

## Methods

If `C` is a quasi-cyclic code, `permute_code` returns a `LinearCode` object.

The most common way to extend a code is to add an extra column to the generator matrix whose values make the sum of the rows zero. This is called an even extension and is the default for `extend(C)`. Alternatively, this new column may be inserted at any index `c` in the matrix, e.g. `extend(C, c)`. In the most general case, one may provide a vector `a` and define the values of the new column to be `-a` dot the row. The standard definition is clearly just the special case that `a` is the all-ones vector.
```@docs
extend
```

Puncturing deletes columns from the generator matrix and then removes any potentially resulting zero rows.
```@docs
puncture
```

Expurgating deletes rows from the generator matrix and then removes any potentially resulting zero columns. This function works directly on the potentially over-complete generator matrix and not on the standard form.
```@docs
expurgate
```

Shortening is expurgating followed by puncturing. This implementation uses the theorem that the dual of code shortened on `L` is equal to the puncture of the dual code on `L`, i.e., `dual(puncture(dual(C), L))`.
```@docs
shorten
```

Augmentation vertically joins the matrix `M` to the bottom of the generator matrix of `C`. This function works directly on the potentially over-complete generator matrix and not on the standard form.
```@docs
augment
```

Lengthening augments the all 1's row and then extends.
```@docs
lengthen
```

```@docs
subcode
```

This function arguments generators of `C1 / C2` to  `C2` until the desired dimenion is reached.
```@docs
subcode_of_dimension_between_codes
```

```@docs
expanded_code
```

The subfield subcode is computed directly via an expansion, whereas the trace code is computed using Delsarte's theorem.
```@docs
subfield_subcode
```

```@docs
trace_code
```

!!! warning "Experimental"
    The next two functions need significantly more testing, but appear to work so far.

```@docs
even_subcode
```

```@docs
doubly_even_subcode
```

## Weight reduction

See the weight reduction tutorial for a more detailed explanation of this function.

```@docs
weight_reduction
```
