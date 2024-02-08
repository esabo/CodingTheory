# Concatenated Codes
## Background
There are at least three different meanings for the term "code concatenation":
1. The concatenation of a code over a finite field $E$ with another code over a subfield $F < E$.
2. The concatenation of two codes over the same field.
3. The generalized concatenation scheme of Blokh and Zyablov.
In the original proposal, there is an $[n_o, k_o, d_o]_{q_o}$ *outer code* $\mathrm{C}_{\mathrm{out}}$ and an $[n_i, k_i, d_i]_{q_i}$ *inner code* $\mathrm{C}_{\mathrm{in}}$ such that $\mathbb{F}_{q_i} < \mathbb{F}_{q_o}$ and $k_i = [\mathbb{F}_{q_o} : \mathbb{F}_{q_i}]$. Each bit of the outer code is expanded to the subfield of the inner code (link). The inner code then encodes $k_i$ bits of the result at a time and the result is concatenated into a single vector. Since the dimension of the inner code is the degree of the field extension, each bit gets expanded into $k_i$ bits, and the inner code encodes each bit of the outer code individually.

The second case is a slight generalization of this. Here, both codes are over the same field and the dimension of the inner code must divide the length of the outer code. As before, the input is first encoded with the outer code and then the inner code encodes $k_i$ bits at a time, concatenating the results. (Using the first case, $k_i$ would have to be 1.)





## Constructors
The constructor for the first two cases examines the input codes and automatically selects the correct procedure. Additionally, skipping one extra step, the constructor will also except an outer code over an extension field for which the second method can then be applied to the inner code and an expanded version of the outer code.
```@docs
concatenate
```

## Attributes
```@docs
inner_code
```

```@docs
outer_code
```

The concatenation type is `:expanded`, `:same`, or `:generalized` depending on which of the three methods above is used. 
```@docs
concatenation_type
```
If the concatenation required expansion, the basis and dual basis used for the expansion are returned via the following commands; otherwise, these are `missing`.
```@docs
expansion_basis
```

```@docs
expansion_dual_basis
```

## Methods
This function accepts valid inputs to both the full concatenated code and the outer code. In the later case, it performs a two-step encoding as described above.
```@docs
encode
```

[È. L. Blokh, V. V. Zyablov, “Coding of Generalized Concatenated Codes”, Probl. Peredachi Inf., 10:3 (1974), 45–50; Problems Inform. Transmission, 10:3 (1974), 218–222]