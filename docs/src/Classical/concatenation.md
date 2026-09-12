# Concatenated Codes

## Background

There are at least three different meanings for the term "code concatenation":

1. the concatenation of a code over a finite field ``E`` with another code over
   a subfield ``F < E``,
2. the concatenation of two codes over the same field, and
3. the generalized concatenation scheme of Blokh and Zyablov.

In the original proposal there is an ``[n_o, k_o, d_o]_{q_o}`` *outer code*
``\mathrm{C}_{\mathrm{out}}`` and an ``[n_i, k_i, d_i]_{q_i}`` *inner code*
``\mathrm{C}_{\mathrm{in}}`` with ``\mathbb{F}_{q_i} < \mathbb{F}_{q_o}`` and
``k_i = [\mathbb{F}_{q_o} : \mathbb{F}_{q_i}]``. Each symbol of the outer code
is expanded to the subfield of the inner code. The inner code then encodes
``k_i`` symbols of the result at a time, and the results are concatenated into
a single vector. Since the dimension of the inner code is the degree of the
field extension, each symbol is expanded into ``k_i`` symbols and the inner
code encodes each symbol of the outer code individually.

The second case is a slight generalization. Both codes are over the same field
and the dimension of the inner code must divide the length of the outer code.
As before, the input is first encoded with the outer code, and then the inner
code encodes ``k_i`` symbols at a time, concatenating the results. Under the
first case ``k_i`` would have to be one.

The constructor examines the input codes and selects the correct procedure
automatically. It also accepts an outer code over an extension field, in which
case the outer code is expanded first and the second method is applied.

The concatenation type is `:expanded`, `:same`, or `:generalized` according to
which of the three methods was used. When the concatenation required expansion,
the basis and dual basis used are available; otherwise they are `missing`.

`encode` accepts valid inputs to both the full concatenated code and the outer
code. In the latter case it performs the two-step encoding described above.

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/concatenation.jl"]
Private = false
```

## References

È. L. Blokh and V. V. Zyablov, "Coding of Generalized Concatenated Codes",
Probl. Peredachi Inf., 10:3 (1974), 45-50; Problems Inform. Transmission,
10:3 (1974), 218-222.
