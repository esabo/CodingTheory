# [Cycles and ACE](@id ldpc-cycles-api)

Message passing is exact on a tree, so the cycles of the Tanner graph are what
limit an LDPC code's iterative decoding performance. This page collects the
tools for measuring them.

Two notions of cycle appear here and they are not interchangeable. A *short*
cycle is one of length at most twice the girth, enumerated per variable node;
these are the cycles that matter for decoding and can be counted cheaply. A
*simple* cycle is any cycle without repeated vertices, and enumerating all of
them is exponential in general, so those functions are for small graphs only.

The approximate cycle extrinsic message degree (ACE) of a cycle is the number
of edges leaving it, counting the degree of each variable node on the cycle
minus two. A cycle with a low ACE value is nearly isolated from the rest of the
graph, so extrinsic information reaches it slowly, which is what makes low-ACE
short cycles the usual culprits in error floors. The ACE spectrum records the
smallest ACE value found at each cycle length, and the distribution functions
summarize per-variable-node ACE values by their mean, median, or mode.

`remove_cycles` greedily edits the parity-check matrix to eliminate cycles, and
the plotting functions require a Makie backend to be loaded.

```@autodocs
Modules = [CodingTheory]
Pages = ["LDPC/cycles.jl"]
Private = false
```
