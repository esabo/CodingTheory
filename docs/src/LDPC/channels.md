# [LDPC Noise Channels](@id ldpc-channels-api)

The channel types record the parameters of the standard binary-input channels
and provide their capacities and log-likelihood ratios. `transmit` draws a
channel realization, which is the input the
[soft-decision decoders](@ref ldpc-decoders-api) expect.

The LDPC ensemble types also live in this file group: an `LDPCEnsemble` is
described by its degree distributions rather than by a specific matrix, and a
`METEnsemble` is the multi-edge-type generalization, where edges are sorted
into classes so that different sockets of a variable node can attach to
different classes of check node. The analysis routines that consume them are in
[Analysis](@ref ldpc-analysis-api).

```@autodocs
Modules = [CodingTheory]
Pages = ["LDPC/channels.jl", "LDPC/types.jl"]
Private = false
```
