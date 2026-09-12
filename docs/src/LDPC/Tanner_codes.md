# Generalized LDPC and Tanner Codes

A Tanner code generalizes an LDPC code by replacing each single parity check
with a short local code imposed on the edges incident to that check vertex.
Taking the local code to be a single parity check recovers an ordinary LDPC
code, while stronger local codes buy better distance at the same degree.

The spectral functions here are what make these codes tractable to reason
about: the Sipser-Spielman bound turns the spectral gap of the underlying graph
into a distance guarantee, so an expander graph yields a code with
provably good distance.

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/Tanner.jl"]
Private = false
```
