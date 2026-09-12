# Quantum LDPC Parameters

A quantum LDPC code is a stabilizer code whose check weights and qubit degrees
stay bounded as the length grows. That is a property of a chosen set of
generators rather than of the code itself, so these functions all describe a
particular presentation: adding a redundant generator changes the degrees
without changing the code.

The accessors come in ``X``, ``Z``, and combined flavors for CSS codes,
mirroring the classical LDPC degree machinery: degree distributions and their
polynomials, maximum and minimum degrees, density, and regularity tests. For
subsystem codes the gauge generators have their own weight functions, since the
gauge group is what is actually measured.

`quantum_LDPC_parameters` collects the headline numbers in one call, which is
usually what you want when comparing families.

```@autodocs
Modules = [CodingTheory]
Pages = ["Quantum/QLDPC.jl"]
Private = false
```
